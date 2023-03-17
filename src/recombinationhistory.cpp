#include"recombinationhistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(BackgroundCosmology *cosmo, double Yp) :
  cosmo(cosmo),
  Yp(Yp){
  
  _H0H0 = cosmo->get_H0()*cosmo->get_H0();
}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}



//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");

  //  set up x-, Xe- and ne-arrays:
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);


  // Calculate recombination history
  // bool saha_regime = true;
  int i = 0;
  // bool saha_regime = ( Xe_current > Xe_saha_limit);
  double Xe_current = 1.;

  while( ( Xe_current > Xe_saha_limit) && i<npts_rec_arrays){
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    Xe_arr[i] = Xe_ne_data.first;
    ne_arr[i] = Xe_ne_data.second;
    Xe_current = Xe_arr[i];

    i++;
  }

  int idx_pee = i-2;

  // The Peebles ODE equation
  ODESolver peebles_Xe_ode;
  ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
    return rhs_peebles_ode(x, Xe, dXedx);
  };

  //  set up and solve Peebles ODE:
  Vector Xe_ini = {Xe_arr[idx_pee]};  // initial condition from Saha regime
  Vector x_arr_rest = Vector(x_array.begin()+idx_pee, x_array.end());
  peebles_Xe_ode.solve(dXedx, x_arr_rest, Xe_ini);
  Vector Xe_arr_peebles = peebles_Xe_ode.get_data_by_component(0);
  
  //  fill array:
  for(int i=idx_pee; i<npts_rec_arrays; i++){
    Xe_arr[i] = Xe_arr_peebles[i-idx_pee];
    ne_arr[i] = Xe_arr[i]*nb_of_x(x_array[i]);
  }

  // spline result:
  Xe_of_x_spline.create(x_array, Xe_arr, "Xe");

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double m_e         = Constants.m_e;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;

  // Fetch cosmological parameters
  const double Tb = cosmo->get_TCMB(x);
  
  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;
  
  double f = 1/nb_of_x(x) * pow(m_e*Tb*k_b/(2*M_PI*_hbhb), 1.5) * exp(-epsilon_0/(k_b*Tb));
  if(f>1e7){ // avoid overflow
    Xe = 1.;
  }
  else{
    Xe = 0.5*f * (sqrt(1+4./f)-1);
  }
  ne = Xe * nb_of_x(x);

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  // const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  // const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double Lambda_2s1s = Constants.Lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;

  // Cosmological parameters
  const double Tb      = cosmo->get_TCMB(x);

  //  define helper variables:
  double eps_kT = epsilon_0/(k_b*Tb); // ϵ_0/(k_b T_b)
  double Hp_a = cosmo->H_of_x(x);     // Hp(x)/a   (= H)
  // double Hp_a = cosmo->Hp_of_x(x)/a;  // Hp(x)/a  

  //  compute ϕ_2:
  double phi2 = 0.448 * log(eps_kT);
  //  compute α^(2):
  double alpha2 = 8*pow(3*M_PI, -0.5) * c*sigma_T * sqrt(eps_kT) * phi2;
  //  compute β:
  double beta = alpha2 * pow((m_e*k_b*Tb / (2*M_PI * _hbhb)), 1.5) * exp(-eps_kT);
  //  compute β^(2):
  double beta2 = 0;
  if(eps_kT < 200) // avoid nan's
    beta2 = beta * exp(0.75*eps_kT);

  //  find n_i's:
  double n_H = (1 - Yp) * nb_of_x(x);
  double n_1s = (1 - X_e) * n_H;

  //  compute Λ_α:
  double fac = (3*epsilon_0/(hbar*c));
  double Lambda_alpha = Hp_a * (fac*fac*fac) / (_8pi*_8pi) / n_1s;
  //  compute C_r:
  double Cr = (Lambda_2s1s + Lambda_alpha) / (Lambda_2s1s + Lambda_alpha + beta2);
  

  //  find rhs:
  double expr = (beta*(1-X_e) - n_H*alpha2*X_e*X_e);
  dXedx[0] = Cr/Hp_a * expr;

  return GSL_SUCCESS;
}





//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("tau");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation 

  int npts = 8001;
  Vector x_array = Utils::linspace(x_start, 0, npts);
  Vector x_array_backward = Utils::linspace(0, x_start, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){
    // Set the derivative for photon optical depth
    dtaudx[0] = -Constants.c*ne_of_x(x)*Constants.sigma_T*exp(x)/cosmo->Hp_of_x(x);

    return GSL_SUCCESS;
  };

  //  Set up + solve the ODE and spline
  Vector tau_0 = {0};  // initial condition
  ODESolver ode_backward; 
  ode_backward.solve(dtaudx, x_array_backward, tau_0);
  Vector sol_backward = ode_backward.get_data_by_component(0);

  Vector tau_arr(npts);
  for(int i=0; i<npts; i++){
    tau_arr[i] = sol_backward[npts-(i+1)];
  }

  tau_of_x_spline.create(x_array, tau_arr, "tau");


  //  Compute gt + derivatives and spline
  Vector gt_arr(npts);
  double xi;
  for(int i=0; i<npts; i++){
    xi = x_array[i];
    gt_arr[i] = exp(-tau_of_x_spline(xi)) * Constants.c*ne_of_x(xi)*Constants.sigma_T*exp(xi)/cosmo->Hp_of_x(xi);
    // gt_arr[i] = - exp(tau_of_x_splmine(xi)) * tau_of_x_spline.deriv_x(xi);
  }
  gt_of_x_spline.create(x_array, gt_arr, "gt");

  Utils::EndTiming("tau");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{
  return tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtaudxx_of_x(double x) const{
  return tau_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::gt_of_x(double x) const{
  return gt_of_x_spline(x);
}

double RecombinationHistory::dgtdx_of_x(double x) const{
  return gt_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgtdxx_of_x(double x) const{
  return gt_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{
  // return exp(log_Xe_of_x_spline(x));
  return Xe_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{
  return nb_of_x(x)*Xe_of_x(x);
}

double RecombinationHistory::nb_of_x(double x) const{
  const double Omegab0 = cosmo->get_Omegab();
  return (3*_H0H0 * Omegab0) / (_8piG * Constants.m_H * exp(3.*x));
}





double RecombinationHistory::get_Yp() const{
  return Yp;
}


//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(OUTPUT_PATH + filename.c_str());
  const int npts       = 5000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtaudxx_of_x(x)     << " ";
    fp << gt_of_x(x)           << " ";
    fp << dgtdx_of_x(x)        << " ";
    fp << ddgtdxx_of_x(x)      << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}
