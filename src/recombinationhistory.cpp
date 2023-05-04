#include"recombinationhistory.h"


//  ----------------------------------
//  CONSTRUCT
//  ----------------------------------
   
RecombinationHistory::RecombinationHistory(BackgroundCosmology *cosmo, double Y_P) :
  cosmo(cosmo),
  Y_P(Y_P){
  
  _H0H0 = cosmo->get_H0()*cosmo->get_H0();
}

void RecombinationHistory::set_Saha_limit(double Xe_Saha_lower_limit){
  Xe_Saha_limit = Xe_Saha_lower_limit;
}



//  ----------------------------------
//  SOLVE
//  ----------------------------------

void RecombinationHistory::solve_number_density_electrons(int nsteps){

  //  set up x-, Xe- and ne-arrays:
  int npts = nsteps+1;
  Vector x_array = Utils::linspace(x_start, x_end, npts);
  Vector Xe_arr(npts);
  Vector ne_arr(npts);


  //  Calculate recombination history
  
  //  SAHA REGIME:

  int i = 0;
  // bool Saha_regime = ( Xe_current > Xe_Saha_limit);
  double Xe_current = 1.;
  while( ( Xe_current > Xe_Saha_limit) && i<npts){
    auto Xe_ne_data = electron_fraction_from_Saha_equation(x_array[i]);

    Xe_arr[i] = Xe_ne_data.first;
    ne_arr[i] = Xe_ne_data.second;
    Xe_current = Xe_arr[i];

    i++;
  }

  //  PEEBLES REGIME:
  
  if(i<npts-1){
    int idx_Pee = i-2;

    //  the Peebles ODE equation:
    ODESolver Peebles_Xe_ode;
    ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
      return rhs_Peebles_ode(x, Xe, dXedx);
    };

    //  set up and solve Peebles ODE:
    Vector Xe_ini = {Xe_arr[idx_Pee]};  // initial condition from Saha regime
    Vector x_arr_rest = Vector(x_array.begin()+idx_Pee, x_array.end());
    Peebles_Xe_ode.solve(dXedx, x_arr_rest, Xe_ini);
    Vector Xe_arr_Peebles = Peebles_Xe_ode.get_data_by_component(0);
    
    //  fill array:
    for(int i=idx_Pee; i<npts; i++){
      Xe_arr[i] = Xe_arr_Peebles[i-idx_Pee];
      ne_arr[i] = Xe_arr[i]*nb_of_x(x_array[i]);
    }

  }

  Vector log_Xe_arr = log(Xe_arr);
  Vector Xe_arr_new = exp(log_Xe_arr);
  // spline result:
  Xe_of_x_spline.create(x_array, Xe_arr_new, "Xe");
  log_Xe_of_x_spline.create(x_array, log_Xe_arr, "logXe");
  

}

void RecombinationHistory::solve_for_optical_depth_tau(int nsteps){
  // Utils::StartTiming("tau");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation 

  int npts = nsteps+1;
  Vector x_array = Utils::linspace(x_start, 0, npts);
  Vector x_array_backward = x_array;// = Utils::linspace(0, x_start, npts);
  std::reverse(x_array_backward.begin(), x_array_backward.end());


  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){
    // Set the derivative for photon optical depth
    dtaudx[0] = -Constants.c*ne_of_x(x)*Constants.sigma_T*exp(x)/cosmo->Hp_of_x(x);
    return GSL_SUCCESS;
  };

  //  Set up + solve the ODE and spline
  Vector tau_0 = {0};  // initial condition
  ODESolver ode_backward; 
  ode_backward.set_accuracy(1e-5, 1e-8, 1e-8);
  ode_backward.solve(dtaudx, x_array_backward, tau_0);
  // Vector sol_backward = ode_backward.get_data_by_component(0);
  Vector tau_arr = ode_backward.get_data_by_component(0);
  std::reverse(tau_arr.begin(), tau_arr.end());

  Vector dtau_arr = ode_backward.get_derivative_data_by_component(0);
  std::reverse(dtau_arr.begin(), dtau_arr.end());

  tau_of_x_spline.create(x_array, tau_arr, "tau");
  dtaudx_of_x_spline.create(x_array, dtau_arr, "dtaudx");


  //  Compute gt  and spline
  Vector gt_arr(npts);
  Vector dgt_arr(npts);
  // gt_arr = -exp(-tau_arr)*dtau_arr;
  double xi, A;
  for(int i=0; i<npts; i++){
    xi = x_array[i];
    A = -exp(-tau_arr[i]);
    gt_arr[i] = A*dtau_arr[i];
    dgt_arr[i] = A*(ddtaudxx_of_x(xi) - dtau_arr[i]*dtau_arr[i]);
  }
  gt_of_x_spline.create(x_array, gt_arr, "gt");
  dgtdx_of_x_spline.create(x_array, dgt_arr, "dgtdx");

}

void RecombinationHistory::solve_for_sound_horizon(int nsteps){

  int npts = nsteps+1;
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  ODEFunction drsdx = [&](double x, const double *rs, double *drsdx){
    drsdx[0] = cs_of_x(x)/cosmo->Hp_of_x(x);
    return GSL_SUCCESS;
  };

  //  Set up + solve the ODE and spline
  Vector rs_0 = {cs_of_x(x_start)/cosmo->Hp_of_x(x_start)};  // initial condition
  ODESolver ode; 
  ode.solve(drsdx, x_array, rs_0);
  Vector rs_arr = ode.get_data_by_component(0);

  rs_of_x_spline.create(x_array, rs_arr, "rs");

}

void RecombinationHistory::solve(bool print_milestones, int nsteps_Xe, int nsteps_tau, int nsteps_rs){
    
  //  compute and spline Xe, ne:
  Utils::StartTiming("Xe");
  solve_number_density_electrons(nsteps_Xe);
  Utils::EndTiming("Xe");

  //  compute and spline tau and g:
  Utils::StartTiming("tau");
  solve_for_optical_depth_tau(nsteps_tau);
  Utils::EndTiming("tau");

  //  compute and spline rs:
  Utils::StartTiming("rs");
  solve_for_sound_horizon(nsteps_rs);
  Utils::EndTiming("rs");


  if(print_milestones){
    milestones(1e6);
  }

}

void RecombinationHistory::milestones(int nsteps){

  //  Locate milestones

  std::pair<double,double>xrange(-10,-4);
  
  double xc = Utils::binary_search_for_value(tau_of_x_spline, 1., xrange);
  double xa = Utils::binary_search_for_value(Xe_of_x_spline, 0.1, xrange);

  std::pair<double,double>xrange_narrow(xc-0.15, xc+0.15);
  double xb = Utils::binary_search_for_value(dgtdx_of_x_spline, 0, xrange_narrow);

  double x_0 = 0.;


  //  Print table

  double conv = 1/Constants.s_per_Gyr*1e6;
  printf("______________________________________________________________________\n");
  printf("                            |      x           z             t       | \n");
  printf("----------------------------------------------------------------------\n");
  printf("tau = 1                     |");
  printf(" %11.7f %11.3f %11.3f ka |\n", xc, cosmo->z_of_x(xc), cosmo->t_of_x(xc)*conv);
  printf("X_e = 0.1                   |");
  printf(" %11.7f %11.3f %11.3f ka |\n", xa, cosmo->z_of_x(xa), cosmo->t_of_x(xa)*conv);
  printf("gt = max{gt}                |");
  printf(" %11.7f %11.3f %11.3f ka |\n", xb, cosmo->z_of_x(xb), cosmo->t_of_x(xb)*conv);
  printf("----------------------------------------------------------------------\n");
  printf("f-o abundance of free electrons today: %11.7e \n", Xe_of_x(x_0));
  printf("sound horizon at decoupling:           %11.3f Mpc\n", rs_of_x(xb)/Constants.Mpc);
  printf("----------------------------------------------------------------------\n");

  // printf("\nUsed %d+1 points in x-array from x=%.1f to x=%.1f.\n\n", nsteps, x_start, x_end);
}



//  ----------------------------------
//  EQUATIONS FOR SOLVING REC. 
//  ----------------------------------

std::pair<double,double> RecombinationHistory::electron_fraction_from_Saha_equation(double x) const{
 
  // physical constants:
  const double k_b    = Constants.k_b;
  const double m_e    = Constants.m_e;
  const double m_H    = Constants.m_H;
  const double eps_0  = Constants.epsilon_0;

  //  fetch cosmological parameters:
  const double Tb = cosmo->get_TCMB(x);  // Tb = Tgamma

  double Ups = eps_0/(k_b*Tb);  // Υ = ϵ_0/(k_b T_b)

  double Xe = 0.0;  // electron fraction
  double ne = 0.0;  // electron number density
  
  double U = 1/nb_of_x(x) * pow(m_e*Tb*k_b/(2*M_PI*_hbhb), 1.5) * exp(-Ups);
  if(U>1e7){ // avoid overflow
    Xe = 1.;
  }
  else{
    Xe = 0.5*U * (sqrt(1+4./U)-1);
  }
  ne = Xe * nb_of_x(x);

  return std::pair<double,double>(Xe, ne);
}




int RecombinationHistory::rhs_Peebles_ode(double x, const double *Xe, double *dXedx){

  const double X_e          = Xe[0];  // current value of X_e

  // physical constants in SI units:
  const double k_b          = Constants.k_b;
  const double c            = Constants.c;
  const double m_e          = Constants.m_e;
  const double hbar         = Constants.hbar;
  const double m_H          = Constants.m_H;
  const double sigma_T      = Constants.sigma_T;
  const double Lambda_2s1s  = Constants.Lambda_2s1s;
  const double eps_0        = Constants.epsilon_0;

  // fetch cosmological parameters:
  const double Tb = cosmo->get_TCMB(x);

  //  define 'helper' variables:
  double Ups = eps_0/(k_b*Tb);    // Υ = ϵ_0/(k_b T_b)
  // double Hp_a = cosmo->Hp_of_x(x)*exp(-x);     // Hp(x)/a   (= H)
  double Hp_a = cosmo->H_of_x(x);          // Hp(x)/a   (= H)

  //  compute ϕ_2:
  double phi2 = 0.448 * log(Ups);
  //  compute α^(2):
  double alpha2 = 8*pow(3*M_PI, -0.5) * c*sigma_T * sqrt(Ups) * phi2;
  //  compute β:
  double beta = alpha2 * pow((m_e*k_b*Tb / (2*M_PI * _hbhb)), 1.5) * exp(-Ups);
  //  compute β^(2):
  double beta2 = 0;
  if(Ups < 200) // avoid nan's
    beta2 = beta * exp(0.75*Ups);

  //  find n_i's:
  double n_H = (1 - Y_P) * nb_of_x(x);
  double n_1s = (1 - X_e) * n_H;

  //  compute Λ_α:
  double fac = (3*eps_0/(hbar*c));
  double Lambda_alpha = Hp_a * (fac*fac*fac) / (_8pi*_8pi) / n_1s;
  //  compute C_r:
  double Lambda = (Lambda_2s1s + Lambda_alpha);
  double Cr = Lambda / (Lambda + beta2);
  

  //  find rhs:
  dXedx[0] = Cr/Hp_a * (beta*(1-X_e) - n_H*alpha2*X_e*X_e);

  return GSL_SUCCESS;
}


//  ----------------------------------
//  GET OPTICAL DEPTH + DERIVATIVES
//  ----------------------------------

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{
  // return -Constants.c*ne_of_x(x)*Constants.sigma_T*exp(x)/cosmo->Hp_of_x(x);
  return dtaudx_of_x_spline(x);
}

double RecombinationHistory::ddtaudxx_of_x(double x) const{
  return dtaudx_of_x_spline.deriv_x(x);
}



//  ----------------------------------
//  GET VISIBILITY + DERIVATIVES
//  ----------------------------------

double RecombinationHistory::gt_of_x(double x) const{
  return gt_of_x_spline(x);
}

double RecombinationHistory::dgtdx_of_x(double x) const{
  return dgtdx_of_x_spline(x);
}

double RecombinationHistory::ddgtdxx_of_x(double x) const{
  return dgtdx_of_x_spline.deriv_x(x);
}



//  ----------------------------------
//  GET ELECTRON FRACTION + NUM. DENSITY
//  ----------------------------------


double RecombinationHistory::Xe_of_x(double x) const{
  // return exp(log_Xe_of_x_spline(x));
  return Xe_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{
  return nb_of_x(x)*Xe_of_x(x);
}



//  ----------------------------------
//  GET MISC. FUNCTIONS
//  ----------------------------------

double RecombinationHistory::nb_of_x(double x) const{
  const double Omegab0 = cosmo->get_Omegab();
  return (3*_H0H0 * Omegab0) / (_8piG * Constants.m_H * exp(3.*x));
}

double RecombinationHistory::R_of_x(double x) const{
  return 0.75 * cosmo->get_Omegab()/cosmo->get_Omegagamma() * exp(x);
}

double RecombinationHistory::cs_of_x(double x) const{
  double U = 1./(3.* (1 + R_of_x(x)));
  return Constants.c * sqrt(U);
}

double RecombinationHistory::rs_of_x(double x) const{
  return rs_of_x_spline(x);
}

double RecombinationHistory::get_Y_P() const{
  return Y_P;
}


//  ----------------------------------
//  HANDLE I/O
//  ----------------------------------


void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Y_P:               " << Y_P             << "\n";
  std::cout << "Saha thresh: Xe > " << Xe_Saha_limit  << "\n";
  std::cout << std::endl;
} 

void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(OUTPUT_PATH + filename.c_str());
  const int npts       = 5000;
  const double x_min   = -12;
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
