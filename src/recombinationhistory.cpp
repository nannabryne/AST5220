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

  // set up x-, Xe- and ne-arrays:
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

    // the Peebles ODE equation:
    ODESolver Peebles_Xe_ode;
    ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
      return rhs_Peebles_ode(x, Xe, dXedx);
    };

    // set up and solve Peebles ODE:
    Vector Xe_ini = {Xe_arr[idx_Pee]};  // initial condition from Saha regime
    Vector x_arr_rest = Vector(x_array.begin()+idx_Pee, x_array.end());
    Peebles_Xe_ode.solve(dXedx, x_arr_rest, Xe_ini);
    Vector Xe_arr_Peebles = Peebles_Xe_ode.get_data_by_component(0);
    
    // fill array:
    for(int i=idx_Pee; i<npts; i++){
      Xe_arr[i] = Xe_arr_Peebles[i-idx_Pee];
      ne_arr[i] = Xe_arr[i]*nb_of_x(x_array[i]);
    }

  }

  Vector log_Xe_arr = log(Xe_arr);
  // spline result:
  log_Xe_of_x_spline.create(x_array, log_Xe_arr, "Xe");

}

void RecombinationHistory::solve_for_optical_depth_tau(int nsteps){
  // Utils::StartTiming("tau");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation 

  int npts = nsteps+1;
  Vector x_array = Utils::linspace(x_start, 0, npts);
  Vector x_array_backward = Utils::linspace(0, x_start, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){
    // Set the derivative for photon optical depth
    dtaudx[0] = dtaudx_of_x(x);
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
    gt_arr[i] =  -exp(-tau_of_x_spline(xi))*dtaudx_of_x(xi);
  }
  gt_of_x_spline.create(x_array, gt_arr, "gt");

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
    
  // compute and spline Xe, ne:
  Utils::StartTiming("Xe");
  solve_number_density_electrons(nsteps_Xe);
  Utils::EndTiming("Xe");

  // compute and spline tau and g:
  Utils::StartTiming("tau");
  solve_for_optical_depth_tau(nsteps_tau);
  Utils::EndTiming("tau");

  // compute and spline rs:
  Utils::StartTiming("rs");
  solve_for_sound_horizon(nsteps_rs);
  Utils::EndTiming("rs");


  if(print_milestones){
    milestones(1e6);
  }

}

void RecombinationHistory::milestones(int nsteps){

  //  Locate milestones
  Vector x_array = Utils::linspace(-10, -4, nsteps+1);
  double gmax = 0., XXmin = 1., tauumin = 1.;
  int i = 0;
  double x_lss = 0., x_rec = 0.;  // last scattering surface, recombination
  double x_lss_alt = 0.;  // alternative definition
  
  
  double g_curr, XX_curr, tauu_curr, x_curr;
  for(i=0; i<nsteps+1; i++){
    x_curr = x_array[i];
    g_curr = gt_of_x(x_curr);
    tauu_curr = abs(tau_of_x(x_curr) - 1.);
    XX_curr = abs(Xe_of_x(x_curr) - 0.1);
    if(g_curr>gmax){
      gmax = g_curr;
      x_lss = x_curr;
    }
    if(tauu_curr<tauumin){
      tauumin = tauu_curr;
      x_lss_alt = x_curr;
    }

    if(XX_curr<XXmin){
      XXmin = XX_curr;
      x_rec = x_curr;
    }
  }

  double x_0 = 0.;

  // x_lss = x_lss_alt;


  //  Print table

  double conv = 1/Constants.s_per_Gyr*1e6;
  printf("______________________________________________________________________\n");
  printf("                            |      x           z             t       | \n");
  printf("----------------------------------------------------------------------\n");
  printf("recombination               |");
  printf(" %11.7f %11.3f %11.3f ka |\n", x_rec, cosmo->z_of_x(x_rec), cosmo->t_of_x(x_rec)*conv);
  printf("last scattering surface     |");
  printf(" %11.7f %11.3f %11.3f ka |\n", x_lss, cosmo->z_of_x(x_lss), cosmo->t_of_x(x_lss)*conv);
  printf("----------------------------------------------------------------------\n");
  printf("f-o abundance of free electrons today: %11.7e \n", Xe_of_x(x_0));
  printf("sound horizon at decoupling:           %11.3f Mpc\n", rs_of_x(x_lss)/Constants.Mpc);
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

  // fetch cosmological parameters:
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
  double Hp_a = cosmo->H_of_x(x);     // Hp(x)/a   (= H)

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
  double Cr = (Lambda_2s1s + Lambda_alpha) / (Lambda_2s1s + Lambda_alpha + beta2);
  

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
  return -Constants.c*ne_of_x(x)*Constants.sigma_T*exp(x)/cosmo->Hp_of_x(x);
}

double RecombinationHistory::ddtaudxx_of_x(double x) const{
  return tau_of_x_spline.deriv_xx(x);
}


//  ----------------------------------
//  GET VISIBILITY + DERIVATIVES
//  ----------------------------------

double RecombinationHistory::gt_of_x(double x) const{
  return gt_of_x_spline(x);
}

double RecombinationHistory::dgtdx_of_x(double x) const{
  return gt_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgtdxx_of_x(double x) const{
  return gt_of_x_spline.deriv_xx(x);
}



//  ----------------------------------
//  GET ELECTRON FRACTION + NUM. DENSITY
//  ----------------------------------


double RecombinationHistory::Xe_of_x(double x) const{
  return exp(log_Xe_of_x_spline(x));
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

double RecombinationHistory::cs_of_x(double x) const{
  double R = 0.75 * cosmo->get_Omegab()/cosmo->get_Omegagamma() * exp(x);
  double U = 1./(3.* (1 + R));
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
