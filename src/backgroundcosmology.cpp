#include "backgroundcosmology.h"

//====================================================
// Constructors
//====================================================

BackgroundCosmology::BackgroundCosmology(
    double h0, 
    double Omegab0, 
    double OmegaCDM0, 
    double Omegak0,
    double Neff, 
    double TCMB0) :
  h0(h0),
  Omegab0(Omegab0),
  OmegaCDM0(OmegaCDM0),
  Omegak0(Omegak0),
  Neff(Neff), 
  TCMB0(TCMB0)
{
 
  
  //  Compute the remaining cosmological parameters 

  H0 = h0 * Constants.H0_over_h;

  double kT = (Constants.k_b*TCMB0);
  double h3 = Constants.hbar*Constants.hbar*Constants.hbar;
  double c5 = Constants.c*Constants.c*Constants.c*Constants.c*Constants.c;
  Omegagamma0 = (M_PI*M_PI)/15. * (kT*kT*kT*kT)/(h3*c5) * (8.*M_PI*Constants.G)/(3.*H0*H0);

  Omeganu0 = Neff * 7./8. * std::pow(4./11., 4./3.) * Omegagamma0;

  OmegaLambda0 = 1 - Omegak0 - (Omegab0 + OmegaCDM0) - (Omegagamma0 + Omeganu0);
}



void BackgroundCosmology::solve_conformal_time(int nsteps){

  Vector x_array = Utils::linspace(x_start, x_end, nsteps);

  // the ODE for deta/dx:
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){
    detadx[0] = Constants.c/Hp_of_x(x);
    return GSL_SUCCESS;
  };

  Vector eta_ini = {Constants.c/Hp_of_x(x_start)};  // initial condition
  // eta_ini = {Constants.c/(H0*std::sqrt(Omegagamma0)) * std::exp(x_start)};
  ODESolver ode; 
  ode.solve(detadx, x_array, eta_ini);
  Vector sol = ode.get_data_by_component(0);

  eta_of_x_spline.create(x_array, sol, "eta");

}

void BackgroundCosmology::solve_cosmic_time(int nsteps){

  Vector x_array = Utils::linspace(x_start, x_end, nsteps);

  // the ODE for dt/dx:
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){
    dtdx[0] = 1./H_of_x(x);
    return GSL_SUCCESS;
  };

  Vector t_ini = {1./H_of_x(x_start)*0.5};  // initial condition
  // t_ini = {std::exp(2*x_start)/(2*H0*std::sqrt(Omegagamma0))};
  ODESolver ode; 
  ode.solve(dtdx, x_array, t_ini);
  Vector sol = ode.get_data_by_component(0);

  t_of_x_spline.create(x_array, sol, "t");

}

void BackgroundCosmology::solve(int nsteps){

  Utils::StartTiming("eta");  
  solve_conformal_time(nsteps);
  Utils::EndTiming("eta");

  Utils::StartTiming("t");  
  solve_cosmic_time(nsteps);
  Utils::EndTiming("t");

}


//====================================================
// Get methods
//====================================================


double BackgroundCosmology::H_of_x(double x) const{

  double expr = 0; 
  expr += (Omegab0 + OmegaCDM0)*std::exp(-3.*x);
  expr += (Omegagamma0 + Omeganu0)*std::exp(-4.*x);
  expr += Omegak0*std::exp(-2.*x);
  expr += OmegaLambda0;

  return H0 * std::sqrt(expr);
}

double BackgroundCosmology::Hp_of_x(double x) const{

  double expr = 0; 
  expr += (Omegab0 + OmegaCDM0)*std::exp(-x);
  expr += (Omegagamma0 + Omeganu0)*std::exp(-2.*x);
  expr += Omegak0;
  expr += OmegaLambda0*std::exp(2.*x);

  return H0 * std::sqrt(expr);
}

double BackgroundCosmology::dHpdx_of_x(double x) const{

  double t1, t2, t3, t4, M0, M1, M2; 
  t1 = (Omegab0 + OmegaCDM0)*std::exp(-x);
  t2 = (Omegagamma0 + Omeganu0)*std::exp(-2.*x);
  t3 = Omegak0;
  t4 = OmegaLambda0*std::exp(2.*x);

  M0 = t1+t2+t3+t4;
  M1 = -t1-2*t2+2*t4;
  // M2 = t1+4*t2+4*t4;

  return H0 * M1/std::sqrt(M0)*0.5;
}

double BackgroundCosmology::ddHpdxx_of_x(double x) const{

  double t1, t2, t3, t4, M0, M1, M2; 
  t1 = (Omegab0 + OmegaCDM0)*std::exp(-x);
  t2 = (Omegagamma0 + Omeganu0)*std::exp(-2.*x);
  t3 = Omegak0;
  t4 = OmegaLambda0*std::exp(2.*x);

  M0 = t1+t2+t3+t4;
  M1 = -t1-2*t2+2*t4;
  M2 = t1+4*t2+4*t4;

  return H0 * 1./(2.*std::sqrt(M0)) * ( M2 - 0.5*(M1*M1)/M0 );
}

double BackgroundCosmology::get_Omegab(double x) const{ 
  if(x == 0.0) return Omegab0;

  double Hx = H_of_x(x);
  return Omegab0 * H0*H0 * std::exp(-3.*x) / (Hx*Hx);
}

double BackgroundCosmology::get_Omegagamma(double x) const{ 
  if(x == 0.0) return Omegagamma0;

  // double Hx = H_of_x(x);
  double Hp_x = Hp_of_x(x);
  // return Omegagamma0 * H0*H0 * std::exp(-4*x) / (Hx*Hx);
  return Omegagamma0 * H0*H0* std::exp(-2.*x) / (Hp_x*Hp_x);
}

double BackgroundCosmology::get_Omeganu(double x) const{ 
  if(x == 0.0) return Omeganu0;

  double Hx = H_of_x(x);
  return Omeganu0 * H0*H0 * std::exp(-4.*x) / (Hx*Hx);
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM0;

  double Hx = H_of_x(x);
  return OmegaCDM0 * H0*H0 * std::exp(-3.*x) / (Hx*Hx);
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda0;

  double Hx = H_of_x(x);
  return OmegaLambda0 * H0*H0 / (Hx*Hx);
}

double BackgroundCosmology::get_Omegak(double x) const{ 
  if(x == 0.0) return Omegak0;

  double Hx = H_of_x(x);
  return Omegak0 * H0*H0 * std::exp(-2.*x) / (Hx*Hx);
}

double BackgroundCosmology::get_OmegaM(double x) const{ 
  if(x == 0.0) return Omegab0 + OmegaCDM0;
  
  return get_Omegab(x) + get_OmegaCDM(x);
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return Omegagamma0 + Omeganu0;
  
  return get_Omegagamma(x) + get_Omeganu(x);
}

double BackgroundCosmology::r_of_Chi(double Chi) const{
  double tol = 1e-3;  // numerical tolerance: smaller values are consideres zero in this context
  double r = Chi;     // radius
  if(Omegak0 < -tol){ // closed Universe (k = +1): Omegak0 < 0 
    double term = std::sqrt(std::abs(Omegak0)) * H0*Chi/Constants.c;
    r *= std::sin(term)/term;
  }
  else if(Omegak0 > +tol){  // open Universe (k = -1): Omegak0 > 0 
    double term = std::sqrt(std::abs(Omegak0)) * H0*Chi/Constants.c;
    r *= std::sinh(term)/term;
  }
  return r;
}

double BackgroundCosmology::Chi_of_x(double x) const{
  return eta_of_x(0)-eta_of_x(x);
}

double BackgroundCosmology::dA_of_x(double x) const{
  double Chi = Chi_of_x(x);
  return std::exp(x)*r_of_Chi(Chi);
}

double BackgroundCosmology::dL_of_x(double x) const{
  double Chi = Chi_of_x(x);
  return r_of_Chi(Chi)/std::exp(x);
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  return dL_of_x(x);
}
double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  return Chi_of_x(x);
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::t_of_x(double x) const{
  return t_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h0; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB0;
  return TCMB0 * exp(-x); 
}



void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "Omegab0:      " << Omegab0      << "\n";
  std::cout << "OmegaCDM0:    " << OmegaCDM0    << "\n";
  std::cout << "OmegaLambda0: " << OmegaLambda0 << "\n";
  std::cout << "Omegak0:      " << Omegak0      << "\n";
  std::cout << "Omeganu0:     " << Omeganu0     << "\n";
  std::cout << "Omegagamma0:  " << Omegagamma0  << "\n";
  std::cout << "Neff:         " << Neff         << "\n";
  std::cout << "h0:           " << h0           << "\n";
  std::cout << "TCMB0:        " << TCMB0        << "\n";
  std::cout << std::endl;
} 


void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -14.0;
  const double x_max =  0.0;
  const int    n_pts =  400;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(OUTPUT_PATH + filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";  // 0
    fp << eta_of_x(x)        << " ";  // 1
    fp << t_of_x(x)          << " ";  // 2
    fp << Hp_of_x(x)         << " ";  // 3
    fp << dHpdx_of_x(x)      << " ";  // 4
    fp << ddHpdxx_of_x(x)    << " ";  // 5
    fp << get_Omegab(x)      << " ";  // 6
    fp << get_OmegaCDM(x)    << " ";  // 7
    fp << get_OmegaLambda(x) << " ";  // 8
    fp << get_Omegagamma(x)  << " ";  // 9
    fp << get_Omeganu(x)     << " ";  // 10
    fp << get_Omegak(x)      << " ";  // 11
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);

}
