#include "backgroundcosmology.h"

//====================================================
// Constructors
//====================================================

BackgroundCosmology::BackgroundCosmology(
    double h, 
    double Omegab, 
    double OmegaCDM, 
    double Omegak,
    double Neff, 
    double TCMB) :
  h(h),
  Omegab(Omegab),
  OmegaCDM(OmegaCDM),
  Omegak(Omegak),
  Neff(Neff), 
  TCMB(TCMB)
{


  /*...................................................
  Compute the remaining cosmological parameters  
  ...................................................*/

  H0 = h * Constants.H0_over_h;

  Omeganu = 0.;

  Omegak = 0.;

  double kT = (Constants.k_b*TCMB);
  double h3 = Constants.hbar*Constants.hbar*Constants.hbar;
  double c5 = Constants.c*Constants.c*Constants.c*Constants.c*Constants.c;
  Omegar = (M_PI*M_PI)/15 * (kT*kT*kT*kT)/(h3*c5) * (8*M_PI*Constants.G)/(3*H0*H0);

  OmegaLambda = 1 - Omegak - (Omegab + OmegaCDM) - (Omegar + Omeganu);
}



//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){

  Vector x_array = Utils::linspace(x_start, x_end, 1e4);
  Utils::StartTiming("eta");  

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){
    detadx[0] = Constants.c/Hp_of_x(x);
    return GSL_SUCCESS;
  };

  Vector eta_ini = {Constants.c/Hp_of_x(x_start)};
  // eta_ini = {Constants.c/(H0*std::sqrt(Omegar)) * std::exp(x_start)};
  ODESolver ode; 
  ode.solve(detadx, x_array, eta_ini);
  Vector sol = ode.get_data_by_component(0);

  eta_of_x_spline.create(x_array, sol, "eta");
  Utils::EndTiming("eta");

}

void BackgroundCosmology::solve_time(){
  Utils::StartTiming("t");  

  // The ODE for dt/dx
  Vector x_array = Utils::linspace(x_start, x_end, 1e4);
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){
    dtdx[0] = 1/H_of_x(x);
    return GSL_SUCCESS;
  };

  Vector t_ini = {1/H_of_x(x_start)*0.5};
  // t_ini = {std::exp(2*x_start)/(2*H0*std::sqrt(Omegar))};
  ODESolver ode; 
  ode.solve(dtdx, x_array, t_ini);
  Vector sol = ode.get_data_by_component(0);

  t_of_x_spline.create(x_array, sol, "t");
  Utils::EndTiming("t");

  
}

//====================================================
// Get methods
//====================================================


double BackgroundCosmology::H_of_x(double x) const{

  double expr = 0; 
  expr += (Omegab + OmegaCDM)*std::exp(-3*x);
  expr += (Omegar + Omeganu)*std::exp(-4*x);
  expr += Omegak*std::exp(-2*x);
  expr += OmegaLambda;

  return H0 * std::sqrt(expr);
}

double BackgroundCosmology::Hp_of_x(double x) const{

  double expr = 0; 
  expr += (Omegab + OmegaCDM)*std::exp(-x);
  expr += (Omegar + Omeganu)*std::exp(-2*x);
  expr += Omegak;
  expr += OmegaLambda*std::exp(2*x);

  return H0 * std::sqrt(expr);
}

double BackgroundCosmology::dHpdx_of_x(double x) const{

  double t1, t2, t3, t4, M0, M1, M2; 
  t1 = (Omegab + OmegaCDM)*std::exp(-x);
  t2 = (Omegar + Omeganu)*std::exp(-2*x);
  t3 = Omegak;
  t4 = OmegaLambda*std::exp(2*x);

  M0 = t1+t2+t3+t4;
  M1 = -t1-2*t2+2*t4;
  // M2 = t1+4*t2+4*t4;

  return H0 * M1/std::sqrt(M0)*0.5;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{

  double t1, t2, t3, t4, M0, M1, M2; 
  t1 = (Omegab + OmegaCDM)*std::exp(-x);
  t2 = (Omegar + Omeganu)*std::exp(-2*x);
  t3 = Omegak;
  t4 = OmegaLambda*std::exp(2*x);

  M0 = t1+t2+t3+t4;
  M1 = -t1-2*t2+2*t4;
  M2 = t1+4*t2+4*t4;

  return H0 * 1/(2*std::sqrt(M0)) * ( M2 - 0.5*(M1*M1)/M0 );
}

double BackgroundCosmology::get_Omegab(double x) const{ 
  if(x == 0.0) return Omegab;

  double Hx = H_of_x(x);
  return Omegab * H0*H0 * std::exp(-3*x) / (Hx*Hx);
}

double BackgroundCosmology::get_Omegar(double x) const{ 
  if(x == 0.0) return Omegar;

  // double Hx = H_of_x(x);
  double Hp_x = Hp_of_x(x);
  // return Omegar * H0*H0 * std::exp(-4*x) / (Hx*Hx);
  return Omegar * H0*H0* std::exp(-2*x) / (Hp_x*Hp_x);
}

double BackgroundCosmology::get_Omeganu(double x) const{ 
  if(x == 0.0) return Omeganu;

  double Hx = H_of_x(x);
  return Omeganu * H0*H0 * std::exp(-4*x) / (Hx*Hx);
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;

  double Hx = H_of_x(x);
  return OmegaCDM * H0*H0 * std::exp(-3*x) / (Hx*Hx);
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;

  double Hx = H_of_x(x);
  return OmegaLambda * H0*H0 / (Hx*Hx);
}

double BackgroundCosmology::get_Omegak(double x) const{ 
  if(x == 0.0) return Omegak;

  double Hx = H_of_x(x);
  return Omegak * H0*H0 * std::exp(-2*x) / (Hx*Hx);
}

double BackgroundCosmology::r_of_Chi(double Chi) const{
  double tol = 1e-12;
  double term = std::sqrt(std::abs(Omegak)) * H0*Chi/Constants.c;
  double r;
  if(Omegak < tol){
    r = Chi * std::sin(term)/term;
  }
  else if(Omegak > tol){
    r = std::sinh(term)/term;
  }
  else{
    r = Chi;
  }
  return r;
}

double BackgroundCosmology::Chi_of_x(double x) const{
  return eta_of_x(0)-eta_of_x(x);
}

double BackgroundCosmology::dA_of_x(double x) const{
  // double a = Hp_of_x(x) / H_of_x(x);
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
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "Omegab:      " << Omegab      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "Omegak:      " << Omegak      << "\n";
  std::cout << "Omeganu:     " << Omeganu     << "\n";
  std::cout << "Omegar:      " << Omegar      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -14.0;
  const double x_max =  0.0;
  const int    n_pts =  400;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  // double age = t_of_x(0) / (365*24*60*60) *1e-9;
  // printf("Age of universe: %.5f Gyrs\n", age);

  std::ofstream fp(OUTPUT_PATH + filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << t_of_x(x)          << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << ddHpddx_of_x(x)    << " ";
    fp << get_Omegab(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_Omegar(x)      << " ";
    fp << get_Omeganu(x)     << " ";
    fp << get_Omegak(x)      << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);



}
