#include "backgroundcosmology.h"


//  ----------------------------------
//  CONSTRUCT
//  ----------------------------------

BackgroundCosmology::BackgroundCosmology(
    double h0, 
    double Omegab0, 
    double Omegac0, 
    double OmegaK0,
    double Neff, 
    double TCMB0) :
  h0(h0),
  Omegab0(Omegab0),
  Omegac0(Omegac0),
  OmegaK0(OmegaK0),
  Neff(Neff), 
  TCMB0(TCMB0){
 
  //  Compute the remaining cosmological parameters 

  H0 = h0 * Constants.H0_over_h;

  double kT = (Constants.k_b*TCMB0);
  double h3 = Constants.hbar*Constants.hbar*Constants.hbar;
  double c5 = Constants.c*Constants.c*Constants.c*Constants.c*Constants.c;
  Omegagamma0 = (M_PI*M_PI)/15. * (kT*kT*kT*kT)/(h3*c5) * (8.*M_PI*Constants.G)/(3.*H0*H0);

  Omeganu0 = Neff * 7./8. * std::pow(4./11., 4./3.) * Omegagamma0;

  OmegaM0 = (Omegac0 + Omegab0);
  OmegaR0 = (Omegagamma0 + Omeganu0);
  OmegaLambda0 = (1 - OmegaK0 - OmegaM0 - OmegaR0);

}


//  ----------------------------------
//  SOLVE
//  ----------------------------------


void BackgroundCosmology::solve_conformal_time(int nsteps){

  Vector x_array = Utils::linspace(x_start, x_end, nsteps+1);

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

  Vector x_array = Utils::linspace(x_start, x_end, nsteps+1);

  // the ODE for dt/dx:
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){
    dtdx[0] = 1./H_of_x(x);
    return GSL_SUCCESS;
  };

  Vector t_ini = {1./H_of_x(x_start) * 0.5};  // initial condition
  // t_ini = {std::exp(2*x_start)/(2*H0*std::sqrt(Omegagamma0))};
  ODESolver ode; 
  ode.solve(dtdx, x_array, t_ini);
  Vector sol = ode.get_data_by_component(0);

  t_of_x_spline.create(x_array, sol, "t");

}

void BackgroundCosmology::solve(bool print_milestones, int nsteps){

  Utils::StartTiming("eta");  
  solve_conformal_time(nsteps);
  Utils::EndTiming("eta");

  Utils::StartTiming("t");  
  solve_cosmic_time(nsteps);
  Utils::EndTiming("t");

  if(print_milestones){
    milestones(nsteps);
  }

}

void BackgroundCosmology::milestones(int nsteps){

  //  Locate milestones

  double x_eq = log(OmegaR0/OmegaM0);
  double x_L = 1./3.*log(OmegaM0/OmegaLambda0);
  // double x_acc = 1./3.*log(OmegaM0/OmegaLambda0*0.5);
  double x_acc = x_L - 1./3.*log(2);
  double x_0 = 0.;
  
  
  //  Print table

  double conv = 1/Constants.s_per_Gyr;
  printf("______________________________________________________________________\n");
  printf("                            |      x           z             t       | \n");
  printf("----------------------------------------------------------------------\n");
  printf("radiation-matter equality   |");
  printf(" %11.7f %11.3e %11.3e ka |\n", x_eq, z_of_x(x_eq), t_of_x(x_eq)*conv*1e6);
  printf("start of acceleration       |");
  printf(" %11.7f %11.7f %11.7f Ga |\n", x_acc, z_of_x(x_acc), t_of_x(x_acc)*conv);
  printf("matter-dark energy equality |");
  printf(" %11.7f %11.7f %11.6f Ga |\n", x_L, z_of_x(x_L), t_of_x(x_L)*conv);
  printf("today                       |");
  printf(" %11.7f %11.7f %11.6f Ga |\n", x_0, z_of_x(x_0), t_of_x(x_0)*conv);
  printf("----------------------------------------------------------------------\n");
  printf("conformal time today: %11.7f c Ga \n", eta_of_x(x_0)/Constants.c*conv);
  printf("----------------------------------------------------------------------\n");
  printf("inverse conformal time at rad.-matter equality: %11.7f Mpc^-1 \n", 1/eta_of_x(x_eq)*Constants.Mpc); // k_eq
  printf("----------------------------------------------------------------------\n");

  printf("\nUsed %d+1 points in x-array from x=%.1f to x=%.1f.\n\n", nsteps, x_start, x_end);
}


//  ----------------------------------
//  COMPUTE HUBBLE + DERIVATIVES
//  ----------------------------------

/*
Comments:
  - Hp_of_x(x) and H_of_x(x) should be computationally efficient! This is the reasoning behind the crazy way of writing the expressions. 
  - The purpose of Xi{m}(x) are to not have to write the expressions over and over again. Some computational cost from this, though only on dHpdx_of_x(x) and ddHpdxx_of_x, so we leave it be.
*/


double BackgroundCosmology::Xi0(double x) const{
  double e_x = exp(-x);   // e^(-x)
  double e_x2 = e_x*e_x;  // e^(-2x)
  return (OmegaM0*e_x + OmegaR0*e_x2 + OmegaK0 + OmegaLambda0/e_x2);
}

double BackgroundCosmology::Xi1(double x) const{
  double e_x = exp(-x);   // e^(-x)
  double e_x2 = e_x*e_x;  // e^(-2x)
  return - OmegaM0*e_x - 2 * OmegaR0*e_x2 + 2 * OmegaLambda0/e_x2 ;
}

double BackgroundCosmology::Xi2(double x) const{
  double e_x = exp(-x);   // e^(-x)
  double e_x2 = e_x*e_x;  // e^(-2x)
  return OmegaM0*e_x + 4 * OmegaR0*e_x2 + 4 * OmegaLambda0/e_x2 ;
}



double BackgroundCosmology::H_of_x(double x) const{
  double e_x = exp(-x);   // e^(-x)
  double e_x2 = e_x*e_x;  // e^(-2x)
  // return H0 * std::sqrt( OmegaM0*std::exp(-3.*x) + OmegaR0*std::exp(-4.*x) + OmegaK0*std::exp(-2.*x) + OmegaLambda0 );
  return H0 * sqrt( OmegaM0*e_x2*e_x + OmegaR0*e_x2*e_x2 + OmegaK0*e_x2 + OmegaLambda0 );
}

double BackgroundCosmology::Hp_of_x(double x) const{
  return H0 * sqrt(Xi0(x));
  // return H0 * std::sqrt( OmegaM0*std::exp(-x) + OmegaR0*std::exp(-2.*x) + OmegaK0 + OmegaLambda0*std::exp(2.*x));
}


double BackgroundCosmology::dHpdx_of_x(double x) const{

  return H0 * Xi1(x)/sqrt(Xi0(x))*0.5;
}

double BackgroundCosmology::ddHpdxx_of_x(double x) const{

  double A0=Xi0(x), A1=Xi1(x), A2=Xi2(x);
  return H0 * 1./(2.*std::sqrt(A0)) * ( A2 - 0.5*(A1*A1)/A0 );
}



//  ----------------------------------
//  GET DENSITY PARAMETERS
//  ----------------------------------



double BackgroundCosmology::get_Omegab(double x) const{ 
  if(x == 0.0) return Omegab0;

  double Hpx = Hp_of_x(x);
  return Omegab0 * H0*H0 * std::exp(-x) / (Hpx*Hpx);
}

double BackgroundCosmology::get_Omegagamma(double x) const{ 
  if(x == 0.0) return Omegagamma0;

  double Hpx = Hp_of_x(x);
  return Omegagamma0 * H0*H0 * std::exp(-2.*x) / (Hpx*Hpx);
}

double BackgroundCosmology::get_Omeganu(double x) const{ 
  if(x == 0.0) return Omeganu0;

  double Hpx = Hp_of_x(x);
  return Omeganu0 * H0*H0 * std::exp(-2.*x) / (Hpx*Hpx);
}

double BackgroundCosmology::get_Omegac(double x) const{ 
  if(x == 0.0) return Omegac0;

  double Hpx = Hp_of_x(x);
  return Omegac0 * H0*H0 * std::exp(-x) / (Hpx*Hpx);
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda0;

  double Hx = H_of_x(x);
  return OmegaLambda0 * H0*H0 / (Hx*Hx);
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK0;

  double Hpx = Hp_of_x(x);
  return OmegaK0 * H0*H0 / (Hpx*Hpx);
}

double BackgroundCosmology::get_OmegaM(double x) const{ 
  if(x == 0.0) return OmegaM0;

  double Hpx = Hp_of_x(x);
  return OmegaM0 * H0*H0 * std::exp(-x) / (Hpx*Hpx);
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return Omegagamma0 + Omeganu0;

  double Hpx = Hp_of_x(x);
  return OmegaR0 * H0*H0 * std::exp(-2.*x) / (Hpx*Hpx);
}




//  ----------------------------------
//  DISTANCE MEASURES
//  ----------------------------------



double BackgroundCosmology::r_of_Chi(double Chi) const{
  double tol = 1e-3;  // numerical tolerance: smaller values are consideres zero in this context
  double r = Chi;     // radius
  if(OmegaK0 < -tol){ // closed Universe (k = +1): OmegaK0 < 0 
    double term = std::sqrt(-OmegaK0) * H0*Chi/Constants.c;
    r *= std::sin(term)/term;
  }
  else if(OmegaK0 > +tol){  // open Universe (k = -1): OmegaK0 > 0 
    double term = std::sqrt(OmegaK0) * H0*Chi/Constants.c;
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



//  ----------------------------------
//  GET SPLINES
//  ----------------------------------



double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::t_of_x(double x) const{
  return t_of_x_spline(x);
}

double BackgroundCosmology::z_of_x(double x) const{
  return std::exp(-x) - 1;
}


//  ----------------------------------
//  GET MISC. PARAMETERS
//  ----------------------------------


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




//  ----------------------------------
//  HANDLE I/O
//  ----------------------------------


void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "Omegab0:      " << Omegab0      << "\n";
  std::cout << "Omegac0:      " << Omegac0      << "\n";
  std::cout << "OmegaLambda0: " << OmegaLambda0 << "\n";
  std::cout << "OmegaK0:      " << OmegaK0      << "\n";
  std::cout << "Omeganu0:     " << Omeganu0     << "\n";
  std::cout << "Omegagamma0:  " << Omegagamma0  << "\n";
  std::cout << "Neff:         " << Neff         << "\n";
  std::cout << "h0:           " << h0           << "\n";
  std::cout << "TCMB0:        " << TCMB0        << "\n";
  std::cout << std::endl;
} 


void BackgroundCosmology::output(const std::string filename) const{
  const int n_pts = 4001;
  double xa=x_start, xb=x_end;
  if(x_start<-20.){xa = -20.;}
  if(x_end>5.){xb = 5.;}
  Vector x_array = Utils::linspace(xa, xb, n_pts);

  std::ofstream fp(OUTPUT_PATH + filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";  // 0
    fp << eta_of_x(x)        << " ";  // 1
    fp << t_of_x(x)          << " ";  // 2
    fp << Hp_of_x(x)         << " ";  // 3
    fp << dHpdx_of_x(x)      << " ";  // 4
    fp << ddHpdxx_of_x(x)    << " ";  // 5
    fp << get_Omegab(x)      << " ";  // 6
    fp << get_Omegac(x)      << " ";  // 7
    fp << get_OmegaLambda(x) << " ";  // 8
    fp << get_Omegagamma(x)  << " ";  // 9
    fp << get_Omeganu(x)     << " ";  // 10
    fp << get_OmegaK(x)      << " ";  // 11
    fp << dL_of_x(x)         << " ";  // 12
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);

}
