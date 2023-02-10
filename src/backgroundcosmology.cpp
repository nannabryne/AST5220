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

  //=============================================================================
  // TODO: Compute Omegar, Omeganu, OmegaLambda, H0, ...
  //=============================================================================
  //...
  //...
  //...
  //...
  //-> H0:
  H0 = h * Constants.H0_over_h;

  // Omeganu:
  Omeganu = 0.;

  // Omegak:
  Omegak = 0.;

  // Omegar:
  double kT = (Constants.k_b*TCMB);
  double h3 = Constants.hbar*Constants.hbar*Constants.hbar;
  double c5 = Constants.c*Constants.c*Constants.c*Constants.c*Constants.c;
  Omegar = (M_PI*M_PI)/15 * (kT*kT*kT*kT)/(h3*c5) * (8*M_PI*Constants.G)/(3*H0*H0);


  // OmegaLambda:
  OmegaLambda = 1 - Omegak - (Omegab + OmegaCDM) - (Omegar + Omeganu);
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
    
  //=============================================================================
  // TODO: Set the range of x and the number of points for the splines
  // For this Utils::linspace(x_start, x_end, npts) is useful
  //=============================================================================
  Vector x_array;

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    //=============================================================================
    // TODO: Set the rhs of the detadx ODE
    //=============================================================================
    //...
    //...

    detadx[0] = 0.0;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set the initial condition, set up the ODE system, solve and make
  // the spline eta_of_x_spline 
  //=============================================================================
  // ...
  // ...
  // ...
  // ...

  Utils::EndTiming("Eta");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...
  // double a = std::exp(x);
  // temporary

  double expr = 0; 
  expr += (Omegab * OmegaCDM)*std::exp(-3*x);
  expr += (Omegar * Omeganu)*std::exp(-4*x);
  expr += Omegak*std::exp(-2*x);
  expr += OmegaLambda;

  return H0 * std::sqrt(expr);
}

double BackgroundCosmology::Hp_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_Omegab(double x) const{ 
  if(x == 0.0) return Omegab;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_Omegar(double x) const{ 
  if(x == 0.0) return Omegar;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_Omeganu(double x) const{ 
  if(x == 0.0) return Omeganu;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_Omegak(double x) const{ 
  if(x == 0.0) return Omegak;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}
double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
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
  const double x_min = -10.0;
  const double x_max =  0.0;
  const int    n_pts =  100;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
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
