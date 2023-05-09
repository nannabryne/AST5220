#include"powerspectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================
  // Vector k_array;
  // Vector log_k_array = log(k_array);


  
  



  //=========================================================================
  // TODO: Make splines for j_ell. 
  // Implement generate_bessel_function_splines
  //=========================================================================
  Utils::StartTiming("besselspline");
  generate_bessel_function_splines();
  Utils::EndTiming("besselspline");

  //=========================================================================
  // TODO: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  double n_sampling = 32;
  double dk = 2.*M_PI/n_sampling * 1./cosmo->eta_of_x(0);
  int n_k = ( k_max - k_min ) / dk;

  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  // Vector k_array = Utils::linspace(k_min, k_max, n_k);
  Vector k_array = exp(log_k_array);

  Utils::StartTiming("lineofsight");
  line_of_sight_integration(k_array);
  Utils::EndTiming("lineofsight");

  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_Cell
  //=========================================================================

  Utils::StartTiming("C_ell");
  auto Cell_TT = solve_for_Cell(log_k_array, ThetaT_ell_of_k_spline, ThetaT_ell_of_k_spline);
  Cell_TT_spline.create(ells, Cell_TT, "Cell_TT_of_ell");
  Utils::EndTiming("C_ell");
  

}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  
  //  make storage for the splines:
  j_ell_splines = std::vector<Spline>(ells.size());

  //  Compute splines for bessel functions j_ell(z)

  const double eta0 = cosmo->eta_of_x(0);   // η(0) - conformal time today

  //  define z-grid:
  const double z_end = k_max * eta0;
  const double z_start = 0.;
  const double n_sampling = 25;
  const double dz = 2*M_PI/n_sampling;
  const int n_z = (int)((z_end - z_start)/dz);
  Vector z_array = Utils::linspace(z_start, z_end, n_z);


  for(size_t iell = 0; iell < ells.size(); iell++){
    const int ell = ells[iell];

    Vector j_ell_arr(n_z);
    for(int iz=0; iz<n_z; iz++){
      j_ell_arr[iz] = Utils::j_ell(ell, z_array[iz]);
    }

    //  make the j_ell_splines[i] spline:
    std::string str = "j_" + std::to_string(ell);
    j_ell_splines[iell].create(z_array, j_ell_arr, str);
  }

}


//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  

  //  make storage for the results:
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));
  
  //  define x-grid:
  const double x_start = -10.;
  const double x_end = 0.;
  const double n_sampling = 500;
  const double dx = 2.*M_PI/n_sampling;
  const int n_x = int( ( x_end - x_start ) / dx );

  const double eta0 = cosmo->eta_of_x(0); // η(0) - conformal time today

  //  For each k, solve the LOS integral for all the ell values

  #pragma omp parallel for
  for(size_t ik=0; ik<k_array.size(); ik++){

    double k = k_array[ik];
    for(int iell=0; iell<ells.size(); iell++){

      int ell = ells[iell];

      //  Solve using the trapezodial rule

      std::function<double(double)> F = [&](double x){
        double z = k * ( eta0 - cosmo->eta_of_x(x) );
        return source_function(x, k) * j_ell_splines[iell](z);
      }; 

      result[iell][ik] = integrate_trapezodial(F, x_start, x_end, dx);

    }

  }
  
  return result;
}




//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int nells      = ells.size();
  
  //  make storage for the splines we are to create:
  ThetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //  Solve for Theta_ell(k) and spline the result

  // Function returning the source function:
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  //  do the line of sight integration
  Vector2D ThetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

  //  spline the result and store it in ThetaT_ell_of_k_spline
  for(int iell=0; iell<nells; iell++){
    std::string name = "Theta_" + std::to_string(ells[iell]);
    ThetaT_ell_of_k_spline[iell].create(k_array, ThetaT_ell_of_k[iell], name);
  }

}




//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================


Vector PowerSpectrum::solve_for_Cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  
  const int nells      = ells.size();

  //  Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  //    ... dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell

  const int n_k = log_k_array.size();
  const double dlogk = log_k_array[1] - log_k_array[0];

  Vector result(nells);

  const double k_p = kpivot_mpc/Constants.Mpc;
  const double fac = 4.*M_PI*A_s * std::pow(k_p, 1.-n_s);

  for(int iell=0; iell<nells; iell++){

    int ell = ells[iell];

    //  Solve using the trapezodial rule

    std::function<double(double)> F = [&](double log_k){
      double k = exp(log_k);
      return std::pow(k, n_s-1.) * f_ell_spline[iell](k)*g_ell_spline[iell](k);
    };

    // double sum = 0;
    // double k = k_min;

    // sum += 0.5 * F(log(k));

    // for(int ik=1; ik<n_k-1; ik++){
    //   sum += F(log_k_array[ik]);
    // }

    // k = k_max;
    // sum += 0.5*F(log(k));
    // sum *= fac;

    result[iell] = fac*integrate_trapezodial(F, log_k_array);
    
  }
  

  return result;
}


double PowerSpectrum::integrate_trapezodial(std::function<double(double)> &F, double z_start, double z_stop, const double dz){

  double sum = 0;
  double z = z_start;

  sum += 0.5 * F(z);

  z += dz;
  while(z<z_stop){
    sum += F(z);
    z += dz;
  }

  z = z_stop;
  sum += 0.5*F(z);

  return sum*dz;

}

double PowerSpectrum::integrate_trapezodial(std::function<double(double)> &F, Vector z_array){

  double sum = 0;
  const int nz = z_array.size();
  const double dz = z_array[1] - z_array[0];
  
  double z = z_array[0];
  sum += 0.5 * F(z);

  for(int iz=1; iz<nz-1; iz++){
    z = z_array[iz];
    sum += F(z);

  }

  z = z_array[nz-1];
  sum += 0.5*F(z);

  return sum*dz;

}




//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  double pofk = 0.0;

  //=============================================================================
  // TODO: Compute the matter power spectrum
  //=============================================================================

  // ...
  // ...
  // ...

  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_Cell_TT(const double ell) const{
  return Cell_TT_spline(ell);
}
double PowerSpectrum::get_Dell(const double ell) const{
  double normfactor = (ell * (ell+1.)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
  return Cell_TT_spline(ell) * normfactor;
}


//====================================================
// Output the Cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(OUTPUT_PATH + filename.c_str());

  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);

  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    // fp << Cell_TT_spline( ell ) * normfactor  << " ";
    fp << get_Dell(ell);

    if(Constants.polarization){
      fp << Cell_EE_spline( ell ) * normfactor  << " ";
      fp << Cell_TE_spline( ell ) * normfactor  << " ";
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}
