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

  //  choose range of k's and the resolution to compute Theta_ell(k)
  const double n_sampling = 32;
  const double dk = 2.*M_PI/n_sampling * 1./cosmo->eta_of_x(0);
  const int n_k = ( k_max - k_min ) / dk + 1;

  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  Vector k_array = exp(log_k_array);

  //  make splines for j_ell
  Utils::StartTiming("besselspline");
  generate_bessel_function_splines();
  Utils::EndTiming("besselspline");


  //  line of sight integration to get Theta_ell(k)
  Utils::StartTiming("lineofsight");
  line_of_sight_integration(k_array);
  Utils::EndTiming("lineofsight");


  //  integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  Utils::StartTiming("C_ell");
  auto Cell_TT = solve_for_Cell(log_k_array, ThetaT_ell_of_k_spline, ThetaT_ell_of_k_spline);
  Cell_TT_spline.create(ells, Cell_TT, "Cell_TT_of_ell");
  Utils::EndTiming("C_ell");



  // (temporary solution)
  //  get the different contributions to Cell
  Utils::StartTiming("C_ell_decomp");

  // const int n_k        = k_array.size();
  const int nells      = ells.size();
  //  make storage for the splines we are to create:
  
  Cell_decomp = std::vector<Spline>(4);

  for(int term=1; term<=4; term++){

    std::cout << term << "..." << std::endl;

    tmp_ThetaT_ell_of_k_spline = std::vector<Spline>(nells);

    // Function returning the source function:
    std::function<double(double,double)> source_function_T = [&](double x, double k){
      return pert->get_Source_T(x,k,term);
    };

    //  do the line of sight integration
    Vector2D ThetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

    //  spline the result and store it in tmp_ThetaT_ell_of_k_spline
    for(int iell=0; iell<nells; iell++){
      tmp_ThetaT_ell_of_k_spline[iell].create(k_array, ThetaT_ell_of_k[iell], "idk");
    }
  
    // solve for Cell:
    auto Cell = solve_for_Cell(log_k_array, tmp_ThetaT_ell_of_k_spline, tmp_ThetaT_ell_of_k_spline);
    Cell_decomp[term-1].create(ells, Cell, "idk");

    std::cout  << "done!" << std::endl;


  }

  Utils::EndTiming("C_ell_decomp");


  


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
  const double n_sampling = 20;
  const double dz = 2*M_PI/n_sampling;
  const int n_z = (int)((z_end - z_start)/dz) + 1;
  Vector z_array = Utils::linspace(z_start, z_end, n_z);
  

  // #pragma omp parallel for schedule(dynamic, 1)
  for(size_t iell = 0; iell < ells.size(); iell++){
    const int ell = ells[iell];
    // Vector z_array = Utils::linspace(z_start, z_end, n_z);

    Vector j_ell_arr(n_z);
    for(int iz=0; iz<n_z; iz++){
      
      j_ell_arr[iz] = Utils::j_ell(ell, z_array[iz]);
    }

    //  make the j_ell_splines[i] spline:
    std::string str = "j_" + std::to_string(int(ell));
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
  const int n_x = int( ( x_end - x_start ) / dx ) + 1;

  const double eta0 = cosmo->eta_of_x(0); // η(0) - conformal time today

  //  For each k, solve the LOS integral for all the ell values

  #pragma omp parallel for schedule(dynamic, 1)
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
  const double dz = z_array[1] - z_array[0];

  return integrate_trapezodial(F, z_array[0], z_array[z_array.size()-1], dz);

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
  //=============================================================================
  // TODO: Compute the matter power spectrum
  //=============================================================================
  const double k = k_mpc / Constants.Mpc;
  const double c = Constants.c;
  const double fac = c*k/cosmo->get_H0();
  const double Hp = cosmo->Hp_of_x(x);

  double DeltaM = fac*fac * 2.* pert->get_Phi(x,k) * exp(x) / (3.*cosmo->get_OmegaM(0));

  double fac2 = 3*Hp/(c*k);
  double Deltab = pert->get_delta_b(x,k) - fac2 * pert->get_u_b(x,k);
  double Deltac = pert->get_delta_c(x,k) - fac2 * pert->get_u_c(x,k);
  // DeltaM = Deltab + Deltac;


  // // ...
  // // ...
  // ...

  double P_R = A_s * pow( k_mpc/kpivot_mpc, n_s-1.) * 2*M_PI*M_PI/ (k_mpc*k_mpc*k_mpc);
  // double P_of_k = DeltaM*DeltaM * primordial_power_spectrum(k) * 2*M_PI*M_PI /(k*k*k);
  double P_of_k = DeltaM*DeltaM * P_R;
  double conv = cosmo->get_h()/Constants.Mpc;
  double h = cosmo->get_h();
  return P_of_k * h*h*h;
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


double PowerSpectrum::get_Dell_comp(const double ell, const int term) const{
  double normfactor = (ell * (ell+1.)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
  return Cell_decomp[term-1](ell) * normfactor;
}



//====================================================
// Output the Cells to file
//====================================================

void PowerSpectrum::output(const std::string filename) const{

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
    fp << get_Dell(ell)                       << " ";
    // fp << get_Dell_comp(ell,1)                     << " ";
    // fp << get_Dell_comp(ell,2)                     << " ";
    // fp << get_Dell_comp(ell,3)                     << " ";
    // fp << get_Dell_comp(ell,4)                     << " ";

    // if(Constants.polarization){
    //   fp << Cell_EE_spline( ell ) * normfactor  << " ";
    //   fp << Cell_TE_spline( ell ) * normfactor  << " ";
    // }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);

}



void PowerSpectrum::output(const std::string filename, 
std::vector<int> & for_ells){

  const int n_ells = for_ells.size();
  std::vector<int> ell_values(n_ells);
  std::vector<int> ell_indices(n_ells);


  for(int i=0; i<n_ells; i++){
    const int ell_target = for_ells[i];
    int ellval=ells[ells.size()-1], ellidx=ells.size()-1;
    int diff = 100;
    int curr_diff, curr_ell;

    for(int iell=0; iell<ells.size(); iell++){
      curr_ell = ells[iell];
      curr_diff = abs(curr_ell - ell_target);
      if(curr_diff == 0){
        ellval = ell_target;
        ellidx = iell;
        break;
      }
      else if(curr_diff < diff){
        diff = curr_diff;
        ellval = curr_ell;
        ellidx = iell;
      }
    }
    ell_values[i] = ellval;
    ell_indices[i] = ellidx;
    if(ellval != ell_target){
      std::cout << "The value of ell that was requested," << ell_target << ", does not exist in 'ells'-array. Using ell = " << ellval << " instead." << std::endl; 
    }

  }

  for_ells = ell_values;

  std::ofstream fp(OUTPUT_PATH + filename.c_str());

  const int npts = 5000;
  const double Mpc = Constants.Mpc;
  const double Mpc_inv = 1./Mpc;

  auto log_k_array = Utils::linspace(log(k_min), log(k_max), npts);
  auto k_array = exp(log_k_array);

  auto print_data = [&] (const double k) {

    double k_SI = k;
    double k_mpc = k_SI*Mpc;

    fp << k_mpc/cosmo->get_h()                  << " ";
    fp << get_matter_power_spectrum(0, k_mpc)   << " ";
    for(int i=0; i<n_ells; i++){
      fp << ThetaT_ell_of_k_spline[ell_indices[i]](k_SI)  << " ";
    }
    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);


}







// void PowerSpectrum::output(int ell, std::string filename) const{
//   const int ell_target = ell;
//   int ellval, ellidx;
//   int diff = 1000;
//   int curr_diff, curr_ell;
//   for(int iell=0; iell<ells.size(); iell++){
//     curr_ell = ells[iell];
//     curr_diff = abs(curr_ell - ell_target);
//     if(curr_diff == 0){
//       ellval = ell_target;
//       ellidx = iell;
//       break;
//     }
//     else if(curr_diff < diff){
//       diff = curr_diff;
//       ellval = curr_ell;
//       ellidx = iell;
//     }
//   }
//   if(ellval != ell_target){
//     std::cout << "The value of ell that was requested does not exist in 'ells'-array. Using ell = " << ellval << " instead." << std::endl; 
//   }


//   std::ofstream fp(OUTPUT_PATH + filename.c_str());

//   const int npts = 5000;
//   const double Mpc = Constants.Mpc;
//   const double Mpc_inv = 1./Mpc;

//   auto log_k_array = Utils::linspace(log(k_min), log(k_max), npts);
//   auto k_array = exp(log_k_array);

//   auto print_data = [&] (const double k) {

//     double k_SI = k;
//     double k_mpc = k_SI*Mpc;

//     fp << k_mpc/cosmo->get_h()                  << " ";
//     fp << get_matter_power_spectrum(0, k_mpc)   << " ";
//     fp << ThetaT_ell_of_k_spline[ellidx](k_SI)  << " ";
//     fp << "\n";
//   };
//   std::for_each(k_array.begin(), k_array.end(), print_data);


// }








