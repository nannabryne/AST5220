#ifndef _POWERSPECTRUM_HEADER
#define _POWERSPECTRUM_HEADER
#ifdef _USEOPENMP
#include <omp.h>
#endif
#include <functional>
#include <utility> 
#include <fstream> 
#include <algorithm>
#include "utils.h"
#include "backgroundcosmology.h"
#include "recombinationhistory.h"
#include "perturbations.h"

using Vector   = std::vector<double>;
using Vector2D = std::vector<Vector>;

class PowerSpectrum {
  private:

    BackgroundCosmology *cosmo = nullptr;
    RecombinationHistory *rec  = nullptr;
    Perturbations *pert        = nullptr;

    // Parameters defining the primordial power-spectrum
    double A_s        = 2.1e-9; // primordial amplitude
    double n_s        = 0.965;  // spectral index
    double kpivot_mpc = 0.05;   // fiducial scale

    // The k-values we compute Theta_ell(k) etc. for
    const int n_k      = 100;
    const double k_min = Constants.k_min;
    const double k_max = Constants.k_max;
    
    // The ells's we will compute Theta_ell and Cell for
    Vector ells{ 
        2,    3,    4,    5,    6,    7,    8,    10,   12,   15,   
        20,   25,   30,   40,   50,   60,   70,   80,   90,   100,  
        120,  140,  160,  180,  200,  225,  250,  275,  300,  350,  
        400,  450,  500,  550,  600,  650,  700,  750,  800,  850,  
        900,  950,  1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 
        1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 
        1900, 1950, 2000};
   
    //=====================================================================
    // [1] Create bessel function splines needed for the LOS integration
    //=====================================================================

    // Splines of bessel-functions for each value of ell in the array above
    std::vector<Spline> j_ell_splines;
    
    // Generate splines of bessel-functions for each ell needed
    // to do the LOS integration
    void generate_bessel_function_splines();
    
    //=====================================================================
    // [2] Do the line of sight integration and spline the result
    //=====================================================================
    
    // Do LOS integration for all ells and all k's in the given k_array
    // and for all the source functions (temperature, polarization, ...)
    void line_of_sight_integration(Vector & k_array);
  
    // Do the line of sight integration for a single quantity
    // for all ells by providing a source_function(x,k) (can be temp, pol, ...)
    Vector2D line_of_sight_integration_single(
        Vector & k_array, 
        std::function<double(double,double)> &source_function);
    
    // Splines of the reusult of the LOS integration
    // Theta_ell(k) and ThetaE_ell(k) for polarization
    std::vector<Spline> ThetaT_ell_of_k_spline;
    std::vector<Spline> ThetaE_ell_of_k_spline;
    
    //=====================================================================
    // [3] Integrate to get power-spectrum
    //=====================================================================
    
    // General method to solve for Cells (allowing for cross-correlations)
    // For auto spectrum (C_TT) then call with f_ell = g_ell = Theta_ell
    // For polarization C_TE call with f_ell = Theta_ell and g_ell = ThetaE_ell
    Vector solve_for_Cell(
        Vector & log_k_array,
        std::vector<Spline> & f_ell, 
        std::vector<Spline> & g_ell);

    // Splines with the power-spectra
    Spline Cell_TT_spline{"Cell_TT_spline"};
    Spline Cell_TE_spline{"Cell_TE_spline"};
    Spline Cell_EE_spline{"Cell_EE_spline"};


    double integrate_trapezodial(std::function<double(double)> &F, double z_start, double z_stop, const double dz);

    double integrate_trapezodial(std::function<double(double)> &F, Vector z_array);


  public:

    // Constructors
    PowerSpectrum() = delete;
    PowerSpectrum(
        BackgroundCosmology *cosmo, 
        RecombinationHistory *rec, 
        Perturbations *pert,
        double A_s,
        double n_s,
        double kpivot_mpc);
    
    // Do all the solving: bessel functions, LOS integration and then compute Cells
    void solve();

    // The dimensionless primordial power-spectrum Delta = 2pi^2/k^3 P(k)
    double primordial_power_spectrum(const double k) const;

    // Get P(k,x) for a given x in units of (Mpc)^3
    double get_matter_power_spectrum(const double x, const double k_mpc) const;

    // Get the quantities we have computed
    double get_Cell_TT(const double ell) const;
    double get_Cell_TE(const double ell) const;
    double get_Cell_EE(const double ell) const;

    double get_Dell(const double ell) const;

    // Output Cells in units of l(l+1)/2pi (muK)^2
    void output(std::string filename) const;
    void output(int ell, std::string filename) const;
};

#endif