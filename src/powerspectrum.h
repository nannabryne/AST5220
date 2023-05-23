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

    BackgroundCosmology *cosmo = nullptr;   // background
    RecombinationHistory *rec  = nullptr;   // recombination history
    Perturbations *pert        = nullptr;   // perturbations

    // Parameters defining the primordial power-spectrum
    double A_s        = 2.1e-9; // primordial amplitude
    double n_s        = 0.965;  // spectral index
    double kpivot_mpc = 0.05;   // fiducial scale

    // The k-values we compute Theta_ell(k) etc. for
    const int n_k        = 200;             // #k's to consider
    const double k_min   = Constants.k_min; // minimum wavenumber k
    const double k_max   = Constants.k_max; // maximum wavenumber k
    
    // The ells's we will compute Theta_ell and Cell for
    Vector ells{ 
        2,    3,    4,    5,    6,    7,    8,    10,   12,   15,   
        20,   25,   30,   40,   50,   60,   70,   80,   90,   100,  
        120,  140,  160,  180,  200,  225,  250,  275,  300,  350,  
        400,  450,  500,  550,  600,  650,  700,  750,  800,  850,  
        900,  950,  1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 
        1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 
        1900, 1950, 2000};

    

    // Normalisation factor without ell's
    const double __normfactor_ish = 1. / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);


    // Splines of bessel-functions for each value of ell in the array above
    std::vector<Spline> j_ell_splines;  // j_ℓ(z) spline

    // Splines of the reusult of the LOS integration
    // Theta_ell(k) and ThetaE_ell(k) for polarization
    std::vector<Spline> ThetaT_ell_of_k_spline;     // Θ_ℓ(0,k) spline

    //  ... has to be a better way....
    std::vector<Spline> tmp_ThetaT_ell_of_k_spline;


    // Splines with the power-spectra
    Spline Cell_TT_spline{"Cell_TT_spline"};
    std::vector<Spline> Cell_decomp;
   

    /*
    Workflow:
        [1] Create bessel function splines needed for the LOS integration
        [2] Do the line of sight integration and spline the result
        [3] Integrate to get power-spectrum
    */

    /**
     * @brief Generate splines of bessel-functions for each ell needed to do the LOS integration.
     */
    void generate_bessel_function_splines();
    
    
    /**
     * @brief Do LOS integration for all ells and all k's in the given k_array.
     * @param k_array values to integrate over
     */
    void line_of_sight_integration(Vector & k_array);
  

    /**
     * @brief Do the line of sight integration for a single quantity for all ells.
     * @param k_array values to integrate over
     * @param source_function source function
     * @return Θ_ℓ(0,k) 
     */
    Vector2D line_of_sight_integration_single(
        Vector & k_array, 
        std::function<double(double,double)> &source_function);
    
    


    /**
     * @brief General method to solve for Cells (allowing for cross-correlations).
     * @param log_k_array the log(k)-values to integrate over
     * @param f_ell transfer function one   ( for us: Θ_ℓ(0,k) )
     * @param g_ell transfer function two   ( for us: Θ_ℓ(0,k) )
     * @return Cell 
     */
    Vector solve_for_Cell(
        Vector & log_k_array,
        std::vector<Spline> & f_ell, 
        std::vector<Spline> & g_ell);

    

    /**
     * @brief Integrate a function F(z) using the trapezoidal rule for a uniform grid z.
     * @param F integrand
     * @param z_start lower limit in integration
     * @param z_stop upper limit in integration
     * @param dz step size in z-dir.
     * @return result 
     */
    double integrate_trapezoidal(std::function<double(double)> &F, double z_start, double z_stop, const double dz);

    /**
     * @brief Integrate a function F(z) using the trapezoidal rule for a uniform grid z.
     * @param F integrand
     * @param z_array the values of z to integrate over (must be linearly spaced)
     * @return result
     */
    double integrate_trapezoidal(std::function<double(double)> &F, Vector z_array);


  public:

    // Constructors:

    PowerSpectrum() = delete;
    /**
     * @brief Create an instance of the PowerSpectrum-class with defined background, recombination history and perturbations.
     * @param cosmo instance of the Backgroundcosmology-class, solved 
     * @param rec instance of the RecombinationHistory-class, solved
     * @param pert instance of the Perturbations-class, solved
     * @param A_s primordial amplitude
     * @param n_s scalar spectral index
     * @param k_p pivot scale
     * @return Perturbations-object
    */
    PowerSpectrum(
        BackgroundCosmology *cosmo, 
        RecombinationHistory *rec, 
        Perturbations *pert,
        double A_s,
        double n_s,
        double kpivot_mpc);
    
    // Do all the solving: bessel functions, LOS integration and then compute Cells
    // Methods:

    /**
     * @brief Do all the solving.
     * @details [1] Bessel functions, [2] LOS integration and then [3] compute Cells
     */
    void solve();

    /**
     * @brief Solve for C(ℓ) by component of the source function. Slow and stupid method...
     */
    void solve_decomposed();



    /**
     * @brief Get the matter power spectrum by component.
     * 
     * @param x time point x = ln(a)
     * @param k_mpc comoving wavenumber in Mpc^(-1)
     * @param comp component key: 1 = CDM, 2 = baryons, 3 = photons
     * @return P_s(k) 
     */
    double get_matter_power_spectrum_comp(const double x, const double k_mpc, int comp) const;

    /**
     * @brief Get the total matter power spectrum P_m(x,k) in units of Mpc^3 h^(-3).
     * @param x time point x = ln(a)
     * @param k_mpc comoving wavenumber in Mpc^(-1)
     * @return P_m(x,k)
     */
    double get_matter_power_spectrum(const double x, const double k_mpc) const;

    /**
     * @brief Get the primordial power spectrum P_R(k) in units of Mpc^3 h^(-3).
     * @param k_mpc comoving wavenumber in Mpc^(-1)
     * @return P_R(k)
     */
    double get_primordial_power_spectrum(const double k_mpc) const;

    /**
     * @brief Get the dimensionless primordial power spectrum Δ^2_R(k).
     * @param k_mpc comoving wavenumber in Mpc^(-1)
     * @return Δ^2_R(k)
     */
    double Delta2_R_of_k(const double k) const;

    /**
     * @brief Get the primordial power spectrum P_R(k).
     * @param k_mpc comoving wavenumber in Mpc^(-1)
     * @return P_R(k)
     */
    double P_R_of_k(const double k) const;


    /**
     * @brief Get the CMB anisotropy spectrum C(ℓ).
     * @param ell multipole moment ℓ
     * @return C(ℓ)
     */
    double get_Cell_TT(const double ell) const;

    /**
     * @brief Get the scaled CMB anisotropy spectrum D(ℓ).
     * @param ell multipole moment ℓ
     * @return D(ℓ) 
     */
    double get_Dell(const double ell) const;

    /**
     * @brief Get the scaled CMB anisotropy spectrum D(ℓ) by component.
     * @param ell multipole moment ℓ
     * @param term component key: 1 = SW, 2 = ISW, 3 = Doppler, 4 = pol
     * @return D(ℓ)^[comp]
     */
    double get_Dell_comp(const double ell, const int term) const;


    /**
     * @brief Print some (not so) useful information about the class.
     */
    void info(); 


    /**
     * @brief Output D(ℓ) in units of (muK)^2.
     * @param filename the scaled CMB anisotropy spectrum
     */
    void output(const std::string filename) const;

    /**
     * @brief Output D(ℓ) in units of (muK)^2 by contributions from the different terms in the source function.
     * @param filename name of file to be stored in output-directory
     */
    void output_decomposed(const std::string filename) const;

    /**
     * @brief Output results for today as functions of wavenumber.
     * @param filename the scaled CMB anisotropy spectrum
     * @param for_ells the ℓ's for which we want Θ_ℓ(0,k)
     */
    void output(const std::string filename, std::vector<int> & for_ells);
};

#endif