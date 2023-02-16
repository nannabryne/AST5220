//====================================================================================================
//
// Method to be used with code template to allow you to do MCMC and get constraints from supernova data
// 
// Things needed to be done to use this:
// * Check that line 109-116 is how you set up your background class and do the solving
// * Edit line 124 with the call to your luminosity function method (mine is called get_luminosity_distance_of_x). 
//   As written here this assumes that the luminosity distance is returned in meters from this method
// * Include this file in Main.cpp: #include "SupernovaFitting.h"
// * Call the function below in main: mcmc_fit_to_supernova_data("supernovadata.txt", "results.txt");
//
// The input [supernovadata] is the path to the file supernovadata.txt. The output file with all the samples is [result_filename]
// This runs for [maxsteps] steps before ending. You can reduce maxsteps or just kill the run if it takes too long
//
// How to analyze the resulting chains:
// * Load the chains and skip the first few hundred samples (the burnin of the chains). E.g. loadtxt(file,skiprows=200) in python
// * Find the minimum chi2 and the corresponding best-fit parameters (you can use np.argmin to get index of the minvalue in python)
// * Select all samples of Omegam and OmegaLambda (computed from Omegam and Omegak) that satisfy chi2 < chi2_min + 3.53 
//   (e.g. Omegam[chi2 < chi2min + 3.53] in python)
// * Scatterplotting these gives you the 1sigma (68.4%) confidence region
// * Find the standard deviation of the samples to get the 1sigma confidence region of the parameters (assuming the posterior is a gaussian)
// * Make and plot a histogram of the samples for the different parameters (Omegam, Omegak, OmegaLambda, H0)
// * You can also compute the mean and standard deviation of the chain values and use this to overplot a gaussian with the same mean and variance for comparison.
//
//====================================================================================================


/*
Based on template by Hans Winther (13/02/22).
  > Changed some variable names to be consistent with my own convention/notation.

INPUT_PATH and OUTPUT_PATH specified in "utils.h"
*/

#include "utils.h"
#include <climits>
#include <random>
#include <array>
#include <vector>



/**
 * @brief Perform Monte Carlo Markov Chain via the Metropolis algorithm to fit the parameters h, Ω_m amd Ω_k
 * @param supernovadata_filename name of file in input-directory containing observational data: {z, d_L[Gpc], err[Gpc]}
 * @param result_filename name of file in output-directory for which to store the result: {χ^2, h, Ω_m, Ω_k}
*/
void mcmc_fit_to_supernova_data(std::string supernovadata_filename, std::string result_filename){

  const int nparam = 3;                     // #parameters we want to fit
  const int maxsteps = 10000;               // maximal #samples to generate
  const int seed = 13092;                   // my random random seed
  const double Gpc = Constants.Mpc * 1000;  // #meters in a Gpc
  
  //  Read observational data from 'supernova_filename'

  std::vector<double> z_obs;    // observed redshift
  std::vector<double> dL_obs;   // observed luminosity distance
  std::vector<double> err_obs;  // error in observed luminosity distance
  
  auto read_data = [&](std::string filename){
    std::string header;
    std::ifstream fp(INPUT_PATH + filename.c_str());
    if(!fp)
      throw std::runtime_error("Error: cannot open file " + filename);
    std::getline(fp, header);
    std::cout << "Reading luminosity data from file:\n";
    while(1){
      //  read line by line:
      double zi, dLi, erri;
      fp >> zi;
      if(fp.eof()) break;
      fp >> dLi;
      fp >> erri;
      z_obs.push_back(zi);
      dL_obs.push_back(dLi);
      err_obs.push_back(erri);
      std::cout << "z: " << zi << " " << dLi << " " << erri << "\n";
    }
    std::cout << "We found n = " << z_obs.size() << " points\n";
  };

  read_data(supernovadata_filename);


  //  Set up random number generator

  std::default_random_engine gen;
  gen.seed(seed);
  std::normal_distribution<double> ndist(0.0,1.0);      // normal distribution with stdv. 1
  std::uniform_real_distribution<double> udist(0.0,1.0);// uniform distribution of numbers between 0 and 1

  //  Define our priors (chi^2 -> Inf if outside these ranges)

  const std::array<double, nparam> prior_high {1.5, 1.0, 1.0};  // max{ h,  Omegam,  Omegak}
  const std::array<double, nparam> prior_low {0.5, 0.0, -1.0};  // min{ h,  Omegam,  Omegak}
  
  //  Set starting point for chain and step size

  std::array<double, nparam> parameters{0.7, 0.25, 0.0};  // { h,  Omegam,  Omegak}
  std::array<double, nparam> stepsize{0.007, 0.05, 0.05}; // {dh, dOmegam, dOmegak}

  /* old params ...*/
  // std::array<double, nparam> parameters{0.7, 0.3, -0.5};
  // std::array<double, nparam> stepsize{0.005, 0.05, 0.05};

  for(int i = 0; i < nparam; i++){
    parameters[i] = prior_low[i] + (prior_high[i]-prior_low[i])*udist(gen);
  }

  std::array<double, nparam> best_parameters = parameters;  // best-fit as we go along in the chain
  double chi2_min = std::numeric_limits<double>::max();     // Chi^2 corresponding to the best-fit
  
  // The chi^2 function
  auto comp_chi2 = [&](std::array<double, nparam> & parameters){
    // Priors: if outside range return huuuuge chi^2
    bool inside_prior = true;
    for(int i = 0; i < nparam; i++){
      if(parameters[i] > prior_high[i]) inside_prior = false;
      if(parameters[i] < prior_low[i]) inside_prior = false;
    }
    if(not inside_prior) return std::numeric_limits<double>::max(); 

    //=========================================================================================
    //======= Here we set up the cosmology class and solve to get the distance functions ======
    //=========================================================================================
    // Set parameters and compute background
    double param_Omegab   = 0.05;           // Not important as we just sample and are sensitive to Omegam = Omegab+OmegaCDM
    double param_Neff     = 0.0;            // Not relevant at late times
    double param_TCMB     = 2.7255;         // Temperature of the CMB
    double param_h        = parameters[0];
    double param_OmegaCDM = parameters[1] - param_Omegab; // OmegaCDM = Omegam - Omegab
    double param_Omegak   = parameters[2];
    BackgroundCosmology cosmo(param_h, param_Omegab, param_OmegaCDM, param_Omegak, param_Neff, param_TCMB);
    cosmo.solve();
    //=========================================================================================

    // Compute chi^2
    double chi2 = 0.0;
    for (size_t i = 0; i < z_obs.size(); i++){
      double x = -std::log(1.0+z_obs[i]);
      //======= Here we call your distance function ======
      double dL = cosmo.get_luminosity_distance_of_x(x) / Gpc; // Luminosity function in units of Gpc
      //==================================================
      chi2 += (dL - dL_obs[i]) * (dL - dL_obs[i]) / (err_obs[i] * err_obs[i]);
    }
    return chi2;
  };

  // Generic Metropolis MCMC algorithm. This requires no changes
  int steps = 0;
  int nsample = 0;
  double oldchi2 = std::numeric_limits<double>::max();
  std::ofstream out(OUTPUT_PATH + result_filename.c_str());
  out << "#          Chi2            h           Omegam           Omegak               Acceptrate\n";
  while(nsample < maxsteps){
    steps++;
    
    // Generate new set of parameters
    std::array<double, nparam> new_parameters;
    for(int i = 0; i < nparam; i++)
      new_parameters[i] = parameters[i] + ndist(gen) * stepsize[i];

    // Compute chi^2
    double chi2 = 0.0;
    try {
      chi2 = comp_chi2(new_parameters);
    } catch(...){
      chi2 = std::numeric_limits<double>::max();
    }

    // Always accept sample if lower chi^2
    bool accept_step = chi2 < oldchi2;
    if(not accept_step){
      // Metropolis step - draw random number and accept if low enough
      if(std::exp(-(chi2-oldchi2)/2.0) > udist(gen)){
        accept_step = true;
      }
    }

    // If the step is accepted change parameters to new parameters and update the old chi^2 value
    if(accept_step){
      nsample++;
      oldchi2 = chi2;
      parameters = new_parameters;

      // Write sample to file and screen
      std::cout << "#          Chi2            h           Omegam           Omegak               Acceptrate\n";
      std::cout << std::setw(15) << chi2 << " ";
      out       << std::setw(15) << chi2 << " ";
      for(int i = 0; i < nparam; i++){
        std::cout << std::setw(15) << parameters[i] << " ";
        out << std::setw(15) << parameters[i] << " ";
      }
      out << "\n";
      std::cout << std::setw(15) << " " << nsample/double(steps)*100.0 << "%\n";
      
      // Record new best-fit 
      if(chi2 < chi2_min){
        chi2_min = chi2;
        best_parameters = parameters;
      }
    }
  }

  // Print best-fit
  std::cout << "Minimum Chi^2 found " << chi2_min << " ";
  for(int i = 0; i < nparam; i++)
    std::cout << best_parameters[i] << " ";
  std::cout << "\n";
}