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
 * @brief Perform Monte Carlo Markov Chain via the Metropolis algorithm to fit the parameters h, Ω_m amd Ω_K
 * @param supernovadata_filename name of file in input-directory containing observational data: {z, d_L[Gpc], err[Gpc]}
 * @param result_filename name of file in output-directory for which to store the result: {χ^2, h, Ω_m, Ω_K}
 * @return BackgroundCosmology object with new parameters
*/
BackgroundCosmology mcmc_fit_to_supernova_data(std::string supernovadata_filename, std::string result_filename){

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

  const std::array<double, nparam> prior_high {1.5, 1.0, 1.0};  // max{h,  Omegam,  OmegaK}
  const std::array<double, nparam> prior_low {0.5, 0.0, -1.0};  // min{h,  Omegam,  OmegaK}
  
  //  Set starting point for chain and step size

  std::array<double, nparam> parameters{0.7, 0.25, 0.0};  // { h,  Omegam,  OmegaK}
  std::array<double, nparam> stepsize{0.007, 0.05, 0.05}; // {dh, dOmegam, dOmegaK}


  for(int i = 0; i < nparam; i++){
    parameters[i] = prior_low[i] + (prior_high[i]-prior_low[i])*udist(gen);
  }

  std::array<double, nparam> best_parameters = parameters;  // best-fit as we go along in the chain
  double chi2_min = std::numeric_limits<double>::max();     // chi^2 corresponding to the best-fit
  
  /** 
   * @brief The chi^2 function
   * @param parameters the parameters we want to fit
  */
  auto comp_chi2 = [&](std::array<double, nparam> & parameters){
    // Priors: if outside range return huuuuge chi^2
    bool inside_prior = true;
    for(int i = 0; i < nparam; i++){
      if(parameters[i] > prior_high[i]) inside_prior = false;
      if(parameters[i] < prior_low[i]) inside_prior = false;
    }
    if(not inside_prior) return std::numeric_limits<double>::max(); 

    
    //  Set parameters and compute background

    double param_Omegab   = 0.05;                         // unimportant as we just sample and are sensitive to Omegam
    double param_Neff     = 0.0;                          // irrelevant at late times
    double param_TCMB     = 2.7255;                       // temperature of the CMB
    double param_h        = parameters[0];                //
    double param_OmegaCDM = parameters[1] - param_Omegab; // OmegaCDM = Omegam - Omegab
    double param_OmegaK   = parameters[2];                //
    
    BackgroundCosmology cosmo(param_h, param_Omegab, param_OmegaCDM, param_OmegaK, param_Neff, param_TCMB);
    //  solve to get the distance function:
    cosmo.solve_conformal_time();


    //  Compute chi^2

    double chi2 = 0.0;
    double num;  // square root of numerator (attempt to avoid unnecessary FLOPs)
    for (size_t i=0; i<z_obs.size(); i++){
      double x = -std::log(1.0+z_obs[i]);   // x = -ln(1+z)
      double dL = cosmo.get_luminosity_distance_of_x(x) / Gpc; // luminosity function in units of Gpc
      num = dL - dL_obs[i];
      chi2 += num*num / (err_obs[i] * err_obs[i]);
    }
    return chi2;
  };

  //  Generic Metropolis MCMC algorithm

  int steps = 0;
  int nsample = 0;
  double oldchi2 = std::numeric_limits<double>::max();
  std::ofstream out(OUTPUT_PATH + result_filename.c_str());
  // out << "#          chi2            h           Omegam           OmegaK               Acceptrate\n";
  out << "#          chi2            h           Omegam           OmegaK          \n";
  
  while(nsample < maxsteps){
    steps++;
    
    // generate new set of parameters:
    std::array<double, nparam> new_parameters;
    for(int i = 0; i < nparam; i++)
      new_parameters[i] = parameters[i] + ndist(gen) * stepsize[i];

    // compute chi^2:
    double chi2 = 0.0;
    try {
      chi2 = comp_chi2(new_parameters);
    } catch(...){
      chi2 = std::numeric_limits<double>::max();
    }

    //  Evaluate step

    bool accept_step = chi2 < oldchi2;  // always accept sample if lower chi^2
    if(not accept_step){
      //  Metropolis step - draw random number and accept if low enough
      if(std::exp(-(chi2-oldchi2)/2.0) > udist(gen)){
        accept_step = true;
      }
    }

    if(accept_step){ 
      // if step is accepted, change parameters to new parameters and update the old chi^2 value
      nsample++;
      oldchi2 = chi2;
      parameters = new_parameters;

      // write sample to file and screen:
      std::cout << "#          chi2            h           Omegam           OmegaK               Acceptrate\n";
      std::cout << std::setw(15) << chi2 << " ";
      out       << std::setw(15) << chi2 << " ";
      for(int i = 0; i < nparam; i++){
        std::cout << std::setw(15) << parameters[i] << " ";
        out << std::setw(15) << parameters[i] << " ";
      }
      out << "\n";
      std::cout << std::setw(15) << " " << nsample/double(steps)*100.0 << "%\n";
      
      // record new best-fit:
      if(chi2 < chi2_min){
        chi2_min = chi2;
        best_parameters = parameters;
      }
    }
  }

  // Print best-fit
  std::cout << "Minimum chi^2 found " << chi2_min << " ";
  for(int i=0; i<nparam; i++)
    std::cout << best_parameters[i] << " ";
  std::cout << "\n";

  double param_Omegab   = 0.05;                         // unimportant as we just sample and are sensitive to Omegam
  double param_Neff     = 0.0;                          // irrelevant at late times
  double param_TCMB     = 2.7255;                       // temperature of the CMB
  double param_h        = best_parameters[0];                //
  double param_OmegaCDM = best_parameters[1] - param_Omegab; // OmegaCDM = Omegam - Omegab
  double param_OmegaK   = best_parameters[2];                //
  
  BackgroundCosmology cosmo(param_h, param_Omegab, param_OmegaCDM, param_OmegaK, param_Neff, param_TCMB);
  return cosmo;
}