#include "utils.h"
#include "backgroundcosmology.h"
#include "supernovafitting.h"




int milestone_one();
int milestone_two();
int milestone_three();
int milestone_four();


/* Goal main-funtion! */
int mani(int narg, char **argv);

/* Temporary main-function */
int main(int narg, char **argv){

    //=========================================================================
    // Parameters
    //=========================================================================

    // Background parameters
    double h           = 0.67;
    double Omegab      = 0.05;
    double OmegaCDM    = 0.267;
    double Omegak      = 0.0;
    double Neff        = 3.046;
    double TCMB        = 2.7255;

    // Recombination parameters
    double Yp          = 0.245;

    // Power-spectrum parameters
    double A_s         = 2.1e-9;
    double n_s         = 0.965;
    double kpivot_mpc  = 0.05;

    // miletone_one();

     // Set up and solve the background
    BackgroundCosmology cosmo(h, Omegab, OmegaCDM, Omegak, Neff, TCMB);
    cosmo.solve();
    cosmo.info();
    
    // Output background evolution quantities
    cosmo.output("cosmology.txt");



    // std::cout << "Success!" << std::endl;

    // std::ofstream f;
    // f.open(OUTPUT_PATH+"dummy.txt");
    // f << "helo\n";
    // f.close();

    return 0;
}


// int milestone_one(){
//     //=========================================================================
//     // Module I
//     //=========================================================================

//     // Set up and solve the background
//     BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
//     cosmo.solve();
//     cosmo.info();
    
//     // Output background evolution quantities
//     cosmo.output("cosmology.txt");

//     // Do the supernova fits. Uncomment when you are ready to run this
//     // Make sure you read the comments on the top of src/SupernovaFitting.h
//     // mcmc_fit_to_supernova_data("data/supernovadata.txt", "results_supernovafitting.txt");

//     return 0;
// }
















int mani(int narg, char **argv){

    std::string keyword;
    int milestone; 

    // //  Option to do milstone I, II, III and IV
    // if(narg < 2){
    //     milestone = 0;
    // }
    // else{
    //     milestone = atoi(argv[1]);
    // }

    // if(milestone == 1){
    //     milestone_one();
    // }
    // // else if(milestone == 2){

    // // }
    // // else if(milestone == 3){
        
    // // }
    // // else if(milestone == 4){
        
    // // }
    // else if(milestone == 0){
    //     milestone_one();
    //     //milestone_two();
    //     //...
    // }

    return 0;

}