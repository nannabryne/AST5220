#include "perturbations.h"

//  ----------------------------------
//  CONSTRUCT
//  ----------------------------------

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//  ----------------------------------
//  SOLVE
//  ----------------------------------

void Perturbations::solve(){

  //  integrate all the perturbation equation and spline the result:
  Utils::StartTiming("integrateperturbation");
  integrate_perturbations();
  Utils::EndTiming("integrateperturbation");

  //  compute source functions and spline the result:
  Utils::StartTiming("source");
  compute_source_functions();
  Utils::EndTiming("source");
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  

  //===================================================================
  // TODO: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================

  //  set up k-array for the k's we are to integrate over:
  Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), n_k));
  

  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  Vector Phi_array(n_x*n_k);

  //  loop over all wavenumbers:  
  for(int ik=0; ik<n_k; ik++){

    //  progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    double k = k_array[ik];   // current value of k

    //  find value to integrate to:
    double x_end_tight = get_tight_coupling_time(k);
    double idk = n_x*(x_end_tight-x_start)/(x_end-x_start);
    int n_x_tight = int(idk);
    Vector x_array_tight = Utils::linspace(x_start, x_end_tight, n_x_tight);
    Vector x_array_full = Utils::linspace(x_end_tight, x_end, n_x-n_x_tight);

    std::cout << "x_end_tight: " <<  x_end_tight << std::endl;
    std::cout << "n_x_tight: " <<  n_x_tight << std::endl;


    //===================================================================
    // TODO: Tight coupling integration
    // Remember to implement the routines:
    // set_ic : The IC at the start
    // rhs_tight_coupling_ode : The dydx for our coupled ODE system
    //===================================================================
    

    //  set up initial conditions in the tight coupling regime:
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system

    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

   
    // Integrate from x_start -> x_end_tight
    // ...
    // ...
    // ...
    // ...
    // ...

    ODESolver ode_tight;
    ode_tight.solve(dydx_tight_coupling, x_array_tight, y_tight_coupling_ini);
    Vector y_tight_coupling = ode_tight.get_final_data();

    std::cout << y_tight_coupling[0] << ", " << y_tight_coupling[2] << std::endl;


    //====i===============================================================
    // TODO: Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_end_tight, k);


    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Integrate from x_end_tight -> x_end
    // ...
    // ...
    // ...
    // ...
    // ...
    ODESolver ode_full;
    ode_full.solve(dydx_full, x_array_full, y_full_ini);
    Vector y_full = ode_full.get_data_by_component(0);


    



    //===================================================================
    // TODO: remember to store the data found from integrating so we can
    // spline it below
    //
    // To compute a 2D spline of a function f(x,k) the data must be given 
    // to the spline routine as a 1D array f_array with the points f(ix, ik) 
    // stored as f_array[ix + n_x * ik]
    // Example:
    // Vector x_array(n_x);
    // Vector k_array(n_k);
    // Vector f(n_x * n_k);
    // Spline2D y_spline;
    // f_spline.create(x_array, k_array, f_array);
    // We can now use the spline as f_spline(x, k)
    //
    // NB: If you use Theta_spline then you have to allocate it first,
    // before using it e.g.
    // Theta_spline = std::vector<Spline2D>(n_ell_Theta);
    //
    //===================================================================
    //...
    //...
    int d = n_x*ik;

    Vector Phi__ = ode_tight.get_data_by_component(Constants.ind_Phi);
    
    // std::cout << Phi__[40] << std::endl;


    for(int i=0; i<n_x_tight; i++){
      // Phi_array[i + d] = y_tight_coupling[Constants.ind_Phi];
      Phi_array[i + d] = ode_tight.get_data_by_component(Constants.ind_Phi)[i];
      // std::cout << Phi_array[i+d] << std::endl;
    }
    for(int i=n_x_tight; i<n_x; i++){
      Phi_array[i + d] = ode_full.get_data_by_component(Constants.ind_Phi)[i];//y_full[Constants.ind_Phi];

    }

    // Phi_array[n_x*ik] = y_full[Constants.ind_Phi];

    // std::cout << Phi_array[n_x*ik] << std::en



  }
  

  //=============================================================================
  // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
  // ...
  // ...
  // ...

  Phi_spline.create(x_array, k_array, Phi_array);
  // std::cout << " ok" << std::endl;

}



//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  Vector y_tc(Constants.n_ell_tot_tc);  // the vector we are going to fill

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)


  const int n_ell_Theta_tc      = Constants.n_ell_Theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  //  reference to the tight coupling quantities:
  double &delta_c      =  y_tc[Constants.ind_deltac_tc];        // δ_c(x,k) 
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];        // δ_b(x,k)
  double &u_c          =  y_tc[Constants.ind_uc_tc];            // u_c(x,k)
  double &u_b          =  y_tc[Constants.ind_ub_tc];            // u_b(x,k)
  double &Phi          =  y_tc[Constants.ind_Phi_tc];           // Φ(x,k) 
  double *Theta        = &y_tc[Constants.ind_start_Theta_tc];   // Θ_ℓ(x,k)
  double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  
  //  fetch cosmological parameters and variables:
  double Hp     = cosmo->Hp_of_x(x);        // Hp(x)
  double dHpdx  = cosmo->dHpdx_of_x(x);     // d/dx[Hp(x)]

  //  fetch recombination variables:
  double R          = rec->R_of_x(x);         // R(x)
  double dtaudx     = rec->dtaudx_of_x(x);    // d/dx[τ(x)]
  double ddtaudxx   = rec->ddtaudxx_of_x(x);  // d^2/dx^2[τ(x)]

  //  fetch perturbation variables:
  double Psi = -2./3.;    // Ψ(x=x_init,k)

  //  define useful factors:
  double ckH      = Constants.c*k/Hp;   // ck/Hp(x)
  double Rplus    = (1+R);              // ( 1 + R )
  double Rminus   = (1-R);              // ( 1 - R )
  double dtau_inv = 1./dtaudx;          // (d/dx[τ(x)])^(-1)


  //  Set initial conditions in the tight coupling regime

  //  compute scalar quantities (Φ, δ, u):

  Phi = - Psi; 

  delta_c = -3./2. * Psi;
  delta_b = delta_c; 

  u_c = -0.5*ckH * Psi;
  u_b = u_c;

  // compute photon temperature perturbations (Θ_ℓ):

  Theta[0] = -0.5 * Psi;
  Theta[1] = ckH/6. * Psi;

  return y_tc;
}


//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  Vector y(Constants.n_ell_tot_full);   // vector to fill
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime


  const int n_ell_Theta         = Constants.n_ell_Theta;
  const int n_ell_Thetap        = Constants.n_ell_Thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime


  const int n_ell_Theta_tc      = Constants.n_ell_Theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  //  reference to the tight coupling quantities:
  const double &delta_c_tc      =  y_tc[Constants.ind_deltac_tc];       // δ_c(x,k) from tight coupling
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];       // δ_b(x,k) from tight coupling
  const double &u_c_tc          =  y_tc[Constants.ind_uc_tc];           // u_c(x,k) from tight coupling
  const double &u_b_tc          =  y_tc[Constants.ind_ub_tc];           // u_b(x,k) from tight coupling
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];          // Φ(x,k) from tight coupling
  const double *Theta_tc        = &y_tc[Constants.ind_start_Theta_tc];  // Θ_ℓ(x,k) from tight coupling
  const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  //  reference to the quantities we are going to set:
  double &delta_c         =  y[Constants.ind_deltac_tc];        // δ_c(x,k) after tight coupling
  double &delta_b         =  y[Constants.ind_deltab_tc];        // δ_b(x,k) after tight coupling
  double &u_c             =  y[Constants.ind_uc_tc];            // u_c(x,k) after tight coupling
  double &u_b             =  y[Constants.ind_ub_tc];            // u_b(x,k) after tight coupling
  double &Phi             =  y[Constants.ind_Phi_tc];           // Φ(x,k) after tight coupling
  double *Theta           = &y[Constants.ind_start_Theta_tc];   // Θ_ℓ(x,k) after tight coupling
  double *Theta_p         = &y[Constants.ind_start_Thetap_tc];
  double *Nu              = &y[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================
  // ...
  // ...
  // ...

  //  fetch constants and variables:
  const double c              = Constants.c;              // c
  const double H0             = cosmo->get_H0();          // H_0
  const double H0H0           = H0*H0;                    // H_0^2
  const double Omegagamma0    = cosmo->get_Omegagamma();  // Ω_γ0


  //  define useful variables:
  double ck     = ck;           // ck
  double a      = exp(x);       // e^x

  double Psi        =  -Phi - 12*H0H0/(ck*ck)*exp(-2*x) * Omegagamma0*Theta[2];  // Ψ(x,k)
  double dtau_inv   = 1./rec->dtaudx_of_x(x); // (d/dx[τ(x)])^(-1)


  //  compute scalar quantities (Φ, δ, u):
  Phi = Phi_tc;

  delta_c = delta_c_tc;
  delta_b = delta_b_tc;

  u_c = u_c_tc;
  u_b = u_b_tc;


  // compute photon temperature perturbations (Θ_ℓ):

  //FIXME
  int ell       = 0;                // ℓ
  int ell_max   = n_ell_Theta - 1;  // ℓ_max

  for(ell=0; ell<2; ell++){
    Theta[ell] = Theta_tc[ell];
  }

  double U = Constants.c*k/cosmo->Hp_of_x(x) * dtau_inv; 
  Theta[2] = -8./15.*U * Theta[1]; //idk
  for(ell=3; ell<=ell_max; ell++){
    Theta[ell] = - ell/(2*ell + 1) * U * Theta[ell-1];
  }

  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

double Perturbations::get_tight_coupling_time(const double k) const{
  

  //=============================================================================
  // TODO: compute and return x for when tight coupling ends
  // Remember all the three conditions in Callin
  //=============================================================================
  // ...
  // ...

  //  Find x for when the tight coupling is over

  double x_rec = -8.3;
  int npts = int((1e6+1)*(x_rec-x_start)/(x_end-x_start));
  Vector x_array = Utils::linspace(x_start, x_rec, npts);
  
  bool one, two, three;
  double ck = Constants.c*k;
  double x, val, lim;

  double x_tight_coupling_end = x_rec;
  
  for(int i=0; i<npts; i++){
    x = x_array[i];
    
    lim = ck/cosmo->Hp_of_x(x);
    if(lim>1)
      lim = 1;
    lim *= 10.;

    val = abs(rec->dtaudx_of_x(x));

    if(val<=10){
      x_tight_coupling_end = x;
      break;
    }

  }


  return x_tight_coupling_end;
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  // ...
  // ...
  Vector k_array;
  Vector x_array;

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)

  Vector ST_array(k_array.size() * x_array.size());   // temperature source
  Vector SE_array(k_array.size() * x_array.size());

  //  fetch constants and variables:
  const double c              = Constants.c;              // c
  const double H0             = cosmo->get_H0();          // H_0
  const double H0H0           = H0*H0;                    // H_0^2
  const double Omegagamma0    = cosmo->get_Omegagamma();  // Ω_γ0
  const double Omegac0        = cosmo->get_Omegac();      // Ω_c0
  const double Omegab0        = cosmo->get_Omegab();      // Ω_b0

  //  compute source functions:
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    double a = exp(x);

    //  fetch all the things we need...:
    const double Hp       = cosmo->Hp_of_x(x);
    const double dHpdx    = cosmo->dHpdx_of_x(x);
    const double tau      = rec->tau_of_x(x);
    const double dtaudx   = rec->dtaudx_of_x(x);
    const double gt       = rec->gt_of_x(x);
    const double dgtdx    = rec->dgtdx_of_x(x);

    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================
      
      //  fetch the rest of the things we need...:
      double ck             = c*k;
      double ck_inv         = 1./ck;
      double ckH            = ck/Hp;
      const double Psi      = get_Psi(x, k);        // Ψ(x,k)
      const double Phi      = get_Phi(x, k);  
      const double Theta0   = get_Theta(x, k, 0);
      const double ub       = get_u_b(x, k);

      double term1, term2, term3;

      term1 = gt * (Theta0 + Psi);


      // OBS!!
      double dubdx = -ub - ckH*Psi + dtaudx ;
      term2 = ck_inv * ( dgtdx*ub*Hp + gt*(dubdx*Hp + ub*dHpdx) );
    
      double dPsidx = - Phi - 12*H0*H0*ck_inv*ck_inv*exp(-2*x)*Omegagamma0*get_Theta(x, k, 2);
      double dPhidx = Psi - 1./3*ckH*ckH*Phi + 0.5*H0H0/(Hp*Hp*a*a) * (Omegac0*get_delta_c(x, k) + Omegab0*get_delta_b(x, k))*a + 4*Omegagamma0*Theta0;
      term3 = exp(-tau)* ( dPsidx - dPhidx );

      // Temperature source
      ST_array[index] = Hp/c*(term1+term2+term3);

      // Polarization source
      if(Constants.polarization){
        SE_array[index] = 0.0;
      }
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }
}





//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)

  const int n_ell_Theta_tc      = Constants.n_ell_Theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

  //  reference to the different quantities in the y-array:
  const double &delta_c         =  y[Constants.ind_deltac_tc];      // δ_c(x,k)
  const double &delta_b         =  y[Constants.ind_deltab_tc];      // δ_b(x,k)
  const double &u_c             =  y[Constants.ind_uc_tc];          // u_c(x,k)
  const double &u_b             =  y[Constants.ind_ub_tc];          // u_b(x,k)
  const double &Phi             =  y[Constants.ind_Phi_tc];         // Φ(x,k)
  const double *Theta           = &y[Constants.ind_start_Theta_tc]; // Θ_ℓ(x,k)
  const double *Nu              = &y[Constants.ind_start_nu_tc];

  //  reference to the quantities we are going to set in the dydx array:
  double &ddelta_cdx      =  dydx[Constants.ind_deltac_tc];       // d/dx[δ_c(x,k)]
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];       // d/dx[δ_b(x,k)]
  double &du_cdx          =  dydx[Constants.ind_uc_tc];           // d/dx[u_c(x,k)]
  double &du_bdx          =  dydx[Constants.ind_ub_tc];           // d/dx[u_b(x,k)]
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];          // d/dx[Φ(x,k)]
  double *dThetadx        = &dydx[Constants.ind_start_Theta_tc];  // d/dx[Θ_ℓ(x,k)]
  double *dNudx           = &dydx[Constants.ind_start_nu_tc];     

  //  fetch cosmological parameters and variables:
  double Hp             = cosmo->Hp_of_x(x);          // Hp(x)
  double dHpdx          = cosmo->dHpdx_of_x(x);       // d/dx[Hp(x)]
  double H0             = cosmo->get_H0();            // H_0
  double Omegagamma0    = cosmo->get_Omegagamma();    // Ω_γ0
  double Omegac0        = cosmo->get_Omegac();        // Ω_c0
  double Omegab0        = cosmo->get_Omegab();        // Ω_b0 
  double eta            = cosmo->eta_of_x(x);         // η(x)

  //  fetch recombination variables::
  double R          = rec->R_of_x(x);           // R(x)
  double dtaudx     = rec->dtaudx_of_x(x);      // d/dx[τ(x)]
  double ddtaudxx   = rec->ddtaudxx_of_x(x);    // d^2/dx^2[τ(x)]

  

  // useful factors:
  double c      = Constants.c;  // c
  double ck     = ck;           // ck
  double ckH    = ck/Hp;        // ck/Hp
  double Rplus  = (1+R);        // ( 1 + R )
  double Rminus = (1-R);        // ( 1 - R )
  double H0H0   = H0*H0;        // H_0^2
  double a      = exp(x);       // e^x

  double expr;

  //  fetch perturbation variables:
  // double Psi = get_Psi(x, k);   // Ψ(x,k)
  double Psi = -Phi - 12*H0H0/(ck*ck*a*a) * Omegagamma0*Theta[2];  // Ψ(x,k)
  

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  
  // reocurring terms:
  double U1 = ( 1 - dHpdx/Hp );             // ( 1 - d/dx[Hp] )
  double U2 = ( -Theta[0] + 2*Theta[2] );   // ( - Θ_0 + 2Θ_2 )


  //  compute scalar quantities (Φ, δ, u):

  expr = (Omegac0*delta_c + Omegab0*delta_b)*a + 4*Omegagamma0*Theta[0];
  dPhidx = Psi - 1./3*ckH*ckH*Phi + 0.5*H0H0/(Hp*Hp*a*a) * expr;
  
  expr = 3*dPhidx;
  ddelta_cdx = ckH*u_c - expr;
  ddelta_bdx = ckH*u_b - expr;
  
  expr = ckH*Psi;
  du_cdx = - u_c - expr;

  double denom = ( dtaudx*Rplus - R*U1 );
  double term1 = ( ddtaudxx*Rplus - dtaudx*Rminus ) * (3*Theta[1] + u_b);
  double term2 = ckH*R * ( Psi + U1*U2 - dThetadx[0]);
  double q = (term1+term2)/denom;

  du_bdx = ( q - R*u_b + ckH*U2)/Rplus - expr;


  //  compute photon temperature perturbations (Θ_ℓ):

  dThetadx[0] = -ckH * Theta[1] - dPhidx;
  dThetadx[1] = 1./3. * ( q - du_bdx);

  // int ell     = 2;                // ℓ
  // int ell_max = n_ell_Theta - 1;  // ℓ_max

  // for(ell; ell<ell_max; ell++){
  //   expr_a = ckH / (2*ell+1) * ( ell*Theta[ell-1] - (ell+1)*ell*Theta[ell+1] );
  //   expr_b = dtaudx * Theta[ell];
  //   dThetadx[ell] = expr_a + expr_b;
  // }
  // dThetadx[2] -= 0.1*dtaudx*Theta[2];
  
  // ell = ell_max;
  // dThetadx[ell] = ckH*dThetadx[ell-1] - (c*(ell+1)/(Hp*eta) + dtaudx)*Theta[ell];

  return GSL_SUCCESS;
}



//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_Theta         = Constants.n_ell_Theta;
  const int n_ell_Thetap        = Constants.n_ell_Thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_c         =  y[Constants.ind_deltac];       // δ_c(x,k)
  const double &delta_b         =  y[Constants.ind_deltab];       // δ_b(x,k)
  const double &u_c             =  y[Constants.ind_uc];           // u_c(x,k)
  const double &u_b             =  y[Constants.ind_ub];           // u_c(x,k)
  const double &Phi             =  y[Constants.ind_Phi];          // Φ(x,k)
  const double *Theta           = &y[Constants.ind_start_Theta];  // Θ_ℓ(x,k)
  const double *Theta_p         = &y[Constants.ind_start_Thetap];
  const double *Nu              = &y[Constants.ind_start_nu];    

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdx      =  dydx[Constants.ind_deltac];        // d/dx[δ_c(x,k)]
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];        // d/dx[δ_b(x,k)]
  double &du_cdx          =  dydx[Constants.ind_uc];            // d/dx[u_c(x,k)]
  double &du_bdx          =  dydx[Constants.ind_ub];            // d/dx[u_b(x,k)]
  double &dPhidx          =  dydx[Constants.ind_Phi];           // d/dx[Φ(x,k)]
  double *dThetadx        = &dydx[Constants.ind_start_Theta];   // d/dx[Θ_ℓ(x,k)]
  double *dTheta_pdx      = &dydx[Constants.ind_start_Thetap];
  double *dNudx           = &dydx[Constants.ind_start_nu];

  // Cosmological parameters and variables 
  double Hp           = cosmo->Hp_of_x(x);          // Hp(x)
  double dHpdx        = cosmo->dHpdx_of_x(x);       // d/dx[Hp(x)]
  double H0           = cosmo->get_H0();            // H_0
  double Omegagamma0  = cosmo->get_Omegagamma();    // Ω_γ0
  double Omegac0      = cosmo->get_Omegac();        // Ω_c0
  double Omegab0      = cosmo->get_Omegab();        // Ω_b0

  // Recombination variables
  double R          = rec->R_of_x(x);           // R(x)
  double dtaudx     = rec->dtaudx_of_x(x);      // d/dx[τ(x)]
  double ddtaudxx   = rec->ddtaudxx_of_x(x);    // d^2/dx^2[τ(x)]

  


  //  useful factors
  double c      = Constants.c;    // c
  double ck     = c*k;            // ck
  double ckH    = ck/Hp;          // ck/Hp
  double Rplus  = (1+R);          // ( 1 + R )
  double Rminus = (1-R);          // ( 1 - R ) 
  double H0H0   = H0*H0;          // (H_0)^2
  double a      = exp(x);         // e^x

  double expr_a, expr_b;

  // Perturbation variables
  // double Psi = get_Psi(x, k);   // Ψ(x,k)
  double Psi = -Phi - 12*H0H0/(ck*ck*a*a) * Omegagamma0*Theta[2];  // Ψ(x,k) ???

  


  //  compute derivatives of scalar quantities (Φ, δ, u):

  expr_a = (Omegac0*delta_c + Omegab0*delta_b)*a + 4*Omegagamma0*Theta[0];
  dPhidx = Psi - 1./3*ckH*ckH*Phi + 0.5*H0H0/(Hp*Hp*a*a) * expr_a;
  
  expr_a = 0.5*ckH; expr_b = 3*dPhidx;
  ddelta_cdx = expr_a*u_c - expr_b;
  ddelta_bdx = expr_a*u_b - expr_b;
  
  expr_b = ckH*Psi;
  du_cdx = - u_c - expr_b;
  du_bdx = - u_b - expr_b + dtaudx * (3*Theta[0]+u_b) / R;


  //  compute derivatives of photon multipoles (Θ_ℓ):

  int ell     = 0;                // ℓ
  int ell_max = n_ell_Theta - 1;  // ℓ_max

  dThetadx[0] = -ckH * Theta[1] + dtaudx*Theta[0] - dPhidx;

  for(ell=1; ell<ell_max; ell++){
    expr_a = ckH / (2*ell+1) * ( ell*Theta[ell-1] - (ell+1)*ell*Theta[ell+1] );
    expr_b = dtaudx * Theta[ell];
    dThetadx[ell] = expr_a + expr_b;
  }
  dThetadx[1] += ckH*Psi + 1./3*dtaudx*u_b;
  dThetadx[2] -= 0.1*dtaudx*Theta[2];
  
  ell = ell_max;
  dThetadx[ell] = ckH*dThetadx[ell-1] - (c*(ell+1)/(Hp*cosmo->eta_of_x(x)) + dtaudx)*Theta[ell];


  return GSL_SUCCESS;
}





//====================================================
// Get methods
//====================================================



//  ----------------------------------
//  GET METHODS
//  ----------------------------------

double Perturbations::get_delta_c(const double x, const double k) const{
  return delta_c_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_u_c(const double x, const double k) const{
  return u_c_spline(x,k);
}
double Perturbations::get_u_b(const double x, const double k) const{
  return u_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltac:         " << Constants.ind_deltac         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_u_c:            " << Constants.ind_uc             << "\n";
  std::cout << "ind_u_b:            " << Constants.ind_ub               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_Theta:    " << Constants.ind_start_Theta      << "\n";
  std::cout << "n_ell_Theta:        " << Constants.n_ell_Theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_Thetap:   " << Constants.ind_start_Thetap   << "\n";
    std::cout << "n_ell_Thetap:       " << Constants.n_ell_Thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltac:         " << Constants.ind_deltac_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_u_c:            " << Constants.ind_uc_tc          << "\n";
  std::cout << "ind_u_b:            " << Constants.ind_ub_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_Theta:    " << Constants.ind_start_Theta_tc   << "\n";
  std::cout << "n_ell_Theta:        " << Constants.n_ell_Theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(OUTPUT_PATH + filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_Pi(x,k)        << " ";
    fp << get_Source_T(x,k)  << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}




/*

\Psi : Ψ
\Phi : Φ
\Theta : Θ
\Omega : Ω



\gamma : γ
\delta : δ
\eta : η
\tau : τ



\ell : ℓ

*/