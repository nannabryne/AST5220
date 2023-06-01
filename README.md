# Microwave music: Finding the cosmic symphony

#### Project in AST5220 Cosmology II (autumn 2023)

The project report is found [here](https://github.com/nannabryne/AST5220/blob/main/tex/main.pdf), or as `tex/main.pdf`. Final draft was finished June 1st.


### ABSTRACT
> **Context.** Project in the master-level subject AST5220 Cosmology II at the University of Oslo (UiO).
> **Aims.** Write a self-contained program to compute theoretical CMB and matter power spectra, starting from fundamental principles.
> **Methods.** By perturbing the background and considering the recombination history, we obtained the necessary quantities for calculating the CMB and matter power spectra. To get there, we used the Boltzmann and Einstein equations and linear perturbation theory,w
amongst other things. We considered the conformal Newtonian gauge, and we neglected effects from neutrinos, helium, reionisation
and polarisation. We used numerical methods s.a. Runge-Kutta and the trapezoidal rule for integrating.
> **Results.** Our code provided satisfactory results to the point where neglected ingredients s.a. neutrinos were bound to play a role. All
results prior to the power spectra essentially substantiate the physics we drew from the latter.
> **Conclusions.** We succeded in building a code for our simple ΛCDM model with results that can be used to understand the governing
physics at large and intermediate scales. However, small-scale physics can be understood better by relatively simple generalisations
and improvements to the code.

### DATA AND FIGURES
The observational data we used are located in `data/input/`. Raw output is found in `data/output/`. Figures are found in `figs/milestone*/`.

## Milestone I: Background cosmology 
*(first draft 05/03/23)*
Relevant text found in Sect. 2 of the report. However, a [copy](https://github.com/nannabryne/AST5220/blob/main/deliverables/report_050323.pdf) of the paper as it was prior to moving on to the next milestone is found as `deliverables/report_050323.pdf`.

### _Products_
> FIGURES: `figs/milestone1/*`

> CODES: `src/backgroundcosmology.*`


## Milestone II: Recombination history
*(second draft 02/04/23)*
Relevant text found in Sect. 3 of the report. However, a [copy](https://github.com/nannabryne/AST5220/blob/main/deliverables/report_020423.pdf) of the paper as it was prior to moving on to the next milestone is found as `deliverables/report_020423.pdf`.

### _Products_

> FIGURES: `figs/milestone2/*`

> CODES: `src/recombinationhistory.*`



## Milestone III: Growth of structure
*(third draft 03/05/23)*
Relevant text found in Sect. 4 of the report. However, a [copy](https://github.com/nannabryne/AST5220/blob/main/deliverables/report_030523.pdf) of the paper as it was prior to moving on to the next milestone is found as `deliverables/report_030523.pdf`. 

### _Products_
> FIGURES: `figs/milestone3/*`

> CODES: `src/perturbations.*`



## Milestone IV: CMB and matter power spectra
*(fourth draft 27/05/23)*
Relevant text found in Sect. 5 of the report. However, a [copy](https://github.com/nannabryne/AST5220/blob/main/deliverables/report_270523.pdf) of the paper as it was prior to revising it one last time is found as `deliverables/report_270523.pdf`. 

### _Products_
> FIGURES: `figs/milestone4/*`

> CODES: `src/powerspectrum.*`




# Code
All code is located in the `src/` directory.
## Prerequisites
For running code in `C++`, the following is needed:
- `GSL`
- `omp`

For running `Python` code, one needs this:
- `numpy`
- `scipy`
- `pandas`
- `matplotlib`
- `seaborn` 

## How to run
Go to the `src/` directory. Create a build directory unless you have not already, by typing:
```
make dir
```
To run all simulations (i.e. run `src/main.cpp`), type:
```
make cmb
```

To produce plots, run:
~~~
make plotsX     # X is the milestone number
~~~


## Plotting source
We produce plots using `Python` and relevant code is found in `src/plotting/`. The file `analysis.py` contains all the data reading and function calls, but plotting scripts are found as `plotsrc_milestone_*.py`.

<!-- # <mark>TO DO:
## Makefile
- [ ] Create a makefile in parent directory (this)
- [x] Fix precompiled headers!
## Plotting
- [x] Clean up plotting code: get sensible structure etc.
- [ ] Decide on style: time to give up on seaborn darkgrid (I think I am ready for a more professional style? Mom come pick me up I am scared)
- [x] Fix units in plots (km/s $\to$ km s $^{-1}$ )
- [ ] Fix colours in density parameters plot
- [ ] Coordinate fonts to be the same as in tex document
- [x] Fix labels in H, dHdx, ...-plots
## Coding
- [x] Structural changes
- [x] (M1) What happened to the luminosity distance????
- [x] Review after feedback from M1
- [x] Review after feedback from M2
- [x] Review after feedback from M3
## Report
- [x] Write figure captions
- [x] Review after feedback from M1
- [x] Review after feedback from M2
- [ ] Review after feedback from M3
## Readme
- [x] Complete list of prerequisites 
- [ ] Structural changes -->
