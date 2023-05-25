# [SOME_TITLE]
### Project in AST5220 Cosmology II (autumn 2023)
The project report is found [here](https://github.com/nannabryne/AST5220/blob/main/tex/main.pdf), or as `tex/main.pdf`.

## Milestone I: Background cosmology 
*(first draft 05/03/23)*
Relevant text found in Sect. 2 of the report. However, a [copy](https://github.com/nannabryne/AST5220/blob/main/deliverables/report_050323.pdf) of the paper as it was prior to moving on to the next milestone is found as `deliverables/report_050323.pdf`.


## Milestone II: Recombination history
*(second draft 02/04/23)*
Relevant text found in Sect. 3 of the report. However, a [copy](https://github.com/nannabryne/AST5220/blob/main/deliverables/report_020423.pdf) of the paper as it was prior to moving on to the next milestone is found as `deliverables/report_020423.pdf`.


## Milestone III: Growth of structure
*(third draft 03/05/23)*
Relevant text found in Sect. 4 of the report. However, a [copy](https://github.com/nannabryne/AST5220/blob/main/deliverables/report_030523.pdf) of the paper as it was prior to moving on to the next milestone is found as `deliverables/report_030523.pdf`.



# Code
All code is located in the `src/` directory.
## Prerequisites <mark>INCOMPLETE</mark>
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
To run all simulations, type:
```
make clean
make all
make run
```
An alias for these commands is:
~~~
make autorun
~~~
To produce plots, run:
~~~
make plotsX     # X is the milestone number
~~~


## Plotting source
We produce plots using `Python` and relevant code is found in `src/plotting/`. The file `analysis.py` ...

# <mark>TO DO:
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
- [ ] Structural changes
- [x] (M1) What happened to the luminosity distance????
- [x] Review after feedback from M1
- [x] Review after feedback from M2
- [ ] Review after feedback from M3
## Report
- [x] Write figure captions
- [ ] Review after feedback from M1
- [ ] Review after feedback from M2
- [ ] Review after feedback from M3
## Readme
- [ ] Complete list of prerequisites 
- [ ] Structural changes
