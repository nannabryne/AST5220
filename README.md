# [SOME_TITLE]
### Project in AST5220 Cosmology II (autumn 2023)
The project report is found [here](https://github.com/nannabryne/AST5220/blob/main/tex/main.pdf), or as `tex/main.pdf`.

## Milestone I: Background cosmology 
*(first draft 05/03/23)*
Relevant text found in Sec. 2 of the report. However, a [copy](https://github.com/nannabryne/AST5220/blob/main/deliverables/report_050323.pdf) of the paper as it was prior to moving on to the next milestone is found as `deliverables/report_050323.pdf`.

# Code
All code is located in the `src/` directory.
## Prerequisites <mark>INCOMPLETE</mark>
For running code in `C++`, the following is needed:
- `GSL`
For running `Python` code, one needs this:
- `numpy`
- `matplotlib`
- `seaborn` 

## How to run
Go to the `src/` directory. To run all simulations, type:
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
make plots
~~~

## Plotting source
We produce plots using `Python` and relevant code is found in `src/plotting/`. 

<mark>**NB**:</mark> This is a work in progress. Currently (03/03/23) it is complete chaos.

# <mark>TO DO:
## Makefile
- [ ] Create a makefile in parent directory (this)
## Plotting
- [ ] Clean up plotting code: get sensible structure etc.
- [ ] Decide on style: time to give up on seaborn darkgrid (I think I am ready for a more professional style? Mom come pick me up I am scared)
- [x] Fix units in plots (km/s $\to$ km s $^{-1}$ )
- [ ] Fix colours in density parameters plot
- [ ] Coordinate fonts to be the same as in tex document
- [ ] Fix labels in H, dHdx, ...-plots
## Coding
- [ ] Structural changes
- [x] (M1) What happened to the luminosity distance????
## Report
- [x] Write figure captions
## Readme
- [ ] Complete list of prerequisites 
- [ ] Structural changes
