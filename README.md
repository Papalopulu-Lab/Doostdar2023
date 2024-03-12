# Doostdar2023
## Mathematical modelling code for the paper titled: _Cell coupling compensates for changes in single-cell Her6 dynamics and provides phenotypic robustness_
**Authors:** Parnian Doostdar<sup>1</sup>, Joshua Hawley<sup>1</sup>, Kunal Chopra, Elli Marinopoulou, Robert Lea, Kiana Arashvand, Veronica Biga*, Nancy Papalopulu* and Ximena Soto Rodriguez*

Division of Developmental Biology and Medicine, School of Medical Sciences, Faculty of Biology, Medicine and Health, The University of Manchester, Oxford Road, Manchester, M13 9PT, UK

**Joint corresponding authors,
<sup>1</sup> joint first authors*

## How to use the code
All code is written in MATLAB. The purpose of the code is to run an optimiser on a multicellular mathematical model of coupled Her6 expression dynamics. Download this folder and run one of the 4 scripts contained within this project. The function of each script is outlined below.

## What each script does
### runOptimiser.m
Runs a pattern search optimiser N-times starting from random initial parameter values, and returns a minimised error parameter set on each run. It optimises on a multicellular model of Her6 expression and minimises for maximal changes in population coefficient of variation when the protein degradation rate is increased.

Toolboxes required to run this script: *Optimization toolbox, Signal toolbox, Statistics toolbox, Wavelet toolbox*

**Note that the optimiser parameter sets identified for two model setups (Model 1 and Model 2 in the paper) have been saved in the Data folder. Therefore you do not need to run this to use the other scripts in this repository because by default they use the workspaces that have already been saved in the Data folder.**

### plotOptimiserResults
This takes workspaces saved in the Data folder which contains the results of optimiser results and analyses and plots the data.

Important variables to check before running: _Boundary, rows, cols, Stochastic, CoupledCells_; do they match the data being used from the optimiser?

### runSingleParameterSet.m
Uses a single parameter set from the optimiser-identified parameter sets in the Data folder, and plots the outputs of a single run of the model.


### plotDifferentDegradationRates.m 
Uses one parameter set from the optimiser parameter sets in the Data folder and simulates the same parameter set over several different protein degradation rates, and then plots these outputs as histograms of Her6 expression distributions.
