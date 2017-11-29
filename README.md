# HPmixed
HPMIXED fits the linear mixed effects models (LME models), y = X*b + Z*u + e, with a simple variance covariance structure (VC: variance components), by solving the Henderson's mixed model equations. 

For current status of the MATLAB toolbox see the HPmixed development available at

- https://github.com/witkovsky/HPmixed

About
=====

HPmixed (the High Performance Mixed Effects Model toolbox) consists of a set of algorithms for fitting the linear mixed models of the form  y = X*b + Z*u + e with a simple variance componets structure by solving the Henderson's mixed model equations. 
                                                                              
The model structure can be specified by the data vector y and the design matrices X and Z. Here, X is assumed to be a full ranked matrix and Z is typically a sparse matrix. 

Alternatively, the model can be specified from the given DATASET and FORMULA by using the developed function hpmixedmodel, which is based on functionality of the LinearMixedModels class (Statistics and Machine Learning Toolbox, MATLAB).

Installation and requirements
=============================

HPmixed was developed with MATLAB Version: 9.2 (R2017a).

To install, you can either clone the directory with Git or download a .zip file. 

## Option 1: Download .zip file

Download a .zip of HPmixed from

- https://github.com/witkovsky/HPmixed/archive/master.zip

After unzipping, you will need to add HPmixed to the MATLAB path. You can do this either (a) by typing
```
addpath(HPmixedRoot), savepath
```
where `HPmixedRoot` is the path to the unzipped directory, (b) by selecting the `HPmixed` directory with the `pathtool` command, or (c) though the File > Set Path... dialog from the MATLAB menubar.

## Option 2: Clone with Git

To clone the CharFunTool repository, first navigate in a terminal to where you want the repository cloned, then type
```
git clone https://github.com/witkovsky/HPmixed.git
```
To use HPmixed in MATLAB, you will need to add the `HPmixed` directory to the MATLAB path as above.


Getting started
===============

We recommend taking a look at the Examples collection. 

To get a taste of what computing with HPmixed is like, try to fit LME to the standard Split-Plot data, see e.g. Stroup (1989). For that, simply type
```
clear
load dsSplitPlotData
formula = 'y ~ A + B + A:B + (1 | Block) + (1 | Block:A)';
model   = hpmixedmodel(SplitPlotData,formula);
lmefit  = hpmixed(model) 
```

HPmixed allows to get more advanced inference from fitted LME. For example, here we illustrate one possible way how to evaluate the LSMEANS contrasts from the Split-Plot data in a broad sense inference space, for the interaction A:B.
```
options.STAT.inference = 'contrasts';
options.STAT.inferenceSpace = 'broad';
[Lambda,options]  = getLambda({'A' 'B' },model,SplitPlotData,options);

options.verbose = false;
STAT   = getStats(Lambda,lmefit,options)
```

HPmixed is efficient for large LME models. For example, consider the modified version of the SAS Proc HPMIXED Example 43.1 Ranking Many Random-Effect Coefficients, see Wang and Tobias (2009).
```
clear
load dsAnimalData
formula = 'Yield ~ Species + (1 | Farm) + (1 | Animal)';
model   = hpmixedmodel(SASAnimalData,formula);

options.ddfMethod = 'Satterthwaite';
lmefit = hpmixed(model,options);

disp(lmefit)

```

License
=======

See `LICENSE.txt` for HPmixed licensing information.
