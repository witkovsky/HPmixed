%% EAXAMPLE: SAS Animal Data / All Cows in Canada
%  Data analysis by Linear Mixed Model with MATLAB function HPMIXED

% Author: Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 06-Jan-2014 13:52:11

%% Data Description
% SAS Proc HPMIXED Example 43.1 Ranking Many Random-Effect Coefficients. 
% In analyzing models with random effects that have many levels, a frequent
% goal is to estimate and rank the predicted values of the coefficients
% corresponding to these levels. For example, in mixed models for animal
% breeding, the predicted coefficient of the random effect for each animal
% is referred to as the estimated breeding value (EBV) and animals with
% relatively high EBVs are chosen for breeding. This example demonstrates
% the use of the HPMIXED procedure for computing EBVs and their precision.
% The data for this problem were generated by simulation. Suppose you are
% considering analyzing EBVs for animals on 100 farms, with about 100
% animals of 5 different species on each farm. The following DATA step
% simulates data with this structure, where about 40 observations of the
% response variable Yield are made per animal.

%% Load dataset SASAnimalData
clear
load dsAnimalData

%% Create the model structure for HPMIXED by using hpmixedmodel
% If MATLAB 2013b + Statistics Toolbox is available
% Otherwise use provided model structure, or construct it manually

formula  = 'Yield ~ Species + (1 | Farm) + (1 | Animal)';
tic;
model = hpmixedmodel(SASAnimalData,formula);
toc

%% Fit the linear mixed model by HPMIXED with limitted  output
options.ddfMethod = 'residual';
lmefit1 = hpmixed(model,options);

disp(lmefit1)

%% Fit the linear mixed model by HPMIXED with complete output
% !! Computing ALL the statistics take considerable longer time
%
opts.verbose = true;
opts.ddfMethod = 'Satterthwaite';
lmefit = hpmixed(model,opts);

disp(lmefit)

%% Fit the linear mixed model by FITLME
tic;
lme = fitlme(SASAnimalData,formula,'FitMethod','REML');
toc

disp(lme)
