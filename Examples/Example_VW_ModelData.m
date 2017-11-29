%% EAXAMPLE: Model Data
%  Data analysis by Linear Mixed Model with MATLAB function HPMIXED

% Author: Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 28-Oct-2014 00:25:56

%% Data Description
% Test model data used to test covergence of lmer procedure in R.
% Minimal working example for different convergences in lme4 and lme4.0

%% Load dataset SplitPlotData 
clear
load dsModelData
%% Create the model design matrices 
%formula  = 'mean1 ~ ambiguity * sdiff + (1 | item) + (-1 + sdiff | item) + (1 | subj) + (-1 + sdiff | subj)';
formula  = 'mean1 ~ ambiguity + sdiff + ambiguity:sdiff + (1 | item) + (-1 + sdiff | item) + (1 | subj) + (-1 + sdiff | subj)';
opts.dummyVarCode = 'reference';
%opts.dummyVarCode = 'effects';
%opts.dummyVarCode = 'full';
model = hpmixedmodel(ModelData,formula,opts);
model.Description = 'ModelData: R Test Example';

%% Fit the linear mixed model by HPMIXED with complete output
%opts.FitMethod = 'None';
opts.FitMethod = 'REML';
opts.verbose = true;
lmefit = hpmixed(model,opts);

disp(lmefit)

%% ANOVA (Type III SS if opts.dummyVarCode = 'effects')

STAT_ANOVA = getAnova(lmefit);
disp(STAT_ANOVA)

%% Statistics for FIXED and RANDOM effects 

STAT_FE = getStats('fixed',lmefit);
disp(STAT_FE)

STAT_RE = getStats('random',lmefit);
disp(STAT_RE)

%% Statistics for FITTED values
STAT_FIT = getStats('fitted',lmefit);
disp(STAT_FIT)

figure
plot(model.y,STAT_FIT.TABLE.Estimate,'o')
grid
xlabel('y observed')
ylabel('y fitted')
title('SplitPlot Data Fitted by HPMIXED')