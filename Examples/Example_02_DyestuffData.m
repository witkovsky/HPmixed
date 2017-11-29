%% EAXAMPLE: Dyestuff Data
%  Data analysis by Linear Mixed Model with MATLAB function HPMIXED

% Author: Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 06-Jan-2014 13:52:11

%% Data Description
% The Dyestuff data frame provides the yield of dyestuff (Naphthalene Black
% 12B) from 5 different preparations from each of 6 different batchs of an
% intermediate product (H-acid). 
%
% Data frames, each with 30 observations on the following 2 variables.
% * Batch - a factor indicating the batch of the intermediate product from
%   which the preparation was created.
% * Yield - the yield of dyestuff from the preparation (grams of standard
%   color).
% 
% The Dyestuff data are described in Davies and Goldsmith (1972) as coming
% from “an investigation to find out how much the variation from batch to
% batch in the quality of an intermediate product (H-acid) contributes to
% the variation in the yield of the dyestuff (Naphthalene Black 12B) made
% from it. In the experiment six samples of the intermediate, representing
% different batches of works manufacture, were obtained, and five
% preparations of the dyestuff were made in the laboratory from each
% sample. The equivalent yield of each preparation as grams of standard
% colour was determined by dye-trial.”
% 
% References
% 
% O.L. Davies and P.L. Goldsmith (eds), Statistical Methods in Research and
% Production, 4th ed., Oliver and Boyd, (1972), section 6.4
% 
% G.E.P. Box and G.C. Tiao, Bayesian Inference in Statistical Analysis,
% Addison-Wesley, (1973), section 5.1.2
% 
% R Example:
% fm1 <- lmer(Yield ~ 1|Batch, Dyestuff)
% 
% Documentation reproduced from R package lme4,
% Version: 1.0-5, Date: 2013-10-25.
% Authors: Douglas Bates, Martin Maechler, Ben Bolker, Steven Walker

%% Load dataset DyestuffData 
clear
load dsDyestuffData
%% Create the model design matrices 
formula  = 'Yield ~ 1 + (1 | Batch)';
%opts.dummyVarCode = 'reference';
opts.dummyVarCode = 'effects';
model = hpmixedmodel(DyestuffData,formula,opts);
model.Description = 'DyestuffData: Dyestuff R Example';

%% Fit the linear mixed model by HPMIXED with limitted  output
opts.verbose = false;
lmefit1 = hpmixed(model,opts);

disp(lmefit1)
%% Fit the linear mixed model by HPMIXED with limitted  output
lmefit = hpmixed(model);

disp(lmefit)
%% EXAMPLE 1: (Statistics for FIXED and RANDOM effects and FITTED values)
STAT_ANOVA = getAnova(lmefit);
disp(STAT_ANOVA)

STAT_FE = getStats('fixed',lmefit);
disp(STAT_FE)

STAT_RE = getStats('random',lmefit);
disp(STAT_RE)

STAT_FIT = getStats('fitted',lmefit);
disp(STAT_FIT)

figure
plot(model.y,STAT_FIT.TABLE.Estimate,'o')
grid
xlabel('y observed')
ylabel('y fitted')
title('Dyestuff Data Fitted by HPMIXED')
%% EXAMPLE (BROAD INFERENCE SPACE for LSMEANS of the factor Batch)

options.STAT.inference = 'lsmeans';
options.STAT.inferenceSpace = 'broad';
[Lambda,options]  = getLambda({'Batch'},model,DyestuffData,options);
options.STAT.alpha = 0.01;
STAT = getStats(Lambda,lmefit,options);
disp(STAT)
%% Fit the linear mixed model by FITLME
tic;
lme = fitlme(DyestuffData,formula,'FitMethod','REML');
toc

disp(lme)