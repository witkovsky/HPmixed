%% EAXAMPLE: Ham Data
%  Data analysis by Linear Mixed Model with MATLAB function HPMIXED

% Author: Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 06-Jan-2014 13:52:11

%% Data Description
% Conjoint study of dry cured ham. One of the purposes of the study was to
% investigate the effect of information given to the consumers measured in
% hedonic liking for the hams. Two of the hams were Spanish and two were
% Norwegian, each origin representing different salt levels and different
% aging time. The information about origin was given in such way that both
% true and false information was given. essentially a 4*2 design with 4
% samples and 2 information levels. A total of 81 Consumers participated in
% the study.
% 
% Consumer factor with 81 levels: 
% * Consumer -  identifying consumers
% * Product - factor with four levels
% * Informedliking - numeric: hedonic liking for the products
% * Information - factor with two levels
% * Gender - factor with two levels (gender)
% * Age - numeric: age of Consumer
% 
% References
% 
% "Alternative methods for combining design variables and consumer
% preference with information about attitudes and demographics in conjoint
% analysis". T. Naes, V.Lengard, S. Bolling Johansen, M. Hersleth.
% 
% R Example:
% m <- lmer(Informed.liking ~ Product*Information*Gender +
% (1|Product:Consumer) , data=ham)
%   
% Documentation reproduced from R package lmerTest,
% Version: 2.0-3, Date: 2012-01-09
% Authors: Author: Alexandra Kuznetsova, Per Bruun Brockhoff, Rune Haubo
% Bojesen Christensen

%% Load dataset HamData
clear
load dsHamData

%% Load dataset HamData and create the model design matrices 
% formula  = 'Informedliking ~ Gender + Information + Product + (1|Consumer)';
% formula  = 'Informedliking ~  Product*Information*Gender + (1|Product:Consumer)';
formula  = 'Informedliking ~  Product*Information*Gender + (1|Consumer)';
% formula  = 'Informedliking ~  Information*Gender + (1|Consumer)';
model = hpmixedmodel(HamData,formula);
model.Description = 'HamData: R lmerTest Example';

%% Fit the linear mixed model by HPMIXED
lmefit = hpmixed(model);

disp(lmefit)

%% EXAMPLE 1: (Statistics for FIXED and RANDOM effects and FITTED values)

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
title('Ham Data Fitted by HPMIXED')

%% EXAMPLE: (NARROW INFERENCE SPACE for LSMEANS of the factor Product)

options.STAT.inference = 'lsmeans';
options.STAT.inferenceSpace = 'narrow';

[Lambda,options]  = getLambda({'Product'},model,HamData,options);

STAT_Product_LSMEANS_NARROW = getStats(Lambda,lmefit,options);
disp(STAT_Product_LSMEANS_NARROW)

%% EXAMPLE: (Broad INFERENCE SPACE for LSMEANS of the factor Product)

options.STAT.inference = 'lsmeans';
options.STAT.inferenceSpace = 'broad';

[Lambda,options]  = getLambda({'Product'},model,HamData,options);

STAT_Product_LSMEANS_BROAD = getStats(Lambda,lmefit,options);
disp(STAT_Product_LSMEANS_BROAD)

%% EXAMPLE: (LSMEANS PAIRWISE COMPARISONS for the factor Product)

options.STAT.inference = 'pairwise';
options.STAT.inferenceSpace = 'narrow';
options.STAT.ddfMethod = 'Satterthwaite';

[Lambda,options]  = getLambda({'Product'},model,HamData,options);

STAT_Product_COMPARISONS = getStats(Lambda,lmefit,options);
disp(STAT_Product_COMPARISONS)

%% Fit the linear mixed model by FITLME
tic;
lme = fitlme(HamData,formula,'FitMethod','REML');
toc

disp(lme)
