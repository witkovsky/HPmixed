%% EAXAMPLE: Split-Plot Data
%  Data analysis by Linear Mixed Model with MATLAB function HPMIXED

% Author: Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 06-Jan-2014 13:52:11

%% Data Description
% Consider the following data from Stroup (1989a), which arise from a
% balanced split-plot design with the whole plots arranged in a randomized
% complete-block design. The variable A is the whole-plot factor, and the
% variable B is the subplot factor. A traditional analysis of these data
% involves the construction of the whole-plot error (A*Block) to test A and
% the pooled residual error (B*Block and A*B*Block) to test B and A*B. 
%
% Documentation reproduced from SAS/STAT(R) 9.22 User's Guide, 
% Example 56.1 Split-Plot Design.

%% Load dataset SplitPlotData 
clear
load dsSplitPlotData

%% Create the model design matrices 
formula  = 'y ~ A + B + A:B + (1 | Block) + (1 | Block:A)';
opts.dummyVarCode = 'reference';
%opts.dummyVarCode = 'effects';
%opts.dummyVarCode = 'full';
model = hpmixedmodel(SplitPlotData,formula,opts);
model.Description = 'SplitPlotData: SAS Example';


%% Fit the linear mixed model by HPMIXED with complete output
%opts.sigma20 = [62.3958   15.3819+25    9.3611]';
%opts.FitMethod = 'None';
opts.FitMethod = 'REML';
opts.verbose = true;
lmefit = hpmixed(model,opts);

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
title('SplitPlot Data Fitted by HPMIXED')

%% EXAMPLE 2: (t-statistics for the LSMEANS CONTRASTS for interaction A:B)

options.STAT.inference = 'contrasts';
options.STAT.inferenceSpace = 'broad';
options.STAT.ddfMethod = 'Satterthwaite';
%options.STAT.ddfMethod = 'None';

[Lambda,options]  = getLambda({'A' 'B' },model,SplitPlotData,options);

STAT_AB_CONTRASTS = getStats(Lambda,lmefit,options);
disp(STAT_AB_CONTRASTS)

%% EXAMPLE 3: (FTEST for the LSMEANS CONTRASTS for interaction A:B)

options.STAT.inference = 'contrasts';
options.STAT.inferenceSpace = 'broad';
options.STAT.ddfMethod = 'FaiCornelius';

[Lambda,options]  = getLambda({'A' 'B' },model,SplitPlotData,options);

STAT_AB_FTEST = getStats(Lambda,lmefit,options);
disp(STAT_AB_FTEST)

%% EXAMPLE 4: (LSMEANS PAIRWISE COMPARISONS for the effect A)

options.STAT.inference = 'pairwise';
options.STAT.inferenceSpace = 'broad';
options.STAT.ddfMethod = 'Satterthwaite';

[Lambda,options]  = getLambda({'A'},model,SplitPlotData,options);

STAT_A_COMPARISONS = getStats(Lambda,lmefit,options);
disp(STAT_A_COMPARISONS)

%% EXAMPLE 5: (Statistics for the LSMEANS for the effect A)

options.STAT.inference = 'lsmeans';
options.STAT.inferenceSpace = 'broad';

[Lambda,options]  = getLambda({'A'},model,SplitPlotData,options);

STAT_A_LSMEANS = getStats(Lambda,lmefit,options);
disp(STAT_A_LSMEANS)

%% EXAMPLE 5A: (Statistics for the LSMEANS for the effect A)

options.STAT.inference = 'contrasts';
options.STAT.inferenceSpace = 'broad';
options.STAT.ddfMethod = 'faicornelius';
[Lambda,options]  = getLambda({'A'},model,SplitPlotData,options);

STAT_A_LSMEANS = getStats(Lambda,lmefit,options);
disp(STAT_A_LSMEANS)

%% EXAMPLE_6: (BROAD INFERENCE SPACE for LSMEANS of the factor Block)

options.STAT.inference = 'lsmeans';
options.STAT.inferenceSpace = 'broad';
[Lambda,options]  = getLambda({'Block'},model,SplitPlotData,options);
options.STAT.alpha = 0.01;
STAT_Block_LSMEANS_BROAD = getStats(Lambda,lmefit,options);
disp(STAT_Block_LSMEANS_BROAD)

%% EXAMPLE 7: (NARROW INFERENCE SPACE for LSMEANS of the factor Block)

options.STAT.inference = 'lsmeans';
options.STAT.inferenceSpace = 'narrow';
options.STAT.alpha = 0.01;

[Lambda,options]  = getLambda({'Block'},model,SplitPlotData,options);

STAT_Block_LSMEANS_NARROW = getStats(Lambda,lmefit,options);
disp(STAT_Block_LSMEANS_NARROW)

%% EXAMPLE 8: (ITERMEDIATE INFERENCE SPACE for LSMEANS of the factor Block)

options.STAT.inference = 'lsmeans';
options.STAT.inferenceSpace = 'intermediate';
options.STAT.includedXcols = (1:size(model.X,2));
options.STAT.includedZcols = model.dimRE(1) + (1:model.dimRE(2));
options.STAT.alpha = 0.01;

[Lambda,options]  = getLambda({'Block'},model,SplitPlotData,options);

STAT_Block_LSMEANS_ITERMEDIATE = getStats(Lambda,lmefit,options);
disp(STAT_Block_LSMEANS_ITERMEDIATE)

%% Comparison with MATLAB function FITLME
tic;
lme = fitlme(SplitPlotData,formula,'FitMethod','REML');
toc
disp(lme)
%% Fit the linear mixed model by FITLME
tic;
lme = fitlme(SplitPlotData,formula,'DummyVarCoding','effects','FitMethod','REML');
toc

disp(lme)

anova(lme,'DFMethod','Satterthwaite')
%% Fit the linear mixed model by MIXED
dim = [4 12];
s20 = [1 1 1];
method = 2;
tic;
[s2,b,u,Is2,C,~,~,loglik] = mixed(model.y,model.X,model.Z,dim,s20,method);
toc
disp(s2)
%% Fit the linear mixed model by MIXEDHD
% dim = [4 12];
% s20 = [1 1 1];
% method = 2;
% lambda = 0.5;
% tic;
% [s2,b,u,Is2,C,~,~,loglik] = mixedHD(model.y,model.X,model.Z,dim,s20,method,lambda);
% toc
% disp(s2)