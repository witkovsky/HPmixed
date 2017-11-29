%% EAXAMPLE: Split-Plot Data / Type III sum of squares
%  Data analysis by Linear Mixed Model with MATLAB function HPMIXED
%  Use MATLAB algorithm (package) HPMIXED:
%  http://www.mathworks.com/matlabcentral/fileexchange/45576-hpmixed

% Author: Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 21-Oct-2014 13:53:35

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

%% Create the model design matrices  / 
%  Use opts.dummyVarCode = 'effects' for Type III SS 
%  Otherwise, we suggest to use opts.dummyVarCode = 'reference', which is
%  in general computationally more efficient
%
formula  = 'y ~ A + B + A:B + (1 | Block) + (1 | Block:A)';
%opts.dummyVarCode = 'full';
%opts.dummyVarCode = 'reference';
opts.dummyVarCode = 'effects';
model = hpmixedmodel(SplitPlotData,formula,opts);
model.Description = 'SplitPlotData: SAS Example';

%% Fit the linear mixed model by HPMIXED with complete output
opts.FitMethod = 'REML';
opts.verbose = true;
lmefit = hpmixed(model,opts);

disp(lmefit)
%% ANOVA Table with Type III SS

STAT_ANOVA = getAnova(lmefit);
disp(STAT_ANOVA)