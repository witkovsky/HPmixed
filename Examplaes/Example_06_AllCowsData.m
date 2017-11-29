%% EAXAMPLE: SAS Animal Data / All Cows in Canada
%  Data analysis by Linear Mixed Model with MATLAB function HPMIXED

% Author: Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 06-Jan-2014 13:52:11

%% Data Description
% Modified version of the SAS Proc HPMIXED Example 43.1 Ranking Many
% Random-Effect Coefficients. See Wang and Tobias (2009).
%
% In analyzing models with random effects that have many levels, a frequent
% goal is to estimate and rank the predicted values of the coefficients
% corresponding to these levels. For example, in mixed models for animal
% breeding, the predicted coefficient of the random effect for each animal
% is referred to as the estimated breeding value (EBV) and animals with
% relatively high EBVs are chosen for breeding. This example demonstrates
% the use of the HPMIXED procedure for computing EBVs and their precision.
% The data for this problem were generated by simulation. Suppose you are
% considering analyzing EBVs for animals on 100 farms, with about 100
% animals of 5 different species on each farm. The AllCowsData dataset 
% simulates data with this structure, where about 40 observations of the
% response variable Yield are made per animal.
%
% Comment from  Wang and Tobias (2009):
% You can also use PROC MIXED and PROC GLIMMIX to compute EBVs, but the
% performance of these general mixed modeling procedures for this
% specialized kind of data and model is quite different from that of PROC
% HPMIXED. Here is the difference: whereas using a desktop-class PC to fit
% this model would take PROC MIXED or PROC GLIMMIX hours (assuming
% sufficient memory is available), the HPMIXED procedure runs in under a
% minute using relatively little memory. The MIXED and GLIMMIX procedures
% are engineered to have good performance properties across a broad class
% of models and analyses, a class much broader than what PROC HPMIXED can
% handle. On the other hand, The HPMIXED procedure can have better
% performance, in terms of both memory and run time, for certain
% specialized models and analyses, of which the current example is one. The
% next section explores this comparison between performance for these three
% procedures in more depth.
%
% REFERENCE
% Tianlin Wang and Randy Tobias: All the Cows in Canada: Massive Mixed
% Modeling with the HPMIXED Procedure in SAS� 9.2. Paper 256-2009, SAS
% Global Forum 2009: Statistics and Data Analysis, SAS Institute Inc., Cary
% NC.

%% Load dataset AllCowsData
clear
load dsAllCowsData

%% Create the model structure manually
% Used formula
formula  = 'Yield ~ Species + Species:Farm + (1 | Animal)';

% Response Variable: y = Yield
y = AllCowsData.Yield;
n = length(y);

% FE design matrix: X. Intercept + Species + Farm:Species
% Sort the levels and take their indices
[SpeciesName,~,jS] = unique(AllCowsData.Species);
[FarmName,~,jF] = unique(AllCowsData.Farm);
[SpeciesFarmName,~,jFS] = unique(AllCowsData.Farm.*AllCowsData.Species);
pS = size(SpeciesName,1);
pF = size(FarmName,1);

% Create the full-ranked design matrix X
rows = (1:n)';
cols = 1;
X = sparse(rows,1,1,n,cols);
cols = pS-1;
X(rows,1+(1:cols)) = sparse(rows(jS~=1),(jS(jS~=1)-1),1,n,cols);
cols = pS*(pF-1);
X(rows,pS+(1:cols)) = sparse(rows(jF~=1),(jFS(jF~=1)-pS),1,n,cols);

% RE design matrix: Z
[AnimalName,~,j] = unique(AllCowsData.Animal);
q = size(AnimalName,1);
Z = sparse(1:n,j,1,n,q);

% Create the model structure
model.Formula.char = formula;
model.y = y;
model.X = X;
model.Z = Z;
model.dimRE = q;
model.Description = 'SAS Animal Data / All Cows in Canada';

%% Fit the linear mixed model by HPMIXED with limitted  output
opts.verbose = false;
opts.ddfMethod = 'residual';
% opts.ddfMethod = 'Satterthwaite';
% opts.isExactGrad = 0;
opts.isExactFIreml = 0;
lmefit = hpmixed(model,opts);

disp(lmefit)
disp(lmefit.ModelInfo)
disp(lmefit.fixedEffects.Statistics.TABLE(1:20,:))
disp(lmefit.randomEffects.Statistics.TABLE(1:20,:))
disp(lmefit.varianceComponents.Statistics.TABLE)

%% The t-statistics for the random effects (Satterthwaite's DF)
tic;
STAT = getStats(Lambda,lmefit);
toc

disp(STAT.TABLE(1:20,:));

%% Sorted STAT TABLE
[~,id] = sort(STAT.TABLE.Lower,'descend');

disp(STAT.TABLE(id(1:20),:));

%% Fit the linear mixed model by FITLME (! Long execution time !)
% tic;
% lme = fitlme(AllCowsData,formula,'FitMethod','REML');
% toc
% 
% disp(lme)