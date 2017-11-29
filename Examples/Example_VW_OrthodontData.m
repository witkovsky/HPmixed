%% EAXAMPLE: Orthodont Data
%  Data analysis by Linear Mixed Model with MATLAB function HPMIXED

% Author: Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 17-Jan-2015 16:46:32

%% Data Description
% Investigators at the University of North Carolina Dental School followed
% the growth of 27 children (16 males, 11 females) from age 8 until age 14.
% Every two years they measured the distance between the pituitary and the
% pterygomaxillary fissure, two points that are easily identified on x-ray
% exposures of the side of the head.
% 
% Source
% Pinheiro, J. C. and Bates, D. M. (2000), Mixed-Effects Models in S and
% S-PLUS, Springer, New York. (Appendix A.17)
% 
% Potthoff, R. F. and Roy, S. N. (1964), “A generalized multivariate
% analysis of variance model useful especially for growth curve problems”,
% Biometrika, 51, 313–326.

% % SAS
% data pr;
%      input Person Gender $ y1 y2 y3 y4;
%      y=y1; Age=8;  output;
%      y=y2; Age=10; output;
%      y=y3; Age=12; output;
%      y=y4; Age=14; output;
%      drop y1-y4;
%      datalines;
%     1   F   21.0    20.0    21.5    23.0
%     2   F   21.0    21.5    24.0    25.5
%     3   F   20.5    24.0    24.5    26.0
%     4   F   23.5    24.5    25.0    26.5
%     5   F   21.5    23.0    22.5    23.5
%     6   F   20.0    21.0    21.0    22.5
%     7   F   21.5    22.5    23.0    25.0
%     8   F   23.0    23.0    23.5    24.0
%     9   F   20.0    21.0    22.0    21.5
%    10   F   16.5    19.0    19.0    19.5
%    11   F   24.5    25.0    28.0    28.0
%    12   M   26.0    25.0    29.0    31.0
%    13   M   21.5    22.5    23.0    26.5
%    14   M   23.0    22.5    24.0    27.5
%    15   M   25.5    27.5    26.5    27.0
%    16   M   20.0    23.5    22.5    26.0
%    17   M   24.5    25.5    27.0    28.5
%    18   M   22.0    22.0    24.5    26.5
%    19   M   24.0    21.5    24.5    25.5
%    20   M   23.0    20.5    31.0    26.0
%    21   M   27.5    28.0    31.0    31.5
%    22   M   23.0    23.0    23.5    25.0
%    23   M   21.5    23.5    24.0    28.0
%    24   M   17.0    24.5    26.0    29.5
%    25   M   22.5    25.5    25.5    26.0
%    26   M   23.0    24.5    26.0    30.0
%    27   M   22.0    21.5    23.5    25.0
%    ;
% 
% /* Model and Solution by SAS */
%    proc mixed data=pr method=reml asycov;
%       class Person Gender;
%       model y = Gender Age Gender*Age / ddfm =kr s;
% 	  random  intercept Age / subject=Person type=un; 
%    run;

%% Asymptotic covariance matrix of variance component estimators
% 
% Sigma = [26.3705   -2.0160    0.1561   -0.6872
%          -2.0160    0.1724   -0.0147    0.0600
%           0.1561   -0.0147    0.0014   -0.0055
%          -0.6872    0.0600   -0.0055    0.1091 ];

%% Load dataset RatData 
clear
load dsOrthodontData

%% Create the model design matrices 
%formula  = 'distance ~  age + Sex + (1 | Subject) + (age-1 | Subject)';
formula  = 'distance ~ age*Sex + (1 | Subject) + (age-1 | Subject)';
%formula  = 'distance ~ age*Sex + (age | Subject)';
%opts.dummyVarCode = 'reference';
opts.dummyVarCode = 'effects';
%opts.dummyVarCode = 'full';
model = hpmixedmodel(OrthodontData,formula,opts);
model.Description = 'OrthodontData: Data set used in Potthoff and Roy (1964)';

%% Fit the linear mixed model by HPMIXED with complete output
opts.FitMethod = 'REML';
opts.verbose = true;
lmefit = hpmixed(model,opts);

disp(lmefit)

%% Statistics for FIXED and RANDOM effects and FITTED values

STAT_ANOVA = getAnova(lmefit);
disp(STAT_ANOVA)

STAT_FE = getStats('fixed',lmefit);
disp(STAT_FE)

STAT_RE = getStats('random',lmefit);
disp(STAT_RE)

STAT_FIT = getStats('fitted',lmefit);
disp(STAT_FIT)

%% Narrow inference
% Compute statistitc for mean gender differences at different ages

options.STAT.inference = 'contrasts';
options.STAT.inferenceSpace = 'narrow';

ds8 = OrthodontData(OrthodontData.age == 8,:);
model8 = hpmixedmodel(ds8,formula,opts);
L8  = getLambda({ 'Sex' },model8,ds8,options);

ds10 = OrthodontData(OrthodontData.age == 10,:);
model10 = hpmixedmodel(ds10,formula,opts);
L10  = getLambda({ 'Sex' },model10,ds10,options);

ds12 = OrthodontData(OrthodontData.age == 12,:);
model12 = hpmixedmodel(ds12,formula,opts);
L12  = getLambda({ 'Sex' },model12,ds12,options);

ds14 = OrthodontData(OrthodontData.age == 14,:);
model14 = hpmixedmodel(ds14,formula,opts);
L14  = getLambda({ 'Sex' },model14,ds14,options);

Lambda = [L8 L10 L12 L14];
%options.STAT.ddfMethod = 'fc';
options.STAT.colnames = {'Females-Males_Age08'; 'Females-Males_Age10' ;...
    'Females-Males_Age12' ;'Females-Males_Age14' };
%
result = getStats(Lambda,lmefit,options);

%% Broad inference
options.STAT.inference = 'contrasts';
options.STAT.inferenceSpace = 'broad';

ds8 = OrthodontData(OrthodontData.age == 8,:);
model8 = hpmixedmodel(ds8,formula,opts);
L8  = getLambda({ 'Sex' },model8,ds8,options);

ds10 = OrthodontData(OrthodontData.age == 10,:);
model10 = hpmixedmodel(ds10,formula,opts);
L10  = getLambda({ 'Sex' },model10,ds10,options);

ds12 = OrthodontData(OrthodontData.age == 12,:);
model12 = hpmixedmodel(ds12,formula,opts);
L12  = getLambda({ 'Sex' },model12,ds12,options);

ds14 = OrthodontData(OrthodontData.age == 14,:);
model14 = hpmixedmodel(ds14,formula,opts);
L14  = getLambda({ 'Sex' },model14,ds14,options);

LambdaBroad = [L8 L10 L12 L14];
%options.STAT.ddfMethod = 'fc';
options.STAT.colnames = {'Females-Males_Age08'; 'Females-Males_Age10' ;...
    'Females-Males_Age12' ;'Females-Males_Age14' };

resultBroad = getStats(LambdaBroad,lmefit,options);

%% Alternatively

LambdaBroad = Lambda;
LambdaBroad(1:end-4,:) = 0;

%options.STAT.ddfMethod = 'fc';
options.STAT.colnames = {'Females-Males_Age08'; 'Females-Males_Age10' ;...
    'Females-Males_Age12' ;'Females-Males_Age14' };
resultBroad = getStats(LambdaBroad,lmefit,options);

%%
figure
plot(model.y,STAT_FIT.TABLE.Estimate,'o')
grid
xlabel('y observed')
ylabel('y fitted')
title('Orthodont Data Fitted by HPMIXED')

%% Fit the linear mixed model by FITLME
tic;
%formula  = 'distance ~ age + Sex + age:Sex + (age | Subject)';
formula  = 'distance ~ age*Sex + (1 | Subject) + (age-1 | Subject)';
lme = fitlme(OrthodontData,formula,'DummyVarCoding','effects','FitMethod','REML');
%lme = fitlme(OrthodontData,formula,'DummyVarCoding','reference','FitMethod','REML');
toc

disp(lme)
anova(lme,'DFMethod','Satterthwaite')
%%

