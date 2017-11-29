%% EAXAMPLE: Rat Data
%  Data analysis by Linear Mixed Model with MATLAB function HPMIXED

% Author: Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 06-Jan-2014 13:52:11

%% Data Description
%    Rat data set used in:
%    Verbeke and Lesaffre (1999), Applied Statistics, 48, 363-375
%    Verbeke and Molenberghs (2000), New-York: Springer-Verlag
%    Verbeke and Molenberghs (2003), Biometrics, 59, 254-262
%    Gelman et al (2005), Biometrics, 61, 74-8
%
%    Variables:
%       (1) obs: observation number
%       (2) treat: treament group ('con': control; 'hig': high dose; 'low': low dose)
%       (3) rat: rat identification number
%       (4) age: age of the rat at the moment the observation is made
%       (5) respons: the response measured

%% Load dataset RatData 
clear
load dsRatData

%% Create the model design matrices 
% t = log(1+(age-45)/10);
%formula  = 'respons ~  t + t:treat +  (1 | rat) + (-1 + t | rat)';
formula  = 'respons ~ t +  t:treat +  (1 | rat)';
opts.dummyVarCode = 'reference';
%opts.dummyVarCode = 'effects';
%opts.dummyVarCode = 'full';
model = hpmixedmodel(RatData,formula,opts);
model.Description = 'RatData: Data set used in Verbeke and Lesaffre (1999)';

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

figure
plot(model.y,STAT_FIT.TABLE.Estimate,'o')
grid
xlabel('y observed')
ylabel('y fitted')
title('Rat Data Fitted by HPMIXED')