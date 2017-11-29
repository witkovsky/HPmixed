%% EAXAMPLE: Toenail Data
%  Data analysis by Linear Mixed Model with MATLAB function HPMIXED

% Author: Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 06-Jan-2014 13:52:11

%% Data Description
% Toenail data set used in:
%    Verbeke and Molenberghs (2000), New-York: Springer-Verlag
%    Verbeke, Lesaffre, and Spiessens (2001), Drug Information Journal, 35, 419-434.
% 
%    Variables:
%       (1) obs: observation number
%       (2) treat: treament group (0: Itraconazol (group B); 1: Lamisil (group A))
%       (3) id: subject identification number
%       (4) time: time at which the observation is taken (months)
%       (5) respons: the response measured (unaffected naillength, mm)
      
%% Load dataset RatData 
clear
load dsToenailData

%% Create the model design matrices 
% t = log(1+(age-45)/10);
%formula  = 'response ~ treat +  time + time:treat +  (1 | id)';
formula  = 'response ~ -1 + treat +  time:treat +  (1 | id)';
%opts.dummyVarCode = 'reference';
%opts.dummyVarCode = 'effects';
opts.dummyVarCode = 'full';
model = hpmixedmodel(ToenailData,formula,opts);
model.Description = 'Toenail Data: Data set used in Verbeke and Molenberghs (2000)';

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
title('Toenail Data Fitted by HPMIXED')