%% EAXAMPLE: Students Data
%  Data analysis by Linear Mixed Model with MATLAB function HPMIXED

% Author: Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 06-Jan-2014 13:52:11

%% Data Description
% vam_data from the R package GPvam.
%
% A simulated data set used to illustrate the functionality of the package.
% The data are simulated according to the VP model, and demonstrate the
% stability of the program in the presence of perfectly correlated future
% year effects.

%% Load dataset RatData 
clear
load dsStudentsData

%% Create the model design matrices 
%formula  = 'y ~ -1 + year +  (1 | teacher) + (1 | student)';
%formula  = 'y ~ -1 + year + contvar + (1 | teacher) + (1 | student)';
formula  = 'y ~ teacher + (1 | student)';
%formula  = 'y ~ year + contvar +  teacher + (1 | student)';
opts.dummyVarCode = 'reference';
%opts.dummyVarCode = 'effects';
%opts.dummyVarCode = 'full';
model = hpmixedmodel(StudentsData,formula,opts);
model.Description = 'StudentsData: vam_data from the R package GPvam)';

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
title('Students Data Fitted by HPMIXED')