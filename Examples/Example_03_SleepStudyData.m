%% EAXAMPLE: Sleep Study Data
%  Data analysis by Linear Mixed Model with MATLAB function HPMIXED

% Author: Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 06-Jan-2014 13:52:11

%% Data Description
% The average reaction time per day for subjects in a sleep deprivation
% study. On day 0 the subjects had their normal amount of sleep. Starting
% that night they were restricted to 3 hours of sleep per night. The
% observations represent the average reaction time on a series of tests
% given each day to each subject.
% 
% A data frame with 180 observations on the following 3 variables.
% * Reaction - Average reaction time (ms)
% * Days - Number of days of sleep deprivation
% * Subject - Subject number on which the observation was made.
% 
% These data are from the study described in Belenky et al. (2003), for the
% sleep-deprived group and for the first 10 days of the study, up to the
% recovery period.
% 
% References
%
% Gregory Belenky, Nancy J. Wesensten, David R. Thorne, Maria L. Thomas,
% Helen C. Sing, Daniel P. Redmond, Michael B. Russo and Thomas J. Balkin
% (2003) Patterns of performance degradation and restoration during sleep
% restriction and subsequent recovery: a sleep dose-response study. Journal
% of Sleep Research 12, 1–12.
%
% R Example:
% fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)
% 
% Documentation reproduced from R package lme4,
% Version: 1.0-5, Date: 2013-10-25.
% Authors: Douglas Bates, Martin Maechler, Ben Bolker, Steven Walker

%% Load dataset SleepStudyData 
clear
load dsSleepStudyData
%% Load dataset SpeepStudyData and create the model design matrices 
formula  = 'Reaction ~ Days + (1 | Subject) + (-1+Days | Subject)';
%opts.dummyVarCode = 'reference';
opts.dummyVarCode = 'effects';
model = hpmixedmodel(SleepStudyData,formula,opts);
model.Description = 'SleepStudyData: R lme4 Example';

%% Fit the linear mixed model by HPMIXED with limitted  output
% Use hpmixedmodel to create model from formula and dataset (requires
% MATLAB 2013b + Statistics Toolbox).
% Otherwise use the provided model structure, or construct the model 
% structure manually
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
title('SleepStudy Data Fitted by HPMIXED')

%% Plot the fitted subjects
% yfit =STAT_FIT.TABLE.Estimate;
% day = SleepStudyData.Days;
% t = linspace(0,9)';
% FE = lmefit.fixedEffects.Estimates;
% RE = lmefit.randomEffects.Estimates;
% plot(day,model.y,'o')
% hold on
% for i = 1:18
% color = rand(3,1);
% plot(day(1:10),model.y((i-1)*10+(1:10)),'o','Color',color,'LineWidth',2)
% line(t,(FE(1)+RE(i))+(FE(2)+RE(18+i))*t,'Color',color,'LineWidth',2)
% end
% hold off

%% Fit the linear mixed model by MIXED
% dim = [18 18];
% s20 = [1 1 1];
% method = 2;
% tic;
% [s2,b,u,Is2,C] = mixed(model.y,model.X,model.Z,dim,s20,method);
% toc
% disp(s2)

%% Fit the linear mixed model by FITLME
% tic;
% lme = fitlme(SleepStudyData,formula,'FitMethod','REML');
% toc
% 
% disp(lme)
%% EXAMPLE: (BROAD INFERENCE SPACE for LSMEANS of the factor Subject)

options.STAT.inference = 'lsmeans';
options.STAT.inferenceSpace = 'broad';
%options.STAT.inferenceSpace = 'narrow';
[Lambda,options]  = getLambda({'Subject'},model,SleepStudyData,options);
options.STAT.alpha = 0.01;
STAT_Subject = getStats(Lambda,lmefit,options);
disp(STAT_Subject)
%% Fit the linear mixed model by FITLME
tic;
lme = fitlme(SleepStudyData,formula,'DummyVarCoding','effects','FitMethod','REML');
toc

disp(lme)

anova(lme,'DFMethod','Satterthwaite')