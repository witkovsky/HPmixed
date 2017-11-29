%% EAXAMPLE: InstEval Data
%  Data analysis by Linear Mixed Model with MATLAB function HPMIXED

% Author: Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 06-Jan-2014 13:52:11

%% Data Description
% University lecture evaluations by students at ETH Zurich, anonymized for
% privacy protection. This is an interesting “medium” sized example of a
% partially nested mixed effect model.
% 
% A data frame with 73421 observations on the following 7 variables.
% * s - a factor with levels 1:2972 denoting individual students.
% * d - a factor with 1128 levels from 1:2160, denoting individual
%   professors or lecturers.
% * studage - an ordered factor with levels 2 < 4 < 6 < 8, denoting
%   student's “age” measured in the semester number the student has been
%   enrolled.
% * lectage - an ordered factor with 6 levels, 1 < 2 < ... < 6, measuring
%   how many semesters back the lecture rated had taken place.
% * service - a binary factor with levels 0 and 1; a lecture is a “service”,
%   if held for a different department than the lecturer's main one.
% * dept - a factor with 14 levels from 1:15, using a random code for the
%   department of the lecture.
% * y - a numeric vector of ratings of lectures by the students, using the
%   discrete scale 1:5, with meanings of ‘poor’ to ‘very good’.
% 
% Each observation is one student's rating for a specific lecture (of one
% lecturer, during one semester in the past).
% 
% The main goal of the survey is to find “the best liked prof”, according
% to the lectures given. Statistical analysis of such data has been the
% basis for a (student) jury selecting the final winners.
% 
% The present data set has been anonymized and slightly simplified on
% purpose.
% 
% R Example:
% system.time(fm <- lmer(y ~ dept*service + (1|s) + (1|d), InstEval))
%   
% Documentation reproduced from R package lme4,
% Version: 1.0-5, Date: 2013-10-25.
% Authors: Douglas Bates, Martin Maechler, Ben Bolker, Steven Walker

%% Load dataset InstEvalData
clear
load dsInstEvalData

%% Create the model design matrices 
formula  = 'y ~ dept*service + (1|s) + (1|d)';

%opts.dummyVarCode = 'reference';
opts.dummyVarCode = 'effects';
model = hpmixedmodel(InstEvalData,formula,opts);
model.Description = 'InstEvalData: InstEval R Example';

%% Fit the linear mixed model by HPMIXED with limitted  output
%opts.verbose = false;
opts.ddfMethod = 'residual';
% opts.ddfMethod = 'Satterthwaite';
opts.isExactGrad = 0;
opts.isExactFIreml = 0;

lmefit = hpmixed(model,opts);

disp(lmefit)
disp(lmefit.ModelInfo)
disp(lmefit.fixedEffects.Statistics.TABLE(1:20,:))
disp(lmefit.randomEffects.Statistics.TABLE(1:20,:))
disp(lmefit.varianceComponents.Statistics.TABLE)

%% ANOVA
STAT_ANOVA = getAnova(lmefit);
disp(STAT_ANOVA)

%% Fit the linear mixed model by FITLME
tic;
lme = fitlme(InstEvalData,formula,'FitMethod','REML');
toc

disp(lme)