function model = hpmixedmodel_Complex(fmri,AMP,SLP,aSLP,TMP,options)
% hpmixedmodel_Complex Creates the complex (very specific hard-wired)
%  linear mixed model for given continuous variables fmri, AMP, SLP, aSLP,
%  TMP, and the nominal variables Subj and session, specified with the
%  model formula:
%   formula = 'fmri ~ 1 + AMP + SLP + aSLP + TMP + Subj + AMP:Subj + ...
%              SLP:Subj + aSLP:Subj + TMP:Subj + (1 | session) + ...
%              (AMP - 1 | session)  + (SLP - 1 | session) + ...
%              (aSLP - 1 | session) + (TMP - 1 | session)';
%
%  Here, the method for creating the dummy variables is based on the value
%  set by options.dummyVarCode (default is options.dummyVarCode = 'effects'
%  alternative choice is options.dummyVarCode = 'reference' or
%  options.dummyVarCode = 'full').    
%
% SYNTAX:
%  model = hpmixedmodel_Complex(fmri,AMP,SLP,aSLP,TMP)
%  or
%  model = hpmixedmodel_Complex(fmri,AMP,SLP,aSLP,TMP,options)
%
% INPUTS:
%  fmri - n-dimensional vector of the observed fmri signal
%  AMP  - n-dimensional vector of the observed AMP signal
%  SLP  - n-dimensional vector of the observed SLP signal
%  aSLP - n-dimensional vector of the observed aSLP signal
%  TMP  - n-dimensional vector of the observed TMP signal
% *options* is a structure with the follwing properties and default values:
%  options.dummyVarCode = 'effects' (also 'reference' or 'full').
% 
% EXAMPLE
% load PainData.mat
% fmri = PainData.fmri;
% AMP = PainData.AMP;
% SLP = PainData.SLP;
% aSLP = PainData.aSLP;
% TMP = PainData.TMP;
% model   = hpmixedmodel_Complex(fmri,AMP,SLP,aSLP,TMP);
% clear options
% options.ddfMethod = 'Satterthwaite';
% lmefit  = hpmixed(model,options);
% STAT_FE = getStats('fixed',lmefit,options);
% disp(STAT_FE.TABLE);
% options.STAT.ddfMethod = 'ChiSquare';
% anovaChi = getAnova(lmefit,options);
% disp(anovaChi.TABLE);
% options.STAT.ddfMethod = 'FaiCornelius';
% anovaF = getAnova(lmefit,options);
% disp(anovaF.TABLE);

% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: '18-Aug-2022 20:03:19

%% CHECK the inputs
narginchk(5,6);

if nargin < 6, options = []; end
if ~isfield(options, 'dummyVarCode'), options.dummyVarCode = 'effects'; end

% Set the fixed parameters specific for the considered model
n = 235600;
p = 190; 
q = 760; 
nSessions = 152; 
nSubjects = 38; 
idxSessions = repmat((1:nSessions),1,n/nSessions); 
auxid = repmat((1:nSubjects),4,1);
auxid = auxid(:)';
idxSubjects = repmat(auxid,1,n/nSessions);

model.Description = 'Complex model';
model.Formula.char = 'fmri ~ AMP + SLP + aSLP + TMP + Subj + AMP:Subj + SLP:Subj+ aSLP:Subj+ TMP:Subj + (1 | session) + (AMP - 1 | session)+ (SLP - 1 | session) + (aSLP - 1 | session)+ (TMP - 1 | session)';
model.Formula.FELinearFormula.TermNames = {'(Intercept)' 'AMP' ...
    'SLP' 'aSLP' 'TMP' 'Subj' 'AMP:Subj' 'SLP:Subj' 'aSLP:Subj' ...
    'TMP:Subj'}';
model.FixedInfo.XCols2Terms = [1 2 3 4 5 6*ones(1,nSubjects-1) ... 
    7*ones(1,nSubjects-1) 8*ones(1,nSubjects-1) 9*ones(1,nSubjects-1) ...
    10*ones(1,nSubjects-1)];
model.FixedInfo.XColNames =  {'(Intercept)'  'AMP'  'SLP'  'aSLP'  'TMP' ...
  'Subj_2' 'Subj_3' 'Subj_4' 'Subj_5' 'Subj_6' 'Subj_7' 'Subj_8' ...
  'Subj_9' 'Subj_10' 'Subj_11' 'Subj_12' 'Subj_13' 'Subj_14' 'Subj_15' ...
  'Subj_16' 'Subj_17' 'Subj_18' 'Subj_19' 'Subj_20' 'Subj_21' 'Subj_22' ...
  'Subj_23' 'Subj_24' 'Subj_25' 'Subj_26' 'Subj_27' 'Subj_28' 'Subj_29' ...
  'Subj_30' 'Subj_31' 'Subj_32' 'Subj_33' 'Subj_34' 'Subj_35' 'Subj_36' ...
  'Subj_37' 'Subj_38' ...
  'AMP:Subj_2' 'AMP:Subj_3' 'AMP:Subj_4' 'AMP:Subj_5' 'AMP:Subj_6' ...
  'AMP:Subj_7' 'AMP:Subj_8' 'AMP:Subj_9' 'AMP:Subj_10' 'AMP:Subj_11' ...
  'AMP:Subj_12' 'AMP:Subj_13' 'AMP:Subj_14' 'AMP:Subj_15' 'AMP:Subj_16' ...
  'AMP:Subj_17' 'AMP:Subj_18' 'AMP:Subj_19' 'AMP:Subj_20' 'AMP:Subj_21' ...
  'AMP:Subj_22' 'AMP:Subj_23' 'AMP:Subj_24' 'AMP:Subj_25' 'AMP:Subj_26' ...
  'AMP:Subj_27' 'AMP:Subj_28' 'AMP:Subj_29' 'AMP:Subj_30' 'AMP:Subj_31' ...
  'AMP:Subj_32' 'AMP:Subj_33' 'AMP:Subj_34' 'AMP:Subj_35' 'AMP:Subj_36' ...
  'AMP:Subj_37' 'AMP:Subj_38' ...
  'SLP:Subj_2' 'SLP:Subj_3' 'SLP:Subj_4' 'SLP:Subj_5' 'SLP:Subj_6' ...
  'SLP:Subj_7' 'SLP:Subj_8' 'SLP:Subj_9' 'SLP:Subj_10' 'SLP:Subj_11' ...
  'SLP:Subj_12' 'SLP:Subj_13' 'SLP:Subj_14' 'SLP:Subj_15' 'SLP:Subj_16' ...
  'SLP:Subj_17' 'SLP:Subj_18' 'SLP:Subj_19' 'SLP:Subj_20' 'SLP:Subj_21' ...
  'SLP:Subj_22' 'SLP:Subj_23' 'SLP:Subj_24' 'SLP:Subj_25' 'SLP:Subj_26' ...
  'SLP:Subj_27' 'SLP:Subj_28' 'SLP:Subj_29' 'SLP:Subj_30' 'SLP:Subj_31' ...
  'SLP:Subj_32' 'SLP:Subj_33' 'SLP:Subj_34' 'SLP:Subj_35' 'SLP:Subj_36' ...
  'SLP:Subj_37' 'SLP:Subj_38' ...
  'aSLP:Subj_2' 'aSLP:Subj_3' 'aSLP:Subj_4' 'aSLP:Subj_5' 'aSLP:Subj_6' ...
  'aSLP:Subj_7' 'aSLP:Subj_8' 'aSLP:Subj_9' 'aSLP:Subj_10' 'aSLP:Subj_11' ...
  'aSLP:Subj_12' 'aSLP:Subj_13' 'aSLP:Subj_14' 'aSLP:Subj_15' 'aSLP:Subj_16' ...
  'aSLP:Subj_17' 'aSLP:Subj_18' 'aSLP:Subj_19' 'aSLP:Subj_20' 'aSLP:Subj_21' ...
  'aSLP:Subj_22' 'aSLP:Subj_23' 'aSLP:Subj_24' 'aSLP:Subj_25' 'aSLP:Subj_26' ...
  'aSLP:Subj_27' 'aSLP:Subj_28' 'aSLP:Subj_29' 'aSLP:Subj_30' 'aSLP:Subj_31' ...
  'aSLP:Subj_32' 'aSLP:Subj_33' 'aSLP:Subj_34' 'aSLP:Subj_35' 'aSLP:Subj_36' ...
  'aSLP:Subj_37' 'aSLP:Subj_38' ...
  'TMP:Subj_2' 'TMP:Subj_3' 'TMP:Subj_4' 'TMP:Subj_5' 'TMP:Subj_6' ...
  'TMP:Subj_7' 'TMP:Subj_8' 'TMP:Subj_9' 'TMP:Subj_10' 'TMP:Subj_11' ...
  'TMP:Subj_12' 'TMP:Subj_13' 'TMP:Subj_14' 'TMP:Subj_15' 'TMP:Subj_16' ...
  'TMP:Subj_17' 'TMP:Subj_18' 'TMP:Subj_19' 'TMP:Subj_20' 'TMP:Subj_21' ...
  'TMP:Subj_22' 'TMP:Subj_23' 'TMP:Subj_24' 'TMP:Subj_25' 'TMP:Subj_26' ...
  'TMP:Subj_27' 'TMP:Subj_28' 'TMP:Subj_29' 'TMP:Subj_30' 'TMP:Subj_31' ...
  'TMP:Subj_32' 'TMP:Subj_33' 'TMP:Subj_34' 'TMP:Subj_35' 'TMP:Subj_36' ...
  'TMP:Subj_37' 'TMP:Subj_38'};
model.dimRE = [nSessions;nSessions;nSessions;nSessions;nSessions];
model.n = n;
model.p = p;
model.q = q;
model.H0 = [];
model.ZXy = [];
model.y2 = [];

%% model.y is optional model output / not necessary for further computation
y = fmri(:);

% optional output
% model.y = y; 

%% model.y2 is required model output / necessary for further computation
y2 = y'*y;

% required output
model.y2 = y2;

%% model.X is optional model output / not necessary for further computation

% HERE, the (nxp)-matrix X (the fixed effects design matrix) consists of
% p = 5 columns / one column for each fixed effect:
% 1. Intercept, 2. AMP, 3. SLP, 4. aSLP, 5. TMP

% X0 is (nx5)-matrix related to the fixed effects [1 AMP SLP aSLP TMP]
X0 = [ones(n,1) AMP(:) SLP(:) aSLP(:) TMP(:)];

% X1 is (nx37)-matrix related to the fixed effects of Subj
X1 = sparse(1:n,idxSubjects,1,n,nSubjects,n); 

% X2 is (nx37)-matrix related to the fixed effects of AMP:Subj
X2 = sparse(1:n,idxSubjects,AMP,n,nSubjects,n); 

% X3 is (nx37)-matrix related to the fixed effects of SLP:Subj
X3 = sparse(1:n,idxSubjects,SLP,n,nSubjects,n); 

% X4 is (nx37)-matrix related to the fixed effects of aSLP:Subj
X4 = sparse(1:n,idxSubjects,aSLP,n,nSubjects,n); 

% X5 is (nx37)-matrix related to the fixed effects of TMP:Subj
X5 = sparse(1:n,idxSubjects,TMP,n,nSubjects,n); 

if strcmp(options.dummyVarCode,'effects')
    X1 = X1 - X1(:,nSubjects); X1(:,end) = [];
    X2 = X2 - X2(:,nSubjects); X2(:,end) = [];
    X3 = X3 - X3(:,nSubjects); X3(:,end) = [];
    X4 = X4 - X4(:,nSubjects); X4(:,end) = [];
    X5 = X5 - X5(:,nSubjects); X5(:,end) = [];
elseif strcmp(options.dummyVarCode,'reference')
    X1 = X1(:,2:end);
    X2 = X2(:,2:end);
    X3 = X3(:,2:end);
    X4 = X4(:,2:end);
    X5 = X5(:,2:end);
% else options.dummyVarCode = 'full'
end

% This is concatenaded (n*190)-matrix X (the design matrix for all fixed
% effects) 
X = [X0 X1 X2 X3 X4 X5];

% optional output
% model.X = X;

%% model.Z is optional model output / not necessary for further computation
% Create the Z matrix

% HERE, the (nxq)-matrix Z (the random effects design matrix) consists of
% q = 5*152 columns / 152 columns for each random effect:
% 1. (1|session), 2. (AMP|session), 3. (SLP|session), 4. (aSLP|session), 5.
% (TMP|session)

% Z1 is (nx152)-matrix related to the 1st random effect (1|session)
Z1 = sparse(1:n,idxSessions,1,n,nSessions,n);

% Z2 is (nx152)-matrix related to the 2nd random effect (AMP|session)
Z2 = sparse(1:n,idxSessions,AMP,n,nSessions,n);

% Z3 is (nx152)-matrix related to the 3rd random effect (SLP|session)
Z3 = sparse(1:n,idxSessions,SLP,n,nSessions,n);

% Z4 is (nx152)-matrix related to the 4th random effect (aSLP|session)
Z4 = sparse(1:n,idxSessions,aSLP,n,nSessions,n);

% Z5 is (nx152)-matrix related to the 5th random effect (TMP|session)
Z5 = sparse(1:n,idxSessions,TMP,n,nSessions,n);

% This is concatenaded (n*760)-matrix Z (the design matrix for all random
% effects)  
Z = [Z1 Z2 Z3 Z4 Z5];

% optional output
% model.Z = Z;

%% model.ZXy is required model output / necessary for further computation
% ZXy is m-dimensional vector of cross-products

ZXy = [Z' * y ; X' * y];

% required output
model.ZXy = ZXy;

%% model.H0 is required model output / necessary for further computation
% H0 is the (mxm)-dimensional matrix of cross-products / matrix of the
% Mixed Model Eqiuations (MME)

H0 = [ triu(Z' * Z) , Z' * X ; sparse(p,q) , triu(X' * X) ];

% required output
model.H0 = H0;