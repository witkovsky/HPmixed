function model = hpmixedmodel_Complex1(fmri,AMP,SLP,aSLP,TMP,n,p,q,nSessions,nSubjects)
% hpmixedmodel_Complex1 Creates the complex (very specific hard-wired)
%  linear mixed model for given continuous variables fmri, AMP, SLP, aSLP,
%  TMP, and the nominal variables Subj and session, specified with the
%  model formula:
%   formula = 'fmri ~ 1 + AMP + SLP + aSLP + TMP + Subj + AMP:Subj + ...
%              SLP:Subj + aSLP:Subj + TMP:Subj + (1 | session) + ...
%              (AMP - 1 | session)  + (SLP - 1 | session) + ...
%              (aSLP - 1 | session) + (TMP - 1 | session)';
%
% Here, the method for creating the dummy variables is based on
% options_simple.dummyVarCode = 'reference' (if used in the function
% hpmixedmodel). 
%
% SYNTAX:
% model = hpmixedmodel_Complex1(fmri,AMP,SLP,aSLP,TMP,n,p,nSessions)
%
% INPUTS:
%  fmri - n-dimensional vector of the observed fmri signal
%  AMP  - n-dimensional vector of the observed AMP signal
%  SLP  - n-dimensional vector of the observed SLP signal
%  aSLP - n-dimensional vector of the observed aSLP signal
%  TMP  - n-dimensional vector of the observed TMP signal
%  n    - length of the vector variables. If empty, the default value is
%         n = 235600; 
%  p    - number of columns of the (nxp)-matrix X (the fixed effects design
%         matrix). If empty, the default value is p = 5 + 5*(37-1) = 190;  
%  q    - number of columns of the (nxq)-matrix Z (the random effects
%         design matrix). If empty, the default value is q = 5* 152 = 760
%         190; 
%  nSessions - number of sessions. If empty, the default value is 
%         nSessions = 152.  
%  nSubjects - number of subjects. If empty, the default value is 
%         nSubjects = 38. 
% 
% EXAMPLE
% load PainData.mat
% fmri = PainData.fmri;
% AMP = PainData.AMP;
% SLP = PainData.SLP;
% aSLP = PainData.aSLP;
% TMP = PainData.TMP;
% n = 235600;
% p = 190;
% q = 760;
% nSessions = 152;
% nSubjects = 38;
% model   = hpmixedmodel_Complex1(fmri,AMP,SLP,aSLP,TMP,n,p,q,nSessions,nSubjects);
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
% Ver.: 14-Aug-2022 10:41:01

%% CHECK the inputs
narginchk(5,10);
if nargin <  10, nSubjects = 38; end
if nargin <  9, nSessions = 152; end
if nargin <  8, q = 760; end
if nargin <  7, p = 190; end
if nargin <  6, n = 235600; end

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
%%  Auxiliary index vectors

idxSubj1  = [1:nSessions:n; 2:nSessions:n; 3:nSessions:n; 4:nSessions:n];
idxSubj1  = idxSubj1(:); 
idxSes1   = 1:nSessions:n;
idxSes1   = idxSes1(:);

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

nXcol = nSubjects-1;
% X0 is (nx5)-matrix related to the fixed effects [1 AMP SLP aSLP TMP]
X0 = [ones(n,1) AMP(:) SLP(:) aSLP(:) TMP(:)];

% X1 is (nx37)-matrix related to the fixed effects of Subj
X1 = sparse(n,nXcol,n);
for i = 1:nXcol
    X1(idxSubj1+i,i) = 1;
end

% X2 is (nx37)-matrix related to the fixed effects of AMP:Subj
X2 = sparse(n,nXcol,n);
for i = 1:nXcol
    X2(idxSubj1+i,i) = AMP(idxSubj1+i);
end

% X3 is (nx37)-matrix related to the fixed effects of SLP:Subj
X3 = sparse(n,nXcol,n);
for i = 1:nXcol
    X3(idxSubj1+i,i) = SLP(idxSubj1+i);
end

% X4 is (nx37)-matrix related to the fixed effects of aSLP:Subj
X4 = sparse(n,nXcol,n);
for i = 1:nXcol
    X4(idxSubj1+i,i) = aSLP(idxSubj1+i);
end

% X5 is (nx37)-matrix related to the fixed effects of TMP:Subj
X5 = sparse(n,nXcol,n);
for i = 1:nXcol
    X5(idxSubj1+i,i) = TMP(idxSubj1+i);
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
Z1 = sparse(idxSes1,1,1,n,nSessions,n);
for i = 2:nSessions
    Z1(idxSes1+i-1,i) = 1;
end

% Z2 is (nx152)-matrix related to the 2nd random effect (AMP|session)
Z2 = sparse(idxSes1,1,AMP(idxSes1),n,nSessions,n);
for i = 2:nSessions
    Z2(idxSes1+i-1,i) = AMP(idxSes1+i-1);
end

% Z3 is (nx152)-matrix related to the 3rd random effect (SLP|session)
Z3 = sparse(idxSes1,1,SLP(idxSes1),n,nSessions,n);
for i = 2:nSessions
    Z3(idxSes1+i-1,i) = SLP(idxSes1+i-1);
end

% Z4 is (nx152)-matrix related to the 4th random effect (aSLP|session)
Z4 = sparse(idxSes1,1,aSLP(idxSes1),n,nSessions,n);
for i = 2:nSessions
    Z4(idxSes1+i-1,i) = aSLP(idxSes1+i-1);
end

% Z5 is (nx152)-matrix related to the 5th random effect (TMP|session)
Z5 = sparse(idxSes1,1,TMP(idxSes1),n,nSessions,n);
for i = 2:nSessions
    Z5(idxSes1+i-1,i) = TMP(idxSes1+i-1);
end

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

H0 = triu(sparse([Z' * Z , Z' * X ; X' * Z , X' * X]));

% required output
model.H0 = H0;
