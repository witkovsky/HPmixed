function model = hpmixedmodel_Simple(fmri,AMP,SLP,aSLP,TMP,n,p,q,nSessions,idxSessions)
% hpmixedmodel_Simple Creates a simple (very specific hard-wired)
%  linear mixed model for given continuous variables fmri, AMP, SLP, aSLP,
%  TMP, and the nominal variable session, specified with the model formula:
%   formula = 'fmri ~ 1 + AMP + SLP + aSLP + TMP + (1 | session) + ... 
%              (AMP - 1 | session)  + (SLP - 1 | session) + ... 
%              (aSLP - 1 | session) + (TMP - 1 | session)'.
%
%  Here, the method for creating the dummy variables is based on
%  options_simple.dummyVarCode = 'reference' (if used in the function
%  hpmixedmodel). 
%
% SYNTAX:
%  model = hpmixedmodel_Simple(fmri,AMP,SLP,aSLP,TMP)
%  or
%  model = hpmixedmodel_Simple(fmri,AMP,SLP,aSLP,TMP,n,p,q,nSessions)
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
%         matrix). If empty, the default value is p = 5;  
%  q    - number of columns of the (nxq)-matrix Z (the random effects
%         design matrix). If empty, the default value is q = 5* 152 = 760
%         190; 
%  nSessions - number of sessions. If empty, the default value is 
%         nSessions = 152.
%  idxSessions - the index vector used to generate sparse Z matrices. If
%         empty the index vector is calculate by the algorithm. 
% 
% EXAMPLE
% load PainData.mat
% fmri = PainData.fmri;
% AMP = PainData.AMP;
% SLP = PainData.SLP;
% aSLP = PainData.aSLP;
% TMP = PainData.TMP;
% n = 235600;
% p = 5;
% q = 760;
% nSessions = 152;
% model   = hpmixedmodel_Simple(fmri,AMP,SLP,aSLP,TMP,n,p,q,nSessions);
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
% Ver.: 16-Aug-2022 09:29:54

%% CHECK the inputs
narginchk(5,10);
if nargin <  6, n = 235600; end
if nargin <  7, p = 5; end
if nargin <  8, q = 760; end
if nargin <  9, nSessions = 152; end
if nargin <  10, idxSessions = repmat((1:nSessions),1,n/nSessions); end

model.Description = 'Simple model';
model.Formula.char = 'fmri ~ AMP + SLP + aSLP + TMP + (1 | session) + (AMP - 1 | session) + (SLP - 1 | session) + (aSLP - 1 | session) + (TMP - 1 | session)';
model.Formula.FELinearFormula.TermNames = {'(Intercept)' 'AMP' ...
    'SLP' 'aSLP' 'TMP'}';
model.FixedInfo.XCols2Terms = [1 2 3 4 5];
model.FixedInfo.XColNames =  {'(Intercept)'  'AMP'  'SLP'  'aSLP'  'TMP'};
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

X = [ones(n,1) AMP(:) SLP(:) aSLP(:) TMP(:)];

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

H0 = [triu(Z' * Z) , Z' * X ; sparse(p,q) , triu(X' * X)];

% required output
model.H0 = H0;