function lmefit = hpmixed(model,options)
%HPMIXED  REML based fitting of the linear mixed effects models (LME models)
%         of the form: y = X*b + Z*u + e, with a simple variance covariance
%         structure (VC: variance components), by solving the Henderson's
%         mixed model equations.
%
%         Here we assume that the design matrix X for the fixed effects is
%         a full-ranked (n x p)-matrix, such that X = [X_0 X_1 ... X_nFE],
%         where X_0 represents the intercept (if included) and X_i, for i =
%         1,...,nFE, are the (n x p_i) design matrices for the fixed
%         effects.
%
%         The assumed variance-covariance matrix Var(y) = V = Z*G*Z' + R,
%         is V = Z_1 * G_1 * Z_1' + ... + Z_nRE * G_nRE * Z_nRE' + R, where
%         G_i = sigma^2_i * eye(q_i), sigma^2_i >= 0, i = 1,...,nRE, and
%         R = sigma^2_{nRE+1} * eye(n), sigma^2_{nRE+1} > 0.
%
% SYNTAX:
% lmefit = hpmixed(model,options)
%
% INPUTS:
% *model* is a structure with the following minimal properties:
% model.y           - n-vector of observations (dependent variable),
% model.X           - (n x p) full-ranked (possibly sparse) design matrix
%                     for the fixed effects,
% model.Z           - (n x q) full-ranked (possibly sparse) design matrix
%                     for random effects, typically Z = [Z_1 Z_2 ... Z_nRE],
% model.dimRE       - nRE-dimensional vector of random effects group sizes.
%
% Alternatively, for large models it can be useful to specify the model by
% the following properties
% model.H0          - H0 is a sparse upper triangular of the (m x m) matrix
%                     [Z X]'*[Z X], which is the constructor for the MMEs,
%                     where m = q + p,
% model.ZXy         - m-vector [Z'*y; X'*y],
% model.y2          - scalar value y2 = y'*y,
% model.dimRE       - nRE-dimensional vector of random effects group sizes,
% model.dimFE       - nFE-dimensional vector of fixed effects group sizes,
% model.n           - n is number of observations, dimension of y,
% model.p           - p is number of fixed effects, columns of X.
%
% Further optional properties of the model structure are:
% model.Description - character string
% model.FixedInfo   - structure with further details on fixed effects
% model.Formula     - structure by classreg.regr.LinearMixedFormula
% model.GroupingInfo- structure with further details on grouping variables
% model.RandomInfo  - structure with further details on random effects
%
% The MODEL structure can be generated from given DATASET and FORMULA by
% using the developed function hpmixedmodel, which is based on
% functionality of the LinearMixedModels class (Statistics Toolbox, MATLAB
% 2013b).
%
% *options* is a structure with the follwing minimal properties:
% options.FitMethod - default (for now the only possible) method is 'REML',
% options.tolerance - default value is sqrt(100*eps),
% options.tolZero   - default value is 1e-12,
% options.verbose   - default value is: true.
%
% EXAMPLE:
% load dsSplitPlotData
% formula  = 'y ~ A + B + A:B + (1 | Block) + (1 | Block:A)';
% model = hpmixedmodel(SplitPlotData,formula);
% lmefit = hpmixed(model)
%
% REFERENCES:
% Witkovský, V. (2012). Estimation, testing, and prediction regions of the
% fixed and random effects by solving the Henderson’s mixed  model
% equations. Measurement Science Review 12(6): 234-248. 
% DOI: 10.2478/v10048-012-0033-6
%
% REMARKS:
%         The Matlab algorithm HPMIXED is an update of the already obsolete
%         version of the Matlab algorithm MIXED, still available at
%         http://www.mathworks.com/matlabcentral/fileexchange/200-mixed
%
%         HPMIXED is an alternative to the standard Matlab algorithms
%         FITLME (class MixedModels), based on different strategy.
% 
%         This is an experimantal (working) version, supposed to be
%         improved and generalized in future.
%
%         HPMIXED can be generalized to the LME models with more general
%         structure of variance-covariance component matrices G_i =
%         G(theta_i), i = 1,...,nRE, with linear parametrization: G_i =
%         {theta_i}_1 * {G_i}_1 + ... + {theta_i}_nREi * {G_i}_nREi, where
%         theta_i are nVC_i-dimensional vector parameters of the
%         variance-covariance components.
%
%         Implementation of the Kenward-Roger approach for DDF (denominator
%         degrees of freedom), as suggested in [1], is another open issue,
%         see also
%         http://www.mathworks.com/matlabcentral/fileexchange/40875-ddfmixed

% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 16-Feb-2014 14:39:28

%% TICTOC start the system-time clock
tic;

%% CHECK the inputs
narginchk(1,2);
if nargin < 2, options = []; end

%% SET the default options
if ~isfield(options, 'FitMethod'), options.FitMethod = 'REML'; end
if ~isfield(options, 'ddfMethod'), options.ddfMethod = 'Satterthwaite'; end
if ~isfield(options, 'alpha'), options.alpha = 0.05; end
if ~isfield(options, 'tolerance'), options.tolerance = 1e4*eps; end
if ~isfield(options, 'tolZero'), options.tolZero = 1e-12; end
if ~isfield(options, 'isExactGrad'), options.isExactGrad = true; end
if ~isfield(options, 'isExactFIreml'), options.isExactFIreml = true; end
if ~isfield(options, 'verbose'), options.verbose = true; end
if ~isfield(options, 'ExactMethodLimitNNZ')
    options.ExactMethodLimitNNZ = 40000;
end

%% CHECK the model structure
if  isfield(model,'H0') && isfield(model,'ZXy') && isfield(model,'y2') ...
        && isfield(model,'dimRE') && isfield(model,'n') ...
        && isfield(model,'p')
    H0 = model.H0;
    ZXy = model.ZXy;
    y2 = model.y2;
    dimRE = model.dimRE;
    n = model.n;
    p = model.p;
    if ~issparse(H0)
        H0 = triu(sparse(H0));
    end
elseif isfield(model,'X') && isfield(model,'Z') && ...
        isfield(model,'y') && isfield(model,'dimRE')
    H0 = triu(sparse([model.Z' * model.Z , model.Z' * model.X ; ...
        model.X' * model.Z , model.X' * model.X]));
    ZXy = [model.Z' * model.y ; model.X' * model.y];
    y2 = model.y' * model.y;
    dimRE = model.dimRE;
    [n,p] = size(model.X);
    model.H0 = H0;
    model.ZXy = ZXy;
    model.y2 = y2;
    model.n = n;
    model.p = p;
else
    error('VW:hpmixed', 'Wrong input arguments')
end

if ~isfield(model, 'Description')
    model.Description = 'User defined model';
end

if ~isfield(model, 'Formula')
    model.Formula.char = 'Not available';
end

%% SET the internal variables
nRE = length(dimRE);
q = sum(dimRE);
nVC = nRE + 1;
m = p + q;

[invperm,H0,ZXy] = findperm(H0,ZXy,m,q);
perm(invperm) = 1:m;
H0DiagInd = 1+(0:m-1)*(m+1);
H0REDiagInd = H0DiagInd(invperm(1:q));
H0REDiag = H0(H0REDiagInd);
nonZerosH0 = nnz(H0);

if options.isExactGrad && options.isExactFIreml ...
        && (nonZerosH0 < options.ExactMethodLimitNNZ)
    tol = options.tolerance;
else
    tol = max(options.tolerance,eps^(1/2));
end

ind = cell(nRE,1);
Gd = cell(nRE,1);
startidx = 0;
for i = 1:(nRE)
    ind{i} = startidx + (1:dimRE(i));
    Gd{i} = sparse(invperm(ind{i}),invperm(ind{i}),ones(1,dimRE(i)),m,m);
    startidx = startidx + dimRE(i);
end

%% OUTPUT structure
lmefit.Title = 'Linear mixed-effects model fit by REML';
lmefit.DataInfo = model.Description;
lmefit.ModelInfo.Formula = model.Formula.char;
lmefit.ModelInfo.FitMethod = options.FitMethod;
lmefit.ModelInfo.I_REML_used = [];
lmefit.ModelInfo.ddfMethod = options.ddfMethod;
lmefit.ModelInfo.observationsNum = n;
lmefit.ModelInfo.fixedEffectsLength = p;
lmefit.ModelInfo.randomEffectsLength = q;
lmefit.ModelInfo.randomEffectsNum = nRE;
lmefit.ModelInfo.randomEffectsDim = dimRE';
lmefit.ModelInfo.varianceComponentsNum = nVC;
lmefit.ModelInfo.MMEmatRows = m;
lmefit.ModelInfo.MMEmatNonzeros = nonZerosH0;
lmefit.ModelInfo.MMEmatSparsity = nonZerosH0/(m*m);
lmefit.ModelInfo.loglik_REML = [];
lmefit.ModelInfo.neg2loglik_REML = [];
lmefit.ModelInfo.AIC = [];
lmefit.ModelInfo.BIC = [];
lmefit.ModelInfo.covergenceCrit = [];
lmefit.ModelInfo.covergenceTol = [];
lmefit.fixedEffects.Estimates = [];
lmefit.fixedEffects.StandardErors = [];
lmefit.fixedEffects.Statistics = [];
lmefit.fixedEffects.Details = [];
lmefit.randomEffects.Estimates = [];
lmefit.randomEffects.StandardErors = [];
lmefit.randomEffects.Statistics = [];
lmefit.randomEffects.Details = [];
lmefit.varianceComponents.Estimates = [];
lmefit.varianceComponents.StandardErors = [];
lmefit.varianceComponents.Statistics = [];
lmefit.varianceComponents.Details = [];
lmefit.Details.Matrices.ubeta = [];
lmefit.Details.Matrices.H_factor = [];
lmefit.Details.Matrices.I_REML = [];
lmefit.Details.Matrices.Perm = perm;
lmefit.Details.Matrices.InvPerm = invperm;
lmefit.Details.Matrices.H0 = H0;
lmefit.Details.Matrices.H0REDiagInd = H0REDiagInd;
lmefit.Details.Matrices.H0REDiag = H0REDiag;
lmefit.Details.Matrices.Gd = Gd;
lmefit.Details.Matrices.ind = ind;
lmefit.Details.Options = options;
if options.verbose
    lmefit.Details.Model = model;
end

%% ALGORITHM / Linear mixed-effects model fit by REML
if ~isfield(model, 'sig20')
    model.sig20 = ones(nRE+1,1);
end
sig2 = model.sig20;
sig20 = zeros(nRE+1,1);

loops = 0;
tol = sqrt(nonZerosH0) * tol;
covergenceCrit = norm((sig2-sig20)/sig2(end));
while  covergenceCrit > tol
    loops = loops + 1;
    sig20 = sig2;
    sig2 = sig20 / sig20(nRE+1);
    
    % fixed and random effects
    H0 = updateH0(sig20,H0,H0REDiagInd,H0REDiag,nRE,ind);
    Hfac = chol(H0);
    ubeta = Hfac \ (Hfac' \ ZXy);
    
    % variance components
    if options.isExactGrad && ...
            (nonZerosH0 < options.ExactMethodLimitNNZ)
        nu = zeros(nRE+1,1);
        nu(nRE+1) = n-p;
        for i = 1:nRE
            trii = sum(nonzeros(Hfac' \ Gd{i}).^2) / sig2(i);
            nu(i) = dimRE(i) - trii;
        end
    else
        nu = sig20.* gradient(@(sig2) -2 * maxloglikFun(sig2,H0,H0REDiag,...
            H0REDiagInd,ind,nRE,n,dimRE',p,sig20),sig20,[],1:nRE);
        nu(nRE+1) = n-p;
    end
    for i = 1:nRE
        ui = ubeta(invperm(ind{i}));
        sig2(i) = (ui' * ui) / nu(i);
    end
    sig2(nRE+1) = (y2 - ubeta' * ZXy) / nu(nRE+1);
    covergenceCrit = norm((sig2-sig20)/sig2(end));
end
Hfac = (round(Hfac/options.tolZero)*options.tolZero)/ sqrt(sig2(nRE+1));

%% SET the first part of the results
lmefit.Details.Matrices.H_factor = Hfac;
lmefit.Details.Matrices.ubeta = ubeta;
lmefit.fixedEffects.Estimates = ubeta(invperm(q+(1:p)));
lmefit.randomEffects.Estimates = ubeta(invperm(1:q));
lmefit.varianceComponents.Estimates = sig2;
lmefit.varianceComponents.Details.nu = nu;
lmefit.ModelInfo.covergenceCrit = covergenceCrit;
lmefit.ModelInfo.covergenceTol = tol;

%% Get maximum value of the REML log-likehood and the Fisher Information
[loglik_REML,AIC,BIC]  = loglikREML(sig2,Hfac,n,p,nRE,dimRE);

if options.isExactFIreml && (nonZerosH0 < options.ExactMethodLimitNNZ)
    I_REML = getFisherInfExact(lmefit);
    I_REML_used = 'Fisher Information Matrix / Exact';
elseif ~options.isExactFIreml
    I_REML = getFisherInfApprox(lmefit);
    I_REML_used = 'Fisher Information Matrix / Approx';
else
    I_REML = getFisherInfApprox(lmefit);
    I_REML_used = 'Fisher Information Matrix / Approx';
end

stdVC =  sqrt(diag(I_REML\eye(size(I_REML))));

lmefit.ModelInfo.I_REML_used = I_REML_used;
lmefit.ModelInfo.loglik_REML = loglik_REML;
lmefit.ModelInfo.neg2loglik_REML = -2*loglik_REML;
lmefit.ModelInfo.AIC = AIC;
lmefit.ModelInfo.BIC = BIC;
lmefit.Details.Matrices.I_REML = I_REML;
lmefit.varianceComponents.StandardErors = stdVC;

%% Basic statistics for fixed and/or random effects
if isfield(model,'FixedInfo')
    lmefit.fixedEffects.Names = model.FixedInfo.XColNames(:);
else
    lmefit.fixedEffects.Names = 'FEcoef_';
end

options.STAT.verbose = false;
options.STAT.ddfMethod = options.ddfMethod;

FE_STAT = getStats('fixed',lmefit,options);
RE_STAT = getStats('random',lmefit,options);
VC_STAT = getStats('variancecomponents',lmefit,options);

lmefit.fixedEffects.Statistics = FE_STAT;
lmefit.randomEffects.Statistics = RE_STAT;
lmefit.fixedEffects.StandardErors = FE_STAT.TABLE.SE;
lmefit.randomEffects.StandardErors = RE_STAT.TABLE.SE;
lmefit.varianceComponents.Statistics = VC_STAT;
lmefit.loops = loops;
lmefit.tictoc = toc;

%% Display the Statistics
if options.verbose
    disp('   MODEL INFORMATION  ')
    disp(lmefit.ModelInfo)
    disp('   FIXED EFFECTS  ')
    disp(lmefit.fixedEffects.Statistics.TABLE)
    disp('   RANDOM EFFECTS ')
    disp(lmefit.randomEffects.Statistics.TABLE)
    disp('   VARIANCE COMPONENTS ')
    disp(lmefit.varianceComponents.Statistics.TABLE)
    disp('   FURTHER Details  ')
    disp(lmefit.Details)
end

end     % END OF THE MAIN FUNCTION HPMIXED


%% FUNCTION findperm
function [perm,H0,ZXy] = findperm(H0,ZXy,m,q)
%function [perm,H0,ZXy] = findperm(H0,ZXy,m,q)
%FINDPERM Find optimum permutation for Cholesky factorization and rearrange
%         (permute) columns and rows of the matrix H0 and the vector ZXy.

% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 06-Oct-2013 18:45:48

ind = 1+(0:q-1)*(m+1);
H0(ind) = H0(ind) + 1;
[~,~,P] = chol(H0);
H0(ind) = H0(ind) - 1;
H0 = triu(P'*(H0 + triu(H0,1)')*P);
ZXy = P'*ZXy;
perm = (1:m)*P';

end
%% Function updateH0
function H0 = updateH0(sig2,H0,H0REDiagInd,H0REDiag,nRE,ind)
%updateH0 - Update profiled MME equations
%
% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 09-Jan-2014 21:55:07

sig2 = sig2 / sig2(nRE+1);
for i = 1:nRE
    H0(H0REDiagInd(ind{i})) = H0REDiag(ind{i}) + 1/sig2(i);
end
end     % END of updateH0

%% FUNCTION loglik_REML
function [loglik_REML,AIC,BIC] = loglikREML(sig2,Hfac,n,p,nRE,dimRE)
% loglikREML Evaluates the maximum of the REML log-likelihood
%            based on REML estimates

% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 05-Nov-2013 22:29:13

const = (n-p) * log(2*pi);
logdetR = n * log(sig2(end));
logdetG =  dimRE' * log(sig2(1:nRE));
logdetH = sum(2*log(diag(Hfac)));
yPysig2 = (n-p);
loglik_REML = (-1/2) * (const + logdetR + logdetG + logdetH + yPysig2);

% Information criteria
D = (nRE + 1);
N = dimRE(1);
AIC = -2*loglik_REML + 2*D;
BIC = -2*loglik_REML + D*log(N);

end     % END of loglik_REML

%% Function maxloglikFun
function fun = maxloglikFun(sig2,H0,H0REDiag,H0REDiagInd,ind,nRE,n,dimRE,p,sig2REML)
%maxloglikFun evaluates the 'REML maximum log-likelihood function'
% fun = -(1/2) * (const+log(det(R))+log(det(G))+log(det(H))+(n-p))
%
% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 09-Jan-2014 21:55:07

%% start maxloglikFun
for inx = 1:nRE
    H0(H0REDiagInd(ind{inx})) = H0REDiag(ind{inx}) + ...
        1/(sig2(inx)/sig2(end));
end

for inx = 1:nRE
    H0(H0REDiagInd(ind{inx})) = H0REDiag(ind{inx}) + ...
        1/(sig2(inx)/sig2REML(end));
end

H0 = H0 / sig2REML(end);
const = (n-p) * log(2*pi);
logdetR = n * log(sig2(end));
logdetG =  dimRE * log(sig2(1:nRE));
logdetH = sum(2*log(diag(chol(H0))));
yPysig2 = (n-p) * sig2REML(end) / sig2(end);
fun = (-1/2) * (const + logdetR + logdetG + logdetH + yPysig2);

end
%% Function GRADIENT
function g = gradient(fun,x,delta,idx)
%gradient - The numerical gradient by finite differences.
%
% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 09-Jan-2014 21:55:07

narginchk(2,4);
p = length(x);
if nargin < 4, idx = []; end
if nargin < 3, delta = []; end
if isempty(idx), idx = 1:p; end
if isempty(delta), delta = eps^(1/3); end

g = zeros(p,1);
for i = idx
    x_minus = x; x_minus(i) = x_minus(i) - delta;
    x_plus = x; x_plus(i) = x_plus(i) + delta;
    g(i) = (fun(x_plus) - fun(x_minus))/2/delta;
end
end     % END of gradient