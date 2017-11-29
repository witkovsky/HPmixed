function FI_REML = getFisherInfApprox(lmefit,options)
% getFisherInfApprox The numerical approximation of the Fisher
%   information matrix FI for the REML estimator of variance components in
%   linear mixed models with linear parametrization of the variance -
%   covariance matrix V = Z'*G*Z + R. The Fisher Information matrix FI_REML
%   is Hessian of the 'REML maximum log-likelihood function', which is
%   evaluated from the fitted linear mixed model(LMEFIT by HPMIXED).
%
% SYNTAX:
% FI_REML = getFisherInfApprox(lmefit,options)
%
% EXAMPLE:
%   load dsSplitPlotData;
%   options.verbose = false;
%   lmefit = hpmixed(model,options);
%   FI_REML = getFisherInfApprox(lmefit)

% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 09-Jan-2014 21:55:07

%% CHECK INPUTS / OUTPUTS
narginchk(1,2);
if nargin < 2, options = []; end

% Default options / Currently only REML
if ~isfield(options, 'FitMethod')
    options.FitMethod = 'REML';
end

if ~isfield(options, 'delta')
    options.delta = eps^(1/4);
end

%% Get the required information from the fitetd model, lmefit
sig2        = lmefit.varianceComponents.Estimates;
H0          = lmefit.Details.Matrices.H0;
H0REDiag    = lmefit.Details.Matrices.H0REDiag;
H0REDiagInd = lmefit.Details.Matrices.H0REDiagInd;
ind         = lmefit.Details.Matrices.ind;
n           = lmefit.ModelInfo.observationsNum;
dimRE       = lmefit.ModelInfo.randomEffectsDim;
nRE         = lmefit.ModelInfo.randomEffectsNum;

%% Calculate the Fisher Information Matrix
delta = options.delta;
FI_REML = hessian(@(sig2) ... 
    maxloglikFun(sig2,H0,H0REDiag,H0REDiagInd,ind,nRE,n,dimRE),sig2,delta);

%% Function maxloglikFun
    function fun = maxloglikFun(sig2,H0,H0REDiag,H0REDiagInd,ind,nRE,n,dimRE)
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

        H0 = H0 / sig2(end);
        logdetR = n * log(sig2(end));
        logdetG =  dimRE * log(sig2(1:nRE));
        logdetH = sum(2*log(diag(chol(H0))));
        fun = (-1/2) * (logdetR + logdetG + logdetH);
    end

%% Function hessian
    function H = hessian(fun,x,delta)
    %hessian - The numerical Hessian by finite differences.
    %
    % (c) Viktor Witkovsky (witkovsky@savba.sk)
    % Ver.: 09-Jan-2014 21:55:07
        
    %% start hessian
        narginchk(2,3);
        if nargin < 3, delta = []; end
        if isempty(delta), delta = eps^(1/4); end
        if isscalar(delta), delta = delta * ones(size(x)); end
        % Correction for small (positive) x
        delta((x-delta) < 0) = x((x-delta) < 0) / 2;        
        
        p = length(x);
        H = zeros(p,p);
        for i = 1:p
            x_minus = x; x_minus(i) = x_minus(i) - delta(i);          
            x_plus = x; x_plus(i) = x_plus(i) + delta(i);          
            H(:,i) = (gradient(fun, x_plus, delta, i) - ...
                gradient(fun, x_minus, delta, i))/2/delta(i);
        end
        H = H + tril(H,-1)';
    end

%% Function gradient
    function g = gradient(fun,x,delta,idx)
    %gradient - The numerical gradient by finite differences.
    %
    % (c) Viktor Witkovsky (witkovsky@savba.sk)
    % Ver.: 09-Jan-2014 21:55:07
        
    %% start gradient
        narginchk(2,4);
        if nargin < 4, idx = []; end
        if nargin < 3, delta = []; end
        if isempty(idx), idx = 1; end
        if isempty(delta), delta = eps^(1/3); end
        if isscalar(delta), delta = delta * ones(size(x)); end
        % Correction for small (positive) x
        delta((x-delta) < 0) = x((x-delta) < 0) / 2;  
        
        p = length(x);
        g = zeros(p,1);
        for i = idx:p
            x_minus = x; x_minus(i) = x_minus(i) - delta(i);          
            x_plus = x; x_plus(i) = x_plus(i) + delta(i);            
            g(i) = (fun(x_plus) - fun(x_minus))/2/delta(i);
        end
    end

end     % END of getFisherInfApprox