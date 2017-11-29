function result = getStats(Lambda,lmefit,options)
%getStats computes the basic t-statistics, p-values, (1-alpha) confidence
% intervals, for the linear functions of the fixed and the random effects
% coefficients, w = Lambda * [u;beta], by using the SATTERTHWAITE's
% denominator degrees of freedom (DDF).
%
% SYNTAX:
% STAT = getStats(Lambda,lmefit,options)
%
% INPUTS:
%  Lambda   - coefficient matrix
%  lmefit   - lmefit, fitted linear mixed effects model by HPMIXED
%  options. - structure of options with the following default values:
%         options.STAT.verbose = true
%         options.STAT.alpha = 0.05
%         options.STAT.Description = ''
%         options.STAT.colnames = 'fit_'
%         options.STAT.effectID = []
%         options.STAT.ddfMethod = 'Satterthwaite'
%
% EXAMPLE 0: (Fit and analysis of the SplitPlotData by using HPMIXED)
%   clear
%   load dsSplitPlotData;
%   formula  = 'y ~ A + B + A:B + (1 | Block) + (1 | Block:A)';
%   model = hpmixedmodel(SplitPlotData,formula);
%   lmefit = hpmixed(model);
%
% EXAMPLE 1: (Get statistics for FIXED and RANDOM effects and FITTED values)
%   STAT_FE = getStats('fixed',lmefit)
%   STAT_RE = getStats('random',lmefit)
%   STAT_FIT = getStats('fitted',lmefit)
%
% EXAMPLE 2: (Get t-statistics for the LSMEANS CONTRASTS for interaction A:B)
%   options.STAT.inference = 'contrasts';
%   options.STAT.inferenceSpace = 'broad';
%   options.STAT.ddfMethod = 'Satterthwaite';
%   [Lambda,options]  = getLambda({'A' 'B' },model,SplitPlotData,options);
%   STAT_AB_CONTRASTS = getStats(Lambda,lmefit,options)
%
% EXAMPLE 3: (FTEST for the LSMEANS CONTRASTS for interaction A:B)
%   options.STAT.inference = 'contrasts';
%   options.STAT.inferenceSpace = 'broad';
%   options.STAT.ddfMethod = 'FaiCornelius';
%   [Lambda,options]  = getLambda({'A' 'B' },model,SplitPlotData,options);
%   STAT_AB_FTEST = getStats(Lambda,lmefit,options)
%
% EXAMPLE 4: (LSMEANS PAIRWISE COMPARISONS for the effect A)
%   options.STAT.inference = 'pairwise';
%   options.STAT.inferenceSpace = 'narrow';
%   options.STAT.ddfMethod = 'Satterthwaite';
%   [Lambda,options]  = getLambda({'A'},model,SplitPlotData,options);
%   STAT_A_COMPARISONS = getStats(Lambda,lmefit,options)
%
% EXAMPLE 5: (Get statistics for the LSMEANS for the effect A)
%   options.STAT.inference = 'lsmeans';
%   options.STAT.inferenceSpace = 'broad';
%   [Lambda,options]  = getLambda({'A'},model,SplitPlotData,options);
%   STAT_A_LSMEANS = getStats(Lambda,lmefit,options)
%
% EXAMPLE_6: (BROAD INFERENCE SPACE for LSMEANS of the factor Block)
%   options.STAT.inference = 'lsmeans';
%   options.STAT.inferenceSpace = 'broad';
%   [Lambda,options]  = getLambda({'Block'},model,SplitPlotData,options);
%   options.STAT.alpha = 0.01;
%   STAT_Block_LSMEANS_BROAD = getStats(Lambda,lmefit,options)
%
% EXAMPLE 7: (NARROW INFERENCE SPACE for LSMEANS of the factor Block)
%   options.STAT.inference = 'lsmeans';
%   options.STAT.inferenceSpace = 'narrow';
%   [Lambda,options]  = getLambda({'Block'},model,SplitPlotData,options);
%   options.STAT.alpha = 0.01;
%   STAT_Block_LSMEANS_NARROW = getStats(Lambda,lmefit,options)
%
% EXAMPLE 8: (ITERMEDIATE INFERENCE SPACE for LSMEANS of the factor Block)
%   options.STAT.inference = 'lsmeans';
%   options.STAT.inferenceSpace = 'intermediate';
%   options.STAT.includedXcols = (1:size(model.X,2));
%   options.STAT.includedZcols = model.dimRE(1) + (1:model.dimRE(2));
%   [Lambda,options]  = getLambda({'Block'},model,SplitPlotData,options);
%   options.STAT.alpha = 0.01;
%   STAT_Block_LSMEANS_ITERMEDIATE = getStats(Lambda,lmefit,options)


% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 12-Jan-2014 11:41:21

%% CHECK INPUTS / OUTPUTS
narginchk(2,3);
if nargin < 3, options = []; end

%% SET the default options
if ~isfield(options, 'STAT'), options.STAT = []; end
if ~isfield(options.STAT, 'verbose'), options.STAT.verbose = true; end
if ~isfield(options.STAT, 'alpha'), options.STAT.alpha = 0.05; end
if ~isfield(options.STAT, 'Description'), options.STAT.Description = ''; end
if ~isfield(options.STAT, 'colnames'), options.STAT.colnames = 'fit_'; end
if ~isfield(options.STAT, 'effectID'), options.STAT.effectID = []; end
if ~isfield(options.STAT, 'ddfMethod'), options.STAT.ddfMethod = 'Satterthwaite'; end
if ~isfield(options.STAT, 'isFitConditional'), options.STAT.isFitConditional = true; end

%% Get the required information from the fittetd model (lmefit)
sig2        = lmefit.varianceComponents.Estimates;
stdVC       = lmefit.varianceComponents.StandardErors;
ubeta       = lmefit.Details.Matrices.ubeta;
I_REML      = lmefit.Details.Matrices.I_REML;
Hfac        = lmefit.Details.Matrices.H_factor;
perm        = lmefit.Details.Matrices.Perm;
invperm     = lmefit.Details.Matrices.InvPerm;
dimRE       = lmefit.ModelInfo.randomEffectsDim;
nRE         = lmefit.ModelInfo.randomEffectsNum;
p           = lmefit.ModelInfo.fixedEffectsLength;
q           = lmefit.ModelInfo.randomEffectsLength;
m           = lmefit.ModelInfo.MMEmatRows;
n           = lmefit.ModelInfo.observationsNum;

isVarianceComponentsTable = false;
%% Set the coefficient matrix Lambda for standard cases
if isempty(Lambda)
    Lambda = speye(m);
    options.STAT.colnames = 'REandFEcoef_';
    options.STAT.Description = 'Estimated RANDOM and FIXED effects';
elseif ischar(Lambda)
    switch lower(Lambda)
        case {'fixed' 'fe'}
            Lambda = sparse(q + (1:p),1:p,1,m,p);
            if isempty(lmefit.fixedEffects.Names)
                options.STAT.colnames = 'FEcoef_';
            else
                options.STAT.colnames = lmefit.fixedEffects.Names;
            end
            options.STAT.Description = 'Estimated FIXED effects';
        case {'random' 're'}
            if isempty(options.STAT.effectID)
                Lambda = sparse(1:q,1:q,1,m,q);
            else
                dim = dimRE(options.STAT.effectID);
                if options.STAT.effectID ==1
                    idx0 = 1;
                else
                    idx0 = dimRE(options.STAT.effectID-1);
                end
                idx = 1:dim;
                Lambda = sparse(idx0+idx,idx,1,m,dim);
            end
            options.STAT.colnames = ...
                ['REcoef ',num2str(options.STAT.effectID),'_'];
            options.STAT.Description = 'Estimated RANDOM effects';
        case {'all' 'coefs' 'coefficients'}
            Lambda = sparse(1:m,1:m,1,m,m);
            options.STAT.colnames = 'REandFEcoef_';
            options.STAT.Description = 'Estimated RANDOM and FIXED effects';
        case {'fitted' 'fit' 'yhat'}
            try
                if options.STAT.isFitConditional
                    Lambda = [lmefit.Details.Model.Z ...
                        sparse(lmefit.Details.Model.X)]';
                else
                    Lambda = [sparse(n,q) ...
                        sparse(lmefit.Details.Model.X)]';
                end
            catch
                error('VW:hpmixed:getStats', ...
                    'MISSING INPUT: lmefit.Details.Model ... ')
            end
            options.STAT.colnames = 'yhat_';
            options.STAT.Description = 'FITTED yhat';
        case {'variancecomponents' 'vc'}
            isVarianceComponentsTable = true;
            options.STAT.ddfMethod = 'none';
            options.STAT.colnames = 'sigma2_';
            options.STAT.Description = 'FITTED variance components';
            STAT = createTTESTtable(sig2,stdVC,Inf,nRE+1,options);            
        otherwise
            error('VW:hpmixed:getStats', ...
                'WRONG INPUT: "Lambda" TRY: getLambda() ...')
    end
end

%% ALGORITHM
mse = [];
if ~isVarianceComponentsTable  
    if isempty(I_REML)        
        I_REML = getFisherInfApprox(lmefit);
        % I_REML_used = 'Fisher Information Matrix / Approx';
    end
    
    switch lower(options.STAT.ddfMethod)
        case {'residual' 'res'}
            [estimate,std,ddf,ndf] = getStatssAndDDF(Lambda,ubeta,sig2,...
                Hfac,perm,n,p,[],[],[],[],[],'residual');
            testType = 'Ttest';
        case {'satterthwaite' 'sat'}
            [estimate,std,ddf,ndf] = getStatssAndDDF(Lambda,ubeta,sig2, ...
                Hfac,perm,n,p,I_REML,m,dimRE,nRE,invperm,'satterthwaite');
            testType = 'Ttest';
        case {'faicornelius' 'fc' 'fai' 'fai-cornelius' 'cornelius'}
            [fstat,~,ddf,ndf,mse] = getStatssAndDDF(Lambda,ubeta,sig2, ...
                Hfac,perm,n,p,I_REML,m,dimRE,nRE,invperm,'faicornelius');
            testType = 'Ftest';
        otherwise
            [estimate,std,ddf,ndf] = getStatssAndDDF(Lambda,ubeta,sig2,...
                Hfac,perm,[],[],[],[],[],[],[],'none');
            testType = 'Ttest';
    end
    
    switch lower(testType)
        case 'ttest'
            STAT = createTTESTtable(estimate,std,ddf,ndf,options);
        case 'ftest'
            [STAT,mse] = createFTESTtable(fstat,mse,ddf,ndf,options);
    end
end


%% RESULT
result.TABLE = STAT;
result.Details.Lambda = Lambda;
result.Details.mse = mse;
result.options = options.STAT;

if options.STAT.verbose
    disp(STAT)
    disp(options.STAT)
end

end     % END of getStats

%% FUNCION getStatsistics
function [estimate,std,ddf,ndf,mse] = getStatssAndDDF(Lambda,ubeta,sig2, ...
    Hfac,perm,n,p,I_REML,m,dimRE,nRE,invperm,method)
% getStatssAndDDF Computes the estimates and the degrees of freedom (DDF) by
%                the required method. Supported DDF methods are: 'none',
%                'residual', 'satterthwaite', 'faicornelius'.

% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 12-Jan-2014 11:41:21

switch method
    case {'faicornelius'}
        ndf = size(Lambda,2);
        mse = Hfac' \ Lambda(perm,:);
        mse = full(mse' * mse);
        estimate = Lambda(perm,:)' * ubeta;
        estimate = (estimate' * (mse \ estimate)) / ndf;
        [U,S] = svd(mse);
        Lambda = Lambda * U;
        sse = diag(S);
        std = [];
    otherwise
        estimate = Lambda(perm,:)' * ubeta;
        sse = full(sum((Hfac' \ Lambda(perm,:)).^2))';
        std = sqrt(sse);
        ndf = size(Lambda,2);
        mse = [];
end

switch method
    case 'none'
        ddf = Inf;
    case 'residual'
        ddf = (n - p);
    case {'satterthwaite' 'faicornelius'}
        ddf = zeros(size(sse));
        SigmaVC = I_REML \ eye(length(sig2));
        idx = cell(nRE,1);
        startidx = 0;
        for ii = 1:nRE
            idx{ii} = startidx + (1:dimRE(ii));
            startidx = startidx + dimRE(ii);
        end
        
        for ii = 1:nRE
            idx{ii} = invperm(idx{ii});
        end
        
        Lambda = Lambda(perm,:);
        SqrtmHatLambda = (Hfac' \ Lambda);
        
        diagSqrtmGinvDer1 = zeros(m,1);
        for ii = 1:nRE
            diagSqrtmGinvDer1(idx{ii}) = - 1 / sig2(ii);
        end
        SqrtmGinvDer1 = spdiags(diagSqrtmGinvDer1,0,m,m);
        for ii = 1:ndf
            SqrtmGinvDer1CLambda = SqrtmGinvDer1 * (Hfac \ SqrtmHatLambda(:,ii));
            %SqrtmGinvDer1CLambda = (SqrtmGinvDer1 * Cfac) * SqrtmHatLambda(:,ii);
            grad = zeros(nRE+1,1);
            gradsum = 0;
            for jj = 1:nRE
                % grad_j = Lambda' * C * GinvDer1_j * C * Lambda
                grad(jj) = sum(nonzeros(SqrtmGinvDer1CLambda(idx{jj})).^2);
                % gradsum = Lambda' * C * Ginv * C * Lambda
                gradsum = gradsum + grad(jj) * sig2(jj);
            end
            % here, sse = Lambda'*C*Lambda
            % sse(ii) = sum(nonzeros(SqrtmHatLambda(:,ii)).^2);
            % grad(nRE+1) = (Lambda'*C*Lambda - Lambda'*C*Ginv*C*Lambda) / sig2(nRE+1)
            grad(nRE+1) = (sse(ii) -  gradsum) / sig2(nRE+1);
            ddf(ii) = 2 * sse(ii)^2 / (grad'*SigmaVC*grad);
        end
        if strcmp(method,'faicornelius')
            E = sum((ddf > 0) .* ddf ./ (ddf - 2));
            ddf = 2 * E ./ (E - ndf);
        end
end
end     % END of getStatssAndDDF

%% Create the T-test STAT table
function STAT = createTTESTtable(estimate,std,ddf,ndf,options)
% createTTESTtable Creates a table with basic t-statistics. The results
%     are presented as a dataset STAT.

% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 12-Jan-2014 11:41:21

%% Set the variables and create the table
alpha = options.STAT.alpha;
if size(options.STAT.colnames,1) == 1 && ndf > 1
    colnames =  strcat({options.STAT.colnames},num2str((1:ndf)','%-d'));
else
    colnames = options.STAT.colnames;
end

stat        = estimate ./ std;
quantile    = tinv(1-alpha/2,ddf);
pval        = 2 * tcdf(-abs(stat),ddf);
low         = estimate - quantile .* std;
upp         = estimate + quantile .* std;

if size(ddf,1) == 1
    ddf = ddf * ones(size(estimate));
    quantile = quantile * ones(size(estimate));
end

VarNames = {'Estimate' 'SE' 'tStat' 'DF' 'Quantile' ...
    'pValue' 'Lower' 'Upper'};
STAT = dataset(estimate, std, stat, ddf, quantile, pval, low, upp, ...
    'VarNames',VarNames,'ObsNames',colnames);

STAT.Properties.Description = options.STAT.Description;
end     % END of createTTESTtable

%% Create the F-test STAT table
function [STAT,mse] = createFTESTtable(fstat,mse,ddf,ndf,options)
% createTTESTtable Creates a table with basic F-statistics. The results
%     are presented as a dataset STAT.

% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 12-Jan-2014 11:41:21

%% Set the variables and create the table
alpha = options.STAT.alpha;
%if ~isfield(options.STAT, 'grpname'), options.STAT.grpname = 'GRP'; end
colnames = options.STAT.colnames;

quantile    = finv(1-alpha,ndf,ddf);
pval        = 1 - fcdf(fstat,ndf,ddf);

if size(ddf,1) == 1
    ddf = ddf * ones(size(fstat));
    quantile = quantile * ones(size(fstat));
end

VarNames = {'Fstat' 'NDF' 'DDF' 'Quantile' 'pValue'};
% STAT = dataset(fstat, ndf, ddf, quantile, pval, ...
%     'VarNames',VarNames,'ObsNames',colnames);

STAT = dataset(fstat, ndf, ddf, quantile, pval, ...
    'VarNames',VarNames);

STAT.Properties.Description = options.STAT.Description;
end     % END of createTTESTtable