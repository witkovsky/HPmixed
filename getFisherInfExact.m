function [FI_REML,Tfac] = getFisherInfExact(lmefit,options)
% getFisherInfExact Computes the 'exact' Fisher information matrix
%   FI_REML for the REML estimator of variance components in linear mixed
%   model with simple variance components (VC) structure.
%   The Fisher Information matrix FI_REML is evaluated from the
%   fitted linear mixed model(LMEFIT by HPMIXED).
%
% SYNTAX:
% FI_REML = getFisherInfExact(lmefit,options)
%
% EXAMPLE:
%   load dsSplitPlotData;
%   options.verbose = false;
%   lmefit = hpmixed(model,options);
%   FI_REML = getFisherInfExact(lmefit)

% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 14-Jan-2014 19:46:30

%% CHECK INPUTS / OUTPUTS
narginchk(1,2);
if nargin < 2, options = []; end

%% Set the default options
if ~isfield(options, 'FImethod'), options.FImethod = 'exact'; end

%% Get the required information from the fitetd model, lmefit
sig2        = lmefit.varianceComponents.Estimates;
Hfac        = lmefit.Details.Matrices.H_factor;
Gd          = lmefit.Details.Matrices.Gd;
n           = lmefit.ModelInfo.observationsNum;
p           = lmefit.ModelInfo.fixedEffectsLength;
q           = lmefit.ModelInfo.randomEffectsLength;
dimRE       = lmefit.ModelInfo.randomEffectsDim;
nRE         = lmefit.ModelInfo.randomEffectsNum;

%% Calculate the exact Fisher Information Matrix for REML of VC    
Tfac = cell(nRE,1);
swptol = 1e-11;
for i = 1:nRE
    Tfac{i} = (round((Gd{i} / Hfac)/swptol)*swptol) / sqrt(sig2(i));
end

trT = zeros(nRE,1);
trTT = zeros(nRE);
for i = 1:nRE
    trT(i) = sum(nonzeros(Tfac{i}).^2);
    trTT(i,i) = blockTrace(Tfac{i});
    for j = (i+1):nRE
        trTT(i,j) = blockTrace(Tfac{i},Tfac{j});
        trTT(j,i) = trTT(i,j);
    end
end

% Evaluate the Fisher Information matrix
FI_REML = zeros(nRE+1);
trcolsum = sum(trTT);
trallsum = sum(trcolsum);
for i = 1:nRE
    FI_REML(i,i) = ((dimRE(i) - 2*trT(i)) + trTT(i,i)) / sig2(i)^2;
    for j = (i+1):nRE
        FI_REML(i,j) = trTT(i,j) / (sig2(i)*sig2(j));
        FI_REML(j,i) = FI_REML(i,j);
    end
    FI_REML(nRE+1,i) = (trT(i)-trcolsum(i)) / (sig2(nRE+1)*sig2(i));
    FI_REML(i,nRE+1) = FI_REML(nRE+1,i);
end
FI_REML(nRE+1,nRE+1) = ((n-p) - q + trallsum) / sig2(nRE+1)^2;
FI_REML = FI_REML/2;

end
%% FUNCTION blockTrace
function tr = blockTrace(T1,T2)
%blockTrace Computes the trace(T1*T2') where T1 and T2 are sparse
% matrices of dimensions (dimT1 x m) and (dimT2 x m)

% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 09-Nov-2013 15:14:26

dimT1 = size(T1,1);
tr = 0;
if nargin < 2
    for i = 1:dimT1
        row = T1(i,:);
        tr = tr + sum(nonzeros(T1*row').^2);
    end
else
    dimT2 = size(T2,1);
    if dimT1 <= dimT2
        for i = 1:dimT1
            row = T1(i,:);
            tr = tr + sum(nonzeros(T2*row').^2);
        end
    else
        for i = 1:dimT2
            row = T2(i,:);
            tr = tr + sum(nonzeros(T1*row').^2);
        end
    end
end
end