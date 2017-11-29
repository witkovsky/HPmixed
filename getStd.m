function [stdErr,Lambda] = getStd(lmefit,Lambda)
% getStd computes the standard errors of the functions of the fixed and
% the random effects, Lambda * [u;beta]
%
% SYNTAX:
% se = getStd(lmefit,Lambda)
%
% INPUTS:
%
% OUTPUTS:
%
% EXAMPLE: (Statistics for the fitted values yfit = X*beta + Z*u)
%   load dsSplitPlotData;
%   Lambda = [model.Z model.X]';
%   se = getStd(lmefit,Lambda)

% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 08-Jan-2014 12:15:14

%% CHECK INPUTS / OUTPUTS
narginchk(1,2);
if nargin < 2, Lambda = []; end

%% Get the required information from the fitetd model, lmefit
Hfac        = lmefit.Details.Matrices.H_factor;
perm        = lmefit.Details.Matrices.Perm;
p           = lmefit.ModelInfo.fixedEffectsLength;
q           = lmefit.ModelInfo.randomEffectsLength;
m           = p + q;

if isempty(Lambda)
    % Random + Fixed effects SE
    Lambda = speye(m);
elseif ischar(Lambda)
    if strcmp(Lambda,'fixed')
        % Fixed effects only
        Lambda = sparse(q + (1:p),1:p,1,m,p);
    elseif strcmp(Lambda,'random')
        % Random effects only
        Lambda = sparse(1:q,1:q,1,m,q);
    elseif strcmp(Lambda,'all')
        % All: Random and Fixed effects
        Lambda = sparse(1:m,1:m,1,m,m);
    end
end

%% Calculate the std for functions Lambda * [u;beta]
stdErr = sqrt(full(sum((Hfac' \ Lambda(perm,:)).^2)))';

end     % end of getStd