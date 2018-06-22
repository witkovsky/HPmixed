function [Lambda,options] = getLambda(varNames,model,dataset,options)
% getLambda creates the coefficient matrix Lambda for LSMEANS and their
% CONTRASTS and PAIRWISE contrasts for computing the basic t-statistics
% resp. F-statistics and their p-values by using the Satterthwaite's
% approximation of the degrees of freedom.
%
% EXAMPLE: (LSMEANS contrasts for interaction A:B)
%   load dsSplitPlotData;
%   formula  = 'y ~ A + B + A:B + (1 | Block) + (1 | Block:A)';
%   model = hpmixedmodel(SplitPlotData,formula);
%
%   options.STAT.inference = 'contrasts';
%   options.STAT.inferenceSpace = 'broad';
%   [Lambda,options]  = getLambda({'A' 'B' },model,SplitPlotData,options);
%
%   options.verbose = false;
%   lmefit = hpmixed(model,options);
%   STAT = getStats(Lambda,lmefit,options)

% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 11-Jan-2014 19:56:44

%% CHECK INPUTS / OUTPUTS
narginchk(3,4);
if nargin < 4, options = []; end

%% SET the default options

if ~isfield(options, 'STAT')
    options.STAT = []; 
end

if ~isfield(options.STAT, 'inference')
    options.STAT.inference = 'LSmeans'; 
end

if ~isfield(options.STAT, 'inferenceSpace') 
    options.STAT.inferenceSpace = 'broad'; 
end

if ~isfield(options.STAT, 'includedZcols')    
    options.STAT.includedZcols = []; 
end

if ~isfield(options.STAT, 'includedXcols')    
    options.STAT.includedXcols = []; 
end

%% Create the coefficient matrix Lambda
[n,p] = size(model.X);
q = size(model.Z,2);
m = p + q;

nVars = length(varNames);
var = dataset.(varNames{1});
grpName = varNames{1};
for i = 2:nVars
    var = var .* dataset.(varNames{i});
    grpName = [grpName,':',varNames{i}];
end
options.STAT.VarDescription = grpName;

numVar = size(var,1);
[levels,~,j] = unique(var);
labels = getlabels(levels);
numLevels = size(levels,1);
Lambda = sparse(1:numVar,j,1,numVar,numLevels);
replications = sum(Lambda)';
Lambda = Lambda * spdiags(1./replications,0,numLevels,numLevels);

switch lower(options.STAT.inferenceSpace)
    case 'broad'
        Lambda = sparse([sparse(n,q) model.X])' * Lambda;
        Description = 'INFERENCE SPACE: Broad, ';
    case 'narrow'
        Lambda = sparse([model.Z model.X])' * Lambda;
        Description = 'INFERENCE SPACE: Narrow, ';
    case {'intermediate' 'inter'}        
        ZX = sparse(n,m);
        ZX(:,options.STAT.includedZcols) = ...
            model.Z(:,options.STAT.includedZcols);
        ZX(:,q + options.STAT.includedXcols) = ...
            model.X(:,options.STAT.includedXcols);
        Lambda = ZX' * Lambda;
        Description = 'INFERENCE SPACE: Intermediate, User INCLUDED columns, ';
    otherwise
        Lambda = sparse([sparse(n,q) model.X])' * Lambda;
        Description = 'INFERENCE SPACE: Broad, ';
end
        
switch lower(options.STAT.inference)
    case {'lsmeans' 'means' 'ls'}
        colNames = cell(numLevels,1);
        for i = 1:numLevels
            colNames{i} = [grpName,'_',labels{i}];
        end
        Description = ['STAT TABLE: LSmeans, ', Description];
    case {'contrasts'}
        contrasts = [ones(numLevels-1,1) (-1)*speye(numLevels-1)];
        Lambda = Lambda * contrasts';
        colNames = cell(numLevels-1,1);
        Description = ['STAT TABLE: LSMEANS contrasts, ', Description];
        for i = 2:numLevels
            colNames{i-1} = [grpName,'_',labels{1},' - ',grpName,'_',labels{i}];
        end
    case {'pairwise' 'pairs' 'comparisons'}
        nPairs = numLevels*(numLevels-1)/2;
        colNames = cell(nPairs,1);
        ind = 0;
        indI = zeros(nPairs,1);
        indJ = zeros(nPairs,1);
        for i = 1:numLevels
            for j = (i+1):numLevels
                ind = ind + 1;
                indI(ind) = i;
                indJ(ind) = j;
                colNames{ind} = [grpName,'_',labels{i},' - ',grpName,'_',labels{j}];
            end
        end
        pairs = sparse(1:nPairs,indI,1,nPairs,numLevels);
        pairs = pairs - sparse(1:nPairs,indJ,1,nPairs,numLevels);
        Lambda = Lambda * pairs';
        Description = ['STAT TABLE: Pairwise comparisons, ', Description];
    otherwise
        colNames = cell(numLevels,1);
        for i = 1:numLevels
            colNames{i} = [grpName,'_',labels{i}];
        end
        Description = ['STAT TABLE: LSmeans, ', Description];  
end

Description = [Description,' ', 'VARIABLE: ',grpName];
options.STAT.Description = Description;
options.STAT.colnames = colNames;
options.STAT.grpname = grpName;

end     % END of getLambda