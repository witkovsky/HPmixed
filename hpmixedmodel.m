function model = hpmixedmodel(ds, formula, options)
%HPMIXEDMODEL  Creates linear mixed effects model structure suitable
%  for the algorithm *hpmixed* from a given dataset and formula.
%
%  Syntax:
%  model = hpmixedmodel(ds, formula)
%
%  EXAMPLE
%  load dsSplitPlotData
%  formula  = 'y ~ A + B + A:B + (1 | Block) + (1 | Block:A)';
%  model = hpmixedmodel(SplitPlotData,formula);
%
%  This is an experimantal version, adapted from MATLABs LinearMixedModel:
%  Based on Statistics Toolbox (version 8.3 (R2013b)):
%  classreg.regr.modelutils.designmatrix
%  classreg.regr.LinearMixedFormula
%
%  Requires Statistics Toolbox (version 8.3 (R2013b)).

% Ver.: 25-Oct-2013 12:18:46

%% CHECK INPUTS / OUTPUTS
narginchk(2,3);

if nargin < 3, options = []; end

if ~isfield(options, 'verbose'), options.verbose = true; end
%if ~isfield(options, 'dummyVarCode'), options.dummyVarCode = 'full'; end
%if ~isfield(options, 'dummyVarCode'), options.dummyVarCode = 'reference'; end
if ~isfield(options, 'dummyVarCode'), options.dummyVarCode = 'effects'; end

F = classreg.regr.LinearMixedFormula(formula,ds.Properties.VarNames);
model.Formula = F;

model.y = ds.(F.ResponseName);
model.X = [];
model.Z = [];
model.dimRE = [];

[includedFormulaVars,DatasetsVarLocation] = ...
    ismember(ds.Properties.VarNames,F.VariableNames);
if ~all(includedFormulaVars)
    ds = ds(:,includedFormulaVars);
end

VarNames = ds.Properties.VarNames;
nVars = length(VarNames);
IsCategorical = zeros(1,nVars);
Range = cell(1,nVars);
for i = 1:nVars
    IsCategorical(i) = iscategorical(eval(['ds.',VarNames{i}]));
    if IsCategorical(i)
        Range{i} = double(getlevels(eval(['ds.',VarNames{i}])));
    else
        Range{i} = [min(eval(['ds.',VarNames{i}])) ...
            max(eval(['ds.',VarNames{i}]))];
    end
end

%% Matrix X by classreg.regr.modelutils.designmatrix
[model.X,~,~,XCols2Terms,XColNames] = ...
    classreg.regr.modelutils.designmatrix(ds,...
    'Model',F.FELinearFormula.Terms(:,DatasetsVarLocation),...
    'PredictorVars',F.FELinearFormula.PredictorNames,...
    'ResponseVar',F.ResponseName,...
    'DummyVarCoding',options.dummyVarCode);

if options.verbose
    model.FixedInfo.XCols2Terms = XCols2Terms;
    model.FixedInfo.XColNames = XColNames;
end

%% Groups
R = length(model.Formula.RELinearFormula);
G = cell(R,1);
GNames = cell(R,1);
for i = 1:R
    intVars = model.Formula.GroupingVariableNames{i};
    [G{i},GNames{i}] = intvars(ds,intVars);
end

Gid = cell(R,1);
GidLevelNames = cell(R,1);
lev = zeros(R,1);
for i = 1:R
    [Gid{i},GidLevelNames{i}] = grp2idx(G{i});
    lev(i) = length(GidLevelNames{i});
end

if options.verbose
    model.GroupingInfo.R = R;
    model.GroupingInfo.G = G;
    model.GroupingInfo.GNames = GNames;
    model.GroupingInfo.Gid = Gid;
    model.GroupingInfo.GidLevelNames = GidLevelNames;
    model.GroupingInfo.lev = lev;
end

%% Matrix Z by classreg.regr.modelutils.designmatrix
numRE = length(F.RELinearFormula);
Z = cell(numRE,1);
ZColNames = cell(numRE,1);
ZColGrps = cell(numRE,1);
q = zeros(numRE,1);
for i = 1:numRE
    [Z{i},~,~,ZColGrps{i},ZColNames{i}] = ...
        classreg.regr.modelutils.designmatrix(ds,...
        'Model',F.RELinearFormula{i}.Terms(:,DatasetsVarLocation),...
        'DummyVarCoding','reference', ...
        'CategoricalVars',logical(IsCategorical(DatasetsVarLocation)), ...
        'CategoricalLevels',Range{DatasetsVarLocation(i)});
    q(i) = size(Z{i},2);
end

if options.verbose
    model.RandomInfo.Z = Z;
    model.RandomInfo.ZColNames = ZColNames;
    model.RandomInfo.ZColGrps = ZColGrps;
    model.RandomInfo.q = q;
    model.RandomInfo.numRE = numRE;
end
model.Z = sparseZ(Z,q,lev,Gid);
model.dimRE = lev;

% Most of the information is in Model.Formula
model.VarInfo.nVars = nVars;
model.VarInfo.VarNames = VarNames;
model.VarInfo.VarNames = VarNames;
model.VarInfo.IsCategorical = IsCategorical;
model.VarInfo.DatasetsVarLocation = DatasetsVarLocation;
model.VarInfo.ResponseName = F.ResponseName;
model.VarInfo.PredictorNames = F.FELinearFormula.PredictorNames;
model.VarInfo.XColNames = XColNames;
model.VarInfo.XCols2Terms = XCols2Terms;
model.VarInfo.Terms = F.FELinearFormula.Terms(:,DatasetsVarLocation);
model.VarInfo.TermNames = model.Formula.FELinearFormula.TermNames;

end

%% Function intvars
function [G,GName] = intvars(ds,intVars)

assert( iscellstr(intVars) );
assert( all(ismember(intVars,ds.Properties.VarNames)) );

k = length(intVars);
if k >= 1
    G = nominal(ds.(intVars{1}));
    GName = intVars{1};
    for i = 2:k
        G = G.*nominal(ds.(intVars{i}));
        GName = [GName,':',intVars{i}];
    end
else
    G = [];
    GName = [];
end

G = removecats(G);
end

%% Function sparseZ
function Zs = sparseZ(Z,q,lev,Gid)

qlev = q .* lev;
R = length(Gid);
N = size(Z{1},1);
Zs = sparse(N,sum(qlev));
for r = 1:R
    for k = 1:lev(r)
        idx = ( Gid{r} == k );
        offset = sum(qlev(1:(r-1))) + (k-1)*q(r);
        Zs(idx, offset + 1 : offset + q(r)) = Z{r}(idx,:);
    end
end
end