function model = hpmixedmodeltable(dataTable, formula, options)
%HPMIXEDMODEL Creates linear mixed-effects model structure suitable
% for the algorithm *hpmixed* from a given data table and formula.
%
% Syntax:
% model = hpmixedmodel(ds, formula)
%
%  EXAMPLE
%  load dsSplitPlotData
%  formula  = 'y ~ A + B + A:B + (1 | Block) + (1 | Block:A)';
%  model = hpmixedmodeltable(SplitPlotDataTable,formula);
%
% This is an experimental version, adapted from MATLAB's LinearMixedModel.

% Check inputs
narginchk(2, 3);
if nargin < 3, options = []; end
if ~isfield(options, 'verbose'), options.verbose = true; end
if ~isfield(options, 'dummyVarCode'), options.dummyVarCode = 'effects'; end

% Formula processing
F = classreg.regr.LinearMixedFormula(formula, dataTable.Properties.VariableNames);
model.Formula = F;

% Response variable
model.y = dataTable.(F.ResponseName);
model.X = [];
model.Z = [];
model.dimRE = [];

% Ensure included formula variables are part of the dataset
[includedFormulaVars, DatasetsVarLocation] = ...
    ismember(dataTable.Properties.VariableNames, F.VariableNames);
if ~all(includedFormulaVars)
    dataTable = dataTable(:, includedFormulaVars);
end

VarNames = dataTable.Properties.VariableNames;
nVars = length(VarNames);

% Determine categorical variables and their levels
IsCategorical = false(1, nVars);
Range = cell(1, nVars);
for i = 1:nVars
    if iscategorical(dataTable.(VarNames{i}))
        IsCategorical(i) = true;
        Range{i} = categories(dataTable.(VarNames{i}));
    else
        Range{i} = [min(dataTable.(VarNames{i})), max(dataTable.(VarNames{i}))];
    end
end

% Fixed effects design matrix
[model.X, ~, ~, XCols2Terms, XColNames] = ...
    classreg.regr.modelutils.designmatrix(dataTable, ...
    'Model', F.FELinearFormula.Terms(:, DatasetsVarLocation), ...
    'PredictorVars', F.FELinearFormula.PredictorNames, ...
    'ResponseVar', F.ResponseName, ...
    'DummyVarCoding', options.dummyVarCode);

if options.verbose
    model.FixedInfo.XCols2Terms = XCols2Terms;
    model.FixedInfo.XColNames = XColNames;
end

% Random effects groups
R = length(model.Formula.RELinearFormula);
G = cell(R, 1);
GNames = cell(R, 1);
for i = 1:R
    intVars = model.Formula.GroupingVariableNames{i};
    [G{i}, GNames{i}] = intvars(dataTable, intVars);
end

% Random effects structure
Gid = cell(R, 1);
GidLevelNames = cell(R, 1);
lev = zeros(R, 1);
for i = 1:R
    [Gid{i}, GidLevelNames{i}] = grp2idx(G{i});
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

% Random effects design matrix
numRE = length(F.RELinearFormula);
Z = cell(numRE, 1);
ZColNames = cell(numRE, 1);
ZColGrps = cell(numRE, 1);
q = zeros(numRE, 1);
for i = 1:numRE
    [Z{i}, ~, ~, ZColGrps{i}, ZColNames{i}] = ...
        classreg.regr.modelutils.designmatrix(dataTable, ...
        'Model', F.RELinearFormula{i}.Terms(:, DatasetsVarLocation), ...
        'DummyVarCoding', 'reference', ...
        'CategoricalVars', logical(IsCategorical(DatasetsVarLocation)), ...
        'CategoricalLevels', Range(DatasetsVarLocation));
    q(i) = size(Z{i}, 2);
end

if options.verbose
    model.RandomInfo.Z = Z;
    model.RandomInfo.ZColNames = ZColNames;
    model.RandomInfo.ZColGrps = ZColGrps;
    model.RandomInfo.q = q;
    model.RandomInfo.numRE = numRE;
end
model.Z = sparseZ(Z, q, lev, Gid);
model.dimRE = lev;

% Additional variable information
model.VarInfo.nVars = nVars;
model.VarInfo.VarNames = VarNames;
model.VarInfo.IsCategorical = IsCategorical;
model.VarInfo.DatasetsVarLocation = DatasetsVarLocation;
model.VarInfo.ResponseName = F.ResponseName;
model.VarInfo.PredictorNames = F.FELinearFormula.PredictorNames;
model.VarInfo.XColNames = XColNames;
model.VarInfo.XCols2Terms = XCols2Terms;
model.VarInfo.Terms = F.FELinearFormula.Terms(:, DatasetsVarLocation);
model.VarInfo.TermNames = F.FELinearFormula.TermNames;

end

%% Helper Function: intvars
function [G, GName] = intvars(ds, intVars)
assert(iscellstr(intVars));
assert(all(ismember(intVars, ds.Properties.VariableNames)));

k = length(intVars);
if k >= 1
    G = categorical(ds.(intVars{1}));
    GName = intVars{1};
    for i = 2:k
        G = G .* categorical(ds.(intVars{i}));
        GName = [GName, ':', intVars{i}];
    end
else
    G = [];
    GName = [];
end

G = removecats(G);
end

%% Helper Function: sparseZ
function Zs = sparseZ(Z, q, lev, Gid)
qlev = q .* lev;
R = length(Gid);
N = size(Z{1}, 1);
Zs = sparse(N, sum(qlev));
for r = 1:R
    for k = 1:lev(r)
        idx = (Gid{r} == k);
        offset = sum(qlev(1:(r-1))) + (k-1)*q(r);
        Zs(idx, offset + 1 : offset + q(r)) = Z{r}(idx, :);
    end
end
end
