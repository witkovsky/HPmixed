function result = getAnova(lmefit,options)
%getAnova Perform hypothesis tests (Type III tests) on fixed effect, based
%         on the lmefit, by using the FAI-CORNELIUS approximate method for
%         denominator degrees of freedom (DDF).
%
% SYNTAX:
% STAT = getAnova(lmefit,options)
%
% INPUTS:
%  lmefit   - lmefit, fitted linear mixed effects model by HPMIXED
%  options. - structure of options with the following default values:
%         options.STAT.verbose = true
%         options.STAT.alpha = 0.05
%         options.STAT.Description = ''
%         options.STAT.colnames = 'fit_'
%         options.STAT.effectID = []
%         options.STAT.ddfMethod = 'Fai-Cornelius' ('ChiSquare')
%
% EXAMPLE 0: (Fit and analysis of the SplitPlotData by using HPMIXED)
%   clear
%   load dsSplitPlotData;
%   formula  = 'y ~ A + B + A:B + (1 | Block) + (1 | Block:A)';
%   model = hpmixedmodel(SplitPlotData,formula);
%   lmefit = hpmixed(model);
%   options.STAT.ddfMethod = 'Fai-Cornelius';
%   anovaF  = getAnova(lmefit,options)
%   options.STAT.ddfMethod = 'ChiSquare';
%   anovaChi = getAnova(lmefit,options)

% (c) Viktor Witkovsky (witkovsky@savba.sk)
% Revision: 14-Aug-2022 16:44:34
% Ver.: 16-Mar-2014 18:00:15

%% CHECK INPUTS / OUTPUTS
narginchk(1,2);
if nargin < 2, options = []; end

%% SET the default options
if ~isfield(options, 'STAT'), options.STAT = []; end
if ~isfield(options.STAT, 'verbose'), options.STAT.verbose = false; end
if ~isfield(options.STAT, 'alpha'), options.STAT.alpha = 0.05; end
if ~isfield(options.STAT, 'Description'), options.STAT.Description = ''; end
if ~isfield(options.STAT, 'colnames'), options.STAT.colnames = 'fit_'; end
if ~isfield(options.STAT, 'effectID'), options.STAT.effectID = []; end
if ~isfield(options.STAT, 'ddfMethod'), options.STAT.ddfMethod = 'Fai-Cornelius'; end

%% Get the required information from the fittetd model (lmefit)
q           = lmefit.ModelInfo.randomEffectsLength;
m           = lmefit.ModelInfo.MMEmatRows;
TermNames   = lmefit.Details.Model.Formula.FELinearFormula.TermNames;
XCols2Terms = lmefit.Details.Model.FixedInfo.XCols2Terms;

%% ALGORITHM
nTerms  = length(TermNames); 
P       = zeros(nTerms,1);
F       = zeros(nTerms,1);
DF1     = zeros(nTerms,1);
DF2     = zeros(nTerms,1);
Lambda  = cell(1,nTerms);

for k = 1:nTerms
    Cols        = find(XCols2Terms == k);
    nCols       = length(Cols);
    Lambda{k}   = sparse(q + Cols,1:nCols,1,m,nCols);
    rowTab      = getStats(Lambda{k},lmefit,options);
    P(k)        = rowTab.TABLE.pValue;
    F(k)        = rowTab.TABLE.Stat;
    DF1(k)      = rowTab.TABLE.NDF;
    DF2(k)      = rowTab.TABLE.DDF;
end

switch lower(options.STAT.ddfMethod)
    case {'faicornelius' 'fc' 'fai' 'fai-cornelius' 'cornelius'}
        result.TABLE = dataset({TermNames,'Term'},{F,'FStat'},...
            {DF1,'DF1'},{DF2,'DF2'},{P,'pValue'});
    case {'chisquare' 'chi' 'chi-square'}
        result.TABLE = dataset({TermNames,'Term'},{F,'ChiSquareStat'},...
            {DF1,'DF1'},{DF2,'DF2'},{P,'pValue'});
end

result.Details.Lambda = Lambda;
result.options = options;

disp(result.TABLE)

end     % END of getAnova