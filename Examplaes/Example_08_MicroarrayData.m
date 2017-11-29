%% EAXAMPLE: Micro Array data
%  Data analysis by Linear Mixed Model with MATLAB function HPMIXED

% Author: Viktor Witkovsky (witkovsky@savba.sk)
% Ver.: 06-Jan-2014 13:52:11

%% Data Description
% SAS 9.3 Example 45.4 Mixed Model Analysis of Microarray Data 
%
% Microarray experiments are an advanced genomic technique used in the
% discovery of new treatments for diseases. Microarray analysis allows for
% the detection of tens of thousands of genes in a single DNA sample. 
% A microarray is a glass slide or membrane that has been spotted or
% "arrayed" with DNA fragments or oligonucleotides representing specific
% genes. The response of the gene detected by a spot is proportional to the
% intensity of fluorescence associated with that spot. These gene responses
% can indicate associations with disease conditions, but they can also be
% affected by systematic biases and different treatments such as sex and
% genotypes. Statistical models for microarray data attempt to assess the
% significance and magnitude of gene effects across treatments while
% adjusting for these systematic biases and to evaluate the significance of
% differences between treatments. There are two statistical approaches
% frequently used in mixed model analysis for microarray data. The first
% approach is to fit multiple gene-specific models to data normalized for
% systematic biases (Wolfinger et al.; 2001; Gibson and Wolfinger; 2004).
% This approach is based on assuming that the biases are independent from
% the gene effects. If this assumption is untenable, then a second approach
% fits a single model that combines both the systematic biases and the gene
% effects (Kerr, Martin, and Churchill; 2000; Churchill; 2002; Littell et
% al.; 2006). When the number of genes is very large, several hundreds to
% tens of thousands, this is an analysis for which the sparse matrix
% approach implemented in the HPMIXED procedure is well suited. 
%
% The SAS statements simulate a microarray experiment with a so-called loop
% design structure, which is commonly used in such studies. There are 500
% genes, each gene occurs in 6 arrays, and each array has 2 dyes.
% A linear mixed model for fitting the log intensity data  from such a
% design is described by Littell et al. (2006).
%
% You can use the SAS Proc HPMIXED procedure with the following statements
% to fit this model:
% proc hpmixed data=microarray;
%    class marray dye trt gene pin dip;
%    model log2i = dye trt gene dye*gene trt*gene pin;
%    random marray marray*gene dip(marray) pin*marray;
%    test trt;
% run;

%% Load dataset dsMicroarrayData
clear
load dsMicroarrayData

%% Create the model structure for HPMIXED by using hpmixedmodel
% If MATLAB 2013b + Statistics Toolbox is available
% Otherwise use provided model structure, or construct it manually

%formula  = 'log2i ~ Trt + (1 | Gene) + (1 | MArray:Dye:Pin:Dip)';
%formula  = 'log2i ~ Dye + Trt + Gene + (1 | MArray) + (1 | MArray:Gene)';
formula  = 'log2i ~   Trt:Gene  +  (1 | MArray)';

tic;
model = hpmixedmodel(MicroArrayData,formula);
toc
model.Description = 'MicroArrayData: SAS PROC HPMIXED Example';

%% Fit the linear mixed model by HPMIXED with limitted  output
opts.verbose = false;
lmefit1 = hpmixed(model,opts);

disp(lmefit1)

%% Fit the linear mixed model by HPMIXED with complete output
opts.verbose = true;
lmefit = hpmixed(model,opts);

disp(lmefit)

%% Fit the linear mixed model by FITLME (! Long execution time !)
% tic;
% lme = fitlme(MicroArrayData,formula,'FitMethod','REML');
% toc
% 
% disp(lme)
