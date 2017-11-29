function [A,B,C] = Design2(N)
%DESIGN2 Creates the design matrices for two-way classification
%        model y_ijk = a_i + b_j + c_ij + e_ijk,
%        given by its incidence matrix N (ixj).
%        The incidence matrix N could have some cells equal to zero.
%
%======================================================================
% Syntax:
%        [A,B,C]=design2(N);
%======================================================================
%        See also:  diagm.m

%======================================================================
% Ver.: 2.0
% Revised 19-Dec-2001 20:31:48
% Copyright (c) 1998-2001 Viktor Witkovsky
%======================================================================
%	BEGIN DESIGN2.M
%======================================================================
[ni,nj] = size(N);
M  = N';
v  = M(:);
v  = v(v>0);
oj = ones(nj,1);
n  = N*oj;
%
A = Diagm(n);
B = Diagm(N(1,:));
for i = 2:ni
    B = [B;Diagm(N(i,:))];
end
C = Diagm(v);
%======================================================================
%	EOF DESIGN2.M
%======================================================================