%--------------------------------------------------------------------------
% BVP_demo3:  
%
% This demo computes the solution of a Boundary Value Problem
%
% This returns the graph of the approximate solution together with the 
% corresponding truncation parameter j and the condition number of the 
% involved linear system and the maximum absolute error. 
%
% Author:  MC De Bonis, V Sagaria
%
% Date last modified: June, 2025
%
% This file is part of the BVP_Real_Axis package Copyright (C) 2025, 
% MC De Bonis, V Sagaria.
%
% The BVP_Real_Axis package is free software: you can redistribute it 
% and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation.
%
% The BVP_Real_Axis package is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with theBVP_Real_Axis package. 
% If not, see <http://www.gnu.org/licenses/>.
% %%
clc
clear
close all

fprintf('Welcome to BVP demo #3\n\n');

% (h(x)+mu*b*a(x)), kind, A(x)=(mu a(x)-1)*exp(x^2/2) and b
h = @(x) (exp(-11*x.^2/(10)).*(-2+4*exp(x.^2/(10)).*(-1+x.^2)));
kind = 2;
A=@(x) (2*exp((2*x.^2)/5)+exp(x.^2/2));
b=4;

% row array of the evaluation points
t=linspace(-2,2,30000);
time=tic;
[fm,C,j,Merr,~] = BVP(kind,h,A,t);
exTime = toc(time);
fprintf('m, truncation index j, condition number and maximum error:\n\n');
fprintf('%d,\t % d,\t  %.4e, \t  %.2e. ',512,j,C,Merr);
fprintf('\n\n');

fprintf('Plot of the solution of BVP #3\n\n');
figure
plot(t,fm,'-k','linewidth',1.5)
set(gca,'fontsize',14)
title('Plot of the solution of BVP #3')
xlabel('t')
xlim([t(1),t(end)])
ylabel('f_{512}(t)')
fprintf('Execution time for evaluating f_{512}(t): %.4f s', exTime);
fprintf('\n\n');