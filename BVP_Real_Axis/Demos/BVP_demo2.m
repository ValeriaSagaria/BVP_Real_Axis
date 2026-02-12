%--------------------------------------------------------------------------
% BVP_demo2:  
%
% This demo computes the solution of a Boundary Value Problem
%
% This reproduces Figure 3 from [2] 
%
% Author:  MC De Bonis, V Sagaria
%
% Date last modified: February, 2026
%
% This file is part of the BVP_Real _Axis package Copyright (C) 2025, 
% MC De Bonis, V Sagaria.
%
% The BVP_Real _Axis package is free software: you can redistribute it 
% and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation.
%
% The BVP_Real _Axis package is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with theBVP_Real _Axis package. 
% If not, see <http://www.gnu.org/licenses/>.
clc
clear
close all

fprintf('Welcome to BVP demo #2 \n \n');
addpath('..');

% h(t), kind and A(x)=(mu a(x)-1)*exp(x^2/2)
h = @(x) (abs (x).^(13/2).*exp(-3*x.^2/2));
kind = 2;
A = @(x) (1./(x.^8+5).*exp(x.^2/2));

% row array of the evaluation points
t=linspace(-6,6);
time=tic;
[fm,C,j,Merr,~] = BVP(kind,h,A,t);
exTime = toc(time);
fprintf('m, truncation index j, condition number and maximum error:\n \n');
fprintf('% d,\t % d,\t  % .4e, \t  % .2e. ',512,j,C,Merr);
fprintf('\n \n');

fprintf('Plot of the solution of BVP #3 \n \n');
figure
plot(t,fm,'-k','linewidth',1.5)
set(gca,'fontsize',14)
title('Plot of the solution of BVP #2')
xlabel('t')
xlim([t(1),t(end)])
ylabel('f_{512}(t)')
fprintf('Execution time for evaluating f_{512}(t): % .4f s', exTime);
fprintf('\n \n');


