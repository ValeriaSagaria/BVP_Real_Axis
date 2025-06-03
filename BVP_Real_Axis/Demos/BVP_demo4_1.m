%--------------------------------------------------------------------------
% BVP_demo4_1:  
%
% This demo computes the solution of a Boundary Value Problem
%
% This reproduces figure ??? from "???", MC De Bonis, V Sagaria, 2025. 
%
% Author:  MC De Bonis, V Sagaria
%
% Date last modified: May, 2025
%
% This returns return both a table containing in each row the values 
% of m, j,cond(A_m) and err_m, for increasing values of m,  and the plot 
% of the trend of the absolute error err_m.
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

fprintf('Welcome to BVP demo #4_1\n\n');

% h(t), kind and AA(x)=(mu a(x)-1)*exp(x^2/2)
h = @(x) (-0.5*exp(-17*x.^2/30)-exp(-x.^2/2)+x.^2.*exp(-x.^2/2));
kind = 2;
A = @(x) (exp(x.^2/2).*(0.5*exp(-x.^2/15)-1));

% expression of exact solution
fe=@(x) (exp(-x.^2/2));

% row array of the evaluation points
t = linspace(-3,3,30000);
time=tic;
l = 9;
M=zeros(l,length(t)+3);
for i=2:l
    m=2^i;
    [fm,C,j] = NystromMethodFIE(m,kind,h,A,t);
    M(i-1,:)=[m j C fm];
end

% evaluation of exact solution
[fm,C,j] = NystromMethodFIE(600,kind,h,A,t);
M(l,:)=[m j C fm]; 

% compute the absolute errors
Merr=zeros(l-1,2);
err=zeros(length(t),l-1);
for k=2:l
    for i=1:length(t)
        err(i,k-1)=abs(M(l,i+3)-M(k-1,i+3));
    end
    Merr(k-1,:)=[2^k max(err(:,k-1))];
end
exTime = toc(time);
fprintf('Table\n\n');
fprintf('m & j & cond(A_m) & err_m \n');
for k=1:l-1
    fprintf('%d \t & %d \t & %.4e \t & %.2e \t ',M(k,1),M(k,2),M(k,3),Merr(k,2))
    fprintf('\n')
end
fprintf('\n')
fprintf('Trend of absolute errors\n\n');
figure(2)
x=M(1:l-1,2);
y=Merr(:,2);
loglog(x,y,'o-k','linewidth',1.5);
set(gca,'fontsize',14)
title('Trend of absolute errors')
xlim([x(1),x(end)])
xlabel('2 j')
ylabel('e_{m}(f)')
xticks(x)
fprintf('Execution time: %.4f s', exTime);
fprintf('\n\n');