%--------------------------------------------------------------------------
% BVP_demo2 _ 1:  
%
% This demo computes the solution of a Boundary Value Problem
%
% This reproduces Table 4 from [2] 
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

fprintf('Welcome to BVP demo #2 _ 1\n \n');
addpath('..');

% h(t), kind and A(x)=(mu a(x)-1)*exp(x^2/2)
h = @(x) (abs (x).^(13/2).*exp(-3*x.^2/2));
kind = 2;
A = @(x) (1./(x.^8+5).*exp(x.^2/2)); 

% row array of the evaluation points
t=linspace(-6,6);
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
fprintf('Table:\n \n');
fprintf('m & j & cond(A_m) & err_m \n');
for k=1:l-1
    fprintf('% d \t & % d \t & % .4e \t & % .2e \t ',M(k,1),M(k,2),M(k,3),Merr(k,2))
    fprintf('\n')
end
fprintf('\n')
fprintf('Trend of absolute errors \n \n');
figure(2)
x=M(1:l-1,2);
y=Merr(:,2);
loglog(x,y,'o-k','linewidth',1.5);
set(gca,'fontsize',14)
title('Trend of absolute errors')
xlim([x(1),x(end)])
xlabel('2 j')
ylabel('e_m (f)')
xticks(x)
fprintf('Execution time: % .4f s', exTime);
fprintf('\n \n');


