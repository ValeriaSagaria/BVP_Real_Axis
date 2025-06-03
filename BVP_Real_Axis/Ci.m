function [c]=Ci(m,t)
%--------------------------------------------------------------------------
% File: Ci.m
%
% Goal: Evaluate the integrals 
%             +oo
%    c_i(t) = | G(x,t)p_i(x) exp(-x^2/2)dx,  k=0,...,m
%           -oo 
%       using a recurrence scheme
%
% Use: [c] = Ci(m,t)
%
% Input:  m - number of Hermite zeros
%         t - row array of evaluation points
%
% Output: c - matrix (c_i(t_j)) 
%             size(c) = (m+1)xlength(t)
%
% Recalls: Hermite.m
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
%--------------------------------------------------------------------------
s = length(t);
p = Hermite(m,t);

c = zeros(m+1,s);
c1i = zeros(m+1,s);
c2i = zeros(m+1,s);

s2 = sqrt(2);
erPlus = erfc((1+t)/s2);
erMinus = erf((t-1)/s2);
sPi = sqrt(pi);
ssPi = sqrt(sPi);
K = -sqrt(exp(1)*sPi)/(2*s2);
et2 = exp(-t.^2/2);

c1i(1,:) = K*(1+erMinus).*(cosh(t)-sinh(t));
c2i(1,:) = K*exp(t).*erPlus;
c(1,:) = c1i(1,:)+c2i(1,:);

c1i(2,:) = (et2/(2*ssPi)+c1i(1,:))*s2;
c2i(2,:) = -(et2/(2*ssPi)+c2i(1,:))*s2;
c(2,:) = (c1i(1,:)-c2i(1,:))*s2;

c1i(3,:) = p(:,2)'.*et2/2+c1i(2,:)+c1i(1,:)/s2;
c2i(3,:) = -p(:,2)'.*et2/2-c2i(2,:)+c2i(1,:)/s2;
c(3,:) = c1i(2,:)-c2i(2,:)+c(1,:)/s2;

for j=3:m
    B = s2/sqrt(j);
    A = B*p(:,j)'.*et2/2;
    C = sqrt((j-1)/j); 
    c1i(j+1,:)=A+B*c1i(j,:)+C*c1i(j-1,:);
    c2i(j+1,:)=-A-B*c2i(j,:)+C*c2i(j-1,:);
    c(j+1,:)=B*(c1i(j,:)-c2i(j,:))+C*c(j-1,:);
end
end