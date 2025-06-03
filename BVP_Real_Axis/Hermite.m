function [p]=Hermite(m,t)
%--------------------------------------------------------------------------
% File: Hermite.m
%
% Goal: Construct the orthonormal Hermite polynomials from degree 0 to 
%       degree m+1 evaluated at the points t  
%
% Use: [p] = Hermite(m,x)
%
% Input:  m - degree of the polynomial
%         t - row array of evaluation points
%
% Output: p - matrix of the orthonormal Hermite polynomials  
%             size(H) = length(t)x(m+1)
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
p = zeros(length(t),m+1);

s4Pi = (pi)^(1/4);
s2 = sqrt(2);

p(:,1) = 1/s4Pi;
p(:,2) = (s2*t)/s4Pi;
p(:,3) = (-2+4*t.^2)/(2*s2*s4Pi);

for j=2:m-1
    p(:,j+2) = (t'*sqrt(2/(j+1))).*p(:,j+1)-sqrt(j/(j+1))*p(:,j);
end
end