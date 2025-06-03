function [H] = Hm(m,h,y,wy,t)
%--------------------------------------------------------------------------
% File: Hm.m
%
% Goal: Approximate the integral
%              +oo
%       H(t) = | G(x,t)h(x) dx, 
%            -oo 
%       using the N-point truncated Gauss-Laguerre quadrature rule
%
% Use: [H] = Hm(m,h,y,wy,t)
%
% Input:  m  - m-point Gauss-Laguerre rule
%         h  - h(x) right-hand side of the BVP
%         y  - row array of the Laguerre zeros y_k, k=1,...,N, N=1024
%         wy - row array of the weights of the N-point Gauss-Laguerre rule
%              with N=1024
%         t  - row array of evaluation points
%
% Output: Hm -  column array of the approximations of H(t) 
%               length(Hm) = length(t)
%
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
H = zeros(s,1);
hMinus = h(t'-y);
hPlus = h(t'+y);
for i = 1:s
    j = 1;         
    while (j<m && abs((hMinus(i,j)+hPlus(i,j))*wy(j))>0.5*10^(-20))
       j = j+1;
    end
    j = j-1;

    H(i) = (hMinus(i,1:j)+hPlus(i,1:j))*wy(1:j)';
    H(i)=-0.5*H(i);
end
end