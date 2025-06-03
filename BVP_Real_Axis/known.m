function [j,x,w,M,Ax,y,wy]=known(m,A,kind)
%--------------------------------------------------------------------------
% File: known.m
%
% Goal: identify the truncation parameter j, evaluate the known function 
%       A(x) and the matrix of the modified moments M_k(x) at the Hermite 
%       zeros
%
% Use: [j,x,w,M,Ax,y,wy]=known(m,A,kind)
%
% Input:  m - number of quadrature knots
%         A - A(x)=(mu a(x)-1)*exp(x^2/2)
%         kind - |1 if the function H(t) has an exact expression
%                |2 if the function H(t) is approximated by Gauss-Laguerre
%                                        quadrature rule
%
% Output: j  - truncation parameter 
%         x  - row array of the Hermite zeros x_k, k=-j,...,j 
%         w  - row array of the Hermite weights w_k, k=-j,...,j
%         M  - matrix of the modifeid moments 
%         Ax - A(x_k)=(mu a(x_k)-1)*exp((x_k)^2/2), k=-j,...,j
%
% Recalls: gaussq.m, Mk.m
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
if kind == 2
    [y,wy]=gaussq(6,1024,0,0,0,[0,0]);
else
    y = 0;
    wy = 0;
end
[x,w]=gaussq(4,m,0,0,0,[0,0]);

Ax = A(x);

M = Mk(m,x,w,x);
j = 1;
ind = floor(m/2);
if abs(Ax(ind)*M(ind,j))<0.5*10^(-16)
    while (j<=m && abs(Ax(ind)*M(ind,j))<0.5*10^(-16))
        j = j+1;
    end
    j = j-1;
    if j~=m
        x = x(j:m-j+1);
        w = w(j:m-j+1);  
        Ax = Ax(j:m-j+1);
        M = M(j:m-j+1,j:m-j+1);
    end
end
j = length(x);
end