function [fm] = NystromInterpolant(m,x,w,z,Ax,kind,h,y,wy,t)
%--------------------------------------------------------------------------
% File: NystromInterpolant.m
%
% Goal: Compute the Nystrom interpolant
%
% Use: [fm] = NystromInterpolant(m,x,w,z,Ax,kind,h,y,wy,t)
%
% Input:  m - number of quadrature knots
%         x - row array of the Hermite zeros x_k, k=-j,...,j
%         w - row array of the weights w_k, k=-j,...,j m-point 
%             Gauss-Hermite rule
%         Ax - A(x_k)=(mu a(x_k)-1)*exp((x_k)^2/2), k=-j,...,j
%         z  - solution of the linear system
%         kind - |1 if the function H(t) of FIE has an exact expression
%                |2 if the function H(t) of FIE is approximated by 
%                                        Gauss-Laguerre quadrature rule
%         h    - |H(t) right-hand side of the FIE if kind == 1
%                |h(x) right-hand side of the BVP if kind == 2
%         y  - row array of the Laguerre zeros y_k, k=1,...,N, N=1024
%         wy - row array of the weights of the N-point Gauss-Laguerre rule
%              with N=1024
%         t  - row array of evaluation points
%
% Output: f - Nystrom interpolant, solution of the BVP
%
% Recalls: Mk.m, Hm.m
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
M = Mk(m,x,w,t);
if kind == 1 
    fm = (z'.*Ax)*M'+ h(t); 
elseif kind == 2 
    fm = (z'.*Ax)*M'+Hm(1024,h,y,wy,t)';  
else 
    disp('WARNING: kind must be 1 or 2')
end
end