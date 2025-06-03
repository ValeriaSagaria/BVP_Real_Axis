function [fm,C,j] = NystromMethodFIE(m,kind,h,A,t)
%--------------------------------------------------------------------------
% File: NystromMethodFIE.m
%
% Goal: Compute the solution of the Fredholm Integral Equation (FIE)
%
%         +oo                         +oo
%   f(t) - |G(x,t)(mu a(x)-1)f(x)dx = |G(x,t)h(x)dx 
%        -oo                         -oo 
%
%
% Use: [fm,C,j,m] = NystromMethodFIE(m,kind,h,A,t)
%
% Input:  m    - number of quadrature knots
%         kind - |1 if the function H(t) of FIE has an exact expression
%                |2 if the function H(t) of FIE is approximated by 
%                                        Gauss-Laguerre quadrature rule
%         h    - |H(t) right-hand side of the FIE if kind == 1
%                |h(x) right-hand side of the BVP if kind == 2
%         A    - A(x)=(mu a(x)-1)*exp(x^2/2)
%         t    - row array of evaluation points
%
% Output: fm   - solution of the BVP
%         C    - condition number of the solved linear system
%         j    - size of the solved linear system
%         
%          
% Recalls: known.m, build.m, NystromInterpolant.m
%
% Author:  MC De Bonis, V Sagaria
%
% Date last modified: May, 2025
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
[j,x,w,M,Ax,y,wy] = known(m,A,kind);

[Am,b] = build(j,x,y,wy,M,Ax,kind,h);

C = cond(Am,inf);
z = Am\b;

fm = NystromInterpolant(m,x,w,z,Ax,kind,h,y,wy,t);

end