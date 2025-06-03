function [fm,C,j,Merr,err] = BVP(kind,h,A,t)
%--------------------------------------------------------------------------
% File: BVP.m
%
% Goal: Compute the solution of the BVP (Boundary Value Problem)
%       
%           |f''(x)-mu a(x)f(x) = h(x)
%           |f(-oo) = f(+oo) =0
%
% Use: [fm,C,j,Merr,err] = BVP(kind,h,A,t)
%
% Input:  kind - |1 if the function H(t) of FIE has an exact expression
%                |2 if the function H(t) of FIE is approximated by 
%                                        Gauss-Laguerre quadrature rule
%         h    - |H(t) right-hand side of the FIE if kind == 1
%                |h(x) right-hand side of the BVP if kind == 2
%         A    - A(x)=(mu a(x)-1)*exp(x^2/2)
%         t    - row array of evaluation points
%
% Output: fm   - solution of the BVP 
%         j    - size of the solved linear system
%         C    - condition number of the solved linear system
%         Merr - maximum of absolute errors
%         err  - row of absolute errors
%         
%
% Recalls: NystromMethodFIE.m
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
[f0,~,~] = NystromMethodFIE(600,kind,h,A,t);
if isnan(f0)
    error('The method is not applicable')
end
[fm,C,j] = NystromMethodFIE(512,kind,h,A,t);
% control of the condition number of the solved linear system 
if C > 1/eps
    sprintf('Warning: Results may be inaccurate. COND = %d',C)
end
% compute the row of absolute errors
err = abs(f0-fm);
% compute the maximum of absolute errors
Merr = max(err);
% control of the convergence of the method using the weigthed absolute error
if Merr > 1
    disp('Error: The method is not convergent');
    fm = [];
end
end