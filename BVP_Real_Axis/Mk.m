function [M] = Mk(m,x,w,t)
%--------------------------------------------------------------------------
% File: M.m
%
% Goal: Build the matrix of the modified moment
%                    +oo
%           M_k(t) = | G(x,t)lmk(x) exp(-x^2/2)dx
%                   -oo 
%
% Use: [M] = Mk(m,x,w,t)
%
% Input:  m - size 
%         x - row array of the Hermite zeros x_k, k=-j,...,j
%         w - row array of the weights w_k, k=-j,...,j m-point 
%             Gauss-Hermite rule
%         t - row array of evaluation points
%
% Output: M - matrix of the modifeid moments evaluated at points t
%             size(M) = length(t)xm 
% 
% Recalls: Ci.m, Hermite.m
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
c = Ci(m-1,t);    
p = Hermite(m-1,x);

M = ((p.*w')*c)';
end