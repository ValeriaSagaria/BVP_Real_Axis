BVP_Real_Axis: Boundary Value Problem on Real Axis Package 

(c) Developed by Maria Carmela De Bonis, Valeria Sagaria

This package contains code implementing the Fredholm Integral Equation (FIE) arising from second-order Boundary Value Problem (BVP) on the real axis. 
In this toolbox we propose a main algorithm, for approximating the solutions of a BVP. Moreover, we create a Graphical User Interface (GUI) MatLab BVP_Real_Axis which allows the user to plot the solution of a BVP as well as the approximate solution at a given point.

The software implements a numerical method developed in the following two very recent papers:

[1] M.C. De Bonis, V. Sagaria, Numerical method for boundary value problems on the real line. Appl. Numer. Math., vol. 200, (2024), 179–194
[2] M.C. De Bonis, V. Sagaria, A new Nyström method for solving boundary value problems on the real axis, submitted 2025

The algorithms implemented here are described in detail in: 

[3] M.C. De Bonis, V. Sagaria, A MatLab package for solving Fredholm Integral Equations arising from second-order Boundary Value Problems on the real axis, submitted 2025


GETTING STARTED
Run

>> BVP_demo1   (Figure (1) (right) from [2])

>> BVP_demo1_1 (Table 1 and Figure 1 (left) from [2])

>> BVP_demo2   (Figure 3 (right) from [2])

>> BVP_demo2_1 (Table  5 and Figure 3 (left) from [2])

>> BVP_demo3   (Graph of the approximate solution)

>> BVP_demo3_1 (Plot of the trend of the absolute error err_m)

>> BVP_demo4   (Graph of the approximate solution)

>> BVP_demo4_1 (Plot of the trend of the absolute error err_m)


as a demonstration of what this package is used for.

-The file gaussq.m is the result of a translation of the fortran procedure http://www.netlib.org/go/gaussq.f

License : 

The BVP_Real_Axis Package is a Matlab library released under the GPL.

The BVP_Real_Axis Package is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The BVP_Real_Axis Package is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the BVP_Real_Axis Package.  
If not, see <http://www.gnu.org/licenses/>.
