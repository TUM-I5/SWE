/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader (bader AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Univ.-Prof._Dr._Michael_Bader)
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * Main section of the doxygen documentation.
 */

/** \mainpage SWE - A Simple Shallow Water Code

SWE is a teaching code that implements simple Finite Volumes models that solve 
the shallow water equations - in a problem setting as it would be used for 
tsunami simulation.  

\section intro The Shallow Water Equations

The shallow water equations describe the behaviour of a fluid, in particular water, 
of a certain (possibly varying) depth <i>h</i> 
in a two-dimensional domain -- imagine, for example, a puddle of water or a shallow pond
(and compare the 1D sketch given below).
The main modelling assumption is that we can neglect effects of flow in vertical 
direction.
The resulting model proved to be useful for the simulation of tsunami propagation
(with appropriate extensions). While an ocean can hardly be considered as 
"shallow" in the usual sense, tsunami waves (in contrast to regular waves induced 
by wind, e.g.) affect the entire water column, such that effects of vertical flow can 
again be neglected. To allow for a non-even sea bottom (as required for accurate 
modelling of tsunamis), we include the elevation <i>b</i> of the sea floor in our model:

\image latex basin_bathy.pdf
\image html basin_bathy.gif

The shallow water equations describe the changes of water depth <i>h</i> and 
horizontal velocities <i>v<sub>x</sub></i> and <i>v<sub>y</sub></i> (in the resp.\ coordinate directions) over time, 
depending on some initial conditions -- in the case of tsunami simulation, these initial 
conditions could, for example, result from an initial elevation of the sea floor caused 
by an earthquake. The respective changes in time can be described via a system of 
partial differential equations:
\f[
\begin{array}{c} 
    \displaystyle
    \frac{\partial h}{\partial t} + \frac{\partial (v_x h)}{\partial x} + \frac{\partial (v_y h)}{\partial y} = 0 \\[1em]
    \displaystyle
    \frac{\partial (h v_x)}{\partial t} + \frac{\partial (h v_x v_x)}{\partial x}  + \frac{\partial (h v_y v_x)}{\partial y}
    + \frac{1}{2} g \frac{\partial (h^2)}{\partial x}  = - gh \frac{\partial b}{\partial x},
\\[1em]
    \displaystyle
    \frac{\partial (h v_y)}{\partial t} + \frac{\partial (h v_x v_y)}{\partial x}  + \frac{\partial (h v_y v_y)}{\partial y}
    + \frac{1}{2} g  \frac{\partial (h^2)}{\partial y} = - gh \frac{\partial b}{\partial y},
\end{array}
\f]
The equation for <i>h</i> is obtained, if we examine the conservation of mass in a control 
volume. The equations for <i>hv<sub>x</sub></i> and <i>hv<sub>y</sub></i> result from conservation of momentum 
(note that <i>h</i> is directly related to the volume, and thus the mass of the water 
-- thus <i>hv<sub>x</sub></i> can be interpreted as a momentum). 

The two terms involving <i>g</i> model a gravity-induced force 
(<i>g</i> being the constant for the gravitational acceleration, <i>g</i> = 9.81 ms<sup>-2</sup>),
which results from the hydrostatic pressure.
The right-hand-side source terms model the effect of an uneven ocean floor
(<i>b</i> obtained from respective bathymetry data).

\subsection finvol Finite Volume Discretisation

The shallow water equations are usually too difficult to be solved exactly -
hence, SWE implements simple discrete models as an approximation. 
As the applied numerical method (typically a Finite Volume discretization) may vary, 
we will stick to the basics at this point.

First, SWE assumes that the unknown functions <i>h(t,x,y)</i>, 
<i>hu(t,x,y) := h(t,x,y) v<sub>x</sub>(t,x,y)</i>, 
<i>hv(t,x,y) := h(t,x,y) v<sub>y</sub>(t,x,y)</i>, 
as well as the given sea bottom level <i>b(x,y)</i>, 
are approximated on a Cartesian mesh of grid cells, as illustrated below. 
In each grid cell, with indices <i>(i,j)</i>, 
the unknowns have constant values 
<i>h<sub>ij</sub></i>, <i>hu<sub>ij</sub></i>, <i>hv<sub>ij</sub></i>, and <i>b<sub>ij</sub></i>:

\image html grid_unknowns.gif
\image latex grid_unknowns.pdf

\subsection fluxes Computing Numerical Fluxes at Edges and Euler Time-Stepping

The details of the numerical schemes are too complicated to be described in this 
overview. Please refer to the accompanying material.
To put it short, we successively perform two main computational steps: 
- we compute so-called <b>numerical fluxes</b> on each edge of the grid (which 
approximate the transfer of mass or momentum between grid cells), 
- based on these numerical fluxes, we then update the unknowns in each cell. 


\section impl Implementation and base class SWE_Block

For the simulation of the shallow water model, we thus require a regular Cartesian 
grid, where each grid cell carries the respective unknowns - water level, momentum 
in x- and y-direction, and bathymetry data. 
The central data structures for Cartesian grid and arrays of unknowns are 
provided with the abstract base class SWE_Block, which has four 2D arrays
SWE_Block::h, SWE_Block::hu, SWE_Block::hv, and SWE_Block::b.
To implement the behaviour of the fluid at boundaries, and also to allow the 
connection of several grid blocks (for parallelization or just to build more 
complicated computational domains),  
each array has an additional layer of so-called <i>ghost cells</i>, 
as illustrated in the following figure:

\image latex ghost_cells.pdf
\image html ghost_cells.gif

\subsection models Parallelisation and Different Models

In each time step, our numerical algorithm will compute the flux terms for each 
edge of the computational domain. To compute the fluxes, we require the values 
of the unknowns in both adjacent cells. At the boundaries of the fluid domain, 
the ghost layer makes sure that we also have two adjacent cells for the 
cell edges on the domain boundary. The values in the ghost layer cells will be 
set to values depending on the values in the adjacent fluid domain. 
We will model three different situations:
\begin{description}
\item{Outflow:} \texttt{h}, \texttt{u}, and \texttt{v} in the ghost cell are
   set to the same value as in the adjacent fluid cell. This models the situation 
   that the unknowns do not change across the domain boundary
   (undisturbed outflow).
\item{Wall:} At a wall, the velocity component normal to the boundary should 
   be $0$, such that no fluid can cross the boundary. To model this case, we 
   set the normal velocity, e.g.\ \texttt{u[0]} at the left boundary, to 
   the negative value of the adjacent cell: \texttt{-u[1]}. 
   The interpolated value at the boundary edge will then be $0$ 
   (\texttt{h} is identical in both cells due to the imposed boundary condition).
   The other two variables are set in the same way as for the outflow situation.  
\item{Connect:} With the connect case, we can connect a domain at two boundaries.
   If we connect the left and right boundary, we will obtain a periodically 
   repeated domain. 
   Here, all ghost values are determined by the values of the 
   unknowns in the fluid cell adjacent to the connected boundary. 
\end{description}

To implement the boundary conditions, the class \texttt{SWE\_Block} contains an 
array of four enum variables, \texttt{boundary[4]} (for left/right/bottom/top 
boundary), that can take the values 
<code>OUTFLOW</code>, <code>WALL</code>, and <code>CONNECT</code>.


\subsection multblocks Multiple Blocks

Via the connect boundary condition, it is also possible to connect several 
Cartesian grid blocks to build a more complicated domain. Figure \ref{fig:connect} 
illustrates the exchange of ghost values for two connected blocks.

\image html connect.gif
\latexonly
\begin{figure}
\begin{center}
  \includegraphics[width=0.8\textwidth]{pics/connect.pdf}
\end{center}
\caption{Exchange of values in ghost layers between two connected \texttt{SWE\_Block}s.}
\label{fig:connect}
\end{figure}
\endlatexonly

To store the neighbour block in case of a <code>CONNECT</code> boundary, <code>SWE_Block</code> 
contains a further array of four pointers, <code>neighbour[4]</code> (for left/right/bottom/top 
boundary), that will store a pointer to the connected adjacent <code>SWE_Block</code>.

The respective block approach can also be exploited for parallelisation: 
the different blocks would then be assigned to the available processors 
(or processor cores) -- each processor (core) works on its share of blocks, 
while the program has to make sure to keep the values in the ghost cells 
up to date (which requires explicit communication in the case of distributed-memory 
computers).

\subsection*{The Time Step Loop}

For each time step, our solver thus performs the following steps -- each step 
is implemented via a separate member function of the class \texttt{SWE\_Block}: 
\begin{enumerate}
\item set the values at the boundaries: \texttt{setBoundaryLayer()};
\item compute the flux terms for all edges: \texttt{computeFluxes()};
\item from the flux terms, compute the in/outflow balance for each cell, 
     and compute the new values of the unknowns for the next time step:
     \texttt{eulerTimestep()}.
\end{enumerate}

\subsection*{Visualisation with ParaView}

The class \texttt{SWE\_Block} contains further methods that will write the
numerical solution into a sequence of files that can be read by the visualisation package 
ParaView (just enter the respective folder from ParaView -- the files will be recognised 
and displayed as one project). 
ParaView allows to visualise the computed time-dependent solution (as ``movie'' or 
in single-step mode). ParaView is pretty self-explanatory for our purposes, but provides 
an online help for further instructions.

\subsection*{Visualisation ``on the fly''}

We also provide a CUDA implementation of the simulation code (requires a computer 
with a CUDA-capable GPU, together with the respective drivers -- visit NVIDIA's 
website on CUDA for details on implementation). 
Apart from the fact that the simulation usually runs a lot faster on the GPU, 
the program is also capable of plotting the computing solution (water surface) 
``on the fly''.

Finally: whoever thinks that they can do a better (faster, ...)
implementation (visualisation, ...) of the provided code is more than welcome to do so!
Feel free to contribute to SWE - for questions, just contact Michael Bader (bader@in.tum.de).

*/
