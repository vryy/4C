/*!
\addtogroup Ale
\brief This module 'Ale' contains all routines and structures necessary
for the 2D and 3D Ale-elements.

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

This module 'Ale' contains all routines and structures necessary
for the 2D and 3D Ale-elements. It provides routines to control the
'academic' pure ale problemtype and to calculate the element stiffness
matrix and the load vectors for both the 2D and 3D elements.
For 2D 4 noded quadrilateral elements are implemented as well as linear
triangles, for 3D a 8 noded brick element is used.

The ALE has been enlarged to treat quadratic 2D elements (quad8 and quad9) as
well. Nevertheless it is somewhat cumbersome to accurately obtain the minimal
Jacobian of a higher order element. So an estimation based on the assumption of
stright sides with central placed mid nodes is used.
It is strongly recommended to ensure that there is no jump in Dirichlet boundary
conditions within the nodes of one side of a higher order element. ie. the nodes
sharing one gline should all take part at the same dbc and include also the end
nodes of a line.

For two-dimensional problems a variety of different elements is available
which include stiffened and nonstiffened pseudo structure elements as well
as elements based on springs and laplacian smoothing. The different ale
elements and the corresponding calculating procedures can be chosen by setting
the appropriate ALE_TYPE within the ALE DYNAMIC block of the input file.
The implemented ALE_TYPEs are the following:

- classic_lin

   Leads to a purely linear treatment of the ale elements. The stiffness
   matrices are calculated and assembled once. The solution for different
   load steps is calculated using different right hand side (rhs) vecotors.
   Within the stiffness calculation it is possible to use the standard 4 (8)
   point gaussian quadrature or a one point quadrature with hourglass
   stabilization.
   It is also possible to disregard the Jacobian determinant integrating
   the element stiffness.

- min_Je_stiff

   Is a pseudo structure treatment of the ale problem as well. The element
   stiffnesses are calculated newly for each time or load step. The actual
   element geometry is taken into account. Further more the element
   stiffness is increased by the inverse of the square of the minimal
   Jacobian determinant found within the element. (For linear quadrilaterals
   the minimal Jacobian determinant can be found at one of the nodes, within
   linear triangles the Jacobian is constant.)
   The ALE_TYPE 'min_Je_stiff' can be used in combination with some initial
   steps set in NUM_INITSTEP. This causes a number of initial steps which are
   not stiffened (but refering to the actual geometry) and which include an
   additional element right hand side that consists of exactly those nodal
   forces which are neccessary to give the unconstraint element a quadratic
   shape after the next step. Initsteps can be useful for unstructured meshes
   of quadrilaterals and can improve the mesh quality.

- two_step

   This calculation follows a paper by Chiandussi et al. 'A simple method
   for automatic update of finite element meshes' in Commun. Numer. Meth.
   Engng. 2000; 16: 1-19.
   The pseudo structure calculation is performed in two steps. A first step
   works with a uniform youngs modulus and evaluates the strain distribution
   which in a second step is used to obtain a non uniform stiffness
   distribution. The elemental stiffness factor within the second step is
   calculated using the square norm of element principal strain criterion, which
   is eq. (13) in the mentioned paper.
   The two_step calculation has been extended to treat quadratic elements as
   well.

- springs

   The ale calculation with springs follows a paper by Farhat et al.:
   'Torsional springs for two-dimensional dynamic unstructured fluid
   meshes'. Comput. Methods Appl. Mech. Engrg. 163 (1998) 231-245.
   The stiffness is obtained from springs connecting every pair of nodes
   within an element as well as torsional springs at the nodes.
   It works for linear triangles and quadrilaterals.
   The spring treatment of quadratic elements is somewhat special. The elemental
   corner nodes are connected by springs just like the linear case. The
   stiffnesses for the edge and mid nodes are determined such that these nodes
   remain on the middle of the respective geometric element.
   Note that this causes non-symmetric stiffness matrices and calls for adequate
   solvers!

- laplace

   Basing on a paper by Loehner et al. 'Improved ale mesh velocities
   for moving bodies' Commun. in Numer. Methods in Engng. 12: 599-608,
   1996 laplacian smoothing is implemented as well. Rather than the
   elasticity equations it is solved

   div( k grad v) = 0,

   for the velocities (or displacement increments) v, where k is a scalar diffusion
   coefficient. In difference to the mentioned
   paper k is not determined from the distance to the moving body (circumventing
   the need for body search algorithms) but rather
   k = 1/(min_det_J)^2, which increases the diffusion coefficient in regions
   where large displacement has already taken place.
   The laplace treatment works for higher order elements without significant
   changes.

There are no loads possible for these elements. Only Dirichlet Boundary
Conditions are considered. They are accounted for as an additional load
vector.

The Ale module also contains control routines for the academic type of a 'pure
ale problem' where Dirichlet boundary conditions are prescribed rather than
obtained from an fluid or structure displacement.

The implementation of the ale problem includes also an ale element quality
monitoring. The QUALITY criterion to monitor can be specified in the input to
- aspect_ratio
- corner_angle
- min_J
- none

All quality criterions are normaised to one for a perfectly shaped element and
decrease to zero for invalid elements. The monitoring plots its results into a
.plt file, which can be used to produce gnuplot plots of the quality behaviour.
The .plt file contains: the actual step, the average quality over all elements,
the standard degression, the minimal element quality and the maximal element
quality.

At the moment quality monitoring is implemented for quad4-elements only.

Quality monitoring has been extended to higher order elements but is in this
case based on estimated values rather than exact ones.
*/
