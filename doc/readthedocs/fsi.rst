old_Fluid-Structure Interaction
===============================

Overview and Problem Definition
-------------------------------

.. figure:: figures/coupled_system_f_s.pdf
   :alt: Coupled FSI problem with fluid and structure domain and a common interface. In the boxes the main variables of each field are given.
   :align: center
   :width: 100%

   Coupled FSI problem with fluid and structure domain and a common
   interface. In the boxes the main variables of each field are
   given.

Time integration
^^^^^^^^^^^^^^^^

Only implicit time integration implemented for structure and fluid

-  Structure – Gen. Alpha.

-  Fluid – one-step-:math:`\theta`, BDF2

They algorithms in ``fsi_fluid(), fsi_structure()`` are copies of
existing single field time integrators.

ALE and Fixed-Grid Methods
--------------------------

Currently, only Arbitrary Lagrangian Eulerian (ALE) methods are
implemented, tested and in use. Fixed grid methods are in preparation
and will be added, once they are ready for general usage.

Partitioned algorithms
----------------------

One has the choice between several iterative methods and sequential
schemes, which however are not generaly applicable for incompressible
flow problems.

Fix-Point Iteration
^^^^^^^^^^^^^^^^^^^

Fixpunkt

::

   fsidyn->ifsi==fsi_iter_stagg_fixed_rel_param || (Verschiebungsrelaxation)
   fsidyn->ifsi==fsi_iter_stagg_AITKEN_rel_param ||
   fsidyn->ifsi==fsi_iter_stagg_steep_desc ||
   fsidyn->ifsi==fsi_iter_stagg_AITKEN_rel_force || (Kraftiteration Uli)
   fsidyn->ifsi==fsi_iter_stagg_steep_desc_force ||

Newton Iterations
^^^^^^^^^^^^^^^^^

(Diplomarbeit Markus Schmidtberger)

::

   fsidyn->ifsi==fsi_iter_stagg_Newton_FD || (Finite Differenzen)
   fsidyn->ifsi==fsi_iter_stagg_Newton_I )   (Fixpunkt ohne Relax)

unused
^^^^^^

::

   fsidyn->ifsi==fsi_iter_stagg_CHEB_rel_param ||

Sequential Algorithms
^^^^^^^^^^^^^^^^^^^^^

Sequential algorithms do not requeire an iteration which make them
cheap. However they have a severe time step restriction and might not be
applicable for interation involving incompressible flow.

::

   basic sequentiel staggered scheme
   sequential staggered scheme with predictor

Coupling for Matching and Non-Matching Grids
--------------------------------------------

Matching Grids
^^^^^^^^^^^^^^

Matching grids means that fluid and structure nodes are located at the
very same position along the interface und can directly exchange values
like velocities or forces.

Calculation of surface forces in 2d and 3d is done either using the
integration

.. math:: \int_\Omega \sigma \vec{n} \mathrm{d}\Omega

or by using consistent nodal forces following Förster

Non-Matching Grids
^^^^^^^^^^^^^^^^^^

Mortar Methods
^^^^^^^^^^^^^^

Mortar (Semesterarbeit Firl, Betr. Förster)

the implemented version compiles, but the test fails (Structur: No
convergence in maxiter steps)

It should be implemented/reorganized again using the Trilinos library
Moertel.

Interpolation
^^^^^^^^^^^^^

Interpolation (Semesterarbeit Florian Henke (Kuettler)) Einschalten per
DEFINES Flag!! nur für ``hex8,hex27,tet4,tet10``

Single field solvers
--------------------

In the following all single field solvers that currwtnly work in the FSI
setting are briefly reviewed.

Fluid
^^^^^

The working horse so far has been the stabilized equal order Q1Q1 and
Q2Q2 finite element fluid solver.

There are also fast Elemente of type ``FLUID3D_FAST Elemente``

New approaches to be included into the FSI algorithms are inf-sup stable
elements and elements using projection methods (Küttler)

Structure
^^^^^^^^^

Wall und Brick Elemente, Shell8

ALE mesh dynamics
^^^^^^^^^^^^^^^^^

The mesh dynamics algorithms are necessary to compute the fluid mesh
movement in the interior of the fluid domain.

It’s possible to run the pure ALE mesh to test ‘mesh materials’
(``ale_dyn_control.c``)

Materials & Element types
"""""""""""""""""""""""""

2D
""

Springs, Laplace, Linear FE, Linear FE mit versteifendem J, 2step
Elementtypen QUAD4, QUAD8, QUAD9, TRI3, TRI6

.. _d-1:

3D
""

Lin. Elastisch Elementtypen HEX8, HEX20, TET4, TET10

Tests
-----

FSI Tests
^^^^^^^^^

::

   ffsi_ow3D.dat
   ffsi_ow3D_usfem_ca.dat
   ffsi_ow3D_usfem_ca_sub2.dat
   ffsi_ow3D_usfem_ca_sub2_sd2.dat
   ffsi_ow3D_usfem.dat
   ffsi_ow3D_usfem_sub2.dat
   fsi_mtr_dc4x4.dat                           --> Mortar Kopplung -->
   fsi_ow32x32.dat
   fsi_ow32x32_force.dat
   fsi_ow32x32_sd.dat
   fsi_ow32x32_usfem_2xaztec.dat
   fsi_ow32x32_usfem.dat
   fsi_ow3D_usfem.dat
   fsi_tank20x10.dat
   fsi_tank3D.dat
   fsi_tank3D_par.dat

ALE tests
^^^^^^^^^

::

   ale2_nln.dat
   ale2_nln_sub4.dat
   ale_3d_hourg.dat
   ale_3d_hourg_sub2.dat
   ale_3d_oll.dat
   ale_3d_sd2.dat
   ale_3d_sd.dat
   ale_3d_umfpack.dat
   ale_3d_umfpack_fastass.dat
