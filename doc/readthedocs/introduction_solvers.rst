Linear solver
================


The heart of a finite element solver is the linear solver, since all equations used to describe the behavior of a structure are finally collected in a single matrix equation of the form

.. math::

    \mathbf{A} x =  b

where the vector :math:`x` is to be calculated based on a known matrix :math:`\mathbf{A}` and the vector :math:`b`. Two numerical methods exist in general to solve this equation.

* Direct solver
* Iterative solver

It is very important to understand that a direct solver simply solves the equation without consideration of the nature of the system, while physics plays a role in iterative methods, and particular, their multigrid preconditioners (see below), so it is important for the latter to comprehend the parameters to be used.
Anyway, the decision on which solver to use is illustrated in the following flow chart:

.. figure:: figures/LinearSolversInBaci_flowchart.png
   :width: 400px
   :align: center

   Flowchart of the decision process for an appropriate linear solver, taken from a presentation held by Max Firmbach in 2022.

At this point, the linear algebra library *Trilinos* is heavily used, which provides a number of packages on different levels:

**Epetra/TPetra**: 
   Sparse Linear Algebra

**Amesos/Amesos2**: 
   Direct Solvers (UMFPACK, Superlu)

**AztecOO/Belos**: 
   Iterative Solvers (CG, GMRES)

**Ifpack**: 
   Preconditioner (ILU, ICHOL, GS)

**ML/MueLu**: 
   Multigrid-Preconditioner


The solvers are called in the solver sections::

   -----------SOLVER 1
   NAME    <arbitrary_name>
   SOLVER  <solver>
   ...
   -----------SOLVER 9

in which all details for a specific solver are defined. The parameter :ref:`SOLVER <solver1_solver>` defines the solver to be used. Most of the solvers require additional parameters. These are explained in the following. 



Direct solver
-------------

In principle, the direct solver identifies :math:`x` by calculating the inverse of :math:`\mathbf{A}`: :math:`x = \mathbf{A}^{-1} b`.

In BACI, we have three different direct solvers included:

   * UMFPACK, using a multifront approach, thus a sequential solver (can only use a single core)
   * SuperLU, an MPI based solver, which may run on many cores and compute nodes.

Compared to iterative solvers, these solvers do not scale well with the numbers of equations, and are therefore not well suited for large systems.
If one has to solve a system with more than 50000 degrees of freedom (approx. equal to number of equations), the iterative solver will be significantly faster.
In addition, the iterative solver is more memory efficient, so it can solve larger system with a computer equipped with low memory.

The benefit of the direct solver is that there are no parameters, which one has to provide, 
since the direct solver does not care about the underlying physics.

Iterative solver
-----------------

The iterative solver can be used for any size of equation systems, but is the more efficient the larger the problem is.
If a good parameter set for the solver is chosen, it scales linearly with the size of the system, either with respect to time or to the number of cores on which the system is solved.

The main drawback is that one has to provide a number of parameters, which are crucial for a fast and correct solution.

Contrary to the direct solver, the matrix inverse is not needed.
Instead, this solution method solves the equation :math:`\mathbf{A} x_k = b`  with an initial guess :math:`x_0 (k=0)` and an iteration

.. math::

   x_{k+1} = \mathbf{P}(x_k, \mathbf{A} x_k, b) \, ,

such that :math:`x_k \rightarrow x \mbox{ for } k \rightarrow \infty`.
Slow progress if :math:`x_0` is not chosen properly. Preconditioner helps by solving
:math:`\mathbf{P}^{-1} \mathbf{A} x = \mathbf{P}^{-1} b`. Ideally :math:`\mathbf{P}^{-1} = \mathbf{A}` (gives the solution for *x*), 
but :math:`\mathbf{P}` should be cheap to calculate.
The larger the problem is, the higher is the benefit of iterative solvers.

Solver interface
^^^^^^^^^^^^^^^^

The current solver is based on Trilinos' **Belos** package, which is the successor of AztecOO. This package provides a bunch of KRYLOV solvers, e.g.

   - CG for symmetric problems, 
   - GMRES (also for unsymmetric problems)
   - BICGSTAB (??)

The parameters are very similar to the AztecOO package, and therefore further information on specific details can be found in the official `AztecOO user guide <https://trilinos.github.io/pdfs/AztecOOUserGuide.pdf>`_. Since an iterative method approximates the correct solution iteratively, 
a residual criterion (here :ref:`AZCONV <solver1_azconv>`), 
a tolerance (:ref:`AZTOL <solver1_aztol>`), 
and the maximum number of iterations :ref:`AZITER <solver1_aziter>` must be given. 
The default of AZCONV (AZ_r0) is reasonable, while the default for the maximum number of iterations (1000) is rather high.

In addition, the parameter :ref:`AZSUB <solver1_azsub>` is important, if the default algorithm, GMRES, is used: 
Krylov solvers build up a subspace of vectors, and they should be rebuilt after a number of iterations. The default of 50 is reasonable.

These parameters are defined in the following way:

::

   -----------SOLVER 1
   SOLVER  Belos
   AZSOLVE [CG|GMRES|BICGSTAB]
   AZCONV  [AZ_r0|AZ_rhs|AZ_Anorm|AZ_noscaled|AZ_sol|AZ_weighted]
   AZITER  <number>



Preconditioners
^^^^^^^^^^^^^^^^

The choice and design of the preconditioner highly affect performance. Within Belos one can choose between 

-	ILU
-	Algebraic Multigrid (AMG) methods

**ILU** (incomplete LU method) 

Perfect scalability is not achieved with this method, but is has the advantage of being less complex.

ILU comes with a single parameter: :ref:`IFPACKGFILL <solver1_ifpackgfill>`.
The default, ``IFPACKGFILL 0``, will not include further elements in the preconditioner P (same sparcity pattern)
setup will be faster, approximation is worse
``IFPACKGFILL = 1 .. n``: The higher the more elements are included, sparcity decreases (a level of 12 might be a full matrix, like a direct solver)
**Remark** One should probably not go beyond 3, maybe start with 0) Only for very strong ill-conditioning one should go towards 3.

::

   -----------SOLVER 1
   AZPREC       ILU
   IFPACKGFILL  [0 .. 12]

**Algebraic Multigrid (AMG) methods**

The current recommendation is the trilinos ML preconditioner, for further information on this, see [Gee07]_.

*Theory:*

The trick is to apply a cheap transfer method to get from the complete system to a smaller one (coarsening/aggregation of the system). The method to be used is given in :ref:`ML_COARSEN <solver1_ml_coarsen>`. In order to get a preconditioner matrix of the same size as the original matrix, of course, the aggregation must be used in the opposite direction afterwards. The default method is UC (uncoupled), which is a good choice as well.

The coarsening reduces the size by a factor :ref:`ML_AGG_SIZE <solver1_ml_agg_size>`. It defines how many lines are comprised to one (good choices are 27 for 3D, 9 for 2D, and 3 for 1D problems).

A smoother is used twice (pre- and post-smoother) for each level of aggregation to reduce the error frequencies in your solution vector. Multiple transfer operations are applied in sequence, since only high frequency components can be tackled by smoothing, while the low frequency errors are still there. The restriction operator restricts the current error to the coarser grid. 
At some point (let say if 10000 dofs are left) the system has a size where one can apply the direct solver. This number is given by :ref:`ML_MAXCOARSESIZE <solver1_ml_maxcoarsesize>`. That is, when the number of remaining dofs is smaller than ML_MAXCOARSESIZE, no more coarsening is conducted. It should be larger than the default of 1000, let say, 5000-10000. Also, the maximum number of coarsenings is given by :ref:`ML_MAXLEVEL <solver1_ml_maxlevel>` (maxlevel should always be high enough). 

One may define three different smoothers:
:ref:`ML_SMOOTHERFINE <solver1_ml_smootherfine>` (for the first / fine level); 
:ref:`Ml_SMOOTHERMED <solver1_ml_smoothermed>` (for all intermediate levels); 
:ref:`ML_SMOOTHERCOARSE <solver1_ml_smoothercoarse>` (probably always a direct solver like UMFPACK). 

While many solvers can be used, five of them are most popular: SGS (symmetric Gauss Seidel), Jacobi, Chebychev, ILU, MLS. Besides that, particularly for the coarsest smoother, a direct solver can be used, as (Umfpack, SuperLU, KLU).

*Chebychev smoother:*
   This is a polynomial smoother. The degree of the polynomial is given by :ref:`ML_SMOTIMES <solver1_ml_smotimes>`. A lower degree is faster (not much), but higher is more accurate; one may use 3, 6 or even 9 [very high])

*Relaxation method (e.g. SGS):*
   For this kind of smoothers, :ref:`ML_SMOTIMES <solver1_ml_smotimes>` is providing the number of sweeps for each smoothening. This one is rather for fluid dynamics problems.

*ILU:*
   Here, :ref:`ML_SMOTIMES <solver1_ml_smotimes>` will be interpreted as the FILL level.

Damping helps with convergence, and it can be appliedto any of the smoothers, see :ref:`ML_DAMPFINE <solver1_ml_dampfine>`, :ref:`ML_DAMPMED <solver1_ml_dampmed>`, :ref:`ML_DAMPCOARSE <solver1_ml_dampcoarse>`. a vlues of 1 cancels damping, 0 means maximum damping. Too much damping increases the iterations, thus, usually it should be between 1 and 0.5. A little bit of damping will probably improve convergence (also from the beginning).

:ref:`ML_PROLONG_SMO <solver1_ml_prolong_smo>` is the main parameter to control the prolongation. Transfer operator from coarse to fine (a tentative prolongator is created by constant interpolation, then try to improve the constant to linear interpolation). 1.33 is a good value for structural, scatra, thermo problems. Different choice would be 1 (maybe for fluids).

.. list-table::
   header-rows: 1

   * - Problem
     - Symmetry
   * - Convection dominated flow
     - nonsymm
   * - elasticity 
     - symm
   * - Contact
     - unsymm


Coupled problems:
^^^^^^^^^^^^^^^^^^

If a :ref:`multiphysics problem <multifieldproblems>` is to be solved, they can be solved sequentially, and the interaction then leads to an iterative procedure, where the influence of one field to the other hopefully converges to a common solution.
From a solver's point of view, the solution is achieved by running the solver of each single field problem independent of others. Therefore, those problems are solved using the methods given above.

If, on the other hand, the interaction of the physical fields is strong, the iterative procedure may converge only slowly (if at all), thus a monolithic solution, solving all degrees of freedom simultaneously is beneficial. This so-called *monolithic solution* will be described in the following:

**Monolithic solution:** all degrees of freedom appear in the linear system. 
Since the stiffness factors of the different physics may be different by orders of magnitude, and the coupling between the physics may again have a different magnitude, the linear system may be particularly ill-conditioned.
On the other hand, **sequential solutions** are handled like single field problems from the solver point of view.

*Monolithic solvers for multiphysics problems:*

Contact with penalty: basically still solid mechanics (probably a bit more ill-conditioned)
Contact with lagrange multipliers and other problems: block structure in the system, and the preconditioner needs to know about it.
If you have a block structure (e.g. TSI monolithic): e.g. Block Gauss Seidel 2 by 2 (BGS2x2)
one has to invert the d
How to connect the different blocks together. In BGS you need to invert the different blocks, you need to provide for preconditioner for each block, which is specified in the other solver sections. You need to define a preconditioner in the structural and the thermo solver.

Future prospects
^^^^^^^^^^^^^^^^

**Switch from ML to MueLu**

   In the near future there will be a major change with respect to the preconditioners. This will affect also the th einput parameters, and even the input style, which will then rather depend on separate solver parameter files than parameters in the .dat file.

Further reading
^^^^^^^^^^^^^^^

.. figure:: figures/TGM_LinearSolversInBaci.pdf
   :width: 400px
   :align: center

   A presentation held by Max Firmbach in 2022 


