/*----------------------------------------------------------------------*/
/*! \file

\brief A service method allowing the application of initial conditions
       for nurbs discretisations.

\maintainer Martin Kronbichler

\level 2
*/
/*----------------------------------------------------------------------*/
#include "drt_apply_nurbs_initial_condition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"



/*----------------------------------------------------------------------*/
/*!
   A service method allowing the application of initial conditions
   for nurbs discretisations. recommended version with allocation
   of a separate solver!
*/
/*----------------------------------------------------------------------*/
void DRT::NURBS::apply_nurbs_initial_condition(DRT::Discretization& dis, FILE* outfile,
    const Teuchos::ParameterList& solverparams, const int startfuncno,
    Teuchos::RCP<Epetra_Vector> initialvals)
{
  // try to cast dis to a nurbs discretisation --- if possible, proceed
  // with setting initial conditions. Otherwise return.
  DRT::NURBS::NurbsDiscretization* nurbsdis = dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&dis);

  if (nurbsdis == NULL)
  {
    return;
  }

  // Owing to experience a very accurate solution has to be enforced here!
  // Thus, we allocate an own solver with VERY strict tolerance!
  Teuchos::ParameterList p(solverparams);
  const double origtol = p.get<double>("AZTOL");
  const double newtol = 1.0e-11;
  p.set("AZTOL", newtol);

  Teuchos::RCP<LINALG::Solver> lssolver = Teuchos::rcp(new LINALG::Solver(p, dis.Comm(), outfile));
  dis.ComputeNullSpaceIfNecessary(lssolver->Params());

  // get the processor ID from the communicator
  const int myrank = dis.Comm().MyPID();

  if (myrank == 0)
    std::cout << "\nSolver tolerance for least squares problem set to " << newtol
              << " (orig: " << origtol << ")";

  apply_nurbs_initial_condition_solve(dis, *lssolver, startfuncno, initialvals);
  return;
}


/*----------------------------------------------------------------------*/
/*!
   A service method allowing the application of initial conditions
                    for nurbs discretisations.
*/
/*----------------------------------------------------------------------*/
void DRT::NURBS::apply_nurbs_initial_condition_solve(DRT::Discretization& dis,
    LINALG::Solver& solver, const int startfuncno, Teuchos::RCP<Epetra_Vector> initialvals)
{
  // try to cast dis to a nurbs discretisation --- if possible, proceed
  // with setting initial conditions. Otherwise return.
  DRT::NURBS::NurbsDiscretization* nurbsdis = dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&dis);

  if (nurbsdis == NULL)
  {
    return;
  }

  // get the knotvector from nurbs discretisation
  Teuchos::RCP<DRT::NURBS::Knotvector> knots = nurbsdis->GetKnotVector();

  // get the processor ID from the communicator
  const int myrank = dis.Comm().MyPID();

  if (myrank == 0)
  {
    printf("\n");
    printf("Setting up least-squares Nurbs approximation of initial field (discretization %s)\n",
        dis.Name().c_str());
  }

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = dis.DofRowMap();

  // -------------------------------------------------------------------
  // create empty mass matrix
  // -------------------------------------------------------------------
  Teuchos::RCP<LINALG::SparseMatrix> massmatrix =
      Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap, 108, false, true));

  // -------------------------------------------------------------------
  // create empty right hand side vector
  // -------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> rhs = LINALG::CreateVector(*dofrowmap, true);

  // -------------------------------------------------------------------
  // call elements to calculate massmatrix and righthandside
  // -------------------------------------------------------------------
  {
    // call elements and assemble
    if (!nurbsdis->Filled()) dserror("FillComplete() was not called");
    if (!nurbsdis->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

    // see what we have for input
    bool assemblemat = massmatrix != Teuchos::null;
    bool assemblevec = rhs != Teuchos::null;

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elemass;
    Epetra_SerialDenseVector elerhs;

    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;

    // loop over column elements
    const int numcolele = nurbsdis->NumMyColElements();

    int every = numcolele / 58;
    // prevent division by zero when dividing by every later on
    if (every < 1) every = 1;

    for (int i = 0; i < numcolele; ++i)
    {
      // first baci progress bar
      if (myrank == 0 && i % every == 0)
      {
        printf(".");
        fflush(0);
      }

      DRT::Element* actele = nurbsdis->lColElement(i);

      // get element location vector, dirichlet flags and ownerships
      lm.clear();
      lmowner.clear();
      actele->LocationVector(*nurbsdis, lm, lmowner, lmstride);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      const int eledim = (int)lm.size();

      if (assemblemat)
      {
        if (elemass.M() != eledim or elemass.N() != eledim)
          elemass.Shape(eledim, eledim);
        else
          memset(elemass.A(), 0, eledim * eledim * sizeof(double));
      }
      if (assemblevec)
      {
        if (elerhs.Length() != eledim)
          elerhs.Size(eledim);
        else
          memset(elerhs.Values(), 0, eledim * sizeof(double));
      }

      {
        int spacedim = -1;

        const DRT::Element::DiscretizationType distype = actele->Shape();
        switch (distype)
        {
          case DRT::Element::nurbs4:
          case DRT::Element::nurbs9:
          {
            spacedim = 2;
            break;
          }
          case DRT::Element::nurbs8:
          case DRT::Element::nurbs27:
          {
            spacedim = 3;
            break;
          }
          default:
            dserror("this method is designed for usage with NurbsDiscretization only");
            break;
        }

        // set element data
        const int iel = actele->NumNode();

        // dofblocks (spacedim or spacedim+1 for fluid problems)
        const int dofblock = eledim / iel;

        // get node coordinates of element
        Epetra_SerialDenseMatrix xyze(spacedim, iel);
        DRT::Node** nodes = actele->Nodes();
        for (int inode = 0; inode < iel; inode++)
        {
          const double* x = nodes[inode]->X();
          for (int dim = 0; dim < spacedim; ++dim)
          {
            xyze(dim, inode) = x[dim];
          }
        }

        // aquire weights from nodes
        Epetra_SerialDenseVector weights(iel);

        for (int inode = 0; inode < iel; ++inode)
        {
          DRT::NURBS::ControlPoint* cp = dynamic_cast<DRT::NURBS::ControlPoint*>(nodes[inode]);

          weights(inode) = cp->W();
        }

        // access elements knot span
        std::vector<Epetra_SerialDenseVector> eleknots(spacedim);

        bool zero_size = false;
        zero_size = knots->GetEleKnots(eleknots, actele->Id());

        // nothing to be done for a zero sized element
        if (zero_size)
        {
          continue;
        }

        Epetra_SerialDenseVector funct(iel);
        Epetra_SerialDenseMatrix xjm(spacedim, spacedim);
        Epetra_SerialDenseMatrix deriv(spacedim, iel);
        Epetra_SerialDenseVector gp(spacedim);
        Epetra_SerialDenseVector position(
            3);  // always three-dimensional coordinates for function evaluation!
        Epetra_SerialDenseVector initialval(dofblock);

        // depending on the spatial dimension, we need a different
        // integration scheme
        switch (spacedim)
        {
          case 2:
          {
            // gaussian points
            const DRT::UTILS::IntegrationPoints2D intpoints(DRT::UTILS::intrule_quad_9point);

            for (int iquad = 0; iquad < intpoints.nquad; ++iquad)
            {
              // set gauss point coordinates
              for (int rr = 0; rr < spacedim; ++rr)
              {
                gp(rr) = intpoints.qxg[iquad][rr];
              }

              DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv(
                  funct, deriv, gp, eleknots, weights, distype);

              // get transposed Jacobian matrix and determinant
              //
              //        +-            -+ T      +-            -+
              //        | dx   dx   dx |        | dx   dy   dz |
              //        | --   --   -- |        | --   --   -- |
              //        | dr   ds   dt |        | dr   dr   dr |
              //        |              |        |              |
              //        | dy   dy   dy |        | dx   dy   dz |
              //        | --   --   -- |   =    | --   --   -- |
              //        | dr   ds   dt |        | ds   ds   ds |
              //        |              |        |              |
              //        | dz   dz   dz |        | dx   dy   dz |
              //        | --   --   -- |        | --   --   -- |
              //        | dr   ds   dt |        | dt   dt   dt |
              //        +-            -+        +-            -+
              //
              // The Jacobian is computed using the formula
              //
              //            +-----
              //   dx_j(r)   \      dN_k(r)
              //   -------  = +     ------- * (x_j)_k
              //    dr_i     /       dr_i       |
              //            +-----    |         |
              //            node k    |         |
              //                  derivative    |
              //                   of shape     |
              //                   function     |
              //                           component of
              //                          node coordinate
              //
              for (int rr = 0; rr < spacedim; ++rr)
              {
                for (int mm = 0; mm < spacedim; ++mm)
                {
                  xjm(rr, mm) = deriv(rr, 0) * xyze(mm, 0);
                  for (int nn = 1; nn < iel; ++nn)
                  {
                    xjm(rr, mm) += deriv(rr, nn) * xyze(mm, nn);
                  }
                }
              }

              // The determinant is computed using Sarrus's rule
              const double det = xjm(0, 0) * xjm(1, 1) - xjm(0, 1) * xjm(1, 0);

              // get real physical coordinates of integration point
              /*
              //              +-----
              //               \
              //    pos (x) =   +      N (x) * x
              //               /        j       j
              //              +-----
              //              node j
              */
              for (int rr = 0; rr < spacedim; ++rr)
              {
                position(rr) = funct(0) * xyze(rr, 0);
                for (int mm = 1; mm < iel; ++mm)
                {
                  position(rr) += funct(mm) * xyze(rr, mm);
                }
              }
              // if spacedim < 3, ensure we define a valid z-coordinate!
              for (int rr = spacedim; rr < 3; ++rr) position(rr) = 0.0;

              for (int rr = 0; rr < dofblock; ++rr)
              {
                // important: position has to have always three components!!
                initialval(rr) = DRT::Problem::Instance()
                                     ->Funct(startfuncno - 1)
                                     .Evaluate(rr, position.Values(), 0.0);
              }


              // check for degenerated elements
              if (det < 1E-16)
              {
                dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f",
                    actele->Id(), det);
              }

              // set total integration factor
              double fac = intpoints.qwgt[iquad] * det;

              for (int vi = 0; vi < iel; ++vi)  // loop rows  (test functions)
              {
                const int fvi = dofblock * vi;

                for (int ui = 0; ui < iel; ++ui)  // loop columns  (test functions)
                {
                  const int fui = dofblock * ui;

                  const double diag = fac * funct(ui) * funct(vi);

                  for (int rr = 0; rr < dofblock; ++rr)
                  {
                    elemass(fvi + rr, fui + rr) += diag;
                  }
                }
                for (int rr = 0; rr < dofblock; ++rr)
                {
                  elerhs(fvi + rr) += fac * funct(vi) * initialval(rr);
                }
              }
            }  // end gaussloop
            break;
          }
          case 3:
          {
            // gaussian points
            const DRT::UTILS::IntegrationPoints3D intpoints(DRT::UTILS::intrule_hex_27point);

            for (int iquad = 0; iquad < intpoints.nquad; ++iquad)
            {
              // set gauss point coordinates
              for (int rr = 0; rr < spacedim; ++rr)
              {
                gp(rr) = intpoints.qxg[iquad][rr];
              }

              DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv(
                  funct, deriv, gp, eleknots, weights, distype);

              // get transposed Jacobian matrix and determinant
              //
              //        +-            -+ T      +-            -+
              //        | dx   dx   dx |        | dx   dy   dz |
              //        | --   --   -- |        | --   --   -- |
              //        | dr   ds   dt |        | dr   dr   dr |
              //        |              |        |              |
              //        | dy   dy   dy |        | dx   dy   dz |
              //        | --   --   -- |   =    | --   --   -- |
              //        | dr   ds   dt |        | ds   ds   ds |
              //        |              |        |              |
              //        | dz   dz   dz |        | dx   dy   dz |
              //        | --   --   -- |        | --   --   -- |
              //        | dr   ds   dt |        | dt   dt   dt |
              //        +-            -+        +-            -+
              //
              // The Jacobian is computed using the formula
              //
              //            +-----
              //   dx_j(r)   \      dN_k(r)
              //   -------  = +     ------- * (x_j)_k
              //    dr_i     /       dr_i       |
              //            +-----    |         |
              //            node k    |         |
              //                  derivative    |
              //                   of shape     |
              //                   function     |
              //                           component of
              //                          node coordinate
              //
              for (int rr = 0; rr < spacedim; ++rr)
              {
                for (int mm = 0; mm < spacedim; ++mm)
                {
                  xjm(rr, mm) = deriv(rr, 0) * xyze(mm, 0);
                  for (int nn = 1; nn < iel; ++nn)
                  {
                    xjm(rr, mm) += deriv(rr, nn) * xyze(mm, nn);
                  }
                }
              }

              // The determinant is computed using Sarrus's rule
              const double det =
                  xjm(0, 0) * xjm(1, 1) * xjm(2, 2) + xjm(2, 0) * xjm(0, 1) * xjm(1, 2) +
                  xjm(0, 2) * (xjm(1, 0) * xjm(2, 1) - xjm(2, 0) * xjm(1, 1)) -
                  xjm(2, 1) * xjm(0, 0) * xjm(1, 2) - xjm(0, 1) * xjm(1, 0) * xjm(2, 2);


              // get real physical coordinates of integration point
              /*
              //              +-----
              //               \
              //    pos (x) =   +      N (x) * x
              //               /        j       j
              //              +-----
              //              node j
              */
              for (int rr = 0; rr < spacedim; ++rr)
              {
                position(rr) = funct(0) * xyze(rr, 0);
                for (int mm = 1; mm < iel; ++mm)
                {
                  position(rr) += funct(mm) * xyze(rr, mm);
                }
              }
              // if spacedim < 3, ensure we define a valid z-coordinate!
              for (int rr = spacedim; rr < 3; ++rr) position(rr) = 0.0;

              for (int rr = 0; rr < dofblock; ++rr)
              {
                // important: position has to have always three components!!
                initialval(rr) = DRT::Problem::Instance()
                                     ->Funct(startfuncno - 1)
                                     .Evaluate(rr, position.Values(), 0.0);
              }

              // check for degenerated elements
              if (det < 0.0)
              {
                dserror(
                    "GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", actele->Id(), det);
              }

              // set total integration factor
              double fac = intpoints.qwgt[iquad] * det;

              for (int vi = 0; vi < iel; ++vi)  // loop rows  (test functions)
              {
                const int fvi = dofblock * vi;

                for (int ui = 0; ui < iel; ++ui)  // loop columns  (test functions)
                {
                  const int fui = dofblock * ui;

                  const double diag = fac * funct(ui) * funct(vi);

                  for (int rr = 0; rr < dofblock; ++rr)
                  {
                    elemass(fvi + rr, fui + rr) += diag;
                  }
                }
                for (int rr = 0; rr < dofblock; ++rr)
                {
                  elerhs(fvi + rr) += fac * funct(vi) * initialval(rr);
                }
              }
            }  // end gaussloop
            break;
          }
          default:
            dserror("expecting two or three-dimensional problems to set the initial conditions\n");
            break;
        }
      }

      int eid = actele->Id();
      if (assemblemat) massmatrix->Assemble(eid, elemass, lm, lmowner);
      if (assemblevec) LINALG::Assemble(*rhs, elerhs, lm, lmowner);
    }  // for (int i=0; i<numcolele; ++i)
  }

  if (myrank == 0)
  {
    printf("\n");
    printf("\n");
    printf("Solving least-squares problem: ");
  }

  // -------------------------------------------------------------------
  // finalize the system matrix
  // -------------------------------------------------------------------
  massmatrix->Complete();

  // -------------------------------------------------------------------
  // solve system
  // -------------------------------------------------------------------

  // always refactor and reset the matrix before a single new solver call
  bool refactor = true;
  bool reset = true;

  initialvals->PutScalar(0.0);

  solver.Solve(massmatrix->EpetraOperator(), initialvals, rhs, refactor, reset);

  // perform resets for solver and matrix
  solver.Reset();
  massmatrix->Reset();

  if (myrank == 0)
  {
    printf("Initial field set.\n\n");
  }

  return;
}
