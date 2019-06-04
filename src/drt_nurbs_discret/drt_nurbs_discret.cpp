/*----------------------------------------------------------------------*/
/*!

\brief a class to manage one nurbs discretization

\level 1

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#include <Epetra_Vector.h>
#include <Epetra_Time.h>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_utils.H"

#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"

#include "drt_nurbs_utils.H"
#include "drt_nurbs_discret.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 05/08|
 *----------------------------------------------------------------------*/
DRT::NURBS::NurbsDiscretization::NurbsDiscretization(
    const std::string name, Teuchos::RCP<Epetra_Comm> comm)
    : DRT::Discretization::Discretization(name, comm), npatches_(0), knots_(Teuchos::null)
{
  return;
}


/*----------------------------------------------------------------------*
 |  add a knotvector to the discretization (public)          gammi 05/08|
 *----------------------------------------------------------------------*/
void DRT::NURBS::NurbsDiscretization::SetKnotVector(Teuchos::RCP<DRT::NURBS::Knotvector> knots)
{
  if (knots == Teuchos::null)
  {
    dserror(
        "You're trying to set an invalid knotvector in the DRT::NURBS::NurbsDiscretization "
        "'%s'. The given know vector is a null vector and can't be set as such.",
        (this->Name()).c_str());
  }

  knots_ = knots;
  filled_ = false;
  return;
}

/*----------------------------------------------------------------------*
 |  get a pointer to knotvector from the discretization         (public)|
 |                                                           gammi 05/08|
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::NURBS::Knotvector> DRT::NURBS::NurbsDiscretization::GetKnotVector()
{
  if (knots_ == Teuchos::null)
  {
    dserror(
        "You're trying to access the NURBS knot vector in the DRT::NURBS::NurbsDiscretization "
        "'%s'. The required knot vector is a null vector and can't be accessed as such.",
        (this->Name()).c_str());
  }
  return knots_;
}


/*----------------------------------------------------------------------*
 |  get a pointer to knotvector from the discretization         (public)|
 |  (const version, read-only)                               gammi 05/08|
 *----------------------------------------------------------------------*/
Teuchos::RCP<const DRT::NURBS::Knotvector> DRT::NURBS::NurbsDiscretization::GetKnotVector() const
{
  if (knots_ == Teuchos::null)
  {
    dserror(
        "You're trying to access the NURBS knot vector in the DRT::NURBS::NurbsDiscretization "
        "'%s'. The required knot vector is a null vector and can't be accessed as such.",
        (this->Name()).c_str());
  }
  return knots_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::DbcNurbs::Evaluate(const DRT::DiscretizationInterface& discret, const double& time,
    const Teuchos::RCP<Epetra_Vector>* systemvectors, Epetra_Vector& toggle,
    Teuchos::RCP<std::set<int>>* dbcgids) const
{
  // --------------------------- Step 1 ---------------------------------------
  DRT::UTILS::Dbc::Evaluate(discret, time, systemvectors, toggle, dbcgids);

  // --------------------------- Step 2 ---------------------------------------
  std::vector<std::string> dbc_cond_names(2, "");
  dbc_cond_names[0] = "Dirichlet";
  dbc_cond_names[1] = "NurbsLSDirichlet";

  std::vector<Teuchos::RCP<DRT::Condition>> conds(0);
  std::vector<Teuchos::RCP<DRT::Condition>> curr_conds(0);
  for (std::vector<std::string>::const_iterator cit_name = dbc_cond_names.begin();
       cit_name != dbc_cond_names.end(); ++cit_name)
  {
    discret.GetCondition(*cit_name, curr_conds);

    conds.reserve(conds.size() + curr_conds.size());
    std::copy(curr_conds.begin(), curr_conds.end(), std::back_inserter(conds));
  }

  ReadDirichletCondition(discret, conds, toggle, dbcgids);

  // --------------------------- Step 3 ---------------------------------------
  conds.clear();
  discret.GetCondition("NurbsLSDirichlet", conds);

  Teuchos::RCP<std::set<int>> dbcgids_nurbs[2] = {Teuchos::null, Teuchos::null};
  dbcgids_nurbs[set_row] = Teuchos::rcp<std::set<int>>(new std::set<int>());
  dbcgids_nurbs[set_col] = Teuchos::rcp<std::set<int>>(new std::set<int>());

  // create a new toggle vector with column layout
  const DRT::NURBS::NurbsDiscretization* discret_nurbs =
      dynamic_cast<const DRT::NURBS::NurbsDiscretization*>(&discret);
  if (not discret_nurbs) dserror("Dynamic cast failed!");

  // build dummy column toggle vector
  Epetra_Vector toggle_col = Epetra_Vector(*discret_nurbs->DofColMap());

  ReadDirichletCondition(discret, conds, toggle_col, dbcgids_nurbs);

  // --------------------------- Step 4 ---------------------------------------
  DoDirichletCondition(discret, conds, time, systemvectors, toggle_col, dbcgids_nurbs);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::UTILS::DbcNurbs::DoDirichletCondition(const DRT::DiscretizationInterface& discret,
    const DRT::Condition& cond, const double& time,
    const Teuchos::RCP<Epetra_Vector>* systemvectors, const Epetra_Vector& toggle,
    const Teuchos::RCP<std::set<int>>* dbcgids) const
{
  // default call
  if (dbcgids[set_col].is_null())
  {
    DRT::UTILS::Dbc::DoDirichletCondition(discret, cond, time, systemvectors, toggle, dbcgids);
    return;
  }

  Epetra_Time timer(discret.Comm());

  // get the processor ID from the communicator
  const int myrank = discret.Comm().MyPID();
  if (myrank == 0) std::cout << "calculating least squares Dirichlet condition in ... ";

  const DRT::NURBS::NurbsDiscretization& nurbs_dis =
      static_cast<const DRT::NURBS::NurbsDiscretization&>(discret);

  // create map extractor to always (re)build dbcmapextractor which is needed later
  Teuchos::RCP<LINALG::MapExtractor> auxdbcmapextractor = Teuchos::rcp(new LINALG::MapExtractor());
  {
    // build map of Dirichlet DOFs
    int nummyelements = 0;
    int* myglobalelements = NULL;
    std::vector<int> dbcgidsv;
    if (dbcgids[set_row]->size() > 0)
    {
      dbcgidsv.reserve(dbcgids[set_row]->size());
      dbcgidsv.assign(dbcgids[set_row]->begin(), dbcgids[set_row]->end());
      nummyelements = dbcgidsv.size();
      myglobalelements = &(dbcgidsv[0]);
    }
    Teuchos::RCP<Epetra_Map> dbcmap = Teuchos::rcp(new Epetra_Map(-1, nummyelements,
        myglobalelements, discret.DofRowMap()->IndexBase(), discret.DofRowMap()->Comm()));
    // build the map extractor of Dirichlet-conditioned and free DOFs
    *auxdbcmapextractor = LINALG::MapExtractor(*(discret.DofRowMap()), dbcmap);
  }

  // column map of all DOFs subjected to a least squares Dirichlet condition
  Teuchos::RCP<Epetra_Map> dbccolmap = Teuchos::null;
  {
    // build map of Dirichlet DOFs
    int nummyelements = 0;
    int* myglobalelements = NULL;
    std::vector<int> dbcgidsv;
    if (dbcgids[set_col]->size() > 0)
    {
      dbcgidsv.reserve(dbcgids[set_col]->size());
      dbcgidsv.assign(dbcgids[set_col]->begin(), dbcgids[set_col]->end());
      nummyelements = dbcgidsv.size();
      myglobalelements = &(dbcgidsv[0]);
    }
    dbccolmap = Teuchos::rcp(new Epetra_Map(-1, nummyelements, myglobalelements,
        nurbs_dis.DofColMap()->IndexBase(), discret.DofRowMap()->Comm()));
  }

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Teuchos::RCP<const Epetra_Map> dofrowmap = auxdbcmapextractor->CondMap();

  if (dofrowmap->NumGlobalElements() == 0) return;  // no dbc gids ->leave

  // read information from condition
  const std::vector<int>* nodeids = cond.Nodes();
  if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");

  const std::vector<int>* funct = cond.Get<std::vector<int>>("funct");
  const std::vector<double>* val = cond.Get<std::vector<double>>("val");


  // determine highest degree of time derivative
  // and first existent system vector to apply DBC to
  unsigned deg = 0;  // highest degree of requested time derivative
  Teuchos::RCP<Epetra_Vector> systemvectoraux = Teuchos::null;  // auxiliar system vector
  if (systemvectors[0] != Teuchos::null)
  {
    deg = 0;
    systemvectoraux = systemvectors[0];
  }
  if (systemvectors[1] != Teuchos::null)
  {
    deg = 1;
    if (systemvectoraux == Teuchos::null) systemvectoraux = systemvectors[1];
  }
  if (systemvectors[2] != Teuchos::null)
  {
    deg = 2;
    if (systemvectoraux == Teuchos::null) systemvectoraux = systemvectors[2];
  }
  dsassert(systemvectoraux != Teuchos::null, "At least one vector must be unequal to null");


  // -------------------------------------------------------------------
  // create empty mass matrix
  // -------------------------------------------------------------------
  Teuchos::RCP<LINALG::SparseMatrix> massmatrix =
      Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap, 108, false, true));

  // -------------------------------------------------------------------
  // create empty right hand side vector
  // -------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> rhs = LINALG::CreateVector(*dofrowmap, true);
  Teuchos::RCP<Epetra_Vector> dbcvector = LINALG::CreateVector(*dofrowmap, true);

  Teuchos::RCP<Epetra_Vector> rhsd = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> dbcvectord = Teuchos::null;
  if (systemvectors[1] != Teuchos::null)
  {
    rhsd = LINALG::CreateVector(*dofrowmap, true);
    dbcvectord = LINALG::CreateVector(*dofrowmap, true);
  }

  Teuchos::RCP<Epetra_Vector> rhsdd = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> dbcvectordd = Teuchos::null;
  if (systemvectors[2] != Teuchos::null)
  {
    rhsdd = LINALG::CreateVector(*dofrowmap, true);
    dbcvectordd = LINALG::CreateVector(*dofrowmap, true);
  }

  const bool assemblevecd = rhsd != Teuchos::null;
  const bool assemblevecdd = rhsdd != Teuchos::null;

  // -------------------------------------------------------------------
  // call elements to calculate massmatrix and righthandside
  // -------------------------------------------------------------------
  {
    // call elements and assemble
    if (!discret.Filled()) dserror("FillComplete() was not called");
    if (!discret.HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

    // see what we have for input
    bool assemblemat = massmatrix != Teuchos::null;
    bool assemblevec = rhs != Teuchos::null;

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elemass;
    // Epetra_SerialDenseVector elerhs;
    std::vector<Epetra_SerialDenseVector> elerhs(deg + 1);

    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;

    std::vector<int> lm_full;
    std::vector<int> lmowner_full;
    std::vector<int> lmstride_full;

    const std::map<int, Teuchos::RCP<DRT::Element>>& geom = cond.Geometry();
    std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator curr;
    for (curr = geom.begin(); curr != geom.end(); ++curr)
    {
      Teuchos::RCP<DRT::Element> actele = curr->second;

      static const int probdim = DRT::Problem::Instance()->NDim();
      const DRT::Element::DiscretizationType distype = actele->Shape();
      const int dim = DRT::UTILS::getDimension(distype);
      const bool isboundary = (dim != probdim);
      const int nen = DRT::UTILS::getNumberOfElementNodes(distype);

      // access elements knot span
      std::vector<Epetra_SerialDenseVector> eleknots(dim);
      Epetra_SerialDenseVector weights(nen);

      bool zero_size = false;
      if (isboundary)
      {
        Teuchos::RCP<DRT::FaceElement> faceele =
            Teuchos::rcp_dynamic_cast<DRT::FaceElement>(actele, true);
        double normalfac = 0.0;
        std::vector<Epetra_SerialDenseVector> pknots(probdim);
        zero_size = DRT::NURBS::GetKnotVectorAndWeightsForNurbsBoundary(actele.get(),
            faceele->FaceMasterNumber(), faceele->ParentElement()->Id(), discret, pknots, eleknots,
            weights, normalfac);
      }
      else
        zero_size = DRT::NURBS::GetMyNurbsKnotsAndWeights(discret, actele.get(), eleknots, weights);

      // nothing to be done for a zero sized element
      if (zero_size)
      {
        continue;
      }

      // get element full location vector, dirichlet flags and ownerships
      lm_full.clear();
      lmowner_full.clear();
      lmstride_full.clear();
      actele->LocationVector(nurbs_dis, lm_full, lmowner_full, lmstride_full);

      // we are only interested in DOFs with dirichlet condition, hence we compare the location
      // vector with the
      // drichlet condition map
      lm.clear();
      lmowner.clear();
      lmstride.clear();
      for (unsigned j = 0; j < lm_full.size(); ++j)
      {
        int gid = lm_full[j];
        if (dbccolmap->MyGID(gid))
        {
          lm.push_back(gid);
          lmowner.push_back(lmowner_full[j]);
          lmstride.push_back(lmstride_full[j]);
        }
      }

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
        for (unsigned i = 0; i < deg + 1; ++i)
        {
          if (elerhs[i].Length() != eledim)
            elerhs[i].Size(eledim);
          else
            memset(elerhs[i].Values(), 0, eledim * sizeof(double));
        }
      }

      if (isboundary) switch (distype)
        {
          case DRT::Element::nurbs2:
            FillMatrixAndRHSForLSDirichletBoundary<DRT::Element::nurbs2>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case DRT::Element::nurbs3:
            FillMatrixAndRHSForLSDirichletBoundary<DRT::Element::nurbs3>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case DRT::Element::nurbs4:
            FillMatrixAndRHSForLSDirichletBoundary<DRT::Element::nurbs4>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case DRT::Element::nurbs9:
            FillMatrixAndRHSForLSDirichletBoundary<DRT::Element::nurbs9>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case DRT::Element::nurbs8:
            FillMatrixAndRHSForLSDirichletBoundary<DRT::Element::nurbs8>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case DRT::Element::nurbs27:
            FillMatrixAndRHSForLSDirichletBoundary<DRT::Element::nurbs27>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          default:
            dserror("invalid element shape for least squares dirichlet evaluation: %s",
                DistypeToString(distype).c_str());
            break;
        }
      else
        switch (distype)
        {
          case DRT::Element::nurbs2:
            FillMatrixAndRHSForLSDirichletDomain<DRT::Element::nurbs2>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case DRT::Element::nurbs3:
            FillMatrixAndRHSForLSDirichletDomain<DRT::Element::nurbs3>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case DRT::Element::nurbs4:
            FillMatrixAndRHSForLSDirichletDomain<DRT::Element::nurbs4>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case DRT::Element::nurbs9:
            FillMatrixAndRHSForLSDirichletDomain<DRT::Element::nurbs9>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case DRT::Element::nurbs8:
            FillMatrixAndRHSForLSDirichletDomain<DRT::Element::nurbs8>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case DRT::Element::nurbs27:
            FillMatrixAndRHSForLSDirichletDomain<DRT::Element::nurbs27>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          default:
            dserror("invalid element shape for least squares dirichlet evaluation: %s",
                DistypeToString(distype).c_str());
            break;
        }

      int eid = actele->Id();
      if (assemblemat) massmatrix->Assemble(eid, elemass, lm, lmowner);
      if (assemblevec) LINALG::Assemble(*rhs, elerhs[0], lm, lmowner);
      if (assemblevecd) LINALG::Assemble(*rhsd, elerhs[1], lm, lmowner);
      if (assemblevecdd) LINALG::Assemble(*rhsdd, elerhs[2], lm, lmowner);
    }
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

  // Owing to experience a very accurate solution has to be enforced here!
  // Thus, we allocate an own solver with VERY strict tolerance!
  // One could think of verifiying an extra solver in the input file...
  Teuchos::ParameterList p = DRT::Problem::Instance()->UMFPACKSolverParams();
  //  const double origtol = p.get<double>("AZTOL");
  //  const double newtol  = 1.0e-11;
  //  p.set("AZTOL",newtol);

  //  if(myrank==0)
  //    cout<<"\nSolver tolerance for least squares problem set to "<<newtol<<"\n";

  Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(
      new LINALG::Solver(p, discret.Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));
  // FixMe actually the const qualifier could stay, if someone adds to each single
  // related ComputeNullSpace routine a "const"....
  const_cast<DRT::DiscretizationInterface&>(discret).ComputeNullSpaceIfNecessary(solver->Params());

  // solve for control point values
  solver->Solve(massmatrix->EpetraOperator(), dbcvector, rhs, refactor, reset);

  // solve for first derivatives in time
  if (assemblevecd) solver->Solve(massmatrix->EpetraOperator(), dbcvectord, rhsd, refactor, reset);

  // solve for second derivatives in time
  if (assemblevecdd)
    solver->Solve(massmatrix->EpetraOperator(), dbcvectordd, rhsdd, refactor, reset);

  // perform resets for solver and matrix
  solver->Reset();
  massmatrix->Reset();

  // insert nodal values to sysvec
  auxdbcmapextractor->InsertCondVector(dbcvector, systemvectors[0]);
  if (assemblevecd) auxdbcmapextractor->InsertCondVector(dbcvectord, systemvectors[1]);
  if (assemblevecdd) auxdbcmapextractor->InsertCondVector(dbcvectordd, systemvectors[2]);

  if (myrank == 0) std::cout << timer.ElapsedTime() << " seconds \n\n";

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                   vuong 08/14|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::UTILS::DbcNurbs::FillMatrixAndRHSForLSDirichletBoundary(Teuchos::RCP<DRT::Element> actele,
    const std::vector<Epetra_SerialDenseVector>* knots, const std::vector<int>& lm,
    const std::vector<int>* funct, const std::vector<double>* val, const unsigned deg,
    const double time, Epetra_SerialDenseMatrix& elemass,
    std::vector<Epetra_SerialDenseVector>& elerhs) const
{
  if (deg + 1 != elerhs.size())
    dserror("given degree of time derivative does not match number or rhs vectors!");

  static const int dim = DRT::UTILS::DisTypeToDim<distype>::dim;

  const int ndbcdofs = (int)lm.size();

  // set element data
  static const int nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  // dofblocks (number of DOFs with Dirichlet condition per node)
  const int dofblock = ndbcdofs / nen;

  // get node coordinates of element
  LINALG::Matrix<dim + 1, nen> xyze;
  DRT::Node** nodes = actele->Nodes();

  for (int inode = 0; inode < nen; inode++)
  {
    const double* x = nodes[inode]->X();
    for (int idim = 0; idim < dim + 1; ++idim)
    {
      xyze(idim, inode) = x[idim];
    }
  }

  // aquire weights from nodes
  Epetra_SerialDenseVector weights(nen);

  for (int inode = 0; inode < nen; ++inode)
  {
    DRT::NURBS::ControlPoint* cp = dynamic_cast<DRT::NURBS::ControlPoint*>(nodes[inode]);

    weights(inode) = cp->W();
  }

  // shape functions
  LINALG::Matrix<nen, 1> shpfunct;
  // coordinates of integration points in parameter space
  LINALG::Matrix<dim, 1> xsi;
  // first derivative of shape functions
  LINALG::Matrix<dim, nen> deriv;
  // coordinates of integration point in physical space
  Epetra_SerialDenseVector position(
      3);  // always three-dimensional coordinates for function evaluation!
  // auxiliary date container for dirichlet evaluation
  std::vector<Epetra_SerialDenseVector> value(deg + 1, dofblock);
  // unit normal on boundary element
  LINALG::Matrix<dim + 1, 1> unitnormal;

  // gaussian points
  const DRT::UTILS::IntPointsAndWeights<dim> intpoints(
      DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // integration factor
    double fac = 0.0;
    double drs = 0.0;

    DRT::UTILS::EvalShapeFuncAtBouIntPoint<distype>(
        shpfunct, deriv, fac, unitnormal, drs, xsi, xyze, intpoints, iquad, knots, &weights, true);

    // get real physical coordinates of integration point
    /*
    //              +-----
    //               \
    //    pos (x) =   +      N (x) * x
    //               /        j       j
    //              +-----
    //              node j
    */
    for (int rr = 0; rr < dim + 1; ++rr)
    {
      position(rr) = shpfunct(0) * xyze(rr, 0);
      for (int mm = 1; mm < nen; ++mm)
      {
        position(rr) += shpfunct(mm) * xyze(rr, mm);
      }
    }
    // if dim < 3, ensure we define a valid z-coordinate!
    for (int rr = dim + 1; rr < 3; ++rr) position(rr) = 0.0;

    for (int rr = 0; rr < dofblock; ++rr)
    {
      // factor given by FUNCTS
      std::vector<double> functimederivfac(deg + 1, 1.0);
      for (unsigned i = 1; i < (deg + 1); ++i) functimederivfac[i] = 0.0;

      const int funct_num = (*funct)[rr];
      if (funct_num > 0)
      {
        // important: position has to have always three components!!
        functimederivfac = DRT::Problem::Instance()
                               ->Funct((*funct)[rr] - 1)
                               .EvaluateTimeDerivative(rr, position.Values(), time, deg);
      }

      // apply factors to Dirichlet value
      for (unsigned i = 0; i < deg + 1; ++i)
      {
        value[i](rr) = (*val)[rr] * functimederivfac[i];
      }
    }

    for (int vi = 0; vi < nen; ++vi)  // loop rows  (test functions)
    {
      const int fvi = dofblock * vi;

      for (int ui = 0; ui < nen; ++ui)  // loop columns  (test functions)
      {
        const int fui = dofblock * ui;

        const double diag = fac * shpfunct(ui) * shpfunct(vi);

        for (int rr = 0; rr < dofblock; ++rr)
        {
          elemass(fvi + rr, fui + rr) += diag;
        }
      }
      for (int rr = 0; rr < dofblock; ++rr)
      {
        for (unsigned i = 0; i < deg + 1; ++i)
          elerhs[i](fvi + rr) += fac * shpfunct(vi) * value[i](rr);
      }
    }
  }  // end gaussloop

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                   vuong 08/14|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::UTILS::DbcNurbs::FillMatrixAndRHSForLSDirichletDomain(Teuchos::RCP<DRT::Element> actele,
    const std::vector<Epetra_SerialDenseVector>* knots, const std::vector<int>& lm,
    const std::vector<int>* funct, const std::vector<double>* val, const unsigned deg,
    const double time, Epetra_SerialDenseMatrix& elemass,
    std::vector<Epetra_SerialDenseVector>& elerhs) const
{
  if (deg + 1 != elerhs.size())
    dserror("given degree of time derivative does not match number or rhs vectors!");

  static const int dim = DRT::UTILS::DisTypeToDim<distype>::dim;
  const int ndbcdofs = (int)lm.size();

  // set element data
  static const int nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  // dofblocks (number of DOFs with Dirichlet condition per node)
  const int dofblock = ndbcdofs / nen;

  // get node coordinates of element
  LINALG::Matrix<dim, nen> xyze;
  DRT::Node** nodes = actele->Nodes();

  for (int inode = 0; inode < nen; inode++)
  {
    const double* x = nodes[inode]->X();
    for (int idim = 0; idim < dim; ++idim)
    {
      xyze(idim, inode) = x[idim];
    }
  }

  // aquire weights from nodes
  Epetra_SerialDenseVector weights(nen);

  for (int inode = 0; inode < nen; ++inode)
  {
    DRT::NURBS::ControlPoint* cp = dynamic_cast<DRT::NURBS::ControlPoint*>(nodes[inode]);

    weights(inode) = cp->W();
  }

  // shape functions
  LINALG::Matrix<nen, 1> shpfunct;
  // coordinates of integration points in parameter space
  LINALG::Matrix<dim, 1> xsi;
  // transposed jacobian "dx/ds"
  LINALG::Matrix<dim, dim> xjm;
  // first derivative of shape functions
  LINALG::Matrix<dim, nen> deriv;
  // coordinates of integration point in physical space
  Epetra_SerialDenseVector position(
      3);  // always three-dimensional coordinates for function evaluation!
  // auxiliary date container for dirichlet evaluation
  std::vector<Epetra_SerialDenseVector> value(deg + 1, dofblock);

  // gaussian points
  const DRT::UTILS::IntPointsAndWeights<dim> intpoints(
      DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    for (int idim = 0; idim < dim; idim++)
    {
      xsi(idim) = intpoints.IP().Point(iquad)[idim];
    }

    // evaluate shape function and derivatevs at integration point
    DRT::NURBS::UTILS::nurbs_get_funct_deriv(shpfunct, deriv, xsi, *knots, weights, distype);

    xjm.MultiplyNT(deriv, xyze);
    double det = xjm.Determinant();

    if (det < 1E-16)
      dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", actele->Id(), det);

    // compute integration factor
    double fac = intpoints.IP().qwgt[iquad] * det;

    // get real physical coordinates of integration point
    /*
    //              +-----
    //               \
    //    pos (x) =   +      N (x) * x
    //               /        j       j
    //              +-----
    //              node j
    */
    for (int rr = 0; rr < dim; ++rr)
    {
      position(rr) = shpfunct(0) * xyze(rr, 0);
      for (int mm = 1; mm < nen; ++mm)
      {
        position(rr) += shpfunct(mm) * xyze(rr, mm);
      }
    }
    // if dim < 3, ensure we define a valid z-coordinate!
    for (int rr = dim; rr < 3; ++rr) position(rr) = 0.0;

    for (int rr = 0; rr < dofblock; ++rr)
    {
      // factor given by FUNCTS
      std::vector<double> functimederivfac(deg + 1, 1.0);
      double functfac = 1.0;

      const int funct_num = (*funct)[rr];
      if (funct_num > 0)
      {
        // important: position has to have always three components!!
        functimederivfac = DRT::Problem::Instance()
                               ->Funct((*funct)[rr] - 1)
                               .EvaluateTimeDerivative(rr, position.Values(), time, deg);

        functfac =
            DRT::Problem::Instance()->Funct((*funct)[rr] - 1).Evaluate(rr, position.Values(), time);
      }

      // apply factors to Dirichlet value
      for (unsigned i = 0; i < deg + 1; ++i)
      {
        value[i](rr) = (*val)[rr] * functfac * functimederivfac[i];
      }
    }

    for (int vi = 0; vi < nen; ++vi)  // loop rows  (test functions)
    {
      const int fvi = dofblock * vi;

      for (int ui = 0; ui < nen; ++ui)  // loop columns  (test functions)
      {
        const int fui = dofblock * ui;

        const double diag = fac * shpfunct(ui) * shpfunct(vi);

        for (int rr = 0; rr < dofblock; ++rr)
        {
          elemass(fvi + rr, fui + rr) += diag;
        }
      }
      for (int rr = 0; rr < dofblock; ++rr)
      {
        for (unsigned i = 0; i < deg + 1; ++i)
          elerhs[i](fvi + rr) += fac * shpfunct(vi) * value[i](rr);
      }
    }
  }  // end gaussloop

  return;
}
