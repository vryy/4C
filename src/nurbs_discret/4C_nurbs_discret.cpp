/*----------------------------------------------------------------------*/
/*! \file

\brief a class to manage one nurbs discretization

\level 1


*/
/*----------------------------------------------------------------------*/

#include "4C_nurbs_discret.hpp"

#include "4C_fem_general_utils_boundary_integration.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_fem_general_utils_nurbs_shapefunctions.hpp"
#include "4C_global_data.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_nurbs_discret_nurbs_utils.hpp"
#include "4C_utils_function.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_Time.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 05/08|
 *----------------------------------------------------------------------*/
Discret::Nurbs::NurbsDiscretization::NurbsDiscretization(
    const std::string name, Teuchos::RCP<Epetra_Comm> comm, const unsigned int n_dim)
    : Discret::Discretization::Discretization(name, comm, n_dim), knots_(Teuchos::null)
{
  return;
}


/*----------------------------------------------------------------------*
 |  add a knotvector to the discretization (public)          gammi 05/08|
 *----------------------------------------------------------------------*/
void Discret::Nurbs::NurbsDiscretization::SetKnotVector(
    Teuchos::RCP<Discret::Nurbs::Knotvector> knots)
{
  if (knots == Teuchos::null)
  {
    FOUR_C_THROW(
        "You're trying to set an invalid knotvector in the Discret::Nurbs::NurbsDiscretization "
        "'%s'. The given know vector is a null vector and can't be set as such.",
        (this->Name()).c_str());
  }

  knots_ = knots;
  filled_ = false;
}

/*----------------------------------------------------------------------*
 |  get a pointer to knotvector from the discretization         (public)|
 |                                                           gammi 05/08|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Discret::Nurbs::Knotvector> Discret::Nurbs::NurbsDiscretization::GetKnotVector()
{
  if (knots_ == Teuchos::null)
  {
    FOUR_C_THROW(
        "You're trying to access the NURBS knot vector in the Discret::Nurbs::NurbsDiscretization "
        "'%s'. The required knot vector is a null vector and can't be accessed as such.",
        (this->Name()).c_str());
  }
  return knots_;
}


/*----------------------------------------------------------------------*
 |  get a pointer to knotvector from the discretization         (public)|
 |  (const version, read-only)                               gammi 05/08|
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Discret::Nurbs::Knotvector> Discret::Nurbs::NurbsDiscretization::GetKnotVector()
    const
{
  if (knots_ == Teuchos::null)
  {
    FOUR_C_THROW(
        "You're trying to access the NURBS knot vector in the Discret::Nurbs::NurbsDiscretization "
        "'%s'. The required knot vector is a null vector and can't be accessed as such.",
        (this->Name()).c_str());
  }
  return knots_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Discret::UTILS::DbcNurbs::evaluate(const Teuchos::ParameterList& params,
    const Discret::Discretization& discret, double time,
    const Teuchos::RCP<Epetra_Vector>* systemvectors, Discret::UTILS::Dbc::DbcInfo& info,
    Teuchos::RCP<std::set<int>>* dbcgids) const
{
  // --------------------------- Step 1 ---------------------------------------
  Discret::UTILS::Dbc::evaluate(params, discret, time, systemvectors, info, dbcgids);

  // --------------------------- Step 2 ---------------------------------------
  std::vector<std::string> dbc_cond_names(2, "");
  dbc_cond_names[0] = "Dirichlet";
  dbc_cond_names[1] = "NurbsLSDirichlet";

  std::vector<Teuchos::RCP<Core::Conditions::Condition>> conds(0);
  std::vector<Teuchos::RCP<Core::Conditions::Condition>> curr_conds(0);
  for (std::vector<std::string>::const_iterator cit_name = dbc_cond_names.begin();
       cit_name != dbc_cond_names.end(); ++cit_name)
  {
    discret.GetCondition(*cit_name, curr_conds);

    conds.reserve(conds.size() + curr_conds.size());
    std::copy(curr_conds.begin(), curr_conds.end(), std::back_inserter(conds));
  }

  Discret::UTILS::Dbc::DbcInfo info2(info.toggle.Map());
  read_dirichlet_condition(params, discret, conds, time, info2, dbcgids);

  // --------------------------- Step 3 ---------------------------------------
  conds.clear();
  discret.GetCondition("NurbsLSDirichlet", conds);

  Teuchos::RCP<std::set<int>> dbcgids_nurbs[2] = {Teuchos::null, Teuchos::null};
  dbcgids_nurbs[set_row] = Teuchos::rcp<std::set<int>>(new std::set<int>());
  dbcgids_nurbs[set_col] = Teuchos::rcp<std::set<int>>(new std::set<int>());

  // create a new toggle vector with column layout
  const Discret::Nurbs::NurbsDiscretization* discret_nurbs =
      dynamic_cast<const Discret::Nurbs::NurbsDiscretization*>(&discret);
  if (not discret_nurbs) FOUR_C_THROW("Dynamic cast failed!");

  // build dummy column toggle vector and auxiliary vectors
  Discret::UTILS::Dbc::DbcInfo info_col(*discret_nurbs->DofColMap());
  read_dirichlet_condition(params, discret, conds, time, info_col, dbcgids_nurbs);

  // --------------------------- Step 4 ---------------------------------------
  do_dirichlet_condition(
      params, discret, conds, time, systemvectors, info_col.toggle, dbcgids_nurbs);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Discret::UTILS::DbcNurbs::do_dirichlet_condition(const Teuchos::ParameterList& params,
    const Discret::Discretization& discret, const Core::Conditions::Condition& cond, double time,
    const Teuchos::RCP<Epetra_Vector>* systemvectors, const Epetra_IntVector& toggle,
    const Teuchos::RCP<std::set<int>>* dbcgids) const
{
  // default call
  if (dbcgids[set_col].is_null())
  {
    Discret::UTILS::Dbc::do_dirichlet_condition(
        params, discret, cond, time, systemvectors, toggle, dbcgids);
    return;
  }

  Teuchos::Time timer("", true);

  // get the processor ID from the communicator
  const int myrank = discret.Comm().MyPID();
  if (myrank == 0) std::cout << "calculating least squares Dirichlet condition in ... ";

  const Discret::Nurbs::NurbsDiscretization& nurbs_dis =
      static_cast<const Discret::Nurbs::NurbsDiscretization&>(discret);

  // create map extractor to always (re)build dbcmapextractor which is needed later
  Teuchos::RCP<Core::LinAlg::MapExtractor> auxdbcmapextractor =
      Teuchos::rcp(new Core::LinAlg::MapExtractor());
  {
    // build map of Dirichlet DOFs
    int nummyelements = 0;
    int* myglobalelements = nullptr;
    std::vector<int> dbcgidsv;
    if (dbcgids[set_row]->size() > 0)
    {
      dbcgidsv.reserve(dbcgids[set_row]->size());
      dbcgidsv.assign(dbcgids[set_row]->begin(), dbcgids[set_row]->end());
      nummyelements = dbcgidsv.size();
      myglobalelements = dbcgidsv.data();
    }
    Teuchos::RCP<Epetra_Map> dbcmap = Teuchos::rcp(new Epetra_Map(-1, nummyelements,
        myglobalelements, discret.dof_row_map()->IndexBase(), discret.dof_row_map()->Comm()));
    // build the map extractor of Dirichlet-conditioned and free DOFs
    *auxdbcmapextractor = Core::LinAlg::MapExtractor(*(discret.dof_row_map()), dbcmap);
  }

  // column map of all DOFs subjected to a least squares Dirichlet condition
  Teuchos::RCP<Epetra_Map> dbccolmap = Teuchos::null;
  {
    // build map of Dirichlet DOFs
    int nummyelements = 0;
    int* myglobalelements = nullptr;
    std::vector<int> dbcgidsv;
    if (dbcgids[set_col]->size() > 0)
    {
      dbcgidsv.reserve(dbcgids[set_col]->size());
      dbcgidsv.assign(dbcgids[set_col]->begin(), dbcgids[set_col]->end());
      nummyelements = dbcgidsv.size();
      myglobalelements = dbcgidsv.data();
    }
    dbccolmap = Teuchos::rcp(new Epetra_Map(-1, nummyelements, myglobalelements,
        nurbs_dis.DofColMap()->IndexBase(), discret.dof_row_map()->Comm()));
  }

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Teuchos::RCP<const Epetra_Map> dofrowmap = auxdbcmapextractor->CondMap();

  if (dofrowmap->NumGlobalElements() == 0) return;  // no dbc gids ->leave

  // read information from condition
  const std::vector<int>* nodeids = cond.GetNodes();
  if (!nodeids) FOUR_C_THROW("Dirichlet condition does not have nodal cloud");

  const auto* funct = &cond.parameters().Get<std::vector<int>>("funct");
  const auto* val = &cond.parameters().Get<std::vector<double>>("val");


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
  FOUR_C_ASSERT(systemvectoraux != Teuchos::null, "At least one vector must be unequal to null");


  // -------------------------------------------------------------------
  // create empty mass matrix
  // -------------------------------------------------------------------
  Teuchos::RCP<Core::LinAlg::SparseMatrix> massmatrix =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dofrowmap, 108, false, true));

  // -------------------------------------------------------------------
  // create empty right hand side vector
  // -------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> rhs = Core::LinAlg::CreateVector(*dofrowmap, true);
  Teuchos::RCP<Epetra_Vector> dbcvector = Core::LinAlg::CreateVector(*dofrowmap, true);

  Teuchos::RCP<Epetra_Vector> rhsd = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> dbcvectord = Teuchos::null;
  if (systemvectors[1] != Teuchos::null)
  {
    rhsd = Core::LinAlg::CreateVector(*dofrowmap, true);
    dbcvectord = Core::LinAlg::CreateVector(*dofrowmap, true);
  }

  Teuchos::RCP<Epetra_Vector> rhsdd = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> dbcvectordd = Teuchos::null;
  if (systemvectors[2] != Teuchos::null)
  {
    rhsdd = Core::LinAlg::CreateVector(*dofrowmap, true);
    dbcvectordd = Core::LinAlg::CreateVector(*dofrowmap, true);
  }

  const bool assemblevecd = rhsd != Teuchos::null;
  const bool assemblevecdd = rhsdd != Teuchos::null;

  // -------------------------------------------------------------------
  // call elements to calculate massmatrix and righthandside
  // -------------------------------------------------------------------
  {
    // call elements and assemble
    if (!discret.Filled()) FOUR_C_THROW("fill_complete() was not called");
    if (!discret.HaveDofs()) FOUR_C_THROW("assign_degrees_of_freedom() was not called");

    // see what we have for input
    bool assemblemat = massmatrix != Teuchos::null;
    bool assemblevec = rhs != Teuchos::null;

    // define element matrices and vectors
    Core::LinAlg::SerialDenseMatrix elemass;
    std::vector<Core::LinAlg::SerialDenseVector> elerhs(deg + 1);

    std::vector<int> lm;
    std::vector<int> lmowner;

    std::vector<int> lm_full;
    std::vector<int> lmowner_full;
    std::vector<int> lmstride_full;

    const std::map<int, Teuchos::RCP<Core::Elements::Element>>& geom = cond.Geometry();
    std::map<int, Teuchos::RCP<Core::Elements::Element>>::const_iterator curr;
    for (curr = geom.begin(); curr != geom.end(); ++curr)
    {
      Teuchos::RCP<Core::Elements::Element> actele = curr->second;

      static const int probdim = Global::Problem::Instance()->NDim();
      const Core::FE::CellType distype = actele->Shape();
      const int dim = Core::FE::getDimension(distype);
      const bool isboundary = (dim != probdim);
      const int nen = Core::FE::getNumberOfElementNodes(distype);

      // access elements knot span
      std::vector<Core::LinAlg::SerialDenseVector> eleknots(dim);
      Core::LinAlg::SerialDenseVector weights(nen);

      bool zero_size = false;
      if (isboundary)
      {
        Teuchos::RCP<Core::Elements::FaceElement> faceele =
            Teuchos::rcp_dynamic_cast<Core::Elements::FaceElement>(actele, true);
        double normalfac = 0.0;
        std::vector<Core::LinAlg::SerialDenseVector> pknots(probdim);
        zero_size = Discret::Nurbs::GetKnotVectorAndWeightsForNurbsBoundary(actele.get(),
            faceele->FaceMasterNumber(), faceele->parent_element()->Id(), discret, pknots, eleknots,
            weights, normalfac);
      }
      else
        zero_size =
            Discret::Nurbs::GetMyNurbsKnotsAndWeights(discret, actele.get(), eleknots, weights);

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
      for (unsigned j = 0; j < lm_full.size(); ++j)
      {
        int gid = lm_full[j];
        if (dbccolmap->MyGID(gid))
        {
          lm.push_back(gid);
          lmowner.push_back(lmowner_full[j]);
        }
      }

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      const int eledim = (int)lm.size();

      if (assemblemat)
      {
        if (elemass.numRows() != eledim or elemass.numCols() != eledim)
          elemass.shape(eledim, eledim);
        else
          elemass.putScalar(0.0);
      }
      if (assemblevec)
      {
        for (unsigned i = 0; i < deg + 1; ++i)
        {
          if (elerhs[i].length() != eledim)
            elerhs[i].size(eledim);
          else
            elerhs[i].putScalar(0.0);
        }
      }

      if (isboundary) switch (distype)
        {
          case Core::FE::CellType::nurbs2:
            fill_matrix_and_rhs_for_ls_dirichlet_boundary<Core::FE::CellType::nurbs2>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case Core::FE::CellType::nurbs3:
            fill_matrix_and_rhs_for_ls_dirichlet_boundary<Core::FE::CellType::nurbs3>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case Core::FE::CellType::nurbs4:
            fill_matrix_and_rhs_for_ls_dirichlet_boundary<Core::FE::CellType::nurbs4>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case Core::FE::CellType::nurbs9:
            fill_matrix_and_rhs_for_ls_dirichlet_boundary<Core::FE::CellType::nurbs9>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case Core::FE::CellType::nurbs8:
            fill_matrix_and_rhs_for_ls_dirichlet_boundary<Core::FE::CellType::nurbs8>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case Core::FE::CellType::nurbs27:
            fill_matrix_and_rhs_for_ls_dirichlet_boundary<Core::FE::CellType::nurbs27>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          default:
            FOUR_C_THROW("invalid element shape for least squares dirichlet evaluation: %s",
                Core::FE::CellTypeToString(distype).c_str());
            break;
        }
      else
        switch (distype)
        {
          case Core::FE::CellType::nurbs2:
            fill_matrix_and_rhs_for_ls_dirichlet_domain<Core::FE::CellType::nurbs2>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case Core::FE::CellType::nurbs3:
            fill_matrix_and_rhs_for_ls_dirichlet_domain<Core::FE::CellType::nurbs3>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case Core::FE::CellType::nurbs4:
            fill_matrix_and_rhs_for_ls_dirichlet_domain<Core::FE::CellType::nurbs4>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case Core::FE::CellType::nurbs9:
            fill_matrix_and_rhs_for_ls_dirichlet_domain<Core::FE::CellType::nurbs9>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case Core::FE::CellType::nurbs8:
            fill_matrix_and_rhs_for_ls_dirichlet_domain<Core::FE::CellType::nurbs8>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          case Core::FE::CellType::nurbs27:
            fill_matrix_and_rhs_for_ls_dirichlet_domain<Core::FE::CellType::nurbs27>(
                actele, &eleknots, lm, funct, val, deg, time, elemass, elerhs);
            break;
          default:
            FOUR_C_THROW("invalid element shape for least squares dirichlet evaluation: %s",
                Core::FE::CellTypeToString(distype).c_str());
            break;
        }

      int eid = actele->Id();
      if (assemblemat) massmatrix->Assemble(eid, elemass, lm, lmowner);
      if (assemblevec) Core::LinAlg::Assemble(*rhs, elerhs[0], lm, lmowner);
      if (assemblevecd) Core::LinAlg::Assemble(*rhsd, elerhs[1], lm, lmowner);
      if (assemblevecdd) Core::LinAlg::Assemble(*rhsdd, elerhs[2], lm, lmowner);
    }
  }
  // -------------------------------------------------------------------
  // finalize the system matrix
  // -------------------------------------------------------------------
  massmatrix->Complete();

  // -------------------------------------------------------------------
  // solve system
  // -------------------------------------------------------------------
  // Owing to experience a very accurate solution has to be enforced here!
  // Thus, we allocate an own solver with VERY strict tolerance!
  // One could think of verifiying an extra solver in the input file...
  Teuchos::ParameterList p = Global::Problem::Instance()->UMFPACKSolverParams();
  //  const double origtol = p.get<double>("AZTOL");
  //  const double newtol  = 1.0e-11;
  //  p.set("AZTOL",newtol);

  //  if(myrank==0)
  //    cout<<"\nSolver tolerance for least squares problem set to "<<newtol<<"\n";

  Teuchos::RCP<Core::LinAlg::Solver> solver =
      Teuchos::rcp(new Core::LinAlg::Solver(p, discret.Comm()));
  // FixMe actually the const qualifier could stay, if someone adds to each single
  // related ComputeNullSpace routine a "const"....
  const_cast<Discret::Discretization&>(discret).compute_null_space_if_necessary(solver->Params());

  // solve for control point values
  // always refactor and reset the matrix before a single new solver call
  Core::LinAlg::SolverParams solver_params;
  solver_params.refactor = true;
  solver_params.reset = true;
  solver->Solve(massmatrix->EpetraOperator(), dbcvector, rhs, solver_params);

  // solve for first derivatives in time
  if (assemblevecd) solver->Solve(massmatrix->EpetraOperator(), dbcvectord, rhsd, solver_params);

  // solve for second derivatives in time
  if (assemblevecdd) solver->Solve(massmatrix->EpetraOperator(), dbcvectordd, rhsdd, solver_params);

  // perform resets for solver and matrix
  solver->Reset();
  massmatrix->Reset();

  // insert nodal values to sysvec
  auxdbcmapextractor->InsertCondVector(dbcvector, systemvectors[0]);
  if (assemblevecd) auxdbcmapextractor->InsertCondVector(dbcvectord, systemvectors[1]);
  if (assemblevecdd) auxdbcmapextractor->InsertCondVector(dbcvectordd, systemvectors[2]);

  if (myrank == 0) std::cout << timer.totalElapsedTime(true) << " seconds \n\n";

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                   vuong 08/14|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::UTILS::DbcNurbs::fill_matrix_and_rhs_for_ls_dirichlet_boundary(
    Teuchos::RCP<Core::Elements::Element> actele,
    const std::vector<Core::LinAlg::SerialDenseVector>* knots, const std::vector<int>& lm,
    const std::vector<int>* funct, const std::vector<double>* val, const unsigned deg,
    const double time, Core::LinAlg::SerialDenseMatrix& elemass,
    std::vector<Core::LinAlg::SerialDenseVector>& elerhs) const
{
  if (deg + 1 != elerhs.size())
    FOUR_C_THROW("given degree of time derivative does not match number or rhs vectors!");

  static const int dim = Core::FE::dim<distype>;

  const int ndbcdofs = (int)lm.size();

  // set element data
  static const int nen = Core::FE::num_nodes<distype>;

  // dofblocks (number of DOFs with Dirichlet condition per node)
  const int dofblock = ndbcdofs / nen;

  // get node coordinates of element
  Core::LinAlg::Matrix<dim + 1, nen> xyze;
  Core::Nodes::Node** nodes = actele->Nodes();

  for (int inode = 0; inode < nen; inode++)
  {
    const auto& x = nodes[inode]->X();
    for (int idim = 0; idim < dim + 1; ++idim)
    {
      xyze(idim, inode) = x[idim];
    }
  }

  // aquire weights from nodes
  Core::LinAlg::SerialDenseVector weights(nen);

  for (int inode = 0; inode < nen; ++inode)
  {
    Discret::Nurbs::ControlPoint* cp = dynamic_cast<Discret::Nurbs::ControlPoint*>(nodes[inode]);

    weights(inode) = cp->W();
  }

  // shape functions
  Core::LinAlg::Matrix<nen, 1> shpfunct;
  // coordinates of integration points in parameter space
  Core::LinAlg::Matrix<dim, 1> xsi;
  // first derivative of shape functions
  Core::LinAlg::Matrix<dim, nen> deriv;
  // coordinates of integration point in physical space
  Core::LinAlg::SerialDenseVector position(
      3);  // always three-dimensional coordinates for function evaluation!
  // auxiliary date container for dirichlet evaluation
  std::vector<Core::LinAlg::SerialDenseVector> value(deg + 1, dofblock);
  // unit normal on boundary element
  Core::LinAlg::Matrix<dim + 1, 1> unitnormal;

  // gaussian points
  const Core::FE::IntPointsAndWeights<dim> intpoints(
      Discret::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // integration factor
    double fac = 0.0;
    double drs = 0.0;

    Core::FE::EvalShapeFuncAtBouIntPoint<distype>(
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
        functimederivfac = Global::Problem::Instance()
                               ->FunctionById<Core::UTILS::FunctionOfSpaceTime>((*funct)[rr] - 1)
                               .evaluate_time_derivative(position.values(), time, deg, rr);
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
template <Core::FE::CellType distype>
void Discret::UTILS::DbcNurbs::fill_matrix_and_rhs_for_ls_dirichlet_domain(
    Teuchos::RCP<Core::Elements::Element> actele,
    const std::vector<Core::LinAlg::SerialDenseVector>* knots, const std::vector<int>& lm,
    const std::vector<int>* funct, const std::vector<double>* val, const unsigned deg,
    const double time, Core::LinAlg::SerialDenseMatrix& elemass,
    std::vector<Core::LinAlg::SerialDenseVector>& elerhs) const
{
  if (deg + 1 != elerhs.size())
    FOUR_C_THROW("given degree of time derivative does not match number or rhs vectors!");

  static const int dim = Core::FE::dim<distype>;
  const int ndbcdofs = (int)lm.size();

  // set element data
  static const int nen = Core::FE::num_nodes<distype>;

  // dofblocks (number of DOFs with Dirichlet condition per node)
  const int dofblock = ndbcdofs / nen;

  // get node coordinates of element
  Core::LinAlg::Matrix<dim, nen> xyze;
  Core::Nodes::Node** nodes = actele->Nodes();

  for (int inode = 0; inode < nen; inode++)
  {
    const auto& x = nodes[inode]->X();
    for (int idim = 0; idim < dim; ++idim)
    {
      xyze(idim, inode) = x[idim];
    }
  }

  // aquire weights from nodes
  Core::LinAlg::SerialDenseVector weights(nen);

  for (int inode = 0; inode < nen; ++inode)
  {
    Discret::Nurbs::ControlPoint* cp = dynamic_cast<Discret::Nurbs::ControlPoint*>(nodes[inode]);

    weights(inode) = cp->W();
  }

  // shape functions
  Core::LinAlg::Matrix<nen, 1> shpfunct;
  // coordinates of integration points in parameter space
  Core::LinAlg::Matrix<dim, 1> xsi;
  // transposed jacobian "dx/ds"
  Core::LinAlg::Matrix<dim, dim> xjm;
  // first derivative of shape functions
  Core::LinAlg::Matrix<dim, nen> deriv;
  // coordinates of integration point in physical space
  Core::LinAlg::SerialDenseVector position(
      3);  // always three-dimensional coordinates for function evaluation!
  // auxiliary date container for dirichlet evaluation
  std::vector<Core::LinAlg::SerialDenseVector> value(deg + 1, dofblock);

  // gaussian points
  const Core::FE::IntPointsAndWeights<dim> intpoints(
      Discret::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    for (int idim = 0; idim < dim; idim++)
    {
      xsi(idim) = intpoints.IP().Point(iquad)[idim];
    }

    // evaluate shape function and derivatevs at integration point
    Core::FE::Nurbs::nurbs_get_funct_deriv(shpfunct, deriv, xsi, *knots, weights, distype);

    xjm.MultiplyNT(deriv, xyze);
    double det = xjm.Determinant();

    if (det < 1E-16)
      FOUR_C_THROW(
          "GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", actele->Id(), det);

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
        functimederivfac = Global::Problem::Instance()
                               ->FunctionById<Core::UTILS::FunctionOfSpaceTime>((*funct)[rr] - 1)
                               .evaluate_time_derivative(position.values(), time, deg, rr);

        functfac = Global::Problem::Instance()
                       ->FunctionById<Core::UTILS::FunctionOfSpaceTime>((*funct)[rr] - 1)
                       .Evaluate(position.values(), time, rr);
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

FOUR_C_NAMESPACE_CLOSE
