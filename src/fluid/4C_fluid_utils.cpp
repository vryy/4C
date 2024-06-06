/*-----------------------------------------------------------*/
/*! \file

\brief utility functions for fluid problems


\level 2

*/
/*-----------------------------------------------------------*/

#include "4C_fluid_utils.hpp"

#include "4C_discretization_dofset.hpp"
#include "4C_discretization_dofset_interface.hpp"
#include "4C_discretization_fem_general_l2_projection.hpp"
#include "4C_discretization_fem_general_utils_superconvergent_patch_recovery.hpp"
#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_implicit_integration.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_linear_solver_method_linalg.hpp"

#include <MLAPI_Aggregation.h>
#include <MLAPI_Workspace.h>
#include <stdio.h>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | constructor                                         Thon/Krank 11/14 |
 *----------------------------------------------------------------------*/
FLD::UTILS::StressManager::StressManager(Teuchos::RCP<Discret::Discretization> discret,
    Teuchos::RCP<Epetra_Vector> dispnp, const bool alefluid, const int numdim)
    : discret_(discret),
      dispnp_(dispnp),
      alefluid_(alefluid),
      numdim_(numdim),
      sep_enr_(Teuchos::null),
      wss_type_(Core::UTILS::IntegralValue<Inpar::FLUID::WSSType>(
          Global::Problem::Instance()->FluidDynamicParams(), "WSS_TYPE")),
      sum_stresses_(Teuchos::null),
      sum_wss_(Teuchos::null),
      sum_dt_stresses_(0.0),
      sum_dt_wss_(0.0),
      isinit_(false)
{
  switch (wss_type_)
  {
    case Inpar::FLUID::wss_standard:
      isinit_ = true;  // in this cases nothing has to be initialized
      break;
    case Inpar::FLUID::wss_aggregation:
      isinit_ = false;  // we do this in InitAggr()
      break;
    case Inpar::FLUID::wss_mean:
      sum_stresses_ = Teuchos::rcp(new Epetra_Vector(*(discret_->dof_row_map()), true)),
      sum_wss_ = Teuchos::rcp(new Epetra_Vector(*(discret_->dof_row_map()), true)), isinit_ = true;
      break;
    default:
      FOUR_C_THROW(
          "There are only the wss calculation types 'standard', 'aggregation' and 'mean'!!");
      break;
  }
}

/*----------------------------------------------------------------------*
 | constructor                                         Thon/Krank 11/14 |
 *----------------------------------------------------------------------*/
void FLD::UTILS::StressManager::InitAggr(Teuchos::RCP<Core::LinAlg::SparseOperator> sysmat)
{
  if (wss_type_ != Inpar::FLUID::wss_aggregation)
    FOUR_C_THROW("One should end up here just in case of aggregated stresses!");

  calc_sep_enr(sysmat);
  if (sep_enr_ == Teuchos::null)
    FOUR_C_THROW("SepEnr matrix has not been build correctly. Strange...");

  isinit_ = true;

  return;
}

/*----------------------------------------------------------------------*
 | update and return WSS vector                        Thon/Krank 11/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::UTILS::StressManager::get_wall_shear_stresses(
    Teuchos::RCP<const Epetra_Vector> trueresidual, const double dt)
{
  Teuchos::RCP<Epetra_Vector> wss = get_wall_shear_stresses_wo_agg(trueresidual);

  switch (wss_type_)
  {
    case Inpar::FLUID::wss_standard:
      // nothing to do
      break;
    case Inpar::FLUID::wss_aggregation:
      wss = aggreagte_stresses(wss);
      break;
    case Inpar::FLUID::wss_mean:
      wss = time_average_wss(wss, dt);
      break;
    default:
      FOUR_C_THROW(
          "There are only the wss calculation types 'standard', 'aggregation' and 'mean'!!");
      break;
  }

  return wss;
}

/*--------------------------------------------------------------------------*
 | return wss vector (without updating the mean stress vector)   Thon 11/14 |
 *--------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::UTILS::StressManager::get_pre_calc_wall_shear_stresses(
    Teuchos::RCP<const Epetra_Vector> trueresidual)
{
  Teuchos::RCP<Epetra_Vector> wss =
      Teuchos::rcp(new Epetra_Vector(*(discret_->dof_row_map()), true));

  switch (wss_type_)
  {
    case Inpar::FLUID::wss_standard:
      wss = get_wall_shear_stresses_wo_agg(trueresidual);
      break;
    case Inpar::FLUID::wss_aggregation:
      wss = get_wall_shear_stresses_wo_agg(trueresidual);
      wss = aggreagte_stresses(wss);
      break;
    case Inpar::FLUID::wss_mean:
      if (sum_dt_wss_ > 0.0)  // iff we have actually calculated some mean wss
        wss->Update(1.0 / sum_dt_wss_, *sum_wss_, 0.0);  // weighted sum of all prior stresses
      break;
    default:
      FOUR_C_THROW(
          "There are only the wss calculation types 'standard', 'aggregation' and 'mean'!!");
      break;
  }

  return wss;
}

/*----------------------------------------------------------------------*
 | return WSS vector always without aggregation        Thon/Krank 11/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::UTILS::StressManager::get_wall_shear_stresses_wo_agg(
    Teuchos::RCP<const Epetra_Vector> trueresidual)
{
  if (not isinit_) FOUR_C_THROW("StressManager not initialized");

  Teuchos::RCP<Epetra_Vector> stresses = calc_stresses(trueresidual);
  // calculate wss from stresses
  Teuchos::RCP<Epetra_Vector> wss = calc_wall_shear_stresses(stresses);

  return wss;
}

/*-----------------------------------------------------------------------------*
 |  calculate traction vector at (Dirichlet) boundary (public) Thon/Krank 07/07|
 *-----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::UTILS::StressManager::GetStressesWOAgg(
    Teuchos::RCP<const Epetra_Vector> trueresidual)
{
  Teuchos::RCP<Epetra_Vector> stresses = calc_stresses(trueresidual);

  return stresses;
}  // FLD::UTILS::StressManager::GetStressesWOAgg()

/*-----------------------------------------------------------------------------*
 |  update and return stress vector                            Thon/Krank 07/07|
 *-----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::UTILS::StressManager::GetStresses(
    Teuchos::RCP<const Epetra_Vector> trueresidual, const double dt)
{
  Teuchos::RCP<Epetra_Vector> stresses = GetStressesWOAgg(trueresidual);

  switch (wss_type_)
  {
    case Inpar::FLUID::wss_standard:
      // nothing to do
      break;
    case Inpar::FLUID::wss_aggregation:
      stresses = aggreagte_stresses(stresses);
      break;
    case Inpar::FLUID::wss_mean:
      stresses = time_average_stresses(stresses, dt);
      break;
    default:
      FOUR_C_THROW(
          "There are only the wss calculation types 'standard', 'aggregation' and 'mean'!!");
      break;
  }

  return stresses;
}  // FLD::UTILS::StressManager::GetStresses()

/*-----------------------------------------------------------------------------*
 | return stress vector (without updating the mean stress vector)   Thon 03/15 |
 *-----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::UTILS::StressManager::GetPreCalcStresses(
    Teuchos::RCP<const Epetra_Vector> trueresidual)
{
  Teuchos::RCP<Epetra_Vector> stresses =
      Teuchos::rcp(new Epetra_Vector(*(discret_->dof_row_map()), true));

  switch (wss_type_)
  {
    case Inpar::FLUID::wss_standard:
      stresses = GetStressesWOAgg(trueresidual);
      break;
    case Inpar::FLUID::wss_aggregation:
      stresses = GetStressesWOAgg(trueresidual);
      stresses = aggreagte_stresses(stresses);
      break;
    case Inpar::FLUID::wss_mean:
      if (sum_dt_stresses_ > 0.0)  // iff we have actually calculated some mean stresses
        stresses->Update(
            1.0 / sum_dt_stresses_, *sum_stresses_, 0.0);  // weighted sum of all prior stresses
      break;
    default:
      FOUR_C_THROW(
          "There are only the wss calculation types 'standard', 'aggregation' and 'mean'!!");
      break;
  }
  return stresses;
}

/*-----------------------------------------------------------------------------*
 |  calculate traction vector at (Dirichlet) boundary (public) Thon/Krank 07/07|
 *-----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::UTILS::StressManager::calc_stresses(
    Teuchos::RCP<const Epetra_Vector> trueresidual)
{
  if (not isinit_) FOUR_C_THROW("StressManager not initialized");
  std::string condstring("FluidStressCalc");
  Teuchos::RCP<Epetra_Vector> integratedshapefunc = integrate_interface_shape(condstring);

  // compute traction values at specified nodes; otherwise do not touch the zero values
  for (int i = 0; i < integratedshapefunc->MyLength(); i++)
  {
    if ((*integratedshapefunc)[i] != 0.0)
    {
      // overwrite integratedshapefunc values with the calculated traction coefficients,
      // which are reconstructed out of the nodal forces (trueresidual_) using the
      // same shape functions on the boundary as for velocity and pressure.
      (*integratedshapefunc)[i] = (*trueresidual)[i] / (*integratedshapefunc)[i];
    }
  }

  return integratedshapefunc;
}  // FLD::UTILS::StressManager::calc_stresses()

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::UTILS::StressManager::integrate_interface_shape(
    std::string condname)
{
  Teuchos::ParameterList eleparams;
  // set action for elements
  eleparams.set<int>("action", FLD::integrate_Shapefunction);

  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  const Epetra_Map* dofrowmap = discret_->dof_row_map();

  // create vector (+ initialization with zeros)
  Teuchos::RCP<Epetra_Vector> integratedshapefunc = Core::LinAlg::CreateVector(*dofrowmap, true);

  // call loop over elements
  discret_->ClearState();
  if (alefluid_)
  {
    discret_->set_state("dispnp", dispnp_);
  }
  discret_->evaluate_condition(eleparams, integratedshapefunc, condname);
  discret_->ClearState();

  return integratedshapefunc;
}

/*----------------------------------------------------------------------*
 |  calculate wall sheer stress from stresses at dirichlet boundary     |
 |                                                     Thon/Krank 11/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::UTILS::StressManager::calc_wall_shear_stresses(
    Teuchos::RCP<Epetra_Vector> stresses)
{
  // -------------------------------------------------------------------
  // first evaluate the normals at the nodes
  // -------------------------------------------------------------------

  Teuchos::ParameterList eleparams;
  // set action for elements
  eleparams.set<int>("action", FLD::ba_calc_node_normal);

  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  const Epetra_Map* dofrowmap = discret_->dof_row_map();

  // vector ndnorm0 with pressure-entries is needed for evaluate_condition
  Teuchos::RCP<Epetra_Vector> ndnorm0 = Core::LinAlg::CreateVector(*dofrowmap, true);

  // call loop over elements, note: normal vectors do not yet have length = 1.0
  discret_->ClearState();  // TODO: (Thon) Do we really have to to this in here?
  if (alefluid_)
  {
    discret_->set_state("dispnp", dispnp_);
  }
  // evaluate the normals of the surface
  discret_->evaluate_condition(eleparams, ndnorm0, "FluidStressCalc");
  discret_->ClearState();

  // -------------------------------------------------------------------
  // normalize the normal vectors
  // -------------------------------------------------------------------
  for (int i = 0; i < ndnorm0->MyLength(); i += numdim_ + 1)
  {
    // calculate the length of the normal
    double L = 0.0;
    for (int j = 0; j < numdim_; j++)
    {
      L += ((*ndnorm0)[i + j]) * ((*ndnorm0)[i + j]);
    }
    L = sqrt(L);

    // normalise the normal vector (if present for the current node)
    if (L > 1e-15)
    {
      for (int j = 0; j < numdim_; j++)
      {
        (*ndnorm0)[i + j] /= L;
      }
    }
  }

  // -------------------------------------------------------------------
  // evaluate the wall shear stress from the traction by removing
  // the normal stresses
  // -------------------------------------------------------------------

  // get traction
  Teuchos::RCP<Epetra_Vector> wss = stresses;

  // loop over all entities within the traction vector
  for (int i = 0; i < ndnorm0->MyLength(); i += numdim_ + 1)
  {
    // evaluate the normal stress = < traction . normal >
    double normal_stress = 0.0;
    for (int j = 0; j < numdim_; j++)
    {
      normal_stress += (*wss)[i + j] * (*ndnorm0)[i + j];
    }

    // subtract the normal stresses from traction
    for (int j = 0; j < numdim_; j++)
    {
      (*wss)[i + j] -= normal_stress * (*ndnorm0)[i + j];
    }
  }

  // -------------------------------------------------------------------
  // return the wall_shear_stress vector
  // -------------------------------------------------------------------
  return wss;
}


/*----------------------------------------------------------------------*
 | smooth stresses/wss via ML-aggregation              Thon/Krank 11/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::UTILS::StressManager::aggreagte_stresses(
    Teuchos::RCP<Epetra_Vector> wss)
{
  if (sep_enr_ == Teuchos::null) FOUR_C_THROW("no scale separation matrix");

  Teuchos::RCP<Epetra_Vector> mean_wss =
      Teuchos::rcp(new Epetra_Vector(*(discret_->dof_row_map()), true));

  // Do the actual aggregation
  sep_enr_->Multiply(false, *wss, *mean_wss);

  return mean_wss;
}

/*----------------------------------------------------------------------*
 | time average stresses                                     Thon 03/15 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::UTILS::StressManager::time_average_stresses(
    Teuchos::RCP<const Epetra_Vector> stresses, const double dt)
{
  sum_stresses_->Update(dt, *stresses, 1.0);  // weighted sum of all prior stresses
  sum_dt_stresses_ += dt;

  Teuchos::RCP<Epetra_Vector> mean_stresses =
      Teuchos::rcp(new Epetra_Vector(*(discret_->dof_row_map()), true));
  mean_stresses->Update(1.0 / sum_dt_stresses_, *sum_stresses_, 0.0);

  return mean_stresses;
}

/*----------------------------------------------------------------------*
 | time average wss                                          Thon 03/15 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FLD::UTILS::StressManager::time_average_wss(
    Teuchos::RCP<const Epetra_Vector> wss, const double dt)
{
  sum_wss_->Update(dt, *wss, 1.0);  // weighted sum of all prior stresses
  sum_dt_wss_ += dt;

  Teuchos::RCP<Epetra_Vector> mean_wss =
      Teuchos::rcp(new Epetra_Vector(*(discret_->dof_row_map()), true));
  mean_wss->Update(1.0 / sum_dt_wss_, *sum_wss_, 0.0);

  return mean_wss;
}

/*----------------------------------------------------------------------*
 | Calculate Aggregation Matrix and set is as member variable SepEnr_   |
 |                                                     Thon/Krank 11/14 |
 *------------------------------------------------- --------------------*/
void FLD::UTILS::StressManager::calc_sep_enr(Teuchos::RCP<Core::LinAlg::SparseOperator> sysmat)
{
  if (wss_type_ == Inpar::FLUID::wss_aggregation)  // iff we have not specified a ML-solver one does
                                                   // not want to smooth the wss
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> sysmat2;

    // Try this:
    sysmat2 = Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(sysmat);
    if (sysmat2 == Teuchos::null)  // if it does not work the fluid matrix probably is a
                                   // BlockSparseMatrix, compare with function use_block_matrix()
      sysmat2 = Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(sysmat)->Merge();
    if (sysmat2 == Teuchos::null)
      FOUR_C_THROW("One of these two dynamic casts should have worked... Sorry!");


    if (discret_->Comm().MyPID() == 0)
      std::cout << "Calculating mean WSS via ML-aggregation:" << std::endl;

    MLAPI::Init();

    int ML_solver =
        (Global::Problem::Instance()->FluidDynamicParams()).get<int>("WSS_ML_AGR_SOLVER");

    if (ML_solver == -1)
      FOUR_C_THROW(
          "If you want to aggregate your stresses you need to specify a WSS_ML_AGR_SOLVER!");

    Teuchos::RCP<Core::LinAlg::Solver> solver = Teuchos::rcp(new Core::LinAlg::Solver(
        Global::Problem::Instance()->SolverParams(ML_solver), discret_->Comm()));

    if (solver == Teuchos::null)
      FOUR_C_THROW(
          "The solver WSS_ML_AGR_SOLVER in the FLUID DYNMAICS section is not a valid solver!");

    Teuchos::ParameterList& mlparams = solver->Params().sublist("ML Parameters");
    // compute the null space,
    discret_->compute_null_space_if_necessary(solver->Params(), true);

    // get nullspace parameters
    double* nullspace = mlparams.get("null space: vectors", (double*)nullptr);
    if (!nullspace) FOUR_C_THROW("No nullspace supplied in parameter list");
    int nsdim = mlparams.get("null space: dimension", 1);
    if (nsdim != 4)
      FOUR_C_THROW("The calculation of mean WSS is only tested for three space dimensions!");

    int lrowdofs = discret_->dof_row_map()->NumMyElements();

    for (int j = 0; j < discret_->NodeRowMap()->NumMyElements(); ++j)
    {
      int gid = discret_->NodeRowMap()->GID(j);

      if (not discret_->NodeRowMap()->MyGID(gid))  // just in case
        FOUR_C_THROW("not on proc");
      {
        Core::Nodes::Node* node = discret_->gNode(gid);
        if (!node) FOUR_C_THROW("Cannot find node");

        int firstglobaldofid = discret_->Dof(node, 0);
        int firstlocaldofid = discret_->dof_row_map()->LID(firstglobaldofid);

        std::vector<Core::Conditions::Condition*> nodedircond;
        node->GetCondition("FluidStressCalc", nodedircond);

        if (not nodedircond.empty())
        {
          // these nodes are wall nodes, so aggregate them
          nullspace[firstlocaldofid] = 1.0;
          nullspace[lrowdofs + firstlocaldofid + 1] = 1.0;
          nullspace[lrowdofs * 2 + firstlocaldofid + 2] = 1.0;
          nullspace[lrowdofs * 3 + firstlocaldofid + 3] = 1.0;
        }
        else
        {
          // set everything to zero
          nullspace[firstlocaldofid] = 0.0;
          nullspace[lrowdofs + firstlocaldofid + 1] = 0.0;
          nullspace[lrowdofs * 2 + firstlocaldofid + 2] = 0.0;
          nullspace[lrowdofs * 3 + firstlocaldofid + 3] = 0.0;
        }

        if (discret_->NumDof(node) > 5)  // in the case of xWall fluid
        {
          nullspace[firstlocaldofid + 4] = 0.0;
          nullspace[lrowdofs + firstlocaldofid + 5] = 0.0;
          nullspace[lrowdofs * 2 + firstlocaldofid + 6] = 0.0;
          nullspace[lrowdofs * 3 + firstlocaldofid + 7] = 0.0;
        }
      }
    }

    // get plain aggregation Ptent
    Teuchos::RCP<Epetra_CrsMatrix> crsPtent;
    MLAPI::GetPtent(*sysmat2->EpetraMatrix(), mlparams, nullspace, crsPtent);
    Core::LinAlg::SparseMatrix Ptent(crsPtent, Core::LinAlg::View);

    // compute scale-separation matrix: S = Ptent*Ptent^T
    sep_enr_ = Core::LinAlg::Multiply(Ptent, false, Ptent, true);
    sep_enr_->Complete();
  }

  return;
}


//----------------------------------------------------------------------*/
//----------------------------------------------------------------------*/
void FLD::UTILS::SetupFluidFluidVelPresSplit(const Discret::Discretization& fluiddis, int ndim,
    const Discret::Discretization& alefluiddis, Core::LinAlg::MapExtractor& extractor,
    Teuchos::RCP<Epetra_Map> fullmap)
{
  std::set<int> veldofset;
  std::set<int> presdofset;

  // for fluid elements
  int numfluidrownodes = fluiddis.NumMyRowNodes();
  for (int i = 0; i < numfluidrownodes; ++i)
  {
    Core::Nodes::Node* fluidnode = fluiddis.lRowNode(i);

    std::vector<int> fluiddof = fluiddis.Dof(0, fluidnode);
    for (unsigned j = 0; j < fluiddof.size(); ++j)
    {
      // test for dof position
      if (j < static_cast<unsigned>(ndim))
      {
        veldofset.insert(fluiddof[j]);
      }
      else
      {
        presdofset.insert(fluiddof[j]);
      }
    }
  }

  // for ale_fluid elements
  int numalefluidrownodes = alefluiddis.NumMyRowNodes();
  for (int i = 0; i < numalefluidrownodes; ++i)
  {
    Core::Nodes::Node* alefluidnode = alefluiddis.lRowNode(i);

    std::vector<int> alefluiddof = alefluiddis.Dof(alefluidnode);
    for (unsigned j = 0; j < alefluiddof.size(); ++j)
    {
      // test for dof position
      if (j < static_cast<unsigned>(ndim))
      {
        veldofset.insert(alefluiddof[j]);
      }
      else
      {
        presdofset.insert(alefluiddof[j]);
      }
    }
  }

  std::vector<int> veldofmapvec;
  veldofmapvec.reserve(veldofset.size());
  veldofmapvec.assign(veldofset.begin(), veldofset.end());
  veldofset.clear();
  Teuchos::RCP<Epetra_Map> velrowmap = Teuchos::rcp(
      new Epetra_Map(-1, veldofmapvec.size(), veldofmapvec.data(), 0, fluiddis.Comm()));
  veldofmapvec.clear();

  std::vector<int> presdofmapvec;
  presdofmapvec.reserve(presdofset.size());
  presdofmapvec.assign(presdofset.begin(), presdofset.end());
  presdofset.clear();
  Teuchos::RCP<Epetra_Map> presrowmap = Teuchos::rcp(
      new Epetra_Map(-1, presdofmapvec.size(), presdofmapvec.data(), 0, alefluiddis.Comm()));
  extractor.Setup(*fullmap, presrowmap, velrowmap);
}



// -------------------------------------------------------------------
// compute forces and moments                          rasthofer 08/13
// -------------------------------------------------------------------
void FLD::UTILS::LiftDrag(const Teuchos::RCP<const Discret::Discretization> dis,
    const Teuchos::RCP<const Epetra_Vector> trueresidual,
    const Teuchos::RCP<const Epetra_Vector> dispnp, const int ndim,
    Teuchos::RCP<std::map<int, std::vector<double>>>& liftdragvals, bool alefluid)
{
  int myrank = dis->Comm().MyPID();

  std::map<const int, std::set<Core::Nodes::Node*>> ldnodemap;
  std::map<const int, const std::vector<double>*> ldcoordmap;
  std::map<const int, const std::vector<double>*> ldaxismap;
  bool axis_for_moment = false;

  // allocate and initialise LiftDrag conditions
  std::vector<Core::Conditions::Condition*> ldconds;
  dis->GetCondition("LIFTDRAG", ldconds);

  // there is an L&D condition if it has a size
  if (ldconds.size())
  {
    // vector with lift&drag forces after communication
    liftdragvals = Teuchos::rcp(new std::map<int, std::vector<double>>);

    for (unsigned i = 0; i < ldconds.size(); ++i)  // loop L&D conditions (i.e. lines in .dat file)
    {
      // get label of present LiftDrag condition
      const int label = ldconds[i]->parameters().Get<int>("label");

      ((*liftdragvals))
          .insert(std::pair<int, std::vector<double>>(label, std::vector<double>(6, 0.0)));
    }

    // prepare output
    if (myrank == 0)
    {
      std::cout << "Lift and drag calculation:"
                << "\n";
      if (ndim == 2)
      {
        std::cout << "lift'n'drag Id      F_x             F_y             M_z :"
                  << "\n";
      }
      if (ndim == 3)
      {
        std::cout << "lift'n'drag Id      F_x             F_y             F_z           ";
        std::cout << "M_x             M_y             M_z :"
                  << "\n";
      }
    }

    // sort data
    for (unsigned i = 0; i < ldconds.size(); ++i)  // loop L&D conditions (i.e. lines in .dat file)
    {
      // get label of present LiftDrag condition
      const int label = ldconds[i]->parameters().Get<int>("label");

      /* get new nodeset for new label OR:
         return pointer to nodeset for known label ... */
      std::set<Core::Nodes::Node*>& nodes = ldnodemap[label];

      // center coordinates to present label
      ldcoordmap[label] = &ldconds[i]->parameters().Get<std::vector<double>>("centerCoord");

      // axis of rotation for present label (only needed for 3D)
      if (ldconds[i]->Type() == Core::Conditions::SurfLIFTDRAG)
      {
        ldaxismap[label] = &ldconds[i]->parameters().Get<std::vector<double>>("axis");
        // get pointer to axis vector (if available)
        const std::vector<double>* axisvecptr = ldaxismap[label];
        if (axisvecptr->size() != 3) FOUR_C_THROW("axis vector has not length 3");
        Core::LinAlg::Matrix<3, 1> axisvec(axisvecptr->data(), false);
        if (axisvec.Norm2() > 1.0e-9) axis_for_moment = true;  // axis has been set
      }

      // get pointer to its nodal Ids
      const std::vector<int>* ids = ldconds[i]->GetNodes();

      /* put all nodes belonging to the L&D line or surface into 'nodes' which are
         associated with the present label */
      for (unsigned j = 0; j < ids->size(); ++j)
      {
        // give me present node Id
        const int node_id = (*ids)[j];
        // put it into nodeset of actual label if node is new and mine
        if (dis->HaveGlobalNode(node_id) && dis->gNode(node_id)->Owner() == myrank)
          nodes.insert(dis->gNode(node_id));
      }
    }  // end loop over conditions


    // now step the label map
    for (std::map<const int, std::set<Core::Nodes::Node*>>::const_iterator labelit =
             ldnodemap.begin();
         labelit != ldnodemap.end(); ++labelit)
    {
      const std::set<Core::Nodes::Node*>& nodes =
          labelit->second;                    // pointer to nodeset of present label
      const int label = labelit->first;       // the present label
      std::vector<double> myforces(3, 0.0);   // vector with lift&drag forces
      std::vector<double> mymoments(3, 0.0);  // vector with lift&drag moments

      // get also pointer to center coordinates
      const std::vector<double>* centerCoordvec = ldcoordmap[label];
      if (centerCoordvec->size() != 3) FOUR_C_THROW("axis vector has not length 3");
      Core::LinAlg::Matrix<3, 1> centerCoord(centerCoordvec->data(), false);

      // loop all nodes within my set
      for (std::set<Core::Nodes::Node*>::const_iterator actnode = nodes.begin();
           actnode != nodes.end(); ++actnode)
      {
        const Core::LinAlg::Matrix<3, 1> x(
            (*actnode)->X().data(), false);  // pointer to nodal coordinates
        const Epetra_BlockMap& rowdofmap = trueresidual->Map();
        const std::vector<int> dof = dis->Dof(*actnode);

        // get nodal forces
        Core::LinAlg::Matrix<3, 1> actforces(true);
        for (int idim = 0; idim < ndim; idim++)
        {
          actforces(idim, 0) = (*trueresidual)[rowdofmap.LID(dof[idim])];
          myforces[idim] += (*trueresidual)[rowdofmap.LID(dof[idim])];
        }
        // z-component remains zero for ndim=2

        // get moment
        Core::LinAlg::Matrix<3, 1> actmoments(true);
        // get vector of point to center point
        Core::LinAlg::Matrix<3, 1> distances;
        distances.Update(1.0, x, -1.0, centerCoord);

        // ALE case: take displacements into account
        if (alefluid)
        {
          if (dispnp == Teuchos::null) FOUR_C_THROW("Displacement expected for ale fluid!");
          for (int idim = 0; idim < ndim; idim++)
          {
            distances(idim, 0) += (*dispnp)[rowdofmap.LID(dof[idim])];
          }
        }

        // calculate nodal angular moment with respect to global coordinate system
        Core::LinAlg::Matrix<3, 1> actmoment_gc(true);
        actmoment_gc(0, 0) =
            distances(1) * actforces(2, 0) - distances(2) * actforces(1, 0);  // zero for 2D
        actmoment_gc(1, 0) =
            distances(2) * actforces(0, 0) - distances(0) * actforces(2, 0);  // zero for 2D
        actmoment_gc(2, 0) = distances(0) * actforces(1, 0) - distances(1) * actforces(0, 0);

        if (axis_for_moment)
        {
          const std::vector<double>* axisvecptr = ldaxismap[label];
          Core::LinAlg::Matrix<3, 1> axisvec(axisvecptr->data(), false);
          double norm = 0.0;
          if (axisvec.Norm2() != 0.0)
          {
            norm = axisvec.Norm2();
            // normed axis vector
            axisvec.Scale(1.0 / norm);
          }
          else
            FOUR_C_THROW("norm==0.0!");
          // projection of moment on given axis
          double mdir = actmoment_gc.Dot(axisvec);

          actmoments(2, 0) = mdir;
        }
        else
        {
          for (int idim = 0; idim < 3; idim++) actmoments(idim, 0) = actmoment_gc(idim, 0);
        }

        for (int idim = 0; idim < 3; idim++) mymoments[idim] += actmoments(idim, 0);
      }  // end: loop over nodes

      // care for the fact that we are (most likely) parallel
      trueresidual->Comm().SumAll(myforces.data(), ((*liftdragvals)[label]).data(), 3);
      trueresidual->Comm().SumAll(mymoments.data(), ((*liftdragvals)[label]).data() + 3, 3);

      // do the output
      if (myrank == 0)
      {
        if (ndim == 2)
        {
          std::cout << "     " << label << "         ";
          std::cout << std::scientific << ((*liftdragvals)[label])[0] << "    ";
          std::cout << std::scientific << ((*liftdragvals)[label])[1] << "    ";
          std::cout << std::scientific << ((*liftdragvals)[label])[5];
          std::cout << "\n";
        }
        if (ndim == 3)
        {
          std::cout << "     " << label << "         ";
          std::cout << std::scientific << ((*liftdragvals)[label])[0] << "    ";
          std::cout << std::scientific << ((*liftdragvals)[label])[1] << "    ";
          std::cout << std::scientific << ((*liftdragvals)[label])[2] << "    ";
          std::cout << std::scientific << ((*liftdragvals)[label])[3] << "    ";
          std::cout << std::scientific << ((*liftdragvals)[label])[4] << "    ";
          std::cout << std::scientific << ((*liftdragvals)[label])[5];
          std::cout << "\n";
        }
      }
    }  // end: loop over L&D labels
    if (myrank == 0)
    {
      std::cout << "\n";
    }
  }

  return;
}

// -------------------------------------------------------------------
// write forces and moments to file                    rasthofer 08/13
// -------------------------------------------------------------------
void FLD::UTILS::WriteLiftDragToFile(
    const double time, const int step, const std::map<int, std::vector<double>>& liftdragvals)
{
  // print to file
  std::ostringstream header;
  header << std::right << std::setw(16) << "Time" << std::right << std::setw(10) << "Step"
         << std::right << std::setw(10) << "Label" << std::right << std::setw(16) << "F_x"
         << std::right << std::setw(16) << "F_y" << std::right << std::setw(16) << "F_z"
         << std::right << std::setw(16) << "M_x" << std::right << std::setw(16) << "M_y"
         << std::right << std::setw(16) << "M_z";


  for (std::map<int, std::vector<double>>::const_iterator liftdragval = liftdragvals.begin();
       liftdragval != liftdragvals.end(); ++liftdragval)
  {
    std::ostringstream s;
    s << std::right << std::setw(16) << std::scientific << time << std::right << std::setw(10)
      << std::scientific << step << std::right << std::setw(10) << std::scientific
      << liftdragval->first << std::right << std::setw(16) << std::scientific
      << liftdragval->second[0] << std::right << std::setw(16) << std::scientific
      << liftdragval->second[1] << std::right << std::setw(16) << std::scientific
      << liftdragval->second[2] << std::right << std::setw(16) << std::scientific
      << liftdragval->second[3] << std::right << std::setw(16) << std::scientific
      << liftdragval->second[4] << std::right << std::setw(16) << std::scientific
      << liftdragval->second[5];

    std::ostringstream slabel;
    slabel << std::setw(3) << std::setfill('0') << liftdragval->first;
    std::ofstream f;
    const std::string fname = Global::Problem::Instance()->OutputControlFile()->FileName() +
                              ".liftdrag_label_" + slabel.str() + ".txt";

    if (step <= 1)
    {
      f.open(fname.c_str(), std::fstream::trunc);
      f << header.str() << std::endl;
    }
    else
    {
      f.open(fname.c_str(), std::fstream::ate | std::fstream::app);
    }

    f << s.str() << "\n";
    f.close();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<int, double> FLD::UTILS::ComputeFlowRates(Discret::Discretization& dis,
    const Teuchos::RCP<Epetra_Vector>& velnp, const std::string& condstring,
    const Inpar::FLUID::PhysicalType physicaltype)
{
  return ComputeFlowRates(dis, velnp, Teuchos::null, Teuchos::null, condstring, physicaltype);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<int, double> FLD::UTILS::ComputeFlowRates(Discret::Discretization& dis,
    const Teuchos::RCP<Epetra_Vector>& velnp, const Teuchos::RCP<Epetra_Vector>& gridv,
    const Teuchos::RCP<Epetra_Vector>& dispnp, const std::string& condstring,
    const Inpar::FLUID::PhysicalType physicaltype)
{
  Teuchos::ParameterList eleparams;
  // set action for elements
  eleparams.set<int>("action", FLD::calc_flowrate);
  eleparams.set<int>("Physical Type", physicaltype);

  // note that the flowrate is not yet divided by the area
  std::map<int, double> volumeflowrateperline;

  // get condition
  std::vector<Core::Conditions::Condition*> conds;
  dis.GetCondition(condstring, conds);

  // each condition is on every proc , but might not have condition elements there
  for (std::vector<Core::Conditions::Condition*>::const_iterator conditer = conds.begin();
       conditer != conds.end(); ++conditer)
  {
    const Core::Conditions::Condition* cond = *conditer;
    const int condID = cond->parameters().Get<int>("ConditionID");

    // get a vector layout from the discretization to construct matching
    // vectors and matrices local <-> global dof numbering
    const Epetra_Map* dofrowmap = dis.dof_row_map();

    // create vector (+ initialization with zeros)
    Teuchos::RCP<Epetra_Vector> flowrates = Core::LinAlg::CreateVector(*dofrowmap, true);

    // call loop over elements
    dis.ClearState();

    dis.set_state("velaf", velnp);
    if (dispnp != Teuchos::null) dis.set_state("dispnp", dispnp);
    if (gridv != Teuchos::null) dis.set_state("gridv", gridv);

    dis.evaluate_condition(eleparams, flowrates, condstring, condID);
    dis.ClearState();

    double local_flowrate = 0.0;
    for (int i = 0; i < dofrowmap->NumMyElements(); i++)
    {
      local_flowrate += ((*flowrates)[i]);
    }

    double flowrate = 0.0;
    dofrowmap->Comm().SumAll(&local_flowrate, &flowrate, 1);

    // if(dofrowmap->Comm().MyPID()==0)
    // std::cout << "gobal flow rate = " << flowrate << "\t condition ID = " << condID << std::endl;

    // ATTENTION: new definition: outflow is positive and inflow is negative
    volumeflowrateperline[condID] = flowrate;
  }
  return volumeflowrateperline;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<int, double> FLD::UTILS::compute_volume(Discret::Discretization& dis,
    const Teuchos::RCP<Epetra_Vector>& velnp, const Teuchos::RCP<Epetra_Vector>& gridv,
    const Teuchos::RCP<Epetra_Vector>& dispnp, const Inpar::FLUID::PhysicalType physicaltype)
{
  Teuchos::ParameterList eleparams;
  // set action for elements
  eleparams.set<int>("action", FLD::calc_volume);
  eleparams.set<int>("Physical Type", physicaltype);

  std::map<int, double> volumeperline;

  // call loop over elements
  dis.ClearState();
  dis.set_state("velnp", velnp);
  if (dispnp != Teuchos::null) dis.set_state("dispnp", dispnp);
  if (gridv != Teuchos::null) dis.set_state("gridv", gridv);

  Teuchos::RCP<Core::LinAlg::SerialDenseVector> volumes =
      Teuchos::rcp(new Core::LinAlg::SerialDenseVector(1));

  // call loop over elements (assemble nothing)
  dis.EvaluateScalars(eleparams, volumes);
  dis.ClearState();

  volumeperline[0] = (*volumes)(0);

  return volumeperline;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<int, Core::LinAlg::Matrix<3, 1>> FLD::UTILS::ComputeSurfaceImpulsRates(
    Discret::Discretization& dis, const Teuchos::RCP<Epetra_Vector> velnp)
{
  Teuchos::ParameterList eleparams;
  // set action for elements
  eleparams.set("action", "calc_impuls_rate");

  std::map<int, Core::LinAlg::Matrix<3, 1>> volumeflowratepersurface;

  // get condition
  std::vector<Core::Conditions::Condition*> conds;
  dis.GetCondition("SurfImpulsRate", conds);

  // collect elements by xfem coupling label
  for (std::vector<Core::Conditions::Condition*>::const_iterator conditer = conds.begin();
       conditer != conds.end(); ++conditer)
  {
    const Core::Conditions::Condition* cond = *conditer;

    const int condID = cond->parameters().Get<int>("ConditionID");

    // create vector (+ initialization with zeros)
    const Epetra_BlockMap mappy = velnp->Map();
    Teuchos::RCP<Epetra_Vector> impulsrates = Teuchos::rcp(new Epetra_Vector(mappy));
    impulsrates->PutScalar(0.0);

    // call loop over elements
    dis.ClearState();
    dis.set_state("velnp", velnp);
    dis.evaluate_condition(eleparams, impulsrates, "SurfImpulsRate", condID);
    dis.ClearState();
    Core::LinAlg::Matrix<3, 1> locflowrate(true);
    for (int inode = 0; inode < dis.NumMyRowNodes(); inode++)
    {
      const Core::Nodes::Node* node = dis.lRowNode(inode);
      static std::vector<int> gdofs(4);
      dis.Dof(node, 0, gdofs);
      for (size_t isd = 0; isd < 3; isd++)
      {
        locflowrate(isd) += (*impulsrates)[dis.DofColMap()->LID(gdofs[isd])];
        //        std::cout << (*impulsrates)[dis.DofColMap()->LID(gdofs[isd])] << std::endl;
      }
    }

    //    Core::LinAlg::Matrix<3,1> flowrate(true);
    //    dofrowmap->Comm().SumAll(&locflowrate(0),&flowrate(0),1);
    //    dofrowmap->Comm().SumAll(&locflowrate(1),&flowrate(1),1);
    //    dofrowmap->Comm().SumAll(&locflowrate(2),&flowrate(2),1);
    //    std::cout << "locflowrate " << locflowrate << std::endl;
    if (volumeflowratepersurface.find(condID) == volumeflowratepersurface.end())
    {
      Core::LinAlg::Matrix<3, 1> tmp(true);
      volumeflowratepersurface.insert(std::make_pair(condID, tmp));
    }
    Core::LinAlg::Matrix<3, 1> tmp = volumeflowratepersurface[condID];
    tmp += locflowrate;
    volumeflowratepersurface[condID] = tmp;
  }

  return volumeflowratepersurface;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::UTILS::WriteDoublesToFile(
    const double time, const int step, const std::map<int, double>& data, const std::string& name)
{
  if (data.empty()) FOUR_C_THROW("data vector is empty");

  // print to file
  std::ostringstream header;
  header << std::right << std::setw(16) << "Time" << std::right << std::setw(10) << "Step"
         << std::right << std::setw(10) << "ID" << std::right << std::setw(16) << name;

  for (std::map<int, double>::const_iterator iter = data.begin(); iter != data.end(); ++iter)
  {
    std::ostringstream s;
    s << std::right << std::setw(16) << std::scientific << time << std::right << std::setw(10)
      << std::scientific << step << std::right << std::setw(10) << std::scientific << iter->first
      << std::right << std::setw(29) << std::setprecision(14) << std::scientific << iter->second;

    std::ostringstream slabel;
    slabel << std::setw(3) << std::setfill('0') << iter->first;
    std::ofstream f;
    const std::string fname = Global::Problem::Instance()->OutputControlFile()->FileName() + "." +
                              name + "_ID_" + slabel.str() + ".txt";

    if (step <= 1)
      f.open(fname.c_str(), std::fstream::trunc);  // f << header.str() << std::endl;
    else
      f.open(fname.c_str(), std::fstream::ate | std::fstream::app);

    f << s.str() << "\n";
    f.close();
  }
}


/*----------------------------------------------------------------------*|
 | project vel gradient and store it in given param list        bk 05/15 |
 *----------------------------------------------------------------------*/
void FLD::UTILS::ProjectGradientAndSetParam(Teuchos::RCP<Discret::Discretization> discret,
    Teuchos::ParameterList& eleparams, Teuchos::RCP<const Epetra_Vector> vel,
    const std::string paraname, bool alefluid)
{
  // project gradient
  Teuchos::RCP<Epetra_MultiVector> projected_velgrad =
      FLD::UTILS::ProjectGradient(discret, vel, alefluid);

  // store multi vector in parameter list after export to col layout
  if (projected_velgrad != Teuchos::null)
    discret->add_multi_vector_to_parameter_list(eleparams, paraname, projected_velgrad);

  return;
}


/*----------------------------------------------------------------------*|
 | Project velocity gradient                                    bk 05/15 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> FLD::UTILS::ProjectGradient(
    Teuchos::RCP<Discret::Discretization> discret, Teuchos::RCP<const Epetra_Vector> vel,
    bool alefluid)
{
  // reconstruction of second derivatives for fluid residual
  Inpar::FLUID::GradientReconstructionMethod recomethod =
      Core::UTILS::IntegralValue<Inpar::FLUID::GradientReconstructionMethod>(
          Global::Problem::Instance()->FluidDynamicParams(), "VELGRAD_PROJ_METHOD");

  const int dim = Global::Problem::Instance()->NDim();
  const int numvec = dim * dim;
  Teuchos::ParameterList params;
  Teuchos::RCP<Epetra_MultiVector> projected_velgrad = Teuchos::null;

  // dependent on the desired projection, just remove this line
  if (not vel->Map().SameAs(*discret->dof_row_map()))
    FOUR_C_THROW("input map is not a dof row map of the fluid");

  switch (recomethod)
  {
    case Inpar::FLUID::gradreco_none:
    {
      // no projection and no parameter in parameter list
    }
    break;
    case Inpar::FLUID::gradreco_spr:
    {
      if (alefluid)
        FOUR_C_THROW(
            "ale fluid is currently not supported everywhere for superconvergent patch recovery, "
            "but it is easy to implement");
      params.set<int>("action", FLD::calc_velgrad_ele_center);
      // project velocity gradient of fluid to nodal level via superconvergent patch recovery
      switch (dim)
      {
        case 3:
          projected_velgrad = Core::FE::compute_superconvergent_patch_recovery<3>(
              *discret, *vel, "vel", numvec, params);
          break;
        case 2:
          projected_velgrad = Core::FE::compute_superconvergent_patch_recovery<2>(
              *discret, *vel, "vel", numvec, params);
          break;
        case 1:
          projected_velgrad = Core::FE::compute_superconvergent_patch_recovery<1>(
              *discret, *vel, "vel", numvec, params);
          break;
        default:
          FOUR_C_THROW("only 1/2/3D implementation available for superconvergent patch recovery");
          break;
      }
    }
    break;
    case Inpar::FLUID::gradreco_l2:
    {
      const int solvernumber =
          Global::Problem::Instance()->FluidDynamicParams().get<int>("VELGRAD_PROJ_SOLVER");
      if (solvernumber < 1) FOUR_C_THROW("you have to specify a VELGRAD_PROJ_SOLVER");
      const auto& solverparams = Global::Problem::Instance()->SolverParams(solvernumber);

      params.set<int>("action", FLD::velgradient_projection);

      // set given state for element evaluation
      discret->ClearState();
      discret->set_state("vel", vel);

      // project velocity gradient of fluid to nodal level via L2 projection
      projected_velgrad =
          Core::FE::compute_nodal_l2_projection(discret, "vel", numvec, params, solverparams);
    }
    break;
    default:
      FOUR_C_THROW("desired projection method not available");
      break;
  }

  return projected_velgrad;
}

FOUR_C_NAMESPACE_CLOSE
