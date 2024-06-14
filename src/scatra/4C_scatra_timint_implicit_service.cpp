/*----------------------------------------------------------------------*/
/*! \file
\brief Service routines of the scalar transport time integration class

\level 1


*/
/*----------------------------------------------------------------------*/

#include "4C_coupling_adapter.hpp"
#include "4C_fem_general_utils_superconvergent_patch_recovery.hpp"
#include "4C_fem_nurbs_discretization.hpp"
#include "4C_fluid_rotsym_periodicbc_utils.hpp"
#include "4C_fluid_turbulence_dyn_smag.hpp"
#include "4C_fluid_turbulence_dyn_vreman.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_gmsh.hpp"
#include "4C_io_runtime_csv_writer.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_timint_meshtying_strategy_base.hpp"
#include "4C_scatra_turbulence_hit_scalar_forcing.hpp"
#include "4C_utils_parameter_list.hpp"

#include <MLAPI_Aggregation.h>
#include <MLAPI_Workspace.h>

FOUR_C_NAMESPACE_OPEN


/*==========================================================================*
 |                                                                          |
 | public:                                                                  |
 |                                                                          |
 *==========================================================================*/

/*==========================================================================*/
//! @name general framework
/*==========================================================================*/

/*----------------------------------------------------------------------*
 | calculate fluxes inside domain and/or on boundary         fang 07/16 |
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::CalcFlux(const bool writetofile)
{
  switch (calcflux_domain_)
  {
    case Inpar::ScaTra::flux_diffusive:
    case Inpar::ScaTra::flux_total:
    {
      flux_domain_ = CalcFluxInDomain();
      break;
    }

    case Inpar::ScaTra::flux_none:
    {
      // do nothing
      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid option for flux calculation inside domain!");
      break;
    }
  }

  switch (calcflux_boundary_)
  {
    case Inpar::ScaTra::flux_convective:
    case Inpar::ScaTra::flux_diffusive:
    case Inpar::ScaTra::flux_total:
    {
      // calculate normal flux vector field only for the user-defined boundary conditions:
      flux_boundary_ = CalcFluxAtBoundary(writetofile);

      break;
    }

    case Inpar::ScaTra::flux_none:
    {
      // do nothing
      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid option for flux calculation on boundary!");
      break;
    }
  }
}


/*----------------------------------------------------------------------*
 | calculate flux vector field inside computational domain   fang 08/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> ScaTra::ScaTraTimIntImpl::CalcFluxInDomain()
{
  // extract dofrowmap from discretization
  const Epetra_Map& dofrowmap = *discret_->dof_row_map();

  // initialize global flux vectors
  Teuchos::RCP<Epetra_MultiVector> flux = Teuchos::rcp(new Epetra_MultiVector(dofrowmap, 3));
  Teuchos::RCP<Epetra_MultiVector> flux_projected =
      Teuchos::rcp(new Epetra_MultiVector(dofrowmap, 3));

  // parameter list for element evaluation
  Teuchos::ParameterList params;

  // set action for elements
  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::calc_flux_domain, params);

  // provide discretization with state vector
  discret_->set_state("phinp", phinp_);

  // evaluate flux vector field inside the whole computational domain (e.g., for visualization of
  // particle path lines)
  discret_->Evaluate(params, Teuchos::null, Teuchos::null, Teuchos::rcp((*flux)(0), false),
      Teuchos::rcp((*flux)(1), false), Teuchos::rcp((*flux)(2), false));

  if (calcflux_domain_lumped_)
  {
    // vector for integrated shape functions
    Teuchos::RCP<Epetra_Vector> integratedshapefcts = Core::LinAlg::CreateVector(dofrowmap);

    // overwrite action for elements
    Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
        "action", ScaTra::Action::integrate_shape_functions, params);

    // integrate shape functions
    Teuchos::RCP<Core::LinAlg::IntSerialDenseVector> dofids =
        Teuchos::rcp(new Core::LinAlg::IntSerialDenseVector(NumDofPerNode()));
    dofids->putScalar(1);  // integrate shape functions for all dofs
    params.set("dofids", dofids);
    discret_->Evaluate(
        params, Teuchos::null, Teuchos::null, integratedshapefcts, Teuchos::null, Teuchos::null);

    // insert values into final flux vector for visualization
    // We do not solve a global, linear system of equations (exact L2 projection with good
    // consistency), but perform mass matrix lumping, i.e., we divide by the values of the
    // integrated shape functions
    if (flux_projected->ReciprocalMultiply(1., *integratedshapefcts, *flux, 0.))
      FOUR_C_THROW("ReciprocalMultiply failed!");
  }

  else
  {
    // solve global, linear system of equations without lumping global mass matrix
    // overwrite action for elements
    Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
        "action", ScaTra::Action::calc_mass_matrix, params);

    // initialize global mass matrix
    Teuchos::RCP<Core::LinAlg::SparseMatrix> massmatrix =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(dofrowmap, 27, false));

    // call loop over elements
    discret_->Evaluate(params, massmatrix, Teuchos::null);

    // finalize global mass matrix
    massmatrix->Complete();

    // compute flux vector field
    Core::LinAlg::SolverParams solver_params;
    solver_params.refactor = true;
    solver_params.reset = true;
    solver_->Solve(massmatrix->EpetraOperator(), flux_projected, flux, solver_params);

    // reset solver
    solver_->Reset();
  }

  return flux_projected;
}  // ScaTra::ScaTraTimIntImpl::CalcFluxInDomain


/*----------------------------------------------------------------------*
 |  calculate mass / heat normal flux at specified boundaries  gjb 06/09|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> ScaTra::ScaTraTimIntImpl::CalcFluxAtBoundary(
    const bool writetofile)
{
  // The normal flux calculation is based on the idea proposed in
  // GRESHO ET AL.,
  // "THE CONSISTENT GALERKIN FEM FOR COMPUTING DERIVED BOUNDARY
  // QUANTITIES IN THERMAL AND/OR FLUIDS PROBLEMS",
  // INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN FLUIDS, VOL. 7, 371-394 (1987)

  // empty vector for (normal) mass or heat flux vectors (always 3D)
  Teuchos::RCP<Epetra_MultiVector> flux =
      Teuchos::rcp(new Epetra_MultiVector(*discret_->dof_row_map(), 3, true));

  // determine the averaged normal vector field for indicated boundaries
  // used for the output of the normal flux as a vector field
  // is computed only once; for ALE formulation recalculation is necessary
  if ((normals_ == Teuchos::null) or isale_)
    normals_ = compute_normal_vectors(std::vector<std::string>(1, "ScaTraFluxCalc"));

  if (calcflux_boundary_ == Inpar::ScaTra::flux_convective)
  {
    // zero out trueresidual vector -> we do not need this info
    trueresidual_->PutScalar(0.0);
  }
  else
  {
    // was the residual already prepared?
    if ((solvtype_ != Inpar::ScaTra::solvertype_nonlinear) and (lastfluxoutputstep_ != step_))
    {
      lastfluxoutputstep_ = step_;

      // For nonlinear problems we already have the actual residual vector
      // from the last convergence test!
      // For linear problems we have to compute this information first, since
      // the residual (w.o. Neumann boundary) has not been computed after the last solve!

      // zero out matrix entries
      sysmat_->Zero();

      // zero out residual vector
      residual_->PutScalar(0.0);

      Teuchos::ParameterList eleparams;
      // action for elements
      Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
          "action", ScaTra::Action::calc_mat_and_rhs, eleparams);

      // other parameters that might be needed by the elements
      // the resulting system has to be solved incrementally
      set_element_time_parameter(true);

      // AVM3 separation for incremental solver: get fine-scale part of scalar
      if (incremental_ and (fssgd_ != Inpar::ScaTra::fssugrdiff_no or
                               turbmodel_ == Inpar::FLUID::multifractal_subgrid_scales))
        av_m3_separation();

      // add element parameters according to time-integration scheme
      // we want to have an incremental solver!!
      add_time_integration_specific_vectors(true);

      // add element parameters according to specific problem
      add_problem_specific_parameters_and_vectors(eleparams);

      {
        // call standard loop over elements
        discret_->Evaluate(
            eleparams, sysmat_, Teuchos::null, residual_, Teuchos::null, Teuchos::null);
      }

      // scaling to get true residual vector for all time integration schemes
      trueresidual_->Update(residual_scaling(), *residual_, 0.0);

      // undo potential changes
      set_element_time_parameter();

    }  // if ((solvtype_!=Inpar::ScaTra::solvertype_nonlinear) && (lastfluxoutputstep_ != step_))
  }    // if (writeflux_==Inpar::ScaTra::flux_convective_boundary)

  // if desired add the convective flux contribution
  // to the trueresidual_ now.
  if (calcflux_boundary_ == Inpar::ScaTra::flux_total or
      calcflux_boundary_ == Inpar::ScaTra::flux_convective)
  {
    // following screen output removed, for the time being
    /*if(myrank_ == 0)
    {
      std::cout << "Convective flux contribution is added to trueresidual_ vector." << std::endl;
      std::cout << "Be sure not to address the same boundary part twice!" << std::endl;
      std::cout << "Two flux calculation boundaries should also not share a common node!" <<
    std::endl;
    }*/

    std::vector<Core::Conditions::Condition*> cond;
    discret_->GetCondition("ScaTraFluxCalc", cond);

    Teuchos::ParameterList params;
    Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
        "action", ScaTra::BoundaryAction::add_convective_mass_flux, params);

    // add element parameters according to time-integration scheme
    add_time_integration_specific_vectors();

    // call loop over boundary elements and add integrated fluxes to trueresidual_
    discret_->evaluate_condition(params, trueresidual_, "ScaTraFluxCalc");
  }

  // vector for effective flux over all defined boundary conditions
  // maximal number of fluxes  (numscal+1 -> ionic species + potential) is generated
  // for OUTPUT standard -> last entry is not used
  std::vector<double> normfluxsum(NumDofPerNode());

  std::vector<Core::Conditions::Condition*> cond;
  discret_->GetCondition("ScaTraFluxCalc", cond);

  // safety check
  if (!cond.size())
    FOUR_C_THROW("Flux output requested without corresponding boundary condition specification!");

  if (myrank_ == 0)
  {
    std::cout << std::endl;
    std::cout << "Normal fluxes at boundary 'ScaTraFluxCalc' on discretization '"
              << discret_->Name() << "':" << std::endl;
    std::cout
        << "+----+-----+-------------------------+------------------+--------------------------+"
        << std::endl;
    std::cout
        << "| ID | DOF | Integral of normal flux | Area of boundary | Mean normal flux density |"
        << std::endl;
  }

  // first, add to all conditions of interest a ConditionID
  for (int condid = 0; condid < (int)cond.size(); condid++)
  {
    // is there already a ConditionID?
    const auto* CondID = cond[condid]->parameters().GetIf<int>("ConditionID");
    if (CondID)
    {
      if ((*CondID) != condid)
        FOUR_C_THROW("Condition 'ScaTraFluxCalc' has non-matching ConditionID!");
    }
    else
    {
      // let's add a ConditionID
      cond[condid]->parameters().Add("ConditionID", condid);
    }
  }

  // now we evaluate the conditions and separate via ConditionID
  for (int icond = 0; icond < static_cast<int>(cond.size()); ++icond)
  {
    // extract dofrowmap associated with current boundary segment
    const Epetra_Map& dofrowmap = *flux_boundary_maps_->Map(icond + 1);

    // extract part of true residual vector associated with current boundary segment
    const Teuchos::RCP<Epetra_Vector> trueresidual_boundary =
        flux_boundary_maps_->ExtractVector(*trueresidual_, icond + 1);

    // initialize vector for nodal values of normal boundary fluxes
    Teuchos::RCP<Epetra_Vector> normalfluxes = Core::LinAlg::CreateVector(dofrowmap);

    // create parameter list for boundary elements
    Teuchos::ParameterList params;

    // initialize variable for value of boundary integral
    double boundaryint(-1.);

    // compute nodal values of normal boundary fluxes
    if (calcflux_boundary_lumped_)
    {
      // lump boundary mass matrix instead of solving a small, linear system of equations
      // calculate integral of shape functions over indicated boundary and its area
      params.set("area", 0.0);
      Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
          "action", ScaTra::BoundaryAction::integrate_shape_functions, params);

      // create vector (+ initialization with zeros)
      const Teuchos::RCP<Epetra_Vector> integratedshapefunc = Core::LinAlg::CreateVector(dofrowmap);

      // call loop over elements
      discret_->evaluate_condition(params, integratedshapefunc, "ScaTraFluxCalc", icond);

      // compute normal boundary fluxes
      for (int idof = 0; idof < dofrowmap.NumMyElements(); ++idof)
        (*normalfluxes)[idof] = (*trueresidual_boundary)[idof] / (*integratedshapefunc)[idof];

      // get value of boundary integral on this processor
      boundaryint = params.get<double>("area");
    }

    else
    {
      // solve small, linear system of equations without lumping boundary mass matrix
      // add action to parameter list
      Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
          "action", ScaTra::BoundaryAction::calc_mass_matrix, params);

      // initialize boundary mass matrix
      Teuchos::RCP<Core::LinAlg::SparseMatrix> massmatrix_boundary =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(dofrowmap, 27, false));

      // call loop over elements
      discret_->evaluate_condition(params, massmatrix_boundary, Teuchos::null, Teuchos::null,
          Teuchos::null, Teuchos::null, "ScaTraFluxCalc", icond);

      // finalize boundary mass matrix
      massmatrix_boundary->Complete();

      // compute normal boundary fluxes
      Core::LinAlg::SolverParams solver_params;
      solver_params.refactor = true;
      solver_params.reset = true;
      solver_->Solve(massmatrix_boundary->EpetraOperator(), normalfluxes, trueresidual_boundary,
          solver_params);

      // reset solver
      solver_->Reset();

      // overwrite action in parameter list
      Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
          "action", ScaTra::BoundaryAction::calc_boundary_integral, params);

      // initialize one-component result vector for value of boundary integral
      Teuchos::RCP<Core::LinAlg::SerialDenseVector> boundaryint_vector =
          Teuchos::rcp(new Core::LinAlg::SerialDenseVector(1));

      // compute value of boundary integral
      discret_->EvaluateScalars(params, boundaryint_vector, "ScaTraFluxCalc", icond);

      // extract value of boundary integral
      boundaryint = (*boundaryint_vector)(0);
    }

    // maximal number of fluxes
    std::vector<double> normfluxintegral(NumDofPerNode());

    // insert values into final flux vector for visualization
    for (int lnodid = 0; lnodid < discret_->NumMyRowNodes(); ++lnodid)
    {
      Core::Nodes::Node* actnode = discret_->lRowNode(lnodid);
      for (int idof = 0; idof < discret_->NumDof(0, actnode); ++idof)
      {
        const int dofgid = discret_->Dof(0, actnode, idof);
        if (dofrowmap.MyGID(dofgid))
        {
          const int doflid = dofrowmap.LID(dofgid);

          // compute integral value for every degree of freedom
          normfluxintegral[idof] += (*trueresidual_boundary)[doflid];

          // care for the slave nodes of rotationally symm. periodic boundary conditions
          double rotangle(0.0);
          bool havetorotate = FLD::IsSlaveNodeOfRotSymPBC(actnode, rotangle);

          // do not insert slave node values here, since they would overwrite the
          // master node values owning the same dof
          // (rotation of slave node vectors is performed later during output)
          if (not havetorotate)
          {
            // for visualization, we plot the normal flux with
            // outward pointing normal vector
            for (int idim = 0; idim < 3; idim++)
            {
              Epetra_Vector* normalcomp = (*normals_)(idim);
              double normalveccomp = (*normalcomp)[lnodid];
              int err =
                  flux->ReplaceGlobalValue(dofgid, idim, (*normalfluxes)[doflid] * normalveccomp);
              if (err != 0) FOUR_C_THROW("Detected error in ReplaceMyValue");
            }
          }
        }
      }
    }

    // care for the parallel case
    std::vector<double> parnormfluxintegral(NumDofPerNode());
    discret_->Comm().SumAll(normfluxintegral.data(), parnormfluxintegral.data(), NumDofPerNode());
    double parboundaryint = 0.0;
    discret_->Comm().SumAll(&boundaryint, &parboundaryint, 1);

    for (int idof = 0; idof < NumDofPerNode(); ++idof)
    {
      // print out results
      if (myrank_ == 0)
      {
        printf("| %2d | %2d  |       %+10.4E       |    %10.4E    |        %+10.4E       |\n",
            icond, idof, parnormfluxintegral[idof], parboundaryint,
            parnormfluxintegral[idof] / parboundaryint);
      }
      normfluxsum[idof] += parnormfluxintegral[idof];
    }

    // statistics section for normfluxintegral
    if (do_boundary_flux_statistics())
    {
      // add current flux value to the sum!
      (*sumnormfluxintegral_)[icond] += parnormfluxintegral[0];  // only first scalar!
      int samstep = step_ - samstart_ + 1;

      // dump every dumperiod steps (i.e., write to screen)
      bool dumpit(false);
      if (dumperiod_ == 0)
      {
        dumpit = true;
      }
      else
      {
        if (samstep % dumperiod_ == 0) dumpit = true;
      }

      if (dumpit)
      {
        double meannormfluxintegral = (*sumnormfluxintegral_)[icond] / samstep;
        // dump statistical results
        if (myrank_ == 0)
        {
          printf("| %2d | Mean normal-flux integral (step %5d -- step %5d) :   %12.5E |\n", icond,
              samstart_, step_, meannormfluxintegral);
        }
      }
    }

    // print out results to file as well (only if really desired)
    if ((myrank_ == 0) and writetofile)
    {
      std::ostringstream temp;
      temp << icond;
      temp << discret_->Name();
      const std::string fname = problem_->OutputControlFile()->file_name() +
                                ".boundaryflux_ScaTraFluxCalc_" + temp.str() + ".txt";

      std::ofstream f;
      if (Step() <= 1)
      {
        f.open(fname.c_str(), std::fstream::trunc);
        f << "#| ID | Step | Time | Area of boundary |";
        for (int idof = 0; idof < NumDofPerNode(); ++idof)
        {
          f << " Integral of normal flux " << idof << " | Mean normal flux density " << idof
            << " |";
        }
        f << "\n";
      }
      else
        f.open(fname.c_str(), std::fstream::ate | std::fstream::app);

      f << icond << " " << Step() << " " << Time() << " " << parboundaryint << " ";
      for (int idof = 0; idof < NumDofPerNode(); ++idof)
      {
        f << parnormfluxintegral[idof] << " " << parnormfluxintegral[idof] / parboundaryint << " ";
      }
      f << "\n";
      f.flush();
      f.close();
    }  // write to file

  }  // loop over condid

  if (myrank_ == 0)
  {
    std::cout
        << "+----+-----+-------------------------+------------------+--------------------------+"
        << std::endl;
  }

  // print out the accumulated normal flux over all indicated boundaries
  // the accumulated normal flux over all indicated boundaries is not printed for numscal+1
  if (myrank_ == 0)
  {
    for (int idof = 0; idof < NumScal(); ++idof)
    {
      printf("| Sum of all normal flux boundary integrals for scalar %d: %+10.5E             |\n",
          idof, normfluxsum[idof]);
    }
    std::cout
        << "+----------------------------------------------------------------------------------+"
        << std::endl;
  }

  return flux;
}  // ScaTra::ScaTraTimIntImpl::CalcFluxAtBoundary


/*--------------------------------------------------------------------*
 | calculate initial time derivatives of state variables   fang 03/15 |
 *--------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::calc_initial_time_derivative()
{
  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + calculate initial time derivative");

  // initial screen output
  if (myrank_ == 0)
  {
    std::cout
        << "SCATRA: calculating initial time derivative of state variables on discretization \""
        << discret_->Name().c_str() << "\" (step " << Step() << ", time " << Time() << ") ... ... ";
  }

  // evaluate Dirichlet and Neumann boundary conditions at time t = 0 to ensure consistent
  // computation of initial time derivative vector Dirichlet boundary conditions should be
  // consistent with initial field
  apply_dirichlet_bc(time_, phinp_, Teuchos::null);
  compute_intermediate_values();
  apply_neumann_bc(neumann_loads_);

  // In most cases, the history vector is entirely zero at this point, since this routine is called
  // before the first time step. However, in levelset simulations with reinitialization, this
  // routine might be called before every single time step. For this reason, the history vector
  // needs to be manually set to zero here and restored at the end of this routine.
  Teuchos::RCP<Epetra_Vector> hist = hist_;
  hist_ = zeros_;

  // In a first step, we assemble the standard global system of equations.
  assemble_mat_and_rhs();

  // In a second step, we need to modify the assembled system of equations, since we want to solve
  // M phidt^0 = f^n - K\phi^n - C(u_n)\phi^n
  // In particular, we need to replace the global system matrix by a global mass matrix,
  // and we need to remove all transient contributions associated with time discretization from the
  // global residual vector.

  // reset global system matrix
  sysmat_->Zero();

  // create and fill parameter list for elements
  Teuchos::ParameterList eleparams;
  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::calc_initial_time_deriv, eleparams);
  add_problem_specific_parameters_and_vectors(eleparams);

  // add state vectors according to time integration scheme
  add_time_integration_specific_vectors();

  // modify global system of equations as explained above
  discret_->Evaluate(eleparams, sysmat_, residual_);

  // apply condensation to global system of equations if necessary
  strategy_->CondenseMatAndRHS(sysmat_, residual_, true);

  // finalize assembly of global mass matrix
  sysmat_->Complete();

  // apply artificial Dirichlet boundary conditions to system of equations to hold initial values of
  // auxiliary degrees of freedom (= non-transported scalars) constant when solving for initial time
  // derivatives of primary degrees of freedom (= transported scalars)
  if (splitter_ != Teuchos::null)
    Core::LinAlg::apply_dirichlet_to_system(
        *sysmat_, *phidtnp_, *residual_, *zeros_, *(splitter_->CondMap()));

  // solve global system of equations for initial time derivative of state variables
  Core::LinAlg::SolverParams solver_params;
  solver_params.refactor = true;
  solver_params.reset = true;
  solver_->Solve(sysmat_->EpetraOperator(), phidtnp_, residual_, solver_params);

  // ToDo: Impose initial time derivatives resulting from Dirichlet boundary conditions.
  // At the moment, the initial time derivatives do not take account of time-dependent
  // Dirichlet boundary conditions. So for example, if you consider a simple setup with
  // one single finite element exhibiting a zero initial scalar field, and if you then
  // ramp up the scalar inside the entire element from zero to one linearly in time via
  // a time-dependent Dirichlet boundary condition, you will currently get zero initial
  // time derivatives (although they should be positive) and, as a consequence, there will
  // be unphysical temporal oscillations in the time derivatives from time step to time
  // step. The following, currently commented line of code would solve this problem, but
  // also cause many existing test cases to fail. There are intricate cases where the
  // imposed initial scalar field is not consistent with Dirichlet boundary conditions,
  // and it is still not clear how to handle these cases. One should have a closer look
  // at this problem and try to activate the following line of code at some point.

  // apply_dirichlet_bc(time_,Teuchos::null,phidtnp_);

  // copy solution
  phidtn_->Update(1.0, *phidtnp_, 0.0);

  // reset global system matrix and its graph, since we solved a very special problem with a special
  // sparsity pattern
  sysmat_->Reset();

  // reset solver
  solver_->Reset();

  // reset true residual vector computed during assembly of the standard global system of equations,
  // since not yet needed
  trueresidual_->PutScalar(0.0);

  // restore history vector as explained above
  hist_ = hist;

  // final screen output
  if (myrank_ == 0) std::cout << "done!" << std::endl;
}  // ScaTra::ScaTraTimIntImpl::calc_initial_time_derivative


/*---------------------------------------------------------------------------*
 | compute nodal density values from nodal concentration values   fang 12/14 |
 *---------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::ComputeDensity()
{
  // loop over all nodes owned by current processor
  for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); ++lnodeid)
  {
    // get current node
    Core::Nodes::Node* lnode = discret_->lRowNode(lnodeid);

    // get associated degrees of freedom
    std::vector<int> nodedofs = discret_->Dof(0, lnode);
    const int numdof = static_cast<int>(nodedofs.size());

    // initialize nodal density value
    double density = 1.;

    // loop over all transported scalars
    for (int k = 0; k < NumScal(); ++k)
    {
      /*
        //                  k=NumScal()-1
        //          /       ----                         \
        //         |        \                            |
        // rho_0 * | 1 +    /      alpha_k * (c_k - c_0) |
        //         |        ----                         |
        //          \       k=0                          /
        //
        // For use of molar mass M_k:  alpha_k = M_k/rho_0
       */

      // global and local dof ID
      const int globaldofid = nodedofs[k];
      const int localdofid = Phiafnp()->Map().LID(globaldofid);
      if (localdofid < 0) FOUR_C_THROW("Local dof ID not found in dof map!");

      // add contribution of scalar k to nodal density value
      density += densific_[k] * ((*(Phiafnp()))[localdofid] - c0_[k]);
    }

    // insert nodal density value into global density vector
    // note that the density vector has been initialized with the dofrowmap of the discretization
    // in case there is more than one dof per node, the nodal density value is inserted into the
    // position of the last dof this way, all nodal density values will be correctly extracted in
    // the fluid algorithm
    const int globaldofid = nodedofs[numdof - 1];
    const int localdofid = Phiafnp()->Map().LID(globaldofid);
    if (localdofid < 0) FOUR_C_THROW("Local dof ID not found in dof map!");

    int err = densafnp_->ReplaceMyValue(localdofid, 0, density);

    if (err) FOUR_C_THROW("Error while inserting nodal density value into global density vector!");
  }  // loop over all local nodes
}  // ScaTra::ScaTraTimIntImpl::ComputeDensity


/*---------------------------------------------------------------------------*
 | compute time derivatives of discrete state variables           fang 01/17 |
 *---------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::compute_time_derivative()
{
  // compute time derivatives of discrete state variables associated with meshtying strategy
  strategy_->compute_time_derivative();
}


/*----------------------------------------------------------------------*
 | output total and mean values of transported scalars       fang 03/16 |
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::output_total_and_mean_scalars(const int num)
{
  if (outputscalars_ != Inpar::ScaTra::outputscalars_none)
  {
    if (outputscalarstrategy_ == Teuchos::null)
      FOUR_C_THROW("output strategy was not initialized!");
    outputscalarstrategy_->output_total_and_mean_scalars(this, num);
  }  // if(outputscalars_)
}  // ScaTra::ScaTraTimIntImpl::output_total_and_mean_scalars


/*--------------------------------------------------------------------------------------------------------*
 | output domain or boundary integrals, i.e., surface areas or volumes of specified nodesets   fang
 06/15 |
 *--------------------------------------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::output_domain_or_boundary_integrals(const std::string& condstring)
{
  // check whether output is applicable
  if ((computeintegrals_ == Inpar::ScaTra::computeintegrals_initial and Step() == 0) or
      computeintegrals_ == Inpar::ScaTra::computeintegrals_repeated)
  {
    if (outputdomainintegralstrategy_ == Teuchos::null)
      FOUR_C_THROW("output strategy was not initialized!");

    outputdomainintegralstrategy_->evaluate_integrals_and_print_results(this, condstring);
  }  // check whether output is applicable
}  // ScaTra::ScaTraTimIntImpl::output_domain_or_boundary_integrals


/*----------------------------------------------------------------------*
 | Evaluate surface/interface permeability for FS3I          Thon 11/14 |
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::SurfacePermeability(
    Teuchos::RCP<Core::LinAlg::SparseOperator> matrix, Teuchos::RCP<Epetra_Vector> rhs)
{
  // time measurement: evaluate condition 'SurfacePermeability'
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'ScaTraCoupling'");

  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
      "action", ScaTra::BoundaryAction::calc_fs3i_surface_permeability, condparams);

  // add element parameters according to time-integration scheme
  add_time_integration_specific_vectors();

  // replace phinp by the mean of phinp if this has been set
  if (mean_conc_ != Teuchos::null)
  {
    // TODO: (thon) this is not a nice way of using the mean concentration instead of phinp!
    // Note: meanconc_ is not cleared by calling 'discret_->ClearState()', hence if you
    // don't want to do this replacement here any more call 'ClearMeanConcentration()'
    discret_->set_state("phinp", mean_conc_);

    if (myrank_ == 0)
      std::cout << "Replacing 'phinp' by 'meanconc_' in the evaluation of the surface permeability"
                << std::endl;
  }

  if (membrane_conc_ == Teuchos::null)
    FOUR_C_THROW("Membrane concentration must already been saved before calling this function!");
  discret_->set_state("MembraneConcentration", membrane_conc_);

  // test if all necessary ingredients had been set
  if (not discret_->HasState(NdsWallShearStress(), "WallShearStress"))
    FOUR_C_THROW(
        "WSS must already been set into one of the secondary dofset before calling this function!");

  // Evaluate condition
  discret_->evaluate_condition(
      condparams, matrix, Teuchos::null, rhs, Teuchos::null, Teuchos::null, "ScaTraCoupling");

  matrix->Complete();
}  // ScaTra::ScaTraTimIntImpl::SurfacePermeability

/*----------------------------------------------------------------------------*
 |  Kedem & Katchalsky second membrane equation for FPS3i       hemmler 07/14 |
 |  see e.g. Kedem, O. T., and A. Katchalsky. "Thermodynamic analysis of the permeability of
 biological membranes to non-electrolytes." Biochimica et biophysica Acta 27 (1958): 229-246.
 *----------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::KedemKatchalsky(
    Teuchos::RCP<Core::LinAlg::SparseOperator> matrix, Teuchos::RCP<Epetra_Vector> rhs)
{
  // time measurement: evaluate condition 'SurfacePermeability'
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'ScaTraCoupling'");

  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
      "action", ScaTra::BoundaryAction::calc_fps3i_surface_permeability, condparams);

  // add element parameters according to time-integration scheme
  add_time_integration_specific_vectors();

  // test if all necessary ingredients for the second Kedem-Katchalsky equations had been set
  if (not discret_->HasState(NdsWallShearStress(), "WallShearStress"))
    FOUR_C_THROW(
        "WSS must already been set into one of the secondary dofset before calling this function!");

  if (not discret_->HasState(NdsPressure(), "Pressure"))
  {
    FOUR_C_THROW(
        "Pressure must already been set into one of the secondary dofset before calling this "
        "function!");
  }

  if (membrane_conc_ == Teuchos::null)
    FOUR_C_THROW("Membrane concentration must already been saved before calling this function!");
  discret_->set_state("MembraneConcentration", membrane_conc_);

  // Evaluate condition
  discret_->evaluate_condition(
      condparams, matrix, Teuchos::null, rhs, Teuchos::null, Teuchos::null, "ScaTraCoupling");

  matrix->Complete();
}  // ScaTra::ScaTraTimIntImpl::KedemKatchalsky


/*==========================================================================*
 |                                                                          |
 | protected:                                                               |
 |                                                                          |
 *==========================================================================*/

/*==========================================================================*/
// general framework
/*==========================================================================*/

/*----------------------------------------------------------------------*
 | add approximation to flux vectors to a parameter list      gjb 05/10 |
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::add_flux_approx_to_parameter_list(Teuchos::ParameterList& p)
{
  Teuchos::RCP<Epetra_MultiVector> flux = CalcFluxInDomain();

  // post_ensight does not support multivectors based on the dofmap
  // for now, I create single vectors that can be handled by the filters

  // get the noderowmap
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  Teuchos::RCP<Epetra_MultiVector> fluxk =
      Teuchos::rcp(new Epetra_MultiVector(*noderowmap, 3, true));
  for (int k = 0; k < NumScal(); ++k)
  {
    std::ostringstream temp;
    temp << k;
    std::string name = "flux_phi_" + temp.str();
    for (int i = 0; i < fluxk->MyLength(); ++i)
    {
      Core::Nodes::Node* actnode = discret_->lRowNode(i);
      int dofgid = discret_->Dof(0, actnode, k);
      fluxk->ReplaceMyValue(i, 0, ((*flux)[0])[(flux->Map()).LID(dofgid)]);
      fluxk->ReplaceMyValue(i, 1, ((*flux)[1])[(flux->Map()).LID(dofgid)]);
      fluxk->ReplaceMyValue(i, 2, ((*flux)[2])[(flux->Map()).LID(dofgid)]);
    }
    discret_->add_multi_vector_to_parameter_list(p, name, fluxk);
  }
}  // ScaTra::ScaTraTimIntImpl::add_flux_approx_to_parameter_list


/*----------------------------------------------------------------------*
 | compute outward pointing unit normal vectors at given b.c.  gjb 01/09|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> ScaTra::ScaTraTimIntImpl::compute_normal_vectors(
    const std::vector<std::string>& condnames)
{
  // create vectors for x,y and z component of average normal vector field
  // get noderowmap of discretization
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  Teuchos::RCP<Epetra_MultiVector> normal =
      Teuchos::rcp(new Epetra_MultiVector(*noderowmap, 3, true));

  // set action for elements
  Teuchos::ParameterList eleparams;
  Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
      "action", ScaTra::BoundaryAction::calc_normal_vectors, eleparams);
  eleparams.set<Teuchos::RCP<Epetra_MultiVector>>("normal vectors", normal);

  // loop over all intended types of conditions
  for (const auto& condname : condnames)
  {
    discret_->evaluate_condition(eleparams, condname);
  }

  // the normal vector field is not properly scaled up to now. We do this here
  int numrownodes = discret_->NumMyRowNodes();
  Epetra_Vector* xcomp = (*normal)(0);
  Epetra_Vector* ycomp = (*normal)(1);
  Epetra_Vector* zcomp = (*normal)(2);
  for (int i = 0; i < numrownodes; ++i)
  {
    double x = (*xcomp)[i];
    double y = (*ycomp)[i];
    double z = (*zcomp)[i];
    double norm = sqrt(x * x + y * y + z * z);
    // form the unit normal vector
    if (norm > 1e-15)
    {
      normal->ReplaceMyValue(i, 0, x / norm);
      normal->ReplaceMyValue(i, 1, y / norm);
      normal->ReplaceMyValue(i, 2, z / norm);
    }
  }

  return normal;
}  // ScaTra:: ScaTraTimIntImpl::compute_normal_vectors


/*------------------------------------------------------------------------------------------------*
 | compute null space information associated with global system matrix if applicable   fang 09/17 |
 *------------------------------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::compute_null_space_if_necessary() const
{
  // extract solver parameters
  Teuchos::ParameterList& solverparams = solver_->Params();

  // compute point-based null space information if applicable
  if (Core::UTILS::IntegralValue<bool>(*params_, "NULLSPACE_POINTBASED"))
  {
    // MueLu preconditioner
    if (solverparams.isSublist("MueLu Parameters"))
    {
      // extract and fill parameter list for MueLu preconditioner
      Teuchos::ParameterList& mllist = solverparams.sublist("MueLu Parameters", true);
      mllist.set("PDE equations", 1);
      mllist.set("null space: dimension", 1);
      mllist.set("null space: type", "pre-computed");
      mllist.set("null space: add default vectors", false);

      Teuchos::RCP<Epetra_MultiVector> nullspace =
          Teuchos::rcp(new Epetra_MultiVector(*(discret_->dof_row_map()), 1, true));
      nullspace->PutScalar(1.0);

      mllist.set<Teuchos::RCP<Epetra_MultiVector>>("nullspace", nullspace);
      mllist.set("null space: vectors", nullspace->Values());
      mllist.set("ML validate parameter list", false);
    }

    else
    {
      FOUR_C_THROW(
          "Point-based null space calculation currently only implemented for MueLu "
          "preconditioner!");
    }
  }

  else
    // compute standard, vector-based null space information if applicable
    discret_->compute_null_space_if_necessary(solverparams, true);
}  // ScaTra::ScaTraTimIntImpl::compute_null_space_if_necessary()


/*----------------------------------------------------------------------*
 | evaluate Neumann inflow boundary condition                  vg 03/09 |
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::compute_neumann_inflow(
    Teuchos::RCP<Core::LinAlg::SparseOperator> matrix, Teuchos::RCP<Epetra_Vector> rhs)
{
  // time measurement: evaluate condition 'Neumann inflow'
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'TransportNeumannInflow'");

  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
      "action", ScaTra::BoundaryAction::calc_Neumann_inflow, condparams);

  // add element parameters and vectors according to time-integration scheme
  add_time_integration_specific_vectors();

  // set thermodynamic pressure in case of loma
  add_problem_specific_parameters_and_vectors(condparams);

  std::string condstring("TransportNeumannInflow");
  discret_->evaluate_condition(
      condparams, matrix, Teuchos::null, rhs, Teuchos::null, Teuchos::null, condstring);
}  // ScaTra::ScaTraTimIntImpl::compute_neumann_inflow

/*----------------------------------------------------------------------*
 | evaluate boundary cond. due to convective heat transfer     vg 10/11 |
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::evaluate_convective_heat_transfer(
    Teuchos::RCP<Core::LinAlg::SparseOperator> matrix, Teuchos::RCP<Epetra_Vector> rhs)
{
  // time measurement: evaluate condition 'TransportThermoConvections'
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'TransportThermoConvections'");

  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
      "action", ScaTra::BoundaryAction::calc_convective_heat_transfer, condparams);

  // add element parameters and vectors according to time-integration scheme
  add_time_integration_specific_vectors();

  std::string condstring("TransportThermoConvections");
  discret_->evaluate_condition(
      condparams, matrix, Teuchos::null, rhs, Teuchos::null, Teuchos::null, condstring);
}  // ScaTra::ScaTraTimIntImpl::evaluate_convective_heat_transfer


/*------------------------------------------------------------------------------------------*
 | output performance statistics associated with linear solver into *.csv file   fang 02/17 |
 *------------------------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::output_lin_solver_stats(
    const Core::LinAlg::Solver& solver,  //!< linear solver
    const double& time,                  //!< solver time maximized over all processors
    const int& step,                     //!< time step
    const int& iteration,                //!< Newton-Raphson iteration number
    const int& size                      //!< size of linear system
)
{
  // extract communicator
  const Epetra_Comm& comm = solver.Comm();

  // write performance statistics to file
  if (comm.MyPID() == 0)
  {
    // set file name
    std::string filename(
        Global::Problem::Instance()->OutputControlFile()->file_name() + ".lin_solver_stats.csv");

    // open file in appropriate mode and write header at beginning
    std::ofstream file;
    if (step == 1 and iteration == 1)
    {
      file.open(filename.c_str(), std::fstream::trunc);
      file << "Step,NewtonRaphsonIteration,SystemSize,SolverIterations,SolverTime" << std::endl;
    }
    else
      file.open(filename.c_str(), std::fstream::app);

    // write results
    file << step << "," << iteration << "," << size << "," << solver.getNumIters() << ","
         << std::setprecision(16) << std::scientific << time << std::endl;

    // close file
    file.close();
  }
}  // ScaTra::ScaTraTimIntImpl::output_lin_solver_stats


/*---------------------------------------------------------------------------------------------*
 | output performance statistics associated with nonlinear solver into *.csv file   fang 04/15 |
 *---------------------------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::output_nonlin_solver_stats(
    const int& iterations,   //!< iteration count of nonlinear solver
    const double& time,      //!< solver time maximized over all processors
    const int& step,         //!< time step
    const Epetra_Comm& comm  //!< communicator
)
{
  // write performance statistics to file
  if (comm.MyPID() == 0)
  {
    // set file name
    std::string filename(
        Global::Problem::Instance()->OutputControlFile()->file_name() + ".nonlin_solver_stats.csv");

    // open file in appropriate mode and write header at beginning
    std::ofstream file;
    if (step == 1)
    {
      file.open(filename.c_str(), std::fstream::trunc);
      file << "Step,NonlinSolverIterations,NonlinSolverTime" << std::endl;
    }
    else
      file.open(filename.c_str(), std::fstream::app);

    // write results
    file << step << "," << iterations << "," << std::setprecision(16) << std::scientific << time
         << std::endl;

    // close file
    file.close();
  }
}  // ScaTra::ScaTraTimIntImpl::output_nonlin_solver_stats


/*----------------------------------------------------------------------*
 | write state vectors to Gmsh postprocessing files        henke   12/09|
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::output_to_gmsh(const int step, const double time) const
{
  // turn on/off screen output for writing process of Gmsh postprocessing file
  const bool screen_out = true;

  // create Gmsh postprocessing file
  const std::string filename =
      Core::IO::Gmsh::GetNewFileNameAndDeleteOldFiles("solution_field_scalar",
          DiscWriter()->output()->file_name(), step, 500, screen_out, discret_->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  //  {
  //    // add 'View' to Gmsh postprocessing file
  //    gmshfilecontent << "View \" " << "Phin \" {" << std::endl;
  //    // draw scalar field 'Phindtp' for every element
  //    Core::IO::Gmsh::ScalarFieldToGmsh(discret_,phin_,gmshfilecontent);
  //    gmshfilecontent << "};" << std::endl;
  //  }
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "Phinp \" {" << std::endl;
    // draw scalar field 'Phinp' for every element
    Core::IO::Gmsh::ScalarFieldToGmsh(discret_, phinp_, gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }
  //  {
  //    // add 'View' to Gmsh postprocessing file
  //    gmshfilecontent << "View \" " << "Phidtn \" {" << std::endl;
  //    // draw scalar field 'Phindtn' for every element
  //    Core::IO::Gmsh::ScalarFieldToGmsh(discret_,phidtn_,gmshfilecontent);
  //    gmshfilecontent << "};" << std::endl;
  //  }
  //  {
  //    // add 'View' to Gmsh postprocessing file
  //    gmshfilecontent << "View \" " << "Phidtnp \" {" << std::endl;
  //    // draw scalar field 'Phindtp' for every element
  //    Core::IO::Gmsh::ScalarFieldToGmsh(discret_,phidtnp_,gmshfilecontent);
  //    gmshfilecontent << "};" << std::endl;
  //  }
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" "
                    << "Convective Velocity \" {" << std::endl;

    // extract convective velocity from discretization
    Teuchos::RCP<const Epetra_Vector> convel =
        discret_->GetState(NdsVel(), "convective velocity field");
    if (convel == Teuchos::null)
      FOUR_C_THROW("Cannot extract convective velocity field from discretization");

    // draw vector field 'Convective Velocity' for every element
    Core::IO::Gmsh::VectorFieldDofBasedToGmsh(discret_, convel, gmshfilecontent, NdsVel());
    gmshfilecontent << "};" << std::endl;
  }
  gmshfilecontent.close();
  if (screen_out) std::cout << " done" << std::endl;
}  // ScaTraTimIntImpl::output_to_gmsh


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::write_restart() const
{
  // output restart information associated with mesh tying strategy
  strategy_->write_restart();
}


/*----------------------------------------------------------------------*
 |  write mass / heat flux vector to BINIO                   gjb   08/08|
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::output_flux(Teuchos::RCP<Epetra_MultiVector> flux,  //!< flux vector
    const std::string& fluxtype  //!< flux type ("domain" or "boundary")
)
{
  // safety checks
  if (flux == Teuchos::null)
    FOUR_C_THROW("Null pointer for flux vector output. Output() called before Update() ??");
  if (fluxtype != "domain" and fluxtype != "boundary")
    FOUR_C_THROW("Unknown flux type. Must be either 'domain' or 'boundary'!");

  // WORK-AROUND FOR NURBS DISCRETIZATIONS
  // using noderowmap is problematic. Thus, we do not add normal vectors
  // to the scalar result field (scalar information is enough anyway)
  auto* nurbsdis = dynamic_cast<Core::FE::Nurbs::NurbsDiscretization*>(&(*discret_));
  if (nurbsdis != nullptr)
  {
    Teuchos::RCP<Epetra_Vector> normalflux = Teuchos::rcp(((*flux)(0)), false);
    output_->write_vector("normalflux", normalflux, Core::IO::dofvector);
    return;  // leave here
  }

  // post_ensight does not support multivectors based on the dofmap
  // for now, I create single vectors that can be handled by the filter

  // get the noderowmap
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  Teuchos::RCP<Epetra_MultiVector> fluxk =
      Teuchos::rcp(new Epetra_MultiVector(*noderowmap, 3, true));
  for (int writefluxid : *writefluxids_)
  {
    std::ostringstream temp;
    temp << writefluxid;
    std::string name = "flux_" + fluxtype + "_phi_" + temp.str();
    for (int i = 0; i < fluxk->MyLength(); ++i)
    {
      Core::Nodes::Node* actnode = discret_->lRowNode(i);
      int dofgid = discret_->Dof(0, actnode, writefluxid - 1);
      // get value for each component of flux vector
      double xvalue = ((*flux)[0])[(flux->Map()).LID(dofgid)];
      double yvalue = ((*flux)[1])[(flux->Map()).LID(dofgid)];
      double zvalue = ((*flux)[2])[(flux->Map()).LID(dofgid)];
      // care for the slave nodes of rotationally symm. periodic boundary conditions
      double rotangle(0.0);  // already converted to radians
      bool havetorotate = FLD::IsSlaveNodeOfRotSymPBC(actnode, rotangle);
      if (havetorotate)
      {
        double xvalue_rot = (xvalue * cos(rotangle)) - (yvalue * sin(rotangle));
        double yvalue_rot = (xvalue * sin(rotangle)) + (yvalue * (cos(rotangle)));
        xvalue = xvalue_rot;
        yvalue = yvalue_rot;
      }
      // insert values
      int err = fluxk->ReplaceMyValue(i, 0, xvalue);
      err += fluxk->ReplaceMyValue(i, 1, yvalue);
      err += fluxk->ReplaceMyValue(i, 2, zvalue);
      if (err != 0) FOUR_C_THROW("Detected error in ReplaceMyValue");
    }
    output_->write_vector(name, fluxk, Core::IO::nodevector);
  }
}  // ScaTraTimIntImpl::OutputFlux


// TODO: SCATRA_ELE_CLEANING: BIOFILM
// This action is not called at element level!!
// Can it be integrated in CalcFlux_domain?
/*----------------------------------------------------------------------*
 |  output of integral reaction                               mc   03/13|
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::OutputIntegrReac(const int num)
{
  if (outintegrreac_)
  {
    if (!discret_->Filled()) FOUR_C_THROW("fill_complete() was not called");
    if (!discret_->HaveDofs()) FOUR_C_THROW("assign_degrees_of_freedom() was not called");

    // set scalar values needed by elements
    discret_->set_state("phinp", phinp_);
    // set action for elements
    Teuchos::ParameterList eleparams;
    Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
        "action", ScaTra::Action::calc_integr_reaction, eleparams);
    Teuchos::RCP<std::vector<double>> myreacnp =
        Teuchos::rcp(new std::vector<double>(NumScal(), 0.0));
    eleparams.set<Teuchos::RCP<std::vector<double>>>("local reaction integral", myreacnp);

    discret_->Evaluate(eleparams);

    myreacnp = eleparams.get<Teuchos::RCP<std::vector<double>>>("local reaction integral");
    // global integral of reaction terms
    std::vector<double> intreacterm(NumScal(), 0.0);
    for (int k = 0; k < NumScal(); ++k)
      phinp_->Map().Comm().SumAll(&((*myreacnp)[k]), &intreacterm[k], 1);

    // print out values
    if (myrank_ == 0)
    {
      for (int k = 0; k < NumScal(); k++)
      {
        std::cout << "Total reaction (r_" << k << "): " << std::setprecision(9) << intreacterm[k]
                  << std::endl;
      }
    }

    // print out results to file as well
    if (myrank_ == 0)
    {
      std::stringstream number;
      number << num;
      const std::string fname =
          problem_->OutputControlFile()->file_name() + number.str() + ".integrreacvalues.txt";

      std::ofstream f;
      if (Step() <= 1)
      {
        f.open(fname.c_str(), std::fstream::trunc);
        f << "#| Step | Time ";
        for (int k = 0; k < NumScal(); k++)
        {
          f << "| Total reaction (r_" << k << ") ";
        }
        f << "\n";
      }
      else
        f.open(fname.c_str(), std::fstream::ate | std::fstream::app);

      f << Step() << " " << Time() << " ";
      for (int k = 0; k < NumScal(); k++)
      {
        f << std::setprecision(9) << intreacterm[k] << " ";
      }
      f << "\n";
      f.flush();
      f.close();
    }

  }  // if(outintegrreac_)
}  // ScaTra::ScaTraTimIntImpl::OutputIntegrReac


/*==========================================================================*/
// ELCH
/*==========================================================================*/



/*----------------------------------------------------------------------*
 |  prepare AVM3-based scale separation                        vg 10/08 |
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::av_m3_preparation()
{
  // time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:            + avm3");

  // create normalized all-scale subgrid-diffusivity matrix
  sysmat_sd_->Zero();

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // action for elements, time factor and stationary flag
  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::calc_subgrid_diffusivity_matrix, eleparams);

  // add element parameters according to time-integration scheme
  add_time_integration_specific_vectors();

  // call loop over elements
  discret_->Evaluate(eleparams, sysmat_sd_, residual_);

  // finalize the normalized all-scale subgrid-diffusivity matrix
  sysmat_sd_->Complete();

  // apply DBC to normalized all-scale subgrid-diffusivity matrix
  Core::LinAlg::apply_dirichlet_to_system(
      *sysmat_sd_, *phinp_, *residual_, *phinp_, *(dbcmaps_->CondMap()));

  // get normalized fine-scale subgrid-diffusivity matrix
  {
    // this is important to have!!!
    // MLAPI::Init() without arguments uses internally MPI_COMM_WOLRD
    MLAPI::Init();

    // extract the ML parameters
    Teuchos::ParameterList& mlparams = solver_->Params().sublist("ML Parameters");
    // remark: we create a new solver with ML preconditioner here, since this allows for also using
    // other solver setups to solve the system of equations get the solver number used form the
    // multifractal subgrid-scale model parameter list
    const int scale_sep_solvernumber =
        extraparams_->sublist("MULTIFRACTAL SUBGRID SCALES").get<int>("ML_SOLVER");
    if (scale_sep_solvernumber != (-1))  // create a dummy solver
    {
      Teuchos::RCP<Core::LinAlg::Solver> solver =
          Teuchos::rcp(new Core::LinAlg::Solver(problem_->SolverParams(scale_sep_solvernumber),
              discret_->Comm(), Global::Problem::Instance()->solver_params_callback(),
              Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
                  Global::Problem::Instance()->IOParams(), "VERBOSITY")));
      // compute the null space,
      discret_->compute_null_space_if_necessary(solver->Params(), true);
      // and, finally, extract the ML parameters
      mlparams = solver->Params().sublist("ML Parameters");
    }

    // get toggle vector for Dirchlet boundary conditions
    const Teuchos::RCP<const Epetra_Vector> dbct = dirichlet_toggle();

    // get nullspace parameters
    double* nullspace = mlparams.get("null space: vectors", (double*)nullptr);
    if (!nullspace) FOUR_C_THROW("No nullspace supplied in parameter list");
    int nsdim = mlparams.get("null space: dimension", 1);

    // modify nullspace to ensure that DBC are fully taken into account
    if (nullspace)
    {
      const int length = sysmat_sd_->OperatorRangeMap().NumMyElements();
      for (int i = 0; i < nsdim; ++i)
        for (int j = 0; j < length; ++j)
          if ((*dbct)[j] != 0.0) nullspace[i * length + j] = 0.0;
    }

    // get plain aggregation Ptent
    Teuchos::RCP<Epetra_CrsMatrix> crsPtent;
    MLAPI::GetPtent(*sysmat_sd_->EpetraMatrix(), mlparams, nullspace, crsPtent);
    Core::LinAlg::SparseMatrix Ptent(crsPtent, Core::LinAlg::View);

    // compute scale-separation matrix: S = I - Ptent*Ptent^T
    Sep_ = Core::LinAlg::Multiply(Ptent, false, Ptent, true);
    Sep_->Scale(-1.0);
    Teuchos::RCP<Epetra_Vector> tmp = Core::LinAlg::CreateVector(Sep_->RowMap(), false);
    tmp->PutScalar(1.0);
    Teuchos::RCP<Epetra_Vector> diag = Core::LinAlg::CreateVector(Sep_->RowMap(), false);
    Sep_->ExtractDiagonalCopy(*diag);
    diag->Update(1.0, *tmp, 1.0);
    Sep_->replace_diagonal_values(*diag);

    // complete scale-separation matrix and check maps
    Sep_->Complete(Sep_->DomainMap(), Sep_->RangeMap());
    if (!Sep_->RowMap().SameAs(sysmat_sd_->RowMap())) FOUR_C_THROW("rowmap not equal");
    if (!Sep_->RangeMap().SameAs(sysmat_sd_->RangeMap())) FOUR_C_THROW("rangemap not equal");
    if (!Sep_->DomainMap().SameAs(sysmat_sd_->DomainMap())) FOUR_C_THROW("domainmap not equal");

    // precomputation of unscaled diffusivity matrix:
    // either two-sided S^T*M*S: multiply M by S from left- and right-hand side
    // or one-sided M*S: only multiply M by S from left-hand side
    if (not incremental_)
    {
      Mnsv_ = Core::LinAlg::Multiply(*sysmat_sd_, false, *Sep_, false);
      // Mnsv_ = Core::LinAlg::Multiply(*Sep_,true,*Mnsv_,false);
    }
  }
}  // ScaTraTimIntImpl::av_m3_preparation

/*----------------------------------------------------------------------*
 |  scaling of AVM3-based subgrid-diffusivity matrix           vg 10/08 |
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::av_m3_scaling(Teuchos::ParameterList& eleparams)
{
  // time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:            + avm3");

  // some necessary definitions
  int ierr;
  double* sgvsqrt = nullptr;
  int length = subgrdiff_->MyLength();

  // square-root of subgrid-viscosity-scaling vector for left and right scaling
  sgvsqrt = (double*)subgrdiff_->Values();
  for (int i = 0; i < length; ++i)
  {
    sgvsqrt[i] = sqrt(sgvsqrt[i]);
    int err = subgrdiff_->ReplaceMyValues(1, &sgvsqrt[i], &i);
    if (err != 0) FOUR_C_THROW("index not found");
  }

  // get unscaled S^T*M*S from Sep
  sysmat_sd_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*Mnsv_));

  // left and right scaling of normalized fine-scale subgrid-viscosity matrix
  ierr = sysmat_sd_->LeftScale(*subgrdiff_);
  if (ierr) FOUR_C_THROW("Epetra_CrsMatrix::LeftScale returned err=%d", ierr);
  ierr = sysmat_sd_->RightScale(*subgrdiff_);
  if (ierr) FOUR_C_THROW("Epetra_CrsMatrix::RightScale returned err=%d", ierr);

  // add the subgrid-viscosity-scaled fine-scale matrix to obtain complete matrix
  Teuchos::RCP<Core::LinAlg::SparseMatrix> sysmat = SystemMatrix();
  sysmat->Add(*sysmat_sd_, false, 1.0, 1.0);

  // set subgrid-diffusivity vector to zero after scaling procedure
  subgrdiff_->PutScalar(0.0);
}  // ScaTraTimIntImpl::AVM3Scaling

/*----------------------------------------------------------------------*
 | construct toggle vector for Dirichlet dofs                  gjb 11/08|
 | assures backward compatibility for avm3 solver; should go away once  |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ScaTra::ScaTraTimIntImpl::dirichlet_toggle()
{
  if (dbcmaps_ == Teuchos::null) FOUR_C_THROW("Dirichlet map has not been allocated");
  Teuchos::RCP<Epetra_Vector> dirichones =
      Core::LinAlg::CreateVector(*(dbcmaps_->CondMap()), false);
  dirichones->PutScalar(1.0);
  Teuchos::RCP<Epetra_Vector> dirichtoggle =
      Core::LinAlg::CreateVector(*(discret_->dof_row_map()), true);
  dbcmaps_->InsertCondVector(dirichones, dirichtoggle);
  return dirichtoggle;
}  // ScaTraTimIntImpl::dirichlet_toggle


/*========================================================================*/
// turbulence and related
/*========================================================================*/


/*----------------------------------------------------------------------*
 | provide access to the dynamic Smagorinsky filter     rasthofer 08/12 |
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::AccessDynSmagFilter(Teuchos::RCP<FLD::DynSmagFilter> dynSmag)
{
  DynSmag_ = Teuchos::rcp(dynSmag.get(), false);

  // access to the dynamic Smagorinsky class is provided
  // by the fluid scatra coupling algorithm
  // therefore, we only have to add some scatra specific parameters here
  if (turbmodel_ == Inpar::FLUID::dynamic_smagorinsky)
  {
    if (myrank_ == 0)
    {
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
      std::cout << "Dynamic Smagorinsky model: provided access for ScaTra       " << std::endl;
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    }
    DynSmag_->AddScatra(discret_);
  }
}

/*----------------------------------------------------------------------*
 | provide access to the dynamic Vreman class               krank 09/13 |
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::AccessVreman(Teuchos::RCP<FLD::Vreman> vrem)
{
  Vrem_ = Teuchos::rcp(vrem.get(), false);

  // access to the dynamic Vreman class is provided
  // by the fluid scatra coupling algorithm
  // therefore, we only have to add some scatra specific parameters here
  if (turbmodel_ == Inpar::FLUID::dynamic_vreman)
  {
    if (myrank_ == 0)
    {
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
      std::cout << "Dynamic Vreman model: provided access for ScaTra            " << std::endl;
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    }
    Vrem_->AddScatra(discret_);
  }
}


/*----------------------------------------------------------------------*
 | read restart data                                         fang 01/17 |
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::read_restart(
    const int step, Teuchos::RCP<Core::IO::InputControl> input)
{
  // read restart data associated with meshtying strategy
  strategy_->read_restart(step, input);
}  // ScaTra::ScaTraTimIntImpl::read_restart()


/*-------------------------------------------------------------------------*
 | calculate mean CsgsB to estimate CsgsD                                  |
 | for multifractal subgrid-scale model                    rasthofer 08/12 |
 *-------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::recompute_mean_csgs_b()
{
  if (Core::UTILS::IntegralValue<int>(
          extraparams_->sublist("MULTIFRACTAL SUBGRID SCALES"), "ADAPT_CSGS_PHI"))
  {
    // mean Cai
    double meanCai = 0.0;

    // variables required for calculation
    // local sums
    double local_sumCai = 0.0;
    double local_sumVol = 0.0;
    // global sums
    double global_sumCai = 0.0;
    double global_sumVol = 0.0;

    // define element matrices and vectors --- dummies
    Core::LinAlg::SerialDenseMatrix emat1;
    Core::LinAlg::SerialDenseMatrix emat2;
    Core::LinAlg::SerialDenseVector evec1;
    Core::LinAlg::SerialDenseVector evec2;
    Core::LinAlg::SerialDenseVector evec3;

    // generate a parameterlist for communication and control
    Teuchos::ParameterList myparams;

    // action for elements
    Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
        "action", ScaTra::Action::calc_mean_Cai, myparams);

    // add element parameters according to time-integration scheme
    add_time_integration_specific_vectors();

    // set thermodynamic pressure in case of loma
    add_problem_specific_parameters_and_vectors(myparams);

    // loop all elements on this proc (excluding ghosted ones)
    for (int nele = 0; nele < discret_->NumMyRowElements(); ++nele)
    {
      // get the element
      Core::Elements::Element* ele = discret_->lRowElement(nele);

      // get element location vector, dirichlet flags and ownerships
      Core::Elements::Element::LocationArray la(discret_->NumDofSets());
      ele->LocationVector(*discret_, la, false);

      // call the element evaluate method to integrate functions
      int err = ele->Evaluate(myparams, *discret_, la, emat1, emat2, evec1, evec2, evec2);
      if (err) FOUR_C_THROW("Proc %d: Element %d returned err=%d", myrank_, ele->Id(), err);

      // get contributions of this element and add it up
      local_sumCai += myparams.get<double>("Cai_int");
      local_sumVol += myparams.get<double>("ele_vol");
    }

    // gather contibutions of all procs
    discret_->Comm().SumAll(&local_sumCai, &global_sumCai, 1);
    discret_->Comm().SumAll(&local_sumVol, &global_sumVol, 1);

    // calculate mean Cai
    meanCai = global_sumCai / global_sumVol;

    // std::cout << "Proc:  " << myrank_ << "  local vol and Cai   "
    //<< local_sumVol << "   " << local_sumCai << "  global vol and Cai   "
    //<< global_sumVol << "   " << global_sumCai << "  mean   " << meanCai << std::endl;

    if (myrank_ == 0)
    {
      std::cout << "\n+----------------------------------------------------------------------------"
                   "----------------+"
                << std::endl;
      std::cout << "Multifractal subgrid scales: adaption of CsgsD from near-wall limit of CsgsB:  "
                << std::setprecision(8) << meanCai << std::endl;
      std::cout << "+------------------------------------------------------------------------------"
                   "--------------+\n"
                << std::endl;
    }

    // set meanCai via pre-evaluate call
    Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
        "action", ScaTra::Action::set_mean_Cai, myparams);
    myparams.set<double>("meanCai", meanCai);
    // call standard loop over elements
    discret_->Evaluate(
        myparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  }
}


/*----------------------------------------------------------------------*
 | calculate intermediate solution                       rasthofer 05/13|
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::calc_intermediate_solution()
{
  if (special_flow_ == "scatra_forced_homogeneous_isotropic_turbulence" and
      extraparams_->sublist("TURBULENCE MODEL").get<std::string>("SCALAR_FORCING") ==
          "isotropic" and
      Core::UTILS::IntegralValue<Inpar::FLUID::ForcingType>(
          extraparams_->sublist("TURBULENCE MODEL"), "FORCING_TYPE") ==
          Inpar::FLUID::linear_compensation_from_intermediate_spectrum)
  {
    bool activate = true;

    if (activate)
    {
      if (homisoturb_forcing_ == Teuchos::null) FOUR_C_THROW("Forcing expected!");

      if (myrank_ == 0)
      {
        std::cout << "+------------------------------------------------------------------------+\n";
        std::cout << "|     calculate intermediate solution\n";
        std::cout << "|" << std::endl;
      }

      // turn off forcing in nonlinear_solve()
      homisoturb_forcing_->ActivateForcing(false);

      // temporary store velnp_ since it will be modified in nonlinear_solve()
      const Epetra_Map* dofrowmap = discret_->dof_row_map();
      Teuchos::RCP<Epetra_Vector> tmp = Core::LinAlg::CreateVector(*dofrowmap, true);
      tmp->Update(1.0, *phinp_, 0.0);

      // compute intermediate solution without forcing
      forcing_->PutScalar(0.0);  // just to be sure
      nonlinear_solve();

      // calculate required forcing
      homisoturb_forcing_->CalculateForcing(step_);

      // reset velnp_
      phinp_->Update(1.0, *tmp, 0.0);

      // recompute intermediate values, since they have been likewise overwritten
      // only for gen.-alpha
      compute_intermediate_values();

      homisoturb_forcing_->ActivateForcing(true);

      if (myrank_ == 0)
      {
        std::cout << "|\n";
        std::cout << "+------------------------------------------------------------------------+\n";
        std::cout << "+------------------------------------------------------------------------+\n";
        std::cout << "|" << std::endl;
      }
    }
    else
      // set force to zero
      forcing_->PutScalar(0.0);
  }
}


/*==========================================================================*/
// Biofilm related
/*==========================================================================*/

/*----------------------------------------------------------------------*/
/* set scatra-fluid displacement vector due to biofilm growth          */
void ScaTra::ScaTraTimIntImpl::SetScFldGrDisp(
    Teuchos::RCP<Epetra_MultiVector> scatra_fluid_growth_disp)
{
  scfldgrdisp_ = scatra_fluid_growth_disp;
}

/*----------------------------------------------------------------------*/
/* set scatra-structure displacement vector due to biofilm growth          */
void ScaTra::ScaTraTimIntImpl::SetScStrGrDisp(
    Teuchos::RCP<Epetra_MultiVector> scatra_struct_growth_disp)
{
  scstrgrdisp_ = scatra_struct_growth_disp;
}

/*----------------------------------------------------------------------*
 | Calculate the reconstructed nodal gradient of phi        winter 04/17|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> ScaTra::ScaTraTimIntImpl::compute_superconvergent_patch_recovery(
    Teuchos::RCP<const Epetra_Vector> state, const std::string& statename, const int numvec,
    Teuchos::ParameterList& eleparams, const int dim)
{
  // Warning, this is only tested so far for 1 scalar field!!!

  // dependent on the desired projection, just remove this line
  if (not state->Map().SameAs(*discret_->dof_row_map()))
    FOUR_C_THROW("input map is not a dof row map of the fluid");

  // set given state for element evaluation
  discret_->set_state(statename, state);

  switch (dim)
  {
    case 3:
      return Core::FE::compute_superconvergent_patch_recovery<3>(
          *discret_, *state, statename, numvec, eleparams);
      break;
    case 2:
      return Core::FE::compute_superconvergent_patch_recovery<2>(
          *discret_, *state, statename, numvec, eleparams);
      break;
    case 1:
      return Core::FE::compute_superconvergent_patch_recovery<1>(
          *discret_, *state, statename, numvec, eleparams);
      break;
    default:
      FOUR_C_THROW("only 1/2/3D implementation available for superconvergent patch recovery");
      break;
  }

  return Teuchos::null;
}

/*--------------------------------------------------------------------------------------------------------------------*
 | convergence check (only for two-way coupled problems, e.g., low-Mach-number flow, ...) |
 *--------------------------------------------------------------------------------------------------------------------*/
bool ScaTra::ScaTraTimIntImpl::convergence_check(int itnum, int itmax, const double ittol)
{
  // declare bool variable for potentially stopping nonlinear iteration
  bool stopnonliniter = false;

  // declare L2-norm of (two) residual, incremental scalar and scalar
  double res1norm_L2(0.0);
  double res2norm_L2(0.0);
  double phi1incnorm_L2(0.0);
  double phi2incnorm_L2(0.0);
  double phi1norm_L2(0.0);
  double phi2norm_L2(0.0);

  // distinguish whether one or two scalars are considered
  if (NumScal() == 2)
  {
    Teuchos::RCP<Epetra_Vector> vec1 = splitter_->ExtractOtherVector(residual_);
    Teuchos::RCP<Epetra_Vector> vec2 = splitter_->ExtractCondVector(residual_);
    vec1->Norm2(&res1norm_L2);
    vec2->Norm2(&res2norm_L2);

    vec1 = splitter_->ExtractOtherVector(increment_);
    vec2 = splitter_->ExtractCondVector(increment_);
    vec1->Norm2(&phi1incnorm_L2);
    vec2->Norm2(&phi2incnorm_L2);

    vec1 = splitter_->ExtractOtherVector(phinp_);
    vec2 = splitter_->ExtractCondVector(phinp_);
    vec1->Norm2(&phi1norm_L2);
    vec2->Norm2(&phi2norm_L2);

    // check for any INF's and NaN's
    if (std::isnan(res1norm_L2) or std::isnan(phi1incnorm_L2) or std::isnan(phi1norm_L2) or
        std::isnan(res2norm_L2) or std::isnan(phi2incnorm_L2) or std::isnan(phi2norm_L2))
      FOUR_C_THROW("At least one of the calculated vector norms is NaN.");

    if (std::isinf(res1norm_L2) or std::isinf(phi1incnorm_L2) or std::isinf(phi1norm_L2) or
        std::isinf(res2norm_L2) or std::isinf(phi2incnorm_L2) or std::isinf(phi2norm_L2))
      FOUR_C_THROW("At least one of the calculated vector norms is INF.");

    // for scalar norm being (close to) zero, set to one
    if (phi1norm_L2 < 1e-5) phi1norm_L2 = 1.0;
    if (phi2norm_L2 < 1e-5) phi2norm_L2 = 1.0;

    if (myrank_ == 0)
    {
      printf(
          "+------------+-------------------+---------------+---------------+---------------+------"
          "---------+\n");
      printf(
          "|- step/max -|- tol      [norm] -|- residual1   -|- scalar1-inc -|- residual2   -|- "
          "scalar2-inc -|\n");
      printf("|  %3d/%3d   | %10.3E[L_2 ]  |  %10.3E   |  %10.3E   |  %10.3E   |  %10.3E   |",
          itnum, itmax, ittol, res1norm_L2, phi1incnorm_L2 / phi1norm_L2, res2norm_L2,
          phi2incnorm_L2 / phi2norm_L2);
      printf("\n");
      printf(
          "+------------+-------------------+---------------+---------------+---------------+------"
          "---------+\n");
    }

    if ((res1norm_L2 <= ittol) and (phi1incnorm_L2 / phi1norm_L2 <= ittol) and
        (res2norm_L2 <= ittol) and (phi2incnorm_L2 / phi2norm_L2 <= ittol))
      stopnonliniter = true;

    // warn if itemax is reached without convergence, but proceed to next timestep
    if ((itnum == itmax) and ((res1norm_L2 > ittol) or (phi1incnorm_L2 / phi1norm_L2 > ittol) or
                                 (res2norm_L2 > ittol) or (phi2incnorm_L2 / phi2norm_L2 > ittol)))
    {
      stopnonliniter = true;
      if (myrank_ == 0)
      {
        printf("|            >>>>>> not converged in itemax steps!             |\n");
        printf("+--------------------------------------------------------------+\n");
      }
    }
  }
  else if (NumScal() == 1)
  {
    // compute L2-norm of residual, incremental scalar and scalar
    residual_->Norm2(&res1norm_L2);
    increment_->Norm2(&phi1incnorm_L2);
    phinp_->Norm2(&phi1norm_L2);

    // check for any INF's and NaN's
    if (std::isnan(res1norm_L2) or std::isnan(phi1incnorm_L2) or std::isnan(phi1norm_L2))
      FOUR_C_THROW("At least one of the calculated vector norms is NaN.");

    if (std::isinf(res1norm_L2) or std::isinf(phi1incnorm_L2) or std::isinf(phi1norm_L2))
      FOUR_C_THROW("At least one of the calculated vector norms is INF.");

    // for scalar norm being (close to) zero, set to one
    if (phi1norm_L2 < 1e-5) phi1norm_L2 = 1.0;

    if (myrank_ == 0)
    {
      printf("+------------+-------------------+--------------+--------------+\n");
      printf("|- step/max -|- tol      [norm] -|- residual   -|- scalar-inc -|\n");
      printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E   |", itnum, itmax, ittol,
          res1norm_L2, phi1incnorm_L2 / phi1norm_L2);
      printf("\n");
      printf("+------------+-------------------+--------------+--------------+\n");
    }

    if ((res1norm_L2 <= ittol) and (phi1incnorm_L2 / phi1norm_L2 <= ittol)) stopnonliniter = true;

    // warn if itemax is reached without convergence, but proceed to next timestep
    if ((itnum == itmax) and ((res1norm_L2 > ittol) or (phi1incnorm_L2 / phi1norm_L2 > ittol)))
    {
      stopnonliniter = true;
      if (myrank_ == 0)
      {
        printf("|            >>>>>> not converged in itemax steps!             |\n");
        printf("+--------------------------------------------------------------+\n");
      }
    }
  }
  else
    FOUR_C_THROW(
        "ScaTra convergence check for number of scalars other than one or two not yet supported!");

  return stopnonliniter;
}  // ScaTra::ScaTraTimIntImplicit::convergence_check


/*----------------------------------------------------------------------------------------------*
 | finite difference check for scalar transport system matrix (for debugging only)   fang 10/14 |
 *----------------------------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::fd_check()
{
  // initial screen output
  if (myrank_ == 0)
    std::cout << std::endl << "FINITE DIFFERENCE CHECK FOR SCATRA SYSTEM MATRIX" << std::endl;

  // make a copy of state variables to undo perturbations later
  Teuchos::RCP<Epetra_Vector> phinp_original = Teuchos::rcp(new Epetra_Vector(*phinp_));

  // make a copy of system matrix as Epetra_CrsMatrix
  Teuchos::RCP<Epetra_CrsMatrix> sysmat_original = Teuchos::null;
  if (Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(sysmat_) != Teuchos::null)
  {
    sysmat_original = (new Core::LinAlg::SparseMatrix(
                           *(Teuchos::rcp_static_cast<Core::LinAlg::SparseMatrix>(sysmat_))))
                          ->EpetraMatrix();
  }
  else if (Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(sysmat_) != Teuchos::null)
  {
    sysmat_original =
        (new Core::LinAlg::SparseMatrix(
             *(Teuchos::rcp_static_cast<Core::LinAlg::BlockSparseMatrixBase>(sysmat_)->Merge())))
            ->EpetraMatrix();
  }
  else
    FOUR_C_THROW("Type of system matrix unknown!");
  sysmat_original->FillComplete();

  // make a copy of system right-hand side vector
  Teuchos::RCP<Epetra_Vector> rhs_original = Teuchos::rcp(new Epetra_Vector(*residual_));

  // initialize counter for system matrix entries with failing finite difference check
  int counter(0);

  // initialize tracking variable for maximum absolute and relative errors
  double maxabserr(0.);
  double maxrelerr(0.);

  for (int colgid = 0; colgid <= sysmat_original->ColMap().MaxAllGID(); ++colgid)
  {
    // check whether current column index is a valid global column index and continue loop if not
    int collid(sysmat_original->ColMap().LID(colgid));
    int maxcollid(-1);
    discret_->Comm().MaxAll(&collid, &maxcollid, 1);
    if (maxcollid < 0) continue;

    // fill state vector with original state variables
    phinp_->Update(1., *phinp_original, 0.);

    // impose perturbation
    if (phinp_->Map().MyGID(colgid))
      if (phinp_->SumIntoGlobalValue(colgid, 0, fdcheckeps_))
        FOUR_C_THROW(
            "Perturbation could not be imposed on state vector for finite difference check!");

    // carry perturbation over to state vectors at intermediate time stages if necessary
    compute_intermediate_values();

    // calculate element right-hand side vector for perturbed state
    assemble_mat_and_rhs();

    // Now we compare the difference between the current entries in the system matrix
    // and their finite difference approximations according to
    // entries ?= (residual_perturbed - residual_original) / epsilon

    // Note that the residual_ vector actually denotes the right-hand side of the linear
    // system of equations, i.e., the negative system residual.
    // To account for errors due to numerical cancellation, we additionally consider
    // entries + residual_original / epsilon ?= residual_perturbed / epsilon

    // Note that we still need to evaluate the first comparison as well. For small entries in the
    // system matrix, the second comparison might yield good agreement in spite of the entries being
    // wrong!
    for (int rowlid = 0; rowlid < discret_->dof_row_map()->NumMyElements(); ++rowlid)
    {
      // get global index of current matrix row
      const int rowgid = sysmat_original->RowMap().GID(rowlid);
      if (rowgid < 0) FOUR_C_THROW("Invalid global ID of matrix row!");

      // get relevant entry in current row of original system matrix
      double entry(0.);
      int length = sysmat_original->NumMyEntries(rowlid);
      int numentries;
      std::vector<double> values(length);
      std::vector<int> indices(length);
      sysmat_original->ExtractMyRowCopy(rowlid, length, numentries, values.data(), indices.data());
      for (int ientry = 0; ientry < length; ++ientry)
      {
        if (sysmat_original->ColMap().GID(indices[ientry]) == colgid)
        {
          entry = values[ientry];
          break;
        }
      }

      // finite difference suggestion (first divide by epsilon and then add for better conditioning)
      const double fdval =
          -(*residual_)[rowlid] / fdcheckeps_ + (*rhs_original)[rowlid] / fdcheckeps_;

      // confirm accuracy of first comparison
      if (abs(fdval) > 1.e-17 and abs(fdval) < 1.e-15)
        FOUR_C_THROW("Finite difference check involves values too close to numerical zero!");

      // absolute and relative errors in first comparison
      const double abserr1 = entry - fdval;
      if (abs(abserr1) > maxabserr) maxabserr = abs(abserr1);
      double relerr1(0.);
      if (abs(entry) > 1.e-17)
        relerr1 = abserr1 / abs(entry);
      else if (abs(fdval) > 1.e-17)
        relerr1 = abserr1 / abs(fdval);
      if (abs(relerr1) > maxrelerr) maxrelerr = abs(relerr1);

      // evaluate first comparison
      if (abs(relerr1) > fdchecktol_)
      {
        std::cout << "sysmat[" << rowgid << "," << colgid << "]:  " << entry << "   ";
        std::cout << "finite difference suggestion:  " << fdval << "   ";
        std::cout << "absolute error:  " << abserr1 << "   ";
        std::cout << "relative error:  " << relerr1 << std::endl;

        counter++;
      }

      // first comparison OK
      else
      {
        // left-hand side in second comparison
        const double left = entry - (*rhs_original)[rowlid] / fdcheckeps_;

        // right-hand side in second comparison
        const double right = -(*residual_)[rowlid] / fdcheckeps_;

        // confirm accuracy of second comparison
        if (abs(right) > 1.e-17 and abs(right) < 1.e-15)
          FOUR_C_THROW("Finite difference check involves values too close to numerical zero!");

        // absolute and relative errors in second comparison
        const double abserr2 = left - right;
        if (abs(abserr2) > maxabserr) maxabserr = abs(abserr2);
        double relerr2(0.);
        if (abs(left) > 1.e-17)
          relerr2 = abserr2 / abs(left);
        else if (abs(right) > 1.e-17)
          relerr2 = abserr2 / abs(right);
        if (abs(relerr2) > maxrelerr) maxrelerr = abs(relerr2);

        // evaluate second comparison
        if (abs(relerr2) > fdchecktol_)
        {
          std::cout << "sysmat[" << rowgid << "," << colgid << "]-rhs[" << rowgid
                    << "]/eps:  " << left << "   ";
          std::cout << "-rhs_perturbed[" << rowgid << "]/eps:  " << right << "   ";
          std::cout << "absolute error:  " << abserr2 << "   ";
          std::cout << "relative error:  " << relerr2 << std::endl;

          counter++;
        }
      }
    }
  }

  // communicate tracking variables
  int counterglobal(0);
  discret_->Comm().SumAll(&counter, &counterglobal, 1);
  double maxabserrglobal(0.);
  discret_->Comm().MaxAll(&maxabserr, &maxabserrglobal, 1);
  double maxrelerrglobal(0.);
  discret_->Comm().MaxAll(&maxrelerr, &maxrelerrglobal, 1);

  // final screen output
  if (myrank_ == 0)
  {
    if (counterglobal)
    {
      printf(
          "--> FAILED AS LISTED ABOVE WITH %d CRITICAL MATRIX ENTRIES IN TOTAL\n\n", counterglobal);
      FOUR_C_THROW("Finite difference check failed for scalar transport system matrix!");
    }
    else
    {
      printf(
          "--> PASSED WITH MAXIMUM ABSOLUTE ERROR %+12.5e AND MAXIMUM RELATIVE ERROR %+12.5e\n\n",
          maxabserrglobal, maxrelerrglobal);
    }
  }

  // undo perturbations of state variables
  phinp_->Update(1., *phinp_original, 0.);
  compute_intermediate_values();

  // recompute system matrix and right-hand side vector based on original state variables
  assemble_mat_and_rhs();
}

/*----------------------------------------------------------------------------------*
 | calculation of relative error with reference to analytical solution   fang 11/16 |
 *----------------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::evaluate_error_compared_to_analytical_sol()
{
  switch (calcerror_)
  {
    case Inpar::ScaTra::calcerror_byfunction:
    case Inpar::ScaTra::calcerror_spherediffusion:
    {
      // create the parameters for the error calculation
      Teuchos::ParameterList eleparams;
      Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
          "action", ScaTra::Action::calc_error, eleparams);
      eleparams.set<int>("calcerrorflag", calcerror_);

      if (calcerror_ == Inpar::ScaTra::calcerror_byfunction)
      {
        const int errorfunctnumber = params_->get<int>("CALCERRORNO");
        if (errorfunctnumber < 1)
          FOUR_C_THROW("invalid value of paramter CALCERRORNO for error function evaluation!");

        eleparams.set<int>("error function number", errorfunctnumber);
      }

      // set vector values needed by elements
      discret_->set_state("phinp", phinp_);

      // get (squared) error values
      Teuchos::RCP<Core::LinAlg::SerialDenseVector> errors =
          Teuchos::rcp(new Core::LinAlg::SerialDenseVector(4 * NumDofPerNode()));
      discret_->EvaluateScalars(eleparams, errors);

      for (int k = 0; k < NumDofPerNode(); ++k)
      {
        if (std::abs((*errors)[k * 4 + 2]) > 1e-14)
          (*relerrors_)[k * 2] = sqrt((*errors)[k * 4]) / sqrt((*errors)[k * 4 + 2]);
        else
          FOUR_C_THROW("Can't compute relative L2 error due to numerical roundoff sensitivity!");
        if (std::abs((*errors)[k * 4 + 3]) > 1e-14)
          (*relerrors_)[k * 2 + 1] = sqrt((*errors)[k * 4 + 1]) / sqrt((*errors)[k * 4 + 3]);
        else
          FOUR_C_THROW("Can't compute relative H1 error due to numerical roundoff sensitivity!");

        if (myrank_ == 0)
        {
          // print last error in a separate file
          //        if ((step_==stepmax_) or (time_==maxtime_))// write results to file
          //        {
          //          std::ostringstream temp;
          //          temp << k;
          //          const std::string simulation = problem_->OutputControlFile()->file_name();
          //          const std::string fname = simulation+"_c"+temp.str()+".relerror";
          //
          //          std::ofstream f;
          //          f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
          //          f << "#| " << simulation << "\n";
          //          f << "#| Step | Time | rel. L2-error scalar | rel. H1-error scalar  |\n";
          //          f << step_ << " " << time_ << " " << (*relerrors_)[k*2] << " " <<
          //          (*relerrors_)[k*2+1] << "\n"; f.flush(); f.close();
          //        }

          std::ostringstream temp;
          temp << k;
          const std::string simulation = problem_->OutputControlFile()->file_name();
          const std::string fname = simulation + "_c" + temp.str() + "_time.relerror";
          std::ofstream f;

          // create new error file and write initial error
          if (step_ == 0)
          {
            f.open(fname.c_str());
            f << "| Step | Time | rel. L2-error | rel. H1-error |" << std::endl;
          }

          // append error of the last time step to the error file
          else
            f.open(fname.c_str(), std::fstream::ate | std::fstream::app);

          f << step_ << " " << time_ << " " << std::setprecision(6) << (*relerrors_)[k * 2] << " "
            << (*relerrors_)[k * 2 + 1] << std::endl;

          f.flush();
          f.close();
        }
      }

      break;
    }

    case Inpar::ScaTra::calcerror_bycondition:
    {
      // extract conditions from discretization
      std::vector<Core::Conditions::Condition*> relerrorconditions;
      discret_->GetCondition("ScatraRelError", relerrorconditions);
      const unsigned ncond = relerrorconditions.size();

      // loop over all conditions
      for (unsigned icond = 0; icond < ncond; ++icond)
      {
        // extract condition ID
        const int condid = relerrorconditions[icond]->parameters().Get<int>("ConditionID");

        // create element parameter list for error calculation
        Teuchos::ParameterList eleparams;
        Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
            "action", ScaTra::Action::calc_error, eleparams);
        eleparams.set<int>("calcerrorflag", Inpar::ScaTra::calcerror_byfunction);
        const int errorfunctnumber = relerrorconditions[icond]->parameters().Get<int>("FunctionID");
        if (errorfunctnumber < 1) FOUR_C_THROW("Invalid function number for error calculation!");
        eleparams.set<int>("error function number", errorfunctnumber);

        // set state vector needed by elements
        discret_->set_state("phinp", phinp_);

        // get (squared) error values
        Teuchos::RCP<Core::LinAlg::SerialDenseVector> errors =
            Teuchos::rcp(new Core::LinAlg::SerialDenseVector(4 * NumDofPerNode()));
        discret_->EvaluateScalars(eleparams, errors, "ScatraRelError", condid);

        // compute index offset
        const int offset = condid * NumDofPerNode() * 2;

        // loop over all degrees of freedom
        for (int k = 0; k < NumDofPerNode(); ++k)
        {
          // compute relative L2 error
          if (std::abs((*errors)[k * 4 + 2]) > 1e-14)
            (*relerrors_)[offset + k * 2] = sqrt((*errors)[k * 4]) / sqrt((*errors)[k * 4 + 2]);
          else
            FOUR_C_THROW("Can't compute relative L2 error due to numerical roundoff sensitivity!");

          // compute relative H1 error
          if (std::abs((*errors)[k * 4 + 3]) > 1e-14)
          {
            (*relerrors_)[offset + k * 2 + 1] =
                sqrt((*errors)[k * 4 + 1]) / sqrt((*errors)[k * 4 + 3]);
          }
          else
            FOUR_C_THROW("Can't compute relative H1 error due to numerical roundoff sensitivity!");
        }
      }

      // print errors to files
      if (myrank_ == 0)
      {
        // loop over all degrees of freedom
        for (int k = 0; k < NumDofPerNode(); ++k)
        {
          // determine name of file associated with current degree of freedom
          std::ostringstream temp;
          temp << k;
          const std::string fname =
              problem_->OutputControlFile()->file_name() + "_dof_" + temp.str() + ".relerror";

          // initialize output file stream
          std::ofstream f;

          // create new error file and write initial error
          if (step_ == 0)
          {
            f.open(fname.c_str());
            f << "| Step | Time |";

            // loop over all conditions
            for (unsigned icond = 0; icond < ncond; ++icond)
            {
              // extract condition ID
              const int condid = relerrorconditions[icond]->parameters().Get<int>("ConditionID");

              // extend headline
              f << " rel. L2-error in domain " << condid << " | rel. H1-error in domain " << condid
                << " |";
            }

            // break line
            f << std::endl;
          }

          // append error of the last time step to the error file
          else
            f.open(fname.c_str(), std::fstream::ate | std::fstream::app);

          // write error
          f << step_ << " " << time_ << std::scientific << std::setprecision(6);

          // loop over all conditions
          for (unsigned icond = 0; icond < ncond; ++icond)
          {
            // extract condition ID
            const int condid = relerrorconditions[icond]->parameters().Get<int>("ConditionID");

            // extend error line
            f << " " << (*relerrors_)[condid * NumDofPerNode() * 2 + k * 2] << " "
              << (*relerrors_)[condid * NumDofPerNode() * 2 + k * 2 + 1];
          }

          // finalize
          f << std::endl;
          f.flush();
          f.close();
        }
      }

      break;
    }

    case Inpar::ScaTra::calcerror_no:
    {
      // do nothing
      break;
    }

    default:
    {
      FOUR_C_THROW("Cannot calculate error. Unknown type of analytical test problem!");
      break;
    }
  }
}  // ScaTra::ScaTraTimIntImpl::evaluate_error_compared_to_analytical_sol


/*------------------------------------------------------------------------------------------------------*
 | do explicit predictor step to obtain better starting value for Newton-Raphson iteration   fang
 01/17 |
 *------------------------------------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::explicit_predictor() const
{
  // predict discrete solution variables associated with meshtying strategy
  strategy_->explicit_predictor();
}


/*--------------------------------------------------------------------------*
 | perform Aitken relaxation                                     fang 08/17 |
 *--------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::perform_aitken_relaxation(
    Epetra_Vector& phinp,  //!< state vector to be relaxed
    const Epetra_Vector&
        phinp_inc_diff  //!< difference between current and previous state vector increments
)
{
  if (solvtype_ == Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro_aitken)
  {
    // compute L2 norm of difference between current and previous increments of macro-scale state
    // vector
    double phinp_inc_diff_L2(0.);
    phinp_inc_diff.Norm2(&phinp_inc_diff_L2);

    // compute dot product between increment of macro-scale state vector and difference between
    // current and previous increments of macro-scale state vector
    double phinp_inc_dot_phinp_inc_diff(0.);
    if (phinp_inc_diff.Dot(*phinp_inc_, &phinp_inc_dot_phinp_inc_diff))
      FOUR_C_THROW("Couldn't compute dot product!");

    // compute Aitken relaxation factor
    if (iternum_outer_ > 1 and phinp_inc_diff_L2 > 1.e-12)
      omega_[0] *= 1 - phinp_inc_dot_phinp_inc_diff / (phinp_inc_diff_L2 * phinp_inc_diff_L2);

    // perform Aitken relaxation
    phinp.Update(omega_[0], *phinp_inc_, 1.);
  }

  else
    FOUR_C_THROW("Invalid Aitken method!");
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::ApplyBCToSystem()
{
  apply_dirichlet_bc(time_, phinp_, Teuchos::null);
  compute_intermediate_values();
  apply_neumann_bc(neumann_loads_);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::evaluate_initial_time_derivative(
    Teuchos::RCP<Core::LinAlg::SparseOperator> matrix, Teuchos::RCP<Epetra_Vector> rhs)
{
  // create and fill parameter list for elements
  Teuchos::ParameterList eleparams;
  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::calc_initial_time_deriv, eleparams);
  add_problem_specific_parameters_and_vectors(eleparams);

  // add state vectors according to time integration scheme
  add_time_integration_specific_vectors();
  // modify global system of equations as explained above
  discret_->Evaluate(eleparams, matrix, rhs);

  // finalize assembly of global mass matrix
  matrix->Complete();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::OutputScalarsStrategyBase::prepare_evaluate(
    const ScaTraTimIntImpl* const scatratimint, Teuchos::ParameterList& eleparams)
{
  const Teuchos::RCP<Core::FE::Discretization>& discret = scatratimint->discret_;

  // add state vector to discretization
  discret->set_state("phinp", scatratimint->phinp_);

  // set action for elements
  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::calc_total_and_mean_scalars, eleparams);
  eleparams.set("inverting", false);
  eleparams.set("calc_grad_phi", output_mean_grad_);
}

/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
void ScaTra::OutputScalarsStrategyBase::print_header_to_screen(const std::string& dis_name)
{
  if (myrank_ == 0)
  {
    // screen output
    std::cout << std::endl;
    std::cout << "Total and mean values of transported scalars from discretization " << dis_name
              << ":" << std::endl;
    std::cout
        << "+-----------+-----------+--------------------+-----------------+-------------------+"
        << std::endl;
    std::cout
        << "| domain ID | scalar ID | total scalar value | domain integral | mean scalar value |"
        << std::endl;
    std::cout
        << "+-----------+-----------+--------------------+-----------------+-------------------+"
        << std::endl;
  }
}

/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
void ScaTra::OutputScalarsStrategyBase::finalize_screen_output()
{
  if (myrank_ == 0)
  {
    std::cout
        << "+-----------+-----------+--------------------+-----------------+-------------------+"
        << std::endl;
  }
}

/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
void ScaTra::OutputScalarsStrategyBase::output_total_and_mean_scalars(
    const ScaTraTimIntImpl* const scatratimint, const int num)
{
  // evaluate domain/condition integrals
  evaluate_integrals(scatratimint);

  // print evaluated data to screen
  print_header_to_screen(scatratimint->discret_->Name());
  print_to_screen();
  finalize_screen_output();

  // write evaluated data to file
  const auto output_data = prepare_csv_output();

  runtime_csvwriter_->WriteDataToFile(scatratimint->Time(), scatratimint->Step(), output_data);
}


/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void ScaTra::OutputScalarsStrategyBase::Init(const ScaTraTimIntImpl* const scatratimint)
{
  myrank_ = scatratimint->myrank_;

  output_mean_grad_ = Core::UTILS::IntegralValue<bool>(
      *scatratimint->ScatraParameterList(), "OUTPUTSCALARSMEANGRAD");

  output_micro_dis_ = scatratimint->MacroScale() and scatratimint->NdsMicro();

  std::string filename("");
  if (scatratimint->ScatraParameterList()->isParameter("output_file_name_discretization") and
      scatratimint->ScatraParameterList()->get<bool>("output_file_name_discretization"))
    filename = scatratimint->discretization()->Name() + "_scalarvalues";
  else
    filename = "scalarvalues";

  runtime_csvwriter_.emplace(myrank_, *scatratimint->DiscWriter()->output(), filename);
  init_strategy_specific(scatratimint);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void ScaTra::OutputScalarsStrategyDomain::init_strategy_specific(
    const ScaTraTimIntImpl* const scatratimint)
{
  if (not scatratimint->scalarhandler_->EqualNumDof())
  {
    FOUR_C_THROW(
        "Output of scalars for entire domain not valid for different numbers of DOFs per node in "
        "ScaTra discretization. \n"
        "Use option `by_condition' for the parameter 'OUTPUTSCALARS' in the SCALAR TRANSPORT "
        "DYNAMIC section instead!");
  }

  numdofpernode_ = scatratimint->scalarhandler_->NumDofPerNode();
  numscal_ = scatratimint->scalarhandler_->NumScal();

  // initialize result vectors associated with entire domain
  totalscalars_[dummy_domain_id_].resize(numdofpernode_, 0.0);
  meanscalars_[dummy_domain_id_].resize(numdofpernode_, 0.0);
  meangradients_[dummy_domain_id_].resize(numscal_, 0.0);
  domainintegral_[dummy_domain_id_] = 0.0;

  FOUR_C_ASSERT(runtime_csvwriter_.has_value(), "internal error: runtime csv writer not created.");

  // register output in csv writer
  runtime_csvwriter_->register_data_vector("Integral of entire domain", 1, 16);
  for (int k = 0; k < numdofpernode_; ++k)
  {
    runtime_csvwriter_->register_data_vector(
        "Total value of scalar " + std::to_string(k + 1) + " in entire domain", 1, 16);
    runtime_csvwriter_->register_data_vector(
        "Mean value of scalar " + std::to_string(k + 1) + " in entire domain", 1, 16);
  }
  if (output_mean_grad_)
  {
    for (int k = 0; k < numscal_; ++k)
    {
      runtime_csvwriter_->register_data_vector(
          "Mean value of gradient of scalar " + std::to_string(k + 1) + " in entire domain", 1, 16);
    }
  }
}

/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
void ScaTra::OutputScalarsStrategyDomain::print_to_screen()
{
  if (myrank_ == 0)
  {
    std::cout
        << "+-----------+-----------+--------------------+-----------------+-------------------+"
        << std::endl;

    for (int k = 0; k < numdofpernode_; ++k)
    {
      std::cout << "|  overall  |    " << std::setw(2) << k + 1 << "     |    " << std::scientific
                << std::setprecision(6) << totalscalars_[dummy_domain_id_][k] << "    |   "
                << domainintegral_[dummy_domain_id_] << "  |    "
                << meanscalars_[dummy_domain_id_][k] << "   |" << std::endl;
    }
    if (output_mean_grad_)
    {
      for (int k = 0; k < numscal_; ++k)
      {
        std::cout << "|  overall  |  grad. " << k + 1 << "  |        ----        |   "
                  << domainintegral_[dummy_domain_id_] << "  |    "
                  << meangradients_[dummy_domain_id_][k] << "   |" << std::endl;
      }
    }
  }
}

/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
std::map<std::string, std::vector<double>> ScaTra::OutputScalarsStrategyDomain::prepare_csv_output()
{
  FOUR_C_ASSERT(runtime_csvwriter_.has_value(), "internal error: runtime csv writer not created.");
  std::map<std::string, std::vector<double>> output_data;
  output_data["Integral of entire domain"] = {domainintegral_[dummy_domain_id_]};

  for (int k = 0; k < numdofpernode_; ++k)
  {
    output_data["Total value of scalar " + std::to_string(k + 1) + " in entire domain"] = {
        totalscalars_[dummy_domain_id_][k]};
    output_data["Mean value of scalar " + std::to_string(k + 1) + " in entire domain"] = {
        meanscalars_[dummy_domain_id_][k]};
  }
  if (output_mean_grad_)
  {
    for (int k = 0; k < numscal_; ++k)
    {
      output_data["Mean value of gradient of scalar " + std::to_string(k + 1) +
                  " in entire domain"] = {meangradients_[dummy_domain_id_][k]};
    }
  }
  return output_data;
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void ScaTra::OutputScalarsStrategyDomain::evaluate_integrals(
    const ScaTraTimIntImpl* const scatratimint)
{
  Teuchos::ParameterList eleparams;

  // fill parameter list and set states in discretization
  prepare_evaluate(scatratimint, eleparams);

  // initialize result vector
  // first components = scalar integrals, last component = domain integral
  const int num_active_scalars = output_mean_grad_ ? numscal_ : 0;
  auto scalars =
      Teuchos::rcp(new Core::LinAlg::SerialDenseVector(numdofpernode_ + 1 + num_active_scalars));

  // perform integration
  scatratimint->discret_->EvaluateScalars(eleparams, scalars);

  // extract domain integral
  domainintegral_[dummy_domain_id_] = (*scalars)[numdofpernode_];

  // compute results
  for (int k = 0; k < numdofpernode_; ++k)
  {
    totalscalars_[dummy_domain_id_][k] = (*scalars)[k];
    meanscalars_[dummy_domain_id_][k] = (*scalars)[k] / domainintegral_[dummy_domain_id_];
  }
  if (output_mean_grad_)
  {
    for (int k = 0; k < numscal_; ++k)
    {
      const double mean_gradient =
          (*scalars)[numdofpernode_ + 1 + k] / domainintegral_[dummy_domain_id_];

      meangradients_[dummy_domain_id_][k] = std::abs(mean_gradient) < 1.0e-10 ? 0.0 : mean_gradient;
    }
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void ScaTra::OutputScalarsStrategyCondition::init_strategy_specific(
    const ScaTraTimIntImpl* const scatratimint)
{
  // extract conditions for calculation of total and mean values of transported scalars
  scatratimint->discret_->GetCondition("TotalAndMeanScalar", conditions_);

  if (conditions_.empty())
  {
    FOUR_C_THROW(
        "No 'TotalAndMeanScalar' conditions defined for scalar output by condition!\n"
        "Change the parameter 'OUTPUTSCALARS' in the SCALAR TRANSPORT DYNAMIC section \n"
        "or include 'DESIGN TOTAL AND MEAN SCALAR' conditions!");
  }

  numdofpernodepercondition_.clear();
  numscalpercondition_.clear();

  // loop over all conditions
  for (auto* condition : conditions_)
  {
    // extract condition ID
    const int condid = condition->parameters().Get<int>("ConditionID");

    // determine the number of dofs on the current condition
    const int numdofpernode = scatratimint->scalarhandler_->num_dof_per_node_in_condition(
        *condition, scatratimint->discret_);
    const int numscalpernode = scatratimint->scalarhandler_->NumScal();

    // save the number of dofs on the current condition
    numdofpernodepercondition_.insert(std::make_pair(condid, numdofpernode));
    numscalpercondition_.insert(std::make_pair(condid, numscalpernode));

    // initialize result vectors associated with current condition
    totalscalars_[condid].resize(numdofpernodepercondition_[condid], 0.0);
    meanscalars_[condid].resize(numdofpernodepercondition_[condid], 0.0);
    meangradients_[condid].resize(numscalpercondition_[condid], 0.0);
    domainintegral_[condid] = 0.0;

    // micro dis has only one dof
    if (output_micro_dis_) micromeanscalars_[condid].resize(1, 0.0);

    FOUR_C_ASSERT(
        runtime_csvwriter_.has_value(), "internal error: runtime csv writer not created.");

    // register all data vectors
    runtime_csvwriter_->register_data_vector("Integral of domain " + std::to_string(condid), 1, 16);
    for (int k = 0; k < numdofpernodepercondition_[condid]; ++k)
    {
      runtime_csvwriter_->register_data_vector(
          "Total value of scalar " + std::to_string(k + 1) + " in domain " + std::to_string(condid),
          1, 16);
      runtime_csvwriter_->register_data_vector(
          "Mean value of scalar " + std::to_string(k + 1) + " in domain " + std::to_string(condid),
          1, 16);
    }
    if (output_mean_grad_)
    {
      for (int k = 0; k < numscalpercondition_[condid]; ++k)
      {
        runtime_csvwriter_->register_data_vector("Mean value of gradient of scalar " +
                                                     std::to_string(k + 1) + " in domain " +
                                                     std::to_string(condid),
            1, 16);
      }
    }
    if (output_micro_dis_)
    {
      runtime_csvwriter_->register_data_vector(
          "Mean value of micro scalar in domain " + std::to_string(condid), 1, 16);
    }
  }
}

/*--------------------------------------------------------------------------------------*
 |  Initialize output class                                            kremheller 11/19 |
 *--------------------------------------------------------------------------------------*/
void ScaTra::OutputDomainIntegralStrategy::Init(const ScaTraTimIntImpl* const scatratimint)
{
  // extract conditions for calculation of total and mean values of transported scalars
  scatratimint->discretization()->GetCondition("DomainIntegral", conditionsdomain_);
  scatratimint->discretization()->GetCondition("BoundaryIntegral", conditionsboundary_);

  if (conditionsdomain_.empty() && conditionsboundary_.empty())
  {
    FOUR_C_THROW(
        "No 'DomainIntegral' or 'BoundaryIntegral' conditions defined for scalar output by "
        "condition!\n"
        "Change the parameter 'COMPUTEINTEGRALS' in the SCALAR TRANSPORT DYNAMIC section \n"
        "or include 'DESIGN DOMAIN INTEGRAL * CONDITIONS' conditions!");
  }

  domainintegralvalues_.resize(conditionsdomain_.size());
  boundaryintegralvalues_.resize(conditionsboundary_.size());
}

/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
void ScaTra::OutputScalarsStrategyCondition::print_to_screen()
{
  if (myrank_ == 0)
  {
    for (auto* condition : conditions_)
    {
      // extract condition ID
      const int condid = condition->parameters().Get<int>("ConditionID");

      // determine the number of dofs on the current condition
      const int numdofpernode = numdofpernodepercondition_[condid];

      std::ostringstream condid_string;
      condid_string << "|    " << std::setw(2) << condid << "     |";
      for (int k = 0; k < numdofpernode; ++k)
      {
        std::cout << condid_string.str() << "    " << std::setw(2) << k + 1 << "     |    "
                  << std::scientific << std::setprecision(6) << totalscalars_[condid][k]
                  << "    |   " << domainintegral_[condid] << "  |    " << meanscalars_[condid][k]
                  << "   |" << std::endl;
      }
      if (output_mean_grad_)
      {
        for (int k = 0; k < numscalpercondition_[condid]; ++k)
        {
          std::cout << condid_string.str() << "  grad. " << k + 1 << "  |        ----        |   "
                    << domainintegral_[condid] << "  |    " << meangradients_[condid][k] << "   |"
                    << std::endl;
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
std::map<std::string, std::vector<double>>
ScaTra::OutputScalarsStrategyCondition::prepare_csv_output()
{
  std::map<std::string, std::vector<double>> output_data;

  for (auto* condition : conditions_)
  {
    // extract condition ID
    const int condid = condition->parameters().Get<int>("ConditionID");

    FOUR_C_ASSERT(
        runtime_csvwriter_.has_value(), "internal error: runtime csv writer not created.");

    output_data["Integral of domain " + std::to_string(condid)] = {domainintegral_[condid]};

    for (int k = 0; k < numdofpernodepercondition_[condid]; ++k)
    {
      output_data["Total value of scalar " + std::to_string(k + 1) + " in domain " +
                  std::to_string(condid)] = {totalscalars_[condid][k]};
      output_data["Mean value of scalar " + std::to_string(k + 1) + " in domain " +
                  std::to_string(condid)] = {meanscalars_[condid][k]};
    }
    if (output_mean_grad_)
    {
      for (int k = 0; k < numscalpercondition_[condid]; ++k)
      {
        output_data["Mean value of gradient of scalar " + std::to_string(k + 1) + " in domain " +
                    std::to_string(condid)] = {meangradients_[condid][k]};
      }
    }
    if (output_micro_dis_)
    {
      output_data["Mean value of micro scalar in domain " + std::to_string(condid)] = {
          micromeanscalars_[condid][0]};
    }
  }
  return output_data;
}

/*--------------------------------------------------------------------------------------*
 |  evaluate domain integrals and print them to file and screen        kremheller 11/19 |
 *--------------------------------------------------------------------------------------*/
void ScaTra::OutputDomainIntegralStrategy::evaluate_integrals_and_print_results(
    const ScaTraTimIntImpl* const scatratimint, const std::string& condstring)
{
  // initialize label
  std::string label;

  // create parameter list
  Teuchos::ParameterList condparams;

  // determine label and element action depending on whether domain or boundary integrals are
  // relevant
  if (condstring == "BoundaryIntegral")
  {
    label = "Boundary";
    Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
        "action", ScaTra::BoundaryAction::calc_boundary_integral, condparams);
  }
  else if (condstring == "DomainIntegral")
  {
    label = "Domain";
    Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
        "action", ScaTra::Action::calc_domain_integral, condparams);
  }
  else
    FOUR_C_THROW("Invalid condition name!");

  // extract conditions for computation of domain or boundary integrals
  std::vector<Core::Conditions::Condition*> conditions;
  Teuchos::RCP<Core::FE::Discretization> discret = scatratimint->discretization();
  const int myrank = discret->Comm().MyPID();
  discret->GetCondition(condstring, conditions);

  // print header
  if (conditions.size() > 0 and myrank == 0)
  {
    std::cout << std::endl << label + " integrals:" << std::endl;
    std::cout << "+----+-------------------------+" << std::endl;
    std::cout << "| ID | value of integral       |" << std::endl;
  }

  // loop over all conditions
  for (int condid = 0; condid < static_cast<int>(conditions.size()); ++condid)
  {
    // initialize one-component result vector for value of current domain or boundary integral
    Teuchos::RCP<Core::LinAlg::SerialDenseVector> integralvalue =
        Teuchos::rcp(new Core::LinAlg::SerialDenseVector(1));

    // compute value of current domain or boundary integral
    discret->EvaluateScalars(condparams, integralvalue, condstring, condid);

    // output results to screen and file
    if (myrank == 0)
    {
      // print results to screen
      std::cout << "| " << std::setw(2) << condid << " |        " << std::setw(6) << std::scientific
                << std::setprecision(3) << (*integralvalue)(0) << "        |" << std::endl;

      // set file name
      const std::string filename(Global::Problem::Instance()->OutputControlFile()->file_name() +
                                 "." + label + "_integrals.csv");

      // open file in appropriate mode and write header at beginning
      std::ofstream file;
      if (scatratimint->Step() == 0 and condid == 0)
      {
        file.open(filename.c_str(), std::fstream::trunc);
        file << "Step,Time," << label << "ID," << label << "Integral" << std::endl;
      }
      else
        file.open(filename.c_str(), std::fstream::app);

      // write value of current domain or boundary integral to file
      file << scatratimint->Step() << "," << scatratimint->Time() << "," << condid << ','
           << std::scientific << std::setprecision(6) << (*integralvalue)(0) << std::endl;

      // close file
      file.close();
    }  // if(myrank_ == 0)

    if (condstring == "BoundaryIntegral")
      boundaryintegralvalues_[condid] = (*integralvalue)(0);
    else if (condstring == "DomainIntegral")
      domainintegralvalues_[condid] = (*integralvalue)(0);
    else
      FOUR_C_THROW("Invalid condition name!");
  }  // loop over all conditions

  // print finish line to screen
  if (conditions.size() and myrank == 0)
    std::cout << "+----+-------------------------+" << std::endl;
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void ScaTra::OutputScalarsStrategyCondition::evaluate_integrals(
    const ScaTraTimIntImpl* const scatratimint)
{
  Teuchos::ParameterList eleparams;

  // fill parameter list and set states in discretization
  prepare_evaluate(scatratimint, eleparams);

  // evaluate scalar integrals and domain integral for each subdomain
  for (auto* condition : conditions_)
  {
    // extract condition ID
    const int condid = condition->parameters().Get<int>("ConditionID");

    // determine the number of dofs on the current condition
    const int numdofpernode = numdofpernodepercondition_[condid];
    const int numscalpernode = numscalpercondition_[condid];

    // initialize result vector
    // first components = scalar integrals, last component = domain integral
    const int num_active_scalars = output_mean_grad_ ? numscalpernode : 0;
    const int num_micro_dis = output_micro_dis_ ? 1 : 0;
    auto scalars = Teuchos::rcp(new Core::LinAlg::SerialDenseVector(
        numdofpernode + 1 + num_active_scalars + num_micro_dis));

    // perform integration
    scatratimint->discret_->EvaluateScalars(eleparams, scalars, "TotalAndMeanScalar", condid);

    // extract domain integral
    domainintegral_[condid] = (*scalars)[numdofpernode];

    // compute results
    for (int k = 0; k < numdofpernode; ++k)
    {
      totalscalars_[condid][k] = (*scalars)[k];
      meanscalars_[condid][k] = (*scalars)[k] / domainintegral_[condid];
    }
    if (output_mean_grad_)
    {
      for (int k = 0; k < numscalpercondition_[condid]; ++k)
      {
        const double mean_gradient = (*scalars)[numdofpernode + 1 + k] / domainintegral_[condid];
        meangradients_[condid][k] = std::abs(mean_gradient) < 1.0e-10 ? 0.0 : mean_gradient;
      }
    }
    if (output_micro_dis_)
      micromeanscalars_[condid][0] = (*scalars)[scalars->length() - 1] / domainintegral_[condid];
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void ScaTra::OutputScalarsStrategyDomainAndCondition::init_strategy_specific(
    const ScaTraTimIntImpl* const scatratimint)
{
  // initialize base classes
  OutputScalarsStrategyCondition::init_strategy_specific(scatratimint);
  OutputScalarsStrategyDomain::init_strategy_specific(scatratimint);
}

/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
void ScaTra::OutputScalarsStrategyDomainAndCondition::print_to_screen()
{
  OutputScalarsStrategyCondition::print_to_screen();
  OutputScalarsStrategyDomain::print_to_screen();
}

/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
std::map<std::string, std::vector<double>>
ScaTra::OutputScalarsStrategyDomainAndCondition::prepare_csv_output()
{
  std::map<std::string, std::vector<double>> output_data;
  auto output_data_condition = OutputScalarsStrategyCondition::prepare_csv_output();
  auto output_data_domain = OutputScalarsStrategyDomain::prepare_csv_output();

  // merge both maps
  output_data.merge(output_data_condition);
  output_data.merge(output_data_domain);

  return output_data;
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void ScaTra::OutputScalarsStrategyDomainAndCondition::evaluate_integrals(
    const ScaTraTimIntImpl* const scatratimint)
{
  // call base classes
  OutputScalarsStrategyCondition::evaluate_integrals(scatratimint);
  OutputScalarsStrategyDomain::evaluate_integrals(scatratimint);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  set up handler class                                   vuong   04/16|
 *----------------------------------------------------------------------*/
void ScaTra::ScalarHandler::Setup(const ScaTraTimIntImpl* const scatratimint)
{
  // save reference to discretization for convenience
  const Teuchos::RCP<Core::FE::Discretization>& discret = scatratimint->discretization();

  // initialize set of all number of dofs per node on this proc
  std::set<int> mynumdofpernode;

  // collect number of DoFs for each node on this proc
  for (int i = 0; i < discret->NumMyRowNodes(); i++)
    mynumdofpernode.insert(discret->NumDof(0, discret->lRowNode(i)));

  // number of different numbers of dofs on all procs
  int maxsize = 0;
  // number of different numbers of dofs on this procs
  int mysize = static_cast<int>(mynumdofpernode.size());
  // communicate
  discret->Comm().MaxAll(&mysize, &maxsize, 1);

  // copy mynumdofpernode into std::vector for communication
  std::vector<int> vecmynumdofpernode(mynumdofpernode.begin(), mynumdofpernode.end());
  vecmynumdofpernode.resize(maxsize, 0);

  // initialize std::vector for communication
  std::vector<int> vecnumdofpernode(maxsize * discret->Comm().NumProc(), 0);

  // communicate
  discret->Comm().GatherAll(vecmynumdofpernode.data(), vecnumdofpernode.data(),
      static_cast<int>(vecmynumdofpernode.size()));

  // copy back into set
  for (int& ndofpernode : vecnumdofpernode)
    if (ndofpernode != 0) numdofpernode_.insert(ndofpernode);

  // check for equal number of Dofs on all nodes in whole discretization
  equalnumdof_ = (numdofpernode_.size() == 1);

  // done
  issetup_ = true;
}

/*-------------------------------------------------------------------------*
|  Determine number of DoFs per node in given condition       vuong   04/16|
 *-------------------------------------------------------------------------*/
int ScaTra::ScalarHandler::num_dof_per_node_in_condition(
    const Core::Conditions::Condition& condition,
    const Teuchos::RCP<const Core::FE::Discretization>& discret) const
{
  check_is_setup();

  // get all nodes in condition
  const std::vector<int>* nodegids = condition.GetNodes();

  // initialize set of all number of dofs per node on all procs
  std::set<int> numdofpernode;
  // initialize set of all number of dofs per node on this proc
  std::set<int> mynumdofpernode;

  // determine all number of dofs per node on this proc
  for (int nodegid : *nodegids)
  {
    // if on this proc
    if (discret->NodeRowMap()->MyGID(nodegid))
    {
      // get node
      Core::Nodes::Node* curnode = discret->gNode(nodegid);
      // save number of dofs in set
      mynumdofpernode.insert(discret->NumDof(0, curnode));
    }
  }

  // number of different numbers of dofs on all procs
  int maxsize = 0;
  // number of different numbers of dofs on this procs
  int mysize = static_cast<int>(mynumdofpernode.size());
  // communicate
  discret->Comm().MaxAll(&mysize, &maxsize, 1);

  // copy mynumdofpernode into std::vector for communication
  std::vector<int> vecmynumdofpernode(mynumdofpernode.begin(), mynumdofpernode.end());
  vecmynumdofpernode.resize(maxsize, 0);

  // initialize std::vector for communication
  std::vector<int> vecnumdofpernode(maxsize * discret->Comm().NumProc(), 0);

  // communicate
  discret->Comm().GatherAll(vecmynumdofpernode.data(), vecnumdofpernode.data(),
      static_cast<int>(vecmynumdofpernode.size()));

  // copy back into set
  for (int& ndofpernode : vecnumdofpernode)
    if (ndofpernode != 0) numdofpernode.insert(ndofpernode);

  // this is the maximum number of dofs per node within the discretization on all procs
  int maxnumdofpernode = *numdofpernode.rbegin();

  if (numdofpernode.size() != 1)
  {
    FOUR_C_THROW(
        "Different number of DOFs within condition. This is not supported. Split the condition "
        "in your input file!");
  }

  return maxnumdofpernode;
}

/*-------------------------------------------------------------------------*
|  Get number of DoFs per node within discretization            vuong   04/16|
 *-------------------------------------------------------------------------*/
int ScaTra::ScalarHandler::NumDofPerNode() const
{
  check_is_setup();

  if (not equalnumdof_)
  {
    FOUR_C_THROW(
        "Number of DOFs per node is not equal for all nodes within the ScaTra discretization!\n"
        "Calling this method is not valid in this case!");
  }
  return *(numdofpernode_.rbegin());
}

/*-----------------------------------------------------------------------------*
 |  check if class is set up                                       rauch 09/16 |
 *-----------------------------------------------------------------------------*/
void ScaTra::ScalarHandler::check_is_setup() const
{
  if (not issetup_) FOUR_C_THROW("ScalarHanlder is not set up. Call Setup() first.");
}

FOUR_C_NAMESPACE_CLOSE
