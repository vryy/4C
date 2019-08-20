/*----------------------------------------------------------------------*/
/*! \file

\brief base level-set algorithm: collection of functions related to reinitialization

    detailed description in header file levelset_algorithm.H

\level 2

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
 *------------------------------------------------------------------------------------------------*/


#include "levelset_algorithm.H"
#include "levelset_intersection_utils.H"
#include "../drt_lib/drt_periodicbc.H"
#include "../drt_geometry/integrationcell.H"
#include "../drt_geometry/geo_utils.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/standardtypes_cpp.H"  // for EPS13 and so on
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"

#include "../drt_io/io_gmsh.H"


/*----------------------------------------------------------------------*
 | algebraic reinitialization via solution of equation  rasthofer 09/13 |
 | pde-based reinitialization according to Sussman 1994                 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::ReinitEq()
{
  if (myrank_ == 0)
    std::cout << "\n---------------------------------------  REINITIALIZATION SOLVER  "
                 "----------------------------\n";

  // -----------------------------------------------------------------
  //            prepare time loop for reinitialization
  // -----------------------------------------------------------------
  // set vectors and all further quantities such that the reinitialization
  // can be solved within the existing framework
  PrepareTimeLoopReinit();

  // -----------------------------------------------------------------
  //            time loop for reinitialization equation
  // -----------------------------------------------------------------
  TimeLoopReinit();

  // -----------------------------------------------------------------
  //            complete reinitialization via equation
  // -----------------------------------------------------------------

  // cleaning of all necessary modifications
  FinishTimeLoopReinit();

  return;
}


/*-------------------------------------------------------------------*
 | set element parameters for reinitialization equation   fang 08/15 |
 *-------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::SetReinitializationElementParameters(
    bool calcinitialtimederivative) const
{
  // create element parameter list
  Teuchos::ParameterList eleparams;

  // set action for elements
  eleparams.set<int>("action", SCATRA::set_lsreinit_scatra_parameter);

  // reinitialization equation is given in convective form
  eleparams.set<int>("convform", INPAR::SCATRA::convform_convective);

  // no ALE intended
  eleparams.set("isale", false);

  // parameters for stabilization, which are the same as for the level-set equation (if turned on)
  eleparams.sublist("stabilization") = params_->sublist("STABILIZATION");

  // set flag for writing the flux vector fields
  eleparams.set<int>("calcflux_domain", calcflux_domain_);

  // set vector containing IDs of scalars for which flux vectors are calculated
  eleparams.set<Teuchos::RCP<std::vector<int>>>("writefluxids", writefluxids_);

  // set level-set reinitialization specific parameters
  eleparams.sublist("REINITIALIZATION") = levelsetparams_->sublist("REINITIALIZATION");

  // turn off stabilization and artificial diffusivity when calculating initial time derivative
  if (calcinitialtimederivative)
  {
    eleparams.sublist("REINITIALIZATION").set<std::string>("STABTYPEREINIT", "no_stabilization");
    eleparams.sublist("REINITIALIZATION").set<std::string>("ARTDIFFREINIT", "no");
  }

  // parameters for finite difference check
  eleparams.set<int>("fdcheck", fdcheck_);
  eleparams.set<double>("fdcheckeps", fdcheckeps_);
  eleparams.set<double>("fdchecktol", fdchecktol_);

  // overwrite some values in general stabilization parameter list by modified values in levelset
  // reinitialization parameter list
  eleparams.sublist("stabilization")
      .set<std::string>("DEFINITION_TAU",
          eleparams.sublist("REINITIALIZATION").get<std::string>("DEFINITION_TAU_REINIT"));
  eleparams.sublist("stabilization")
      .set<std::string>(
          "STABTYPE", eleparams.sublist("REINITIALIZATION").get<std::string>("STABTYPEREINIT"));
  eleparams.sublist("stabilization").set<std::string>("SUGRVEL", "no");
  eleparams.sublist("stabilization")
      .set<std::string>("DEFINITION_ASSGD",
          eleparams.sublist("REINITIALIZATION").get<std::string>("DEFINITION_ARTDIFFREINIT"));

  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  return;
}


/*----------------------------------------------------------------------*
 | set time parameters for reinitialization equation    rasthofer 12/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::SetReinitializationElementTimeParameters()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action", SCATRA::set_time_parameter);

  eleparams.set<bool>("using generalized-alpha time integration", false);
  eleparams.set<bool>("using stationary formulation", false);
  // reinitialization equation only implemented incrementally, since it is nonlinear
  eleparams.set<bool>("incremental solver", true);

  eleparams.set<double>("time-step length", dtau_);
  eleparams.set<double>("total time", dtau_ * pseudostep_);
  eleparams.set<double>("time factor", thetareinit_ * dtau_);
  eleparams.set<double>("alpha_F", 1.0);

  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  return;
}


/*----------------------------------------------------------------------*
 | prepare internal time loop for reinitialization equation             |
 |                                                      rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::PrepareTimeLoopReinit()
{
  // set switch flag to true to active reinitialization specific parts
  switchreinit_ = true;

  // initial or start phi of reinitialization process
  initialphireinit_->Update(1.0, *phinp_, 0.0);
  phin_->Update(1.0, *phinp_, 0.0);

  // set internal step counter to zero
  pseudostep_ = 0;

  // set time-integration parameters for reinitialization equation
  SetReinitializationElementTimeParameters();
  // set element parameters for reinitialization equation
  SetReinitializationElementParameters();

  return;
}


/*----------------------------------------------------------------------*
 | internal time loop for reinitialization equation     rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::TimeLoopReinit()
{
  // time measurement: time loop
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:  + reinitialization time loop");

  //       e.g., steady state of interface nodal values
  //             integrated gradient norm
  bool converged = false;
  while (pseudostep_ < pseudostepmax_ and not converged)
  {
    // -------------------------------------------------------------------
    //                  prepare time step
    // -------------------------------------------------------------------
    PrepareTimeStepReinit();

    // -------------------------------------------------------------------
    //                  solve nonlinear equation
    // -------------------------------------------------------------------
    SolveReinit();

    // -------------------------------------------------------------------
    //                  interface correction
    // -------------------------------------------------------------------
    if (reinitcorrector_) CorrectionReinit();

    // -------------------------------------------------------------------
    //                        check for convergence
    // -------------------------------------------------------------------
    converged = ConvergenceCheckReinit();

    // -------------------------------------------------------------------
    //                        update solution
    // -------------------------------------------------------------------
    UpdateReinit();
  }

  return;
}


/*----------------------------------------------------------------------*
 | clean internal time loop for reinitialization equation               |
 |                                                      rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::FinishTimeLoopReinit()
{
  // reset quantites that may have been overwritten
  // reset internal step counter
  pseudostep_ = 0;

  // reset time-integration parameters for element evaluation
  SetElementTimeParameter();
  // reset general parameters for element evaluation
  SetElementGeneralParameters();
  SetElementTurbulenceParameters();

  //  {
  //    // turn on/off screen output for writing process of Gmsh postprocessing file
  //  const bool screen_out = true;
  //
  //  // create Gmsh postprocessing file
  //  const std::string filename =
  //  IO::GMSH::GetNewFileNameAndDeleteOldFiles("reinitialization_scalar", step_, 500, screen_out,
  //  discret_->Comm().MyPID()); std::ofstream gmshfilecontent(filename.c_str());
  //
  //  {
  //    // add 'View' to Gmsh postprocessing file
  //    gmshfilecontent << "View \" " << "Inital Phi \" {" << std::endl;
  //    // draw scalar field 'Phinp' for every element
  //    IO::GMSH::ScalarFieldToGmsh(discret_,initialphireinit_,gmshfilecontent);
  //    gmshfilecontent << "};" << std::endl;
  //  }
  //
  //  {
  //    // add 'View' to Gmsh postprocessing file
  //    gmshfilecontent << "View \" " << "Final Phi \" {" << std::endl;
  //    // draw scalar field 'Phinp' for every element
  //    IO::GMSH::ScalarFieldToGmsh(discret_,phinp_,gmshfilecontent);
  //    gmshfilecontent << "};" << std::endl;
  //  }
  //
  //  {
  //    // add 'View' to Gmsh postprocessing file
  //    gmshfilecontent << "View \" " << "Reinit Velocity \" {" << std::endl;
  //    // draw vector field 'Convective Velocity' for every element
  //    IO::GMSH::VectorFieldNodeBasedToGmsh(discret_,reinitvel_,gmshfilecontent);
  //    gmshfilecontent << "};" << std::endl;
  //  }
  //  gmshfilecontent.close();
  //  if (screen_out) std::cout << " done" << std::endl;
  //  }
  //  dserror("ENDE");

  return;
}


/*----------------------------------------------------------------------*
 | convergence check for reinitialization equation      rasthofer 03/14 |
 *----------------------------------------------------------------------*/
bool SCATRA::LevelSetAlgorithm::ConvergenceCheckReinit()
{
  bool abortreinitloop = false;

  if (reinit_tol_ > 0.0)
  {
    if (myrank_ == 0)
      std::cout << "## WARNING: convergence criterion for reinitialization equation not yet "
                   "carefully checked"
                << std::endl;

    // stop criterion according to Sussman et al 1994
    //         sum_(nodes A with abs(phi_n) < alpha) abs(phi_A_n+1 -phi_A_n)
    //  err = --------------------------------------------------------------- < dtau*h^2
    //                      sum_(nodes A with abs(phi_n) < alpha) 1

    double local_sum = 0.0;
    int local_num_nodes = 0;

    for (int inode = 0; inode < discret_->NumMyRowNodes(); inode++)
    {
      if (std::abs((*phin_)[inode]) < reinitbandwidth_)
      {
        local_sum += std::abs((*phinp_)[inode] - (*phin_)[inode]);
        local_num_nodes += 1;
      }
    }

    // communicate sums
    double global_sum = 0.0;
    int global_num_nodes = 0;
    discret_->Comm().SumAll(&local_sum, &global_sum, 1);
    discret_->Comm().SumAll(&local_num_nodes, &global_num_nodes, 1);

    // compute current error in band
    const double err = global_sum / ((double)global_num_nodes);

    if (myrank_ == 0)
      std::cout << "Convergence Check Reinitialization: Err  " << err << "  Tol  " << reinit_tol_
                << " Number of nodes in band  " << global_num_nodes << std::endl;

    if (err <
        reinit_tol_)  //(dtau_*char_ele_length*char_ele_length) suggested by Sussman et al 1994
      abortreinitloop = true;
  }

  return abortreinitloop;
}


/*----------------------------------------------------------------------*
 | setup the variables to do a new reinitialization time step           |
 |                                                      rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::PrepareTimeStepReinit()
{
  // prepare first time step
  if (pseudostep_ == 0)
  {
    // TODO: Restructure enforcement of Dirichlet boundary conditions on phin_
    ApplyDirichletBC(time_, phin_, Teuchos::null);
    CalcInitialTimeDerivative();
  }

  // increment time and step
  pseudostep_ += 1;

  // -------------------------------------------------------------------
  // set part of the rhs vector belonging to the old timestep
  // -------------------------------------------------------------------
  SetOldPartOfRighthandside();
  SetReinitializationElementTimeParameters();

  // -------------------------------------------------------------------
  // compute node-based velocity field
  // -------------------------------------------------------------------
  // TODO phin oder phinp: fuer phinp in AddProblemSpecificParametersAndVectors
  //     muss vor oder gleich zu Beginn von AssembleMatAndRHS gerufen werden
#ifdef USE_PHIN_FOR_VEL
  if (useprojectedreinitvel_ == INPAR::SCATRA::vel_reinit_node_based) CalcNodeBasedReinitVel();
#endif

  return;
}


/*----------------------------------------------------------------------*
 | calculate node-based velocity field via L2-projection                |
 | for reinitialization                                 rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::CalcNodeBasedReinitVel()
{
  // loop all space dimensions,
  // since assembler can only deal with one dof per node here
  for (int idim = 0; idim < 3; idim++)
  {
    // define vector for velocity component
    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    Teuchos::RCP<Epetra_Vector> velcomp = LINALG::CreateVector(*dofrowmap, true);
    velcomp->PutScalar(0.0);

    if (lsdim_ == INPAR::SCATRA::ls_3D or (lsdim_ == INPAR::SCATRA::ls_2Dx and idim != 0) or
        (lsdim_ == INPAR::SCATRA::ls_2Dy and idim != 1) or
        (lsdim_ == INPAR::SCATRA::ls_2Dz and idim != 2))
    {
      // zero out matrix and rhs entries
      sysmat_->Zero();
      residual_->PutScalar(0.0);

      // create the parameters for the discretization
      Teuchos::ParameterList eleparams;

      // parameters for the elements
      // action
      eleparams.set<int>("action", SCATRA::calc_node_based_reinit_velocity);
      // set current spatial direction
      // we have to loop the dimensions, since we merely have one dof per node here
      eleparams.set<int>("direction", idim);
      // activate reinitialization calculation routines
      eleparams.set<bool>("solve reinit eq", true);

      discret_->ClearState();  // TODO Caution if called from NonlinearSolve
      // set initial phi, i.e., solution of level-set equation
      discret_->SetState("phizero", initialphireinit_);

      switch (reinitaction_)
      {
        case INPAR::SCATRA::reinitaction_sussman:
        {
          // set phin as phi used for velocity
          // note:read as phinp in SysmatNodalVel()
#ifdef USE_PHIN_FOR_VEL
          discret_->SetState("phinp", phin_);
#else
          discret_->SetState("phinp", phinp_);
#endif
          break;
        }
        case INPAR::SCATRA::reinitaction_ellipticeq:
        {
          discret_->SetState("phinp", phinp_);
          break;
        }
        default:
        {
          dserror("Unknown reinitialization method for projection!");
          exit(EXIT_FAILURE);
        }
      }
      // call loop over elements
      discret_->Evaluate(eleparams, sysmat_, residual_);
      discret_->ClearState();

      // finalize the complete matrix
      sysmat_->Complete();

      // solve for velocity component
      solver_->Solve(sysmat_->EpetraOperator(), velcomp, residual_, true, true);

      SystemMatrix()->Reset();
      // reset the solver as well
      solver_->Reset();

      // TODO: add simple ILU/SGS/... for projection case, such that the standard system can be
      // solved by efficient AMG methods
    }

    // loop over all local nodes of scatra discretization
    for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
    {
      // store velocity in reinitialization velocity
      const double val = (*velcomp)[lnodeid];
      ((*nb_grad_val_)(idim))->ReplaceMyValues(1, &val, &lnodeid);
    }
  }

#if 0
  {
    // turn on/off screen output for writing process of Gmsh postprocessing file
  const bool screen_out = true;

  // create Gmsh postprocessing file
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("reinitialization_scalar", step_, 500, screen_out, discret_->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Inital Phi \" {" << std::endl;
    // draw scalar field 'Phinp' for every element
    IO::GMSH::ScalarFieldToGmsh(discret_,initialphireinit_,gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Projection \" {" << std::endl;
    // draw vector field 'Convective Velocity' for every element
    IO::GMSH::VectorFieldNodeBasedToGmsh(discret_,nb_grad_val_,gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }
  gmshfilecontent.close();
  if (screen_out) std::cout << " done" << std::endl;
  }
  dserror("ENDE");
#endif
  return;
}


/*----------------------------------------------------------------------*
 | contains call of nonlinear solver for reinitialization equation      |
 |                                                      rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::SolveReinit()
{
  // we simply call the NonlinearSolve (since the reinitialization equation is
  // indeed nonlinear), and all the rest concerning the correct action type and
  // parameters is handled via the switchreinit_-flag in the concrete time-integration
  // schemes for level-set problems
  NonlinearSolve();

  return;
}


/*----------------------------------------------------------------------*
 | correction step according to Sussman & Fatemi 1999   rasthofer 12/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::CorrectionReinit()
{
  if (myrank_ == 0) std::cout << "\n---------------  Correction projection\n";

  // this correction step should force the interface to stay fixed during reinitialization
  // according to Sussman & Fatemi 1999, it is given as
  //
  //    phinp_final = phinp_reinit + dtau * penaltyparameter * deriv H(phizero) ||nabla phizero||,
  //
  // which is used here in form of a projection

  // zero out matrix and rhs entries !
  sysmat_->Zero();
  residual_->PutScalar(0.0);

  // generate a parameterlist for communication and control
  Teuchos::ParameterList eleparams;
  // action for elements
  eleparams.set<int>("action", SCATRA::calc_mat_and_rhs_lsreinit_correction_step);
  eleparams.set<bool>("solve reinit eq", true);

  // set state vectors
  discret_->ClearState();
  discret_->SetState("phizero", initialphireinit_);
  discret_->SetState("phinp", phinp_);


  // call loop over elements
  discret_->Evaluate(eleparams, sysmat_, residual_);
  discret_->ClearState();

  // residual_->Print(std::cout);

  // finalize the complete matrix
  sysmat_->Complete();

  // solve for corrected phinp
  solver_->Solve(sysmat_->EpetraOperator(), phinp_, residual_, true, true);

  // phinp_->Print(std::cout);

  SystemMatrix()->Reset();
  // reset the solver as well
  solver_->Reset();

  return;
}


///*----------------------------------------------------------------------*
// | convergence phi                                      rasthofer 12/13 |
// *----------------------------------------------------------------------*/
// bool SCATRA::LevelSetAlgorithm::ReinitSteadyState()
//{
//  bool steady_state = false;
//
//  // abort reinitialization time loop if gradient of phi is smaller
//  // than prescribed tolerance
//  if (actgraderr <= gradtol)
//    abortreinitloop = true;
//
//  return steady_state;
//}
//#endif


#if 0
/*----------------------------------------------------------------------*
 | update solution after reinitialization step          rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::UpdateReinit(
  double& oldgraderr,
  bool&   stoploop)
{
  // allows for aborting reinitialization loop if reinitialization does not further improve the gradient
  // TODO name nicht gut
  INPAR::SCATRA::ReInitialStationaryCheck reinit_stationary_check
      = INPAR::SCATRA::reinit_stationarycheck_L1normintegrated;
  //TODO: get from parameters
  //  = DRT::INPUT::IntegralValue<INPAR::SCATRA::ReInitialStationaryCheck>(combustdynreinit_->sublist("COMBUSTION PDE REINITIALIZATION"),"STATIONARY_CHECK");

  double actgraderr = EvaluateGradientNormError();

  if((actgraderr >= oldgraderr) and
      reinit_stationary_check == INPAR::SCATRA::reinit_stationarycheck_L1normintegrated)
  {
     if (myrank_ == 0)
        std::cout << "---------------- no further gradient improvements -------------" << std::endl;
     // the reinitialization step did not further improve the gradient
     stoploop = true;
     // and keep the phinp from the last reinitialization step
     phinp_->Update(1.0,*phin_,0.0);

     return;
  }
  else
  {
     // update gradient for convergence check
     oldgraderr = actgraderr;
     // update solution
     UpdatePhiReinit();
     stoploop = false;
  }

  return;
}
#endif


#if 0
/*----------------------------------------------------------------------*
 |  calculate error in relative gradient norm           rasthofer 09/13 |
 *----------------------------------------------------------------------*/
double SCATRA::LevelSetAlgorithm::EvaluateGradientNormError()
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;

  // parameters for the elements
  p.set<int>("action",SCATRA::calc_error_reinit);
  p.set<double>("L1 integrated gradient error", 0.0);
  // volume of element for which gradient error is evaluated
  // usually band around interface
  // TODO: include band width: check current version
  p.set<double>("volume_error_band", 0.0);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);

  // get (squared) error values
  discret_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  discret_->ClearState();

  // get procs local errors
  double locL1gradienterr   = p.get<double>("L1 integrated gradient error");
  double locvolume          = p.get<double>("volume_error_band");

  // initialize global errors
  double L1gradient_err     = 0.0;
  double volume             = 0.0;
  double rel_gradient_err   = 0.0;

  // sum over all processors
  discret_->Comm().SumAll(&locL1gradienterr,&L1gradient_err,1);
  discret_->Comm().SumAll(&locvolume,&volume,1);

  // calculate relative gradient error
  if(fabs(volume) > 1e-014) rel_gradient_err = L1gradient_err / volume;
  else dserror("volume is smaller than 1e-14: check your band with");

  // print norms to screen
  if (myrank_ == 0)
  {
    printf("\nConvergence check for reinitialization:\n");
    printf("absolute gradient error || (||grad(phi)||-1.0) ||_L2(Omega) %15.8e\n"
           "relative gradient error                                     %15.8e unsing band of volume %15.8e\n\n",
           L1gradient_err,rel_gradient_err,volume);
  }

  return rel_gradient_err;
}
#endif


/*----------------------------------------------------------------------*
 | geometric reinitialization via distance to interface rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::ReinitGeo(const std::map<int, GEO::BoundaryIntCells>& interface)
{
  if (myrank_ == 0)
    std::cout << "---  reinitializing level-set field by computing distance to interface ..."
              << std::flush;

  // set switch flag to true to active reinitialization specific parts
  switchreinit_ = true;

  // map holding pbc nodes (masters or slaves) <pbc node id, distance to flame front>
  std::map<int, double> pbcnodes;

  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // determine the number of nodes per element
  int numnodesperele = 0;
  if (discret_->NumMyRowElements() <= 0)
    dserror("This discretization does not have any row elements.");
  switch (discret_->lRowElement(0)->Shape())
  {
    case DRT::Element::hex8:
      numnodesperele = 8;
      break;
    case DRT::Element::hex20:
      numnodesperele = 20;
      std::cout << "Warning, the fast signed distance reinitialization has not been tested with "
                   "hex20 elements!"
                << std::endl;
      break;
    case DRT::Element::hex27:
      numnodesperele = 27;
      std::cout << "Warning, the fast signed distance reinitialization has not been tested with "
                   "hex27 elements!"
                << std::endl;
      break;
    default:
    {
      dserror(
          "The fast signed distance reinitialization only supports hex8, hex20 and hex27 "
          "elements.");
      break;
    }
  }

  //========================================================================
  // get the following information about the pbc
  // - planenormaldirection e.g. (1,0,0)
  // - minimum in planenormaldirection
  // - maximum in planenormaldirection
  //========================================================================
  // std::vector<DRT::Condition*>* surfacepbcs = pbc_->ReturnSurfacePBCs();
  // get periodic surface boundary conditions
  std::vector<DRT::Condition*> surfacepbcs;
  discret_->GetCondition("SurfacePeriodic", surfacepbcs);
  if (surfacepbcs.empty()) discret_->GetCondition("LinePeriodic", surfacepbcs);

  std::vector<int> planenormal(0);
  std::vector<double> globalmins(0);
  std::vector<double> globalmaxs(0);

  for (size_t i = 0; i < surfacepbcs.size(); ++i)
  {
    const std::string* ismaster =
        surfacepbcs[i]->Get<std::string>("Is slave periodic boundary condition");
    if (*ismaster == "Master")
    {
      const int masterid = surfacepbcs[i]->GetInt("Id of periodic boundary condition");
      std::vector<int> nodeids(*(surfacepbcs[i]->Nodes()));
      for (size_t j = 0; j < surfacepbcs.size(); ++j)
      {
        const int slaveid = surfacepbcs[j]->GetInt("Id of periodic boundary condition");
        if (masterid == slaveid)
        {
          const std::string* isslave =
              surfacepbcs[j]->Get<std::string>("Is slave periodic boundary condition");
          if (*isslave == "Slave")
          {
            const std::vector<int>* slavenodeids = surfacepbcs[j]->Nodes();
            // append slave node Ids to node Ids for the complete condition
            for (size_t k = 0; k < slavenodeids->size(); ++k)
              nodeids.push_back(slavenodeids->at(k));
          }
        }
      }

      // Get normal direction of pbc plane
      const std::string* pbcplane =
          surfacepbcs[i]->Get<std::string>("degrees of freedom for the pbc plane");
      if (*pbcplane == "yz")
        planenormal.push_back(0);
      else if (*pbcplane == "xz")
        planenormal.push_back(1);
      else if (*pbcplane == "xy")
        planenormal.push_back(2);
      else
        dserror("A PBC condition could not provide a plane normal.");

      double min = +10e19;
      double max = -10e19;
      for (size_t j = 0; j < nodeids.size(); ++j)
      {
        const int gid = nodeids[j];
        const int lid = discret_->NodeRowMap()->LID(gid);
        if (lid < 0) continue;
        const DRT::Node* lnode = discret_->lRowNode(lid);
        const double* coord = lnode->X();
        if (coord[planenormal.back()] < min) min = coord[planenormal.back()];
        if (coord[planenormal.back()] > max) max = coord[planenormal.back()];
      }
      globalmins.resize(planenormal.size());
      globalmaxs.resize(planenormal.size());
      discret_->Comm().MinAll(&min, &(globalmins.back()), 1);
      discret_->Comm().MaxAll(&max, &(globalmaxs.back()), 1);
    }
  }  // end loop over all surfacepbcs


  //=======================================================================
  // Create a vector of eleGIDs and a vector of those eles' node coords and
  // redundantly store it on each proc
  //=======================================================================
  std::vector<int> allcuteleids;
  std::vector<double> allnodecoords;
  {
    // Here we simply take the eleids from the boundaryIntCells map, which leads to our list of cut
    // elements also there is no distribution necessary, as this map is already stored on every proc
    for (std::map<int, GEO::BoundaryIntCells>::const_iterator elepatches = interface.begin();
         elepatches != interface.end(); ++elepatches)
      allcuteleids.push_back(elepatches->first);

    // our local nodecoords
    std::vector<double> nodecoords(3 * numnodesperele * (allcuteleids.size()), 0.0);
    allnodecoords.resize(nodecoords.size(), 0.0);

    // write the node coordinates of every cut rownode of this proc into nodecoords
    for (size_t ivec = 0; ivec < allcuteleids.size(); ++ivec)
    {
      int elegid = allcuteleids[ivec];
      int elelid = discret_->ElementRowMap()->LID(elegid);
      if (elelid >= 0)
      {
        const int coordbase = 3 * numnodesperele * ivec;
        const DRT::Element* ele = discret_->lRowElement(elelid);
        const DRT::Node* const* nodes = ele->Nodes();
        for (int inode = 0; inode < ele->NumNode(); ++inode)
        {
          const int nodecoordbase = coordbase + 3 * inode;
          nodecoords[nodecoordbase + 0] = nodes[inode]->X()[0];
          nodecoords[nodecoordbase + 1] = nodes[inode]->X()[1];
          nodecoords[nodecoordbase + 2] = nodes[inode]->X()[2];
        }
      }
    }

    discret_->Comm().SumAll(&(nodecoords[0]), &(allnodecoords[0]), (int)nodecoords.size());
  }

  //================================================================
  // loop all row nodes on the processor
  // those nodes will receive new phi values
  //================================================================
  for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); ++lnodeid)
  {
    // get the processor local node
    const DRT::Node* lnode = discret_->lRowNode(lnodeid);

    // get the dof associated with this node
    const int dofgid =
        discret_->Dof(0, lnode, 0);  // since this is a scalar field the dof is always 0
    int doflid = dofrowmap->LID(dofgid);
    if (doflid < 0) dserror("Proc %d: Cannot find dof gid=%d in Epetra_Vector", myrank_, dofgid);

    // get physical coordinates of this node
    LINALG::Matrix<3, 1> nodecoord(false);
    nodecoord(0) = lnode->X()[0];
    nodecoord(1) = lnode->X()[1];
    nodecoord(2) = lnode->X()[2];

    //=======================================================================================
    // Build a list< pair< int eleGID, double distance > >
    // the distance is based on the distance between the current node and the closest node of
    // the cut element. This guarantees an estimated distance <= the real distance
    //=======================================================================================
    std::list<std::pair<int, double>> eledistance;

    {
      // loop all cut elements
      for (size_t ieleid = 0; ieleid < allcuteleids.size(); ++ieleid)
      {
        const size_t coordbase = 3 * numnodesperele * ieleid;
        double distance = 1.0e19;

        // loop all cut element's nodes
        for (int inode = 0; inode < numnodesperele; ++inode)
        {
          const int nodecoordbase = coordbase + 3 * inode;
          LINALG::Matrix<3, 1> delta(false);
          delta(0) = allnodecoords[nodecoordbase + 0];
          delta(1) = allnodecoords[nodecoordbase + 1];
          delta(2) = allnodecoords[nodecoordbase + 2];

          delta.Update(1.0, nodecoord, -1.0);

          // take care of PBCs
          for (size_t ipbc = 0; ipbc < planenormal.size(); ++ipbc)
          {
            const double fulllength = (globalmaxs[ipbc] - globalmins[ipbc]);
            if (delta(planenormal[ipbc]) >= fulllength / 2.0)
              delta(planenormal[ipbc]) = fulllength - delta(planenormal[ipbc]);
            else if (delta(planenormal[ipbc]) <= -fulllength / 2.0)
              delta(planenormal[ipbc]) = delta(planenormal[ipbc]) + fulllength;
          }
          const double thisdistance =
              sqrt(delta(0) * delta(0) + delta(1) * delta(1) + delta(2) * delta(2));

          if (thisdistance < distance) distance = thisdistance;
        }

        std::pair<int, double> thispair;
        thispair.first = allcuteleids[ieleid];
        thispair.second = distance;
        eledistance.push_back(thispair);
      }
    }
    if (eledistance.empty()) dserror("No intersected elements available! G-function correct?");


    //==================================================================
    // sort the the vector in ascending order by the estimated distance
    //==================================================================
    // this is the STL sorting, which is pretty fast
    eledistance.sort(MyComparePairs);

    //--------------------------------------------------------------------------------
    // if a reinitbandwith is used the nodes not within the band will be set to the
    // estimated distance for all others the actual distance will be determined
    //--------------------------------------------------------------------------------
    if (!reinitband_ or (reinitband_ and fabs(eledistance.front().second) <= reinitbandwidth_))
    {
      //========================================================================
      // + update the eledistance vector with the real distance to the interface
      //   starting with the closest estimated element.
      // + Sort the vector by distance after every iteration.
      // + if the distance of the first element in the vector does not change
      //   any more, we have found the shortest distance
      //========================================================================
      int oldeleid = -1;
      while (oldeleid !=
             eledistance.front()
                 .first)  // this is just a safety check. usually loop should abort earlier.
      {
        oldeleid = eledistance.front().first;

        // the minimal distance, if all element patches and the PBCs are considered
        double pbcmindist = 1.0e19;

        // get patches belonging to first entry
        std::map<int, GEO::BoundaryIntCells>::const_iterator elepatches =
            interface.find(eledistance.front().first);
        if (elepatches == interface.end())
          dserror("Could not find the boundary integration cells belonging to Element %d.",
              eledistance.front().first);

        // number of flamefront patches for this element
        const std::vector<GEO::BoundaryIntCell> patches = elepatches->second;
        const int numpatch = patches.size();

        //--------------------------------------------------------------------
        // due to the PBCs the node might actually be closer to the
        // interface then would be calculated if one only considered
        // the actual position of the node. In order to find the
        // smallest distance the node is copied along all PBC directions
        //
        //   +------------------+ - - - - - - - - - -+
        //   +             II   +
        //   +   x        I  I  +    y               +
        //   +             II   +
        //   +------------------+ - - - - - - - - - -+
        //         original           copy
        //
        //   x: current node
        //   y: copy of current node
        //   I: interface
        //   +: pbc
        //--------------------------------------------------------------------
        if (planenormal.size() > 3)
          dserror(
              "Sorry, but currently a maximum of three periodic boundary conditions are supported "
              "by the combustion reinitializer.");

        // since there is no stl pow(INT, INT) function, we calculate it manually
        size_t looplimit = 1;
        for (size_t i = 0; i < planenormal.size(); ++i) looplimit *= 2;

        for (size_t ipbc = 0; ipbc < looplimit; ++ipbc)
        {
          double mindist = 1.0e19;
          LINALG::Matrix<3, 1> tmpcoord(nodecoord);

          // determine which pbcs have to be applied
          //
          // loopcounter | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
          // ------------+---+---+---+---+---+---+---+---+
          //  first PBC  |     x       x       x       x
          // second PBC  |         x   x           x   x
          //  third PBC  |                 x   x   x   x
          //
          // this is equivalent to the binary representation
          // of the size_t
          if (ipbc & 0x01)
          {
            const double pbclength = globalmaxs[0] - globalmins[0];
            if (nodecoord(0) < globalmins[0] or nodecoord(0) > globalmaxs[0]) continue;
            if (nodecoord(planenormal[0]) > globalmins[0] + pbclength / 2.0)
              tmpcoord(planenormal[0]) = nodecoord(planenormal[0]) - pbclength;
            else
              tmpcoord(planenormal[0]) = nodecoord(planenormal[0]) + pbclength;
          }
          if (ipbc & 0x02)
          {
            const double pbclength = globalmaxs[1] - globalmins[1];
            if (nodecoord(1) < globalmins[1] or nodecoord(1) > globalmaxs[1]) continue;
            if (nodecoord(planenormal[1]) > globalmins[1] + pbclength / 2.0)
              tmpcoord(planenormal[1]) = nodecoord(planenormal[1]) - pbclength;
            else
              tmpcoord(planenormal[1]) = nodecoord(planenormal[1]) + pbclength;
          }
          if (ipbc & 0x04)
          {
            const double pbclength = globalmaxs[2] - globalmins[2];
            if (nodecoord(2) < globalmins[2] or nodecoord(2) > globalmaxs[2]) continue;
            if (nodecoord(planenormal[2]) > globalmins[2] + pbclength / 2.0)
              tmpcoord(planenormal[2]) = nodecoord(planenormal[2]) - pbclength;
            else
              tmpcoord(planenormal[2]) = nodecoord(planenormal[2]) + pbclength;
          }

          //-----------------------------------------
          // loop flame front patches of this element
          //-----------------------------------------
          for (int ipatch = 0; ipatch < numpatch; ++ipatch)
          {
            // get a single patch from group of flamefront patches
            const GEO::BoundaryIntCell patch = patches[ipatch];

            // only triangles and quadrangles are allowed as flame front patches (boundary cells)
            if (!(patch.Shape() == DRT::Element::tri3 or patch.Shape() == DRT::Element::quad4))
            {
              dserror("invalid type of boundary integration cell for reinitialization");
            }

            // get coordinates of vertices defining flame front patch
            const LINALG::SerialDenseMatrix& patchcoord = patch.CellNodalPosXYZ();

            // compute normal vector to flame front patch
            LINALG::Matrix<3, 1> normal(true);
            ComputeNormalVectorToInterface(patch, patchcoord, normal);

            //-----------------------------------------
            // find flame front patches facing the node
            //-----------------------------------------
            // boolean indicating if facing patch was found
            bool facenode = false;
            // distance to the facing patch
            double patchdist = 1.0e19;  // default value
            // check if this patch faces the node
            FindFacingPatchProjCellSpace(tmpcoord, patch, patchcoord, normal, facenode, patchdist);

            // a facing patch was found
            if (facenode == true)
            {
              // overwrite smallest distance if computed patch distance is smaller
              if (fabs(patchdist) < fabs(mindist))
              {
                // if G-value at the node is negative, the minimal distance has to be negative
                if ((*phinp_)[doflid] < 0.0)
                  mindist = -patchdist;
                else
                  mindist = patchdist;
              }
            }

            //-------------------------------------------------------------
            // compute smallest distance to edges of this flame front patch
            //-------------------------------------------------------------
            // distance to the patch edge
            double edgedist = 1.0e19;
            ComputeDistanceToEdge(tmpcoord, patch, patchcoord, edgedist);

            if (fabs(edgedist) < fabs(mindist))
            {
              // if G-value at the node is negative, the minimal distance has to be negative
              if ((*phinp_)[doflid] < 0.0)
                mindist = -edgedist;
              else
                mindist = edgedist;
            }

            //----------------------------------------------------------------
            // compute smallest distance to vertices of this flame front patch
            //----------------------------------------------------------------
            // distance to the patch vertex
            double vertexdist = 1.0e19;
            ComputeDistanceToPatch(tmpcoord, patch, patchcoord, vertexdist);

            if (fabs(vertexdist) < fabs(mindist))
            {
              // if G-value at the node is negative, the minimal distance has to be negative
              if ((*phinp_)[doflid] < 0.0)
                mindist = -vertexdist;
              else
                mindist = vertexdist;
            }
          }  // loop over flamefront patches

          if (fabs(mindist) < fabs(pbcmindist))
          {
            pbcmindist = mindist;
          }
        }  // loop over PBCs

        // store the new distance, which is >= the estimated distance
        eledistance.front().second = pbcmindist;


        //==============================================================
        // sort the the vector in ascending order by the distance
        //==============================================================
        // here we use the fact, that everything is already sorted but the first list item
        std::pair<int, double> tmppair = eledistance.front();
        std::list<std::pair<int, double>>::iterator insertiter = eledistance.begin();
        ++insertiter;

        int loopcount = 0;
        // find the place where the item must be inserted
        // while (fabs(tmppair.second) > fabs(insertiter->second) and insertiter !=
        // eledistance.end())
        while (insertiter != eledistance.end() and fabs(tmppair.second) > fabs(insertiter->second))
        {
          insertiter++;
          loopcount++;
        }

        // this removes the item from the front and inserts it where insertiter points to
        eledistance.splice(insertiter, eledistance, eledistance.begin());

        // if item was inserted at the beginnig of the list, it must be the shortest distance
        // possible and we can stop checking the other elements' distances
        if (loopcount == 0) break;

      }  // loop over eledistance
    }
    // if outside the reinit band
    else
    {
      // correct the sign of estimated distance
      if ((*phinp_)[doflid] < 0.0) eledistance.front().second = -eledistance.front().second;
    }

    int err = phinp_->ReplaceMyValues(1, &(eledistance.front().second), &doflid);
    if (err) dserror("this did not work");
  }

  if (myrank_ == 0) std::cout << " done" << std::endl;

  return;
}


/*--------------------------------------------------------------------- -----------------*
 | find a facing flame front patch by projection of node into boundary cell space        |
 |                                                                           henke 12/09 |
 *----------------------------------------------------------------------  -------------- */
void SCATRA::LevelSetAlgorithm::FindFacingPatchProjCellSpace(const LINALG::Matrix<3, 1>& node,
    const GEO::BoundaryIntCell& patch, const LINALG::SerialDenseMatrix& patchcoord,
    const LINALG::Matrix<3, 1>& normal, bool& facenode, double& patchdist)
{
  // indicator
  facenode = false;

  static LINALG::Matrix<2, 1> eta(true);
  double alpha = 0.0;

  //-------------------------------------------------------
  // perform Newton-Raphson method to project node on patch
  //-------------------------------------------------------
  bool converged = false;
  switch (patch.Shape())
  {
    case DRT::Element::tri3:
    {
      converged =
          ProjectNodeOnPatch<DRT::Element::tri3>(node, patch, patchcoord, normal, eta, alpha);
      break;
    }
    case DRT::Element::quad4:
    {
      converged =
          ProjectNodeOnPatch<DRT::Element::quad4>(node, patch, patchcoord, normal, eta, alpha);
      break;
    }
    default:
    {
      dserror("unknown type of boundary integration cell");
      break;
    }
  }

  // Newton iteration converged
  //  std::cout << "Newton iteration converged in " << iter << " steps!" << std::endl;

  //----------------------------------------------------
  // check if projection lies within boundary cell space
  //----------------------------------------------------
  // remark: - tolerance has to be of same order as the tolerance that coordinates of projected
  // nodes
  //           differ from an exact position on edges of patches (e.g. 1.0E-7 ~ 1.0E-8 -> 1.0E-6)
  //         - if this is not the case, the level set function can become tilted, since valid
  //           patches are ignored
  double TOL = 1e-6;

  switch (patch.Shape())
  {
    case DRT::Element::tri3:
    {
      // criteria for tri3 patch
      if ((eta(0) > -TOL) and (eta(0) < 1.0 + TOL) and (eta(1) > -TOL) and (eta(1) < 1.0 + TOL) and
          (1.0 - eta(0) - eta(1) > -TOL) and (1.0 - eta(0) - eta(1) < 1.0 + TOL) and converged)
      {
        facenode = true;
        patchdist = fabs(alpha);
        //      std::cout << "facing patch found (tri3 patch)! coordinates eta(0): " << eta(0) << "
        //      eta(1) " << eta(1) << std::endl;
      }
      break;
    }
    case DRT::Element::quad4:
    {
      // criteria for quad4 patch
      if ((eta(0) > -1.0 - TOL) and (eta(0) < 1.0 + TOL) and (eta(1) > -1.0 - TOL) and
          (eta(1) < 1.0 + TOL) and converged)
      {
        facenode = true;
        patchdist = fabs(alpha);
        //      std::cout << "facing patch found (quad4 patch)!" << std::endl;
      }
      break;
    }
    default:
    {
      dserror("unknown type of boundary integration cell");
      break;
    }
  }
  //  if (!converged)
  //  {
  //    std::cout << "node x component " << node(0,0) << std::endl;
  //    std::cout << "node y component " << node(1,0) << std::endl;
  //    std::cout << "node z component " << node(2,0) << std::endl;
  //    std::cout << "eta1 " << eta(0) << std::endl;
  //    std::cout << "eta2 " << eta(1) << std::endl;
  //    std::cout << "alpha " << alpha << std::endl;
  //    std::cout << "patch vertices x component " << patchcoord(0,0) << " " << patchcoord(0,1) << "
  //    " << patchcoord(0,2) << std::endl; std::cout << "patch vertices y component " <<
  //    patchcoord(1,0) << " " << patchcoord(1,1) << " " << patchcoord(1,2) << std::endl; std::cout
  //    << "patch vertices z component " << patchcoord(2,0) << " " << patchcoord(2,1) << " " <<
  //    patchcoord(2,2) << std::endl;
  //  }

  return;
}


/*---------------------------------------------------------------------------------------*
 | compute distance to edge of patch                                         henke 08/09 |
 *-------------------------------------------------------------------------------------- */
void SCATRA::LevelSetAlgorithm::ComputeDistanceToEdge(const LINALG::Matrix<3, 1>& node,
    const GEO::BoundaryIntCell& patch, const LINALG::SerialDenseMatrix& patchcoord,
    double& edgedist)
{
  // set temporary edgedist to large value
  double edgedisttmp = edgedist;

  // number of vertices of flame front patch (3 for tri3; 4 for quad4)
  const size_t numvertices = patchcoord.N();

  // current vertex of the patch (first vertex)
  static LINALG::Matrix<3, 1> vertex1(true);
  // current next vertex of the patch (second vertex)
  static LINALG::Matrix<3, 1> vertex2(true);
  // distance vector from first vertex to node
  static LINALG::Matrix<3, 1> vertex1tonode(true);
  // distance vector from first vertex to second vertex
  static LINALG::Matrix<3, 1> vertex1tovertex2(true);

  // compute distance to all vertices of patch
  for (size_t ivert = 0; ivert < numvertices; ++ivert)
  {
    // vertex1 of flame front patch
    vertex1(0) = patchcoord(0, ivert);
    vertex1(1) = patchcoord(1, ivert);
    vertex1(2) = patchcoord(2, ivert);

    if (ivert < (numvertices - 1))
    {
      vertex2(0) = patchcoord(0, ivert + 1);
      vertex2(1) = patchcoord(1, ivert + 1);
      vertex2(2) = patchcoord(2, ivert + 1);
    }
    else if (ivert == (numvertices - 1))
    {
      vertex2(0) = patchcoord(0, 0);
      vertex2(1) = patchcoord(1, 0);
      vertex2(2) = patchcoord(2, 0);
    }

    // compute distance vector from node to current first
    vertex1tonode.Update(1.0, node, -1.0, vertex1);
    // compute distance vector from current second first vertex to current frist vertex (edge)
    vertex1tovertex2.Update(1.0, vertex2, -1.0, vertex1);
    double normvertex1tovertex2 = vertex1tovertex2.Norm2();
    // normalize vector
    vertex1tovertex2.Scale(1.0 / normvertex1tovertex2);

    // scalar product of vertex1tonode and the normed vertex1tovertex2
    double lotfusspointdist = vertex1tovertex2.Dot(vertex1tonode);

    if ((lotfusspointdist >= 0.0) and
        (lotfusspointdist <= normvertex1tovertex2))  // lotfusspoint on edge
    {
      LINALG::Matrix<3, 1> lotfusspoint(true);
      lotfusspoint.Update(1.0, vertex1, lotfusspointdist, vertex1tovertex2);
      LINALG::Matrix<3, 1> nodetolotfusspoint(true);
      nodetolotfusspoint.Update(1.0, lotfusspoint, -1.0, node);

      // determine length of vector from node to lot fuss point
      edgedisttmp = nodetolotfusspoint.Norm2();
      if (edgedisttmp < edgedist) edgedist = edgedisttmp;
    }
  }

  return;
}


/*---------------------------------------------- ----------------------------------------*
 | compute distance to vertex of patch                                       henke 08/09 |
 *-------------------------------------------------------------------------------------- */
void SCATRA::LevelSetAlgorithm::ComputeDistanceToPatch(const LINALG::Matrix<3, 1>& node,
    const GEO::BoundaryIntCell& patch, const LINALG::SerialDenseMatrix& patchcoord,
    double& vertexdist)
{
  // set temporary vertexdist to large value
  double vertexdisttmp = vertexdist;

  // number of vertices of flame front patch (3 for tri3; 4 for quad4)
  const size_t numvertices = patchcoord.N();

  // current vertex of the patch
  static LINALG::Matrix<3, 1> vertex(true);
  // distance vector from patch to node
  static LINALG::Matrix<3, 1> dist(true);

  // compute distance to all vertices of patch
  for (size_t ivert = 0; ivert < numvertices; ++ivert)
  {
    // vertex of flame front patch
    vertex(0) = patchcoord(0, ivert);
    vertex(1) = patchcoord(1, ivert);
    vertex(2) = patchcoord(2, ivert);

    // compute distance vector from flame front to node
    dist.Update(1.0, node, -1.0, vertex);

    // compute L2-norm of distance vector
    vertexdisttmp = dist.Norm2();
    if (vertexdisttmp < vertexdist) vertexdist = vertexdisttmp;
  }

  return;
}


/*-------------------------------------------------------------------------------------*
 | compute normal vector to interface patch                                henke 08/09 |
 *------------------------------------------------- ---------------------------------- */
void SCATRA::LevelSetAlgorithm::ComputeNormalVectorToInterface(const GEO::BoundaryIntCell& patch,
    const LINALG::SerialDenseMatrix& patchcoord, LINALG::Matrix<3, 1>& normal)
{
  // first point of flame front patch
  LINALG::Matrix<3, 1> point1;
  point1(0) = patchcoord(0, 0);
  point1(1) = patchcoord(1, 0);
  point1(2) = patchcoord(2, 0);

  // second point of flame front patch
  LINALG::Matrix<3, 1> point2;
  point2(0) = patchcoord(0, 1);
  point2(1) = patchcoord(1, 1);
  point2(2) = patchcoord(2, 1);

  // first edge of flame front patch
  LINALG::Matrix<3, 1> edge1;
  edge1.Update(1.0, point2, -1.0, point1);

  // third point of flame front patch
  point2(0) = patchcoord(0, 2);
  point2(1) = patchcoord(1, 2);
  point2(2) = patchcoord(2, 2);

  // second edge of flame front patch (if patch is triangle; if not: edge 2 is secant of polygon)
  LINALG::Matrix<3, 1> edge2;
  edge2.Update(1.0, point2, -1.0, point1);

  // compute normal vector of patch (cross product: edge1 x edge2)
  // remark: normal vector points into unburnt domain (G<0)
  normal(0) = (edge1(1) * edge2(2) - edge1(2) * edge2(1));
  normal(1) = (edge1(2) * edge2(0) - edge1(0) * edge2(2));
  normal(2) = (edge1(0) * edge2(1) - edge1(1) * edge2(0));

  //  const Epetra_Comm& comm = scatra_.Discretization()->Comm();
  //  std::cout << "proc " << comm.MyPID() << " normal " <<  normal << std::endl;
  //  std::cout << "proc " << comm.MyPID() << " patch " <<  patchcoord << std::endl;

  // compute unit (normed) normal vector
  double norm = sqrt(normal(0) * normal(0) + normal(1) * normal(1) + normal(2) * normal(2));
  if (norm == 0.0) dserror("norm of normal vector is zero!");
  normal.Scale(1.0 / norm);

  return;
}


/*-------------------------------------------------------------------------------------*
 | project node into the boundary cell space                               henke 08/09 |
 *------------------------------------------------- ---------------------------------- */
template <DRT::Element::DiscretizationType DISTYPE>
bool SCATRA::LevelSetAlgorithm::ProjectNodeOnPatch(const LINALG::Matrix<3, 1>& node,
    const GEO::BoundaryIntCell& patch, const LINALG::SerialDenseMatrix& patchcoord,
    const LINALG::Matrix<3, 1>& normal, LINALG::Matrix<2, 1>& eta, double& alpha)
{
  // indicator for convergence of Newton-Raphson scheme
  bool converged = false;
  // number space dimensions for 3d combustion problems
  const size_t nsd = 3;
  // here, a triangular boundary integration cell is assumed (numvertices = 3)
  const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // get coordinates of vertices of flame front patch
  // remark: here we only get a view (bool true) on the SerialDenseMatrix returned by
  // CellNodalPosXYZ()
  LINALG::Matrix<nsd, numvertices> patchcoordfix(patchcoord.A(), true);

  static LINALG::Matrix<numvertices, 1> funct(true);
  static LINALG::Matrix<2, numvertices> deriv(true);
  static LINALG::Matrix<nsd, 1> projX(true);
  static LINALG::Matrix<nsd, 2> gradprojX(true);

  //----------------------------------
  // start values for iterative scheme
  //----------------------------------
  // start position (barycenter of triangular boundary cell)
  eta(0) = 1.0 / 3.0;
  eta(1) = 1.0 / 3.0;
  // auxiliary variable
  // remark: third unknown to close system of equations; arbitrary value
  alpha = 0.0;

  // function F (system of equations)
  static LINALG::Matrix<nsd, 1> f(true);
  // gradient of function F (dF/deta(0), dF/deta(1), dF/dalpha)
  static LINALG::Matrix<nsd, nsd> gradf(true);
  // increment in Newton iteration (unknown to be solved for)
  static LINALG::Matrix<nsd, 1> incr(true);

  // maximum number Newton iterations
  size_t maxiter = 3;
  // convergence tolerance
  double conv = 0.0;

  //------------------------------------------------------
  // Newton-Raphson loop for non-linear projection problem
  //------------------------------------------------------
  for (size_t iter = 0; iter < maxiter; ++iter)
  {
    // evaluate shape functions in boundary cell space at current position \eta_1,\eta_2 on the
    // patch
    funct.Clear();
    DRT::UTILS::shape_function_2D(funct, eta(0), eta(1), patch.Shape());
    // evaluate derivatives of shape functions in boundary cell space at current position
    // \eta_1,\eta_2 on the patch
    deriv.Clear();
    DRT::UTILS::shape_function_2D_deriv1(deriv, eta(0), eta(1), patch.Shape());

    // evaluate projection X of node P at current position \eta_1,\eta_2 on the patch
    // projX(i,j) = patchcoord(i,k)*funct(k,1)
    projX.Clear();
    projX.MultiplyNN(patchcoordfix, funct);

    // evaluate gradient of projection X of node P at current position \eta_1,\eta_2 on the patch
    // gradprojX(i,j) = patchcoord(i,k)*deriv(j,k)
    gradprojX.Clear();
    gradprojX.MultiplyNT(patchcoordfix, deriv);

    //---------------------------------------------------
    // build system of equations F and its gradient gradF
    //---------------------------------------------------
    // TODO documentaton missing
    f.Clear();
    gradf.Clear();
    incr.Clear();
    for (size_t icoord = 0; icoord < nsd; ++icoord)
    {
      // evaluate function f
      f(icoord) = projX(icoord) + alpha * normal(icoord) - node(icoord);
      // evaluate gradient of function at current position on patch in boundary cell space
      gradf(icoord, 0) = gradprojX(icoord, 0);
      gradf(icoord, 1) = gradprojX(icoord, 1);
      gradf(icoord, 2) = normal(icoord);
    }

    // check convergence
    conv = sqrt(f(0) * f(0) + f(1) * f(1) + f(2) * f(2));
    // std::cout << "iteration " << iter << ": -> |f|=" << conv << std::endl;
    if (conv <= 1.0E-12) break;

    //----------------------------------------------------
    // solve linear system of equations: gradF * incr = -F
    //----------------------------------------------------
    // F = F*-1.0
    f.Scale(-1.0);
    // solve A.X=B
    LINALG::FixedSizeSerialDenseSolver<nsd, nsd, 1> solver;
    solver.SetMatrix(gradf);               // set A=gradF
    solver.SetVectors(incr, f);            // set X=incr, B=F
    solver.FactorWithEquilibration(true);  // "some easy type of preconditioning" (Michael)
    int err2 = solver.Factor();            // ?
    int err = solver.Solve();              // incr = gradF^-1.F
    if ((err != 0) || (err2 != 0))
      dserror("solving linear system in Newton-Raphson method for projection failed");

    // update eta and alpha
    eta(0) += incr(0);
    eta(1) += incr(1);
    alpha += incr(2);
    // std::cout << "solution vector: component 1: " << eta(0) << " component 2: " << eta(1) << "
    // alpha: " << alpha << std::endl;
  }
  // change sign to preserve sign of G-function
  alpha = -alpha;

  // Newton iteration unconverged
  if (conv > 1.0E-12)
  {
    alpha = 7777.7;
    //        std::cout << "projection did not converge" << std::endl;
    // dserror("projection did not converge!");
  }
  else
  {
    converged = true;
    // std::cout << "convergence criterion " << conv << std::endl;
    // std::cout << "solution vector: component 1: " << eta(0) << " component 2: " << eta(1) << "
    // alpha: " << alpha << std::endl;
  }

  return converged;
}


/*------------------------------------------------------------------------------------------------*
 | correct the volume of the minus domain after reinitialization                  rasthofer 07/11 |
 |                                                                                    DA wichmann |
 | Idea: shift level-set so that volume is conserved                                              |
 *------------------------------------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::CorrectVolume()
{
  double volminus = 0.0;
  double volplus = 0.0;
  double surface = 0.0;
  std::map<int, GEO::BoundaryIntCells> interface;
  interface.clear();
  // reconstruct interface and calculate volumes, etc ...
  SCATRA::LEVELSET::Intersection intersect;
  intersect.CaptureZeroLevelSet(phinp_, discret_, volminus, volplus, surface, interface);

  const double voldelta = initvolminus_ - volminus;
  if (myrank_ == 0)
    IO::cout << "Correcting volume of minus(-) domain by " << voldelta << " ... " << IO::endl;

  // This is a guess on how thick a layer needs to be added to the surface of the minus domain.
  // Due to $ \grad \phi \approx 1 $ this also happens to be the value that needs to be subtracted
  // of all phis. To make sure that \grad \phi really is close to 1 this function should only be
  // called after a reinitialization.
  const double thickness = -voldelta / surface;

  Teuchos::RCP<Epetra_Vector> one = Teuchos::rcp(new Epetra_Vector(phin_->Map()));
  one->PutScalar(1.0);

  // update phi
  phinp_->Update(thickness, *one, 1.0);

  if (myrank_ == 0) IO::cout << "done" << IO::endl;

  return;
}

/*----------------------------------------------------------------------*
 | elliptic reinitialization                            rasthofer 09/14 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::ReinitElliptic(std::map<int, GEO::BoundaryIntCells>& interface)
{
  // store interface
  interface_eleq_ = Teuchos::rcp(new std::map<int, GEO::BoundaryIntCells>(interface));

  // call the executing method
  ReinitializeWithEllipticEquation();
}

/*----------------------------------------------------------------------*
 | elliptic reinitialization                            rasthofer 09/14 |
 *----------------------------------------------------------------------*/
void SCATRA::LevelSetAlgorithm::ReinitializeWithEllipticEquation()
{
  //-------------------------------------------------
  // preparations
  //-------------------------------------------------

  // set switch flag to true to activate reinitialization specific parts
  switchreinit_ = true;

  // set element parameters for reinitialization equation
  SetReinitializationElementParameters();

  // vector for initial phi (solution of level-set equation) of reinitialization process
  // this vector is only initialized: currently function CalcNodeBasedReinitVel() is also
  // used to compute nodal level-set gradients, and this function expects that initialphireinit_ has
  // been set although it is not used for the present purposes
  initialphireinit_ = LINALG::CreateVector(*(discret_->DofRowMap()), true);

  //-------------------------------------------------
  // solve
  //-------------------------------------------------

  // we simply call the LinearSolve (since the elliptic reinitialization equation is
  // indeed nonlinear), and all the rest concerning the correct action type and
  // parameters is handled via the switchreinit_-flag in the concrete time-integration
  // schemes for level-set problems

  // some preparations
  Teuchos::RCP<Epetra_Vector> phinmloc = Teuchos::rcp(new Epetra_Vector(*phinp_));
  Teuchos::RCP<Epetra_Vector> inc = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap()), true));
  int step = 0;
  bool not_conv = true;

  while (not_conv)
  {
    step += 1;

    //-----------------------------
    // compute node-based gradient
    //-----------------------------
    if (projection_) CalcNodeBasedReinitVel();

    //-----------------------------
    // setup and solve system
    //-----------------------------
    // caution: we can only use LinearSolve together with linear-full strategy here
    LinearSolve();

    //-----------------------------
    // check convergence
    //-----------------------------
    inc->Update(1.0, *phinp_, -1.0, *phinmloc, 0.0);
    double norm = 0.0;
    inc->Norm2(&norm);

    if (myrank_ == 0)
      std::cout << "STEP:  " << step << "/" << pseudostepmax_
                << "  -- inc norm L2:  " << std::setprecision(3) << std::scientific << norm
                << std::endl;

    if (reinit_tol_ > 0.0)
    {
      if (norm < reinit_tol_ or step >= pseudostepmax_) not_conv = false;
    }
    else
    {
      if (step >= pseudostepmax_) not_conv = false;
    }

    phinmloc->Update(1.0, *phinp_, 0.0);
  }

  //-------------------------------------------------
  // finish
  //-------------------------------------------------

  // reset time-integration parameters for element evaluation
  // SetElementTimeParameter(); -> have not been modified
  // reset general parameters for element evaluation
  SetElementGeneralParameters();
  SetElementTurbulenceParameters();

  // clear variables
  interface_eleq_ = Teuchos::null;
  initialphireinit_ = Teuchos::null;
  if (projection_ == true) nb_grad_val_->PutScalar(0.0);

  return;
}
