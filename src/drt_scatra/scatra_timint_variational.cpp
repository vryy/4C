/*----------------------------------------------------------------------*/
/*! \file

\brief  Time integration for variational formulation problems

\level 2

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/
#include "scatra_timint_variational.H"
// To read initial conditions from input file
#include "../drt_lib/drt_globalproblem.H"
// To extract the material
#include "../drt_lib/drt_node.H"
#include "../drt_mat/scatra_mat_var_chemdiffusion.H"
#include "../drt_nurbs_discret/drt_apply_nurbs_initial_condition.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
// For preconditioning
#include "../drt_scatra/scatra_timint_meshtying_strategy_std_var_chemdiff.H"
// To output analytic solution
#include "../drt_io/io.H"
// To use SCATRA::calc_error when computing analytic solution
#include "../drt_scatra_ele/scatra_ele_action.H"
// For initialization of splitter
#include "../drt_fluid/fluid_utils.H"  // for splitter



/*----------------------------------------------------------------------*
 |  Constructor (public)                                   deanda 10/17 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntVariational::TimIntVariational(
    Teuchos::RCP<DRT::Discretization> actdis,          //!< discretization
    Teuchos::RCP<LINALG::Solver> solver,               //!< linear solver
    Teuchos::RCP<Teuchos::ParameterList> params,       //!< parameter list
    Teuchos::RCP<Teuchos::ParameterList> extraparams,  //!< supplementary parameter list
    Teuchos::RCP<IO::DiscretizationWriter> output,     //!< output writer
    const int probnum                                  //!< global problem number
    )
    : ScaTraTimIntImpl(actdis, solver, params, extraparams, output, probnum),
      varparams_(params),
      phi0_(Teuchos::null)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                             deanda 10/17 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVariational::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  ScaTraTimIntImpl::Init();

  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                             deanda 10/17 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVariational::Setup()
{
  // call Setup()-functions of base classes
  // note: this order is important
  // back-up of the initial condition
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  phi0_ = LINALG::CreateVector(*dofrowmap, true);
  ScaTraTimIntImpl::Setup();

  if ((calcerror_ == INPAR::SCATRA::calcerror_AnalyticSeries && nsd_ == 1) ||
      (calcerror_ == INPAR::SCATRA::calcerror_byfunction))
    phiAnalytic_ = LINALG::CreateVector(*dofrowmap, true);

  return;
}

/*----------------------------------------------------------------------*
 |  initialize splitter for preconditioning                deanda 10/17 |
 *----------------------------------------------------------------------*/

void SCATRA::TimIntVariational::SetupSplitter()
{
  // set up the concentration-ch.potential splitter
  // prepare sets for concentration and potential dofs
  std::set<int> conddofset;
  std::set<int> otherdofset;

  // fill sets
  for (int inode = 0; inode < discret_->NumMyRowNodes(); ++inode)
  {
    std::vector<int> dofs = discret_->Dof(0, discret_->lRowNode(inode));
    for (unsigned idof = 0; idof < dofs.size(); ++idof)
      if (idof < (unsigned)NumScal())
        otherdofset.insert(dofs[idof]);
      else
        conddofset.insert(dofs[idof]);
  }

  // transform sets to maps
  std::vector<int> conddofmapvec(conddofset.begin(), conddofset.end());
  const Teuchos::RCP<const Epetra_Map> conddofmap = Teuchos::rcp(
      new Epetra_Map(-1, conddofmapvec.size(), &conddofmapvec[0], 0, discret_->Comm()));
  std::vector<int> otherdofmapvec(otherdofset.begin(), otherdofset.end());
  const Teuchos::RCP<const Epetra_Map> otherdofmap = Teuchos::rcp(
      new Epetra_Map(-1, otherdofmapvec.size(), &otherdofmapvec[0], 0, discret_->Comm()));

  // set up concentration-potential splitter
  splitter_ =
      Teuchos::rcp(new LINALG::MapExtractor(*discret_->DofRowMap(), conddofmap, otherdofmap));

  //  FLD::UTILS::SetupFluidSplit(*discret_,NumScal(),*splitter_);
}
/*----------------------------------------------------------------------*
 |  adds initial state for reference                       deanda 09/17 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVariational::AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver)
{
  // call base class routine
  ScaTraTimIntImpl::AddTimeIntegrationSpecificVectors(forcedincrementalsolver);

  discret_->SetState("phi0", phi0_);

  if ((calcerror_ == INPAR::SCATRA::calcerror_AnalyticSeries && nsd_ == 1) ||
      (calcerror_ == INPAR::SCATRA::calcerror_byfunction))
    discret_->SetState("phiAnalytic", phiAnalytic_);

  return;
}

/*----------------------------------------------------------------------*
| Destructor dtor (public)                                 deanda 10/17 |
*-----------------------------------------------------------------------*/
SCATRA::TimIntVariational::~TimIntVariational() { return; }


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                         deanda 10/17 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVariational::Update(const int num)
{
  ScaTraTimIntImpl::Update(num);

  return;
}


/*----------------------------------------------------------------------*
 | explicit predictor for nonlinear solver                 deanda 10/17 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVariational::ExplicitPredictor() const
{
  // call base class routine
  ScaTraTimIntImpl::ExplicitPredictor();

  return;
}

/*----------------------------------------------------------------------*
 |  set initial field for phi                              deanda 10/17 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVariational::SetInitialField(
    const INPAR::SCATRA::InitialField init, const int startfuncno)
{
  switch (init)
  {
    case INPAR::SCATRA::initfield_zero_field:
    {
      phi0_->PutScalar(0.0);
      break;
    }
    case INPAR::SCATRA::initialfield_algebraic_field_dependence:
    {
      // Loop through NodeRowMap
      for (int n = 0; n < discret_->NodeRowMap()->NumMyElements(); ++n)
      {
        // Needed variables
        DRT::Node* lnode = discret_->lRowNode(n);
        std::vector<int> node_dofs_gids = discret_->Dof(0, lnode);
        int numdofs = node_dofs_gids.size();
        const int indexOffsetIndependentVariable = numdofs / 2;

        // Safety check-up
        if (node_dofs_gids.size() % 2) dserror("Number of degrees of freedom must be pair!");

        // Gets the material properties needed to work-out the dependence
        DRT::Element** elelist = lnode->Elements();
        // get material from first (arbitrary!) element adjacent to this node
        const Teuchos::RCP<MAT::Material> matlistptr = elelist[0]->Material();
        dsassert(matlistptr->MaterialType() == INPAR::MAT::m_var_chemdiffusion,
            "material is not of type m_var_chemdiffusion");
        const Teuchos::RCP<const MAT::ScatraMatVarChemDiffusion>& actmat =
            Teuchos::rcp_dynamic_cast<const MAT::ScatraMatVarChemDiffusion>(matlistptr);

        // Gets the local id of the d.o.f. for each node
        std::vector<int> node_dof_lids(node_dofs_gids.size());
        for (int dof = 0; dof < numdofs; ++dof)
          node_dof_lids.at(dof) = discret_->DofRowMap()->LID(node_dofs_gids.at(dof));


        // Iterates through the Independent d.o.f. and solves for the dependent d.o.f.
        for (int loc_id = indexOffsetIndependentVariable; loc_id < numdofs; ++loc_id)
        {
          // Assigns value to the independent field from a given function in the input file
          double initialval = problem_->Funct(startfuncno - 1).Evaluate(loc_id, lnode->X(), time_);
          int err = phin_->ReplaceMyValues(1, &initialval, &node_dof_lids.at(loc_id));
          if (err != 0) dserror("d.o.f. not on processor");

          // Assigns value to the dependent field from a given function in the material parameters
          initialval = actmat->ComputeDualInternalEnergy(initialval, actmat->RefMu(),
              actmat->RefC(), actmat->RefTemp() * actmat->GasConstant(), 1);
          int indexDependentVariable = node_dof_lids.at(loc_id - indexOffsetIndependentVariable);
          err = phin_->ReplaceMyValues(1, &initialval, &indexDependentVariable);
          if (err != 0) dserror("d.o.f. not on processor");
        }
      }

      // for NURBS discretizations we have to solve a least squares problem,
      // with high accuracy! (do nothing for Lagrangian polynomials)
      const Teuchos::ParameterList& scatradyn = problem_->ScalarTransportDynamicParams();
      const int lstsolver = scatradyn.get<int>("LINEAR_SOLVER");

      DRT::NURBS::NurbsDiscretization* nurbsdis =
          dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*discret_));
      if (nurbsdis != NULL)
      {
        if (lstsolver == (-1))
          dserror(
              "no linear solver defined for least square NURBS problem. Please set LINEAR_SOLVER "
              "in SCALAR TRANSPORT DYNAMIC to a valid number! Note: this solver block is misused "
              "for the least square problem. Maybe one should add a separate parameter for this.");

        DRT::NURBS::apply_nurbs_initial_condition(
            *discret_, errfile_, problem_->SolverParams(lstsolver), startfuncno, phin_);
      }

      // initialize also the solution vector. These values are a pretty good guess for the
      // solution after the first time step (much better than starting with a zero vector)
      phinp_->Update(1.0, *phin_, 0.0);

      // saves the initial condition for future use
      phi0_->Update(1.0, *phin_, 0.0);

      break;
    }
    default:
    {
      ScaTraTimIntImpl::SetInitialField(init, startfuncno);
      break;
    }
  }  // switch(init)



}  // ScaTraTimIntImpl::SetInitialField


/*----------------------------------------------------------------------------------------*
 | create scalar manager                                                      deanda 10/17 |
 *----------------------------------------------------------------------------------------*/
void SCATRA::TimIntVariational::CreateScalarHandler()
{
  scalarhandler_ = Teuchos::rcp(new ScalarHandlerVar());

  return;
}  // TimIntVariational::CreateScalarHandler()


/*----------------------------------------------------------------------*
 |  calculate error compared to analytical solution         deanda 10/17 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVariational::EvaluateErrorComparedToAnalyticalSol()
{
  switch (calcerror_)
  {
    case INPAR::SCATRA::calcerror_no:  // do nothing (the usual case)
      break;
    case INPAR::SCATRA::calcerror_AnalyticSeries:
    {
      //   References:
      //   John Crank
      //   "Mathematics of Diffusion"
      //    Plane sheet problem
      //   1975, Ed II, 22-24 Eq(2.54)

      if (nsd_ != 1) dserror("The analytic solution has not been defined for this dimension");

      // Needed parameters at the element level for the analytical solution
      // Extracts imposed Dirichlet condition
      Teuchos::RCP<Epetra_Vector> phinpDirichlet = dbcmaps_->ExtractCondVector(phinp_);
      double Cext = phinpDirichlet->operator[](0);  // Special condition since is Dirichlet only at
                                                    // one point Gets length of the domain
      double L = 0;
      Teuchos::RCP<DRT::Discretization> discret = Discretization();
      for (int i = 0; i < discret->NodeRowMap()->NumMyElements(); ++i)
      {
        // Get the coordinates of the mesh
        int gid = discret->NodeRowMap()->GID(i);
        DRT::Node* node = discret->gNode(gid);
        if (node->X()[0] > L) L = node->X()[0];
      }

      // Verify all the parameters were specified in the dat file
      if (!params_->isParameter("CALCERRORNO"))
        dserror(
            "You forgot to specify the number of iterations to be taken by the series. Specify it "
            "using the parameter CALCERRORNO in the datfile");

      // Set the parameters for the error calculation in the element
      Teuchos::ParameterList eleparams;
      const int errorfunctnumber = params_->get<int>("CALCERRORNO");
      if (errorfunctnumber < 1)
        dserror(
            "The value for stopping parameter CALCERRORNO for the series evaluation must be "
            "specified in the datfile and be positive!");
      eleparams.set<int>("error function number", errorfunctnumber);
      eleparams.set<int>("action", SCATRA::calc_error);
      eleparams.set("total time", time_);
      eleparams.set<int>("calcerrorflag", calcerror_);
      eleparams.set<double>("Dirichlet_values", Cext);
      eleparams.set<double>("length_domain", L);

      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("phinp", phinp_);

      // get (squared) error values
      Teuchos::RCP<Epetra_SerialDenseVector> errors = Teuchos::rcp(new Epetra_SerialDenseVector(4));

      discret_->EvaluateScalars(eleparams, errors);
      discret_->ClearState();

      double error_DiffConcentration = 0.0;
      double error_DiffChemPot = 0.0;
      double L2_AnalyticConcentration = 0.0;
      double L2_AnalyticChemPot = 0.0;

      // for the L2 norm, we need the square root
      if (NumScal() < 3)
      {
        error_DiffConcentration = sqrt((*errors)[0]);
        error_DiffChemPot = sqrt((*errors)[1]);
        L2_AnalyticConcentration = sqrt((*errors)[2]);
        L2_AnalyticChemPot = sqrt((*errors)[3]);
      }
      else
        dserror(
            "The analytical solution of the plane sheet has only being coded so far for one "
            "species");

      if (myrank_ == 0 && time_ != 0)
      {
        printf("\nRelative L2_err for plane sheet using series of error functions (time = %f):\n",
            time_);
        printf(" concentration %15.8e\n chemical potential %15.8e\n \n",
            error_DiffConcentration / L2_AnalyticConcentration,
            error_DiffChemPot / L2_AnalyticChemPot);
      }
      if (step_ == 0)
        SaveError2File(0.0, 0.0);
      else
      {
        SaveError2File(error_DiffConcentration / L2_AnalyticConcentration,
            error_DiffChemPot / L2_AnalyticChemPot);
        (*relerrors_)[0] =
            error_DiffConcentration /
            L2_AnalyticConcentration;  // TODO:Update to n species, see for
                                       // referenceSCATRA::ScaTraTimIntImpl::EvaluateErrorComparedToAnalyticalSol()
        (*relerrors_)[1] = error_DiffChemPot / L2_AnalyticChemPot;
      }
      break;
    }
    case INPAR::SCATRA::calcerror_byfunction:
    {
      // create the parameters for the error calculation
      Teuchos::ParameterList eleparams;
      eleparams.set<int>("action", SCATRA::calc_error);
      eleparams.set<int>("calcerrorflag", calcerror_);

      // Get Error function number and past it to the elements
      const int errorfunctnumber = params_->get<int>("CALCERRORNO");
      if (errorfunctnumber < 1)
        dserror("invalid value of parameter CALCERRORNO for error function evaluation!");

      eleparams.set<int>("error function number", errorfunctnumber);

      // provide displacement field in case of ALE
      if (isale_) eleparams.set<int>("ndsdisp", nds_disp_);

      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("phinp", phinp_);

      // get (squared) error values
      Teuchos::RCP<Epetra_SerialDenseVector> errors =
          Teuchos::rcp(new Epetra_SerialDenseVector(4 * NumDofPerNode()));
      discret_->EvaluateScalars(eleparams, errors);
      discret_->ClearState();

      // error parameters
      double error_DiffConcentration = 0.0;
      double error_DiffChemPot = 0.0;
      double L2_AnalyticConcentration = 0.0;
      double L2_AnalyticChemPot = 0.0;

      // for the L2 norm, we need the square root
      error_DiffConcentration = sqrt((*errors)[0]);
      error_DiffChemPot = sqrt((*errors)[1]);
      L2_AnalyticConcentration = sqrt((*errors)[2]);
      L2_AnalyticChemPot = sqrt((*errors)[3]);

      if (myrank_ == 0 && time_ != 0)
      {
        printf("\nRelative L2_error (time = %f):\n", time_);
        printf(" concentration %15.8e\n chemical potential %15.8e\n \n",
            error_DiffConcentration / L2_AnalyticConcentration,
            error_DiffChemPot / L2_AnalyticChemPot);
      }
      if (step_ == 0)
        SaveError2File(0.0, 0.0);
      else
      {
        SaveError2File(error_DiffConcentration / L2_AnalyticConcentration,
            error_DiffChemPot / L2_AnalyticChemPot);
        (*relerrors_)[0] =
            error_DiffConcentration /
            L2_AnalyticConcentration;  // TODO:Update to n species, see for
                                       // referenceSCATRA::ScaTraTimIntImpl::EvaluateErrorComparedToAnalyticalSol()
        (*relerrors_)[1] = error_DiffChemPot / L2_AnalyticChemPot;
      }
      break;
    }
    default:
    {
      // call base class routine
      ScaTraTimIntImpl::EvaluateErrorComparedToAnalyticalSol();
      break;
    }
  }
  return;
}  // SCATRA::TimIntVariational::EvaluateErrorComparedToAnalyticalSol


/*-------------------------------------------------------------------------*
 | valid parameters for the diffusion-conduction formulation  deanda 10/17 |
 *-------------------------------------------------------------------------*/
void SCATRA::TimIntVariational::ValidParameterDiffCond()
{
  if (myrank_ == 0)
  {
    if (DRT::INPUT::IntegralValue<int>(varparams_->sublist("VARIATIONAL"), "BLOCKPRECOND") == true)
      dserror("Block preconditioner is not supported so far!!");
  }

  return;
}
/*----------------------------------------------------------------------*
 |  set dirichlet bc for dependent fields                   deanda 09/17 |
  brief: This part is needed since in the variational formulation only the boundary conditions
  are prescribed for the chemical potential. Based on that when the concentration is also defined
  at the nodes it must be used the Internal energy to relate the concentration at the dirichlet
  boundary to that of the chemical potential
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVariational::ApplyDirichletBC(
    const double time, Teuchos::RCP<Epetra_Vector> phinp, Teuchos::RCP<Epetra_Vector> phidt)
{
  ScaTraTimIntImpl::ApplyDirichletBC(time, phinp, phidt);

  // Extract nodes with Dirichlet condition imposed
  Teuchos::RCP<Epetra_Vector> phinpDirichlet = dbcmaps_->ExtractCondVector(phinp);

  // Safety checks
  if (phinpDirichlet->MyLength() % 2)
    dserror("Number of degrees of freedom with Dirichlet condition must be pair!");
  const int numdofpernode = NumDofPerNode();
  if (numdofpernode != 2)
    dserror("Boundary condition only implemented for two degrees of freedom at the node yet!");

  for (int i = 0; i < phinpDirichlet->MyLength(); i = i + numdofpernode)
  {
    // Locates gid of dependent dof and gets value of independent dof
    int gid_dependent_field = phinpDirichlet->Map().GID(i);
    double value_independent_field = phinpDirichlet->operator[](i + 1);

    // Gets the material parameters and functions needed to define the dependent field
    DRT::Node* lnode = discret_->lRowNode(phinpDirichlet->Map().LID(gid_dependent_field));
    DRT::Element** elelist = lnode->Elements();
    const Teuchos::RCP<MAT::Material> matlistptr = elelist[0]->Material();
    dsassert(matlistptr->MaterialType() == INPAR::MAT::m_var_chemdiffusion,
        "material is not of type m_var_chemdiffusion");
    const Teuchos::RCP<const MAT::ScatraMatVarChemDiffusion>& actmat =
        Teuchos::rcp_dynamic_cast<const MAT::ScatraMatVarChemDiffusion>(matlistptr);

    // Computes dependent field
    double value_dependent_field = actmat->ComputeDualInternalEnergy(value_independent_field,
        actmat->RefMu(), actmat->RefC(), actmat->RefTemp() * actmat->GasConstant(), 1);


    int lid_in_global = phinp->Map().LID(gid_dependent_field);
    if (lid_in_global == -1) dserror("didn't find gid in global vector on this proc");
    phinpDirichlet->operator[](i) = value_dependent_field;
    phinp->operator[](gid_dependent_field) = value_dependent_field;
  }


  return;
}  // SCATRA::TimIntVariational::ApplyDirichletBC

/*----------------------------------------------------------------------*
 | add parameters depending on the problem                  deanda 10/17 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVariational::AddProblemSpecificParametersAndVectors(
    Teuchos::ParameterList& params  //!< parameter list
)
{
  return;
}

/*----------------------------------------------------------------------------------------*
 | initialize meshtying strategy (including standard case without meshtying) deanda 10/17 |
 *----------------------------------------------------------------------------------------*/
void SCATRA::TimIntVariational::CreateMeshtyingStrategy()
{
  // standard case without meshtying
  strategy_ = Teuchos::rcp(new MeshtyingStrategyStdVar(this));

  return;
}  // SCATRA::TimIntVariational::CreateMeshtyingStrategy()

/*----------------------------------------------------------------------*
 |  write current state to BINIO                          deanda   10/17|
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVariational::OutputState()
{
  ScaTraTimIntImpl::OutputState();

  if (!DRT::INPUT::IntegralValue<int>(params_->sublist("VARIATIONAL"), "ANALYTIC2PARAVIEW")) return;
  switch (calcerror_)
  {
    case INPAR::SCATRA::calcerror_byfunction:
    case INPAR::SCATRA::calcerror_AnalyticSeries:
    {
      if (time_ == 0) *phiAnalytic_ = *phinp_;

      if (nsd_ != 1 && calcerror_ == INPAR::SCATRA::calcerror_AnalyticSeries)
        dserror("An analytic series solution is not defined for this dimension ");

      const int numdofpernode = NumDofPerNode();
      Teuchos::RCP<DRT::Discretization> discret = Discretization();
      double L = 0;  // Saves length of the domain when needed

      // Extracts imposed Dirichlet condition
      Teuchos::RCP<Epetra_Vector> phinpDirichlet = dbcmaps_->ExtractCondVector(phinp_);
      double Cext =
          phinpDirichlet->operator[](0);  // Special condition since is Dirichlet only at one point

      if (calcerror_ == INPAR::SCATRA::calcerror_AnalyticSeries)
      {
        for (int i = 0; i < discret->NodeRowMap()->NumMyElements(); ++i)
        {
          // Get the coordinates of the mesh
          int gid = discret->NodeRowMap()->GID(i);
          DRT::Node* node = discret->gNode(gid);
          if (node->X()[0] > L) L = node->X()[0];
        }
      }

      // Get the coordinates of the mesh
      for (int i = 0; i < discret->NodeRowMap()->NumMyElements(); ++i)
      {
        // Get the coordinates of the mesh
        int gid = discret->NodeRowMap()->GID(i);
        DRT::Node* node = discret->gNode(gid);
        double my_x_coord = node->X()[0];
        double my_y_coord = node->X()[1];
        double my_z_coord = node->X()[2];

        // Extracts material parameters
        DRT::Element** elelist = node->Elements();
        const Teuchos::RCP<MAT::Material> matlistptr = elelist[0]->Material();
        dsassert(matlistptr->MaterialType() == INPAR::MAT::m_var_chemdiffusion,
            "material is not of type m_var_chemdiffusion");
        const Teuchos::RCP<const MAT::ScatraMatVarChemDiffusion>& actmat =
            Teuchos::rcp_dynamic_cast<const MAT::ScatraMatVarChemDiffusion>(matlistptr);


        double c_analytic(0);   // Analytic solution concentration
        double mu_analytic(0);  // Analytic solution chemical potential

        switch (calcerror_)
        {
          case INPAR::SCATRA::calcerror_AnalyticSeries:  // Computes concentration
          {
            const int series_end =
                params_->get<int>("CALCERRORNO");  // Number of iterations for analytic series

            if (time_ != 0)
            {
              AnalyticSolution_SeriesErrorFnt(c_analytic, series_end, my_x_coord, time_, L,
                  actmat->Diffusivity(0), actmat->RefC(), Cext);
              PostProcess_ChemPot(mu_analytic, c_analytic, node);
            }
            break;
          }  // case INPAR::SCATRA::calcerror_AnalyticSeries:
          case INPAR::SCATRA::calcerror_byfunction:
          {
            // function evaluation requires a 3D position vector!!
            double position[3] = {0.0, 0.0, 0.0};
            position[0] = my_x_coord;
            position[1] = my_y_coord;
            position[2] = my_z_coord;

            // gets function to evaluate
            const int errorfunctno = params_->get<int>("CALCERRORNO");
            if (errorfunctno < 1)
              dserror("invalid value of parameter CALCERRORNO for error function evaluation!");

            // Ensures we only output one species solution (since it is only code like this so far)
            int k = 0;  // TODO: Remove when updating for multple species
            if (numdofpernode != 2) dserror("The output is coded only for one species");

            // Computes concentration
            c_analytic =
                DRT::Problem::Instance()->Funct(errorfunctno - 1).Evaluate(k, position, time_);
            PostProcess_ChemPot(mu_analytic, c_analytic, node);
            break;
          }
          default:
            break;
        }  // switch(calcerror_)

        phiAnalytic_->operator[](i* numdofpernode) = c_analytic;
        phiAnalytic_->operator[](i* numdofpernode + 1) = mu_analytic;


      }  // For computation for all nodes

      output_->WriteVector("phiAnalytic", phiAnalytic_);
      break;
    }  // case
    default:
      break;
  }  // switch(calcerror_)

  return;
}  // ScaTraTimIntImpl::OutputState

void SCATRA::TimIntVariational::SaveError2File(
    const double RelError_Conc,    // !< Relative error for Concentration
    const double RelError_ChemPot  // !< Relative error for Chemical potential
)
{
  if (myrank_ == 0)
  {
    // save to file
    std::ostringstream temp;
    temp << NumDofPerNode() / 2;
    const std::string simulation = problem_->OutputControlFile()->FileName();
    const std::string fname = simulation + ".relerror";

    std::ofstream f;

    // create new error file and write initial error
    if (step_ == 0)
    {
      f.open(fname.c_str());
      f << "| Step | Time | rel. L2-error concentration | rel. L2-error chemical potential |"
        << std::endl;
    }

    // append error of the last time step to the error file
    else
      f.open(fname.c_str(), std::fstream::ate | std::fstream::app);

    f << step_ << " " << time_ << " " << std::setprecision(6) << RelError_Conc << " "
      << RelError_ChemPot << std::endl;

    f.flush();
    f.close();
  }

  return;
}  // ScaTraTimIntImpl::SaveError2File

/*----------------------------------------------------------------------*
 |  Analytic solution functions                           deanda   10/17|
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVariational::AnalyticSolution_SeriesErrorFnt(double& Conc,  //!< Concentration
    const int series_end,  //!< stopping number for the series
    const double x,        //!< evaluation position in x
    const double t,        //!< evaluation time
    const double L,        //!< length of the domain
    const double D,        //!< Diffusion parameter
    const double c_0,      //!< initial/reference concentration
    const double Cext      //!< concentration at Dirichlet nodes
)
{
  //   References:
  //   John Crank
  //   "Mathematics of Diffusion"
  //    Plane sheet problem
  //   1975, Ed II, 22-24 Eq(2.54)

  // Series parameters
  double sum1 = 0;
  double sum2 = 0;

  for (int m = 0; m < series_end; ++m)
  {
    sum1 = sum1 + pow(-1, m) * erfc(((2 * m + 1) * L - x) / (2 * sqrt(D * t)));
    sum2 = sum2 + pow(-1, m) * erfc(((2 * m + 1) * L + x) / (2 * sqrt(D * t)));
  }  // Series loop

  Conc = (Cext - c_0) * (sum1 + sum2) + c_0;

}  // ScaTraTimIntImpl::AnalyticSolution_SeriesErrorFnt

void SCATRA::TimIntVariational::PostProcess_ChemPot(double& ChemPot,  //!< chemical potential
    const double conct,                                               //!< concentration
    DRT::Node* node)  //!< Node: to extract corresponding post-process relation
{
  // Extracts material parameters
  DRT::Element** elelist = node->Elements();
  const Teuchos::RCP<MAT::Material> matlistptr = elelist[0]->Material();
  dsassert(matlistptr->MaterialType() == INPAR::MAT::m_var_chemdiffusion,
      "material is not of type m_var_chemdiffusion");
  const Teuchos::RCP<const MAT::ScatraMatVarChemDiffusion>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::ScatraMatVarChemDiffusion>(matlistptr);

  // Uses conjugate relation to compute the chemical potential form the concentration
  ChemPot = actmat->ComputeInternalEnergy(
      conct, actmat->RefMu(), actmat->RefC(), actmat->RefTemp() * actmat->GasConstant(), 1);

}  // TimIntVariational::PostProcess_ChemPot


/*************************************************************************
 *******          ScalarHandlerVar            ******
 ************************************************************************/

SCATRA::ScalarHandlerVar::ScalarHandlerVar() : numscal_() { return; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScalarHandlerVar::Setup(const ScaTraTimIntImpl* const scatratimint)
{
  // call base class
  ScalarHandler::Setup(scatratimint);

  // Safety check-up
  if (NumDofPerNode() % 2) dserror("Number of degrees of freedom must be pair!");

  numscal_.clear();
  numscal_.insert(NumDofPerNode() / 2);

  return;
}

/*-------------------------------------------------------------------------*
|  Determine number of Scalars per node in given condition     deanda 10/17 |
 *-------------------------------------------------------------------------*/
int SCATRA::ScalarHandlerVar::NumScalInCondition(
    const DRT::Condition& condition, const Teuchos::RCP<const DRT::Discretization>& discret) const
{
  CheckIsSetup();
  // for now only equal dof numbers are supported
  if (not equalnumdof_)
    dserror(
        "Different number of DOFs per node within ScaTra discretization! This is not supported for "
        "Elch!");

  return NumScal();
}
