/*----------------------------------------------------------------------*/
/*!

\brief  Time integration for variational formulation of chemical diffusion,
    it inherits directly from scatra_timint_ost.

\level 2

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/
#include "scatra_timint_var_chemdiffusion.H"
// To read initial conditions from input file
#include "../drt_lib/drt_globalproblem.H"
// To extract the material
#include "../drt_lib/drt_node.H"
#include "../drt_mat/scatra_mat_var_chemdiffusion.H"
#include "../drt_nurbs_discret/drt_apply_nurbs_initial_condition.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                     deanda 09/17|
 *----------------------------------------------------------------------*/
SCATRA::TimIntVarChemDiffusionOST::TimIntVarChemDiffusionOST(
    Teuchos::RCP<DRT::Discretization> actdis,          //!< discretization
    Teuchos::RCP<LINALG::Solver> solver,               //!< linear solver
    Teuchos::RCP<Teuchos::ParameterList> params,       //!< parameter list
    Teuchos::RCP<Teuchos::ParameterList> extraparams,  //!< supplementary parameter list
    Teuchos::RCP<IO::DiscretizationWriter> output,     //!< output writer
    const int probnum                                  //!< global problem number
    )
    : ScaTraTimIntImpl(actdis, solver, params, extraparams, output, probnum),
      TimIntOneStepTheta(actdis, solver, params, extraparams, output, probnum),
      phi0_(Teuchos::null),
      semImplicitFunctional_(
          DRT::INPUT::IntegralValue<int>(params->sublist("VARIATIONAL"), "SEMIMPLICITFUNCTIONAL"))
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                              deanda 09/17|
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVarChemDiffusionOST::Init()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntOneStepTheta::Init();

  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                             deanda 09/17|
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVarChemDiffusionOST::Setup()
{
  // call Setup()-functions of base classes
  // note: this order is important
  // back-up of the initial condition

  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  phi0_ = LINALG::CreateVector(*dofrowmap, true);
  TimIntOneStepTheta::Setup();
  return;
}

void SCATRA::TimIntVarChemDiffusionOST::AddTimeIntegrationSpecificVectors(
    bool forcedincrementalsolver)
{
  // call base class routine
  ScaTraTimIntImpl::AddTimeIntegrationSpecificVectors(forcedincrementalsolver);
  TimIntOneStepTheta::AddTimeIntegrationSpecificVectors(forcedincrementalsolver);

  discret_->SetState("phi0", phi0_);

  return;
}

/*----------------------------------------------------------------------*
| Destructor dtor (public)                                  deanda 09/17|
*-----------------------------------------------------------------------*/
SCATRA::TimIntVarChemDiffusionOST::~TimIntVarChemDiffusionOST() { return; }


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                          deanda 09/17|
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVarChemDiffusionOST::Update(const int num)
{
  TimIntOneStepTheta::Update(num);

  return;
}


/*----------------------------------------------------------------------*
 | explicit predictor for nonlinear solver                   deanda 09/17|
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVarChemDiffusionOST::ExplicitPredictor() const
{
  // call base class routine
  TimIntOneStepTheta::ExplicitPredictor();

  return;
}

/*----------------------------------------------------------------------*
 |  set initial field for phi                               deanda 09/17 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVarChemDiffusionOST::SetInitialField(
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

/*----------------------------------------------------------------------*
 |  set dirichlet bc for dependent fields                   deanda 09/17 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVarChemDiffusionOST::ApplyDirichletBC(
    const double time, Teuchos::RCP<Epetra_Vector> phinp, Teuchos::RCP<Epetra_Vector> phidt)
{
  ScaTraTimIntImpl::ApplyDirichletBC(time, phinp, phidt);


  /*!
    Algebraic dependence between dof with Dirichlet condition
    brief: This part is needed since in the variational formulation only the boundary conditions are
    prescribed for the chemical potential. Based on that when the concentration is also defined at
    the nodes it must be used the Internal energy to relate the concentration at the dirichlet
    boundary to that of the chemical potential
   */
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

    double value_dependent_field = actmat->ComputeDualInternalEnergy(value_independent_field,
        actmat->RefMu(), actmat->RefC(), actmat->RefTemp() * actmat->GasConstant(), 1);


    int lid_in_global = phinp->Map().LID(gid_dependent_field);
    if (lid_in_global == -1) dserror("didn't find gid in global vector on this proc");
    phinpDirichlet->operator[](i) = value_dependent_field;
    phinp->operator[](gid_dependent_field) = value_dependent_field;
  }

  return;
}  // SCATRA::ScaTraTimIntImpl::ApplyDirichletBC

/*----------------------------------------------------------------------*
 | add parameters depending on the problem                   deanda 09/17|
 *----------------------------------------------------------------------*/
void SCATRA::TimIntVarChemDiffusionOST::AddProblemSpecificParametersAndVectors(
    Teuchos::ParameterList& params  //!< parameter list
)
{
  // Defines the type of scheme to use in the variational functional: Implicit or semi-implicit
  params.set<bool>("Is_semImplicit_Functional", semImplicitFunctional_);
  return;
}
