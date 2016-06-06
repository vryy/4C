/*---------------------------------------------------------------------*/
/*!
\file str_model_evaluator_contact.cpp

\brief Evaluation and assembly of all contact terms

\maintainer Michael Hiermeier

\date Feb 3, 2016

\level 3

*/
/*---------------------------------------------------------------------*/


#include "str_model_evaluator_contact.H"
#include "str_model_evaluator_data.H"
#include "str_timint_base.H"
#include "str_integrator.H"
#include "str_dbc.H"
#include "str_utils.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparsematrix.H"

#include "../drt_contact/contact_poro_lagrange_strategy.H"
#include "../drt_contact/contact_strategy_factory.H"
#include "../drt_contact_aug/contact_augmented_strategy.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::Contact::Contact()
    : strategy_ptr_(Teuchos::null)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::Setup()
{
  CheckInit();
  eval_contact_ptr_ = EvalData().ContactPtr();
  // ---------------------------------------------------------------------
  // create the contact factory
  // ---------------------------------------------------------------------
  Teuchos::RCP<CONTACT::STRATEGY::Factory> factory_ptr = Teuchos::null;
  factory_ptr =  Teuchos::rcp(new CONTACT::STRATEGY::Factory());
  factory_ptr->Init(GStatePtr());
  factory_ptr->Setup();

  // check the problem dimension
  factory_ptr->CheckDimension();

  // create some local variables (later to be stored in strategy)
  std::vector<Teuchos::RCP<CONTACT::CoInterface> > interfaces;
  Teuchos::ParameterList cparams;

  // read and check contact input parameters
  factory_ptr->ReadAndCheckInput(cparams);

  // check for FillComplete of discretization
  if (!Discret().Filled())
    dserror("Discretization is not fillcomplete");
  // ---------------------------------------------------------------------
  // build the contact interfaces
  // ---------------------------------------------------------------------
  // FixMe Would be great, if we get rid of these poro parameters...
  bool poroslave = false;
  bool poromaster = false;
  factory_ptr->BuildInterfaces(cparams, interfaces,poroslave,poromaster);

  // ---------------------------------------------------------------------
  // build the solver strategy object
  // ---------------------------------------------------------------------
  strategy_ptr_ = factory_ptr->BuildStrategy(cparams,poroslave,
      poromaster,DofOffset(),interfaces);

  // build the search tree
  factory_ptr->BuildSearchTree(interfaces);

  // print final screen output
  factory_ptr->Print(interfaces,strategy_ptr_,cparams);

  // ---------------------------------------------------------------------
  // final touches to the contact strategy
  // ---------------------------------------------------------------------
  strategy_ptr_->StoreDirichletStatus(Int().GetDbc().GetDBCMapExtractor());
  strategy_ptr_->SetState(MORTAR::state_new_displacement,
      Int().GetDbc().GetZeros());
  strategy_ptr_->SaveReferenceState(Int().GetDbc().GetZerosPtr());
  strategy_ptr_->EvaluateReferenceState(Int().GetDbc().GetZerosPtr());
  strategy_ptr_->Inttime_init();
  strategy_ptr_->RedistributeContact(GState().GetDisN());
  strategy_ptr_->InitBinStrategyforTimestep(GState().GetVelN());

  CheckPseudo2D();

  // ---------------------------------------------------------------------
  // print current algorithm info
  // ---------------------------------------------------------------------
  PrintBanner();

  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::CheckPseudo2D() const
{
  // print messages for multifield problems (e.g FSI)
  const PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();
  if ((probtype != prb_structure) and (GState().GetMyRank()==0))
  {
    // warnings
  #ifdef CONTACTPSEUDO2D
    std::cout << "WARNING: The flag CONTACTPSEUDO2D is switched on. If this "
         << "is a real 3D problem, switch it off!" << std::endl;
  #else
    std::cout << "WARNING: The flag CONTACTPSEUDO2D is switched off. If this "
         << "is a 2D problem modeled pseudo-3D, switch it on!" << std::endl;
  #endif // #ifdef CONTACTPSEUDO2D
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::PrintBanner() const
{
  if (GState().GetMyRank()!=0)
    return;

  // some parameters
  const Teuchos::ParameterList&   smortar   = DRT::Problem::Instance()->MortarCouplingParams();
  const Teuchos::ParameterList&   scontact  = DRT::Problem::Instance()->ContactDynamicParams();
  INPAR::MORTAR::ShapeFcn         shapefcn  = DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(smortar,"LM_SHAPEFCN");
  INPAR::CONTACT::SolvingStrategy soltype   = DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(scontact,"STRATEGY");
  INPAR::CONTACT::SystemType      systype   = DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(scontact,"SYSTEM");
  INPAR::MORTAR::AlgorithmType    algorithm = DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(smortar,"ALGORITHM");
  bool smoothing = DRT::INPUT::IntegralValue<int>(scontact,"DISCR_SMOOTHING");

  if(smoothing)
  {
    if(soltype == INPAR::CONTACT::solution_lagmult)
    {
      std::cout << "================================================================\n";
      std::cout << "========= !!! EXPERIMENTAL VERSION  !!!     ====================\n";
      std::cout << "================================================================\n\n";
      std::cout << "================================================================\n";
      std::cout << "===== Interface smoothing approach with     ====================\n";
      std::cout << "===== Standard Lagrange multiplier strategy ====================\n";
      std::cout << "===== (Saddle point formulation) ===============================\n";
      std::cout << "================================================================\n\n";
    }
    else if (INPAR::CONTACT::solution_penalty)
    {
      std::cout << "================================================================\n";
      std::cout << "========= !!! EXPERIMENTAL VERSION  !!!     ====================\n";
      std::cout << "================================================================\n\n";
      std::cout << "================================================================\n";
      std::cout << "===== Interface smoothing approach with     ====================\n";
      std::cout << "===== Standard Penalty strategy             ====================\n";
      std::cout << "===== (Pure displacement formulation)===========================\n";
      std::cout << "================================================================\n\n";
    }
    else
      dserror("ERROR: Invalid system type for contact/meshtying interface smoothing");
  }
  else
  {
    if(algorithm == INPAR::MORTAR::algorithm_mortar)
    {
      // saddle point formulation
      if (systype == INPAR::CONTACT::system_saddlepoint)
      {
        if (soltype == INPAR::CONTACT::solution_lagmult && shapefcn == INPAR::MORTAR::shape_standard)
        {
          std::cout << "================================================================\n";
          std::cout << "===== Standard Lagrange multiplier strategy ====================\n";
          std::cout << "===== (Saddle point formulation) ===============================\n";
          std::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_lagmult && shapefcn == INPAR::MORTAR::shape_dual)
        {
          std::cout << "================================================================\n";
          std::cout << "===== Dual Lagrange multiplier strategy ========================\n";
          std::cout << "===== (Saddle point formulation) ===============================\n";
          std::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_lagmult && shapefcn == INPAR::MORTAR::shape_petrovgalerkin)
        {
          std::cout << "================================================================\n";
          std::cout << "===== Petrov-Galerkin Lagrange multiplier strategy =============\n";
          std::cout << "===== (Saddle point formulation) ===============================\n";
          std::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_penalty && shapefcn == INPAR::MORTAR::shape_standard)
        {
          std::cout << "================================================================\n";
          std::cout << "===== Standard Penalty strategy ================================\n";
          std::cout << "===== (Pure displacement formulation) ==========================\n";
          std::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_penalty && shapefcn == INPAR::MORTAR::shape_dual)
        {
          std::cout << "================================================================\n";
          std::cout << "===== Dual Penalty strategy ====================================\n";
          std::cout << "===== (Pure displacement formulation) ==========================\n";
          std::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_uzawa && shapefcn == INPAR::MORTAR::shape_standard)
        {
          std::cout << "================================================================\n";
          std::cout << "===== Uzawa Augmented Lagrange strategy ========================\n";
          std::cout << "===== (Pure displacement formulation) ==========================\n";
          std::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_uzawa && shapefcn == INPAR::MORTAR::shape_dual)
        {
          std::cout << "================================================================\n";
          std::cout << "===== Dual Uzawa Augmented Lagrange strategy ===================\n";
          std::cout << "===== (Pure displacement formulation) ==========================\n";
          std::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_augmented && shapefcn == INPAR::MORTAR::shape_standard)
        {
          std::cout << "================================================================\n";
          std::cout << "===== Augmented Lagrange strategy ==============================\n";
          std::cout << "===== (Saddle point formulation) ===============================\n";
          std::cout << "================================================================\n\n";
        }
        else dserror("ERROR: Invalid strategy or shape function type for contact/meshtying");
      }

      // condensed formulation
      else if (systype == INPAR::CONTACT::system_condensed || systype == INPAR::CONTACT::system_condensed_lagmult)
      {
        if (soltype == INPAR::CONTACT::solution_lagmult && shapefcn == INPAR::MORTAR::shape_dual)
        {
          std::cout << "================================================================\n";
          std::cout << "===== Dual Lagrange multiplier strategy ========================\n";
          std::cout << "===== (Condensed formulation) ==================================\n";
          std::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_lagmult && shapefcn == INPAR::MORTAR::shape_petrovgalerkin)
        {
          std::cout << "================================================================\n";
          std::cout << "===== Petrov-Galerkin Lagrange multiplier strategy =============\n";
          std::cout << "===== (Condensed formulation) ==================================\n";
          std::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_penalty && shapefcn == INPAR::MORTAR::shape_standard)
        {
          std::cout << "================================================================\n";
          std::cout << "===== Standard Penalty strategy ================================\n";
          std::cout << "===== (Pure displacement formulation) ==========================\n";
          std::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_penalty && shapefcn == INPAR::MORTAR::shape_dual)
        {
          std::cout << "================================================================\n";
          std::cout << "===== Dual Penalty strategy ====================================\n";
          std::cout << "===== (Pure displacement formulation) ==========================\n";
          std::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_uzawa && shapefcn == INPAR::MORTAR::shape_standard)
        {
          std::cout << "================================================================\n";
          std::cout << "===== Uzawa Augmented Lagrange strategy ========================\n";
          std::cout << "===== (Pure displacement formulation) ==========================\n";
          std::cout << "================================================================\n\n";
        }
        else if (soltype == INPAR::CONTACT::solution_uzawa && shapefcn == INPAR::MORTAR::shape_dual)
        {
          std::cout << "================================================================\n";
          std::cout << "===== Dual Uzawa Augmented Lagrange strategy ===================\n";
          std::cout << "===== (Pure displacement formulation) ==========================\n";
          std::cout << "================================================================\n\n";
        }
        else dserror("ERROR: Invalid strategy or shape function type for contact/meshtying");
      }
    }
    else if(algorithm == INPAR::MORTAR::algorithm_nts)
    {
      std::cout << "================================================================\n";
      std::cout << "===== Node-To-Segment approach =================================\n";
      std::cout << "================================================================\n\n";
    }
    else if(algorithm == INPAR::MORTAR::algorithm_gpts)
    {
      std::cout << "================================================================\n";
      std::cout << "===== Gauss-Point-To-Segment approach ==========================\n";
      std::cout << "================================================================\n\n";
    }
    // invalid system type
    else
      dserror("ERROR: Invalid system type for contact/meshtying");
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::Reset(const Epetra_Vector& x)
{
  CheckInitSetup();
  // get the displacement DoFs
  Strategy().Reset(*GState().GetDisNp());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::Reset(const Epetra_Vector& x,
    LINALG::SparseOperator& jac)
{
  CheckInitSetup();
  Reset(x);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Contact::ApplyForce(const Epetra_Vector& x,
    Epetra_Vector& f)
{
  CheckInitSetup();
  Reset(x);
  bool ok = true;
  // --- evaluate contact contributions ---------------------------------
  EvalContact().SetActionType(MORTAR::eval_force);
  Strategy().Evaluate(EvalData().Contact());
  // --- assemble right-hand-side ---------------------------------------
  AssembleRhs(f);

  return ok;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Contact::ApplyStiff(const Epetra_Vector& x,
    LINALG::SparseOperator& jac)
{
  CheckInitSetup();
  Reset(x,jac);
  bool ok = true;
  // --- evaluate contact contributions ---------------------------------
  EvalContact().SetActionType(MORTAR::eval_force_stiff);
  Strategy().Evaluate(EvalData().Contact());
  // --- assemble jacobian matrix ---------------------------------------
  AssembleJacobian(jac);

  return ok;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Contact::ApplyForceStiff(
    const Epetra_Vector& x,
    Epetra_Vector& f,
    LINALG::SparseOperator& jac)
{
  CheckInitSetup();
  Reset(x,jac);
  bool ok = true;
  // --- evaluate contact contributions ---------------------------------
  EvalContact().SetActionType(MORTAR::eval_force_stiff);
  Strategy().Evaluate(EvalData().Contact());
  // --- assemble right-hand-side ---------------------------------------
  AssembleRhs(f);
  // --- assemble jacobian matrix ---------------------------------------
  AssembleJacobian(jac);

  return ok;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::AssembleRhs(Epetra_Vector& f) const
{
  Teuchos::RCP<const Epetra_Vector> block_vec_ptr = Teuchos::null;
  if (Strategy().IsCondensedSystem())
  {
    block_vec_ptr = Strategy().GetCondensedRhsPtr();
    // if there are no active contact contributions, we can skip this...
    if (block_vec_ptr.is_null()) return;
    STR::AssembleVector(1.0,f,1.0,*block_vec_ptr);
  }
  else if (Strategy().IsSaddlePointSystem())
  {
    block_vec_ptr = Strategy().GetRhsBlockPtr(STR::block_displ);
    // if there are no active contact contributions, we can skip this...
    if (block_vec_ptr.is_null()) return;
    STR::AssembleVector(1.0,f,1.0,*block_vec_ptr);
    block_vec_ptr = Strategy().GetRhsBlockPtr(STR::block_constraint);
    if (block_vec_ptr.is_null())
      dserror("The constraint vector is a NULL pointer, although \n"
          "the structural part indicates, that contact contributions \n"
          "are present!");
    STR::AssembleVector(1.0,f,1.0,*block_vec_ptr);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::AssembleJacobian(
    LINALG::SparseOperator& jac) const
{
  Teuchos::RCP<const LINALG::SparseMatrix> block_ptr = Teuchos::null;
  // ---------------------------------------------------------------------
  // condensed system of equations
  // ---------------------------------------------------------------------
  if (Strategy().IsCondensedSystem())
  {
    block_ptr = Strategy().GetCondensedMatrixBlockPtr();
    // if there are no active contact contributions, we can skip this...
    if (block_ptr.is_null()) return;
    Teuchos::RCP<LINALG::SparseMatrix> jac_dd =
        GState().ExtractDisplBlock(jac);
    jac_dd->Add(*block_ptr,false,1.0,1.0);
  }
  // ---------------------------------------------------------------------
  // saddle-point system of equations or no contact contributions
  // ---------------------------------------------------------------------
  else if (Strategy().SystemType() == INPAR::CONTACT::system_saddlepoint)
  {
    // --- Kdd - block ---------------------------------------------------
    block_ptr = Strategy().GetMatrixBlockPtr(STR::block_displ_displ);
    /* if there are no active contact contributions, we put a identity
     * matrix at the (lm,lm)-block */
    if (block_ptr.is_null())
    {
      Teuchos::RCP<Epetra_Vector> ones =
          Teuchos::rcp(new Epetra_Vector(
              GState().BlockMap(INPAR::STR::model_contact),false));
      ones->PutScalar(1.0);
      block_ptr = Teuchos::rcp(new LINALG::SparseMatrix(*ones));
      GState().AssignModelBlock(jac,*block_ptr,INPAR::STR::model_contact,
          STR::block_lm_lm);
      // done ...
      return;
    }
    Teuchos::RCP<LINALG::SparseMatrix> jac_dd =
        GState().ExtractModelBlock(jac,INPAR::STR::model_contact,
            STR::block_displ_displ);
    jac_dd->Add(*block_ptr,false,1.0,1.0);
    // reset the block pointer, just to be on the safe side
    block_ptr = Teuchos::null;
    // --- Kdz - block ---------------------------------------------------
    block_ptr = Strategy().GetMatrixBlockPtr(STR::block_displ_lm);
    GState().AssignModelBlock(jac,*block_ptr,INPAR::STR::model_contact,
        STR::block_displ_lm);
    // reset the block pointer, just to be on the safe side
    block_ptr = Teuchos::null;
    // --- Kzd - block ---------------------------------------------------
    block_ptr = Strategy().GetMatrixBlockPtr(STR::block_lm_displ);
    GState().AssignModelBlock(jac,*block_ptr,INPAR::STR::model_contact,
        STR::block_lm_displ);
    // reset the block pointer, just to be on the safe side
    block_ptr = Teuchos::null;
    // --- Kzz - block ---------------------------------------------------
    block_ptr = Strategy().GetMatrixBlockPtr(STR::block_lm_lm);
    GState().AssignModelBlock(jac,*block_ptr,INPAR::STR::model_contact,
        STR::block_lm_lm);
    // reset the block pointer, just to be on the safe side
    block_ptr = Teuchos::null;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::WriteRestart(
    IO::DiscretizationWriter& iowriter,
    const bool& forced_writerestart) const
{
  // quantities to be written for restart
  std::map<std::string,Teuchos::RCP<Epetra_Vector> > restart_vectors;

  Strategy().DoWriteRestart(restart_vectors,forced_writerestart);

  // export restart information for contact to problem dof row map
  Teuchos::RCP<const Epetra_Map> problemdofs = Strategy().ProblemDofs();
  Teuchos::RCP<Epetra_Vector> lagrmultoldexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  LINALG::Export(*(Strategy().GetLagrMultN(true)), *lagrmultoldexp);
  iowriter.WriteVector("lagrmultold", lagrmultoldexp);

  // write all vectors specified by used strategy
  for (std::map<std::string,Teuchos::RCP<Epetra_Vector> >::const_iterator p=restart_vectors.begin();
      p!=restart_vectors.end();++p)
  {
    Teuchos::RCP<Epetra_Vector> expVec = Teuchos::rcp(new Epetra_Vector(*problemdofs));
    LINALG::Export(*p->second,*expVec);
    iowriter.WriteVector(p->first,expVec);
  }

  // since the global OutputStepState() routine is not called, if the
  // restart is written, we have to do it here manually.
  OutputStepState(iowriter);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::ReadRestart(
    IO::DiscretizationReader& ioreader)
{
  // reader strategy specific stuff
  Strategy().DoReadRestart(ioreader,GState().GetDisN(),
      EvalData().ContactPtr());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::UpdateStepState()
{
  /* Note: DisN() and DisNp() have the same value at this stage, since
   * we call the structural model evaluator always in first place! */
  strategy_ptr_->Update(GState().GetDisN());
  PostUpdateStepState();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::PostUpdateStepState()
{
  strategy_ptr_->Inttime_init();
  strategy_ptr_->RedistributeContact(GState().GetDisN());
  strategy_ptr_->InitBinStrategyforTimestep(GState().GetVelN());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::UpdateStepElement()
{
  /* empty */
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::RecoverState(
    const Epetra_Vector& xold,
    const Epetra_Vector& dir,
    const Epetra_Vector& xnew)
{
  CheckInitSetup();
  Reset(xnew);
  std::vector<Teuchos::RCP<const Epetra_Vector> > eval_vec(3,Teuchos::null);
  eval_vec[0] = Teuchos::rcpFromRef(xold);
  eval_vec[1] = Teuchos::rcpFromRef(dir);
  eval_vec[2] = Teuchos::rcpFromRef(xnew);
  EvalContact().SetActionType(MORTAR::eval_recover);
  Strategy().Evaluate(EvalData().Contact(),eval_vec);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::DetermineStressStrain()
{
  // evaluate contact tractions
  Strategy().OutputStresses();

  if (Strategy().WeightedWear())
  {
    /* *******************************************************************
     * We do not compute the non-weighted wear here. we just write
     * the output. the non-weighted wear will be used as dirichlet-b.
     * for the ale problem. n.w.wear will be called in
     * stru_ale_algorithm.cpp and computed in strategy.OutputWear();
     *                                                         farah 06/13
     * ******************************************************************/
    // evaluate wear (not weighted)
    Strategy().ContactWear()->PutScalar(0.0);
    Strategy().OutputWear();
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::DetermineEnergy()
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Contact::OutputStepState(
    IO::DiscretizationWriter& iowriter) const
{
  // *********************************************************************
  // active contact set and slip set
  // *********************************************************************

  // evaluate active set and slip set
  Teuchos::RCP<Epetra_Vector> activeset =
      Teuchos::rcp( new Epetra_Vector(*Strategy().ActiveRowNodes()));
  activeset->PutScalar(1.0);
  if (Strategy().Friction())
  {
    Teuchos::RCP<Epetra_Vector> slipset =
        Teuchos::rcp(new Epetra_Vector(*Strategy().SlipRowNodes()));
    slipset->PutScalar(1.0);
    Teuchos::RCP<Epetra_Vector> slipsetexp =
        Teuchos::rcp(new Epetra_Vector(*Strategy().ActiveRowNodes()));
    LINALG::Export(*slipset, *slipsetexp);
    activeset->Update(1.0, *slipsetexp, 1.0);
  }

  // export to problem node row map
  Teuchos::RCP<const Epetra_Map> problemnodes = Strategy().ProblemNodes();
  Teuchos::RCP<Epetra_Vector> activesetexp = Teuchos::rcp( new Epetra_Vector(*problemnodes));
  LINALG::Export(*activeset, *activesetexp);

  if (Strategy().WearBothDiscrete())
  {
    Teuchos::RCP<Epetra_Vector> mactiveset = Teuchos::rcp(new Epetra_Vector(*Strategy().MasterActiveNodes()));
    mactiveset->PutScalar(1.0);
    Teuchos::RCP<Epetra_Vector> slipset = Teuchos::rcp(new Epetra_Vector(*Strategy().MasterSlipNodes()));
    slipset->PutScalar(1.0);
    Teuchos::RCP<Epetra_Vector> slipsetexp = Teuchos::rcp( new Epetra_Vector(*Strategy().MasterActiveNodes()));
    LINALG::Export(*slipset, *slipsetexp);
    mactiveset->Update(1.0, *slipsetexp, 1.0);

    Teuchos::RCP<Epetra_Vector> mactivesetexp = Teuchos::rcp( new Epetra_Vector(*problemnodes));
    LINALG::Export(*mactiveset, *mactivesetexp);
    activesetexp->Update(1.0, *mactivesetexp, 1.0);
  }

  iowriter.WriteVector("activeset", activesetexp);

  // *********************************************************************
  // contact tractions
  // *********************************************************************

  // export to problem dof row map
  Teuchos::RCP<const Epetra_Map> problemdofs = Strategy().ProblemDofs();

  // normal direction
  Teuchos::RCP<const Epetra_Vector> normalstresses = Strategy().ContactNorStress();
  Teuchos::RCP<Epetra_Vector> normalstressesexp = Teuchos::rcp( new Epetra_Vector(*problemdofs));
  LINALG::Export(*normalstresses, *normalstressesexp);

  // tangential plane
  Teuchos::RCP<const Epetra_Vector> tangentialstresses = Strategy().ContactTanStress();
  Teuchos::RCP<Epetra_Vector> tangentialstressesexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  LINALG::Export(*tangentialstresses, *tangentialstressesexp);

  // write to output
  // contact tractions in normal and tangential direction
  iowriter.WriteVector("norcontactstress", normalstressesexp);
  iowriter.WriteVector("tancontactstress", tangentialstressesexp);

  #ifdef CONTACTFORCEOUTPUT
  dserror("Untested in the new structural framework!");
  // *********************************************************************
  // contact forces on slave non master side,
  // in normal and tangential direction
  // *********************************************************************
  // vectors for contact forces
  Teuchos::RCP<Epetra_Vector> fcslavenor = Teuchos::rcp(new Epetra_Vector(Strategy().DMatrix()->RowMap()));
  Teuchos::RCP<Epetra_Vector> fcslavetan = Teuchos::rcp(new Epetra_Vector(Strategy().DMatrix()->RowMap()));
  Teuchos::RCP<Epetra_Vector> fcmasternor = Teuchos::rcp(new Epetra_Vector(Strategy().MMatrix()->DomainMap()));
  Teuchos::RCP<Epetra_Vector> fcmastertan = Teuchos::rcp(new Epetra_Vector(Strategy().MMatrix()->DomainMap()));

  // vectors with problem dof row map
  Teuchos::RCP<Epetra_Vector> fcslavenorexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  Teuchos::RCP<Epetra_Vector> fcslavetanexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  Teuchos::RCP<Epetra_Vector> fcmasternorexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));
  Teuchos::RCP<Epetra_Vector> fcmastertanexp = Teuchos::rcp(new Epetra_Vector(*problemdofs));

  // multiplication
  Strategy().DMatrix()->Multiply(true, *normalstresses, *fcslavenor);
  Strategy().DMatrix()->Multiply(true, *tangentialstresses, *fcslavetan);
  Strategy().MMatrix()->Multiply(true, *normalstresses, *fcmasternor);
  Strategy().MMatrix()->Multiply(true, *tangentialstresses, *fcmastertan);

  #ifdef MASTERNODESINCONTACT
  //BEGIN: to output the global ID's of the master nodes in contact - devaal 02.2011

  int dim = DRT::Problem::Instance()->NDim();

  if (dim == 2)
  dserror("Only working for 3D");

  std::vector<int> lnid, gnid;

  //std::cout << "MasterNor" << fcmasternor->MyLength() << std::endl;

  for (int i=0; i<fcmasternor->MyLength(); i=i+3)
  {

    //check if master node in contact
    if (sqrt(((*fcmasternor)[i])*((*fcmasternor)[i])+((*fcmasternor)[i+1])*((*fcmasternor)[i+1])+((*fcmasternor)[i+2])*((*fcmasternor)[i]+2)) > 0.00001)
    {
      lnid.push_back((fcmasternor->Map()).GID(i)/3);
    }
  }

  // we want to gather data from on all procs
  std::vector<int> allproc(Comm().NumProc());
  for (int i=0; i<Comm().NumProc(); ++i) allproc[i] = i;

  // communicate all data to proc 0
  LINALG::Gather<int>(lnid,gnid,(int)allproc.size(),&allproc[0],Comm());

  //std::cout << " size of gnid:" << gnid.size() << std::endl;

  ////////////////
  ///// attempt at obtaining the nid and relative displacement u of master nodes in contact - devaal
  // define my own interface
  MORTAR::StrategyBase& myStrategy = strategy;
  CoAbstractStrategy& myContactStrategy = dynamic_cast<CoAbstractStrategy&>(myStrategy);

  std::vector<Teuchos::RCP<CONTACT::CoInterface> > myInterface = Strategy().ContactInterfaces();

  //check interface size - just doing this now for a single interface

  if (myInterface.size() != 1)
  dserror("Interface size should be 1");

  std::cout << "OUTPUT OF MASTER NODE IN CONTACT" << std::endl;
  //std::cout << "Master_node_in_contact x_dis y_dis z_dis" << std::endl;
  for (int i=0; i<(int)gnid.size(); ++i)
  {
    int myGid = gnid[i];
    std::cout << gnid[i] << std::endl; // << " " << myUx << " " << myUy << " " << myUz << std::endl;
  }

  #endif  //MASTERNODESINCONTACT: to output the global ID's of the master nodes in contact
  // export
  LINALG::Export(*fcslavenor,*fcslavenorexp);
  LINALG::Export(*fcslavetan,*fcslavetanexp);
  LINALG::Export(*fcmasternor,*fcmasternorexp);
  LINALG::Export(*fcmastertan,*fcmastertanexp);

  // contact forces on slave and master side
  iowriter.WriteVector("norslaveforce",fcslavenorexp);
  iowriter.WriteVector("tanslaveforce",fcslavetanexp);
  iowriter.WriteVector("normasterforce",fcmasternorexp);
  iowriter.WriteVector("tanmasterforce",fcmastertanexp);

  #ifdef CONTACTEXPORT
  // export averaged node forces to xxx.force
  double resultnor[fcslavenor->NumVectors()];
  double resulttan[fcslavetan->NumVectors()];
  fcslavenor->Norm2(resultnor);
  fcslavetan->Norm2(resulttan);

  if(Comm().MyPID()==0)
  {
    std::cout << "resultnor= " << resultnor[0] << std::endl;
    std::cout << "resulttan= " << resulttan[0] << std::endl;

    FILE* MyFile = NULL;
    std::ostringstream filename;
    const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix();
    filename << filebase << ".force";
    MyFile = fopen(filename.str().c_str(), "at+");
    if (MyFile)
    {
      //fprintf(MyFile,valuename.c_str());
      fprintf(MyFile, "%g\t",resultnor[0]);
      fprintf(MyFile, "%g\n",resulttan[0]);
      fclose(MyFile);
    }
    else
    dserror("ERROR: File for Output could not be opened.");
  }
  #endif //CONTACTEXPORT
  #endif //CONTACTFORCEOUTPUT

  // Evaluate the interface forces for the augmented Lagrange formulation
  const CONTACT::AugmentedLagrangeStrategy* aug_strat_ptr =
      dynamic_cast<const CONTACT::AugmentedLagrangeStrategy*>(&Strategy());
  if (aug_strat_ptr != NULL)
  {
    Teuchos::RCP<Epetra_Vector> augfs_lm = Teuchos::rcp(new Epetra_Vector(*problemdofs));
    Teuchos::RCP<Epetra_Vector> augfs_g  = Teuchos::rcp(new Epetra_Vector(*problemdofs));
    Teuchos::RCP<Epetra_Vector> augfm_lm = Teuchos::rcp(new Epetra_Vector(*problemdofs));
    Teuchos::RCP<Epetra_Vector> augfm_g  = Teuchos::rcp(new Epetra_Vector(*problemdofs));

    // evaluate augmented contact forces
    aug_strat_ptr->AugForces(*augfs_lm,*augfs_g,*augfm_lm,*augfm_g);

    // contact forces on slave and master side
    iowriter.WriteVector("norslaveforcelm",augfs_lm);
    iowriter.WriteVector("norslaveforceg" ,augfs_g);
    iowriter.WriteVector("normasterforcelm",augfm_lm);
    iowriter.WriteVector("normasterforceg" ,augfm_g);
  }

  // *********************************************************************
  // wear with internal state variable approach
  // *********************************************************************
  if (Strategy().WeightedWear())
  {
    // write output
    Teuchos::RCP<const Epetra_Vector> wearoutput    = Strategy().ContactWear();
    Teuchos::RCP<Epetra_Vector> wearoutputexp = Teuchos::rcp( new Epetra_Vector(*problemdofs));
    LINALG::Export(*wearoutput, *wearoutputexp);
    iowriter.WriteVector("wear", wearoutputexp);
  }

  // *********************************************************************
  // poro contact
  // *********************************************************************
  if (Strategy().HasPoroNoPenetration())
  {
    //output of poro no penetration lagrange multiplier!
    const CONTACT::PoroLagrangeStrategy& poro_strategy =
        dynamic_cast<const CONTACT::PoroLagrangeStrategy&>(Strategy());
    Teuchos::RCP<const Epetra_Vector> lambdaout     = poro_strategy.LambdaNoPen();
    Teuchos::RCP<Epetra_Vector> lambdaoutexp  = Teuchos::rcp(new Epetra_Vector(*problemdofs));
    LINALG::Export(*lambdaout, *lambdaoutexp);
    iowriter.WriteVector("poronopen_lambda",lambdaoutexp);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::ContactData& STR::MODELEVALUATOR::Contact::EvalContact()
{
  CheckInitSetup();
  return *eval_contact_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Teuchos::RCP<CONTACT::CoAbstractStrategy>&
    STR::MODELEVALUATOR::Contact::StrategyPtr()
{
  CheckInitSetup();
  return strategy_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::CoAbstractStrategy& STR::MODELEVALUATOR::Contact::Strategy()
{
  CheckInitSetup();
  return *strategy_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const CONTACT::CoAbstractStrategy& STR::MODELEVALUATOR::Contact::Strategy()
    const
{
  CheckInitSetup();
  return *strategy_ptr_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::Contact::
    GetBlockDofRowMapPtr() const
{
  return Strategy().LMDoFRowMapPtr(false);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Contact::
    GetCurrentSolutionPtr() const
{
  Teuchos::RCP<Epetra_Vector> curr_lm_ptr =
      Teuchos::rcp(new Epetra_Vector(*Strategy().GetLagrMultNp(false)));
  if (not curr_lm_ptr.is_null())
    curr_lm_ptr->ReplaceMap(*GetBlockDofRowMapPtr());
  return curr_lm_ptr;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Contact::
    GetLastTimeStepSolutionPtr() const
{
  Teuchos::RCP<Epetra_Vector> old_lm_ptr =
      Teuchos::rcp(new Epetra_Vector(*Strategy().GetLagrMultN(false)));
  if (not old_lm_ptr.is_null())
    old_lm_ptr->ReplaceMap(*GetBlockDofRowMapPtr());
  return old_lm_ptr;
}
