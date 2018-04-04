/*-----------------------------------------------------------*/
/*!
\file str_timint_noxinterface.cpp

\brief Concrete implementation of the Jacobian, Required and
       Preconditioner %NOX::NLN interfaces.

\maintainer Michael Hiermeier

\date Nov 27, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_timint_noxinterface.H"
#include "str_timint_base.H"
#include "str_dbc.H"
#include "str_impl_generic.H"
#include "str_utils.H"

#include "../solver_nonlin_nox/nox_nln_aux.H"
#include "../solver_nonlin_nox/nox_nln_constraint_group.H"

#include "../linalg/linalg_sparseoperator.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_io/io_pstream.H"

#include <NOX_Epetra_Vector.H>

#include <boost/algorithm/string/predicate.hpp>  // case insensitive string compare

#include "../linalg/linalg_sparsematrix.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::NoxInterface::NoxInterface()
    : isinit_(false),
      issetup_(false),
      gstate_ptr_(Teuchos::null),
      implint_ptr_(Teuchos::null),
      dbc_ptr_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::NoxInterface::Init(
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
    const Teuchos::RCP<STR::IMPLICIT::Generic>& implint_ptr,
    const Teuchos::RCP<STR::Dbc>& dbc_ptr,
    const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr)
{
  // reset the setup flag
  issetup_ = false;

  gstate_ptr_  = gstate_ptr;
  timint_ptr_  = timint_ptr;
  implint_ptr_ = implint_ptr;
  dbc_ptr_     = dbc_ptr;

  // set the initialization flag
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::NoxInterface::Setup()
{
  CheckInit();

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::NoxInterface::CheckInit() const
{
  if (not IsInit())
    dserror("Call Init() first!");
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::NoxInterface::CheckInitSetup() const
{
  if (not IsInit() or not IsSetup())
    dserror("Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::IMPLICIT::Generic& STR::TIMINT::NoxInterface::ImplInt()
{
  CheckInitSetup();
  return *implint_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TIMINT::NoxInterface::computeF(const Epetra_Vector& x,
    Epetra_Vector& F, const FillType fillFlag)
{
  CheckInitSetup();

  if (not implint_ptr_->ApplyForce(x,F))
    return false;

  /* Apply the DBC on the right hand side, since we need the Dirichlet free
   * right hand side inside NOX for the convergence check, etc.               */
  Teuchos::RCP<Epetra_Vector> rhs_ptr = Teuchos::rcp(&F,false);
  dbc_ptr_->ApplyDirichletToRhs(rhs_ptr);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TIMINT::NoxInterface::computeJacobian(const Epetra_Vector& x,
    Epetra_Operator& jac)
{
  CheckInitSetup();

  LINALG::SparseOperator* jac_ptr =
      dynamic_cast<LINALG::SparseOperator*>(&jac);
  if (jac_ptr == NULL)
    dserror("Dynamic cast failed.");

  if (not implint_ptr_->ApplyStiff(x,*jac_ptr))
    return false;

  /* We do not consider the jacobian DBC at this point. The Dirichlet conditions
   * are applied inside the NOX::NLN::LinearSystem::applyJacobianInverse()
   * routine, instead. See the runPreApplyJacobianInverse() implementation
   * for more information.                               hiermeier 01/15/2016 */

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TIMINT::NoxInterface::computeFandJacobian(const Epetra_Vector& x,
    Epetra_Vector& rhs, Epetra_Operator& jac)
{
  CheckInitSetup();

  LINALG::SparseOperator* jac_ptr =
      dynamic_cast<LINALG::SparseOperator*>(&jac);
  if (jac_ptr==NULL)
    dserror("Dynamic cast failed!");

  if (not implint_ptr_->ApplyForceStiff(x,rhs,*jac_ptr))
    return false;

  /* Apply the DBC on the right hand side, since we need the Dirchilet free
   * right hand side inside NOX for the convergence check, etc.               */
  Teuchos::RCP<Epetra_Vector> rhs_ptr = Teuchos::rcp(&rhs,false);
  dbc_ptr_->ApplyDirichletToRhs(rhs_ptr);

  /* We do not consider the jacobian DBC at this point. The Dirichlet conditions
   * are applied inside the NOX::NLN::LinearSystem::applyJacobianInverse()
   * routine, instead. See the runPreApplyJacobianInverse() implementation
   * for more information.                               hiermeier 01/15/2016 */

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TIMINT::NoxInterface::computeCorrectionSystem(
    const enum NOX::NLN::CorrectionType type,
    const NOX::Abstract::Group& grp,
    const Epetra_Vector& x,
    Epetra_Vector& rhs,
    Epetra_Operator& jac )
{
  CheckInitSetup();

  LINALG::SparseOperator* jac_ptr =
      dynamic_cast<LINALG::SparseOperator*>(&jac);
  if (jac_ptr==NULL)
    dserror("Dynamic cast failed!");

  std::vector<INPAR::STR::ModelType> constraint_models;
  FindConstraintModels( &grp, constraint_models );

  if ( not implint_ptr_->ApplyCorrectionSystem(type,constraint_models,x,rhs,
      *jac_ptr) )
    return false;

  /* Apply the DBC on the right hand side, since we need the Dirchilet free
   * right hand side inside NOX for the convergence check, etc.               */
  Teuchos::RCP<Epetra_Vector> rhs_ptr = Teuchos::rcpFromRef(rhs);
  dbc_ptr_->ApplyDirichletToRhs(rhs_ptr);

   return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TIMINT::NoxInterface::computePreconditioner(
    const Epetra_Vector& x,
    Epetra_Operator& M, Teuchos::ParameterList* precParams)
{
  CheckInitSetup();
  // currently not supported
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TIMINT::NoxInterface::GetPrimaryRHSNorms(
    const Epetra_Vector& F,
    const NOX::NLN::StatusTest::QuantityType& checkquantity,
    const NOX::Abstract::Vector::NormType& type,
    const bool& isscaled) const
{
  CheckInitSetup();
  double rhsnorm = -1.0;

  // convert the given quantity type to a model type
  const INPAR::STR::ModelType mt =
      STR::NLN::ConvertQuantityType2ModelType(checkquantity);

  switch (checkquantity)
  {
    case NOX::NLN::StatusTest::quantity_structure:
    case NOX::NLN::StatusTest::quantity_cardiovascular0d:
    {
      // export the model specific solution if necessary
      Teuchos::RCP<Epetra_Vector> rhs_ptr =
          gstate_ptr_->ExtractModelEntries(mt,F);

      // remove entries specific to element technology
      gstate_ptr_->RemoveElementTechnologies(rhs_ptr);

      implint_ptr_->RemoveCondensedContributionsFromRhs(*rhs_ptr);

      rhsnorm = CalculateNorm(rhs_ptr,type,isscaled);

      break;
    }
    case NOX::NLN::StatusTest::quantity_pressure:
    {
      // export the model specific solution if necessary
      Teuchos::RCP<Epetra_Vector> rhs_ptr =
          gstate_ptr_->ExtractModelEntries(mt,F);

      // extract entries specific to element technology
      gstate_ptr_->ExtractElementTechnologies(NOX::NLN::StatusTest::quantity_pressure,rhs_ptr);

      rhsnorm = CalculateNorm(rhs_ptr,type,isscaled);

      break;
    }
    default:
    {
      /* Nothing to do. Functionality is supposed to be extended. */
      break;
    }
  }

  return rhsnorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TIMINT::NoxInterface::GetPrimarySolutionUpdateRMS(
    const Epetra_Vector& xnew, const Epetra_Vector& xold,
    const double& atol, const double& rtol,
    const NOX::NLN::StatusTest::QuantityType& checkquantity,
    const bool& disable_implicit_weighting) const
{
  CheckInitSetup();

  double rms = -1.0;

  // convert the given quantity type to a model type
  const INPAR::STR::ModelType mt =
      STR::NLN::ConvertQuantityType2ModelType(checkquantity);

  switch (checkquantity)
  {
    case NOX::NLN::StatusTest::quantity_structure:
    case NOX::NLN::StatusTest::quantity_cardiovascular0d:
    {
      // export the displacement solution if necessary
      Teuchos::RCP<Epetra_Vector> model_incr_ptr =
          Teuchos::rcp(new Epetra_Vector(*gstate_ptr_->ExtractModelEntries(mt,xold)));
      Teuchos::RCP<Epetra_Vector> model_xnew_ptr =
          gstate_ptr_->ExtractModelEntries(mt,xnew);

      // remove entries specific to element technology
      gstate_ptr_->RemoveElementTechnologies(model_incr_ptr);
      gstate_ptr_->RemoveElementTechnologies(model_xnew_ptr);

      model_incr_ptr->Update(1.0,*model_xnew_ptr,-1.0);
      rms = NOX::NLN::AUX::RootMeanSquareNorm(atol,rtol,model_xnew_ptr,model_incr_ptr,
          disable_implicit_weighting);

      break;
    }
    case NOX::NLN::StatusTest::quantity_pressure:
    {
      // export the displacement solution if necessary
      Teuchos::RCP<Epetra_Vector> model_incr_ptr =
          Teuchos::rcp(new Epetra_Vector(*gstate_ptr_->ExtractModelEntries(mt,xold)));
      Teuchos::RCP<Epetra_Vector> model_xnew_ptr =
          gstate_ptr_->ExtractModelEntries(mt,xnew);

      // extract entries specific to element technology
      gstate_ptr_->ExtractElementTechnologies(
          NOX::NLN::StatusTest::quantity_pressure,model_incr_ptr);
      gstate_ptr_->ExtractElementTechnologies(
          NOX::NLN::StatusTest::quantity_pressure,model_xnew_ptr);

      model_incr_ptr->Update(1.0,*model_xnew_ptr,-1.0);
      rms = NOX::NLN::AUX::RootMeanSquareNorm(atol,rtol,model_xnew_ptr,model_incr_ptr,
          disable_implicit_weighting);
      break;
    }
    case NOX::NLN::StatusTest::quantity_eas:
    case NOX::NLN::StatusTest::quantity_plasticity:
    {
      rms = implint_ptr_->GetCondensedSolutionUpdateRMS(checkquantity);
      break;
    }
    default:
    {
      /* Nothing to do. Functionality is supposed to be extended. */
      break;
    }
  }

  return rms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TIMINT::NoxInterface::GetPrimarySolutionUpdateNorms(
    const Epetra_Vector& xnew, const Epetra_Vector& xold,
    const NOX::NLN::StatusTest::QuantityType& checkquantity,
    const NOX::Abstract::Vector::NormType& type,
    const bool& isscaled) const
{
  CheckInitSetup();

  double updatenorm = -1.0;

  // convert the given quantity type to a model type
  const INPAR::STR::ModelType mt =
      STR::NLN::ConvertQuantityType2ModelType(checkquantity);

  switch (checkquantity)
  {
    case NOX::NLN::StatusTest::quantity_structure:
    case NOX::NLN::StatusTest::quantity_cardiovascular0d:
    {
      // export the displacement solution if necessary
      Teuchos::RCP<Epetra_Vector> model_incr_ptr =
          gstate_ptr_->ExtractModelEntries(mt,xold);
      Teuchos::RCP<Epetra_Vector> model_xnew_ptr =
          gstate_ptr_->ExtractModelEntries(mt,xnew);

      // remove entries specific to element technology
      gstate_ptr_->RemoveElementTechnologies(model_incr_ptr);
      gstate_ptr_->RemoveElementTechnologies(model_xnew_ptr);

      model_incr_ptr->Update(1.0,*model_xnew_ptr,-1.0);
      updatenorm = CalculateNorm(model_incr_ptr,type,isscaled);

      break;
    }
    case NOX::NLN::StatusTest::quantity_pressure:
    {
      // export the displacement solution if necessary
      Teuchos::RCP<Epetra_Vector> model_incr_ptr =
          gstate_ptr_->ExtractModelEntries(mt,xold);
      Teuchos::RCP<Epetra_Vector> model_xnew_ptr =
          gstate_ptr_->ExtractModelEntries(mt,xnew);

      // extract entries specific to element technology
      gstate_ptr_->ExtractElementTechnologies(
          NOX::NLN::StatusTest::quantity_pressure,model_incr_ptr);
      gstate_ptr_->ExtractElementTechnologies(
          NOX::NLN::StatusTest::quantity_pressure,model_xnew_ptr);

      model_incr_ptr->Update(1.0,*model_xnew_ptr,-1.0);
      updatenorm = CalculateNorm(model_incr_ptr,type,isscaled);

      break;
    }
    case NOX::NLN::StatusTest::quantity_eas:
    case NOX::NLN::StatusTest::quantity_plasticity:
    {
      // get the update norm of the condensed quantities
      updatenorm =
          implint_ptr_->GetCondensedUpdateNorm(checkquantity);
      // do the scaling if desired
      if (isscaled)
      {
        int gdofnumber = implint_ptr_->GetCondensedDofNumber(checkquantity);
        updatenorm /= static_cast<double>(gdofnumber);
      }
      break;
    }
    default:
    {
      /* Nothing to do. Functionality is supposed to be extended. */
      break;
    }
  }

  return updatenorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TIMINT::NoxInterface::GetPreviousPrimarySolutionNorms(
    const Epetra_Vector& xold,
    const NOX::NLN::StatusTest::QuantityType& checkquantity,
    const NOX::Abstract::Vector::NormType& type,
    const bool& isscaled) const
{
  CheckInitSetup();

  double xoldnorm = -1.0;

  // convert the given quantity type to a model type
  const INPAR::STR::ModelType mt =
      STR::NLN::ConvertQuantityType2ModelType(checkquantity);

  switch (checkquantity)
  {
    case NOX::NLN::StatusTest::quantity_structure:
    case NOX::NLN::StatusTest::quantity_cardiovascular0d:
    {
      // export the displacement solution if necessary
      Teuchos::RCP<Epetra_Vector> model_xold_ptr =
          gstate_ptr_->ExtractModelEntries(mt,xold);

      // remove entries specific to element technology
      gstate_ptr_->RemoveElementTechnologies(model_xold_ptr);

      xoldnorm = CalculateNorm(model_xold_ptr,type,isscaled);

      break;
    }
    case NOX::NLN::StatusTest::quantity_pressure:
    {
      // export the displacement solution if necessary
      Teuchos::RCP<Epetra_Vector> model_xold_ptr =
          gstate_ptr_->ExtractModelEntries(mt,xold);

      // extract entries specific to element technology
      gstate_ptr_->ExtractElementTechnologies(
          NOX::NLN::StatusTest::quantity_pressure,model_xold_ptr);

      xoldnorm = CalculateNorm(model_xold_ptr,type,isscaled);

      break;
    }
    case NOX::NLN::StatusTest::quantity_eas:
    case NOX::NLN::StatusTest::quantity_plasticity:
    {
      // get the update norm of the condensed quantities
      xoldnorm =
          implint_ptr_->GetCondensedPreviousSolNorm(checkquantity);
      if (isscaled)
      {
        int gdofnumber = implint_ptr_->GetCondensedDofNumber(checkquantity);
        xoldnorm /= static_cast<double>(gdofnumber);
      }
      break;
    }
    default:
    {
      /* Nothing to do. Functionality is supposed to be extended. */
      break;
    }
  }

  return xoldnorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TIMINT::NoxInterface::CalculateNorm(
    Teuchos::RCP<Epetra_Vector> quantity,
    const NOX::Abstract::Vector::NormType type,
    const bool isscaled ) const
{
  Teuchos::RCP<const NOX::Epetra::Vector> quantity_nox =
      Teuchos::rcp(new NOX::Epetra::Vector(quantity,
          NOX::Epetra::Vector::CreateView));

  double norm = quantity_nox->norm(type);
  // do the scaling if desired
  if (isscaled)
    norm /= static_cast<double>(quantity_nox->length());

  return norm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TIMINT::NoxInterface::GetModelValue(
    const Epetra_Vector& x,
    const Epetra_Vector& F,
    const enum NOX::NLN::MeritFunction::MeritFctName merit_func_type ) const
{
  CheckInitSetup();

  double omval = 0.0;

  switch ( merit_func_type )
  {
    case NOX::NLN::MeritFunction::mrtfct_lagrangian_active:
    case NOX::NLN::MeritFunction::mrtfct_lagrangian:
    case NOX::NLN::MeritFunction::mrtfct_energy:
    {
      // get the current displacement vector
      Teuchos::RCP<const Epetra_Vector> disnp =
          gstate_ptr_->ExtractDisplEntries(x);

      Teuchos::ParameterList p;
      // parameter needed by the elements
      p.set("action", "calc_struct_energy");

      Teuchos::RCP<DRT::DiscretizationInterface> discret_ptr =
          gstate_ptr_->GetMutableDiscret();
      // set vector values needed by elements
      discret_ptr->ClearState();
      discret_ptr->SetState("displacement", disnp);
      // get internal structural energy
      Teuchos::RCP<Epetra_SerialDenseVector> energy
        = Teuchos::rcp(new Epetra_SerialDenseVector(1));
      discret_ptr->EvaluateScalars(p, energy);
      discret_ptr->ClearState();

      omval = (*energy)(0);

      IO::cout << __LINE__ << " - " << __FUNCTION__ << "\n";
      IO::cout << "--- Structural energy\n";
      IO::cout << "energy = " << std::setprecision(15) << std::scientific << omval << "\n";

      break;
    }
    case NOX::NLN::MeritFunction::mrtfct_infeasibility_two_norm:
    case NOX::NLN::MeritFunction::mrtfct_infeasibility_two_norm_active:
    {
      // do nothing in the primary field
      break;
    }
    default:
    {
      dserror("There is no objective model value for %s | %d.",
          NOX::NLN::MeritFunction::MeritFuncName2String( merit_func_type ).c_str(),
          merit_func_type );
      exit( EXIT_FAILURE );
    }
  }

  return omval;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TIMINT::NoxInterface::GetLinearizedModelTerms(
    const NOX::Abstract::Group* group,
    const Epetra_Vector& dir,
    const enum NOX::NLN::MeritFunction::MeritFctName mf_type,
    const enum NOX::NLN::MeritFunction::LinOrder linorder,
    const enum NOX::NLN::MeritFunction::LinType lintype ) const
{
  switch ( mf_type )
  {
    case NOX::NLN::MeritFunction::mrtfct_lagrangian:
    case NOX::NLN::MeritFunction::mrtfct_lagrangian_active:
    {
      return GetLinearizedEnergyModelTerms( group, dir, linorder, lintype );
    }
    case NOX::NLN::MeritFunction::mrtfct_infeasibility_two_norm:
    case NOX::NLN::MeritFunction::mrtfct_infeasibility_two_norm_active:
      return 0.0;
    default:
    {
      dserror("There is no linearization for the objective model %s | %d.",
          NOX::NLN::MeritFunction::MeritFuncName2String(mf_type).c_str(),
          mf_type );
      exit( EXIT_FAILURE );
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TIMINT::NoxInterface::GetLinearizedEnergyModelTerms(
    const NOX::Abstract::Group* group,
    const Epetra_Vector& dir,
    const enum NOX::NLN::MeritFunction::LinOrder linorder,
    const enum NOX::NLN::MeritFunction::LinType lintype ) const
{
  double lin_val = 0.0;

  switch ( linorder )
  {
    case NOX::NLN::MeritFunction::linorder_first:
    case NOX::NLN::MeritFunction::linorder_all:
    {
      switch( lintype )
      {
        case NOX::NLN::MeritFunction::lin_wrt_all_dofs:
        case NOX::NLN::MeritFunction::lin_wrt_primary_dofs:
        {
          Epetra_Vector str_gradient( dir.Map(), true );

          std::vector<INPAR::STR::ModelType> constraint_models;
          FindConstraintModels( group, constraint_models );

          // assemble the force and exclude all constraint models
          implint_ptr_->AssembleForce( str_gradient, &constraint_models );
          str_gradient.Dot( dir, &lin_val );

          break;
        }
        default:
        {
          /* do nothing, there are only primary dofs */
          break;
        }
      }

      break;
    }
    default:
    {
      /* do nothing, there are no high order terms */
      break;
    }
  }

  return lin_val;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::NoxInterface::FindConstraintModels(
    const NOX::Abstract::Group* grp,
    std::vector<INPAR::STR::ModelType>& constraint_models ) const
{
  const NOX::NLN::CONSTRAINT::Group* constr_grp =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>( grp );

  // direct return if this is no constraint problem
  if ( not constr_grp )
    return;

  // find the constraint model types
  const auto& imap = constr_grp->GetConstrInterfaces();
  constraint_models.reserve( imap.size() );

  for ( auto cit = imap.begin(); cit != imap.end(); ++cit )
  {
    const enum NOX::NLN::SolutionType soltype = cit->first;
    const enum INPAR::STR::ModelType mtype = STR::NLN::ConvertSolType2ModelType( soltype );

    constraint_models.push_back( mtype );
  }
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TIMINT::NoxInterface::CalcRefNormForce()
{
  CheckInitSetup();
  const NOX::Epetra::Vector::NormType& nox_normtype =
      timint_ptr_->GetDataSDyn().GetNoxNormType();
  return implint_ptr_->CalcRefNormForce(nox_normtype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix>
STR::TIMINT::NoxInterface::CalcJacobianContributionsFromElementLevelForPTC()
{
  CheckInitSetup();
  Teuchos::RCP<LINALG::SparseMatrix> scalingMatrixOpPtr =
      Teuchos::rcp(new LINALG::SparseMatrix(*gstate_ptr_->DofRowMap(),81,true,true));
  implint_ptr_->ComputeJacobianContributionsFromElementLevelForPTC(scalingMatrixOpPtr);

  return scalingMatrixOpPtr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::NoxInterface::CreateBackupState(
    const Epetra_Vector& dir )
{
  CheckInitSetup();
  implint_ptr_->CreateBackupState(dir);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::NoxInterface::RecoverFromBackupState()
{
  CheckInitSetup();
  implint_ptr_->RecoverFromBackupState();
}
