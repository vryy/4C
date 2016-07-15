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

#include "../solver_nonlin_nox/nox_nln_aux.H"

#include "../linalg/linalg_sparseoperator.H"
#include "../drt_lib/drt_discret.H"

#include <NOX_Epetra_Vector.H>

#include <boost/algorithm/string/predicate.hpp>  // case insensitive string compare

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

  /* Apply the DBC on the right hand side, since we need the Dirchilet free
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
bool STR::TIMINT::NoxInterface::computePreconditioner(
    const Epetra_Vector& x,
    Epetra_Operator& M, Teuchos::ParameterList* precParams)
{
  CheckInitSetup();
  // currently not supported
  // ToDo add the scaled thickness conditioning (STC) approach here
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

  switch (checkquantity)
  {
    case NOX::NLN::StatusTest::quantity_structure:
    {
      // export the displacement solution if necessary
      Teuchos::RCP<Epetra_Vector> rhs_ptr =
          gstate_ptr_->ExportDisplEntries(F);

      // transform to a NOX::Epetra::Vector
      Teuchos::RCP<const NOX::Epetra::Vector> rhs_nox_ptr =
          Teuchos::rcp(new NOX::Epetra::Vector(rhs_ptr,
              NOX::Epetra::Vector::CreateView));
      rhsnorm = rhs_nox_ptr->norm(type);
      // do the scaling if desired
      if (isscaled)
        rhsnorm /= static_cast<double>(rhs_nox_ptr->length());
      break;
    }
    case NOX::NLN::StatusTest::quantity_cardiovascular0d:
    {
      // export the solution if necessary
      Teuchos::RCP<Epetra_Vector> rhs0d_ptr =
          gstate_ptr_->ExportModelEntries(INPAR::STR::model_cardiovascular0d,F);

      // transform to a NOX::Epetra::Vector
      Teuchos::RCP<const NOX::Epetra::Vector> rhs0d_nox_ptr =
          Teuchos::rcp(new NOX::Epetra::Vector(rhs0d_ptr,
              NOX::Epetra::Vector::CreateView));
      rhsnorm = rhs0d_nox_ptr->norm(type);
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

  switch (checkquantity)
  {
    case NOX::NLN::StatusTest::quantity_structure:
    {
      // export the displacement solution if necessary
      Teuchos::RCP<Epetra_Vector> disincr_ptr =
          Teuchos::rcp(new Epetra_Vector(*gstate_ptr_->ExportDisplEntries(xold)));
      Teuchos::RCP<const Epetra_Vector> disnew_ptr =
          gstate_ptr_->ExportDisplEntries(xnew);

      disincr_ptr->Update(1.0,*disnew_ptr,-1.0);
      rms = NOX::NLN::AUX::RootMeanSquareNorm(atol,rtol,disnew_ptr,disincr_ptr,
          disable_implicit_weighting);
      break;
    }
    case NOX::NLN::StatusTest::quantity_eas:
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

  switch (checkquantity)
  {
    case NOX::NLN::StatusTest::quantity_structure:
    {
      // export the displacement solution if necessary
      Teuchos::RCP<Epetra_Vector> disincr_ptr =
          gstate_ptr_->ExportDisplEntries(xold);
      Teuchos::RCP<const Epetra_Vector> disnew_ptr =
          gstate_ptr_->ExportDisplEntries(xnew);

      disincr_ptr->Update(1.0,*disnew_ptr,-1.0);
      Teuchos::RCP<const NOX::Epetra::Vector> disincr_nox_ptr =
          Teuchos::rcp(new NOX::Epetra::Vector(disincr_ptr,NOX::Epetra::Vector::CreateView));

      updatenorm = disincr_nox_ptr->norm(type);
      // do the scaling if desired
      if (isscaled)
        updatenorm /= static_cast<double>(disincr_nox_ptr->length());

      break;
    }
    case NOX::NLN::StatusTest::quantity_cardiovascular0d:
    {
      // export the solution if necessary
      Teuchos::RCP<Epetra_Vector> cv0dincr_ptr =
          gstate_ptr_->ExportModelEntries(INPAR::STR::model_cardiovascular0d,xold);
      Teuchos::RCP<const Epetra_Vector> cv0dnew_ptr =
          gstate_ptr_->ExportModelEntries(INPAR::STR::model_cardiovascular0d,xnew);

      cv0dincr_ptr->Update(1.0,*cv0dnew_ptr,-1.0);
      Teuchos::RCP<const NOX::Epetra::Vector> cv0dincr_nox_ptr =
          Teuchos::rcp(new NOX::Epetra::Vector(cv0dincr_ptr,NOX::Epetra::Vector::CreateView));

      updatenorm = cv0dincr_nox_ptr->norm(type);

      break;
    }
    case NOX::NLN::StatusTest::quantity_eas:
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

  switch (checkquantity)
  {
    case NOX::NLN::StatusTest::quantity_structure:
    {
      // export the displacement solution if necessary
      Teuchos::RCP<Epetra_Vector> disold_ptr =
          gstate_ptr_->ExportDisplEntries(xold);

      Teuchos::RCP<const NOX::Epetra::Vector> disold_nox_ptr =
          Teuchos::rcp(new NOX::Epetra::Vector(disold_ptr,NOX::Epetra::Vector::CreateView));

      xoldnorm = disold_nox_ptr->norm(type);
      // do the scaling if desired
      if (isscaled)
        xoldnorm /= static_cast<double>(disold_nox_ptr->length());

      break;
    }
    case NOX::NLN::StatusTest::quantity_cardiovascular0d:
    {
      // export the solution if necessary
      Teuchos::RCP<Epetra_Vector> cv0dold_ptr =
          gstate_ptr_->ExportModelEntries(INPAR::STR::model_cardiovascular0d,xold);

      Teuchos::RCP<const NOX::Epetra::Vector> cv0dold_nox_ptr =
          Teuchos::rcp(new NOX::Epetra::Vector(cv0dold_ptr,NOX::Epetra::Vector::CreateView));

      xoldnorm = cv0dold_nox_ptr->norm(type);

      break;
    }
    case NOX::NLN::StatusTest::quantity_eas:
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
double STR::TIMINT::NoxInterface::GetObjectiveModelValue(
    const Epetra_Vector& x,
    const Epetra_Vector& F,
    const std::string& name) const
{
  CheckInitSetup();

  double omval = 0.0;

  if (boost::iequals(name,"energy"))
  {
    // get the current displacement vector
    Teuchos::RCP<const Epetra_Vector> disnp =
        gstate_ptr_->ExportDisplEntries(x);

    Teuchos::ParameterList p;
    // parameter needed by the elements
    p.set("action", "calc_struct_energy");

    Teuchos::RCP<DRT::Discretization> discret_ptr =
        Teuchos::rcp_const_cast<DRT::Discretization>(gstate_ptr_->GetDiscret());
    // set vector values needed by elements
    discret_ptr->ClearState();
    discret_ptr->SetState("displacement", disnp);
    // get internal structural energy
    Teuchos::RCP<Epetra_SerialDenseVector> energy
      = Teuchos::rcp(new Epetra_SerialDenseVector(1));
    discret_ptr->EvaluateScalars(p, energy);
    discret_ptr->ClearState();

    omval = (*energy)(0);
  }
  else
    dserror("There is no objective model value with the name \"%s\".",
        name.c_str());

  return omval;
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
