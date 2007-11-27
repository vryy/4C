
#ifdef CCADISCRET

#include <Thyra_DefaultProductVector.hpp>
#include <Thyra_DefaultBlockedLinearOp.hpp>

#include "mfsi_nox_thyra_group.H"
#include "mfsi_algorithm.H"

MFSI::NOX_Thyra_Group::NOX_Thyra_Group(const NOX::Thyra::Vector& initial_guess,
                                       const Algorithm* mfsi)
  : NOX::Thyra::Group(initial_guess, Teuchos::rcp(mfsi,false)),
    mfsi_(mfsi)
{
  // we know we already have the first linear system calculated

  Thyra::DefaultProductVector<double>& f = dynamic_cast<Thyra::DefaultProductVector<double>&>(*f_vec_->getThyraRCPVector());
  mfsi->SetupRHS(f);

  Teuchos::RCP<Thyra::DefaultBlockedLinearOp<double> > mat = mfsi_->SysMat();
  mfsi->SetupSysMat(*mat);

  is_valid_f_ = true;
  is_valid_jacobian_ = true;
}


NOX::Abstract::Group::ReturnType MFSI::NOX_Thyra_Group::computeF()
{
  NOX::Abstract::Group::ReturnType ret = NOX::Thyra::Group::computeF();

  if (ret==NOX::Abstract::Group::Ok and not this->isJacobian())
  {
    Teuchos::RCP<Thyra::DefaultBlockedLinearOp<double> > mat = mfsi_->SysMat();
    mfsi_->SetupSysMat(*mat);
    is_valid_jacobian_ = true;
  }

  return ret;
}


NOX::Abstract::Group::ReturnType MFSI::NOX_Thyra_Group::computeJacobian()
{
  NOX::Abstract::Group::ReturnType ret = NOX::Thyra::Group::computeJacobian();

  if (ret==NOX::Abstract::Group::Ok and not this->isF())
  {
    Thyra::DefaultProductVector<double>& f = dynamic_cast<Thyra::DefaultProductVector<double>&>(*f_vec_->getThyraRCPVector());
    mfsi_->SetupRHS(f);
    is_valid_f_ = true;
  }

  return ret;
}

#endif
