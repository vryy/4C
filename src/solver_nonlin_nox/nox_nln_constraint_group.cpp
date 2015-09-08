/*-----------------------------------------------------------*/
/*!
\file nox_nln_constraint_group.cpp

\maintainer Michael Hiermeier

\date Jun 29, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_CONSTRAINT_GROUP_CPP_
#define NOX_NLN_CONSTRAINT_GROUP_CPP_

#include "nox_nln_constraint_group.H"
#include "nox_nln_interface_required.H"

#include <NOX_Epetra_LinearSystem.H>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::CONSTRAINT::Group::Group
(
  Teuchos::ParameterList& printParams,
  const Teuchos::RCP<NOX::Epetra::Interface::Required>& i,
  const NOX::Epetra::Vector& x,
  const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys,
  const Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required>& iConstr
) : NOX::NLN::Group(printParams,i,x,linSys),
  userConstraintInterfacePtr_(iConstr)
{
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::CONSTRAINT::Group::Group(const NOX::NLN::CONSTRAINT::Group& source, CopyType type) :
  NOX::NLN::Group(source,type),
  userConstraintInterfacePtr_(source.userConstraintInterfacePtr_)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Abstract::Group> NOX::NLN::CONSTRAINT::Group::clone(CopyType type) const
{
  Teuchos::RCP<NOX::Abstract::Group> newgrp =
    Teuchos::rcp(new NOX::NLN::CONSTRAINT::Group(*this, type));
  return newgrp;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::CONSTRAINT::Group::computeX(
    const NOX::Abstract::Group& grp,
    const NOX::Abstract::Vector& d,
    double step)
{
  // Cast to appropriate type, then call the "native" computeX
  const NOX::NLN::CONSTRAINT::Group* constrNlnGrp =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group*> (&grp);
  if (constrNlnGrp==NULL)
    throwError("computeX","dyn_cast to nox_nln_constraint_group failed!");

  const NOX::Epetra::Vector& epetrad =
    dynamic_cast<const NOX::Epetra::Vector&> (d);
  computeX(*constrNlnGrp, epetrad, step);
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::CONSTRAINT::Group::computeX(
    const NOX::NLN::CONSTRAINT::Group& grp,
    const NOX::Epetra::Vector& d,
    double step)
{
  // -------------------------------------
  // Update the external variables
  // -------------------------------------
  // We have to remove the const state, because we want to change class variables
  // of the constraint interface class.
  Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> iConstr =
      Teuchos::rcp_const_cast<NOX::NLN::CONSTRAINT::Interface::Required>(GetConstrPtr());
  iConstr->SetLagrangeMultiplier(grp.xVector,d,step);
  // -------------------------------------
  // Update the nox internal variables
  // and the primary variables
  // -------------------------------------
  NOX::NLN::Group::computeX(grp,d,step);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const std::vector<double> > NOX::NLN::CONSTRAINT::Group::GetRHSNorms(
    const std::vector<NOX::NLN::StatusTest::QuantityType>& chQ,
    const std::vector<NOX::Abstract::Vector::NormType>& type,
    Teuchos::RCP<const std::vector<NOX::StatusTest::NormF::ScaleType> > scale
    ) const
{
  if (scale.is_null())
    scale = Teuchos::rcp(new std::vector<NOX::StatusTest::NormF::ScaleType>(chQ.size(),
        NOX::StatusTest::NormF::Unscaled));

  Teuchos::RCP<std::vector<double> > norms = Teuchos::rcp(new std::vector<double>(0));

  double rval = -1.0;
  for (std::size_t i=0;i<chQ.size();++i)
  {
    rval = GetNlnReqInterfacePtr()->GetPrimaryRHSNorms(chQ[i],type[i],
        (*scale)[i]==NOX::StatusTest::NormF::Scaled);
    if (rval>=0.0)
    {
      norms->push_back(rval);
      continue;
    }

    rval = userConstraintInterfacePtr_->GetConstraintRHSNorms(chQ[i],type[i],
        (*scale)[i]==NOX::StatusTest::NormF::Scaled);
    if (rval>=0.0)
      norms->push_back(rval);
    else
    {
      std::ostringstream msg;
      msg << "The desired quantity"
          " for the \"NormF\" Status Test could not be found! (enum="
          << chQ[i] << " | " << NOX::NLN::StatusTest::QuantityType2String(chQ[i])
          << ")" << std::endl;
      throwError("GetRhsNorms",msg.str());
    }
  }

  return norms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<double> > NOX::NLN::CONSTRAINT::Group::GetSolutionUpdateRMS(
    const std::vector<double>& aTol, const std::vector<double>& rTol,
    const std::vector<NOX::NLN::StatusTest::QuantityType>& chQ,
    const std::vector<bool>& disable_implicit_weighting) const
{
  Teuchos::RCP<std::vector<double> > rms = Teuchos::rcp(new std::vector<double>(0));

  double rval = -1.0;
  for (std::size_t i=0;i<chQ.size();++i)
  {
    rval = GetNlnReqInterfacePtr()->GetPrimarySolutionUpdateRMS(
        aTol[i],rTol[i],chQ[i],disable_implicit_weighting[i]);
    if (rval>=0.0)
    {
      rms->push_back(rval);
      continue;
    }

    rval = userConstraintInterfacePtr_->GetLagrangeMultiplierUpdateRMS(
        aTol[i],rTol[i],chQ[i],disable_implicit_weighting[i]);
    if (rval>=0)
      rms->push_back(rval);
    else
    {
      std::ostringstream msg;
      msg << "The desired quantity"
          " for the \"NormWRMS\" Status Test could not be found! (enum="
          << chQ[i] << " | " << NOX::NLN::StatusTest::QuantityType2String(chQ[i])
          << ")" << std::endl;
      throwError("GetSolutionUpdateRMS",msg.str());
    }
  }

  return rms;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::CONSTRAINT::Group::throwError(
    const std::string& functionName,
    const std::string& errorMsg) const
{
  if (utils.isPrintType(NOX::Utils::Error)) {
    utils.out() << "ERROR - NOX::NLN::CONSTRAINT::Group::" << functionName
     << " - " << errorMsg << std::endl;
  }
  throw "NOX Error";
}

#endif /* NOX_NLN_CONSTRAINT_GROUP_CPP_ */
