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
#include "nox_nln_aux.H"

#include <NOX_Epetra_LinearSystem.H>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::CONSTRAINT::Group::Group(Teuchos::ParameterList& printParams,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& i,
    const NOX::Epetra::Vector& x,
    const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys,
    const std::map<NOX::NLN::SolutionType,Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> >& iConstr)
    : NOX::NLN::Group(printParams,i,x,linSys),
      userConstraintInterfaces_(iConstr)
{
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::CONSTRAINT::Group::Group(const NOX::NLN::CONSTRAINT::Group& source,
    CopyType type)
    : NOX::NLN::Group(source,type),
      userConstraintInterfaces_(source.userConstraintInterfaces_)
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
const std::map<NOX::NLN::SolutionType,Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> >&
    NOX::NLN::CONSTRAINT::Group::GetConstrInterfaces() const
{
  return userConstraintInterfaces_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const NOX::NLN::CONSTRAINT::Interface::Required>
    NOX::NLN::CONSTRAINT::Group::GetConstraintInterfacePtr(
    const NOX::NLN::SolutionType& soltype) const
{
  return GetConstraintInterfacePtr(soltype,true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const NOX::NLN::CONSTRAINT::Interface::Required>
    NOX::NLN::CONSTRAINT::Group::GetConstraintInterfacePtr(
    const NOX::NLN::SolutionType& soltype,
    const bool& errflag) const
{
  Teuchos::RCP<const NOX::NLN::CONSTRAINT::Interface::Required> constrptr =
      Teuchos::null;
  std::map<NOX::NLN::SolutionType,Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::
           Required> >::const_iterator it;
  it = userConstraintInterfaces_.find(soltype);
  if (errflag and it == userConstraintInterfaces_.end())
  {
    std::ostringstream msg;
    msg << "The given NOX::NLN::SolutionType \""
        << NOX::NLN::SolutionType2String(soltype)
        << "\" could not be found!";
    throwError("GetConstraintInterfacePtr",msg.str());
  }
  else if (it != userConstraintInterfaces_.end())
    constrptr = it->second;

  return constrptr;
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
  std::map<NOX::NLN::SolutionType,Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> >::iterator it;
  for (it=userConstraintInterfaces_.begin();it!=userConstraintInterfaces_.end();++it)
  {
    it->second->SetLagrangeMultiplier(grp.xVector,d,step);
  }
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
    const std::vector<NOX::Abstract::Vector::NormType>& type,
    const std::vector<NOX::NLN::StatusTest::QuantityType>& chQ,
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

    enum NOX::NLN::SolutionType soltype =
        NOX::NLN::AUX::ConvertQuantityType2SolutionType(chQ[i]);
    Teuchos::RCP<const NOX::NLN::CONSTRAINT::Interface::Required> constrptr =
        GetConstraintInterfacePtr(soltype);
    rval = constrptr->GetConstraintRHSNorms(chQ[i],type[i],
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

    enum NOX::NLN::SolutionType soltype =
        NOX::NLN::AUX::ConvertQuantityType2SolutionType(chQ[i]);
    Teuchos::RCP<const NOX::NLN::CONSTRAINT::Interface::Required> constrptr =
        GetConstraintInterfacePtr(soltype);

    rval = constrptr->GetLagrangeMultiplierUpdateRMS(
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
enum NOX::StatusTest::StatusType
    NOX::NLN::CONSTRAINT::Group::GetActiveSetStatus(
    const enum NOX::NLN::StatusTest::QuantityType& qtype,
    int& activesetsize,
    int& cyclesize) const
{
  enum NOX::NLN::SolutionType soltype =
      NOX::NLN::AUX::ConvertQuantityType2SolutionType(qtype);

  return GetConstraintInterfacePtr(soltype)->GetActiveSetStatus(
      activesetsize,cyclesize);
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
