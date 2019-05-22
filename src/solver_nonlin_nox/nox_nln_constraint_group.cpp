/*-----------------------------------------------------------*/
/*!
\file nox_nln_constraint_group.cpp

\brief %NOX::NLN implementation of a %NOX::Epetra::Group
       to handle constrained problems.

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_constraint_group.H"
#include "nox_nln_interface_required.H"
#include "nox_nln_aux.H"

#include <NOX_Epetra_LinearSystem.H>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::CONSTRAINT::Group::Group(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& grpOptionParams,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& i, const NOX::Epetra::Vector& x,
    const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys,
    const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr)
    : NOX::Epetra::Group(printParams, i, x, linSys),
      NOX::NLN::Group(printParams, grpOptionParams, i, x, linSys),
      userConstraintInterfaces_(iConstr)
{
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::CONSTRAINT::Group::Group(const NOX::NLN::CONSTRAINT::Group& source, CopyType type)
    : NOX::Epetra::Group(source, type),
      NOX::NLN::Group(source, type),
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
NOX::Abstract::Group& NOX::NLN::CONSTRAINT::Group::operator=(const NOX::Epetra::Group& source)
{
  NOX::NLN::Group::operator=(source);

  userConstraintInterfaces_ =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group&>(source).userConstraintInterfaces_;

  return *this;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const NOX::NLN::CONSTRAINT::ReqInterfaceMap& NOX::NLN::CONSTRAINT::Group::GetConstrInterfaces()
    const
{
  return userConstraintInterfaces_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const NOX::NLN::CONSTRAINT::Interface::Required>
NOX::NLN::CONSTRAINT::Group::GetConstraintInterfacePtr(const NOX::NLN::SolutionType& soltype) const
{
  return GetConstraintInterfacePtr(soltype, true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const NOX::NLN::CONSTRAINT::Interface::Required>
NOX::NLN::CONSTRAINT::Group::GetConstraintInterfacePtr(
    const NOX::NLN::SolutionType& soltype, const bool& errflag) const
{
  Teuchos::RCP<const NOX::NLN::CONSTRAINT::Interface::Required> constrptr = Teuchos::null;

  ReqInterfaceMap::const_iterator it = userConstraintInterfaces_.find(soltype);
  if (errflag and it == userConstraintInterfaces_.end())
  {
    std::ostringstream msg;
    msg << "The given NOX::NLN::SolutionType \"" << NOX::NLN::SolutionType2String(soltype)
        << "\" could not be found!";
    throwError("GetConstraintInterfacePtr", msg.str());
  }
  else if (it != userConstraintInterfaces_.end())
    constrptr = it->second;

  return constrptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::CONSTRAINT::Group::GetModelValue(
    const enum NOX::NLN::MeritFunction::MeritFctName mf_type) const
{
  double mrtFctVal = NOX::NLN::Group::GetModelValue(mf_type);

  // constraint contributions
  double constr_contr = 0.0;
  for (const auto& constr_iter : userConstraintInterfaces_)
    constr_contr += constr_iter.second->GetModelValue(mf_type);

  return mrtFctVal + constr_contr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::CONSTRAINT::Group::GetLinearizedModelTerms(const NOX::Abstract::Vector& dir,
    const enum NOX::NLN::MeritFunction::MeritFctName mf_type,
    const enum NOX::NLN::MeritFunction::LinOrder linorder,
    const enum NOX::NLN::MeritFunction::LinType lintype) const
{
  // contributions of the primary field
  double linVal = NOX::NLN::Group::GetLinearizedModelTerms(dir, mf_type, linorder, lintype);

  const NOX::Epetra::Vector& dir_nox_epetra = dynamic_cast<const NOX::Epetra::Vector&>(dir);
  const Epetra_Vector& dir_epetra = dir_nox_epetra.getEpetraVector();

  // constraint contributions
  for (const auto& constr_iter : userConstraintInterfaces_)
    linVal += constr_iter.second->GetLinearizedModelTerms(dir_epetra, mf_type, linorder, lintype);

  return linVal;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const std::vector<double>> NOX::NLN::CONSTRAINT::Group::GetRHSNorms(
    const std::vector<NOX::Abstract::Vector::NormType>& type,
    const std::vector<NOX::NLN::StatusTest::QuantityType>& chQ,
    Teuchos::RCP<const std::vector<NOX::StatusTest::NormF::ScaleType>> scale) const
{
  if (scale.is_null())
    scale = Teuchos::rcp(new std::vector<NOX::StatusTest::NormF::ScaleType>(
        chQ.size(), NOX::StatusTest::NormF::Unscaled));

  Teuchos::RCP<std::vector<double>> norms = Teuchos::rcp(new std::vector<double>(0));

  double rval = -1.0;
  for (std::size_t i = 0; i < chQ.size(); ++i)
  {
    rval = GetNlnReqInterfacePtr()->GetPrimaryRHSNorms(RHSVector.getEpetraVector(), chQ[i], type[i],
        (*scale)[i] == NOX::StatusTest::NormF::Scaled);
    if (rval >= 0.0)
    {
      norms->push_back(rval);
      continue;
    }

    // avoid the execution of this block for NaN and Inf return values
    if (rval < 0.0)
    {
      enum NOX::NLN::SolutionType soltype = NOX::NLN::AUX::ConvertQuantityType2SolutionType(chQ[i]);
      Teuchos::RCP<const NOX::NLN::CONSTRAINT::Interface::Required> constrptr =
          GetConstraintInterfacePtr(soltype);
      rval = constrptr->GetConstraintRHSNorms(RHSVector.getEpetraVector(), chQ[i], type[i],
          (*scale)[i] == NOX::StatusTest::NormF::Scaled);
    }

    if (rval >= 0.0)
      norms->push_back(rval);
    else if (rval < 0.0)
    {
      std::ostringstream msg;
      msg << "The desired quantity"
             " for the \"NormF\" Status Test could not be found! (enum="
          << chQ[i] << " | " << NOX::NLN::StatusTest::QuantityType2String(chQ[i])
          << " | return value=" << rval << ")" << std::endl;
      dserror(msg.str());
    }
    else
    {
      dserror("The norm value %e for quantity %s is not valid!", rval,
          NOX::NLN::StatusTest::QuantityType2String(chQ[i]).c_str());
    }
  }

  return norms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<double>> NOX::NLN::CONSTRAINT::Group::GetSolutionUpdateNorms(
    const NOX::Abstract::Vector& xOld, const std::vector<NOX::Abstract::Vector::NormType>& type,
    const std::vector<StatusTest::QuantityType>& chQ,
    Teuchos::RCP<const std::vector<StatusTest::NormUpdate::ScaleType>> scale) const
{
  const NOX::Epetra::Vector& xOldEpetra = dynamic_cast<const NOX::Epetra::Vector&>(xOld);
  if (scale.is_null())
    scale = Teuchos::rcp(new std::vector<StatusTest::NormUpdate::ScaleType>(
        chQ.size(), StatusTest::NormUpdate::Unscaled));

  Teuchos::RCP<std::vector<double>> norms = Teuchos::rcp(new std::vector<double>(0));

  double rval = -1.0;
  for (std::size_t i = 0; i < chQ.size(); ++i)
  {
    rval = GetNlnReqInterfacePtr()->GetPrimarySolutionUpdateNorms(xVector.getEpetraVector(),
        xOldEpetra.getEpetraVector(), chQ[i], type[i],
        (*scale)[i] == StatusTest::NormUpdate::Scaled);
    if (rval >= 0.0)
    {
      norms->push_back(rval);
      continue;
    }
    enum NOX::NLN::SolutionType soltype = NOX::NLN::AUX::ConvertQuantityType2SolutionType(chQ[i]);
    Teuchos::RCP<const NOX::NLN::CONSTRAINT::Interface::Required> constrptr =
        GetConstraintInterfacePtr(soltype);
    rval = constrptr->GetLagrangeMultiplierUpdateNorms(xVector.getEpetraVector(),
        xOldEpetra.getEpetraVector(), chQ[i], type[i],
        (*scale)[i] == StatusTest::NormUpdate::Scaled);

    if (rval >= 0.0)
      norms->push_back(rval);
    else
    {
      std::ostringstream msg;
      msg << "The desired quantity"
             " for the \"NormIncr\" Status Test could not be found! (enum="
          << chQ[i] << " | " << NOX::NLN::StatusTest::QuantityType2String(chQ[i]) << ")"
          << std::endl;
      throwError("GetSolutionUpdateNorms", msg.str());
    }
  }

  return norms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<double>> NOX::NLN::CONSTRAINT::Group::GetPreviousSolutionNorms(
    const NOX::Abstract::Vector& xOld, const std::vector<NOX::Abstract::Vector::NormType>& type,
    const std::vector<StatusTest::QuantityType>& chQ,
    Teuchos::RCP<const std::vector<StatusTest::NormUpdate::ScaleType>> scale) const
{
  const NOX::Epetra::Vector& xOldEpetra = dynamic_cast<const NOX::Epetra::Vector&>(xOld);
  if (scale.is_null())
    scale = Teuchos::rcp(new std::vector<StatusTest::NormUpdate::ScaleType>(
        chQ.size(), StatusTest::NormUpdate::Unscaled));

  Teuchos::RCP<std::vector<double>> norms = Teuchos::rcp(new std::vector<double>(0));

  double rval = -1.0;
  for (std::size_t i = 0; i < chQ.size(); ++i)
  {
    rval = GetNlnReqInterfacePtr()->GetPreviousPrimarySolutionNorms(xOldEpetra.getEpetraVector(),
        chQ[i], type[i], (*scale)[i] == StatusTest::NormUpdate::Scaled);
    if (rval >= 0.0)
    {
      norms->push_back(rval);
      continue;
    }
    enum NOX::NLN::SolutionType soltype = NOX::NLN::AUX::ConvertQuantityType2SolutionType(chQ[i]);
    Teuchos::RCP<const NOX::NLN::CONSTRAINT::Interface::Required> constrptr =
        GetConstraintInterfacePtr(soltype);
    rval = constrptr->GetPreviousLagrangeMultiplierNorms(xOldEpetra.getEpetraVector(), chQ[i],
        type[i], (*scale)[i] == StatusTest::NormUpdate::Scaled);

    if (rval >= 0.0)
      norms->push_back(rval);
    else
    {
      std::ostringstream msg;
      msg << "The desired quantity"
             " for the \"NormUpdate\" Status Test could not be found! (enum="
          << chQ[i] << " | " << NOX::NLN::StatusTest::QuantityType2String(chQ[i]) << ")"
          << std::endl;
      throwError("GetPreviousSolutionNorms", msg.str());
    }
  }

  return norms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<double>> NOX::NLN::CONSTRAINT::Group::GetSolutionUpdateRMS(
    const NOX::Abstract::Vector& xOld, const std::vector<double>& aTol,
    const std::vector<double>& rTol, const std::vector<NOX::NLN::StatusTest::QuantityType>& chQ,
    const std::vector<bool>& disable_implicit_weighting) const
{
  const NOX::Epetra::Vector& xOldEpetra = dynamic_cast<const NOX::Epetra::Vector&>(xOld);
  Teuchos::RCP<std::vector<double>> rms = Teuchos::rcp(new std::vector<double>(0));

  double rval = -1.0;
  for (std::size_t i = 0; i < chQ.size(); ++i)
  {
    rval = GetNlnReqInterfacePtr()->GetPrimarySolutionUpdateRMS(xVector.getEpetraVector(),
        xOldEpetra.getEpetraVector(), aTol[i], rTol[i], chQ[i], disable_implicit_weighting[i]);
    if (rval >= 0.0)
    {
      rms->push_back(rval);
      continue;
    }

    enum NOX::NLN::SolutionType soltype = NOX::NLN::AUX::ConvertQuantityType2SolutionType(chQ[i]);
    Teuchos::RCP<const NOX::NLN::CONSTRAINT::Interface::Required> constrptr =
        GetConstraintInterfacePtr(soltype);

    rval = constrptr->GetLagrangeMultiplierUpdateRMS(xVector.getEpetraVector(),
        xOldEpetra.getEpetraVector(), aTol[i], rTol[i], chQ[i], disable_implicit_weighting[i]);
    if (rval >= 0)
      rms->push_back(rval);
    else
    {
      std::ostringstream msg;
      msg << "The desired quantity"
             " for the \"NormWRMS\" Status Test could not be found! (enum="
          << chQ[i] << " | " << NOX::NLN::StatusTest::QuantityType2String(chQ[i]) << ")"
          << std::endl;
      throwError("GetSolutionUpdateRMS", msg.str());
    }
  }

  return rms;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> NOX::NLN::CONSTRAINT::Group::GetCurrentActiveSetMap(
    const enum NOX::NLN::StatusTest::QuantityType& qt) const
{
  enum NOX::NLN::SolutionType soltype = NOX::NLN::AUX::ConvertQuantityType2SolutionType(qt);

  return GetConstraintInterfacePtr(soltype)->GetCurrentActiveSetMap(qt);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> NOX::NLN::CONSTRAINT::Group::GetOldActiveSetMap(
    const enum NOX::NLN::StatusTest::QuantityType& qt) const
{
  enum NOX::NLN::SolutionType soltype = NOX::NLN::AUX::ConvertQuantityType2SolutionType(qt);

  return GetConstraintInterfacePtr(soltype)->GetOldActiveSetMap(qt);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
enum NOX::StatusTest::StatusType NOX::NLN::CONSTRAINT::Group::GetActiveSetInfo(
    const enum NOX::NLN::StatusTest::QuantityType& qt, int& activeset_size) const
{
  enum NOX::NLN::SolutionType soltype = NOX::NLN::AUX::ConvertQuantityType2SolutionType(qt);

  return GetConstraintInterfacePtr(soltype)->GetActiveSetInfo(qt, activeset_size);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::CONSTRAINT::Group::throwError(
    const std::string& functionName, const std::string& errorMsg) const
{
  if (utils.isPrintType(NOX::Utils::Error))
  {
    utils.out() << "ERROR - NOX::NLN::CONSTRAINT::Group::" << functionName << " - " << errorMsg
                << std::endl;
  }
  throw "NOX Error";
}
