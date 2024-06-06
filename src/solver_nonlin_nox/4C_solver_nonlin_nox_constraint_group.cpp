/*-----------------------------------------------------------*/
/*! \file

\brief %NOX::NLN implementation of a %::NOX::Epetra::Group
       to handle constrained problems.



\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_solver_nonlin_nox_constraint_group.hpp"

#include "4C_solver_nonlin_nox_aux.hpp"
#include "4C_solver_nonlin_nox_interface_required.hpp"

#include <NOX_Epetra_LinearSystem.H>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::CONSTRAINT::Group::Group(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& grpOptionParams,
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& i, const ::NOX::Epetra::Vector& x,
    const Teuchos::RCP<::NOX::Epetra::LinearSystem>& linSys,
    const NOX::Nln::CONSTRAINT::ReqInterfaceMap& iConstr)
    : ::NOX::Epetra::Group(printParams, i, x, linSys),
      NOX::Nln::Group(printParams, grpOptionParams, i, x, linSys),
      user_constraint_interfaces_(iConstr)
{
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::CONSTRAINT::Group::Group(const NOX::Nln::CONSTRAINT::Group& source, ::NOX::CopyType type)
    : ::NOX::Epetra::Group(source, type),
      NOX::Nln::Group(source, type),
      user_constraint_interfaces_(source.user_constraint_interfaces_)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Abstract::Group> NOX::Nln::CONSTRAINT::Group::clone(::NOX::CopyType type) const
{
  Teuchos::RCP<::NOX::Abstract::Group> newgrp =
      Teuchos::rcp(new NOX::Nln::CONSTRAINT::Group(*this, type));
  return newgrp;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group& NOX::Nln::CONSTRAINT::Group::operator=(const ::NOX::Epetra::Group& source)
{
  NOX::Nln::Group::operator=(source);

  user_constraint_interfaces_ =
      dynamic_cast<const NOX::Nln::CONSTRAINT::Group&>(source).user_constraint_interfaces_;

  return *this;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const NOX::Nln::CONSTRAINT::ReqInterfaceMap& NOX::Nln::CONSTRAINT::Group::GetConstrInterfaces()
    const
{
  return user_constraint_interfaces_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const NOX::Nln::CONSTRAINT::Interface::Required>
NOX::Nln::CONSTRAINT::Group::get_constraint_interface_ptr(
    const NOX::Nln::SolutionType soltype) const
{
  return get_constraint_interface_ptr(soltype, true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const NOX::Nln::CONSTRAINT::Interface::Required>
NOX::Nln::CONSTRAINT::Group::get_constraint_interface_ptr(
    const NOX::Nln::SolutionType soltype, const bool errflag) const
{
  Teuchos::RCP<const NOX::Nln::CONSTRAINT::Interface::Required> constrptr = Teuchos::null;

  ReqInterfaceMap::const_iterator it = user_constraint_interfaces_.find(soltype);
  if (errflag and it == user_constraint_interfaces_.end())
  {
    std::ostringstream msg;
    msg << "The given NOX::Nln::SolutionType \"" << NOX::Nln::SolutionType2String(soltype)
        << "\" could not be found!";
    throw_error("get_constraint_interface_ptr", msg.str());
  }
  else if (it != user_constraint_interfaces_.end())
    constrptr = it->second;

  return constrptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::Nln::CONSTRAINT::Group::GetModelValue(
    const enum NOX::Nln::MeritFunction::MeritFctName merit_func_type) const
{
  double mrtFctVal = NOX::Nln::Group::GetModelValue(merit_func_type);

  // constraint contributions
  double constr_contr = 0.0;
  for (const auto& constr_iter : user_constraint_interfaces_)
    constr_contr += constr_iter.second->GetModelValue(merit_func_type);

  return mrtFctVal + constr_contr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::Nln::CONSTRAINT::Group::get_linearized_model_terms(const ::NOX::Abstract::Vector& dir,
    const enum NOX::Nln::MeritFunction::MeritFctName merit_func_type,
    const enum NOX::Nln::MeritFunction::LinOrder linorder,
    const enum NOX::Nln::MeritFunction::LinType lintype) const
{
  // contributions of the primary field
  double linVal =
      NOX::Nln::Group::get_linearized_model_terms(dir, merit_func_type, linorder, lintype);

  const ::NOX::Epetra::Vector& dir_nox_epetra = dynamic_cast<const ::NOX::Epetra::Vector&>(dir);
  const Epetra_Vector& dir_epetra = dir_nox_epetra.getEpetraVector();

  // constraint contributions
  for (const auto& constr_iter : user_constraint_interfaces_)
    linVal += constr_iter.second->get_linearized_model_terms(
        dir_epetra, merit_func_type, linorder, lintype);

  return linVal;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const std::vector<double>> NOX::Nln::CONSTRAINT::Group::GetRHSNorms(
    const std::vector<::NOX::Abstract::Vector::NormType>& type,
    const std::vector<NOX::Nln::StatusTest::QuantityType>& chQ,
    Teuchos::RCP<const std::vector<::NOX::StatusTest::NormF::ScaleType>> scale) const
{
  if (scale.is_null())
    scale = Teuchos::rcp(new std::vector<::NOX::StatusTest::NormF::ScaleType>(
        chQ.size(), ::NOX::StatusTest::NormF::Unscaled));

  Teuchos::RCP<std::vector<double>> norms = Teuchos::rcp(new std::vector<double>(0));

  double rval = -1.0;
  for (std::size_t i = 0; i < chQ.size(); ++i)
  {
    rval = get_nln_req_interface_ptr()->GetPrimaryRHSNorms(RHSVector.getEpetraVector(), chQ[i],
        type[i], (*scale)[i] == ::NOX::StatusTest::NormF::Scaled);
    if (rval >= 0.0)
    {
      norms->push_back(rval);
      continue;
    }

    // avoid the execution of this block for NaN and Inf return values
    if (rval < 0.0)
    {
      enum NOX::Nln::SolutionType soltype = NOX::Nln::Aux::ConvertQuantityType2SolutionType(chQ[i]);
      Teuchos::RCP<const NOX::Nln::CONSTRAINT::Interface::Required> constrptr =
          get_constraint_interface_ptr(soltype, false);
      if (constrptr != Teuchos::null)
        rval = constrptr->get_constraint_rhs_norms(RHSVector.getEpetraVector(), chQ[i], type[i],
            (*scale)[i] == ::NOX::StatusTest::NormF::Scaled);
      else
        rval = 0.0;
    }

    if (rval >= 0.0)
      norms->push_back(rval);
    else if (rval < 0.0)
    {
      std::ostringstream msg;
      msg << "The desired quantity"
             " for the \"NormF\" Status Test could not be found! (enum="
          << chQ[i] << " | " << NOX::Nln::StatusTest::QuantityType2String(chQ[i])
          << " | return value=" << rval << ")" << std::endl;
      FOUR_C_THROW(msg.str());
    }
    else
    {
      FOUR_C_THROW("The norm value %e for quantity %s is not valid!", rval,
          NOX::Nln::StatusTest::QuantityType2String(chQ[i]).c_str());
    }
  }

  return norms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<double>> NOX::Nln::CONSTRAINT::Group::get_solution_update_norms(
    const ::NOX::Abstract::Vector& xOld, const std::vector<::NOX::Abstract::Vector::NormType>& type,
    const std::vector<StatusTest::QuantityType>& chQ,
    Teuchos::RCP<const std::vector<StatusTest::NormUpdate::ScaleType>> scale) const
{
  const ::NOX::Epetra::Vector& xOldEpetra = dynamic_cast<const ::NOX::Epetra::Vector&>(xOld);
  if (scale.is_null())
    scale = Teuchos::rcp(new std::vector<StatusTest::NormUpdate::ScaleType>(
        chQ.size(), StatusTest::NormUpdate::Unscaled));

  Teuchos::RCP<std::vector<double>> norms = Teuchos::rcp(new std::vector<double>(0));

  double rval = -1.0;
  for (std::size_t i = 0; i < chQ.size(); ++i)
  {
    rval = get_nln_req_interface_ptr()->get_primary_solution_update_norms(xVector.getEpetraVector(),
        xOldEpetra.getEpetraVector(), chQ[i], type[i],
        (*scale)[i] == StatusTest::NormUpdate::Scaled);
    if (rval >= 0.0)
    {
      norms->push_back(rval);
      continue;
    }
    enum NOX::Nln::SolutionType soltype = NOX::Nln::Aux::ConvertQuantityType2SolutionType(chQ[i]);
    Teuchos::RCP<const NOX::Nln::CONSTRAINT::Interface::Required> constrptr =
        get_constraint_interface_ptr(soltype);
    rval = constrptr->get_lagrange_multiplier_update_norms(xVector.getEpetraVector(),
        xOldEpetra.getEpetraVector(), chQ[i], type[i],
        (*scale)[i] == StatusTest::NormUpdate::Scaled);

    if (rval >= 0.0)
      norms->push_back(rval);
    else
    {
      std::ostringstream msg;
      msg << "The desired quantity"
             " for the \"NormIncr\" Status Test could not be found! (enum="
          << chQ[i] << " | " << NOX::Nln::StatusTest::QuantityType2String(chQ[i]) << ")"
          << std::endl;
      throw_error("get_solution_update_norms", msg.str());
    }
  }

  return norms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<double>> NOX::Nln::CONSTRAINT::Group::get_previous_solution_norms(
    const ::NOX::Abstract::Vector& xOld, const std::vector<::NOX::Abstract::Vector::NormType>& type,
    const std::vector<StatusTest::QuantityType>& chQ,
    Teuchos::RCP<const std::vector<StatusTest::NormUpdate::ScaleType>> scale) const
{
  const ::NOX::Epetra::Vector& xOldEpetra = dynamic_cast<const ::NOX::Epetra::Vector&>(xOld);
  if (scale.is_null())
    scale = Teuchos::rcp(new std::vector<StatusTest::NormUpdate::ScaleType>(
        chQ.size(), StatusTest::NormUpdate::Unscaled));

  Teuchos::RCP<std::vector<double>> norms = Teuchos::rcp(new std::vector<double>(0));

  double rval = -1.0;
  for (std::size_t i = 0; i < chQ.size(); ++i)
  {
    rval = get_nln_req_interface_ptr()->get_previous_primary_solution_norms(
        xOldEpetra.getEpetraVector(), chQ[i], type[i],
        (*scale)[i] == StatusTest::NormUpdate::Scaled);
    if (rval >= 0.0)
    {
      norms->push_back(rval);
      continue;
    }
    enum NOX::Nln::SolutionType soltype = NOX::Nln::Aux::ConvertQuantityType2SolutionType(chQ[i]);
    Teuchos::RCP<const NOX::Nln::CONSTRAINT::Interface::Required> constrptr =
        get_constraint_interface_ptr(soltype);
    rval = constrptr->get_previous_lagrange_multiplier_norms(xOldEpetra.getEpetraVector(), chQ[i],
        type[i], (*scale)[i] == StatusTest::NormUpdate::Scaled);

    if (rval >= 0.0)
      norms->push_back(rval);
    else
    {
      std::ostringstream msg;
      msg << "The desired quantity"
             " for the \"NormUpdate\" Status Test could not be found! (enum="
          << chQ[i] << " | " << NOX::Nln::StatusTest::QuantityType2String(chQ[i]) << ")"
          << std::endl;
      throw_error("get_previous_solution_norms", msg.str());
    }
  }

  return norms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::vector<double>> NOX::Nln::CONSTRAINT::Group::get_solution_update_rms(
    const ::NOX::Abstract::Vector& xOld, const std::vector<double>& aTol,
    const std::vector<double>& rTol, const std::vector<NOX::Nln::StatusTest::QuantityType>& chQ,
    const std::vector<bool>& disable_implicit_weighting) const
{
  const ::NOX::Epetra::Vector& xOldEpetra = dynamic_cast<const ::NOX::Epetra::Vector&>(xOld);
  Teuchos::RCP<std::vector<double>> rms = Teuchos::rcp(new std::vector<double>(0));

  double rval = -1.0;
  for (std::size_t i = 0; i < chQ.size(); ++i)
  {
    rval = get_nln_req_interface_ptr()->get_primary_solution_update_rms(xVector.getEpetraVector(),
        xOldEpetra.getEpetraVector(), aTol[i], rTol[i], chQ[i], disable_implicit_weighting[i]);
    if (rval >= 0.0)
    {
      rms->push_back(rval);
      continue;
    }

    enum NOX::Nln::SolutionType soltype = NOX::Nln::Aux::ConvertQuantityType2SolutionType(chQ[i]);
    Teuchos::RCP<const NOX::Nln::CONSTRAINT::Interface::Required> constrptr =
        get_constraint_interface_ptr(soltype);

    rval = constrptr->get_lagrange_multiplier_update_rms(xVector.getEpetraVector(),
        xOldEpetra.getEpetraVector(), aTol[i], rTol[i], chQ[i], disable_implicit_weighting[i]);
    if (rval >= 0)
      rms->push_back(rval);
    else
    {
      std::ostringstream msg;
      msg << "The desired quantity"
             " for the \"NormWRMS\" Status Test could not be found! (enum="
          << chQ[i] << " | " << NOX::Nln::StatusTest::QuantityType2String(chQ[i]) << ")"
          << std::endl;
      throw_error("get_solution_update_rms", msg.str());
    }
  }

  return rms;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> NOX::Nln::CONSTRAINT::Group::get_current_active_set_map(
    const enum NOX::Nln::StatusTest::QuantityType& qt) const
{
  enum NOX::Nln::SolutionType soltype = NOX::Nln::Aux::ConvertQuantityType2SolutionType(qt);

  return get_constraint_interface_ptr(soltype)->get_current_active_set_map(qt);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> NOX::Nln::CONSTRAINT::Group::GetOldActiveSetMap(
    const enum NOX::Nln::StatusTest::QuantityType& qt) const
{
  enum NOX::Nln::SolutionType soltype = NOX::Nln::Aux::ConvertQuantityType2SolutionType(qt);

  return get_constraint_interface_ptr(soltype)->GetOldActiveSetMap(qt);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
enum ::NOX::StatusTest::StatusType NOX::Nln::CONSTRAINT::Group::GetActiveSetInfo(
    const enum NOX::Nln::StatusTest::QuantityType& qt, int& activeset_size) const
{
  enum NOX::Nln::SolutionType soltype = NOX::Nln::Aux::ConvertQuantityType2SolutionType(qt);

  return get_constraint_interface_ptr(soltype)->GetActiveSetInfo(qt, activeset_size);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::Nln::CONSTRAINT::Group::throw_error(
    const std::string& functionName, const std::string& errorMsg) const
{
  if (utils.isPrintType(::NOX::Utils::Error))
  {
    utils.out() << "ERROR - NOX::Nln::CONSTRAINT::Group::" << functionName << " - " << errorMsg
                << std::endl;
  }
  throw "NOX Error";
}

FOUR_C_NAMESPACE_CLOSE
