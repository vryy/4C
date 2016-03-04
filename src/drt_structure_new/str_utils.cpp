/*-----------------------------------------------------------*/
/*!
\file str_utils.cpp

\maintainer Michael Hiermeier

\date Dec 2, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_utils.H"

#include "../drt_lib/drt_dserror.H"

#include "../linalg/linalg_utils.H"

#include <Epetra_Vector.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::Abstract::Vector::NormType STR::NLN::Convert2NoxNormType(
    const enum INPAR::STR::VectorNorm& normtype)
{
  enum NOX::Abstract::Vector::NormType nox_normtype =
      NOX::Abstract::Vector::TwoNorm;

  switch (normtype)
  {
    case INPAR::STR::norm_l2:
      nox_normtype = NOX::Abstract::Vector::TwoNorm;
      break;
    case INPAR::STR::norm_l1:
      nox_normtype = NOX::Abstract::Vector::OneNorm;
      break;
    case INPAR::STR::norm_inf:
      nox_normtype = NOX::Abstract::Vector::MaxNorm;
      break;
    case INPAR::STR::norm_rms:
    case INPAR::STR::norm_vague:
    default:
      dserror("Unknown conversion for the given vector norm type: \" %s \"!",
          INPAR::STR::VectorNormString(normtype).c_str());
      break;
  } // switch case normtype

  return nox_normtype;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::ConvertModelType2SolType(
    std::vector<enum NOX::NLN::SolutionType>& soltypes,
    std::map<enum NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& slinsolvers,
    const std::set<enum INPAR::STR::ModelType>& modeltypes,
    const std::map<enum INPAR::STR::ModelType,Teuchos::RCP<LINALG::Solver> >& mlinsolvers)
{
  // initialize the vector and/or force the length to zero
  if (soltypes.size()>0)
  {
    soltypes.clear();
    slinsolvers.clear();
  }

  // pre-set the vector size
  soltypes.reserve(modeltypes.size());

  // The strings of the different enums have to fit!
  std::set<enum INPAR::STR::ModelType>::const_iterator mt_iter;
  for (mt_iter=modeltypes.begin();mt_iter!=modeltypes.end();++mt_iter)
  {
    const std::string name = INPAR::STR::ModelTypeString(*mt_iter);
    const enum NOX::NLN::SolutionType soltype =
        NOX::NLN::String2SolutionType(name);

    // check if the corresponding enum could be found.
    if (soltype == NOX::NLN::sol_unknown)
      dserror("The corresponding solution-type was not found. "
          "Given string: %s", name.c_str());

    soltypes.push_back(soltype);
    // copy the linsolver pointers into the new map
    if (mlinsolvers.find(*mt_iter) != mlinsolvers.end())
      slinsolvers[soltype] = mlinsolvers.at(*mt_iter);
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::NLN::GlobalData::OptimizationProblemType
    STR::NLN::OptimizationType(
    const std::vector<enum NOX::NLN::SolutionType>& soltypes)
{
  enum NOX::NLN::GlobalData::OptimizationProblemType opttype =
      NOX::NLN::GlobalData::opt_unconstrained;
  std::vector<enum NOX::NLN::SolutionType>::const_iterator st_iter;

  for (st_iter=soltypes.begin();st_iter!=soltypes.end();++st_iter)
  {
    switch (*st_iter)
    {
      // -----------------------------------
      // Inequality constraint
      // active set decision necessary
      // saddle point structure or condensed
      // -----------------------------------
      case NOX::NLN::sol_contact:
        return NOX::NLN::GlobalData::opt_inequality_constrained;
        break;
      // -----------------------------------
      // Equality constraint
      // no active set decision necessary
      // saddle point structure or condensed
      // -----------------------------------
      case NOX::NLN::sol_meshtying:
      case NOX::NLN::sol_lag_pen_constraint:
      case NOX::NLN::sol_windkessel:
        opttype = NOX::NLN::GlobalData::opt_equality_constrained;
        break;
      // -----------------------------------
      // Unconstrained problem
      // pure structural problem
      // no saddle point structure
      // -----------------------------------
      case NOX::NLN::sol_structure:
      case NOX::NLN::sol_springdashpot:
      default:
        break;
    }
  }

  return opttype;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::CreateConstraintInterfaces(
    std::map<enum NOX::NLN::SolutionType,Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> >& iconstr,
    const std::vector<enum NOX::NLN::SolutionType>& soltypes)
{
  if (iconstr.size()>0)
    iconstr.clear();

  std::vector<enum NOX::NLN::SolutionType>::const_iterator st_iter;
  for (st_iter=soltypes.begin();st_iter!=soltypes.end();++st_iter)
  {
    switch (*st_iter)
    {
      case NOX::NLN::sol_contact:
//        iconstr_[NOX::NLN::sol_contact] = itimint.GetContactManager();
//        break;
      case NOX::NLN::sol_windkessel:
//        iconstr_[NOX::NLN::sol_windkessel] = itimint.GetWindkesselManager();
//        break;
      case NOX::NLN::sol_lag_pen_constraint:
//        iconstr_[NOX::NLN::sol_windkessel] = itimint.GetLagPenConstrManager();
//        break;
      default:
        break;
    }
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::AssembleForce(Epetra_Vector& target,
    const Epetra_Vector& source, const double& scalar_source)
{
  Teuchos::RCP<const Epetra_Vector> source_ptr;
  // brief check if a export is necessary
  if (not target.Map().PointSameAs(source.Map()))
  {
    // brief check if the source map is a sub-map of the target map
    if (target.Map().NumMyElements() < source.Map().NumMyElements())
      dserror("The target map should be greater or equal than the "
          "source map on each single processor, since the source map "
          "is supposed to be a sub-map of the target map.");

    Teuchos::RCP<Epetra_Vector> source_exp_ptr =
        Teuchos::rcp(new Epetra_Vector(target.Map()));
    LINALG::Export(source,*source_exp_ptr);
    source_ptr = source_exp_ptr;
  }
  else
    source_ptr = Teuchos::rcp(&source,false);

  if (target.Update(scalar_source,*source_ptr,1.0))
    dserror("Update failed");

  return;
}
