/*-----------------------------------------------------------*/
/*!
\file nox_nln_aux.cpp

\maintainer Michael Hiermeier

\date Jul 31, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_aux.H"
#include "nox_nln_statustest_normf.H"
#include "nox_nln_statustest_combo.H"

#include <Epetra_Vector.h>

#include <NOX_Abstract_ImplicitWeighting.H>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::AUX::RootMeanSquareNorm(
    const double& atol,
    const double& rtol,
    Teuchos::RCP<const Epetra_Vector> xnew,
    Teuchos::RCP<const Epetra_Vector> xincr,
    const bool& disable_implicit_weighting)
{
  double rval = 0.0;

  // calculate the old iterate (k-1)
  Teuchos::RCP<Epetra_Vector> v = Teuchos::rcp(new Epetra_Vector(*xnew));
  v->Update(-1.0,*xincr,1.0);

  // new auxiliary vector
  Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(xnew->Map(),false));

  // create the weighting factor u = RTOL |x^(k-1)| + ATOL
  u->PutScalar(1.0);
  u->Update(rtol,*v,atol);

  // v = xincr/u (elementwise)
  v->ReciprocalMultiply(1.0,*u,*xincr,0);

  // Turn off implicit scaling of norm if the vector supports it
  // ToDo Check if this makes any sense for pure Epetra_Vectors
  Teuchos::RCP<NOX::Abstract::ImplicitWeighting> iw_v;
  iw_v = Teuchos::rcp_dynamic_cast<NOX::Abstract::ImplicitWeighting>(v,false);
  bool saved_status = false;
  if (Teuchos::nonnull(iw_v) and disable_implicit_weighting)
  {
    saved_status = iw_v->getImplicitWeighting();
    iw_v->setImplicitWeighting(false);
  }

  // rval = sqrt (v * v / N)
  v->Norm2(&rval);
  rval /= std::sqrt(static_cast<double>(v->GlobalLength()));

  // Set the implicit scaling back to original value
  if (nonnull(iw_v) && disable_implicit_weighting)
    iw_v->setImplicitWeighting(saved_status);

  return rval;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::AUX::GetNormFClassVariable(
    const NOX::StatusTest::Generic& test,
    const NOX::NLN::StatusTest::QuantityType& qType,
    const std::string& classVariableName)
{
  // try to cast the given test to a NOX_StatusTest_Combo
  const NOX::NLN::StatusTest::Combo* comboTest =
      dynamic_cast<const NOX::NLN::StatusTest::Combo*>(&test);

  // if it is no combo test, we just have to check for the desired type
  if (comboTest == 0)
  {
    const NOX::NLN::StatusTest::NormF* normFTest =
          dynamic_cast<const NOX::NLN::StatusTest::NormF*>(&test);

    // no normF StatusTest...
    if (normFTest==0)
      return -1.0;
    // yeah we found one...
    else
    {
      // look for the right quantity and get the desired class variable value
      if (classVariableName == "NormF")
        return normFTest->GetNormF(qType);
      else if (classVariableName == "TrueTolerance")
        return normFTest->GetTrueTolerance(qType);
      else if (classVariableName == "SpecifiedTolerance")
        return normFTest->GetSpecifiedTolerance(qType);
      else if (classVariableName == "InitialTolerance")
        return normFTest->GetInitialTolerance(qType);
    }
  }
  // if the nox_nln_statustest_combo Test cast was successful
  else
  {
    const std::vector<Teuchos::RCP<NOX::StatusTest::Generic> >& tests =
        comboTest->GetTestVector();
    double ret = -1.0;
    for (std::size_t i=0;i<tests.size();++i)
    {
      // recursive function call
      ret = GetNormFClassVariable(*(tests[i]),qType,classVariableName);
      if (ret!=-1.0)
        return ret;
    }
  }

  // default return
  return -1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::NLN::SolutionType NOX::NLN::AUX::ConvertQuantityType2SolutionType(
      const enum NOX::NLN::StatusTest::QuantityType& qtype)
{
  enum NOX::NLN::SolutionType soltype = NOX::NLN::sol_unknown;
  switch(qtype)
  {
    case NOX::NLN::StatusTest::quantity_structure:
    case NOX::NLN::StatusTest::quantity_eas:
    case NOX::NLN::StatusTest::quantity_pressure:
      soltype = NOX::NLN::sol_structure;
      break;
    case NOX::NLN::StatusTest::quantity_lag_pen_constraint:
      soltype = NOX::NLN::sol_lag_pen_constraint;
      break;
    case NOX::NLN::StatusTest::quantity_contact:
      soltype = NOX::NLN::sol_contact;
      break;
    case NOX::NLN::StatusTest::quantity_meshtying:
      soltype = NOX::NLN::sol_meshtying;
      break;
    case NOX::NLN::StatusTest::quantity_windkessel:
      soltype = NOX::NLN::sol_windkessel;
      break;
    case NOX::NLN::StatusTest::quantity_plasticity:
    case NOX::NLN::StatusTest::quantity_unknown:
    default:
      dserror("Unknown conversion for the quantity type \"%s\".",
          NOX::NLN::StatusTest::QuantityType2String(qtype).c_str());
  }
  // return the corresponding solution type
  return soltype;
}
