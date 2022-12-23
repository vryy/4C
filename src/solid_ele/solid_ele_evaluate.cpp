/*----------------------------------------------------------------------*/
/*! \file

\brief A C++ wrapper for the solid element

This file contains the element-specific evaluation routines such as
Evaluate(...), EvaluateNeumann(...), etc.

\level 1
*/

#include "elements_paramsinterface.H"
#include "str_elements_paramsinterface.H"
#include "dserror.H"
#include "solid_ele.H"
#include "solid_ele_factory.H"
#include "solid_ele_interface.H"

namespace
{
  void LumpMatrix(Epetra_SerialDenseMatrix& matrix)
  {
    dsassert(matrix.N() == matrix.M(), "The provided mass matrix is not a square matrix!");

    // we assume mass is a square matrix
    for (int c = 0; c < matrix.N(); ++c)  // parse columns
    {
      double d = 0.0;
      for (int r = 0; r < matrix.M(); ++r)  // parse rows
      {
        d += matrix(r, c);  // accumulate row entries
        matrix(r, c) = 0.0;
      }
      matrix(c, c) = d;  // apply sum of row entries on diagonal
    }
  }
}  // namespace

int DRT::ELEMENTS::Solid::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  if (!material_post_setup_)
  {
    DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->MaterialPostSetup(*this);
    material_post_setup_ = true;
  }

  // get ptr to interface to time integration
  SetParamsInterfacePtr(params);

  const ELEMENTS::ActionType action = std::invoke(
      [&]()
      {
        if (IsParamsInterface())
          return ParamsInterface().GetActionType();
        else
          return String2ActionType(params.get<std::string>("action", "none"));
      });

  switch (action)
  {
    case DRT::ELEMENTS::struct_calc_nlnstiff:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->EvaluateNonlinearForceStiffnessMass(
          *this, discretization, lm, params, &elevec1, &elemat1, nullptr);
      return 0;
    }
    case struct_calc_internalforce:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->EvaluateNonlinearForceStiffnessMass(
          *this, discretization, lm, params, &elevec1, nullptr, nullptr);
      return 0;
    }
    case struct_calc_nlnstiffmass:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->EvaluateNonlinearForceStiffnessMass(
          *this, discretization, lm, params, &elevec1, &elemat1, &elemat2);
      return 0;
    }
    case struct_calc_nlnstifflmass:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->EvaluateNonlinearForceStiffnessMass(
          *this, discretization, lm, params, &elevec1, &elemat1, &elemat2);
      LumpMatrix(elemat2);
      return 0;
    }
    case DRT::ELEMENTS::struct_calc_update_istep:
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->Update(*this, discretization, lm, params);
      return 0;
    case DRT::ELEMENTS::struct_calc_recover:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->Recover(*this, discretization, lm, params);
      return 0;
    }
    case struct_calc_stress:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->CalculateStress(
          *this, discretization, lm, params);
      return 0;
    }
    case struct_postprocess_stress:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->PostProcessStressStrain(
          *this, discretization, lm, params);
      return 0;
    }
    case struct_init_gauss_point_data_output:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->InitializeGaussPointDataOutput(*this);
      return 0;
    }
    case struct_gauss_point_data_output:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->EvaluateGaussPointDataOutput(*this);
      return 0;
    }
    case DRT::ELEMENTS::struct_calc_predict:
      // do nothing for now
      return 0;
    default:
      dserror("The element action %s is not yet implemented for the new solid elements",
          ActionType2String(action).c_str());
  }

  return 0;
}
int DRT::ELEMENTS::Solid::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  dserror("not implemented");
  return 1;
}
