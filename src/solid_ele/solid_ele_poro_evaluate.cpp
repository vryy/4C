/*----------------------------------------------------------------------*/
/*! \file

\brief A C++ wrapper for the solid-poro element

This file contains the element-specific evaluation routines such as
Evaluate(...), EvaluateNeumann(...), etc.

\level 1
*/

#include "solid_ele_poro.H"
#include "solid_ele_poro_factory.H"
#include "solid_ele.H"
#include "solid_ele_factory.H"
#include "solid_ele_calc_lib.H"


int DRT::ELEMENTS::SolidPoro::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  if (!material_post_setup_)
  {
    DRT::ELEMENTS::SolidFactory::ProvideImpl(
        this, GetEleTech(), GetEleKinematicType(), GetEAStype())
        ->MaterialPostSetup(*this, this->StructPoroMaterial());
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
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, GetEleTech(), GetEleKinematicType(), GetEAStype())
          ->EvaluateNonlinearForceStiffnessMass(*this, this->StructPoroMaterial(), discretization,
              la[0].lm_, params, &elevec1, &elemat1, nullptr);

      if (la.Size() > 2 and this->NumMaterial() > 1)
      {
        if (discretization.HasState(1, "porofluid"))
        {
          DRT::ELEMENTS::SolidPoroFactory::ProvideImpl(this, this->GetElePoroType())
              ->EvaluateNonlinearForceStiffness(*this, this->StructPoroMaterial(),
                  this->FluidPoroMultiMaterial(), this->GetEleKinematicType(), discretization, la,
                  params, &elevec1, &elemat1);
        }
      }
      return 0;
    }
    case struct_calc_internalforce:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, GetEleTech(), GetEleKinematicType(), GetEAStype())
          ->EvaluateNonlinearForceStiffnessMass(*this, this->StructPoroMaterial(), discretization,
              la[0].lm_, params, &elevec1, nullptr, nullptr);

      if (la.Size() > 2 and this->NumMaterial() > 1)
      {
        if (discretization.HasState(1, "porofluid"))
        {
          DRT::ELEMENTS::SolidPoroFactory::ProvideImpl(this, this->GetElePoroType())
              ->EvaluateNonlinearForceStiffness(*this, this->StructPoroMaterial(),
                  this->FluidPoroMultiMaterial(), this->GetEleKinematicType(), discretization, la,
                  params, &elevec1, nullptr);
        }
      }
      return 0;
    }
    case struct_calc_nlnstiffmass:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, GetEleTech(), GetEleKinematicType(), GetEAStype())
          ->EvaluateNonlinearForceStiffnessMass(*this, this->StructPoroMaterial(), discretization,
              la[0].lm_, params, &elevec1, &elemat1, &elemat2);

      // we skip this evaluation if the coupling is not setup yet, i.e.
      // if the secondary dofset or the secondary material was not set
      // this can happen during setup of the time integrator or restart
      // there might be a better way. For instance do not evaluate
      // before the setup of the multiphysics problem is completed.
      if (la.Size() > 2 and this->NumMaterial() > 1)
      {
        if (discretization.HasState(1, "porofluid"))
        {
          DRT::ELEMENTS::SolidPoroFactory::ProvideImpl(this, this->GetElePoroType())
              ->EvaluateNonlinearForceStiffness(*this, this->StructPoroMaterial(),
                  this->FluidPoroMultiMaterial(), this->GetEleKinematicType(), discretization, la,
                  params, &elevec1, &elemat1);
        }
      }
      return 0;
    }
    case struct_calc_nlnstifflmass:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, GetEleTech(), GetEleKinematicType(), GetEAStype())
          ->EvaluateNonlinearForceStiffnessMass(*this, this->StructPoroMaterial(), discretization,
              la[0].lm_, params, &elevec1, &elemat1, &elemat2);
      DRT::ELEMENTS::LumpMatrix(elemat2);
      return 0;
    }
    case struct_poro_calc_scatracoupling:
    {
      // no coupling -> return
      return 0;
    }
    case struct_poro_calc_fluidcoupling:
    {
      if (la.Size() > 2)
      {
        if (discretization.HasState(1, "porofluid"))
        {
          DRT::ELEMENTS::SolidPoroFactory::ProvideImpl(this, this->GetElePoroType())
              ->CouplingPoroelast(*this, this->StructPoroMaterial(), this->FluidPoroMultiMaterial(),
                  this->GetEleKinematicType(), discretization, la, params, elemat1);
        }
      }
      return 0;
    }
    case DRT::ELEMENTS::struct_calc_update_istep:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, GetEleTech(), GetEleKinematicType(), GetEAStype())
          ->Update(*this, SolidPoroMaterial(), discretization, la[0].lm_, params);
      return 0;
    }
    case DRT::ELEMENTS::struct_calc_recover:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, GetEleTech(), GetEleKinematicType(), GetEAStype())
          ->Recover(*this, discretization, la[0].lm_, params);
      return 0;
    }
    case struct_calc_stress:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, GetEleTech(), GetEleKinematicType(), GetEAStype())
          ->CalculateStress(*this, this->StructPoroMaterial(),
              StressIO{GetIOStressType(*this, params), GetMutableStressData(*this, params)},
              StrainIO{GetIOStrainType(*this, params), GetMutableStrainData(*this, params)},
              discretization, la[0].lm_, params);

      if (la.Size() > 2)
      {
        if (discretization.HasState(1, "porofluid"))
        {
          DRT::ELEMENTS::SolidPoroFactory::ProvideImpl(this, this->GetElePoroType())
              ->CouplingStress(*this, discretization, la[0].lm_, params);
        }
      }
      return 0;
    }
    case struct_init_gauss_point_data_output:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, GetEleTech(), GetEleKinematicType(), GetEAStype())
          ->InitializeGaussPointDataOutput(*this, SolidPoroMaterial(),
              *ParamsInterface().MutableGaussPointDataOutputManagerPtr());
      return 0;
    }
    case struct_gauss_point_data_output:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, GetEleTech(), GetEleKinematicType(), GetEAStype())
          ->EvaluateGaussPointDataOutput(*this, SolidPoroMaterial(),
              *ParamsInterface().MutableGaussPointDataOutputManagerPtr());
      return 0;
    }
    case DRT::ELEMENTS::struct_calc_predict:
    {
      // do nothing for now
      return 0;
    }
    default:
      dserror("The element action %s is not yet implemented for the new solid elements",
          ActionType2String(action).c_str());
      return 0;
  }
}

int DRT::ELEMENTS::SolidPoro::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  dserror("not implemented");
  return 1;
}