/*-----------------------------------------------------------*/
/*! \file

\brief evaluation routines for the fluid poro element


\level 2

*/
/*-----------------------------------------------------------*/


#include "baci_fluid_ele_action.hpp"
#include "baci_fluid_ele_factory.hpp"
#include "baci_fluid_ele_interface.hpp"
#include "baci_fluid_ele_parameter_poro.hpp"
#include "baci_fluid_ele_poro.hpp"
#include "baci_global_data.hpp"

BACI_NAMESPACE_OPEN

void DRT::ELEMENTS::FluidPoroEleType::PreEvaluate(DRT::Discretization& dis,
    Teuchos::ParameterList& p, Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix2,
    Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
    Teuchos::RCP<Epetra_Vector> systemvector3)
{
  const auto action = CORE::UTILS::GetAsEnum<FLD::Action>(p, "action");

  // poro specific actions
  if (action == FLD::set_poro_parameter)
  {
    DRT::ELEMENTS::FluidEleParameterPoro* fldpara =
        DRT::ELEMENTS::FluidEleParameterPoro::Instance();
    fldpara->SetElementPoroParameter(p, dis.Comm().MyPID());
  }
  else
  {
    // call standard fluid type
    FluidType::PreEvaluate(
        dis, p, systemmatrix1, systemmatrix2, systemvector1, systemvector2, systemvector3);
  }
}

int DRT::ELEMENTS::FluidPoro::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  // get the action required
  const auto act = CORE::UTILS::GetAsEnum<FLD::Action>(params, "action");

  // get material
  Teuchos::RCP<MAT::Material> mat = Material();

  // switch between different physical types as used below
  std::string impltype = "poro";
  switch (params.get<int>("Physical Type", INPAR::FLUID::physicaltype_undefined))
  {
    case INPAR::FLUID::poro:
      impltype = "poro";
      break;
    case INPAR::FLUID::poro_p1:
      impltype = "poro_p1";
      break;
    default:
      dserror("invalid physical type for porous fluid!");
      break;
  }

  switch (act)
  {
    //-----------------------------------------------------------------------
    // standard implementation enabling time-integration schemes such as
    // one-step-theta, BDF2, and generalized-alpha (n+alpha_F and n+1)
    //-----------------------------------------------------------------------
    case FLD::calc_fluid_systemmat_and_residual:
    {
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), impltype)
          ->Evaluate(
              this, discretization, lm, params, mat, elemat1, elemat2, elevec1, elevec2, elevec3);
    }
    //-----------------------------------------------------------------------
    // standard implementation enabling time-integration schemes such as
    // one-step-theta, BDF2, and generalized-alpha (n+alpha_F and n+1)
    // for the particular case of porous flow
    //-----------------------------------------------------------------------
    case FLD::calc_porousflow_fluid_coupling:
    {
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), impltype)
          ->Evaluate(this, discretization, lm, params, mat, elemat1, elemat2, elevec1, elevec2,
              elevec3, true);
    }
    //-----------------------------------------------------------------------
    // standard implementation enabling time-integration schemes such as
    // one-step-theta, BDF2, and generalized-alpha (n+alpha_F and n+1)
    // for evaluation of off-diagonal matrix block for monolithic
    // porous flow with scalar transport
    //-----------------------------------------------------------------------
    case FLD::calc_poroscatra_mono_odblock:
    {
      switch (params.get<int>("Physical Type", INPAR::FLUID::physicaltype_undefined))
      {
        case INPAR::FLUID::poro:
        {
          // no coupling
          return 0;
        }
        default:
          dserror("Invalid physical type for monolithic poroelasticity with scalar transport\n");
          break;
      }
    }
    break;
    case FLD::calc_volume:
    {
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), impltype)
          ->EvaluateService(
              this, params, mat, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);
    }
    case FLD::set_poro_parameter:
      break;
    default:
      // call evaluate of standard fluid
      return Fluid::Evaluate(
          params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);
  }

  return 0;
}

BACI_NAMESPACE_CLOSE
