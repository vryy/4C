/*-----------------------------------------------------------*/
/*! \file

\brief evaluation routines for the fluid poro element


\level 2

*/
/*-----------------------------------------------------------*/


#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_ele_factory.hpp"
#include "4C_fluid_ele_interface.hpp"
#include "4C_fluid_ele_parameter_poro.hpp"
#include "4C_fluid_ele_poro.hpp"
#include "4C_global_data.hpp"

FOUR_C_NAMESPACE_OPEN

void Discret::ELEMENTS::FluidPoroEleType::pre_evaluate(Core::FE::Discretization& dis,
    Teuchos::ParameterList& p, Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix1,
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix2,
    Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
    Teuchos::RCP<Epetra_Vector> systemvector3)
{
  const auto action = Core::UTILS::GetAsEnum<FLD::Action>(p, "action");

  // poro specific actions
  if (action == FLD::set_poro_parameter)
  {
    Discret::ELEMENTS::FluidEleParameterPoro* fldpara =
        Discret::ELEMENTS::FluidEleParameterPoro::instance();
    fldpara->set_element_poro_parameter(p, dis.get_comm().MyPID());
  }
  else
  {
    // call standard fluid type
    FluidType::pre_evaluate(
        dis, p, systemmatrix1, systemmatrix2, systemvector1, systemvector2, systemvector3);
  }
}

int Discret::ELEMENTS::FluidPoro::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  // get the action required
  const auto act = Core::UTILS::GetAsEnum<FLD::Action>(params, "action");

  // get material
  Teuchos::RCP<Core::Mat::Material> mat = material();

  // switch between different physical types as used below
  std::string impltype = "poro";
  switch (params.get<int>("Physical Type", Inpar::FLUID::physicaltype_undefined))
  {
    case Inpar::FLUID::poro:
      impltype = "poro";
      break;
    case Inpar::FLUID::poro_p1:
      impltype = "poro_p1";
      break;
    default:
      FOUR_C_THROW("invalid physical type for porous fluid!");
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
      return Discret::ELEMENTS::FluidFactory::provide_impl(shape(), impltype)
          ->evaluate(
              this, discretization, lm, params, mat, elemat1, elemat2, elevec1, elevec2, elevec3);
    }
    //-----------------------------------------------------------------------
    // standard implementation enabling time-integration schemes such as
    // one-step-theta, BDF2, and generalized-alpha (n+alpha_F and n+1)
    // for the particular case of porous flow
    //-----------------------------------------------------------------------
    case FLD::calc_porousflow_fluid_coupling:
    {
      return Discret::ELEMENTS::FluidFactory::provide_impl(shape(), impltype)
          ->evaluate(this, discretization, lm, params, mat, elemat1, elemat2, elevec1, elevec2,
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
      switch (params.get<int>("Physical Type", Inpar::FLUID::physicaltype_undefined))
      {
        case Inpar::FLUID::poro:
        {
          // no coupling
          return 0;
        }
        default:
          FOUR_C_THROW(
              "Invalid physical type for monolithic poroelasticity with scalar transport\n");
          break;
      }
    }
    break;
    case FLD::calc_volume:
    {
      return Discret::ELEMENTS::FluidFactory::provide_impl(shape(), impltype)
          ->evaluate_service(
              this, params, mat, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);
    }
    case FLD::set_poro_parameter:
      break;
    default:
      // call evaluate of standard fluid
      return Fluid::evaluate(
          params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);
  }

  return 0;
}

FOUR_C_NAMESPACE_CLOSE
