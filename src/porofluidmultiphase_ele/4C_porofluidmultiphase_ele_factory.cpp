/*----------------------------------------------------------------------*/
/*! \file
 \brief factory class providing the implementations of the porofluidmultiphase
        element evaluation routines

   \level 3

 *----------------------------------------------------------------------*/



#include "4C_porofluidmultiphase_ele_factory.hpp"

#include "4C_global_data.hpp"
#include "4C_porofluidmultiphase_ele_calc.hpp"

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------------*
 | provide the implementation of evaluation class      (public) vuong 08/16 |
 *--------------------------------------------------------------------------*/
Discret::ELEMENTS::PoroFluidMultiPhaseEleInterface*
Discret::ELEMENTS::PoroFluidMultiPhaseFactory::ProvideImpl(
    Core::FE::CellType distype, const int numdofpernode, const std::string& disname)
{
  // -------------------------------------- number of degrees of freedom
  // number of degrees of freedom
  static const int ndim = Global::Problem::Instance()->NDim();

  switch (distype)
  {
    case Core::FE::CellType::quad4:
    {
      if (ndim == 2)
        return define_problem_type<Core::FE::CellType::quad4>(numdofpernode, disname);
      else
        FOUR_C_THROW("invalid problem dimension for quad4 porofluidmultiphase element!");
      break;
    }
      //  case Core::FE::CellType::quad8:
      //  {
      //    if(ndim==2)
      //      return
      //      define_problem_type<Core::FE::CellType::quad8>(numdofpernode,disname);
      //    else
      //      FOUR_C_THROW("invalid problem dimension for quad8 porofluidmultiphase element!");
      //    break;
      //  }
    case Core::FE::CellType::quad9:
    {
      if (ndim == 2)
        return define_problem_type<Core::FE::CellType::quad9>(numdofpernode, disname);
      else
        FOUR_C_THROW("invalid problem dimension for quad9 porofluidmultiphase element!");
      break;
    }
    case Core::FE::CellType::tri3:
    {
      if (ndim == 2)
        return define_problem_type<Core::FE::CellType::tri3>(numdofpernode, disname);
      else
        FOUR_C_THROW("invalid problem dimension for tri3 porofluidmultiphase element!");
      break;
    }
      //  case Core::FE::CellType::tri6:
      //  {
      //    if(ndim==2)
      //      return
      //      define_problem_type<Core::FE::CellType::tri6>(numdofpernode,disname);
      //    else
      //      FOUR_C_THROW("invalid problem dimension for tri6 porofluidmultiphase element!");
      //    break;
      //  }
    case Core::FE::CellType::line2:
    {
      if (ndim == 1)
        return define_problem_type<Core::FE::CellType::line2>(numdofpernode, disname);
      else
        FOUR_C_THROW("invalid problem dimension for line2 porofluidmultiphase element!");
      break;
    }
    case Core::FE::CellType::line3:
    {
      if (ndim == 1)
        return define_problem_type<Core::FE::CellType::line3>(numdofpernode, disname);
      else
        FOUR_C_THROW("invalid problem dimension for line3 porofluidmultiphase element!");
      break;
    }
    case Core::FE::CellType::hex8:
    {
      if (ndim == 3)
        return define_problem_type<Core::FE::CellType::hex8>(numdofpernode, disname);
      else
        FOUR_C_THROW("invalid problem dimension for hex8 porofluidmultiphase element!");
      break;
    }
    case Core::FE::CellType::hex27:
    {
      if (ndim == 3)
        return define_problem_type<Core::FE::CellType::hex27>(numdofpernode, disname);
      else
        FOUR_C_THROW("invalid problem dimension for hex27 porofluidmultiphase element!");
      break;
    }
    case Core::FE::CellType::tet4:
    {
      if (ndim == 3)
        return define_problem_type<Core::FE::CellType::tet4>(numdofpernode, disname);
      else
        FOUR_C_THROW("invalid problem dimension for tet4 porofluidmultiphase element!");
      break;
    }
    case Core::FE::CellType::tet10:
    {
      if (ndim == 3)
        return define_problem_type<Core::FE::CellType::tet10>(numdofpernode, disname);
      else
        FOUR_C_THROW("invalid problem dimension for tet10 porofluidmultiphase element!");
      break;
    }
    default:
      FOUR_C_THROW("Element shape %s not activated. Just do it.",
          Core::FE::CellTypeToString(distype).c_str());
      break;
  }
  return nullptr;
}

/*--------------------------------------------------------------------------*
 | provide the implementation of evaluation class      (public) vuong 08/16 |
 *--------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::PoroFluidMultiPhaseEleInterface*
Discret::ELEMENTS::PoroFluidMultiPhaseFactory::define_problem_type(
    const int numdofpernode, const std::string& disname)
{
  return Discret::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::Instance(numdofpernode, disname);
}

FOUR_C_NAMESPACE_CLOSE
