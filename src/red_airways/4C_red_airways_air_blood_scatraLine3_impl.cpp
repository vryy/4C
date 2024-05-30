/*---------------------------------------------------------------------*/
/*! \file

\brief Incomplete! - Purpose: Internal implementation of RedAirBloodScatraLine3 element


\level 3

*/
/*---------------------------------------------------------------------*/



#include "4C_red_airways_air_blood_scatraLine3_impl.hpp"

#include "4C_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_discretization_fem_general_utils_gder2.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_newtonianfluid.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_red_airways_evaluation_data.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAirBloodScatraLine3ImplInterface*
DRT::ELEMENTS::RedAirBloodScatraLine3ImplInterface::Impl(
    DRT::ELEMENTS::RedAirBloodScatraLine3* red_acinus)
{
  switch (red_acinus->Shape())
  {
    case CORE::FE::CellType::line3:
    {
      static RedAirBloodScatraLine3Impl<CORE::FE::CellType::line3>* acinus;
      if (acinus == nullptr)
      {
        acinus = new RedAirBloodScatraLine3Impl<CORE::FE::CellType::line3>;
      }
      return acinus;
    }
    default:
      FOUR_C_THROW(
          "shape %d (%d nodes) not supported", red_acinus->Shape(), red_acinus->num_node());
  }
  return nullptr;
}



/*----------------------------------------------------------------------*
  | constructor (public)                                    ismail 01/10 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::RedAirBloodScatraLine3Impl()
{
}

/*----------------------------------------------------------------------*
 | evaluate (public)                                       ismail 01/10 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
int DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::Evaluate(RedAirBloodScatraLine3* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra, Teuchos::RCP<CORE::MAT::Material> mat)
{
  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::Initial(RedAirBloodScatraLine3* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
    Teuchos::RCP<const CORE::MAT::Material> material)
{
  REDAIRWAYS::EvaluationData& evaluation_data = REDAIRWAYS::EvaluationData::get();

  //--------------------------------------------------------------------
  // get the generation numbers
  //--------------------------------------------------------------------
  //  if(myrank == ele->Owner())
  {
    int gid = ele->Id();
    double val = -2.0;
    evaluation_data.generations->ReplaceGlobalValues(1, &val, &gid);
  }

}  // RedAirBloodScatraLine3Impl::Initial

/*----------------------------------------------------------------------*
 |  Get the coupled the values on the coupling interface    ismail 07/10|
 |  of the 3D/reduced-D problem                                         |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::GetCoupledValues(
    RedAirBloodScatraLine3* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    Teuchos::RCP<CORE::MAT::Material> material)
{
}

FOUR_C_NAMESPACE_CLOSE
