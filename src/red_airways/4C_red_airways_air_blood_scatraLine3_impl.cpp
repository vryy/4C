/*---------------------------------------------------------------------*/
/*! \file

\brief Incomplete! - Purpose: Internal implementation of RedAirBloodScatraLine3 element


\level 3

*/
/*---------------------------------------------------------------------*/



#include "4C_red_airways_air_blood_scatraLine3_impl.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_general_utils_gder2.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_newtonianfluid.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_red_airways_evaluation_data.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::RedAirBloodScatraLine3ImplInterface*
Discret::ELEMENTS::RedAirBloodScatraLine3ImplInterface::Impl(
    Discret::ELEMENTS::RedAirBloodScatraLine3* red_acinus)
{
  switch (red_acinus->Shape())
  {
    case Core::FE::CellType::line3:
    {
      static RedAirBloodScatraLine3Impl<Core::FE::CellType::line3>* acinus;
      if (acinus == nullptr)
      {
        acinus = new RedAirBloodScatraLine3Impl<Core::FE::CellType::line3>;
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
template <Core::FE::CellType distype>
Discret::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::RedAirBloodScatraLine3Impl()
{
}

/*----------------------------------------------------------------------*
 | evaluate (public)                                       ismail 01/10 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::evaluate(RedAirBloodScatraLine3* ele,
    Teuchos::ParameterList& params, Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
    Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
    Core::LinAlg::SerialDenseVector& elevec1_epetra,
    Core::LinAlg::SerialDenseVector& elevec2_epetra,
    Core::LinAlg::SerialDenseVector& elevec3_epetra, Teuchos::RCP<Core::Mat::Material> mat)
{
  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::Initial(RedAirBloodScatraLine3* ele,
    Teuchos::ParameterList& params, Core::FE::Discretization& discretization, std::vector<int>& lm,
    Teuchos::RCP<const Core::Mat::Material> material)
{
  ReducedLung::EvaluationData& evaluation_data = ReducedLung::EvaluationData::get();

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
template <Core::FE::CellType distype>
void Discret::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::GetCoupledValues(
    RedAirBloodScatraLine3* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Teuchos::RCP<Core::Mat::Material> material)
{
}

FOUR_C_NAMESPACE_CLOSE
