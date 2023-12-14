/*---------------------------------------------------------------------*/
/*! \file

\brief Incomplete! - Purpose: Internal implementation of RedAirBloodScatraLine3 element


\level 3

*/
/*---------------------------------------------------------------------*/



#include "baci_red_airways_air_blood_scatraLine3_impl.H"

#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"
#include "baci_discretization_fem_general_utils_gder2.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_utils.H"
#include "baci_mat_air_0d_O2_saturation.H"
#include "baci_mat_hemoglobin_0d_O2_saturation.H"
#include "baci_mat_list.H"
#include "baci_mat_newtonianfluid.H"
#include "baci_mat_par_bundle.H"
#include "baci_red_airways_evaluation_data.h"
#include "baci_utils_function.H"

BACI_NAMESPACE_OPEN

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
      dserror("shape %d (%d nodes) not supported", red_acinus->Shape(), red_acinus->NumNode());
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
    CORE::LINALG::SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat)
{
  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::Initial(RedAirBloodScatraLine3* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
    Teuchos::RCP<const MAT::Material> material)
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
    DRT::Discretization& discretization, std::vector<int>& lm, Teuchos::RCP<MAT::Material> material)
{
}

BACI_NAMESPACE_CLOSE
