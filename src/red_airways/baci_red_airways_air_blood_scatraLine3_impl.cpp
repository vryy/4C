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
#include "baci_lib_function.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_utils.H"
#include "baci_mat_air_0d_O2_saturation.H"
#include "baci_mat_hemoglobin_0d_O2_saturation.H"
#include "baci_mat_list.H"
#include "baci_mat_newtonianfluid.H"
#include "baci_mat_par_bundle.H"
#include "baci_red_airways_evaluation_data.h"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAirBloodScatraLine3ImplInterface*
DRT::ELEMENTS::RedAirBloodScatraLine3ImplInterface::Impl(
    DRT::ELEMENTS::RedAirBloodScatraLine3* red_acinus)
{
  switch (red_acinus->Shape())
  {
    case DRT::Element::DiscretizationType::line3:
    {
      static RedAirBloodScatraLine3Impl<DRT::Element::DiscretizationType::line3>* acinus;
      if (acinus == nullptr)
      {
        acinus = new RedAirBloodScatraLine3Impl<DRT::Element::DiscretizationType::line3>;
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
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::RedAirBloodScatraLine3Impl()
{
}

/*----------------------------------------------------------------------*
 | evaluate (public)                                       ismail 01/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
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
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::Initial(RedAirBloodScatraLine3* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
    Teuchos::RCP<const MAT::Material> material)
{
  DRT::REDAIRWAYS::EvaluationData& evaluation_data = DRT::REDAIRWAYS::EvaluationData::get();

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
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::Sysmat(RedAirBloodScatraLine3* ele,
    CORE::LINALG::SerialDenseVector& epnp, CORE::LINALG::SerialDenseVector& epn,
    CORE::LINALG::SerialDenseVector& epnm, CORE::LINALG::SerialDenseMatrix& sysmat,
    CORE::LINALG::SerialDenseVector& rhs, Teuchos::RCP<const MAT::Material> material,
    Teuchos::ParameterList& params, double time, double dt)
{
}

/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 01/10|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::EvaluateTerminalBC(
    RedAirBloodScatraLine3* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, CORE::LINALG::SerialDenseVector& rhs,
    Teuchos::RCP<MAT::Material> material)
{
}


/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 01/10|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::CalcFlowRates(RedAirBloodScatraLine3* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    CORE::LINALG::SerialDenseVector& elevec1,  // a_volumenp,
    CORE::LINALG::SerialDenseVector& elevec2,  // a_volume_strain_np,
    std::vector<int>& lm, Teuchos::RCP<MAT::Material> material)

{
}

/*----------------------------------------------------------------------*
 |  Get the coupled the values on the coupling interface    ismail 07/10|
 |  of the 3D/reduced-D problem                                         |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::GetCoupledValues(
    RedAirBloodScatraLine3* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Teuchos::RCP<MAT::Material> material)
{
}

/*----------------------------------------------------------------------*
 |  solve the transport of O2 from air to blood             ismail 02/13|
 |                                                                      |
 | Example of use (1):                                                  |
 |--------------------                                                  |
 |        [RED_AIRWAY element]          [RED_ACINUS element]            |
 |             |                             |                          |
 |             |                             |                          |
 |  (node1)    V     (node2)        (node2)  V   (node1)                |
 |     o======>>>=======o              o============o                   |
 |     |       ^        |              |                                |
 |     |       |        |              |                                |
 |     |(flow direction)|              |                                |
 |     |                |              |                                |
 |     V                V              |                                |
 |     o=====           o===           |                                |
 |  (node1)  =====   (node3)==         |                                |
 |    or 3      ^ ==== or 1   ==       |                                |
 |              |     ======    ===    |                                |
 |              |           =====  ==  V                                |
 |              |                ======o                                |
 |              |                   (node2)                             |
 |              |                                                       |
 |              |                                                       |
 |    [RED_AIR_BLOOD_SCATRA_LINE3 element]                              |
 |                                                                      |
 | Example of use (2):                                                  |
 |--------------------                                                  |
 |        [RED_AIRWAY element]          [RED_ACINUS element]            |
 |             |                             |                          |
 |             |                             |                          |
 |  (node1)    V     (node2)        (node2)  V   (node1)                |
 |     o======>>>=======o              o============o                   |
 |     |       ^        |              |                                |
 |     |       |        |              |                                |
 |     |(flow direction)|              |                                |
 |     |                |              |                                |
 |     V                V              |                                |
 |     o================o==============o                                |
 |  (node1)     ^    (node2)        (node3)                             |
 |    or 3      |      or 1           or 1                              |
 |              |                                                       |
 |              |                                                       |
 |    [RED_AIR_BLOOD_SCATRA_LINE3 element]                              |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::RedAirBloodScatraLine3Impl<distype>::SolveBloodAirTransport(
    RedAirBloodScatraLine3* ele, CORE::LINALG::SerialDenseVector& dscatra,
    CORE::LINALG::SerialDenseVector& dvo2, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Teuchos::RCP<MAT::Material> material)
{
}
