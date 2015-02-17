/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_std.cpp

\brief evaluation of scatra boundary terms at integration points

 <pre>
   Maintainer: Moritz Thon & Andre Hemmler
               thon@mhpc.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-10364
 </pre>
 */
/*----------------------------------------------------------------------*/

#include "scatra_ele_boundary_calc_poro.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_action.H"
#include "scatra_ele.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"
#include "../drt_geometry/position_array.H"

/*----------------------------------------------------------------------*
 |  Singleton access method                               hemmler 07/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype> * DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype>::Instance(const int numdofpernode, const int numscal, bool create )
{
  static ScaTraEleBoundaryCalcPoro<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleBoundaryCalcPoro<distype>(numdofpernode, numscal);
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}

/*----------------------------------------------------------------------*
 |  Clean up                                              hemmler 07/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
    Instance( 0, 0, false );
}


/*----------------------------------------------------------------------*
 |  Private constructor                                   hemmler 07/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype>::ScaTraEleBoundaryCalcPoro(const int numdofpernode, const int numscal)
  : DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::ScaTraEleBoundaryCalc(numdofpernode,numscal)
{
  return;
}

/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<distype>::EvaluateAction(
    DRT::ELEMENTS::TransportBoundary*   ele,
    Teuchos::ParameterList&             params,
    DRT::Discretization&                discretization,
    SCATRA::BoundaryAction              action,
    std::vector<int>&                   lm,
    Epetra_SerialDenseMatrix&           elemat1_epetra,
    Epetra_SerialDenseMatrix&           elemat2_epetra,
    Epetra_SerialDenseVector&           elevec1_epetra,
    Epetra_SerialDenseVector&           elevec2_epetra,
    Epetra_SerialDenseVector&           elevec3_epetra
    )
{
  // determine and evaluate action
  switch(action)
  {
  case SCATRA::bd_calc_fps3i_surface_permeability:
  case SCATRA::bd_calc_fs3i_surface_permeability:
  case SCATRA::bd_calc_Neumann:
  {
    my::EvaluateAction(
        ele,
        params,
        discretization,
        action,
        lm,
        elemat1_epetra,
        elemat2_epetra,
        elevec1_epetra,
        elevec2_epetra,
        elevec3_epetra
        );
    break;
  }
  default:
  {
    dserror("Invalid action parameter!");
    break;
  }
  } // switch(action)

  return 0;
}


// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcPoro<DRT::Element::nurbs9>;
