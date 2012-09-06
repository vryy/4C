/*!-------------------------------------------------------------------------*\
 * \file meshfree_scatra_cell_evaluate.cpp
 *
 * \brief evaluation of a meshfree scatra cell
 *
 * <pre>
 * Maintainer: Keijo Nissen (nis)
 *             nissen@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15253
 * </pre>
 *
\*--------------------------------------------------------------------------*/

#include "meshfree_scatra_cell.H"       // eh klar...
#include "../drt_inpar/inpar_scatra.H"  // enums ScaTraType and FluxType
#include "meshfree_scatra_impl.H"       // DRT::ELEMENTS::MeshfreeScaTraImplInterface

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             nis Mar12 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::MeshfreeTransport::Evaluate(
    ParameterList&            params,
    DRT::Discretization&      discretization,
    vector<int>&              lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  // the type of scalar transport problem has to be provided for all actions!
  const INPAR::SCATRA::ScaTraType scatratype = DRT::INPUT::get<INPAR::SCATRA::ScaTraType>(params, "scatratype");
  if (scatratype == INPAR::SCATRA::scatratype_undefined)
    dserror("Element parameter SCATRATYPE has not been set!");

  switch (distype_) {
    /*----------------------------------------------------------------------*
       Meshfree discretisation integrated on a backgroundmesh:
       ========================================================
       all physics-related stuff is included in the implementation class that
       can be used in principle inside any element. Elements themself are only
       used for non-overlapping domain decomposition and the creation of Gauss
       points. Since no information like basis function values can be stored
       in an impl-class approach, this implementation serves the purposes of
       validation. For efficient computations, non-meshbased schemes should be
       implemented and used.
     *----------------------------------------------------------------------*/

  case hex8:
  case tet4:
  case quad4:
  case tri3:
  case line2:

    return DRT::ELEMENTS::MeshfreeScaTraImplInterface::Impl(this,scatratype)->Evaluate(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elemat2,
        elevec1,
        elevec2,
        elevec3
        );
    break;

    /*----------------------------------------------------------------------*
       Meshfree discretisation with meshfree integration:
       ========================================================
       to be implemented
     *----------------------------------------------------------------------*/

  default:
    dserror("No element evaluate implemented for this distype.");
    return(-1);
  }
}

/*----------------------------------------------------------------------*
 |  do nothing (public)                                       nis Mar12 |
 |                                                                      |
 |  Not decided on an implementation of Neumann BC for meshfree schemes |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::MeshfreeTransport::EvaluateNeumann(ParameterList& params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    vector<int>&              lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}
