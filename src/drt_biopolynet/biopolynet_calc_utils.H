/*-----------------------------------------------------------*/
/*!
\file biopolynet_calc_utils.H

\brief utils for biopolymer network business

\maintainer Jonas Eichinger, Maximilian Grill

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BIOPOLYNET_CALC_UTILS_H_
#define BIOPOLYNET_CALC_UTILS_H_

#include "Teuchos_RCP.hpp"
#include "../linalg/linalg_fixedsizematrix.H"

// forward declaration
class Epetra_FEVector;
class Epetra_Vector;
namespace LINALG {
  class SerialDenseVector;
  class SerialDenseMatrix;
  class SparseMatrix;
}
namespace DRT {
  class Element;
  class Discretization;
}
namespace GEO
{
namespace MESHFREE
{
  class BoundingBox;
}
}
namespace BEAMINTERACTION {
class CrosslinkingParams;
class BeamToBeamLinkage;
}

namespace BIOPOLYNET
{
namespace UTILS
{

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void PeriodicBoundaryConsistentDisVector(
    Teuchos::RCP<Epetra_Vector> dis,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& pbb,
    Teuchos::RCP<DRT::Discretization> const&        discret);

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::vector<int> Permutation(const int& number);

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void GetCurrentElementDis(
    DRT::Discretization const& discret,
    DRT::Element const* ele,
    Teuchos::RCP<Epetra_Vector> const& ia_discolnp,
    std::vector<double>& eledisp);

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void GetPosAndTriadOfBindingSpot(
    DRT::Element* ele,
    Teuchos::RCP<Epetra_Vector> const& ia_discolnp,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& pbb,
    int const locbspotnum,
    LINALG::Matrix<3,1>& bspotpos,
    LINALG::Matrix<3,3>& bspottriad,
    std::vector<double>& eledisp);

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void GetPosAndTriadOfBindingSpot(
    DRT::Discretization const& discret,
    DRT::Element* ele,
    Teuchos::RCP<Epetra_Vector> const& ia_discolnp,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& pbb,
    int const locbspotnum,
    LINALG::Matrix<3,1>& bspotpos,
    LINALG::Matrix<3,3>& bspottriad);

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void FEAssembleEleForceStiffIntoSystemVectorMatrix(
    const DRT::Discretization&         discret,
    std::vector<int> const&            elegid,
    std::vector< LINALG::SerialDenseVector > const& elevec,
    std::vector< std::vector< LINALG::SerialDenseMatrix > > const& elemat,
    Teuchos::RCP<Epetra_FEVector>      fe_sysvec,
    Teuchos::RCP<LINALG::SparseMatrix> fe_sysmat);

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void AssembleCenterlineDofForceStiffIntoElementForceStiff(
    DRT::Discretization const&                                     discret,
    std::vector<int> const&                                        elegid,
    std::vector< LINALG::SerialDenseVector > const&                eleforce_centerlineDOFs,
    std::vector< std::vector< LINALG::SerialDenseMatrix > > const& elestiff_centerlineDOFs,
    std::vector< LINALG::SerialDenseVector >*                      eleforce,
    std::vector< std::vector< LINALG::SerialDenseMatrix > >*       elestiff);

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ApplyBpotForceToParentElements(
    DRT::Discretization const&                             discret,
    const Teuchos::RCP<Epetra_Vector>                      disp_np_col,
    const Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> elepairptr,
    std::vector< LINALG::SerialDenseVector > const&        bspotforce,
    std::vector< LINALG::SerialDenseVector >&              eleforce);

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ApplyBpotStiffToParentElements(
    DRT::Discretization const&                                     discret,
    const Teuchos::RCP<Epetra_Vector>                              disp_np_col,
    const Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage>         elepairptr,
    std::vector< std::vector< LINALG::SerialDenseMatrix > > const& bspotstiff,
    std::vector< std::vector< LINALG::SerialDenseMatrix > >&       elestiff);

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ApplyBpotForceStiffToParentElements(
    DRT::Discretization const&                                     discret,
    const Teuchos::RCP<Epetra_Vector>                              disp_np_col,
    const Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage>         elepairptr,
    std::vector< LINALG::SerialDenseVector > const&                bspotforce,
    std::vector< std::vector< LINALG::SerialDenseMatrix > > const& bspotstiff,
    std::vector< LINALG::SerialDenseVector >&                      eleforce,
    std::vector< std::vector< LINALG::SerialDenseMatrix > >&       elestiff);

}
}

#endif
