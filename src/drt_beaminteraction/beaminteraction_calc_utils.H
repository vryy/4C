/*-----------------------------------------------------------*/
/*!
\file beaminteraction_calc_utils.H

\brief utils for biopolymer network business

\maintainer Jonas Eichinger, Maximilian Grill

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BEAMINTERACTION_CALC_UTILS_H_
#define BEAMINTERACTION_CALC_UTILS_H_

#include "Teuchos_RCP.hpp"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_mapextractor.H"


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

namespace BEAMINTERACTION
{
namespace UTILS
{

/// specific MultiMapExtractor to handle different types of element during beam interaction
class MapExtractor : public LINALG::MultiMapExtractor
{
public:

  enum
  {
    beam   = 0,
    sphere = 1,
    solid  = 2
  };

  MAP_EXTRACTOR_VECTOR_METHODS( Beam, beam )
  MAP_EXTRACTOR_VECTOR_METHODS( Sphere, sphere )
  MAP_EXTRACTOR_VECTOR_METHODS( Solid, solid )
};

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void PeriodicBoundaryConsistentDisVector(
    Teuchos::RCP<Epetra_Vector> dis,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& pbb,
    Teuchos::RCP<DRT::Discretization> const&        discret);

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void UpdateCellsPositionRandomly(
    Teuchos::RCP<Epetra_Vector> dis,
    Teuchos::RCP<DRT::Discretization> const& discret,
    int const step);

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
inline int CalculateNumberOfBeamElementsFromNumberOfNodesOnFilament(
    int const numnodes,
    int const numnodesperele)
{
  // from: nodesperfil = nodesperele + ( numele - 1 ) * ( nodesperele - 1 )
  return ( (numnodes - numnodesperele ) / ( numnodesperele - 1.0 ) ) + 1.0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::vector<int> Permutation( int const number );

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void GetCurrentElementDis(
    DRT::Discretization const& discret,
    DRT::Element const* ele,
    Teuchos::RCP<Epetra_Vector> const& ia_discolnp,
    std::vector<double>& eledisp);

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void GetCurrentUnshiftedElementDis(
    DRT::Discretization const& discret,
    DRT::Element const* ele,
    Teuchos::RCP<Epetra_Vector> const& ia_discolnp,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& pbb,
    std::vector<double>& eledisp);

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void GetPosAndTriadOfBindingSpot(
    DRT::Element* ele,
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
    const DRT::Discretization&                                     discret,
    std::vector<int> const&                                        elegid,
    std::vector< LINALG::SerialDenseVector > const&                elevec,
    std::vector< std::vector< LINALG::SerialDenseMatrix > > const& elemat,
    Teuchos::RCP<Epetra_FEVector>                                  fe_sysvec,
    Teuchos::RCP<LINALG::SparseMatrix>                             fe_sysmat);

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
void ExtractPosDofVecAbsoluteValues(
    DRT::Discretization const& discret,
    DRT::Element const* ele,
    Teuchos::RCP<Epetra_Vector> const& ia_discolnp,
    std::vector<double>& element_posdofvec_absolutevalues);

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ApplyBindingSpotForceToParentElements(
    DRT::Discretization const&                             discret,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> const&        pbb,
    const Teuchos::RCP<Epetra_Vector>                      disp_np_col,
    const Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> elepairptr,
    std::vector< LINALG::SerialDenseVector > const&        bspotforce,
    std::vector< LINALG::SerialDenseVector >&              eleforce);

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ApplyBindingSpotStiffToParentElements(
    DRT::Discretization const&                                     discret,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> const&                pbb,
    const Teuchos::RCP<Epetra_Vector>                              disp_np_col,
    const Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage>         elepairptr,
    std::vector< std::vector< LINALG::SerialDenseMatrix > > const& bspotstiff,
    std::vector< std::vector< LINALG::SerialDenseMatrix > >&       elestiff);

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ApplyBindingSpotForceStiffToParentElements(
    DRT::Discretization const&                                     discret,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> const&                pbb,
    const Teuchos::RCP<Epetra_Vector>                              disp_np_col,
    const Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage>         elepairptr,
    std::vector< LINALG::SerialDenseVector > const&                bspotforce,
    std::vector< std::vector< LINALG::SerialDenseMatrix > > const& bspotstiff,
    std::vector< LINALG::SerialDenseVector >&                      eleforce,
    std::vector< std::vector< LINALG::SerialDenseMatrix > >&       elestiff);

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void SetupEleTypeMapExtractor(
    Teuchos::RCP<DRT::Discretization> const& discret,
    Teuchos::RCP<LINALG::MultiMapExtractor>& eletypeextractor);

}
}

#endif
