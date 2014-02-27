/*!-----------------------------------------------------------------------------------------------*
\file xfem_fluidwizard.H

\brief class that provides the interface that bridges the cut libraries and fluid part of XFEM

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef XFEM_FLUIDWIZARD_H
#define XFEM_FLUIDWIZARD_H


#include <Teuchos_RCP.hpp>
#include "../drt_geometry/geo_utils.H"
#include "physics.H"


class Epetra_Vector;
class Epetra_IntVector;

namespace DRT
{
  class Discretization;
  class Element;
}

namespace GEO
{
  class CutWizard;

  namespace CUT
  {
    class ElementHandle;
    class SideHandle;
    class Side;
    class Node;
  }
}

namespace Teuchos
{
  class ParameterList;
}

namespace XFEM
{
  class FieldEnr;
  class ElementAnsatz;
  class FluidDofSet;

class FluidWizard
{
public:

  /*!
  \brief The interface that bridges the cut libraries and fluid part of XFEM
  */
  FluidWizard( DRT::Discretization & backdis, DRT::Discretization & cutterdis )
    : backdis_( backdis ),
      cutterdis_( cutterdis )
  {
  }

  /*!
  \brief Cut routine for the new XFEM framework (XFSI and XFLUIDFLUID)
  */
  void Cut(   bool include_inner,             //!< perform cut within the structure
              const Epetra_Vector & idispcol, //!< col vector holding interface displacements
              std::string VCellgausstype,     //!< Gauss point generation method for Volumecell
              std::string BCellgausstype,     //!< Gauss point generation method for Boundarycell
              bool parallel,                  //!< use parallel cut algorithms
              bool gmsh_output,               //!< print write gmsh output for cut
              bool positions = true,          //!< set inside and outside point, facet and volumecell positions);
              bool tetcellsonly = false,      //!< generate only tet cells
              bool screenoutput = true,       //!< print screen output
              bool cutinrefconf = false       //!< do not try to perform coordinate update
  );

  //! create a new fluid dofset and reserve a dofs for a maximal number of dofsets
  Teuchos::RCP<FluidDofSet> DofSet(int maxNumMyReservedDofs);

  //! get the cut wizard
  GEO::CutWizard& CutWizard();

  //! get the elementhandle created within the cut
  GEO::CUT::ElementHandle* GetElement( DRT::Element * ele );

  //! get the sidehandle created within the cut
  GEO::CUT::SideHandle* GetSide( std::vector<int> & nodeids );

  //! get the sidehandle created within the cut
  GEO::CUT::SideHandle* GetSideHandle( int sid );

  //! get the sidehandle for cut sides
  GEO::CUT::SideHandle * GetCutSide( int sid, int mi );

  //! get the node created within the cut
  GEO::CUT::Node * GetNode( int nid );

  //! set the nodes representing crack tip in FSI with crack structure simulations
  void setCrackTipNodes( std::map<int, LINALG::Matrix<3,1> > & tip ){ tip_nodes_ = tip; }

private:

  /// cut wizard
  Teuchos::RCP<GEO::CutWizard> cut_;

  /// background discretization
  DRT::Discretization & backdis_;

  /// cutter discretization
  DRT::Discretization & cutterdis_;

  /// nodes of the crack tip when simulating FSI with crack
  std::map<int, LINALG::Matrix<3,1> > tip_nodes_;


};

}

#endif
