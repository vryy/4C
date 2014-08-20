/*!-----------------------------------------------------------------------------------------------*
\file geo_intersection.H

\brief class that provides the common functionality for a mesh cut based on a level set field or on surface meshes

<pre>
Maintainer: Benedikt Schott and Magnus Winter
            schott@lnm.mw.tum.de, winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef GEO_INTERSECTION_H
#define GEO_INTERSECTION_H

#include <Teuchos_RCP.hpp>

#include "../drt_inpar/inpar_cut.H"

class Epetra_Map;
class Epetra_SerialDenseMatrix;

namespace DRT
{
  class Discretization;
  class Element;
}

namespace GEO
{
  namespace CUT
  {
    class ParentIntersection;
    class SideHandle;
    class ElementHandle;
    class Node;
  }

  /// contains the cut, and shared functionality between the level set and mesh cut.
  class CutWizard
  {
  public:

    /*!
    \brief Constructor
    */
    CutWizard( const DRT::Discretization & dis);

    /*========================================================================*/
    //! @name Set routines
    /*========================================================================*/

    /*!
    \brief Set the value of positions in the cut_option.H
    */
    void SetFindPositions( bool positions );

    /*========================================================================*/
    //! @name Get routines
    /*========================================================================*/

    /*!
    \brief Get this side (not from cut meshes) from the cut libraries
     */
    GEO::CUT::SideHandle * GetSide( std::vector<int>& nodeids );

    /*!
    \brief Get this side (not from cut meshes) from the cut libraries
     */
    GEO::CUT::SideHandle * GetSide( int sid );

    /*!
    \brief Get this element from the cut libraries
     */
    GEO::CUT::ElementHandle * GetElement( DRT::Element * ele );

    /*!
    \brief Get this node from the cut libraries
     */
    GEO::CUT::Node * GetNode( int nid );

//    virtual Teuchos::RCP<Epetra_CrsGraph> MatrixGraph( const CutDofSet & dofset, const Epetra_Map & dbcmap );

    /*========================================================================*/
    //! @name Output routines
    /*========================================================================*/

    /*!
    \brief Print the number of volumecells and boundarycells generated over the whole mesh during the cut
     */
    void PrintCellStats();

    /*!
    \brief Write the DOF details of the nodes
     */
    void DumpGmshNumDOFSets( bool include_inner );

    /*!
    \brief Write volumecell output in GMSH format throughout the domain
     */
    void DumpGmshVolumeCells( bool include_inner );

    /*!
    \brief Write the integrationcells and boundarycells in GMSH format throughout the domain
     */
    void DumpGmshIntegrationCells();

  protected:

    /*========================================================================*/
    //! @name protected class variables
    /*========================================================================*/

    const DRT::Discretization & dis_;                            ///< discretization
    int myrank_;                                                 ///< my processor Id
    Teuchos::RCP<GEO::CUT::ParentIntersection> mesh_;            ///< parent intersection object

  };
}

#endif
