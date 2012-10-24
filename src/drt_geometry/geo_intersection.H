/*!-----------------------------------------------------------------------------------------------*
\file geo_intersection.H

\brief class that provides to set up a mesh cut based on a level set field or on further surface meshes

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef GEO_INTERSECTION_H
#define GEO_INTERSECTION_H


#include <Teuchos_RCP.hpp>


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
    class MeshIntersection;
    class SideHandle;
    class ElementHandle;
    class Side;
    class Node;
  }

  /// contains the cut
  class CutWizard
  {
  public:

    /*!
    \brief Constructor
    */
    CutWizard( DRT::Discretization & dis, bool include_inner, int numcutmesh=1 );

    /*!
    \brief Set the value of positions in the cut_option.H
    */
    void SetFindPositions( bool positions );

    /*!
    \brief Add this cut element to the cut libraries
    */
    void AddCutSide( int mi, DRT::Element * ele, const Epetra_SerialDenseMatrix & xyze );

    /*!
    \brief Add this background mesh element to the cut libraries
    */
    void AddElement( DRT::Element * ele );

    /*!
    \brief Get this side (not from cut meshes) from the cut libraries
     */
    GEO::CUT::Side * GetSide( std::vector<int>& nodeids );

    /*!
    \brief Get this side (not from cut meshes) from the cut libraries
     */
    GEO::CUT::SideHandle * GetSide( int sid );

    /*!
    \brief Get this side (not from cut meshes) from the cut libraries
     */
    GEO::CUT::SideHandle * GetCutSide( int sid, int mi );

    /*!
    \brief Get this element from the cut libraries
     */
    GEO::CUT::ElementHandle * GetElement( DRT::Element * ele );

    /*!
    \brief Get this node from the cut libraries
     */
    GEO::CUT::Node * GetNode( int nid );

    /*!
    \brief Cut routine for parallel framework in XFSI and XFLUIDFLUID
    */
    void CutParallel( bool include_inner, std::string VCellgausstype, std::string BCellgausstype );

    /*!
    \brief Routine for finding node positions and computing volumecell dofsets in a parallel way
    */
    void CutParallel_FindPositionDofSets(bool include_inner, bool communicate);

    /*!
    \brief Cut routine for standard non-parallel framework which is used only for cut_test
    */
    void Cut( bool include_inner, std::string VCellgausstype, std::string BCellgausstype );

//     Teuchos::RCP<Epetra_CrsGraph> MatrixGraph( const CutDofSet & dofset, const Epetra_Map & dbcmap );

    /*!
    \brief get the mesh intersection object
     */
    GEO::CUT::MeshIntersection & Mesh() { return *mesh_; }

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

  private:

    DRT::Discretization & dis_;                             ///< discretization
    int myrank_;                                            ///< my processor Id
    Teuchos::RCP<GEO::CUT::MeshIntersection> mesh_;         ///< mesh intersection object
    bool include_inner_;                                    ///< include the inner part
  };
}

#endif
