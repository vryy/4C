/*!-----------------------------------------------------------------------------------------------*
\file geo_meshintersection.H

\brief class that provides to set up a mesh cut based on surface meshes

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef GEO_MESHINTERSECTION_H
#define GEO_MESHINTERSECTION_H

#include "../drt_geometry/geo_intersection.H"


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
  }

  /// contains the cut
  class CutWizardMesh : public CutWizard
  {
  public:

    /*!
    \brief Constructor
    */
    CutWizardMesh( DRT::Discretization & dis, int numcutmesh=1 );

    /*========================================================================*/
    //! @name Add functionality for elements
    /*========================================================================*/

    /*!
    \brief Add this cut element to the cut libraries
    */
    void AddCutSide( int mi, DRT::Element * ele, const Epetra_SerialDenseMatrix & xyze );

    /*!
    \brief Add this background mesh element to the cut libraries
    */
    void AddElement( DRT::Element * ele, const Epetra_SerialDenseMatrix & xyze );

    /*!
    \brief Build the bounding volume tree for the collision detection in the context of the selfcut
    */
    void BuildBVTree();

    /*!
    \brief Build the static search tree for the collision detection
    */
    void BuildStaticSearchTree();

    /*========================================================================*/
    //! @name Cut functionality routines
    /*========================================================================*/

    /*!
    \brief Cut routine for standard non-parallel framework which is used only for cut_test
    */
    void Cut(
        bool include_inner,
        INPAR::CUT::VCellGaussPts VCellgausstype,  //!< Gauss point generation method for Volumecell
        INPAR::CUT::BCellGaussPts BCellgausstype,  //!< Gauss point generation method for Boundarycell
        bool screenoutput=true
    );


    /*========================================================================*/
    //! @name Get routines
    /*========================================================================*/

    /*!
    \brief Get this side (not from cut meshes) from the cut libraries
     */
    GEO::CUT::SideHandle * GetCutSide( int sid, int mi );

    /*!
    \brief get the mesh intersection object
     */
    Teuchos::RCP<GEO::CUT::MeshIntersection> Mesh() {return meshintersection_; }


  protected:

    void FindNodePositions();

  private:


    /*========================================================================*/
    //! @name private class variables
    /*========================================================================*/

    Teuchos::RCP<GEO::CUT::MeshIntersection> meshintersection_;         ///< mesh intersection object
  };
}

#endif
