/*!-----------------------------------------------------------------------------------------------*
\file cut_cutwizard.H

\brief class that provides the common functionality for a mesh cut based on a level set field or on surface meshes

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef CUT_CUTWIZARD_H
#define CUT_CUTWIZARD_H

#include <Teuchos_RCP.hpp>

#include "../drt_inpar/inpar_cut.H"

class Epetra_Vector;
class Epetra_SerialDenseMatrix;

namespace DRT
{
  class Discretization;
  class DiscretizationXFEM;
  class Element;
}

namespace GEO
{
  namespace CUT
  {
    class CombIntersection;
    class ElementHandle;
    class Node;
    class SideHandle;
  }

  /// contains the cut, and shared functionality between the level set and mesh cut.
  class CutWizardNEW
  {
  public:


    /*========================================================================*/
    //! @name Constructor and Destructor
    /*========================================================================*/

    /*!
    \brief Constructor
    */
    CutWizardNEW( Teuchos::RCP<DRT::DiscretizationXFEM> xdis, Teuchos::RCP<DRT::Discretization> cutterdis );

    /*!
    \brief Destructor
    */
    virtual ~CutWizardNEW(){};

    //@}


    /*========================================================================*/
    //! @name Setters
    /*========================================================================*/

    //! set options and flags used during the cut
    void SetOptions(
        INPAR::CUT::VCellGaussPts VCellgausstype,   //!< Gauss point generation method for Volumecell
        INPAR::CUT::BCellGaussPts BCellgausstype,   //!< Gauss point generation method for Boundarycell
        bool gmsh_output,                           //!< print write gmsh output for cut
        bool positions,                             //!< set inside and outside point, facet and volumecell positions
        bool tetcellsonly,                          //!< generate only tet cells
        bool screenoutput                           //!< print screen output
    );

    //! set displacement and level-set vectors used during the cut
    void SetState(
        Teuchos::RCP<const Epetra_Vector> back_disp_col,      //!< col vector holding background ALE displacements for backdis
        Teuchos::RCP<const Epetra_Vector> cutter_disp_col,    //!< col vector holding interface displacements for cutterdis
        Teuchos::RCP<const Epetra_Vector> back_levelset_col   //!< col vector holding nodal level-set values based on backdis
    );

    //! set crack tip nodes for crack applications
    void setCrackTipNodes( std::map<int, LINALG::Matrix<3,1> > & tip );

    //@}


    /*========================================================================*/
    //! @name main Cut call
    /*========================================================================*/

    void Cut(
        bool include_inner  //!< perform cut in the interior of the cutting mesh
    );


    /*========================================================================*/
    //! @name Getters
    /*========================================================================*/

    //! Get this side (not from cut meshes) (faces of background elements) from the cut libraries
    GEO::CUT::SideHandle * GetSide( std::vector<int>& nodeids );

    //! Get this side (not from cut meshes) from the cut libraries
    GEO::CUT::SideHandle * GetSide( int sid );

    //! Get this element from the cut libraries
    GEO::CUT::ElementHandle * GetElement( DRT::Element * ele );

    //! Get this node from the cut libraries
    GEO::CUT::Node * GetNode( int nid );

    //! Get the sidehandle for cutting sides
    GEO::CUT::SideHandle * GetMeshCuttingSide( int sid, int mi );

    //! is there a level-set side with the given sid?
    bool HasLSCuttingSide( int sid );


  private:



    /*========================================================================*/
    //! @name Add functionality for elements and cutting sides
    /*========================================================================*/

    //! add all cutting sides (mesh and level-set sides)
    void AddCuttingSides();

    //! add level-set cutting side
    void AddLSCuttingSide();

    //! add all cutting sides from the cut-discretization
    void AddMeshCuttingSide();

    //! add elements from the background discretization
    void AddBackgroundElements();

    //! Add this cuttin side element with given global coordiantes to the intersection class
    void AddMeshCuttingSide( int mi, DRT::Element * ele, const Epetra_SerialDenseMatrix & xyze );

    //! Add this background mesh element to the intersection class
    void AddElement( DRT::Element * ele, const Epetra_SerialDenseMatrix & xyze, double* myphinp = NULL, bool ls_only_plus_domain = false );

    //@}


    /*========================================================================*/
    //! @name Major steps to prepare the cut, to perform it and to do postprocessing
    /*========================================================================*/

    //! prepare the cut, add background elements and cutting sides
    void Prepare();

    //! perform the actual cut, the intersection
    void Run_Cut(
        bool include_inner  //!< perform cut in the interior of the cutting mesh
    );

    //! routine for finding node positions and computing volume-cell dofsets in a parallel way
    void FindPositionDofSets(bool include_inner, bool screenoutput);

    //! write statistics and output to screen and files
    void Output(bool include_inner);

    //@}


    /*========================================================================*/
    //! @name Output routines
    /*========================================================================*/

    //! Print the number of volumecells and boundarycells generated over the whole mesh during the cut
    void PrintCellStats();

    //! Write the DOF details of the nodes
    void DumpGmshNumDOFSets( bool include_inner );

    //! Write volumecell output in GMSH format throughout the domain
    void DumpGmshVolumeCells( bool include_inner );

    //! Write the integrationcells and boundarycells in GMSH format throughout the domain
    void DumpGmshIntegrationCells();

    //@}


    //---------------------------------discretizations----------------------------

    //! @name discretizations
    Teuchos::RCP<DRT::DiscretizationXFEM> backdis_;              ///< background discretization
    Teuchos::RCP<DRT::Discretization> cutterdis_;                ///< cutter discretization

    int myrank_;                                                 ///< my processor Id
    //@}


    //---------------------------------main intersection class----------------------------
    //! @name main intersection class and flags
    Teuchos::RCP<GEO::CUT::CombIntersection> intersection_;      ///< combined intersection object which handles cutting mesh sides and a level-set side

    bool do_mesh_intersection_;                                  ///< flag to perform intersection with mesh sides
    bool do_levelset_intersection_;                              ///< flag to perform intersection with a level-set side
    //@}

    //---------------------------------state vectors ----------------------------

    //! @name state vectors holding displacements and level-set values
    Teuchos::RCP<const Epetra_Vector> back_disp_col_;      ///< col vector holding background ALE displacements for backdis
    Teuchos::RCP<const Epetra_Vector> cutter_disp_col_;    ///< col vector holding interface displacements for cutterdis
    Teuchos::RCP<const Epetra_Vector> back_levelset_col_;  ///< col vector holding nodal level-set values based on backdis
    //@}


    //---------------------------------Options ----------------------------

    //! @name Options
    INPAR::CUT::VCellGaussPts VCellgausstype_;             ///< integration type for volume-cells
    INPAR::CUT::BCellGaussPts BCellgausstype_;             ///< integration type for boundary-cells
    bool gmsh_output_;                                     ///< write gmsh output?
    bool tetcellsonly_;                                    ///< enforce to create tetrahedral integration cells exclusively
    bool screenoutput_;                                    ///< write output to screen
    bool lsv_only_plus_domain_;                            ///< consider only plus domain of level-set field as physical field
    //@}

    //---------------------------------Crack specific stuff ----------------------------

    //! @name Options
    std::map<int, LINALG::Matrix<3,1> > tip_nodes_;        ///< nodes of the crack tip when simulating FSI with crack
    bool is_crack_;                                        ///< are crack tipnodes involved?
    //@}

  };

}

#endif
