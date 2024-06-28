/*----------------------------------------------------------------------*/
/*! \file

\brief class that provides the common functionality for a mesh cut based on a level set field or on
surface meshes

\level 3
*------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_CUTWIZARD_HPP
#define FOUR_C_CUT_CUTWIZARD_HPP

#include "4C_config.hpp"

#include "4C_fem_general_node.hpp"
#include "4C_inpar_cut.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Elements
{
  class Element;
}

namespace Core::LinAlg
{
  class SerialDenseMatrix;
}

namespace XFEM
{
  class ConditionManager;
}

namespace Core::Geo
{
  namespace Cut
  {
    class CombIntersection;
    class ElementHandle;
    class Node;
    class SideHandle;
  }  // namespace Cut

  /// contains the cut, and shared functionality between the level set and mesh cut.
  class CutWizard
  {
   public:
    /*------------------------------------------------------------------------*/
    /*! \brief Container class for the background mesh object
     *
     *  \author hiermeier \date 01/17 */
    class BackMesh
    {
     public:
      /// constructor
      explicit BackMesh(
          const Teuchos::RCP<Core::FE::Discretization>& backdis, Core::Geo::CutWizard* wizard)
          : wizard_(wizard),
            back_discret_(backdis),
            back_disp_col_(Teuchos::null),
            back_levelset_col_(Teuchos::null)
      {
        if (backdis.is_null()) FOUR_C_THROW("null pointer to background dis, invalid!");
      }

      virtual ~BackMesh() = default;

      void init(const Teuchos::RCP<const Epetra_Vector>& back_disp_col,
          const Teuchos::RCP<const Epetra_Vector>& back_levelset_col);

      const Teuchos::RCP<Core::FE::Discretization>& GetPtr() { return back_discret_; }

      Core::FE::Discretization& get() { return *back_discret_; }

      const Core::FE::Discretization& get() const { return *back_discret_; }

      virtual int NumMyColElements() const;

      virtual const Core::Elements::Element* lColElement(int lid) const;

      inline bool IsBackDisp() const { return (not back_disp_col_.is_null()); }

      const Epetra_Vector& BackDispCol() const
      {
        if (not IsBackDisp())
          FOUR_C_THROW("The background displacement was not initialized correctly!");

        return *back_disp_col_;
      }

      inline bool IsLevelSet() const { return (not back_levelset_col_.is_null()); }

      const Epetra_Vector& BackLevelSetCol() const
      {
        if (not IsLevelSet())
          FOUR_C_THROW("No level-set values set for the background discretization!");

        return *back_levelset_col_;
      }


     protected:
      Core::Geo::CutWizard* wizard_;

     private:
      /// background discretization
      Teuchos::RCP<Core::FE::Discretization> back_discret_;

      /// col vector holding background ALE displacements for backdis
      Teuchos::RCP<const Epetra_Vector> back_disp_col_;

      /// col vector holding nodal level-set values based on backdis
      Teuchos::RCP<const Epetra_Vector> back_levelset_col_;
    };

    /*------------------------------------------------------------------------*/
    /*!
     * \brief Container class for a certain cutting mesh objects
     */
    class CutterMesh
    {
     public:
      //! ctor
      CutterMesh(Teuchos::RCP<Core::FE::Discretization> cutterdis,
          Teuchos::RCP<const Epetra_Vector> cutter_disp_col, const int start_ele_gid)
          : cutterdis_(cutterdis), cutter_disp_col_(cutter_disp_col), start_ele_gid_(start_ele_gid)
      {
      }

      //---------------------------------discretization-----------------------------

      //! @name cutter discretization
      Teuchos::RCP<Core::FE::Discretization> cutterdis_;  ///< cutter discretization
      //@}

      //---------------------------------state vectors ----------------------------

      //! @name state vectors holding displacements
      Teuchos::RCP<const Epetra_Vector>
          cutter_disp_col_;  ///< col vector holding interface displacements for cutterdis
      //@}

      //!
      int start_ele_gid_;
    };

    /*========================================================================*/
    //! @name Constructor and Destructor
    /*========================================================================*/

    /**
     * \brief Constructor.
     *
     * Create CutWizard object with the given background discretization. The optional
     * function @p global_dof_indices can be used to retrieve the global dof indices from
     * the background discretization.
     */
    CutWizard(const Teuchos::RCP<Core::FE::Discretization>& backdis,
        std::function<void(const Core::Nodes::Node& node, std::vector<int>& lm)>
            global_dof_indices = nullptr);


    /*!
    \brief Destructor
    */
    virtual ~CutWizard() = default;

    //@}

    /*========================================================================*/
    //! @name Setters
    /*========================================================================*/

    //! set options and flags used during the cut
    void SetOptions(Inpar::Cut::NodalDofSetStrategy
                        nodal_dofset_strategy,     //!< strategy for nodal dofset management
        Inpar::Cut::VCellGaussPts VCellgausstype,  //!< Gauss point generation method for Volumecell
        Inpar::Cut::BCellGaussPts
            BCellgausstype,  //!< Gauss point generation method for Boundarycell
        bool gmsh_output,    //!< print write gmsh output for cut
        bool positions,      //!< set inside and outside point, facet and volumecell positions
        bool tetcellsonly,   //!< generate only tet cells
        bool screenoutput    //!< print screen output
    );

    virtual void SetBackgroundState(
        Teuchos::RCP<const Epetra_Vector>
            back_disp_col,  //!< col vector holding background ALE displacements for backdis
        Teuchos::RCP<const Epetra_Vector>
            back_levelset_col,  //!< col vector holding nodal level-set values based on backdis
        int level_set_sid       //!< global id for level-set side
    );

    void AddCutterState(const int mc_idx, Teuchos::RCP<Core::FE::Discretization> cutter_dis,
        Teuchos::RCP<const Epetra_Vector> cutter_disp_col);

    void AddCutterState(const int mc_idx, Teuchos::RCP<Core::FE::Discretization> cutter_dis,
        Teuchos::RCP<const Epetra_Vector> cutter_disp_col, const int start_ele_gid);

    // Find marked background-boundary sides.
    //  Extract these sides and create boundary cell for these!
    void set_marked_condition_sides(
        // const int mc_idx,
        Teuchos::RCP<Core::FE::Discretization> cutter_dis,
        // Teuchos::RCP<const Epetra_Vector> cutter_disp_col,
        const int start_ele_gid);

    //@}

    /*========================================================================*/
    //! @name main Cut call
    /*========================================================================*/

    //! prepare the cut, add background elements and cutting sides
    void Prepare();

    void Cut(bool include_inner  //!< perform cut in the interior of the cutting mesh
    );

    /*========================================================================*/
    //! @name Accessors
    /*========================================================================*/

    //! Get this side (not from cut meshes) (faces of background elements) from the cut libraries
    Core::Geo::Cut::SideHandle* get_side(std::vector<int>& nodeids);

    //! Get this side (not from cut meshes) from the cut libraries
    Core::Geo::Cut::SideHandle* get_side(int sid);

    //! Get this side from cut meshes from the cut libraries
    Core::Geo::Cut::SideHandle* GetCutSide(int sid);

    //! Get this element from the cut libraries by element id
    Core::Geo::Cut::ElementHandle* GetElement(const int eleid) const;

    //! Get this element from the cut libraries by element pointer
    Core::Geo::Cut::ElementHandle* GetElement(const Core::Elements::Element* ele) const;

    //! Get this node from the cut libraries
    Core::Geo::Cut::Node* GetNode(int nid);

    //! Get the sidehandle for cutting sides
    Core::Geo::Cut::SideHandle* GetMeshCuttingSide(int sid, int mi);

    //! is there a level-set side with the given sid?
    bool HasLSCuttingSide(int sid);

    //! update the coordinates of the cut boundary cells
    void update_boundary_cell_coords(Teuchos::RCP<Core::FE::Discretization> cutterdis,
        Teuchos::RCP<const Epetra_Vector> cutter_disp_col, const int start_ele_gid);

    //! Cubaturedegree for creating of integrationpoints on boundarycells
    int get_bc_cubaturedegree() const;

   protected:
    /** \brief hidden constructor for derived classes only
     *
     *  \author hiermeier \date 01/17 */
    CutWizard(const Epetra_Comm& comm);

    Teuchos::RCP<BackMesh>& back_mesh_ptr() { return back_mesh_; }

    Teuchos::RCP<const BackMesh> back_mesh_ptr() const { return back_mesh_.getConst(); }

    virtual void get_physical_nodal_coordinates(
        const Core::Elements::Element* element, Core::LinAlg::SerialDenseMatrix& xyze) const;

    Core::Geo::Cut::CombIntersection& intersection()
    {
      if (intersection_.is_null()) FOUR_C_THROW("nullptr pointer!");

      return *intersection_;
    }

   private:
    /*========================================================================*/
    //! @name Add functionality for elements and cutting sides
    /*========================================================================*/

    //! add all cutting sides (mesh and level-set sides)
    void add_cutting_sides();

    //! add level-set cutting side
    void add_ls_cutting_side();

    //! add all cutting sides from the cut-discretization
    void add_mesh_cutting_side();

    //! add elements from the background discretization
    void add_background_elements();

    //! Add all cutting side elements of given cutter discretization with given displacement field
    //! to the intersection class
    void add_mesh_cutting_side(Teuchos::RCP<Core::FE::Discretization> cutterdis,
        Teuchos::RCP<const Epetra_Vector> cutter_disp_col,
        const int start_ele_gid = 0  ///< global start index for element id numbering
    );

    //! Add this cutting side element with given global coordinates to the intersection class
    void add_mesh_cutting_side(int mi, Core::Elements::Element* ele,
        const Core::LinAlg::SerialDenseMatrix& xyze, const int start_ele_gid);

    //! Add this background mesh element to the intersection class
    void add_element(const Core::Elements::Element* ele,
        const Core::LinAlg::SerialDenseMatrix& xyze, double* myphinp = nullptr,
        bool lsv_only_plus_domain = false);

    //@}


    /*========================================================================*/
    //! @name Major steps to prepare the cut, to perform it and to do postprocessing
    /*========================================================================*/

    //! perform the actual cut, the intersection
    void run_cut(bool include_inner  //!< perform cut in the interior of the cutting mesh
    );

    //! routine for finding node positions and computing volume-cell dofsets in a parallel way
    void find_position_dof_sets(bool include_inner);

    //! write statistics and output to screen and files
    void output(bool include_inner);

    //! Check that cut is initialized correctly
    bool safety_checks(bool is_prepare_cut_call);

    //@}

    /*========================================================================*/
    //! @name Output routines
    /*========================================================================*/

    /*! Print the number of volumecells and boundarycells generated over the
     *  whole mesh during the cut */
    void print_cell_stats();

    //! Write the DOF details of the nodes
    void dump_gmsh_num_dof_sets(bool include_inner);

    //! Write volumecell output in GMSH format throughout the domain
    void dump_gmsh_volume_cells(bool include_inner);

    //! Write the integrationcells and boundarycells in GMSH format throughout the domain
    void dump_gmsh_integration_cells();

    //@}

    //---------------------------------discretizations----------------------------

    //! @name meshes
    Teuchos::RCP<BackMesh> back_mesh_;
    std::function<void(const Core::Nodes::Node& node, std::vector<int>& lm)> global_dof_indices_;
    std::map<int, Teuchos::RCP<CutterMesh>> cutter_meshes_;
    const Epetra_Comm& comm_;
    int myrank_;  ///< my processor Id
    //@}

    //---------------------------------main intersection class----------------------------
    //! @name main intersection class and flags
    Teuchos::RCP<Core::Geo::Cut::CombIntersection>
        intersection_;  ///< combined intersection object which handles cutting mesh sides and a
                        ///< level-set side

    bool do_mesh_intersection_;      ///< flag to perform intersection with mesh sides
    bool do_levelset_intersection_;  ///< flag to perform intersection with a level-set side
    //@}

    //---------------------------------state vectors ----------------------------

    //! @name state vectors holding displacements and level-set values
    int level_set_sid_;
    //@}

    //---------------------------------Options ----------------------------

    //! @name Options
    Inpar::Cut::VCellGaussPts v_cellgausstype_;  ///< integration type for volume-cells
    Inpar::Cut::BCellGaussPts b_cellgausstype_;  ///< integration type for boundary-cells
    bool gmsh_output_;                           ///< write gmsh output?
    bool tetcellsonly_;          ///< enforce to create tetrahedral integration cells exclusively
    bool screenoutput_;          ///< write output to screen
    bool lsv_only_plus_domain_;  ///< consider only plus domain of level-set field as physical field
    //@}

    //--------------------------------- Initialization flags ----------------------------

    //! @name Flags whether wizard is initialized correctly
    bool is_set_options_;
    bool is_cut_prepare_performed_;
    //@}

  };  // class CutWizard
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
