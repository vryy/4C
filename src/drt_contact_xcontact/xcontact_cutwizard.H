/*----------------------------------------------------------------------*/
/*! \file
\brief Cutwizard of the xcontact algorithm ( level-set cut of the structural
       XFEM discretization )


\level 3
*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_XCONTACT_XCONTACT_CUTWIZARD_H_
#define SRC_DRT_CONTACT_XCONTACT_XCONTACT_CUTWIZARD_H_

// DEBUG FLAGS
//#define DEBUG_XCONTACT_INTERSECTION

#include "../headers/pairedvector.H"
#include "../drt_lib/drt_element.H"
#include "../drt_cut/cut_utils.H"

#include "../drt_cut/cut_cutwizard.H"

// forward declarations
class Epetra_Vector;
class Epetra_Map;
class Epetra_Comm;

namespace LINALG
{
  class SerialDenseMatrix;
}  // namespace LINALG

namespace DRT
{
  class Discretization;
  class DiscretizationInterface;
}  // namespace DRT

namespace GEO
{
  class BoundaryIntCell;
  typedef std::vector<Teuchos::RCP<GEO::BoundaryIntCell>> BoundaryIntCellPtrs;
  namespace CUT
  {
    class ElementHandle;
    class LevelSetIntersection;
    class CombIntersection;
  }  // namespace CUT
}  // namespace GEO

namespace XCONTACT
{
  class CutWizard : public GEO::CutWizard
  {
   protected:
    typedef GEN::pairedvector<const DRT::Element *, std::vector<const DRT::FaceElement *>>
        EleFacePair;

    /*! \brief Container class for the background mesh object
     *
     *  \author hiermeier \date 01/17 */
    class BackMesh : public GEO::CutWizard::BackMesh
    {
     public:
      /// constructor
      BackMesh(const Teuchos::RCP<DRT::Discretization> &backdis, XCONTACT::CutWizard *wizard)
          : GEO::CutWizard::BackMesh(backdis, wizard)
      {
        /* left blank */
      }

      /// destructor
      virtual ~BackMesh(){};

      /// returns the number of considered column elements
      virtual int NumMyColElements() const;

      /// returns a pointer to the col element with the given local ID
      virtual const DRT::Element *lColElement(int lid) const;

     private:
      /// call-back to the wrapping cut wizard class (read and write)
      XCONTACT::CutWizard *Wizard() { return static_cast<XCONTACT::CutWizard *>(wizard_); }

      /// call-back to the wrapping cut wizard class (read-only)
      const XCONTACT::CutWizard *Wizard() const
      {
        return static_cast<XCONTACT::CutWizard *>(wizard_);
      }
    };  // class BackMesh

   public:
    /// constructor
    CutWizard(const Teuchos::RCP<DRT::DiscretizationInterface> &backdis);

    /// destructor
    virtual ~CutWizard(){};

    void SetBackgroundState(Teuchos::RCP<const Epetra_Vector> back_disp_col,
        Teuchos::RCP<const Epetra_Vector> back_levelset_row, int level_set_sid);

    /** \brief Set cut wizard options for the xcontact formulation
     *
     *  \param gen_bcells_position   (in) : generate boundary cells of facets with the given
     * position \param nodal_dofset_strategy (in) : strategy for nodal dofset management \param
     * VCellgausstype        (in) : Gauss point generation method for Volumecell \param
     * BCellgausstype        (in) : Gauss point generation method for Boundarycell \param
     * gmsh_output           (in) : print write gmsh output for cut \param positions (in) : set
     * inside and outside point, facet and volumecell positions \param tetcellsonly          (in) :
     * generate only tet cells \param screenoutput          (in) : print screen output
     *
     *  \author hiermeier \date 01/17 */
    void SetOptions(INPAR::CUT::BoundaryCellPosition gen_bcells_position,
        INPAR::CUT::NodalDofSetStrategy nodal_dofset_strategy,
        INPAR::CUT::VCellGaussPts VCellgausstype, INPAR::CUT::BCellGaussPts BCellgausstype,
        bool gmsh_output, bool positions, bool tetcellsonly, bool screenoutput);

    /** \brief access the created boundary integration cells on the contact interface
     *
     *  \param face_ele_id (in) : face element ID
     *
     *  \return This method returns a pointer to the boundery integration cell vector
     *          or a NULL pointer, if there are no boundary integration cells for the
     *          given face element ID.
     *
     *  \author hiermeier \date 01/17 */
    const GEO::BoundaryIntCellPtrs *GetBoundaryIntegrationCells(int face_ele_id) const;

   protected:
    /// do stuff after the Run_Cut routine is finished
    virtual void Post_Run_Cut(bool include_inner);

    virtual void GetPhysicalNodalCoordinates(
        const DRT::Element *ele, LINALG::SerialDenseMatrix &xyze) const;

   private:
    /// access the conditioned parent-face element paired vector
    inline const EleFacePair &CondColParentFaceElementPair() const
    {
      return cond_col_parent_face_element_pair_;
    }

    /// access the face elements corresponding to the xfem parent element
    inline const std::vector<const DRT::FaceElement *> *CondColParentFaceElementPair(
        const DRT::Element *parent_ele) const
    {
      EleFacePair::const_iterator cit = cond_col_parent_face_element_pair_.find(parent_ele);

      if (cit != cond_col_parent_face_element_pair_.end()) return (&cit->second);

      return NULL;
    }

    /// collect the cut background elements
    void GetCutBackgroundElements(std::vector<const DRT::Element *> &cut_elements) const;

    /** \brief Calculate the local coordinates of the cell vertices in the
     *  parent element parameter space coordinate system
     *
     *  \param parent_element (in)  : parent element ( cut element or element face )
     *  \param xyz_cell       (in)  : physical / global coordinates of the cell vertices
     *  \param cell_type      (in)  : cell shape type
     *  \param rs_cell        (out) : local coordinates of the cell vertices
     *
     *  \return \c TRUE, if all the cell vertices are within the parent element.
     *
     *  \author hiermeier \date 01/17 */
    bool GetNodalLocalCoordinatesInParentElement(const DRT::Element *parent_element,
        const Epetra_SerialDenseMatrix &xyz_cell, DRT::Element::DiscretizationType celltype,
        LINALG::SerialDenseMatrix &rs_bc) const;

    void AddToCondBoundaryIntCells(const Teuchos::RCP<GEO::BoundaryIntCell> &bicell,
        int face_ele_id, std::map<int, GEO::BoundaryIntCellPtrs> &cond_bicells_map) const;

    void CreateCondBoundaryIntCells(const DRT::Element *cut_element,
        const GEO::CUT::BoundaryCell &bcell,
        std::map<int, GEO::BoundaryIntCellPtrs> &cond_bicells_map) const;

    void CreateCondBoundaryIntegrationCells(const std::vector<const DRT::Element *> &cut_elements,
        std::map<int, GEO::BoundaryIntCellPtrs> &bicells_map) const;

    void RemoveNonInterfaceNonStandardNodalDofSets();

    void CheckBoundaryCellType(DRT::Element::DiscretizationType distype_bc) const;

    /** \brief Expand level-set values along the element edges perpendicular to
     *  the slave interface in the parameter space
     *
     *  \param xdiscret        (in) : xfem discretization
     *  \param levelset_values (in) : nodal level set values in row layout
     *  \param desired_conds   (in) : desired conditions specifying the interface elements
     *  \param cond_col_parent_face_element_pair (out) : pair of parent and face elements at
     *                                condition interface
     *  \param cond_face_node_col_map (out) : column map of all conditioned nodes on the interface
     *
     *  \return expanded level-set value vector
     *
     *  \author hiermeier \date 01/17 */
    Teuchos::RCP<Epetra_Vector> ExpandLevelSetValues(const DRT::Discretization &xdiscret,
        const Epetra_Vector &levelset_values, const std::vector<std::string> &desired_conds,
        EleFacePair &cond_col_parent_face_element_pair,
        Teuchos::RCP<Epetra_Map> &cond_face_node_col_map) const;

    /** \brief Expand level-set values along the element edges perpendicular to
     *  the slave interface in the parameter space
     *
     *  \param xdiscret        (in) : xfem discretization
     *  \param levelset_values (in) : nodal level set values in row layout
     *  \param desired_cond    (in) : desired condition specifying the interface elements
     *  \param cond_col_parent_face_element_pair (out) : pair of parent and face elements at
     *                                condition interface
     *  \param cond_face_node_col_map (out) : column map of all conditioned nodes on the interface
     *
     *  \return expanded level-set value vector
     *
     *  \author hiermeier \date 01/17 */
    inline Teuchos::RCP<Epetra_Vector> ExpandLevelSetValues(const DRT::Discretization &xdiscret,
        const Epetra_Vector &levelset_values, const std::string &desired_cond,
        EleFacePair &cond_col_parent_face_element_pair,
        Teuchos::RCP<Epetra_Map> &cond_face_node_col_map) const
    {
      const std::vector<std::string> desired_conds(1, desired_cond);
      return ExpandLevelSetValues(xdiscret, levelset_values, desired_conds,
          cond_col_parent_face_element_pair, cond_face_node_col_map);
    }

    /** \brief project the level-set values along the given line
     *
     *  \param line_node_gids       (in) : node GIDs of the current line element
     *  \param num_line_nodes       (in) : number of nodes of the line
     *  \param common_line_node_lid (in) : LID of the line node on the face element
     *                                     (holding the level-set value)
     *  \param expanded_level_set_values (out) : vector which will contain the expanded
     *                                     level-set values.
     *
     *  \author hiermeier \date 01/17 */
    void ProjectLevelSetValuesAlongLine(const int *line_node_gids, unsigned num_line_nodes,
        int common_line_node_lid, Epetra_Vector &expanded_level_set_values) const;

    void CollectConditionSideInfo(const std::map<int, Teuchos::RCP<DRT::Element>> &cond_elements,
        const Epetra_Comm &comm, Teuchos::RCP<Epetra_Map> &cond_node_col_map,
        EleFacePair &cond_col_parent_face_element_pair,
        Teuchos::RCP<Epetra_Map> &cond_face_node_col_map) const;

    const Epetra_Map &CondFaceNodeColMap() const
    {
      if (cond_face_node_col_map_.is_null()) dserror("NULL pointer!");

      return *cond_face_node_col_map_;
    }

   public:
    /// debug output to check the cut algorithm
    void Cut_Debug(const DRT::Element *cut_ele) const;

   private:
    std::map<int, GEO::BoundaryIntCellPtrs> cond_bicells_map_;

    EleFacePair cond_col_parent_face_element_pair_;

    Teuchos::RCP<Epetra_Map> cond_face_node_col_map_;
  };
}  // namespace XCONTACT



#endif /* SRC_DRT_CONTACT_XCONTACT_XCONTACT_CUTWIZARD_H_ */
