/*----------------------------------------------------------------------*/
/*! \file

\brief Element classes that represent faces, i.e. surface elements.

\level 1
*/
// End doxygen header.


#ifndef FOUR_C_GEOMETRY_PAIR_ELEMENT_FACES_HPP
#define FOUR_C_GEOMETRY_PAIR_ELEMENT_FACES_HPP


#include "4C_config.hpp"

#include "4C_geometry_pair_element.hpp"
#include "4C_lib_element.hpp"

#include <unordered_map>

namespace
{
  class GeometryPairLineToSurfacePatchTest;
}

FOUR_C_NAMESPACE_OPEN

// Forward declarations.
namespace GEOMETRYPAIR
{
  struct ConnectedFace;
}  // namespace GEOMETRYPAIR
namespace INPAR
{
  namespace GEOMETRYPAIR
  {
    enum class SurfaceNormals;
  }
}  // namespace INPAR


namespace GEOMETRYPAIR
{
  /**
   * \brief This structure "converts" a DRT face element type to the underlying volume element.
   */
  template <CORE::FE::CellType discretization>
  struct FaceDiscretizationTypeToVolumeElement
  {
    using volume_type_ = void;
  };
  template <>
  struct FaceDiscretizationTypeToVolumeElement<CORE::FE::CellType::quad4>
  {
    using volume_type_ = t_hex8;
  };
  template <>
  struct FaceDiscretizationTypeToVolumeElement<CORE::FE::CellType::quad8>
  {
    using volume_type_ = t_hex20;
  };
  template <>
  struct FaceDiscretizationTypeToVolumeElement<CORE::FE::CellType::quad9>
  {
    using volume_type_ = t_hex27;
  };

  /**
   * \brief Utility structure to represent the connection of faces to a patch.
   */
  struct ConnectedFace
  {
    //! This vector stores the LIDs of the nodes of this connected face corresponding to the nodes
    //! of the patch.
    std::vector<int> my_node_patch_lid_;

    //! This map represents the link between the nodes of a connected face and the patch face. The
    //! key is the LID of the node on the connected face and the value is the LID of the node on the
    //! patch face.
    std::map<int, int> node_lid_map_;
  };

  /**
   * \brief Base, non templated class for an object that represents a surface element.
   */
  class FaceElement
  {
    friend GeometryPairLineToSurfacePatchTest;

   public:
    /**
     * \brief Constructor.
     * @param face_element (in) Pointer to the DRT face element.
     */
    FaceElement(const Teuchos::RCP<const DRT::FaceElement>& face_element)
        : drt_face_element_(face_element), part_of_pair_(false), patch_dof_gid_(){};

    /**
     * \brief Destructor.
     */
    virtual ~FaceElement() = default;

    /**
     * \brief Get the RCP to the DRT face element.
     * @return RCP to the DRT face element.
     */
    const DRT::FaceElement* GetDrtFaceElement() const { return drt_face_element_.getRawPtr(); }

    /**
     * \brief Setup the object. Has to be implemented in derived class.
     * @param discret (in) Pointer to the discretization.
     * @param face_elements (in) Vector with all face elements in the surface condition.
     */
    virtual void Setup(const Teuchos::RCP<const DRT::Discretization>& discret,
        const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements) = 0;

    /**
     * \brief Set the needed displacement vectors for this face. Has to be implemented in derived
     * class.
     *
     * @param displacement (in) Current displacement vector.
     * @param face_elements (in) Map with all the faces in this condition.
     */
    virtual void set_state(const Teuchos::RCP<const Epetra_Vector>& displacement,
        const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements) = 0;

    /**
     * \brief Calculate the averaged normals on the nodes of this face. Has to be implemented in
     * derived class.
     *
     * @param face_elements (in) Vector with all face elements in the surface condition.
     */
    virtual void calculate_averaged_reference_normals(
        const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements){};

    /**
     * \brief Return the position on this face as a double.
     *
     * @param xi (in) Parameter coordinate on the face.
     * @param r (out) Position on the face.
     * @param reference (in) If the reference position or the current position should be returned.
     */
    virtual void evaluate_face_position_double(const CORE::LINALG::Matrix<2, 1, double>& xi,
        CORE::LINALG::Matrix<3, 1, double>& r, bool reference = false) const = 0;

    /**
     * \brief Return a normal on the element.
     *
     * If an averaged normal is requested and the element does not have averaged normals (e.g.
     * nurbs) a zero vector is returned.
     *
     * @param xi (in) Parameter coordinate on the face.
     * @param n (out) Normal on the face.
     * @param reference (in) If the reference position or the current position should be returned.
     * @param averaged_normal (in) If an averaged normal or an geometrical normal should be
     * returned.
     */
    virtual void evaluate_face_normal_double(const CORE::LINALG::Matrix<2, 1, double>& xi,
        CORE::LINALG::Matrix<3, 1, double>& n, const bool reference,
        const bool averaged_normal) const = 0;

    /**
     * \brief Return if the face is part of a pair.
     * @return True if the face is part of a pair, false otherwise.
     */
    bool IsPartOfPair() const { return part_of_pair_; }

    /**
     * \brief Set the part_of_pair_ flag.
     * @param part_of_pair (in) Value to set.
     */
    void SetPartOfPair(bool part_of_pair) { part_of_pair_ = part_of_pair; }

    /**
     * \brief Get the local to global indices for this surface patch.
     * @return Local to global indices for this surface patch.
     */
    const std::vector<int>& GetPatchGID() const { return patch_dof_gid_; }

   protected:
    //! Pointer to the drt face element.
    Teuchos::RCP<const DRT::FaceElement> drt_face_element_;

    //! Flag if this face element is part of a contact pair, i.e. if it has evaluate it's averaged
    //! normals.
    bool part_of_pair_;

    //! Global DOF IDs of this face patch.
    std::vector<int> patch_dof_gid_;
  };


  /**
   * \brief An object that represents a surface element which only depends on a single face. For
   * example nurbs faces since they are C1 continuous.
   *
   * @tparam surface Type of surface element.
   * @tparam scalar_type Scalar type for FAD evaluations.
   */
  template <typename surface, typename scalar_type>
  class FaceElementTemplate : public FaceElement
  {
   public:
    //! Shortcut to the type of this templated object.
    using my_type = FaceElementTemplate<surface, scalar_type>;

   public:
    /**
     * \brief Constructor (derived).
     */
    FaceElementTemplate(const Teuchos::RCP<const DRT::FaceElement>& face_element)
        : FaceElement(face_element), n_dof_other_element_(0){};


    /**
     * \brief Get the face GIDs and set the reference configuration.
     *
     * @param discret (in) Pointer to the discretization.
     * @param face_elements (in) Vector with all face elements in the surface condition.
     */
    void Setup(const Teuchos::RCP<const DRT::Discretization>& discret,
        const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
        override;

    /**
     * \brief Set the needed displacement vectors for this face (derived).
     */
    void set_state(const Teuchos::RCP<const Epetra_Vector>& displacement,
        const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
        override;

    /**
     * \brief Return the reference element data of this face.
     */
    const auto& get_face_reference_element_data() const { return face_reference_position_; }

    /**
     * \brief Return the current element data of this face.
     */
    const auto& GetFaceElementData() const { return face_position_; }

    /**
     * \brief Return the reference normals on this face.
     */
    virtual const CORE::LINALG::Matrix<3 * surface::n_nodes_, 1, double>* GetReferenceNormals()
        const
    {
      return nullptr;
    }

    /**
     * \brief Return the current normals on this face.
     */
    virtual const CORE::LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* GetCurrentNormals()
        const
    {
      return nullptr;
    }

    /**
     * \brief Return the position on this face as a double. (derived)
     */
    void evaluate_face_position_double(const CORE::LINALG::Matrix<2, 1, double>& xi,
        CORE::LINALG::Matrix<3, 1, double>& r, bool reference = false) const override;

    /**
     * \brief Return a normal on the element. (derived)
     */
    void evaluate_face_normal_double(const CORE::LINALG::Matrix<2, 1, double>& xi,
        CORE::LINALG::Matrix<3, 1, double>& n, const bool reference,
        const bool averaged_normal) const override;

    /**
     * \brief Return the number of DOFs of the interacting element.
     */
    [[nodiscard]] unsigned int get_number_of_dof_other_element() const
    {
      return n_dof_other_element_;
    }

    /**
     * \brief Return the number of DOFs of the interacting element.
     */
    void set_number_of_dof_other_element(const unsigned int n_dof_other_element)
    {
      if (n_dof_other_element_ == 0)
      {
        n_dof_other_element_ = n_dof_other_element;
      }
      else if (n_dof_other_element_ != n_dof_other_element)
      {
        FOUR_C_THROW(
            "FaceElementTemplate only allows other elements with the same number of DOFs. You "
            "already have an other element with %d DOFs and now want to set another one with %d "
            "DOFs",
            n_dof_other_element_, n_dof_other_element);
      }
    }

   protected:
    //! Reference position.
    ElementData<surface, double> face_reference_position_;

    //! Current position.
    ElementData<surface, scalar_type> face_position_;

    //! Number of DOFs used for the element that will be interacting with this face. This is
    //! required to correctly set the FAD types. For now, a single face element can only be used in
    //! interactions with other elements of the same type.
    unsigned int n_dof_other_element_;
  };


  /**
   * \brief An object that represents a surface element and stores averaged normal information as
   * well as position information of the face.
   *
   * @tparam surface Type of surface element.
   * @tparam scalar_type Scalar type for FAD evaluations.
   */
  template <typename surface, typename scalar_type>
  class FaceElementPatchTemplate : public FaceElementTemplate<surface, scalar_type>
  {
    friend GeometryPairLineToSurfacePatchTest;

   public:
    //! Shortcut to the type of this templated object.
    using my_type = FaceElementPatchTemplate<surface, scalar_type>;

    //! Shortcut to the base class.
    using base_class = FaceElementTemplate<surface, scalar_type>;

   public:
    /**
     * \brief Constructor (derived).
     * @param evaluate_current_normals (in) If the current normals should be evaluated.
     */
    FaceElementPatchTemplate(const Teuchos::RCP<const DRT::FaceElement>& face_element,
        const bool evaluate_current_normals)
        : base_class(face_element),
          connected_faces_(),
          evaluate_current_normals_(evaluate_current_normals)
    {
    }

    /**
     * \brief Set the patch information of the patch connected to this face element.
     *
     * Find the connected faces of this patch. Order the nodes, starting with the nodes of the main
     * face element and then add the nodes of the connected faces. Also create a vector with the
     * GIDs of this patch.
     *
     * @param discret (in) Pointer to the discretization.
     * @param face_elements (in) Vector with all face elements in the surface condition.
     */
    void Setup(const Teuchos::RCP<const DRT::Discretization>& discret,
        const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
        override;

    /**
     * \brief Set the needed displacement vectors for this face (derived).
     */
    void set_state(const Teuchos::RCP<const Epetra_Vector>& displacement,
        const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
        override;

    /**
     * \brief Calculate the averaged normals on the nodes of this face (derived).
     */
    void calculate_averaged_reference_normals(
        const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
        override;

    /**
     * \brief Return a normal on the element. (derived)
     */
    void evaluate_face_normal_double(const CORE::LINALG::Matrix<2, 1, double>& xi,
        CORE::LINALG::Matrix<3, 1, double>& n, const bool reference,
        const bool averaged_normal) const override;

   private:
    /**
     * \brief Average normals at the face nodes.
     * @param normals (in) Sum of all normals at a each node.
     * @param averaged_normals (out) Averaged normals, already in the vector format needed for this
     * element.
     */
    template <typename T>
    void average_nodal_normals(
        CORE::LINALG::Matrix<surface::n_nodes_, 1, CORE::LINALG::Matrix<3, 1, T>>& normals,
        CORE::LINALG::Matrix<3 * surface::n_nodes_, 1, T>& averaged_normals) const;

   protected:
    //! Store the relevant information of the connected faces. The key of this map is the GID of the
    //! connected face element.
    std::map<int, ConnectedFace> connected_faces_;

    //! If the current normals should be evaluated. If they should not be evaluated, then no
    //! information about the patch is required.
    bool evaluate_current_normals_;
  };


  /**
   * \brief Class to handle extended volume coupling.
   */
  template <typename surface, typename scalar_type, typename volume>
  class FaceElementTemplateExtendedVolume : public FaceElementTemplate<surface, scalar_type>
  {
   public:
    //! Shortcut to the type of this templated object.
    using my_type = FaceElementTemplateExtendedVolume<surface, scalar_type, volume>;

    //! Shortcut to the base class.
    using base_class = FaceElementTemplate<surface, scalar_type>;

   public:
    /**
     * \brief Constructor (derived).
     */
    FaceElementTemplateExtendedVolume(const Teuchos::RCP<const DRT::FaceElement>& face_element)
        : base_class(face_element),
          surface_dof_lid_map_(true),
          face_to_volume_coordinate_axis_map_(true),
          face_to_volume_coordinate_axis_factor_(true),
          third_direction_(-1),
          third_direction_factor_(0){};

    /**
     * \brief Get the face GIDs and set the reference configuration.
     *
     * @param discret (in) Pointer to the discretization.
     * @param face_elements (in) Vector with all face elements in the surface condition.
     */
    void Setup(const Teuchos::RCP<const DRT::Discretization>& discret,
        const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
        override;

    /**
     * \brief Set the needed displacement vectors for this face (derived).
     */
    void set_state(const Teuchos::RCP<const Epetra_Vector>& displacement,
        const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
        override;

    /**
     * \brief Calculate the nodal normals for this face element.
     * @tparam scalar_type_normal Scalar type for normals
     * @param volume_position (in) Face position vector of the volume.
     * @param surface_position (in) Face position vector of the surface.
     * @param normals (out) Normals on the nodes.
     */
    template <typename scalar_type_normal>
    void CalculateNormals(
        const GEOMETRYPAIR::ElementData<volume, scalar_type_normal>& volume_position,
        const GEOMETRYPAIR::ElementData<surface, scalar_type_normal>& surface_position,
        CORE::LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type_normal>& normals) const;

    /**
     * \brief Return a normal on the element. (derived)
     */
    void evaluate_face_normal_double(const CORE::LINALG::Matrix<2, 1, double>& xi,
        CORE::LINALG::Matrix<3, 1, double>& n, const bool reference,
        const bool averaged_normal) const override;

   protected:
    /**
     * \brief Convert the face element parameter coordinates to the volume element parameter
     * coordinates.
     */
    void xi_face_to_xi_volume(const CORE::LINALG::Matrix<2, 1, double>& xi_face,
        CORE::LINALG::Matrix<3, 1, double>& xi_volume) const;

   protected:
    //! Map between surface DOFs and volume DOFs.
    CORE::LINALG::Matrix<surface::n_dof_, 1, int> surface_dof_lid_map_;

    //! Reference position.
    GEOMETRYPAIR::ElementData<volume, double> volume_reference_position_;

    //! Current position.
    GEOMETRYPAIR::ElementData<volume, scalar_type> volume_position_;

    //! Map the face coordinate axis to the volume coordinate axis.
    CORE::LINALG::Matrix<2, 1, int> face_to_volume_coordinate_axis_map_;

    //! Factor for the transformed coordinate axis.
    CORE::LINALG::Matrix<2, 1, int> face_to_volume_coordinate_axis_factor_;

    //! Third, i.e. normal direction.
    int third_direction_;

    //! Factor for normal direction.
    int third_direction_factor_;
  };


  /**
   * \brief Create the templated version of the face element.
   * @param drt_face_element (in) Pointer to the DRT face element.
   * @param fad_order (in) Order of the created FAD type (0 means double).
   * @param surface_normal_strategy (in) strategy to be used for surface normals.
   * @return RCP to the created GEOMETRYPAIR FaceElement.
   */
  Teuchos::RCP<FaceElement> FaceElementFactory(
      const Teuchos::RCP<const DRT::FaceElement>& drt_face_element, const int fad_order,
      const INPAR::GEOMETRYPAIR::SurfaceNormals surface_normal_strategy);

}  // namespace GEOMETRYPAIR

FOUR_C_NAMESPACE_CLOSE

#endif
