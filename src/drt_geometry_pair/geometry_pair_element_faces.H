/*----------------------------------------------------------------------*/
/*! \file

\brief Element classes that represent faces, i.e. surface elements.

\level 1
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_ELEMENT_FACES_H_
#define GEOMETRY_PAIR_ELEMENT_FACES_H_


#include "geometry_pair_element_functions.H"

#include "../drt_lib/drt_element.H"

#include <unordered_map>


// Forward declarations.
namespace GEOMETRYPAIR
{
  struct ConnectedFace;
  class GeometryPairLineToSurfacePatch_TestSuite;
}  // namespace GEOMETRYPAIR


namespace GEOMETRYPAIR
{
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
    friend GEOMETRYPAIR::GeometryPairLineToSurfacePatch_TestSuite;

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
    virtual ~FaceElement(){};

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
    virtual void SetState(const Teuchos::RCP<const Epetra_Vector>& displacement,
        const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements) = 0;

    /**
     * \brief Calculate the averaged normals on the nodes of this face. Has to be implemented in
     * derived class.
     *
     * @param face_elements (in) Vector with all face elements in the surface condition.
     */
    virtual void CalculateAveragedReferenceNormals(
        const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements){};

    /**
     * \brief Return the position on this face as a double.
     *
     * @param xi (in) Parameter coordinate on the face.
     * @param r (out) Position on the face.
     * @param reference (in) If the reference position or the current position should be returned.
     */
    virtual void EvaluateFacePositionDouble(const LINALG::Matrix<2, 1, double>& xi,
        LINALG::Matrix<3, 1, double>& r, bool reference = false) const = 0;

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
    virtual void EvaluateFaceNormalDouble(const LINALG::Matrix<2, 1, double>& xi,
        LINALG::Matrix<3, 1, double>& n, const bool reference,
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
        : FaceElement(face_element), face_reference_position_(true), face_position_(true){};

    /**
     * \brief Destructor.
     */
    virtual ~FaceElementTemplate(){};

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
    void SetState(const Teuchos::RCP<const Epetra_Vector>& displacement,
        const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
        override;

    /**
     * \brief Return the reference position of this face.
     */
    const LINALG::Matrix<surface::n_dof_, 1, double>& GetFaceReferencePosition() const
    {
      return face_reference_position_;
    }

    /**
     * \brief Return the current position of this face.
     */
    const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& GetFacePosition() const
    {
      return face_position_;
    }

    /**
     * \brief Return the reference normals on this face.
     */
    virtual const LINALG::Matrix<3 * surface::n_nodes_, 1, double>* GetReferenceNormals() const
    {
      return nullptr;
    }

    /**
     * \brief Return the current normals on this face.
     */
    virtual const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* GetCurrentNormals() const
    {
      return nullptr;
    }

    /**
     * \brief Return the position on this face as a double. (derived)
     */
    void EvaluateFacePositionDouble(const LINALG::Matrix<2, 1, double>& xi,
        LINALG::Matrix<3, 1, double>& r, bool reference = false) const override;

    /**
     * \brief Return a normal on the element. (derived)
     */
    void EvaluateFaceNormalDouble(const LINALG::Matrix<2, 1, double>& xi,
        LINALG::Matrix<3, 1, double>& n, const bool reference,
        const bool averaged_normal) const override;

   protected:
    //! Reference position.
    LINALG::Matrix<surface::n_dof_, 1, double> face_reference_position_;

    //! Current position.
    LINALG::Matrix<surface::n_dof_, 1, scalar_type> face_position_;

    //! Number of DOFs used for a beam element.
    static const unsigned int n_beam_dof_ = 12;
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
    friend GEOMETRYPAIR::GeometryPairLineToSurfacePatch_TestSuite;

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
          evaluate_current_normals_(evaluate_current_normals),
          reference_normals_(true),
          current_normals_(true)
    {
    }

    /**
     * \brief Destructor.
     */
    ~FaceElementPatchTemplate(){};

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
    void SetState(const Teuchos::RCP<const Epetra_Vector>& displacement,
        const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
        override;

    /**
     * \brief Calculate the averaged normals on the nodes of this face (derived).
     */
    void CalculateAveragedReferenceNormals(
        const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements)
        override;

    /**
     * \brief Return the reference normals on this face.
     */
    const LINALG::Matrix<3 * surface::n_nodes_, 1, double>* GetReferenceNormals() const override
    {
      return &reference_normals_;
    }

    /**
     * \brief Return the current normals on this face.
     */
    const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* GetCurrentNormals() const override
    {
      if (evaluate_current_normals_)
        return &current_normals_;
      else
        return nullptr;
    }

    /**
     * \brief Return a normal on the element. (derived)
     */
    void EvaluateFaceNormalDouble(const LINALG::Matrix<2, 1, double>& xi,
        LINALG::Matrix<3, 1, double>& n, const bool reference,
        const bool averaged_normal) const override;

   private:
    /**
     * \brief Average normals at the face nodes.
     * @param normals (in) Sum of all normals at a each node.
     * @param averaged_normals (out) Averaged normals, already in the vector format needed for this
     * element.
     */
    template <typename T>
    void AverageNodalNormals(LINALG::Matrix<surface::n_nodes_, 1, LINALG::Matrix<3, 1, T>>& normals,
        LINALG::Matrix<3 * surface::n_nodes_, 1, T>& averaged_normals) const;

   protected:
    //! Store the relevant information of the connected faces. The key of this map is the GID of the
    //! connected face element.
    std::map<int, ConnectedFace> connected_faces_;

    //! If the current normals should be evaluated. If they should not be evaluated, then no
    //! information about the patch is required.
    bool evaluate_current_normals_;

    //! Reference normals.
    LINALG::Matrix<3 * surface::n_nodes_, 1, double> reference_normals_;

    //! Current normals.
    LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type> current_normals_;
  };


  /**
   * \brief Create the templated version of the face element.
   * @param drt_face_element (in) Pointer to the DRT face element.
   * @param is_fad (in) If the created face element should consist of fad types.
   * @return RCP to the created GEOMETRYPAIR FaceElement.
   */
  Teuchos::RCP<FaceElement> FaceElementFactory(
      const Teuchos::RCP<const DRT::FaceElement>& drt_face_element, const bool is_fad);

}  // namespace GEOMETRYPAIR

#endif
