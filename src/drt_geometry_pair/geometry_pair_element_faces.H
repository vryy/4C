/*----------------------------------------------------------------------*/
/*! \file

\brief Element classes that represent faces, i.e. surface elements.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_ELEMENT_FACES_H_
#define GEOMETRY_PAIR_ELEMENT_FACES_H_


#include "geometry_pair_element_functions.H"

#include "../drt_lib/drt_element.H"

#include <unordered_map>


namespace GEOMETRYPAIR
{
  /**
   * \brief Base, non templated class for an object that represents a surface element.
   */
  class FaceElement
  {
   public:
    /**
     * \brief Constructor.
     * @param face_element (in) Pointer to the DRT face element.
     */
    FaceElement(const Teuchos::RCP<const DRT::FaceElement>& face_element)
        : drt_face_element_(face_element), part_of_pair_(false){};

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
     * \brief Set the reference position of this element. Has to be implemented in derived class.
     * @param face_elements (in) Vector with all face elements in the surface condition.
     */
    virtual void Setup() = 0;

    /**
     * \brief Set the needed displacement vectors for this face. Has to be implemented in derived
     * class.
     *
     * @param discret (in) Discretization.
     * @param displacement (in) Current displacement vector.
     */
    virtual void SetState(const Teuchos::RCP<const DRT::Discretization>& discret,
        const Teuchos::RCP<const Epetra_Vector>& displacement) = 0;

    /**
     * \brief Calculate the averaged normals on the nodes of this face. Has to be implemented in
     * derived class.
     *
     * @param face_elements (in) Vector with all face elements in the surface condition.
     */
    virtual void CalculateAveragedNormals(
        const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements) = 0;

    /**
     * \brief Calculate the nodal normal at the node with GID node_gid. If this node is not part of
     * this face, false is returned. (Has to be implemented in derived class)
     *
     * @param node_gid (in) GID of the node at which the normals should be calculated.
     * @param reference_normal (out) Normal at the node in the reference configuration.
     * @param current_normal (out) Normal at the node in the current configuration.
     * @return True if node if part of this face, false otherwise.
     */
    virtual bool EvaluateNodalNormal(const int node_gid,
        LINALG::Matrix<3, 1, double>& reference_normal,
        LINALG::Matrix<3, 1, double>& current_normal) const = 0;

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
     * \brief Return a normal on this face as a double.
     *
     * @param xi (in) Parameter coordinate on the face.
     * @param n (out) Normal on the face.
     * @param reference (in) If the reference position or the current position should be returned.
     */
    virtual void EvaluateFaceNormalDouble(const LINALG::Matrix<2, 1, double>& xi,
        LINALG::Matrix<3, 1, double>& n, bool reference = false) const = 0;

    /**
     * \brief Return a normal of the averaged normal field on this face as a double.
     *
     * @param xi (in) Parameter coordinate on the face.
     * @param n (out) Normal on the face.
     * @param reference (in) If the reference position or the current position should be returned.
     */
    virtual void EvaluateFaceAveragedNormalDouble(const LINALG::Matrix<2, 1, double>& xi,
        LINALG::Matrix<3, 1, double>& n, bool reference = false) const = 0;

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

   protected:
    /**
     * \brief Get the face elements connected to this face element (within the surface condition).
     * @param face_elements (in) Map with all face elements of the surface condition.
     * @param connected_faces (out) Vector with the connected face elements, not including this
     * face.
     */
    void GetConnectedFaces(
        const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements,
        std::vector<Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& connected_faces) const;

    /**
     * \brief Get the local to global indices for this surface patch.
     * @param connected_faces (in) Vector with the connected face elements.
     * @param ltg (out) Local to global indices for this surface patch.
     */
    void GetPatchLocalToGlobalIndices(
        const std::vector<Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& connected_faces,
        std::vector<int>& ltg) const;

   protected:
    //! Pointer to the drt face element.
    Teuchos::RCP<const DRT::FaceElement> drt_face_element_;

    //! Flag if this face element is part of a contact pair, i.e. if it has evaluate it's averaged
    //! normals.
    bool part_of_pair_;
  };


  /**
   * \brief An object that represents a surface element and stores averaged normal information as
   * well as position information of the face.
   *
   * @tparam surface Type of surface element.
   */
  template <typename surface>
  class FaceElementTemplate : public FaceElement
  {
   public:
    /**
     * \brief Constructor (derived).
     */
    FaceElementTemplate(const Teuchos::RCP<const DRT::FaceElement>& face_element)
        : FaceElement(face_element),
          face_reference_position_(true),
          face_position_(true),
          reference_normals_(true),
          current_normals_(true),
          face_dof_gid_()
    {
    }

    /**
     * \brief Destructor.
     */
    ~FaceElementTemplate(){};

    /**
     * \brief Set the reference position of this element (derived).
     */
    void Setup() override;

    /**
     * \brief Set the needed displacement vectors for this face (derived).
     */
    void SetState(const Teuchos::RCP<const DRT::Discretization>& discret,
        const Teuchos::RCP<const Epetra_Vector>& displacement) override;

    /**
     * \brief Calculate the averaged normals on the nodes of this face (derived).
     */
    void CalculateAveragedNormals(
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
    const LINALG::Matrix<surface::n_dof_, 1, double>& GetFacePosition() const
    {
      return face_position_;
    }

    /**
     * \brief Return the reference normals on this face.
     */
    const LINALG::Matrix<3 * surface::n_nodes_, 1, double>* GetReferenceNormals() const
    {
      return &reference_normals_;
    }

    /**
     * \brief Return the current normals on this face.
     */
    const LINALG::Matrix<3 * surface::n_nodes_, 1, double>* GetCurrentNormals() const
    {
      return &current_normals_;
    }

    /**
     * \brief Return the position on this face as a double.
     *
     * @param xi (in) Parameter coordinate on the face.
     * @param r (out) Position on the face.
     * @param reference (in) If the reference position or the current position should be returned.
     */
    void EvaluateFacePositionDouble(const LINALG::Matrix<2, 1, double>& xi,
        LINALG::Matrix<3, 1, double>& r, bool reference = false) const override;

    /**
     * \brief Return a normal on this face as a double.
     *
     * @param xi (in) Parameter coordinate on the face.
     * @param n (out) Normal on the face.
     * @param reference (in) If the reference position or the current position should be returned.
     */
    void EvaluateFaceNormalDouble(const LINALG::Matrix<2, 1, double>& xi,
        LINALG::Matrix<3, 1, double>& n, bool reference = false) const override;

    /**
     * \brief Return a normal of the averaged normal field on this face as a double.
     *
     * @param xi (in) Parameter coordinate on the face.
     * @param n (out) Normal on the face.
     * @param reference (in) If the reference position or the current position should be returned.
     */
    void EvaluateFaceAveragedNormalDouble(const LINALG::Matrix<2, 1, double>& xi,
        LINALG::Matrix<3, 1, double>& n, bool reference = false) const override;

   private:
    /**
     * \brief Calculate the nodal normal at the node with GID node_gid. If this node is not part of
     * this face, false is returned. (derived)
     */
    bool EvaluateNodalNormal(const int node_gid, LINALG::Matrix<3, 1, double>& reference_normal,
        LINALG::Matrix<3, 1, double>& current_normal) const override;

   private:
    //! Reference position.
    LINALG::Matrix<surface::n_dof_, 1, double> face_reference_position_;

    //! Current position.
    LINALG::Matrix<surface::n_dof_, 1, double> face_position_;

    //! Reference normals.
    LINALG::Matrix<3 * surface::n_nodes_, 1, double> reference_normals_;

    //! Current normals.
    LINALG::Matrix<3 * surface::n_nodes_, 1, double> current_normals_;

    //! Global DOF IDs of this face element.
    std::vector<int> face_dof_gid_;
  };


  /**
   * A few methods of the FaceElementTemplate class have to be specialized for nurbs9 faces, since
   * they do not need averaged nodal normals due to their C1-continuous discretization.
   *
   * Full explicit template specializations have to be inline.
   * https://stackoverflow.com/questions/4445654/multiple-definition-of-template-specialization-when-using-different-objects/4445772
   *
   * Due to a gcc bug we have to perform the template specialization inside the namespace.
   * (https://stackoverflow.com/questions/25594644/warning-specialization-of-template-in-different-namespace)
   */
  template <>
  inline void FaceElementTemplate<t_nurbs9>::CalculateAveragedNormals(
      const std::unordered_map<int, Teuchos::RCP<FaceElement>>& face_elements){};
  template <>
  inline const LINALG::Matrix<3 * t_nurbs9::n_nodes_, 1, double>*
  FaceElementTemplate<t_nurbs9>::GetReferenceNormals() const
  {
    return nullptr;
  };
  template <>
  inline const LINALG::Matrix<3 * t_nurbs9::n_nodes_, 1, double>*
  FaceElementTemplate<t_nurbs9>::GetCurrentNormals() const
  {
    return nullptr;
  };


  /**
   * \brief Create the templated version of the face element.
   * @param drt_face_element (in) Pointer to the DRT face element.
   * @return RCP to the created GEOMETRYPAIR FaceElement.
   */
  Teuchos::RCP<FaceElement> FaceElementFactory(
      const Teuchos::RCP<const DRT::FaceElement>& drt_face_element);

}  // namespace GEOMETRYPAIR

#endif
