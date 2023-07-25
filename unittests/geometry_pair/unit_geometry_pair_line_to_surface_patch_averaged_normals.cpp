/*----------------------------------------------------------------------*/
/*! \file

\brief Unit tests for line to surface geometry pairs that rely on averaged nodal normals.

\level 1
*/
// End doxygen header.


#include <gtest/gtest.h>

#include <Epetra_SerialComm.h>
#include "baci_geometry_pair_element_faces.H"
#include "baci_geometry_pair_element.H"
#include "baci_geometry_pair_scalar_types.H"
#include "unit_geometry_pair_line_to_surface_patch_geometry.H"
#include "unit_geometry_pair_line_to_surface_patch_results.H"


using namespace GEOMETRYPAIR;

namespace
{
  /**
   * \brief Class to test the surface patch functionality of the geometry pairs.
   */
  class GeometryPairLineToSurfacePatchTest : public ::testing::Test
  {
   protected:
    /**
     * \brief Set up the testing environment, is called before each test.
     */
    GeometryPairLineToSurfacePatchTest()
    {
      Teuchos::RCP<Epetra_SerialComm> comm =
          Teuchos::rcp<Epetra_SerialComm>(new Epetra_SerialComm());
      discret_ = Teuchos::rcp(new DRT::Discretization("unit_test", comm));
    }

    /**
     * \brief Return a reference to the connected faces of a face element.
     */
    template <typename A>
    std::map<int, ConnectedFace>& GetConnectedFaces(A& face_element)
    {
      return face_element->connected_faces_;
    }

    /**
     * \brief Get the number of dofs for the beam.
     */
    template <typename A>
    unsigned int GetNBeamDof(A& face_element)
    {
      return face_element->n_beam_dof_;
    }

    //! Pointer to the discretization object that holds the geometry for the tests.
    Teuchos::RCP<DRT::Discretization> discret_;
  };

  /**
   * \brief Test the evaluation of averaged normals on a patch of hex8/quad4 elements.
   */
  TEST_F(GeometryPairLineToSurfacePatchTest, TestSurfacePatchAveragedNormalsQuad4)
  {
    // Define the type of the face elements.
    using surface = GEOMETRYPAIR::t_quad4;
    using scalar_type = GEOMETRYPAIR::line_to_surface_patch_scalar_type;
    using face_element_type = GEOMETRYPAIR::FaceElementPatchTemplate<surface, scalar_type>;

    // Tolerance for the result tests.
    const double eps = 1e-12;

    // Fill the discretization object with the geometry.
    std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>> face_elements_map;
    XtestSurfacePatchQuad4<face_element_type>(discret_, face_elements_map);

    // Load the result vectors.
    std::vector<double> reference_normals, current_normals, position;
    std::vector<std::vector<double>> current_normals_derivative, position_derivative;
    std::vector<std::vector<std::vector<double>>> current_normals_derivative_2,
        position_derivative_2;
    XtestSurfacePatchQuad4Results(reference_normals, current_normals, current_normals_derivative,
        current_normals_derivative_2, position, position_derivative, position_derivative_2);

    // Face element that will be analyzed.
    const unsigned int investigated_face_element_volume_id = 14;
    Teuchos::RCP<face_element_type> face_element = Teuchos::rcp_dynamic_cast<face_element_type>(
        face_elements_map[investigated_face_element_volume_id]);

    // Offset in the derivatives for the beam dof.
    const unsigned int beam_dof_offset = GetNBeamDof(face_element);

    // Setup all face elements and get the patch information.
    for (auto& face_element_map_iterator : face_elements_map)
      face_element_map_iterator.second->Setup(discret_, face_elements_map);

    {
      // Check if the GID are correct.
      std::vector<int> patch_dof_gid_reference = {126, 127, 128, 111, 112, 113, 117, 118, 119, 129,
          130, 131, 120, 121, 122, 102, 103, 104, 99, 100, 101, 108, 109, 110, 114, 115, 116};
      EXPECT_EQ(face_element->GetPatchGID().size(), patch_dof_gid_reference.size());
      for (unsigned int i = 0; i < face_element->GetPatchGID().size(); i++)
        EXPECT_EQ(face_element->GetPatchGID()[i], patch_dof_gid_reference[i]);

      // Check if the local node ID map of the connected faces to the main face could be found.
      EXPECT_EQ(GetConnectedFaces(face_element).size(), 3);

      EXPECT_EQ(GetConnectedFaces(face_element)[10].node_lid_map_.size(), 1);
      EXPECT_EQ(GetConnectedFaces(face_element)[10].node_lid_map_[3], 1);
      EXPECT_EQ(GetConnectedFaces(face_element)[10].my_node_patch_lid_.size(), 4);
      EXPECT_EQ(GetConnectedFaces(face_element)[10].my_node_patch_lid_[0], 5);
      EXPECT_EQ(GetConnectedFaces(face_element)[10].my_node_patch_lid_[1], 6);
      EXPECT_EQ(GetConnectedFaces(face_element)[10].my_node_patch_lid_[2], 7);
      EXPECT_EQ(GetConnectedFaces(face_element)[10].my_node_patch_lid_[3], 1);

      EXPECT_EQ(GetConnectedFaces(face_element)[11].node_lid_map_.size(), 2);
      EXPECT_EQ(GetConnectedFaces(face_element)[11].node_lid_map_[0], 1);
      EXPECT_EQ(GetConnectedFaces(face_element)[11].node_lid_map_[3], 2);
      EXPECT_EQ(GetConnectedFaces(face_element)[11].my_node_patch_lid_.size(), 4);
      EXPECT_EQ(GetConnectedFaces(face_element)[11].my_node_patch_lid_[0], 1);
      EXPECT_EQ(GetConnectedFaces(face_element)[11].my_node_patch_lid_[1], 7);
      EXPECT_EQ(GetConnectedFaces(face_element)[11].my_node_patch_lid_[2], 8);
      EXPECT_EQ(GetConnectedFaces(face_element)[11].my_node_patch_lid_[3], 2);

      EXPECT_EQ(GetConnectedFaces(face_element)[13].node_lid_map_.size(), 2);
      EXPECT_EQ(GetConnectedFaces(face_element)[13].node_lid_map_[2], 1);
      EXPECT_EQ(GetConnectedFaces(face_element)[13].node_lid_map_[3], 0);
      EXPECT_EQ(GetConnectedFaces(face_element)[13].my_node_patch_lid_.size(), 4);
      EXPECT_EQ(GetConnectedFaces(face_element)[13].my_node_patch_lid_[0], 4);
      EXPECT_EQ(GetConnectedFaces(face_element)[13].my_node_patch_lid_[1], 5);
      EXPECT_EQ(GetConnectedFaces(face_element)[13].my_node_patch_lid_[2], 1);
      EXPECT_EQ(GetConnectedFaces(face_element)[13].my_node_patch_lid_[3], 0);
    }

    // Calculate the averaged reference normals on the face.
    face_element->CalculateAveragedReferenceNormals(face_elements_map);
    {
      for (unsigned int i = 0; i < reference_normals.size(); i++)
        EXPECT_NEAR((*face_element->GetReferenceNormals())(i), reference_normals[i], eps);
    }

    // Set the state in the face element, here also the FAD variables for each patch are set.
    auto gid_map = Teuchos::rcp(new Epetra_Map(
        discret_->NumGlobalNodes() * 3, discret_->NumGlobalNodes() * 3, 0, discret_->Comm()));
    auto displacement_vector = Teuchos::rcp(new Epetra_Vector(*gid_map));
    for (int i = 0; i < displacement_vector->GlobalLength(); i++)
      (*displacement_vector)[i] = i * 0.01;
    face_element->SetState(displacement_vector, face_elements_map);
    {
      // Check the values of the averaged normals.
      for (unsigned int i_dof = 0; i_dof < 3 * surface::n_nodes_; i_dof++)
      {
        EXPECT_NEAR(CORE::FADUTILS::CastToDouble((*face_element->GetCurrentNormals())(i_dof)),
            current_normals[i_dof], eps);
        for (unsigned int i_der = 0; i_der < face_element->GetPatchGID().size(); i_der++)
        {
          EXPECT_NEAR(CORE::FADUTILS::CastToDouble(
                          (*face_element->GetCurrentNormals())(i_dof).dx(beam_dof_offset + i_der)),
              current_normals_derivative[i_dof][i_der], eps);
          for (unsigned int i_der_2 = 0; i_der_2 < face_element->GetPatchGID().size(); i_der_2++)
          {
            EXPECT_NEAR(CORE::FADUTILS::CastToDouble((*face_element->GetCurrentNormals())(i_dof)
                                                         .dx(beam_dof_offset + i_der)
                                                         .dx(beam_dof_offset + i_der_2)),
                current_normals_derivative_2[i_dof][i_der][i_der_2], eps);
          }
        }
      }

      // Check an surface position on the element.
      CORE::LINALG::Matrix<3, 1, double> xi;
      xi(0) = 0.2;
      xi(1) = -0.8;
      xi(2) = 0.69;
      CORE::LINALG::Matrix<3, 1, scalar_type> r;
      GEOMETRYPAIR::EvaluateSurfacePosition<surface>(xi, face_element->GetFacePosition(), r,
          face_element->GetDrtFaceElement(), face_element->GetCurrentNormals());
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      {
        EXPECT_NEAR(CORE::FADUTILS::CastToDouble(r(i_dim)), position[i_dim], eps);
        for (unsigned int i_der = 0; i_der < face_element->GetPatchGID().size(); i_der++)
        {
          EXPECT_NEAR(CORE::FADUTILS::CastToDouble(r(i_dim).dx(beam_dof_offset + i_der)),
              position_derivative[i_dim][i_der], eps);
          for (unsigned int i_der_2 = 0; i_der_2 < face_element->GetPatchGID().size(); i_der_2++)
            EXPECT_NEAR(CORE::FADUTILS::CastToDouble(
                            r(i_dim).dx(beam_dof_offset + i_der).dx(beam_dof_offset + i_der_2)),
                position_derivative_2[i_dim][i_der][i_der_2], eps);
        }
      }
    }
  }

}  // namespace
