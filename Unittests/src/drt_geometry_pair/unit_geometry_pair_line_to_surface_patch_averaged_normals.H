/*----------------------------------------------------------------------*/
/*! \file

\brief Unit tests for line to surface geometry pairs that rely on averaged nodal normals.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef UNIT_GEOMETRY_PAIR_LINE_TO_SURFACE_PATCH_AVERAGED_NORMALS_H_
#define UNIT_GEOMETRY_PAIR_LINE_TO_SURFACE_PATCH_AVERAGED_NORMALS_H_


#include <cxxtest/TestSuite.h>

#include <Epetra_SerialComm.h>
#include "src/drt_geometry_pair/geometry_pair_element_faces.H"
#include "src/drt_geometry_pair/geometry_pair_element.H"
#include "src/drt_geometry_pair/geometry_pair_scalar_types.H"
#include "unit_geometry_pair_line_to_surface_patch_geometry.H"
#include "unit_geometry_pair_line_to_surface_patch_results.H"


// Forward declarations.
namespace GEOMETRYPAIR
{
  class GeometryPairLineToSurfacePatch_TestSuite;
}


/**
 * \brief Class to test the surface patch functionality of the geometry pairs.
 */
class GEOMETRYPAIR::GeometryPairLineToSurfacePatch_TestSuite : public CxxTest::TestSuite
{
 public:
  /**
   * \brief Set up the testing environment, is called before each test.
   */
  void setUp()
  {
    Teuchos::RCP<Epetra_SerialComm> comm = Teuchos::rcp<Epetra_SerialComm>(new Epetra_SerialComm());
    discret_ = Teuchos::rcp(new DRT::Discretization("unit_test", comm));
  }

  /**
   * \brief Delete pointers and other class variables, is called after each test.
   */
  void tearDown() { discret_ = Teuchos::null; }

  /**
   * \brief Test the evaluation of averaged normals on a patch of hex8/quad4 elements.
   */
  void test_surface_patch_averaged_normals_quad4()
  {
    // Define the type of the face elements.
    using surface = GEOMETRYPAIR::t_quad4;
    using scalar_type = GEOMETRYPAIR::line_to_surface_patch_scalar_type;
    using face_element_type = GEOMETRYPAIR::FaceElementPatchTemplate<surface, scalar_type>;

    // Offset in the derivatives for the beam dof.
    const unsigned int beam_dof_offset = face_element_type::n_beam_dof_;

    // Tolerance for the result tests.
    const double eps = 1e-12;

    // Fill the discretization object with the geometry.
    std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>> face_elements_map;
    xtest_surface_patch_quad4<face_element_type>(discret_, face_elements_map);

    // Load the result vectors.
    std::vector<double> reference_normals, current_normals, position;
    std::vector<std::vector<double>> current_normals_derivative, position_derivative;
    std::vector<std::vector<std::vector<double>>> current_normals_derivative_2,
        position_derivative_2;
    xtest_surface_patch_quad4_results(reference_normals, current_normals,
        current_normals_derivative, current_normals_derivative_2, position, position_derivative,
        position_derivative_2);

    // Face element that will be analyzed.
    const unsigned int investigated_face_element_volume_id = 14;
    Teuchos::RCP<face_element_type> face_element = Teuchos::rcp_dynamic_cast<face_element_type>(
        face_elements_map[investigated_face_element_volume_id]);

    // Setup all face elements and get the patch information.
    for (auto& face_element_map_iterator : face_elements_map)
      face_element_map_iterator.second->Setup(discret_, face_elements_map);

    {
      // Check if the GID are correct.
      std::vector<int> patch_dof_gid_reference = {126, 127, 128, 111, 112, 113, 117, 118, 119, 129,
          130, 131, 120, 121, 122, 102, 103, 104, 99, 100, 101, 108, 109, 110, 114, 115, 116};
      TS_ASSERT_EQUALS(face_element->patch_dof_gid_.size(), patch_dof_gid_reference.size());
      for (unsigned int i = 0; i < face_element->patch_dof_gid_.size(); i++)
        TS_ASSERT_EQUALS(face_element->patch_dof_gid_[i], patch_dof_gid_reference[i]);

      // Check if the local node ID map of the connected faces to the main face could be found.
      TS_ASSERT_EQUALS(face_element->connected_faces_.size(), 3);

      TS_ASSERT_EQUALS(face_element->connected_faces_[10].node_lid_map_.size(), 1);
      TS_ASSERT_EQUALS(face_element->connected_faces_[10].node_lid_map_[3], 1);
      TS_ASSERT_EQUALS(face_element->connected_faces_[10].my_node_patch_lid_.size(), 4);
      TS_ASSERT_EQUALS(face_element->connected_faces_[10].my_node_patch_lid_[0], 5);
      TS_ASSERT_EQUALS(face_element->connected_faces_[10].my_node_patch_lid_[1], 6);
      TS_ASSERT_EQUALS(face_element->connected_faces_[10].my_node_patch_lid_[2], 7);
      TS_ASSERT_EQUALS(face_element->connected_faces_[10].my_node_patch_lid_[3], 1);

      TS_ASSERT_EQUALS(face_element->connected_faces_[11].node_lid_map_.size(), 2);
      TS_ASSERT_EQUALS(face_element->connected_faces_[11].node_lid_map_[0], 1);
      TS_ASSERT_EQUALS(face_element->connected_faces_[11].node_lid_map_[3], 2);
      TS_ASSERT_EQUALS(face_element->connected_faces_[11].my_node_patch_lid_.size(), 4);
      TS_ASSERT_EQUALS(face_element->connected_faces_[11].my_node_patch_lid_[0], 1);
      TS_ASSERT_EQUALS(face_element->connected_faces_[11].my_node_patch_lid_[1], 7);
      TS_ASSERT_EQUALS(face_element->connected_faces_[11].my_node_patch_lid_[2], 8);
      TS_ASSERT_EQUALS(face_element->connected_faces_[11].my_node_patch_lid_[3], 2);

      TS_ASSERT_EQUALS(face_element->connected_faces_[13].node_lid_map_.size(), 2);
      TS_ASSERT_EQUALS(face_element->connected_faces_[13].node_lid_map_[2], 1);
      TS_ASSERT_EQUALS(face_element->connected_faces_[13].node_lid_map_[3], 0);
      TS_ASSERT_EQUALS(face_element->connected_faces_[13].my_node_patch_lid_.size(), 4);
      TS_ASSERT_EQUALS(face_element->connected_faces_[13].my_node_patch_lid_[0], 4);
      TS_ASSERT_EQUALS(face_element->connected_faces_[13].my_node_patch_lid_[1], 5);
      TS_ASSERT_EQUALS(face_element->connected_faces_[13].my_node_patch_lid_[2], 1);
      TS_ASSERT_EQUALS(face_element->connected_faces_[13].my_node_patch_lid_[3], 0);
    }

    // Calculate the averaged reference normals on the face.
    face_element->CalculateAveragedReferenceNormals(face_elements_map);
    {
      for (unsigned int i = 0; i < reference_normals.size(); i++)
        TS_ASSERT_DELTA((*face_element->GetReferenceNormals())(i), reference_normals[i], eps);
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
        TS_ASSERT_DELTA(FADUTILS::CastToDouble((*face_element->GetCurrentNormals())(i_dof)),
            current_normals[i_dof], eps);
        for (unsigned int i_der = 0; i_der < face_element->patch_dof_gid_.size(); i_der++)
        {
          TS_ASSERT_DELTA(FADUTILS::CastToDouble((*face_element->GetCurrentNormals())(i_dof).dx(
                              beam_dof_offset + i_der)),
              current_normals_derivative[i_dof][i_der], eps);
          for (unsigned int i_der_2 = 0; i_der_2 < face_element->patch_dof_gid_.size(); i_der_2++)
          {
            TS_ASSERT_DELTA(FADUTILS::CastToDouble((*face_element->GetCurrentNormals())(i_dof)
                                                       .dx(beam_dof_offset + i_der)
                                                       .dx(beam_dof_offset + i_der_2)),
                current_normals_derivative_2[i_dof][i_der][i_der_2], eps);
          }
        }
      }

      // Check an surface position on the element.
      LINALG::Matrix<3, 1, double> xi;
      xi(0) = 0.2;
      xi(1) = -0.8;
      xi(2) = 0.69;
      LINALG::Matrix<3, 1, scalar_type> r;
      GEOMETRYPAIR::EvaluateSurfacePosition<surface>(xi, face_element->GetFacePosition(), r,
          face_element->GetDrtFaceElement(), face_element->GetCurrentNormals());
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      {
        TS_ASSERT_DELTA(FADUTILS::CastToDouble(r(i_dim)), position[i_dim], eps);
        for (unsigned int i_der = 0; i_der < face_element->patch_dof_gid_.size(); i_der++)
        {
          TS_ASSERT_DELTA(FADUTILS::CastToDouble(r(i_dim).dx(beam_dof_offset + i_der)),
              position_derivative[i_dim][i_der], eps);
          for (unsigned int i_der_2 = 0; i_der_2 < face_element->patch_dof_gid_.size(); i_der_2++)
            TS_ASSERT_DELTA(FADUTILS::CastToDouble(
                                r(i_dim).dx(beam_dof_offset + i_der).dx(beam_dof_offset + i_der_2)),
                position_derivative_2[i_dim][i_der][i_der_2], eps);
        }
      }
    }
  }

 private:
  //! Pointer to the discretization object that holds the geometry for the tests.
  Teuchos::RCP<DRT::Discretization> discret_;
};

#endif
