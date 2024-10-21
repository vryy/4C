#ifndef FOUR_C_PRE_EXODUS_VALIDATE_HPP
#define FOUR_C_PRE_EXODUS_VALIDATE_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_linalg_serialdensematrix.hpp"

#include <Teuchos_RCP.hpp>

#include <string>

FOUR_C_NAMESPACE_OPEN

namespace EXODUS
{
  // forward declaration
  class Mesh;
  class ElementBlock;

  //! validate a given datfile
  void validate_input_file(const Teuchos::RCP<Epetra_Comm> comm, const std::string datfile);

  //! Check Elements for positive Jacobian and otherwise 'rewind' them
  void validate_mesh_element_jacobians(EXODUS::Mesh& mymesh);

  //! Check for positive Jacobian for Element of distype and otherwise 'rewind' them
  void validate_element_jacobian(
      EXODUS::Mesh& mymesh, const Core::FE::CellType distype, EXODUS::ElementBlock&);

  //! Check Elements of distype with full gauss integration rule for positive det at all gps and
  //! return number of negative dets
  int validate_element_jacobian_fullgp(
      Mesh& mymesh, const Core::FE::CellType distype, ElementBlock& eb);

  //! Check one element for positive Jacobi determinant
  bool positive_ele(const int& eleid, const std::vector<int>& nodes, const EXODUS::Mesh& mymesh,
      const Core::LinAlg::SerialDenseMatrix& deriv);
  int ele_sane_sign(
      const std::vector<int>& nodes, const std::map<int, std::vector<double>>& nodecoords);

  //! Rewind one Element
  std::vector<int> rewind_ele(std::vector<int> old_nodeids, const Core::FE::CellType distype);

}  // namespace EXODUS

FOUR_C_NAMESPACE_CLOSE

#endif
