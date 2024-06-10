/*----------------------------------------------------------------------*/
/*! \file

\brief validate a given .dat-file


\level 1

Validate a given 4C input file (after all preprocessing steps)

*/

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
  void ValidateInputFile(const Teuchos::RCP<Epetra_Comm> comm, const std::string datfile);

  //! Check Elements for positive Jacobian and otherwise 'rewind' them
  void ValidateMeshElementJacobians(EXODUS::Mesh& mymesh);

  //! Check for positive Jacobian for Element of distype and otherwise 'rewind' them
  void ValidateElementJacobian(
      EXODUS::Mesh& mymesh, const Core::FE::CellType distype, Teuchos::RCP<EXODUS::ElementBlock>);

  //! Check Elements of distype with full gauss integration rule for positive det at all gps and
  //! return number of negative dets
  int ValidateElementJacobian_fullgp(
      Mesh& mymesh, const Core::FE::CellType distype, Teuchos::RCP<ElementBlock> eb);

  //! Check one element for positive Jacobi determinant
  bool PositiveEle(const int& eleid, const std::vector<int>& nodes, const EXODUS::Mesh& mymesh,
      const Core::LinAlg::SerialDenseMatrix& deriv);
  int EleSaneSign(
      const std::vector<int>& nodes, const std::map<int, std::vector<double>>& nodecoords);

  //! Rewind one Element
  std::vector<int> RewindEle(std::vector<int> old_nodeids, const Core::FE::CellType distype);

}  // namespace EXODUS

FOUR_C_NAMESPACE_CLOSE

#endif
