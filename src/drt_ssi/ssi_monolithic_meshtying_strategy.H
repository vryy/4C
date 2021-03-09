/*----------------------------------------------------------------------*/
/*! \file
\brief Mesh tying strategy for monolithic SSI

\level 2

 */
/*----------------------------------------------------------------------*/
#ifndef SSI_MONOLITHIC_MESHTYING_STRATEGY_H
#define SSI_MONOLITHIC_MESHTYING_STRATEGY_H

#include "Teuchos_RCPDecl.hpp"
#include "vector"

#include "../linalg/linalg_sparseoperator.H"

namespace LINALG
{
  class SparseMatrix;
  class SparseOperator;
}  // namespace LINALG

namespace SSI
{
  class SSIMono;

  class MeshtyingStrategyBase
  {
   public:
    //! constructor
    explicit MeshtyingStrategyBase(const SSI::SSIMono& ssi_mono);

    /*!
     * @brief apply mesh tying to structure right hand side vector
     *
     * @param[in] structure_rhs  structure right hand side vector without mesh tying contributions
     * @return structure right hand side vector including mesh tying contributions
     */
    Epetra_Vector ApplyMeshtyingToStructureRHS(Teuchos::RCP<const Epetra_Vector> structure_rhs);

   private:
    //! monolithic algorithm for scalar-structure interaction
    const SSI::SSIMono& SSIMono() const { return ssi_mono_; }

    //! indicates, if meshtying includes three entities
    bool Meshtying3DomainIntersection() const { return meshtying_3_domain_intersection_; }

    //! indicates, if meshtying includes three entities
    const bool meshtying_3_domain_intersection_;

    //! monolithic algorithm for scalar-structure interaction
    const SSI::SSIMono& ssi_mono_;
  };

  //! build specific meshtying strategy
  Teuchos::RCP<SSI::MeshtyingStrategyBase> BuildMeshtyingStrategy(const SSI::SSIMono& ssi_mono);
}  // namespace SSI
#endif  // #ifndef SSI_MONOLITHIC_MESHTYING_STRATEGY_H