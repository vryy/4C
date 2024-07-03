/*-----------------------------------------------------------*/
/*! \file

\brief


\level 3

*/
/*-----------------------------------------------------------*/
#ifndef FOUR_C_STRUCTURE_NEW_NLN_LINEARSYSTEM_SCALING_HPP
#define FOUR_C_STRUCTURE_NEW_NLN_LINEARSYSTEM_SCALING_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"

#include <Epetra_Map.h>
#include <NOX_Epetra_Scaling.H>

namespace NOX
{
  namespace Epetra
  {
    class Scaling;
  }
}  // namespace NOX

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Solid
{
  namespace TimeInt
  {
    class BaseDataSDyn;
    class BaseDataGlobalState;
  }  // namespace TimeInt
}  // namespace Solid

namespace Core::LinAlg
{
  class SparseMatrix;
}

namespace Solid
{
  namespace Nln
  {
    namespace LinSystem
    {
      class StcScaling : public ::NOX::Epetra::Scaling
      {
       public:
        //! Constructor.
        StcScaling(const Solid::TimeInt::BaseDataSDyn& DataSDyn,
            Solid::TimeInt::BaseDataGlobalState& GState);

        //! Scales the linear system.
        void scaleLinearSystem(Epetra_LinearProblem& problem) override;

        //! Remove the scaling from the linear system.
        void unscaleLinearSystem(Epetra_LinearProblem& problem) override;

       private:
        //! stiffness matrix after scaling
        Teuchos::RCP<Core::LinAlg::SparseMatrix> stiff_scaled_;

        //! scale thickness of shells
        const enum Inpar::Solid::StcScale stcscale_;

        //! number of layers for multilayered case
        const int stclayer_;

        //! scaling matrix for STC
        Teuchos::RCP<Core::LinAlg::SparseMatrix> stcmat_;
      };
    }  // namespace LinSystem
  }    // namespace Nln
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif