/*-----------------------------------------------------------*/
/*! \file

\brief


\level 3

*/
/*-----------------------------------------------------------*/
#ifndef FOUR_C_STRUCTURE_NEW_NLN_LINEARSYSTEM_SCALING_HPP
#define FOUR_C_STRUCTURE_NEW_NLN_LINEARSYSTEM_SCALING_HPP

#include "baci_config.hpp"

#include "baci_inpar_structure.hpp"

#include <Epetra_Map.h>
#include <NOX_Epetra_Scaling.H>

namespace NOX
{
  namespace Epetra
  {
    class Scaling;
  }
}  // namespace NOX

BACI_NAMESPACE_OPEN

// forward declarations
namespace STR
{
  namespace TIMINT
  {
    class BaseDataSDyn;
    class BaseDataGlobalState;
  }  // namespace TIMINT
}  // namespace STR

namespace CORE::LINALG
{
  class SparseMatrix;
}

namespace STR
{
  namespace NLN
  {
    namespace LinSystem
    {
      class StcScaling : public ::NOX::Epetra::Scaling
      {
       public:
        //! Constructor.
        StcScaling(
            const STR::TIMINT::BaseDataSDyn& DataSDyn, STR::TIMINT::BaseDataGlobalState& GState);


        //! Scales the linear system.
        void scaleLinearSystem(Epetra_LinearProblem& problem) override;

        //! Remove the scaling from the linear system.
        void unscaleLinearSystem(Epetra_LinearProblem& problem) override;

       private:
        //! stiffness matrix after scaling
        Teuchos::RCP<CORE::LINALG::SparseMatrix> stiff_scaled_;

        //! scale thickness of shells
        const enum INPAR::STR::STC_Scale stcscale_;

        //! number of layers for multilayered case
        const int stclayer_;

        //! scaling matrix for STC
        Teuchos::RCP<CORE::LINALG::SparseMatrix> stcmat_;
      };
    }  // namespace LinSystem
  }    // namespace NLN
}  // namespace STR

BACI_NAMESPACE_CLOSE

#endif