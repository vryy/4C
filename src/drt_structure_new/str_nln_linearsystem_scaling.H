/*-----------------------------------------------------------*/
/*! \file

\brief

\maintainer Amadeus Gebauer

\level 3

*/
/*-----------------------------------------------------------*/

#include <NOX_Epetra_Scaling.H>
#include "Epetra_Map.h"
#include "../drt_inpar/inpar_structure.H"

// forward declarations
namespace STR
{
  namespace TIMINT
  {
    class BaseDataSDyn;
    class BaseDataGlobalState;
  }  // namespace TIMINT
}  // namespace STR

namespace LINALG
{
  class SparseMatrix;
}

namespace NOX
{
  namespace Epetra
  {
    class Scaling;
  }
}  // namespace NOX

namespace STR
{
  namespace NLN
  {
    namespace LinSystem
    {
      class StcScaling : public NOX::Epetra::Scaling
      {
       public:
        //! Constructor.
        StcScaling(STR::TIMINT::BaseDataSDyn DataSDyn, STR::TIMINT::BaseDataGlobalState GState);

        //! Virtual destructor
        virtual ~StcScaling(){};

        //! Scales the linear system.
        void scaleLinearSystem(Epetra_LinearProblem& problem);

        //! Remove the scaling from the linear system.
        void unscaleLinearSystem(Epetra_LinearProblem& problem);

       private:
        //! stiffness matrix after scaling
        Teuchos::RCP<LINALG::SparseMatrix> stiff_scaled_;

        //! scale thickness of shells
        const enum INPAR::STR::STC_Scale stcscale_;

        //! number of layers for multilayered case
        const int stclayer_;

        //! scaling matrix for STC
        Teuchos::RCP<LINALG::SparseMatrix> stcmat_;
      };
    }  // namespace LinSystem
  }    // namespace NLN
}  // namespace STR
