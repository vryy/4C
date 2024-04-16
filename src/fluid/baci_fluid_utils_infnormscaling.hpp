/*-----------------------------------------------------------*/
/*! \file

\brief infnorm-scaling utility class for preconditioning of fluid problems


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_UTILS_INFNORMSCALING_HPP
#define FOUR_C_FLUID_UTILS_INFNORMSCALING_HPP


#include "baci_config.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

// forward declarations
namespace CORE::LINALG
{
  class MapExtractor;
  class SparseOperator;
}  // namespace CORE::LINALG

namespace FLD
{
  namespace UTILS
  {
    class FluidInfNormScaling
    {
     public:
      //! constructor
      FluidInfNormScaling(CORE::LINALG::MapExtractor& mapextractor);

      //! destructor
      virtual ~FluidInfNormScaling() = default;

      //! perform infnorm-scaling of linear system
      void ScaleSystem(Teuchos::RCP<CORE::LINALG::SparseOperator> matrix, Epetra_Vector& b);

      //! perform un-scaling of solution (and the system, just to be on the safe side)
      void UnscaleSolution(
          Teuchos::RCP<CORE::LINALG::SparseOperator> matrix, Epetra_Vector& x, Epetra_Vector& b);

     private:
      //! processor id
      const int myrank_;

      //! Extractor for splitting of velocity and pressure dofs
      CORE::LINALG::MapExtractor& velpressplitter_;

      Teuchos::RCP<Epetra_Vector> srowsum_;
      Teuchos::RCP<Epetra_Vector> scolsum_;
      Teuchos::RCP<Epetra_Vector> prowsum_;
      Teuchos::RCP<Epetra_Vector> pcolsum_;

      // flags
      const bool leftscale_momentum_;
      const bool leftscale_continuity_;

    };  // class FluidInfNormScaling

  }  // namespace UTILS
}  // namespace FLD

BACI_NAMESPACE_CLOSE

#endif
