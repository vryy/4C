/*!----------------------------------------------------------------------
\file fluid_utils_infnormscaling.H
\brief infnorm-scaling utility class for preconditioning of fluid problems

\level 3

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#ifndef FLUID_UTILS_INFNORMSCALING_H
#define FLUID_UTILS_INFNORMSCALING_H


#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>

// forward declarations
namespace LINALG
{
  class MapExtractor;
  class SparseOperator;
}  // namespace LINALG

namespace FLD
{
  namespace UTILS
  {
    class FluidInfNormScaling
    {
     public:
      //! constructor
      FluidInfNormScaling(LINALG::MapExtractor& mapextractor);

      //! destructor
      virtual ~FluidInfNormScaling(){};

      //! perform infnorm-scaling of linear system
      void ScaleSystem(Teuchos::RCP<LINALG::SparseOperator> matrix, Epetra_Vector& b);

      //! perform un-scaling of solution (and the system, just to be on the safe side)
      void UnscaleSolution(
          Teuchos::RCP<LINALG::SparseOperator> matrix, Epetra_Vector& x, Epetra_Vector& b);

     private:
      //! processor id
      const int myrank_;

      //! Extractor for splitting of velocity and pressure dofs
      LINALG::MapExtractor& velpressplitter_;

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

#endif  // #ifndef FLUID_UTILS_H
