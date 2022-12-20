/*-----------------------------------------------------------*/
/*! \file

\brief Evaluate calls for magno-hydrodynamics


\level 2

*/
/*-----------------------------------------------------------*/


#ifndef FLUID_MHD_EVALUATE_H
#define FLUID_MHD_EVALUATE_H

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>


// forward declarations
namespace DRT
{
  class Discretization;
}

namespace LINALG
{
  class SparseMatrix;
}

namespace FLD
{
  class FluidMHDEvaluate
  {
   public:
    /*!

    \brief Standard Constructor (public)

    o trigger numbering of dofs (including periodic boundaries)

    o create empty system matrix and right hand side

    o allocate distributed solution vectors

    o generate maps to split velocity and pressure dofs

    o set alpha_F and alpha_M

    o calculate and set gamma parameter for generalised alpha method
      according to the second order necessary condition

    */
    FluidMHDEvaluate(Teuchos::RCP<DRT::Discretization> actdis);

    /*!
   \brief Destructor

   */
    virtual ~FluidMHDEvaluate();


    void BoundaryElementLoop(Teuchos::ParameterList& mhdbcparams,
        Teuchos::RCP<Epetra_Vector> velaf_, Teuchos::RCP<Epetra_Vector> velnp_,
        Teuchos::RCP<Epetra_Vector> residual_, Teuchos::RCP<LINALG::SparseMatrix> sysmat_);


   private:
    //! the parent discretization
    Teuchos::RCP<DRT::Discretization> pdiscret_;
    //! the boundary discretisation
    Teuchos::RCP<DRT::Discretization> bnd_discret_;
    //! dofrowmap for
    Teuchos::RCP<Epetra_Map> subdofrowmap_;
    //! (small) system FEmatrix for boundary discretisation
    Teuchos::RCP<LINALG::SparseMatrix> bndmat_;
  };

}  // namespace FLD
#endif  // not FLUID_MHD_EVALUATE
