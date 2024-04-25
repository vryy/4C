/*----------------------------------------------------------------------*/
/*! \file

\brief Methods for spring and dashpot constraints / boundary conditions:

\level 2


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CONSTRAINT_SPRINGDASHPOT_MANAGER_HPP
#define FOUR_C_CONSTRAINT_SPRINGDASHPOT_MANAGER_HPP

#include "4C_config.hpp"

#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace CORE::LINALG
{
  class SparseMatrix;
}  // namespace CORE::LINALG

namespace IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}  // namespace IO

namespace CONSTRAINTS
{
  class SpringDashpot;

  class SpringDashpotManager
  {
   public:
    /*!
      \brief Constructor
    */
    SpringDashpotManager(Teuchos::RCP<DRT::Discretization> dis);

    /*!
     \brief Return if there are spring dashpots
    */
    bool HaveSpringDashpot() const { return havespringdashpot_; };

    //! add contribution of spring dashpot BC to residual vector and stiffness matrix
    void StiffnessAndInternalForces(Teuchos::RCP<CORE::LINALG::SparseMatrix> stiff,
        Teuchos::RCP<Epetra_Vector> fint, Teuchos::RCP<Epetra_Vector> disn,
        Teuchos::RCP<Epetra_Vector> veln, Teuchos::ParameterList parlist);

    //! update for each new time step
    void Update();

    //! output of gap, normal, and nodal stiffness
    void Output(Teuchos::RCP<IO::DiscretizationWriter> output,
        Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<Epetra_Vector> disp);

    //! output of prestressing offset for restart
    void OutputRestart(Teuchos::RCP<IO::DiscretizationWriter> output,
        Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<Epetra_Vector> disp);

    /*!
     \brief Read restart information
    */
    void ReadRestart(IO::DiscretizationReader& reader, const double& time);

    //! reset spring after having done a MULF prestressing update (mhv 12/2015)
    void ResetPrestress(Teuchos::RCP<Epetra_Vector> disold);

   private:
    Teuchos::RCP<DRT::Discretization> actdisc_;         ///< standard discretization
    std::vector<Teuchos::RCP<SpringDashpot>> springs_;  ///< all spring dashpot instances

    bool havespringdashpot_;  ///< are there any spring dashpot BCs at all?
    int n_conds_;             ///< number of spring dashpot conditions
  };                          // class
}  // namespace CONSTRAINTS
FOUR_C_NAMESPACE_CLOSE

#endif
