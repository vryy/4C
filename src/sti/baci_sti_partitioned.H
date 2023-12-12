/*----------------------------------------------------------------------*/
/*! \file

\brief partitioned coupling algorithm for scatra-thermo interaction

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_STI_PARTITIONED_H
#define BACI_STI_PARTITIONED_H

#include "baci_inpar_sti.H"
#include "baci_sti_algorithm.H"

BACI_NAMESPACE_OPEN

namespace STI
{
  //! partitioned coupling algorithm for scatra-thermo interaction
  class Partitioned : public Algorithm
  {
   public:
    //! constructor
    explicit Partitioned(const Epetra_Comm& comm,  //! communicator
        const Teuchos::ParameterList& stidyn,      //! parameter list for scatra-thermo interaction
        const Teuchos::ParameterList&
            scatradyn,  //! scalar transport parameter list for scatra and thermo fields
        const Teuchos::ParameterList&
            solverparams_scatra,  //! solver parameter list for scatra field
        const Teuchos::ParameterList&
            solverparams_thermo  //! solver parameter list for thermo field
    );

   private:
    //! convergence check for iterative staggered TSI solver
    bool ExitOuterCoupling() const;

    //! evaluate time step using outer coupling iteration
    void Solve() override;

    //! evaluate time step using one-way coupling iteration
    void SolveOneWay() const;

    //! evaluate time step using two-way coupling iteration
    void SolveTwoWay();

    //! type of coupling between scatra and thermo fields
    const INPAR::STI::CouplingType couplingtype_;

    //! maximum value of Aitken relaxation parameter
    const double omegamax_;
  };  // class Partitioned : public Algorithm
}  // namespace STI
BACI_NAMESPACE_CLOSE

#endif
