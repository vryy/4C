/*----------------------------------------------------------------------*/
/*!
\file sti_partitioned.H

\brief partitioned coupling algorithm for scatra-thermo interaction

\level 2

\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
*/
/*----------------------------------------------------------------------*/
#ifndef STI_PARTITIONED_H
#define STI_PARTITIONED_H

#include "sti_algorithm.H"

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
    void Solve();

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
#endif
