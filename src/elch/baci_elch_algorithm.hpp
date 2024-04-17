/*----------------------------------------------------------------------*/
/*! \file

\brief Basis of all ELCH algorithms

\level 2

*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_ELCH_ALGORITHM_HPP
#define FOUR_C_ELCH_ALGORITHM_HPP

#include "baci_config.hpp"

#include "baci_scatra_algorithm.hpp"

FOUR_C_NAMESPACE_OPEN

namespace ELCH
{
  /// ELCH algorithm base
  /*!

    Base class of ELCH algorithms. Derives from ScaTraAlgorithm.

    \author gjb
    \date 03/08
   */
  class Algorithm : public SCATRA::ScaTraAlgorithm
  {
   public:
    /// constructor
    explicit Algorithm(const Epetra_Comm& comm,     ///< communicator
        const Teuchos::ParameterList& elchcontrol,  ///< elch parameter list
        const Teuchos::ParameterList& scatradyn,    ///< scatra parameter list
        const Teuchos::ParameterList& fdyn,         ///< fluid parameter list
        const Teuchos::ParameterList& solverparams  ///< solver parameter list
    );


   protected:
    /// provide information about initial field
    void PrepareTimeLoop() override;

    /// print scatra solver type to screen
    void PrintScaTraSolver() override;

    /// convergence check for natural convection solver
    bool ConvergenceCheck(int natconvitnum, int natconvitmax, double natconvittol) override;
  };
}  // namespace ELCH

FOUR_C_NAMESPACE_CLOSE

#endif
