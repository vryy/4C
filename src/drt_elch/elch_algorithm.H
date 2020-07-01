/*----------------------------------------------------------------------*/
/*! \file

\brief Basis of all ELCH algorithms

\level 2

*/
/*----------------------------------------------------------------------*/


#ifndef ELCH_ALGORITHM_H
#define ELCH_ALGORITHM_H

#include "../drt_scatra/scatra_algorithm.H"

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

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

   protected:
    /// provide information about initial field
    void PrepareTimeLoop();

    /// print scatra solver type to screen
    void PrintScaTraSolver();

    /// convergence check for natural convection solver
    bool ConvergenceCheck(int natconvitnum, int natconvitmax, double natconvittol);
  };
}  // namespace ELCH

#endif  // ELCH_ALGORITHM_H
