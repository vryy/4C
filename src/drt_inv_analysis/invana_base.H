/*----------------------------------------------------------------------*/
/*! \file
\brief inverse analysis base class

\level 3

\maintainer Sebastian Brandstaeter
*/
/*----------------------------------------------------------------------*/
#ifndef INVANA_BASE_H
#define INVANA_BASE_H

#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"

// forward declarations
namespace DRT
{
  class ResultTest;
  class Discretization;
}  // namespace DRT

namespace IO
{
  class OutputControl;
  class InputControl;
}  // namespace IO

namespace INVANA
{
  class ObjectiveFunct;
  class MatParManager;
  class RegularizationBase;
  class InitialGuess;

  /*!
  \brief A base class for the inverse analysis
  \author kehl

  This class is a base class to setup an optimzation problem. It's main purpose is to hold
  all base class pointers to all the functional classes used to make up an optimization
  problem, i.e.:\n
  - the discretization (more for convenience)\n
  - an objective function\n
  - a parameterization of the material parameters to be estimated\n
  - a regularization to smooth the objective function\n
 */

  class InvanaBase
  {
   public:
    //! standard constructor
    InvanaBase();

    //! destructor
    virtual ~InvanaBase(){};

    /*! Initialize class members.
     *
     *  Setup of the neccessary elements to be found in every
     *  standard inverse analysis procedure
     */
    void Init(Teuchos::RCP<DRT::Discretization> discret,  ///< the underlying discretization
        Teuchos::RCP<INVANA::ObjectiveFunct> objfunct,    ///< an objective function
        Teuchos::RCP<INVANA::MatParManager> matman,       ///< the material parameterization
        Teuchos::RCP<INVANA::RegularizationBase> regman,  ///< a regularization strategy
        Teuchos::RCP<INVANA::InitialGuess> initguess      ///< initial mean and covariance
    );

    /*! Setup IO
     *
     * Writing simulation output in baci for every field is handled by the
     * discretization (at least for structure i have no idea about fluid).
     * In the inverse analysis we actually want no output or at least dont
     * write forward problem output into the invana problem control file. So
     * I create a new control file bla_forward and pass it over to the the
     * discretization which can then use it in the field output.
     *
     * Also input for the forward problem (in case of restarts should not be
     * read from the invana problem inputfile but from a forward problem
     * inputfile.
     *
     * I keep the BINIO flag for the forward problem output as specified in the
     * IO section.
     *
     * This is here in the base class so the the INVANA::OptimizerBase can
     * rely on this behavior.
     *
     * \author kehl \date 06/16
     */
    void SetupIO();

    //! Setup of the derived classes.
    virtual void Setup() = 0;

    /*! \brief Value and gradient of the functional
      \param In
      sol - the current solution of the optimization problem
      \param Out
      val - the value of the objective function
      \param Out
      optprob - the gradient of the objective function wrt to the optimization parameters

      returns 0 if successfull (conditioned on the erros tested in
      the forward problem evaluator)
    */
    virtual int Evaluate(const Epetra_MultiVector& sol,  ///< input to the problem
        double* val,                                     ///> value of the functional to be filled
        Teuchos::RCP<Epetra_MultiVector> gradient        ///> gradient of the functional
        ) = 0;

    //! @name General access
    //@{

    //! a pointer to the discretization
    Teuchos::RCP<DRT::Discretization> Discret() { return discret_; }

    //! a pointer to the objective function
    Teuchos::RCP<INVANA::ObjectiveFunct> ObjectiveFunct() { return objfunct_; }

    //! a pointer to the material parameter manager
    Teuchos::RCP<INVANA::MatParManager> Matman() { return matman_; }

    //! a pointer to the regularization of the objective function
    Teuchos::RCP<INVANA::RegularizationBase> Regman() { return regman_; }

    //! a pointer to the initial guess object
    Teuchos::RCP<INVANA::InitialGuess> InitialGuess() { return initguess_; }

    //! a reference to the programs communicator
    const Epetra_Comm& Comm();

    //! the unique layout of the optimization parameters
    Teuchos::RCP<Epetra_Map> VectorRowLayout();

    //! inputcontrolfile
    Teuchos::RCP<IO::InputControl>& InputControl();

    //! Restart of the forward problem
    int FPRestart() { return fprestart_; }
    //@}

    //! @name Special acess
    //@{

    /*! \brief Get primal variable access
     *
     * For some optimization problems (at least the ones where pdes are involved)
     * one might be interested in the solution of the pde (the primal variable).
     * E.g., when the solution of the pde under the optimized conditions is of
     * interest, or when prediction based on the optimization result is desired.
     *
     * \author kehl \date 09/2016
     */
    virtual const Epetra_MultiVector& GetPrimalVariables() const = 0;
    //@}

   private:
    //! a reference to the discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! objective function to be optimized
    Teuchos::RCP<INVANA::ObjectiveFunct> objfunct_;

    //! handler for calls to the parametrization
    Teuchos::RCP<INVANA::MatParManager> matman_;

    //! handler for calls to the regularization
    Teuchos::RCP<INVANA::RegularizationBase> regman_;

    //! initial guess
    Teuchos::RCP<INVANA::InitialGuess> initguess_;

    //! input control
    Teuchos::RCP<IO::InputControl> inputcontrol_;

    //! forward problem restart step
    int fprestart_;

    //! initialization bool
    bool isinit_;
  };
}  // namespace INVANA

#endif /*INVANA_BASE_H_*/
