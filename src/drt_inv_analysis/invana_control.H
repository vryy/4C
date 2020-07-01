/*----------------------------------------------------------------------------*/
/*! \file
\brief Control object to handle solution of the inverse analysis

\level 3

*/

/*----------------------------------------------------------------------------*/

#ifndef INVANA_CONTROL_H
#define INVANA_CONTROL_H

/*----------------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace INVANA
{
  class InvanaBase;
  class OptimizerBase;
}  // namespace INVANA

namespace Teuchos
{
  class ParameterList;
}

namespace DRT
{
  class ResultTest;
}

namespace IO
{
  class DiscretizationReader;
}

class Epetra_MultiVector;

/*----------------------------------------------------------------------------*/
namespace INVANA
{
  /*! \class InvanaControl
   *  \brief Class to control an inverse problem and its solution
   *
   *  TODO: Documentation
   *
   *  \author kehl \date 12/2015
   */
  class InvanaControl
  {
   public:
    //! @name Construction/destruction
    //@{

    //! Constructor (empty)
    InvanaControl();

    //! Destructor
    virtual ~InvanaControl() {}

    //@}

    //! @name Initialization
    //@{

    /*! \brief Initialization
     *
     *  Initialization of the nonlinear inverse problem
     *  and the algorithm used to solve it
     *
     *  \author kehl \date 12/2015
     */
    void Init(const Teuchos::ParameterList& invp);

    //@}

    //! @name Solution
    //@{

    /*! \brief Start the solution
     *
     *  Call the solver object to solve the inverse problem
     *
     *  \author kehl \date 12/2015
     */
    int Solve(int restart);

    //! @name Access to the solution
    //@{

    //! Get the gradient
    const Epetra_MultiVector& GetGradient() { return *f_; }

    //! Get functional value
    const double& GetValue() { return val_; }
    //@}

    //! add invana specific result test to the global problem field tests
    Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

    //!@name Public Access wrapper
    //@{

    //! the inverse problem pointer
    Teuchos::RCP<INVANA::InvanaBase> InverseProblem() { return invprob_; }
    //@}

   private:
    //!@name Internal Access wrapper
    //@{

    //! The INVANA optimizer
    Teuchos::RCP<INVANA::OptimizerBase> InvanaOpti() { return invanaopt_; }

    //@}

    //! @name What to be solved and how
    //{@

    //! the inverse problem
    Teuchos::RCP<INVANA::InvanaBase> invprob_;

    //! an inverse analysis internal solver
    Teuchos::RCP<INVANA::OptimizerBase> invanaopt_;

    //@}

    //! @name Input
    //@{

    //! input reader
    Teuchos::RCP<IO::DiscretizationReader> input_;

    //@}

    //! @name Solution of the inverse problem
    //@{

    //! solution using INVANA
    void invanasolve(int restart);

    //@}

    //! @name solution variables and storage
    //{@

    //! the current solution
    Teuchos::RCP<Epetra_MultiVector> x_;

    //! the resulting nonlinear function
    Teuchos::RCP<Epetra_MultiVector> f_;

    //! the functional value
    double val_;
    //@}

  };  // class InvanaControl
}  // namespace INVANA

#endif /* INVANA_CONTROL_H */
