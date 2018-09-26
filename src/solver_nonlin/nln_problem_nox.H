/*----------------------------------------------------------------------------*/
/*!
\file nln_problem_nox.H

\brief Field problem using NOX interface to BACI

\level 3

\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------------*/

#ifndef NLN_PROBLEM_NOX_H
#define NLN_PROBLEM_NOX_H

/*----------------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_VerboseObject.hpp>

// baci
#include "nln_problem.H"

#include "../drt_io/io.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */

class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Operator;

namespace LINALG
{
  class SparseOperator;
}  // namespace LINALG

namespace NOX
{
  namespace Abstract
  {
    class Group;
  }  // namespace Abstract
}  // namespace NOX

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \class NlnProblemNox
   *  \brief Field problem using NOX interface to BACI
   *
   *  Interface class between the nonlinear solver and BACI. We utilize the
   *  NOX::Abstract::Group() to access the time integration within BACI for
   *  Evaluate() / computeF().
   *
   *  NLNSOL::NlnProblem acts on the DofRowMap() of the (fine) level, only. For
   *  coarse level representations, we use NLNSOL::NlnProblemCoarseLevel().
   *
   *  \author mayr.mt \date 10/2015
   */
  class NlnProblemNox : public NlnProblem
  {
   public:
    //! @name Construction / Destruction
    //@{

    //! constructor (empty)
    NlnProblemNox();

    //! Destructor
    virtual ~NlnProblemNox(){};

    //@}

    //! @name Setup
    //@{

    //! Setup of nonlinear problem
    virtual void Setup();

    /*! \brief Set actual model evaluator
     *
     *  Extract Teuchos::RCP to the actual BACI model evaluator form
     *  Teuchos::ParameterList and store it.
     *
     *  \return Bool: true if successful, false otherwise
     *
     *  \author mayr.mt \date 10/2015
     */
    virtual bool SetModelEvaluator();

    /*! \brief Set Jacobian operator
     *
     *  Extract Jacobian operator form Teuchos::ParameterList and store it.
     *
     *  \return Bool: true if successful, false otherwise
     *
     *  \author mayr.mt \date 10/2015
     */
    virtual bool SetJacobianOperator();

    //@}

    //! @name Actions on the nonlinear problem at hand
    //@{

    /*! \brief Compute nonlinear residual \f$F(x)\f$
     *
     *  Set the most recent solution 'x' in the NOX group #noxgrp_. Afterwards,
     *  the marching time integrator is asked to evaluate the residual via the
     *  computeF()-call.
     *
     *  Updating the solution variables in the marching time integrator is done
     *  during computeF().
     *
     *  Return the most recent residual in 'f' as a descent direction.
     *
     *  \note Solvers implemented in the NLNSOL framework require a residual
     *  that points into descent direction.
     *
     *  \note One can define the 'mechanical' residual as
     *  \f$F(x) = f_{int} - f_{ext}\f$ or as \f$F(x) = f_{ext} - f_{int}\f$.
     *  This choice is arbitrary, since both formulations result in exactly the
     *  same solution \f$x\f$. In many Finite Element Codes as well as in baci,
     *  often the first version is utilized which produces an ascent direction.
     *  For a descent direction, one has to use the latter formula.
     *  To be consistent with mathematical literature, the residual needs
     *  to be a descent direction. Hence, we need to switch the sign of the
     *  'mechanical' baci-residual.
     *
     *  \author mayr.mt \date 09/2014
     */
    virtual void ComputeF(const Epetra_MultiVector& x,  ///< current solution
        Epetra_MultiVector& f                           ///< current residual (to be filled)
        ) const;

    /*! \brief Compute nonlinear residual \f$F(x)\f$
     *
     *  \author mayr.mt \date 10/2015
     */
    virtual Teuchos::RCP<const Epetra_MultiVector> ComputeF(
        const Epetra_MultiVector& x  ///< current solution
        ) const;

    /*! \brief Evaluate Jacobian (linearization of residual)
     *
     *  Updating the solution variables in the marching time integrator has
     *  already been done during computeF().
     *
     *  \author mayr.mt \date 06/2014
     */
    virtual void ComputeJacobian() const;

    //@}

    //! @name Access to underlying NOX group elements
    //@{

    //! get Jacobian as Epetra_Operator
    virtual Teuchos::RCP<Epetra_Operator> GetJacobianOperator() const;

    //@}

    //! @name Utilities
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::NlnProblemNox"; }

    //@}

   protected:
   private:
    //! NOX group
    virtual NOX::Abstract::Group& NOXGroup() const;

    //! NOX group for access to time integrator via call-back routines
    Teuchos::RCP<NOX::Abstract::Group> noxgrp_;

    /*! \brief Pointer to Jacobian matrix of underlying problem
     *
     *  This needs to be a LINALG::SparseOperator since this is the base class
     *  for all matrix types in BACI. We can always extract an Epetra_Operator
     *  when we want to call a Trilinos package/solver.
     *
     *  \warning Access to time integration WIHTOUT using the NOX group
     *  #noxgrp_. Seems to be inconsistent since all other accesses to time
     *  integration are done via the NOX group #noxgrp_. However, usage of NOX
     *  group might go away at some point.
     */
    Teuchos::RCP<LINALG::SparseOperator> jac_;

  };  // class NlnProblemNox
}  // namespace NLNSOL

#endif /* NLN_PROBLEM_NOXH */
