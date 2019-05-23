/*----------------------------------------------------------------------------*/
/*!

\brief Field problem base class for different types of BACI problems

\maintainer Matthias Mayr

\level 3
*/

/*----------------------------------------------------------------------------*/

#ifndef NLN_PROBLEM_H
#define NLN_PROBLEM_H

/*----------------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

// baci
#include "nln_problem_base.H"

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

namespace NLNSOL
{
  namespace UTILS
  {
    class DebugWriterBase;
  }  // namespace UTILS
}  // namespace NLNSOL

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
  /*! \class NlnProblem
   *  \brief Field problem base class for different types of BACI problems
   *
   *  This is a base class for field problems and acts as common entry to all
   *  derived field problem classes. Derived classes wrap BACI time integrators
   *  or NOX::Abstract::Groups that allow to evaluate the physical
   *  problem at hand.
   *
   *  Derive from NLNSOL::NlnProblem to implement the interface to a certain
   *  time integrator in BACI.
   *
   *  \author mayr.mt \date 11/2015
   */
  class NlnProblem : public NlnProblemBase
  {
   public:
    //! @name Construction / Destruction
    //@{

    //! constructor (empty)
    NlnProblem();

    //! Destructor
    virtual ~NlnProblem(){};

    //@}

    //! @name Setup
    //@{

    //! Setup of nonlinear problem
    virtual void Setup() = 0;

    /*! \brief Set actual model evaluator
     *
     *  Extract Teuchos::RCP to the actual BACI model evaluator form
     *  Teuchos::ParameterList and store it.
     *
     *  \return Bool: true if successful, false otherwise
     *
     *  \author mayr.mt \date 10/2015
     */
    virtual bool SetModelEvaluator() = 0;

    /*! \brief Set Jacobian operator
     *
     *  Extract Jacobian operator form Teuchos::ParameterList and store it.
     *
     *  \return Bool: true if successful, false otherwise
     *
     *  \author mayr.mt \date 10/2015
     */
    virtual bool SetJacobianOperator() = 0;

    /*! \brief Set coarse level residual modifications for FAS
     *
     *  \note Only needed for coarse levels. Implemented here only for purpose
     *  of inheritance.
     *
     *  \sa NLNSOL::NlnProblemCoarseLevel::SetupResidualModification
     */
    virtual void SetupResidualModification(
        Teuchos::RCP<const Epetra_MultiVector> xbar,  ///< \f$\bar{x}\f$
        Teuchos::RCP<const Epetra_MultiVector> fbar   ///< \f$\bar{F}\f$
    )
    {
      dserror(
          "Not implemented in base class. This routine is needed only in "
          "case of a coarse level problem within FAS-AMG.");
    }

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
        ) const = 0;

    /*! \brief Compute nonlinear residual \f$F(x)\f$
     *
     *  \author mayr.mt \date 10/2015
     */
    virtual Teuchos::RCP<const Epetra_MultiVector> ComputeF(
        const Epetra_MultiVector& x  ///< current solution
        ) const = 0;

    /*! \brief Evaluate Jacobian (linearization of residual)
     *
     *  Updating the solution variables in the marching time integrator has
     *  already been done during computeF().
     *
     *  \author mayr.mt \date 06/2014
     */
    virtual void ComputeJacobian() const = 0;

    //@}

    //! @name Get attributes
    //@{

    /*! \brief Level ID of level where this nonlinear problem lives on
     *
     *  Since this is the fine level problem, the level ID is always 0.
     *
     *  \author mayr.mt \date 03/2015
     */
    virtual int LevelID() const { return 0; }

    //@}

    //! Debug output of data
    //@{

    /*! \brief Write Teuchos::RCP<Epetra_MultiVector> to output
     *
     *  \author mayr.mt \date 06/2015
     */
    virtual void WriteVector(
        Teuchos::RCP<const Epetra_MultiVector> vec,  ///< vector with result data (to be outputted)
        const std::string& description,              ///< description of output quantity
        const IO::VectorType vt = IO::dofvector      ///< vector type
        ) const;

    /*! \brief Wrapper to write plain Epetra_MultiVector to output
     *
     *  \author mayr.mt \date 06/2015
     */
    virtual void WriteVector(
        const Epetra_MultiVector& vec,           ///< vector with result data (to be outputted)
        const std::string& description,          ///< description of output quantity
        const IO::VectorType vt = IO::dofvector  ///< vector type
        ) const;

    //@}

    //! @name Access to underlying NOX group elements
    //@{

    //! get Jacobian as Epetra_Operator
    virtual Teuchos::RCP<Epetra_Operator> GetJacobianOperator() const = 0;

    //@}

    //! @name Utilities
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::NlnProblem"; }

    //@}

   protected:
   private:
  };  // class NlnProblem
}  // namespace NLNSOL

#endif /* NLN_PROBLEM_H */
