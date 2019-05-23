/*----------------------------------------------------------------------------*/
/*!

\brief Coarse level interface of nonlinear solver to BACI

\maintainer Matthias Mayr

\level 3
*/

/*----------------------------------------------------------------------------*/

#ifndef NLN_PROBLEM_COARSELEVEL_H
#define NLN_PROBLEM_COARSELEVEL_H

/*----------------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

// baci
#include "nln_problem_base.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */

class Epetra_MultiVector;

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

// baci
namespace NLNSOL
{
  class NlnProblem;
}  // namespace NLNSOL

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \class NlnProblemCoarseLevel
   *  \brief Coarse level interface of nonlinear solver to BACI
   *
   *  This is a coarse-level-wrapper around the NLNSOL::NlnProblem class. It
   *  handles access to the nonlinear problem on coarse levels of a FAS-AMG
   *  procedure.
   *
   *  Access to respective fine level routines (ComputeF(), ComputeJacobian(),
   *  ConvergenceCheck(), ...) is provided where necessary level transfer
   *  operations are performed.
   *
   *  \sa NlnProblem
   *  \sa FAS::NlnAMGHierarchy
   *
   *  \author mayr.mt \date 06/2014
   */
  class NlnProblemCoarseLevel : public NlnProblemBase
  {
   public:
    //! @name Construction / Destruction
    //@{

    //! constructor (empty)
    NlnProblemCoarseLevel();

    //! Destructor
    virtual ~NlnProblemCoarseLevel(){};

    //@}

    //! @name Setup
    //@{

    /*! \brief Setup of nonlinear problem
     *
     *  On the coarse level, we just need to store the access to the
     *  NLNSOL::FAS::AMGHierarchy and the level ID #levelid_ of the current
     *  level. All other Setup stuff is done by NLNSOL::NlnProblem.
     *
     *  It is desirable to initialize the residual corrections #fhat_ and #fbar_
     *  during Setup. However, the dof row map of the specific level is not
     *  accessible at that time, yet, so we can't do so. This has to be done
     *  later by calling SetupResidualModification() with zeroed vectors.
     *
     *  \sa NLNSOL::NlnProblem::Setup
     *  \sa SetupResidualModification
     */
    virtual void Setup();

    /*! \brief Set actual model evaluator
     *
     *  Extract Teuchos::RCP to the actual BACI model evaluator form
     *  Teuchos::ParameterList store it.
     *
     *  \note For the coarse level problem, the model evaluator is of type
     *  NLNSOL::NlnProblem.
     *
     *  \return Bool: true if successful, false otherwise
     *
     *  \sa NLNSOL::NlnProblem()
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
    virtual bool SetJacobianOperator()
    { /* nothing to do */
      return true;
    }

    /*! \brief Set coarse level residual modifications for FAS
     *
     *  Two coarse level residual modifications are needed for FAS:
     *  <ul>
     *  <li> \f$\hat{F}_c = R\cdot F_f(x_f)\f$: restriction of the current fine
     *       level residual </li>
     *  <li> \f$\bar{F}_c = R\cdot F_f(P\bar{x})\f$: variationally evaluated
     *       residual of the restricted fine level solution
     *       \f$\bar{x} = R\cdot x_f\f$ </li>
     *  </ul>
     *
     *  They are computed and stored for later use in ComputeF().
     */
    virtual void SetupResidualModification(
        Teuchos::RCP<const Epetra_MultiVector> xbar,  ///< \f$\bar{x}\f$
        Teuchos::RCP<const Epetra_MultiVector> fbar   ///< \f$\bar{F}\f$
    );

    //@}

    //! @name Actions on the nonlinear problem at hand
    //@{

    /*! \brief Compute the modified residual on a coarse level
     *
     *  First, we evaluate the plain residual \f$F\f$ in a variational way
     *  using ComputPlainF(). Due to FAS, we need to apply coarse level residual
     *  corrections \f$\hat{F}\f$ and \f$\bar{F}\f$ afterwards. The modified
     *  residual finally reads
     *  \f[
     *    F - \hat{F} + \bar{F} \rightarrow 0
     *  \f]
     *
     *  The current coarse level ID is stored internally. Thus, the
     *  ComputeF()-call looks exactly like in the fine level case.
     *
     *  \author mayr.mt \date 06/2014
     */
    virtual void ComputeF(const Epetra_MultiVector& xc,  ///< current solution on coarse level
        Epetra_MultiVector& fc  ///< current residual on coarse level (to be filled)
        ) const;

    /*! \brief Compute the modified residual on a coarse level
     *
     *  \author mayr.mt \date 06/2014
     */
    virtual Teuchos::RCP<const Epetra_MultiVector> ComputeF(
        const Epetra_MultiVector& xc  ///< current solution on coarse level
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

    //! @name Attribute access routines
    //@{

    /*! \brief Level ID of level where this nonlinear problem lives on
     *
     *  Since this is the coarse level problem, the level ID is #levelid_ > 0.
     *
     *  \return ID of current level
     *
     *  \author mayr.mt \date 03/2015
     */
    virtual int LevelID() const { return levelid_; }

    //@}

    //! @name Access to underlying NOX group elements
    //@{

    //! get Jacobian as Epetra_Operator
    virtual Teuchos::RCP<Epetra_Operator> GetJacobianOperator() const;

    //@}

    //! Debug output of data
    //@{

    /*! \brief Write Teuchos::RCP<Epetra_MultiVector> to output
     *
     *  The debug output writer lives on the fine level discretization, only.
     *  Thus, we first prolongate the coarse level vector \p vec to the fine
     *  level and then output it using the fine level discretization writer.
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

    //! @name Utilities
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::NlnProblemCoarseLevel"; }

    //@}

   protected:
   private:
    /*! \brief Compute the plain residual on a coarse level
     *
     *  This is an ComputeF() in a variational way. Hence, we compute
     *  \f$RF(Px)\f$, i.e. we first prolongate the current coarse level solution
     *  to the fine level, then call ComputeF() on the fine level using the fine
     *  level nonlinear problem NLNSOL::NlnProblem, and restrict the
     *  residual to the current coarse level, afterwards.
     *
     *  \return residual on coarse level
     *
     *  \author mayr.mt \date 10/2015
     */
    virtual Teuchos::RCP<const Epetra_MultiVector> ComputePlainF(
        const Epetra_MultiVector& xc  ///< current solution on coarse level
        ) const;

    //! @name Access to members
    //@{

    //! Access to the hierarchy of multigrid levels
    virtual const NLNSOL::FAS::AMGHierarchy& Hierarchy() const;

    //@}

    //! Field problem to be evaluated
    Teuchos::RCP<NLNSOL::NlnProblem> nlnproblem_;

    //! AMG-FAS object with level hierarchy
    Teuchos::RCP<const NLNSOL::FAS::AMGHierarchy> hierarchy_;

    /*! \brief Residual correction \f$\hat{F}_c\f$ on coarse level
     *
     *  \f$\bar{F}_c = R\cdot F_f(P\bar{x})\f$: variationally evaluated
     *  residual of the restricted fine level solution
     *  \f$\bar{x} = R\cdot x_f\f$
     */
    Teuchos::RCP<const Epetra_MultiVector> fhat_;

    /*! \brief Residual correction \f$\bar{F}_c\f$ on coarse level
     *
     *  \f$\hat{F}_c = R\cdot F_f(x_f)\f$: restriction of the current fine
     *  level residual
     */
    Teuchos::RCP<const Epetra_MultiVector> fbar_;

    //! coarse level ID
    int levelid_;

  };  // class NlnProblemCoarseLevel
}  // namespace NLNSOL

#endif /* NLN_PROBLEM_COARSELEVEL_H */
