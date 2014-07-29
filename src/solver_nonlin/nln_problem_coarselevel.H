/*----------------------------------------------------------------------------*/
/*!
\file nln_problem_coarselevel.H

\brief Coarse level interface of nonlinear solver to BACI

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

#ifndef NLN_PROBLEM_COARSELEVEL_H
#define NLN_PROBLEM_COARSELEVEL_H

/*----------------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

// baci
#include "nln_problem.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */

class Epetra_Comm;
class Epetra_MultiVector;

namespace Teuchos
{
  class ParameterList;
} // namespace Teuchos

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \brief Coarse level interface of nonlinear solver to BACI
   *
   *  This is a coarse-level-wrapper around the NLNSOL::NlnProblem class. It
   *  handles prolongation to fine and restriction to coarse level and calls the
   *  respective fine level routines (Evaluate(), ConvergenceCheck(), ...).
   *
   *  \sa NLNSOL::NlnProblem()
   *  \sa NLNSOL::FAS::NlnAMGHierarchy()
   *
   *  \author mayr.mt \date 06/2014
   */
  class NlnProblemCoarseLevel : public NlnProblem
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

    //! Setup of nonlinear problem
    virtual void Setup();

    //! Set coarse level residual corrections for FAS
    virtual void SetFHatFBar(Teuchos::RCP<Epetra_MultiVector> fhat, ///< \f$\hat{F}\f$
        Teuchos::RCP<Epetra_MultiVector> fbar ///< \f$\bar{F}\f$
        );

    //@}

    //! @name Actions on the nonlinear problem at hand
    //@{

    /*! \brief Evaluate the residual on a coarse level
     *
     *  This is an Evaluate() in a variational way. Hence, we compute \f$Rf(Px)\f$,
     *  i.e. we first prolongate the current coarse level solution to the fine
     *  level, then call Evaluate() on the fine level using the fine level
     *  nonlinear problem NLNSOL::NlnProblem, and finally restrict the residual
     *  to the current coarse level.
     *
     *  The current coarse level ID is stored internally. Thus, the Evaluate()-call
     *  looks exactly like in the fine level case.
     *
     *  \author mayr.mt \date 06/2014
     */
    virtual void Evaluate(const Epetra_MultiVector& xc, ///< current solution on coarse level
        Epetra_MultiVector& fc  ///< current residual on coarse level (to be filled)
        ) const;

    //@}

    //! @name Utilities
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::NlnProblemCoarseLevel"; }

    //@}

  protected:

  private:

    //! @name Access to members
    //@{

    //! Access to the hierarchy of multigrid levels
    virtual const NLNSOL::FAS::AMGHierarchy& Hierarchy() const;

    //@}

    //! AMG-FAS object with level hierarchy
    Teuchos::RCP<const NLNSOL::FAS::AMGHierarchy> hierarchy_;

    //! residual correction fhat on coarse level
    Teuchos::RCP<Epetra_MultiVector> fhat_;

    //! residual correction fbar on coarse level
    Teuchos::RCP<Epetra_MultiVector> fbar_;

    //! coarse level ID
    int levelid_;

  }; // class NlnProblemCoarseLevel
} // namespace NLNSOL

#endif /* NLN_PROBLEM_COARSELEVEL_H */
