/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_fas.H

\brief AMG-FAS as nonlinear preconditioner

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

#ifndef NLN_OPERATOR_FAS_H
#define NLN_OPERATOR_FAS_H

/*----------------------------------------------------------------------------*/
/* headers */

// standard
#include <vector>
#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>

// baci
#include "nln_operator_base.H"

#include "../drt_inpar/inpar_solver_nonlin.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */

class Epetra_Comm;
class Epetra_MultiVector;

namespace NLNSOL
{
  namespace FAS
  {
    class AMGHierarchy;
    class NlnLevel;
  }  // namespace FAS
}  // namespace NLNSOL

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

/*----------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \class NlnOperatorFas
   *  \brief Nonlinear Algebraic Multigrid with Full Approximation Scheme.
   *
   *  The FAS operator holds a NLNSOL::FAS::AMGHierarchy which stores all the
   *  multigrid levels. The algorithm that loops over all levels as solver or
   *  preconditioner is implemented in ApplyInverse(). Currently, just a simple
   *  V-cycle() is available.
   *
   *  To evaluate the problem in the time integrator,
   *  Hierarchy()->NlnLevel(level)->NlnProblem()->Evaluate(x,f) is used.
   *
   *  All level smoothers and level specific stuff is stored in the
   *  NLNSOL::FAS::AMGHierarchy in a collection of NLNSOL::FAS::NlnLevel
   *  objects.
   *
   *  <h3>References:</h3>
   *  <ul>
   *  <li>[1] Briggs WL, Van Henson E and McCormick SF: A Multigrid Tutorial,
   *          SIAM, ed. 2, 2000 </li>
   *  <li>[2] Gee MW and Tuminaro RS: Nonlinear algebraic multigrid for
   *          constrained solid mechanics problems, Sandia National
   *          Laboratories, Technical Report Nr. SAND2006-2256, 2006 </li>
   *  <li>[3] Wiesner TA, Gee MW, Prokopenko A and Hu JJ: The MueLu Tutorial,
   *          Sandia National Laboratories, Technical Report, 2014
   *          </li>
   *  <li>[4] Prokopenko A, Hu JJ, Wiesner TA, Siefert CM and Tuminaro RS:
   *          MueLu User's Guide 1.0, Sandia National Laboratories, Technical
   *          Report Nr. SAND2014-18874, 2014 </li>
   *  </ul>
   *
   *  \sa NlnOperatorBase
   *
   *  \author mayr.mt \date 06/2014
   */
  class NlnOperatorFas : public NlnOperatorBase
  {
   public:
    //! @name Construction
    //@{

    //! Constructor (empty)
    NlnOperatorFas();

    //! Destructor
    virtual ~NlnOperatorFas() {}

    //@}

    //! @name Setup (public)
    //@{

    //! Setup of the preconditioner object
    virtual void Setup();

    /*! \brief Rebuild preconditioner
     *
     *  To rebuild the FAS preconditioner, we basically need to refresh the
     *  coarse level operators and, thus, need to recompute the RAPs.
     *
     *  \author mayr.mt \date 03/2015
     */
    virtual void RebuildPrec();

    /*! \brief Rebuild preconditioner (as const routine)
     *
     *  \sa RebuildPrec()
     *
     *  \note Can be called from a function, that is not allowed to change
     *  members.
     *
     *  \author mayr.mt \date 03/2015
     */
    virtual void RebuildPrecConst() const;

    //@}

    //! @name Mathematical functions
    //@{

    /*! \brief Apply the preconditioner, i.e. do iterations on FAS-Cycle
     *
     *  Do iterations on FAS cycle, where the actual cycling routine is called
     *  recursively. ApplyIncerse only contains the fine level entry to the
     *  cycling, but no actual computation. This is all done in the respective
     *  cycling routines.
     *
     *  \return Integer error code, set to 0 if successful.
     *
     *  \sa Cycle, VCycle, WCycle
     *
     *  \author mayr.mt \date 03/2014
     */
    virtual int ApplyInverse(
        const Epetra_MultiVector& f_do_not_use,  ///< residual (do not use, just due to inheritance)
        Epetra_MultiVector& x                    ///< solution
        ) const;

    //@}

    //! @name Attribute access functions
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::NlnOperatorFas"; }

    //! Returns a character string describing the operator briefly
    virtual const char* LabelShort() const { return "Fas"; }

    //@}

   protected:
   private:
    //! Compute step length, i.e. line search parameter
    virtual void ComputeStepLength(const Epetra_MultiVector& x,  ///< current solution
        const Epetra_MultiVector& f,                             ///< residual at current solution
        const Epetra_MultiVector& inc,                           ///< search direction
        const double fnorm2,                                     ///< current L2-norm of residual
        double& lsparam,  ///< line search parameter (to be determined)
        bool& suffdecr    ///< achieved sufficient decrease? (to be determined)
        ) const
    {
      dserror("There is no line search in FAS.");
      return;
    }

    //! @name Multigrid cycles
    //@{

    //! Entry to all kind of multigrid cycles
    virtual int Cycle(Teuchos::RCP<Epetra_MultiVector> xbar  ///< solution from finer level
        ) const;

    //! V-Cycle (recursive definition)
    virtual int VCycle(Teuchos::RCP<Epetra_MultiVector> xbar,  ///< solution from finer level
        const int level                                        ///< level ID
        ) const;

    //! V-Cycle with 2 levels (explicit definition)
    virtual int VCycleTwoLevel(Teuchos::RCP<Epetra_MultiVector> x  ///< solution on level 0
        ) const;

    //! V-Cycle with 3 levels (explicit definition)
    virtual int VCycleThreeLevel(Teuchos::RCP<Epetra_MultiVector> x  ///< solution on level 0
        ) const;

    //! W-Cycle // ToDo (mayr) implement W-Cycle
    virtual int WCycle() const;

    //@}

    //! @name Core routines of a multigrid cycle
    //@{

    //! Do the presmoothing
    virtual int PreSmoothing(Epetra_MultiVector& x,  ///< solution from finer level
        const int level                              ///< level ID
        ) const;

    //! Do the postsmoothing
    virtual int PostSmoothing(Epetra_MultiVector& x,  ///< solution from finer level
        const int level                               ///< level ID
        ) const;

    //! Do the coarse level solve
    virtual int CoarseLevelSolve(Epetra_MultiVector& x,  ///< solution from finer level
        const int level                                  ///< level ID
        ) const;

    //@}

    //! Access to hierarchy of multigrid levels
    const Teuchos::RCP<const NLNSOL::FAS::AMGHierarchy> Hierarchy() const;

    //! Hierarchy of multigrid levels
    mutable Teuchos::RCP<NLNSOL::FAS::AMGHierarchy> hierarchy_;

    //! Type of multigrid cycle
    INPAR::NLNSOL::FAS::CycleType cycletype_;

  };  // class NlnOperatorFas
}  // namespace NLNSOL

#endif /* NLN_OPERATOR_FAS_H */
