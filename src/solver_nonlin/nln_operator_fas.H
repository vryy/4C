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
  } // namespace FAS
} // namespace NLNSOL

namespace Teuchos
{
  class ParameterList;
} // namespace Teuchos

/*----------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \brief Nonlinear Algebraic Multigrid with Full Approximation Scheme.
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
   *  <li>[1] Briggs, W. L., Van Henson, E. and McCormick, S. F: A Multigrid
   *          Tutorial, SIAM, ed. 2, 2000 </li>
   *  <li>[2] Gee, M. W and Tuminaro, R. S.: Nonlinear algebraic multigrid for
   *          constrained solid mechanics problems, Sandia National
   *          Laboratories, Technical Report Nr. SAND2006-2256, 2006 </li>
   *  <li>[3] Wiesner, T. A., Gee, M. W, Prokopenko, A. and Hu, J.: The MueLu
   *          Tutorial, Sandia National Laboratories, Technical Report, 2014
   *          </li>
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
    virtual ~NlnOperatorFas(){}

    //@}

    //! @name Setup (public)
    //@{

    //! Setup of the preconditioner object
    virtual void Setup();

    //@}

    //! @name Mathematical functions
    //@{

    /*! \brief Apply the preconditioner, i.e. do one step of Quasi-Newton
     *
     * \return Integer error code, set to 0 if successful.
     *
     * \author mayr.mt \date 03/2014
     */
    virtual int ApplyInverse(const Epetra_MultiVector& f, ///< residual
        Epetra_MultiVector& x ///< solution
        )  const;

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
    virtual const double ComputeStepLength(
        const Epetra_MultiVector& x, ///< current solution
        const Epetra_MultiVector& f, ///< residual at current solution
        const Epetra_MultiVector& inc, ///< search direction
        double fnorm2 ///< current L2-norm of residual
        ) const { dserror("There is no line search in FAS."); return 0.0; }

    //! @name Multigrid cycles
    //@{

    //! V-Cycle (recursive definition)
    virtual void VCycle(const Epetra_MultiVector& f, ///< residual
                        Epetra_MultiVector& x, ///< solution
                        const int level ///< level ID
                        ) const;

    //! W-Cycle
    virtual void WCycle() const { dserror("Not implemented, yet."); }

    //@}

    //! Access to hierarchy of multigrid levels
    const Teuchos::RCP<const NLNSOL::FAS::AMGHierarchy> Hierarchy() const;

    //! Hierarchy of multigrid levels
    Teuchos::RCP<NLNSOL::FAS::AMGHierarchy> hierarchy_;

  }; // class NlnOperatorFas
} // namespace NLNSOL

#endif /* NLN_OPERATOR_FAS_H */
