/*----------------------------------------------------------------------*/
/*!
\file nln_operator_fas.H

\brief AMG-FAS as nonlinear preconditioner

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef NLN_OPERATOR_FAS_H
#define NLN_OPERATOR_FAS_H

/*----------------------------------------------------------------------*/
/* headers */

// standard
#include <vector>
#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>

// baci
#include "nln_operator.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*/
/* forward declarations */

class Epetra_Comm;
class Epetra_MultiVector;

namespace NLNSOL
{
  namespace FAS
  {
    class AMGHierarchy;
    class NlnLevel;
  }
}

namespace Teuchos
{
  class ParameterList;
}

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
   *  NLNSOL::FAS::AMGHierarchy in a collection of NLNSOL::FAS::NlnLevel objects.
   *
   *  \sa NLNSOL::NlnOperator()
   *
   *  \author mayr.mt \date 06/2014
   */
  class NlnOperatorFas : public NlnOperator
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
    virtual const char* Label() const { return "NlnOperatorFas"; }

    //@}

  protected:

  private:

    //! @name Multigrid cycles
    //@{

    //! V-Cycle (recursive definition)
    virtual void VCycle(const Epetra_MultiVector& f,  ///< residual
                        Epetra_MultiVector& x,        ///< solution
                        const int level               ///< level ID
                        ) const;

    //! W-Cycle
    virtual void WCycle() const { dserror("Not implemented, yet."); }

    //@}

    //! Access to hierarchy of multigrid levels
    const Teuchos::RCP<const NLNSOL::FAS::AMGHierarchy> Hierarchy() const;

    //! Hierarchy of multigrid levels
    Teuchos::RCP<NLNSOL::FAS::AMGHierarchy> hierarchy_;

  };
}

#endif /* NLN_OPERATOR_FAS_H */
