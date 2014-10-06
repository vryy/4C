/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_dinverse.H

\brief Use inverse of diagonal of Jacobian as preconditioner

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

#ifndef NLN_OPERATOR_D_INVERSE_H
#define NLN_OPERATOR_D_INVERSE_H

/*----------------------------------------------------------------------------*/
/* headers */

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
  class LineSearchBase;
} // namespace NLNSOL

namespace Teuchos
{
  class ParameterList;
} // namespace Teuchos

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \brief Use inverse of diagonal of Jacobian as preconditioner
   *
   *  \sa NLNSOL::NlnOperatorBase()
   *
   *  \author mayr.mt \date 10/2014
   */
  class NlnOperatorDInverse : public NlnOperatorBase
  {
  public:

    //! @name Construction
    //@{

    //! Constructor (empty)
    NlnOperatorDInverse();

    //! Destructor
    virtual ~NlnOperatorDInverse(){}

    //@}

    //! @name Setup
    //@{

    //! Setup of the preconditioner object
    virtual void Setup();

    //@}

    //! @name Mathematical functions
    //@{

    /*! \brief Apply the preconditioner
     *
     * \return Integer error code, set to 0 if successful.
     *
     * \author mayr.mt \date 08/2014
     */
    virtual int ApplyInverse(
        const Epetra_MultiVector& f, ///< residual (dummy vector to satisfy inheritance from Epetra_Operator)
        Epetra_MultiVector& x ///< solution
        )  const;

    //! @name Attribute access functions
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::NlnOperatorDInverse"; }

    //@}

  protected:

  private:

    //! @name Setup
    //@{

    //! Setup of line search
    virtual void SetupLineSearch();

    //@}

    //! @name Search Direction
    //@{

    /*! \brief Compute search direction
     *
     *  The new search direction \f$p\f$ is computed via
     *  \[ p = - D^{-1} f \]
     *  with \f$D = diag(J)\f$ denoting the diagonal of the Jacobian matrix.
     *
     *  \author mayr.mt \date 10/2014
     */
    virtual const int ComputeSearchDirection(
        const Epetra_MultiVector& rhs, ///< residual
        Epetra_MultiVector& inc ///< increment (to be filled)
        ) const;

    //! Compute line search parameter (=step length)
    virtual const double ComputeStepLength(
        const Epetra_MultiVector& x, ///< current solution
        const Epetra_MultiVector& inc, ///< search direction
        double fnorm2 ///< current L2-norm of residual
        ) const;

    //@}

    //! @name Necessary sub-algorithms
    //@{

    //! Line search object
    Teuchos::RCP<NLNSOL::LineSearchBase> linesearch_;

    //@}

  }; // class NlnOperatorDInverse
} // namespace NLNSOL

#endif /* NLN_OPERATOR_D_INVERSE_H */
