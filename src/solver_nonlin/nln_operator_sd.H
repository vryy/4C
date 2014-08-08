/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_sd.H

\brief Steepest-descent algorithm

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

#ifndef NLN_OPERATOR_SD_H
#define NLN_OPERATOR_SD_H

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
  /*! \brief Steepest-descent algorithm as nonlinear preconditioner
   *
   *  \sa NLNSOL::NlnOperatorBase()
   *
   *  \author mayr.mt \date 08/2014
   */
  class NlnOperatorSD : public NlnOperatorBase
  {
  public:

    //! @name Construction
    //@{

    //! Constructor (empty)
    NlnOperatorSD();

    //! Destructor
    virtual ~NlnOperatorSD(){}

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
    virtual int ApplyInverse(const Epetra_MultiVector& f, ///< residual (dummy vector to satisfy inheritance from Epetra_Operator)
        Epetra_MultiVector& x ///< solution
        )  const;

    //! @name Attribute access functions
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::NlnOperatorSD"; }

    //@}

  protected:

  private:

    //! @name Setup
    //@{

    //! Setup of line search
    virtual void SetupLineSearch();

    //@}

    //! @name Access routines
    //@{

    //! Get max allowed number of iterations
    virtual const int GetMaxIter() const { return maxiter_; }

    //@}

    //! @name Search Direction
    //@{

    virtual const int ComputeSearchDirection(
        const Epetra_MultiVector& rhs, ///< residual
        Epetra_MultiVector& inc ///< increment (to be filled)
        ) const;

    //@}

    //! @name Necessary sub-algorithms
    //@{

    //! line search object
    Teuchos::RCP<NLNSOL::LineSearchBase> linesearch_;

    //@}

    //! @name Algorithmic parameters
    //@{

    //! max number of iterations
    int maxiter_;

    //@}

  }; // class NlnOperatorSD
} // namespace NLNSOL

#endif /* NLN_OPERATOR_SD_H */
