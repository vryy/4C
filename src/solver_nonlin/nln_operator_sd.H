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
}  // namespace NLNSOL

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \class NlnOperatorSD
   *  \brief Steepest-descent algorithm as nonlinear preconditioner
   *
   *  When solving nonlinear equations, steepest descent means that the residual
   *  is chosen as the search direction.
   *
   *  \sa NlnOperatorBase
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
    virtual ~NlnOperatorSD() {}

    //@}

    //! @name Setup
    //@{

    //! Setup of the preconditioner object
    virtual void Setup();

    /*! \brief Rebuild preconditioner
     *
     *  No action required.
     *
     *  \author mayr.mt \date 03/2015
     */
    virtual void RebuildPrec(){/* no action required */};

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
        const Epetra_MultiVector&
            f,  ///< residual (dummy vector to satisfy inheritance from Epetra_Operator)
        Epetra_MultiVector& x  ///< solution
        ) const;

    //! @name Attribute access functions
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::NlnOperatorSD"; }

    //! Returns a character string describing the operator briefly
    virtual const char* LabelShort() const { return "SD"; }

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

    //! Compute search direction
    virtual int ComputeSearchDirection(const Epetra_MultiVector& rhs,  ///< residual
        Epetra_MultiVector& inc                                        ///< increment (to be filled)
        ) const;

    //! Compute line search parameter (=step length)
    virtual void ComputeStepLength(const Epetra_MultiVector& x,  ///< current solution
        const Epetra_MultiVector& f,                             ///< residual at current solution
        const Epetra_MultiVector& inc,                           ///< search direction
        const double fnorm2,                                     ///< current L2-norm of residual
        double& lsparam,  ///< line search parameter (to be determined)
        bool& suffdecr    ///< achieved sufficient decrease? (to be determined)
        ) const;

    //@}

    //! @name Necessary sub-algorithms
    //@{

    //! Line search object
    Teuchos::RCP<NLNSOL::LineSearchBase> linesearch_;

    //@}

  };  // class NlnOperatorSD
}  // namespace NLNSOL

#endif /* NLN_OPERATOR_SD_H */
