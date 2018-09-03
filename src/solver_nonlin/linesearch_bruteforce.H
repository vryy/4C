/*----------------------------------------------------------------------------*/
/*!
\file linesearch_bruteforce.H

\brief Line search based on a brute force backtracking strategy

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

#ifndef LINESEARCH_BRUTEFORCE_H
#define LINESEARCH_BRUTEFORCE_H

/*----------------------------------------------------------------------------*/
/* headers */

// baci
#include "linesearch_base.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \class LineSearchBruteForce
   *  \brief Line search based on a brute force backtracking strategy
   *
   *  The step length parameter is varied by constant increments in an interval
   *  \f$\alpha\in]0,1]\f$. Then, the \f$alpha\f$ that results in the minimal
   *  residual norm is chosen as the final step length.
   *
   *  This algorithms seems to be quite accurate, though it is very expensive.
   *
   *  \author mayr.mt \date 01/2015
   */
  class LineSearchBruteForce : public NLNSOL::LineSearchBase
  {
   public:
    //! @name Construction / Destruction
    //@{

    //! constructor (empty)
    LineSearchBruteForce();

    //! Destructor
    virtual ~LineSearchBruteForce(){};

    //@}

    //! @name Setup
    //@{

    //! Setup the line search object
    virtual void Setup();

    //@}

    //! @name Computation
    //@{

    /*! \brief Compute the step length / line search parameter
     *
     *  During brute force backtracking, we simply try a lot of line search
     *  parameters in the interval \f$]0.0, 1.0]\f$ and finally select the one
     *  with the minimal residual norm. This seems to be quite accurate, though
     *  very costly.
     *
     *  The sufficient decrease condition is checked.
     *
     *  \warning This algorithm is brute force and, thus, not very efficient,
     *  but seems to be very stable and accurate.
     *
     *  \note Use this rather for testing than for real computations.
     *
     *  \sa NLNSOL::LineSearchBase::IsSufficientDecrease()
     *
     *  \author mayr.mt \date 11/2013
     */
    virtual void ComputeLSParam(double& lsparam,  ///< line search parameter (to be determined)
        bool& suffdecr  ///< achieved sufficient decrease? (to be determined)
        ) const;

    //@}

    //! @name Description of this class
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::LineSearchBruteForce"; }

    //! Returns a character string describing the operator briefly
    virtual const char* LabelShort() const { return "LS_BruteForce"; }

    //@}

   protected:
   private:
    //! increment for trial line search parameters
    double trialstepsize_;

  };  // class LineSearchBruteForce
}  // namespace NLNSOL



#endif /* LINESEARCH_BRUTEFORCE_H */
