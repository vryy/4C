/*----------------------------------------------------------------------------*/
/*!
\file linesearch_backtracking.H

\brief Line search based on a simple backtracking strategy

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

#ifndef LINESEARCH_BACKTRACKING_H
#define LINESEARCH_BACKTRACKING_H

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
  /*! \class LineSearchBacktracking
   *  \brief Line search based on a simple backtracking strategy
   *
   *  We start with a trial step length of \f$alpha = 1.0\f$. if the sufficient
   *  decrease condition is not met, the trial step length is halved.
   *
   *  <h3> References </h3>
   *  <ul>
   *  <li>[1] Kelley CT: Iterative Methods for Linear and Nonlinear Equations,
   *          SIAM, Philadelphia (1995)
   *  </ul>
   *
   *  \author mayr.mt \date 01/2015
   */
  class LineSearchBacktracking : public NLNSOL::LineSearchBase
  {
   public:
    //! @name Construction / Destruction
    //@{

    //! constructor (empty)
    LineSearchBacktracking();

    //! Destructor
    virtual ~LineSearchBacktracking(){};

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
     *  During backtracking, we start with a line search parameter
     *  \f$\alpha = 1.0]\f$. It will be halved by \f$2\f$ until the sufficient
     *  decrease condition is met.
     *
     *  The sufficient decrease condition is checked.
     *
     *  \sa NLNSOL::LineSearchBase::IsSufficientDecrease()
     *
     *  \author mayr.mt \date 01/2015
     */
    virtual void ComputeLSParam(double& lsparam,  ///< line search parameter (to be determined)
        bool& suffdecr  ///< achieved sufficient decrease? (to be determined)
        ) const;

    //@}

    //! @name Description of this class
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::LineSearchBacktracking"; }

    //! Returns a character string describing the operator briefly
    virtual const char* LabelShort() const { return "LS_Backtracking"; }

    //@}

   protected:
   private:
    //! maximum number of backtracking steps
    int itermax_;

    //! factor for step size reduction per iteration
    double backtrackfac_;

  };  // class LineSearchBacktracking
}  // namespace NLNSOL



#endif /* LINESEARCH_BACKTRACKING_H */
