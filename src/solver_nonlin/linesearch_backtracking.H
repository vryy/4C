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
  /*!
   * \brief Line search based on a simple backtracking strategy
   *
   * \author mayr.mt
   * \date 05/2014
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
     *  During backtracking, we simply try a lot of line search parameters in the
     *  interval \f$]0.0, 1.0]\f$ and finally select the one with the minimal
     *  residual norm.
     *
     *  The sufficient decrease condition is checked.
     *
     *  \warning This algorithm has notion of brute force, is not very
     *  efficient, but seems to be very stable and accurate.
     *
     *  \note Use this rather for testing than for real computations.
     *
     *  \sa NLNSOL::LineSearchBase::IsSufficientDecrease()
     *
     *  \author mayr.mt
     *  \date 11/2013
     */
    virtual const double ComputeLSParam() const;

    //@}

  protected:

  private:

    //! increment for trial line search parameters
    double trialstepsize_;

  }; // class LineSearchBacktracking
} // namespace NLNSOL



#endif /* LINESEARCH_BACKTRACKING_H */
