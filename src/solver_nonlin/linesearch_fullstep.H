/*----------------------------------------------------------------------------*/
/*!
\file linesearch_fullstep.H

\brief Line search with fixed full step

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

#ifndef LINESEARCH_FULLSTEP_H_
#define LINESEARCH_FULLSTEP_H_

/*----------------------------------------------------------------------*/
/* headers */

// baci
#include "linesearch_base.H"

/*----------------------------------------------------------------------*/
/* forward declarations */

/*----------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*!
   * \brief Line search with fixed full step
   *
   * \author mayr.mt
   * \date 05/2014
   */
  class LineSearchFullStep : public NLNSOL::LineSearchBase
  {
  public:

    //! @name Construction / Destruction
    //@{

    //! constructor (empty)
    LineSearchFullStep();

    //! Destructor
    virtual ~LineSearchFullStep(){};

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
     *  We always do a full step, here. So just do some initial checks to be
     *  consistent with the other line search algorithms and then return 1.0.
     *
     *  We do not care about the sufficient decrease condition, here.
     *
     *  \author mayr.mt \date 06/2014
     */
    virtual const double ComputeLSParam() const;

    //@}

  protected:

  private:

  };
}



#endif /* LINESEARCH_FULLSTEP_H_ */
