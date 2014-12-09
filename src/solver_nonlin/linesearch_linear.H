/*----------------------------------------------------------------------------*/
/*!
\file linesearch_linear.H

\brief Line search based on linearity assumption

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

#ifndef LINESEARCH_LINEAR_H
#define LINESEARCH_LINEAR_H

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
   * \brief Line search based on linearity assumption
   *
   * \author mayr.mt \date 10/2014
   */
  class LineSearchLinear : public NLNSOL::LineSearchBase
  {
  public:

    //! @name Construction / Destruction
    //@{

    //! constructor (empty)
    LineSearchLinear();

    //! Destructor
    virtual ~LineSearchLinear(){};

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
     *  We use the analytical formula to compute the exact and optimal line
     *  search parameter for the case of a linear problem:
     *
     *  \[\alpha = \]
     *
     *  \author mayr.mt
     *  \date 10/2014
     */
    virtual const double ComputeLSParam() const;

    //@}

  protected:

  private:

  }; // class LineSearchLinear
} // namespace NLNSOL



#endif /* LINESEARCH_LINEAR_H */
