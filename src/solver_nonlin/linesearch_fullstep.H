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

#ifndef LINESEARCH_FULLSTEP_H
#define LINESEARCH_FULLSTEP_H

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
  /*! \class LineSearchFullStep
   *  \brief Line search with fixed full step
   *
   *  We do not compute the step length. We always take the full step,
   *  i.e. step length \f$\alpha = 1.0\f$.
   *
   * \author mayr.mt \date 05/2014
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
     *  We do not care about the sufficient decrease condition, here. However,
     *  we evaluate it since the calling nonlinear operator might care.
     *
     *  \author mayr.mt \date 06/2014
     */
    virtual void ComputeLSParam(double& lsparam,  ///< line search parameter (to be determined)
        bool& suffdecr  ///< achieved sufficient decrease? (to be determined)
        ) const;

    /*! \brief Compute the step length / line search parameter
     *
     *  We always do a full step, here. So just do some initial checks to be
     *  consistent with the other line search algorithms and then return 1.0.
     *
     *  We do not care about the sufficient decrease condition, here.
     *
     *  \author mayr.mt \date 06/2014
     */
    virtual void ComputeLSParam(double& lsparam  ///< line search parameter (to be determined)
        ) const;

    //@}

    //! @name Description of this class
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::LineSearchFullStep"; }

    //! Returns a character string describing the operator briefly
    virtual const char* LabelShort() const { return "LS_FullStep"; }

    //@}

   protected:
   private:
  };  // class LineSearchFullStep
}  // namespace NLNSOL



#endif /* LINESEARCH_FULLSTEP_H */
