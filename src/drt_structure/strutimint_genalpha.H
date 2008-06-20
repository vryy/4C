/*----------------------------------------------------------------------*/
/*!
\file strutimint.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STRUTIMINT_GENALPHA_H
#define STRUTIMINT_GENALPHA_H

/*----------------------------------------------------------------------*/
/* headers */
#include "strutimint_impl.H"


/*======================================================================*/
/*!
\brief Generalised-alpha time integration
\author bborn
*/
class StruTimIntGenAlpha : public StruTimIntImpl
{

public:

  //! Constructor
  StruTimIntGenAlpha
  (
    const Teuchos::ParameterList& sdynparams,
    const Teuchos::ParameterList& genalphaparams,
    DRT::Discretization& actis,
    LINALG::Solver& solver,
    IO::DiscretizationWriter& output
  );

  //! Destructor
  // ....

protected:

  //! equal operator is NOT wanted
  StruTimIntGenAlpha operator=(const StruTimIntGenAlpha& old);

  //! copy constructor is NOT wanted
  StruTimIntGenAlpha(const StruTimIntGenAlpha& old);

  //! @name parameters
  //@{
  double beta_;  //!< factor (0,1/2]
  double gamma_;  //!< factor (0,1]
  double alphaf_;  //!< factor [0,1)
  double alpham_;  //!< factor [0,1)
  //@}

};  // end class StruTimIntGenAlpha

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRUTIMINT_GENALPHA_H
#endif  // #ifdef CCADISCRET
