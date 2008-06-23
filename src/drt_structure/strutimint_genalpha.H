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

  //! Mid-average type of internal forces
  enum MidAverageEnum
  {
    midavg_vague = 0,  //!< undefined mid-averaging type
    midavg_imrlike,  //!< alphaf-mid-averaging is done IMR-like, i.e. 
                     //!< \f$F_{int,m}\f$ 
                     //!< \f$= F_{int}(D_m)\f$
                     //!< \f$= F_{int}( (1-\alpha_f)*D_{n+1} + \alpha_f*D_n )\f$
                     //!< (IMR means implicit mid-point rule.)
    midavg_trlike    //!< alphaf-mid-averaging is done TR-like, i.e.
                     //!< \f$F_{int,m}\f$
                     //!< \f$= (1-\alpha_f)*F_{int,n+1} + \alpha_f*F_{int,n}\f$
                     //!<  \f$= (1-\alpha_f)*F_{int}(D_{n+1}) + \alpha_f*F_{int}(D_n)\f$
                     //!<  (TR means trapezoidal rule.)
  };

  //! convert input string to enum for mid-average type
  enum MidAverageEnum MidAvgStringToEnum
  (
    const std::string instring  //!< input string
  );

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

  //! @name set-up
  //@{
  enum MidAverageEnum midavg_;
  //@}

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
