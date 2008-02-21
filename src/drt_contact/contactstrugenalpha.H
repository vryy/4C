/*!----------------------------------------------------------------------
\file contactstrugenalpha.H

\class ContactStruGenAlpha

\brief Generalized Alpha time integration for structural problems with contact

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef CONTACTSTRUGENALPHA_H
#define CONTACTSTRUGENALPHA_H

#include "../drt_structure/strugenalpha.H"

/*!
\class ContactStruGenAlpha

\brief Generalized Alpha time integration for structural problems with contact

\author popp (popp@lnm.mw.tum.de)
*/
class ContactStruGenAlpha : public virtual StruGenAlpha
{
public:

  /*!
  \brief Standard Constructor

  */
  explicit ContactStruGenAlpha(ParameterList& params,
                               DRT::Discretization& dis,
                               LINALG::Solver& solver,
                               IO::DiscretizationWriter& output);

  /*!
  \brief Destructor

  */
  virtual ~ContactStruGenAlpha() {};

  /*!
  \brief Do time integration

  */
  void Integrate();

  /*!
  \brief do constant predictor step

  */
  void ConstantPredictor();

  /*!
  \brief do consistent constant predictor step

  */
  void ConsistentPredictor();

  /*!
  \brief do full Newton iteration

  */
  void FullNewton();

  /*!
  \brief do pseudo transient continuation nonlinear iteration

  Pseudo transient continuation is a variant of full newton that has a
  larger convergence radius than newton and is therefore more stable and/or
  can do larger time steps

  */
  virtual void PTC();
    
  /*!
  \brief update configuration and output to file/screen

  */
  void UpdateandOutput();

  /*!
  \brief restart calculation from a given step

  */
  void ReadRestart(int step);

protected:

  // don't want = operator and cctor
  ContactStruGenAlpha operator = (const ContactStruGenAlpha& old);
  ContactStruGenAlpha(const ContactStruGenAlpha& old);

  RCP<CONTACT::Manager>         contactmanager_;

}; // class ContactStruGenAlpha

#endif  // #ifndef CONTACTSTRUGENALPHA_H
#endif  // #ifdef CCADISCRET
