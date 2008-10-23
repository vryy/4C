/*!----------------------------------------------------------------------
\file contactstrugenalpha.H

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich
              
Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed, 
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de) 
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de                   

-------------------------------------------------------------------------
<\pre>

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

namespace CONTACT
{

// forward declarations

/*!
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

  PRIMAL-DUAL ACTIVE SET STRATEGY is interpreted as a fixed-point
  method: The search for the correct active set (=contact nonlinearity)
  is represented by a fixed-point approach, whereas the large deformation
  linearization (=geometrical nonlinearity) is treated by a standard
  Newton scheme. This yields TWO nested iteration loops!

  */
  void FullNewton();

  /*!
  \brief do full Newton iteration with line search

  */
  void FullNewtonLineSearch();
    
  /*!
  \brief do semi-smooth full Newton iteration

  PRIMAL-DUAL ACTIVE SET STRATEGY is interpreted as a semi-smooth Newton
  method: The search for the correct active set (=contact nonlinearity)
  and the large deformation linearization (=geometrical nonlinearity) are
  merged into one semi-smooth Newton method and can therefore besolved
  within ONE iteration loop!
  
  */
  void SemiSmoothNewton();
  
  /*!
  \brief do semi-smooth full Newton iteration with line search

  */
  void SemiSmoothNewtonLineSearch();
    
  /*!
  \brief update configuration and output to file/screen

  */
  virtual void UpdateandOutput();

  /// update configuration
  virtual void Update();

  /// output to file/screen
  virtual void Output();

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
} // namespace CONTACT

#endif  // #ifndef CONTACTSTRUGENALPHA_H
#endif  // #ifdef CCADISCRET
