/*!----------------------------------------------------------------------
\file contactstrugenalpha.H

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef BEAM3CONTACTSTRUGENALPHA_H
#define BEAM3CONTACTSTRUGENALPHA_H

#include "../drt_structure/strugenalpha.H"

#ifdef D_BEAM3
#include "beam3contact_manager.H"
#endif // #ifdef D_BEAM3

namespace CONTACT
{

// forward declarations

/*!
\brief Generalized Alpha time integration for structural problems with 3D beam contact

\author popp (popp@lnm.mw.tum.de)
*/
class Beam3ContactStruGenAlpha : public virtual StruGenAlpha
{
public:

  /*!
  \brief Standard Constructor

  */
  explicit Beam3ContactStruGenAlpha(ParameterList& params,
                               DRT::Discretization& dis,
                               LINALG::Solver& solver,
                               IO::DiscretizationWriter& output);

  /*!
  \brief Destructor

  */
  virtual ~Beam3ContactStruGenAlpha() {};

protected:

  // don't want = operator and cctor
  Beam3ContactStruGenAlpha operator = (const Beam3ContactStruGenAlpha& old);
  Beam3ContactStruGenAlpha(const Beam3ContactStruGenAlpha& old);

  // only possible to do this if D_BEAM3 activated
#ifdef D_BEAM3
  RCP<CONTACT::Beam3cmanager>    beamcmanager_;
#endif // #ifdef D_BEAM3
  
}; // class Beam3ContactStruGenAlpha
} // namespace CONTACT

#endif  // #ifndef BEAM3CONTACTSTRUGENALPHA_H
#endif  // #ifdef CCADISCRET
