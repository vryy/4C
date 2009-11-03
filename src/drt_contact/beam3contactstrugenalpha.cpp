/*!----------------------------------------------------------------------
\file contactstrugenalpha.cpp
\brief Generalized Alpha time integration for structural problems with 3D beam contact

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET


#include "beam3contactstrugenalpha.H"
#include <iostream>

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/09|
 *----------------------------------------------------------------------*/
CONTACT::Beam3ContactStruGenAlpha::Beam3ContactStruGenAlpha(
  ParameterList& params,  DRT::Discretization& dis,
  LINALG::Solver& solver, IO::DiscretizationWriter& output) :
StruGenAlpha(params,dis,solver,output)
{
  // print welcome message
  if (!myrank_)
  {
    cout << "\n*******************************";
    cout << "\n* Welcome to 3D BEAM CONTACT! *";
    cout << "\n*******************************\n" << endl;
  }
  
  // only possible to do this if D_BEAM3 activated
#ifdef D_BEAM3
  // -------------------------------------------------------------------
  // check again whether we have beam contact and create beam3cmanager
  // -------------------------------------------------------------------
  // Check for beam contact
  const Teuchos::ParameterList& scontact = DRT::Problem::Instance()->StructuralContactParams();
  if (Teuchos::getIntegralValue<INPAR::CONTACT::ContactType>(scontact,"CONTACT") == INPAR::CONTACT::contact_beams)
    beamcmanager_ = rcp(new CONTACT::Beam3cmanager(dis));
  else
    dserror("ERROR: How did you arrive here...???");
#else
  dserror("ERROR: Beam3 contact time integration called without D_BEAM3 activated...???");
#endif //#ifdef D_BEAM3
  
  return;
} // Beam3ContactStruGenAlpha::Beam3ContactStruGenAlpha

#endif  // #ifdef CCADISCRET
