/*----------------------------------------------------------------------*/
/*!
\file thr_contact.cpp
\brief Thermal contact routines for

<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                                     06/11 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*
 | headers                                                    mgit 10/10 |
 *----------------------------------------------------------------------*/
#include <sstream>

#include "thrtimint.H"
#include "thrtimint_impl.H"
#include "thr_aux.H"
#include "thr_contact.H"

#include "../drt_mortar/mortar_manager_base.H"
#include "../drt_mortar/mortar_utils.H"
#include "../drt_contact/meshtying_manager.H"
#include "../drt_contact/contact_manager.H"
#include "../drt_contact/contact_interface.H"
#include "../drt_contact/contact_abstract_strategy.H"
#include "../drt_contact/contact_node.H"
#include "../drt_contact/friction_node.H"


/*----------------------------------------------------------------------*
 | constructor                                               mgit 06/11 |
 *----------------------------------------------------------------------*/
THR::ThermoContactMan::ThermoContactMan(Teuchos::RCP<MORTAR::ManagerBase> cmtman,
                                        Teuchos::RCP<DRT::Discretization> discretstruct,
                                        Teuchos::RCP<DRT::Discretization> discretthermo):
cmtman_(cmtman),
discretstruct_(discretstruct),
discretthermo_(discretthermo)
{
  // done so far
  return;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
