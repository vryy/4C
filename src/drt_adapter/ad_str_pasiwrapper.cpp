/*!----------------------------------------------------------------------
\file ad_str_pasiwrapper.cpp

\brief structural adapter for PASI problems

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                               sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
#include "ad_str_pasiwrapper.H"

/*----------------------------------------------------------------------*
 | pasi adapter                                          sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
ADAPTER::PASIStructureWrapper::PASIStructureWrapper(Teuchos::RCP<Structure> structure)
: StructureWrapper(structure)
{

} // ADAPTER::PASIStructureWrapper::PASIStructureWrapper()
