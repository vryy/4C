/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_interface.cpp

\brief definition of map of instances

<pre>
Maintainer: Andreas Rauch
            rauch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 -15240
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_ele_interface.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/

std::map<int,std::map<int,DRT::ELEMENTS::FluidEleInterface* >* > DRT::ELEMENTS::FluidEleInterface::instances_;
