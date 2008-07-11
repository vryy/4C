/*!----------------------------------------------------------------------
\file combust_create_gfunction.H

\brief create G-function discretization based on fluid discretization

	The G-function field is a convection-diffusion field. Based on 
	the fluid discretization, the G-function discretization is 
	created. BACI knows two discretizations.
	
<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef COMBUST_CREATE_GFUNCTION_H
#define COMBUST_CREATE_GFUNCTION_H

namespace COMBUST
{

/*----------------------------------------------------------------------*
 * create a G-function discretization from the fluid        henke 06/08 * 
 *----------------------------------------------------------------------*/
void CreateGfuncDiscretization(int disnumff,int disnumgff);

} //namespace COMBUST

#endif // #ifdef COMBUST_CREATE_GFUNCTION_H
#endif // #ifdef CCADISCRET
