/*----------------------------------------------------------------------*/
/*!
\file elch_create_condif.H

\brief Create convection-diffusion discretization based on fluid discretization

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ELCH_CREATE_CONDIF_H
#define ELCH_CREATE_CONDIF_H

namespace ELCH
{
/// automatically create a matching convection-diffusion discretization 
/// from the fluid discretization
void CreateConDifDiscretization(int disnumff,int disnumscatra);

}

#endif
#endif
