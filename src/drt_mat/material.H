/*----------------------------------------------------------------------*/
/*!
\file material.H

\brief base class for handling of history variables in Gauss points
arising from complex material behaviour

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef MATERIAL_H
#define MATERIAL_H

#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

class Material
{
public:

  virtual ~Material();

};

#endif
#endif
#endif
