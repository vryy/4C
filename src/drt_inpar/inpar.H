/*----------------------------------------------------------------------*/
/*!
\file inpar.H

\brief General input parameters

<pre>
Maintainer: Karl-Robert Wichmann
            wichmann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_H
#define INPAR_H

/*----------------------------------------------------------------------*/
namespace INPAR
{

  //! Geometry reading specification
  enum GeometryType
  {
    geometry_full,
    geometry_box,
    geometry_file
  };

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_H
