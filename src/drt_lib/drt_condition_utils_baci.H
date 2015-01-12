/*----------------------------------------------------------------------*/
/*!
\file drt_condition_utils_baci.H

\brief

<pre>
Maintainer: Christian Roth
            roth@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef DRT_CONDITION_UTILS_BACI_H
#define DRT_CONDITION_UTILS_BACI_H


#include <Teuchos_RCP.hpp>


namespace DRT
{
  class Element;

  namespace UTILS
  {

    /// Writes debug output for volumetric coupling
    void WriteBoundarySurfacesVolumeCoupling(
        std::map< std::vector<int>, Teuchos::RCP<DRT::Element> > surfmap, //map containing the surfaces of the boundary condition
        int condID,                                                       //condition ID
        int numproc,                                                      //number of procs
        int mypid                                                         //current proc
        );

  }
}

#endif
