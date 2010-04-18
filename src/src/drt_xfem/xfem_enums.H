/*----------------------------------------------------------------------*/
/*!
\file xfem_enums.H

\brief basic enums used in xfem routines

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef XFEM_ENUMS_H
#define XFEM_ENUMS_H

namespace XFEM
{

  //! Allows to distinguish between normal assembly and xfem assembly
  enum AssemblyType
  {
    standard_assembly,    ///< fixed layout of element matrices
    xfem_assembly         ///< variable layout of element matrices - ask dofmanagers for layout information
  };

}

#endif
#endif
