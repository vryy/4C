/*----------------------------------------------------------------------*/
/*!
\file xfem_utils_combust.H

\brief basic enums used in xfem routines

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
*/
/*----------------------------------------------------------------------*/


#ifndef XFEM_UTILS_COMBUST_H
#define XFEM_UTILS_COMBUST_H


namespace XFEM
{
  //! global Position of degree of freedom in Epetra vector
  typedef int DofGID;

  //! Allows to distinguish between normal assembly and xfem assembly
  enum AssemblyType
  {
    standard_assembly,    ///< fixed layout of element matrices
    xfem_assembly         ///< variable layout of element matrices - ask dofmanagers for layout information
  };

  //! identify interface side due to phivalue value
  inline bool plusDomain(double phiValue)
  {
    double TOL = 1.0e-14;
    if (phiValue<-TOL)
      return false;
    else
      return true;
  };
}

#endif
