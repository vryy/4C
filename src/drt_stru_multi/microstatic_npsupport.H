/*!----------------------------------------------------------------------
\file microstatic_npsupport.H

<pre>
Maintainer: Lena Yoshihara
            yoshihara@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/

#ifndef STRU_MICRO_STAT_NPSUPPORT_H
#define STRU_MICRO_STAT_NPSUPPORT_H

/*!
\brief Namespace for structural multi-scale related algorithms

\author yoshihara (yoshihara@lnm.mw.tum.de)
*/
namespace STRUMULTI
{
  // std::endless loop for supporting procs in multi scale problems
  void np_support_drt();


}  // namespace STRUMULTI
#endif
