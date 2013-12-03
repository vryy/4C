/*----------------------------------------------------------------------*/
/*!
\file adapter_structure_timint.cpp

\brief Structure field adapter for time step size adaptivity within monolithic FSI

<pre>
Maintainer: Matthias Mayr
            mayr@lnm.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "ad_str_fsi_timint_adaptive.H"


/*======================================================================*/
/* constructor */
ADAPTER::StructureFSITimIntAda::StructureFSITimIntAda(
  Teuchos::RCP<STR::TimAda> sta,
  Teuchos::RCP<Structure> sti
)
: FSIStructureWrapper(sti),
  StructureTimIntAda(sta, sti)
{

}
