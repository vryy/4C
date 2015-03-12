/*!----------------------------------------------------------------------
\file inpar_cell.H

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#ifndef INPAR_CELL_H
#define INPAR_CELL_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* The coupling methods for Cell Migration. */
/*----------------------------------------------------------------------*/
typedef enum _CELL_COUPLING
{
  cell_coupling_undefined=0,
  cell_basic_sequ_stagg=1,
  cell_iter_stagg_fixed_rel_param=2,
  cell_iter_stagg_AITKEN_rel_param=3,
} CELL_COUPLING;


/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_FSI_H
