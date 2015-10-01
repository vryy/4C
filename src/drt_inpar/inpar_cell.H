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

namespace INPAR
{

namespace CELL
{

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

typedef enum _CELL_MIGRATION_TYPE
{
  cell_migration_undefined=0,
  cell_migration_ameboid=1,
  cell_migration_proteolytic=2
} CELL_MIGRATION_TYPE;

typedef enum _SEGREGATION_BY_TYPE
{
  segregation_by_undefined=0,
  segregation_by_dirichlet=1,
  segregation_by_neumann=2
} SEGREGATION_BY_TYPE;

typedef enum _SEGREGATION_TYPE
{
  segregation_undefined=0,
  segregation_volumetric=1,
  segregation_surface=2
} SEGREGATION_TYPE;

typedef enum _SEGREGATION_LAW_TYPE
{
  segregation_law_undefined=0,
  segregation_law_constant=1,
  segregation_law_forcedependent=2
} SEGREGATION_LAW_TYPE;

/// set the cell parameters
void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

} // CELL
} // INPAR
/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_CELL_H
