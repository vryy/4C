/*----------------------------------------------------------------------*/
/*!
\file inpar_cell.H

\brief input parameters and conditions for cell migration

\maintainer Andreas Rauch
            rauch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 -15240

\level 2
*/
/*----------------------------------------------------------------------*/

#ifndef INPAR_CELL_H
#define INPAR_CELL_H

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT { namespace INPUT { class ConditionDefinition; } }

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

typedef enum _SIM_TYPE
{
  sim_type_pureFSI=0,
  sim_type_pureAdhesion=1,
  sim_type_pureConfinement=2,
  sim_type_pureProtrusionFormation=3,
  sim_type_pureContraction=4,
  sim_type_pureEndoExocytosis=5,
  sim_type_multiphysics=6
} SIM_TYPE;

typedef enum _ART_FLUID_MOD
{
  mod_none=0,
  mod_simplified=1,
  mod_full=2
} ART_FLUID_MOD;

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

typedef enum _GROWTH_COUP_TYPE
{
  coup_growth_undefined=0,
  coup_growth_growth=1,
  coup_growth_ssi=2
} GROWTH_COUP_TYPE;

typedef enum _EXOCYTOSIS_TYPE
{
  exo_surface=0,
  exo_volume=1
} EXOCYTOSIS_TYPE;

/// set the cell parameters
void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

/// set specific cell conditions
void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);

} // CELL
} // INPAR
/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_CELL_H
