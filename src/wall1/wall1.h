#include "wall1_calc.h"

/*----------------------------------------------------------------------*
 | type of 2D problem                                       al 01/02    |
 *----------------------------------------------------------------------*/
typedef enum _WALL_TYPE
{
                       plain_strain, 
                       plain_stress,
                       rotat_symmet  
} WALL_TYPE;
/*----------------------------------------------------------------------*
 | wall                                                     al 01/02    |
 *----------------------------------------------------------------------*/
typedef struct _WALL1
{
enum _WALL_TYPE    wtype;       /* type of 2D problem, see above */

int           nGP[4];

double        thick;

W1_ELE_WA     *elewa;

} WALL1;
