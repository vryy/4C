/*!---------------------------------------------------------------------
\file
\brief The 2D fluid element

<pre>
Maintainer: Thomas Hettich
            hettich@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hettich/
            0771 - 685-6575
</pre>

------------------------------------------------------------------------*/

/*!---------------------------------------------------------------------
\brief fluid2                                                 

<pre>                                                         he  12/02 

In this structure all variables used for element evaluation by the 2D
fluid element fluid2_tu are stored.

</pre>

--------------------------------------------------------------------------*/
typedef struct _FLUID2_TU
{

INT                nGP[2];   /*!< number of gaussian points in rs direct. */
INT                is_ale;   /*!< flag whether there is ale to me or not  */
struct _ELEMENT   *my_ale;   /*!< pointer to my ale ele, otherwise NULL   */

/*--------------------------------- element sizes for stability parameter */
DOUBLE             strom;    /*!< streamlenght=elemtlenght                */
DOUBLE             strom_dc; /*!< streamlenghtfor DISC. CAPT.             */

} FLUID2_TU;

