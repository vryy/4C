/*!---------------------------------------------------------------------
\file
\brief The 2D fluid element

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

int                ntyp;     /*!< flag for element type: 1=quad; 2=tri    */
int                nGP[2];   /*!< number of gaussian points in rs direct. */
int                is_ale;   /*!< flag whether there is ale to me or not  */
struct _ELEMENT   *my_ale;   /*!< pointer to my ale ele, otherwise NULL   */

/*--------------------------------- element sizes for stability parameter */
double             strom;    /*!< streamlenght=elemtlenght                */
double             strom_dc; /*!< streamlenghtfor DISC. CAPT.             */

} FLUID2_TU;

