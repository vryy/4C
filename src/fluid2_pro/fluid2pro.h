/*!---------------------------------------------------------------------
\file
\brief The 2D fluid element used for projection methods

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/

/*!---------------------------------------------------------------------
\brief fluid2

<pre>                                                         genk 03/02

In this structure all variables used for element evaluation by the 2D
fluid element fluid2_pro are stored.

</pre>

--------------------------------------------------------------------------*/
typedef struct _FLUID2_PRO
{

INT                ntyp;     /*!< flag for element type: 1=quad; 2=tri    */
INT                nGP[2];   /*!< number of gaussian points in rs direct. */
INT                is_ale;   /*!< flag whether there is ale to me or not  */
struct _ELEMENT   *my_ale;   /*!< pointer to my ale ele, otherwise NULL   */
enum   _DISMODE    dm;

} FLUID2_PRO;

