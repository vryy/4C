/*----------------------------------------------------------------------*
 | fluid3                                                 m.gee 8/00    |
 *----------------------------------------------------------------------*/
typedef struct _FLUID3
{
int                nGP[3];   /* number of gaussian points in rst directions */
int                nGP_tri;  /* number of gaussian points for triangle elements */

int                is_ale;   /* flag whether there is ale to me or not */
struct _ELEMENT   *my_ale;   /* pointer to my ale element, otherwise NULL */
} FLUID3;



/*----------------------------------------------------------------------*
 |  f3_inpele.c                                          m.gee 11/01    |
 *----------------------------------------------------------------------*/
void f3inp(ELEMENT *ele);
