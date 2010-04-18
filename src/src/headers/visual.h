/*!---------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

---------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | control structure                                      m.gee 5/02    |
 | for nonlinear structural dynamics                                    |
 *----------------------------------------------------------------------*/
typedef struct _VISUAL_DATA
{
INT                          step;

long int                     time;

long int                     handle_of_node_handles;
INT                          node_fdim;
INT                          node_sdim;
long int                    *node_handles;

} VISUAL_DATA;
