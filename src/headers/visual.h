/*----------------------------------------------------------------------*
 | control structure                                      m.gee 5/02    |
 | for nonlinear structural dynamics                                    |           
 *----------------------------------------------------------------------*/
typedef struct _VISUAL_DATA
{
int                          step;     

long int                     time;

long int                     handle_of_node_handles;
int                          node_fdim;
int                          node_sdim;
long int                    *node_handles;

} VISUAL_DATA;
