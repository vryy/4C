/*----------------------------------------------------------------------*
 | control structure                                      m.gee 5/02    |
 | for nonlinear structural dynamics                                    |           
 *----------------------------------------------------------------------*/
typedef struct _RESTART_DYNSTRUCT
{
int                          step;     
struct _STRUCT_DYNAMIC       sdyn;
struct _STRUCT_DYN_CALC      dynvar;

long int                     dist_vec_rhs[10];
long int                     dist_vec_sol[10];
long int                     dist_vec_dispi[10];
long int                     dist_vec_vel[10];
long int                     dist_vec_acc[10];
long int                     dist_vec_fie[10];
long int                     dist_vec_work[10];

long int                     intforce;
long int                     dirich;

long int                     handle_of_node_handles;
int                          node_fdim;
int                          node_sdim;
long int                   **node_handles;

long int                     handle_of_ele_handles;
int                          ele_fdim;
int                          ele_sdim;
long int                   **ele_handles;

} RESTART_DYNSTRUCT;
