/*----------------------------------------------------------------------*
 | control structure                                      m.gee 5/02    |
 | for nonlinear structural dynamics                                    |           
 *----------------------------------------------------------------------*/
typedef struct _RESTART_DYNSTRUCT
{
INT                          step;     
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
INT                          node_fdim;
INT                          node_sdim;
long int                   **node_handles;

long int                     handle_of_ele_handles;
INT                          ele_fdim;
INT                          ele_sdim;
long int                   **ele_handles;

} RESTART_DYNSTRUCT;

/*----------------------------------------------------------------------*
 | control structure                                      m.gee 5/02    |
 | for nonlinear structural dynamics                                    |           
 *----------------------------------------------------------------------*/
typedef struct _RESTART_STATSTRUCT
{
INT                          step;     
struct _STATIC_VAR           statvar;
struct _STANLN               nln_data;

long int                     dist_vec_rhs[10];
long int                     dist_vec_sol[10];
long int                     dist_vec_dispi[10];

long int                     handle_of_node_handles;
INT                          node_fdim;
INT                          node_sdim;
long int                   **node_handles;

long int                     handle_of_ele_handles;
INT                          ele_fdim;
INT                          ele_sdim;
long int                   **ele_handles;

} RESTART_STATSTRUCT;

