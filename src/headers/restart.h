/*----------------------------------------------------------------------*
 | control structure                                      m.gee 5/02    |
 | for nonlinear structural dynamics                                    |           
 *----------------------------------------------------------------------*/
typedef struct _RESTART_DYNSTRUCT
{
int                          step;     
struct _STRUCT_DYNAMIC       sdyn;
struct _STRUCT_DYN_CALC      dynvar;

int                          dist_vec_rhs[10];
int                          dist_vec_sol[10];
int                          dist_vec_dispi[10];
int                          dist_vec_vel[10];
int                          dist_vec_acc[10];
int                          dist_vec_fie[10];
int                          dist_vec_work[10];

int                          intforce;
int                          dirich;

int                          handle_of_node_handles;
struct _ARRAY                node_handles;

int                          handle_of_ele_handles;
struct _ARRAY                ele_handles;

} RESTART_DYNSTRUCT;
