/*!---------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/

#ifndef GLOBAL_PROTOTYPES_H
#define GLOBAL_PROTOTYPES_H



/*----------------------------------------------------------------------*
  |  main_ccarat.c                                      m.gee 11/01    |
 *----------------------------------------------------------------------*/
INT main(
    INT            argc,
    char          *argv[]);


/*----------------------------------------------------------------------*
  |  global_ass_dof.c                                   m.gee 11/01    |
 *----------------------------------------------------------------------*/
void assign_dof(
    FIELD         *actfield);

void init_dof_discretization(
  DISCRET       *actdis
  );

void assign_dof_discretization(
  DISCRET       *actdis
  );

void init_dof_discretization_xfem(
  DISCRET       *actdis
  );

void assign_dof_discretization_xfem(
  DISCRET       *actdis
  );

/*----------------------------------------------------------------------*
  |  global_ass_dof_ndis.c                               genk 08/02    |
 *----------------------------------------------------------------------*/
void assign_dof_ndis(
    FIELD         *actfield);


/*----------------------------------------------------------------------*
  |  global_cal_control.c                               m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntacal(void);


/*----------------------------------------------------------------------*
  |  global_cal_control.c                               genk 10/03     |
 *----------------------------------------------------------------------*/
void global_result_test(void);


/*----------------------------------------------------------------------*
  |  global_dyn_control.c                               m.gee 11/01    |
 *----------------------------------------------------------------------*/
void caldyn(void);


/*----------------------------------------------------------------------*
  | global_control.c                                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntam(
    INT            argc,
    char          *argv[]);


/*----------------------------------------------------------------------*
  | global_init_control.c                               m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntaini(
    INT            argc,
    char          *argv[]);


/*----------------------------------------------------------------------*
  | global_inp_control.c                                m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntainp(void);


/*----------------------------------------------------------------------*
  | global_locsys.c                                       genk 01/03    |
 *----------------------------------------------------------------------*/
void locsys_inherit_to_node(void);
void locsys_trans(ELEMENT *ele, DOUBLE **estif1, DOUBLE **estif2,
                                DOUBLE *vec1,    DOUBLE *vec2);
void locsys_trans_sol(FIELD *actfield, INT idis, INT array,
                      INT place, INT flag);
void locsys_trans_sol_dirich(FIELD *actfield, INT idis, INT array,
                             INT place, INT flag);
void locsys_trans_nodval(ELEMENT *actele, DOUBLE *val, INT numdf,
                         INT iloccsys, INT flag);

/*----------------------------------------------------------------------*
  | global_monitoring.c                                  genk 01/03    |
 *----------------------------------------------------------------------*/
void monitoring(
                  FIELD         *actfield,
                  INT            numf,
                  INT            actpos,
                  DOUBLE         time
               );


/*----------------------------------------------------------------------*
  | global_open_files.c                                 m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntadev(
    INT            argc,
    char          *argv[]);


/*----------------------------------------------------------------------*
  |  global_node_find.c                                 m.gee 11/01    |
 *----------------------------------------------------------------------*/
void iscouple_find_node_comp(
    NODE          *actnode,
    FIELD         *searchfield,
    NODE         **partnernode,
    INT            coupleID,
    INT            dof);

void iscouple_find_node_comp_discretization(
    NODE          *actnode,
    DISCRET       *searchdis,
    NODE         **partnernode,
    INT            coupleID,
    INT            dof);

void cheque_distance(
    DOUBLE        *x1,
    DOUBLE        *x2,
    DOUBLE         tol,
    INT           *ierr);

void find_assign_coupset(
    FIELD         *actfield,
    INT            coupleID,
    INT           *counter);

void find_assign_coupset_discretization(
    DISCRET       *actdis,
    INT            coupleID,
    INT           *counter);

/*----------------------------------------------------------------------*
  | global_timecurve.c                                  m.gee 02/02    |
 *----------------------------------------------------------------------*/
void dyn_init_curve(
    INT            actcurve,
    INT            nstep,
    DOUBLE         dt,
    DOUBLE         maxtime);

void dyn_facfromcurve(
    INT            actcurve,
    DOUBLE         T,
    DOUBLE        *fac);

DOUBLE dyn_facexplcurve(
    INT            actcurve,
    DOUBLE         T);


/*----------------------------------------------------------------------*
  | global_calelm.c                                     m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calelm(
    FIELD         *actfield,     /* active field */
    SOLVAR        *actsolv,      /* active SOLVAR */
    PARTITION     *actpart,      /* my partition of this field */
    INTRA         *actintra,     /* my intra-communicator */
    INT            sysarray1,    /* number of first sparse system matrix */
    INT            sysarray2,    /* number of secnd system matrix, if present, else -1 */
    CONTAINER     *container,    /*!< contains variables defined in container.h */
    CALC_ACTION   *action);       /* calculation option passed to element routines */

void calinit(
    FIELD         *actfield,   /* the active physical field */
    PARTITION     *actpart,    /* my partition of this field */
    CALC_ACTION   *action,
    CONTAINER     *container); /*!< contains variables defined in container.h */

void calreduce(
    FIELD         *actfield,    /* the active field */
    PARTITION     *actpart,     /* my partition of this field */
    INTRA         *actintra,    /* the field's intra-communicator */
    CALC_ACTION   *action,      /* action for element routines */
    CONTAINER     *container);  /*!< contains variables defined in container.h */


/*----------------------------------------------------------------------*
  | global_calrhs.c                                     m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calrhs(
    FIELD         *actfield,     /* the active field */
    SOLVAR        *actsolv,      /* the active SOLVAR */
    PARTITION     *actpart,      /* my partition of this field */
    INTRA         *actintra,     /* the field's intra-communicator */
    INT            actsysarray,  /* the active sparse array */
    DIST_VECTOR   *rhs1,         /* 2 dist. vectors for rhs */
    CALC_ACTION   *action,       /* action to be passed to element routines */
    CONTAINER     *container);   /*!< contains variables defined in container.h */

void rhs_point_neum(
    DOUBLE        *rhs,
    INT            dimrhs,
    PARTITION     *actpart);


/*----------------------------------------------------------------------*
  | global_check_max.c                                        mn 04/04 |
 *----------------------------------------------------------------------*/
#ifdef CHECK_MAX
void check_max_sizes(void);
#endif




#endif

