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

#include "standardtypes.h"
#include "solution.h"

/*----------------------------------------------------------------------*
 |  main_ccarat.c                                        m.gee 11/01    |
 *----------------------------------------------------------------------*/
INT main(INT argc, char *argv[]);
/*----------------------------------------------------------------------*
 |  global_ass_dof.c                                     m.gee 11/01    |
 *----------------------------------------------------------------------*/
void assign_dof(FIELD *actfield);
/*----------------------------------------------------------------------*
 |  global_ass_dof_ndis.c                                 genk 08/02    |
 *----------------------------------------------------------------------*/
void assign_dof_ndis(FIELD *actfield);
/*----------------------------------------------------------------------*
 |  global_cal_control.c                                 m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntacal(void);
/*----------------------------------------------------------------------*
 |  global_cal_control.c                                 genk 10/03     |
 *----------------------------------------------------------------------*/
void global_result_test(void); 
/*----------------------------------------------------------------------*
 |  cal_dyn_control.c                                    m.gee 11/01    |
 *----------------------------------------------------------------------*/
void caldyn(void);
/*----------------------------------------------------------------------*
 | cal_static_control.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calsta(void);
void stalin(void);
/*----------------------------------------------------------------------*
 | cal_static_service.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calstatserv_findcontroldof(FIELD     *actfield,
                                INT        control_node_global,
                                INT        control_dof,
                                NODE     **node,
                                INT       *cdof); 
/*----------------------------------------------------------------------*
 | cal_static_service.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void calstatserv_findreldofs(FIELD     *actfield,
                             INT       *reldisnode_ID,
                             INT       *reldis_dof,
                             INT        num_reldis,
                             INT       *reldof); 
/*----------------------------------------------------------------------*
 | cal_static_service.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void get_stepsize(INT         kstep,
                  STATIC_VAR *statvar); 
/*----------------------------------------------------------------------*
 | global_control.c                                      m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntam(INT argc, char *argv[]);
/*----------------------------------------------------------------------*
 | global_visual.c                                       genk  07/02    |
 *----------------------------------------------------------------------*/
/*!---------------------------------------------------------------------                                         
\brief call of Visualisation tools

<pre>                                                         genk 07/02       

This routine checks the type of problem and based on the program  options
a visualisation tool is called.
At the moment implemented:
VISUAL2

</pre>  
\return void                                                                       

------------------------------------------------------------------------*/
void ntavisual(void); 
/*----------------------------------------------------------------------*
 | global_init_control.c                                 m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntaini(INT argc, char *argv[]);
/*----------------------------------------------------------------------*
 | global_inp_control.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntainp(void);
/*----------------------------------------------------------------------*
 | global_mask_matrices.c                                m.gee 11/01    |
 *----------------------------------------------------------------------*/
void mask_global_matrices(void);
/*----------------------------------------------------------------------*
 | global_monitoring.c                                    genk 01/03    |
 *----------------------------------------------------------------------*/
void monitoring(FIELD *actfield,INT numf, INT actpos, INT actstep, DOUBLE time); 
/*----------------------------------------------------------------------*
 |  machine_hpux.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntadev(INT argc, char *argv[]);
/*----------------------------------------------------------------------*
 |  map_node_find.c                                      m.gee 11/01    |
 *----------------------------------------------------------------------*/
void iscouple_find_node_comp(NODE  *actnode, 
                                FIELD *searchfield, 
                                NODE **partnernode,
                                INT    coupleID,
                                INT    dof);
void cheque_distance(DOUBLE *x1, DOUBLE *x2, DOUBLE tol, INT *ierr);
void find_assign_coupset(FIELD *actfield, 
                            INT    coupleID, 
                            INT   *counter);
/*----------------------------------------------------------------------*
 |  dyn_timecurve.c                                      m.gee 02/02    |
 *----------------------------------------------------------------------*/
void dyn_init_curve(INT actcurve,
                   INT    nstep,
                   DOUBLE dt,
                   DOUBLE maxtime);
void dyn_facfromcurve(INT actcurve,
                   DOUBLE T,
                   DOUBLE *fac);
DOUBLE dyn_facexplcurve(INT actcurve,
                      DOUBLE T);		   
/*----------------------------------------------------------------------*
 |  out_global.c                                         m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_general(void);
void out_sol(FIELD *actfield, PARTITION *actpart, INTRA *actintra, 
             INT step, INT place);

void out_fluidmf(FIELD *fluidfield);
void out_fsi(FIELD *fluidfield);
void out_fluidtu(FIELD *actfield, INTRA *actintra, INT step, INT place);
/*----------------------------------------------------------------------*
 |  out_gid_sol.c                                        m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_sol_init(void);
void out_gid_domains(FIELD *actfield);
void out_gid_sol(char string[], FIELD *actfield, INTRA  *actintra, INT step,
                 INT place, DOUBLE time);
void out_fsi(FIELD *fluidfield);
void out_gid_sol_fsi(FIELD *fluidfield, FIELD *structfield);
/*----------------------------------------------------------------------*
 |  out_gid_soldyn.c                                     m.gee 5/03     |
 *----------------------------------------------------------------------*/
void out_gid_soldyn(char string[], FIELD *actfield, INTRA  *actintra, INT step,
                   INT place, DOUBLE totaltime);
/*----------------------------------------------------------------------*
 |  out_gid_msh.c                                        m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_msh(void);
void out_gid_allcoords(FILE *out);
/*----------------------------------------------------------------------*
 |  out_monitor.c                                         genk 01/03    |
 *----------------------------------------------------------------------*/
void out_monitor(FIELD *actfield, INT numf);
void out_area(ARRAY totarea_a); 

/*----------------------------------------------------------------------*
 |  out_checkfilesize.c                                   genk 08/03    |
 *----------------------------------------------------------------------*/
void out_checkfilesize(INT opt); 
/*----------------------------------------------------------------------*
 |  out_plt.c                                            chfoe 01/04    |
 *----------------------------------------------------------------------*/
void plot_liftdrag(DOUBLE time, DOUBLE *liftdrag);
void plot_lte(	DOUBLE  time, 
                INT     step, 
                DOUBLE  norm, 
                DOUBLE  dt, 
                INT     itnum);
void plot_ale_quality(FIELD *field,INT step, INTRA *actintra, 
                      PARTITION *actpart);
/*----------------------------------------------------------------------*
 |  inherit_insidedesign.c                                  m.gee 3/02  |
 *----------------------------------------------------------------------*/
void inherit_dirich_coup_indesign(void);
/*----------------------------------------------------------------------*
 |  inherit_design_dis.c                                    m.gee 3/02  |
 *----------------------------------------------------------------------*/
void inherit_design_dis_dirichlet(DISCRET *actdis);
void inherit_design_dis_couple(DISCRET *actdis);
void inherit_design_dis_fsicouple(DISCRET *actdis);
void inherit_design_dis_freesurf(DISCRET *actdis);
void inherit_design_dis_neum(DISCRET *actdis);
/*----------------------------------------------------------------------*
 |  inherit_design_ele.c                                   chfoe 01/04  |
 *----------------------------------------------------------------------*/
void inherit_design_ele(DISCRET *actdis);
/*----------------------------------------------------------------------*
 |  input_conditions.c                                  m.gee 11/01     |
 *----------------------------------------------------------------------*/
void inp_conditions(void);
/*----------------------------------------------------------------------*
 |  input_control_global.c                                  m.gee 11/01 |
 *----------------------------------------------------------------------*/
void inpctr(void);
void inpctrprob(void);
void inpctrdyn(void);
void inpctrstat(void);
void inpctreig(void);
void inpctr_dyn_struct(STRUCT_DYNAMIC *sdyn);
void inpctr_dyn_ale(ALE_DYNAMIC *adyn);
void inpctr_eig_struct(ALLEIG *alleig);
/*!---------------------------------------------------------------------                                         
\brief input of the FLUID DYNAMIC block in the input-file

<pre>                                                         genk 03/02

In this routine the data in the FLUID DYNAMIC block of the input file
are read and stored in fdyn	       

</pre>
\param  *fdyn 	  FLUID_DATA       (o)	   
\return void                                                                       

------------------------------------------------------------------------*/
void inpctr_dyn_fluid(FLUID_DYNAMIC *fdyn);
/*!---------------------------------------------------------------------                                         
\brief input of the FSI DYNAMIC block in the input-file

<pre>                                                         genk 09/02

In this routine the data in the FSI DYNAMIC block of the input file
are read and stored in fsidyn	       

</pre>
\param  *fsidyn 	  FSI_DATA       (o)	   
\return void                                                                       

------------------------------------------------------------------------*/
void inpctr_dyn_fsi(FSI_DYNAMIC *fsidyn);
/*----------------------------------------------------------------------*
 |  input_ctr_head.c                                  m.gee 11/01       |
 *----------------------------------------------------------------------*/
void inpctrhed(void);
void inptrace(void);
/*----------------------------------------------------------------------*
 |  input_curves.c                                  m.gee 11/01         |
 *----------------------------------------------------------------------*/
void inp_cond_curve(void);
void inp_read_curve(char *string);

/*----------------------------------------------------------------------*
    input_funct.c
 *----------------------------------------------------------------------*/
void inp_cond_funct(void);


/*----------------------------------------------------------------------*
 |  input_design.c                                  m.gee 11/01         |
 *----------------------------------------------------------------------*/
void inpdesign(void);
void inp_dnode(void);
void read_1_dnode(DNODE *dnode, INT readId);
void inp_dline(void);
void read_1_dline(DLINE *dline, INT readId);
void inp_dsurface(void);
void read_1_dsurf(DSURF *dsurf, INT readId);
void inp_dvolume(void);
void read_1_dvol(DVOL *dvol, INT readId);
void inp_designsize(void);
/*----------------------------------------------------------------------*
 |  input_design_top.c                                  m.gee 11/01     |
 *---------------------------------------------------------------------*/
void inpdesign_topology_design(void);
void inpdesign_topology_fe(void);
/*----------------------------------------------------------------------*
 |  input_material.c                                  m.gee 11/01       |
 *----------------------------------------------------------------------*/
void inp_material(void);
void inp_multimat(void);
/*----------------------------------------------------------------------*
 |  input_mesh.c                                  m.gee 11/01           |
 *----------------------------------------------------------------------*/
void inpfield(void);
void inp_assign_nodes(DISCRET *actdis);
void inpdis(FIELD *actfield);
void inpnodes(void);
void inp_struct_field(FIELD *structfield);
void inp_fluid_field(FIELD *fluidfield);
void inp_ale_field(FIELD *alefield);
/*----------------------------------------------------------------------*
 |  input_monitor                                         genk 01/03    |
 *----------------------------------------------------------------------*/
void inp_monitor(void);
/*----------------------------------------------------------------------*
 |  input_resultdescr                                       uk 05/04    |
 *----------------------------------------------------------------------*/
void inp_resultdescr(void);
/*----------------------------------------------------------------------*
 |  input_topology.c                                  m.gee 11/01       |
 *----------------------------------------------------------------------*/
void inp_topology(DISCRET *actdis);
void inp_detailed_topology(DISCRET   *actdis);
/*----------------------------------------------------------------------*
 |  math1.c                                               m.gee 11/01   |
 *----------------------------------------------------------------------*/
void math_array_copy(DOUBLE **from, INT n, INT m, DOUBLE **to);
void math_inv3(DOUBLE **a, DOUBLE *det);
void math_tran(DOUBLE **a, INT n);
void math_unvc(DOUBLE *enorm,DOUBLE *vec, INT n);
void math_matvecdense(DOUBLE  *r,
                         DOUBLE **A,
                         DOUBLE  *b,
                         INT      ni,
                         INT      nk,
                         INT      init,
                         DOUBLE   factor);
void math_mattrnvecdense(DOUBLE  *r,
                         DOUBLE **A,
                         DOUBLE  *b,
                         INT      ni,
                         INT      nk,
                         INT      init,
                         DOUBLE   factor);
void math_matmatdense(DOUBLE **R,
                         DOUBLE **A,
                         DOUBLE **B,
                         INT      ni,
                         INT      nk,
                         INT      nj,
                         INT      init,
                         DOUBLE   factor);
void math_mattrnmatdense(DOUBLE **R,
                            DOUBLE **A,
                            DOUBLE **B,
                            INT      ni,
                            INT      nk,
                            INT      nj,
                            INT      init,
                            DOUBLE   factor);
void math_matmattrndense(DOUBLE **R,
                            DOUBLE **A,
                            DOUBLE **B,
                            INT      ni,
                            INT      nk,
                            INT      nj,
                            INT      init,
                            DOUBLE   factor);
void math_sym_inv(DOUBLE **A, INT dim);
void math_unsym_inv(DOUBLE **A, INT dimr, INT dimc);
void math_sppr(DOUBLE *spat, DOUBLE *a, DOUBLE *b, DOUBLE *c);
void math_addab(DOUBLE **a, DOUBLE **b, INT dim1, INT dim2, DOUBLE fact);
/*!---------------------------------------------------------------------                                         
\brief extract digits from integer number

<pre>                                                         genk 04/02		     
</pre>   
\param  num	 INT   (i)    integer number
\param *it	 INT   (o)    integer on position "thousand"
\param *ih       INT   (o)    integer on position "hundred"
\param *id       INT   (o)    integer on position "ten"
\param *id       INT   (o)    integer on position "one"
\return void 

------------------------------------------------------------------------*/
void math_intextract(
                    INT num,    
                    INT *it,    
		    INT *ih,    
		    INT *id,    
		    INT *io     
	            );

/*----------------------------------------------------------------------*
 |  sort_find.c                                          m.gee 11/01    |
 *----------------------------------------------------------------------*/
INT cmp_int(const void *a, const void *b );
DOUBLE cmp_double(const void *a, const void *b );
void mg_sort(INT list[], INT N, INT list2[], DOUBLE list3[]);
INT quick_find(INT key, INT list[], INT length, INT shift, INT bins[]);
void init_quick_find(INT list[], INT length, INT *shift, INT *bins);
INT find_index(INT key, INT list[], INT length);
/*----------------------------------------------------------------------*
 |  par_assignmesh.c                                     m.gee 11/01    |
 *----------------------------------------------------------------------*/
void part_assignfield(void);
/*----------------------------------------------------------------------*
 |  par_initmetis.c                                      m.gee 11/01    |
 *----------------------------------------------------------------------*/
void part_fields(void);
/*----------------------------------------------------------------------*
 |  par_make_comm.c                                      m.gee 11/01    |
 *----------------------------------------------------------------------*/
void create_communicators(void);
/*----------------------------------------------------------------------*
 |  pss_am.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | shift a pointer by a certain distance                      m.gee 3/02|
 | NOTE:                                                                |
 | This routine is usefull when data is moved in the memory (e.g.) by a |
 | call to REALLOC, and there are pointers which used to point to this  |
 | data. After moving the data these pointers do no longer point to the |
 | correct adresses, but they can be moved to the new location of the   |
 | data by using this function.                                         |
 | The correct functionality of this routine is dependent on the proper |
 | definition of the data type PTRSIZE amde in definitions.h:           |
 |#ifdef SIXTYFOUR                a 64 bit pointer is of size long int  |
 |typedef long int PTRSIZE;                                             |
 |#else                                a 32 bit pointer is of size INT  |
 |typedef INT PTRSIZE;                                                  |
 |#endif                                                                |
 | This function takes as argument                                      |
 | - the adress of the pointer to be moved                              |
 | - the relative difference between the new and the old adress of the  |
 |   data which was pointed to                                          |
 | EXAMPLE:                                                             |
 | i=5;                                                                 |
 | j=12;                                                                |
 | iptr = &i;                                                           |
 | now move iptr from i to j:                                           |
 | ShiftPointer((void**)&iptr,(PTRSIZE)&j-(PTRSIZE)&i);                 |
 | result:                                                              |
 | *ptr = 12 now                                                        |
 | WARNING:                                                             |
 | This is an extremely power and helpfull little routine, but when using
 | this you should be REALLY sure what you are doing, 'cause misuse     |
 | leads to a very nasty type of error, which can be hardly detected by |
 | debugging!                                                           |
 |                                                                      |
 |                                                                      |
 | This routine appears courtesy of                                     |
 | Dr. Ralf Diekmann, Hilti AG, Fuerstentum Liechtenstein               |
 *----------------------------------------------------------------------*/
void ShiftPointer(void **ptr, PTRSIZE diff);
/*----------------------------------------------------------------------*
 | redefinition of malloc DEBUG version                       m.gee 2/02|
 | bhaves exactly like malloc conform to ansi c standard                |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | redefinition of malloc FAST version                        m.gee 2/02|
 | bhaves exactly like malloc conform to ansi c standard                |
 *----------------------------------------------------------------------*/
void *CCAMALLOC(INT size);
/*----------------------------------------------------------------------*
 | redefinition of calloc DEBUG version                   m.gee 2/02    |
 | bhaves exactly like calloc conform to ansi c standard                |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | redefinition of calloc FAST version                    m.gee 2/02    |
 | bhaves exactly like calloc conform to ansi c standard                |
 *----------------------------------------------------------------------*/
void *CCACALLOC(INT num, INT size);
/*----------------------------------------------------------------------*
 | redefinition of realloc DEBUG version                  m.gee 2/02    |
 | bhaves exactly like realloc conform to ansi c standard               |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | redefinition of realloc FAST version                  m.gee 2/02     |
 | bhaves exactly like realloc conform to ansi c standard               |
 *----------------------------------------------------------------------*/
void *CCAREALLOC(void *oldptr, INT size);
/*----------------------------------------------------------------------*
 | redefinition of free DEBUG version                     m.gee 2/02    |
 | bhaves exactly like free conform to ansi c standard                  |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | redefinition of free FAST version                      m.gee 2/02    |
 | bhaves exactly like free conform to ansi c standard                  |
 *----------------------------------------------------------------------*/
void *CCAFREE(void *oldptr);
/*----------------------------------------------------------------------*
 | define array                                           m.gee 8/00    |
 | allocate a 1 or 2 - D vector of type INT or DOUBLE in the structure  |
 | ARRAY *a. See also am.h for the structure ARRAY                      |
 | char *namstr  (input)   name of array                                |
 | ARRAY *a      (input)   adress of structure ARRAY the vector lives in|
 | INT fdim      (input)   first dimension of 2D vector                 |
 |                         dimension of 1D vector                       |
 | INT sdim      (input)   scnd dimension of 2D vector                  |
 | char typstr[] (input)   type of array to allocate                    |
 |               ="IV"     allocate integer vector in a->a.iv           |
 |               ="IA"     allocate DOUBLE array  in a->a.ia            |
 |               ="DV"     allocate integer vector in a->a.dv           |
 |               ="DA"     allocate DOUBLE array  in a->a.da            |
 | return value:                                                        |
 | void pointer to allocated memory                                     |
 *----------------------------------------------------------------------*/
void* amdef(char *namstr,ARRAY *a,INT fdim, INT sdim, char typstr[]);
/*----------------------------------------------------------------------*
 | redefine array                                         m.gee 8/00    |
 | changes an already allocated array a in dimensions                   |
 | a typecast of the values in the array is not possible                |
 | (no INT to DOUBLE or vice versa transformation)                      |
 |                                                                      |
 | if the new dimension of an the array is larger then the old one,     |
 | the values inside the array are kept, the new entries are initialized|
 | with zero                                                            |
 |                                                                      |
 | if the new dimension of the  array is smaller then the old one,      |
 | values outside the new dimension are dropped                         |
 |                                                                      |
 | this routine handles every dimension (fdim and sdim) separately,     |
 | that means, it does NOT behave like a fortran array                  |
 |                                                                      |
 | a cast from iv to ia and from dv to da and vice versa is allowed     |
 |                                                                      |
 | ARRAY *a      (input)   adress of structure ARRAY the vector lives in|
 | INT newfdim   (input)   new first dimension of 2D vector             |
 |                         new dimension of 1D vector                   |
 | INT newsdim   (input)   new scnd dimension of 2D vector              |
 | char newtypstr[] (input)   type of array to allocate                 |
 |               ="IV"     allocate integer vector in a->a.iv           |
 |               ="IA"     allocate DOUBLE array  in a->a.ia            |
 |               ="DV"     allocate integer vector in a->a.dv           |
 |               ="DA"     allocate DOUBLE array  in a->a.da            |
 | return value:                                                        |
 | void pointer to allocated memory                                     |
 *----------------------------------------------------------------------*/
void* amredef(ARRAY *a,INT newfdim, INT newsdim, char newtypstr[]);
/*----------------------------------------------------------------------*
 | delete         array                                   m.gee 8/00    |
 | frees the vector or array located in the ARRAY                       |
 | *array (input) the adress of the structure holding the vector/array  |
 *----------------------------------------------------------------------*/
void  amdel(ARRAY *array);
/*----------------------------------------------------------------------*
 | initialize an array by zero                            m.gee 8/00    |
 | initializes the content of the ARRAY array to zero                   |
 | put 0 to integer fields, 0.0 to DOUBLE fields                        |
 | ARRAY *array (input) adress of the ARRAY array                       |
 *----------------------------------------------------------------------*/
void  amzero(ARRAY *array);
/*----------------------------------------------------------------------*
 | multiply an array by a scalar                          m.gee 6/01    |
 | scales the contents of a field by a given value                      |
 | ARRAY *array (input) adress of the ARRAY array                       |
 | *value (input) adress of the scaling parameter, this may be of type  |
 |                INT* or DOUBLE* and must be casted to void* in the    |
 |                parameter list                                        |
 |                example: amscal(&val,(void*)(&ione));                 |        
 *----------------------------------------------------------------------*/
void  amscal(ARRAY *array, void *value);
/*----------------------------------------------------------------------*
 | initialize an array by value                           m.gee 6/01    |
 | inits the contents of a field by a given value                       |
 | ARRAY *array (input) adress of the ARRAY array                       |
 | *value (input) adress of the initvalues, this may be of type         |
 |                INT* or DOUBLE* and must be casted to void* in the    |
 |                parameter list                                        |
 |                example: aminit(&val,(void*)(&ione));                 |        
 *----------------------------------------------------------------------*/
void  aminit(ARRAY *array, void *value);
/*----------------------------------------------------------------------*
 | allocate and make a copy of ARRAY                      m.gee 8/00    |
 | alocates a new field in array_to of the same size and type as        |
 | in array_from. Then copies the contents from array_from to array_to  |
 | ARRAY *array_from (input) the adress of the field one wants to copy  |
 | ARRAY *array_to   (in/out) the adress of the structure where the new |
 |                            field is allocated in                     |
 | return value: adress of the new field                                |
 *----------------------------------------------------------------------*/
void* am_alloc_copy(ARRAY *array_from, ARRAY *array_to);
/*----------------------------------------------------------------------*
 | make a copy of ARRAY,                                  m.gee 6/01    |
 | Copies the contents from array_from to array_to                      |
 | ARRAY *array_from (input) the adress of the field one wants to copy  |
 | ARRAY *array_to   (in/out) the adress of the field values are copied |
 |                            to                                        |
 | user must provide fields of matching type and size!                  |
 | return value: array_to->a.iv                                         |
 *----------------------------------------------------------------------*/
void* amcopy(ARRAY *array_from, ARRAY *array_to);
/*----------------------------------------------------------------------*
 |                                                        m.gee 2/02    |
 | makes array_to += array_from * factor                                |
 | if init==1 array_to is initialized to zero                           |
 | user must provide matching array-types and sufficient space in       |
 | array_to                                                             |
 | ARRAY *array_to   (output) adress of structure values are added to   |
 | ARRAY *array_from (input)  adress of structure values are taken from |
 | DOUBLE factor     (input)  scaling factor, must be casted to         |
 |                            DOUBLE in the call to this routine,       |
 |                            but also operates on integer fields       |
 |                            (this is a bit dirty I know....)          |
 | INT init          (input)  flag                                      |
 | ==1 array_to is initialized to zero                                  |
 | else values are assembled to array_to                                |
 *----------------------------------------------------------------------*/
void  amadd(ARRAY *array_to, ARRAY *array_from, DOUBLE factor, INT init);
/*----------------------------------------------------------------------*
 | define 4D array                                       m.gee 12/01    |
 | similar to amdef, but for 3D and 4D fields                           |
 | the field is NOT initialized                                         |
 |                                                                      |
 | char *namstr  (input)   name of array                                |
 | ARRAY4D *a    (input)   adress of structure ARRAY4D                  |
 | INT fdim      (input)   first dimension of 3D or 4D array            |
 | INT sdim      (input)   scnd dimension of 3D or 4D array             |
 | INT tdim      (input)   third dimension of 3D or 4D array            |
 | INT fodim     (input)   fourth dimension of 4D array,                |
 |                         ==0 of array is 3D                           |
 | char typstr[] (input)   type of field to allocate                    |
 |                         ="I3" 3D integer field                       |
 |                         ="I4" 4D integer field                       |
 |                         ="D3" 3D DOUBLE  field                       |
 |                         ="D4" 4D DOUBLE  field                       |
 | return value:                                                        |
 | void pointer to allocated memory                                     |
 *----------------------------------------------------------------------*/
void* am4def(char *namstr, ARRAY4D *a, INT fdim, INT sdim, INT tdim, 
             INT fodim, char typstr[]);
/*----------------------------------------------------------------------*
 | delete 4dimensional array                             m.gee 12/01    |
 | frees all field memory in array                                      |
 | ARRAY4D *array    (input)                adress of array-structure   |
 | no return value                                                      |
 *----------------------------------------------------------------------*/
void  am4del(ARRAY4D *array);
/*----------------------------------------------------------------------*
 | initialize an 4D array by zero                           m.gee 12/01 |
 | see head of amzero                                                   |
 *----------------------------------------------------------------------*/
void  am4zero(ARRAY4D *array);
/*----------------------------------------------------------------------*
 | initialize a 4D array by value                         m.gee 6/01    |
 | see head of aminit                                                   |
 *----------------------------------------------------------------------*/
void  am4init(ARRAY4D *array, void *value);
/*----------------------------------------------------------------------*
 | allocate and make a copy of ARRAY4D                     m.gee 12/01  |
 | see head of am_alloc_copy                                            |
 *----------------------------------------------------------------------*/
void* am4_alloc_copy(ARRAY4D *array_from, ARRAY4D *array_to);
/*----------------------------------------------------------------------*
 | make a copy of ARRAY4D                                  m.gee 12/01  |
 | the given arrays must be of equal type and size                      |
 | see head of amcopy                                                   |
 *----------------------------------------------------------------------*/
void* am4copy(ARRAY4D *array_from, ARRAY4D *array_to);
/*----------------------------------------------------------------------*
 | redefine of ARRAY4D                                     m.gee 12/01  |
 | NOTE:  - a dimension cast like it is possible with ARRAY doesn't work|
 |        - 4-dimensional arrays with fodim=0 are not allowed           |
 |        - remember, this routine can be VERY expensive !              |
 |        - aditional space in the redefined array is set to zero       |
 |          (unlike in am4def which does NOT initialize)                |
 |        usage similar to amredef                                      |
 *----------------------------------------------------------------------*/
void* am4redef(ARRAY4D *array, INT newfdim, INT newsdim, INT newtdim, 
               INT newfodim);
/*----------------------------------------------------------------------*
 |  pss_ds.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | report the amount of actual allocated memory           m.gee 2/02    |
 | only works if DEBUG is defined                                       |
 *----------------------------------------------------------------------*/
void dsmemreport(void);
/*----------------------------------------------------------------------*
 |                                                        m.gee 3/02    |
 | this routine does nothing if the boolean criterium is true           |
 | but aborts the programm, if it is not                                |
 | The routine is empty in an optimezed (not DEBUG) compilation and     |
 | can therefor be excessively used to develop a secure code, without   |
 | making it slow when running as fast-exe                              |
 | parameter list                                                       |
 | test - (input) boolean criterium                                     |
 | string - (input) error message                                       |
 *----------------------------------------------------------------------*/
void dsassert(INT test, char string[]);
/*----------------------------------------------------------------------*
 | report an error and stop program                       m.gee 8/00    |
 | prints error message string to console and *.err                     |
 | prints call tree, if DEBUG was defined                               |
 | aborts parallel and sequentiell programm                             |
 *----------------------------------------------------------------------*/
void dserror(char string[], ...);
/*----------------------------------------------------------------------*
 | collects warnings during running process                  ck 07/03   |
 | and writes them to the screen at the end                             |
 *----------------------------------------------------------------------*/
void dswarning(INT task, INT warning);
/*----------------------------------------------------------------------*
 | Initialize bugtracing systems                          m.gee 8/00    |
 *----------------------------------------------------------------------*/
void dsinit(void);
/*----------------------------------------------------------------------*
 | report intrance to routine named string to the tracing system  m.gee |
 | is empty is DEBUG is not defined                               8/00  |
 *----------------------------------------------------------------------*/
void dstrc_enter(char string[]);
/*----------------------------------------------------------------------*
 | report exit to routine to the tracing system           m.gee 8/00    |
 | is empty is DEBUG is not defined                                     |
 *----------------------------------------------------------------------*/
void dstrc_exit(void);
/*----------------------------------------------------------------------*
 | report a new DOUBLE array to the bugtracing system     m.gee 8/00    |    
 | this routine is called by the am-system only !                       | 
 *----------------------------------------------------------------------*/
void dsreportarray(void *array, INT typ);
/*----------------------------------------------------------------------*
 | report a new DOUBLE array to the bugtracing system     m.gee 8/00    |    
 | this routine is called by the am-system only !                       |
 *----------------------------------------------------------------------*/
void dsdeletearray(void *array, INT typ);
/*----------------------------------------------------------------------*
 | write a report about all arrays to the .err file       m.gee 8/00    |    
 | does nothing if DEBUG is not defined                                 |
 | writes a list of all ARRAY and ARRAY4D structure generated           |
 | by the am-System to the *.err file                                   |
 *----------------------------------------------------------------------*/
void dstrace_to_err(void);
/*----------------------------------------------------------------------*
 | routine to initialise the cpu - time                  genk 05/02     |
 *----------------------------------------------------------------------*/
void ds_cputime_init(void);
/*----------------------------------------------------------------------*
 | routine to meassure the cpu - time                    genk 05/02     |
 *----------------------------------------------------------------------*/
DOUBLE ds_cputime(void);
/*----------------------------------------------------------------------*
 |  pss_fr.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | init file reading system                               m.gee 8/00    |
 *----------------------------------------------------------------------*/
void frinit(void);
/*----------------------------------------------------------------------*
 | rewind the copy of the input file                      m.gee 8/00    |
 *----------------------------------------------------------------------*/
void frrewind(void);
/*----------------------------------------------------------------------*
 | find a character string                                m.gee 8/00    |
 | char string[] (input) character string to search for copy of inputfil|
 |                       terminates programm if not found               |
 *----------------------------------------------------------------------*/
INT frfind(char string[]);
/*----------------------------------------------------------------------*
 |                                                        m.gee 8/00    |
 | sets a pointer to the next line in thecopy of input file on all procs|
 *----------------------------------------------------------------------*/
void frread(void);
/*----------------------------------------------------------------------*
 | reads n integers from input_file                       m.gee 4/01    |
 | char string[] (input) keyword to search for in actual line           |
 | INT *var      (output) adress of field to hold values read           |
 | INT num       (input)  number of values to read                      |
 | INT *ierr     (output) =0 keyword not found / =1 values read         |
 *----------------------------------------------------------------------*/
void frint_n(char string[],INT *var,INT num, INT *ierr);
/*----------------------------------------------------------------------*
 | reads an integer from input_file                       m.gee 8/00    |
 | starting in allfiles.actrow                                          |
 | char string[] (input) keyword to search for in actual line           |
 | INT *var      (output)adress of variable to read to                  |
 | INT *ierr     (output) flag to indicate success                      |
 | ierr=0 keyword not found                                             |
 | ierr=1 integer read                                                  |
 *----------------------------------------------------------------------*/
void frint(char string[],INT *var, INT *ierr);
/*----------------------------------------------------------------------*
 | reads a DOUBLE from input file line allfiles.line      m.gee 8/00    |
 | see frint_n                                                          |
 | ierr=0 kenner not found                                              |
 | ierr=1 DOUBLE  read                                                  |
 *----------------------------------------------------------------------*/
void frdouble_n(char string[],DOUBLE *var,INT num, INT *ierr);
/*----------------------------------------------------------------------*
 | reads a DOUBLE from input file line allfiles.line      m.gee 8/00    |
 | see frint                                                            |
 | ierr=0 kenner not found                                              |
 | ierr=1 DOUBLE  read                                                  |
 *----------------------------------------------------------------------*/
void frdouble(char string[],DOUBLE *var, INT *ierr);
/*----------------------------------------------------------------------*
 | reads a charstring from input file                                   |
 | user must assure, that the given charpointer space is long enough    |   
 | to hold the string                                                   |
 | ierr=0 keyword not found on line                                     |
 | ierr=1 char string read                                 m.gee 8/00   |
 *----------------------------------------------------------------------*/
void frchar(char string[],char *var, INT *ierr);
/*----------------------------------------------------------------------*
 | checks for a keyword in actual line                                  |
 | ierr=0 not found                                                     |
 | ierr=1 found                                            m.gee 8/00   |
 | char string[] (input) character string to check actual line for      |
 *----------------------------------------------------------------------*/
void frchk(char string[], INT *ierr);
/*----------------------------------------------------------------------*
 | close and delete input file copy                        m.gee 4/01   |
 *----------------------------------------------------------------------*/
void frend(void);
/*----------------------------------------------------------------------*
 |  pss_pss.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void pss_write(char          *name, 
               INT            fdim, 
               INT            sdim,
               INT            byte,
               const void    *startaddress,
               long int      *handle,
               FILE          *out, 
               INT           *ierr);
void pss_write_array(const ARRAY *array, 
                     long int    *handle,
                     FILE        *out, 
                     INT         *ierr);
void pss_read_name(char      *name, 
                   INT       *fdim, 
                   INT       *sdim,
                   INT       *byte,
                   void      *ziel,
                   long int  *handle, 
                   FILE      *in,
                   INT       *ierr);
void pss_read_name_handle(char       *name, 
                          INT	     *fdim, 
                          INT	     *sdim,
                          INT	     *byte,
                          void       *ziel, 
                          long int   *handle,
                          FILE       *in, 
                          INT	     *ierr);
void pss_read_array_name(char       *name, 
                         ARRAY      *array,
                         long int   *handle,
                         FILE       *in,
                         INT        *ierr);
void pss_read_array_name_handle(char       *name, 
                                ARRAY	   *array,
                                long int   *handle,
                                FILE       *in,
                                INT	   *ierr);
void pss_read_array_handle(ARRAY      *array,
                           long int   *handle,
                           FILE       *in,
                           INT        *ierr);
void pss_chck(char       *name,
              long int   *handle, 
              FILE       *in,
              INT        *ierr);
void pss_chck_handle(char       *name,
                     long int   *handle, 
                     FILE       *in,
                     INT        *ierr);
void pss_getdims_name(char       *name, 
                      INT	 *fdim,
                      INT	 *sdim,
                      INT	 *byte,
                      long int	 *handle,
                      FILE       *in,
                      INT	 *ierr);
void pss_getdims_name_handle(char       *name, 
                             INT	*fdim,
                             INT	*sdim,
                             INT	*byte,
                             long int	*handle,
                             FILE       *in,
                             INT	*ierr);
void pss_status_to_err(FILE *inout);
/*----------------------------------------------------------------------*
 | pss_visual.c                                            genk 12/02   |
 *----------------------------------------------------------------------*/
void visual_writepss(FIELD  *actfield,  
                       INT     ntsteps,  
		       ARRAY  *time_a	 
		    );
void visual_readpss(FIELD   *actfield, 
                      INT     *ntsteps,  
		      ARRAY   *time_a	 
		     );
		     
/*----------------------------------------------------------------------*
 | ccarat_visual2.c                                        genk 07/02   |
 *----------------------------------------------------------------------*/
void vis2caf(INT numf, INT numa, INT nums);
void v2movie(void);
void v2cell(FIELD *actfield);
void v2_init(
             char *titl, INT *iopt, INT *cmncol, char *cmfile, INT *cmunit,
	     INT *xypix, float *xymin, float *xymax,
	     INT *nkeys, INT *ikeys, INT *tkeys, INT *fkeys, float **flims,
	     INT *mnode, INT *mptri, INT *mpptri,
	     INT *mface, INT *mpface, INT *medge, INT *mpedge
	    );
/*----------------------------------------------------------------------*
 | ccarat_visual3.c                                        genk 01/04   |
 *----------------------------------------------------------------------*/
void vis3caf(INT numff, INT numaf, INT numsf);

/*----------------------------------------------------------------------*
 | ccarat_visual3.c                                        genk 01/04   |
 *----------------------------------------------------------------------*/
void vis3caf(INT numff, INT numaf, INT numsf);

/*----------------------------------------------------------------------*
 | visual_readflaviares.c                                  genk 01/04   |
 *----------------------------------------------------------------------*/
void visual_readflaviares(FIELD   *actfield, 
                          INT     *ntsteps,  
		          ARRAY   *time_a,
			  ARRAY   *step_a,
			  INT     *FIRSTSTEP,
			  INT     *LASTSTEP,
			  INT     *DSTEP	 
		         )  ;
/*!---------------------------------------------------------------------                                         
\brief input of optimization data 

<pre>                                                          al  05/01      
</pre>  		 
\return void                                                                       

------------------------------------------------------------------------*/
void inpctropt(void);
/*!---------------------------------------------------------------------                                         
\brief control execution of optimization

<pre>                                                          al  05/01      
</pre>  		 
\return void                                                                       

------------------------------------------------------------------------*/
void caloptmain(void);


/* ====================================================================
 * file: ps_perf.c
 * ==================================================================== */
/*!---------------------------------------------------------------------
  \brief routine to meassure the performance                                              

  <pre>                                                        mn 01/04 
  Gets the current cpu tics.
  </pre>
  \return void                                                

  ------------------------------------------------------------------------*/
DOUBLE perf_cpu ();
/*!---------------------------------------------------------------------
  \brief routine to meassure the performance                                              

  <pre>                                                        mn 01/04 
  Gets the current system time.
  </pre>
  \return void                                                

  ------------------------------------------------------------------------*/
DOUBLE perf_time ();
/*!---------------------------------------------------------------------
  \brief routine to meassure the performance                                              

  <pre>                                                        mn 01/04 
  Initializes all counters with 0.
  </pre>
  \return void                                                

  ------------------------------------------------------------------------*/
void perf_init_all ();
/*!---------------------------------------------------------------------
  \brief routine to meassure the performance

  <pre>                                                        mn 01/04 
  Initializes one counters with 0.
  </pre>
  \return void                                                

  ------------------------------------------------------------------------*/
void perf_init (INT index);
/*!---------------------------------------------------------------------
  \brief routine to meassure the performance

  <pre>                                                        mn 01/04 
  start of the region for one timer.
  </pre>
  \return void                                                

  ------------------------------------------------------------------------*/
void perf_begin (INT index);
/*!---------------------------------------------------------------------
  \brief routine to meassure the performance

  <pre>                                                        mn 01/04 
  end of the region for one timer.
  </pre>
  \return void                                                


  ------------------------------------------------------------------------*/
void perf_end (INT index);
/*!---------------------------------------------------------------------
  \brief routine to meassure the performance

  <pre>                                                        mn 01/04 
  Print the results for one timer.
  </pre>
  \return void                                                

  ------------------------------------------------------------------------*/
void perf_print (INT index, char string[], INT bezug, INT ops);
/*!---------------------------------------------------------------------
  \brief routine to meassure the performance

  <pre>                                                        mn 01/04 
  Print the results for all timers.
  </pre>
  \return void                                                

  ------------------------------------------------------------------------*/
void perf_out ();


#ifdef CHECK_MAX
/* ====================================================================
 * file: check_max_sizes.c
 * ==================================================================== */
/*!----------------------------------------------------------------------
\brief check the values of the max sizes

<pre>                                                              mn 04/04
This routine determines the optimal values for maxele, maxnod, maxdofpernode
and maxgauss and compares those to the given values of the respective defines.
</pre>

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: --- 

*----------------------------------------------------------------------*/
void check_max_sizes(
    );
#endif

