/*----------------------------------------------------------------------*
 |  main_ccarat.c                                        m.gee 11/01    |
 *----------------------------------------------------------------------*/
int main(int argc, char *argv[]);
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
                                int        control_node_global,
                                int        control_dof,
                                NODE     **node,
                                int       *cdof); 
/*----------------------------------------------------------------------*
 | global_control.c                                      m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntam(int argc, char *argv[]);
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
void ntaini(int argc, char *argv[]);
/*----------------------------------------------------------------------*
 | global_inp_control.c                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntainp(void);
/*----------------------------------------------------------------------*
 | global_mask_matrices.c                                m.gee 11/01    |
 *----------------------------------------------------------------------*/
void mask_global_matrices(void);
/*----------------------------------------------------------------------*
 |  machine_hpux.c                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ntadev(int argc, char *argv[]);
/*----------------------------------------------------------------------*
 |  restart_control.c                                    m.gee 02/02    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |  map_node_find.c                                      m.gee 11/01    |
 *----------------------------------------------------------------------*/
void iscouple_find_node_comp(NODE  *actnode, 
                                FIELD *searchfield, 
                                NODE **partnernode,
                                int    coupleID,
                                int    dof);
void cheque_distance(double *x1, double *x2, double tol, int *ierr);
void find_assign_coupset(FIELD *actfield, 
                            int    coupleID, 
                            int   *counter);
/*----------------------------------------------------------------------*
 |  dyn_timecurve.c                                      m.gee 02/02    |
 *----------------------------------------------------------------------*/
void dyn_init_curve(int actcurve,
                   int    nstep,
                   double dt,
                   double maxtime);
void dyn_facfromcurve(int actcurve,
                   double T,
                   double *fac);
double dyn_facexplcurve(int actcurve,
                      double T);		   
/*----------------------------------------------------------------------*
 |  out_global.c                                         m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_general(void);
void out_sol(FIELD *actfield, PARTITION *actpart, INTRA *actintra, 
             int step, int place);
/*----------------------------------------------------------------------*
 |  out_gid_sol.c                                        m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_sol_init(void);
void out_gid_domains(FIELD *actfield);
void out_gid_sol(char string[], FIELD *actfield, INTRA  *actintra, int step,
                 int place);
/*----------------------------------------------------------------------*
 |  out_gid_msh.c                                        m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_msh(void);
void out_gid_allcoords(FILE *out);
/*----------------------------------------------------------------------*
 |  inherit_insidedesign.c                                  m.gee 3/02  |
 *----------------------------------------------------------------------*/
void inherit_dirich_coup_indesign(void);
/*----------------------------------------------------------------------*
 |  inherit_design_dis.c                                    m.gee 3/02  |
 *----------------------------------------------------------------------*/
void inherit_design_dis_dirichlet(DISCRET *actdis);
void inherit_design_dis_couple(DISCRET *actdis);
void inherit_design_dis_neum(DISCRET *actdis);
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
void inpctr_dyn_struct(STRUCT_DYNAMIC *sdyn);
/*!---------------------------------------------------------------------                                         
\brief input of the FLUID DYNAMIC block in the input-file

<pre>                                                         genk 03/02

In this routine the data in the FLUID DYNAMIC block of the input file
are read and stored in fdyn	       

</pre>
\param  *data 	  FLUID_DATA       (o)	   
\return void                                                                       

------------------------------------------------------------------------*/
void inpctr_dyn_fluid(FLUID_DYNAMIC *fdyn);
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
 |  input_design.c                                  m.gee 11/01         |
 *----------------------------------------------------------------------*/
void inpdesign(void);
void inp_dnode(void);
void read_1_dnode(DNODE *dnode, int readId);
void inp_dline(void);
void read_1_dline(DLINE *dline, int readId);
void inp_dsurface(void);
void read_1_dsurf(DSURF *dsurf, int readId);
void inp_dvolume(void);
void read_1_dvol(DVOL *dvol, int readId);
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
 |  input_topology.c                                  m.gee 11/01       |
 *----------------------------------------------------------------------*/
void inp_topology(DISCRET *actdis);
void inp_detailed_topology(DISCRET   *actdis);
/*----------------------------------------------------------------------*
 |  math1.c                                               m.gee 11/01   |
 *----------------------------------------------------------------------*/
void math_array_copy(double **from, int n, int m, double **to);
void math_inv3(double **a, double *det);
void math_tran(double **a, int n);
void math_unvc(double *enorm,double *vec, int n);
void math_matvecdense(double  *r,
                         double **A,
                         double  *b,
                         int      ni,
                         int      nk,
                         int      init,
                         double   factor);
void math_mattrnvecdense(double  *r,
                         double **A,
                         double  *b,
                         int      ni,
                         int      nk,
                         int      init,
                         double   factor);
void math_matmatdense(double **R,
                         double **A,
                         double **B,
                         int      ni,
                         int      nk,
                         int      nj,
                         int      init,
                         double   factor);
void math_mattrnmatdense(double **R,
                            double **A,
                            double **B,
                            int      ni,
                            int      nk,
                            int      nj,
                            int      init,
                            double   factor);
void math_matmattrndense(double **R,
                            double **A,
                            double **B,
                            int      ni,
                            int      nk,
                            int      nj,
                            int      init,
                            double   factor);
void math_sym_inv(double **A, int dim);
void math_unsym_inv(double **A, int dimr, int dimc);
void math_sppr(double *spat, double *a, double *b, double *c);
void math_addab(double **a, double **b, int dim1, int dim2);
/*!---------------------------------------------------------------------                                         
\brief extract digits from integer number

<pre>                                                         genk 04/02		     
</pre>   
\param  num	 int   (i)    integer number
\param *it	 int   (o)    integer on position "thousand"
\param *ih       int   (o)    integer on position "hundred"
\param *id       int   (o)    integer on position "ten"
\param *id       int   (o)    integer on position "one"
\return void 

------------------------------------------------------------------------*/
void math_intextract(
                    int num,    
                    int *it,    
		    int *ih,    
		    int *id,    
		    int *io     
	            );

/*----------------------------------------------------------------------*
 |  sort_find.c                                          m.gee 11/01    |
 *----------------------------------------------------------------------*/
int cmp_int(const void *a, const void *b );
double cmp_double(const void *a, const void *b );
void mg_sort(int list[], int N, int list2[], double list3[]);
int quick_find(int key, int list[], int length, int shift, int bins[]);
void init_quick_find(int list[], int length, int *shift, int *bins);
int find_index(int key, int list[], int length);
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
 |#else                                a 32 bit pointer is of size int  |
 |typedef int PTRSIZE;                                                  |
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
void *CCAMALLOC(int size);
/*----------------------------------------------------------------------*
 | redefinition of calloc DEBUG version                   m.gee 2/02    |
 | bhaves exactly like calloc conform to ansi c standard                |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | redefinition of calloc FAST version                    m.gee 2/02    |
 | bhaves exactly like calloc conform to ansi c standard                |
 *----------------------------------------------------------------------*/
void *CCACALLOC(int num, int size);
/*----------------------------------------------------------------------*
 | redefinition of realloc DEBUG version                  m.gee 2/02    |
 | bhaves exactly like realloc conform to ansi c standard               |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | redefinition of realloc FAST version                  m.gee 2/02     |
 | bhaves exactly like realloc conform to ansi c standard               |
 *----------------------------------------------------------------------*/
void *CCAREALLOC(void *oldptr, int size);
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
 | int fdim      (input)   first dimension of 2D vector                 |
 |                         dimension of 1D vector                       |
 | int sdim      (input)   scnd dimension of 2D vector                  |
 | char typstr[] (input)   type of array to allocate                    |
 |               ="IV"     allocate integer vector in a->a.iv           |
 |               ="IA"     allocate double array  in a->a.ia            |
 |               ="DV"     allocate integer vector in a->a.dv           |
 |               ="DA"     allocate double array  in a->a.da            |
 | return value:                                                        |
 | void pointer to allocated memory                                     |
 *----------------------------------------------------------------------*/
void* amdef(char *namstr,ARRAY *a,int fdim, int sdim, char typstr[]);
/*----------------------------------------------------------------------*
 | redefine array                                         m.gee 8/00    |
 | changes an already allocated array a in dimensions                   |
 | a typecast of the values in the array is not possible                |
 | (no int to double or vice versa transformation)                      |
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
 | int newfdim   (input)   new first dimension of 2D vector             |
 |                         new dimension of 1D vector                   |
 | int newsdim   (input)   new scnd dimension of 2D vector              |
 | char newtypstr[] (input)   type of array to allocate                 |
 |               ="IV"     allocate integer vector in a->a.iv           |
 |               ="IA"     allocate double array  in a->a.ia            |
 |               ="DV"     allocate integer vector in a->a.dv           |
 |               ="DA"     allocate double array  in a->a.da            |
 | return value:                                                        |
 | void pointer to allocated memory                                     |
 *----------------------------------------------------------------------*/
void* amredef(ARRAY *a,int newfdim, int newsdim, char newtypstr[]);
/*----------------------------------------------------------------------*
 | delete         array                                   m.gee 8/00    |
 | frees the vector or array located in the ARRAY                       |
 | *array (input) the adress of the structure holding the vector/array  |
 *----------------------------------------------------------------------*/
void  amdel(ARRAY *array);
/*----------------------------------------------------------------------*
 | initialize an array by zero                            m.gee 8/00    |
 | initializes the content of the ARRAY array to zero                   |
 | put 0 to integer fields, 0.0 to double fields                        |
 | ARRAY *array (input) adress of the ARRAY array                       |
 *----------------------------------------------------------------------*/
void  amzero(ARRAY *array);
/*----------------------------------------------------------------------*
 | multiply an array by a scalar                          m.gee 6/01    |
 | scales the contents of a field by a given value                      |
 | ARRAY *array (input) adress of the ARRAY array                       |
 | *value (input) adress of the scaling parameter, this may be of type  |
 |                int* or double* and must be casted to void* in the    |
 |                parameter list                                        |
 |                example: amscal(&val,(void*)(&ione));                 |        
 *----------------------------------------------------------------------*/
void  amscal(ARRAY *array, void *value);
/*----------------------------------------------------------------------*
 | initialize an array by value                           m.gee 6/01    |
 | inits the contents of a field by a given value                       |
 | ARRAY *array (input) adress of the ARRAY array                       |
 | *value (input) adress of the initvalues, this may be of type         |
 |                int* or double* and must be casted to void* in the    |
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
void  amadd(ARRAY *array_to, ARRAY *array_from, double factor, int init);
/*----------------------------------------------------------------------*
 | define 4D array                                       m.gee 12/01    |
 | similar to amdef, but for 3D and 4D fields                           |
 | the field is NOT initialized                                         |
 |                                                                      |
 | char *namstr  (input)   name of array                                |
 | ARRAY4D *a    (input)   adress of structure ARRAY4D                  |
 | int fdim      (input)   first dimension of 3D or 4D array            |
 | int sdim      (input)   scnd dimension of 3D or 4D array             |
 | int tdim      (input)   third dimension of 3D or 4D array            |
 | int fodim     (input)   fourth dimension of 4D array,                |
 |                         ==0 of array is 3D                           |
 | char typstr[] (input)   type of field to allocate                    |
 |                         ="I3" 3D integer field                       |
 |                         ="I4" 4D integer field                       |
 |                         ="D3" 3D double  field                       |
 |                         ="D4" 4D double  field                       |
 | return value:                                                        |
 | void pointer to allocated memory                                     |
 *----------------------------------------------------------------------*/
void* am4def(char *namstr, ARRAY4D *a, int fdim, int sdim, int tdim, 
             int fodim, char typstr[]);
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
void* am4redef(ARRAY4D *array, int newfdim, int newsdim, int newtdim, 
               int newfodim);
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
 | true - (input) boolean criterium                                     |
 | string - (input) error message                                       |
 *----------------------------------------------------------------------*/
void dsassert(int true, char string[]);
/*----------------------------------------------------------------------*
 | report an error and stop program                       m.gee 8/00    |
 | prints error message string to console and *.err                     |
 | prints call tree, if DEBUG was defined                               |
 | aborts parallel and sequentiell programm                             |
 *----------------------------------------------------------------------*/
void dserror(char string[]);
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
 | report a new double array to the bugtracing system     m.gee 8/00    |    
 | this routine is called by the am-system only !                       | 
 *----------------------------------------------------------------------*/
void dsreportarray(void *array, int typ);
/*----------------------------------------------------------------------*
 | report a new double array to the bugtracing system     m.gee 8/00    |    
 | this routine is called by the am-system only !                       |
 *----------------------------------------------------------------------*/
void dsdeletearray(void *array, int typ);
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
double ds_cputime(void);
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
void frfind(char string[]);
/*----------------------------------------------------------------------*
 |                                                        m.gee 8/00    |
 | sets a pointer to the next line in thecopy of input file on all procs|
 *----------------------------------------------------------------------*/
void frread(void);
/*----------------------------------------------------------------------*
 | reads n integers from input_file                       m.gee 4/01    |
 | char string[] (input) keyword to search for in actual line           |
 | int *var      (output) adress of field to hold values read           |
 | int num       (input)  number of values to read                      |
 | int *ierr     (output) =0 keyword not found / =1 values read         |
 *----------------------------------------------------------------------*/
void frint_n(char string[],int *var,int num, int *ierr);
/*----------------------------------------------------------------------*
 | reads an integer from input_file                       m.gee 8/00    |
 | starting in allfiles.actrow                                          |
 | char string[] (input) keyword to search for in actual line           |
 | int *var      (output)adress of variable to read to                  |
 | int *ierr     (output) flag to indicate success                      |
 | ierr=0 keyword not found                                             |
 | ierr=1 integer read                                                  |
 *----------------------------------------------------------------------*/
void frint(char string[],int *var, int *ierr);
/*----------------------------------------------------------------------*
 | reads a double from input file line allfiles.line      m.gee 8/00    |
 | see frint_n                                                          |
 | ierr=0 kenner not found                                              |
 | ierr=1 double  read                                                  |
 *----------------------------------------------------------------------*/
void frdouble_n(char string[],double *var,int num, int *ierr);
/*----------------------------------------------------------------------*
 | reads a double from input file line allfiles.line      m.gee 8/00    |
 | see frint                                                            |
 | ierr=0 kenner not found                                              |
 | ierr=1 double  read                                                  |
 *----------------------------------------------------------------------*/
void frdouble(char string[],double *var, int *ierr);
/*----------------------------------------------------------------------*
 | reads a charstring from input file                                   |
 | user must assure, that the given charpointer space is long enough    |   
 | to hold the string                                                   |
 | ierr=0 keyword not found on line                                     |
 | ierr=1 char string read                                 m.gee 8/00   |
 *----------------------------------------------------------------------*/
void frchar(char string[],char *var, int *ierr);
/*----------------------------------------------------------------------*
 | checks for a keyword in actual line                                  |
 | ierr=0 not found                                                     |
 | ierr=1 found                                            m.gee 8/00   |
 | char string[] (input) character string to check actual line for      |
 *----------------------------------------------------------------------*/
void frchk(char string[], int *ierr);
/*----------------------------------------------------------------------*
 | close and delete input file copy                        m.gee 4/01   |
 *----------------------------------------------------------------------*/
void frend(void);
/*----------------------------------------------------------------------*
 |  pss_pss.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void pss_write(char          *name, 
               int            fdim, 
               int            sdim,
               int            byte,
               const void    *startaddress,
               long int      *handle,
               FILE          *out, 
               int           *ierr);
void pss_write_array(const ARRAY *array, 
                     long int    *handle,
                     FILE        *out, 
                     int         *ierr);
void pss_read_name(char      *name, 
                   int       *fdim, 
                   int       *sdim,
                   int       *byte,
                   void      *ziel,
                   long int  *handle, 
                   FILE      *in,
                   int       *ierr);
void pss_read_name_handle(char       *name, 
                          int	     *fdim, 
                          int	     *sdim,
                          int	     *byte,
                          void       *ziel, 
                          long int   *handle,
                          FILE       *in, 
                          int	     *ierr);
void pss_read_array_name(char       *name, 
                         ARRAY      *array,
                         long int   *handle,
                         FILE       *in,
                         int        *ierr);
void pss_read_array_name_handle(char       *name, 
                                ARRAY	   *array,
                                long int   *handle,
                                FILE       *in,
                                int	   *ierr);
void pss_read_array_handle(ARRAY      *array,
                           long int   *handle,
                           FILE       *in,
                           int        *ierr);
void pss_chck(char       *name,
              long int   *handle, 
              FILE       *in,
              int        *ierr);
void pss_chck_handle(char       *name,
                     long int   *handle, 
                     FILE       *in,
                     int        *ierr);
void pss_getdims_name(char       *name, 
                      int	 *fdim,
                      int	 *sdim,
                      int	 *byte,
                      long int	 *handle,
                      FILE       *in,
                      int	 *ierr);
void pss_getdims_name_handle(char       *name, 
                             int	*fdim,
                             int	*sdim,
                             int	*byte,
                             long int	*handle,
                             FILE       *in,
                             int	*ierr);
void pss_status_to_err(FILE *inout);
/*----------------------------------------------------------------------*
 | ccarat_visual2.c                                        genk 07/02   |
 *----------------------------------------------------------------------*/
/*!---------------------------------------------------------------------                                         
\brief call of visual2 for fluid 

<pre>                                                         genk 07/02      
</pre>  
\param  numf   int      (i)       actual number of fluid field
\return void                                                                       

------------------------------------------------------------------------*/
void vis2caf(int numf);
/*!---------------------------------------------------------------------                                         
\brief  compute the array WCELL for use in QAT2V2 

<pre>                                                         genk 07/02      
</pre>  
\param  *actfield     FIELD    (i)  actual field		 
\return void                                                                       
\warning QAT2V2 requires FORTRAN numbering, so increase node-numbers by one!!

------------------------------------------------------------------------*/
void v2cell(FIELD *actfield);
/*!---------------------------------------------------------------------                                         
\brief dummy routine 

<pre>                                                         genk 07/02      

since VISUAL2 is called by a fortran routine, this one is necessary, if
VIS2-routines are not compiled into the program

</pre>  		 
\return void                                                                       

------------------------------------------------------------------------*/
void v2_init(void);
