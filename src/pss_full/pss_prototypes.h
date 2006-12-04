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

#ifndef PSS_PROTOTYPES_H
#define PSS_PROTOTYPES_H


/*----------------------------------------------------------------------*
 |  pss_am.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
void ShiftPointer(
    void    **ptr,
    PTRSIZE   diff);

/*----------------------------------------------------------------------*
 | redefinition of malloc DEBUG version                       m.gee 2/02|
 | bhaves exactly like malloc conform to ansi c standard                |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | redefinition of malloc FAST version                        m.gee 2/02|
 | bhaves exactly like malloc conform to ansi c standard                |
 *----------------------------------------------------------------------*/
void *CCAMALLOC(
    INT   size);

/*----------------------------------------------------------------------*
 | redefinition of calloc DEBUG version                   m.gee 2/02    |
 | bhaves exactly like calloc conform to ansi c standard                |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | redefinition of calloc FAST version                    m.gee 2/02    |
 | bhaves exactly like calloc conform to ansi c standard                |
 *----------------------------------------------------------------------*/
void *CCACALLOC(
    INT    num,
    INT    size);

/*----------------------------------------------------------------------*
 | redefinition of realloc DEBUG version                  m.gee 2/02    |
 | bhaves exactly like realloc conform to ansi c standard               |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | redefinition of realloc FAST version                  m.gee 2/02     |
 | bhaves exactly like realloc conform to ansi c standard               |
 *----------------------------------------------------------------------*/
void *CCAREALLOC(
    void    *oldptr,
    INT      size);

/*----------------------------------------------------------------------*
 | redefinition of free DEBUG version                     m.gee 2/02    |
 | bhaves exactly like free conform to ansi c standard                  |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | redefinition of free FAST version                      m.gee 2/02    |
 | bhaves exactly like free conform to ansi c standard                  |
 *----------------------------------------------------------------------*/
void *CCAFREE(
    void     *oldptr);

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
void* amdef(
    char     *namstr,
    ARRAY    *a,
    INT       fdim,
    INT       sdim,
    char      typstr[]);

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
void* amredef(
    ARRAY     *a,
    INT        newfdim,
    INT        newsdim,
    char       newtypstr[]);

/*----------------------------------------------------------------------*
 | delete         array                                   m.gee 8/00    |
 | frees the vector or array located in the ARRAY                       |
 | *array (input) the adress of the structure holding the vector/array  |
 *----------------------------------------------------------------------*/
void  amdel(
    ARRAY     *array);

/*----------------------------------------------------------------------*
 | initialize an array by zero                            m.gee 8/00    |
 | initializes the content of the ARRAY array to zero                   |
 | put 0 to integer fields, 0.0 to DOUBLE fields                        |
 | ARRAY *array (input) adress of the ARRAY array                       |
 *----------------------------------------------------------------------*/
void  amzero(
    ARRAY     *array);

/*----------------------------------------------------------------------*
 | multiply an array by a scalar                          m.gee 6/01    |
 | scales the contents of a field by a given value                      |
 | ARRAY *array (input) adress of the ARRAY array                       |
 | *value (input) adress of the scaling parameter, this may be of type  |
 |                INT* or DOUBLE* and must be casted to void* in the    |
 |                parameter list                                        |
 |                example: amscal(&val,(void*)(&ione));                 |
 *----------------------------------------------------------------------*/
void  amscal(
    ARRAY     *array,
    void      *value);

/*----------------------------------------------------------------------*
 | initialize an array by value                           m.gee 6/01    |
 | inits the contents of a field by a given value                       |
 | ARRAY *array (input) adress of the ARRAY array                       |
 | *value (input) adress of the initvalues, this may be of type         |
 |                INT* or DOUBLE* and must be casted to void* in the    |
 |                parameter list                                        |
 |                example: aminit(&val,(void*)(&ione));                 |
 *----------------------------------------------------------------------*/
void  aminit(
    ARRAY     *array,
    void      *value);

/*----------------------------------------------------------------------*
 | allocate and make a copy of ARRAY                      m.gee 8/00    |
 | alocates a new field in array_to of the same size and type as        |
 | in array_from. Then copies the contents from array_from to array_to  |
 | ARRAY *array_from (input) the adress of the field one wants to copy  |
 | ARRAY *array_to   (in/out) the adress of the structure where the new |
 |                            field is allocated in                     |
 | return value: adress of the new field                                |
 *----------------------------------------------------------------------*/
void* am_alloc_copy(
    ARRAY     *array_from,
    ARRAY     *array_to);

/*----------------------------------------------------------------------*
 | make a copy of ARRAY,                                  m.gee 6/01    |
 | Copies the contents from array_from to array_to                      |
 | ARRAY *array_from (input) the adress of the field one wants to copy  |
 | ARRAY *array_to   (in/out) the adress of the field values are copied |
 |                            to                                        |
 | user must provide fields of matching type and size!                  |
 | return value: array_to->a.iv                                         |
 *----------------------------------------------------------------------*/
void* amcopy(
    ARRAY     *array_from,
    ARRAY     *array_to);

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
void  amadd(
    ARRAY      *array_to,
    ARRAY      *array_from,
    DOUBLE      factor,
    INT         init);

void amprint(FILE* err,ARRAY *a,INT fdim, INT sdim);

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
void* am4def(
    char     *namstr,
    ARRAY4D  *a,
    INT       fdim,
    INT       sdim,
    INT       tdim,
    INT       fodim,
    char      typstr[]);

/*----------------------------------------------------------------------*
 | delete 4dimensional array                             m.gee 12/01    |
 | frees all field memory in array                                      |
 | ARRAY4D *array    (input)                adress of array-structure   |
 | no return value                                                      |
 *----------------------------------------------------------------------*/
void  am4del(
    ARRAY4D    *array);

/*----------------------------------------------------------------------*
 | initialize an 4D array by zero                           m.gee 12/01 |
 | see head of amzero                                                   |
 *----------------------------------------------------------------------*/
void  am4zero(
    ARRAY4D    *array);

/*----------------------------------------------------------------------*
 | initialize a 4D array by value                         m.gee 6/01    |
 | see head of aminit                                                   |
 *----------------------------------------------------------------------*/
void  am4init(
    ARRAY4D     *array,
    void        *value);

/*----------------------------------------------------------------------*
 | allocate and make a copy of ARRAY4D                     m.gee 12/01  |
 | see head of am_alloc_copy                                            |
 *----------------------------------------------------------------------*/
void* am4_alloc_copy(
    ARRAY4D      *array_from,
    ARRAY4D      *array_to);

/*----------------------------------------------------------------------*
 | make a copy of ARRAY4D                                  m.gee 12/01  |
 | the given arrays must be of equal type and size                      |
 | see head of amcopy                                                   |
 *----------------------------------------------------------------------*/
void* am4copy(
    ARRAY4D     *array_from,
    ARRAY4D     *array_to);

/*----------------------------------------------------------------------*
 | redefine of ARRAY4D                                     m.gee 12/01  |
 | NOTE:  - a dimension cast like it is possible with ARRAY doesn't work|
 |        - 4-dimensional arrays with fodim=0 are not allowed           |
 |        - remember, this routine can be VERY expensive !              |
 |        - aditional space in the redefined array is set to zero       |
 |          (unlike in am4def which does NOT initialize)                |
 |        usage similar to amredef                                      |
 *----------------------------------------------------------------------*/
void* am4redef(
    ARRAY4D     *array,
    INT          newfdim,
    INT          newsdim,
    INT          newtdim,
    INT          newfodim);


void amprint(FILE* err,ARRAY *a,INT fdim, INT sdim);


/*----------------------------------------------------------------------*
 |  pss_ds.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | report the amount of actual allocated memory           m.gee 2/02    |
 | only works if DEBUG is defined                                       |
 *----------------------------------------------------------------------*/
void dsmemreport( void);

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
void dsassert_func(
    char*               file,
    INT                 line,
    INT                 test,
    char                string[]
    );

/* to make sure we do not have a redefinition here as its also defined in dserror.H */
#ifndef DSERROR_H
#define DSERROR_H

#define dsassert(test, string) dsassert_func(__FILE__, __LINE__, test, string)

/*!
  \brief always call dslatest before dserror so that we have the right
  position.

  The point here is that the macro does not accept arguments. This way
  the dserror name is replaced by a dslatest call and the arguments
  the user supplied to dserror apply to the function dserror_func.
  However, the two function calls are not separated by a semicolon (;)
  but by a comma (,) and thus the whole line is still one
  statement. In particular this dserror macro can be used in if
  clauses without braces ({}).
 */
#define dserror dslatest(__FILE__, __LINE__), dserror_func

#endif 

/*----------------------------------------------------------------------*
 | report an error and stop program                       m.gee 8/00    |
 | prints error message string to console and *.err                     |
 | prints call tree, if DEBUG was defined                               |
 | aborts parallel and sequentiell programm                             |
 *----------------------------------------------------------------------*/
void dserror_func(char *string, ...);


/*!
  \brief set a file name and line number to be used by the next
  dserror_func
 */
void dslatest(char* file, INT line);


/*----------------------------------------------------------------------*
 | collects warnings during running process                  ck 07/03   |
 | and writes them to the screen at the end                             |
 *----------------------------------------------------------------------*/
void dswarning(
    INT     task,
    INT     warning);

/*----------------------------------------------------------------------*
 | Initialize bugtracing systems                          m.gee 8/00    |
 *----------------------------------------------------------------------*/
void dsinit(void);

/*----------------------------------------------------------------------*
 | report intrance to routine named string to the tracing system  m.gee |
 | is empty is DEBUG is not defined                               8/00  |
 *----------------------------------------------------------------------*/
void dstrc_enter(
    char      string[]);

/*----------------------------------------------------------------------*
 | report exit to routine to the tracing system           m.gee 8/00    |
 | is empty is DEBUG is not defined                                     |
 *----------------------------------------------------------------------*/
void dstrc_exit(void);

void dstrc_whereami(void);

/*----------------------------------------------------------------------*
 | report a new DOUBLE array to the bugtracing system     m.gee 8/00    |
 | this routine is called by the am-system only !                       |
 *----------------------------------------------------------------------*/
void dsreportarray(
    void     *array,
    INT       typ);

/*----------------------------------------------------------------------*
 | report a new DOUBLE array to the bugtracing system     m.gee 8/00    |
 | this routine is called by the am-system only !                       |
 *----------------------------------------------------------------------*/
void dsdeletearray(
    void      *array,
    INT        typ);

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
INT frfind(
    char      string[]);

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
void frint_n(
    char       string[],
    INT       *var,
    INT        num,
    INT       *ierr);

/*----------------------------------------------------------------------*
 | reads an integer from input_file                       m.gee 8/00    |
 | starting in allfiles.actrow                                          |
 | char string[] (input) keyword to search for in actual line           |
 | INT *var      (output)adress of variable to read to                  |
 | INT *ierr     (output) flag to indicate success                      |
 | ierr=0 keyword not found                                             |
 | ierr=1 integer read                                                  |
 *----------------------------------------------------------------------*/
void frint(
    char      string[],
    INT      *var,
    INT      *ierr);

/*----------------------------------------------------------------------*
 | reads a DOUBLE from input file line allfiles.line      m.gee 8/00    |
 | see frint_n                                                          |
 | ierr=0 kenner not found                                              |
 | ierr=1 DOUBLE  read                                                  |
 *----------------------------------------------------------------------*/
void frdouble_n(
    char      string[],
    DOUBLE   *var,
    INT       num,
    INT      *ierr);

/*----------------------------------------------------------------------*
 | reads a DOUBLE from input file line allfiles.line      m.gee 8/00    |
 | see frint                                                            |
 | ierr=0 kenner not found                                              |
 | ierr=1 DOUBLE  read                                                  |
 *----------------------------------------------------------------------*/
void frdouble(
    char       string[],
    DOUBLE    *var,
    INT       *ierr);

/*----------------------------------------------------------------------*
 | reads a charstring from input file                                   |
 | user must assure, that the given charpointer space is long enough    |
 | to hold the string                                                   |
 | ierr=0 keyword not found on line                                     |
 | ierr=1 char string read                                 m.gee 8/00   |
 *----------------------------------------------------------------------*/
void frchar(
    char       string[],
    char      *var,
    INT       *ierr);

/*----------------------------------------------------------------------*
 | checks for a keyword in actual line                                  |
 | ierr=0 not found                                                     |
 | ierr=1 found                                            m.gee 8/00   |
 | char string[] (input) character string to check actual line for      |
 *----------------------------------------------------------------------*/
void frchk(
    char      string[],
    INT      *ierr);

/*----------------------------------------------------------------------*
 | close and delete input file copy                        m.gee 4/01   |
 *----------------------------------------------------------------------*/
void frend(void);

void frword(char string[],char *var, INT *ierr);

/* Compare two words */
INT frwordcmp(CHAR* p1, CHAR* p2);

INT frcheckyes(CHAR* p);
INT frcheckno(CHAR* p);
INT frreadyes(CHAR* key, INT* flag);




/*----------------------------------------------------------------------*
 |  pss_pss.c                                            m.gee 11/01    |
 *----------------------------------------------------------------------*/
void pss_write(
    char          *name,
    INT            fdim,
    INT            sdim,
    INT            byte,
    const void    *startaddress,
    long int      *handle,
    FILE          *out,
    INT           *ierr);

void pss_write_array(
    const ARRAY *array,
    long int    *handle,
    FILE        *out,
    INT         *ierr);

void pss_read_name(
    char      *name,
    INT       *fdim,
    INT       *sdim,
    INT       *byte,
    void      *ziel,
    long int  *handle,
    FILE      *in,
    INT       *ierr);

void pss_read_name_handle(
    char       *name,
    INT	     *fdim,
    INT	     *sdim,
    INT	     *byte,
    void       *ziel,
    long int   *handle,
    FILE       *in,
    INT	     *ierr);

void pss_read_array_name(
    char       *name,
    ARRAY      *array,
    long int   *handle,
    FILE       *in,
    INT        *ierr);

void pss_read_array_name_handle(
    char       *name,
    ARRAY	   *array,
    long int   *handle,
    FILE       *in,
    INT	   *ierr);

void pss_read_array_handle(
    ARRAY      *array,
    long int   *handle,
    FILE       *in,
    INT        *ierr);

void pss_chck(
    char       *name,
    long int   *handle,
    FILE       *in,
    INT        *ierr);

void pss_chck_handle(
    char       *name,
    long int   *handle,
    FILE       *in,
    INT        *ierr);

void pss_getdims_name(
    char       *name,
    INT	 *fdim,
    INT	 *sdim,
    INT	 *byte,
    long int	 *handle,
    FILE       *in,
    INT	 *ierr);

void pss_getdims_name_handle(
    char       *name,
    INT	*fdim,
    INT	*sdim,
    INT	*byte,
    long int	*handle,
    FILE       *in,
    INT	*ierr);

void pss_status_to_err(
    FILE *inout);



/*----------------------------------------------------------------------*
 | pss_visual.c                                            genk 12/02   |
 *----------------------------------------------------------------------*/
void visual_writepss(
    FIELD  *actfield,
    INT     ntsteps,
    ARRAY  *time_a);
void visual_readpss(
    FIELD   *actfield,
    INT     *ntsteps,
    ARRAY   *time_a);



/* ====================================================================
 * file: ps_perf.c
 * ==================================================================== */
#ifdef PERF
/*!---------------------------------------------------------------------
  \brief routine to meassure the performance

  <pre>                                                        mn 01/04
  Initializes all counters with 0.
  </pre>
  \return void
  ------------------------------------------------------------------------*/
void perf_init_all (void);


/*!---------------------------------------------------------------------
  \brief routine to meassure the performance

  <pre>                                                        mn 01/04
  Initializes one counters with 0.
  </pre>
  \return void
  ------------------------------------------------------------------------*/
void perf_init (
    INT      index);


/*!---------------------------------------------------------------------
  \brief routine to meassure the performance

  <pre>                                                        mn 01/04
  start of the region for one timer.
  </pre>
  \return void
  ------------------------------------------------------------------------*/
void perf_begin (
    INT      index);


/*!---------------------------------------------------------------------
  \brief routine to meassure the performance

  <pre>                                                        mn 01/04
  end of the region for one timer.
  </pre>
  \return void
  ------------------------------------------------------------------------*/
void perf_end (
    INT       index);



/*!---------------------------------------------------------------------
  \brief routine to meassure the performance

  <pre>                                                        mn 01/04
  Print the results for all timers.
  </pre>
  \return void
  ------------------------------------------------------------------------*/
void perf_out (void);

#else
#define perf_init_all()
#define perf_init()
#define perf_begin(index)
#define perf_end(index)
#define perf_out()
#endif

/*!---------------------------------------------------------------------
  \brief routine to meassure the performance

  <pre>                                                        mn 01/04
  Gets the current system time.
  </pre>
  \return void
  ------------------------------------------------------------------------*/
DOUBLE perf_time (void);


#endif

