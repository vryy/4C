/*!---------------------------------------------------------------------
\brief Prototype classes

\maintainer Martin Kronbichler

\level 1

---------------------------------------------------------------------*/

#ifndef PSS_PROTOTYPES_H
#define PSS_PROTOTYPES_H


/*----------------------------------------------------------------------*
 |  pss_am.c                                             m.gee 11/01    |
 *----------------------------------------------------------------------*/
/* void ShiftPointer(
    void    **ptr,
    PTRSIZE   diff);
*/
/*----------------------------------------------------------------------*
 | redefinition of malloc DEBUG version                       m.gee 2/02|
 | bhaves exactly like malloc conform to ansi c standard                |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | redefinition of malloc FAST version                        m.gee 2/02|
 | bhaves exactly like malloc conform to ansi c standard                |
 *----------------------------------------------------------------------*/
void *CCAMALLOC(unsigned size);

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
void *amdef(char *namstr, ARRAY *a, INT fdim, INT sdim, char typstr[]);

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
void *amredef(ARRAY *a, INT newfdim, INT newsdim, char newtypstr[]);

/*----------------------------------------------------------------------*
 | delete         array                                   m.gee 8/00    |
 | frees the vector or array located in the ARRAY                       |
 | *array (input) the adress of the structure holding the vector/array  |
 *----------------------------------------------------------------------*/
void amdel(ARRAY *array);

/*----------------------------------------------------------------------*
 | initialize an array by zero                            m.gee 8/00    |
 | initializes the content of the ARRAY array to zero                   |
 | put 0 to integer fields, 0.0 to DOUBLE fields                        |
 | ARRAY *array (input) adress of the ARRAY array                       |
 *----------------------------------------------------------------------*/
void amzero(ARRAY *array);

/*----------------------------------------------------------------------*
 | multiply an array by a scalar                          m.gee 6/01    |
 | scales the contents of a field by a given value                      |
 | ARRAY *array (input) adress of the ARRAY array                       |
 | *value (input) adress of the scaling parameter, this may be of type  |
 |                INT* or DOUBLE* and must be casted to void* in the    |
 |                parameter list                                        |
 |                example: amscal(&val,(void*)(&ione));                 |
 *----------------------------------------------------------------------*/
void amscal(ARRAY *array, void *value);

/*----------------------------------------------------------------------*
 | initialize an array by value                           m.gee 6/01    |
 | inits the contents of a field by a given value                       |
 | ARRAY *array (input) adress of the ARRAY array                       |
 | *value (input) adress of the initvalues, this may be of type         |
 |                INT* or DOUBLE* and must be casted to void* in the    |
 |                parameter list                                        |
 |                example: aminit(&val,(void*)(&ione));                 |
 *----------------------------------------------------------------------*/
void aminit(ARRAY *array, void *value);

/*----------------------------------------------------------------------*
 | allocate and make a copy of ARRAY                      m.gee 8/00    |
 | alocates a new field in array_to of the same size and type as        |
 | in array_from. Then copies the contents from array_from to array_to  |
 | ARRAY *array_from (input) the adress of the field one wants to copy  |
 | ARRAY *array_to   (in/out) the adress of the structure where the new |
 |                            field is allocated in                     |
 | return value: adress of the new field                                |
 *----------------------------------------------------------------------*/
void *am_alloc_copy(ARRAY *array_from, ARRAY *array_to);

/*----------------------------------------------------------------------*
 | make a copy of ARRAY,                                  m.gee 6/01    |
 | Copies the contents from array_from to array_to                      |
 | ARRAY *array_from (input) the adress of the field one wants to copy  |
 | ARRAY *array_to   (in/out) the adress of the field values are copied |
 |                            to                                        |
 | user must provide fields of matching type and size!                  |
 | return value: array_to->a.iv                                         |
 *----------------------------------------------------------------------*/
void *amcopy(ARRAY *array_from, ARRAY *array_to);

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
void amadd(ARRAY *array_to, ARRAY *array_from, DOUBLE factor, INT init);

void amprint(FILE *err, ARRAY *a, INT fdim, INT sdim);

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
void *am4def(char *namstr, ARRAY4D *a, INT fdim, INT sdim, INT tdim, INT fodim, char typstr[]);

/*----------------------------------------------------------------------*
 | delete 4dimensional array                             m.gee 12/01    |
 | frees all field memory in array                                      |
 | ARRAY4D *array    (input)                adress of array-structure   |
 | no return value                                                      |
 *----------------------------------------------------------------------*/
void am4del(ARRAY4D *array);

/*----------------------------------------------------------------------*
 | initialize an 4D array by zero                           m.gee 12/01 |
 | see head of amzero                                                   |
 *----------------------------------------------------------------------*/
void am4zero(ARRAY4D *array);

/*----------------------------------------------------------------------*
 | initialize a 4D array by value                         m.gee 6/01    |
 | see head of aminit                                                   |
 *----------------------------------------------------------------------*/
void am4init(ARRAY4D *array, void *value);

/*----------------------------------------------------------------------*
 | allocate and make a copy of ARRAY4D                     m.gee 12/01  |
 | see head of am_alloc_copy                                            |
 *----------------------------------------------------------------------*/
void *am4_alloc_copy(ARRAY4D *array_from, ARRAY4D *array_to);

/*----------------------------------------------------------------------*
 | make a copy of ARRAY4D                                  m.gee 12/01  |
 | the given arrays must be of equal type and size                      |
 | see head of amcopy                                                   |
 *----------------------------------------------------------------------*/
void *am4copy(ARRAY4D *array_from, ARRAY4D *array_to);

/*----------------------------------------------------------------------*
 | redefine of ARRAY4D                                     m.gee 12/01  |
 | NOTE:  - a dimension cast like it is possible with ARRAY doesn't work|
 |        - 4-dimensional arrays with fodim=0 are not allowed           |
 |        - remember, this routine can be VERY expensive !              |
 |        - aditional space in the redefined array is set to zero       |
 |          (unlike in am4def which does NOT initialize)                |
 |        usage similar to amredef                                      |
 *----------------------------------------------------------------------*/
void *am4redef(ARRAY4D *array, INT newfdim, INT newsdim, INT newtdim, INT newfodim);


void amprint(FILE *err, ARRAY *a, INT fdim, INT sdim);


#endif
