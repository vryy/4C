static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b), (dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b), (iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))

#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "ctype.h"
#include "definitions.h"
#include "am.h"

void amdel(ARRAY *array);
void amzero(ARRAY *array);
void* am_alloc_copy(ARRAY *array_from, ARRAY *array_to);
/*!---------------------------------------------------------------------
\brief define a 1 or 2 dimensional array in structure ARRAY                                              

<pre>                                                        m.gee 8/00 
allocate a 1 or 2 - D vector of type INT or DOUBLE in the structure  
memory is allocted and pointed to in the following style:

1D arrays:
ptr[0 1 2..........................................]

2d arrays (e.g. ptr[3][5]):
ptr[0][0 1 2 3 4 5 6 7 8 9 10 11 12 13 14]
                |          |
   [1]----------           |
                           |
   [2]---------------------

This means, that the SECOND (or last) indize is continous in contrast to
fortran, where the first indize is continous. Fortran routines called on
am-allocated 2D arrays therefore operate on the transpose of the array.
</pre>
\param namstr   char*   (i)   name of array                                
\param a        ARRAY*  (i)   adress of structure ARRAY the vector lives in
\param fdim     int     (i)   first dimension of 2D vector dimension of 1D vector                       
\param sdim     int     (i)   scnd dimension of 2D vector                  
\param typstr   char[]  (i)   type of array to allocate                    
              ="IV"     allocate integer vector in a->a.iv           
              ="IA"     allocate double array  in a->a.ia            
              ="DV"     allocate integer vector in a->a.dv           
              ="DA"     allocate double array  in a->a.da            
\return void pointer to allocated memory                                               
\sa am4def()                                    

------------------------------------------------------------------------*/
void* amdef(char *namstr,ARRAY *a,int fdim, int sdim, char typstr[])
{
register int i=0;
#ifdef DEBUG 
dstrc_enter("amdef");
#endif
/*----------------------------------------------------------------------*/
strncpy(a->name,namstr,9);
a->fdim = fdim;
a->sdim = sdim;
if (strncmp("DA",typstr,2)==0) { a->Typ=cca_DA; goto next;}
if (strncmp("IA",typstr,2)==0) { a->Typ=cca_IA; goto next;}
if (strncmp("DV",typstr,2)==0) { a->Typ=cca_DV; goto next;}
if (strncmp("IV",typstr,2)==0) { a->Typ=cca_IV; goto next;}
next:
switch (a->Typ)
{
case cca_DA: /* -------------------------------------------double array */
a->a.da    = (double**)malloc((fdim*sizeof(double*)));
if (!(a->a.da))    dserror("Allocation of memory failed");
a->a.da[0] = (double*) malloc((fdim*sdim*sizeof(double)));
if (!(a->a.da[0])) dserror("Allocation of memory failed");
for (i=1; i<fdim; i++) a->a.da[i] = &(a->a.da[0][i*sdim]);
break;


case cca_IA: /* ------------------------------------------integer array */
a->a.ia    = (int**)malloc((fdim*sizeof(int*)));
if (!(a->a.ia))    dserror("Allocation of memory failed");
a->a.ia[0] = (int*) malloc((fdim*sdim*sizeof(int)));
if (!(a->a.ia[0])) dserror("Allocation of memory failed");
for (i=1; i<fdim; i++) a->a.ia[i] = &(a->a.ia[0][i*sdim]);
break;

case cca_DV: /* ------------------------------------------double vector */
a->a.dv = (double*)malloc((fdim*sdim*sizeof(double)));
if (!(a->a.dv)) dserror("Allocation of memory failed");
break;

case cca_IV: /* -----------------------------------------integer vector */
a->a.iv = (int*)malloc((fdim*sdim*sizeof(int)));
if (!(a->a.iv)) dserror("Allocation of memory failed");
break;

default:
dserror("Unknown type of array given");
}
/*----------------------------------------------------------------------*/
return((void*)(a->a.iv));
} /* end of amdef */





/*!---------------------------------------------------------------------

\brief redefine a 1 or 2 dimensional array in structure ARRAY                                              


<pre>                                                        m.gee 8/00 
changes an already allocated array a in dimensions
a typecast of the values in the array is not possible              
(no int to double or vice versa transformation)                    

a cast from iv to ia and from dv to da and vice versa is allowed     

if the new dimension of an the array is larger then the old one,     
the values inside the array are kept, the new entries are initialized
with zero                                                            
                                                                     
if the new dimension of the  array is smaller then the old one,      
values outside the new dimension are dropped                         
                                                                     
this routine handles every dimension (fdim and sdim) separately,     
that means, it does NOT behave like a fortran array. If a dimesnion gets 
smaller, it does NOT brak the spare values to the next line of the other
dimension.
it does NOT: array[3][3] ->redefine-> array[9] all values kept
it does    : array[3][3] ->redefine-> array[9] 
             values array[0][0..2] kept in array[0..2]
             values array[1..2][0..2] dropped, array[3..8]=0.0
</pre>
\param a         ARRAY* (i) adress of structure ARRAY the vector lives in
\param newfdim   int    (i) new first dimension of 2D vector dimension of 1D vector                       
\param newsdim   int    (i) new scnd dimension of 2D vector                  
\param newtypstr char[] (i) type the array shall be reallocated to                    
                 ="IV"  convert existing iv or ia to iv              
                 ="IA"  convert existing iv or ia to ia 
                 ="DV"  convert existing dv or da to dv              
                 ="DA"  convert existing dv or da to da              
\warning This routine is very expensive
\return void pointer to reallocated and reorganized memory                                               
\sa am4redef()                                    

------------------------------------------------------------------------*/
void* amredef(ARRAY *a,int newfdim, int newsdim, char newtypstr[])      
{
int i, j;
enum {amredefvoid, newIV, newIA, newDV, newDA} newtyp = amredefvoid;
int size1,size2;
ARRAY copyarray;                    

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("amredef");
#endif
/*--------------------------------------- find out the new typ of array */
if (strncmp("IV",newtypstr,2)==0){ newtyp =  newIV; goto next; }
if (strncmp("IA",newtypstr,2)==0){ newtyp =  newIA; goto next; }
if (strncmp("DV",newtypstr,2)==0){ newtyp =  newDV; goto next; }
if (strncmp("DA",newtypstr,2)==0){ newtyp =  newDA; goto next; }
next:
switch(newtyp) /*------------- what is the new array typ supposed to be */
{
case newIV:/*---------------------------------------------new typ is IV */ 
   switch(a->Typ) /*------------------------- what is the old array typ */
   {

   case cca_IV: /*------------------------------------- conversion IV to IV */
      am_alloc_copy(a,&copyarray);
      amdel(a);
      amdef(copyarray.name,a,newfdim,newsdim,"IV");
      amzero(a);
      size1 = IMIN(newfdim*newsdim,copyarray.fdim*copyarray.sdim);
      for (i=0; i<size1; a->a.iv[i++]=copyarray.a.iv[i]);
      amdel(&copyarray);
   goto end;
   
   case cca_IA: /*------------------------------------- conversion IA to IV */
      am_alloc_copy(a,&copyarray);
      amdel(a);
      amdef(copyarray.name,a,newfdim,newsdim,"IV");
      amzero(a);
      if (newfdim==1) {
          size1 = IMIN(newsdim,copyarray.sdim);
          for (i=0; i<size1; a->a.iv[i++]=copyarray.a.ia[0][i]);}  
      else {
          size1 = IMIN(newfdim,copyarray.fdim); 
          for (i=0; i<size1; a->a.iv[i++]=copyarray.a.ia[i][0]);}  
      amdel(&copyarray);
   goto end;

   default:
      dserror("conversion from integer to double or vice versa not allowed");
   goto end;
   }

case newIA:/*---------------------------------------------new typ is IA */ 
   switch(a->Typ) /*------------------------- what is the old array typ */
   {
   case cca_IV: /*------------------------------------- conversion IV to IA */
      am_alloc_copy(a,&copyarray);
      amdel(a);
      amdef(copyarray.name,a,newfdim,newsdim,"IA");
      amzero(a);
      if (copyarray.fdim==1){
         size1 = IMIN(newsdim,copyarray.sdim);
         for (i=0; i<size1; a->a.ia[0][i++]=copyarray.a.iv[i]);}  
      else{
         size1 = IMIN(newfdim,copyarray.fdim);
         for (i=0; i<size1; a->a.ia[i++][0]=copyarray.a.iv[i]);}  
      amdel(&copyarray);   
   goto end;

   case cca_IA: /*--------------------------------------conversion IA to IA */
      am_alloc_copy(a,&copyarray);
      amdel(a);
      amdef(copyarray.name,a,newfdim,newsdim,"IA");
      amzero(a);
      size1 = IMIN(newfdim,copyarray.fdim);
      size2 = IMIN(newsdim,copyarray.sdim);
      for (i=0; i<size1; i++)
      {
         for (j=0; j<size2; j++)
         {
            a->a.ia[i][j] = copyarray.a.ia[i][j];
         }
      }
      amdel(&copyarray);   
   goto end;
   default:
      dserror("conversion from integer to double or vice versa not allowed");
   goto end;
   }

case newDV:/*---------------------------------------------new typ is DV */ 
   switch(a->Typ) /*------------------------- what is the old array typ */
   {
   case cca_DV:/*---------------------------------------conversion DV to DV */
      am_alloc_copy(a,&copyarray);
      amdel(a);
      amdef(copyarray.name,a,newfdim,newsdim,"DV");
      amzero(a);
      size1 = IMIN(newfdim*newsdim,copyarray.fdim*copyarray.sdim);
      for (i=0; i<size1; a->a.dv[i++]=copyarray.a.dv[i]);
      amdel(&copyarray);
   goto end;

   case cca_DA:/*---------------------------------------conversion DA to DV */
      am_alloc_copy(a,&copyarray);
      amdel(a);
      amdef(copyarray.name,a,newfdim,newsdim,"DV");
      amzero(a);
      if (newfdim==1){
          size1 = IMIN(newsdim,copyarray.sdim);
          for (i=0; i<size1; a->a.dv[i++]=copyarray.a.da[0][i]);}  
      else{
          size1 = IMIN(newfdim,copyarray.fdim); 
          for (i=0; i<size1; a->a.dv[i++]=copyarray.a.da[i][0]);}  
      amdel(&copyarray);
   goto end;
   default:
      dserror("conversion from integer to double or vice versa not allowed");
   goto end;
   }

case newDA:/*---------------------------------------------new typ is DA */ 
   switch(a->Typ) /*------------------------- what is the old array typ */
   {
   case cca_DV:/*---------------------------------------conversion DV to DA */
      am_alloc_copy(a,&copyarray);
      amdel(a);
      amdef(copyarray.name,a,newfdim,newsdim,"DA");
      amzero(a);
      if (copyarray.fdim==1){
         size1 = IMIN(newsdim,copyarray.sdim);
         for (i=0; i<size1; a->a.da[0][i++]=copyarray.a.dv[i]);}  
      else{
         size1 = IMIN(newfdim,copyarray.fdim);
         for (i=0; i<size1; a->a.da[i++][0]=copyarray.a.dv[i]);}  
      amdel(&copyarray);   
   goto end;

   case cca_DA:/*---------------------------------------conversion DA to DA */
      am_alloc_copy(a,&copyarray);
      amdel(a);
      amdef(copyarray.name,a,newfdim,newsdim,"DA");
      amzero(a);
      size1 = IMIN(newfdim,copyarray.fdim);
      size2 = IMIN(newsdim,copyarray.sdim);
      for (i=0; i<size1; i++)
      {
         for (j=0; j<size2; j++)
         {
            a->a.da[i][j] = copyarray.a.da[i][j];
         }
      }
      amdel(&copyarray);   
   goto end;
   default:
      dserror("conversion from integer to double or vice versa not allowed");
   goto end;
   }

default:
   dserror("the new array typ is unknown");
goto end;
}
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return((void*)(a->a.iv));
} /* end of amredef */






/*!----------------------------------------------------------------------
\brief delete array allocated via am system
                                                      
<pre>                                                         m.gee 8/00    
frees the vector or array located in the ARRAY
sets 
array->name     = "DELETED"
array->Typ      = cca_XX
array->a.iv     = NULL;
array->mytracer = NULL
</pre>
\param array  ARRAY* (i/o) adress of structure ARRAY the vector lives in
\return void
\sa amdef() , amredef() , am_alloc_copy()

*----------------------------------------------------------------------*/
void amdel(ARRAY *array)
{
int i=0;
int size;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("amdel");
#endif
/*-------------------------------------------------delete name of array */
strncpy(array->name,"DELETED",9);
/*-------------------------------------------------------free the space */
switch(array->Typ)
{
case cca_DA:
   if (array->sdim)             free(array->a.da[0]);
   if (array->fdim) free(array->a.da);
break;
case cca_IA:
   if (array->sdim)             free(array->a.ia[0]);
   if (array->fdim) free(array->a.ia);
break;
case cca_DV:
   size = array->fdim * array->sdim;
   if (size) free(array->a.dv);
break;
case cca_IV:
   size = array->fdim * array->sdim;
   if (size) free(array->a.iv);
break;
default:
dserror("Unknown type of array given");
}
/*---------------------------------------------------deletes dimensions */
array->fdim=0;
array->sdim=0;
/*----------------------------------------------------- delete the type */
array->Typ = cca_XX;
/*------------------------- delete the array from the bugtracing system */
return;
} /* end of amdel */





/*!----------------------------------------------------------------------
\brief initialize an array by zero
                                                      
<pre>                                                         m.gee 8/00    
initializes the content of the ARRAY array to zero
put 0 to integer fields, 0.0 to double fields
</pre>
\param array  ARRAY* (i/o) adress of structure ARRAY the vector lives in
\return void
\sa aminit() , amscal()

*----------------------------------------------------------------------*/
void amzero(ARRAY *array)
{
register int i;
int          dim;
int         *iptr;
double      *dptr;
#ifdef DEBUG 
dstrc_enter("amzero");
#endif
/*----------------------------------------------------------------------*/
dim = (array->fdim) * (array->sdim);
switch (array->Typ)
{
case cca_DA:
   dptr = array->a.da[0];
   for (i=0; i<dim; i++) *(dptr++) = 0.0;
   break; 
case cca_DV:
   dptr = array->a.dv;
   for (i=0; i<dim; i++) *(dptr++) = 0.0;
   break;
case cca_IA:
   iptr = array->a.ia[0];
   for (i=0; i<dim; i++) *(iptr++) = 0; 
   break; 
case cca_IV:
   iptr = array->a.iv;
   for (i=0; i<dim; i++) *(iptr++) = 0; 
   break; 
default:
   dserror("Unknown type of array given");
}  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of amzero */





/*!----------------------------------------------------------------------
\brief multiply an array by a scalar
                                                      
<pre>                                                         m.gee 6/01    
scales the contents of a field by a given value 
to avoid warnings in compilation, the call has to take the form
amscal(&val,(void*)(&int_one));
or
amscal(&val,(void*)(&double_one));
depending on whether val holds an int or double array
</pre>
\param array  ARRAY* (i/o) adress of structure ARRAY the vector lives in
\param value  void*  (i)   adress of scaling value casted to void
\return void
\sa aminit() , amzero()

*----------------------------------------------------------------------*/
void amscal(ARRAY *array, void *value)
{
register int     i;
int              dim;
int             *ivalue;
double          *dvalue;
int             *iptr;
double          *dptr;
#ifdef DEBUG 
dstrc_enter("amscal");
#endif
/*----------------------------------------------------------------------*/
dim = (array->fdim) * (array->sdim);
switch (array->Typ)
{
case cca_DA:
   dvalue = (double*)value;
   dptr   = array->a.da[0];
   for (i=0; i<dim; i++) *(dptr++) *= (*dvalue);
   break; 
case cca_DV:
   dvalue = (double*)value;
   dptr = array->a.dv;
   for (i=0; i<dim; i++) *(dptr++) *= (*dvalue);
   break;
case cca_IA:
   ivalue = (int*)value;
   iptr = array->a.ia[0];
   for (i=0; i<dim; i++) *(iptr++) *= (*ivalue);
   break; 
case cca_IV:
   ivalue = (int*)value;
   iptr = array->a.iv;
   for (i=0; i<dim; i++) *(iptr++) *= (*ivalue); 
   break; 
default:
   dserror("Unknown type of array given");
}  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of amscal */


/*!----------------------------------------------------------------------
\brief initialize an array by value
                                                      
<pre>                                                         m.gee 6/01    
sets the contents of a field to a given value
to avoid warnings in compilation, the call has to take the form
aminit(&val,(void*)(&int_one));
or
aminit(&val,(void*)(&double_one));
depending on whether val holds an int or double array
</pre>
\param array  ARRAY* (i/o) adress of structure ARRAY the vector lives in
\param value  void*  (i)   adress of value casted to void
\return void
\sa amscal() , amzero() , am4zero()

*----------------------------------------------------------------------*/
void aminit(ARRAY *array, void *value)
{
register int     i;
int              dim;
int             *ivalue;
double          *dvalue;
int             *iptr;
double          *dptr;
#ifdef DEBUG 
dstrc_enter("aminit");
#endif
/*----------------------------------------------------------------------*/
dim = (array->fdim) * (array->sdim);
switch (array->Typ)
{
case cca_DA:
   dvalue = (double*)value;
   dptr   = array->a.da[0];
   for (i=0; i<dim; i++) *(dptr++) = *dvalue;
   break; 
case cca_DV:
   dvalue = (double*)value;
   dptr = array->a.dv;
   for (i=0; i<dim; i++) *(dptr++) = *dvalue;
   break;
case cca_IA:
   ivalue = (int*)value;
   iptr = array->a.ia[0];
   for (i=0; i<dim; i++) *(iptr++) = *ivalue;
   break; 
case cca_IV:
   ivalue = (int*)value;
   iptr = array->a.iv;
   for (i=0; i<dim; i++) *(iptr++) = *ivalue; 
   break; 
default:
   dserror("Unknown type of array given");
}  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of aminit */



/*!----------------------------------------------------------------------
\brief allocate and make a copy of ARRAY *array_from
                                                      
<pre>                                                         m.gee 8/00    
alocates a new field in array_to of the same size and type as        
in array_from. Then copies the contents from array_from to array_to  .
user must provide an previously NOT allocted structure array_to
</pre>
\param array_from  ARRAY* (i)   adress of existing structure to be copied
\param array_to    ARRAY* (i/o) adress of non-existing structure to be copied to
\return void*      startadress of new array
\sa amcopy() , am4_alloc_copy()

*----------------------------------------------------------------------*/
void* am_alloc_copy(ARRAY *array_from, ARRAY *array_to)
{
register int i;
int          dim;
int         *iptr_from, *iptr_to;
double      *dptr_from, *dptr_to;
#ifdef DEBUG 
dstrc_enter("am_alloc_copy");
#endif
/*----------------------------------------------------------------------*/
dim = array_from->fdim * array_from->sdim;
switch (array_from->Typ)
{
case cca_DA:
   amdef(array_from->name,array_to,array_from->fdim,array_from->sdim,"DA");
   dptr_from = array_from->a.da[0];
   dptr_to   = array_to->a.da[0];
   for (i=0; i<dim; i++) *(dptr_to++) = *(dptr_from++);
   break; 
case cca_DV:
   amdef(array_from->name,array_to,array_from->fdim,array_from->sdim,"DV");
   dptr_from = array_from->a.dv;
   dptr_to   = array_to->a.dv;
   for (i=0; i<dim; i++) *(dptr_to++) = *(dptr_from++);
   break;
case cca_IA:
   amdef(array_from->name,array_to,array_from->fdim,array_from->sdim,"IA");
   iptr_from = array_from->a.ia[0];
   iptr_to   = array_to->a.ia[0];
   for (i=0; i<dim; i++) *(iptr_to++) = *(iptr_from++);
   break; 
case cca_IV:
   amdef(array_from->name,array_to,array_from->fdim,array_from->sdim,"IV");
   iptr_from = array_from->a.iv;
   iptr_to   = array_to->a.iv;
   for (i=0; i<dim; i++) *(iptr_to++) = *(iptr_from++);
   break; 
default:
   dserror("Unknown type of array given");
}  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return((void*)(array_to->a.iv));
} /* end of am_alloc_copy */





/*!----------------------------------------------------------------------
\brief make a copy of ARRAY *array_from to *array_to
                                                      
<pre>                                                         m.gee 6/01    
Copies the contents from array_from to array_to                      
user must provide exisiting arrays of matching type and size!                  
</pre>
\param array_from  ARRAY* (i)   adress of existing structure to be copied from
\param array_to    ARRAY* (i/o) adress of existing structure to be copied to
\return void*      array_to->a.iv
\sa am_alloc_copy() , am4copy()

*----------------------------------------------------------------------*/
void* amcopy(ARRAY *array_from, ARRAY *array_to)
{
register int i;
int          dim1, dim2;
int         *iptr_from, *iptr_to;
double      *dptr_from, *dptr_to;
#ifdef DEBUG 
dstrc_enter("amcopy");
#endif
/*----------------------------------------------------------------------*/
dim1 = array_from->fdim * array_from->sdim;
dim2 = array_to->fdim   * array_to->sdim;
if (dim1 != dim2) 
   dserror("mismatching dimensions, cannot copy ARRAYs");
if (array_from->Typ != array_to->Typ)
   dserror("mismatching typ of ARRAYs, cannot copy ARRAYs");
/*----------------------------------------------------------------------*/
switch (array_from->Typ)
{
case cca_DA:
   dptr_from = array_from->a.da[0];
   dptr_to   = array_to->a.da[0];
   for (i=0; i<dim1; i++) *(dptr_to++) = *(dptr_from++);
   break; 
case cca_DV:
   dptr_from = array_from->a.dv;
   dptr_to   = array_to->a.dv;
   for (i=0; i<dim1; i++) *(dptr_to++) = *(dptr_from++);
   break;
case cca_IA:
   iptr_from = array_from->a.ia[0];
   iptr_to   = array_to->a.ia[0];
   for (i=0; i<dim1; i++) *(iptr_to++) = *(iptr_from++);
   break; 
case cca_IV:
   iptr_from = array_from->a.iv;
   iptr_to   = array_to->a.iv;
   for (i=0; i<dim1; i++) *(iptr_to++) = *(iptr_from++);
   break; 
default:
   dserror("Unknown type of array given");
}  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return((void*)(array_to->a.iv));
} /* end of amcopy */




/*!----------------------------------------------------------------------
\brief makes array_to = or += array_from * factor
                                                      
<pre>                                                         m.gee 2/02    
if init==1 array_to = array_from * factor  
if init==0 array_to = array_to + array_from * factor        
user must provide matching array-types and sufficient space in array_to                                                             
</pre>
\param array_from  ARRAY* (i)   adress of existing structure to be added from
\param array_to    ARRAY* (i/o) adress of existing structure to be added to
\param factor      double (i)   scaling factor of array_from
\param init        int    (i)   init flag
\return void
\warning In the case of int-arrays round_off takes place !
\sa aminit() , amscal() , amcopy() , am_alloc_copy()

*----------------------------------------------------------------------*/
void amadd(ARRAY *array_to, ARRAY *array_from, double factor, int init)
{
register int i,j;
int          fdim, sdim;
double     **dafrom,**dato;
double      *dvfrom, *dvto;
int        **iafrom,**iato;
int         *ivfrom, *ivto;
#ifdef DEBUG 
dstrc_enter("amadd");
#endif
/*----------------------------------------------------------------------*/
if (array_to->Typ != array_from->Typ) dserror("Mismatch of Types");
/*----------------------------------------------------------------------*/
if (init==1) amzero(array_to);
/*----------------------------------------------------------------------*/
switch (array_from->Typ)
{
case cca_DA:
   fdim = DMIN(array_to->fdim,array_from->fdim);
   sdim = DMIN(array_to->sdim,array_from->sdim);
   dafrom = array_from->a.da;
   dato   = array_to->a.da;
   for (i=0; i<fdim; i++)
   for (j=0; j<sdim; j++)
   dato[i][j] += dafrom[i][j] * factor;
break; 
case cca_DV:
   fdim   = array_to->fdim * array_to->sdim;
   sdim   = array_from->fdim * array_from->sdim;
   fdim   = DMIN(fdim,sdim);
   dvfrom = array_from->a.dv;
   dvto   = array_to->a.dv;
   for (i=0; i<fdim; i++)
   dvto[i] += dvfrom[i] * factor;
break;
case cca_IA:
   fdim = DMIN(array_to->fdim,array_from->fdim);
   sdim = DMIN(array_to->sdim,array_from->sdim);
   iafrom = array_from->a.ia;
   iato   = array_to->a.ia;
   for (i=0; i<fdim; i++)
   for (j=0; j<sdim; j++)
   iato[i][j] += (int)((double)iafrom[i][j] * factor);
break; 
case cca_IV:
   fdim   = array_to->fdim * array_to->sdim;
   sdim   = array_from->fdim * array_from->sdim;
   fdim   = DMIN(fdim,sdim);
   ivfrom = array_from->a.iv;
   ivto   = array_to->a.iv;
   for (i=0; i<fdim; i++)
   ivto[i] += (int)((double)ivfrom[i] * factor);
break; 
default:
   dserror("Unknown type of array given");
}  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of amadd */




/*!---------------------------------------------------------------------

\brief define a 3 or 4 dimensional array in structure ARRAY4D                                              

<pre>                                                        m.gee 12/01
allocate a 3 or 4 - D vector of type INT or DOUBLE in the structure  
In the case of 3D arrays, the fourth dimension fodim must be zero.

memory is allocted and pointed to in the following style:
3d arrays (e.g. ptr[2][3][3]):
ptr[0][0][0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17]
      [1]--------
      [2]--------------
   [1]-
      [0]--------------------
      [1]----------------------------
      [2]-------------------------------------

(same logic for 4d arrays)

This means, that the LAST indize is continous in contrast to
fortran, where the first indize is continous. Fortran routines called on
am-allocated 3/4D arrays therefore operate on the transpose of the array.
</pre>
\param namstr    char*    (i)   name of array                                
\param a         ARRAY4D* (i/o) adress of structure ARRAY4D the vector lives in
\param fdim      int      (i)   first dimension of 2D vector dimension of 1D vector                       
\param sdim      int      (i)   scnd dimension of 2D vector                  
\param tdim      int      (i)   thrd dimension of 2D vector                  
\param fodim     int      (i)   fourth dimension of 2D vector                  
\param typstr    char[]   (i)   type of array to allocate                    
                 ="I3"     allocate integer vector in a->a.i3           
                 ="I4"     allocate double array  in a->a.i4            
                 ="D3"     allocate integer vector in a->a.d3           
                 ="D4"     allocate double array  in a->a.d4            
\return void pointer to allocated memory                                               
\sa amdef() , am4redef                                   

------------------------------------------------------------------------*/
void* am4def(char    *namstr,
             ARRAY4D *a,
             int      fdim, 
             int      sdim,
             int      tdim,
             int      fodim, 
             char     typstr[])
{
register int i;
int          endloop;
#ifdef DEBUG 
dstrc_enter("am4def");
#endif
/*----------------------------------------------------------------------*/
strncpy(a->name,namstr,9);
a->fdim =  fdim;
a->sdim =  sdim;
a->tdim =  tdim;
a->fodim = fodim;
if (strncmp("D3",typstr,2)==0) { a->Typ=cca_D3; goto next;}
if (strncmp("D4",typstr,2)==0) { a->Typ=cca_D4; goto next;}
if (strncmp("I3",typstr,2)==0) { a->Typ=cca_I3; goto next;}
if (strncmp("I4",typstr,2)==0) { a->Typ=cca_I4; goto next;}
next:
switch (a->Typ)
{
case cca_D3: /* --------------------------------------------double D3 array */
if (fodim != 0) dserror("Illegal fourth dimension in call to am4def");
a->a.d3       = (double***)malloc((fdim*sizeof(double**)));
if (!(a->a.d3))       dserror("Allocation of memory failed");
a->a.d3[0]    = (double**) malloc((fdim*sdim*sizeof(double*)));
if (!(a->a.d3[0]))    dserror("Allocation of memory failed");
a->a.d3[0][0] = (double*)  malloc((fdim*sdim*tdim*sizeof(double)));
if (!(a->a.d3[0][0])) dserror("Allocation of memory failed");

for (i=1; i<fdim; i++)    a->a.d3[i]    = &(a->a.d3[0][i*sdim]);
endloop=fdim*sdim;
for (i=1; i<endloop; i++) a->a.d3[0][i] = &(a->a.d3[0][0][i*tdim]);
break;

case cca_I3: /* ----------------------------------------------int I3 array */
if (fodim != 0) dserror("Illegal fourth dimension in call to am4def");
a->a.i3       = (int***)malloc((fdim*sizeof(int**)));
if (!(a->a.i3))       dserror("Allocation of memory failed");
a->a.i3[0]    = (int**) malloc((fdim*sdim*sizeof(int*)));
if (!(a->a.i3[0]))    dserror("Allocation of memory failed");
a->a.i3[0][0] = (int*)  malloc((fdim*sdim*tdim*sizeof(int)));
if (!(a->a.i3[0][0])) dserror("Allocation of memory failed");

for (i=1; i<fdim; i++)    a->a.i3[i]    = &(a->a.i3[0][i*sdim]);
endloop=fdim*sdim;
for (i=1; i<endloop; i++) a->a.i3[0][i] = &(a->a.i3[0][0][i*tdim]);
break;

case cca_D4: /* --------------------------------------------double D4 array */
a->a.d4          = (double****)malloc((fdim*sizeof(double***)));
if (!(a->a.d4))          dserror("Allocation of memory failed");
a->a.d4[0]       = (double***) malloc((fdim*sdim*sizeof(double**)));
if (!(a->a.d4[0]))       dserror("Allocation of memory failed");
a->a.d4[0][0]    = (double**)  malloc((fdim*sdim*tdim*sizeof(double*)));
if (!(a->a.d4[0][0]))    dserror("Allocation of memory failed");
a->a.d4[0][0][0] = (double*)   malloc((fdim*sdim*tdim*fodim*sizeof(double)));
if (!(a->a.d4[0][0][0])) dserror("Allocation of memory failed");

for (i=1; i<fdim; i++)    a->a.d4[i]       = &(a->a.d4[0][i*sdim]);
endloop=fdim*sdim;
for (i=1; i<endloop; i++) a->a.d4[0][i]    = &(a->a.d4[0][0][i*tdim]);
endloop=fdim*sdim*tdim;
for (i=1; i<endloop; i++) a->a.d4[0][0][i] = &(a->a.d4[0][0][0][i*fodim]);
break;

case cca_I4: /* ----------------------------------------------int I4 array */
a->a.i4          = (int****)malloc((fdim*sizeof(int***)));
if (!(a->a.i4))          dserror("Allocation of memory failed");
a->a.i4[0]       = (int***) malloc((fdim*sdim*sizeof(int**)));
if (!(a->a.i4[0]))       dserror("Allocation of memory failed");
a->a.i4[0][0]    = (int**)  malloc((fdim*sdim*tdim*sizeof(int*)));
if (!(a->a.i4[0][0]))    dserror("Allocation of memory failed");
a->a.i4[0][0][0] = (int*)   malloc((fdim*sdim*tdim*fodim*sizeof(int)));
if (!(a->a.i4[0][0][0])) dserror("Allocation of memory failed");

for (i=1; i<fdim; i++)    a->a.i4[i]       = &(a->a.i4[0][i*sdim]);
endloop=fdim*sdim;
for (i=1; i<endloop; i++) a->a.i4[0][i]    = &(a->a.i4[0][0][i*tdim]);
endloop=fdim*sdim*tdim;
for (i=1; i<endloop; i++) a->a.i4[0][0][i] = &(a->a.i4[0][0][0][i*fodim]);
break;

default:
dserror("Unknown type of array given");
}
/*------------------- make report about new array to bugtraceing system */
return((void*)(a->a.d3));
} /* end of am4def */




/*!----------------------------------------------------------------------
\brief delete 3/4D array
                                                      
<pre>                                                         m.gee 12/01    
frees all field memory in array
sets 
array->name     = "DELETED"
array->a.iv     = NULL;
array->mytracer = NULL
</pre>
\param array  ARRAY4D* (i/o) adress of structure ARRAY the vector lives in
\return void
\sa am4def() , am4redef() , amdel()

*----------------------------------------------------------------------*/
void am4del(ARRAY4D *array)
{
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("am4del");
#endif
/*-------------------------------------------------delete name of array */
strncpy(array->name,"DELETED",9);
/*---------------------------------------------------deletes dimensions */
array->fdim=0;
array->sdim=0;
array->tdim=0;
array->fodim=0;
/*-------------------------------------------------------free the space */
switch(array->Typ)
{
case cca_D3:
   free(array->a.d3[0][0]);
   free(array->a.d3[0]);
   free(array->a.d3);
break;
case cca_I3:
   free(array->a.i3[0][0]);
   free(array->a.i3[0]);
   free(array->a.i3);
break;
case cca_D4:
   free(array->a.d4[0][0][0]);
   free(array->a.d4[0][0]);
   free(array->a.d4[0]);
   free(array->a.d4);
break;
case cca_I4:
   free(array->a.i4[0][0][0]);
   free(array->a.i4[0][0]);
   free(array->a.i4[0]);
   free(array->a.i4);
break;
default:
dserror("Unknown type of array given");
}
/*----------------------------------------------------- delete the type */
array->Typ = cca_XX4D;
/*------------------------- delete the array from the bugtracing system */
return;
} /* end of am4del */



/*!----------------------------------------------------------------------
\brief initialize an ARRAY4D array by zero
                                                      
<pre>                                                         m.gee 12/01    
initializes the content of the ARRAY4D array to zero
put 0 to integer fields, 0.0 to double fields
</pre>
\param array  ARRAY* (i/o) adress of structure ARRAY the vector lives in
\return void
\sa amzero()

*----------------------------------------------------------------------*/
void am4zero(ARRAY4D *array)
{
register int i;
int          dim;
int         *iptr;
double      *dptr;
#ifdef DEBUG 
dstrc_enter("am4zero");
#endif
/*----------------------------------------------------------------------*/
switch (array->Typ)
{
case cca_D3:
   dim = (array->fdim) * (array->sdim) * (array->tdim);
   dptr = array->a.d3[0][0];
   for (i=0; i<dim; i++) *(dptr++) = 0.0;
break; 
case cca_D4:
   dim = (array->fdim) * (array->sdim) * (array->tdim) * (array->fodim);
   dptr = array->a.d4[0][0][0];
   for (i=0; i<dim; i++) *(dptr++) = 0.0;
break;
case cca_I3:
   dim = (array->fdim) * (array->sdim) * (array->tdim);
   iptr = array->a.i3[0][0];
   for (i=0; i<dim; i++) *(iptr++) = 0; 
break; 
case cca_I4:
   dim = (array->fdim) * (array->sdim) * (array->tdim) * (array->fodim);
   iptr = array->a.i4[0][0][0];
   for (i=0; i<dim; i++) *(iptr++) = 0; 
break; 
default:
   dserror("Unknown type of array given");
}  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of am4zero */




/*!----------------------------------------------------------------------
\brief initialize an ARRAY4D array by value
                                                      
<pre>                                                         m.gee 6/01    
sets the contents of a field to a given value
to avoid warnings in compilation, the call has to take the form
am4init(&val,(void*)(&int_one));
or
am4init(&val,(void*)(&double_one));
depending on whether val holds an int or double array
</pre>
\param array  ARRAY4D* (i/o) adress of structure ARRAY4D the vector lives in
\param value  void*  (i)   adress of value casted to void
\return void
\sa aminit() , am4zero()

*----------------------------------------------------------------------*/
void am4init(ARRAY4D *array, void *value)
{
register int     i;
int              dim;
int             *ivalue;
double          *dvalue;
int             *iptr;
double          *dptr;
#ifdef DEBUG 
dstrc_enter("am4init");
#endif
/*----------------------------------------------------------------------*/
switch (array->Typ)
{
case cca_D3:
   dim = (array->fdim) * (array->sdim) * (array->tdim);
   dvalue = (double*)value;
   dptr   = array->a.d3[0][0];
   for (i=0; i<dim; i++) *(dptr++) = *dvalue;
break; 
case cca_D4:
   dim = (array->fdim) * (array->sdim) * (array->tdim) * (array->fodim);
   dvalue = (double*)value;
   dptr = array->a.d4[0][0][0];
   for (i=0; i<dim; i++) *(dptr++) = *dvalue;
break;
case cca_I3:
   dim = (array->fdim) * (array->sdim) * (array->tdim);
   ivalue = (int*)value;
   iptr = array->a.i3[0][0];
   for (i=0; i<dim; i++) *(iptr++) = *ivalue;
break; 
case cca_I4:
   dim = (array->fdim) * (array->sdim) * (array->tdim) * (array->fodim);
   ivalue = (int*)value;
   iptr = array->a.i4[0][0][0];
   for (i=0; i<dim; i++) *(iptr++) = *ivalue; 
break; 
default:
   dserror("Unknown type of array given");
}  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of am4init */



/*!----------------------------------------------------------------------
\brief allocate and make a copy of ARRAY4D *array_from
                                                      
<pre>                                                         m.gee 12/01    
alocates a new field in array_to of the same size and type as        
in array_from. Then copies the contents from array_from to array_to  .
user must provide an previously NOT allocted structure array_to
</pre>
\param array_from  ARRAY4D* (i)   adress of existing structure to be copied
\param array_to    ARRAY4D* (i/o) adress of non-existing structure to be copied to
\return void*      startadress of new array
\sa am_alloc_copy() , am4copy() , amcopy()

*----------------------------------------------------------------------*/
void* am4_alloc_copy(ARRAY4D *array_from, ARRAY4D *array_to)
{
register int i;
int          dim;
int         *iptr_from, *iptr_to;
double      *dptr_from, *dptr_to;
#ifdef DEBUG 
dstrc_enter("am4_alloc_copy");
#endif
/*----------------------------------------------------------------------*/
switch (array_from->Typ)
{
case cca_D3:
   dim = array_from->fdim * array_from->sdim * array_from->tdim;
   am4def(array_from->name,
          array_to,
          array_from->fdim,
          array_from->sdim,
          array_from->tdim,
          0,
          "D3");
   dptr_from = array_from->a.d3[0][0];
   dptr_to   = array_to->a.d3[0][0];
   for (i=0; i<dim; i++) *(dptr_to++) = *(dptr_from++);
break; 
case cca_D4:
   dim = array_from->fdim * array_from->sdim * array_from->tdim * array_from->fodim;
   am4def(array_from->name,
          array_to,
          array_from->fdim,
          array_from->sdim,
          array_from->tdim,
          array_from->fodim,
          "D4");
   dptr_from = array_from->a.d4[0][0][0];
   dptr_to   = array_to->a.d4[0][0][0];
   for (i=0; i<dim; i++) *(dptr_to++) = *(dptr_from++);
break;
case cca_I3:
   dim = array_from->fdim * array_from->sdim * array_from->tdim;
   am4def(array_from->name,
          array_to,
          array_from->fdim,
          array_from->sdim,
          array_from->tdim,
          0,
          "I3");
   iptr_from = array_from->a.i3[0][0];
   iptr_to   = array_to->a.i3[0][0];
   for (i=0; i<dim; i++) *(iptr_to++) = *(iptr_from++);
break; 
case cca_I4:
   dim = array_from->fdim * array_from->sdim * array_from->tdim * array_from->fodim;
   am4def(array_from->name,
          array_to,
          array_from->fdim,
          array_from->sdim,
          array_from->tdim,
          array_from->fodim,
          "I4");
   iptr_from = array_from->a.i4[0][0][0];
   iptr_to   = array_to->a.i4[0][0][0];
   for (i=0; i<dim; i++) *(iptr_to++) = *(iptr_from++);
break; 
default:
   dserror("Unknown type of array given");
}  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return((void*)(array_to->a.d3));
} /* end of am4_alloc_copy */




/*!----------------------------------------------------------------------
\brief make a copy of ARRAY4D *array_from to *array_to
                                                      
<pre>                                                         m.gee 12/01    
Copies the contents from array_from to array_to                      
user must provide exisiting arrays of matching type and size!                  
</pre>
\param array_from  ARRAY4D* (i)   adress of existing structure to be copied from
\param array_to    ARRAY4D* (i/o) adress of existing structure to be copied to
\return void*      array_to->a.i3
\sa am_alloc_copy() , am4_alloc_copy() , amcopy()

*----------------------------------------------------------------------*/
void* am4copy(ARRAY4D *array_from, ARRAY4D *array_to)
{
register int i;
int          dim,dimnew;
int         *iptr_from, *iptr_to;
double      *dptr_from, *dptr_to;
#ifdef DEBUG 
dstrc_enter("am4copy");
#endif
/*--------------------------------------------------------- check types */
if (array_from->Typ != array_to->Typ)
   dserror("mismatching typ of ARRAY4Ds, cannot copy ARRAY4Ds");
/*----------------------------------------------------------------------*/
switch (array_from->Typ)
{
case cca_D3:
   dim    = array_from->fdim * array_from->sdim * array_from->tdim;
   dimnew = array_to->fdim * array_to->sdim * array_to->tdim;
   if (dim != dimnew) dserror("mismatching dimensions, cannot copy ARRAY4D");
   dptr_from = array_from->a.d3[0][0];
   dptr_to   = array_to->a.d3[0][0];
   for (i=0; i<dim; i++) *(dptr_to++) = *(dptr_from++);
break; 
case cca_D4:
   dim    = array_from->fdim * array_from->sdim * array_from->tdim * array_from->fodim;
   dimnew = array_to->fdim * array_to->sdim * array_to->tdim * array_to->fodim;
   if (dim != dimnew) dserror("mismatching dimensions, cannot copy ARRAY4D");
   dptr_from = array_from->a.d4[0][0][0];
   dptr_to   = array_to->a.d4[0][0][0];
   for (i=0; i<dim; i++) *(dptr_to++) = *(dptr_from++);
break;
case cca_I3:
   dim    = array_from->fdim * array_from->sdim * array_from->tdim;
   dimnew = array_to->fdim * array_to->sdim * array_to->tdim;
   if (dim != dimnew) dserror("mismatching dimensions, cannot copy ARRAY4D");
   iptr_from = array_from->a.i3[0][0];
   iptr_to   = array_to->a.i3[0][0];
   for (i=0; i<dim; i++) *(iptr_to++) = *(iptr_from++);
break; 
case cca_I4:
   dim = array_from->fdim * array_from->sdim * array_from->tdim * array_from->fodim;
   dimnew = array_to->fdim * array_to->sdim * array_to->tdim * array_to->fodim;
   if (dim != dimnew) dserror("mismatching dimensions, cannot copy ARRAY4D");
   iptr_from = array_from->a.i4[0][0][0];
   iptr_to   = array_to->a.i4[0][0][0];
   for (i=0; i<dim; i++) *(iptr_to++) = *(iptr_from++);
break; 
default:
   dserror("Unknown type of array given");
}  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return((void*)(array_to->a.d3));
} /* end of am4copy */




/*!---------------------------------------------------------------------

\brief redefine a 3 or 4 dimensional array in structure ARRAY4D                                              


<pre>                                                        m.gee 12/01 
changes an already allocated array a in dimensions
a typecast of the values in the array is not possible              
(no int to double or vice versa transformation)                    

a cast between 3D and 4D arrays is NOT allowed     

if the new dimension of an the array is larger then the old one,     
the values inside the array are kept, the new entries are initialized
with zero                                                            
                                                                     
if the new dimension of the  array is smaller then the old one,      
values outside the new dimension are dropped                         
                                                                     
this routine handles every dimension (fdim and sdim) separately,     
that means, it does NOT behave like a fortran array. If a dimension gets 
smaller, it does NOT break the spare values to the next line of the other
dimension.

Additional memory in in redefined array is set to zero

</pre>
\param a         ARRAY4D* (i/o) adress of structure ARRAY the vector lives in
\param newfdim   int      (i)   new first dimension of 3/4D array dimension                       
\param newsdim   int      (i)   new scnd dimension of 3/4D array                  
\param newtdim   int      (i)   new thrd dimension of 3/4D array                  
\param newfodim  int      (i)   new fourth dimension of 3/4D array                  
\warning This routine can be very expensive
\return void pointer to reallocated and reorganized memory                                               
\sa amredef()                                    

------------------------------------------------------------------------*/
void* am4redef(ARRAY4D *array, 
               int newfdim, 
               int newsdim, 
               int newtdim,
               int newfodim)
{
register int i,j,k,l;
int          size1,size2,size3,size4;
ARRAY4D copyarray;
#ifdef DEBUG 
dstrc_enter("am4redef");
#endif
/*----------------------------------------------------------------------*/
switch (array->Typ)
{
case cca_D3:
   am4_alloc_copy(array,&copyarray);
   am4del(array);
   am4def(copyarray.name,array,newfdim,newsdim,newtdim,0,"D3");
   am4zero(array);
   size1 = IMIN(copyarray.fdim,newfdim);
   size2 = IMIN(copyarray.sdim,newsdim);
   size3 = IMIN(copyarray.tdim,newtdim);
   for (i=0; i<size1; i++)
   for (j=0; j<size2; j++)
   for (k=0; k<size3; k++)
   array->a.d3[i][j][k] = copyarray.a.d3[i][j][k];
   am4del(&copyarray);
break;
case cca_D4:
   am4_alloc_copy(array,&copyarray);
   am4del(array);
   am4def(copyarray.name,array,newfdim,newsdim,newtdim,newfodim,"D4");
   am4zero(array);
   size1 = IMIN(copyarray.fdim,newfdim);
   size2 = IMIN(copyarray.sdim,newsdim);
   size3 = IMIN(copyarray.tdim,newtdim);
   size4 = IMIN(copyarray.fodim,newfodim);
   for (i=0; i<size1; i++)
   for (j=0; j<size2; j++)
   for (k=0; k<size3; k++)
   for (l=0; l<size4; l++)
   array->a.d4[i][j][k][l] = copyarray.a.d4[i][j][k][l];
   am4del(&copyarray);
break;
case cca_I3:
   am4_alloc_copy(array,&copyarray);
   am4del(array);
   am4def(copyarray.name,array,newfdim,newsdim,newtdim,0,"I3");
   am4zero(array);
   size1 = IMIN(copyarray.fdim,newfdim);
   size2 = IMIN(copyarray.sdim,newsdim);
   size3 = IMIN(copyarray.tdim,newtdim);
   for (i=0; i<size1; i++)
   for (j=0; j<size2; j++)
   for (k=0; k<size3; k++)
   array->a.i3[i][j][k] = copyarray.a.i3[i][j][k];
   am4del(&copyarray);
break;
case cca_I4:
   am4_alloc_copy(array,&copyarray);
   am4del(array);
   am4def(copyarray.name,array,newfdim,newsdim,newtdim,newfodim,"I4");
   am4zero(array);
   size1 = IMIN(copyarray.fdim,newfdim);
   size2 = IMIN(copyarray.sdim,newsdim);
   size3 = IMIN(copyarray.tdim,newtdim);
   size4 = IMIN(copyarray.fodim,newfodim);
   for (i=0; i<size1; i++)
   for (j=0; j<size2; j++)
   for (k=0; k<size3; k++)
   for (l=0; l<size4; l++)
   array->a.i4[i][j][k][l] = copyarray.a.i4[i][j][k][l];
   am4del(&copyarray);
break;
default:
   dserror("Unknown type of array given");
break;
}  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return((void*)(array->a.d3));
} /* end of am4redef */


/*! @} (documentation module close)*/

int dserror(char string[])
{
printf("%s\n",string);
fflush(stdout);
exit(EXIT_FAILURE);
return;
}
