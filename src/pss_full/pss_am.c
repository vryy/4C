#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | tracing variables                                                    |
 | defined in pss_ds.c                                                  |
 *----------------------------------------------------------------------*/
#ifdef DEBUG
extern struct _TRACE         trace;
#endif



/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | global variable  num_byte_allocated                                  |
 | long int num_byte_allocated                                          |
 *----------------------------------------------------------------------*/
#ifdef DEBUG
long int num_byte_allocated;
#endif
/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | in debug mode all memory is allocated one double bigger then         |
 | asked for, so the size of what was allocated can be stored           |
 *----------------------------------------------------------------------*/
#ifdef DEBUG
#define DWORD  (sizeof(double)) /* Be sure it's not smaller than double */
#define MYBYTE (sizeof(unsigned char)) /*  exactly one byte (hopefully) */ 
#endif


/*----------------------------------------------------------------------*
 | redefinition of malloc DEBUG version                       m.gee 2/02|
 *----------------------------------------------------------------------*/
#ifdef DEBUG
void *MALLOC(int size)
{
char *buf;

buf = (char*)malloc((size_t)(size+DWORD));
if (buf)
{
   ((int*)buf)[0] = size;
   num_byte_allocated += (size+DWORD);
   return (void*)(buf+DWORD);
}
else return (void*)(NULL);
}/* end of MALLOC */
/*----------------------------------------------------------------------*
 | redefinition of malloc FAST version                        m.gee 2/02|
 *----------------------------------------------------------------------*/
#else
void *MALLOC(int size)
{
void *buf;
buf = malloc((size_t)size);
return (void*)(buf);
}/* end of MALLOC */
#endif



/*----------------------------------------------------------------------*
 | redefinition of calloc DEBUG version                   m.gee 2/02    |
 *----------------------------------------------------------------------*/
#ifdef DEBUG
void *CALLOC(int num, int size)
{
char *buf;

buf = (char*)calloc((size_t)(num*size+DWORD),(size_t)MYBYTE);
if (buf)
{
   ((int*)buf)[0] = size*num;
   num_byte_allocated += (num*size+DWORD);
   return (void*)(buf+DWORD);
}
else return (void*)(NULL);
}/* end of CALLOC */
/*----------------------------------------------------------------------*
 | redefinition of calloc FAST version                    m.gee 2/02    |
 *----------------------------------------------------------------------*/
#else
void *CALLOC(int num, int size)
{
void *buf;
buf = calloc((size_t)num,(size_t)size);
return (void *)(buf);
}/* end of CALLOC */
#endif



/*----------------------------------------------------------------------*
 | redefinition of realloc DEBUG version                  m.gee 2/02    |
 *----------------------------------------------------------------------*/
#ifdef DEBUG
void *REALLOC(void *oldptr, int size)
{
int   n;
char *buf = ((char*)oldptr) - DWORD;
if (!oldptr) dserror("Tried to realloc NULL pointer");
if (!buf)    dserror("Tried to realloc NULL-DWORD pointer");
n = ((int*)buf)[0];
num_byte_allocated -= (n+DWORD);

buf = (char*)realloc(buf,(size_t)(size+DWORD));
if (buf)
{
   ((int*)buf)[0] = size;
   num_byte_allocated += (size+DWORD);
   return (void*)(buf+DWORD);
}
else return (void*)(NULL);
}/* end of REALLOC */
/*----------------------------------------------------------------------*
 | redefinition of realloc FAST version                  m.gee 2/02    |
 *----------------------------------------------------------------------*/
#else
void *REALLOC(void *oldptr, int size)
{
void *buf;
buf = realloc(oldptr,(size_t)size);
return (void*)(buf);
}/* end of REALLOC */
#endif



/*----------------------------------------------------------------------*
 | redefinition of free DEBUG version                     m.gee 2/02    |
 *----------------------------------------------------------------------*/
#ifdef DEBUG
void *FREE(void *oldptr)
{
int   n;
char *p = ((char*)oldptr) - DWORD;
if (!oldptr) dserror("Tried to free NULL pointer");
if (!p)      dserror("Tried to free NULL-DWORD pointer");
n = ((int*)p)[0];
*((int*) p) = 0;
num_byte_allocated -= (n+DWORD);
free(p);
return (oldptr=NULL);
}/* end of FREE */
/*----------------------------------------------------------------------*
 | redefinition of free FAST version                      m.gee 2/02    |
 *----------------------------------------------------------------------*/
#else
void *FREE(void *oldptr)
{
free(oldptr);
return (oldptr=NULL);
}/* end of FREE */
#endif



/*----------------------------------------------------------------------*
 | define array                                           m.gee 8/00    |
 *----------------------------------------------------------------------*/
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
if (strncmp("DA",typstr,2)==0) { a->Typ=DA; goto next;}
if (strncmp("IA",typstr,2)==0) { a->Typ=IA; goto next;}
if (strncmp("DV",typstr,2)==0) { a->Typ=DV; goto next;}
if (strncmp("IV",typstr,2)==0) { a->Typ=IV; goto next;}
next:
switch (a->Typ)
{
case DA: /* -----------------------------------------------double array */
a->a.da    = (double**)MALLOC((fdim*sizeof(double*)));
if (!(a->a.da))    dserror("Allocation of memory failed");
a->a.da[0] = (double*) MALLOC((fdim*sdim*sizeof(double)));
if (!(a->a.da[0])) dserror("Allocation of memory failed");
for (i=1; i<fdim; i++) a->a.da[i] = &(a->a.da[0][i*sdim]);
break;


case IA: /* ----------------------------------------------integer array */
a->a.ia    = (int**)MALLOC((fdim*sizeof(int*)));
if (!(a->a.ia))    dserror("Allocation of memory failed");
a->a.ia[0] = (int*) MALLOC((fdim*sdim*sizeof(int)));
if (!(a->a.ia[0])) dserror("Allocation of memory failed");
for (i=1; i<fdim; i++) a->a.ia[i] = &(a->a.ia[0][i*sdim]);
break;

case DV: /* ----------------------------------------------double vector */
a->a.dv = (double*)MALLOC((fdim*sdim*sizeof(double)));
if (!(a->a.dv)) dserror("Allocation of memory failed");
break;

case IV: /* ---------------------------------------------integer vector */
a->a.iv = (int*)MALLOC((fdim*sdim*sizeof(int)));
if (!(a->a.iv)) dserror("Allocation of memory failed");
break;

default:
dserror("Unknown type of array given");
}
/*------------------- make report about new array to bugtraceing system */
#ifdef DEBUG 
if (trace.trace_on==1)
{
dstracereport(a);
}
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
return((void*)(a->a.iv));
} /* end of amdef */



/*----------------------------------------------------------------------*
 | redefine array                                         m.gee 8/00    |
   amredef has to be rewritten in the style amzero and aminit are
 *----------------------------------------------------------------------*/
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

   case IV: /*------------------------------------- conversion IV to IV */
      am_alloc_copy(a,&copyarray);
      amdel(a);
      amdef(copyarray.name,a,newfdim,newsdim,"IV");
      amzero(a);
      size1 = IMIN(newfdim*newsdim,copyarray.fdim*copyarray.sdim);
      for (i=0; i<size1; a->a.iv[i++]=copyarray.a.iv[i]);
      amdel(&copyarray);
   goto end;
   
   case IA: /*------------------------------------- conversion IA to IV */
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
   case IV: /*------------------------------------- conversion IV to IA */
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

   case IA: /*--------------------------------------conversion IA to IA */
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
   case DV:/*---------------------------------------conversion DV to DV */
      am_alloc_copy(a,&copyarray);
      amdel(a);
      amdef(copyarray.name,a,newfdim,newsdim,"DV");
      amzero(a);
      size1 = IMIN(newfdim*newsdim,copyarray.fdim*copyarray.sdim);
      for (i=0; i<size1; a->a.dv[i++]=copyarray.a.dv[i]);
      amdel(&copyarray);
   goto end;

   case DA:/*---------------------------------------conversion DA to DV */
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
   case DV:/*---------------------------------------conversion DV to DA */
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

   case DA:/*---------------------------------------conversion DA to DA */
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




/*----------------------------------------------------------------------*
 | delete         array                                   m.gee 8/00    |
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
case DA:
   if (array->sdim)             FREE(array->a.da[0]);
   if (array->fdim) array->a.da=FREE(array->a.da);
break;
case IA:
   if (array->sdim)             FREE(array->a.ia[0]);
   if (array->fdim) array->a.ia=FREE(array->a.ia);
break;
case DV:
   size = array->fdim * array->sdim;
   if (size) array->a.dv=FREE(array->a.dv);
break;
case IV:
   size = array->fdim * array->sdim;
   if (size) array->a.iv=FREE(array->a.iv);
break;
default:
dserror("Unknown type of array given");
}
/*---------------------------------------------------deletes dimensions */
array->fdim=0;
array->sdim=0;
/*----------------------------------------------------- delete the type */
array->Typ = XX;
/*------------------------- delete the array from the bugtracing system */
#ifdef DEBUG 
if (trace.trace_on==1)
{
    trace.arrays[array->place_in_trace]=NULL;
    trace.num_arrays--;
    array->place_in_trace=0;
}                     
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
return;
} /* end of amdel */





/*----------------------------------------------------------------------*
 | initialize an array by zero                            m.gee 8/00    |
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
case DA:
   dptr = array->a.da[0];
   for (i=0; i<dim; i++) *(dptr++) = 0.0;
   break; 
case DV:
   dptr = array->a.dv;
   for (i=0; i<dim; i++) *(dptr++) = 0.0;
   break;
case IA:
   iptr = array->a.ia[0];
   for (i=0; i<dim; i++) *(iptr++) = 0; 
   break; 
case IV:
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



/*----------------------------------------------------------------------*
 | initialize an array by value                           m.gee 6/01    |
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
case DA:
   dvalue = (double*)value;
   dptr   = array->a.da[0];
   for (i=0; i<dim; i++) *(dptr++) = *dvalue;
   break; 
case DV:
   dvalue = (double*)value;
   dptr = array->a.dv;
   for (i=0; i<dim; i++) *(dptr++) = *dvalue;
   break;
case IA:
   ivalue = (int*)value;
   iptr = array->a.ia[0];
   for (i=0; i<dim; i++) *(iptr++) = *ivalue;
   break; 
case IV:
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



/*----------------------------------------------------------------------*
 | allocate and make a copy of ARRAY                      m.gee 8/00    |
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
case DA:
   amdef(array_from->name,array_to,array_from->fdim,array_from->sdim,"DA");
   dptr_from = array_from->a.da[0];
   dptr_to   = array_to->a.da[0];
   for (i=0; i<dim; i++) *(dptr_to++) = *(dptr_from++);
   break; 
case DV:
   amdef(array_from->name,array_to,array_from->fdim,array_from->sdim,"DV");
   dptr_from = array_from->a.dv;
   dptr_to   = array_to->a.dv;
   for (i=0; i<dim; i++) *(dptr_to++) = *(dptr_from++);
   break;
case IA:
   amdef(array_from->name,array_to,array_from->fdim,array_from->sdim,"IA");
   iptr_from = array_from->a.ia[0];
   iptr_to   = array_to->a.ia[0];
   for (i=0; i<dim; i++) *(iptr_to++) = *(iptr_from++);
   break; 
case IV:
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



/*----------------------------------------------------------------------*
 | make a copy of ARRAY,                                  m.gee 6/01    |
 | user must provide sufficient space & typ                             |
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
case DA:
   dptr_from = array_from->a.da[0];
   dptr_to   = array_to->a.da[0];
   for (i=0; i<dim1; i++) *(dptr_to++) = *(dptr_from++);
   break; 
case DV:
   dptr_from = array_from->a.dv;
   dptr_to   = array_to->a.dv;
   for (i=0; i<dim1; i++) *(dptr_to++) = *(dptr_from++);
   break;
case IA:
   iptr_from = array_from->a.ia[0];
   iptr_to   = array_to->a.ia[0];
   for (i=0; i<dim1; i++) *(iptr_to++) = *(iptr_from++);
   break; 
case IV:
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



/*----------------------------------------------------------------------*
 | define 4D array                                       m.gee 12/01    |
 *----------------------------------------------------------------------*/
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
if (strncmp("D3",typstr,2)==0) { a->Typ=D3; goto next;}
if (strncmp("D4",typstr,2)==0) { a->Typ=D4; goto next;}
if (strncmp("I3",typstr,2)==0) { a->Typ=I3; goto next;}
if (strncmp("I4",typstr,2)==0) { a->Typ=I4; goto next;}
next:
switch (a->Typ)
{
case D3: /* --------------------------------------------double D3 array */
if (fodim != 0) dserror("Illegal fourth dimension in call to am4def");
a->a.d3       = (double***)MALLOC((fdim*sizeof(double**)));
if (!(a->a.d3))       dserror("Allocation of memory failed");
a->a.d3[0]    = (double**) MALLOC((fdim*sdim*sizeof(double*)));
if (!(a->a.d3[0]))    dserror("Allocation of memory failed");
a->a.d3[0][0] = (double*)  MALLOC((fdim*sdim*tdim*sizeof(double)));
if (!(a->a.d3[0][0])) dserror("Allocation of memory failed");

for (i=1; i<fdim; i++)    a->a.d3[i]    = &(a->a.d3[0][i*sdim]);
endloop=fdim*sdim;
for (i=1; i<endloop; i++) a->a.d3[0][i] = &(a->a.d3[0][0][i*tdim]);
break;

case I3: /* ----------------------------------------------int I3 array */
if (fodim != 0) dserror("Illegal fourth dimension in call to am4def");
a->a.i3       = (int***)MALLOC((fdim*sizeof(int**)));
if (!(a->a.i3))       dserror("Allocation of memory failed");
a->a.i3[0]    = (int**) MALLOC((fdim*sdim*sizeof(int*)));
if (!(a->a.i3[0]))    dserror("Allocation of memory failed");
a->a.i3[0][0] = (int*)  MALLOC((fdim*sdim*tdim*sizeof(int)));
if (!(a->a.i3[0][0])) dserror("Allocation of memory failed");

for (i=1; i<fdim; i++)    a->a.i3[i]    = &(a->a.i3[0][i*sdim]);
endloop=fdim*sdim;
for (i=1; i<endloop; i++) a->a.i3[0][i] = &(a->a.i3[0][0][i*tdim]);
break;

case D4: /* --------------------------------------------double D4 array */
a->a.d4          = (double****)MALLOC((fdim*sizeof(double***)));
if (!(a->a.d4))          dserror("Allocation of memory failed");
a->a.d4[0]       = (double***) MALLOC((fdim*sdim*sizeof(double**)));
if (!(a->a.d4[0]))       dserror("Allocation of memory failed");
a->a.d4[0][0]    = (double**)  MALLOC((fdim*sdim*tdim*sizeof(double*)));
if (!(a->a.d4[0][0]))    dserror("Allocation of memory failed");
a->a.d4[0][0][0] = (double*)   MALLOC((fdim*sdim*tdim*fodim*sizeof(double)));
if (!(a->a.d4[0][0][0])) dserror("Allocation of memory failed");

for (i=1; i<fdim; i++)    a->a.d4[i]       = &(a->a.d4[0][i*sdim]);
endloop=fdim*sdim;
for (i=1; i<endloop; i++) a->a.d4[0][i]    = &(a->a.d4[0][0][i*tdim]);
endloop=fdim*sdim*tdim;
for (i=1; i<endloop; i++) a->a.d4[0][0][i] = &(a->a.d4[0][0][0][i*fodim]);
break;

case I4: /* ----------------------------------------------int I4 array */
a->a.i4          = (int****)MALLOC((fdim*sizeof(int***)));
if (!(a->a.i4))          dserror("Allocation of memory failed");
a->a.i4[0]       = (int***) MALLOC((fdim*sdim*sizeof(int**)));
if (!(a->a.i4[0]))       dserror("Allocation of memory failed");
a->a.i4[0][0]    = (int**)  MALLOC((fdim*sdim*tdim*sizeof(int*)));
if (!(a->a.i4[0][0]))    dserror("Allocation of memory failed");
a->a.i4[0][0][0] = (int*)   MALLOC((fdim*sdim*tdim*fodim*sizeof(int)));
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
#ifdef DEBUG 
/* ARRAY4D's cannot be traced yet m.gee
if (trace.trace_on==1)
{
dstracereport(a);
} 
*/
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
return((void*)(a->a.d3));
} /* end of am4def */

/*----------------------------------------------------------------------*
 | delete 4dimensional array                             m.gee 12/01    |
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
case D3:
   FREE(array->a.d3[0][0]);
   FREE(array->a.d3[0]);
   array->a.d3=FREE(array->a.d3);
break;
case I3:
   FREE(array->a.i3[0][0]);
   FREE(array->a.i3[0]);
   array->a.i3=FREE(array->a.i3);
break;
case D4:
   FREE(array->a.d4[0][0][0]);
   FREE(array->a.d4[0][0]);
   FREE(array->a.d4[0]);
   array->a.d4=FREE(array->a.d4);
break;
case I4:
   FREE(array->a.i4[0][0][0]);
   FREE(array->a.i4[0][0]);
   FREE(array->a.i4[0]);
   array->a.i4=FREE(array->a.i4);
break;
default:
dserror("Unknown type of array given");
}
/*----------------------------------------------------- delete the type */
array->Typ = XX4D;
/*------------------------- delete the array from the bugtracing system */
#ifdef DEBUG 
/* ARRAY4D cannot be traced yet m.gee 
if (trace.trace_on==1)
{
    trace.arrays[array->place_in_trace]=NULL;
    trace.num_arrays--;
    array->place_in_trace=0;
}                     */
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
return;
} /* end of am4del */

/*----------------------------------------------------------------------*
 | initialize an 4D array by zero                           m.gee 12/01 |
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
case D3:
   dim = (array->fdim) * (array->sdim) * (array->tdim);
   dptr = array->a.d3[0][0];
   for (i=0; i<dim; i++) *(dptr++) = 0.0;
break; 
case D4:
   dim = (array->fdim) * (array->sdim) * (array->tdim) * (array->fodim);
   dptr = array->a.d4[0][0][0];
   for (i=0; i<dim; i++) *(dptr++) = 0.0;
break;
case I3:
   dim = (array->fdim) * (array->sdim) * (array->tdim);
   iptr = array->a.i3[0][0];
   for (i=0; i<dim; i++) *(iptr++) = 0; 
break; 
case I4:
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

/*----------------------------------------------------------------------*
 | initialize an array by value                           m.gee 6/01    |
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
case D3:
   dim = (array->fdim) * (array->sdim) * (array->tdim);
   dvalue = (double*)value;
   dptr   = array->a.d3[0][0];
   for (i=0; i<dim; i++) *(dptr++) = *dvalue;
break; 
case D4:
   dim = (array->fdim) * (array->sdim) * (array->tdim) * (array->fodim);
   dvalue = (double*)value;
   dptr = array->a.d4[0][0][0];
   for (i=0; i<dim; i++) *(dptr++) = *dvalue;
break;
case I3:
   dim = (array->fdim) * (array->sdim) * (array->tdim);
   ivalue = (int*)value;
   iptr = array->a.i3[0][0];
   for (i=0; i<dim; i++) *(iptr++) = *ivalue;
break; 
case I4:
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


/*----------------------------------------------------------------------*
 | allocate and make a copy of ARRAY4D                     m.gee 12/01  |
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
case D3:
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
case D4:
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
case I3:
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
case I4:
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


/*----------------------------------------------------------------------*
 | make a copy of ARRAY4D                                  m.gee 12/01  |
 | the given arrays must be of equal type and size                      |
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
case D3:
   dim    = array_from->fdim * array_from->sdim * array_from->tdim;
   dimnew = array_to->fdim * array_to->sdim * array_to->tdim;
   if (dim != dimnew) dserror("mismatching dimensions, cannot copy ARRAY4D");
   dptr_from = array_from->a.d3[0][0];
   dptr_to   = array_to->a.d3[0][0];
   for (i=0; i<dim; i++) *(dptr_to++) = *(dptr_from++);
break; 
case D4:
   dim    = array_from->fdim * array_from->sdim * array_from->tdim * array_from->fodim;
   dimnew = array_to->fdim * array_to->sdim * array_to->tdim * array_to->fodim;
   if (dim != dimnew) dserror("mismatching dimensions, cannot copy ARRAY4D");
   dptr_from = array_from->a.d4[0][0][0];
   dptr_to   = array_to->a.d4[0][0][0];
   for (i=0; i<dim; i++) *(dptr_to++) = *(dptr_from++);
break;
case I3:
   dim    = array_from->fdim * array_from->sdim * array_from->tdim;
   dimnew = array_to->fdim * array_to->sdim * array_to->tdim;
   if (dim != dimnew) dserror("mismatching dimensions, cannot copy ARRAY4D");
   iptr_from = array_from->a.i3[0][0];
   iptr_to   = array_to->a.i3[0][0];
   for (i=0; i<dim; i++) *(iptr_to++) = *(iptr_from++);
break; 
case I4:
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




/*----------------------------------------------------------------------*
 | redefine of ARRAY4D                                     m.gee 12/01  |
 | NOTE:  - a type cast like it is possible with ARRAY doesn't work here|
 |        - 4-dimensional arrays with fodim=0 are not allowed           |
 |        - remember, this routine can be VERY expensive !              |
 |        - aditional space in the redefined array is set to zero       |
 |          (unlike in am4def which does NOT initialize)                |
 *----------------------------------------------------------------------*/
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
case D3:
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
case D4:
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
case I3:
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
case I4:
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




