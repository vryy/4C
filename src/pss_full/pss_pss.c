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

#include "../headers/standardtypes.h"
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 | write a record to pss-file                            m.gee 11/00    |
 |                                                                      |
 | ierr=0     error in writing                                          |
 | ierr=1     writing successfull                                       |
 | ierr=2     record with this name already exists                      |
 |                                                                      |
 | name          name of record                                     (I) |
 | fdim          first dimension of record                          (I) |
 | sdim          scnd dimension of record                           (I) |
 | byte          size of one object inside record                   (I) |
 | startaddress  starting address of record to write                (I) |
 | ierr          error handler (see above)                          (O) |
 |                                                                      |
 | Attention: a record to write has to be contigous in memory !!        |
 |                                                                      |
 | How does a record on pss-file look like:                             |
 |                                                                      |
 | INT           number_chars_of_name                                   |
 | char          name[]                                                 |
 | INT           handle                                                 |
 | INT           fdim                                                   |
 | INT           sdim                                                   |
 | INT           byte                                                   |
 | void          record[]                                               |
 |                                                                      |
 |                                                                      |
 *----------------------------------------------------------------------*/
void pss_write(char          *name, 
               INT            fdim, 
               INT            sdim,
               INT            byte,
               const void    *startaddress,
               long int      *handle,
               FILE          *out, 
               INT           *ierr)
{
INT          name_size=0;
INT          write_error=0;
INT          dimensions[3];
#ifdef DEBUG 
dstrc_enter("pss_write");
#endif

*ierr=0;
/*----------------------------------------------------- calc name_size */
name_size=strlen(name);

/*----------------------------------------------------- get the handle */
*handle = ftell(out);

/*---------------------------------------------------- right name_size */
write_error = fwrite(&name_size,sizeof(INT),1,out);
if (write_error!=1) dserror("Error writing pss-file");

/*----------------------------------------------- write name of record */
write_error = fwrite(name,sizeof(char),name_size,out);
if (write_error!=name_size) dserror("Error writing pss-file");

/*--------------------------------------------------------write handle */
write_error = fwrite(handle,sizeof(long int),1,out);
if (write_error!=1) dserror("Error writing pss-file");

/*--------------------------------------------------- write dimensions */
dimensions[0]=fdim;
dimensions[1]=sdim;
dimensions[2]=byte;

write_error = fwrite(dimensions,sizeof(INT),3,out);
if (write_error!=3) dserror("Error writing pss-file");

/*--------------------------------------------------- write the record */
write_error = fwrite(startaddress,byte,fdim*sdim,out);
if (write_error!=fdim*sdim) dserror("Error writing pss-file");
else *ierr=1;

fflush(out);

#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of pss_write */




/*----------------------------------------------------------------------*
 | write the content of ARRAY to pss-file                m.gee 11/00    |
 |                                                                      |
 | ierr=0     error in writing                                          |
 | ierr=1     writing successfull                                       |
 | ierr=2     record with this name already exists                      |
 |                                                                      |
 | const ARRAY *array (input) adress of array to be written             |
 | INT         *handle(output) unique handle returned by the pss-system |
 | INT         *ierr  (output) success flag                             |
 |                                                                      |
 *----------------------------------------------------------------------*/
void pss_write_array(const ARRAY *array, 
                     long int    *handle,
                     FILE        *out, 
                     INT         *ierr)
{
INT           name_size=0;
INT           write_error=0;
INT           dimensions[3];
#ifdef DEBUG 
dstrc_enter("pss_write_array");
#endif

*ierr=0;
/*----------------------------------------------------- calc name_size */
name_size=strlen(array->name);

if (name_size>9) name_size=9;
/*----------------------------------------------------- get the handle */
*handle = ftell(out);

/*---------------------------------------------------- right name_size */
write_error = fwrite(&name_size,sizeof(INT),1,out);
if (write_error!=1) dserror("Error writing pss-file");

/*----------------------------------------------- write name of record */
write_error = fwrite(array->name,sizeof(char),name_size,out);
if (write_error!=name_size) dserror("Error writing pss-file");

/*------------------------------------------- make handle and write it */
write_error = fwrite(handle,sizeof(long int),1,out);
if (write_error!=1) dserror("Error writing pss-file");

/*--------------------------------------------------- write dimensions */
dimensions[0]=array->fdim;
dimensions[1]=array->sdim;
switch (array->Typ)
{
  case cca_DA:
    dimensions[2]=sizeof(DOUBLE);
    break;
  case cca_DV:
    dimensions[2]=sizeof(DOUBLE);
    break;
  case cca_IA:
    dimensions[2]=sizeof(INT);
    break;
  case cca_IV:
    dimensions[2]=sizeof(INT);
    break;
  default:
    dserror("Wrong array type!");
    break;
}
write_error = fwrite(dimensions,sizeof(INT),3,out);
if (write_error!=3) dserror("Error writing pss-file");

/*--------------------------------------------------- write the record */
switch (array->Typ)
{
   case cca_DA:
      write_error = fwrite(
                           array->a.da[0],
                           sizeof(DOUBLE),
                           (dimensions[0]*dimensions[1]),
                           out
                          );
      if (write_error != dimensions[0]*dimensions[1]) 
      dserror("Error writing pss-file");
      else *ierr=1;
   break;
   case cca_DV:
      write_error = fwrite(
                           array->a.dv,
                           sizeof(DOUBLE),
                           (dimensions[0]*dimensions[1]),
                           out
                          );
      if (write_error != dimensions[0]*dimensions[1]) 
      dserror("Error writing pss-file");
      else *ierr=1;
   break;
   case cca_IA:
      write_error = fwrite(
                           array->a.ia[0],
                           sizeof(INT),
                           (dimensions[0]*dimensions[1]),
                           out
                          );
      if (write_error != dimensions[0]*dimensions[1]) 
      dserror("Error writing pss-file");
      else *ierr=1;
   break;
   case cca_IV:
      write_error = fwrite(
                           array->a.iv,
                           sizeof(INT),
                           (dimensions[0]*dimensions[1]),
                           out
                          );
      if (write_error != dimensions[0]*dimensions[1]) 
      dserror("Error writing pss-file");
      else *ierr=1;
   break;
  default:
    dserror("Wrong array type!");
    break;
}
fflush(out);

#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of pss_write_array */









/*----------------------------------------------------------------------*
 | read a record from pss-file                           m.gee 06/01    |
 |                                                                      |
 | reads the first record in pss with name NAME, doesn't care for       |
 | the handle, but returns it                                           |                       
 | ierr=0     error in reading                                          |
 | ierr=1     reading successfull                                       |
 | ierr=2     record with this name doesn't exists                      |
 |                                                                      |
 | char *name (input) name of the record to read                        |
 | INT  *fdim (output) first dimension of record                        |
 | INT  *sdim (output) scnd dimension of record                         |
 | INT  *byte (output) size in byte of one entry in record              |
 | 
 |                                                                      |
 *----------------------------------------------------------------------*/
void pss_read_name(char      *name, 
                   INT       *fdim, 
                   INT       *sdim,
                   INT       *byte,
                   void      *ziel,
                   long int  *handle, 
                   FILE      *in,
                   INT       *ierr)
{
long int      cur_pos = 0;
long int      offset = 0;
INT           err=0;
INT           foundit=0;
INT           name_size=0;
INT           sizename=0;
INT           dimensions[3];
char          test_name[200];
long int      handle_dummy;
#ifdef DEBUG 
dstrc_enter("pss_read_name");
#endif

*ierr=0;
/*------------------------------- save current position of file-pointer */
cur_pos = ftell(in);

/*----------------------------------------------------- calc name_size */
name_size=strlen(name);

/*--------------------------------------------------------- rewind file */
rewind(in);

/*------------------------------- loop all records till record is found */
do
{
/*--------------------------------------------------- read size of name */
err=fread(&sizename,sizeof(INT),1,in);
/*----------------------------------------------- check for end of file */
if ( feof(in) != 0)
{
   foundit=0;
   err = fseek(in,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=2;
   goto end;    
}
/*------------------------------------------------------ read test_name */ 
err=fread(test_name,sizeof(char),sizename,in);
if (err != sizename) dserror("error reading pss-file");

/*------------------------------------ check the test_name against name */
if ( strncmp(name,test_name,strlen(name)) != 0 )
{
/*-------------------- not the right record, so jump it and start again */
   err=fread(&handle_dummy,sizeof(long int),1,in);
   if (err != 1) dserror("error reading pss-file");
   err=fread(dimensions,sizeof(INT),3,in);
   if (err != 3) dserror("error reading pss-file");
   offset = dimensions[0]*dimensions[1]*dimensions[2];
   err = fseek(in,offset,SEEK_CUR);
   if (err == -1) dserror("error reading pss-file");
}
else
{
   foundit=1;
/*--------------------------------------------------------- read handle */
   err=fread(handle,sizeof(long int),1,in);
   if (err != 1) dserror("error reading pss-file");
/*--------------------------------------- read the dimensions of record */
   err=fread(dimensions,sizeof(INT),3,in);
   if (err != 3) dserror("error reading pss-file");
   *fdim=dimensions[0];
   *sdim=dimensions[1];
   *byte=dimensions[2];
/*--------------------------------------------------------- read record */
   err=fread(ziel,(*byte),(*fdim)*(*sdim),in);
   if (err != (*fdim)*(*sdim)) dserror("error reading pss-file");
/*---------------------------------- reset the position of file-pointer */
   err = fseek(in,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=1;    
}
} while (foundit != 1);

end:
err = fseek(in,cur_pos,SEEK_SET);
if (err == -1) dserror("error reading pss-file");
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of pss_read_name */





/*----------------------------------------------------------------------*
 | read a record from pss-file                           m.gee 06/01    |
 |                                                                      |
 | reads the first record in pss with name NAME and correct handle      |
 | ierr=0     error in reading                                          |
 | ierr=1     reading successfull                                       |
 | ierr=2     record with this name doesn't exists                      |
 |                                                                      |
 |                                                                      |
 *----------------------------------------------------------------------*/
void pss_read_name_handle(char       *name, 
                          INT	     *fdim, 
                          INT	     *sdim,
                          INT	     *byte,
                          void       *ziel, 
                          long int   *handle,
                          FILE       *in, 
                          INT	     *ierr)
{
long int      cur_pos=0;
INT           err=0;
INT           foundit=0;
INT           name_size=0;
INT           sizename=0;
INT           dimensions[3];
long int      handle_dummy;
char          test_name[200];
#ifdef DEBUG 
dstrc_enter("pss_read_name_handle");
#endif

*ierr=0;
/*------------------------------- save current position of file-pointer */
cur_pos = ftell(in);

/*------------------------------------------------------ calc name_size */
name_size=strlen(name);

/*--------------------------------- set file pointer to handle position */
err = fseek(in,*handle,SEEK_SET);
if (err == -1) dserror("error reading pss-file");

/*--------------------------------------------------- read size of name */
err = fread(&sizename,sizeof(INT),1,in);
if (err != 1) dserror("error reading pss-file");

/*------------------------------------------------------ read test_name */ 
err=fread(test_name,sizeof(char),sizename,in);
if (err != sizename) dserror("error reading pss-file");

/*--------------------------------------------------- read handle_dummy */
err=fread(&handle_dummy,sizeof(long int),1,in);
if (err != 1) dserror("error reading pss-file");

/*---- check the test_name against name and handle_dummy against handle */
if ( strncmp(name,test_name,strlen(name)) != 0 || handle_dummy != *handle)
{
/*------------------------- not the right record, so return with ierr=2 */
   *ierr=2;
   goto end;
}
else
{
   foundit=1;
/*--------------------------------------- read the dimensions of record */
   err=fread(dimensions,sizeof(INT),3,in);
   if (err != 3) dserror("error reading pss-file");
   *fdim=dimensions[0];
   *sdim=dimensions[1];
   *byte=dimensions[2];
/*--------------------------------------------------------- read record */
   err=fread(ziel,(*byte),(*fdim)*(*sdim),in);
   if (err != (*fdim)*(*sdim)) dserror("error reading pss-file");
/*---------------------------------- reset the position of file-pointer */
   err = fseek(in,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=1;    
}

end:
err = fseek(in,cur_pos,SEEK_SET);
if (err == -1) dserror("error reading pss-file");
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of pss_read_name_handle */








/*----------------------------------------------------------------------*
 | read an array from pss-file                           m.gee 11/00    |
 |                                                                      |
 | routine reads the first ARRAY with name NAME, doesn't care for       |
 | the handle but returns it                                            |
 | ierr=0     error in reading                                          |
 | ierr=1     reading successfull                                       |
 | ierr=2     record with this name doesn't exists                      |
 |                                                                      |
 |  NOTE: the structure for the array must be provided by user,         |
 |        space must be allocated in correct amount to hold record      |
 |        before calling this routine                                   |
 |                                                                      |
 *----------------------------------------------------------------------*/
void pss_read_array_name(char       *name, 
                         ARRAY      *array,
                         long int   *handle,
                         FILE       *in,
                         INT        *ierr)
{
long int     cur_pos=0;
long int     offset=0;
INT          err=0;
INT          foundit=0;
INT          name_size=0; 
INT          sizename=0;
INT          dimensions[3];
long int     handle_dummy;
char         test_name[200];
#ifdef DEBUG 
dstrc_enter("pss_read_array_name");
#endif

*ierr=0;
/*------------------------------- save current position of file-pointer */
cur_pos = ftell(in);

/*----------------------------------------------------- calc name_size */
name_size=strlen(name);

/*--------------------------------------------------------- rewind file */
rewind(in);

/*------------------------------- loop all records till record is found */
do
{
/*--------------------------------------------------- read size of name */
err=fread(&sizename,sizeof(INT),1,in);
/*----------------------------------------------- check for end of file */
if ( feof(in) != 0)
{
   foundit=0;
   err = fseek(in,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=2;
   goto end;    
}
/*------------------------------------------------------ read test_name */ 
err=fread(test_name,sizeof(char),sizename,in);
if (err != sizename) dserror("error reading pss-file");
/*--------------------------------------------------------- read handle */
err=fread(&handle_dummy,sizeof(long int),1,in);
if (err != 1) dserror("error reading pss-file");
/*------------------------------------ check the test_name against name */
if ( strncmp(name,test_name,strlen(name)) != 0 )
{
/*-------------------- not the right record, so jump it and start again */
   err=fread(dimensions,sizeof(INT),3,in);
   if (err != 3) dserror("error reading pss-file");
   offset = dimensions[0]*dimensions[1]*dimensions[2];
   err = fseek(in,offset,SEEK_CUR);
   if (err == -1) dserror("error reading pss-file");
}
else
{
   foundit=1;
/*--------------------------------------------------------- read handle */
   *handle = handle_dummy;
/*--------------------------------------- read the dimensions of record */
   err=fread(dimensions,sizeof(INT),3,in);
   if (err != 3) dserror("error reading pss-file");
   if ((unsigned)dimensions[2] == sizeof(INT)) /* it's a INT-record */
   {
       if ( array->Typ == cca_IV )/* it's an INT-vector */
       {
          err=fread(&(array->a.iv[0]),dimensions[2],(dimensions[0]*dimensions[1]),in);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }
       else /* it's an INT-array */
       {
          err=fread(&(array->a.ia[0][0]),dimensions[2],(dimensions[0]*dimensions[1]),in);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }
   }
   else if ((unsigned)dimensions[2] == sizeof(DOUBLE))/* it's a DOUBLE-record */
   {
       if ( array->Typ == cca_DV )/* it's a DOUBLE-vector */
       {
          err=fread(&(array->a.dv[0]),dimensions[2],(dimensions[0]*dimensions[1]),in);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }
       else /* it's a DOUBLE-array */
       {
          err=fread(&(array->a.da[0][0]),dimensions[2],(dimensions[0]*dimensions[1]),in);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }  
   }
   else /* the size in byte does not fit to INT or DOUBLE, so it's unknown */
   {
      dserror("error reading array from pss-file: unexpected byte size");
   }

/*---------------------------------- reset the position of file-pointer */
   err = fseek(in,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=1;    
}
} while (foundit != 1);


end:
err = fseek(in,cur_pos,SEEK_SET);
if (err == -1) dserror("error reading pss-file");
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of pss_read_array_name */





/*----------------------------------------------------------------------*
 | read an array from pss-file                           m.gee 11/00    |
 |                                                                      |
 | routine reads the first ARRAY with name NAME and handle HANDLE       |
 | ierr=0     error in reading                                          |
 | ierr=1     reading successfull                                       |
 | ierr=2     record with this name doesn't exists                      |
 |                                                                      |
 |  NOTE: the structure for the array must be provided by user,         |
 |        space must be allocated in correct amount to hold record      |
 |        before calling this routine                                   |
 |                                                                      |
 *----------------------------------------------------------------------*/
void pss_read_array_name_handle(char       *name, 
                                ARRAY	   *array,
                                long int   *handle,
                                FILE       *in,
                                INT	   *ierr)
{
long int    cur_pos = 0;
INT         err=0; 
INT         foundit=0;
unsigned int name_size=0; 
INT         sizename=0;
INT         dimensions[3];
long int    handle_dummy;
char        test_name[200];
#ifdef DEBUG 
dstrc_enter("pss_read_array_name_handle");
#endif

*ierr=0;
/*------------------------------- save current position of file-pointer */
cur_pos = ftell(in);

/*------------------------------------------------------ calc name_size */
name_size=strlen(name);

/*----------------------------------------- set file to handle position */
err = fseek(in,*handle,SEEK_SET);
if (err == -1) dserror("error reading pss-file");

/*--------------------------------------------------- read size of name */
err=fread(&sizename,sizeof(INT),1,in);
/*------------------------------------------------------ read test_name */ 
err=fread(test_name,sizeof(char),sizename,in);
if (err != sizename) dserror("error reading pss-file");
/*--------------------------------------------------------- read handle */
err=fread(&handle_dummy,sizeof(long int),1,in);
if (err != 1) dserror("error reading pss-file");
/*---- check the test_name against name and handle_dummy against handle */
if ( strncmp(name,test_name,strlen(name)) != 0 || handle_dummy != *handle)
{
/*-------------------------- not the right record, so return with error */
   *ierr=2;
   goto ende;
}
else
{
   foundit=1;
/*--------------------------------------- read the dimensions of record */
   err=fread(dimensions,sizeof(INT),3,in);
   if (err != 3) dserror("error reading pss-file");
   if ((unsigned)dimensions[2] == sizeof(INT)) /* it's a INT-record */
   {
       if ( array->Typ == cca_IV )/* it's an INT-vector */
       {
          err=fread(&(array->a.iv[0]),dimensions[2],(dimensions[0]*dimensions[1]),in);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }
       else /* it's an INT-array */
       {
          err=fread(&(array->a.ia[0][0]),dimensions[2],(dimensions[0]*dimensions[1]),in);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }
   }
   else if ((unsigned)dimensions[2] == sizeof(DOUBLE))/* it's a DOUBLE-record */
   {
       if ( array->Typ == cca_DV )/* it's a DOUBLE-vector */
       {
          err=fread(&(array->a.dv[0]),dimensions[2],(dimensions[0]*dimensions[1]),in);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }
       else /* it's a DOUBLE-array */
       {
          err=fread(&(array->a.da[0][0]),dimensions[2],(dimensions[0]*dimensions[1]),in);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }  
   }
   else /* the size in byte does not fit to INT or DOUBLE, so it's unknown */
   {
      dserror("error reading array from pss-file: unexpected byte size");
   }

/*---------------------------------- reset the position of file-pointer */
   err = fseek(in,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=1;    
}


ende:
err = fseek(in,cur_pos,SEEK_SET);
if (err == -1) dserror("error reading pss-file");
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of pss_read_array_name_handle */




/*----------------------------------------------------------------------*
 | read an array from pss-file                           m.gee 11/00    |
 |                                                                      |
 | routine reads the first ARRAY with name NAME and handle HANDLE       |
 | ierr=0     error in reading                                          |
 | ierr=1     reading successfull                                       |
 | ierr=2     record with this name doesn't exists                      |
 |                                                                      |
 |  NOTE: the structure for the array must be provided by user,         |
 |        space must be allocated in correct amount to hold record      |
 |        before calling this routine                                   |
 |                                                                      |
 *----------------------------------------------------------------------*/
void pss_read_array_handle(ARRAY      *array,
                           long int   *handle,
                           FILE       *in,
                           INT        *ierr)
{
long int    cur_pos = 0;
INT         err=0; 
INT         foundit=0;
unsigned int name_size=0; 
INT         sizename=0;
INT         dimensions[3];
long int    handle_dummy;
char        test_name[200];
char        *name = "         "; 
#ifdef DEBUG 
dstrc_enter("pss_read_array_handle");
#endif

*ierr=0;
/*------------------------------- save current position of file-pointer */
cur_pos = ftell(in);

/*------------------------------------------------------ calc name_size */
name_size=strlen(name);

/*----------------------------------------- set file to handle position */
err = fseek(in,*handle,SEEK_SET);
/*--------------------------------------------------- read size of name */
err=fread(&sizename,sizeof(INT),1,in);
/*------------------------------------------------------ read test_name */ 
err=fread(test_name,sizeof(char),sizename,in);
if (err != sizename) dserror("error reading pss-file");
/*--------------------------------------------------------- read handle */
err=fread(&handle_dummy,sizeof(long int),1,in);
if (err != 1) dserror("error reading pss-file");
/*---------------------------------- check  handle_dummy against handle */
if (handle_dummy != *handle)
{
/*-------------------- not the right record, so jump it and start again */
   *ierr=2;
   goto end;
}
else
{
   foundit=1;
/*------------------------------------------ copy the name of the array */
   sizename = IMIN(sizename,9);
   strncpy(name,test_name,sizename);
/*--------------------------------------- read the dimensions of record */
   err=fread(dimensions,sizeof(INT),3,in);
   if (err != 3) dserror("error reading pss-file");
   if ((unsigned)dimensions[2] == sizeof(INT)) /* it's a INT-record */
   {
       if ( array->Typ == cca_IV )/* it's an INT-vector */
       {
          err=fread(&(array->a.iv[0]),dimensions[2],(dimensions[0]*dimensions[1]),in);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }
       else /* it's an INT-array */
       {
          err=fread(&(array->a.ia[0][0]),dimensions[2],(dimensions[0]*dimensions[1]),in);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }
   }
   else if ((unsigned)dimensions[2] == sizeof(DOUBLE))/* it's a DOUBLE-record */
   {
       if ( array->Typ == cca_DV )/* it's a DOUBLE-vector */
       {
          err=fread(&(array->a.dv[0]),dimensions[2],(dimensions[0]*dimensions[1]),in);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }
       else /* it's a DOUBLE-array */
       {
          err=fread(&(array->a.da[0][0]),dimensions[2],(dimensions[0]*dimensions[1]),in);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }  
   }
   else /* the size in byte does not fit to INT or DOUBLE, so it's unknown */
   {
      dserror("error reading array from pss-file: unexpected byte size");
   }
   *ierr=1;    
}


end:
err = fseek(in,cur_pos,SEEK_SET);
if (err == -1) dserror("error reading pss-file");
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of pss_read_array_handle */







/*----------------------------------------------------------------------*
 | check  whether record exists in pss-file              m.gee 11/00    |
 |                                                                      |
 | ierr=0     record doesn't exist                                      |
 | ierr=1     record exists                                             |
 |                                                                      |
 *----------------------------------------------------------------------*/
void pss_chck(char       *name,
              long int   *handle, 
              FILE       *in,
              INT        *ierr)
{
long int    cur_pos=0;
long int    offset=0;
INT         foundit=0;
INT         sizename=0;
INT         err=0;
char        test_name[200];
INT         dimensions[3];
long int    handle_dummy;
#ifdef DEBUG 
dstrc_enter("pss_chck");
#endif

/*------------------------------- save current position of file-pointer */
cur_pos = ftell(in);

/*--------------------------------------------------------- rewind file */
rewind(in);

/*------------------------------- loop all records till record is found */
do
{
/*--------------------------------------------------- read size of name */
err=fread(&sizename,sizeof(INT),1,in);
/*----------------------------------------------- check for end of file */
if ( feof(in) != 0)
{
   foundit=0;
   err = fseek(in,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=0;
   *handle=-1;
   goto end;    
}
/*------------------------------------------------------ read test_name */ 
err=fread(test_name,sizeof(char),sizename,in);
if (err != sizename) dserror("error reading pss-file");

/*------------------------------------ check the test_name against name */
if ( strncmp(name,test_name,strlen(name)) != 0 )
{
/*-------------------- not the right record, so jump it and start again */
   err = fread(&handle_dummy,sizeof(long int),1,in);
   if (err != 1) dserror("error reading pss-file");
   err=fread(dimensions,sizeof(INT),3,in);
   if (err != 3) dserror("error reading pss-file");
   offset = dimensions[0]*dimensions[1]*dimensions[2];
   err = fseek(in,offset,SEEK_CUR);
   if (err == -1) dserror("error reading pss-file");
}
else
{
   foundit=1;
/*----------------------------------------------------- read the handle */   
   err = fread(handle,sizeof(long int),1,in);
   if (err != 1) dserror("error reading pss-file");
/*---------------------------------- reset the position of file-pointer */
   err = fseek(in,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=1;    
}
} while (foundit != 1);
end:
err = fseek(in,cur_pos,SEEK_SET);
if (err == -1) dserror("error reading pss-file");
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of pss_chck */




/*----------------------------------------------------------------------*
 | check  whether record and handle exists in pss-file   m.gee 06/01    |
 |                                                                      |
 | ierr=0     record doesn't exist                                      |
 | ierr=1     record exists                                             |
 |                                                                      |
 *----------------------------------------------------------------------*/
void pss_chck_handle(char       *name,
                     long int   *handle, 
                     FILE       *in,
                     INT        *ierr)
{
long int    cur_pos=0;
INT         foundit=0;
INT         sizename=0;
INT         err=0;
char        test_name[200];
long int    handle_dummy;
#ifdef DEBUG 
dstrc_enter("pss_chck_handle");
#endif

/*------------------------------- save current position of file-pointer */
cur_pos = ftell(in);

/*----------------------------------------- set file to handle position */
err = fseek(in,*handle,SEEK_SET);

/*--------------------------------------------------- read size of name */
err=fread(&sizename,sizeof(INT),1,in);
/*------------------------------------------------------ read test_name */ 
err=fread(test_name,sizeof(char),sizename,in);
if (err != sizename) dserror("error reading pss-file");
/*----------------------------------------------------- read the handle */   
err = fread(&handle_dummy,sizeof(long int),1,in);
if (err != 1) dserror("error reading pss-file");
/*---- check the test_name against name and handle_dummy against handle */
if ( strncmp(name,test_name,strlen(name)) != 0 && handle_dummy != *handle)
{
/*-------------------- not the right record, so jump it and start again */
   *ierr=0;
   goto end;
}
else
{
   foundit=1;
/*---------------------------------- reset the position of file-pointer */
   err = fseek(in,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=1;    
}

end:
err = fseek(in,cur_pos,SEEK_SET);
if (err == -1) dserror("error reading pss-file");
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of pss_chck_handle */





/*----------------------------------------------------------------------*
 | returns the dimensions of record name                 m.gee 11/00    |
 |                                                                      |
 | takes first record with name NAME and neglects the handle but        |                            
 | returns it                                                           |
 *----------------------------------------------------------------------*/
void pss_getdims_name(char       *name, 
                      INT	 *fdim,
                      INT	 *sdim,
                      INT	 *byte,
                      long int	 *handle,
                      FILE       *in,
                      INT	 *ierr)
{
long int      cur_pos=0;
long int      offset=0;
INT           foundit=0;
INT           sizename=0;
INT           err=0;
char          test_name[200];
INT           dimensions[3];
long int      handle_dummy;
#ifdef DEBUG 
dstrc_enter("pss_getdims_name");
#endif

/*------------------------------- save current position of file-pointer */
cur_pos = ftell(in);

/*--------------------------------------------------------- rewind file */
rewind(in);
*ierr=0;
/*------------------------------- loop all records till record is found */
do
{
/*--------------------------------------------------- read size of name */
err=fread(&sizename,sizeof(INT),1,in);
/*----------------------------------------------- check for end of file */
if ( feof(in) != 0)
{
   foundit=0;
   err = fseek(in,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=1;
   goto end;    
}
/*------------------------------------------------------ read test_name */ 
err=fread(test_name,sizeof(char),sizename,in);
if (err != sizename) dserror("error reading pss-file");
/*----------------------------------------------- read the handle_dummy */
err=fread(&handle_dummy,sizeof(long int),1,in);
if (err != 1) dserror("error reading pss-file");
/*------------------------------------ check the test_name against name */
if ( strncmp(name,test_name,strlen(name)) != 0 )
{
/*-------------------- not the right record, so jump it and start again */
   err=fread(dimensions,sizeof(INT),3,in);
   if (err != 3) dserror("error reading pss-file");
   offset = dimensions[0]*dimensions[1]*dimensions[2];
   err = fseek(in,offset,SEEK_CUR);
   if (err == -1) dserror("error reading pss-file");
}
else
{
   foundit=1;
   err=fread(dimensions,sizeof(INT),3,in);
   if (err != 3) dserror("error reading pss-file");
   *fdim=dimensions[0];
   *sdim=dimensions[1];
   *byte=dimensions[2];
   *handle = handle_dummy;
/*---------------------------------- reset the position of file-pointer */
   err = fseek(in,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=1;    
}
} while (foundit != 1);

end:
err = fseek(in,cur_pos,SEEK_SET);
if (err == -1) dserror("error reading pss-file");
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of pss_getdims_name */




/*----------------------------------------------------------------------*
 | returns the dimensions of record name                 m.gee 06/01    |
 |                                                                      |
 | takes first record with name NAME and handle HANDLE                  |                            
 |                                                                      |
 *----------------------------------------------------------------------*/
void pss_getdims_name_handle(char       *name, 
                             INT	*fdim,
                             INT	*sdim,
                             INT	*byte,
                             long int	*handle,
                             FILE       *in,
                             INT	*ierr)
{
long int      cur_pos=0;
INT           foundit=0;
INT           sizename=0;
INT           err=0;
char          test_name[200];
INT           dimensions[3];
long int           handle_dummy;
#ifdef DEBUG 
dstrc_enter("pss_getdims_name_handle");
#endif

/*------------------------------- save current position of file-pointer */
*ierr=0;
cur_pos = ftell(in);
/*----------------------------------------- set file to handle position */
err = fseek(in,*handle,SEEK_SET);
if (err == -1) dserror("error reading pss-file");
/*--------------------------------------------------- read size of name */
err=fread(&sizename,sizeof(INT),1,in);
/*------------------------------------------------------ read test_name */ 
err=fread(test_name,sizeof(char),sizename,in);
if (err != sizename) dserror("error reading pss-file");
/*----------------------------------------------- read the handle_dummy */
err=fread(&handle_dummy,sizeof(long int),1,in);
if (err != 1) dserror("error reading pss-file");
/*---- check the test_name against name and handle_dummy against handle */
if ( strncmp(name,test_name,strlen(name)) != 0 || handle_dummy != *handle)
{
/*-------------------------- not the right record, so return with error */
   *ierr=2;
   goto end;
}
else
{
   foundit=1;
   err=fread(dimensions,sizeof(INT),3,in);
   if (err != 3) dserror("error reading pss-file");
   *fdim=dimensions[0];
   *sdim=dimensions[1];
   *byte=dimensions[2];
/*---------------------------------- reset the position of file-pointer */
   err = fseek(in,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=1;    
}

end:
err = fseek(in,cur_pos,SEEK_SET);
if (err == -1) dserror("error reading pss-file");
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of pss_getdims_name_handle */






/*----------------------------------------------------------------------*
 | prints status report about all records                m.gee 11/00    |
 | on the pss-file to err-file                                          |
 |                                                                      |
 |                                                                      |
 *----------------------------------------------------------------------*/
void pss_status_to_err(FILE *inout)
{
long int cur_pos=0;
long int current=0;
long int endoffile=0;
long int offset=0;
INT      sizename=0;
INT      err=0;
char     test_name[200];
char     writename[20];
INT      dimensions[3];
INT      counter=0;
long int handle;
#ifdef DEBUG 
dstrc_enter("pss_status_to_err");
#endif

fprintf(allfiles.out_err,"===========================================\n");
fprintf(allfiles.out_err,"pss-status - record report about pss-file\n");
fprintf(allfiles.out_err,"===========================================\n");
/*------------------------------- save current position of file-pointer */
cur_pos = ftell(inout);
fseek(inout,0,SEEK_END);
endoffile = ftell(inout);
/*--------------------------------------------------------- rewind file */
rewind(inout);

/*------------------------------------ loop all records till end of file*/
do
{
/*--------------------------------------------------- read size of name */
   current=ftell(inout);
   if (current==endoffile) goto end;
   err=fread(&sizename,sizeof(INT),1,inout);
   if (err != 1) goto end;
/*------------------------------------------------------ read test_name */ 
   err=fread(test_name,sizeof(char),sizename,inout);
   if (err != sizename) goto end;
/*--------------------------------------------------------- read handle */   
   err=fread(&handle,sizeof(long int),1,inout);
   if (err != 1) goto end;
/*----------------------------------------------------- read dimensions */
   err=fread(dimensions,sizeof(INT),3,inout);
   if (err != 3) goto end;
/*----------------------------------------------------- jump the record */
   offset = dimensions[0]*dimensions[1]*dimensions[2];
   err = fseek(inout,offset,SEEK_CUR);
   if (feof(inout) == 1) goto end;
   if (err == -1) goto end;
   counter++;
/*-------------------------------------- write report about this record */
   strncpy(writename,"                    ",20);
   strncpy(writename,test_name,sizename);
   fprintf(allfiles.out_err,"RECORD No %5d: DIMENSIONS: %6d x %6d x %6d BYTE NAME= %s HANDLE=%d\n",
        counter,
        dimensions[0],
        dimensions[1],
        dimensions[2],
        writename,
        handle);
   strncpy(writename,"                    ",20);

} while (1);
end:
fprintf(allfiles.out_err,"END OF FILE AT: %d BYTE\n",endoffile);
fprintf(allfiles.out_err,"BYTES READABLE: %d BYTE\n",current);
if ((current-endoffile)==0)
fprintf(allfiles.out_err,"PSS-FILE O.K. AND COMPLETE\n");
else
{
fprintf(allfiles.out_err,"WARNING: PSS-FILE MAY BE DAMAGED\n");
fprintf(allfiles.out_err,"WARNING: ERRORS OCCURED WHILE READING\n");
}
fprintf(allfiles.out_err,"===========================================\n");
fprintf(allfiles.out_err,"pss-status - number of records: %d END\n",counter);
fprintf(allfiles.out_err,"===========================================\n");

/*---------------------------------- reset the position of file-pointer */
err = fseek(inout,cur_pos,SEEK_SET);
if (err == -1) dserror("error reading pss-file");
fflush(allfiles.out_err);

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of pss_status_to_err */



