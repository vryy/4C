#include "../headers/standardtypes.h"

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
 | int           number_chars_of_name                                   |
 | char          name[]                                                 |
 | int           handle                                                 |
 | int           fdim                                                   |
 | int           sdim                                                   |
 | int           byte                                                   |
 | void          record[]                                               |
 |                                                                      |
 |                                                                      |
 *----------------------------------------------------------------------*/
void pss_write(char          *name, 
                  int            fdim, 
                  int            sdim,
                  int            byte,
                  const void    *startaddress,
                  int           *handle, 
                  int           *ierr)
{
int          i;
int          err=0;
int          name_size=0;
int          write_error=0;
int          dimensions[3];
#ifdef DEBUG 
dstrc_enter("pss_write");
#endif

*ierr=0;
/*----------------------------------------------------- calc name_size */
name_size=strlen(name);

/*---------------------------------------------------- right name_size */
write_error = fwrite(&name_size,sizeof(int),1,allfiles.out_pss);
if (write_error!=1) dserror("Error writing pss-file");

/*----------------------------------------------- write name of record */
write_error = fwrite(name,sizeof(char),name_size,allfiles.out_pss);
if (write_error!=name_size) dserror("Error writing pss-file");

/*------------------------------------------- make handle and write it */
*handle = allfiles.pss_counter;
(allfiles.pss_counter)++;
write_error = fwrite(handle,sizeof(int),1,allfiles.out_pss);
if (write_error!=1) dserror("Error writing pss-file");

/*--------------------------------------------------- write dimensions */
dimensions[0]=fdim;
dimensions[1]=sdim;
dimensions[2]=byte;

write_error = fwrite(dimensions,sizeof(int),3,allfiles.out_pss);
if (write_error!=3) dserror("Error writing pss-file");

/*--------------------------------------------------- write the record */
write_error = fwrite(startaddress,byte,fdim*sdim,allfiles.out_pss);
if (write_error!=fdim*sdim) dserror("Error writing pss-file");
else *ierr=1;

fflush(allfiles.out_pss);

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
 |                                                                      |
 *----------------------------------------------------------------------*/
void pss_write_array(const ARRAY *array, 
                        int         *handle, 
                        int         *ierr)
{
int           i;
int           err=0;
int           name_size=0;
int           write_error=0;
int           dimensions[3];
#ifdef DEBUG 
dstrc_enter("pss_write_array");
#endif

*ierr=0;
/*----------------------------------------------------- calc name_size */
name_size=strlen(array->name);

if (name_size>9) name_size=9;

/*---------------------------------------------------- right name_size */
write_error = fwrite(&name_size,sizeof(int),1,allfiles.out_pss);
if (write_error!=1) dserror("Error writing pss-file");

/*----------------------------------------------- write name of record */
write_error = fwrite(array->name,sizeof(char),name_size,allfiles.out_pss);
if (write_error!=name_size) dserror("Error writing pss-file");

/*------------------------------------------- make handle and write it */
*handle = allfiles.pss_counter;
(allfiles.pss_counter)++;
write_error = fwrite(handle,sizeof(int),1,allfiles.out_pss);
if (write_error!=1) dserror("Error writing pss-file");

/*--------------------------------------------------- write dimensions */
dimensions[0]=array->fdim;
dimensions[1]=array->sdim;
switch (array->Typ)
{
   case DA:
   dimensions[2]=sizeof(double);
   break;
   case DV:
   dimensions[2]=sizeof(double);
   break;
   case IA:
   dimensions[2]=sizeof(int);
   break;
   case IV:
   dimensions[2]=sizeof(int);
   break;
}
write_error = fwrite(dimensions,sizeof(int),3,allfiles.out_pss);
if (write_error!=3) dserror("Error writing pss-file");

/*--------------------------------------------------- write the record */
switch (array->Typ)
{
   case DA:
      write_error = fwrite(
                           array->a.da[0],
                           sizeof(double),
                           (dimensions[0]*dimensions[1]),
                           allfiles.out_pss
                          );
      if (write_error != dimensions[0]*dimensions[1]) 
      dserror("Error writing pss-file");
      else *ierr=1;
   break;
   case DV:
      write_error = fwrite(
                           array->a.dv,
                           sizeof(double),
                           (dimensions[0]*dimensions[1]),
                           allfiles.out_pss
                          );
      if (write_error != dimensions[0]*dimensions[1]) 
      dserror("Error writing pss-file");
      else *ierr=1;
   break;
   case IA:
      write_error = fwrite(
                           array->a.ia[0],
                           sizeof(int),
                           (dimensions[0]*dimensions[1]),
                           allfiles.out_pss
                          );
      if (write_error != dimensions[0]*dimensions[1]) 
      dserror("Error writing pss-file");
      else *ierr=1;
   break;
   case IV:
      write_error = fwrite(
                           array->a.iv,
                           sizeof(int),
                           (dimensions[0]*dimensions[1]),
                           allfiles.out_pss
                          );
      if (write_error != dimensions[0]*dimensions[1]) 
      dserror("Error writing pss-file");
      else *ierr=1;
   break;
}
fflush(allfiles.out_pss);

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
 |                                                                      |
 *----------------------------------------------------------------------*/
void pss_read_name(char       *name, 
                      int        *fdim, 
                      int        *sdim,
                      int        *byte,
                      void       *ziel,
                      int        *handle, 
                      int        *ierr)
{
long int      cur_pos = 0;
long int      offset = 0;
int           i;
int           err=0;
int           foundit=0;
int           name_size=0;
int           sizename=0;
int           write_error=0;
int           dimensions[3];
char          test_name[200];
int           handle_dummy;
#ifdef DEBUG 
dstrc_enter("pss_read_name");
#endif

*ierr=0;
/*------------------------------- save current position of file-pointer */
cur_pos = ftell(allfiles.out_pss);

/*----------------------------------------------------- calc name_size */
name_size=strlen(name);

/*---------------------------------- check whether name already exists */
pss_chck(name,&handle_dummy,&err);
if (err!=1)
{
   *ierr=2;
   goto end;
}
/*--------------------------------------------------------- rewind file */
rewind(allfiles.out_pss);

/*------------------------------- loop all records till record is found */
do
{
/*--------------------------------------------------- read size of name */
err=fread(&sizename,sizeof(int),1,allfiles.out_pss);
/*----------------------------------------------- check for end of file */
if ( feof(allfiles.out_pss) != 0)
{
   foundit=0;
   err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=2;
   goto end;    
}
/*------------------------------------------------------ read test_name */ 
err=fread(test_name,sizeof(char),sizename,allfiles.out_pss);
if (err != sizename) dserror("error reading pss-file");

/*------------------------------------ check the test_name against name */
if ( strncmp(name,test_name,strlen(name)) != 0 )
{
/*-------------------- not the right record, so jump it and start again */
   err=fread(&handle_dummy,sizeof(int),1,allfiles.out_pss);
   if (err != 1) dserror("error reading pss-file");
   err=fread(dimensions,sizeof(int),3,allfiles.out_pss);
   if (err != 3) dserror("error reading pss-file");
   offset = dimensions[0]*dimensions[1]*dimensions[2];
   err = fseek(allfiles.out_pss,offset,SEEK_CUR);
   if (err == -1) dserror("error reading pss-file");
}
else
{
   foundit=1;
/*--------------------------------------------------------- read handle */
   err=fread(handle,sizeof(int),1,allfiles.out_pss);
   if (err != 1) dserror("error reading pss-file");
/*--------------------------------------- read the dimensions of record */
   err=fread(dimensions,sizeof(int),3,allfiles.out_pss);
   if (err != 3) dserror("error reading pss-file");
   *fdim=dimensions[0];
   *sdim=dimensions[1];
   *byte=dimensions[2];
/*--------------------------------------------------------- read record */
   err=fread(ziel,(*byte),(*fdim)*(*sdim),allfiles.out_pss);
   if (err != (*fdim)*(*sdim)) dserror("error reading pss-file");
/*---------------------------------- reset the position of file-pointer */
   err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=1;    
}
} while (foundit != 1);

end:
err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
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
                             int        *fdim, 
                             int        *sdim,
                             int        *byte,
                             void       *ziel, 
                             int        *handle, 
                             int        *ierr)
{
long int      cur_pos=0;
long int      offset=0;
int           i;
int           err=0;
int           foundit=0;
int           name_size=0;
int           sizename=0;
int           write_error=0;
int           dimensions[3];
int           handle_dummy;
char          test_name[200];
#ifdef DEBUG 
dstrc_enter("pss_read_name_handle");
#endif

*ierr=0;
/*------------------------------- save current position of file-pointer */
cur_pos = ftell(allfiles.out_pss);

/*----------------------------------------------------- calc name_size */
name_size=strlen(name);

/*------------------------------- check whether name and handle exists */
handle_dummy = *handle;
pss_chck_handle(name,&handle_dummy,&err);
if (err!=1)
{
   *ierr=2;
   goto end;
}
/*--------------------------------------------------------- rewind file */
rewind(allfiles.out_pss);

/*------------------------------- loop all records till record is found */
do
{
/*--------------------------------------------------- read size of name */
err=fread(&sizename,sizeof(int),1,allfiles.out_pss);
/*----------------------------------------------- check for end of file */
if ( feof(allfiles.out_pss) != 0)
{
   foundit=0;
   err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=2;
   goto end;    
}
/*------------------------------------------------------ read test_name */ 
err=fread(test_name,sizeof(char),sizename,allfiles.out_pss);
if (err != sizename) dserror("error reading pss-file");
/*--------------------------------------------------- read handle_dummy */
   err=fread(&handle_dummy,sizeof(int),1,allfiles.out_pss);
   if (err != 1) dserror("error reading pss-file");
/*---- check the test_name against name and handle_dummy against handle */
if ( strncmp(name,test_name,strlen(name)) != 0 || handle_dummy != *handle)
{
/*-------------------- not the right record, so jump it and start again */
   err=fread(dimensions,sizeof(int),3,allfiles.out_pss);
   if (err != 3) dserror("error reading pss-file");
   offset = dimensions[0]*dimensions[1]*dimensions[2];
   err = fseek(allfiles.out_pss,offset,SEEK_CUR);
   if (err == -1) dserror("error reading pss-file");
}
else
{
   foundit=1;
/*--------------------------------------- read the dimensions of record */
   err=fread(dimensions,sizeof(int),3,allfiles.out_pss);
   if (err != 3) dserror("error reading pss-file");
   *fdim=dimensions[0];
   *sdim=dimensions[1];
   *byte=dimensions[2];
/*--------------------------------------------------------- read record */
   err=fread(ziel,(*byte),(*fdim)*(*sdim),allfiles.out_pss);
   if (err != (*fdim)*(*sdim)) dserror("error reading pss-file");
/*---------------------------------- reset the position of file-pointer */
   err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=1;    
}
} while (foundit != 1);

end:
err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
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
 |        space is allocated in correct amount to hold record           |
 |        by this routine                                               |
 |                                                                      |
 *----------------------------------------------------------------------*/
void pss_read_array_name(char       *name, 
                            ARRAY      *array,
                            int        *handle,
                            int        *ierr)
{
long int     cur_pos=0;
long int     offset=0;
int          i;
int          err=0;
int          foundit=0;
int          name_size=0; 
int          sizename=0;
int          write_error=0;
int          dimensions[3];
int          handle_dummy;
char         test_name[200];
#ifdef DEBUG 
dstrc_enter("pss_read_array_name");
#endif

*ierr=0;
/*------------------------------- save current position of file-pointer */
cur_pos = ftell(allfiles.out_pss);

/*----------------------------------------------------- calc name_size */
name_size=strlen(name);

/*---------------------------------- check whether name already exists */
pss_chck(name,&handle_dummy,&err);
if (err!=1)
{
   *ierr=2;
   goto end;
}
/*--------------------------------------------------------- rewind file */
rewind(allfiles.out_pss);

/*------------------------------- loop all records till record is found */
do
{
/*--------------------------------------------------- read size of name */
err=fread(&sizename,sizeof(int),1,allfiles.out_pss);
/*----------------------------------------------- check for end of file */
if ( feof(allfiles.out_pss) != 0)
{
   foundit=0;
   err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=2;
   goto end;    
}
/*------------------------------------------------------ read test_name */ 
err=fread(test_name,sizeof(char),sizename,allfiles.out_pss);
if (err != sizename) dserror("error reading pss-file");
/*--------------------------------------------------------- read handle */
err=fread(&handle_dummy,sizeof(int),1,allfiles.out_pss);
if (err != 1) dserror("error reading pss-file");
/*------------------------------------ check the test_name against name */
if ( strncmp(name,test_name,strlen(name)) != 0 )
{
/*-------------------- not the right record, so jump it and start again */
   err=fread(dimensions,sizeof(int),3,allfiles.out_pss);
   if (err != 3) dserror("error reading pss-file");
   offset = dimensions[0]*dimensions[1]*dimensions[2];
   err = fseek(allfiles.out_pss,offset,SEEK_CUR);
   if (err == -1) dserror("error reading pss-file");
}
else
{
   foundit=1;
/*--------------------------------------------------------- read handle */
   *handle = handle_dummy;
/*--------------------------------------- read the dimensions of record */
   err=fread(dimensions,sizeof(int),3,allfiles.out_pss);
   if (err != 3) dserror("error reading pss-file");
   if (dimensions[2] == sizeof(int)) /* it's a int-record */
   {
       if ( (dimensions[0]==1)||(dimensions[1]==1) )/* it's an int-vector */
       {
          amdef(name,array,dimensions[0],dimensions[1],"IV");
          err=fread(&(array->a.iv[0]),dimensions[2],(dimensions[0]*dimensions[1]),allfiles.out_pss);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }
       else /* it's an int-array */
       {
          amdef(name,array,dimensions[0],dimensions[1],"IA");
          err=fread(&(array->a.ia[0][0]),dimensions[2],(dimensions[0]*dimensions[1]),allfiles.out_pss);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }
   }
   else if (dimensions[2] == sizeof(double))/* it's a double-record */
   {
       if ( (dimensions[0]==1)||(dimensions[1]==1) )/* it's a double-vector */
       {
          amdef(name,array,dimensions[0],dimensions[1],"DV");
          err=fread(&(array->a.dv[0]),dimensions[2],(dimensions[0]*dimensions[1]),allfiles.out_pss);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }
       else /* it's a double-array */
       {
          amdef(name,array,dimensions[0],dimensions[1],"DA");
          err=fread(&(array->a.da[0][0]),dimensions[2],(dimensions[0]*dimensions[1]),allfiles.out_pss);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }  
   }
   else /* the size in byte does not fit to int or double, so it's unknown */
   {
      dserror("error reading array from pss-file: unexpected byte size");
   }

/*---------------------------------- reset the position of file-pointer */
   err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=1;    
}
} while (foundit != 1);


end:
err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
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
 |        space is allocated in correct amount to hold record           |
 |        by this routine                                               |
 |                                                                      |
 *----------------------------------------------------------------------*/
void pss_read_array_name_handle(char       *name, 
                                   ARRAY      *array,
                                   int        *handle,
                                   int        *ierr)
{
long int    cur_pos = 0;
long int    offset = 0;
int         i; 
int         err=0; 
int         foundit=0;
int         name_size=0; 
int         sizename=0;
int         write_error=0;
int         dimensions[3];
int         handle_dummy;
char        test_name[200];
#ifdef DEBUG 
dstrc_enter("pss_read_array_name_handle");
#endif

*ierr=0;
/*------------------------------- save current position of file-pointer */
cur_pos = ftell(allfiles.out_pss);

/*----------------------------------------------------- calc name_size */
name_size=strlen(name);

/*---------------------------------- check whether name already exists */
pss_chck(name,&handle_dummy,&err);
if (err!=1)
{
   *ierr=2;
   goto end;
}
/*--------------------------------------------------------- rewind file */
rewind(allfiles.out_pss);

/*------------------------------- loop all records till record is found */
do
{
/*--------------------------------------------------- read size of name */
err=fread(&sizename,sizeof(int),1,allfiles.out_pss);
/*----------------------------------------------- check for end of file */
if ( feof(allfiles.out_pss) != 0)
{
   foundit=0;
   err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=2;
   goto end;    
}
/*------------------------------------------------------ read test_name */ 
err=fread(test_name,sizeof(char),sizename,allfiles.out_pss);
if (err != sizename) dserror("error reading pss-file");
/*--------------------------------------------------------- read handle */
err=fread(&handle_dummy,sizeof(int),1,allfiles.out_pss);
if (err != 1) dserror("error reading pss-file");
/*---- check the test_name against name and handle_dummy against handle */
if ( strncmp(name,test_name,strlen(name)) != 0 || handle_dummy != *handle)
{
/*-------------------- not the right record, so jump it and start again */
   err=fread(dimensions,sizeof(int),3,allfiles.out_pss);
   if (err != 3) dserror("error reading pss-file");
   offset = dimensions[0]*dimensions[1]*dimensions[2];
   err = fseek(allfiles.out_pss,offset,SEEK_CUR);
   if (err == -1) dserror("error reading pss-file");
}
else
{
   foundit=1;
/*--------------------------------------- read the dimensions of record */
   err=fread(dimensions,sizeof(int),3,allfiles.out_pss);
   if (err != 3) dserror("error reading pss-file");
   if (dimensions[2] == sizeof(int)) /* it's a int-record */
   {
       if ( (dimensions[0]==1)||(dimensions[1]==1) )/* it's an int-vector */
       {
          amdef(name,array,dimensions[0],dimensions[1],"IV");
          err=fread(&(array->a.iv[0]),dimensions[2],(dimensions[0]*dimensions[1]),allfiles.out_pss);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }
       else /* it's an int-array */
       {
          amdef(name,array,dimensions[0],dimensions[1],"IA");
          err=fread(&(array->a.ia[0][0]),dimensions[2],(dimensions[0]*dimensions[1]),allfiles.out_pss);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }
   }
   else if (dimensions[2] == sizeof(double))/* it's a double-record */
   {
       if ( (dimensions[0]==1)||(dimensions[1]==1) )/* it's a double-vector */
       {
          amdef(name,array,dimensions[0],dimensions[1],"DV");
          err=fread(&(array->a.dv[0]),dimensions[2],(dimensions[0]*dimensions[1]),allfiles.out_pss);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }
       else /* it's a double-array */
       {
          amdef(name,array,dimensions[0],dimensions[1],"DA");
          err=fread(&(array->a.da[0][0]),dimensions[2],(dimensions[0]*dimensions[1]),allfiles.out_pss);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }  
   }
   else /* the size in byte does not fit to int or double, so it's unknown */
   {
      dserror("error reading array from pss-file: unexpected byte size");
   }

/*---------------------------------- reset the position of file-pointer */
   err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=1;    
}
} while (foundit != 1);


end:
err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
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
 |        space is allocated in correct amount to hold record           |
 |        by this routine                                               |
 |                                                                      |
 *----------------------------------------------------------------------*/
void pss_read_array_handle(ARRAY      *array,
                              int        *handle,
                              int        *ierr)
{
long int    cur_pos = 0;
long int    offset = 0;
int         i; 
int         err=0; 
int         foundit=0;
int         name_size=0; 
int         sizename=0;
int         write_error=0;
int         dimensions[3];
int         handle_dummy;
char        test_name[200];
char        *name = "         "; 
#ifdef DEBUG 
dstrc_enter("pss_read_array_handle");
#endif

*ierr=0;
/*------------------------------- save current position of file-pointer */
cur_pos = ftell(allfiles.out_pss);

/*----------------------------------------------------- calc name_size */
name_size=strlen(name);

/*--------------------------------------------------------- rewind file */
rewind(allfiles.out_pss);

/*------------------------------- loop all records till record is found */
do
{
/*--------------------------------------------------- read size of name */
err=fread(&sizename,sizeof(int),1,allfiles.out_pss);
/*----------------------------------------------- check for end of file */
if ( feof(allfiles.out_pss) != 0)
{
   foundit=0;
   err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=2;
   goto end;    
}
/*------------------------------------------------------ read test_name */ 
err=fread(test_name,sizeof(char),sizename,allfiles.out_pss);
if (err != sizename) dserror("error reading pss-file");
/*--------------------------------------------------------- read handle */
err=fread(&handle_dummy,sizeof(int),1,allfiles.out_pss);
if (err != 1) dserror("error reading pss-file");
/*---------------------------------- check  handle_dummy against handle */
if (handle_dummy != *handle)
{
/*-------------------- not the right record, so jump it and start again */
   err=fread(dimensions,sizeof(int),3,allfiles.out_pss);
   if (err != 3) dserror("error reading pss-file");
   offset = dimensions[0]*dimensions[1]*dimensions[2];
   err = fseek(allfiles.out_pss,offset,SEEK_CUR);
   if (err == -1) dserror("error reading pss-file");
}
else
{
   foundit=1;
/*------------------------------------------ copy the name of the array */
   sizename = IMIN(sizename,9);
   strncpy(name,test_name,sizename);
/*--------------------------------------- read the dimensions of record */
   err=fread(dimensions,sizeof(int),3,allfiles.out_pss);
   if (err != 3) dserror("error reading pss-file");
   if (dimensions[2] == sizeof(int)) /* it's a int-record */
   {
       if ( (dimensions[0]==1)||(dimensions[1]==1) )/* it's an int-vector */
       {
          amdef(name,array,dimensions[0],dimensions[1],"IV");
          err=fread(&(array->a.iv[0]),dimensions[2],(dimensions[0]*dimensions[1]),allfiles.out_pss);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }
       else /* it's an int-array */
       {
          amdef(name,array,dimensions[0],dimensions[1],"IA");
          err=fread(&(array->a.ia[0][0]),dimensions[2],(dimensions[0]*dimensions[1]),allfiles.out_pss);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }
   }
   else if (dimensions[2] == sizeof(double))/* it's a double-record */
   {
       if ( (dimensions[0]==1)||(dimensions[1]==1) )/* it's a double-vector */
       {
          amdef(name,array,dimensions[0],dimensions[1],"DV");
          err=fread(&(array->a.dv[0]),dimensions[2],(dimensions[0]*dimensions[1]),allfiles.out_pss);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }
       else /* it's a double-array */
       {
          amdef(name,array,dimensions[0],dimensions[1],"DA");
          err=fread(&(array->a.da[0][0]),dimensions[2],(dimensions[0]*dimensions[1]),allfiles.out_pss);
          if (err != (dimensions[0]*dimensions[1]))
          dserror("error reading pss-file");
       }  
   }
   else /* the size in byte does not fit to int or double, so it's unknown */
   {
      dserror("error reading array from pss-file: unexpected byte size");
   }
   *ierr=1;    
}
} while (foundit != 1);


end:
err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
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
                 int        *handle, 
                 int        *ierr)
{
long int    cur_pos=0;
long int    offset=0;
int         foundit=0;
int         sizename=0;
int         err=0;
char        test_name[200];
int         dimensions[3];
int         handle_dummy;
#ifdef DEBUG 
dstrc_enter("pss_chck");
#endif

/*------------------------------- save current position of file-pointer */
cur_pos = ftell(allfiles.out_pss);

/*--------------------------------------------------------- rewind file */
rewind(allfiles.out_pss);

/*------------------------------- loop all records till record is found */
do
{
/*--------------------------------------------------- read size of name */
err=fread(&sizename,sizeof(int),1,allfiles.out_pss);
/*----------------------------------------------- check for end of file */
if ( feof(allfiles.out_pss) != 0)
{
   foundit=0;
   err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=0;
   *handle=-1;
   goto end;    
}
/*------------------------------------------------------ read test_name */ 
err=fread(test_name,sizeof(char),sizename,allfiles.out_pss);
if (err != sizename) dserror("error reading pss-file");

/*------------------------------------ check the test_name against name */
if ( strncmp(name,test_name,strlen(name)) != 0 )
{
/*-------------------- not the right record, so jump it and start again */
   err = fread(&handle_dummy,sizeof(int),1,allfiles.out_pss);
   if (err != 1) dserror("error reading pss-file");
   err=fread(dimensions,sizeof(int),3,allfiles.out_pss);
   if (err != 3) dserror("error reading pss-file");
   offset = dimensions[0]*dimensions[1]*dimensions[2];
   err = fseek(allfiles.out_pss,offset,SEEK_CUR);
   if (err == -1) dserror("error reading pss-file");
}
else
{
   foundit=1;
/*----------------------------------------------------- read the handle */   
   err = fread(handle,sizeof(int),1,allfiles.out_pss);
   if (err != 1) dserror("error reading pss-file");
/*---------------------------------- reset the position of file-pointer */
   err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=1;    
}
} while (foundit != 1);
end:
err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
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
                        int        *handle, 
                        int        *ierr)
{
long int    cur_pos=0;
long int    offset=0;
int         foundit=0;
int         sizename=0;
int         err=0;
char        test_name[200];
int         dimensions[3];
int         handle_dummy;
#ifdef DEBUG 
dstrc_enter("pss_chck_handle");
#endif

/*------------------------------- save current position of file-pointer */
cur_pos = ftell(allfiles.out_pss);

/*--------------------------------------------------------- rewind file */
rewind(allfiles.out_pss);

/*------------------------------- loop all records till record is found */
do
{
/*--------------------------------------------------- read size of name */
err=fread(&sizename,sizeof(int),1,allfiles.out_pss);
/*----------------------------------------------- check for end of file */
if ( feof(allfiles.out_pss) != 0)
{
   foundit=0;
   err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=0;
   *handle=-1;
   goto end;    
}
/*------------------------------------------------------ read test_name */ 
err=fread(test_name,sizeof(char),sizename,allfiles.out_pss);
if (err != sizename) dserror("error reading pss-file");
/*----------------------------------------------------- read the handle */   
err = fread(&handle_dummy,sizeof(int),1,allfiles.out_pss);
if (err != 1) dserror("error reading pss-file");
/*---- check the test_name against name and handle_dummy against handle */
if ( strncmp(name,test_name,strlen(name)) != 0 && handle_dummy != *handle)
{
/*-------------------- not the right record, so jump it and start again */
   err=fread(dimensions,sizeof(int),3,allfiles.out_pss);
   if (err != 3) dserror("error reading pss-file");
   offset = dimensions[0]*dimensions[1]*dimensions[2];
   err = fseek(allfiles.out_pss,offset,SEEK_CUR);
   if (err == -1) dserror("error reading pss-file");
}
else
{
   foundit=1;
/*---------------------------------- reset the position of file-pointer */
   err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=1;    
}
} while (foundit != 1);
end:
err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
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
                         int        *fdim,
                         int        *sdim,
                         int        *byte,
                         int        *handle,
                         int        *ierr)
{
long int      cur_pos=0;
long int      offset=0;
int           foundit=0;
int           sizename=0;
int           err=0;
char          test_name[200];
int           dimensions[3];
int           handle_dummy;
#ifdef DEBUG 
dstrc_enter("pss_getdims_name");
#endif

/*------------------------------- save current position of file-pointer */
cur_pos = ftell(allfiles.out_pss);

/*--------------------------------------------------------- rewind file */
rewind(allfiles.out_pss);

/*------------------------------- loop all records till record is found */
do
{
/*--------------------------------------------------- read size of name */
err=fread(&sizename,sizeof(int),1,allfiles.out_pss);
/*----------------------------------------------- check for end of file */
if ( feof(allfiles.out_pss) != 0)
{
   foundit=0;
   err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=1;
   goto end;    
}
/*------------------------------------------------------ read test_name */ 
err=fread(test_name,sizeof(char),sizename,allfiles.out_pss);
if (err != sizename) dserror("error reading pss-file");
/*----------------------------------------------- read the handle_dummy */
err=fread(&handle_dummy,sizeof(int),1,allfiles.out_pss);
if (err != 1) dserror("error reading pss-file");
/*------------------------------------ check the test_name against name */
if ( strncmp(name,test_name,strlen(name)) != 0 )
{
/*-------------------- not the right record, so jump it and start again */
   err=fread(dimensions,sizeof(int),3,allfiles.out_pss);
   if (err != 3) dserror("error reading pss-file");
   offset = dimensions[0]*dimensions[1]*dimensions[2];
   err = fseek(allfiles.out_pss,offset,SEEK_CUR);
   if (err == -1) dserror("error reading pss-file");
}
else
{
   foundit=1;
   err=fread(dimensions,sizeof(int),3,allfiles.out_pss);
   if (err != 3) dserror("error reading pss-file");
   *fdim=dimensions[0];
   *sdim=dimensions[1];
   *byte=dimensions[2];
   *handle = handle_dummy;
/*---------------------------------- reset the position of file-pointer */
   err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=0;    
}
} while (foundit != 1);

end:
err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
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
                                int        *fdim,
                                int        *sdim,
                                int        *byte,
                                int        *handle,
                                int        *ierr)
{
long int      cur_pos=0;
long int      offset=0;
int           foundit=0;
int           sizename=0;
int           err=0;
char          test_name[200];
int           dimensions[3];
int           handle_dummy;
#ifdef DEBUG 
dstrc_enter("pss_getdims_name_handle");
#endif

/*------------------------------- save current position of file-pointer */
cur_pos = ftell(allfiles.out_pss);

/*--------------------------------------------------------- rewind file */
rewind(allfiles.out_pss);

/*------------------------------- loop all records till record is found */
do
{
/*--------------------------------------------------- read size of name */
err=fread(&sizename,sizeof(int),1,allfiles.out_pss);
/*----------------------------------------------- check for end of file */
if ( feof(allfiles.out_pss) != 0)
{
   foundit=0;
   err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=1;
   goto end;    
}
/*------------------------------------------------------ read test_name */ 
err=fread(test_name,sizeof(char),sizename,allfiles.out_pss);
if (err != sizename) dserror("error reading pss-file");
/*----------------------------------------------- read the handle_dummy */
err=fread(&handle_dummy,sizeof(int),1,allfiles.out_pss);
if (err != 1) dserror("error reading pss-file");
/*----- check the test_name against name and danle_dummy against handle */
if ( strncmp(name,test_name,strlen(name)) != 0 || handle_dummy != *handle)
{
/*-------------------- not the right record, so jump it and start again */
   err=fread(dimensions,sizeof(int),3,allfiles.out_pss);
   if (err != 3) dserror("error reading pss-file");
   offset = dimensions[0]*dimensions[1]*dimensions[2];
   err = fseek(allfiles.out_pss,offset,SEEK_CUR);
   if (err == -1) dserror("error reading pss-file");
}
else
{
   foundit=1;
   err=fread(dimensions,sizeof(int),3,allfiles.out_pss);
   if (err != 3) dserror("error reading pss-file");
   *fdim=dimensions[0];
   *sdim=dimensions[1];
   *byte=dimensions[2];
/*---------------------------------- reset the position of file-pointer */
   err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
   if (err == -1) dserror("error reading pss-file");
   *ierr=0;    
}
} while (foundit != 1);

end:
err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
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
void pss_status_to_err()
{
long int cur_pos=0;
long int current=0;
long int endoffile=0;
long int offset=0;
int      sizename=0;
int      err=0;
char     test_name[200];
char     writename[20];
int      dimensions[3];
int      counter=0;
int      handle;
#ifdef DEBUG 
dstrc_enter("pss_status_to_err");
#endif

fprintf(allfiles.out_err,"===========================================\n");
fprintf(allfiles.out_err,"pss-status - record report about pss-file\n");
fprintf(allfiles.out_err,"===========================================\n");
/*------------------------------- save current position of file-pointer */
cur_pos = ftell(allfiles.out_pss);
fseek(allfiles.out_pss,0,SEEK_END);
endoffile = ftell(allfiles.out_pss);
/*--------------------------------------------------------- rewind file */
rewind(allfiles.out_pss);

/*------------------------------------ loop all records till end of file*/
do
{
/*--------------------------------------------------- read size of name */
   current=ftell(allfiles.out_pss);
   if (current==endoffile) goto end;
   err=fread(&sizename,sizeof(int),1,allfiles.out_pss);
   if (err != 1) goto end;
/*------------------------------------------------------ read test_name */ 
   err=fread(test_name,sizeof(char),sizename,allfiles.out_pss);
   if (err != sizename) goto end;
/*--------------------------------------------------------- read handle */   
   err=fread(&handle,sizeof(int),1,allfiles.out_pss);
   if (err != 1) goto end;
/*----------------------------------------------------- read dimensions */
   err=fread(dimensions,sizeof(int),3,allfiles.out_pss);
   if (err != 3) goto end;
/*----------------------------------------------------- jump the record */
   offset = dimensions[0]*dimensions[1]*dimensions[2];
   err = fseek(allfiles.out_pss,offset,SEEK_CUR);
   if (feof(allfiles.out_pss) == 1) goto end;
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
err = fseek(allfiles.out_pss,cur_pos,SEEK_SET);
if (err == -1) dserror("error reading pss-file");
fflush(allfiles.out_err);

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of pss_status_to_err */
