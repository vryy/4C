#include "../headers/standardtypes.h"

/*----------------------------------------------------------------------*
 | init data bunkers                                      m.gee 6/01    |
 *----------------------------------------------------------------------*/
void db_init()
{
#ifdef DEBUG 
dstrc_enter("db_init");
#endif
/*----------------------------------------------------------------------*/
db(1,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of db_init */





/*----------------------------------------------------------------------*
 | create a data bunker                                   m.gee 6/01    |
 *----------------------------------------------------------------------*/
void db_create(int *bunker_Id)
{
int i=0;

#ifdef DEBUG 
dstrc_enter("db_create");
#endif
/*----------------------------------------------------------------------*/
db(2,bunker_Id,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of db_create */





/*----------------------------------------------------------------------*
 | init the contents of one bunker                        m.gee 6/01    |
 *----------------------------------------------------------------------*/
void db_content_init(BUNKER *bptr)
{
int i;
#ifdef DEBUG 
dstrc_enter("db_content_init");
#endif
/*----------------------------------------------------------------------*/
bptr->used=0;
bptr->size=1000;
bptr->shelf          = (ARRAY**)calloc(bptr->size,sizeof(ARRAY*));
bptr->access_handles = (int*)   calloc(bptr->size,sizeof(int));
if (bptr->shelf==NULL || bptr->access_handles==NULL) 
   dserror("Allocation of BUNKER memory failed");
for (i=0; i<bptr->size; i++)
{
   bptr->shelf[i]=NULL;
   bptr->access_handles[i]=-1;
}  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of db_content_init */





/*----------------------------------------------------------------------*
 | create a bunker account for a certain bunker           m.gee 6/01    |
 *----------------------------------------------------------------------*/
void db_create_access(int bunker_Id, int *access)
{
int i=0;

#ifdef DEBUG 
dstrc_enter("db_create_access");
#endif
/*----------------------------------------------------------------------*/
db(3,&bunker_Id,access,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of db_create_access */





/*----------------------------------------------------------------------*
 | enlarge bunker by addaccess                            m.gee 6/01    |
 *----------------------------------------------------------------------*/
void db_bunker_enlarge(BUNKER *bptr, int addaccess)
{
int i;
#ifdef DEBUG 
dstrc_enter("db_bunker_enlarge");
#endif
/*----------------------------------------------------------------------*/
bptr->access_handles = (int*)realloc((void*)(bptr->access_handles),
                                      (bptr->size+addaccess)*sizeof(int));
bptr->shelf = (ARRAY**)realloc((void*)(bptr->shelf),
                               (bptr->size+addaccess)*sizeof(ARRAY*));
if (bptr->access_handles==NULL || bptr->shelf==NULL)
   dserror("Enlargement of BUNKER failed");                                                               
for (i=bptr->size; i<bptr->size+addaccess; i++)
{
   bptr->access_handles[i]=-1;
   bptr->shelf[i]=NULL;
}
bptr->size += addaccess;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of db_bunker_enlarge */





/*----------------------------------------------------------------------*
 | put data to a new entry                                m.gee 6/01    |
 *----------------------------------------------------------------------*/
void db_put_dat_new(int   bunker_Id,
                       int   access,
                       int  *handle,
                       char  string[],
                       void *startadress,
                       int   fdim,
                       int   sdim,
                       char  typ[])
{
int i;
long int size;
int *iptr;
double *dptr;
ARRAY tmp;
#ifdef DEBUG 
dstrc_enter("db_put_dat_new");
#endif
/*----------------------------------------------------------------------*/
amdef(string,&tmp,fdim,sdim,typ);
switch(tmp.Typ)
{
case DA:
dptr = (double*)startadress;
size = fdim*sdim;
for (i=0; i<size; i++)
{
   tmp.a.da[0][i] = dptr[i];
}
break;
case DV:
dptr = startadress;
size = fdim*sdim;
for (i=0; i<size; i++)
{
   tmp.a.dv[i] = dptr[i];
}
break;
case IA:
iptr = (int*)startadress;
size = fdim*sdim;
for (i=0; i<size; i++)
{
   tmp.a.ia[0][i] = iptr[i];
}
break;
case IV:
iptr = (int*)startadress;
size = fdim*sdim;
for (i=0; i<size; i++)
{
   tmp.a.iv[i] = iptr[i];
}
break;
}
db(4,&bunker_Id,&access,handle,&tmp,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
amdel(&tmp);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of db_put_dat_new */





/*----------------------------------------------------------------------*
 | put array to a new entry                               m.gee 6/01    |
 *----------------------------------------------------------------------*/
void db_put_array_new(int    bunker_Id,
                         int    access,
                         int   *handle,
                         ARRAY *array)
{
int i;
#ifdef DEBUG 
dstrc_enter("db_put_array_new");
#endif
/*----------------------------------------------------------------------*/
db(4,&bunker_Id,&access,handle,array,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of db_put_array_new */





/*----------------------------------------------------------------------*
 | put array to an existing handle                        m.gee 6/01    |
 *----------------------------------------------------------------------*/
void db_put_array_overwrite(int    bunker_Id,
                               int    access,
                               int    handle,
                               ARRAY *array)
{
int i;
#ifdef DEBUG 
dstrc_enter("db_put_array_overwrite");
#endif
/*----------------------------------------------------------------------*/
db(7,&bunker_Id,&access,&handle,array,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of db_put_array_overwrite */





/*----------------------------------------------------------------------*
 | probe info about an array                              m.gee 6/01    |
 |  ierr=0 everything o.k.                                              |
 |  ierr=1 bunker doesn't exist                                         |
 |  ierr=2 access doesn't exist                                         |
 |  ierr=3 handle doesn't exist                                         |
 |  ierr=4 ARRAY is somehow damaged                                     |
 |                                                                      |
 *----------------------------------------------------------------------*/
void db_probe_array(int bunker_Id,
                       int  access,
                       int  handle,
                       char name[],
                       int *fdim,
                       int *sdim,
                       char typ[],
                       int *ierr)
{
int i;
#ifdef DEBUG 
dstrc_enter("db_probe_array");
#endif
/*----------------------------------------------------------------------*/
db(5,&bunker_Id,&access,&handle,NULL,fdim,sdim,name,typ,NULL,NULL,ierr);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of db_probe_array */





/*----------------------------------------------------------------------*
 | probe existence about an array                         m.gee 6/01    |
 |  ierr=0 everything o.k.                                              |
 |  ierr=1 bunker doesn't exist                                         |
 |  ierr=2 access doesn't exist                                         |
 |  ierr=3 handle doesn't exist                                         |
 |                                                                      |
 *----------------------------------------------------------------------*/
void db_check_existence(BUNKER *bunker,
                           int     bunker_Id,
                           int     numbunker,
                           int     access,
                           int     handle,
                           int    *ierr)
{
BUNKER *bptr;
#ifdef DEBUG 
dstrc_enter("db_check_existence");
#endif
/*----------------------------------------------------------------------*/
   *ierr= 0;
   if (bunker_Id >= numbunker)
   {
      *ierr=1;
      goto end;
   }
   bptr = &(bunker[bunker_Id]);
   if (access >= bptr->used)
   {
     *ierr=2;
     goto end;
   }
   if (handle >= bptr->access_handles[access])
   {
      *ierr=3;
      goto end;
   }
   if (bptr->shelf[access][handle].Typ == XX)
   {
      *ierr=3;
      goto end;
   }
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of db_check_existence */





/*----------------------------------------------------------------------*
 | allocate an ARRAY and fill data from the bunker to it  m.gee 6/01    |
 |  ierr=0 everything o.k.                                              |
 |  ierr=1 bunker doesn't exist                                         |
 |  ierr=2 access doesn't exist                                         |
 |  ierr=3 handle doesn't exist                                         |
 |                                                                      |
 *----------------------------------------------------------------------*/
void db_get_array_alloc(int    bunker_Id,
                           int    access,
                           int    handle,
                           ARRAY *array,
                           int   *ierr)
{
int i;
char    name[9];
char    typ[3];
int     fdim,sdim;
int    *iptr;
double *dptr;
#ifdef DEBUG 
dstrc_enter("db_get_array_alloc");
#endif
/*----------------------------------------------------------------------*/
db(5,&bunker_Id,&access,&handle,NULL,&fdim,&sdim,name,typ,NULL,NULL,ierr);
if (*ierr!=0) goto end;
amdef(name,array,fdim,sdim,typ);
switch(array->Typ)
{
   case DA:
   dptr = array->a.da[0];
   db(6,&bunker_Id,&access,&handle,NULL,&fdim,&sdim,NULL,NULL,dptr,NULL,ierr);
   break;
   case DV:
   dptr = array->a.dv;
   db(6,&bunker_Id,&access,&handle,NULL,&fdim,&sdim,NULL,NULL,dptr,NULL,ierr);
   break;
   case IA:
   iptr = array->a.ia[0];
   db(6,&bunker_Id,&access,&handle,NULL,&fdim,&sdim,NULL,NULL,iptr,NULL,ierr);
   break;
   case IV:
   iptr = array->a.iv;
   db(6,&bunker_Id,&access,&handle,NULL,&fdim,&sdim,NULL,NULL,iptr,NULL,ierr);
   break;
}
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of db_get_array_alloc */





/*----------------------------------------------------------------------*
 | get data from bunker to an existing array              m.gee 6/01    |
 |  user must ensure enough allocated space and typ of arary            |
 |  ierr=1 bunker doesn't exist                                         |
 |  ierr=2 access doesn't exist                                         |
 |  ierr=3 handle doesn't exist                                         |
 |                                                                      |
 *----------------------------------------------------------------------*/
void db_get_array(int    bunker_Id,
                     int    access,
                     int    handle,
                     ARRAY *array,
                     int   *ierr)
{
int     i;
char    name[9];
char    typ[3];
int     fdim,sdim;
int    *iptr;
double *dptr;
#ifdef DEBUG 
dstrc_enter("db_get_array");
#endif
/*----------------------------------------------------------------------*/
switch(array->Typ)
{
   case DA:
   dptr = array->a.da[0];
   db(6,&bunker_Id,&access,&handle,NULL,&fdim,&sdim,NULL,NULL,dptr,NULL,ierr);
   break;
   case DV:
   dptr = array->a.dv;
   db(6,&bunker_Id,&access,&handle,NULL,&fdim,&sdim,NULL,NULL,dptr,NULL,ierr);
   break;
   case IA:
   iptr = array->a.ia[0];
   db(6,&bunker_Id,&access,&handle,NULL,&fdim,&sdim,NULL,NULL,iptr,NULL,ierr);
   break;
   case IV:
   iptr = array->a.iv;
   db(6,&bunker_Id,&access,&handle,NULL,&fdim,&sdim,NULL,NULL,iptr,NULL,ierr);
   break;
}
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of db_get_array */





/*----------------------------------------------------------------------*
 | get data                                               m.gee 6/01    |
 |  ierr=0 everything o.k.                                              |
 |  ierr=1 bunker doesn't exist                                         |
 |  ierr=2 access doesn't exist                                         |
 |  ierr=3 handle doesn't exist                                         |
 |                                                                      |
 *----------------------------------------------------------------------*/
void db_get_dat(int bunker_Id,
                   int  access,
                   int  handle,
                   void *data,
                   int  *fdim,
                   int  *sdim,
                   int  *ierr)
{
int i;
#ifdef DEBUG 
dstrc_enter("db_get_dat");
#endif
/*----------------------------------------------------------------------*/
db(6,&bunker_Id,&access,&handle,NULL,fdim,sdim,NULL,NULL,data,NULL,ierr);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of db_get_dat */





/*----------------------------------------------------------------------*
 | delete an entry, return -1 for handle                  m.gee 6/01    |
 |                                                                      |
 *----------------------------------------------------------------------*/
void db_del_entry(int  bunker_Id,
                     int  access,
                     int *handle)
{
int i;
#ifdef DEBUG 
dstrc_enter("db_del_entry");
#endif
/*----------------------------------------------------------------------*/
db(8,&bunker_Id,&access,handle,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of db_del_entry */





/*----------------------------------------------------------------------*
 | delete an access, return -1 for access                 m.gee 6/01    |
 | also delete all arrays which have been associated with this access   |
 |                                                                      |
 *----------------------------------------------------------------------*/
void db_del_access(int  bunker_Id,
                      int *access)
{
int i;
#ifdef DEBUG 
dstrc_enter("db_del_access");
#endif
/*----------------------------------------------------------------------*/
db(9,&bunker_Id,access,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of db_del_access */





/*----------------------------------------------------------------------*
 | delete a bunker                                        m.gee 6/01    |
 | also delete all arrays which have been associated with this bunker   |
 |                                                                      |
 *----------------------------------------------------------------------*/
void db_del_bunker(int *bunker_Id)
{
int i;
#ifdef DEBUG 
dstrc_enter("db_del_bunker");
#endif
/*----------------------------------------------------------------------*/
db(10,bunker_Id,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of db_del_bunker */







/*----------------------------------------------------------------------*
 | write a bunker to pss-file and return a pss_handle     m.gee 6/01    |
 | name of the record on pss-file is bunker                             |
 | with pss_handle and name it can be read from pss-file again          |
 |                                                                      |
 *----------------------------------------------------------------------*/
void db_bunker_to_pss(int bunker_Id, int *pss_handle)
{
int i;
#ifdef DEBUG 
dstrc_enter("db_bunker_to_pss");
#endif
/*----------------------------------------------------------------------*/
db(11,&bunker_Id,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,pss_handle,NULL);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of db_bunker_to_pss */





/*----------------------------------------------------------------------*
 | read a bunker from pss-file                            m.gee 6/01    |
 | name of the record on pss-file is bunker                             |
 | with pss_handle and name it can be read from pss-file again          |
 |                                                                      |
 *----------------------------------------------------------------------*/
void db_bunker_from_pss(int bunker_Id, int pss_handle, int *err)
{
int i;
#ifdef DEBUG 
dstrc_enter("db_bunker_from_pss");
#endif
/*----------------------------------------------------------------------*/
db(12,&bunker_Id,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,&pss_handle,err);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of db_bunker_from_pss */
