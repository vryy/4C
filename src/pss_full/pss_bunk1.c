#include "../headers/standardtypes.h"

/*----------------------------------------------------------------------*
 | data bunker                                            m.gee 6/01    |
 *----------------------------------------------------------------------*/
void db(
           int      option,
           int     *bunker_Id,
           int     *access,
           int     *handle,
           ARRAY   *array,
           int     *fdim,
           int     *sdim,
           char     name[],
           char     typ[],
           void    *data,
           int     *pss_handle,
           int     *ierr
          )
{
int            i,j;
int            minusone=-1;
int            one=1;
int            byte;
int           *iptr;
double        *dptr;
int            size;
static int     numbunker;
static BUNKER *bunker;
BUNKER        *bptr;
ARRAY         *aptr;
ARRAY          access_handles_picture;
ARRAY          shelf_picture;
int            found;
int            pss_err;
#ifdef DEBUG 
dstrc_enter("db");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------------- bunker initialisation */
/* option==1 is called from

void db_init();

*/
switch (option)
{
/*------------------------------------------------------ init DB-System */
/* option==1 is called from

void db_init();

*/
case 1:
   numbunker=0;
   bunker=NULL;
   bptr=NULL;
break;
/*------------------------ create a new bunker and return the bunker_Id */
/* option==2 is called from

void db_create(int *bunker_Id);

*/
case 2:
   if (numbunker==0)
   {
      bunker = (BUNKER*)calloc(1,sizeof(BUNKER));
      if (bunker==NULL) dserror("Creation of BUNKER failed");
      numbunker=1;
      *bunker_Id = 0;
      bunker[0].bunker_Id=0;
      bptr = &(bunker[0]);
      db_content_init(bptr);
   }
   else
   {
      numbunker++;
      bunker = (BUNKER*)realloc((void*)bunker,(size_t)(numbunker*sizeof(BUNKER)));
      if (bunker==NULL) dserror("Creation of BUNKER failed");
      *bunker_Id = numbunker-1;
      bptr = &(bunker[numbunker-1]);
      bptr->bunker_Id=numbunker-1;
      db_content_init(bptr);
   }
break;
/*---------------------------- create access_number to a certain bunker */
/* option==3 is called from

void db_create_access(int bunker_Id, int *access);

*/
case 3:
   bptr = &(bunker[*bunker_Id]);
   *access = bptr->used;
   bptr->used += 1;

   if (bptr->used == bptr->size) db_bunker_enlarge(bptr,1000);

   bptr->shelf[*access]          = (ARRAY*)calloc(3,sizeof(ARRAY));
   if (bptr->shelf[*access]==NULL) dserror("Allocation of BUNKER memory failed");
   bptr->access_handles[*access] = 3;
   for (i=0; i<3; i++) 
   {
      bptr->shelf[*access][i].Typ = XX;
   }
break;
/*----------------------------------------------- put a new array entry */
/* option==4 is called from

void db_put_dat_new(int   bunker_Id,
                       int   access,
                       int  *handle,
                       char  string[],
                       void *startadress,
                       int   fdim,
                       int   sdim,
                       char  typ[]);

void db_put_array_new(int    bunker_Id,
                         int    access,
                         int   *handle,
                         ARRAY *array);

*/
case 4:
   bptr = &(bunker[*bunker_Id]);
/* find a free space */
   found=0;
   for (i=0; i<bptr->access_handles[*access]; i++)
   {
      if (bptr->shelf[*access][i].Typ==XX)
      {
         found=1;
         break;
      }
   }
   if (found==0) /* no free space, so enlarge by three */
   {
      bptr->shelf[*access] = 
      (ARRAY*)realloc((void*)(bptr->shelf[*access]),(bptr->access_handles[*access]+3)*sizeof(ARRAY));
      if (bptr->shelf[*access]==NULL) dserror("Enlargement of bunker failed");
      for (i=bptr->access_handles[*access]; i<bptr->access_handles[*access]+3; i++)
      {
         bptr->shelf[*access][i].Typ=XX;
      }
      *handle = bptr->access_handles[*access];
      bptr->access_handles[*access] += 3;
   }
   else
   {
      *handle = i;
   }
   am_alloc_copy(array,&(bptr->shelf[*access][*handle]));
break;
/*--------------------- probe existence and values of an existing array */
/* option==5 is called from

void db_probe_array(int bunker_Id,
                       int  access,
                       int  handle,
                       char name[],
                       int *fdim,
                       int *sdim,
                       char typ[],
                       int *ierr);

*/
case 5:
   db_check_existence(bunker,*bunker_Id,numbunker,*access,*handle,ierr);
   if (*ierr != 0) break;
   bptr = &(bunker[*bunker_Id]);
   aptr = &(bptr->shelf[*access][*handle]);
   strcpy(name,aptr->name);
   *fdim = aptr->fdim;
   *sdim = aptr->sdim;
   switch(aptr->Typ)
   {
      case DA:
      strcpy(typ,"DA");
      break;
      case DV:
      strcpy(typ,"DV");
      break;
      case IA:
      strcpy(typ,"IA");
      break;
      case IV:
      strcpy(typ,"IV");
      break;
      default:
      strcpy(typ,"XX");
      *ierr=4;
      break;
   }
break;
/*------------------------------------------------------------ get data */
/* option==6 is called from

void db_get_dat(int bunker_Id,
                   int  access,
                   int  handle,
                   void *data,
                   int  *fdim,
                   int  *sdim,
                   int  *ierr);

void db_get_array_alloc(int    bunker_Id,
                           int    access,
                           int    handle,
                           ARRAY *array,
                           int   *ierr);

void db_get_array(int    bunker_Id,
                           int    access,
                           int    handle,
                           ARRAY *array,
                           int   *ierr);

*/
case 6:
   db_check_existence(bunker,*bunker_Id,numbunker,*access,*handle,ierr);
   if (*ierr != 0) break;
   bptr = &(bunker[*bunker_Id]);
   aptr = &(bptr->shelf[*access][*handle]);
   *fdim = aptr->fdim;
   *sdim = aptr->sdim;
   size = (aptr->fdim)*(aptr->sdim);
   switch(aptr->Typ)
   {
      case DA:
      dptr = (double*)data;
      for (i=0; i<size; i++) dptr[i] = aptr->a.da[0][i];
      break;
      case DV:
      dptr = (double*)data;
      for (i=0; i<size; i++) dptr[i] = aptr->a.dv[i];
      break;
      case IA:
      iptr = (int*)data;
      for (i=0; i<size; i++) iptr[i] = aptr->a.ia[0][i];
      break;
      case IV:
      iptr = (int*)data;
      for (i=0; i<size; i++) iptr[i] = aptr->a.iv[i];
      break;
   }
break;
/*----------------------------------- overwrite an existing array entry */
/* option==7 is called from

void db_put_array_overwrite(int    bunker_Id,
                               int    access,
                               int    handle,
                               ARRAY *array);

*/
case 7:
   bptr = &(bunker[*bunker_Id]);
   amdel(&(bptr->shelf[*access][*handle]));
   am_alloc_copy(array,&(bptr->shelf[*access][*handle]));
break;
/*---------------------------------------------------- delete an array */
/* option==8 is called from

void db_del_entry(int  bunker_Id,
                     int  access,
                     int *handle);

*/
case 8:
   bptr = &(bunker[*bunker_Id]);
   amdel(&(bptr->shelf[*access][*handle]));
   *handle = -1;
break;
/*--------------------------------delete an access and all his arrays */
/* option==9 is called from

void db_del_access(int  bunker_Id,
                      int *access);

*/
case 9:
   bptr = &(bunker[*bunker_Id]);
   for (i=0; i<bptr->access_handles[*access]; i++)
   {
      if (bptr->shelf[*access][i].Typ != XX)
      amdel(&(bptr->shelf[*access][i]));
   }
   free(bptr->shelf[*access]);
   bptr->shelf[*access]=NULL;
   bptr->access_handles[*access]=-1;
   *access=-1;
break;
/*--------------------------------delete a bunker and all arrays     */
/* option==10 is called from

void db_del_bunker(int *bunker_Id);

*/
case 10:
   bptr = &(bunker[*bunker_Id]);
   for (i=0; i<(bptr->used); i++)
   {
      if (bptr->shelf[i]!=NULL)
      {
         for (j=0; j<bptr->access_handles[i]; j++)
         {
            if (bptr->shelf[i][j].Typ != XX)
            {
               amdel(&(bptr->shelf[i][j]));
            }
         }
         free(bptr->shelf[i]);
      }
   }
   free(bptr->shelf);
   free(bptr->access_handles);
   bptr->shelf=NULL;
   bptr->access_handles=NULL;
   bptr->size=0;
   bptr->used=0;
break;
/*--------------------------------------------write bunker to pss-file  */
/* option==11 is called from

void db_bunker_to_pss(int bunker_Id, int *pss_handle);

*/
case 11:
   bptr = &(bunker[*bunker_Id]);
   /*-------------------------- size of access_handles_picture is known */
   amdef("ahp_pic",&access_handles_picture,bptr->used,1,"IV");
   size=0;
   for (i=0; i<bptr->used; i++)
   {
      access_handles_picture.a.iv[i] = bptr->access_handles[i];
      if (size<bptr->access_handles[i]) size = bptr->access_handles[i];
   } 
   /*--------------------- now the max. size of the shelf is also known */
   amdef("shelf_p",&shelf_picture,bptr->used,size,"IA");
   aminit(&shelf_picture,&minusone);
   /*------------------------- write access_handles_picture to pss-file */
   pss_write_array(
                      &access_handles_picture,
                      &(bptr->ahp_pss_handle),
                      &pss_err
                     );
   if (pss_err != 1) dserror("Writing of data bunker failed");
   /*------------------------ now loop the shelf and write the contents */
   for (i=0; i<bptr->used; i++)
   {
      if (bptr->shelf[i]!=NULL || bptr->access_handles[i] != -1)
      {
         for (j=0; j<bptr->access_handles[i]; j++)
         {
            if (bptr->shelf[i][j].Typ != XX)
            {
               pss_write_array(
                                  &(bptr->shelf[i][j]),
                                  &(shelf_picture.a.ia[i][j]),
                                  &pss_err
                                 );
               if (pss_err != 1) dserror("Writing of data bunker failed");
            }
         }
      }
   }
   /*-------------------------------------- now write the shelf_picture */
   pss_write_array(
                      &shelf_picture,
                      &(bptr->sp_pss_handle),
                      &pss_err
                     );
   if (pss_err != 1) dserror("Writing of data bunker failed");
   /*---------------------------------- now delete the temporary arrays */
   amdel(&access_handles_picture);
   amdel(&shelf_picture);
   /*---------------------------- now write the bunker structure itself */                  
   pss_write(
                "bunker",
                1,
                1,
                sizeof(BUNKER),
                bptr,
                pss_handle,
                &pss_err
               );
   if (pss_err != 1) dserror("Writing of data bunker failed");

break;
/*--------------------------------------read a bunker from pss_file     */
/* option==12 is called from

void db_bunker_from_pss(int bunker_Id, int pss_handle, int *err);

*/
case 12:
   if (numbunker > *bunker_Id)/*-------- array of bunkers large enough */
   {
      bptr = &(bunker[*bunker_Id]);
      /*-------------------------- check whether this bunker is in use */
      if (bptr->shelf != NULL) /*---------------------- yes, delete it */ 
      {
          /*------------------------------ this is nice, call myself ! */
          db(10,bunker_Id,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
          bptr = &(bunker[*bunker_Id]);
      }
   }
   else /*-------------------------- array of bunkers not large enough */
   {
      bunker = (BUNKER*)realloc((void*)bunker,(size_t)((*bunker_Id+1)*sizeof(BUNKER)));
      if (bunker==NULL) dserror("Creation of BUNKER failed");
      for (i=numbunker; i<(*bunker_Id+1); i++)
      {
         bptr = &(bunker[i]);
         bptr->shelf=NULL;
         bptr->access_handles=NULL;
         bptr->size=0;
         bptr->used=0;
      }
      numbunker = *bunker_Id + 1;
      bptr = &(bunker[*bunker_Id]);
   }
   /*------------ now we got the right bunker, so read BUNKER from pss */
   pss_read_name_handle(
                           "bunker",
                           &one,
                           &one,
                           &byte,
                           bptr,
                           pss_handle,
                           &pss_err
                          );
   if (pss_err!=1) dserror("Reading of bunker from pss failed");
   /* we now can read the arrays access_handles_picture & shelf_picture */
   pss_read_array_name_handle(
                                 "ahp_pic",
                                 &access_handles_picture,
                                 &(bptr->ahp_pss_handle),
                                 &pss_err
                                );
   if (pss_err!=1) dserror("Reading of bunker from pss failed");                              
   pss_read_array_name_handle(
                                 "shelf_p",
                                 &shelf_picture,
                                 &(bptr->sp_pss_handle),
                                 &pss_err
                                );
   if (pss_err!=1) dserror("Reading of bunker from pss failed");                              
   /*----------- we now can restore the arrays access_handles and shelf */
   bptr->shelf          = (ARRAY**)calloc(bptr->size,sizeof(ARRAY*));
   bptr->access_handles = (int*)   calloc(bptr->size,sizeof(int));
   if (bptr->shelf==NULL || bptr->access_handles==NULL) 
   dserror("Allocation of BUNKER memory failed");
   /*------------------------------------------------ do access_handles */
   for (i=0; i<bptr->used; i++) 
      bptr->access_handles[i] = access_handles_picture.a.iv[i];
   for (i=bptr->used; i<bptr->size; i++)
      bptr->access_handles[i] = -1;
   /*------------------------------------------- now allocate the shelf */
   for (i=0; i<bptr->used; i++) 
   {
      if (bptr->access_handles[i]==-1) continue;
      bptr->shelf[i] = (ARRAY*)calloc(bptr->access_handles[i],sizeof(ARRAY));
      if (!(bptr->shelf[i])) dserror("Allocation of bunker shelf failed");
   }
   /*---------------------------------------------- read shelf from pss */
   for (i=0; i<bptr->used; i++)
   {
      for (j=0; j<bptr->access_handles[i]; j++)
      {
         if (shelf_picture.a.ia[i][j]==-1) continue;
         else
         {
            found = shelf_picture.a.ia[i][j];
            pss_read_array_handle(
                                     &(bptr->shelf[i][j]),
                                     &found,
                                     &pss_err
                                    );
            if (pss_err!=1) dserror("Reading of bunker from pss failed");     
            
         }
      }
   } 
   /*-------------------------------------- delete the temporary arrays */
   amdel(&access_handles_picture);
   amdel(&shelf_picture);

break;
/*--------------------------------------------no specified action found */
default:
   dserror("BUNKER action unclear");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of db */
