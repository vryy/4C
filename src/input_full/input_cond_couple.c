#include "../headers/standardtypes.h"

/*----------------------------------------------------------------------*
 | input of coupling conditions                           m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_couple(FIELD *field)
{
int  ierr, ierr2;
int  i,j;
int  nodeid;
int  coupleID;
int  is_structure=0, is_fluid=0, is_ale=0;
char *colpointer;
int  *couple_dof;
int  *couple_geom;
char  buffer[50];
NODE  *actnode;
#ifdef DEBUG 
dstrc_enter("inp_couple");
#endif
/*----------------------------------------------------------------------*/
couple_dof  = (int*)calloc(genprob.numdf,sizeof(int));
couple_geom = (int*)calloc(genprob.numdf,sizeof(int));
if (couple_dof==NULL || couple_geom==NULL)
dserror("Allocation of temporary fields failed");
/*------------------------------------------------- read the conditions */
frfind("---COUPLING CONDITIONS");
start:
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   is_structure=is_fluid=is_ale=0;
   frchar("FIELD",buffer,&ierr);
   if (strcmp(buffer,"structure")==0 && field->fieldtyp != structure) goto start;
   if (strcmp(buffer,"fluid")==0     && field->fieldtyp != fluid)     goto start;
   if (strcmp(buffer,"ale")==0       && field->fieldtyp != ale)       goto start;
   if (strcmp(buffer,"structure")==0) is_structure=1;
   if (strcmp(buffer,"fluid")==0) is_fluid=1;
   if (strcmp(buffer,"ale")==0) is_ale=1;
   nodeid = strtol(allfiles.actplace,&colpointer,10);
   nodeid--;
   for (i=0; i<field->numnp; i++)
   {
      if (field->node[i].Id == nodeid)
      {
         actnode = &(field->node[i]);
         break;
      }
   }
   if (actnode->c==NULL)
   {
      actnode->c=(COND_NODE*)calloc(1,sizeof(COND_NODE));
      if (actnode->c==NULL) dserror("Allocation of COND_NODE failed"); 
      actnode->c->isneum=0;
      actnode->c->isdirich=0;
      actnode->c->iscoupled=0;
      actnode->c->fsi_iscoupled=0;
      actnode->c->neum_onoff.a.iv=NULL;
      actnode->c->neum_val.a.iv=NULL;
      actnode->c->dirich_onoff.a.iv=NULL;
      actnode->c->dirich_val.a.iv=NULL;
      actnode->c->couple.a.iv=NULL;
   }
   for (i=0; i<genprob.numdf; i++) {couple_dof[i]=0; couple_geom[i]=0;}
/*----------------------- check whether there already exists a coupling */      
   if (actnode->c->iscoupled==0)
   {
      /* all geostationary couplings are written in line 0
         all dof           couplings are written in line 1
         in the line the coupleID is written on the place of 
         the dof to couple */ 
      frchk("GEOMCOUP",&ierr);
      frchk("DOFCOUP",&ierr2);
      ierr+=ierr2;
      if (ierr!=0) 
      {
         actnode->c->iscoupled=1;
         amdef("couple",&(actnode->c->couple),genprob.numdf,4,"IA");
         amzero(&(actnode->c->couple));
      }
   }
/*--------------------------------------------------- read the coupleID */   
   ierr=0;
   frint("COUP_ID",&coupleID,&ierr);
   if (ierr!=1) dserror("Cannot read coupling ID");
/*----------------------------------------- read geostationary coupling */   
   frchk("GEOMCOUP",&ierr);
   if (ierr==1)
   {
      frint_n("GEOMCOUP",couple_geom,genprob.numdf,&ierr);
      for (i=0; i<genprob.numdf; i++)
      {
         if (couple_geom[i]!=0) actnode->c->couple.a.ia[i][0]=coupleID;
      }
   }
/*------------------------------------------------ read the dofcoupling */
   frchk("DOFCOUP",&ierr);
   if (ierr==1)
   {
      frint_n("DOFCOUP",couple_dof,genprob.numdf,&ierr);
      for (i=0; i<genprob.numdf; i++)
      {
         if (couple_dof[i]!=0) actnode->c->couple.a.ia[i][1]=coupleID;
      }
   }
/*----------------------------------------------- read the FSI coupling */
   frchk("FSI",&ierr);
   if (ierr==1) 
   {
      actnode->c->fsi_iscoupled=coupleID;
   }
/*------------------------------------------------------ read next line */
   frread();
}
/*----------------------------------------------------------------------*/
free(couple_dof); free(couple_geom);
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_couple */



