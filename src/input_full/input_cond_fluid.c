#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 | input of fluidal nodal conditions                      m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_cond_nodal_fluid(FIELD *field)
{
int  ierr;
int  i,j;
int  nodeid;
int    itmp[50];
double dtmp[50];
double eps=1.0E-08;
char *colpointer;
char buffer[50];
NODE *actnode;
#ifdef DEBUG 
dstrc_enter("inp_cond_nodal_fluid");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------ init the condition pointer */
for (i=0; i<field->dis[0].numnp; i++) field->dis[0].node[i].c==NULL;
/*------------------------------------------------- read the conditions */
frfind("--FLUID NODAL CONDITIONS");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   nodeid = strtol(allfiles.actplace,&colpointer,10);
   nodeid--;
   for (j=0; j<field->dis[0].numnp; j++)
   {
      if (field->dis[0].node[j].Id == nodeid) 
      {
         actnode = &(field->dis[0].node[j]);
         break;
      }
   }
   if (actnode->c==NULL)
   {
      actnode->c=(COND_NODE*)CALLOC(1,sizeof(COND_NODE));
      if (actnode->c==NULL) dserror("Allocation of   COND_NODE failed"); 
      actnode->c->isneum=0;
      actnode->c->isdirich=0;
      actnode->c->iscoupled=0;
      actnode->c->neum_onoff.a.iv  =NULL;
      actnode->c->neum_val.a.iv    =NULL;
      actnode->c->dirich_onoff.a.iv=NULL;
      actnode->c->dirich_val.a.iv  =NULL;
      actnode->c->couple.a.iv      =NULL;
   }
/*--------------------------------------------------- read curve number */   
/* be carefull, at the moment only one curve possible (in GID 5!)       */
   frchar("CURVE",buffer,&ierr);
   if (ierr!=1) dserror("cannot read COND_FLUID_NODE from file");

   if (strncmp(buffer,"none",4)==0) actnode->c->curve=0;
   else
   {
      frint("CURVE",&(actnode->c->curve),&ierr);
      if (ierr!=1) dserror("cannot read COND_FLUID_NODE from file");
   }
/*----------------------------------------------------- read the values */   
   for (i=0; i<50; i++) { itmp[i]=0; dtmp[i]=0.0; }

   if (actnode->c->dirich_onoff.a.iv==NULL)
   amdef("dirich1",&(actnode->c->dirich_onoff),genprob.numdf,1,"IV");
   amzero(&(actnode->c->dirich_onoff));

   if (actnode->c->dirich_val.a.iv==NULL)
   amdef("dirich2",&(actnode->c->dirich_val),genprob.numdf,1,"DV");
   amzero(&(actnode->c->dirich_val));
   
   actnode->c->isdirich=1;
   frint_n("ONOFF",&(itmp[0]),4,&ierr);
   frdouble_n("VAL",&(dtmp[0]),4,&ierr);
   
   for (j=0; j<4; j++)
   {
      if (actnode->c->dirich_onoff.a.iv[j]==0) 
          actnode->c->dirich_onoff.a.iv[j] = itmp[j];
      if (FABS(actnode->c->dirich_val.a.dv[j])<eps) 
          actnode->c->dirich_val.a.dv[j] = dtmp[j];
   }
/*----------------------------------------------------------------------*/
   frread();
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_cond_nodal_fluid */





/*----------------------------------------------------------------------*
 | input of ale nodal conditions                          m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_cond_nodal_ale(FIELD *field)
{
int  ierr;
int  i,j;
int  nodeid;
char *colpointer;
char buffer[50];
int    itmp[50];
double dtmp[50];
double eps=1.0E-08;
NODE *actnode;
#ifdef DEBUG 
dstrc_enter("inp_cond_nodal_ale");
#endif
/*------------------------------------------ init the condition pointer */
for (i=0; i<field->dis[0].numnp; i++) field->dis[0].node[i].c==NULL;
/*------------------------------------------------- read the conditions */
frfind("--ALE NODAL CONDITIONS");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   nodeid = strtol(allfiles.actplace,&colpointer,10);
   nodeid--;
   for (j=0; j<field->dis[0].numnp; j++)
   {
      if (field->dis[0].node[j].Id == nodeid) 
      {
         actnode = &(field->dis[0].node[j]);
         break;
      }
   }
   if (actnode->c==NULL)
   {
      actnode->c=(COND_NODE*)CALLOC(1,sizeof(COND_NODE));
      if (actnode->c==NULL) dserror("Allocation of   COND_NODE failed"); 
      actnode->c->dirich_onoff.a.iv=NULL;
      actnode->c->dirich_val.a.iv  =NULL;
      actnode->c->couple.a.iv      =NULL;
   }
/*----------------------------------------------------- read the values */   
   for (i=0; i<50; i++) { itmp[i]=0; dtmp[i]=0.0; }
   if (actnode->c->dirich_onoff.a.iv==0)
   amdef("dirich1",&(actnode->c->dirich_onoff),genprob.numdf,1,"IV");
   amzero(&(actnode->c->dirich_onoff));
   if (actnode->c->dirich_val.a.iv==0)
   amdef("dirich2",&(actnode->c->dirich_val),genprob.numdf,1,"DV");
   amzero(&(actnode->c->dirich_val));

   actnode->c->isdirich=1;
   frint_n("ONOFF",&(itmp[0]),3,&ierr);
   frdouble_n("VAL",&(dtmp[0]),3,&ierr);

   for (j=0; j<3; j++)
   {
      if (actnode->c->dirich_onoff.a.iv[j]==0) 
          actnode->c->dirich_onoff.a.iv[j] = itmp[j];
      if (FABS(actnode->c->dirich_val.a.dv[j])<eps) 
          actnode->c->dirich_val.a.dv[j] = dtmp[j];
   }
   
   frread();
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_cond_nodal_ale */
