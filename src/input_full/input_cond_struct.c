#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 | input of structural nodal conditions                   m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_cond_nodal_struct(FIELD *field)
{
int    ierr;
int    i,j;
int    nodeid;
char  *colpointer;
char   buffer[50];
int    itmp[50];
double dtmp[50];
double eps=1.0E-08;
NODE  *actnode;
#ifdef DEBUG 
dstrc_enter("inp_cond_nodal_struct");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------ init the condition pointer */
for (i=0; i<field->numnp; i++) field->node[i].c==NULL;
/*------------------------------------------------- read the conditions */
frfind("--STRUCT NODAL CONDITIONS");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   nodeid = strtol(allfiles.actplace,&colpointer,10);
   nodeid--;
   for (j=0; j<field->numnp; j++)
   {
      if (field->node[j].Id == nodeid) 
      {
         actnode = &(field->node[j]);
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
      actnode->c->neum_onoff.a.iv=NULL;
      actnode->c->neum_val.a.iv=NULL;
      actnode->c->dirich_onoff.a.iv=NULL;
      actnode->c->dirich_val.a.iv=NULL;
      actnode->c->couple.a.iv=NULL;
   }
/*--------------------------------------------------- read curve number */   
/* be carefull, at the moment only one curve possible (in GID 5!)       */
   frchar("CURVE",buffer,&ierr);
   if (ierr!=1) dserror("cannot read COND_STRUCT_NODE from file");
   if (strncmp(buffer,"none",4)==0) actnode->c->curve=0;
   else
   {
      frint("CURVE",&(actnode->c->curve),&ierr);
      if (ierr!=1) dserror("cannot read COND_STRUCT_NODE from file");
   }
/*------------- read whether its a neum condition or a dirich condition */    
   for (i=0; i<50; i++) { itmp[i]=0; dtmp[i]=0.0; }
   frchk("Neum",&ierr);
   if (ierr==1) 
   {
      actnode->c->isneum=1;
      if (actnode->c->neum_onoff.a.iv==NULL)
      {
         amdef("neum1",&(actnode->c->neum_onoff),genprob.numdf,1,"IV");
         amzero(&(actnode->c->neum_onoff));
         amdef("neum2",&(actnode->c->neum_val)  ,genprob.numdf,1,"DV");
         amzero(&(actnode->c->neum_val));
      }
      frint_n("ONOFF" ,&(itmp[0]),genprob.numdf,&ierr);
      frdouble_n("VAL",&(dtmp[0]),genprob.numdf,&ierr);
      
      for (j=0; j<genprob.numdf; j++)
      {
         if (actnode->c->neum_onoff.a.iv[j]==0) 
             actnode->c->neum_onoff.a.iv[j] = itmp[j];
      }
      for (j=0; j<genprob.numdf; j++)
      {
         if (FABS(actnode->c->neum_val.a.dv[j])<eps)
             actnode->c->neum_val.a.dv[j]   = dtmp[j];
      }
   }
   frchk("Dirich",&ierr);
   if (ierr==1) 
   {
      actnode->c->isdirich=1;
      if (actnode->c->dirich_onoff.a.iv==NULL)
      {
         amdef("dirich1",&(actnode->c->dirich_onoff),genprob.numdf,1,"IV");
         amzero(&(actnode->c->dirich_onoff));
         amdef("dirich2",&(actnode->c->dirich_val)  ,genprob.numdf,1,"DV");
         amzero(&(actnode->c->dirich_val));
      }
      frint_n("ONOFF" ,&(itmp[0]),genprob.numdf,&ierr);
      frdouble_n("VAL",&(dtmp[0]),genprob.numdf,&ierr);
      
      for (j=0; j<genprob.numdf; j++)
      {
         if (actnode->c->dirich_onoff.a.iv[j]==0) 
             actnode->c->dirich_onoff.a.iv[j] = itmp[j];
      }
      for (j=0; j<genprob.numdf; j++)
      {
         if (FABS(actnode->c->dirich_val.a.dv[j])<eps)
             actnode->c->dirich_val.a.dv[j]   = dtmp[j];
      }
   }
   
   frread();
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_cond_nodal_struct */



/*----------------------------------------------------------------------*
 | input of structural element conditions                 m.gee 5/01    |
 *----------------------------------------------------------------------*/
void inp_cond_ele_struct(FIELD *field)
{
int    ierr;
int    i,j;
int    eleid;
char   *colpointer;
char   buffer[50];
int    ibuff[50];
double dbuff[50];
ELEMENT  *actele;
#ifdef DEBUG 
dstrc_enter("inp_cond_ele_struct");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------ init the condition pointer */
for (i=0; i<field->numele; i++) field->element[i].c==NULL;
/*------------------------------------------------------- start reading */
frfind("--STRUCTURAL ELEMENT CONDITIONS");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
/*------------------------ read element number and find correct element */
   eleid = strtol(allfiles.actplace,&colpointer,10);
   eleid--;
   for (j=0; j<field->numele; j++)
   {
      if (field->element[j].Id == eleid) 
      {
         actele = &(field->element[j]);
         break;
      }
   }
/*----------------------- allocate a condition structure to the element */   
   if (actele->c==NULL)
   {
      actele->c = (COND_ELEMENT*)calloc(1,sizeof(COND_ELEMENT));
      if (actele->c==0) dserror("Allocation of COND_ELEMENT failed");
      actele->c->neum_onoff.a.iv=NULL;
      actele->c->neum_val.a.iv=NULL;
   }
/*---------------------------------------------------------- read curve */
   frchar("CURVE",buffer,&ierr);
   if (strncmp(buffer,"none",4)==0) actele->c->curve=0;
   else
   {
      frint("CURVE",&(actele->c->curve),&ierr);
      if (ierr!=1) dserror("Reading of CURVE in COND_ELEMENT failed");
   }
/*----------------------------------------------- read typ of condition */
/*   frchar("TYP",actele->c->condtyp,&ierr);*/
   frchar("TYP",buffer,&ierr);
   if (ierr!=1) dserror("Reading of condtyp in COND_ELEMENT failed");
   if (strncmp(buffer,"Live",4)==0) actele->c->condtyp = ne_live;
   if (strncmp(buffer,"Dead",4)==0) actele->c->condtyp = ne_dead;
/*---------------------------------------------------------- read ONOFF */
   for (j=0; j<genprob.numdf; j++){ ibuff[j]=0; dbuff[j]=0.0;}
   frint_n("ONOFF",ibuff,genprob.numdf,&ierr);
   frdouble_n("VAL",dbuff,genprob.numdf,&ierr);
/* element condition at the moment can only be of type neumann, so check
   for any switched dofs and allocate ARRAY neum_onoff and neum_val */   
   if (actele->c->isneum!=1)
   {
      for (j=0; j<genprob.numdf; j++)
      {
         if (ibuff[j]==1) 
         {
            actele->c->isneum = 1;
            break;
         }
         else
         {
            actele->c->isneum = 0;
         }
      }
   }
   if (actele->c->isneum==1)
   {
      if (actele->c->neum_onoff.a.iv==NULL)
      {
         amdef("neum1",&(actele->c->neum_onoff),genprob.numdf,1,"IV");
         amzero(&(actele->c->neum_onoff));
         amdef("neum2",&(actele->c->neum_val),genprob.numdf,1,"DV");
         amzero(&(actele->c->neum_val));
      }
      for (j=0; j<genprob.numdf; j++)
      {
         if (ibuff[j]==1)
         {
            actele->c->neum_onoff.a.iv[j]=1;
            actele->c->neum_val.a.dv[j]+=dbuff[j];
         }
      }
   }
   frread();
} 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_cond_ele_struct */



