#include "../headers/standardtypes.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure which holds all file pointers                              |
 *----------------------------------------------------------------------*/
struct _FILES           allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 *----------------------------------------------------------------------*/
extern struct _DYNAMIC  *dyn;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate static variables if needed                       |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _STATIC_VAR  *statvar;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 | struct _FIELD         *field;                                        |
 *----------------------------------------------------------------------*/
extern struct _FIELD       *field;      



/*----------------------------------------------------------------------*
 | input of control information                           m.gee 8/00    |
 *----------------------------------------------------------------------*/
void inpctr()
{
#ifdef DEBUG 
dstrc_enter("inpctr");
#endif

/*---------------------------------------------- input of problem title */
   inpctrhed();
/*--------------------------------------- input of general problem data */
   inpctrprob(); 
/*---------------------------------------- input of general solver data */
   if (genprob.probtyp == prb_fsi)
   {
      if (genprob.numfld!=3) dserror("numfld != 3 for FSI");
      
      solv = (SOLVAR*)CALLOC(genprob.numfld,sizeof(SOLVAR));
      if (!solv) dserror("Allocation of SOLVAR failed");
      
      solv[0].fieldtyp = structure;
      inpctrsol(&(solv[0]));

      solv[1].fieldtyp = fluid;
      inpctrsol(&(solv[1]));

      solv[2].fieldtyp = ale;
      inpctrsol(&(solv[2]));
   }
   if (genprob.probtyp == prb_structure)
   {
      if (genprob.numfld!=1) dserror("numfld != 1 for Structural Problem");
      
      solv = (SOLVAR*)CALLOC(genprob.numfld,sizeof(SOLVAR));
      if (!solv) dserror("Allocation of SOLVAR failed");
      
      solv[0].fieldtyp = structure;
      inpctrsol(&(solv[0]));
   }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpctr */




/*----------------------------------------------------------------------*
 | input of general problem data                          m.gee 2/01    |
 *----------------------------------------------------------------------*/
void inpctrprob()
{
int  ierr;
int  i;
char buffer[50];
#ifdef DEBUG 
dstrc_enter("inpctrprob");
#endif

frfind("-PROBLEM SIZE");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   frint("ELEMENTS", &(genprob.nele),&ierr);
   frint("NODES",    &(genprob.nnode),&ierr);
   frint("DIM",      &(genprob.ndim),&ierr);
   frint("MATERIALS",&(genprob.nmat),&ierr);
   frint("NUMDF",    &(genprob.numdf),&ierr);

   frread();
}
frrewind();


frfind("-PROBLEM TYP");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   frchar("PROBLEMTYP",buffer,   &ierr);
   if (ierr==1)
   {
      if (strncmp("Structure"                  ,buffer, 9)==0) genprob.probtyp = prb_structure;
      if (strncmp("Fluid"                      ,buffer, 5)==0) genprob.probtyp = prb_fluid;
      if (strncmp("Fluid_Structure_Interaction",buffer,24)==0) genprob.probtyp = prb_fsi;
      if (strncmp("Optimisation"               ,buffer,12)==0) genprob.probtyp = prb_opt;
   }

   frchar("TIMETYP"   ,buffer,            &ierr);
   if (ierr==1)
   {
      if (strncmp("Static" ,buffer,6)==0) genprob.timetyp=time_static;
      if (strncmp("Dynamic",buffer,7)==0) genprob.timetyp=time_dynamic;
   }

   frchar("DESIGN"    ,buffer,            &ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"on",2)==0) genprob.design=1;
      else                           genprob.design=0;
   }
   frint("RESTART"    ,&(genprob.restart),&ierr);
   
   frint("NUMFIELD",&(genprob.numfld),&ierr);

   frread();
}
frrewind();


frfind("---IO");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   frchar("STRUCT_DISP_FILE",buffer,&ierr);
   if (ierr)
   {
      if (strncmp(buffer,"Yes",3)==0) ioflags.struct_disp_file=1;
      else                            ioflags.struct_disp_file=0;
   }
   frchar("STRUCT_STRESS_FILE",buffer,&ierr);
   if (ierr)
   {
      if (strncmp(buffer,"Yes",3)==0) ioflags.struct_stress_file=1;
      else                            ioflags.struct_stress_file=0;
   }
   frchar("STRUCT_DISP_GID",buffer,&ierr);
   if (ierr)
   {
      if (strncmp(buffer,"Yes",3)==0) ioflags.struct_disp_gid=1;
      else                            ioflags.struct_disp_gid=0;
   }
   frchar("STRUCT_STRESS_GID",buffer,&ierr);
   if (ierr)
   {
      if (strncmp(buffer,"Yes",3)==0) ioflags.struct_stress_gid=1;
      else                            ioflags.struct_stress_gid=0;
   }
   

   frread();
}
frrewind();


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpctrprob */








/*----------------------------------------------------------------------*
 | input of static problem data                           m.gee 2/01    |
 *----------------------------------------------------------------------*/
void inpctrstat()
{
int  ierr;
int  i;
char buffer[50];
#ifdef DEBUG 
dstrc_enter("inpctrstat");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------- allocate a structure STATIC */
statvar = (STATIC_VAR*)CALLOC(1,sizeof(STATIC_VAR));
if (!statvar) dserror("Allocation of STATIC failed");
/*------------------------------------------------------- start reading */
frfind("-STATIC");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   frchk("GEOLINEAR",&ierr);
   if (ierr==1) {statvar->geolinear=1;statvar->geononlinear=0;}
   frchk("GEONONLINEAR",&ierr);
   if (ierr==1) {statvar->geononlinear=1;statvar->geolinear=0;}
   /*------------------------- read for typ of pathfollowing technique */   
   frchar("NEWTONRAPHSO",buffer,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"Displacement_Control",20)==0) 
      statvar->nr_controltyp = control_disp;
      if (strncmp(buffer,"Load_Control",12)==0)
      statvar->nr_controltyp = control_load;
      if (strncmp(buffer,"Arc_Control",11)==0)
      statvar->nr_controltyp = control_arc;
      if (strncmp(buffer,"none",4)==0)
      statvar->nr_controltyp = control_none;
   }
   /*---------------------------------------- read for arcscaling flag */
   frchar("IARC"   ,buffer,&ierr);
   if (ierr)
   {
      if (strncmp(buffer,"Yes",3)==0) statvar->iarc=1;
      else                            statvar->iarc=0;
   }  
   /*-------------------------------------------------- read variables */
   frint("NUMSTEP",&(statvar->nstep)  ,&ierr);
   frint("MAXITER",&(statvar->maxiter),&ierr);

   frdouble("TOLRESID",&(statvar->tolresid),&ierr);
   frdouble("TOLDISP" ,&(statvar->toldisp) ,&ierr);
   frdouble("STEPSIZE",&(statvar->stepsize),&ierr);
   frdouble("ARCSCL"  ,&(statvar->arcscl)  ,&ierr);
   /*-------------------------------------------------------------------*/
   frread();
}
frrewind();
/*------------------------------ in geononlinear case find control node */
if (statvar->geononlinear)
{
   /*---------------------------------------------------- start reading */
   frfind("-CONTROL NODE");
   frread();
   while(strncmp(allfiles.actplace,"------",6)!=0)
   {
      frint("NODE",&(statvar->control_node_global),&ierr);
      if (ierr) statvar->control_node_global--;
      frint("DOF",&(statvar->control_dof),&ierr);
      /*----------------------------------------------------------------*/
      frread();
   }
frrewind();
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpctrstat */




/*----------------------------------------------------------------------*
 | input of dynamic problem data                          m.gee 2/01    |
 *----------------------------------------------------------------------*/
void inpctrdyn()
{
int    ierr;
int    i;
char   buffer[50];
FIELD *actfield;
#ifdef DEBUG 
dstrc_enter("inpctrdyn");
#endif
/*----------------------------------------------------------------------*/
dyn = (DYNAMIC*)CALLOC(genprob.numfld,sizeof(DYNAMIC));
if (!dyn) dserror("Allocation of DYNAMIC failed");
/*----------------------------------------------------------------------*/
for (i=0; i<genprob.numfld; i++)
{
   actfield = &(field[i]);
   switch(actfield->fieldtyp)
   {
   case fluid:
   break;
   case ale:
   break;
   case structure:
      dyn[i].sdyn = (STRUCT_DYNAMIC*)CALLOC(1,sizeof(STRUCT_DYNAMIC));
      if (!dyn[i].sdyn) dserror("Allocation of DYNAMIC failed");
      inpctr_dyn_struct(dyn[i].sdyn);
   break;
   case none:
      dserror("Cannot find type of field");
   break;
   default:
      dserror("Cannot find type of field");
   break;
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpctrdyn */



/*----------------------------------------------------------------------*
 | input of dynamic problem data  for field structure     m.gee 2/02    |
 *----------------------------------------------------------------------*/
int inpctr_dyn_struct(STRUCT_DYNAMIC *sdyn)
{
int    ierr;
int    i;
char   buffer[50];
#ifdef DEBUG 
dstrc_enter("inpctr_dyn_struct");
#endif

frfind("-STRUCTURAL DYNAMIC");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
/*--------------read chars */
   frchar("DYNAMICTYP",sdyn->dyntyp,&ierr);
   frchar("DAMPING"   ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"Yes",3)==0 ||
          strncmp(buffer,"yes",3)==0 ||
          strncmp(buffer,"YES",3)==0 )
          
          sdyn->damp=1;
      else
          sdyn->damp=0;
   }
   frchar("ITERATION"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"full",3)==0 ||
          strncmp(buffer,"Full",3)==0 ||
          strncmp(buffer,"FULL",3)==0 )
          
          sdyn->iter=1;
      else
          sdyn->iter=0;
   }

/*--------------read int */
   frint("NUMSTEP",&(sdyn->nstep)  ,&ierr);
   frint("MAXITER",&(sdyn->maxiter),&ierr);
   
/*--------------read double */
   frdouble("TIMESTEP",&(sdyn->dt)     ,&ierr);
   frdouble("MAXTIME" ,&(sdyn->maxtime),&ierr);
   frdouble("BETA"    ,&(sdyn->beta)   ,&ierr);
   frdouble("GAMMA"   ,&(sdyn->gamma)  ,&ierr);
   frdouble("ALPHA_M" ,&(sdyn->alpha_m),&ierr);
   frdouble("ALPHA_F" ,&(sdyn->alpha_f),&ierr);
   frdouble("M_DAMP"  ,&(sdyn->m_damp),&ierr);
   frdouble("K_DAMP"  ,&(sdyn->k_damp),&ierr);

   frread();
}
frrewind();
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpctr_dyn_struct */
