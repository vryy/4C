/*!----------------------------------------------------------------------
\file
\brief input of control information and general problem data

------------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
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
 | ALLDYNA       *alldyn;                                               |
 *----------------------------------------------------------------------*/
extern ALLDYNA        *alldyn;   
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
      
      solv = (SOLVAR*)CCACALLOC(genprob.numfld,sizeof(SOLVAR));
      if (!solv) dserror("Allocation of SOLVAR failed");
      
      solv[genprob.numsf].fieldtyp = structure;
      inpctrsol(&(solv[genprob.numsf]));

      solv[genprob.numff].fieldtyp = fluid;
      inpctrsol(&(solv[genprob.numff]));

      solv[genprob.numaf].fieldtyp = ale;
      inpctrsol(&(solv[genprob.numaf]));
   }
   if (genprob.probtyp == prb_structure)
   {
      if (genprob.numfld!=1) dserror("numfld != 1 for Structural Problem");
      
      solv = (SOLVAR*)CCACALLOC(genprob.numfld,sizeof(SOLVAR));
      if (!solv) dserror("Allocation of SOLVAR failed");
      
      solv[genprob.numsf].fieldtyp = structure;
      inpctrsol(&(solv[genprob.numsf]));
   }
   if (genprob.probtyp == prb_opt)
   {
      if (genprob.numfld!=1) dserror("numfld != 1 for Structural Problem");
      
      solv = (SOLVAR*)CCACALLOC(genprob.numfld,sizeof(SOLVAR));
      if (!solv) dserror("Allocation of SOLVAR failed");
      
      solv[0].fieldtyp = structure;
      inpctrsol(&(solv[0]));
   }
   if (genprob.probtyp == prb_fluid)   
   {
      solv = (SOLVAR*)CCACALLOC(genprob.numfld,sizeof(SOLVAR));
      if (!solv) dserror("Allocation of SOLVAR failed");
      
      solv[genprob.numff].fieldtyp = fluid;
      inpctrsol(&(solv[genprob.numff]));
      
      if (genprob.numfld==2)
      {
        solv[genprob.numaf].fieldtyp = ale; 
	inpctrsol(&(solv[genprob.numaf]));
      }
                
   }
   if (genprob.probtyp == prb_ale)
   {
      if (genprob.numfld!=1) dserror("numfld != 1 for Ale Problem");

      solv = (SOLVAR*)CCACALLOC(genprob.numfld,sizeof(SOLVAR));
      if (!solv) dserror("Allocation of SOLVAR failed");

      solv[genprob.numaf].fieldtyp = ale;
      inpctrsol(&(solv[genprob.numaf]));
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

/*------------------------------------------ set initial field numbers */
genprob.numsf=-1;
genprob.numff=-1;
genprob.numaf=-1;

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
/*------------------------------------------------------ check values */
if (genprob.nmat<=0)
   dserror ("No Material defined!");
/*------------------------------  default value for multidis flag = 0 */
genprob.multidis=0;
/*-------------------------------------- default value for monitoring */
ioflags.monitor=0;
/*----------------------------------------- defualt values for output */
ioflags.struct_disp_file   =0;
ioflags.struct_stress_file =0;
ioflags.struct_disp_gid    =0;
ioflags.struct_stress_gid  =0;
ioflags.fluid_sol_gid      =0;
ioflags.fluid_sol_file     =0; 
ioflags.fluid_vis_file     =0;
ioflags.ale_disp_file      =0;
ioflags.ale_disp_gid       =0;

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
      if (strncmp("Ale"                        ,buffer, 3)==0) genprob.probtyp = prb_ale;
   }

   frchar("TIMETYP"   ,buffer,            &ierr);
   if (ierr==1)
   {
      if (strncmp("Static" ,buffer,6)==0) genprob.timetyp=time_static;
      if (strncmp("Dynamic",buffer,7)==0) genprob.timetyp=time_dynamic;
   }

   frchar("MULTIDIS"   ,buffer,            &ierr);
   if (ierr==1)
   {
      if (strncmp("yes" ,buffer,3)==0) genprob.multidis=1;
      if (strncmp("YES" ,buffer,3)==0) genprob.multidis=1;
      if (strncmp("Yes" ,buffer,3)==0) genprob.multidis=1;
   }      

   frint("RESTART"    ,&(genprob.restart),&ierr);
   
   frint("NUMFIELD",&(genprob.numfld),&ierr);

   frread();
}
frrewind();

/*----------------- set field numbers depending on problem type and numfld */
if (genprob.probtyp==prb_fsi)
{
   genprob.numsf=0;
   genprob.numff=1;
   genprob.numaf=2;
}
if (genprob.probtyp==prb_fluid)
{
   genprob.numff=0;
   if (genprob.numfld==2) genprob.numaf=1;
}
if (genprob.probtyp==prb_ale) genprob.numaf=0;
if (genprob.probtyp==prb_structure) genprob.numsf=0;


frfind("---IO");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   frchar("STRUCT_DISP_FILE",buffer,&ierr);
   if (ierr)
   {
      if (strncmp(buffer,"yes",3)==0) ioflags.struct_disp_file=1;
      if (strncmp(buffer,"YES",3)==0) ioflags.struct_disp_file=1;
      if (strncmp(buffer,"Yes",3)==0) ioflags.struct_disp_file=1;
   }
   frchar("STRUCT_STRESS_FILE",buffer,&ierr);
   if (ierr)
   {
      if (strncmp(buffer,"yes",3)==0) ioflags.struct_stress_file=1;
      if (strncmp(buffer,"YES",3)==0) ioflags.struct_stress_file=1;
      if (strncmp(buffer,"Yes",3)==0) ioflags.struct_stress_file=1;
   }
   frchar("STRUCT_DISP_GID",buffer,&ierr);
   if (ierr)
   {
      if (strncmp(buffer,"yes",3)==0) ioflags.struct_disp_gid=1;
      if (strncmp(buffer,"YES",3)==0) ioflags.struct_disp_gid=1;
      if (strncmp(buffer,"Yes",3)==0) ioflags.struct_disp_gid=1;
   }
   frchar("STRUCT_STRESS_GID",buffer,&ierr);
   if (ierr)
   {
      if (strncmp(buffer,"yes",3)==0) ioflags.struct_disp_file=1;
      if (strncmp(buffer,"YES",3)==0) ioflags.struct_disp_file=1;
      if (strncmp(buffer,"Yes",3)==0) ioflags.struct_disp_file=1;
   }
   frchar("FLUID_SOL_GID",buffer,&ierr);
   if (ierr)
   {
      if (strncmp(buffer,"yes",3)==0) ioflags.fluid_sol_gid=1;
      if (strncmp(buffer,"YES",3)==0) ioflags.fluid_sol_gid=1;
      if (strncmp(buffer,"Yes",3)==0) ioflags.fluid_sol_gid=1;
   }     
   frchar("FLUID_SOL_FILE",buffer,&ierr);
   if (ierr)
   {
      if (strncmp(buffer,"yes",3)==0) ioflags.fluid_sol_file=1;
      if (strncmp(buffer,"YES",3)==0) ioflags.fluid_sol_file=1;
      if (strncmp(buffer,"Yes",3)==0) ioflags.fluid_sol_file=1;
   }     
   frchar("FLUID_VIS_FILE",buffer,&ierr);
   if (ierr)
   {
      if (strncmp(buffer,"yes",3)==0) ioflags.fluid_vis_file=1;
      if (strncmp(buffer,"YES",3)==0) ioflags.fluid_vis_file=1;
      if (strncmp(buffer,"Yes",3)==0) ioflags.fluid_vis_file=1;
   } 
   frchar("ALE_DISP_FILE",buffer,&ierr);
   if (ierr)
   {
      if (strncmp(buffer,"yes",3)==0) ioflags.ale_disp_file=1;
      if (strncmp(buffer,"YES",3)==0) ioflags.ale_disp_file=1;
      if (strncmp(buffer,"Yes",3)==0) ioflags.ale_disp_file=1;
   }
   frchar("ALE_DISP_GID",buffer,&ierr);
   if (ierr)
   {
      if (strncmp(buffer,"yes",3)==0) ioflags.ale_disp_gid=1;
      if (strncmp(buffer,"YES",3)==0) ioflags.ale_disp_gid=1;
      if (strncmp(buffer,"Yes",3)==0) ioflags.ale_disp_gid=1;
   }
   frchar("MONITOR"   ,buffer,            &ierr);
   if (ierr==1)
   {
      if (strncmp("yes" ,buffer,3)==0) ioflags.monitor=1;
      if (strncmp("YES" ,buffer,3)==0) ioflags.monitor=1;
      if (strncmp("Yes" ,buffer,3)==0) ioflags.monitor=1;
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
statvar = (STATIC_VAR*)CCACALLOC(1,sizeof(STATIC_VAR));
if (!statvar) dserror("Allocation of STATIC failed");
/*----------------------------- default: results written every loadstep */
statvar->resevry_disp=1;
statvar->resevry_stress=1;
statvar->resevery_restart=1;
/*------------------------------------------------------- start reading */
frfind("-STATIC");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   frchk("LINEAR",&ierr);
   if (ierr==1) {statvar->linear=1;statvar->nonlinear=0;}
   frchk("NONLINEAR",&ierr);
   if (ierr==1) {statvar->nonlinear=1;statvar->linear=0;}
   /*------------------------- read for typ of kinematic sh 03/03 */   
   frchar("KINTYP",buffer,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"Geo_Lin",7)==0) 
      statvar->kintyp = geo_linear;
      if (strncmp(buffer,"Upd_Lagr",8)==0)
      statvar->kintyp = upd_lagr;
      if (strncmp(buffer,"Tot_Lagr",8)==0)
      statvar->kintyp = tot_lagr;
   }
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
   /*------------------------------ read for sign-changing-by-csp flag */
   frchar("SIGNCHCSP"   ,buffer,&ierr);
   if (ierr)
   {
      if (strncmp(buffer,"Yes",3)==0) statvar->signchcsp=1;
      else                            statvar->signchcsp=0;
   }  
   /*-------------------------------------------------- read variables */
   frint("NUMSTEP",&(statvar->nstep)  ,&ierr);
   frint("MAXITER",&(statvar->maxiter),&ierr);
   frint("RESEVRYDISP",&(statvar->resevry_disp)  ,&ierr);
   frint("RESEVRYSTRS",&(statvar->resevry_stress),&ierr);
   frint("RESTARTEVRY",&(statvar->resevery_restart),&ierr);

   frdouble("TOLRESID",&(statvar->tolresid),&ierr);
   frdouble("TOLDISP" ,&(statvar->toldisp) ,&ierr);
   frdouble("STEPSIZE",&(statvar->stepsize),&ierr);
   frdouble("ARCSCL"  ,&(statvar->arcscl)  ,&ierr);
   /*-------------------------------------------------------------------*/
   frread();
}
frrewind();
/*--------------------------------- in nonlinear case find control node */
if (statvar->nonlinear)
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
 |                                                         genk 3/01    |
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
if (genprob.probtyp==prb_fsi || 
   (genprob.probtyp==prb_fluid &&  genprob.numfld>1))
   alldyn = (ALLDYNA*)CCACALLOC((genprob.numfld+1),sizeof(ALLDYNA));
else
   alldyn = (ALLDYNA*)CCACALLOC(genprob.numfld,sizeof(ALLDYNA));
if (!alldyn) dserror("Allocation of ALLDYNA failed");
/*----------------------------------------------------------------------*/
for (i=0; i<genprob.numfld; i++)
{
   actfield = &(field[i]);
   switch(actfield->fieldtyp)
   {
   case fluid:
#ifdef D_FLUID   
      alldyn[i].fdyn = (FLUID_DYNAMIC*)CCACALLOC(1,sizeof(FLUID_DYNAMIC));
      if (!alldyn[i].fdyn) dserror("Allocation of FLUID_DYNAMIC failed");
      inpctr_dyn_fluid(alldyn[i].fdyn);
#else
      dserror("General FLUID problem not defined in Makefile!!!");
#endif            
   break;
   case ale:
      alldyn[i].sdyn = (STRUCT_DYNAMIC*)CCACALLOC(1,sizeof(STRUCT_DYNAMIC));
      if (!alldyn[i].sdyn) dserror("Allocation of STRUCT_DYNAMIC failed");
      inpctr_dyn_struct(alldyn[i].sdyn);
   break;
   case structure:
      alldyn[i].sdyn = (STRUCT_DYNAMIC*)CCACALLOC(1,sizeof(STRUCT_DYNAMIC));
      if (!alldyn[i].sdyn) dserror("Allocation of STRUCT_DYNAMIC failed");
      inpctr_dyn_struct(alldyn[i].sdyn);
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
if (genprob.probtyp==prb_fsi || 
   (genprob.probtyp==prb_fluid &&  genprob.numfld>1))
{
#ifdef D_FSI
   alldyn[i].fsidyn = (FSI_DYNAMIC*)CCACALLOC(1,sizeof(FSI_DYNAMIC));
   inpctr_dyn_fsi(alldyn[i].fsidyn);
#else
   dserror("General FSI problem not defined in Makefile!!!");
#endif 
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
void inpctr_dyn_struct(STRUCT_DYNAMIC *sdyn)
{
int    ierr;
int    i;
char   buffer[50];
#ifdef DEBUG 
dstrc_enter("inpctr_dyn_struct");
#endif

sdyn->updevry_disp=1;
sdyn->updevry_stress=1;
sdyn->res_write_evry=1;
sdyn->eigen=0;

frfind("-STRUCTURAL DYNAMIC");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
/*--------------read chars */
   frchar("DYNAMICTYP",buffer,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"Centr_Diff",10)==0) sdyn->Typ = centr_diff;
      if (strncmp(buffer,"Gen_Alfa",8)==0)    sdyn->Typ = gen_alfa;
   }
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

/*-----------------read int */
   frint("EIGEN"      ,&(sdyn->eigen)         ,&ierr);
   frint("NUMSTEP"    ,&(sdyn->nstep)         ,&ierr);
   frint("MAXITER"    ,&(sdyn->maxiter)       ,&ierr);
   frint("RESEVRYDISP",&(sdyn->updevry_disp)  ,&ierr);
   frint("RESEVRYSTRS",&(sdyn->updevry_stress),&ierr);
   frint("RESTARTEVRY",&(sdyn->res_write_evry),&ierr);
   frint("CONTACT"    ,&(sdyn->contact)       ,&ierr);
   
/*--------------read double */
   frdouble("TIMESTEP",&(sdyn->dt)     ,&ierr);
   frdouble("MAXTIME" ,&(sdyn->maxtime),&ierr);
   frdouble("BETA"    ,&(sdyn->beta)   ,&ierr);
   frdouble("GAMMA"   ,&(sdyn->gamma)  ,&ierr);
   frdouble("ALPHA_M" ,&(sdyn->alpha_m),&ierr);
   frdouble("ALPHA_F" ,&(sdyn->alpha_f),&ierr);
   frdouble("M_DAMP"  ,&(sdyn->m_damp) ,&ierr);
   frdouble("K_DAMP"  ,&(sdyn->k_damp) ,&ierr);
   frdouble("TOLDISP" ,&(sdyn->toldisp),&ierr);

/*------ read time adaption */
   frint("TIMEADAPT"  ,&(sdyn->timeadapt),&ierr);
   frint("ITWANT"     ,&(sdyn->itwant)   ,&ierr);
   frdouble("MAXDT"   ,&(sdyn->maxdt)    ,&ierr);
   frdouble("RESULTDT",&(sdyn->resultdt) ,&ierr);

   frread();
}
frrewind();
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpctr_dyn_struct */

#ifdef D_FLUID

/*!--------------------------------------------------------------------- 
\brief input of the FLUID DYNAMIC block in the input-file

<pre>                                                         genk 03/02

In this routine the data in the FLUID DYNAMIC block of the input file
are read and stored in fdyn	       

</pre>
\param  *fdyn 	  FLUID_DATA       (o)	   
\return void                                                                       

------------------------------------------------------------------------*/
void inpctr_dyn_fluid(FLUID_DYNAMIC *fdyn)
{
int    ierr;
int    i;
int    thetafound=0;
char   buffer[50];
#ifdef DEBUG 
dstrc_enter("inpctr_dyn_fluid");
#endif

/*-------------------------------------------------- set default values */
fdyn->iop=2;
fdyn->freesurf=0;
fdyn->surftens=0;
fdyn->iops=2;
fdyn->ite=1;
fdyn->itchk=1;	 
fdyn->itnorm=2;
fdyn->stchk=0;
fdyn->init=0;
fdyn->viscstr=0;
fdyn->numdf=3;  
fdyn->numcont=0;
fdyn->uppss=1;  
fdyn->upout=1;  
fdyn->nstep=1;
fdyn->stchk=5;  
fdyn->nums=0;
fdyn->init=-1;
fdyn->iprerhs=1;
fdyn->itemax=3;
fdyn->dt=0.01;    
fdyn->maxtime=1000.0;
fdyn->alpha=1.0;
fdyn->gamma=1.0; 
fdyn->theta=0.5;
fdyn->ittol=EPS6; 
fdyn->sttol=EPS6; 
 
frfind("-FLUID DYNAMIC");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
/*--------------read chars */
   frchar("DYNAMICTYP",fdyn->dyntyp,&ierr);
   frchar("TIMEINTEGR"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"Stationary",10)==0) 
         fdyn->iop=0;
      else if (strncmp(buffer,"Semi_Impl_One_Step",18)==0) 
         fdyn->iop=2;                  
      else if (strncmp(buffer,"Semi_Impl_Two_Step",18)==0) 
         fdyn->iop=3;
      else if (strncmp(buffer,"One_Step_Theta",14)==0) 
         fdyn->iop=4;
      else if (strncmp(buffer,"Fract_Step_Theta",16)==0) 
         fdyn->iop=5;      
      else
         dserror("TIMEINTEGR-Type unknown");
   }
   frchar("STARTINGALGO"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"Semi_Impl_One_Step",18)==0) 
         fdyn->iops=2;                  
      else if (strncmp(buffer,"Semi_Impl_Two_Step",18)==0) 
         fdyn->iops=3;
      else if (strncmp(buffer,"One_Step_Theta",14)==0) 
         fdyn->iops=4;
      else if (strncmp(buffer,"Fract_Step_Theta",16)==0) 
         fdyn->iops=5;      
      else
         dserror("STARTINGALGO unknown");
   }   
   frchar("NONLINITER"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"No",2)==0 ||
          strncmp(buffer,"NO",2)==0 ||
	  strncmp(buffer,"no",2)==0   ) 
         fdyn->ite=0;
      else if (strncmp(buffer,"fixed_point_like",16)==0) 
         fdyn->ite=1;
      else if (strncmp(buffer,"Newton",6)==0) 
         fdyn->ite=2;                  
      else if (strncmp(buffer,"fixed_point",11)==0) 
         fdyn->ite=3;
      else
         dserror("NONLINITER-Type unknown!");     
   }
   frchar("CONVCHECK"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"No",2)==0 ||
          strncmp(buffer,"NO",2)==0 ||
	  strncmp(buffer,"no",2)==0    )
      {  
         fdyn->itchk=0;
	 fdyn->itnorm=-1;      	  
      }
      else
      {   
	 fdyn->itchk=1;	       
         if (strncmp(buffer,"L_infinity_norm",15)==0) 
            fdyn->itnorm=0;
         else if (strncmp(buffer,"L_1_norm",8)==0) 
            fdyn->itnorm=1;                  
         else if (strncmp(buffer,"L_2_norm",8)==0) 
            fdyn->itnorm=2;
         else
            dserror("Norm for CONVCHECK unknown");   
      }  
   }   
   frchar("STEADYCHECK"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"No",2)==0 ||
          strncmp(buffer,"NO",2)==0 ||
	  strncmp(buffer,"no",2)==0    )
      { 
         fdyn->stchk=0;
	 fdyn->stnorm=-1;	  
      }
      else
      {
         fdyn->stchk=-1;	       
         if (strncmp(buffer,"L_infinity_norm",15)==0) 
            fdyn->stnorm=0;
         else if (strncmp(buffer,"L_1_norm",8)==0) 
            fdyn->stnorm=1;                  
         else if (strncmp(buffer,"L_2_norm",8)==0) 
            fdyn->stnorm=2;
         else 	  
	    dserror("Norm for CONVCHECK unknown");   
      }
   }
   frchar("INITIALFIELD"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"zero_field",10)==0) 
         fdyn->init=0;
      else if (strncmp(buffer,"field_from_file",15)==0) 
         fdyn->init=1;
      else if (strncmp(buffer,"field_by_function",17)==0) 
         fdyn->init=-1;                      
      else
         dserror("INITIALFIELD unknown!");
   }   
   frchar("VISCSTRESS"   ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"yes",3)==0 ||
          strncmp(buffer,"YES",3)==0 ||                 
	  strncmp(buffer,"Yes",3)==0    )
         fdyn->viscstr=1;
      else if (strncmp(buffer,"No",2)==0 ||
               strncmp(buffer,"NO",2)==0 ||
	       strncmp(buffer,"no",2)==0   )
         fdyn->viscstr=0;
      else
         dserror("VISCSTRESS unknown!\n");	 
   }
   frchar("FREESURFACE"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"No",2)==0 ||
          strncmp(buffer,"NO",2)==0 ||
	  strncmp(buffer,"no",2)==0   ) 
         fdyn->freesurf=0;
      else if (strncmp(buffer,"explicit",8)==0) 
         fdyn->freesurf=1;
      else if (strncmp(buffer,"implicit",8)==0) 
         fdyn->freesurf=2;                  
      else
         dserror("Parameter FREESURFACE unknown!\n");     
   }
   frchar("SURFTENSION"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"yes",3)==0 ||
          strncmp(buffer,"YES",3)==0 ||                 
	  strncmp(buffer,"Yes",3)==0    )
         fdyn->surftens=1;
      else if (strncmp(buffer,"No",2)==0 ||
               strncmp(buffer,"NO",2)==0 ||
	       strncmp(buffer,"no",2)==0   )
         fdyn->surftens=0;
      else
         dserror("SURFTENSION unknown!\n");	 
   }
   frchar("CHECKAREA"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"yes",3)==0 ||
          strncmp(buffer,"YES",3)==0 ||                 
	  strncmp(buffer,"Yes",3)==0    )
         fdyn->checkarea=1;
      else if (strncmp(buffer,"No",2)==0 ||
               strncmp(buffer,"NO",2)==0 ||
	       strncmp(buffer,"no",2)==0   )
         fdyn->checkarea=0;
      else
         dserror("CHECKAREA unknown!\n");	 
   }

/*--------------read int */
   frint("NUMDF"     ,&(fdyn->numdf)  ,&ierr);
   frint("NUMCONT"   ,&(fdyn->numcont),&ierr);
   frint("UPPSS"     ,&(fdyn->uppss)  ,&ierr);
   frint("UPOUT"     ,&(fdyn->upout)  ,&ierr);
   frint("NUMSTEP"   ,&(fdyn->nstep)  ,&ierr);
   frint("STEADYSTEP",&i              ,&ierr);   
   if (ierr==1)
   {
      if (fdyn->stchk==-1)
          fdyn->stchk=i;
   }
   frint("NUMSTASTEPS"  ,&(fdyn->nums),&ierr);
   frint("STARTFUNCNO"  ,&i             ,&ierr);
   if (ierr==1)
   {
      if (fdyn->init==-1)
          fdyn->init=i;
   }
   frint("IPRERHS",&(fdyn->iprerhs),&ierr); 
   frint("ITEMAX" ,&(fdyn->itemax) ,&ierr);  
/*--------------read double */
   frdouble("TIMESTEP" ,&(fdyn->dt)     ,&ierr);
   frdouble("MAXTIME"  ,&(fdyn->maxtime),&ierr);
   frdouble("ALPHA"    ,&(fdyn->alpha)  ,&ierr);
   frdouble("GAMMA"    ,&(fdyn->gamma)  ,&ierr);
   if (thetafound==0)
   {
      frdouble("THETA"    ,&(fdyn->theta)  ,&ierr);
      if (ierr==1) thetafound++;
   }
   frdouble("CONVTOL"  ,&(fdyn->ittol)  ,&ierr);
   frdouble("STEADYTOL",&(fdyn->sttol) ,&ierr);
   frdouble("START_THETA",&(fdyn->thetas),&ierr);

   frread();
}
frrewind();
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpctr_dyn_fluid */

#endif

#ifdef D_FSI

/*!--------------------------------------------------------------------- 
\brief input of the FSI DYNAMIC block in the input-file

<pre>                                                         genk 09/02

In this routine the data in the FSI DYNAMIC block of the input file
are read and stored in fsidyn	       

</pre>
\param  *fsidyn 	  FSI_DATA       (o)	   
\return void                                                                       

------------------------------------------------------------------------*/
void inpctr_dyn_fsi(FSI_DYNAMIC *fsidyn)
{
int    ierr;
int    i;
int    length;
char   buffer[50];
#ifdef DEBUG 
dstrc_enter("inpctr_dyn_fsi");
#endif

/*-------------------------------------------------- set default values */
fsidyn->ifsi=1;
fsidyn->ipre=1;
fsidyn->inrmfsi=1;
fsidyn->ichecke=0;
fsidyn->inest=0;
fsidyn->ichopt=0;
fsidyn->iait=0;
fsidyn->itechapp=0;
fsidyn->ichmax=1;
fsidyn->isdmax=1;  
fsidyn->nstep=1;   
fsidyn->itemax=1;  
fsidyn->iale=1;    
fsidyn->uppss=1;
fsidyn->dt=0.1;	
fsidyn->maxtime=1000.0;
fsidyn->entol=EPS6;  
fsidyn->relax=1.0;  
fsidyn->convtol=EPS6;
    
frfind("-FSI DYNAMIC");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
/*--------------read chars */
   frchar("COUPALGO"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      length = strlen(buffer);
      if (strncmp(buffer,"basic_sequ_stagg",16)==0 && length==16) 
         fsidyn->ifsi=1;
      else if (strncmp(buffer,"sequ_stagg_pred",15)==0 && length==15) 
         fsidyn->ifsi=2;                 
      else if (strncmp(buffer,"sequ_stagg_shift",15)==0 && length==15)
         fsidyn->ifsi=3;  
      else if (strncmp(buffer,"iter_stagg_fixed_rel_param",26)==0 && length==26) 
         fsidyn->ifsi=4;
      else if (strncmp(buffer,"iter_stagg_AITKEN_rel_param",27)==0 && length==27) 
         fsidyn->ifsi=5;      
      else if (strncmp(buffer,"iter_stagg_steep_desc",21)==0 && length==21) 
         fsidyn->ifsi=6;
      else if (strncmp(buffer,"iter_stagg_CHEB_rel_param",25)==0 && length==25) 
         fsidyn->ifsi=7;      
      else
         dserror("Coupling Algorithm COUPALGO unknown");
   }
   frchar("PREDICTOR"  ,buffer    ,&ierr);
   if (ierr==1)
   {            
      length = strlen(buffer);
      if (strncmp(buffer,"d(n)",4)==0 && length==4)
         fsidyn->ipre=1;
      else if (strncmp(buffer,"d(n)+dt*(1.5*v(n)-0.5*v(n-1))",29)==0 && length==29)
         fsidyn->ipre=2;
      else if (strncmp(buffer,"d(n)+dt*v(n)",12)==0 && length==12)
            fsidyn->ipre=3;	 
      else if (strncmp(buffer,"d(n)+dt*v(n)+0.5*dt^2*a(n)",26)==0 && length==26)
            fsidyn->ipre=4;	 
      else dserror("PREDICTOR unknown!\n"); 
   }   
   frchar("CONVCRIT"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      length = strlen(buffer);
      if (strncmp(buffer,"||g(i)||:sqrt(neq)",18)==0 && length==18) 
         fsidyn->inrmfsi=1;
      else if (strncmp(buffer,"||g(i)||:||g(0)||",17)==0 && length==17) 
         fsidyn->inrmfsi=2;
      else
         dserror("Convergence criterion CONVCRIT unknown!");     
   }
   frchar("ENERGYCHECK"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      length = strlen(buffer);
      if (strncmp(buffer,"No",2)==0 ||
          strncmp(buffer,"NO",2)==0 ||
	  strncmp(buffer,"no",2)==0 && length==2) 
         fsidyn->ichecke=0;
      else if (strncmp(buffer,"yes",3)==0 ||
               strncmp(buffer,"YES",3)==0 ||                 
	       strncmp(buffer,"Yes",3)==0 && length==3) 
         fsidyn->ichecke=1;
      else
         dserror("Parameter for ENERGYCHECK unknown!");     
   }

   frchar("NESTEDITER"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      length = strlen(buffer);
      if (strncmp(buffer,"No",2)==0 ||
          strncmp(buffer,"NO",2)==0 ||
	  strncmp(buffer,"no",2)==0 && length==2) 
         fsidyn->inest=0;
      else if (strncmp(buffer,"fluid+structure",15)==0 && length==15) 
         fsidyn->inest=1;
      else if (strncmp(buffer,"structure",9)==0 && length==9) 
         fsidyn->inest=2;
      else
         dserror("Parameter unknown: NESTEDITER");     
   }      
   frchar("ICHOPT"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      length = strlen(buffer);
      if (strncmp(buffer,"no",2)==0 && length==2) 
         fsidyn->ichopt=0;
      else if (strncmp(buffer,"yes",3)==0 && length==3) 
         fsidyn->ichopt=1;
      else
         dserror("Parameter unknown: ICHOPT");     
   }   
   frchar("IALE"    ,buffer    ,&ierr);
   if (ierr==1)
   {
      length = strlen(buffer);
      if (strncmp(buffer,"Pseudo_Structure",16)==0 && length==16)
         fsidyn->iale=1;
      else
         dserror("Parameter unknown: IALE");
   }
/*--------------read int */
   frint("ITECHAPP"  ,&(fsidyn->itechapp),&ierr);
   frint("ICHMAX"    ,&(fsidyn->ichmax)  ,&ierr);
   frint("ISDMAX"    ,&(fsidyn->isdmax)  ,&ierr);
   frint("NUMSTEP"   ,&(fsidyn->nstep)   ,&ierr);
   frint("ITEMAX"    ,&(fsidyn->itemax)  ,&ierr);
   frint("UPPSS"     ,&(fsidyn->uppss)   ,&ierr);
/*--------------read double */
   frdouble("TIMESTEP"   ,&(fsidyn->dt)     ,&ierr);
   frdouble("MAXTIME"    ,&(fsidyn->maxtime),&ierr);
   frdouble("TOLENCHECK" ,&(fsidyn->entol)  ,&ierr);
   frdouble("RELAX"      ,&(fsidyn->relax)  ,&ierr);
   frdouble("CONVTOL"    ,&(fsidyn->convtol),&ierr);      
   frread();
}
frrewind();
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpctr_dyn_fsi */

#endif
