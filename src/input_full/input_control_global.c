/*!----------------------------------------------------------------------
\file
\brief input of control information and general problem data

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
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
 |                                                          al 08/02    |
 | pointer to allocate eigensolution variables                          |
 | dedfined in global_control.c                                         |
 | struct _ALLEIG       *alleig;                                        |
 *----------------------------------------------------------------------*/
extern ALLEIG              *alleig;   
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


#ifdef WALLCONTACT
extern struct _WALL_CONTACT contact;
#endif

#ifdef D_LS
void inpctr_dyn_ls(LS_DYNAMIC *lsdyn);
#endif

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
      
      solv[genprob.numsf].fieldtyp = structure;
      inpctrsol(&(solv[genprob.numsf]));
   }
   if (genprob.probtyp == prb_opt)
   {
      if (genprob.numfld!=1) dserror("numfld != 1 for Structural Problem");
      
      solv = (SOLVAR*)CCACALLOC(genprob.numfld,sizeof(SOLVAR));
      
      solv[0].fieldtyp = structure;
      inpctrsol(&(solv[0]));
   }
   if (genprob.probtyp == prb_fluid)   
   {
      solv = (SOLVAR*)CCACALLOC(genprob.numfld,sizeof(SOLVAR));
      
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

      solv[genprob.numaf].fieldtyp = ale;
      inpctrsol(&(solv[genprob.numaf]));
   }
#ifdef D_LS
   if (genprob.probtyp == prb_twophase)
   {
     if (genprob.numfld!=2) dserror("numfld != 2 for Two Phase Fluid Flow Problem");
     
     /* initialize solver object */
     solv = (SOLVAR*)CCACALLOC(genprob.numfld,sizeof(SOLVAR));

     /* initialize solver for the fluid field */
     solv[genprob.numff].fieldtyp = fluid;
     inpctrsol(&(solv[genprob.numff]));
     
     /* initialize solver for the levelset field */
     solv[genprob.numls].fieldtyp = levelset;
     inpctrsol(&(solv[genprob.numls]));
   }
#endif
   
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
INT  ierr;
INT  restart;
char buffer[50];
#ifdef DEBUG 
dstrc_enter("inpctrprob");
#endif

/*------------------------------------------ set initial field numbers */
genprob.numsf=-1;
genprob.numff=-1;
genprob.numaf=-1;
genprob.numls=-1;
 
if (frfind("-PROBLEM SIZE")==0) dserror("frfind: PROBLEM SIZE not in input file");
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
/*-------------------------------------- default value for monitoring */
ioflags.monitor=0;
/*----------------------------------------- defualt values for output */
ioflags.struct_disp_file      =0;
ioflags.struct_stress_file    =0;
ioflags.struct_disp_gid       =0;
ioflags.struct_stress_gid     =0;
ioflags.struct_stress_gid_smo =0;
ioflags.fluid_sol_gid         =0;
ioflags.fluid_stress_gid      =0;
ioflags.fluid_sol_file        =0; 
ioflags.fluid_vis_file        =0;
ioflags.ale_disp_file         =0;
ioflags.ale_disp_gid          =0;
ioflags.relative_displ        =0;


if (frfind("-PROBLEM TYP")==0) dserror("frfind: PROBLEM TYP not in input file");
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
    if (strncmp("Two_Phase_Fluid_Flow"       ,buffer,20)==0) genprob.probtyp = prb_twophase;    
  }
  
#ifdef D_XFEM
  frchar("ENRONOFF"   ,buffer    ,&ierr);
  if (ierr==1)
  {
    if (strncmp(buffer,"yes",3)==0)
      genprob.xfem_on_off=1;
    else
      genprob.xfem_on_off=0;
  }
#endif  
  
  frchar("TIMETYP"   ,buffer,            &ierr);
  if (ierr==1)
  {
    if (strncmp("Static" ,buffer,6)==0) genprob.timetyp=time_static;
    if (strncmp("Dynamic",buffer,7)==0) genprob.timetyp=time_dynamic;
  }
  frint("RESTART"    ,&(restart),&ierr);
  if (ierr==1)
  {
     if (genprob.restart==0 && restart>0)
     dserror("Restart defined in input file but not as program argument!\n");
     genprob.restart=restart;
  }

  frint("NUMFIELD",&(genprob.numfld),&ierr);

  frint("GRADERW",&(genprob.graderw),&ierr);
  frread();
}

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
#ifdef D_LS
if (genprob.probtyp==prb_twophase)
{
  /* redefine the number of materials */
  genprob.nmat     = 3;
  /* set field numbers */
  genprob.numfld   = 2;
  genprob.numff    = 0;
  genprob.numls    = 1;  
}
#endif

if (frfind("---IO")==1)
{
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
      if (strncmp(buffer,"yes",3)==0) ioflags.struct_stress_gid=1;
      if (strncmp(buffer,"YES",3)==0) ioflags.struct_stress_gid=1;
      if (strncmp(buffer,"Yes",3)==0) ioflags.struct_stress_gid=1;
    }
    frchar("STRUCT_STRESS_GID_SMO",buffer,&ierr);
    if (ierr)
    {
      if (strncmp(buffer,"yes",3)==0) ioflags.struct_stress_gid_smo=1;
      if (strncmp(buffer,"YES",3)==0) ioflags.struct_stress_gid_smo=1;
      if (strncmp(buffer,"Yes",3)==0) ioflags.struct_stress_gid_smo=1;
    }
    frchar("FLUID_STRESS_GID",buffer,&ierr);
    if (ierr)
    {
      if (strncmp(buffer,"yes",3)==0) ioflags.fluid_stress_gid=1;
      if (strncmp(buffer,"YES",3)==0) ioflags.fluid_stress_gid=1;
      if (strncmp(buffer,"Yes",3)==0) ioflags.fluid_stress_gid=1;
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
    frint("RELATIVE_DISP_NUM",&(ioflags.relative_displ),&ierr);
    frread();
  }
}

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
INT  ierr;
char buffer[50];
INT  counter;
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
statvar->isrelstepsize=0;
/*------------------------------------------------------- start reading */
if (frfind("-STATIC")==1)
{
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
    frint("GRADERW",&(statvar->graderw),&ierr);

    frdouble("TOLRESID",&(statvar->tolresid),&ierr);
    frdouble("TOLDISP" ,&(statvar->toldisp) ,&ierr);
    frdouble("STEPSIZE",&(statvar->stepsize),&ierr);
    frdouble("ARCSCL"  ,&(statvar->arcscl)  ,&ierr);
    /*-------------------------------------------------------------------*/
    frchar("VARSTEPSI",buffer,&ierr);
    if (ierr)
    {
      if (strncmp(buffer,"yes",3)==0) statvar->isrelstepsize=1;
      if (strncmp(buffer,"YES",3)==0) statvar->isrelstepsize=1;
      if (strncmp(buffer,"Yes",3)==0) statvar->isrelstepsize=1;
    }
    /*-------------------------------------------------------------------*/
    frread();
  }
}
/*--------------------------------- in nonlinear case find control node */
if (statvar->nonlinear)
{
  /*---------------------------------------------------- start reading */
  if (frfind("-CONTROL NODE")==1)
  {
    frread();
    while(strncmp(allfiles.actplace,"------",6)!=0)
    {
      frint("NODE",&(statvar->control_node_global),&ierr);
      if (ierr) statvar->control_node_global--;
      frint("DOF",&(statvar->control_dof),&ierr);
      frread();
    }
  }
  frrewind();
}
/*----------------------------------------------------------------*/
if (ioflags.relative_displ>0)
{
  if (frfind("-RELATIVE DISPLACEMENT NODES")==1)
  {
    frread();
    counter=0;
    while(strncmp(allfiles.actplace,"------",6)!=0)
    {
      frint("NODE",&(statvar->reldisnode_ID[counter]),&ierr);
      if (ierr) statvar->reldisnode_ID[counter]--;
      frint("DOF",&(statvar->reldis_dof[counter]),&ierr);
      counter++;
      frread();
    }
  }
}
/*----------------------------------------------------------------*/
if (statvar->isrelstepsize==1)
{
  if (frfind("-VARIABLE STEP SIZES")==1)
  {
    frread();
    counter=0;
    while(strncmp(allfiles.actplace,"------",6)!=0)
    {
      frint("NSTEP",&(statvar->actstep[counter]),&ierr);
      frdouble("ACTSTEPSIZE",&(statvar->actstepsize[counter]),&ierr);
      counter++;
      frread();
    }
    statvar->numcurve = counter;
    if(counter>19)  dserror("not more than 20 different stepsizes!");
  }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpctrstat */


/*----------------------------------------------------------------------*
 | input of data for eigensolution                           al 8/02    |
 *----------------------------------------------------------------------*/
void inpctreig()
{
INT    i;
FIELD *actfield;
#ifdef DEBUG 
dstrc_enter("inpctreig");
#endif
/*----------------------------------------------------------------------*/
alleig = (ALLEIG*)CCACALLOC(1,sizeof(ALLEIG));
if (!alleig) dserror("Allocation of ALLEIG failed");
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
   case levelset:
   break;
   case structure:
     inpctr_eig_struct(alleig);
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
} /* end of inpctreig */

/*----------------------------------------------------------------------*
 | input of eigensolution problem data                       al 8/02    |
 *----------------------------------------------------------------------*/
void inpctr_eig_struct(ALLEIG *alleig)
{ 

INT    ierr;
char   buffer[50];
#ifdef DEBUG 
dstrc_enter("inpctr_eig_struct");
#endif

if (frfind("--EIGENVALUE ANALYSIS")==0) goto end;
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
/*--------------read chars */
   frchar("SOLTYP"   ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"SUBSPACE",8)==0)
          
          alleig->soltyp=subspace;
      else
          alleig->soltyp=eig_none;
   }
/*--------------read INT */
   frint("STURM ",&(alleig->sturm ),&ierr);
   frint("SUBTYP",&(alleig->subtyp),&ierr);
   frint("IFSH"  ,&(alleig->ifsh  ),&ierr);
   frint("ILMP"  ,&(alleig->ilmp  ),&ierr);
   frint("RANGE" ,&(alleig->range ),&ierr);
   frint("NUMVEC",&(alleig->numvec),&ierr);
   frint("NROOT" ,&(alleig->nroot ),&ierr);
   frint("ITEMAX",&(alleig->itemax),&ierr);
   frint("IFCTR" ,&(alleig->ifctr ),&ierr);
/*--------------read DOUBLE */
   frdouble("TOLEIG",&(alleig->toleig),&ierr);
   frdouble("SHIFT" ,&(alleig->shift ),&ierr);
   frdouble("BOULO" ,&(alleig->boulo ),&ierr);
   frdouble("BOUUP" ,&(alleig->bouup ),&ierr);
   frread();
}
frrewind();

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpctr_eig_struct */


/*----------------------------------------------------------------------*
 | input of dynamic problem data                          m.gee 2/01    |
 |                                                         genk 3/01    |
 *----------------------------------------------------------------------*/
void inpctrdyn()
{
INT    i;
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
         inpctr_dyn_fluid(alldyn[i].fdyn);
#else
         dserror("General FLUID problem not defined in Makefile!!!");
#endif            
         break;
       case ale:
#ifdef D_ALE
         alldyn[i].adyn = (ALE_DYNAMIC*)CCACALLOC(1,sizeof(ALE_DYNAMIC));
         inpctr_dyn_ale(alldyn[i].adyn);
#else
         dserror("General ALE problem not defined in Makefile!!!");
#endif
         break;
       case structure:
         alldyn[i].sdyn = (STRUCT_DYNAMIC*)CCACALLOC(1,sizeof(STRUCT_DYNAMIC));
         inpctr_dyn_struct(alldyn[i].sdyn);
         break;
#ifdef D_LS
       case levelset:
         alldyn[i].lsdyn = (LS_DYNAMIC*)CCACALLOC(1,sizeof(LS_DYNAMIC));
         inpctr_dyn_ls(alldyn[i].lsdyn);
         break;	
#endif	
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
INT    ierr;
char   buffer[50];
#ifdef DEBUG 
dstrc_enter("inpctr_dyn_struct");
#endif

sdyn->updevry_disp=1;
sdyn->updevry_stress=1;
sdyn->res_write_evry=1;
sdyn->eigen=0;

if (frfind("-STRUCTURAL DYNAMIC")==0) goto end;
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
/*--------------read chars */
   frchar("DYNAMICTYP",buffer,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"Centr_Diff",10)==0) sdyn->Typ = centr_diff;
      if (strncmp(buffer,"Gen_Alfa",8)==0)    sdyn->Typ = gen_alfa;
      if (strncmp(buffer,"Gen_EMM",7)==0)    sdyn->Typ = Gen_EMM;
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

/*-----------------read INT */
   frint("EIGEN"      ,&(sdyn->eigen)         ,&ierr);
   frint("NUMSTEP"    ,&(sdyn->nstep)         ,&ierr);
   frint("MAXITER"    ,&(sdyn->maxiter)       ,&ierr);
   frint("RESEVRYDISP",&(sdyn->updevry_disp)  ,&ierr);
   frint("RESEVRYSTRS",&(sdyn->updevry_stress),&ierr);
   frint("RESTARTEVRY",&(sdyn->res_write_evry),&ierr);
   frint("CONTACT"    ,&(sdyn->contact)       ,&ierr);
#ifdef WALLCONTACT
   frint("CET_flag",   &(contact.CET_flag)    ,&ierr);
   frint("FR_flag",    &(contact.FR_flag)     ,&ierr);
#endif
/*--------------read DOUBLE */
   frdouble("TIMESTEP",&(sdyn->dt)     ,&ierr);
#ifdef WALLCONTACT
   if (ierr) contact.dt = sdyn->dt;
#endif
   frdouble("MAXTIME" ,&(sdyn->maxtime),&ierr);
   frdouble("BETA"    ,&(sdyn->beta)   ,&ierr);
   frdouble("GAMMA"   ,&(sdyn->gamma)  ,&ierr);
   frdouble("ALPHA_M" ,&(sdyn->alpha_m),&ierr);
   frdouble("ALPHA_F" ,&(sdyn->alpha_f),&ierr);
#ifdef GEMM   
   frdouble("XSI"     ,&(sdyn->xsi)    ,&ierr);
#endif
#ifdef WALLCONTACT
   frdouble("NPP"     ,&(contact.n_pen_par) , &ierr);
   frdouble("TPP"     ,&(contact.t_pen_par) , &ierr);
   frdouble("FR_COEF" ,&(contact.fr_coef)   , &ierr);
#endif
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

end:
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
INT    ierr;
INT    i;
INT    thetafound=0;
char   buffer[50];
#ifdef DEBUG 
dstrc_enter("inpctr_dyn_fluid");
#endif

/*-------------------------------------------------- set default values */
fdyn->dyntyp=0;
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
fdyn->liftdrag=0;
fdyn->numdf=genprob.ndim+1;  
fdyn->numcont=0;
fdyn->uppss=1;  
fdyn->upout=1;  
fdyn->upres=1;  
fdyn->res_write_evry=20;
fdyn->nstep=1;
fdyn->stchk=5;  
fdyn->nums=0;
fdyn->iprerhs=1;
fdyn->itemax=3;
fdyn->dt=0.01;    
fdyn->maxtime=1000.0;
fdyn->theta=0.5;
fdyn->ittol=EPS6; 
fdyn->sttol=EPS6; 
fdyn->turbu=0;
fdyn->adaptive=0;	/* default: no adaptive time stepping */
fdyn->time_rhs=1;	/* default: build time rhs as W.A. Wall describes */
fdyn->lte = 1.0e-03;	/* default: local truncation error for adapt. time*/
fdyn->max_dt = 1.0;	/* default: maximal time step size to 1.0 */
fdyn->min_dt = 0.0;	/* default: minimal time step size to 0.0 */
fdyn->mlfem = 0;	/* default: no multi-level algorithm */

 
if (frfind("-FLUID DYNAMIC")==0) goto end;
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
/*--------------read chars */
   frchar("DYNAMICTYP",buffer,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"Projection_Method",17)==0) 
         fdyn->dyntyp=1;
      else if (strncmp(buffer,"Nlin_Time_Int",13)==0)
         fdyn->dyntyp=0;
      else
         dserror("DYNAMICTYP unknown");      
   }
   frchar("TIMEINTEGR"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"Stationary",10)==0) 
         fdyn->iop=0; 
      else if (strncmp(buffer,"Gen_Alfa",8)==0) 
         fdyn->iop=1; 
      else if (strncmp(buffer,"Gen_Alpha",9)==0) 
         fdyn->iop=1;     
      else if (strncmp(buffer,"One_Step_Theta",14)==0) 
         fdyn->iop=4;    
      else if (strncmp(buffer,"BDF2",4)==0) 
         fdyn->iop=7;  
      else
         dserror("TIMEINTEGR-Type unknown");
   }
   frchar("STARTINGALGO"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"One_Step_Theta",14)==0) 
         fdyn->iops=4;
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
      else if (strncmp(buffer,"BELTRAMI-FLOW",13)==0) 
         fdyn->init=8;
      else if (strncmp(buffer,"KIM-MOIN-FLOW",13)==0) 
         fdyn->init=9;
      else if (strncmp(buffer,"BREAKING-DAM",12)==0) 
         fdyn->init=10;
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
   frchar("LIFTDRAG"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"yes",3)==0 ||
          strncmp(buffer,"YES",3)==0 ||                 
	  strncmp(buffer,"Yes",3)==0    )
         fdyn->liftdrag=1;
      else if (strncmp(buffer,"No",2)==0 ||
               strncmp(buffer,"NO",2)==0 ||
	       strncmp(buffer,"no",2)==0   )
         fdyn->liftdrag=0;
      else
         dserror("LIFTDRAG unknown!\n");	 
   }
   frchar("TURBULENCE"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"No",2)==0 ||
          strncmp(buffer,"NO",2)==0 ||
	  strncmp(buffer,"no",2)==0   ) 
         fdyn->turbu=0;
      else if (strncmp(buffer,"algebraic",9)==0) 
         fdyn->turbu=1;
      else if (strncmp(buffer,"kappa-eps",9)==0) 
         fdyn->turbu=2;                      
      else if (strncmp(buffer,"kappa-omega",11)==0) 
         fdyn->turbu=3;                      
      else
         dserror("TURBULENCE unknown!");
   }   
   frchar("DISC_CAPT"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"No",2)==0 ||
          strncmp(buffer,"NO",2)==0 ||
	  strncmp(buffer,"no",2)==0   ) 
         fdyn->dis_capt=0;
      else if (strncmp(buffer,"Yes",3)==0) 
         fdyn->dis_capt=1;
      else
         dserror("DISC_CAPT unknown!");
   }
   frchar("ADAPT_TIME"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"No",2)==0 ||
          strncmp(buffer,"NO",2)==0 ||
	  strncmp(buffer,"no",2)==0   ) 
         fdyn->adaptive=0;
      else if (strncmp(buffer,"Yes",3)==0 ||
               strncmp(buffer,"YES",3)==0 ||
	       strncmp(buffer,"yes",3)==0   ) 
         fdyn->adaptive=1;
      else
         dserror("ADAPT_TIME can not be read (yes/no)!");     
   }
   frchar("TIME_RHS"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"mass",4)==0) 
         fdyn->time_rhs=0;
      else if (strncmp(buffer,"classic",7)==0) 
         fdyn->time_rhs=1;
      else
         dserror("TIME_RHS unknown!");     
   }
   frchar("CONVECTERM"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"convective",10)==0) 
         fdyn->conte=0;
      else if (strncmp(buffer,"divergence",10)==0) 
         fdyn->conte=1;                      
      else if (strncmp(buffer,"skew_symmetric",14)==0) 
         fdyn->conte=2;                      
      else
         dserror("Unknown convective term!");
   }   
   frchar("VISCTERM"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"conventional",12)==0) 
         fdyn->vite=0;
      else if (strncmp(buffer,"stress_divergence",17)==0) 
         fdyn->vite=1;                      
      else
         dserror("Unknown viscous term!");
   }   
   frchar("SUBGRIDVISC"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"no",2)==0) 
         fdyn->sgvisc=0;
      else if (strncmp(buffer,"artificial",10)==0) 
         fdyn->sgvisc=1;                      
      else if (strncmp(buffer,"Smagorinsky",11)==0) 
         fdyn->sgvisc=2;                      
      else
         dserror("Unknown subgrid viscosity!");
   }   

/*--------------read INT */
   frint("NUMCONT"    ,&(fdyn->numcont)       ,&ierr);
   frint("UPPSS"      ,&(fdyn->uppss)         ,&ierr);
   frint("UPOUT"      ,&(fdyn->upout)         ,&ierr);
   frint("UPRES"      ,&(fdyn->upres)         ,&ierr);
   frint("RESSTEP"    ,&(fdyn->resstep)       ,&ierr);
   frint("RESTARTEVRY",&(fdyn->res_write_evry),&ierr);
   frint("NUMSTEP"    ,&(fdyn->nstep)         ,&ierr);
   frint("STEADYSTEP" ,&i                     ,&ierr);   
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
   
/*--------------read DOUBLE */
   frdouble("TIMESTEP" ,&(fdyn->dt)     ,&ierr);
   frdouble("MAXTIME"  ,&(fdyn->maxtime),&ierr);
   frdouble("ALPHA_M"  ,&(fdyn->alpha_m),&ierr);
   frdouble("ALPHA_F"  ,&(fdyn->alpha_f),&ierr);
   if (thetafound==0)
   {
      frdouble("THETA"    ,&(fdyn->theta)  ,&ierr);
      if (ierr==1) thetafound++;
   }
   frdouble("CONVTOL"     ,&(fdyn->ittol)  ,&ierr);
   frdouble("STEADYTOL"   ,&(fdyn->sttol) ,&ierr);
   frdouble("START_THETA" ,&(fdyn->thetas),&ierr);
   frdouble("INT_LENGHT"  ,&(fdyn->lenght),&ierr);
   frdouble("ROUGHTNESS"  ,&(fdyn->rought),&ierr);
   frdouble("SC_COORD_X"  ,&(fdyn->coord_scale[0]),&ierr);
   frdouble("SC_COORD_Y"  ,&(fdyn->coord_scale[1]),&ierr);
   frdouble("MAX_DT"      ,&(fdyn->max_dt),&ierr);
   frdouble("MIN_DT"      ,&(fdyn->min_dt),&ierr);
   frdouble("LOC_TRUN_ERR",&(fdyn->lte)   ,&ierr);
   frdouble("SMAGCONST",&(fdyn->smagcon),&ierr);

   frread();
}
frrewind();

/*----------------------------------------------------------------------*/
if (frfind("-MULTILEVEL FLUID DYNAMIC")==0) goto end;
frread();
/* allocate fluid mulitlevel variables */
fdyn->mlvar = (FLUID_DYN_ML*)CCACALLOC(1,sizeof(FLUID_DYN_ML));

while(strncmp(allfiles.actplace,"------",6)!=0)
{
/*--------------read chars */
   frchar("MULEVELFEM"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"no",2)==0) 
         fdyn->mlfem=0;
      else if (strncmp(buffer,"yes",3)==0) 
         fdyn->mlfem=1;
      else
         dserror("unknown Multilevel");
   }
   frchar("TRANSTERM"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"no_approx",9)==0) 
         fdyn->mlvar->transterm=0;
      else if (strncmp(buffer,"neglect_sst",11)==0) 
         fdyn->mlvar->transterm=1;
      else if (strncmp(buffer,"neglect_lst",11)==0) 
         fdyn->mlvar->transterm=2;
      else if (strncmp(buffer,"neglect_both",12)==0) 
         fdyn->mlvar->transterm=3;
      else
         dserror("unknown treatment of transient term");
   }
   frchar("QUASTATBUB"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"no",2)==0) 
         fdyn->mlvar->quastabub=0;
      else if (strncmp(buffer,"yes",3)==0) 
         fdyn->mlvar->quastabub=1;
      else
         dserror("unknown quasi-static bubbles");
   }
   frchar("CONVECVEL"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"no_approx",9)==0) 
         fdyn->mlvar->convel=0;
      else if (strncmp(buffer,"ls_approx",9)==0) 
         fdyn->mlvar->convel=1;
      else
         dserror("unknown convective velocity");
   }	 
   frchar("SMSGVISC"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"no",2)==0) 
         fdyn->mlvar->smsgvi=0;
      else if (strncmp(buffer,"artificial",10)==0) 
         fdyn->mlvar->smsgvi=1;                      
      else if (strncmp(buffer,"Smagorinsky",11)==0) 
         fdyn->mlvar->smsgvi=2;                      
      else if (strncmp(buffer,"dynamic",7)==0) 
         fdyn->mlvar->smsgvi=3;                      
      else
         dserror("Unknown submesh subgrid viscosity!");
   }   
   frchar("SMSTABI"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"no",2)==0) 
         fdyn->mlvar->smstabi=0;
      else if (strncmp(buffer,"yes",3)==0) 
         fdyn->mlvar->smstabi=1;
      else
         dserror("unknown submesh stabilization");
   }
   frchar("SMSTABDO"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"USFEM_instat",12)==0) 
         fdyn->mlvar->smstado=-1;
      else if (strncmp(buffer,"USFEM_CDR",9)==0) 
         fdyn->mlvar->smstado=-2;
      else if (strncmp(buffer,"USFEM_stat",10)==0) 
         fdyn->mlvar->smstado=-3;
      else if (strncmp(buffer,"GLS-_instat",11)==0) 
         fdyn->mlvar->smstado=1;
      else if (strncmp(buffer,"GLS-_CDR",8)==0) 
         fdyn->mlvar->smstado=2;
      else if (strncmp(buffer,"GLS-_stat",9)==0) 
         fdyn->mlvar->smstado=3;
      else
         dserror("unknown differential operator for submesh stabilization");
   }
   frchar("SMSTABNORM"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"L_1",3)==0) 
         fdyn->mlvar->smstano=1;
      else if (strncmp(buffer,"L_2",3)==0) 
         fdyn->mlvar->smstano=2;
      else
         dserror("unknown submesh stabilization norm");
   }
   frchar("SMSTABMK"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"-1",2)==0) 
         fdyn->mlvar->smstamk=-1;
      else if (strncmp(buffer,"0",1)==0) 
         fdyn->mlvar->smstamk=0;
      else if (strncmp(buffer,"1",1)==0) 
         fdyn->mlvar->smstamk=1;
      else
         dserror("unknown submesh mk");
   }
   frchar("SMSTABNINTHS"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"at_center",9)==0) 
         fdyn->mlvar->smstani=0;
      else if (strncmp(buffer,"every_intpt",11)==0) 
         fdyn->mlvar->smstani=1;
      else
         dserror("unknown submesh number of int. points for elesize");
   }
   frchar("SMORDER"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"linear",6)==0) 
         fdyn->mlvar->smorder=1;
      else if (strncmp(buffer,"quadratic",9)==0) 
         fdyn->mlvar->smorder=2;
      else
         dserror("unknown submesh element order");
   }
   frchar("SMNONUNIF"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"no",2)==0) 
         fdyn->mlvar->smnunif=0;
      else if (strncmp(buffer,"shishkin_type",13)==0) 
         fdyn->mlvar->smnunif=1;
      else if (strncmp(buffer,"quadratic",9)==0) 
         fdyn->mlvar->smnunif=2;
      else if (strncmp(buffer,"cubic",5)==0) 
         fdyn->mlvar->smnunif=3;
      else
         dserror("unknown sub-submesh element order");
   }
   frchar("SSMORD"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"linear",6)==0) 
         fdyn->mlvar->ssmorder=1;
      else if (strncmp(buffer,"quadratic",9)==0) 
         fdyn->mlvar->ssmorder=2;
      else
         dserror("unknown sub-submesh element order");
   }
   frchar("SSMNNUN"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"no",2)==0) 
         fdyn->mlvar->ssmnunif=0;
      else if (strncmp(buffer,"shishkin_type",13)==0) 
         fdyn->mlvar->ssmnunif=1;
      else if (strncmp(buffer,"quadratic",9)==0) 
         fdyn->mlvar->ssmnunif=2;
      else if (strncmp(buffer,"cubic",5)==0) 
         fdyn->mlvar->ssmnunif=3;
      else
         dserror("unknown sub-submesh element order");
   }
/*--------------read int */
   frint("SMELESIZE"      ,&(fdyn->mlvar->smesize)  ,&ierr);
   frint("SMSTABPAR"      ,&(fdyn->mlvar->smstapa)  ,&ierr);
   frint("SMELEMENTS"     ,&(fdyn->mlvar->smelenum) ,&ierr);
   frint("SMNUMGP"        ,&(fdyn->mlvar->smnumgp)  ,&ierr);
   frint("SSMELEM"        ,&(fdyn->mlvar->ssmelenum),&ierr);
   frint("SSMNGP"         ,&(fdyn->mlvar->ssmnumgp) ,&ierr);
/*--------------read double */
   frdouble("SMSMAGCON"   ,&(fdyn->mlvar->smsmagcon),&ierr);

   frread();
}
frrewind();
/*----------------------------------------------------------------------*/

end:
/*-------------------------------------------------- plausibility check */
if (fdyn->freesurf>0)
#ifndef D_FSI
dserror("Freesurf requires FSI functions to compile in!\n");
#endif

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
INT    ierr;
INT    i;
INT    length;
INT    mod_res_write;
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
fsidyn->upres=1;
fsidyn->res_write_evry=1;
fsidyn->dt=ONE/TEN;	
fsidyn->maxtime=1000.0;
fsidyn->entol=EPS6;  
fsidyn->relax=ONE;  
fsidyn->convtol=EPS6;
    
if (frfind("-FSI DYNAMIC")==0) goto end;
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
/*--------------read INT */
   frint("ITECHAPP"   ,&(fsidyn->itechapp)         ,&ierr);
   frint("ICHMAX"     ,&(fsidyn->ichmax)           ,&ierr);
   frint("ISDMAX"     ,&(fsidyn->isdmax)           ,&ierr);
   frint("NUMSTEP"    ,&(fsidyn->nstep)            ,&ierr);
   frint("ITEMAX"     ,&(fsidyn->itemax)           ,&ierr);
   frint("UPPSS"      ,&(fsidyn->uppss)            ,&ierr);
   frint("UPRES"      ,&(fsidyn->upres)            ,&ierr);
   frint("RESTARTEVRY",&(fsidyn->res_write_evry)   ,&ierr);
/*--------------read DOUBLE */
   frdouble("TIMESTEP"    ,&(fsidyn->dt)     ,&ierr);
   frdouble("MAXTIME"     ,&(fsidyn->maxtime),&ierr);
   frdouble("TOLENCHECK"  ,&(fsidyn->entol)  ,&ierr);
   frdouble("RELAX"       ,&(fsidyn->relax)  ,&ierr);
   frdouble("CONVTOL"     ,&(fsidyn->convtol),&ierr);  
   frread();
}
frrewind();

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpctr_dyn_fsi */

#endif


/*!---------------------------------------------------------------------                                         
\brief input of the ALE DYNAMIC block in the input-file

<pre>                                                           ck 12/02

In this routine the data in the ALE DYNAMIC block of the input file
are read and stored in adyn	       

</pre>
\param  *adyn 	  ALE_DATA       (o)	   
\return void                                                                       

------------------------------------------------------------------------*/
#ifdef D_ALE
void inpctr_dyn_ale(ALE_DYNAMIC *adyn)
{
INT    ierr;
char buffer[50];

#ifdef DEBUG 
dstrc_enter("inpctr_dyn_ale");
#endif
/*---------------------------------------------------- some defaults ---*/
adyn->num_initstep = 0;
adyn->step = 0;
adyn->nstep = 10;
adyn->dt = 0.1;
adyn->maxtime = 1000.0;
adyn->updevry_disp = 20;
adyn->measure_quality = no_quality;
adyn->typ = classic_lin;

if (frfind("-ALE DYNAMIC")==0) goto end;
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
/*--------------read chars */
   frchar("ALE_TYPE",buffer,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"classic_lin",11)==0) adyn->typ = classic_lin;
      else if (strncmp(buffer,"min_Je_stiff",12)==0) adyn->typ = min_Je_stiff;
      else if (strncmp(buffer,"two_step",8)==0) adyn->typ = two_step;
      else if (strncmp(buffer,"springs",7)==0) adyn->typ = springs;
      else if (strncmp(buffer,"laplace",7)==0) adyn->typ = laplace;
      else dserror("unknown ALE_TYPE");
   }
   frchar("QUALITY",buffer,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"aspect_ratio",12)==0) adyn->measure_quality = aspect_ratio;
      else if (strncmp(buffer,"ASPECT_RATIO",12)==0) adyn->measure_quality = aspect_ratio;
      else if (strncmp(buffer,"corner_angle",12)==0) adyn->measure_quality = corner_angle;
      else if (strncmp(buffer,"CORNER_ANGLE",12)==0) adyn->measure_quality = corner_angle;
      else if (strncmp(buffer,"min_J",5)==0) adyn->measure_quality = min_detF;
      else if (strncmp(buffer,"MIN_J",5)==0) adyn->measure_quality = min_detF;
      else if (strncmp(buffer,"none",4)==0) adyn->measure_quality = no_quality;
      else if (strncmp(buffer,"NONE",4)==0) adyn->measure_quality = no_quality;
      else dserror("unknown QUALITY in ALE DYNAMIC");
   }
/*--------------read INT */
   frint("NUMSTEP" ,   &(adyn->nstep) ,       &ierr);
   frint("RESEVRYDISP",&(adyn->updevry_disp), &ierr);
   frint("NUM_INITSTEP",&(adyn->num_initstep), &ierr);
/*--------------read DOUBLE */
   frdouble("TIMESTEP", &(adyn->dt)     , &ierr);
   frdouble("MAXTIME" , &(adyn->maxtime), &ierr);
   
   frread();
}
frrewind();

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpctr_dyn_ale */
#endif


/*!---------------------------------------------------------------------                                         
\brief input of the LS DYNAMIC block in the input-file

<pre>                                                        irhan 04/04

In this routine the data in the LS DYNAMIC block of the input file
are read and stored in adyn	       

</pre>
\param  *lsdyn->lsdata 	  LS_GEN_DATA       (o)	   
\return void                                                                       

------------------------------------------------------------------------*/

#ifdef D_LS
void inpctr_dyn_ls(LS_DYNAMIC *lsdyn)
{
  INT ierr;
  char buffer[50];
  
#ifdef DEBUG 
  dstrc_enter("inpctr_dyn_ls");
#endif
/*----------------------------------------------------------------------*/  
  
  /* allocate memory for lsdata */
  lsdyn->lsdata = (LS_GEN_DATA*)CCACALLOC(1,sizeof(LS_GEN_DATA));  
  /* set some parameters */
  lsdyn->lsdata->boundary_on_off = 0;
  lsdyn->lsdata->reinitflag = 0;
  lsdyn->lsdata->print_on_off = 0;  
  /* read in lsdyn */
  if (frfind("-LEVELSET DYNAMIC")==0) dserror("frfind: LEVELSET DYNAMIC not in input file");
  frread();
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /* read INT */
    frint("NSTEP", &(lsdyn->nstep), &ierr);
    frint("NFSTEP", &(lsdyn->nfstep), &ierr);
    frint("NFREINIT", &(lsdyn->nfreinit), &ierr);    
    frint("NFREPRN", &(lsdyn->nfreprn), &ierr);    
    frint("ITEMAX", &(lsdyn->itemax), &ierr);	
    frint("INIT", &(lsdyn->init), &ierr);
    frint("IOP", &(lsdyn->iop), &ierr);
    frint("ITE", &(lsdyn->ite), &ierr);
    frint("ITCHK", &(lsdyn->itchk), &ierr);
    frint("STCHK", &(lsdyn->stchk), &ierr);
    frint("ITNRM", &(lsdyn->itnrm), &ierr);
    /* read DOUBLE */
    frdouble("DT", &(lsdyn->dt)  , &ierr);
    frdouble("MAXTIME", &(lsdyn->maxtime), &ierr);
    frdouble("ITTOL", &(lsdyn->ittol) , &ierr);
    frdouble("STTOL", &(lsdyn->sttol) , &ierr);
    frdouble("THETA", &(lsdyn->theta) , &ierr);
    frread();
  }
  frrewind();

  /* read in lsdata */
  if (frfind("-LEVELSET GENERAL DATA")==0) dserror("frfind: LEVELSET GENERAL DATA not in input file");
  frread();
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /* read INT */
    frint("FLAGVEL", &(lsdyn->lsdata->flag_vel), &ierr);
    frint("NUMLAYER", &(lsdyn->lsdata->numlayer), &ierr);
    frint("NUMREINIT", &(lsdyn->lsdata->numreinit), &ierr);
    frint("NUMCIRC", &(lsdyn->lsdata->numcirc), &ierr);
    frint("NUMLINE", &(lsdyn->lsdata->numline), &ierr);    
    /* read DOUBLE */
    frdouble("XC1", &(lsdyn->lsdata->xc1), &ierr);
    frdouble("YC1", &(lsdyn->lsdata->yc1), &ierr);
    frdouble("RAD1", &(lsdyn->lsdata->rad1), &ierr);
    frdouble("XC2", &(lsdyn->lsdata->xc2), &ierr);
    frdouble("YC2", &(lsdyn->lsdata->yc2), &ierr);
    frdouble("RAD2", &(lsdyn->lsdata->rad2), &ierr);
    frdouble("XC3", &(lsdyn->lsdata->xc3), &ierr);
    frdouble("YC3", &(lsdyn->lsdata->yc3), &ierr);
    frdouble("RAD3", &(lsdyn->lsdata->rad3), &ierr);
    frdouble("XS1", &(lsdyn->lsdata->xs1), &ierr);
    frdouble("YS1", &(lsdyn->lsdata->ys1), &ierr);
    frdouble("XE1", &(lsdyn->lsdata->xe1), &ierr);
    frdouble("YE1", &(lsdyn->lsdata->ye1), &ierr);
    frdouble("XS2", &(lsdyn->lsdata->xs2), &ierr);
    frdouble("YS2", &(lsdyn->lsdata->ys2), &ierr);
    frdouble("XE2", &(lsdyn->lsdata->xe2), &ierr);
    frdouble("YE2", &(lsdyn->lsdata->ye2), &ierr);
    frdouble("XS3", &(lsdyn->lsdata->xs3), &ierr);
    frdouble("YS3", &(lsdyn->lsdata->ys3), &ierr);
    frdouble("XE3", &(lsdyn->lsdata->xe3), &ierr);
    frdouble("YE3", &(lsdyn->lsdata->ye3), &ierr);
    frdouble("RDT", &(lsdyn->lsdata->rdt), &ierr);    
    frdouble("EPSILON", &(lsdyn->lsdata->epsilon), &ierr);
    /* read character */
    frchar("BOUNDARYONOFF"   ,buffer    ,&ierr);
    if (ierr==1)
    {
      if (strncmp(buffer,"on",2)==0)
        lsdyn->lsdata->boundary_on_off = 1;
      else if (strncmp(buffer,"off",3)==0)
        lsdyn->lsdata->boundary_on_off = 0;
      else
        dserror("BOUNDARYONOFF unknown!\n");	 
    }
    frchar("SETVEL"   ,buffer    ,&ierr);
    if (ierr==1)
    {
      if (strncmp(buffer,"by_user",7)==0)
        lsdyn->lsdata->setvel = 1;
      else if (strncmp(buffer,"by_fluid",8)==0)
        lsdyn->lsdata->setvel = 2;
      else
        dserror("SETVEL unknown!\n");	 
    }    
    frchar("ISSTAB"   ,buffer    ,&ierr);
    if (ierr==1)
    {
      if (strncmp(buffer,"yes",3)==0)
        lsdyn->lsdata->isstab = 1;
      else if (strncmp(buffer,"no",2)==0)
        lsdyn->lsdata->isstab = 0;
      else
        dserror("ISSTAB unknown!\n");	 
    }
    frchar("LOCALIZATION"   ,buffer    ,&ierr);
    if (ierr==1)
    {
      if (strncmp(buffer,"on",2)==0)
        lsdyn->lsdata->localization = 1;
      else if (strncmp(buffer,"off",3)==0)
        lsdyn->lsdata->localization = 0;
      else
        dserror("LOCALIZATION unknown!\n");	 
    }
    frchar("REINITIALIZATION"   ,buffer    ,&ierr);
    if (ierr==1)
    {
      if (strncmp(buffer,"on",2)==0)
        lsdyn->lsdata->reinitialization = 1;
      else if (strncmp(buffer,"off",3)==0)
        lsdyn->lsdata->reinitialization = 0;
      else
        dserror("REINITIALIZATION unknown!\n");	 
    }
    frchar("ALGO"   ,buffer    ,&ierr);
    if (ierr==1)
    {
      if (strncmp(buffer,"explicit",8)==0)
        lsdyn->lsdata->algo = 0;
      else if (strncmp(buffer,"implicit",8)==0)
        lsdyn->lsdata->algo = 1;
      else
        dserror("ALGO unknown!\n");	 
    }
    frchar("ISSHARP"   ,buffer    ,&ierr);
    if (ierr==1)
    {
      if (strncmp(buffer,"yes",3)==0)
        lsdyn->lsdata->is_sharp = 1;
      else if (strncmp(buffer,"no",2)==0)
        lsdyn->lsdata->is_sharp = 0;
      else
        dserror("ISSHARP unknown!\n");	 
    }
    frread();
  }
  frrewind();
 end:
  
/*----------------------------------------------------------------------*/  
#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
} /* end of inpctr_dyn_ls */
#endif
