/*!----------------------------------------------------------------------
\file
\brief contains the routine 'inpctropt' which reads
       optimization control data

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#ifdef D_OPTIM                   /* include optimization code to ccarat */
/*----------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
#include "../headers/optimization.h"
/*----------------------------------------------------------------------*
 |                                                          al 06/01    |
 | optimization data                                                    |
 *----------------------------------------------------------------------*/
extern struct _OPTI *opt;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 | input of optimization                                a.lipka 5/01    |
 *----------------------------------------------------------------------*/
void inpctropt()
{
int    ierr;
int    i, idum, iloop;
int    ccdv=0; /* n.of.variables */
char   buffer[50]; 
OSNLP *nlp;
OSFSD *fsd;
#ifdef DEBUG 
dstrc_enter("inpctropt");
#endif
/*----------------------------------------------------------------------*/
opt = (OPTI*)calloc(1,sizeof(OPTI));
if (opt==NULL) dserror("Allocation of OPTIMIZATION failed");
/*------------------------------------------------- initialize data ----*/
  opt->opttype   = ot_topology_optimization;
  opt->numiter   = 0;
  opt->strategy  = os_none;
  opt->objective = oj_none;
/*----------------------------------------------------------------------*/
frfind("--OPTIMIZATION");
frread();
i=0;
while(strncmp(allfiles.actplace,"------",6)!=0)
{ /* loop optimization data */
/*------------------------------------ type of optimization strategy ---*/
  frchar(  "OPT_TYPE"   ,buffer          ,&ierr);
  if (ierr==1)
  {
     if (strncmp(buffer,"Shape",5)==0) opt->opttype=ot_shape_optimization;
  }
/*---------------------------------- total number of iteration steps ---*/
  frint(   "OPT_NUMITER"    ,&(opt->numiter) ,&ierr);
/*------------------------------------ type of optimization strategy ---*/
  frchar(  "OPT_STRATEGY"   ,buffer          ,&ierr);
  if (ierr==1)
  {
     if (strncmp(buffer,"FSD",3)==0 )
     {    
         opt->strategy=os_fsd;
         opt->strat.fsd = (OSFSD*)calloc(1,sizeof(OSFSD));
         if (opt->strat.fsd==NULL) 
                             dserror("Allocation of OPTIMIZATION failed");
         /*---------------------- initialize  data for fsd-stragedy ----*/
         fsd = opt->strat.fsd;
         fsd->numiter = 0;   
         fsd->fgrad   = implicit;
         fsd->grdinc  = 1.0E-07;   /* forward difference step increment */
         fsd->acc     = 1.0;  
         fsd->alpha   = 0.0;  
         fsd->beta    = 0.0;  
         fsd->delta   = 0.0;  
         fsd->gamma   = 0.0;  
     }
     else if (strncmp(buffer,"NLP",3)==0 )
     {    
         opt->strategy=os_nlp;
         opt->strat.nlp = (OSNLP*)calloc(1,sizeof(OSNLP));
         if (opt->strat.nlp==NULL) 
                             dserror("Allocation of OPTIMIZATION failed");
         /*---------------------- initialize  data for nlp-stragedy ----*/
         nlp = opt->strat.nlp;
         nlp->fgrad   = implicit;
         nlp->grdinc  = 1.0E-07;   /* forward difference step increment */
         nlp->sclvar  = no;
         nlp->sclfun  = 1.0;
         nlp->sclcon  = 1.0;
         nlp->scbou   = 1.0;
         nlp->numiter = 0;
         nlp->nlpldl  = 1;
         nlp->nlluen  = 1;
         nlp->maxfun  = 100;
         nlp->lise    = 1;
         nlp->amue    = 0.1;
         nlp->merit   = 1;
     }
  }
/*---------------------------------------------------- nlp-stragedy ----*/
  if (opt->strategy==os_nlp)
  {
    nlp = opt->strat.nlp;
/*-------------------------------------- read data for nlp-stragedy ----*/
/*--------------------------------- evaluation of function gradients ---*/
    frchar(  "NLP_GRAD"   ,buffer          ,&ierr);
    if (ierr==1)
    {
      if (strncmp(buffer,"Explicit",8)==0 )         nlp->fgrad=explicit;
      else if (strncmp(buffer,"Variational",11)==0 )nlp->fgrad=variational;
    }
/*------------------------------------ scaling factor for objective  ---*/
    frdouble("NLP_RELINC"     ,&(nlp->grdinc)  ,&ierr);
/*------------------------------------ type of scaling the variables ---*/
    frchar(  "NLP_SCLVAR"   ,buffer          ,&ierr);
    if (ierr==1)
    {
      if (strncmp(buffer,"Initial",8)==0 )         nlp->sclvar=initial;
      else if (strncmp(buffer,"Hessian",11)==0 )   nlp->sclvar=hessian;
    }
/*------------------------------------ scaling factor for objective  ---*/
    frdouble("NLP_SCLFUN"     ,&(nlp->sclfun)  ,&ierr);
/*--------------------------------- scaling factors for constraints  ---*/
    frdouble("NLP_SCLCON"     ,&(nlp->sclcon)  ,&ierr);
/*------------------------ variable for autom. scaling of objective  ---*/
    frdouble("NLP_SCBOU"      ,&(nlp->scbou)   ,&ierr);
/*------------------------ variable for autom. scaling of objective  ---*/
    frdouble("NLP_ACC"        ,&(nlp->acc)     ,&ierr);
/*------------------------------ number of iteration steps with nlp  ---*/
    frint(   "NLP_NUMITER"    ,&(nlp->numiter) ,&ierr);
/*------------------------ flag for control of nonlinear programming ---*/
    frint(   "NLP_LDL"        ,&(nlp->nlpldl)  ,&ierr);
/*--------------------------- lueunberger selfscaling bfgs parameter ---*/
    frint(   "NLP_LUEN"       ,&(nlp->nlluen)  ,&ierr);
/*------------------------ max. n. o. function calls in line search  ---*/
    frint(   "NLP_MAXFUN"     ,&(nlp->maxfun)  ,&ierr);
/*----------------------- flag for control of nonlinear programming  ---*/
    frint(   "NLP_LISE"       ,&(nlp->lise  )  ,&ierr);
/*------------------------------- factor for line search adjustment  ---*/
    frdouble("NLP_AMUE"       ,&(nlp->amue)    ,&ierr);
/*--------------------------- aug. lagrangian / l1-penalty function  ---*/
    frint(   "NLP_MERIT"      ,&(nlp->merit  ) ,&ierr);
/*----------------------------------------------------------------------*/
  }
/*---------------------------------------------------- nlp-strategy ----*/

/*---------------------------------------------------- fsd-stragedy ----*/
  if (opt->strategy==os_fsd)
  {
    fsd = opt->strat.fsd;
/*-------------------------------------- read data for fsd-stragedy ----*/
/*--------------------------------- evaluation of function gradients ---*/
    frchar(  "FSD_GRAD"   ,buffer          ,&ierr);
    if (ierr==1)
    {
      if (strncmp(buffer,"Explicit",8)==0 )         fsd->fgrad=explicit;
      else if (strncmp(buffer,"Variational",11)==0 )fsd->fgrad=variational;
    }
/*------------------------------ number of iteration steps with fsd  ---*/
    frint(   "FSD_NUMITER"    ,&(fsd->numiter ) ,&ierr);
/*--------------------------------------- erforderliche genauigkeit  ---*/
    frdouble("FSD_ACC"        ,&(fsd->acc     ) ,&ierr);
/*--------------------- genauigkeit fuer die gleichheitsrestriktion  ---*/
    frdouble("FSD_ALPHA"      ,&(fsd->alpha   ) ,&ierr);
/*--------------------------------------------- schreitweitenfaktor  ---*/
    frdouble("FSD_BETA"       ,&(fsd->beta    ) ,&ierr);
/*----------------------------------------- schrittweitenbegrenzung  ---*/
    frdouble("FSD_DELTA"      ,&(fsd->delta   ) ,&ierr);
/*------------------------------------- gleichheitsrestriktionswert  ---*/
    frdouble("FSD_GAMMA"      ,&(fsd->gamma   ) ,&ierr);
/*----------------------------------------------------------------------*/
  }
/*---------------------------------------------------- fsd-strategy ----*/

/*---------------------- smoothing of objectives or densities or ... ---*/
    frchar(  "OPT_SMO"   ,buffer          ,&ierr);
    if (ierr==1)
    {
      if (strncmp(buffer,"ON",2)==0 )
      {
        opt->optsmooth=sm_on;
        /*----*/
        frread(); /* read line */
        frchar(  "SMO_TYPE"   ,buffer          ,&ierr);
        if (ierr==1)
        {
          if (strncmp(buffer,"gradient",8)==0 ) opt->smoothtype = sm_grad;
        }
        /*----*/
        frread(); /* read line */
        frdouble("SMO_ERAD"   ,&(opt->smoothrad) ,&ierr);
        frread(); /* read line */
        frdouble("SMO_EXPO"   ,&(opt->smoothexp) ,&ierr);
        /*----*/
        
      }
    }
/*---------------------- smoothing of objectives or densities or ... ---*/

/*-------------------------------------------------------- objective ---*/
  frchar(  "OPT_OBJ"   ,buffer          ,&ierr);
  if (ierr==1)
  {
          if (strncmp(buffer,"Stiffness",9)==0) opt->objective=oj_strain_energy;
     else if (strncmp(buffer,"Frequency",9)==0) opt->objective=oj_frequency;
  }
/*---------- input of all variable and restriction optimization data ---*/
  frint(   "OV_VAR",&(iloop) ,&ierr); /* number of variables */
  if (ierr==1)
  {
    opt->numvar = iloop; /* store n.of.variables    */
    /* allocate memory */
    opt->ovar = (OBTV*)calloc(opt->numvar,sizeof(OBTV));
    if (opt->ovar==NULL) dserror("Allocation for OPT.VARIABLES failed");
    
    for (i=0;i<iloop;i++)
    {
      frread(); /* read line */
      frint(   "ELEofMAT",&(idum) ,&ierr); 
      if(ierr==1)
      { /* ELEofMAT 1 TYPE DENS  DL 1.0E-1  DU 1.0E0 SCL 1. */
        opt->ovar[ccdv].objId  =  idum;
        opt->ovar[ccdv].ovatt  =  eleofmat;
        opt->ovar[ccdv].numvar =  1;
        frdouble("DL"   ,&(opt->ovar[ccdv].blower) ,&ierr); 
        frdouble("DU"   ,&(opt->ovar[ccdv].bupper) ,&ierr); 
        frdouble("SCL"  ,&(opt->ovar[ccdv].scva  ) ,&ierr); 
        ccdv++;
      }
    }
  }
/*-------------------------------------- global equality constraints ---*/
/*
OC_ECO 1                  : global equality constraints
              DFACE 1 TYP Volume  VAL 10.27899 CON 1 */   
  frint(   "OC_ECO",&(iloop) ,&ierr); /* number of constraints */
  if(ierr==1) 
  {
    opt->numeqc = iloop;           /* store  global equality constraints */
    
    /* allocate memory */
    opt->oeqc = (OEQC*)calloc(opt->numeqc,sizeof(OEQC));
    if (opt->oeqc==NULL) dserror("Allocation for GLOB.EQU.CON. failed");
    
    for (i=0; i<iloop; i++)
    {
      frread(); /* read line */
      /* */
      frint(   "DFACE",&(idum) ,&ierr); /* type of design object */ 
      if (ierr==1)
      {
        opt->oeqc[i].objId  = idum;
        opt->oeqc[i].oeqc_object = eqdface;
        
        frchar(  "TYP"   ,buffer          ,&ierr); /* type of constraint */
        if (strncmp(buffer,"Volume",6)==0 ) opt->oeqc[i].oeqc_type = volume;
        frdouble("VAL"   ,&(opt->oeqc[i].val) ,&ierr); /* initial ? */

      }
      /* */
      frint(   "DVOLU",&(idum) ,&ierr); /* type of design object */ 
      if (ierr==1)
      { /* DVOLU 1 TYP WEIGHT  SCL 1. */
        opt->oeqc[i].objId  = idum;
        opt->oeqc[i].oeqc_object = eqdvolu;
        
        frchar(  "TYP"   ,buffer          ,&ierr); /* type of constraint */
        if (strncmp(buffer,"WEIGHT",6)==0 ) opt->oeqc[i].oeqc_type = mass;
        frdouble("SCL"   ,&(opt->oeqc[i].scl) ,&ierr); /* scaling factor */

      }
    }
  }
/*----------------------------------------------------------------------*/
  frread();
} /* loop optimization data */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpctropt */

/*----------------------------------------------------------------------*/
#endif /* stop including optimization code to ccarat :*/
/*----------------------------------------------------------------------*/

/*! @} (documentation module close)*/
