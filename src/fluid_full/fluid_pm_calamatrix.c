/*!----------------------------------------------------------------------
\file
\brief projection method algorithm for fluid: calculation of A-matrix

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2_PRO
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fluid_prototypes.h"
#include "fluid_pm_prototypes.h"
static INT     veldis=0;           /* flag, if vel. discr. is used     */
static INT     predis=1;           /* flag, if pres. discr. is used    */
extern struct _ARRAY lmass_global; /* element lumped mass matrix       */
extern struct _ARRAY gradopr_global; /*element gradient operator       */
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
enum _CALC_ACTION calc_action[MAXFIELD];

/*!---------------------------------------------------------------------
\brief projection method: calculation of A-matrix

<pre>                                                        genk 11/02
control routine to calculate the A-matrix, gradient and reduced
gradient operator for solving the pressure problem of
the projection method.

According to Gresho it's necessary to compute the A-matrix with all
existing dofs (no boundary conditions)!

If we are parallel we have to renumber the dofs, since the function
mlpcg_precond_PtKP requires a very special dof numbering.
Afterwards the elements are called and the gradient matrix C, and the
lumped mass matrix is calculated.
Both are stored in CSR-format
Lumped mass is then copied to a redundant vector.
The matrix matrix matrix product is called: A = C_trans * M_l * C.
The gradient matrix is copied to OLL format and the reduced gradient
matrix is also created. Finally The A-matrix is copied to OLL-format.

</pre>

\param  *actfield         FIELD	         (i) the actual field
\param  *actintra         INTRA	         (i) the actual intra communicator
\param  *actpart          PARTITION      (i) the actual partition
\param  *actsolv          SOLVAR         (i) the actual solver
\param  *amatrix_oll      SPARSE_ARRAY   (o) A-matrix OLL format
\param  *fdyn             FLUID_DYNAMIC  (-) fluid dynamic parameters
\param  *gradmatrix_oll   OLL            (o) gradient matrix (OLL-format)
\param  *rgradmatrix_oll  OLL            (o) reduced gradient matrix (OLL-format)
\param  *lmass            DOUBLE         (o) lumped mass matrix (redundant)
\param  *numeq_total      INT            (i) total number of eq.
\param  *numeq            INT            (i) number of eq. on this proc
\param  *numeq_full       INT            (i) number of eq. on this proc	for full system
\param  *numeq_total_full INT            (i) total number of eqs for full system
\param  *action           CALC_ACTION    (i) action
\param  *container        CONTAINER      (-) container
\return void

------------------------------------------------------------------------*/
void fluid_pm_calamatrix(
                          FIELD           *actfield,
			  INTRA           *actintra,
			  PARTITION       *actpart,
			  SOLVAR          *actsolv,
			  SPARSE_ARRAY    *amatrix_oll,
			  FLUID_DYNAMIC   *fdyn,
			  OLL             *gradmatrix_oll,
			  OLL             *rgradmatrix_oll,
			  DOUBLE          *lmass,
			  INT             *numeq_total,
			  INT             *numeq,
			  INT             *numeq_full,
	 		  INT             *numeq_total_full,
			  CALC_ACTION     *action,
			  CONTAINER       *container
			 )
{
#ifdef MLPCG
INT     firstdof[2];     /* first dofnumber on this proc (vel&pres)    */
INT     lastdof[2];      /* last dofnumber on this proc  (vel&pres)    */
INT     nnz_guess;       /* nnz guess for opening CSR-matrix           */
INT     numnp_vel;       /* number of velocity nodes                   */
DOUBLE  t0,t1;           /* time measuring                             */
DBCSR  *work_csr;        /* working matrix for the matrix product      */
DBCSR  *amatrix_csr;     /* A-matrix in CSR format                     */
DBCSR  *lumpedmass_csr;  /* lumped mass matrix in CSR-format           */
DBCSR  *gradmatrix_csr;  /* gradient matrix in CSR-format              */
#ifdef PARALLEL
INT     mone=-1;
ARRAY   saveveldof_a;    /* vel dofs before renumbering                */
ARRAY   savepredof_a;    /* pre dofs before renumbering                */
ARRAY   newveldof_a;     /* vel dofs after renumbering                 */
ARRAY   newpredof_a;     /* pre dofs after renumbering                 */
ARRAY   velpointer_a;    /* index vector to old pressure dofs	       */
ARRAY   prepointer_a;	 /* index vector to old velocity dofs	       */
INT    *saveveldof;      /* old velocity dofs                          */
INT    *savepredof;      /* old pressure dofs                          */
INT    *newveldof;       /* new velocity dofs                          */
INT    *newpredof;       /* new pressure dofs                          */
INT    *prepointer;      /* index vector to old pressure dofs          */
INT    *velpointer;	 /* index vector to old velocity dofs          */
#endif

#ifdef DEBUG
dstrc_enter("fluid_pm_calamatrix");
#endif

/*-------------------------- allocate the lumped mass matrix CSR-format */
lumpedmass_csr=(DBCSR*)CCACALLOC(1,sizeof(DBCSR));

/*-------------------------------- allocate a working matrix CSR-format */
work_csr=(DBCSR*)CCACALLOC(1,sizeof(DBCSR));

/*--------------------------------- allocate gradient matrix CSR-format */
gradmatrix_csr=(DBCSR*)CCACALLOC(1,sizeof(DBCSR));

/*------------------------------------ allocate the A-matrix CSR-format */
amatrix_csr=(DBCSR*)CCACALLOC(1,sizeof(DBCSR));

numnp_vel = actfield->dis[veldis].numnp;

#ifdef PARALLEL
/* pressure arrays are also allocated with numnp of veldis, since the
   reference is created by the global node Id                           */
saveveldof = amdef("saveveldof",&saveveldof_a,numeq_total_full[veldis],1,"IV");
savepredof = amdef("savepredof",&savepredof_a,numnp_vel,1,"IV");
newveldof  = amdef("newveldof" ,&newveldof_a ,numeq_total_full[veldis],1,"IV");
newpredof  = amdef("newpredof" ,&newpredof_a ,numnp_vel,1,"IV");
velpointer = amdef("velpointer",&velpointer_a,numeq_total_full[veldis],1,"IV");
prepointer = amdef("prepointer",&prepointer_a,numeq_total_full[predis],1,"IV");
aminit(&saveveldof_a,&mone);
aminit(&savepredof_a,&mone);
aminit(&newveldof_a,&mone);
aminit(&newpredof_a,&mone);
#endif

/*----------------------------------------------------- set some values */
lumpedmass_csr->numeq	     = numeq_full[veldis];
lumpedmass_csr->numeq_total  = numeq_total_full[veldis];
gradmatrix_csr->numeq        = numeq_full[veldis];
gradmatrix_csr->numeq_total  = numeq_total_full[veldis];
work_csr->numeq              = numeq_full[predis];
work_csr->numeq_total        = numeq_total_full[predis];
amatrix_csr->numeq	     = numeq_full[predis];
amatrix_csr->numeq_total     = numeq_total_full[predis];

/*--------------------------------------------------- renumber the dofs */
fluid_pm_dofrenumber(actfield,
#ifdef PARALLEL
                    saveveldof,savepredof,newveldof,newpredof,
                    velpointer,prepointer,
#endif
		    actintra,actpart,
                    amatrix_csr,lumpedmass_csr,numeq_total_full,numeq_full,
		    numnp_vel,firstdof,lastdof,0);
#ifdef PARALLEL
amdel(&newveldof_a);
amdel(&newpredof_a);
#endif

/*------------------------------------------ open the lumpedmass matrix */
nnz_guess = 3*numeq_total_full[veldis];
mlpcg_csr_open(lumpedmass_csr,firstdof[veldis],lastdof[veldis],
               lumpedmass_csr->numeq_total,nnz_guess,actintra);
mlpcg_csr_zero(lumpedmass_csr,actintra);

/*------------------------------------------------- open the gradmatrix */
mlpcg_csr_open(gradmatrix_csr,firstdof[veldis],lastdof[veldis],
               gradmatrix_csr->numeq_total,2*nnz_guess,actintra);
mlpcg_csr_zero(gradmatrix_csr,actintra);

/*--------------------------------------------- open the working matrix */
mlpcg_csr_open(work_csr,firstdof[predis],lastdof[predis],
               work_csr->numeq_total,2*nnz_guess,actintra);
mlpcg_csr_zero(work_csr,actintra);

/*--------------------------------------------------- open the A-matrix */
nnz_guess = 3*numeq_total_full[predis];
mlpcg_csr_open(amatrix_csr,firstdof[predis],lastdof[predis],
               amatrix_csr->numeq_total,nnz_guess,actintra);
mlpcg_csr_zero(amatrix_csr,actintra);

/*------------------------------ call the elements and calculate lumped
                                 mass matrix and the gradient matrices  */
if (par.myrank==0) printf("  -> calling elements ...");
container->nii = 0;
container->gradmatrix  = gradmatrix_csr;
container->lumpedmass  = lumpedmass_csr;
container->actndis = veldis;
fdyn->pro_calmat = 1;
fdyn->pro_calrhs = 0;
fdyn->pro_calveln = 1;
fdyn->pro_mvv = 1;
fdyn->pro_kvv = 0;
fdyn->pro_gra = 1;
fdyn->pro_lum = 1;
t0 = ds_cputime();
calelm(actfield,actsolv,actpart,actintra,-1,-1,container,action);
t1 = ds_cputime()-t0;
if (par.myrank==0) printf("             {%10.3E}\n",t1);

/*---------------------------------- close grad- and lumped mass matrix */
mlpcg_csr_close(lumpedmass_csr);
mlpcg_csr_close(gradmatrix_csr);

/*------------------------------------------ invert reduced lumped mass */
fluid_pm_lumpedmass(lumpedmass_csr,actintra,
                    lmass,
#ifdef PARALLEL
		    velpointer,saveveldof,
#endif
		    numeq_total);

/*----------------------------------------------- perform A = CT*ML-1*C */
if (par.myrank==0) printf("  -> evaluate the matrix product ...");
t0 = ds_cputime();
mlpcg_precond_PtKP(gradmatrix_csr,lumpedmass_csr,amatrix_csr,
                   work_csr,NULL,0,actintra);
t1 = ds_cputime()-t0;
if (par.myrank==0) printf("  {%10.3E}\n",t1);

/*------------------------------------------------------ close A-matrix */
mlpcg_csr_close(amatrix_csr);

/*---------------------- copy Amatrix CSR-format to Amatrix OLL and
                         and reduce the pressure boundary conditions    */
if (par.myrank==0) printf("  -> copying matrices ...");
t0 = ds_cputime();
fluid_pm_redcpmat(
#ifdef PARALLEL
		 savepredof,prepointer,savepredof,prepointer,
#endif
                 numeq_total[predis],numeq_total[predis],
		 numnp_vel,numnp_vel,
		 amatrix_csr,amatrix_oll->oll);

/*----------------- copy gradmatrix CSR-format to gradmatrix OLL format */
fluid_pm_redcpmat(
#ifdef PARALLEL
		 saveveldof,velpointer,savepredof,prepointer,
#endif
                 numeq_total_full[veldis],numeq_total_full[predis],
		 numeq_total_full[veldis],numnp_vel,
		 gradmatrix_csr,gradmatrix_oll);

/*--------- copy gradmatrix CSR-format to reduced gradmatrix OLL format */
fluid_pm_redcpmat(
#ifdef PARALLEL
		 saveveldof,velpointer,savepredof,prepointer,
#endif
                 numeq_total[veldis],numeq_total[predis],
                 numeq_total_full[veldis],numnp_vel,
		 gradmatrix_csr,rgradmatrix_oll);
t1 = ds_cputime()-t0;
if (par.myrank==0) printf("             {%10.3E}\n",t1);

/*--------------------------------------------------- renumber the dofs */
fluid_pm_dofrenumber(actfield,
#ifdef PARALLEL
                    saveveldof,savepredof,NULL,NULL,
                    velpointer,prepointer,
#endif
		    actintra,actpart,
                    NULL,NULL,NULL,NULL,
		    numnp_vel,NULL,NULL,1);

/*------------------------------------------ destroy temporary matrices */
mlpcg_csr_destroy(amatrix_csr);
mlpcg_csr_destroy(lumpedmass_csr);
mlpcg_csr_destroy(gradmatrix_csr);
mlpcg_csr_destroy(work_csr);
CCAFREE(amatrix_csr);
CCAFREE(lumpedmass_csr);
CCAFREE(work_csr);
CCAFREE(gradmatrix_csr);
#ifdef PARALLEL
amdel(&saveveldof_a);
amdel(&savepredof_a);
amdel(&velpointer_a);
amdel(&prepointer_a);
#endif

/*----------------------------------------------------------------------*/
#else
  dserror("MLPCG needed for fluid projection method, but not defined!!");
#endif

#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of fluid_promethod_cal_amatrix */


/*!---------------------------------------------------------------------
\brief renumbering the dofs

<pre>                                                        genk 11/02

control routine to renumber the dofs

</pre>

\param  *actfield         FIELD	      (i)     the actual field
\param  *saveveldof       INT         (i/o)   actual velocity dofs
\param  *savepredof       INT         (i/o)   actual pressure dofs
\param  *newveldof        INT         (i/o)   velocity dofs after renumbering
\param  *newpredof        INT         (i/o)   pressure dofs after renumbering
\param  *velpointer       INT         (i/o)   pointer to old vel dofs
\param  *velpointer       INT         (i/o)   pointer to old pre dofs
\param  *actintra         INTRA	      (-)     the actual intra communicator
\param  *actpart          PARTITION   (-)     the actual partition
\param  *amatrix_csr      DBCSR       (-)     gradient matrix (CSR-format)
\param  *lumpedmass       DBCSR       (-)     lumped mass matrix (CSR-format)
\param  *numeq_total_full INT         (i)     total number of eq.
\param  *numeq_full       INT         (i)     number of eq. on this proc
\param   numnp_vel        INT         (i)     number of vel nodes
\param  *firsdof          INT         (o)     first dofnumbers on this proc
\param  *lastdof          INT         (o)     last dofnumbers on this proc
\param   option           INT         (i)     evaluation flag
\return void
\warning this is all in progress

------------------------------------------------------------------------*/
void fluid_pm_dofrenumber(FIELD      *actfield,
#ifdef PARALLEL
                          INT        *saveveldof,
			  INT        *savepredof,
			  INT        *newveldof,
			  INT        *newpredof,
                          INT        *velpointer,
			  INT        *prepointer,
#endif
			  INTRA      *actintra,
			  PARTITION  *actpart,
			  DBCSR      *amatrix_csr,
			  DBCSR      *lumpedmass,
			  INT        *numeq_total_full,
			  INT        *numeq_full,
			  INT         numnp_vel,
			  INT        *firstdof,
			  INT        *lastdof,
			  INT         option
			)
{
#ifdef PARALLEL
INT            Id;            /* global node Id                         */
INT            i,j,k;         /* simply some counters                   */
INT            fd,ld;         /* temporary first and last dof           */
INT            imyrank;       /* actual proc number                     */
INT            inproc;        /* number of procs                        */
INT            actdof;        /* actual dof                             */
DISCRET       *actdis;        /* actual discretisation                  */
PARTDISCRET   *actpdis;       /* actual pdis                            */
NODE          *actnode;       /* actual node                            */
#endif

#ifdef DEBUG
dstrc_enter("fluid_pmdofrenumber");
#endif


switch (option)
{
case 0: /* define new node numbers and save the old ones                */
   /*-------------------- save old dofs for vel and pres discretisation */
#ifdef PARALLEL
   actdis = &(actfield->dis[veldis]);
   for (i=0;i<actdis->numnp;i++)
   {
      actnode = &(actdis->node[i]);
      Id = actnode->Id;
      saveveldof[2*Id]   = actnode->dof[0];
      saveveldof[2*Id+1] = actnode->dof[1];
   }
   actdis = &(actfield->dis[predis]);
   for (i=0;i<actdis->numnp;i++)
   {
      actnode = &(actdis->node[i]);
      savepredof[actnode->Id] = actnode->dof[0];
   }

   /*----------------------- renumber the dofs for both discretisations */
   if (par.nprocs>1)
   {
      imyrank     = actintra->intra_rank;
      inproc      = actintra->intra_nprocs;
      actpdis     = &(actpart->pdis[veldis]);
      fluid_pm_newdofs(imyrank,inproc,actfield,actpdis,actintra,
                       lumpedmass,veldis);
      actpdis     = &(actpart->pdis[predis]);
      fluid_pm_newdofs(imyrank,inproc,actfield,actpdis,actintra,
                       amatrix_csr,predis);
   }

   /*-------------------- save new dofs for vel and pres discretisation */
   actdis = &(actfield->dis[veldis]);
   for (i=0;i<actdis->numnp;i++)
   {
      actnode = &(actdis->node[i]);
      Id = actnode->Id;
      newveldof[2*Id]   = actnode->dof[0];
      newveldof[2*Id+1] = actnode->dof[1];
   }
   actdis = &(actfield->dis[predis]);
   for (i=0;i<actdis->numnp;i++)
   {
      actnode = &(actdis->node[i]);
      newpredof[actnode->Id] = actnode->dof[0];
   }

   /*-------------------------------------------- create pointer vector */
   for (i=0;i<numeq_total_full[veldis];i++)
   {
      actdof = newveldof[i];
      if (actdof==-1) continue;
      velpointer[actdof] = i;
   }
   for (i=0;i<numnp_vel;i++)
   {
      actdof = newpredof[i];
      if (actdof==-1) continue;
      prepointer[actdof] = i;
   }
   /*---------------------------------------- find firstdof and lastdof */
   for (i=0;i<2;i++)
   {
      fd = numeq_total_full[i]+1;
      ld = -1;
      actpdis = &(actpart->pdis[i]);
      for (j=0;j<actpdis->numnp;j++)
      {
         actnode = actpdis->node[j];
         for (k=0;k<actnode->numdf;k++)
         {
            if (actnode->dof[k] > ld) ld = actnode->dof[k];
            if (actnode->dof[k] < fd) fd = actnode->dof[k];
         }
      }
      firstdof[i] = fd;
      lastdof[i]  = ld;
   }
#else
   firstdof[veldis]=0;
   firstdof[predis]=0;
   lastdof[veldis]=numeq_total_full[veldis]-1;
   lastdof[predis]=numeq_total_full[predis]-1;
#endif
break;
case 1: /* write old dof numbers to the nodes again */
#ifdef PARALLEL
   actdis = &(actfield->dis[veldis]);
   for (i=0;i<actdis->numnp;i++)
   {
      actnode = &(actdis->node[i]);
      actdof = actnode->dof[0];
      fd = velpointer[actdof];
      dsassert(fd>=0,"problems while renumbering dofs!\n");
      actnode->dof[0] = saveveldof[fd];
      actdof = actnode->dof[1];
      fd = velpointer[actdof];
      dsassert(fd>=0,"problems while renumbering dofs!\n");
      actnode->dof[1] = saveveldof[fd];
   }
   actdis = &(actfield->dis[predis]);
   for (i=0;i<actdis->numnp;i++)
   {
      actnode = &(actdis->node[i]);
      actdof = actnode->dof[0];
      fd = prepointer[actdof];
      dsassert(fd>=0,"problems while renumbering dofs!\n");
      actnode->dof[0] = savepredof[fd];
   }
#endif
break;
default:
   dserror("option out of range: don't know what to do\n");
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of fluid_pmdofrenumber */

/*!---------------------------------------------------------------------
\brief assembly of lumped mass and gradient matrix

<pre>                                                        genk 11/02

In this routine the lumped mass and the gradient matrix are assembled.
Both matrices are in CSR-format and the mask for the matrices is
created during the assembly.
This is a common assembly routine. However, instead of adding up the
values on our own, we use here the routine  mlpcg_csr_addentry(),
which blows up the CSR-matrix if necessary.

</pre>

\param  *container    CONTAINER      container                       (-)
\param  *actvele      ELEMENT        actual element (vel. discr.)    (-)
\param  *actpele      ELEMENT        actual element (pre. discr.)    (-)
\param  *actintra     INTRA	     the actual intra communicator   (-)
\warning coupled dofs not possible!!!
\warning lmass_global and gradopr_global are external global ARRAYS
\return void
------------------------------------------------------------------------*/
void assemble_fluid_amatrix(
                             CONTAINER *container,
                             ELEMENT   *actvele,
                             ELEMENT   *actpele,
			     INTRA     *actintra
			   )
{
#ifdef MLPCG
INT     i,j;                  /* simply some counters                   */
INT     ndv,ndp;              /* number of dofs per element             */
INT     counter;              /* a counter                              */
INT     numnpv,numnpp;        /* number of nodes per element            */
INT     lmv[MAXDOFPERELE];    /* vel location matrix                    */
INT     lmp[MAXDOFPERELE];    /* pre location matrix                    */
INT     rindex;               /* row index                              */
INT     cindex;               /* column index                           */
#ifdef PARALLEL
INT     ownerv[MAXDOFPERELE]; /* owner vel dofs                         */
#endif
INT     myrank;               /* proc number                            */
DOUBLE  val;                  /* value to assemble                      */
DOUBLE **elmass;              /* element mass matrix                    */
DOUBLE **egraop;              /* element gradient operator              */
DBCSR   *lumpedmass;          /* lumped mass matrix in CSR-format       */
DBCSR   *gradmatrix;          /* gradient operator in CSR-format        */
NODE    *actnode;             /* actual node                            */


#ifdef DEBUG
dstrc_enter("assemble_fluid_amatrix");
#endif


/*---------------------------------------- set some values and pointers */
/*  ATTENTION: lmass_global and gradopr_global are external global
               ARRAYS in this file                                      */
myrank       = actintra->intra_rank;
elmass       = lmass_global.a.da;
egraop       = gradopr_global.a.da;
lumpedmass   = container->lumpedmass;
gradmatrix   = container->gradmatrix;
numnpv       = actvele->numnp;
numnpp       = actpele->numnp;
ndv          = numnpv*actvele->node[0]->numdf;
ndp          = numnpp*actpele->node[0]->numdf;

/*----------------------------------------------- make location vectors */
counter=0;
for (i=0;i<numnpv;i++)
{
   actnode = actvele->node[i];
   for (j=0;j<actnode->numdf;j++)
   {
     lmv[counter] = actnode->dof[j];
#ifdef PARALLEL
     ownerv[counter] = actvele->node[i]->proc;
#endif
     counter++;
   }
}
dsassert(counter==ndv,"assemblage failed due to wrong dof numbering");
counter=0;
for (i=0;i<numnpp;i++)
{
   actnode = actpele->node[i];
   for (j=0;j<actnode->numdf;j++)
   {
     lmp[counter] = actnode->dof[j];
     counter++;
   }
}
dsassert(counter==ndp,"assemblage failed due to wrong dof numbering");

/*-----------------------------------------------------------------------*
 | assemble the lumped mass matrix:					 |
 | the lumped mass matrix is a						 |
 | diagonal matrix, so one loop is sufficient				 |
 | WARNING: COUPLED DOFS NOT POSSIBLE!!!!!!!!!!!!!!!!!!!		 |
 *-----------------------------------------------------------------------*/
for (j=0;j<ndv;j++)
{
   /*-------------------------------------------- loop only my own rows */
#ifdef PARALLEL
   if (ownerv[j]!=myrank) continue;
#endif
   rindex = lmv[j];
   val = elmass[j][j];
   mlpcg_csr_addentry(lumpedmass, val, rindex, rindex, actintra);
}

/*-----------------------------------------------------------------------*
 | assemble the gradient matrix:					 |
 | the gradient matrix is a rectangular matrix				 |
 | with size veldofs x presdofs 			                 |
 | WARNING: COUPLED DOFS NOT POSSIBLE!!!!!!!!!!!!!!!!!!!		 |
 *-----------------------------------------------------------------------*/

/*----------------------------------------- loop over the element row i */
for (i=0;i<ndv;i++)
{
  /*-------------------------------------------- loop only my own rows */
#ifdef PARALLEL
   if (ownerv[i]!=myrank) continue;
#endif
   rindex = lmv[i];
   /*----------------------------------- loop over the element column j */
   for (j=0;j<ndp;j++)
   {
      cindex = lmp[j];
      val = egraop[i][j];
      mlpcg_csr_addentry(gradmatrix, val, rindex, cindex, actintra);
   }
}

/*----------------------------------------------------------------------*/
#else
  dserror("MLPCG needed for fluid projection method, but not defined!!");
#endif

#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of assemble_fluid_amatrix */


/*!---------------------------------------------------------------------
\brief inverse of lumped mass

<pre>                                                        genk 10/03

In this function the lumped mass matrix (CSR) is inverted. Afterwards
non supported values are copied to a redundant vector.

</pre>
\param *lumpedmass_CSR  DBCSR   (i/o)  lumped mass matrix CSR-format
\param *actintra        INTRA   (i)    intra-communicator
\param *lmass           DOUBLE  (o)    lumped mass matrix (redundant vector)
\param *velpointer      INT     (i)    index vector to old vel dofs
\param *saveveldof      INT     (i)    old vel dofs
\param *numeq_total     INT     (i)    total number of equations
\return void

------------------------------------------------------------------------*/
void fluid_pm_lumpedmass(DBCSR *lumpedmass_csr,INTRA *actintra,
                         DOUBLE *lmass,
#ifdef PARALLEL
			 INT *velpointer, INT *saveveldof,
#endif
			 INT *numeq_total)
{
INT         i;                /* simply some counters                   */
INT         numeq;            /* number of equations on this proc       */
DOUBLE     *val;              /* matrix entries                         */
#ifdef PARALLEL
INT         j,k;              /* simply some counters                   */
INT        *update;           /* dofs updated on this proc              */
INT         andof,aodof,fd;   /* actual dof numbers                     */
ARRAY       lmsend_a;
DOUBLE     *lmsend;           /* receive buffer for lumped mass         */
#endif

#ifdef DEBUG
dstrc_enter("fluid_pm_lumpedmass");
#endif

/*--------------------------------------- invert full lumped mass matrix*/
numeq=lumpedmass_csr->numeq;
val = lumpedmass_csr->a.a.dv;
for (i=0;i<numeq;i++) val[i] = 1/val[i];

/*------------- store inverse of lumped mass matrix in redundant vector */
#ifdef PARALLEL
/*----------------------------------------------------------------------*
 | loop all entries on this proc      			        	|
 | get old dofnumber                                                    |
 | skip supported dofs                                                  |
 | store actual value in lmsend (old dofnumbers!!!)			|
 | allreduce to lmass                                                   |
 *----------------------------------------------------------------------*/
lmsend=amdef("lmsend",&lmsend_a,numeq_total[veldis],1,"DV");
amzero(&lmsend_a);
update = lumpedmass_csr->update.a.iv;
for (j=0;j<numeq;j++)
{
   andof=update[j];
   fd   = velpointer[andof];
   aodof=saveveldof[fd];
   if (aodof>=numeq_total[veldis]) continue;
   lmsend[aodof] = val[j];
}
/*------------------------------------ allreduce the lumped mass vector */
MPI_Allreduce(lmsend,lmass,numeq_total[veldis],
              MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
amdel(&lmsend_a);
#else
for (i=0;i<numeq_total[veldis];i++) lmass[i]=val[i];
#endif

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
}/*end of void fluid_pm_lumpedmass*/

/*!---------------------------------------------------------------------
\brief copying/reducing CSR to OLL

<pre>                                                        genk 10/03

Up to now, we calculated the A-matrix / gradmatrix including all
existing pressure dofs. However the system during the timeloop
has to be solved without supported dofs.
Therefore the matrices are reduced in this function. In addition they
are copied from the CSR to OLL format, so that one can use any
solver for symmetric systems for getting the pressure solution.

</pre>
\param  *saverowdof        INT     (i)   row index of dofs OLL-format
\param  *rowpointer        INT     (i)   index vector to row dofs
\param  *saverowdof        INT     (i)   column index of dofs OLL-format
\param  *colpointer        INT     (i)   index vector to col dofs
\param   numeq_totalr      INT     (i)   total number of row equations
\param   numeq_totalc      INT     (i)   total number of column equations
\param   numeq_totalr_full INT     (i)   total number of row equations (full)
\param   numeq_totalc_full INT     (i)   total number of column equations (full)
\param  *matrix_csr        DBCSR   (i)   matrix CSR format (full)
\param  *matrix_oll        OLL     (o)   matrix OLL-format (reduced)
\return void

------------------------------------------------------------------------*/
void fluid_pm_redcpmat(
#ifdef PARALLEL
		       INT *saverowdof,      INT *rowpointer,
                       INT *savecoldof,      INT *colpointer,
#endif
		       INT numeq_totalr,     INT numeq_totalc,
		       INT numeq_totalr_full,INT numeq_totalc_full,
		       DBCSR *matrix_csr,    OLL *matrix_oll)
{
INT      i,j;             /* simply some counters                       */
INT     *ia,*ja,*update;  /* matrix indices                             */
INT      numeq;           /* number of equations on this proc           */
INT      rindex,cindex;   /* row/column index                           */
INT      counter;
INT     *lm;              /* column location vector                     */
DOUBLE  *a,*val;          /* matrix entries                             */
ARRAY    lm_a;
ARRAY    val_a;

#ifdef PARALLEL
INT      k;               /* simply some counters                       */
#endif

#ifdef DEBUG
dstrc_enter("fluid_pmredcpam");
#endif

/*---------------------------------------------------------- initialise */
ja        = matrix_csr->ja.a.iv;
ia        = matrix_csr->ia.a.iv;
a         = matrix_csr->a.a.dv;
update    = matrix_csr->update.a.iv;
numeq     = matrix_csr->numeq;

/*------------------------------------------- allocate temporary arrays */
lm  = amdef("lm",&lm_a,numeq_totalc_full,1,"IV");
val = amdef("val",&val_a,numeq_totalc_full,1,"DV");

/*---------------------------------------------------- loop rows of CSR */
for (i=0;i<numeq;i++)
{
   rindex = update[i];
#ifdef PARALLEL
   /*------------------------------------ find old number of actual dof */
   k = rowpointer[rindex];
   rindex = saverowdof[k];
#endif
   /*------------------------------------------ check for supported dof */
   if (rindex>=numeq_totalr) continue;
   /*---------------------------------------------- loop columns of CSR */
   counter=0;
   for (j=ia[i];j<ia[i+1];j++)
   {
      cindex = ja[j];
#ifdef PARALLEL
      /*--------------------------------- find old number of actual dof */
      k = colpointer[cindex];
      cindex = savecoldof[k];
#endif
      /*--------------------------------------- check for supported dof */
      if (cindex>=numeq_totalc) continue;
      lm[counter] = cindex;
      val[counter] = a[j];
      counter++;
   }
   /*-------------------------------------------- add row to OLL matrix */
   oll_addrow(matrix_oll,rindex,lm,val,counter);
}

amdel(&lm_a);
amdel(&val_a);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
}/*end of void fluid_pmredcpam*/

/*!---------------------------------------------------------------------
\brief renumber dofs

<pre>                                                        genk 10/03

This routine renumbers the dofs in the nodes in a sense, that partition's
dofs are numbers first followed by the boundary dofs coupled to other
processor's equations.
In the one processor or sequentiell case this is not done

In this case, the A-matrix is calculated using ALL existing dofs, so
the supported dofs must not be excluded!!!

</pre>
\param myrank     INT         (i)   this processors intra rank
\param nproc      INT         (i)   number of processors in this intra-communicator
\param actfield   FIELD*      (i)   catual physical field (structure)
\param actpart    PARTITION*  (i)   this processors partition
\param actintra   INTRA       (i)   the intra-communicator of this field
\warning this routine renumbers the dofs and does NOT support coupling!!!!
\sa mlpcg_renumberdofs
\return void

------------------------------------------------------------------------*/
void fluid_pm_newdofs(INT            myrank,
                      INT            nproc,
                      FIELD         *actfield,
                      PARTDISCRET   *actpdiscret,
                      INTRA         *actintra,
                      DBCSR         *bdcsr,
		      INT            dis)
{
#ifdef PARALLEL
INT            i,j,k;
ELEMENT       *actele;
NODE          *actnode;
ARRAY          dofflag_a;
INT           *dofflag;
ARRAY          newdof_a;
INT           *newdof;
ARRAY          newdofrecv_a;
INT           *newdofrecv;
INT            locsize_send[MAXPROC];
INT            locsize[MAXPROC];
INT            startdof;
INT            savedof;
INT            actdof;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("fluid_pm_newdofs");
#endif

/*----------------------------------------------------------------------*/
dofflag = amdef("tmp",&dofflag_a,bdcsr->numeq_total,1,"IV");
amzero(&dofflag_a);
newdof = amdef("tmp",&newdof_a,bdcsr->numeq_total,1,"IV");
amzero(&newdof_a);
newdofrecv = amdef("tmp",&newdofrecv_a,bdcsr->numeq_total,1,"IV");
amzero(&newdofrecv_a);

/*----------------------------------------------------------------------*/
for (i=0; i<nproc; i++) locsize_send[i]=0;
locsize_send[myrank]=bdcsr->numeq;
MPI_Allreduce(locsize_send,locsize,nproc,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
startdof=0;
for (i=0; i<myrank; i++) startdof += locsize[i];
savedof = startdof;

/*---------------------------------------- make flags for boundary dofs */
for (i=0; i<actpdiscret->bou_numele; i++)
{
   actele = actpdiscret->bou_element[i];
   for (j=0; j<actele->numnp; j++)
   {
      actnode = actele->node[j];
      if (actnode->proc != myrank) continue;
      for (k=0; k<actnode->numdf; k++)
      {
         actdof = actnode->dof[k];
         if (actdof>=bdcsr->numeq_total) continue;
         dofflag[actdof]=1;
      }
   }
}

/*-------- give ascending numbers to the dofs, start with internal dofs */
for (i=0; i<actpdiscret->numnp; i++)
{
   actnode = actpdiscret->node[i];
   dsassert(actnode->proc==myrank,"Partitioning got mixed up");
   for (j=0; j<actnode->numdf; j++)
   {
      actdof = actnode->dof[j];
      if (dofflag[actdof]==1) continue;
      newdof[actdof] = startdof;
      startdof++;
   }
}

/* save the first dof which is interproc-coupled in the DBCSR_ROOT matrix */
bdcsr->firstcoupledof = startdof;

/*------ now number the dofs, which have interproc off-diagonal entries */
for (i=0; i<actpdiscret->numnp; i++)
{
   actnode = actpdiscret->node[i];
   for (j=0; j<actnode->numdf; j++)
   {
      actdof = actnode->dof[j];
      if (dofflag[actdof]!=1) continue;
      newdof[actdof] = startdof;
      startdof++;
   }
}

/*------------------------------------------------- check local dof sum */
if (startdof-savedof != bdcsr->numeq)
dserror("Local number of equations wrong");

/*--------------------------------------- allreduce the new dof numbers */
MPI_Allreduce(newdof,newdofrecv,bdcsr->numeq_total,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);

/*----------------------- now put new dofnumbers to the node structures */
for (i=0; i<actfield->dis[dis].numnp; i++)
{
   actnode = &(actfield->dis[dis].node[i]);
   for (j=0; j<actnode->numdf; j++)
   {
      actdof = actnode->dof[j];
      actnode->dof[j] = newdofrecv[actdof];
   }
}

/*----------------------------------------------------------------------*/
amdel(&dofflag_a);
amdel(&newdof_a);
amdel(&newdofrecv_a);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
#endif
return;
} /* end of fluid_pm_newdofs */
#endif
/*! @} (documentation module close)*/
