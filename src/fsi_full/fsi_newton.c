/*!-----------------------------------------------
\file
\brief functions for iter_stagg_Newton_ methods

<pre>
Maintainer: Markus Schmidberger
            schmidbe@in.tum.de
<\pre>
--------------------------------------------------*/

#ifndef CCADISCRET

#ifdef D_FSI
#ifdef FSI_NEWTONCOUPLING

#include "../headers/standardtypes.h"
#include "fsi_prototypes.h"


/*------------------------------------
 |fluid dynamic variables            |
 ------------------------------------*/
static FSI_DYNAMIC *fsidyn;


/*!----------------------------------------------------------------------
\brief file pointers
<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*--------------------------------------------------------------*
 | pointer to allocate dynamic variables if needed              |
 | dedfined in global_control.c                                 |
 | ALLDYNA               *alldyn;                               |
 *--------------------------------------------------------------*/
extern ALLDYNA      *alldyn;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*---------------------------------------------------------------
 | Initialisation of GMRES SOLVER,                              |
 | included by libraries -lewduni -lifcore -lirc                |
 ---------------------------------------------------------------*/
void ewdgmres_(int* icode, double* bg, double* dg, double* wkg1, double* wkg2, double* z, double* v, double* hh, double* yy, double* cc, double* ss, int* ninner, int* nouter, int* ndc, double* etolr, int* its, int* iit, char* aid);

FLUID_DYNAMIC  *fdyn;     /* fluid dynamic variables * */

/*!---------------------------------------------------------------------*
\brief Solve linear System:  R' delta d = -R [rS=-R]

<pre>
Solve linear System with GMRES nach Behr
  Problem: Approximation of Jacobimatrix
  a) R'=I
  b) R' with finite Differences

see DA: Markus Schmidberger
</pre>
------------------------------------------------------------------------*/
double *fsi_newton_SolveSystem(
    DOUBLE             *rS,
    FIELD              *structfield,
    FIELD              *fluidfield,
    FIELD              *alefield,
    FSI_STRUCT_WORK    *struct_work,
    FSI_FLUID_WORK     *fluid_work,
    FSI_ALE_WORK       *ale_work,
    INT                s_disnum_calc,
    INT                s_disnum_io,
    INT                f_disnum_calc,
    INT                f_disnum_io,
    INT                a_disnum_calc,
    INT                a_disnum_io,
    INT                itnum
  )
{
  INT 		       i, j;           /* some counters                                 */
  INT  		       numnp_total;   /* total number of nodes                         */
  INT    	       numdf_total;   /* total number of dofs                          */
  INT    	       numdf,dof;     /* actual number of dofs, actual dof             */
  INT  		       *sid;           /* structural interface dofs */
  NODE  	       *actsnode;      /* actual struct node */
  ARRAY_POSITION       *ipos;

  INT 		       icode=0, its=0;    /* Reverse communication switch, Counters */
  INT 		       iit=0, wert;
  INT 		       ninner, nouter;    /* Size of Krylov-Raum, Number of outer */
  DOUBLE 	       *dg,*wkg1, *wkg2, *z, *v, *hh, *yy, *cc, *ss; /*  Pointers to solver array */
  static  ARRAY        dg_a, wkg1_a, wkg2_a, z_a, v_a, hh_a, yy_a, cc_a, ss_a; /* Arrays for solver */
  DOUBLE 	       etolr;          /* Convergence criterium */
  CHAR 		       aid;             /* String for output */

  DOUBLE        alpha = 1e-4;            /* increment for finite Differenc */
  ARRAY         dtilde_a, zw_a;   /* Array : dtilde */
  DOUBLE        *dtilde, *zw;     /* Vector for dtilde */

  static double norm;

  INT    numff = genprob.numff;
  fdyn          = alldyn[numff].fdyn;
  double ittolSave = fdyn->ittol;

  numnp_total = structfield->dis[s_disnum_calc].numnp; /* Number of nodes in this structurefield */
  fsidyn = alldyn[3].fsidyn; /* Pointer to structure FSI-Dynamic */
  sid = fsidyn->sid.a.iv; /* Array with Structure-Interface DOFs */
  numdf_total = fsidyn->sid.fdim; /* Number of structure-Interface DOFs */
  ipos = &(structfield->dis[s_disnum_calc].ipos); /* meaning of this discretization's node array entries */

  printf("Solving linear System with GMRES");
  switch (fsidyn->ifsi)
  {

  case fsi_iter_stagg_Newton_I:
  {
    ninner = 1 ; nouter = 50 ; etolr= fsidyn->convtol;
    printf("(GMRESTol: %e)...\n", etolr);
    printf("# START #########################################################\n");
    printf("#################################################################\n");
    break;
  }

  case fsi_iter_stagg_Newton_FD:
  {
    /* set variables */
    ninner = floor(sqrt(sqrt(fsidyn->numsid))) ; nouter = 50 ; etolr= fsidyn->convtol;

    /* ALPHA */
    if (itnum == 0) {
      norm=0.0;
      for (i=0;i<numdf_total;i++) {
        norm=norm + (rS[i]*rS[i]);
      }
      norm=sqrt(norm);
    }
    alpha=norm;

    /* alpha=1e-4; */


    printf("(GMRES-Tol: %e, GMRES-Ninner: %d, Alpha:"YELLOW" %e"END_COLOR")...\n", etolr,ninner, alpha);
    printf("# START #########################################################\n");
    printf("#################################################################\n");


    /*loop nodes and reading displacement dtilde at interface */
    dtilde = amdef("dtilde",&dtilde_a,numdf_total,1,"DV");
    amzero(&dtilde_a);
    for (i=0;i<numnp_total;i++)
    {
      actsnode  = &(structfield->dis[s_disnum_calc].node[i]);
      numdf = actsnode->numdf;
      for (j=0;j<numdf;j++)
      {
        dof = actsnode->dof[j];
        dsassert(dof<numdf_total,"dofnumber not valid!\n");
        if (sid[dof]==0) continue;
        dtilde[dof] = actsnode->sol_mf.a.da[ipos->mf_dispnp][j];
      }
    }
    break;
  }
  default:
    dserror("Ups, no Newton Method");
  }

  aid = '|';
  dg = amdef("dg",&dg_a,numdf_total,1,"DV");
  wert=1; /* Initial solution guess */
  aminit(&dg_a, rS);
  wkg1 = amdef("wkg1",&wkg1_a,numdf_total,1,"DV");
  wkg2 = amdef("wkg2", &wkg2_a,numdf_total,1,"DV");
  z = amdef("z",&z_a,numdf_total*ninner, 1,"DV");
  v = amdef("v",&v_a,numdf_total*(ninner+1), 1,"DV");
  hh = amdef("hh",&hh_a,(ninner+1)*ninner, 1,"DV");
  yy = amdef("yy",&yy_a,ninner+1, 1,"DV");
  cc = amdef("cc",&cc_a,ninner,1, "DV");
  ss = amdef("ss",&ss_a,ninner,1,"DV");

  /* GMRES see M.Behr */
  gmres:
  ewdgmres_(&icode, rS, dg, wkg1, wkg2, z, v, hh, yy, cc, ss, &ninner, &nouter, &numdf_total, &etolr, &its, &iit, &aid);
  if (icode == 1)
  {
    /* Left-Prekonditioning onto wkg1, saving in wkg2 */
    for (j=0;j<numdf_total;j++)
      wkg2[j]=wkg1[j];
      goto gmres;
  }
  if (icode == 2)
  {
    /* Matrix-Vector-Product from A to wkg1, saving in wkg2 */
    switch (fsidyn->ifsi)
    {

    case fsi_iter_stagg_Newton_I:
      /*-----------------------------------------------------
      | Approximation by Identity Matrix                    |
      | == Block-Gauß-Seidel                                |
      ------------------------------------------------------*/
    {
      for (j=0;j<numdf_total;j++)
        wkg2[j]=wkg1[j];
      goto gmres;
    }

    case fsi_iter_stagg_Newton_FD:
      /*----------------------------------------------------
       | Approximation by Finite Differenc                  |
       ----------------------------------------------------*/
    {
      solserv_sol_copy(structfield, s_disnum_calc, node_array_sol_mf, node_array_sol_mf, 8,ipos->mf_dispnp);
      solserv_sol_copy(structfield, s_disnum_calc, node_array_sol_mf, node_array_sol_mf, 9,ipos->mf_reldisp);
      if (itnum > 0) solserv_sol_copy(structfield,s_disnum_calc,node_array_sol,node_array_sol,11, 9);
      if (itnum == 0) solserv_sol_copy(structfield,s_disnum_calc,node_array_sol,node_array_sol,12, 1);

      /*Calculation of S(F(d+alpha*wkg1)) */
      /* Verschiebung vorgeben, passend für ale Löser, d.h. in
       * struct_sol_mf[0]*/
      for (i=0;i<numnp_total;i++)
        {
          actsnode  = &(structfield->dis[s_disnum_calc].node[i]);
          numdf = actsnode->numdf;
          for (j=0;j<numdf;j++)
          {
            dof = actsnode->dof[j];
            dsassert(dof<numdf_total,"dofnumber not valid!\n");
            if (sid[dof]==0) continue;
            actsnode->sol_mf.a.da[ipos->mf_dispnp][j] += alpha*wkg1[dof];
          }
        }


      fsi_ale_calc(ale_work,alefield,a_disnum_calc,a_disnum_io,structfield,s_disnum_calc);

      /* one step for fluidsolver, problems with itermax=1 */
      fdyn->ittol = 10;
      fsi_fluid_calc(fluid_work,fluidfield,f_disnum_calc,f_disnum_io,alefield,a_disnum_calc);
      fdyn->ittol = ittolSave;

      /*itnum=1 da nur eine Iteration und immer wieder neu gestartet wird*/
      fsi_struct_calc(struct_work,structfield,s_disnum_calc,s_disnum_io,1,fluidfield,f_disnum_calc);


      /*Reading of the displacement*/
      zw = amdef("zw",&zw_a,numdf_total,1,"DV");
      amzero(&zw_a);
      for (i=0;i<numnp_total;i++)
      {
        actsnode  = &(structfield->dis[s_disnum_calc].node[i]);
        numdf = actsnode->numdf;
        for (j=0;j<numdf;j++)
        {
          dof = actsnode->dof[j];
          dsassert(dof<numdf_total,"dofnumber not valid!\n");
          if (sid[dof]==0) continue;
          zw[dof] = actsnode->sol_mf.a.da[ipos->mf_dispnp][j];
        }
      }

      /* Calculation of Finite Differenz */
      for (j=0;j<numdf_total;j++) {
        wkg2[j]=(dtilde[j]-zw[j]+alpha*wkg1[j])/alpha;;
      }
      amdel(&zw_a);
      goto gmres;
    }
    default:
      dserror("Ups, no Newton Method");
    }
  }

  /*free memory*/
  amdel(&wkg1_a);amdel(&wkg2_a);amdel(&z_a);amdel(&v_a);amdel(&hh_a);
  amdel(&yy_a);amdel(&cc_a);amdel(&ss_a);

  if (fsidyn->ifsi == fsi_iter_stagg_Newton_FD) amdel(&dtilde_a);

  printf("# ENDE ##########################################################\n");
  printf("#################################################################\n\n");


  return(dg);
}






/*!---------------------------------------------------------------------*
\brief Calculation of the right Side

<pre>

see DA: Markus Schmidberger
</pre>
------------------------------------------------------------------------*/
double *fsi_newton_rightSide(
  FIELD           *structfield,
  INT             s_disnum_calc
  )
{
  INT     i,j;           /* some counters                                 */
  INT     numnp_total;   /* total number of nodes                         */
  INT     numdf_total;   /* total number of dofs                          */
  INT     numdf,dof;     /* actual number of dofs, actual dof             */
  INT    *sid;           /* structural interface dofs */
  NODE   *actsnode;      /* actual struct node */
  ARRAY_POSITION  *ipos;
  static double *rS;/* Vector for right Side */
  static ARRAY   rS_a;   /* Array : right Side */

  numnp_total = structfield->dis[s_disnum_calc].numnp; /* Number of nodes in this structurefield */
  fsidyn = alldyn[3].fsidyn; /* Pointer to structure FSI-Dynamic */
  sid = fsidyn->sid.a.iv; /* Array with Structure-Interface DOFs */
  numdf_total = fsidyn->sid.fdim; /* Number of structure-Interface DOFs */
  ipos = &(structfield->dis[s_disnum_calc].ipos); /* meaning of this discretization's node array entries */

  rS = amdef("rS",&rS_a,numdf_total,1,"DV"); /* define vektor in structure ARRAY */
  amzero(&rS_a); /* assignment of the array */

  /* loop over all structure-nodes and calculate residuum at the interface */
  for (i=0;i<numnp_total;i++)
  {
    actsnode  = &(structfield->dis[s_disnum_calc].node[i]); /* actual node */
    /* check for coupling nodes */
    numdf = actsnode->numdf; /* Number of degrees of freedom (DOF) */
    /* loop over dofs */
    for (j=0;j<numdf;j++)
    {
      dof = actsnode->dof[j];/* DOF numbers*/
      dsassert(dof<numdf_total,"dofnumber not valid!\n");
      if (sid[dof]==0) continue;
      rS[dof] = actsnode->sol_mf.a.da[ipos->mf_dispnp][j] - actsnode->sol_mf.a.da[ipos->mf_reldisp][j]; /* dtilde-d */
    }
  }
  return(rS);
}





/*!---------------------------------------------------------------------*
\brief  Calculation of new displacement and writing into sol_mf[0][j]
<pre>

see DA: Markus Schmidberger
</pre>
------------------------------------------------------------------------*/
void fsi_newton_final(
  DOUBLE          *dg,
  FIELD           *structfield,
  INT             s_disnum_calc
  )
{
  INT     i,j;           /* some counters                                 */
  INT     numnp_total;   /* total number of nodes                         */
  INT     numdf_total;   /* total number of dofs                          */
  INT     numdf,dof;     /* actual number of dofs, actual dof             */
  INT    *sid;           /* structural interface dofs */
  NODE   *actsnode;      /* actual struct node */
  ARRAY_POSITION  *ipos;

  numnp_total = structfield->dis[s_disnum_calc].numnp;
  fsidyn = alldyn[3].fsidyn;
  sid = fsidyn->sid.a.iv;
  numdf_total = fsidyn->sid.fdim;
  ipos = &(structfield->dis[s_disnum_calc].ipos);

  for (i=0;i<numnp_total;i++)
    {
      actsnode  = &(structfield->dis[s_disnum_calc].node[i]);
      numdf = actsnode->numdf;
      for (j=0;j<numdf;j++)
      {
        dof = actsnode->dof[j];
        dsassert(dof<numdf_total,"dofnumber not valid!\n");
        if (sid[dof]==0) continue;
        actsnode->sol_mf.a.da[ipos->mf_dispnp][j] = actsnode->sol_mf.a.da[8][j] + dg[dof] ;
      }
    }
}


#endif
#endif
#endif
