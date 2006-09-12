/*-----------------------------------------------------------------------*/
/*!
\file
\brief coupling conditions for fsi problems


<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

 */
/*-----------------------------------------------------------------------*/


/*!
\addtogroup FSI
*//*! @{ (documentation module open)*/


#ifdef D_FSI
#ifndef FSI_NONMATCH

#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../fluid2/fluid2.h"
#include "../fluid3/fluid3.h"
#include "fsi_prototypes.h"


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;



/*-----------------------------------------------------------------------*/
/*!
  \brief create fsi coupling conditions

  In this routine the fsi-coupling conditions from the input file are
  evaluated and transformed into Neumann conditions for structure
  and into Dirichlet conditions for fluid and ale.

  \return void

  \author genk
  \date   09/02

 */
/*-----------------------------------------------------------------------*/
void fsi_createfsicoup()
{

  INT i,j;                                  /* simply some counters       */
  INT hasdirich,hascouple,hasfsi,hasneum;   /* different flags            */

  FIELDTYP  fieldtyp;

  DNODE    *actdnode;
  DLINE    *actdline;
  DSURF    *actdsurf;


#ifdef DEBUG
  dstrc_enter("fsi_createfsicoup");
#endif


  /* loop dsurfs */
  /* ----------- */
  for (i=0; i<design->ndsurf; i++)
  {
    hasdirich=0;
    hascouple=0;
    hasfsi   =0;
    hasneum  =0;
    actdsurf = &(design->dsurf[i]);


    /* check for conditions */
    if (actdsurf->dirich!=NULL) hasdirich++;
    if (actdsurf->couple!=NULL) hascouple++;
    if (actdsurf->fsicouple!=NULL) hasfsi++;
    if (actdsurf->neum!=NULL) hasneum++;



    /* surface is a fsi coupling surface */
    if (hasfsi!=0)
    {

      fieldtyp=actdsurf->fsicouple->fieldtyp;
      switch (fieldtyp)
      {

        case structure:
          dsassert(hasdirich==0,
              "dirich- and fsi-coupling condition defined on same DSURF\n");
          dsassert(hascouple==0,
              "coupling- and fsi-coupling condition defined on same DSURF\n");

          if (hasneum==0)
          {
            actdsurf->neum = (NEUM_CONDITION*)CCACALLOC(1,sizeof(NEUM_CONDITION));
            actdsurf->neum->neum_type = neum_FSI;
          }
          else if (hasneum>0 && actdsurf->neum->neum_type==neum_live)
          {
            /* combination of live load and FSI load */
            actdsurf->neum->neum_type = neum_live_FSI;
          }
          else if (hasneum>0 && actdsurf->neum->neum_type==neum_orthopressure)
          {
            /* combination of live load and FSI load */
            actdsurf->neum->neum_type = neum_opres_FSI;
          }
          else
            dserror("neumann- and fsi-couping condition defined on same DSRUF\n");
          break;


        case fluid:
          dsassert(hasdirich==0,
              "dirich- and fsi-coupling condition defined on same DSURF\n");
          dsassert(hascouple==0,
              "coupling- and fsi-coupling condition defined on same DSURF\n");
          dsassert(hasneum==0,
              "neumann- and fsi-coupling condition defined on same DSURF\n");


          /* allocate space for a dirichlet condition in this dsurf */
          actdsurf->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
          amdef("onoff",&(actdsurf->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
          amdef("val",&(actdsurf->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
          amdef("curve",&(actdsurf->dirich->curve),MAXDOFPERNODE,1,"IV");
          amdef("function",&(actdsurf->dirich->funct),MAXDOFPERNODE,1,"IV");
          amzero(&(actdsurf->dirich->dirich_onoff));
          amzero(&(actdsurf->dirich->dirich_val));
          amzero(&(actdsurf->dirich->curve));
          amzero(&(actdsurf->dirich->funct));

          /* initialise for fsi-coupling */
          for (j=0;j<genprob.ndim;j++)
            actdsurf->dirich->dirich_onoff.a.iv[j] = 1;
          actdsurf->dirich->dirich_type=dirich_FSI;


          /* also create the dirichlet condition for the ale */
          if (genprob.create_ale ==1)
          {
            /* allocate space for a dirichlet condition in this dsurf */
            actdsurf->ale_dirich =
              (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
            amdef("onoff",    &(actdsurf->ale_dirich->dirich_onoff),MAXDOFPERNODE, 1,"IV");
            amdef("val",      &(actdsurf->ale_dirich->dirich_val),  MAXDOFPERNODE, 1,"DV");
            amdef("curve",    &(actdsurf->ale_dirich->curve),       MAXDOFPERNODE, 1,"IV");
            amdef("function", &(actdsurf->ale_dirich->funct),       MAXDOFPERNODE, 1,"IV");
            amzero(&(actdsurf->ale_dirich->dirich_onoff));
            amzero(&(actdsurf->ale_dirich->dirich_val));
            amzero(&(actdsurf->ale_dirich->funct));
            amzero(&(actdsurf->ale_dirich->curve));

            /* initialise for fsi-coupling */
            for (j=0;j<genprob.ndim;j++)
              actdsurf->ale_dirich->dirich_onoff.a.iv[j] = 1;
            actdsurf->ale_dirich->dirich_type=dirich_FSI;
          }
          break;


        case ale:
          dsassert(hasdirich==0,
              "dirich- and fsi-coupling condition defined on same DSURF\n");
          dsassert(hascouple==0,
              "coupling- and fsi-coupling condition defined on same DSURF\n");
          dsassert(hasneum==0,
              "neumann- and fsi-coupling condition defined on same DSURF\n");

          /* allocate space for a dirichlet condition in this dsurf */
          actdsurf->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
          amdef("onoff",&(actdsurf->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
          amdef("val",&(actdsurf->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
          amdef("curve",&(actdsurf->dirich->curve),MAXDOFPERNODE,1,"IV");
          amdef("function",&(actdsurf->dirich->funct),MAXDOFPERNODE,1,"IV");
          amzero(&(actdsurf->dirich->dirich_onoff));
          amzero(&(actdsurf->dirich->dirich_val));
          amzero(&(actdsurf->dirich->funct));
          amzero(&(actdsurf->dirich->curve));

          /* initialise for fsi-coupling */
          for (j=0;j<genprob.ndim;j++)
            actdsurf->dirich->dirich_onoff.a.iv[j] = 1;
          actdsurf->dirich->dirich_type=dirich_FSI;
          break;


        default:
          dserror("fieldtyp unknown!\n");

      } /* switch(fieldtyp) */

    }  /* if (hasfsi!=0) */


    /* surface is not a fsi coupling surface */
    else
    {

#if 0
      /* surface is an exterior surface */
      if (actdsurf->ndvol == 1 )
      {

        /* create dirichlet condition for ale */
        actdsurf->ale_dirich =
          (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
        amdef("onoff",    &(actdsurf->ale_dirich->dirich_onoff),MAXDOFPERNODE, 1,"IV");
        amdef("val",      &(actdsurf->ale_dirich->dirich_val),  MAXDOFPERNODE, 1,"DV");
        amdef("curve",    &(actdsurf->ale_dirich->curve),       MAXDOFPERNODE, 1,"IV");
        amdef("function", &(actdsurf->ale_dirich->funct),       MAXDOFPERNODE, 1,"IV");
        amzero(&(actdsurf->ale_dirich->dirich_onoff));
        amzero(&(actdsurf->ale_dirich->dirich_val));
        amzero(&(actdsurf->ale_dirich->funct));
        amzero(&(actdsurf->ale_dirich->curve));

        /* initialise for ale boundary */
        for (j=0;j<genprob.ndim;j++)
          actdsurf->ale_dirich->dirich_onoff.a.iv[j] = 1;
        actdsurf->ale_dirich->dirich_type = dirich_none;
      }
#endif


    }  /* else */

  }  /* for (i=0; i<design->ndsurf; i++) */



  /* loop dlines */
  /* ----------- */
  for (i=0; i<design->ndline; i++)
  {
    hasdirich=0;
    hascouple=0;
    hasfsi   =0;
    hasneum  =0;
    actdline = &(design->dline[i]);

    /* check for conditions */
    if (actdline->dirich!=NULL) hasdirich++;
    if (actdline->couple!=NULL) hascouple++;
    if (actdline->fsicouple!=NULL) hasfsi++;
    if (actdline->neum!=NULL) hasneum++;


    /* surface is a fsi coupling surface */
    if (hasfsi!=0)
    {

      fieldtyp=actdline->fsicouple->fieldtyp;
      switch (fieldtyp)
      {

        case structure:
          if (hasdirich!=0) dswarning(1,7);
          dsassert(hasdirich==0,
              "dirich- and fsi-coupling condition defined on same DLINE\n");
          dsassert(hascouple==0,
              "coupling- and fsi-coupling condition defined on same DLINE\n");

          actdline->neum = (NEUM_CONDITION*)CCACALLOC(1,sizeof(NEUM_CONDITION));
          if (!actdline->neum) dserror("Allocation of memory failed");
          actdline->neum->neum_type = neum_FSI;
          break;


        case fluid:
          dsassert(hasdirich==0,
              "dirich- and fsi-coupling condition defined on same DLINE\n");
          dsassert(hascouple==0,
              "coupling- and fsi-coupling condition defined on same DLINE\n");
          dsassert(hasneum==0,
              "neumann- and fsi-coupling condition defined on same DLINE\n");

          /* allocate space for a dirichlet condition in this dline */
          actdline->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
          amdef("onoff",&(actdline->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
          amdef("val",&(actdline->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
          amdef("curve",&(actdline->dirich->curve),MAXDOFPERNODE,1,"IV");
          amdef("function",&(actdline->dirich->funct),MAXDOFPERNODE,1,"IV");
          amzero(&(actdline->dirich->dirich_onoff));
          amzero(&(actdline->dirich->dirich_val));
          amzero(&(actdline->dirich->curve));
          amzero(&(actdline->dirich->funct));

          /* initialise for fsi-coupling */
          for (j=0;j<genprob.ndim;j++)
            actdline->dirich->dirich_onoff.a.iv[j] = 1;
          actdline->dirich->dirich_type=dirich_FSI;


          /* also create the dirichlet condition for the ale */
          if (genprob.create_ale ==1)
          {
            /* allocate space for a dirichlet condition in this dsurf */
            actdline->ale_dirich =
              (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
            amdef("onoff",    &(actdline->ale_dirich->dirich_onoff), MAXDOFPERNODE,1,"IV");
            amdef("val",      &(actdline->ale_dirich->dirich_val),   MAXDOFPERNODE,1,"DV");
            amdef("curve",    &(actdline->ale_dirich->curve),        MAXDOFPERNODE,1,"IV");
            amdef("function", &(actdline->ale_dirich->funct),        MAXDOFPERNODE,1,"IV");
            amzero(&(actdline->ale_dirich->dirich_onoff));
            amzero(&(actdline->ale_dirich->dirich_val));
            amzero(&(actdline->ale_dirich->funct));
            amzero(&(actdline->ale_dirich->curve));

            /* initialise for fsi-coupling */
            for (j=0;j<genprob.ndim;j++)
              actdline->ale_dirich->dirich_onoff.a.iv[j] = 1;
            actdline->ale_dirich->dirich_type=dirich_FSI;
          }
          break;


        case ale:
          dsassert(hasdirich==0,
              "dirich- and fsi-coupling condition defined on same DLINE\n");
          dsassert(hascouple==0,
              "coupling- and fsi-coupling condition defined on same DLINE\n");
          dsassert(hasneum==0,
              "neumann- and fsi-coupling condition defined on same DLINE\n");

          /* allocate space for a dirichlet condition in this dline */
          actdline->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
          amdef("onoff",&(actdline->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
          amdef("val",&(actdline->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
          amdef("curve",&(actdline->dirich->curve),MAXDOFPERNODE,1,"IV");
          amdef("function",&(actdline->dirich->funct),MAXDOFPERNODE,1,"IV");
          amzero(&(actdline->dirich->dirich_onoff));
          amzero(&(actdline->dirich->dirich_val));
          amzero(&(actdline->dirich->curve));
          amzero(&(actdline->dirich->funct));

          /* initialise for fsi-coupling */
          for (j=0;j<genprob.ndim;j++)
            actdline->dirich->dirich_onoff.a.iv[j] = 1;
          actdline->dirich->dirich_type=dirich_FSI;
          break;


        default:
          dserror("fieldtyp unknown!\n");

      } /* switch(fieldtyp) */

    }  /* if (hasfsi!=0) */


    /* line is not a fsi coupling line */
    else
    {

      /* line is an exterior line */
      if (actdline->ndsurf == 1 )
      {

        /* create dirichlet condition for ale */
        actdline->ale_dirich =
          (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
        amdef("onoff",    &(actdline->ale_dirich->dirich_onoff),MAXDOFPERNODE, 1,"IV");
        amdef("val",      &(actdline->ale_dirich->dirich_val),  MAXDOFPERNODE, 1,"DV");
        amdef("curve",    &(actdline->ale_dirich->curve),       MAXDOFPERNODE, 1,"IV");
        amdef("function", &(actdline->ale_dirich->funct),       MAXDOFPERNODE, 1,"IV");
        amzero(&(actdline->ale_dirich->dirich_onoff));
        amzero(&(actdline->ale_dirich->dirich_val));
        amzero(&(actdline->ale_dirich->funct));
        amzero(&(actdline->ale_dirich->curve));

        /* initialise for ale boundary */
        for (j=0;j<genprob.ndim;j++)
          actdline->ale_dirich->dirich_onoff.a.iv[j] = 1;
        actdline->ale_dirich->dirich_type = dirich_none;
      }


    }  /* else */

  }  /* for (i=0; i<design->ndline; i++) */



  /* loop dnodes */
  /* ----------- */
  for (i=0; i<design->ndnode; i++)
  {
    hasdirich=0;
    hascouple=0;
    hasfsi   =0;
    hasneum  =0;
    actdnode = &(design->dnode[i]);

    /* check for conditions */
    if (actdnode->dirich!=NULL) hasdirich++;
    if (actdnode->couple!=NULL) hascouple++;
    if (actdnode->fsicouple!=NULL) hasfsi++;
    if (actdnode->neum!=NULL) hasneum++;
    if (hasfsi==0) continue;


    fieldtyp=actdnode->fsicouple->fieldtyp;
    switch (fieldtyp)
    {

      case structure:
        dsassert(hasdirich==0,
            "dirich- and fsi-coupling condition defined on same DNODE\n");
        dsassert(hascouple==0,
            "coupling- and fsi-coupling condition defined on same DNODE\n");

        actdnode->neum = (NEUM_CONDITION*)CCACALLOC(1,sizeof(NEUM_CONDITION));
        if (!actdnode->neum) dserror("Allocation of memory failed");
        actdnode->neum->neum_type = neum_FSI;
        break;


      case fluid:
        dsassert(hasdirich==0,
            "dirich- and fsi-coupling condition defined on same DNODE\n");
        dsassert(hascouple==0,
            "coupling- and fsi-coupling condition defined on same DNODE\n");
        dsassert(hasneum==0,
            "neumann- and fsi-coupling condition defined on same DNODE\n");

        /* allocate space for a dirichlet condition in this dnode */
        actdnode->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
        amdef("onoff",&(actdnode->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
        amdef("val",&(actdnode->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
        amdef("curve",&(actdnode->dirich->curve),MAXDOFPERNODE,1,"IV");
        amdef("function",&(actdnode->dirich->funct),MAXDOFPERNODE,1,"IV");
        amzero(&(actdnode->dirich->dirich_onoff));
        amzero(&(actdnode->dirich->dirich_val));
        amzero(&(actdnode->dirich->curve));
        amzero(&(actdnode->dirich->funct));

        /* initialise for fsi-coupling */
        for (j=0;j<genprob.ndim;j++)
          actdnode->dirich->dirich_onoff.a.iv[j] = 1;
        actdnode->dirich->dirich_type=dirich_FSI;


        /* also create the dirichlet condition for the ale */
        if (genprob.create_ale ==1)
        {
          /* allocate space for a dirichlet condition in this dsurf */
          actdnode->ale_dirich =
            (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
          amdef("onoff",   &(actdnode->ale_dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
          amdef("val",     &(actdnode->ale_dirich->dirich_val),  MAXDOFPERNODE,1,"DV");
          amdef("curve",   &(actdnode->ale_dirich->curve),       MAXDOFPERNODE,1,"IV");
          amdef("function",&(actdnode->ale_dirich->funct),       MAXDOFPERNODE,1,"IV");
          amzero(&(actdnode->ale_dirich->dirich_onoff));
          amzero(&(actdnode->ale_dirich->dirich_val));
          amzero(&(actdnode->ale_dirich->funct));
          amzero(&(actdnode->ale_dirich->curve));

          /* initialise for fsi-coupling */
          for (j=0;j<genprob.ndim;j++)
            actdnode->ale_dirich->dirich_onoff.a.iv[j] = 1;
          actdnode->ale_dirich->dirich_type=dirich_FSI;
        }
        break;


      case ale:
        dsassert(hasdirich==0,
            "dirich- and fsi-coupling condition defined on same DNODE\n");
        dsassert(hascouple==0,
            "coupling- and fsi-coupling condition defined on same DNODE\n");
        dsassert(hasneum==0,
            "neumann- and fsi-coupling condition defined on same DNODE\n");

        /* allocate space for a dirichlet condition in this dnode */
        actdnode->dirich = (DIRICH_CONDITION*)CCACALLOC(1,sizeof(DIRICH_CONDITION));
        amdef("onoff",&(actdnode->dirich->dirich_onoff),MAXDOFPERNODE,1,"IV");
        amdef("val",&(actdnode->dirich->dirich_val),MAXDOFPERNODE,1,"DV");
        amdef("curve",&(actdnode->dirich->curve),MAXDOFPERNODE,1,"IV");
        amdef("function",&(actdnode->dirich->funct),MAXDOFPERNODE,1,"IV");
        amzero(&(actdnode->dirich->dirich_onoff));
        amzero(&(actdnode->dirich->dirich_val));
        amzero(&(actdnode->dirich->curve));
        amzero(&(actdnode->dirich->funct));

        /* initialise for fsi-coupling */
        for (j=0;j<genprob.ndim;j++)
          actdnode->dirich->dirich_onoff.a.iv[j] = 1;
        actdnode->dirich->dirich_type=dirich_FSI;
        break;


      default:
        dserror("fieldtyp unknown!\n");


    } /* switch(fieldtyp) */

  }  /* for (i=0; i<design->ndnode; i++) */


#ifdef DEBUG
  dstrc_exit();
#endif


  return;
} /* end of fsi_creatcoup*/





/*-----------------------------------------------------------------------*/
/*!
  \brief initialise fsi coupling conditions

  create mulitfield solution history and set pointers to the corresponding
  nodes of the other fields

  \param *structfield   FIELD         (i)      structure field
  \param *fluidfield    FIELD         (i)      fluid field
  \param *alefield      FIELD         (i)      ale field

  \return void

  \author genk
  \date   09/02

 */
/*-----------------------------------------------------------------------*/
void fsi_initcoupling(
    FIELD              *structfield,
    INT                 disnum_s,
    FIELD              *fluidfield,
    INT                 disnum_f,
    FIELD              *alefield,
    INT                 disnum_a
    )

{

  INT     numfnp, numsnp, numanp;              /* number of nodes         */
  INT     numdf;                               /* number of dofs          */
  INT     numc=0;                              /* number of columns in mf */
  INT     numaf,numsf,numff;
  INT     i,j;                                 /* simply some counters    */
  INT     ierr;                                /* flag                    */
  INT     dim;                                 /* dimension of problem    */
  INT     sfound,afound;                       /* flag                    */
  INT     is_ale;

  DOUBLE  tol=EPS4;                            /* tolerance for node dist */

  NODE   *actfnode=NULL, *actsnode=NULL, *actanode=NULL;  /* actual nodes */
  GNODE  *actfgnode=NULL,*actsgnode=NULL,*actagnode=NULL; /* actual gnodes*/

  ELEMENT *actele;

  ARRAY   aindex_a, sindex_a;
  INT    *aindex, *sindex;

  FSI_DYNAMIC    *fsidyn;


#ifdef DEBUG
  dstrc_enter("fsi_initcoupling");
#endif


  fsidyn= alldyn[3].fsidyn;


  /* find number of nodes in different fields */
  numsnp  = structfield->dis[disnum_s].numnp;
  numfnp  = fluidfield->dis[disnum_f].numnp;
  numanp  = alefield->dis[disnum_a].numnp;
  numaf   = genprob.numaf;
  numff   = genprob.numff;
  numsf   = genprob.numsf;
  dim     = genprob.ndim;

#ifdef PERF
perf_begin(53);
#endif



  /* allocate space for mulitfield solution history:
   * -----------------------------------------------
   *
   * the following data are necassary:
   * ale needs displacements of the structural fsi coupling nodes
   * fluid needs velocity of all ale nodes
   * structure needs stresses of the fluid fsi coupling nodes
   */


  /* multifield solution history of fluid nodes:
   * -------------------------------------------
   *   actfnode->sol_mf.a.da[0][i]: velocity transfered to ale
   *   actfnode->sol_mf.a.da[1][i]: stresses transfered to structure
   *     2D: SIGMA11,SIGMA22,SIGMA12
   *     3D: SIGMA11,SIGMA22,SIGMA33,SIGMA12,SIGMA13,SIGMA23
   */


  /* loop fluid nodes */
  for (i=0;i<numfnp;i++)
  {
    actfnode  = &(fluidfield->dis[disnum_f].node[i]);
    actfgnode = actfnode->gnode;
    numdf = actfnode->numdf;

    if (numdf==3 || numdf==5) numc=numdf;
    else if (numdf==4 || numdf==7) numc=IMAX(6,numdf);
    else  dserror("number of fluid dofs not possible!\n");

    amdef("sol_mf",&actfnode->sol_mf,3,numc,"DA");
    amzero(&(actfnode->sol_mf));

    actfgnode->mfcpnode=(NODE**)CCACALLOC(3,sizeof(NODE*));

    if (actfgnode->mfcpnode==NULL)
      dserror("Allocation of coupling node pointers failed");

    for (j=0;j<3;j++) actfgnode->mfcpnode[j]=NULL;

  } /* end of loop over fluid nodes */



  /* multifield solution history of structural nodes:
   * ------------------------------------------------
   * actsnode->sol_mf.a.da[0][i]: actual displacements transfered to ale
   * actsnode->sol_mf.a.da[1][i]: displacements of old iteration step
   * actsnode->sol_mf.a.da[2][i]: displacements of old time step
   * actsnode->sol_mf.a.da[3][i]: dispi
   * actsnode->sol_mf.a.da[4][i]: couplingforces at the end of time step
   * actsnode->sol_mf.a.da[5][i]: coupling forces at the beginning of
   * time step
   * actsnode->sol_mf.a.da[6][i]: used in Gradientenverf
   * actsnode->sol_mf.a.da[7][i]: used in Newton
   * actsnode->sol_mf.a.da[8][i]: used in Newton
   * [9]
   */

  /* loop structure nodes */
  for (i=0;i<numsnp;i++)
  {
    actsnode  = &(structfield->dis[disnum_s].node[i]);
    actsgnode = actsnode->gnode;
    numdf = actsnode->numdf;

    amdef("sol_mf",&actsnode->sol_mf,10,numdf,"DA");
    amzero(&(actsnode->sol_mf));

    actsgnode->mfcpnode=(NODE**)CCACALLOC(3,sizeof(NODE*));

    if (actsgnode->mfcpnode==NULL)
      dserror("Allocation of coupling node pointers failed");

    for (j=0;j<3;j++) actsgnode->mfcpnode[j]=NULL;

  } /* end of loop over struct nodes */



  /* multifield solution history of ale nodes:
   * -----------------------------------------
   * actanode->sol_mf.a.da[0][i]: displacements at (n)
   * actanode->sol_mf.a.da[1][i]: displacements at (n+1)
   */

  /* loop ale nodes */
  for (i=0;i<numanp;i++)
  {
    actanode  = &(alefield->dis[disnum_a].node[i]);
    actagnode = actanode->gnode;
    numdf = actanode->numdf;

    amdef("sol_mf",&actanode->sol_mf,4,numdf,"DA");
    amzero(&(actanode->sol_mf));

    actagnode->mfcpnode=(NODE**)CCACALLOC(3,sizeof(NODE*));

    if (actagnode->mfcpnode==NULL)
      dserror("Allocation of coupling node pointers failed");

    for (j=0;j<3;j++) actagnode->mfcpnode[j]=NULL;

  } /* end of loop over ale nodes */


#ifdef PERF
perf_end(53);
#endif





  /* find fsi coupled nodes
   * and set ptrs to the corresponding nodes of the other fields */

#ifdef PERF
perf_begin(54);
#endif



  /* create and initialsise index arrays */
  aindex = amdef("aindex",&aindex_a,numanp,1,"IV");
  sindex = amdef("sindex",&sindex_a,numsnp,1,"IV");
  for (i=0;i<numanp;i++) aindex[i]=1;
  for (i=0;i<numsnp;i++) sindex[i]=1;



  /* loop fluid nodes  */
  for (i=0;i<numfnp;i++)
  {
    actfnode  = &(fluidfield->dis[disnum_f].node[i]);
    actfgnode = actfnode->gnode;
    is_ale=0;

    for (j=0;j<actfnode->numele;j++)
    {

      actele= actfnode->element[j];

      switch (actele->eltyp)
      {

#ifdef D_FLUID2
        case el_fluid2:
          if (actele->e.f2->is_ale>0) is_ale++;
          break;
#endif


#ifdef D_FLUID3
        case el_fluid3:
          if (actele->e.f3->is_ale>0) is_ale++;
          break;
#endif


#ifdef D_FLUID3_F
        case el_fluid3_fast:
          if (actele->e.f3->is_ale>0) is_ale++;
          break;
#endif

        default:
          dserror("eltyp unknow\n");

      }  /* switch (actele->eltyp) */

      if (is_ale>0) break;

    }  /* for (j=0;j<actfnode->numele;j++) */


    sfound=0;
    afound=0;


#ifdef PERF
perf_begin(56);
#endif


    /* loop struct nodes at interface and find corresponding node */
    if (actfgnode->fsicouple != NULL)
      for (j=0;j<numsnp;j++)
      {

        if (sindex[j]==0) continue;

        actsnode   = &(structfield->dis[disnum_s].node[j]);

        cheque_distance(&(actfnode->x[0]),&(actsnode->x[0]),tol,&ierr);
        if(ierr==0) continue;

        sfound++;
        actsgnode = actsnode->gnode;
        sindex[j]=0;

        break;
      } /* end of loop over structnodes */


#ifdef PERF
    perf_end(56);
#endif





#ifdef PERF
    perf_begin(57);
#endif

    if (genprob.create_ale != 1 || fluidfield->subdivide > 0)
    {
      /* loop ale nodes and find corresponding node */
      if (is_ale>0)
        for (j=0;j<numanp;j++)
        {

          if (aindex[j]==0) continue;

          actanode  = &(alefield->dis[disnum_a].node[j]);

          cheque_distance(&(actfnode->x[0]),&(actanode->x[0]),tol,&ierr);
          if(ierr==0) continue;

          afound++;
          actagnode = actanode->gnode;
          aindex[j]=0;

          break;

        } /* end of loop over ale nodes */
    }

    else
    {
      if (is_ale>0)
      {
        for (j=0; j<actfnode->element[0]->numnp; j++)
        {
          if (actfnode->Id == actfnode->element[0]->node[j]->Id )
          {
            switch (actfnode->element[0]->eltyp)
            {
#ifdef D_FLUID3
              case el_fluid3:
              case el_fluid3_fast:
                actanode = actfnode->element[0]->e.f3->ale_ele->node[j];
                break;
#endif

#ifdef D_FLUID2
              case el_fluid2:
                actanode = actfnode->element[0]->e.f2->ale_ele->node[j];
                break;
#endif

              default:
                dserror("");
                break;
            }  /* switch (actfnode->element[0]->eltyp) */

            break;  /* break the loop over the nodes of the element */

          }  /* if (actfnode->Id == actfnode->element[0]->node[j]->Id ) */
        }


#ifdef DEBUG
        cheque_distance(&(actfnode->x[0]),&(actanode->x[0]),tol,&ierr);
        if(ierr==0) dserror("Wrong node matching for fluid <-> ale!!!");
#endif

        afound++;
        actagnode = actanode->gnode;
      }
    }

#ifdef PERF
perf_end(57);
#endif




    /* set pointers to corresponding nodes */
    if(sfound>0)   actfgnode->mfcpnode[numsf]=actsnode;
                   actfgnode->mfcpnode[numff]=actfnode;
    if(afound>0)   actfgnode->mfcpnode[numaf]=actanode;
    if(sfound>0)   actsgnode->mfcpnode[numsf]=actsnode;
    if(sfound>0)   actsgnode->mfcpnode[numff]=actfnode;
    if(sfound>0)   actsgnode->mfcpnode[numaf]=actanode;
    if(sfound>0)   actagnode->mfcpnode[numsf]=actsnode;
    if(afound>0)   actagnode->mfcpnode[numff]=actfnode;
    if(afound>0)   actagnode->mfcpnode[numaf]=actanode;

  }  /* for (i=0;i<numfnp;i++) */


#ifdef PERF
perf_end(54);
#endif


  amdel(&aindex_a);
  amdel(&sindex_a);


  if (genprob.visual>0) goto end;
  if (fsidyn->coupmethod == 0) goto end;

#ifdef PERF
perf_begin(55);
#endif



  /* plausibility checks */
  /* ------------------- */
  for (i=0;i<numfnp;i++)
  {
    actfnode  = &(fluidfield->dis[disnum_f].node[i]);
    actfgnode = actfnode->gnode;

    if (actfgnode->fsicouple==NULL) continue;

    actsnode = actfgnode->mfcpnode[numsf];
    actanode = actfgnode->mfcpnode[numaf];
    actsgnode = actsnode->gnode;
    actagnode = actanode->gnode;


    if (actsgnode == NULL)
      dserror("No structure node for fluid node: #%d", actfnode->Id);


    if (actagnode == NULL)
      dserror("No ale node for fluid node: #%d", actfnode->Id);


    /* check locsys */
    if (actfnode->locsysId!=0)
      dserror("No locsys at FSI coupling node: #%d", actfnode->Id);
    if (actsnode->locsysId!=0)
      dserror("No locsys at FSI coupling node: #%d", actsnode->Id);
    if (actanode->locsysId!=0)
      dserror("No locsys at FSI coupling node: #%d", actanode->Id);


    /* check coupling conditions */
    if (actsgnode->fsicouple==NULL)
    {
      dserror("FSI Coupling Condition Fluid-Struct not the same: struct node #%d",
          actsnode->Id);
    }

    if (actagnode->fsicouple==NULL)
      dserror("FSI Coupling Condition Fluid-Ale not the same: ale node #%d",
          actanode->Id);


    /* check coupling Ids */
    if (actfgnode->fsicouple->fsi_coupleId!=actsgnode->fsicouple->fsi_coupleId)
      dserror("FSI Coupling Condition Fluid-Struct: wrong coupleId, fluid #%d, struct #%d",
          actfnode->Id, actsnode->Id);

    if (actfgnode->fsicouple->fsi_coupleId!=actagnode->fsicouple->fsi_coupleId)
      dserror("FSI Coupling Condition Fluid-Ale: wrong coupleId, fluid #%d, ale #%d",
          actfnode->Id, actsnode->Id);


    /* check mesh */
    if (actfgnode->fsicouple->fsi_mesh!=actsgnode->fsicouple->fsi_mesh)
      dserror("FSI Coupling Condition Fluid-Struct: wrong mesh, fluid #%d, struct #%d",
          actfnode->Id, actsnode->Id);

    if (actfgnode->fsicouple->fsi_mesh!=actagnode->fsicouple->fsi_mesh)
      dserror("FSI Coupling Condition Fluid-Ale: wrong mesh, fluid #%d, ale #%d",
          actfnode->Id, actanode->Id);


    /* check dirich conds of fluid node */
    if (actfgnode->dirich==NULL)
      dserror("No dirich condition for fsi-coupling fluid node #%d",actfnode->Id);
    if(actfgnode->dirich->dirich_type!=dirich_FSI) continue;
    if (actfnode->numdf!=dim+1)
      dserror("numdf not possible for fluid node #%d at FSI interface!", actfnode->Id);
    for (j=0;j<dim;j++)
      if (actfgnode->dirich->dirich_onoff.a.iv[j]!=1)
        dserror("onoff(%d)=%d at fluid node #%d",
            j,actfgnode->dirich->dirich_onoff.a.iv[j],actfnode->Id);
    if (actfgnode->dirich->dirich_onoff.a.iv[dim]!=0)
      dserror("onoff(%d)=%d at fluid node #%d",
          dim,actfgnode->dirich->dirich_onoff.a.iv[j],actfnode->Id);

  }  /* for (i=0;i<numfnp;i++) */


  for (i=0;i<numanp;i++)
  {
    actanode  = &(alefield->dis[disnum_a].node[i]);
    actagnode = actanode->gnode;

    if (actagnode->fsicouple==NULL) continue;

    if (actagnode->dirich==NULL)
      dserror("No dirich condition for fsi-coupling ale node #%d",actanode->Id);

    if(actagnode->dirich->dirich_type!=dirich_FSI) continue;

    for(j=0;j<actanode->numdf;j++)
      if(actagnode->dirich->dirich_onoff.a.iv[j]!=1)
        dserror("wrong onoff() at ale node #%d",actanode->Id);

    for (j=actanode->numdf;j<MAXDOFPERNODE;j++)
      actagnode->dirich->dirich_onoff.a.iv[j]=0;

  }  /* for (i=0;i<numanp;i++) */



#ifdef PERF
perf_end(55);
#endif

  /* print out coupling */
  if (genprob.visual==0)
    out_fsi(fluidfield);

end:


#ifdef DEBUG
  dstrc_exit();
#endif


  return;
} /* end of fsi_initcoup*/







/*-----------------------------------------------------------------------*/
/*!
  \brief determine structural fsi interface dofs

  create array sid (structural interface dofs)

  \param *structfield   FIELD         (i)      structure field
  \param *fsidyn        FSI_DYNAMIC   (i)

  \return void

  \author genk
  \date   01/03

 */
/*-----------------------------------------------------------------------*/
void fsi_struct_intdofs(
    FIELD              *structfield,
    INT                 disnum
    )

{

  INT    i,j;                    /* some counters                         */
  INT    numnp_total;            /* total number of structure dofs        */
  INT    dof;                    /* actual dof                            */
  INT    counter=0;
  INT    numaf;
  INT   *sid;                    /* structural interface dofs             */

  NODE  *actsnode, *actanode;    /* actual nodes                          */
  GNODE *actsgnode, *actagnode;  /* actual gnodes                         */

  FSI_DYNAMIC *fsidyn;


#ifdef DEBUG
  dstrc_enter("fsi_struct_intdofs");
#endif


  fsidyn = alldyn[3].fsidyn;

  numnp_total = structfield->dis[disnum].numnp;
  numaf       = genprob.numaf;

  sid = amdef("sid",&fsidyn->sid,structfield->dis[disnum].numdf,1,"IV");
  amzero(&fsidyn->sid);


  /* loop structure nodes */
  for (i=0;i<numnp_total;i++)
  {
    actsnode  = &(structfield->dis[disnum].node[i]);
    actsgnode = actsnode->gnode;

    actanode  = actsgnode->mfcpnode[numaf];
    if (actanode == NULL) continue;
    actagnode = actanode->gnode;

    /* check for coupling nodes */
    if(actagnode->dirich == NULL)
      dserror("no dirich condition for coupled ALE node #%d",actanode->Id);

    if (actagnode->dirich->dirich_type != dirich_FSI) continue;

    for (j=0;j<actanode->numdf;j++)
    {
      dof = actsnode->dof[j];
      sid[dof]=1;
      counter++;
    }

  } /* end of loop over nodes */

  fsidyn->numsid=counter;


#ifdef DEBUG
  dstrc_exit();
#endif


  return;

} /* end of fsi_struct_intdofs*/



#endif
#endif  /* ifdef D_FSI */

/*! @} (documentation module close)*/

