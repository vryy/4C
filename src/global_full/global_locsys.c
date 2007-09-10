/*!----------------------------------------------------------------------
\file global_locsys.c
\brief all routines for local co-systems

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127

General explanation:
   Local co-systems may be used for defining tangential or normal
   DBCs along inclined or even curved Lines/Surfaces, e.g.
   slip dirichlet condition for fluids along a curved wall.
   In such a case the DBCs  cannot be prescribed in the XYZ-system
   any more.

   One defines a new co-system xyz* for the specific design object.
   The DBCs in the inpute file are then given in this new co-system.
   The locsys information is inherited from the Dobjects to the nodes
   and all elements belonging to one of this node get this
   information, too. (locsys_inherit_to_node).

Bringing in the locsys information into the solution:
   All nodal sol_arrays are still in the XYZ co-system. The element
   matrix and the ele-RHS are also build in the XYZ-system.
   This means before applying the dirichlet conditions one has to
   tranform the used sol-array to the xyz* system, then the prescribed
   values are applied and afterwards one has to transform back the
   sol array. (locsys_trans_sol_dirich)

   After building the element stiffness matrix and the ele-RHS, one has
   to transform them to the xyz* system (before assembling!)
   (locsys_trans).

   When calculating the RHS due to dirichlet conditions (condensation)
   the prescribed values have to be in the xyz* co-system.

   In the resulting solution the values can be in different co-sys.
   After copying back them to the sol-arrays at the nodes the solution
   at the nodes with a locsys have to be transformed to the XYZ system
   (locsys_trans_sol)

Final Remarks:
   To summerise: in the nodal sol-arrays the values are stored in the
   XYZ co-system. Temporarily one transforms them (all or just the
   dirichlet values) to xyz*. After the calculation one immediatetely
   transforms them back.

   The whole procedure looks quite complicated. However I did not
   see another possiblity for the used fluid-solver. It might
   be simpler, when only working with struct elements ...

</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2.h"
#include "../fluid3/fluid3.h"
#include "../ale2/ale2.h"
#include "../ale3/ale3.h"
#ifdef D_SOLID3
#include "../solid3/solid3.h"
#endif
/*!----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | INT                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
extern INT            numcurve;
extern struct _CURVE *curve;
/*----------------------------------------------------------------------*
 |                                                        genk 04/04    |
 | number of local co-ordinate systems                                  |
 | vector of structures of local co-ordinate systems                    |
 | defined in input_locsys.c                                          |
 *----------------------------------------------------------------------*/
extern INT            numlocsys;
extern struct _LOCSYS *locsys;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

/*-------------------------------------- static variables for this file */
static ARRAY     trans_a;
static DOUBLE  **trans;          /* transformation matrix */
static ARRAY     workm_a;
static DOUBLE  **workm;          /* working matrix for temp. result */
static ARRAY     workv_a;
static DOUBLE   *workv;          /* working vector for temp. result */
static ARRAY     nodalwork_a;
static DOUBLE   *nodalwork;      /* nodal working vector */

/*!---------------------------------------------------------------------
\brief inherit local co-ordinate system to elements

<pre>                                                         genk 04/04

local co-ordinate systems are assigned to design elements (points, lines
surfs, vols). They are NOT inherited to their lower design elements.

create local co-ordinate system on node level

</pre>
\return void

------------------------------------------------------------------------*/
void locsys_inherit_to_node()
{
INT      i,j,k,l;
INT      locsysId=-1;
NODE     *actnode;
GNODE    *actgnode;
ELEMENT  *actele;

#ifdef DEBUG
dstrc_enter("locsys_inherit_to_ele");
#endif


/*-------------------------- allocate transformation and working arrays
  since there may be internally defined locsys, these arrays have to be
  allocated for ALL cases                                               */
trans    =amdef("trans"    ,&trans_a    ,MAXDOFPERELE,MAXDOFPERELE,"DA");
workm    =amdef("workm"    ,&workm_a    ,MAXDOFPERELE,MAXDOFPERELE,"DA");
workv    =amdef("workv"    ,&workv_a    ,MAXDOFPERELE,1           ,"DV");
nodalwork=amdef("nodalwork",&nodalwork_a,MAXDOFPERELE,1           ,"DV");

if (numlocsys==0) goto end;

/*--------------------------------------------------------- loop fields */
for (i=0; i<genprob.numfld; i++)
{
   /*--------------------------------------------- loop discretisations */
   for (j=0; j<field[i].ndis; j++)
   {
      /*---------------------------------------------------- loop nodes */
      for (k=0; k<field[i].dis[j].numnp; k++)
      {
         actnode = &(field[i].dis[j].node[k]);
         actgnode = actnode->gnode;
         locsysId = 0;  /* initialise */
         /*------ local co-cordinate system defined by design condition */
         switch(actgnode->ondesigntyp)
         {
         case ondnode:
           /* inherit of design point */
           if (actgnode->d.dnode->locsysId > 0)
           {
             locsysId = actgnode->d.dnode->locsysId;
           }
#if 0
/* #ifdef LOCALSYSTEMS_ST */
           /* On the one hand, thee following deactivated inheritances are
            * a nice idea to avoid assignments of local systems to design
            * entities. On the other hand, the inheritance is not transparent
            * to the user such that a inherently applied local system affects
            */
           /* inherit of design lines */
           if (locsysId == 0)
           {
             INT idl = 0;  /* design line index */
             /* loop all adjacent design lines and give its locsysId to
              * the node */
             while ( (idl<actgnode->d.dnode->ndline) && (locsysId == 0) )
             {
               if (actgnode->d.dnode->dline[idl]->locsysId > 0)
               {
                 locsysId = actgnode->d.dnode->dline[idl]->locsysId;
               }
               idl++;
             }
           }
           /* inherit of design surfaces */
           if (locsysId == 0)
           {
             INT idl = 0; /* design line index */
             /* loop all adjacent design lines to loop their design surfaces */
             while ( (idl<actgnode->d.dnode->ndline) && (locsysId == 0) )
             {
               DLINE* dline = actgnode->d.dnode->dline[idl];  /* curr. design line */
               INT ids = 0;  /* design surface index */
               /* loop all adjacent surfaces */
               while ( (ids<dline->ndsurf) && (locsysId == 0) )
               {
                 if (dline->dsurf[ids]->locsysId > 0)
                 {
                   locsysId = dline->dsurf[ids]->locsysId;
                 }
                 ids++;
               }
               idl++;
             }
           }
           /* inherit of design volumes */
           if (locsysId == 0)
           {
             INT idl = 0; /* design line index */
             /* loop all adjacent design lines to loop their design surfaces */
             while ( (idl<actgnode->d.dnode->ndline) && (locsysId == 0) )
             {
               DLINE* dline = actgnode->d.dnode->dline[idl];   /* curr. design line */
               INT ids = 0;  /* design surface index */
               /* loop all adjacent surfaces */
               while ( (ids<dline->ndsurf) && (locsysId == 0) )
               {
                 DSURF* dsurf = dline->dsurf[ids];  /* curr. design surface */
                 INT idv = 0;  /* design volume index */
                 /* loop all adjacent volumes */
                 while ( (idv<dsurf->ndvol) && (locsysId == 0) )
                 {
                   if (dsurf->dvol[idv]->locsysId > 0)
                   {
                     locsysId = dsurf->dvol[idv]->locsysId;
                   }
                   idv++;
                 }
                 ids++;
               }
               idl++;
             }
           }
#endif
           break;
         case ondline:
           /* inherit of design lines */
           if (actgnode->d.dline->locsysId > 0)
           {
             locsysId = actgnode->d.dline->locsysId;
           }
#if 0
/* #ifdef LOCALSYSTEMS_ST */
           /* inherit of design surfaces */
           if (locsysId == 0)
           {
             DLINE* dline = actgnode->d.dline;  /* curr. design line */
             INT ids = 0;  /* design surface index */
             /* loop all adjacent surfaces */
             while ( (ids<dline->ndsurf) && (locsysId == 0) )
             {
               if (dline->dsurf[ids]->locsysId > 0)
               {
                 locsysId = dline->dsurf[ids]->locsysId;
               }
               ids++;
             }
           }
           /* inherit of design volumes */
           if (locsysId == 0)
           {
             DLINE* dline = actgnode->d.dline;   /* curr. design line */
             INT ids = 0;  /* design surface index */
             /* loop all adjacent surfaces */
             while ( (ids<dline->ndsurf) && (locsysId == 0) )
             {
               DSURF* dsurf = dline->dsurf[ids];  /* curr. design surface */
               INT idv = 0;  /* design volume index */
               /* loop all adjacent volumes */
               while ( (idv<dsurf->ndvol) && (locsysId == 0) )
               {
                 if (dsurf->dvol[idv]->locsysId > 0)
                 {
                   locsysId = dsurf->dvol[idv]->locsysId;
                 }
                 idv++;
               }
               ids++;
             }
           }
#endif
           break;
         case ondsurf:
           /* inherit of design surface */
           if (actgnode->d.dsurf->locsysId > 0)
           {
             locsysId = actgnode->d.dsurf->locsysId;
           }
#if 0
/* #ifdef LOCALSYSTEMS_ST */
           /* inherit of design volumes */
           if (locsysId == 0)
           {
             DSURF* dsurf = actgnode->d.dsurf;  /* curr. design surface */
             INT idv = 0;  /* design volume index */
             /* loop all adjacent volumes */
             while ( (idv<dsurf->ndvol) && (locsysId == 0) )
             {
               if (dsurf->dvol[idv]->locsysId > 0)
               {
                 locsysId = dsurf->dvol[idv]->locsysId;
               }
               idv++;
             }
           }
#endif
           break;
         case ondvol:
           /* inherit of design volume */
           if (actgnode->d.dvol->locsysId > 0)
           {
             locsysId = actgnode->d.dvol->locsysId;
           }
           break;
         case ondnothing:
           dserror("GNODE not owned by any design object");
           break;
         default:
           dserror("Cannot create locsys on element level");
           break;
         }
         /* set NODE local system ID */
         actnode->locsysId=locsysId;
         /* forward "NODE has local system ID" to its associated elements */
         if (locsysId>0)
         {
            for (l=0;l<actnode->numele;l++)
            {
               actele=actnode->element[l];
               actele->locsys=locsys_yes;
            } /* end loop over nodes */
         } /* endif (locsysId>0) */
      } /* end loop over nodes */
   } /* end loop over discretisations */
} /* end loop over fields */

/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of locsys_create */

/*!---------------------------------------------------------------------
\brief transform element stiffness matrix

<pre>                                                         genk 04/04

the element stiffness matrix of the actual element is transformed from
the global XYZ - co-ordinate system to the local one.
This names are a bit confusing, since the "local" co-ordinate system is
the one we are solving in, so we better introduce new names:

estif  = stiffness matrix in the given XYZ cartesian co-system
estif* = stiffness matrix in the alternative xyz* co-system defined
         in the input file

eload  = load vector in the given XYZ cartesian co-system
eload* = load vector in the alternative xyz* co-system

trans  = transformation matrix containing the direction cosini between
         XYZ and xyz*

</pre>
\param   *ele        ELEMENT     (i)      actual element
\param  **estif1     DOUBLE      (i/o)    element matrix
\param  **estif2     DOUBLE      (i/o)    element matrix
\param   *vec1       DOUBLE      (i/o)    element RHS
\param   *vec2       DOUBLE      (i/o)    element RHS
\sa locsys_trans_nodval

\return void

------------------------------------------------------------------------*/
void locsys_trans(ELEMENT *ele, DOUBLE **estif1, DOUBLE **estif2,
                                DOUBLE *vec1,    DOUBLE *vec2)
{
INT      i,j;        /* simply some counters */
INT      nd;         /* counter for total number of element dofs */
INT      ilocsys;    /* index of locsys */
INT      numdf;      /* number of dofs at actual node */
NODE     *actnode;   /* actual node */

#ifdef DEBUG
dstrc_enter("locsys_trans");
#endif
#ifdef PERF
  perf_begin(22);
#endif

/*------------------------------------------------- initialise matrices */
amzero(&trans_a);

/*------------------------------------------ fill transformation matrix */
nd=0;
switch (ele->eltyp)
{
#ifdef D_FLUID2
case el_fluid2:
   dsassert (ele->e.f2->fs_on<=2,
             "no local co-ordinate system on free surface allowed!\n");
   /*----------------------------------------------- loop element nodes */
   for (i=0;i<ele->numnp;i++)
   {
      actnode=ele->node[i];
      numdf=actnode->numdf;
      ilocsys=actnode->locsysId-1;
      if(ilocsys>=0)
      {
         dsassert(ilocsys<numlocsys,"locsysId not existent!\n");
         if (numdf<4)
         {
            trans[nd][nd]     = locsys[ilocsys].lXx;
            trans[nd+1][nd]   = locsys[ilocsys].lXy;
            trans[nd][nd+1]   = locsys[ilocsys].lYx;
            trans[nd+1][nd+1] = locsys[ilocsys].lYy;
            trans[nd+2][nd+2] = ONE;
         }
         else if (numdf==4)
         {
            dserror("transformation for fluid node with 4 dofs not implemented!\n");
         }
         else if (numdf==5) /*  node at free surf. w/ five dofs
                            [vel, vel, pre, velg, velg]                 */
         {
            trans[nd][nd]     = locsys[ilocsys].lXx;
            trans[nd+1][nd]   = locsys[ilocsys].lXy;
            trans[nd][nd+1]   = locsys[ilocsys].lYx;
            trans[nd+1][nd+1] = locsys[ilocsys].lYy;
            trans[nd+2][nd+2] = ONE;
            trans[nd+3][nd+3] = locsys[ilocsys].lXx;
            trans[nd+4][nd+3] = locsys[ilocsys].lXy;
            trans[nd+3][nd+4] = locsys[ilocsys].lYx;
            trans[nd+4][nd+4] = locsys[ilocsys].lYy;
         }
         else
            dserror("transformation not possible!\n");
      }
      else
      {
         for (j=0;j<numdf;j++) trans[nd+j][nd+j] = ONE;
      }
      nd+=numdf;
   } /* end loop over element nodes */
break;
#endif
#ifdef D_FLUID3
case el_fluid3:
   dsassert (ele->e.f3->fs_on<=2,
             "no local co-ordinate system on free surface allowed!\n");
   /*----------------------------------------------- loop element nodes */
   for (i=0;i<ele->numnp;i++)
   {
      actnode=ele->node[i];
      numdf=actnode->numdf;
      ilocsys=actnode->locsysId-1;
      if(ilocsys>=0)
      {
         dsassert(ilocsys<numlocsys,"locsysId not existent!\n");
         if (numdf<5)
         {
            trans[nd][nd]     = locsys[ilocsys].lXx;
            trans[nd+1][nd]   = locsys[ilocsys].lXy;
            trans[nd+2][nd]   = locsys[ilocsys].lXz;
            trans[nd][nd+1]   = locsys[ilocsys].lYx;
            trans[nd+1][nd+1] = locsys[ilocsys].lYy;
            trans[nd+2][nd+1] = locsys[ilocsys].lYz;
            trans[nd][nd+2]   = locsys[ilocsys].lZx;
            trans[nd+1][nd+2] = locsys[ilocsys].lZy;
            trans[nd+2][nd+2] = locsys[ilocsys].lZz;
            trans[nd+3][nd+3] = ONE;
         }
         else if (numdf==5)
         {
            dserror("transformation for fluid node with 5 dofs not implemented!\n");
         }
         else if (numdf==7) /*  node at free surf. w/ five dofs
                          [vel, vel, vel, pre, velg, velg, velg]        */
         {
            trans[nd][nd]     = locsys[ilocsys].lXx;
            trans[nd+1][nd]   = locsys[ilocsys].lXy;
            trans[nd+2][nd]   = locsys[ilocsys].lXz;
            trans[nd][nd+1]   = locsys[ilocsys].lYx;
            trans[nd+1][nd+1] = locsys[ilocsys].lYy;
            trans[nd+2][nd+1] = locsys[ilocsys].lYz;
            trans[nd][nd+2]   = locsys[ilocsys].lZx;
            trans[nd+1][nd+2] = locsys[ilocsys].lZy;
            trans[nd+2][nd+2] = locsys[ilocsys].lZz;
            trans[nd+3][nd+3] = ONE;
            trans[nd+4][nd+4] = locsys[ilocsys].lXx;
            trans[nd+5][nd+4] = locsys[ilocsys].lXy;
            trans[nd+6][nd+4] = locsys[ilocsys].lXz;
            trans[nd+4][nd+5] = locsys[ilocsys].lYx;
            trans[nd+5][nd+5] = locsys[ilocsys].lYy;
            trans[nd+6][nd+5] = locsys[ilocsys].lYz;
            trans[nd+4][nd+6] = locsys[ilocsys].lZx;
            trans[nd+5][nd+6] = locsys[ilocsys].lZy;
            trans[nd+6][nd+6] = locsys[ilocsys].lZz;
         }
         else
            dserror("transformation not possible!\n");
      }
      else
      {
         for (j=0;j<numdf;j++) trans[nd+j][nd+j] = ONE;
      }
      nd+=numdf;
   } /* end loop over element nodes */
break;
#endif
#ifdef D_ALE
case el_ale2:
   /*----------------------------------------------- loop element nodes */
   for (i=0;i<ele->numnp;i++)
   {
      actnode=ele->node[i];
      numdf=actnode->numdf;
      dsassert(numdf==2,
            "numdf of ale2-ele not possible to combine with locsys!\n");
      ilocsys=actnode->locsysId-1;
      if(ilocsys>=0)
      {
         dsassert(ilocsys<numlocsys,"locsysId not existent!\n");
         trans[nd][nd]     = locsys[ilocsys].lXx;
         trans[nd+1][nd]   = locsys[ilocsys].lXy;
         trans[nd][nd+1]   = locsys[ilocsys].lYx;
         trans[nd+1][nd+1] = locsys[ilocsys].lYy;
      }
      else
      {
         for (j=0;j<numdf;j++) trans[nd+j][nd+j] = ONE;
      }
      nd+=numdf;
   } /* end loop over element nodes */
break;
#endif
default: dserror("no transformation implemented for this kind of element!\n");
} /* end switch (ele->eltyp */

/*------ perform the transformation: estif* = trans * estif * trans^t --*/
if (estif1!=NULL)
{
   /* workm = estif1 * trans^t */
   math_matmattrndense(workm,estif1,trans,nd,nd,nd,0,ONE);
   /* estif1* = trans * workm */
   math_matmatdense(estif1,trans,workm,nd,nd,nd,0,ONE);
} /* endif (estif1!=NULL) */

/*------ perform the transformation: estif* = trans * estif * trans^t --*/
if (estif2!=NULL)
{
   /* workm = estif2 * trans^t */
   math_matmattrndense(workm,estif2,trans,nd,nd,nd,0,ONE);
   /* estif2* = trans * workm */
   math_matmatdense(estif2,trans,workm,nd,nd,nd,0,ONE);
} /* endif (estif2!=NULL) */

/*------------------ perform the transformation: eload* = trans * eload */
if (vec1!=NULL)
{
   /* workv = trans * vec1 */
   math_matvecdense(workv,trans,vec1,nd,nd,0,ONE);
   /* copy result to vec1 */
   for(i=0;i<nd;i++) vec1[i]=workv[i];
} /* endif (vec1!=NULL) */

/*------------------ perform the transformation: eload* = trans * eload */
if (vec2!=NULL)
{
   /* workv = trans * vec2 */
   math_matvecdense(workv,trans,vec2,nd,nd,0,ONE);
   /* copy result to vec2 */
   for(i=0;i<nd;i++) vec2[i]=workv[i];
} /* endif (vec2!=NULL) */

/*----------------------------------------------------------------------*/
#ifdef PERF
  perf_end(22);
#endif
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of locsys_trans_ele */


/*======================================================================*/
/*!
\brief Transform element stiffness matrix components which are
       subjected to DirichletBCs defined in local systems
       Different local systems with a single element can be treated.

<pre>
The XYZ-oriented element stiffness matrix of the current element
is transformed to the local x'y'z' system.
This names are a bit confusing, since the "local" co-ordinate system is
the one we are solving in, so we better introduce new names:

estif  = stiffness matrix in the given XYZ cartesian co-system
estif' = stiffness matrix in the alternative x'y'z' co-system defined
         in the input file

eload  = load vector in the given XYZ cartesian co-system
eload' = load vector in the alternative x'y'z' co-system

trans  = transformation matrix containing the direction cosines between
         XYZ and x'y'z'
</pre>

\param   *ele        ELEMENT     (i)      actual element
\param  **estif1     DOUBLE      (i/o)    element matrix
\param  **estif2     DOUBLE      (i/o)    element matrix
\param   *vec1       DOUBLE      (i/o)    element RHS
\param   *vec2       DOUBLE      (i/o)    element RHS
\return void

\author bborn
\date 06/07
*/
#ifdef LOCALSYSTEMS_ST
void locsys_trans_equant_dirich(ELEMENT* ele,
                                ARRAY* estif_global,
                                ARRAY* emass_global,
                                ARRAY* eforce_global,
                                DOUBLE* eforce)
{
  INT nd;  /* total number of element dofs */

#ifdef DEBUG
  dstrc_enter("locsys_trans_equant_dirich");
#endif

  /*--------------------------------------------------------------------*/
  /* initialise matrix */
  amzero(&trans_a);

  /*--------------------------------------------------------------------*/
  /* make rotation matrix */
  switch (ele->eltyp)
  {
#ifdef D_WALL1
  case el_wall1:
  {
    nd = 2 * ele->numnp;  /* total number of element DOFs */
    INT inod;  /* node index */
    for (inod=0; inod<ele->numnp; inod++)
    {
      NODE* actnode = ele->node[inod];
      const INT xdof = 2*inod + 0;  /* index of X-DOF at node */
      const INT ydof = xdof + 1;  /* index of Y-DOF at node */
      if ( (actnode->locsysId > 0) && (actnode->gnode->dirich != NULL) )
      {
        const INT ilocsys = actnode->locsysId - 1;  /* 0-based system ID */
        /* rotation matrix from (XYZ) to (x'y'z') */
        trans[xdof][xdof] = locsys[ilocsys].lXx;
        trans[ydof][xdof] = locsys[ilocsys].lXy;
        trans[xdof][ydof] = locsys[ilocsys].lYx;
        trans[ydof][ydof] = locsys[ilocsys].lYy;
      }
      else
      {
        /* identity matrix */
        trans[xdof][xdof] = 1.0;
        trans[ydof][xdof] = 0.0;
        trans[xdof][ydof] = 0.0;
        trans[ydof][ydof] = 1.0;
      }
    }
    break;
  }
#endif
#ifdef D_BRICK1
  case el_brick1:
  {
    nd = 3 * ele->numnp;  /* total number of element DOFs */
    INT inod;  /* node index */
    for (inod=0; inod<ele->numnp; inod++)
    {
      NODE* actnode = ele->node[inod];
      const INT xdof = 3*inod + 0;  /* index of X-DOF at node */
      const INT ydof = xdof + 1;  /* index of Y-DOF at node */
      const INT zdof = xdof + 2;  /* index of Z-DOF at node */
      if ( (actnode->locsysId > 0) && (actnode->gnode->dirich != NULL) )
      {
        const INT ilocsys = actnode->locsysId - 1;  /* 0-based system ID */
        /* rotation matrix from (XYZ) to (x'y'z') */
        trans[xdof][xdof] = locsys[ilocsys].lXx;
        trans[ydof][xdof] = locsys[ilocsys].lXy;
        trans[zdof][xdof] = locsys[ilocsys].lXz;
        trans[xdof][ydof] = locsys[ilocsys].lYx;
        trans[ydof][ydof] = locsys[ilocsys].lYy;
        trans[zdof][ydof] = locsys[ilocsys].lYz;
        trans[xdof][zdof] = locsys[ilocsys].lZx;
        trans[ydof][zdof] = locsys[ilocsys].lZy;
        trans[zdof][zdof] = locsys[ilocsys].lZz;
      }
      else
      {
        /* identity matrix */
        trans[xdof][xdof] = 1.0;
        trans[ydof][xdof] = 0.0;
        trans[zdof][xdof] = 0.0;
        trans[xdof][ydof] = 0.0;
        trans[ydof][ydof] = 1.0;
        trans[zdof][ydof] = 0.0;
        trans[xdof][zdof] = 0.0;
        trans[ydof][zdof] = 0.0;
        trans[zdof][zdof] = 1.0;
      }
    }
    break;
  }
#endif
#ifdef D_SOLID3
  case el_solid3:
  {
    nd = ele->numnp*NUMDOF_SOLID3;  /* total number of element DOFs */
    INT inod;  /* node index */
    for (inod=0; inod<ele->numnp; inod++)
    {
      NODE* actnode = ele->node[inod];
      const INT xdof = NUMDOF_SOLID3*inod + 0;  /* index of X-DOF at node */
      const INT ydof = xdof + 1;  /* index of Y-DOF at node */
      const INT zdof = xdof + 2;  /* index of Z-DOF at node */
      if ( (actnode->locsysId > 0) && (actnode->gnode->dirich != NULL) )
      {
        const INT ilocsys = actnode->locsysId - 1;  /* 0-based system ID */
        /* rotation matrix from (XYZ) to (x'y'z') */
        trans[xdof][xdof] = locsys[ilocsys].lXx;
        trans[ydof][xdof] = locsys[ilocsys].lXy;
        trans[zdof][xdof] = locsys[ilocsys].lXz;
        trans[xdof][ydof] = locsys[ilocsys].lYx;
        trans[ydof][ydof] = locsys[ilocsys].lYy;
        trans[zdof][ydof] = locsys[ilocsys].lYz;
        trans[xdof][zdof] = locsys[ilocsys].lZx;
        trans[ydof][zdof] = locsys[ilocsys].lZy;
        trans[zdof][zdof] = locsys[ilocsys].lZz;
      }
      else
      {
        /* identity matrix */
        trans[xdof][xdof] = 1.0;
        trans[ydof][xdof] = 0.0;
        trans[zdof][xdof] = 0.0;
        trans[xdof][ydof] = 0.0;
        trans[ydof][ydof] = 1.0;
        trans[zdof][ydof] = 0.0;
        trans[xdof][zdof] = 0.0;
        trans[ydof][zdof] = 0.0;
        trans[zdof][zdof] = 1.0;
      }
    }
    break;
  }
#endif
  default:
    dserror("no transformation implemented for this kind of element!\n");
  } /* end switch (ele->eltyp */

  /*--------------------------------------------------------------------*/
  /* perform the transformation: estif* = trans * estif * trans^t */
  if (estif_global != NULL)
  {
    DOUBLE** estif1 = estif_global->a.da;
    /* workm = estif1 * trans^t */
    math_matmattrndense(workm, estif1, trans, nd, nd, nd, 0, ONE);
    /* estif1* = trans * workm */
    math_matmatdense(estif1, trans, workm, nd, nd, nd, 0, ONE);
  }

  /* perform the transformation: estif* = trans * estif * trans^t */
  if (emass_global != NULL)
  {
    DOUBLE** estif2 = emass_global->a.da;
    /* workm = estif2 * trans^t */
    math_matmattrndense(workm, estif2, trans, nd, nd, nd, 0, ONE);
    /* estif2* = trans * workm */
    math_matmatdense(estif2, trans, workm, nd, nd, nd, 0, ONE);
  }

  /* perform the transformation: eload* = trans * eload */
  if (eforce_global != NULL)
  {
    DOUBLE* vec1 = eforce_global->a.dv;
    /* workv = trans * vec1 */
    math_matvecdense(workv, trans, vec1, nd, nd, 0, ONE);
    /* copy result to vec1 */
    INT i = 0;
    for(i=0; i<nd; i++)
    {
      vec1[i] = workv[i];
    }
  }

  /* perform the transformation: eload* = trans * eload */
  if (eforce != NULL)
  {
    DOUBLE* vec2 = eforce;
    /* workv = trans * vec2 */
    math_matvecdense(workv, trans, vec2, nd, nd, 0, ONE);
    /* copy result to vec2 */
    INT i = 0;
    for(i=0; i<nd; i++)
    {
      vec2[i] = workv[i];
    }
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}
#endif

/*!---------------------------------------------------------------------
\brief transform solution to global co-ordinate system

<pre>                                                         genk 04/04

the solution of the actual node may have to be transformed between
the global XYZ - co-ordinate system to the local one.
This names are a bit confusing, since the "local" co-ordinate system is
the one we are solving in, so we better introduce new names:

sol  = sol vector in the given XYZ cartesian co-system
sol* = sol vector in the alternative xyz* co-system

array   index of the array,        0 = sol
                                   1 = sol_increment
                                   2 = sol_residual
                                   3 = sol_mf

flag = 1 : transform sol in xyz* to XYZ
flag = 0 : transform sol in XYZ to xyz*
HINT: For this flag the type LOCSYS_TRF_KIND is defined in locsys.h

</pre>
\param   *actfield      FIELD       (i)   actual field
\param    dis           INT         (i/o) index for dis
\param    array         INT         (i)   index for nodal array
\param    place         INT         (i)   place in nodal array
\param    flag          INT         (i)   just a flag
\sa  locsys_trans_nodval()

\return void

------------------------------------------------------------------------*/
void locsys_trans_sol(FIELD *actfield, INT idis, INT array,
                      INT place, INT flag)
{
INT      i,j;             /* simply some counters */
INT      iloccsys;        /* index of locsys */
INT      numdf;           /* number of dofs at node */
INT      numnp_total;     /* total number of nodes in field */
DOUBLE   **nodalsol=NULL; /* pointer to nodal array */
NODE     *actnode;        /* actual node */
ELEMENT  *actele;         /* actual element */

#ifdef DEBUG
dstrc_enter("locsys_trans_sol");
#endif

if (numlocsys==0) goto end;

numnp_total=actfield->dis[idis].numnp;

for (i=0;i<numnp_total;i++)
{
   actnode=&(actfield->dis[idis].node[i]);
   /*----- any element can be used to find the local co-ordinate system */
   iloccsys=actnode->locsysId-1;
   if (iloccsys<0) continue; /* no locsys for this node */
   /*----- any element can be used to find the local co-ordinate system */
   actele=actnode->element[0];
   numdf=actnode->numdf;
   /*----------------------------------- fill the local solution vector */
   switch (array)
   {
   case 0:  nodalsol=actnode->sol.a.da;           break;
   case 1:  nodalsol=actnode->sol_increment.a.da; break;
   case 2:  nodalsol=actnode->sol_residual.a.da;  break;
   case 3:  nodalsol=actnode->sol_mf.a.da;        break;
   default: dserror("index out of range!\n");     break;
   } /* end switch (array) */

   /* copy values to working vector */
   for (j=0;j<numdf;j++) nodalwork[j]=nodalsol[place][j];
   /* transform vector */
   locsys_trans_nodval(actele,nodalwork,numdf,iloccsys,flag);
   /* copy result back to nodal sol-field */
   for (j=0;j<numdf;j++) nodalsol[place][j]=nodalwork[j];
} /* end loop numnp_total */

/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of locsys_trans_sol */

/*!---------------------------------------------------------------------
\brief transform solution to global co-ordinate system

<pre>                                                         genk 04/04

the solution of the actual node may have to be transformed between
the global XYZ - co-ordinate system to the local one.
Here only nodes with a dirichlet condition are tranformed.
This names are a bit confusing, since the "local" co-ordinate system is
the one we are solving in, so we better introduce new names:

sol  = sol vector in the given XYZ cartesian co-system
sol* = sol vector in the alternative xyz* co-system

array   index of the array,        0 = sol
                                   1 = sol_increment
                                   2 = sol_residual
                                   3 = sol_mf

flag = 1 : transform sol in xyz* to XYZ
flag = 0 : transform sol in XYZ to xyz*
HINT: For this flag the type LOCSYS_TRF_KIND is defined in locsys.h

</pre>
\param   *actfield      FIELD       (i)   actual field
\param    dis           INT         (i/o) index for dis
\param    array         INT         (i)   index for nodal array
\param    place         INT         (i)   place in nodal array
\param    flag          INT         (i)   just a flag
\sa  locsys_trans_nodval()

\return void

------------------------------------------------------------------------*/
void locsys_trans_sol_dirich(FIELD *actfield, INT idis, INT array,
                             INT place, INT flag)
{
INT      i,j;           /* simply some counters */
INT      iloccsys;      /* index of locsys */
INT      numdf;         /* number of dofs at node */
INT      numnp_total;   /* total number of nodes in field */
DOUBLE   **nodalsol;    /* pointer to nodal array */
NODE     *actnode;      /* actual node */
GNODE    *actgnode;     /* actual gnode */
ELEMENT  *actele;       /* actual element */

#ifdef DEBUG
dstrc_enter("locsys_trans_sol_dirich");
#endif

if (numlocsys==0) goto end;

numnp_total=actfield->dis[idis].numnp;

for (i=0;i<numnp_total;i++)
{
   actnode=&(actfield->dis[idis].node[i]);
   actgnode=actnode->gnode;
   if (actgnode->dirich==NULL) continue;
   /*----- any element can be used to find the local co-ordinate system */
   iloccsys=actnode->locsysId-1;
   if (iloccsys<0) continue; /* no locsys for this node */
   /*----- any element can be used to find the local co-ordinate system */
   actele=actnode->element[0];
   numdf=actnode->numdf;
   /*----------------------------------- fill the local solution vector */
   switch (array)
   {
     case 0:
       nodalsol=actnode->sol.a.da;           break;
     case 1:
       nodalsol=actnode->sol_increment.a.da; break;
     case 2:
       nodalsol=actnode->sol_residual.a.da;  break;
     case 3:
       nodalsol=actnode->sol_mf.a.da;        break;
     default:
       nodalsol=NULL;
       dserror("index out of range!\n");
       break;
   } /* end switch (array) */

   /* copy values to working vector */
   for (j=0;j<numdf;j++) nodalwork[j]=nodalsol[place][j];
   /* transform vector */
   locsys_trans_nodval(actele,nodalwork,numdf,iloccsys,flag);
   /* copy result back to nodal sol-field */
   for (j=0;j<numdf;j++) nodalsol[place][j]=nodalwork[j];
} /* end loop numnp_total */

/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of locsys_trans_sol */

/*!---------------------------------------------------------------------
\brief transform solution to global co-ordinate system

<pre>                                                         genk 04/04

the solution of the actual node is transformed from
the global XYZ - co-ordinate system to the local one.
This names are a bit confusing, since the "local" co-ordinate system is
the one we are solving in, so we better introduce new names:

    Transformation of displacements (3D):

     | (Dx*) |   | cos(Xx*)   cos(Yx*)   cos(Zx*) | | (DX) |
     | (Dy*) | = | cos(Xy*)   cos(Yy*)   cos(Zy*) | | (DY) |
     | (Dz*) |   | cos(Xz*)   cos(Yz*)   cos(Zz*) | | (DZ) |

      val*     =                  T                   val

cos(Yx*) =  cosine of angle between the given Y-base vector
            and the alternative x*-base vector

val   Displ./Vel. vector in the given XYZ cartesian co-system
val*  Displ./Vel. vector in the alternative xyz* co-system
T     Transformation matrix

flag = 1 : transform val in xyz* to XYZ
flag = 0 : transform val in XYZ to xyz*
HINT: For this flag the type LOCSYS_TRF_KIND is defined in locsys.h

</pre>
\param   *actele     ELEMENT        (i)   actual element
\param   *val        DOUBLE         (i/o) values to transform
\param    numdf      INT            (i)   number of dofs
\param    ilocsys    INT            (i)   index of locsys
\param    flag       INT            (i)   just a flag

\sa inp_read_locsys()
\return void

------------------------------------------------------------------------*/
void locsys_trans_nodval(ELEMENT *actele, DOUBLE *val, INT numdf,
                         INT iloccsys, INT flag)
{
INT      i,j;           /* simply some counters */
LOCSYS   *actlocsys;    /* the actual local co-system xzy* */

#ifdef DEBUG
dstrc_enter("locsys_trans_nodval");
#endif
#ifdef PERF
  perf_begin(22);
#endif
if (numlocsys==0) goto end;

/*------------------------------------------------- initialise matrices */
for(i=0;i<numdf;i++)
   for(j=0;j<numdf;j++)
      trans[i][j]=ZERO;

dsassert(iloccsys>=0,"locsysId out of range!\n");
dsassert(iloccsys<numlocsys,"locsysId not existent!\n");
actlocsys=&(locsys[iloccsys]);

/*-------------------------------- fill the nodal transformation matrix */
switch (actele->eltyp)
{
#ifdef D_FLUID2
case el_fluid2:
   if (numdf<4)
   {
      trans[0][0] = actlocsys->lXx;
      trans[1][0] = actlocsys->lXy;
      trans[0][1] = actlocsys->lYx;
      trans[1][1] = actlocsys->lYy;
      /*-- don't transform pressure dof, so reduce number of nodal dofs */
      numdf--;
   }
   else if (numdf==4)
   {
      dserror("transformation for fluid node with 4 dofs not implemented!\n");
   }
   else if (numdf==5) /*  node at free surf. w/ five dofs
                          [vel, vel, pre, velg, velg]                   */
   {
      trans[0][0] = actlocsys->lXx;
      trans[1][0] = actlocsys->lXy;
      trans[0][1] = actlocsys->lYx;
      trans[1][1] = actlocsys->lYy;
      trans[2][2] = ONE;
      trans[3][3] = actlocsys->lXx;
      trans[4][3] = actlocsys->lXy;
      trans[3][4] = actlocsys->lYx;
      trans[4][4] = actlocsys->lYy;
   }
   else
      dserror("transformation not possible!\n");
break;
#endif
#ifdef D_FLUID3
case el_fluid3:
   if (numdf<5)
   {
      trans[0][0] = actlocsys->lXx;
      trans[1][0] = actlocsys->lXy;
      trans[2][0] = actlocsys->lXz;
      trans[0][1] = actlocsys->lYx;
      trans[1][1] = actlocsys->lYy;
      trans[2][1] = actlocsys->lYz;
      trans[0][2] = actlocsys->lZx;
      trans[1][2] = actlocsys->lZy;
      trans[2][2] = actlocsys->lZz;
      trans[3][3] = ONE;
   }
   else if (numdf==5)
   {
      dserror("transformation for fluid node with 5 dofs not implemented!\n");
   }
   else if (numdf==7) /*  node at free surf. w/ five dofs
                    [vel, vel, vel, pre, velg, velg, velg]        */
   {
      trans[0][0] = actlocsys->lXx;
      trans[1][0] = actlocsys->lXy;
      trans[2][0] = actlocsys->lXz;
      trans[0][1] = actlocsys->lYx;
      trans[1][1] = actlocsys->lYy;
      trans[2][1] = actlocsys->lYz;
      trans[0][2] = actlocsys->lZx;
      trans[1][2] = actlocsys->lZy;
      trans[2][2] = actlocsys->lZz;
      trans[3][3] = ONE;
      trans[4][4] = actlocsys->lXx;
      trans[5][4] = actlocsys->lXy;
      trans[6][4] = actlocsys->lXz;
      trans[4][5] = actlocsys->lYx;
      trans[5][5] = actlocsys->lYy;
      trans[6][5] = actlocsys->lYz;
      trans[4][6] = actlocsys->lZx;
      trans[5][6] = actlocsys->lZy;
      trans[6][6] = actlocsys->lZz;
   }
   else
      dserror("transformation not possible!\n");
break;
#endif
#ifdef D_ALE
case el_ale2:
   trans[0][0] = actlocsys->lXx;
   trans[1][0] = actlocsys->lXy;
   trans[0][1] = actlocsys->lYx;
   trans[1][1] = actlocsys->lYy;
break;
#endif
#ifdef D_WALL1
case el_wall1:
   trans[0][0] = actlocsys->lXx;
   trans[1][0] = actlocsys->lXy;
   trans[0][1] = actlocsys->lYx;
   trans[1][1] = actlocsys->lYy;
break;
#endif
#ifdef D_BRICK1
case el_brick1:
   trans[0][0] = actlocsys->lXx;
   trans[1][0] = actlocsys->lXy;
   trans[2][0] = actlocsys->lXz;
   trans[0][1] = actlocsys->lYx;
   trans[1][1] = actlocsys->lYy;
   trans[2][1] = actlocsys->lYz;
   trans[0][2] = actlocsys->lZx;
   trans[1][2] = actlocsys->lZy;
   trans[2][2] = actlocsys->lZz;
break;
#endif
#ifdef D_SOLID3
case el_solid3:
   trans[0][0] = actlocsys->lXx;
   trans[1][0] = actlocsys->lXy;
   trans[2][0] = actlocsys->lXz;
   trans[0][1] = actlocsys->lYx;
   trans[1][1] = actlocsys->lYy;
   trans[2][1] = actlocsys->lYz;
   trans[0][2] = actlocsys->lZx;
   trans[1][2] = actlocsys->lZy;
   trans[2][2] = actlocsys->lZz;
break;
#endif
default: dserror("no transformation implemented for this kind of element!\n");
} /* end switch (actele->eltyp) */

switch (flag)
{
case 0:/* transformation: val* = trans * val */
   math_matvecdense(workv,trans,val,numdf,numdf,0,ONE);
break;
case 1: /* transformation: val = trans^t * val* */
   math_mattrnvecdense(workv,trans,val,numdf,numdf,0,ONE);
break;
default:
   dserror("flag out of range!\n");
} /* end switch (flag) */

/* copy result back to val */
for (j=0;j<numdf;j++) val[j]=workv[j];

/*----------------------------------------------------------------------*/
end:
#ifdef PERF
  perf_end(22);
#endif
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of locsys_trans_nodval */
#endif
