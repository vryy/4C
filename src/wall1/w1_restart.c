/*!----------------------------------------------------------------------
\file
\brief contains the routines
 - w1_write_restart:  writes all the element data that is needed for a restart
 - w1_read_restart:   reads all the element data that is needed for a restart


<pre>
Maintainer: Stefan Hartmann
            hartmann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hartmann/
            0771 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;


/*!----------------------------------------------------------------------
\brief write the data needed to restart this element

<pre>                                                           sh 03/04
This routine writes the data needed to restart the shell9 element.
</pre>
\param  ELEMENT  *actele   (i) actual element
\param  MATERIAL *mat      (i)  the material structure
\param  INT       nhandle  (i) size of handles
\param  long int *handles  ( ) unique handle returned by the pss-system
\param  INT       init    (i) flag for initializing arrays

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: wall1()   [w1_main.c]

*----------------------------------------------------------------------*/
void w1_write_restart(ELEMENT *actele, MATERIAL  *mat, INT nhandle, long int *handles, INT init)
{
INT j,n,k,ierr;
INT size_j;
static ARRAY   elewares_a;  static DOUBLE **elewares;     /* array for elewa */
FILE *out;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("w1_write_restart");
#endif
/*---------------------------------------------------------- init phase */
if (init==1)
{
   /*array for elewares_a
     fdim = 1
     sdim = max number of possible entries in ipwa per GP (30 at the moment)
            * max number of GP (3x3 = 9 at the moment) */

    /*  pl_mises */
    /*  DOUBLE       epstn;      equivalent strain                            */
    /*  INT          yip;        stress state: 1=elastic 2=plastic            */
    /*  DOUBLE       sig[4];     stresses                               [4]   */
    /*  DOUBLE       eps[4]      strains                                [4]   */
    /*  DOUBLE       qn[4]       backstress                             [4]   */
    /*      ( == 14 )                                                         */

    /*  pl_mises_3D */
    /*  DOUBLE       epstn;      equivalent strain                            */
    /*  INT          yip;        stress state: 1=elastic 2=plastic            */
    /*  DOUBLE       sig[4];     stresses                               [4]   */
    /*  DOUBLE       eps[4]      strains                                [4]   */
    /*  DOUBLE       qn[4]       backstress                             [4]   */
    /*  DOUBLE       sigc[4];    stresses condensed                     [4]   */
    /*  DOUBLE       sigi[4];    stresses for condensation              [4]   */
    /*  DOUBLE       epsi[4]     strains for condensation               [4]   */
    /*  DOUBLE       di[4]       part of constitutive matrix for cond.  [4]   */
    /*      ( == 30 )                                                         */

    /*  pl_epc3D */
    /*  INT          yip;        stress state: 1=elastic 2=plastic            */
    /*  DOUBLE       dlam[0];    equivalent strain (tension)                  */
    /*  DOUBLE       dlam[1];    equivalent strain (compression)              */
    /*  DOUBLE       sig[4];     stresses                               [4]   */
    /*  DOUBLE       eps[4]      strains                                [4]   */
    /*  DOUBLE       sigc[4];    stresses condensed                     [4]   */
    /*  DOUBLE       grad[4];    total gradient                         [4]   */
    /*  DOUBLE       sigi[4];    stresses for condensation              [4]   */
    /*  DOUBLE       epsi[4]     strains for condensation               [4]   */
    /*  DOUBLE       di[4]       part of constitutive matrix for cond.  [4]   */
    /*      ( == 31 )                                                         */

   elewares     = amdef("elewares"  ,&elewares_a,1,31*9,"DA");

   goto end;
}
/*-------------------------------------------------------- uninit phase */
else if (init==-1)
{
   amdel(&elewares_a);

   goto end;
}
/*----------------------------------------------------------------------*/
/*
   the wall element has to write
   - elewa for material nonlinearity
*/
/*----------------------------------------------------------------------*/
out = allfiles.out_pss;
/*----------------------------------------------------------------------*/
handles[0]=0;

/*----------------------------------------------------------------------*/
/*--------- now check for ELEWA --->  material nonlinearity ------------*/
/*check if there is a working array for this element*/
if (!actele->e.w1->elewa) goto end; /*no elewa to this element*/

/*if there is a working array, write all the data into the elewares_a */
/*number of gausspoints to this element*/
size_j = actele->e.w1->nGP[0]*actele->e.w1->nGP[1];

n = 0;

/* pl_mises */
if     (mat->mattyp == m_pl_mises )
{
  for (k=0; k<size_j; k++)/*read for every gausspoint*/
  {
    elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].epstn;
    n++;
    elewares[0][n] = (double) actele->e.w1->elewa[0].ipwa[k].yip;
    n++;
    for (j=0; j<4; j++)
    {
      elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].sig[j];
      n++;
      elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].eps[j];
      n++;
      elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].qn[j];
      n++;
    }
  }
} /*end pl_mises */

/* pl_mises_3D */
else if(mat->mattyp == m_pl_mises_3D )
{
  for (k=0; k<size_j; k++)/*read for every gausspoint*/
  {
    elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].epstn;
    n++;
    elewares[0][n] = (double) actele->e.w1->elewa[0].ipwa[k].yip;
    n++;
    for (j=0; j<4; j++)
    {
      elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].sig[j];
      n++;
      elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].eps[j];
      n++;
      elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].qn[j];
      n++;
      elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].sigc[j];
      n++;
      elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].sigi[j];
      n++;
      elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].epsi[j];
      n++;
      elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].di[j];
      n++;
    }
  }
} /*end pl_mises_3D */

/* pl_epc_3D */
else if(mat[actele->mat-1].mattyp == m_pl_epc3D )
{
  for (k=0; k<size_j; k++)/*read for every gausspoint*/
  {
    elewares[0][n] = (double) actele->e.w1->elewa[0].ipwa[k].yip;
    n++;
    elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].dlam[0];
    n++;
    elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].dlam[1];
    n++;
    for (j=0; j<4; j++)
    {
      elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].sig[j];
      n++;
      elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].eps[j];
      n++;
      elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].sigc[j];
      n++;
      elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].grad[j];
      n++;
      elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].sigi[j];
      n++;
      elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].epsi[j];
      n++;
      elewares[0][n] = actele->e.w1->elewa[0].ipwa[k].di[  j];
      n++;
    }
  }
} /*end pl_epc3D */

else printf("WARNING: Unknown mattyp in 'w1_write_restart'\n");

handles[0]+=1;
if (handles[0]+1 > nhandle) dserror("Handle range too small for element in 'w1_write_restart' ");
/*-------------------------------------------------- write elewares */
pss_write_array(&elewares_a,&(handles[1]),out,&ierr);
if (ierr != 1) dserror("Error writing restart data in 'w1_write_restart' ");


/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_write_restart */

/*!----------------------------------------------------------------------
\brief read the data needed to restart this element

<pre>                                                            sh 03/04
This routine reads the data needed to restart the wall1 element.
</pre>
\param  ELEMENT  *actele   (i) actual element
\param  MATERIAL *mat      (i)  the material structure
\param  INT       nhandle  (i) size of handles
\param  long int *handles  ( ) unique handle returned by the pss-system
\param  INT       init    (i) flag for initializing arrays

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: wall1()   [w1_main.c]

*----------------------------------------------------------------------*/
void w1_read_restart(ELEMENT *actele, MATERIAL  *mat, long int *handles, INT init)
{
INT j,n,k,ierr;
INT size_j;
static ARRAY   elewares_a;  static DOUBLE **elewares;     /* array for elewa */
INT dims[3];
FILE *in;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("w1_read_restart");
#endif
/*---------------------------------------------------------- init phase */
if (init==1)
{
   /*for dimensions see w1_write_restart */

   elewares     = amdef("elewares"  ,&elewares_a,1,31*9,"DA");

   goto end;
}
/*-------------------------------------------------------- uninit phase */
else if (init==-1)
{
   amdel(&elewares_a);

   goto end;
}
/*----------------------------------------------------------------------*/
/*
   the wall element has to read
   - elewa for material nonlinearity
*/
/*----------------------------------------------------------------------*/
in = allfiles.in_pss;
if (!in) dserror("There is no restart input file open in 'w1_read_restart'");
/*----------------------------------------------------------------------*/
/*--------- now check for ELEWA --->  material nonlinearity ------------*/
/*check if there is a working array for this element*/
if (!actele->e.w1->elewa) goto end; /*no elewa to this element*/

/*if there is a working array, read all the data into the elewares_a */
/*------------------------------------------------- check dimensions elewares_a*/
pss_getdims_name_handle(elewares_a.name,&dims[0],&dims[1],&dims[2],&handles[1],in,&ierr);
if (ierr != 1) dserror("Cannot read restart data in 'w1_read_restart'");
if (elewares_a.fdim != dims[0] ||
    elewares_a.sdim != dims[1])
    dserror("Mismatch in reading element restart data in 'w1_read_restart'");
/*-------------------------------------------------------- read elewares_a */
pss_read_array_name_handle(elewares_a.name,&elewares_a,&handles[1],in,&ierr);
if (ierr != 1) dserror("Cannot read restart data in 'w1_read_restart'");


/*----- now write the restart data back to elewa -------------------------*/

/*number of gausspoints to this element*/
size_j = actele->e.w1->nGP[0]*actele->e.w1->nGP[1];

n = 0;

/* pl_mises */
if     (mat->mattyp == m_pl_mises )
{
  for (k=0; k<size_j; k++)/*read for every gausspoint*/
  {
    actele->e.w1->elewa[0].ipwa[k].epstn   = elewares[0][n];
    n++;
    actele->e.w1->elewa[0].ipwa[k].yip     = (int) elewares[0][n];
    n++;
    for (j=0; j<4; j++)
    {
      actele->e.w1->elewa[0].ipwa[k].sig[j]    = elewares[0][n];
      n++;
      actele->e.w1->elewa[0].ipwa[k].eps[j]    = elewares[0][n];
      n++;
      actele->e.w1->elewa[0].ipwa[k].qn[j]     = elewares[0][n];
      n++;
    }
  }
} /*end pl_mises */

/* pl_mises_3D */
else if(mat->mattyp == m_pl_mises_3D )
{
  for (k=0; k<size_j; k++)/*read for every gausspoint*/
  {
    actele->e.w1->elewa[0].ipwa[k].epstn   = elewares[0][n];
    n++;
    actele->e.w1->elewa[0].ipwa[k].yip     = (int) elewares[0][n];
    n++;
    for (j=0; j<4; j++)
    {
      actele->e.w1->elewa[0].ipwa[k].sig[j]    = elewares[0][n];
      n++;
      actele->e.w1->elewa[0].ipwa[k].eps[j]    = elewares[0][n];
      n++;
      actele->e.w1->elewa[0].ipwa[k].qn[j]     = elewares[0][n];
      n++;
      actele->e.w1->elewa[0].ipwa[k].sigc[j]   = elewares[0][n];
      n++;
      actele->e.w1->elewa[0].ipwa[k].sigi[j]   = elewares[0][n];
      n++;
      actele->e.w1->elewa[0].ipwa[k].epsi[j]   = elewares[0][n];
      n++;
      actele->e.w1->elewa[0].ipwa[k].di[  j]   = elewares[0][n];
      n++;
    }
  }
} /*end pl_mises_3D */

/* pl_epc_3D */
else if(mat[actele->mat-1].mattyp == m_pl_epc3D )
{
  for (k=0; k<size_j; k++)/*read for every gausspoint*/
  {
    actele->e.w1->elewa[0].ipwa[k].yip     = (int) elewares[0][n];
    n++;
    actele->e.w1->elewa[0].ipwa[k].dlam[0] = elewares[0][n];
    n++;
    actele->e.w1->elewa[0].ipwa[k].dlam[1] = elewares[0][n];
    n++;
    for (j=0; j<4; j++)
    {
      actele->e.w1->elewa[0].ipwa[k].sig[j]    = elewares[0][n];
      n++;
      actele->e.w1->elewa[0].ipwa[k].eps[j]    = elewares[0][n];
      n++;
      actele->e.w1->elewa[0].ipwa[k].sigc[j]   = elewares[0][n];
      n++;
      actele->e.w1->elewa[0].ipwa[k].grad[j]   = elewares[0][n];
      n++;
      actele->e.w1->elewa[0].ipwa[k].sigi[j]   = elewares[0][n];
      n++;
      actele->e.w1->elewa[0].ipwa[k].epsi[j]   = elewares[0][n];
      n++;
      actele->e.w1->elewa[0].ipwa[k].di[  j]   = elewares[0][n];
      n++;
    }
  }
} /*end pl_epc3D */

else printf("WARNING: Unknown mattyp in 'w1_write_restart'\n");

/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_read_restart */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
