/*!-----------------------------------------------------------------------
\file
\brief contains the routines 
 - if_write_restart:  writes all the element data that is needed for a restart
 - if_read_restart:   reads all the element data that is needed for a restart
<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>

*-----------------------------------------------------------------------*/
#ifdef D_INTERF
#include "../headers/standardtypes.h"
#include "interf.h"
#include "interf_prototypes.h"

/*! 
\addtogroup INTERF
*/
/*! @{ (documentation module open)*/

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

<pre>                                                           ah 06/04
This routine writes the data needed to restart the shell9 element. 
</pre>
\param  ELEMENT  *actele   (i) actual element
\param  INT       nhandle  (i) size of handles
\param  long int *handles  ( ) unique handle returned by the pss-system 
\param  INT       init    (i) flag for initializing arrays

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: interf()   [if_main.c]

*----------------------------------------------------------------------*/

void if_write_restart(ELEMENT   *actele, 
                      MATERIAL  *mat,
                      INT        nhandle, 
                      long int  *handles, 
                      INT        init)
{
INT n,k,ierr;
INT size_j;
static ARRAY   elewares_a;  static DOUBLE **elewares;     /* array for elewa */
FILE *out;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("if_write_restart");
#endif
/*--------------------------------------------------------------------- */
/*---------------------------------------------------------- init phase */
if (init==1)
{
   /*array for elewares_a  
     fdim = 1
     sdim = sum of all possible entries in ipwa per GP (10 at the moment) 
            * max number of GP (3 at the moment) */

    /*  DOUBLE    Tt;          tangential stresses                      */
    /*  DOUBLE    Tn;          normal stresses                          */
    /*  DOUBLE    dt;          tangential damage                        */
    /*  DOUBLE    dn;          normal damage                            */
    /*  INT       yip;         stress state: 1=elastic 2=plastic        */
    /*  DOUBLE   jump_ut_pl;   tangential plastic displacement jump     */
    /*  DOUBLE   Q[2][2];      material tangente                        */
    /*      ( == 10 )                                                       */


   elewares     = amdef("elewares"  ,&elewares_a,1,10*3,"DA");    

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
   the interface element has to write
   - elewa for material nonlinearity
*/
/*----------------------------------------------------------------------*/
out = allfiles.out_pss;
/*----------------------------------------------------------------------*/
handles[0]=0;
/*----------------------------------------------------------------------*/
/*--------- now check for ELEWA --->  material nonlinearity ------------*/
/*check if there is a working array for this element*/
if (!actele->e.interf->elewa) goto end; /*no elewa to this element*/

/*if there is a working array, write all the data into the elewares_a */
/*number of gausspoints to this element*/     
size_j = actele->e.interf->nGP;

n = 0;
for (k=0; k<size_j; k++)/*write for every gausspoint -> elewares*/
{
  elewares[0][n] = actele->e.interf->elewa[0].ipwa[k].Tt;
  n++;
  elewares[0][n] = actele->e.interf->elewa[0].ipwa[k].Tn;
  n++;
  elewares[0][n] = actele->e.interf->elewa[0].ipwa[k].dt;
  n++;
  elewares[0][n] = actele->e.interf->elewa[0].ipwa[k].dn;
  n++;
  elewares[0][n] = actele->e.interf->elewa[0].ipwa[k].yip;
  n++;
  elewares[0][n] = actele->e.interf->elewa[0].ipwa[k].jump_ut_pl;
  n++;
  elewares[0][n] = actele->e.interf->elewa[0].ipwa[k].Q[0][0];
  n++;
  elewares[0][n] = actele->e.interf->elewa[0].ipwa[k].Q[0][1];
  n++;
  elewares[0][n] = actele->e.interf->elewa[0].ipwa[k].Q[1][0];
  n++;
  elewares[0][n] = actele->e.interf->elewa[0].ipwa[k].Q[1][1];
  n++;
}

handles[0]+=1;
if (handles[0]+1 > nhandle) dserror("Handle range too small for element");
/*-------------------------------------------------- write elewares */
pss_write_array(&elewares_a,&(handles[1]),out,&ierr);
if (ierr != 1) dserror("Error writing restart data");

/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of if_write_restart */



/*!----------------------------------------------------------------------
\brief read the data needed to restart this element                                       

<pre>                                                            ah 06/04
This routine reads the data needed to restart the wall1 element. 
</pre>
\param  ELEMENT  *actele   (i) actual element
\param  INT       nhandle  (i) size of handles
\param  long int *handles  ( ) unique handle returned by the pss-system 
\param  INT       init    (i) flag for initializing arrays

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: wall1()   [w1_main.c]

*----------------------------------------------------------------------*/
void if_read_restart(ELEMENT  *actele, 
                     MATERIAL *mat, 
                     long int *handles, 
                     INT       init)
{
INT n,k,ierr;
INT size_j;
static ARRAY   elewares_a;  static DOUBLE **elewares;     /* array for elewa */
INT dims[3];
FILE *in;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("if_read_restart");
#endif
/*---------------------------------------------------------- init phase */
if (init==1)
{
   /*for dimensions see w1_write_restart */
               
   elewares     = amdef("elewares"  ,&elewares_a,1,10*3,"DA");    

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
   the interface element has to read
   - elewa for material nonlinearity
*/
/*----------------------------------------------------------------------*/
in = allfiles.in_pss;
if (!in) dserror("There is no restart input file open");
/*--------- now check for ELEWA --->  material nonlinearity ------------*/
/*check if there is a working array for this element*/
if (!actele->e.interf->elewa) goto end; /*no elewa to this element*/

/*if there is a working array, read all the data into the elewares_a */
/*------------------------------------------------- check dimensions elewares_a*/
pss_getdims_name_handle(elewares_a.name,&dims[0],&dims[1],&dims[2],&handles[1],in,&ierr);
if (ierr != 1) dserror("Cannot read restart data");
if (elewares_a.fdim != dims[0] ||
    elewares_a.sdim != dims[1])
    dserror("Mismatch in reading element restart data");
/*-------------------------------------------------------- read elewares_a */
pss_read_array_name_handle(elewares_a.name,&elewares_a,&handles[1],in,&ierr);
if (ierr != 1) dserror("Cannot read restart data");  


/*----- now write the restart data back to elewa -------------------------*/

/*number of gausspoints to this element*/     
size_j = actele->e.interf->nGP;
n = 0;
for (k=0; k<size_j; k++)/*read for every gausspoint*/
{
  actele->e.interf->elewa[0].ipwa[k].Tt         = elewares[0][n];
  n++;
  actele->e.interf->elewa[0].ipwa[k].Tn         = elewares[0][n];
  n++;
  actele->e.interf->elewa[0].ipwa[k].dt         = elewares[0][n];
  n++;
  actele->e.interf->elewa[0].ipwa[k].dn         = elewares[0][n];
  n++;
  actele->e.interf->elewa[0].ipwa[k].yip        = elewares[0][n];
  n++;
  actele->e.interf->elewa[0].ipwa[k].jump_ut_pl = elewares[0][n];
  n++;
  actele->e.interf->elewa[0].ipwa[k].Q[0][0]    = elewares[0][n];
  n++;
  actele->e.interf->elewa[0].ipwa[k].Q[0][1]    = elewares[0][n];
  n++;
  actele->e.interf->elewa[0].ipwa[k].Q[1][0]    = elewares[0][n];
  n++;
  actele->e.interf->elewa[0].ipwa[k].Q[1][1]    = elewares[0][n];
  n++;
}


/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of if_read_restart */

/*----------------------------------------------------------------------*/
#endif /*D_INTERF*/
/*! @} (documentation module close)*/

