/*!----------------------------------------------------------------------
\file
\brief contains the routines
 - s9_write_restart:  writes all the element data that is needed for a restart
 - s9_read_restart:   reads all the element data that is needed for a restart


<pre>
Maintainer: Stefan Hartmann
            hartmann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hartmann/
            0771 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
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
\brief vector of multilayer material law

<pre>                                                            sh 10/02
This structure struct _MULTIMAT  *multimat is defined in global_control.c
and the type is in materials.h                                                  
It holds all information about the layered section data
</pre>
*----------------------------------------------------------------------*/
extern struct _MULTIMAT  *multimat;


/*!----------------------------------------------------------------------
\brief write the data needed to restart this element                                       

<pre>                     m.gee 5/02              modified by    sh 02/03
This routine writes the data needed to restart the shell9 element. 
</pre>
\param  ELEMENT  *actele   (i) actual element
\param  INT       nhandle  (i) size of handles
\param  long int *handles  ( ) unique handle returned by the pss-system 
\param  INT       init     (i) flag for initializing arrays

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: shell9()   [s9_main.c]

*----------------------------------------------------------------------*/
void s9_write_restart(ELEMENT *actele, INT nhandle, long int *handles, INT init)
{
INT j,n,k,kl,ml,ierr;
INT size_j;
INT num_klay;
INT num_mlay;
INT actlay;

MULTIMAT      *actmultimat;                               /* material of actual material layer */
static ARRAY   elewares_a;  static DOUBLE **elewares;     /* array for elewa */
FILE *out;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s9_write_restart");
#endif
/*---------------------------------------------------------- init phase */
if (init==1)
{
   /*array for every layer elewares_a
     fdim = MAXLAY_SHELL9 (total number of layers)

     sdim = max number of possible entries in ipwa per GP (35 at the moment) 
            * max number of GP (3x3x2 = 18 at the moment) */

    /*  pl_mises or pl_dp*/         
    /*  DOUBLE       epstn;        equivalent strain                          */
    /*  INT          yip;          stress state: 1=elastic 2=plastic          */
    /*  DOUBLE       sig[6];       stresses                           [6]     */
    /*  DOUBLE       eps[6]        strains                            [6]     */
    /*  DOUBLE       qn[6]         backstress                         [6]     */
    /*      ( == 20 )                                                         */

    /*  pl_epc */      
    /*  INT          yip;          stress state: 1=elastic 2=plastic          */
    /*  DOUBLE       kappa_t;      equivalent strain (tension)                */
    /*  DOUBLE       kappa_c;      equivalent strain (compression)            */
    /*  DOUBLE       sig[6];       stresses                            [6]    */
    /*  DOUBLE       eps[6]        strains                             [6]    */
    /*      ( == 15 )                                                         */

    /*  pl_hoff */      
    /*  INT          yip;          stress state: 1=elastic 2=plastic          */
    /*  DOUBLE       dhard;        equivalent strain                          */
    /*  DOUBLE       sig[6];       stresses                             [6]   */
    /*  DOUBLE       eps[6]        strains                              [6]   */
    /*  DOUBLE       dkappa[6];                                         [6]   */
    /*  DOUBLE       gamma[6];                                          [6]   */
    /*  DOUBLE       rkappa[9];                                         [9]   */
    /*      ( == 35 )                                                         */

   elewares     = amdef("elewares"  ,&elewares_a,MAXLAY_SHELL9,35*18,"DA");    

   goto end;
}
/*-------------------------------------------------------- uninit phase */
else if (init==-1)
{
   amdel(&elewares_a);   

   goto end;
}
/*----------------------------------------------------------------------*/
num_klay = actele->e.s9->num_klay;
/*
   the shell element has to write
   - eas parameters if necessary
   - elewa for material nonlinearity
*/
/*----------------------------------------------------------------------*/
out = allfiles.out_pss;
/*----------------------------------------------------------------------*/
handles[0]=0;
if (actele->e.s9->nhyb)
{
   handles[0]+=4;
   if (handles[0]+1 > nhandle) dserror("Handle range too small for element");
   /*------------------------------------------------------- write alfa */
   pss_write_array(&(actele->e.s9->alfa),&(handles[1]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   /*------------------------------------------------ write Dtildinv */
   pss_write_array(&(actele->e.s9->Dtildinv),&(handles[2]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   /*------------------------------------------------------ write L */
   pss_write_array(&(actele->e.s9->L),&(handles[3]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
   /*------------------------------------------------------ write Lt */
/*   pss_write_array(&(actele->e.s9->Lt),&(handles[3]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");*/
   /*-------------------------------------------------- write Rtilde */
   pss_write_array(&(actele->e.s9->Rtilde),&(handles[4]),out,&ierr);
   if (ierr != 1) dserror("Error writing restart data");
}
/*----------------------------------------------------------------------*/
/*--------- now check for ELEWA --->  material nonlinearity ------------*/
/*check if there is a working array for this element*/
if (!actele->e.s9->elewa) goto end; /*no elewa to this element*/

/*if there is a working array, write all the data into the elewares_a */
actlay = 0;
for (kl=0; kl<num_klay; kl++) /*loop all kinematic layers*/
{
   num_mlay = actele->e.s9->kinlay[kl].num_mlay;
   for (ml=0; ml<num_mlay; ml++) /*loop all material layers*/
   {
        /*check if there is a ipwa for this layer */
        if (!actele->e.s9->elewa[actlay].ipwa) goto next; /*no ipwa to this layer*/
        
        actmultimat = &(multimat[actele->e.s9->kinlay[kl].mmatID[ml]-1]);
        
        /*number of gausspoints in one layer*/     
        size_j = actele->e.s9->nGP[0]*actele->e.s9->nGP[1]*actele->e.s9->nGP[2];

        n = 0;

        /* pl_mises or pl_dp */
        if (actmultimat->mattyp == m_pl_mises ||
            actmultimat->mattyp == m_pl_dp    )
        {
           for (k=0; k<size_j; k++)/*write for every gausspoint -> elewares*/
           {
             elewares[actlay][n] = actele->e.s9->elewa[actlay].ipwa[k].epstn;
             n++;
             elewares[actlay][n] = (double) actele->e.s9->elewa[actlay].ipwa[k].yip;
             n++;

             for (j=0; j<6; j++)
             {
               elewares[actlay][n] = actele->e.s9->elewa[actlay].ipwa[k].sig[j];
               n++;
               elewares[actlay][n] = actele->e.s9->elewa[actlay].ipwa[k].eps[j];
               n++;
               elewares[actlay][n] = actele->e.s9->elewa[actlay].ipwa[k].qn[j];
               n++;
             }
           }
        } /*end pl_mises of pl_dp*/

        /* pl_epc */
        else if (actmultimat->mattyp == m_pl_epc   )
        {
           for (k=0; k<size_j; k++)/*write for every gausspoint -> elewares*/
           {
             elewares[actlay][n] = (double) actele->e.s9->elewa[actlay].ipwa[k].yip;
             n++;
             elewares[actlay][n] = actele->e.s9->elewa[actlay].ipwa[k].kappa_t;
             n++;
             elewares[actlay][n] = actele->e.s9->elewa[actlay].ipwa[k].kappa_c;
             n++;

             for (j=0; j<6; j++)
             {
               elewares[actlay][n] = actele->e.s9->elewa[actlay].ipwa[k].sig[j];
               n++;
               elewares[actlay][n] = actele->e.s9->elewa[actlay].ipwa[k].eps[j];
               n++;
             }
           }
        } /*end pl_epc */

        /* pl_hoff */
        else if (actmultimat->mattyp ==  m_pl_hoff)
        {
           for (k=0; k<size_j; k++)/*write for every gausspoint -> elewares*/
           {
             elewares[actlay][n] = (double) actele->e.s9->elewa[actlay].ipwa[k].yip;
             n++;
             elewares[actlay][n] = actele->e.s9->elewa[actlay].ipwa[k].dhard;
             n++;

             for (j=0; j<6; j++)
             {
               elewares[actlay][n] = actele->e.s9->elewa[actlay].ipwa[k].sig[j];
               n++;
               elewares[actlay][n] = actele->e.s9->elewa[actlay].ipwa[k].eps[j];
               n++;
               elewares[actlay][n] = actele->e.s9->elewa[actlay].ipwa[k].dkappa[j];
               n++;
               elewares[actlay][n] = actele->e.s9->elewa[actlay].ipwa[k].gamma[j];
               n++;
             }
             for (j=0; j<9; j++)
             {
               elewares[actlay][n] = actele->e.s9->elewa[actlay].ipwa[k].rkappa[j];
               n++;
             }
           }
        } /*end pl_hoff*/

        next:
        actlay += 1;

   }/*end loop material layers*/
}/*end loop kinematic layers*/

handles[0]+=1;
if (handles[0]+1 > nhandle) dserror("Handle range too small for element");
/*-------------------------------------------------- write elewares */
pss_write_array(&elewares_a,&(handles[5]),out,&ierr);
if (ierr != 1) dserror("Error writing restart data");


/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_write_restart */

/*!----------------------------------------------------------------------
\brief read the data needed to restart this element                                       

<pre>                     m.gee 5/02              modified by    sh 02/03
This routine reads the data needed to restart the shell9 element. 
</pre>
\param  ELEMENT  *actele   (i) actual element
\param  INT       nhandle  (i) size of handles
\param  long int *handles  ( ) unique handle returned by the pss-system 
\param  INT       init     (i) flag for initializing arrays

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: shell9()   [s9_main.c]

*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                        m.gee 5/02    |
 | read the data needed to restart this element                         |
 | modified from shell8                                    sh 02/03     |
 *----------------------------------------------------------------------*/
void s9_read_restart(ELEMENT *actele, long int *handles, INT init)
{
INT j,n,k,kl,ml,ierr;
INT size_j;
INT num_klay;
INT num_mlay;
INT actlay;

MULTIMAT      *actmultimat;                               /* material of actual material layer */
static ARRAY   elewares_a;  static DOUBLE **elewares;     /* array for elewa */
INT dims[3];
FILE *in;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s9_read_restart");
#endif
/*---------------------------------------------------------- init phase */
if (init==1)
{
   /*for dimensions see s9_write_restart */
               
   elewares     = amdef("elewares"  ,&elewares_a,MAXLAY_SHELL9,35*18,"DA");    

   goto end;
}
/*-------------------------------------------------------- uninit phase */
else if (init==-1)
{
   amdel(&elewares_a);   

   goto end;
}
num_klay = actele->e.s9->num_klay;
/*----------------------------------------------------------------------*/
/*
   the shell element has to read
   - eas parameters if necessary
   - elewa for material nonlinearity
*/
/*----------------------------------------------------------------------*/
in = allfiles.in_pss;
if (!in) dserror("There is no restart input file open");
/*----------------------------------------------------------------------*/
if (actele->e.s9->nhyb)
{
   /*------------------------------------------------- check dimensions alfa*/
   pss_getdims_name_handle(actele->e.s9->alfa.name,&dims[0],&dims[1],&dims[2],&handles[1],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if (actele->e.s9->alfa.fdim != dims[0] ||
       actele->e.s9->alfa.sdim != dims[1])
       dserror("Mismatch in reading element restart data");
   /*-------------------------------------------------------- read alfa */
   pss_read_array_name_handle(actele->e.s9->alfa.name,&(actele->e.s9->alfa),&handles[1],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");  

   /*------------------------------------------------- check dimensions Dtildinv*/
   pss_getdims_name_handle(actele->e.s9->Dtildinv.name,&dims[0],&dims[1],&dims[2],&handles[2],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if (actele->e.s9->Dtildinv.fdim != dims[0] ||
       actele->e.s9->Dtildinv.sdim != dims[1])
       dserror("Mismatch in reading element restart data");
   /*---------------------------------------------------- read Dtildinv */
   pss_read_array_name_handle(actele->e.s9->Dtildinv.name,&(actele->e.s9->Dtildinv),&handles[2],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");  

   /*------------------------------------------------- check dimensions L*/
   pss_getdims_name_handle(actele->e.s9->L.name,&dims[0],&dims[1],&dims[2],&handles[3],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if (actele->e.s9->L.fdim != dims[0] ||
       actele->e.s9->L.sdim != dims[1])
       dserror("Mismatch in reading element restart data");
   /*---------------------------------------------------------- read L */
   pss_read_array_name_handle(actele->e.s9->L.name,&(actele->e.s9->L),&handles[3],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");  

   /*------------------------------------------------- check dimensions Lt*/
/*   pss_getdims_name_handle(actele->e.s9->Lt.name,&dims[0],&dims[1],&dims[2],&handles[3],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if (actele->e.s9->Lt.fdim != dims[0] ||
       actele->e.s9->Lt.sdim != dims[1])
       dserror("Mismatch in reading element restart data");*/
   /*---------------------------------------------------------- read Lt */
/*   pss_read_array_name_handle(actele->e.s9->Lt.name,&(actele->e.s9->Lt),&handles[3],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data"); */ 

   /*------------------------------------------------- check dimensions Rtilde*/
   pss_getdims_name_handle(actele->e.s9->Rtilde.name,&dims[0],&dims[1],&dims[2],&handles[4],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");
   if (actele->e.s9->Rtilde.fdim != dims[0] ||
       actele->e.s9->Rtilde.sdim != dims[1])
       dserror("Mismatch in reading element restart data");
   /*------------------------------------------------------ read Rtilde */
   pss_read_array_name_handle(actele->e.s9->Rtilde.name,&(actele->e.s9->Rtilde),&handles[4],in,&ierr);
   if (ierr != 1) dserror("Cannot read restart data");  
}
/*----------------------------------------------------------------------*/
/*--------- now check for ELEWA --->  material nonlinearity ------------*/
/*check if there is a working array for this element*/
if (!actele->e.s9->elewa) goto end; /*no elewa to this element*/

/*if there is a working array, read all the data into the elewares_a */
/*------------------------------------------------- check dimensions elewares_a*/
pss_getdims_name_handle(elewares_a.name,&dims[0],&dims[1],&dims[2],&handles[5],in,&ierr);
if (ierr != 1) dserror("Cannot read restart data");
if (elewares_a.fdim != dims[0] ||
    elewares_a.sdim != dims[1])
    dserror("Mismatch in reading element restart data");
/*-------------------------------------------------------- read elewares_a */
pss_read_array_name_handle(elewares_a.name,&elewares_a,&handles[5],in,&ierr);
if (ierr != 1) dserror("Cannot read restart data");  


/*----- now write the restart data back to elewa -------------------------*/
actlay = 0;
for (kl=0; kl<num_klay; kl++) /*loop all kinematic layers*/
{
   num_mlay = actele->e.s9->kinlay[kl].num_mlay;
   for (ml=0; ml<num_mlay; ml++) /*loop all material layers*/
   {
    /*check if there is a ipwa for this layer */
        if (!actele->e.s9->elewa[actlay].ipwa) goto next; /*no ipwa to this layer*/
        
        actmultimat = &(multimat[actele->e.s9->kinlay[kl].mmatID[ml]-1]);
        
        /*number of gausspoints in one layer*/     
        size_j = actele->e.s9->nGP[0]*actele->e.s9->nGP[1]*actele->e.s9->nGP[2];

        n = 0;

        /* pl_mises or pl_dp */
        if (actmultimat->mattyp == m_pl_mises ||
            actmultimat->mattyp == m_pl_dp    )
        {
           for (k=0; k<size_j; k++)/*read for every gausspoint*/
           {
             actele->e.s9->elewa[actlay].ipwa[k].epstn   = elewares[actlay][n];
             n++;
             actele->e.s9->elewa[actlay].ipwa[k].yip     = (int) elewares[actlay][n];
             n++;

             for (j=0; j<6; j++)
             {
               actele->e.s9->elewa[actlay].ipwa[k].sig[j]    = elewares[actlay][n];
               n++;
               actele->e.s9->elewa[actlay].ipwa[k].eps[j]    = elewares[actlay][n];
               n++;
               actele->e.s9->elewa[actlay].ipwa[k].qn[j]     = elewares[actlay][n];
               n++;
             }
           }
        } /*end pl_mises of pl_dp*/

        /* pl_epc */
        else if (actmultimat->mattyp == m_pl_epc   )
        {
           for (k=0; k<size_j; k++)/*read for every gausspoint*/
           {
             actele->e.s9->elewa[actlay].ipwa[k].yip     = (int) elewares[actlay][n];
             n++;
             actele->e.s9->elewa[actlay].ipwa[k].kappa_t = elewares[actlay][n];
             n++;
             actele->e.s9->elewa[actlay].ipwa[k].kappa_c = elewares[actlay][n];
             n++;

             for (j=0; j<6; j++)
             {
               actele->e.s9->elewa[actlay].ipwa[k].sig[j]    = elewares[actlay][n];
               n++;
               actele->e.s9->elewa[actlay].ipwa[k].eps[j]    = elewares[actlay][n];
               n++;
             }
           }
        } /*end pl_epc */

        /* pl_hoff */
        else if (actmultimat->mattyp ==  m_pl_hoff)
        {
           for (k=0; k<size_j; k++)/*read for every gausspoint*/
           {
             actele->e.s9->elewa[actlay].ipwa[k].yip     = (int) elewares[actlay][n];
             n++;
             actele->e.s9->elewa[actlay].ipwa[k].dhard   = elewares[actlay][n];
             n++;

             for (j=0; j<6; j++)
             {
               actele->e.s9->elewa[actlay].ipwa[k].sig[j]    = elewares[actlay][n];
               n++;
               actele->e.s9->elewa[actlay].ipwa[k].eps[j]    = elewares[actlay][n];
               n++;
               actele->e.s9->elewa[actlay].ipwa[k].dkappa[j] = elewares[actlay][n];
               n++;
               actele->e.s9->elewa[actlay].ipwa[k].gamma[j]  = elewares[actlay][n];
               n++;
             }
             for (j=0; j<9; j++)
             {
               actele->e.s9->elewa[actlay].ipwa[k].rkappa[j] = elewares[actlay][n];
               n++;
             }
           }
        } /*end pl_hoff*/

        next:
        actlay += 1;

   }/*end loop material layers*/
}/*end loop kinematic layers*/

/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_read_restart */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
