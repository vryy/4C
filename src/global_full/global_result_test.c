/*!----------------------------------------------------------------------
\file
\brief Testing of results

------------------------------------------------------------------------*/
#ifdef RESULTTEST
#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
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
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;                      
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*!---------------------------------------------------------------------                                         
\brief testing of result 

<pre>                                                         genk 10/03      

Before checking in the latest version it's necessery to check the whole
program. In this context it seems to be useful to check the numerical
results, too.
Based on the problem titel one can add verified results to this function.

</pre>  
\return void                                                                       
\warning Titel must not be longer than one line!

------------------------------------------------------------------------*/
void global_result_test() 
{
FIELD  *alefield,*structfield,*fluidfield;
NODE   *actnode;
DOUBLE  actresult, givenresult;
FILE   *err = allfiles.out_err;

#ifdef DEBUG 
dstrc_enter("global_result_test");
#endif

if (genprob.numff>-1) fluidfield  = &(field[genprob.numff]);
if (genprob.numsf>-1) structfield = &(field[genprob.numsf]);
if (genprob.numaf>-1) alefield    = &(field[genprob.numaf]);
   
switch (genprob.probtyp)
{
case prb_fluid:
   /*------------------------ check results of f2_drivencavity20x20.dat */
   if(strstr(allfiles.title[0],"dRIVEN_cAVITY_iNPUT_tESTING") != NULL)
   {
      /* check result of node with global Id 172 (center) at T = 1.00 s */
      printf("\nChecking results ...\n");
      actnode = &(fluidfield->dis[0].node[172]);
      /*---------------------------------------------------- check velx */
      actresult   = actnode->sol.a.da[0][0];
      givenresult = -0.01144271046249421;
      fprintf(err,"actual = %24.16lf, given = %24.16lf\n",actresult,givenresult);
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: velx not correct!");
      /*---------------------------------------------------- check vely */
      actresult   = actnode->sol.a.da[0][1];
      givenresult = 0.0001027305751481696;
      fprintf(err,"actual = %24.16lf, given = %24.16lf\n",actresult,givenresult);
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: vely not correct!");
      /*---------------------------------------------------- check pres */
      actresult   = actnode->sol.a.da[0][2];
      givenresult = 0.001862187203092821;
      fprintf(err,"actual = %24.16lf, given = %24.16lf\n",actresult,givenresult);
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: pressure not correct!");
   }
   /*----------------------------- check results of f2_freeosz20x20.dat */
   else if(strstr(allfiles.title[0],"freeosz_iNPUT_tESTING") != NULL)
   {
      /* check result of node with global Id 172 (center) at T = 0.01 s */
      printf("\nChecking results ...\n");
      actnode = &(fluidfield->dis[0].node[172]);
      /*---------------------------------------------------- check velx */
      actresult   = actnode->sol.a.da[0][0];
      givenresult = 0.06307272483594084;
      fprintf(err,"actual = %24.16lf, given = %24.16lf\n",actresult,givenresult);
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: velx not correct!");
      /*---------------------------------------------------- check vely */
      actresult   = actnode->sol.a.da[0][1];
      givenresult = -0.00031073969819614;
      fprintf(err,"actual = %24.16lf, given = %24.16lf\n",actresult,givenresult);
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: vely not correct!");
      /*---------------------------------------------------- check pres */
      actresult   = actnode->sol.a.da[0][2];
      givenresult = 489.9264285641417;
      fprintf(err,"actual = %24.16lf, given = %24.16lf\n",actresult,givenresult);
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: pressure not correct!");
   }
   /*--------------------------- check results of f2pro_channelflow.dat */
   else if(strstr(allfiles.title[0],"f2pro_channelflow_iNPUT_tESTING") != NULL)
   {
      /*---------- check result of node with global Id 379 at T = 0.5 s */
      printf("\nChecking results ...\n");
      actnode = &(fluidfield->dis[0].node[379]);
      /*---------------------------------------------------- check velx */
      actresult   = actnode->sol.a.da[0][0];
      givenresult = 0.09770648919466481;
      fprintf(err,"actual = %24.16lf, given = %24.16lf\n",actresult,givenresult);
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: velx not correct!");
      /*---------------------------------------------------- check vely */
      actresult   = actnode->sol.a.da[0][1];
      givenresult = 0.0007473681392594219;
      fprintf(err,"actual = %24.16lf, given = %24.16lf\n",actresult,givenresult);
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: vely not correct!"); 
      /*---------------------------------------------------- check pres */
      actresult   = actnode->sol.a.da[0][2];
      givenresult = 0.02664107119660589;
      fprintf(err,"actual = %24.16lf, given = %24.16lf\n",actresult,givenresult);
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: pressure not correct!");
   }
break;
case prb_fsi:
   /*------------------------------- check results of fsi_tank20x10.dat */
   if(strstr(allfiles.title[0],"fsi_swfsload_nieder_iNPUT_tESTING") != NULL)
   {
      /*--- check result of fluid node with global Id 237 at T = 0.25 s */
      printf("\nChecking results ...\n");
      actnode = &(fluidfield->dis[0].node[111]);
      /*---------------------------------------------------- check vely */
      actresult   = actnode->sol.a.da[0][1];
      givenresult = -0.0002498596571028341;
      fprintf(err,"actual = %24.16lf, given = %24.16lf\n",actresult,givenresult);
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: vely not correct!");
      /*---------------------------------------------------- check pres */
      actresult   = actnode->sol.a.da[0][2];
      givenresult = 244.996442716119;
      fprintf(err,"actual = %24.16lf, given = %24.16lf\n",actresult,givenresult);
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: pressure not correct!");
   }
   /*-------------------------------- check results of  fsi_ow32x32.dat */
   if(strstr(allfiles.title[0],"Ueberstroemter_Hohlraum_iNPUT_tESTING") != NULL)
   {
      /*-- check result of fluid node with global Id 2154 at T = 0.25 s */
      printf("\nChecking results ...\n");
      actnode = &(fluidfield->dis[0].node[1052]);
      /*---------------------------------------------------- check vely */
      actresult   = actnode->sol.a.da[0][1];
      givenresult = -0.0001088850381827999;
      fprintf(err,"actual = %24.16lf, given = %24.16lf\n",actresult,givenresult);
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: vely not correct!");
      /*---------------------------------------------------- check pres */
      actresult   = actnode->sol.a.da[0][2];
      givenresult = 0.0002240443391469128;
      fprintf(err,"actual = %24.16lf, given = %24.16lf\n",actresult,givenresult);
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: pressure not correct!");
      /*-- check result of struct node with global Id 1781 at T = 0.5 s */
      actnode = &(structfield->dis[0].node[16]);
      /*------------------------------------------------------ check dy */
      actresult   = actnode->sol.a.da[0][1];
      givenresult = -4.573526776331859e-05;
      fprintf(err,"actual = %24.16lf, given = %24.16lf\n",actresult,givenresult);
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: dy not correct!");	 
   }
break; 
case prb_structure:
break;
case prb_opt:
break;
case prb_ale:
break;
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of V2GRID */
#endif
