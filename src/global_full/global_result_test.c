/*!----------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

*----------------------------------------------------------------------*/
#ifdef RESULTTEST
#include "../headers/standardtypes.h"
#include "../axishell/axishell.h"
#include "../shell9/shell9.h"
#include "../fluid_full/fluid_prototypes.h"
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
#ifdef D_FLUID
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
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: velx is NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: velx not correct!");
      /*---------------------------------------------------- check vely */
      actresult   = actnode->sol.a.da[0][1];
      givenresult = 0.0001027305751481696;
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: vely is NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: vely not correct!");
      /*---------------------------------------------------- check pres */
      actresult   = actnode->sol.a.da[0][2];
      givenresult = 0.001862187203092821;
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: pressure is NAN!");
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
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: velx is NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: velx not correct!");
      /*---------------------------------------------------- check vely */
      actresult   = actnode->sol.a.da[0][1];
      givenresult = -0.00031073969819614;
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: vely is NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: vely not correct!");
      /*---------------------------------------------------- check pres */
      actresult   = actnode->sol.a.da[0][2];
      givenresult = 489.9264285641417;
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: pressure is NAN!");
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
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: velx is NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: velx not correct!");
      /*---------------------------------------------------- check vely */
      actresult   = actnode->sol.a.da[0][1];
      givenresult = 0.0007473681392594219;
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: vely is NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: vely not correct!"); 
      /*---------------------------------------------------- check pres */
      actresult   = actnode->sol.a.da[0][2];
      givenresult = 0.02664107119660589;
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: pressure is NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: pressure not correct!");
   }
   else if(strstr(allfiles.title[0],"BELTRAMI-FLOW") != NULL)
   {
      printf("\nCalculating Errors for Beltrami-Flow ...\n");

      fluid_cal_error(fluidfield,1);
   }
   else if(strstr(allfiles.title[0],"KIM-MOIN-FLOW") != NULL)
   {
      printf("\nCalculating Errors for Kim-Moin-Flow ...\n");

      fluid_cal_error(fluidfield,2);
   }
break;
#endif
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
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: vely is NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: vely not correct!");
      /*---------------------------------------------------- check pres */
      actresult   = actnode->sol.a.da[0][2];
      givenresult = 244.996442716119;
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: pressure is NAN!");
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
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: vely is NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: vely not correct!");
      /*---------------------------------------------------- check pres */
      actresult   = actnode->sol.a.da[0][2];
      givenresult = 0.0002240443391469128;
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: pressure is NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: pressure not correct!");
      /*-- check result of struct node with global Id 1781 at T = 0.5 s */
      actnode = &(structfield->dis[0].node[16]);
      /*------------------------------------------------------ check dy */
      actresult   = actnode->sol.a.da[0][1];
      givenresult = -4.573526776331859e-05;
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: dy is NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: dy not correct!");	 
   }
break; 
case prb_structure:
#ifdef D_AXISHELL
   /*-------------------------------- check results of  axishell.dat */
   if(strstr(allfiles.title[0],"AXISHELL_TEST_EXAMPLE") != NULL)
   {
      /*-- check result of fluid node with global Id 123 */
      printf("\nChecking results for AXISHELL_TEST_EXAMPLE ...\n");
      /*-------------------------------------------- check displacement */
      actresult   = structfield->dis[0].node[123].sol.a.da[0][0];
      givenresult = -0.02778858368483676;
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: vely is NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: vely not correct!");
      /*-------------------------------------------------- check stress */  
      /*actresult   = structfield->dis[0].element[123].e.saxi.stress_GP;*/
      actresult   = structfield->dis[0].element[123].e.saxi->stress_GP.a.d3[0][2][0];
      givenresult = -169.4900467515912;
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: pressure is NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: pressure not correct!");
   }
#endif
#ifdef D_SHELL9
   /*-------------------------- check results of  shell9_el_lay_eans.dat */
   if(strstr(allfiles.title[0],"SHELL9_7P_el_lay") != NULL)
   {
      /*-- check result of node with global Id 6 */
      printf("\nChecking results for SHELL9_7P_el_lay ...\n");
      /*-------------------------------------------- check displacement */
      actresult   = structfield->dis[0].node[6].sol.a.da[0][1];
      givenresult = 0.004085462018049922;
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: displacement is NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: displacement not correct!");
      /*-------------------------------------------------- check stresses */
      actresult   = structfield->dis[0].element[3].e.s9->stresses.a.d3[0][2][0];
      givenresult = 1.323600311248807;
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: stresses are NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: stresses not correct!");
   }
   /*-------------------------- check results of  shell9_kreuz_easnl.dat */
   if(strstr(allfiles.title[0],"SHELL9_ortho_geoNL") != NULL)
   {
      /*-- check result of node with global Id 48 after 3rd step*/
      printf("\nChecking results for SHELL9_ortho_geoNL ...\n");
      /*-------------------------------------------- check displacement */
      actresult   = structfield->dis[0].node[48].sol.a.da[0][1];
      givenresult = -0.01117777927876936;
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: displacement is NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: displacement not correct!");
      /*-------------------------------------------------- check stresses */
      actresult   = structfield->dis[0].element[8].e.s9->stresses.a.d3[0][0][3];
      givenresult = 95.62657876539751;
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: stresses are NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: stresses not correct!");

      actresult = sqrt(-actresult);


      
   }
   /*-------------------------- check results of  shell9_matnl.dat */
   if(strstr(allfiles.title[0],"SHELL9_scordelis_matNL") != NULL)
   {
      /*-- check result of node with global Id 11 after 3rd step*/
      printf("\nChecking results for SHELL9_scordelis_matNL ...\n");
      /*-------------------------------------------- check displacement */
      actresult   = structfield->dis[0].node[11].sol.a.da[0][2];
      givenresult = 0.006214781874804086;
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: displacement is NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: displacement not correct!");
      /*-------------------------------------------------- check stresses */
      actresult   = structfield->dis[0].element[11].e.s9->stresses.a.d3[0][0][1];
      givenresult = 0.3065420533487162;
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: stresses are NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: stresses not correct!");
   }
#endif /*D_SHELL9*/
break;
case prb_opt:
break;
case prb_ale:
#ifdef D_ALE
   /*-------------------------------- check results of  ale_3d.dat */
   if(strstr(allfiles.title[0],"PURE_ALE_TEST_EXAMPLE") != NULL)
   {
      /*-- check result of fluid node with global Id 123 */
      printf("\nChecking results for PURE_ALE_TEST_EXAMPLE ...\n");
      /*-------------------------------------------- check displacement */
      actresult   = alefield->dis[0].node[123].sol.a.da[0][1];
      givenresult = -0.5152044455551235;
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (!(FABS(actresult-givenresult)==FABS(actresult-givenresult)))
         dserror("RESULTCHECK: vely is NAN!");
      if (FABS(actresult-givenresult)>EPS6)
         dserror("RESULTCHECK: vely not correct!");
   }
   /*-------------------------------- check results of ale2_nln.dat */
   if(strstr(allfiles.title[0],"mOVING hOLE, nonlinear ale2-test") != NULL)
   {
      /*-- check result of ale node with global Id 123 */
      printf("\nChecking results for nonlinear ale2-test example ...\n");
      /*-------------------------------------------- check displacement */
      actresult   = alefield->dis[0].node[123].sol.a.da[0][1];
      givenresult = 1.6274434734043111;
      fprintf(err,"result test\n");
      fprintf(err,"actual = %24.16f, given = %24.16f\n",actresult,givenresult);
      if (FABS(actresult/givenresult-ONE)>EPS6)
         dserror("RESULTCHECK: disp-y not correct!");
   }
#endif
break;
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of global_result_test */
#endif
