#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
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
 |  print out turbulence statistics spatially averaged over             |   
 |  homogeneous directions for fluid                        gravem 02/03|
 *----------------------------------------------------------------------*/
void out_turbstat(FIELD     *actfield, 
                  PARTITION *actpart, 
		  INTRA     *actintra,
		  DOUBLE   **statistic,
		  INT       *indicator,
		  INT        numdco,
		  INT        numnhd,
		  INT        numsnhd,
		  INT        statimst,
		  INT        runtimst)
{
INT        j,k,l;
FILE      *out = allfiles.out_out;
NODE      *actnode;
INT        myrank;
INT        nprocs;
INT        imyrank;
INT        inprocs;
INT	   numnp;
#ifdef DEBUG 
dstrc_enter("out_turbstat");
#endif

numnp = actfield->dis[0].numnp;
/*----------------------------------------------------------------------*/
myrank = par.myrank;
nprocs = par.nprocs;
imyrank= actintra->intra_rank;
inprocs= actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
if (imyrank==0 && myrank==0)
{
/*-------------------------------------------------------- print header */
/*----------------------------------------------------------------------*/
fprintf(out,"================================================================================\n");
fprintf(out,"TURBULENCE STATISTICS (SPATIALLY AVERAGED OVER HOMOGENEOUS DIRECTIONS)\n");
fprintf(out,"(started at timestep %6d and run for %6d timesteps)\n",statimst,runtimst);
fprintf(out,"================================================================================\n");
fprintf(out,"Mean Values of Velocity (x/y/z), Pressure and Vorticity\n"); 
fprintf(out,"================================================================================\n");
/*------------------------------------ print (homogenized) nodal values */
for (j=0; j<numdco; j++)
{
  for (k=0; k<numnp; k++)
  {
    if (indicator[k]==j) 
    { 
      actnode = &(actfield->dis[0].node[k]);
      goto identified1;
    }  
  }
  identified1:
  if (numsnhd!=0)
  {
    if (j<numsnhd) fprintf(out,"Y-COORDINATE %10.3E    ",actnode->x[1]);
    else           fprintf(out,"X-COORDINATE %10.3E    ",actnode->x[0]);
  }
  else fprintf(out,"COORDINATE %10.3E    ",actnode->x[numnhd]);
  fprintf(out,"%20.7E ",statistic[j][0]);
  fprintf(out,"%20.7E ",statistic[j][4]);
  fprintf(out,"%20.7E ",statistic[j][8]);
  fprintf(out,"%20.7E ",statistic[j][12]);
  fprintf(out,"%20.7E ",statistic[j][15]);
  fprintf(out,"\n");
}
fprintf(out,"________________________________________________________________________________\n\n");
fprintf(out,"\n");
fprintf(out,"\n");

fprintf(out,"================================================================================\n");
fprintf(out,"Stresses xy(12): Viscous Shear Stresses, Reynolds Stresses, Total Stresses\n"); 
fprintf(out,"================================================================================\n");
/*------------------------------------ print (homogenized) nodal values */
for (j=0; j<numdco; j++)
{
  for (k=0; k<numnp; k++)
  {
    if (indicator[k]==j) 
    { 
      actnode = &(actfield->dis[0].node[k]);
      goto identified2;
    }  
  }
  identified2:
  if (numsnhd!=0)
  {
    if (j<numsnhd) fprintf(out,"Y-COORDINATE %10.3E    ",actnode->x[1]);
    else           fprintf(out,"X-COORDINATE %10.3E    ",actnode->x[0]);
  }
  else fprintf(out,"COORDINATE %10.3E    ",actnode->x[numnhd]);
  fprintf(out,"%20.7E ",statistic[j][17]);
  fprintf(out,"%20.7E ",statistic[j][18]);
  fprintf(out,"%20.7E ",statistic[j][19]);
  fprintf(out,"\n");
}
fprintf(out,"________________________________________________________________________________\n\n");
fprintf(out,"\n");
fprintf(out,"\n");

fprintf(out,"================================================================================\n");
fprintf(out,"Root-mean-square Velocities (x/y/z), Pressure and Vorticity\n"); 
fprintf(out,"================================================================================\n");
/*------------------------------------ print (homogenized) nodal values */
for (j=0; j<numdco; j++)
{
  for (k=0; k<numnp; k++)
  {
    if (indicator[k]==j) 
    { 
      actnode = &(actfield->dis[0].node[k]);
      goto identified3;
    }  
  }
  identified3:
  if (numsnhd!=0)
  {
    if (j<numsnhd) fprintf(out,"Y-COORDINATE %10.3E    ",actnode->x[1]);
    else           fprintf(out,"X-COORDINATE %10.3E    ",actnode->x[0]);
  }
  else fprintf(out,"COORDINATE %10.3E    ",actnode->x[numnhd]);
  fprintf(out,"%20.7E ",statistic[j][20]);
  fprintf(out,"%20.7E ",statistic[j][21]);
  fprintf(out,"%20.7E ",statistic[j][22]);
  fprintf(out,"%20.7E ",statistic[j][23]);
  fprintf(out,"%20.7E ",statistic[j][24]);
  fprintf(out,"\n");
}
fprintf(out,"________________________________________________________________________________\n\n");
fprintf(out,"\n");
fprintf(out,"\n");

fprintf(out,"================================================================================\n");
fprintf(out,"Skewness of Velocities (x/y/z)\n"); 
fprintf(out,"================================================================================\n");
/*------------------------------------ print (homogenized) nodal values */
for (j=0; j<numdco; j++)
{
  for (k=0; k<numnp; k++)
  {
    if (indicator[k]==j) 
    { 
      actnode = &(actfield->dis[0].node[k]);
      goto identified4;
    }  
  }
  identified4:
  if (numsnhd!=0)
  {
    if (j<numsnhd) fprintf(out,"Y-COORDINATE %10.3E    ",actnode->x[1]);
    else           fprintf(out,"X-COORDINATE %10.3E    ",actnode->x[0]);
  }
  else fprintf(out,"COORDINATE %10.3E    ",actnode->x[numnhd]);
  fprintf(out,"%20.7E ",statistic[j][25]);
  fprintf(out,"%20.7E ",statistic[j][26]);
  fprintf(out,"%20.7E ",statistic[j][27]);
  fprintf(out,"\n");
}
fprintf(out,"________________________________________________________________________________\n\n");
fprintf(out,"\n");
fprintf(out,"\n");

fprintf(out,"================================================================================\n");
fprintf(out,"Flatness of Velocities (x/y/z)\n"); 
fprintf(out,"================================================================================\n");
/*------------------------------------ print (homogenized) nodal values */
for (j=0; j<numdco; j++)
{
  for (k=0; k<numnp; k++)
  {
    if (indicator[k]==j) 
    { 
      actnode = &(actfield->dis[0].node[k]);
      goto identified5;
    }  
  }
  identified5:
  if (numsnhd!=0)
  {
    if (j<numsnhd) fprintf(out,"Y-COORDINATE %10.3E    ",actnode->x[1]);
    else           fprintf(out,"X-COORDINATE %10.3E    ",actnode->x[0]);
  }
  else fprintf(out,"COORDINATE %10.3E    ",actnode->x[numnhd]);
  fprintf(out,"%20.7E ",statistic[j][28]);
  fprintf(out,"%20.7E ",statistic[j][29]);
  fprintf(out,"%20.7E ",statistic[j][30]);
  fprintf(out,"\n");
}
fprintf(out,"________________________________________________________________________________\n\n");
fprintf(out,"\n");
fprintf(out,"\n");

/*----------------------------------------------------------------------*/
} /* end of if (myrank==0 && imyrank==0) */
/*----------------------------------------------------------------------*/
if (myrank==0 && imyrank==0) fflush(out);
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_turbstat */

/*----------------------------------------------------------------------*
 |  print out wall shear stress and flowrate with the respective        |   
 |  temporally averaged values for turbulent channel flow   gravem 04/03|
 *----------------------------------------------------------------------*/
void out_washstr(FLUID_DYNAMIC *fdyn,
                 FLUID_DYN_CALC *dynvar,
                 FIELD     *actfield, 
                 PARTITION *actpart, 
		 INTRA     *actintra)
{
FILE      *out = allfiles.out_out;
INT        myrank;
INT        imyrank;
DOUBLE     flowrate,flowrate1;
DOUBLE     avwastr,avflrate;
#ifdef DEBUG 
dstrc_enter("out_washstr");
#endif

/*----------------------------------------------------------------------*/
dynvar->avwastr+=dynvar->washstr;
avwastr=dynvar->avwastr/fdyn->step;
flowrate1=dynvar->flowrate;
#ifdef PARALLEL 
MPI_Allreduce(&flowrate1,&flowrate,1,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
flowrate=flowrate/2.8574;
#else
flowrate=flowrate1/2.8574;
#endif
dynvar->flowrate=flowrate;
dynvar->avflrate+=flowrate;
avflrate=dynvar->avflrate/fdyn->step;

myrank = par.myrank;
imyrank= actintra->intra_rank;
/*----------------------------------------------------------------------*/
if (imyrank==0 && myrank==0)
{
  fprintf(out,"wall shear stress in step %5d at time %6.2f: %10.9f \n",
          fdyn->step,fdyn->time,dynvar->washstr);
  fprintf(out,"flow rate in step %5d at time %6.2f: %10.9f \n",
          fdyn->step,fdyn->time,flowrate);
  fprintf(out,"averaged wall shear stress in step %5d at time %6.2f: %10.9f \n",
          fdyn->step,fdyn->time,avwastr);
  fprintf(out,"averaged flow rate in step %5d at time %6.2f: %10.9f \n",
          fdyn->step,fdyn->time,avflrate);
  fprintf(out,"\n");
} /* end of if (myrank==0 && imyrank==0) */
/*----------------------------------------------------------------------*/
if (myrank==0 && imyrank==0) fflush(out);

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_washstr */

/*----------------------------------------------------------------------*
 |  print out energy (and vorticity thickness as well as enstrophy)     |   
 |  for turb. 3-D driv. cav. (and 2-D mixing layer)         gravem 04/03|
 *----------------------------------------------------------------------*/
void out_energy(FLUID_DYNAMIC *fdyn,
                FLUID_DYN_CALC *dynvar,
                FIELD	  *actfield, 
                PARTITION *actpart, 
		INTRA	  *actintra)
{
FILE      *out = allfiles.out_out;
INT        myrank;
INT        imyrank;
DOUBLE     energy,energy1;
DOUBLE     enstrophy,enstrophy1;
#ifdef DEBUG 
dstrc_enter("out_energy");
#endif

/*----------------------------------------------------------------------*/
energy1=dynvar->energy;
enstrophy1=dynvar->enstrophy;
#ifdef PARALLEL 
MPI_Allreduce(&energy1,&energy,1,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
MPI_Allreduce(&enstrophy1,&enstrophy,1,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#else
energy=energy1;
enstrophy=enstrophy1;
#endif
dynvar->energy=energy;
dynvar->enstrophy=enstrophy;

myrank = par.myrank;
imyrank= actintra->intra_rank;
/*----------------------------------------------------------------------*/
if (imyrank==0 && myrank==0)
{
  if (fdyn->turbstat==1) 
  {
    if (fdyn->step==0)
    fprintf(out,"step / time / relative vorticity thickness / energy / enstrophy \n");
    fprintf(out,"%5d  %8.5f  %10.9f  %10.9f  %10.9f \n",fdyn->step,fdyn->time,
            dynvar->vorthick,dynvar->energy,dynvar->enstrophy);
  }
  else if (fdyn->turbstat==4) 
  {
    if (fdyn->step==0) fprintf(out,"step / time / energy \n");
    fprintf(out,"%5d  %8.5f  %10.9f \n",fdyn->step,fdyn->time,dynvar->energy);
  }  	      
} /* end of if (myrank==0 && imyrank==0) */
/*----------------------------------------------------------------------*/
if (myrank==0 && imyrank==0) fflush(out);

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_energy */


 
