/*!----------------------------------------------------------------------
\file
\brief reading fluid data from fluid_start_data

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_FLUID
#include "../headers/standardtypes.h"
/*---------------------------------------------------------- prototypes */
void inp_fluid_frfind(char string[]);
/*---------------------------------------------------- global variables */
FILE *in;
char line [500];
INT eof;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
/*----------------------------------------------------------------------*
 | variables needed for parallel comp.                    m.gee 8/00    |
 *----------------------------------------------------------------------*/
extern struct _PAR   par;
/*!---------------------------------------------------------------------
\brief input from file 'fluid_start-data'

<pre>                                                         genk 07/02

in this routine the inital fluid data are read form 'fluid_start_data'
and stored in the array 'start'

</pre>
\return void
------------------------------------------------------------------------*/
void inp_fluid_start_data( FIELD   *actfield,
                           ARRAY_POSITION    *ipos,
                           FLUID_DYNAMIC *fdyn
			 )
{
INT     irc=1;                    /* flag for file opening              */
INT     i;                        /* counters                           */
INT     numnp=0; 	          /* number of nodes                    */
INT     numdf=0;                  /* number dofs per node               */
INT     num;                      /* actual number of node during input */
INT     stepin=-1;
INT     mone=-1;
INT     datastep = 0;
INT     step     = 0;
INT     foundstep;
DOUBLE  dens;
DOUBLE  time     = 0.0;
ARRAY   pre_a;
DOUBLE *pre;
ARRAY   velx_a;
DOUBLE *velx;
ARRAY   vely_a;
DOUBLE *vely;
ARRAY   velz_a;
DOUBLE *velz     = NULL;
ARRAY   globloc;                   /* global - local node Ids            */
NODE  *actnode;
ELEMENT *actele;
char   *foundit = NULL;
char   *end;

#ifdef DEBUG
dstrc_enter("inp_fluid_start_data");
#endif

/* ------------------------------------------ open file fluid_start.data */
if (par.myrank==0)
{
   in = fopen("fluid_start.data","r");
   if (in == NULL) irc=0;
}
#ifdef PARALLEL
MPI_Bcast(&irc,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
if (irc==0)
{
   if (par.myrank==0)
   printf("opening of file fluid_start.data failed\n");
#ifdef PARALLEL
   MPI_Finalize();
#else
   exit(1);
#endif
}
/*------------------------------------------------- store / check values */
numnp = actfield->dis[0].numnp;
numdf = fdyn->numdf;
/*------------------------------------- allocate vector for storing data */
velx = amdef("velx",&velx_a,numnp,1,"DV");
vely = amdef("vely",&vely_a,numnp,1,"DV");
pre  = amdef("pre" ,&pre_a ,numnp,1,"DV");
if (numdf==4)
velz = amdef("velz",&velz_a,numnp,1,"DV");
/*------------------------------------- determine global local node Ids */
amdef("globloc",&globloc,genprob.nnode,1,"IV");
aminit(&globloc,&mone);
for (i=0;i<numnp;i++)
{
   actnode = &(actfield->dis[0].node[i]);
   globloc.a.iv[actnode->Id] = actnode->Id_loc;
}

/*------------------------------------------ read initial data from file */
if (par.myrank==0)
{
  if (fdyn->resstep==-1)
  {
     while (eof==0)
     {
        inp_fluid_frfind("# RESULT velocity on FIELD fluid");
        stepin++;
        if (fgets(line,499,in)==NULL && eof==0)
           dserror("An error occured reading a line from fluid_start.data");
     }
     datastep=stepin-1;
  }
  else if (fdyn->resstep>0)
  {
     foundstep=0;
     while (eof==0 && foundstep==0)
     {
        inp_fluid_frfind("# RESULT velocity on FIELD fluid");
	stepin++;
        inp_fluid_frfind("# STEP");
        foundit=strpbrk(line,"-.1234567890");
        step = strtod(foundit,&end);
        if (step==fdyn->resstep) foundstep++;
        if (fgets(line,499,in)==NULL && eof==0)
           dserror("An error occured reading a line from fluid_start.data");
     }
     if (foundstep==1)
        datastep=stepin;
     else
        dserror("Restart step not in file fluid_start.data\n");
  }
  else
     dserror("An error occured reading a line from fluid_start.data");
  rewind(in);
  if (fgets(line,499,in)==NULL)
     dserror("An error occured reading a line from fluid_start.data");
  stepin=-1;
  eof=0;
  while (stepin<datastep)
  {
     inp_fluid_frfind("# RESULT velocity on FIELD fluid");
     stepin++;
     if (eof==1)
        dserror("Cannot read from fluid_start.data: stepin non-existent!\n");
     if (fgets(line,499,in)==NULL)
        dserror("An error occured reading a line from fluid_start.data");
  }
  /*--------------------------------------------- find & read time & step*/
  inp_fluid_frfind("# TIME");
  foundit=strpbrk(line,"-.1234567890");
  time = strtod(foundit,&end);
  inp_fluid_frfind("# STEP");
  foundit=strpbrk(line,"-.1234567890");
  step = strtod(foundit,&end);
   /*-------------------------------------- find & read velocity results */
  inp_fluid_frfind("VALUES");
  for (i=0;i<numnp;i++)
  {
     if (fgets(line,499,in)==NULL)
  	dserror("An error occured reading a line from fluid_start.data");
     foundit = strstr(line," ");
     foundit=strpbrk(foundit,"-.1234567890");
     /*-------------------------------------------- read global node Id */
     num = strtod(foundit,&end)-1;
     /*---------------------------------------- determine local node Id */
     num = globloc.a.iv[num];
     if (num<0) dserror("node number not valid!\n");
     foundit = strstr(foundit," ");
     foundit=strpbrk(foundit,"-.1234567890");
     velx[num] = strtod(foundit,&end);
     foundit = strstr(foundit," ");
     foundit=strpbrk(foundit,"-.1234567890");
     vely[num] = strtod(foundit,&end);
     if (numdf==4)
     {
        foundit = strstr(foundit," ");
        foundit=strpbrk(foundit,"-.1234567890");
        velz[num] = strtod(foundit,&end);
     }
  }
  /*---------------------------------------------- plausibility checkes */
  if (fgets(line,499,in)==NULL)
     dserror("An error occured reading a line from fluid_start.data\n");
  if (strstr(line,"END VALUES") == NULL)
     dserror("Number of Fluid nodes not correct in fluid_start.datan");
  /*------------------------------------- find  & read pressure results */
  inp_fluid_frfind("# RESULT pressure on FIELD fluid");
  inp_fluid_frfind("VALUES");
  for (i=0;i<numnp;i++)
  {
     if (fgets(line,499,in)==NULL)
  	dserror("An error occured reading a line from fluid_start.data");
     foundit = strstr(line," ");
     foundit=strpbrk(foundit,"-.1234567890");
     /*-------------------------------------------- read global node Id */
     num = strtod(foundit,&end)-1;
     /*---------------------------------------- determine local node Id */
     num = globloc.a.iv[num];
     if (num<0) dserror("node number not valid!\n");
     foundit = strstr(foundit," ");
     foundit=strpbrk(foundit,"-.1234567890");
     pre[num] = strtod(foundit,&end);
  }
  /*---------------------------------------------- plausibility checkes */
  if (fgets(line,499,in)==NULL)
     dserror("An error occured reading a line from fluid_start.data\n");
  if (strstr(line,"END VALUES") == NULL)
     dserror("Number of Fluid nodes not correct in fluid_start.data\n");

/*----------------------------------------- close file fluid_start.data */
   fclose(in);
   printf("initial field read from    fluid_start.data\n\n");
}

/*----------------------------------------------- broadcast array start */
#ifdef PARALLEL
MPI_Bcast(velx,numnp,MPI_DOUBLE,0,MPI_COMM_WORLD);
MPI_Bcast(vely,numnp,MPI_DOUBLE,0,MPI_COMM_WORLD);
if (numdf==4)
MPI_Bcast(velz,numnp,MPI_DOUBLE,0,MPI_COMM_WORLD);
MPI_Bcast(pre ,numnp,MPI_DOUBLE,0,MPI_COMM_WORLD);
MPI_Bcast(&time,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
MPI_Bcast(&step,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
/*------------------------------------------ copy values to the nodes  */
fdyn->acttime=time;
fdyn->step=step;
actele = &(actfield->dis[0].element[0]);
dens  = mat[actele->mat-1].m.fluid->density;
for (i=0;i<numnp;i++)
{
   actnode=&(actfield->dis[0].node[i]);
   num = actnode->Id_loc;
   actnode->sol.a.da[0][0] = velx[num];
   actnode->sol_increment.a.da[ipos->velnm][0] = velx[num];
   actnode->sol_increment.a.da[ipos->veln][0] = velx[num];
   actnode->sol_increment.a.da[ipos->velnp][0] = velx[num];
   actnode->sol.a.da[0][1] = vely[num];
   actnode->sol_increment.a.da[ipos->velnm][1] = vely[num];
   actnode->sol_increment.a.da[ipos->veln][1] = vely[num];
   actnode->sol_increment.a.da[ipos->velnp][1] = vely[num];
   if(numdf==4)
   {
      actnode->sol.a.da[0][2] = velz[num];
      actnode->sol_increment.a.da[ipos->velnm][2] = velz[num];
      actnode->sol_increment.a.da[ipos->veln][2] = velz[num];
      actnode->sol_increment.a.da[ipos->velnp][2] = velz[num];
      actnode->sol.a.da[0][3] = pre[num]/dens;
      actnode->sol_increment.a.da[ipos->velnm][3] = pre[num]/dens;
      actnode->sol_increment.a.da[ipos->veln][3] = pre[num]/dens;
      actnode->sol_increment.a.da[ipos->velnp][3] = pre[num]/dens;
   }
   else
   {
      actnode->sol.a.da[0][2] = pre[num];
      actnode->sol_increment.a.da[ipos->velnm][2] = pre[num]/dens;
      actnode->sol_increment.a.da[ipos->veln][2] = pre[num]/dens;
      actnode->sol_increment.a.da[ipos->velnp][2] = pre[num]/dens;
   }
}
/*--------------------------------------------- delete temporary arrays */
amdel(&velx_a);
amdel(&vely_a);
if (numdf==4) amdel(&velz_a);
amdel(&pre_a);
amdel(&globloc);
/*----------------------------------------------------------------------*/

#ifdef DEBUG
dstrc_exit();
#endif
return;
}
/* end of inp_fluid_start_data */

/*!---------------------------------------------------------------------
\brief find a character string in fluid_start.data

<pre>                                                        genk 07/03
searches for a given character string in fluid_start.data
</pre>
\param string   char[]   (i)   string to search for in input
\return void

------------------------------------------------------------------------*/
void inp_fluid_frfind(char string[])
{

#ifdef DEBUG
dstrc_enter("inp_fluid_frfind");
#endif

while ( strstr(line,string) == NULL )
{
   if (fgets(line,499,in)==NULL)
   {
      if (feof(in)!=0)
      {
         eof = 1;
         goto end;
      }
      else
      dserror("An error occured reading a line from fluid_start.data");
   }
}
end:

#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inp_fluid_frfind */
#endif
#endif
