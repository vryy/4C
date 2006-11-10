/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design                                           |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;

/*----------------------------------------------------------------------*
 | input of design                                        m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inpdesign()
{
#ifdef DEBUG
dstrc_enter("inpctrdesign");
#endif

design = (DESIGN*)CCACALLOC(1,sizeof(DESIGN));
/*----------------------- input of design size (number of dnodes etc... */
inp_designsize();
/*------------------------------------------------input of design nodes */
inp_dnode();
/*------------------------------------------------input of design lines */
inp_dline();
/*---------------------------------------------input of design surfaces */
inp_dsurface();
/*----------------------------------------------input of design volumes */
inp_dvolume();
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inpctrdesign */




/*----------------------------------------------------------------------*
 | input of design size                                   m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_designsize()
{
INT  i,ierr;
INT  numdnode;
INT  numdline;
INT  numdsurf;
INT  numdvol;
#ifdef DEBUG
dstrc_enter("inp_designsize");
#endif

/*-------------------------------------- count number of design objects */
if (frfind("--DESIGN DESCRIPTION")==0)
  dserror("frfind: DESIGN DESCRIPTION is not in input file");

frread();
frint("NDPOINT",&numdnode,&ierr);
if (!ierr) dserror("Cannot read design");
frread();
frint("NDLINE",&numdline,&ierr);
if (!ierr) dserror("Cannot read design");
frread();
frint("NDSURF",&numdsurf,&ierr);
if (!ierr) dserror("Cannot read design");
frread();
frint("NDVOL",&numdvol,&ierr);
if (!ierr) dserror("Cannot read design");
/*-----------------------------------------------put sizes to structure */
design->ndnode = numdnode;
design->ndline = numdline;
design->ndsurf = numdsurf;
design->ndvol  = numdvol;
/*----------------------------------------------allocate design vectors */
design->dnode = (DNODE*)CCACALLOC(design->ndnode,sizeof(DNODE));
for (i=0; i<design->ndnode; i++) design->dnode[i].Id = i;

design->dline = (DLINE*)CCACALLOC(design->ndline,sizeof(DLINE));
for (i=0; i<design->ndline; i++) design->dline[i].Id = i;

design->dsurf = (DSURF*)CCACALLOC(design->ndsurf,sizeof(DSURF));
for (i=0; i<design->ndsurf; i++) design->dsurf[i].Id = i;

design->dvol = (DVOL*)CCACALLOC(design->ndvol,sizeof(DVOL));
for (i=0; i<design->ndvol; i++) design->dvol[i].Id = i;
/*-------------------------------------------------------------- rewind */
frrewind();
/*----------------------------------------------------------------------*/

#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inp_designsize */






/*----------------------------------------------------------------------*
 | input of design nodes                                  m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_dnode()
{
INT    i;
INT    readID;
DNODE *actdnode;
#ifdef DEBUG
dstrc_enter("inp_dnode");
#endif
/*-------------------------------------------------------------- rewind */
frrewind();
/*------------------------ now read the description of the design nodes */
if (frfind("--DESIGN POINTS")==0)
  dserror("frfind: DESIGN POINTS is not in input file");

frread();
for (i=0; i<design->ndnode; i++)
{
   actdnode = &(design->dnode[i]);
   /*------------------------------------- GID writes fortran style Ids */
   readID   = actdnode->Id+1;
   /*---------------------------------------------- find the word POINT */
   read_1_dnode(actdnode,readID);
}
/*----------------------------------------------------------------------*/

#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inp_dnode */

/*----------------------------------------------------------------------*
 | input of one design nodes                              m.gee 1/02    |
 *----------------------------------------------------------------------*/
void read_1_dnode(DNODE *dnode, INT readId)
{
INT    i,ierr;
#ifdef DEBUG
dstrc_enter("read_1_dnode");
#endif
/*----------------------------------------------------------------------*/
ierr=0;
frchk("POINT",&ierr);
while (!ierr) {frread(); frchk("POINT",&ierr);}

frchk("Num:",&ierr);
while (!ierr) {frread(); frchk("Num:",&ierr);}
frint("Num:",&i,&ierr);
if (!ierr) dserror("Cannot read DNODE");

if (i != readId) dserror("DNODEs got mixed up");

frint("conditions:",&(dnode->ncond),&ierr);
if (!ierr) dserror("Cannot read DNODE");

frchk("Coord:",&ierr);
while (!ierr) {frread(); frchk("Coord:",&ierr);}
frdouble_n("Coord:",dnode->x,3,&ierr);
if (!ierr) dserror("Cannot read DNODE");

frchk("END POINT",&ierr);
while (!ierr) {frread(); frchk("END POINT",&ierr);}
frread();
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of read_1_dnode */








/*----------------------------------------------------------------------*
 | input of design lines                                  m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_dline()
{
INT    i;
DLINE *actdline;
INT    readID;
#ifdef DEBUG
dstrc_enter("inp_dline");
#endif

/*------------------------ now read the description of the design lines */
if (frfind("--DESIGN LINES")==0)
  dserror("frfind: DESIGN LINES is not in input file");

frread();

for (i=0; i<design->ndline; i++)
{
   actdline = &(design->dline[i]);
   /*------------------------------------- GID writes fortran style Ids */
   readID   = actdline->Id+1;
   /*---------------------------------------------- find the word POINT */
   read_1_dline(actdline,readID);
}
/*----------------------------------------------------------------------*/

#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inp_dline */






/*----------------------------------------------------------------------*
 | input of one design line                               m.gee 1/02    |
 *----------------------------------------------------------------------*/
void read_1_dline(
    DLINE              *dline,
    INT                 readId
    )

{

  INT    i,ierr;


#ifdef DEBUG
  dstrc_enter("read_1_dline");
#endif


  ierr=0;


  /* read line typ */
  frchk("LINE",&ierr);
  while (!ierr) {frread(); frchk("LINE",&ierr);}
  frchk("STLINE",&ierr);
  if (ierr==1) dline->typ = stline;
  frchk("NURBLINE",&ierr);
  if (ierr==1) dline->typ = nurbline;
  frchk("ARCLINE",&ierr);
  if (ierr==1) dline->typ = arcline;


  /* read line number */
  frchk("Num:",&ierr);
  while (!ierr) {frread(); frchk("Num:",&ierr);}
  frint("Num:",&i,&ierr);
  if (!ierr) dserror("Cannot read DLINE");
  if (i != readId) dserror("DLINEs got mixed up");


  /* read number of conditions to this line */
  frint("conditions:",&(dline->ncond),&ierr);
  if (!ierr) dserror("Cannot read DLINE");


  /* read dpoints connected to this line */
  frchk("Points:",&ierr);
  while (!ierr) {frread(); frchk("Points:",&ierr);}
  frint_n("Points:",dline->my_dnodeId,2,&ierr);
  if (!ierr) dserror("Cannot read DLINE");
  dline->my_dnodeId[0]--;
  dline->my_dnodeId[1]--;
  dline->ndnode=2;


  /* read arcline properties */
  if(dline->typ == arcline)
  {
    dline->props.arcline = (ARCLINE*)CCACALLOC(1,sizeof(ARCLINE));
    if (dline->props.arcline==NULL) dserror("Allocation of dline failed");

    /* find center point */
    frchk("2D center",&ierr);
    while (!ierr) {frread(); frchk("2D center",&ierr);}

    /* read center point */
    frdouble_n("2D center=",&(dline->props.arcline->center[0]),2,&ierr);
    if (!ierr) dserror("Cannot read DLINE");

    /* read radius of this arcline */
    frdouble("radius=",&(dline->props.arcline->radius),&ierr);
    if (!ierr) dserror("Cannot read DLINE");

    /* read initang of this arcline */
    frdouble("initang=",&(dline->props.arcline->initang),&ierr);
    if (!ierr) dserror("Cannot read DLINE");

    /* read endang of this arcline */
    frdouble("endang=",&(dline->props.arcline->endang),&ierr);
    if (!ierr) dserror("Cannot read DLINE");


    /* find m[1,1] */
    frchk("m[1,1]",&ierr);
    while (!ierr) {frread(); frchk("m[1,1]",&ierr);}
    /* read first line */
    frdouble("m[1,1]=",&(dline->props.arcline->trans_matrix[0][0]),&ierr);
    if (!ierr) dserror("Cannot read DLINE");
    frdouble("m[1,2]=",&(dline->props.arcline->trans_matrix[0][1]),&ierr);
    if (!ierr) dserror("Cannot read DLINE");
    frdouble("m[1,3]=",&(dline->props.arcline->trans_matrix[0][2]),&ierr);
    if (!ierr) dserror("Cannot read DLINE");
    frdouble("m[1,4]=",&(dline->props.arcline->trans_matrix[0][3]),&ierr);
    if (!ierr) dserror("Cannot read DLINE");

    /* find m[2,1] */
    frchk("m[2,1]",&ierr);
    while (!ierr) {frread(); frchk("m[2,1]",&ierr);}
    /* read first line */
    frdouble("m[2,1]=",&(dline->props.arcline->trans_matrix[1][0]),&ierr);
    if (!ierr) dserror("Cannot read DLINE");
    frdouble("m[2,2]=",&(dline->props.arcline->trans_matrix[1][1]),&ierr);
    if (!ierr) dserror("Cannot read DLINE");
    frdouble("m[2,3]=",&(dline->props.arcline->trans_matrix[1][2]),&ierr);
    if (!ierr) dserror("Cannot read DLINE");
    frdouble("m[2,4]=",&(dline->props.arcline->trans_matrix[1][3]),&ierr);
    if (!ierr) dserror("Cannot read DLINE");

    /* find m[3,1] */
    frchk("m[3,1]",&ierr);
    while (!ierr) {frread(); frchk("m[3,1]",&ierr);}
    /* read first line */
    frdouble("m[3,1]=",&(dline->props.arcline->trans_matrix[2][0]),&ierr);
    if (!ierr) dserror("Cannot read DLINE");
    frdouble("m[3,2]=",&(dline->props.arcline->trans_matrix[2][1]),&ierr);
    if (!ierr) dserror("Cannot read DLINE");
    frdouble("m[3,3]=",&(dline->props.arcline->trans_matrix[2][2]),&ierr);
    if (!ierr) dserror("Cannot read DLINE");
    frdouble("m[3,4]=",&(dline->props.arcline->trans_matrix[2][3]),&ierr);
    if (!ierr) dserror("Cannot read DLINE");

    /* find m[4,1] */
    frchk("m[4,1]",&ierr);
    while (!ierr) {frread(); frchk("m[4,1]",&ierr);}
    /* read first line */
    frdouble("m[4,1]=",&(dline->props.arcline->trans_matrix[3][0]),&ierr);
    if (!ierr) dserror("Cannot read DLINE");
    frdouble("m[4,2]=",&(dline->props.arcline->trans_matrix[3][1]),&ierr);
    if (!ierr) dserror("Cannot read DLINE");
    frdouble("m[4,3]=",&(dline->props.arcline->trans_matrix[3][2]),&ierr);
    if (!ierr) dserror("Cannot read DLINE");
    frdouble("m[4,4]=",&(dline->props.arcline->trans_matrix[3][3]),&ierr);
    if (!ierr) dserror("Cannot read DLINE");


    /* find total length */
    frchk("TotalLength",&ierr);
    while (!ierr) {frread(); frchk("TotalLength",&ierr);}

    /* read total length */
    frdouble("TotalLength=",&(dline->props.arcline->total_length),&ierr);
    if (!ierr) dserror("Cannot read DLINE");
  }


  /* read stline properties */
  if(dline->typ == stline)
  {
    dline->props.stline = (STLINE*)CCACALLOC(1,sizeof(STLINE));
    if (dline->props.stline==NULL) dserror("Allocation of dline failed");
  }


  /* read nurbline properties */
  if(dline->typ == nurbline)
  {
    dline->props.nurbline = (NURBLINE*)CCACALLOC(1,sizeof(NURBLINE));
    if (dline->props.nurbline==NULL) dserror("Allocation of dline failed");

    /* find "Number of Control Points" */
    frchk("Number of Control Points",&ierr);
    while (!ierr) {frread(); frchk("Number of Control Points",&ierr);}

    /* read "Number of Control Points" */
    frint("Number of Control Points=",&(dline->props.nurbline->num_cp),&ierr);
    if (!ierr) dserror("Cannot read DLINE");

    /* allocate cp */
    dline->props.nurbline->cp =
      (DOUBLE**)CCAMALLOC((dline->props.nurbline->num_cp*sizeof(DOUBLE*)));

    dline->props.nurbline->cp[0] =
      (DOUBLE*) CCAMALLOC((dline->props.nurbline->num_cp*3*sizeof(DOUBLE)));

    for (i=1; i<dline->props.nurbline->num_cp; i++)
      dline->props.nurbline->cp[i] = &(dline->props.nurbline->cp[0][i*3]);


    /* read degree */
    frint("Degree=",&(dline->props.nurbline->degree),&ierr);
    if (!ierr) dserror("Cannot read DLINE");


    /* read cp's */
    for (i=0; i<dline->props.nurbline->num_cp; i++)
    {
      /* read next line */
      frread();
      frdouble_n("coords:",dline->props.nurbline->cp[i],3,&ierr);
      if (!ierr) dserror("Cannot read DLINE");
    }

    /* find "Number of knots" */
    frchk("Number of knots",&ierr);
    while (!ierr) {frread(); frchk("Number of knots",&ierr);}

    /* read "Number of knots" */
    frint("Number of knots=",&(dline->props.nurbline->num_knots),&ierr);
    if (!ierr) dserror("Cannot read DLINE");

    /* allocate knots */
    dline->props.nurbline->knots =
      (DOUBLE*)CCAMALLOC((dline->props.nurbline->num_knots*sizeof(DOUBLE)));


    /* read knots */
    for (i=0; i<dline->props.nurbline->num_knots; i++)
    {
      /* read next line */
      frread();
      frdouble("value=",&(dline->props.nurbline->knots[i]),&ierr);
      if (!ierr) dserror("Cannot read DLINE");
    }


    /* allocate weights */
    dline->props.nurbline->weights =
      (DOUBLE*)CCAMALLOC((dline->props.nurbline->num_cp*sizeof(DOUBLE)));


    /* find "ational" */
    frchk("ational",&ierr);
    while (!ierr) {frread(); frchk("ational",&ierr);}


    frchk("Rational weights:",&ierr);
    if (ierr==1)
    {
      dline->props.nurbline->rational = 1;

      /* read weights */
      for (i=0; i<dline->props.nurbline->num_cp; i++)
      {
        /* read next line */
        frread();
        frdouble("",&(dline->props.nurbline->weights[i]),&ierr);
        if (!ierr) dserror("Cannot read DLINE");
      }

    }  /* if (ierr==1) */

    frchk("Non rational",&ierr);
    if (ierr==1)
    {
      dline->props.nurbline->rational = 0;

      /* read weights */
      for (i=0; i<dline->props.nurbline->num_cp; i++)
      {
        dline->props.nurbline->weights[i] = 1.0;
      }
    }


    /* find total length */
    frchk("TotalLength",&ierr);
    while (!ierr) {frread(); frchk("TotalLength",&ierr);}

    /* read total length */
    frdouble("TotalLength=",&(dline->props.nurbline->total_length),&ierr);
    if (!ierr) dserror("Cannot read DLINE");
  }


  frchk("END",&ierr);
  while (!ierr) {frread(); frchk("END",&ierr);}
  frread();


#ifdef DEBUG
  dstrc_exit();
#endif

  return;

} /* end of read_1_dline */




/*----------------------------------------------------------------------*
 | input of design surfaces                               m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_dsurface()
{
INT    i;
DSURF *actdsurf;
INT    readId;
#ifdef DEBUG
dstrc_enter("inp_dsurface");
#endif
/*------------------------ now read the description of the design lines */
if (frfind("--DESIGN SURFACES")==0)
  dserror("frfind: DESIGN SURFACES is not in input file");

frread();
for (i=0; i<design->ndsurf; i++)
{
   actdsurf = &(design->dsurf[i]);
   /*------------------------------------- GID writes fortran style Ids */
   readId   = actdsurf->Id+1;
   /*---------------------------------------------- find the word POINT */
   read_1_dsurf(actdsurf,readId);
}
/*----------------------------------------------------------------------*/

#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inp_dsurface */



/*----------------------------------------------------------------------*
 | input of one design surface                            m.gee 1/02    |
 *----------------------------------------------------------------------*/
void read_1_dsurf(
    DSURF              *dsurf,
    INT                 readId
    )

{

  INT    i,j,ierr,ierr2;
  char   buffer[100];
  INT    endloop;


#ifdef DEBUG
  dstrc_enter("read_1_dsurf");
#endif


  ierr=0;
  ierr2=0;
  frchk("NURBSURFACE",&ierr);
  frchk("MESHSURFACE",&ierr2);
  while (!ierr && !ierr2)
  {
    frread();
    frchk("NURBSURFACE",&ierr);
    frchk("MESHSURFACE",&ierr2);
  }

  if (ierr != 0)
  {
    dsurf->typ = nurbsurf;

    /* read surf number */
    frchk("Num:",&ierr);
    while (!ierr) {frread(); frchk("Num:",&ierr);}
    frint("Num:",&i,&ierr);
    if (!ierr) dserror("Cannot read DSURF");

    if (i != readId) dserror("DSURF got mixed up");

    /* read number of conditions to this surf */
    frint("conditions:",&(dsurf->ncond),&ierr);
    if (!ierr) dserror("Cannot read DSURF");


    /* read dlines connected to this surf */
    frchk("NumLines:",&ierr);
    while (!ierr) {frread(); frchk("NumLines:",&ierr);}
    frint("NumLines:",&(dsurf->ndline),&ierr);
    if (!ierr) dserror("Cannot read DSURF");
    amdef("my_dlineIDs",&(dsurf->my_dlineId),dsurf->ndline,2,"IA");

    frchk("Line:",&ierr);
    while (!ierr) {frread(); frchk("Line:",&ierr);}
    for (i=0; i<dsurf->ndline; i++)
    {
      frint("Line:",&(dsurf->my_dlineId.a.ia[i][0]),&ierr);
      if (!ierr) dserror("Cannot read DSURF");
      dsurf->my_dlineId.a.ia[i][0]--;
      frchar("Orientation:",buffer,&ierr);
      if (!ierr) dserror("Cannot read DSURF");
      if (strncmp("SAME1ST",buffer,7)==0) dsurf->my_dlineId.a.ia[i][1]=0;
      else                                dsurf->my_dlineId.a.ia[i][1]=1;
      frread();
    }


    /* allocate nurbsurf porperties */
    dsurf->props.nurbsurf = (NURBSURF*)CCACALLOC(1,sizeof(NURBSURF));
    if (dsurf->props.nurbsurf==NULL) dserror("Allocation of dsurf failed");

    /* find "Number of Control Points" */
    frchk("Number of Control Points",&ierr);
    while (!ierr) {frread(); frchk("Number of Control Points",&ierr);}

    /* read "Number of Control Points" */
    frint_n("Number of Control Points=",dsurf->props.nurbsurf->num_cp,2,&ierr);
    if (!ierr) dserror("Cannot read DSURF");

    /* read degree */
    frint_n("Degree=",dsurf->props.nurbsurf->degree,2,&ierr);
    if (!ierr) dserror("Cannot read DSURF");



    /* ---------------- *
     *  Control points  *
     * ---------------- */

    /* allocate cp */
    dsurf->props.nurbsurf->cp       =
      (DOUBLE***)
      CCAMALLOC((dsurf->props.nurbsurf->num_cp[0]*sizeof(DOUBLE**)));

    dsurf->props.nurbsurf->cp[0]    =
      (DOUBLE**)
      CCAMALLOC((dsurf->props.nurbsurf->num_cp[0]*dsurf->props.nurbsurf->num_cp[1]*sizeof(DOUBLE*)));

    dsurf->props.nurbsurf->cp[0][0] =
      (DOUBLE*)
      CCAMALLOC((dsurf->props.nurbsurf->num_cp[0]*dsurf->props.nurbsurf->num_cp[1]*3*sizeof(DOUBLE)));

    for (i=1; i<dsurf->props.nurbsurf->num_cp[0]; i++)
      dsurf->props.nurbsurf->cp[i]    =
        &(dsurf->props.nurbsurf->cp[0][i*dsurf->props.nurbsurf->num_cp[1]]);

    endloop=dsurf->props.nurbsurf->num_cp[0]*dsurf->props.nurbsurf->num_cp[1];
    for (i=1; i<endloop; i++)
      dsurf->props.nurbsurf->cp[0][i] = &(dsurf->props.nurbsurf->cp[0][0][i*3]);


    /* read cp's */
    for (j=0; j<dsurf->props.nurbsurf->num_cp[1]; j++)
    {
      for (i=0; i<dsurf->props.nurbsurf->num_cp[0]; i++)
      {
        /* read next line */
        frread();
        frdouble_n("coords:",dsurf->props.nurbsurf->cp[i][j],3,&ierr);
        if (!ierr) dserror("Cannot read DSURF");
      }
    }



    /* ---------------------- *
     *  Knots in U direction  *
     * ---------------------- */

    /* find "Number of knots in U" */
    frchk("Number of knots in U",&ierr);
    while (!ierr) {frread(); frchk("Number of knots in U",&ierr);}

    /* read "Number of knots in U" */
    frint("Number of knots in U=",&(dsurf->props.nurbsurf->num_knots[0]),&ierr);
    if (!ierr) dserror("Cannot read DSURF");

    /* allocate knots */
    dsurf->props.nurbsurf->knots[0] =
      (DOUBLE*)CCAMALLOC((dsurf->props.nurbsurf->num_knots[0]*sizeof(DOUBLE)));

    /* read knots */
    for (i=0; i<dsurf->props.nurbsurf->num_knots[0]; i++)
    {
      /* read next line */
      frread();
      frdouble("value=",&(dsurf->props.nurbsurf->knots[0][i]),&ierr);
      if (!ierr) dserror("Cannot read DSURF");
    }



    /* ---------------------- *
     *  Knots in V direction  *
     * ---------------------- */

    /* find "Number of knots in V" */
    frchk("Number of knots in V",&ierr);
    while (!ierr) {frread(); frchk("Number of knots in V",&ierr);}

    /* read "Number of knots in V" */
    frint("Number of knots in V=",&(dsurf->props.nurbsurf->num_knots[1]),&ierr);
    if (!ierr) dserror("Cannot read DSURF");

    /* allocate knots */
    dsurf->props.nurbsurf->knots[1] =
      (DOUBLE*)CCAMALLOC((dsurf->props.nurbsurf->num_knots[1]*sizeof(DOUBLE)));

    /* read knots */
    for (i=0; i<dsurf->props.nurbsurf->num_knots[1]; i++)
    {
      /* read next line */
      frread();
      frdouble("value=",&(dsurf->props.nurbsurf->knots[1][i]),&ierr);
      if (!ierr) dserror("Cannot read DSURF");
    }




    /* ------------------ *
     *  Rational Weights  *
     * ------------------ */

    /* allocate weights */
    dsurf->props.nurbsurf->weights =
      (DOUBLE**)
      CCAMALLOC((dsurf->props.nurbsurf->num_cp[0]*sizeof(DOUBLE*)));

    dsurf->props.nurbsurf->weights[0] =
      (DOUBLE*)
      CCAMALLOC((dsurf->props.nurbsurf->num_cp[0]*dsurf->props.nurbsurf->num_cp[1]*sizeof(DOUBLE)));

    for (i=1; i<dsurf->props.nurbsurf->num_cp[0]; i++)
      dsurf->props.nurbsurf->weights[i] =
        &(dsurf->props.nurbsurf->weights[0][i*dsurf->props.nurbsurf->num_cp[1]]);


    /* find "ational" */
    frchk("ational",&ierr);
    while (!ierr) {frread(); frchk("ational",&ierr);}


    frchk("Rational weights:",&ierr);
    if (ierr==1)
    {
      dsurf->props.nurbsurf->rational = 1;

      /* read weights */
      for (j=0; j<dsurf->props.nurbsurf->num_cp[1]; j++)
      {
        for (i=0; i<dsurf->props.nurbsurf->num_cp[0]; i++)
        {
          /* read next line */
          frread();
          frdouble("",&(dsurf->props.nurbsurf->weights[i][j]),&ierr);
          if (!ierr) dserror("Cannot read DSURF");
        }
      }
    }


    frchk("Non rational",&ierr);
    if (ierr==1)
    {
      dsurf->props.nurbsurf->rational = 0;

      for (i=0; i<dsurf->props.nurbsurf->num_cp[0]; i++)
        for (j=0; j<dsurf->props.nurbsurf->num_cp[1]; j++)
          dsurf->props.nurbsurf->weights[i][j] = 1.0;

    }


    /* find "IsTrimmed" */
    frchk("IsTrimmed",&ierr);
    while (!ierr) {frread(); frchk("IsTrimmed",&ierr);}

    /* read "IsTrimmed" */
    frint("IsTrimmed:",&(dsurf->props.nurbsurf->istrimmed),&ierr);
    if (!ierr) dserror("Cannot read DSURF");



    /* find "Center:" */
    frchk("Center:",&ierr);
    while (!ierr) {frread(); frchk("Center:",&ierr);}

    /* read "Center:" */
    frdouble("Center:",dsurf->props.nurbsurf->center,&ierr);
    if (!ierr) dserror("Cannot read DSURF");


    /* find "Normal:" */
    frchk("Normal:",&ierr);
    while (!ierr) {frread(); frchk("Normal:",&ierr);}

    /* read "Normal:" */
    frdouble("Normal:",dsurf->props.nurbsurf->normal,&ierr);
    if (!ierr) dserror("Cannot read DSURF");



    frchk("END NURBSURFACE",&ierr);
    while (!ierr) {frread(); frchk("END NURBSURFACE",&ierr);}
    frread();
  }
  else if (ierr2 != 0)
  {
    dsurf->typ = meshsurf;

    /* read surf number */
    frchk("Num:",&ierr);
    while (!ierr) {frread(); frchk("Num:",&ierr);}
    frint("Num:",&i,&ierr);
    if (!ierr) dserror("Cannot read DSURF");

    if (i != readId) dserror("DSURF got mixed up");

    /* read number of conditions to this surf */
    frint("conditions:",&(dsurf->ncond),&ierr);
    if (!ierr) dserror("Cannot read DSURF");


    /* read dlines connected to this surf */
    frchk("NumLines:",&ierr);
    while (!ierr) {frread(); frchk("NumLines:",&ierr);}
    frint("NumLines:",&(dsurf->ndline),&ierr);
    if (!ierr) dserror("Cannot read DSURF");
    amdef("my_dlineIDs",&(dsurf->my_dlineId),dsurf->ndline,2,"IA");

    frchk("Line:",&ierr);
    while (!ierr) {frread(); frchk("Line:",&ierr);}
    for (i=0; i<dsurf->ndline; i++)
    {
      frint("Line:",&(dsurf->my_dlineId.a.ia[i][0]),&ierr);
      if (!ierr) dserror("Cannot read DSURF");
      dsurf->my_dlineId.a.ia[i][0]--;
      frchar("Orientation:",buffer,&ierr);
      if (!ierr) dserror("Cannot read DSURF");
      if (strncmp("SAME1ST",buffer,7)==0) dsurf->my_dlineId.a.ia[i][1]=0;
      else                                dsurf->my_dlineId.a.ia[i][1]=1;
      frread();
    }

    /* "Flag: 0" is supposed to be the last line in each meshsurf entry. */
    frchk("Flag: 0",&ierr);
    while (!ierr) {frread(); frchk("Flag: 0",&ierr);}
    frread();
  }
  else
    dserror("confused while reading design surfaces");

  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of read_1_dsurf */





/*----------------------------------------------------------------------*
 | input of design volumes                                m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_dvolume()
{
INT   i;
DVOL *actdvol;
INT   readId;
#ifdef DEBUG
dstrc_enter("inp_dvolume");
#endif

/*---------------------- now read the description of the design volumes */
if (frfind("--DESIGN VOLUMES")==0)
  dserror("frfind: DESIGN VOLUMES is not in input file");

frread();
for (i=0; i<design->ndvol; i++)
{
   actdvol = &(design->dvol[i]);
   /*------------------------------------- GID writes fortran style Ids */
   readId   = actdvol->Id+1;
   /*---------------------------------------------- find the word POINT */
   read_1_dvol(actdvol,readId);
}
/*----------------------------------------------------------------------*/

#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inp_dvolume */

/*----------------------------------------------------------------------*
 | input of one design volume                             m.gee 1/02    |
 *----------------------------------------------------------------------*/
void read_1_dvol(DVOL *dvol, INT readId)
{
INT    i,ierr;
char   buffer[100];
#ifdef DEBUG
dstrc_enter("read_1_dvol");
#endif
/*----------------------------------------------------------------------*/
ierr=0;
frchk("VOLUME",&ierr);
while (!ierr) {frread(); frchk("VOLUME",&ierr);}

frchk("Num:",&ierr);
while (!ierr) {frread(); frchk("Num:",&ierr);}
frint("Num:",&i,&ierr);
if (!ierr) dserror("Cannot read DVOL");

if (i != readId) dserror("DVOLs got mixed up");

frint("conditions:",&(dvol->ncond),&ierr);
if (!ierr) dserror("Cannot read DVOL");

frchk("NumSurfaces:",&ierr);
while (!ierr) {frread(); frchk("NumSurfaces:",&ierr);}
frint("NumSurfaces:",&(dvol->ndsurf),&ierr);
if (!ierr) dserror("Cannot read DVOL");
amdef("my_dsurfIDs",&(dvol->my_dsurfId),dvol->ndsurf,2,"IA");

frchk("Surface:",&ierr);
while (!ierr) {frread(); frchk("Surface:",&ierr);}
for (i=0; i<dvol->ndsurf; i++)
{
   frint("Surface:",&(dvol->my_dsurfId.a.ia[i][0]),&ierr);
   if (!ierr) dserror("Cannot read DVOL");
   dvol->my_dsurfId.a.ia[i][0]--;
   frchar("Orientation:",buffer,&ierr);
   if (!ierr) dserror("Cannot read DVOL");
   if (strncmp("SAME1ST",buffer,7)==0) dvol->my_dsurfId.a.ia[i][1]=0;
   else                                dvol->my_dsurfId.a.ia[i][1]=1;
   frread();
}
frchk("END VOLUME",&ierr);
while (!ierr) {frread(); frchk("END VOLUME",&ierr);}
frread();
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of read_1_dvol */
