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

/*----------------------------------------------------------------------*
 | prototypes of functions callable only in this file                   |
 *----------------------------------------------------------------------*/
static void inp_read_funct(char *string, INT id);


extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 *----------------------------------------------------------------------*/
extern struct _DYNAMIC *dyn;   
/*----------------------------------------------------------------------*
 |                                                          mn 02/04    |
 | number of spatual load functions   numcurve                          |
 | vector of structures of functions                                    |
 | defined in input_funct.c                                             |
 | INT                 numfunct;                                        |
 | struct _FUNCT      *funct;                                           |
 *----------------------------------------------------------------------*/
INT                 numfunct;
struct _FUNCT      *funct;

/*----------------------------------------------------------------------*
 | input of curves                                          mn 02/04    |
 *----------------------------------------------------------------------*/
void inp_cond_funct()
{
INT  ierr;
INT  i;
INT  counter;
char buffer[50];
#ifdef DEBUG 
dstrc_enter("inp_cond_funct");
#endif
/*----------------------------------------------------------------------*/
/* count the number of different functions (max=6) */
numfunct=0;
if (frfind("--FUNCT1")==1)
{
  frread();
  frchk("---",&ierr);
  if (ierr==0) (numfunct)++;
}

if (frfind("--FUNCT2")==1)
{
  frread();
  frchk("---",&ierr);
  if (ierr==0) (numfunct)++;
}

if (frfind("--FUNCT3")==1)
{
  frread();
  frchk("---",&ierr);
  if (ierr==0) (numfunct)++;
}

if (frfind("--FUNCT4")==1)
{
  frread();
  frchk("---",&ierr);
  if (ierr==0) (numfunct)++;
}

if (frfind("--FUNCT5")==1)
{
  frread();
  frchk("---",&ierr);
  if (ierr==0) (numfunct)++;
}

if (frfind("--FUNCT6")==1)
{
  frread();
  frchk("---",&ierr);
  if (ierr==0) (numfunct)++;
}


/*------------------------------------------------- allocate the curves */
funct = (FUNCT*)CCACALLOC(numfunct,sizeof(FUNCT));
/*----------------------------------------------------- read the curves */
inp_read_funct("--FUNCT1",0);
inp_read_funct("--FUNCT2",1);
inp_read_funct("--FUNCT3",2);
inp_read_funct("--FUNCT4",3);
inp_read_funct("--FUNCT5",4);
inp_read_funct("--FUNCT6",5);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_cond_funct */







/*----------------------------------------------------------------------*
 | input of funct                                           mn 02/04    |
 *----------------------------------------------------------------------*/
void inp_read_funct(char *string, INT id)
{
INT      ierr;
INT      counter=0;
char     buffer[50];
FUNCT   *actfunct;
DOUBLE   tmp[20];
DOUBLE   length;
#ifdef DEBUG 
dstrc_enter("inp_read_funct");
#endif
/*----------------------------------------------------------------------*/
/* check whether there is info on this funct */
if (frfind(string)==0) goto end;
frread();
frchk("---",&ierr);
if (ierr==1) goto end;
/* set actcurve to correct curev */
actfunct=&(funct[id]);

if (frfind(string)==0) goto end;
   frread();

/* read funct Id */
   frint("FUNCT",&(actfunct->Id),&ierr);
   if (ierr!=1) dserror("cannot read FUNCT");

/* read typ of funct */   
   frchk("LINEAR",&ierr);
   if (ierr==1) 
   {
      actfunct->functtyp = funct_lin;
      actfunct->typ.funct_lin = (FUNCT_LIN*)CCACALLOC(1,sizeof(FUNCT_LIN));
      /* read 8 double values */
      frdouble_n("LINEAR",&(tmp[0]),8,&ierr);
      actfunct->typ.funct_lin->x1[0] = tmp[0];
      actfunct->typ.funct_lin->x1[1] = tmp[1];
      actfunct->typ.funct_lin->x1[2] = tmp[2];
      actfunct->typ.funct_lin->val1  = tmp[3];
      actfunct->typ.funct_lin->x2[0] = tmp[4];
      actfunct->typ.funct_lin->x2[1] = tmp[5];
      actfunct->typ.funct_lin->x2[2] = tmp[6];
      actfunct->typ.funct_lin->val2  = tmp[7];
      /* calculate slope and offset */
      actfunct->typ.funct_lin->b     = tmp[3];
      actfunct->typ.funct_lin->m     = (tmp[7]-tmp[3]);
      actfunct->typ.funct_lin->length= sqrt ( (tmp[4]-tmp[0])*(tmp[4]-tmp[0]) +
                                              (tmp[5]-tmp[1])*(tmp[5]-tmp[1]) +
                                              (tmp[6]-tmp[2])*(tmp[6]-tmp[2]) );
   }

   frchk("QUADRATIC",&ierr);
   if (ierr==1) 
   {
      actfunct->functtyp = funct_qua;
      actfunct->typ.funct_qua = (FUNCT_QUA*)CCACALLOC(1,sizeof(FUNCT_QUA));
      /* read 6 double values */
      frdouble_n("QUADRATIC",&(tmp[0]),6,&ierr);
      actfunct->typ.funct_qua->x1[0] = tmp[0];
      actfunct->typ.funct_qua->x1[1] = tmp[1];
      actfunct->typ.funct_qua->x1[2] = tmp[2];
      actfunct->typ.funct_qua->x2[0] = tmp[3];
      actfunct->typ.funct_qua->x2[1] = tmp[4];
      actfunct->typ.funct_qua->x2[2] = tmp[5];
      /* calculate slope and offset */
      actfunct->typ.funct_qua->length= sqrt ( (tmp[3]-tmp[0])*(tmp[3]-tmp[0]) +
                                              (tmp[4]-tmp[1])*(tmp[4]-tmp[1]) +
                                              (tmp[5]-tmp[2])*(tmp[5]-tmp[2]) );
   }

   frchk("BELTRAMI",&ierr);
   if (ierr==1) 
   {
      actfunct->functtyp = funct_bel;
      actfunct->typ.funct_bel = NULL;
   }

   frchk("KIM-MOIN",&ierr);
   if (ierr==1) 
   {
      actfunct->functtyp = funct_kim;
      actfunct->typ.funct_kim = NULL;
   }

   frchk("CYLINDER_3D",&ierr);
   if (ierr==1) 
   {
      actfunct->functtyp = funct_cyl;
      actfunct->typ.funct_cyl = (FUNCT_CYL*)CCACALLOC(1,sizeof(FUNCT_CYL));
      /* read 1 double values */
      frdouble_n("CYLINDER_3D",&(tmp[0]),1,&ierr);
      actfunct->typ.funct_cyl->um = tmp[0];
   }
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_read_funct */
