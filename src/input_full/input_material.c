#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 | input of materials                                     m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_material()
{
int  ierr;
int  i;
char *colpointer;
char buffer[50];
#ifdef DEBUG 
dstrc_enter("inp_material");
#endif
/*----------------------------------------------------------------------*/
mat = (MATERIAL*)calloc(genprob.nmat,sizeof(MATERIAL));
if (mat==NULL) dserror("Allocation of MATERIAL failed");
/*----------------------------------------------------------------------*/
frfind("--MATERIALS");
frread();
i=0;
while(strncmp(allfiles.actplace,"------",6)!=0)
{
   if (i==genprob.nmat) dserror("number of materials incorrect");

   frint("MAT",&(mat[i].Id),&ierr);

   frchk("MAT_fluid",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_fluid;
      mat[i].m.fluid = (FLUID*)calloc(1,sizeof(FLUID));
      if (mat[i].m.fluid==NULL) dserror("Alloocation of FLUID material failed");
      frdouble("VISCOSITY",&(mat[i].m.fluid->viscosity),&ierr);
      frdouble("DENSITY"  ,&(mat[i].m.fluid->density)  ,&ierr);
   }
   frchk("MAT_Struct_Linear_Elastic",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_lin_el;
      mat[i].m.lin_el = (LINEAR_ELASTIC*)calloc(1,sizeof(LINEAR_ELASTIC));
      if (mat[i].m.lin_el==NULL) dserror("Alloocation of LINEAR_ELASTIC material failed");
      frdouble("YOUNG",&(mat[i].m.lin_el->youngs)        ,&ierr);
      frdouble("NUE"  ,&(mat[i].m.lin_el->possionratio)  ,&ierr);
   }
   frchk("MAT_Struct_NeoHooke",&ierr);
   if (ierr==1)
   {
      mat[i].mattyp = m_neohooke;
      mat[i].m.neohooke = (NEO_HOOKE*)calloc(1,sizeof(NEO_HOOKE));
      if (mat[i].m.neohooke==NULL) dserror("Alloocation of NEO_HOOKE material failed");
      frdouble("YOUNG",&(mat[i].m.neohooke->youngs)        ,&ierr);
      frdouble("NUE"  ,&(mat[i].m.neohooke->possionratio)  ,&ierr);
   }
   i++;
/*----------------------------------------------------------------------*/
   frread();
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_material */
