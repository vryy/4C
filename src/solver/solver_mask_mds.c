#include "../headers/standardtypes.h"
#include "../headers/solution.h"
int cmp_int(const void *a, const void *b );
double cmp_double(const void *a, const void *b );


/*----------------------------------------------------------------------*
 |                                                            al  10/01 |
 |  calculate the mask of a column pointer, row index sparse  matrix    |
 *----------------------------------------------------------------------*/
void mask_mds(FIELD        *actfield, 
              PARTITION    *actpart, 
              SOLVAR       *actsolv,
              INTRA        *actintra, 
              ML_ARRAY_MDS *mds)
{
int       i,j,k,l;
int       numeq;
int       nnz;
int     **dof_connect;
#ifdef DEBUG 
dstrc_enter("mask_mds");
#endif
/*------------------------------------------- put total size of problem */
 mds->numeq = actfield->numeq;
/*---------------------------------------------- allocate vector colstr */
 amdef("colstr1",&(mds->colstr),numeq,1,"IV");
 amzero(&(mds->colstr));
/*---------------------------------- calculate dof connectivity list ---*/
 dof_connect = (int**)CALLOC(mds->numeq,sizeof(int*));
 if (!dof_connect) dserror("Allocation of dof_connect failed");

 dofconnectivity(actfield,dof_connect,&nnz);

 mds->nnz  = nnz;
/*--------------------------------------- allocate colstr and rowind ---*/
amdef("rowind" ,&(mds->rowind)  ,(mds->nnz  +1),1,"IV");
amdef("colstr" ,&(mds->colstr)  ,(mds->numeq+1),1,"IV"); 
/*amdef("value"  ,&(mds->value )  ,(mds->nnz  +1),1,"DV");*/ 
/*amzero(&(mds->value));*/

mds_make_colstr_rowind(actsolv, mds,dof_connect, mds->numeq);

/*---------------------------------------- delete the array dof_connect */
 for (i=0; i<mds->numeq; i++)
{
   if (dof_connect[i]) FREE(dof_connect[i]);
}
FREE(dof_connect);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mask_mds */







