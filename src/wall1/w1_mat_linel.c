#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 | constitutive matrix - linear elastic - 2D              al    9/01    |
 | plane stress, plane strain, rotational symmetry                      |
 *----------------------------------------------------------------------*/
void w1_mat_linel(LINEAR_ELASTIC *mat, double **d)
{
double e1, e2, e3, a1, b1, c1;
double ym,pv;/*------------------------------------------ mat constants */
#ifdef DEBUG 
dstrc_enter("w1_mat_linel");
#endif
/*----------------------------------------------------------------------*/
ym  = mat->youngs;
pv  = mat->possionratio;
/*----------------------------------------------------- plane stress ---*/

e1=ym/(1.0 - pv*pv);
e2=pv*e1;
e3=e1*(1.0 - pv)/2.0;

d[0][0]=e1;
d[0][1]=e2;
d[0][2]=0.0;
d[1][0]=e2;
d[1][1]=e1;
d[1][2]=0.0;
d[2][0]=0.0;
d[2][1]=0.0;
d[2][2]=e3;

/*-------------------------------- plane strain, rotational symmetry ---*/
/*
        c1=ym/(1.0+pv);
        b1=c1*pv/(1.0-2.0*pv);
        a1=b1+c1;

        d(1,1)=a1;
        d(1,2)=b1;
        d(2,1)=b1;
        d(2,2)=a1;
        d(3,3)=c1/2.0;

        d(1,4)=b1;
        d(2,4)=b1;
        d(4,1)=b1;
        d(4,2)=b1;
        d(4,4)=a1;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_mat_linel */
