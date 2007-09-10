/*!--------------------------------------------------------------------
 * \file
 * \brief contains the routine
 * c1_mat_ogden_uncoupled:	calculates the principal stretches and
 * the eigenvectors of right Cauchy Green
 * and calls the decoupled ogden material
 * law at the materials folder and gives
 * back the stress and constitutive tensor
 * c1_calc_eigenval_eigenvec_jacobi:	spectral decomposition of a
 * 3x3 tensor following Numerical
 * Recip 11.1
 * c1_rotation		helper function for spectral decomposition
 * <pre>
 * Maintainer:		Robert Metzke
 * metzke@lnm.mw.tum.de
 * http://www.lnm.mw.tum.de/Members/metzke
 * 089 - 289 15244
 * </pre>
 *
 *-------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_BRICK1
#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"
#ifdef D_MAT
#include "../materials/mat_prototypes.h"
#endif
/*!----------------------------------------------------------------------
 * \brief consitutive matrix for rubber like elastic orthotropic material
 *
 * <pre>                                                            rm 05/07
 * This routine prepares the data for the decoupled Ogden material in the
 * materials folder (D_MAT). In order to work you need to use D_MAT flag in
 * your configure file.
 * history:
 * compressible ogden-material                               m.gee 6/03
 * split in volumetric and deviatoric strains
 * adapted from shell 8 to brick                               rm 03/07
 * </pre>
 * \param  DOUBLE    *mat	        (i)  material information
 * \param  DOUBLE    *disd         (i)  displacement derivatives
 * \param  DOUBLE    *stress       (o)  stress tensor
 * \param  DOUBLE    **d           (o)  constitutive matrix
 *
 * \warning There is nothing special to this routine
 * \return void
 * \sa c1_calc_eigenval_eigenvec_jacobi,mat_el_ogden_decoupled
 *----------------------------------------------------------------------*/
void c1_mat_ogden_decoupled(
        COMPOGDEN *mat,
        DOUBLE *disd,
        DOUBLE *stress,
        DOUBLE **d) {
    INT i, j, k;
    static DOUBLE **FT;
    static ARRAY FT_a;
    static DOUBLE **CG;
    static ARRAY CG_a;
    static DOUBLE **N;
    static ARRAY N_a;
    static DOUBLE *lam;
    static ARRAY lam_a;
    static DOUBLE *lam2;
    static ARRAY lam2_a;
    static DOUBLE **PK2;
    static ARRAY PK2_a;

    if (FT==NULL) {
        FT = amdef("FT", &FT_a, 3, 3, "DA");
        CG = amdef("CG", &CG_a, 3, 3, "DA");
        N = amdef("N", &N_a, 3, 3, "DA");
        lam = amdef("lam", &lam_a, 3, 1, "DV");
        lam2 = amdef("lam2", &lam2_a, 3, 1, "DV");
        PK2 = amdef("PK2", &PK2_a, 3, 3, "DA");
    }

    amzero(&FT_a);
    amzero(&CG_a);
    amzero(&N_a);
    amzero(&lam_a);
    amzero(&lam2_a);
    amzero(&PK2_a);

    DOUBLE  C[3][3][3][3];

#ifdef DEBUG
dstrc_enter("c1_mat_ogden_uncoupled");
#endif
/*----------------------------------- Deformation Gradient, transposed */
FT[0][0]=disd[0]+1.0;
FT[1][1]=disd[1]+1.0;
FT[2][2]=disd[2]+1.0;
FT[0][1]=disd[3];
FT[1][0]=disd[4];
FT[1][2]=disd[5];
FT[2][1]=disd[6];
FT[0][2]=disd[7];
FT[2][0]=disd[8];
/*------------------------------------------- Right Cauchy Green Tensor*/
for (k=0; k<3; k++) {
    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            CG[i][k]+=(FT[i][j]*FT[k][j]); } } }
/*------------------------------------ make spectral decoposition of CG */
c1_calc_eigenval_eigenvec_jacobi(CG, lam2, N);
/*-------------------------------------------- make principal stretches */
lam[0] = sqrt(lam2[0]);
lam[1] = sqrt(lam2[1]);
lam[2] = sqrt(lam2[2]);
#if 1
for (i=0; i<3; i++) mat->l[i] = lam[i];
#endif
/*------------------------------ call ogden material in material folder */
#ifdef D_MAT
mat_el_ogden_decoupled(mat, lam, N, PK2, C);
#else
dserror("Please use D_MAT in your configure file, mate!");
#endif
/*------------------------------------------------------------ Stresses */
stress[0] = PK2[0][0];
stress[1] = PK2[1][1];
stress[2] = PK2[2][2];
stress[3] = PK2[0][1];
stress[4] = PK2[1][2];
stress[5] = PK2[0][2];
/*---------------------------------------------------Constitutive Matrix*/
d[0][0]=C[0][0][0][0];
d[0][1]=C[0][0][1][1];
d[0][2]=C[0][0][2][2];
d[0][3]=C[0][0][1][0];
d[0][4]=C[0][0][2][1];
d[0][5]=C[0][0][2][0];

d[1][0]=C[1][1][0][0];
d[1][1]=C[1][1][1][1];
d[1][2]=C[1][1][2][2];
d[1][3]=C[1][1][1][0];
d[1][4]=C[1][1][2][1];
d[1][5]=C[1][1][2][0];

d[2][0]=C[2][2][0][0];
d[2][1]=C[2][2][1][1];
d[2][2]=C[2][2][2][2];
d[2][3]=C[2][2][1][0];
d[2][4]=C[2][2][2][1];
d[2][5]=C[2][2][2][0];

d[3][0]=C[1][0][0][0];
d[3][1]=C[1][0][1][1];
d[3][2]=C[1][0][2][2];
d[3][3]=C[1][0][1][0];
d[3][4]=C[1][0][2][1];
d[3][5]=C[1][0][2][0];

d[4][0]=C[2][1][0][0];
d[4][1]=C[2][1][1][1];
d[4][2]=C[2][1][2][2];
d[4][3]=C[2][1][1][0];
d[4][4]=C[2][1][2][1];
d[4][5]=C[2][1][2][0];

d[5][0]=C[2][0][0][0];
d[5][1]=C[2][0][1][1];
d[5][2]=C[2][0][2][2];
d[5][3]=C[2][0][1][0];
d[5][4]=C[2][0][2][1];
d[5][5]=C[2][0][2][0];

#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of c1_mat_ogden_uncoupled */



/*----------------------------------------------------------------------*
 * eigenvalues and eigenvectors of a matrix, applying the jacobi method |
 * taken from Numerical Recipes and specialized to 3x3 case    rm 03.07 |
 *----------------------------------------------------------------------*/
void c1_calc_eigenval_eigenvec_jacobi(DOUBLE **C, DOUBLE *d, DOUBLE **V) {
#ifdef DEBUG
    dstrc_enter("c1_calc_eigenval_eigenvec_jacobi");
#endif
INT n=3;
INT j, iq, ip, i;
DOUBLE tresh, theta, tau, t, sm, s, h, g, c;
/*DOUBLE C_loc[3][3];*/
DOUBLE **C_loc;
ARRAY C_loc_a;

C_loc = amdef("C_loc", &C_loc_a, 3, 3, "DA");
DOUBLE b[n];
DOUBLE z[n];
/*------------------------------------------------------ Einheitsmatrix */
for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
        C_loc[i][j]=C[i][j];
    }
}
for (ip=0; ip<n; ip++) {
    for (iq=0; iq<n; iq++) {
        V[ip][iq]=0.0;
        V[ip][ip]=1.0;
    }
}
/*---------------------- vectors b,d = diagonal elements of A, z = NULL */
for (ip=0;ip<n;ip++) {
    b[ip]=d[ip]=C_loc[ip][ip];
    z[ip]=0.0;
}
/*------------------------------------------------------- 50 Iterationen*/
for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=0;ip<n-1;ip++) {
        for (iq=ip+1;iq<n;iq++) {
            sm += fabs(C_loc[ip][iq]);
        }
    }
    if (sm == 0.0) {
        return;
    }
    if (i < 4)
        tresh=0.2*sm/(n*n);
    else
        tresh=0.0;
    for (ip=0;ip<n-1;ip++) {
        for (iq=ip+1;iq<n;iq++) {
            g=100.0 * fabs(C_loc[ip][iq]);
            if (i > 4 && (fabs(d[ip])+g) == fabs(d[ip]) && (fabs(d[iq])+g) == fabs(d[iq]))
                C_loc[ip][iq]=0.0;
            else if (fabs(C_loc[ip][iq]) > tresh){
                h=d[iq]-d[ip];
                if ((fabs(h)+g) == fabs(h))
                    t=(C_loc[ip][iq])/h;
                else {
                    theta=0.5f*h/(C_loc[ip][iq]);
                    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                    if (theta < 0.0f)
                        t = -t;
                }
                c=1.0/sqrt(1+t*t);
                s=t*c;
                tau=s/(1.0+c);
                h=t*C[ip][iq];
                z[ip] -= h;
                z[iq] += h;
                d[ip] -= h;
                d[iq] += h;
                C_loc[ip][iq]=0.0;
                for (j=0;j<=ip-1;j++) {
                    c1_rotation(C_loc, j, ip, j, iq, tau, s);
                }
                for (j=ip+1;j<=iq-1;j++) {
                    c1_rotation(C_loc, ip, j, j, iq, tau, s);
                }
                for (j=iq+1;j<n;j++) {
                    c1_rotation(C_loc, ip, j, iq, j, tau, s);
                }
                for (j=0;j<n;j++) {
                    c1_rotation(V, j, ip, j, iq, tau, s);
                }
            }
        }
    }
    for (ip=0;ip<n;ip++) {
        b[ip] += z[ip];
        d[ip]=b[ip];
        z[ip]=0.0;
    }
}
}


void c1_rotation(DOUBLE **C, INT i, INT j, INT k, INT l, DOUBLE tau, DOUBLE s) {
#ifdef DEBUG
    dstrc_enter("c1_rotation");
#endif
DOUBLE g, h;
g=C[i][j];
h=C[k][l];
C[i][j]=g-s*(h+g*tau);
C[k][l]=h+s*(g-h*tau);
#ifdef DEBUG
dstrc_exit();
#endif
return;
}

#endif /*D_BRICK1 */
/*! @} (documentation module close)*/
#endif
