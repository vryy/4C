/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 | EAS-Ansaetze                                           m.gee 7/01    |
 *----------------------------------------------------------------------*/
void s8_eas(const INT    nhyb,
               const DOUBLE e1,
               const DOUBLE e2,
               const INT    iel,
               const INT   *eas,
               DOUBLE     **P)
{
INT       place_P=0;
const INT nrr=0;
const INT nss=3;
const INT nrs=1;
const INT mrr=6;
const INT mss=9;
const INT mrs=7;
const INT qr=2;
const INT qs=4;
const INT qt=5;
const INT sr=8;
const INT ss=10;
const INT st=11;
DOUBLE    e1e2;
DOUBLE    e1e1;
DOUBLE    e2e2;
DOUBLE    e1e1e2;
DOUBLE    e1e2e2;
DOUBLE    e1e1e2e2;
#ifdef DEBUG
dstrc_enter("s8_eas");
#endif
/*----------------------------------------------------------------------*/
   e1e2     = e1*e2;
   e1e1     = e1*e1;
   e2e2     = e2*e2;
   e1e1e2   = e1*e1e2;
   e1e2e2   = e1e2*e2;
   e1e1e2e2 = e1e2*e1e2;
/*--------------------------------------------------- nine node element */
if (iel>4)
{
/*----------------------------------------------------------------------
      MEMBRAN: E11,E12,E22 KONSTANT
  ----------------------------------------------------------------------*/
   switch(eas[0])
   {
   case 0:
   break;
   case 7:
      P[nrr][place_P]   = e2-3.0*e1e1e2;
      P[nrr][place_P+1] = e2e2-3.0*e1e1e2e2;
      P[nss][place_P+2] = e1-3.0*e1e2e2;
      P[nss][place_P+3] = e1e1-3.0*e1e1e2e2;
      P[nrs][place_P+4] = e2-3.0*e1e1e2;
      P[nrs][place_P+5] = e1-3.0*e1e2e2;
      P[nrs][place_P+6] = 1.0-3.0*(e1e1+e2e2)+9.0*e1e1e2e2;
      place_P+=7;
   break;
   case 9:
      P[nrr][place_P]   = 1.0-3.0*e1e1;;
      P[nrr][place_P+1] = e2-3.0*e1e1e2;;
      P[nss][place_P+2] = 1.0-3.0*e2e2;
      P[nss][place_P+3] = e1-3.0*e1e2e2;
      P[nrs][place_P+4] = 1.0-3.0*e1e1e2;
      P[nrs][place_P+5] = 1.0-3.0*e1e2e2;
      P[nrs][place_P+6] = e2-3.0*e1e1e2;
      P[nrs][place_P+7] = e1-3.0*e1e2e2;
      P[nrs][place_P+8] = 1.0-3.0*(e1e1+e2e2)+9.0*e1e1e2e2;
      place_P+=9;
   break;
   case 11:
      P[nrr][place_P]    = 1.0-3.0*e1e1;
      P[nrr][place_P+1]  = e2-3.0*e1e1e2;
      P[nrr][place_P+2]  = e2e2-3.0*e1e1e2e2;
      P[nss][place_P+3]  = 1.0-3.0*e2e2;
      P[nss][place_P+4]  = e1-3.0*e1e2e2;
      P[nss][place_P+5]  = e1e1-3.0*e1e1e2e2;
      P[nrs][place_P+6]  = 1.0-3.0*e1e1;
      P[nrs][place_P+7]  = 1.0-3.0*e2e2;
      P[nrs][place_P+8]  = e2-3.0*e1e1e2;
      P[nrs][place_P+9]  = e1-3.0*e1e2e2;
      P[nrs][place_P+10] = 1.0-3.0*(e1e1+e2e2)+9.0*e1e1e2e2;
      place_P+=11;
   break;
   default:
      dserror("eas: MEMBRAN: E11,E12,E22 KONSTANT other then 0,7,9,11");
   }
/*----------------------------------------------------------------------
      BIEGUNG: E11,E12,E22 LINEAR
  ----------------------------------------------------------------------*/
   switch(eas[1])
   {
   case 0:
   break;
   case 9:
      P[mrr][place_P]   = 1.0-3.0*e1e1;
      P[mrr][place_P+1] = e2-3.0*e1e1e2;
      P[mss][place_P+2] = 1.0-3.0*e2e2;
      P[mss][place_P+3] = e1-3.0*e1e2e2;
      P[mrs][place_P+4] = 1.0-3.0*e1e1;
      P[mrs][place_P+5] = 1.0-3.0*e2e2;
      P[mrs][place_P+6] = e2-3.0*e1e1e2;
      P[mrs][place_P+7] = e1-3.0*e1e2e2;
      P[mrs][place_P+8] = 1.0-3.0*(e1e1+e2e2)+9.0*e1e1e2e2;
      place_P+=9;
   break;
   case 11:
      P[mrr][place_P]    = 1.0-3.0*e1e1;
      P[mrr][place_P+1]  = e2-3.0*e1e1e2;
      P[mrr][place_P+2]  = e2e2-3.0*e1e1e2e2;
      P[mss][place_P+3]  = 1.0-3.0*e2e2;
      P[mss][place_P+4]  = e1-3.0*e1e2e2;
      P[mss][place_P+5]  = e1e1-3.0*e1e1e2e2;
      P[mrs][place_P+6]  = 1.0-3.0*e1e1;
      P[mrs][place_P+7]  = 1.0-3.0*e2e2;
      P[mrs][place_P+8]  = e2-3.0*e1e1e2;
      P[mrs][place_P+9]  = e1-3.0*e1e2e2;
      P[mrs][place_P+10] = 1.0-3.0*(e1e1+e2e2)+9.0*e1e1e2e2;
      place_P+=11;
   break;
   default:
      dserror("eas: BIEGUNG: E11,E12,E22 LINEAR other then 0,9,11");
   }
/*----------------------------------------------------------------------
      DICKENRICHTUNG: E33 LINEAR (--> 7P - FORMULIERUNG)
  ----------------------------------------------------------------------*/
   switch(eas[2])
   {
   case 0:
   break;
   case 1:
      P[st][place_P] = 1.0;
      place_P+=1;
   break;
   case 3:
      P[st][place_P]   = 1.0;
      P[st][place_P+1] = e1;
      P[st][place_P+2] = e2;
      place_P+=3;
   break;
   case 4:
      P[st][place_P]   = 1.0;
      P[st][place_P+1] = e1;
      P[st][place_P+2] = e2;
      P[st][place_P+3] = e1e2;
      place_P+=4;
   break;
   case 6:
      P[st][place_P]   = 1.0;
      P[st][place_P+1] = e1;
      P[st][place_P+2] = e2;
      P[st][place_P+3] = e1e2;
      P[st][place_P+4] = e1e1;
      P[st][place_P+5] = e2e2;
      place_P+=6;
   break;
   case 8:
      P[st][place_P]   = 1.0;
      P[st][place_P+1] = e1;
      P[st][place_P+2] = e2;
      P[st][place_P+3] = e1e2;
      P[st][place_P+4] = e1e1;
      P[st][place_P+5] = e2e2;
      P[st][place_P+6] = e1e1e2;
      P[st][place_P+7] = e1e2e2;
      place_P+=8;
   break;
   case 9:
      P[st][place_P]   = 1.0;
      P[st][place_P+1] = e1;
      P[st][place_P+2] = e2;
      P[st][place_P+3] = e1e2;
      P[st][place_P+4] = 1.0-3.0*e1e1;
      P[st][place_P+5] = 1.0-3.0*e2e2;
      P[st][place_P+6] = e1e1e2;
      P[st][place_P+7] = e1e2e2;
      P[st][place_P+8] = 1.0-9.0*e1e1e2e2;
      place_P+=9;
   break;
   default:
      dserror("eas: DICKENRICHTUNG: E33 LINEAR other then 0,1,3,4,5,8,9");
   }
/*----------------------------------------------------------------------
      QUERSCHUB: E13,E23 KONSTANT
  ----------------------------------------------------------------------*/
   switch(eas[3])
   {
   case 0:
   break;
   case 2:
      P[qr][place_P]   = e2-3.0*e1e1e2;
      P[qs][place_P+1] = e1-3.0*e1e2e2;
      place_P+=2;
   break;
   case 4:
      P[qr][place_P]   = 1.0-3.0*e1e1;
      P[qr][place_P+1] = e2-3.0*e1e1e2;
      P[qs][place_P+2] = 1.0-3.0*e2e2;
      P[qs][place_P+3] = e1-3.0*e1e2e2;
      place_P+=4;
   break;
   case 6:
      P[qr][place_P]   = 1.0-3.0*e1e1;
      P[qr][place_P+1] = e2-3.0*e1e1e2;
      P[qr][place_P+2] = e2e2-3.0*e1e1e2e2;
      P[qs][place_P+3] = 1.0-3.0*e2e2;
      P[qs][place_P+4] = e1-3.0*e1e2e2;
      P[qs][place_P+5] = e1e1-3.0*e1e1e2e2;
      place_P+=6;
   break;
   default:
      dserror("eas: QUERSCHUB: E13,E23 KONSTANT other then 0,2,4,6");
   }
/*----------------------------------------------------------------------
      QUERSCHUB: E13,E23 LINEAR
  ----------------------------------------------------------------------*/
   switch(eas[4])
   {
   case 0:
   break;
   case 2:
      P[sr][place_P]   = e1e1;
      P[ss][place_P+1] = e2e2;
      place_P+=2;
   break;
   case 4:
      P[sr][place_P]   = e1e1;
      P[sr][place_P+1] = e1e1e2e2;
      P[ss][place_P+2] = e2e2;
      P[ss][place_P+3] = e1e1e2e2;
      place_P+=4;
   break;
   case 6:
      P[sr][place_P]   = e1e1;
      P[sr][place_P+1] = e1e1e2;
      P[sr][place_P+2] = e1e1e2e2;
      P[ss][place_P+3] = e2e2;
      P[ss][place_P+4] = e1e2e2;
      P[ss][place_P+5] = e1e1e2e2;
      place_P+=6;
   break;
   default:
      dserror("eas: QUERSCHUB: E13,E23 LINEAR other then 0,2,4,6");
   }
   /*------------------------------------------------ check correctness */
   if (place_P != nhyb) dserror("wrong parameter count");
}
/*--------------------------------------------------- four node element */
else if (iel==4)
{
/*----------------------------------------------------------------------
      MEMBRAN: E11,E12,E22 KONSTANT
  ----------------------------------------------------------------------*/
   switch(eas[0])
   {
   case 0:
   break;
   case 1:
      P[nss][place_P]=e2;
      place_P+=1;
   break;
   case 2:
      P[nrs][place_P]  =e1;
      P[nrs][place_P+1]=e2;
      place_P+=2;
   break;
   case 3:
      P[nrs][place_P]  =e1;
      P[nrs][place_P+1]=e2;
      P[nrs][place_P+2]=e1e2;
      place_P+=3;
   break;
   case 4:
      P[nrr][place_P]  =e1;
      P[nss][place_P+1]=e2;
      P[nrs][place_P+2]=e1;
      P[nrs][place_P+3]=e2;
      place_P+=4;
   break;
   case 5:
      P[nrr][place_P]  =e1;
      P[nss][place_P+1]=e2;
      P[nrs][place_P+2]=e1;
      P[nrs][place_P+3]=e2;
      P[nrs][place_P+4]=e1e2;
      place_P+=5;
   break;
   case 7:
      P[nrr][place_P]  =e1;
      P[nss][place_P+1]=e2;
      P[nrs][place_P+2]=e1;
      P[nrs][place_P+3]=e2;
      P[nrr][place_P+4]=e1e2;
      P[nss][place_P+5]=e1e2;
      P[nrs][place_P+6]=e1e2;
      place_P+=7;
   break;
   }
/*----------------------------------------------------------------------
      BIEGUNG: E11,E12,E22 LINEAR
  ----------------------------------------------------------------------*/
   switch(eas[1])
   {
   case 0:
   break;
   case 4:
      P[mrr][place_P]  =e1;
      P[mss][place_P+1]=e2;
      P[mrs][place_P+2]=e1;
      P[mrs][place_P+3]=e2;
      place_P+=4;
   break;
   case 5:
      P[mrr][place_P]  =e1;
      P[mss][place_P+1]=e2;
      P[mrs][place_P+2]=e1;
      P[mrs][place_P+3]=e2;
      P[mrs][place_P+4]=e1e2;
      place_P+=5;
   break;
   case 7:
      P[mrr][place_P]  =e1;
      P[mss][place_P+1]=e2;
      P[mrs][place_P+2]=e1;
      P[mrs][place_P+3]=e2;
      P[mrr][place_P+4]=e1e2;
      P[mss][place_P+5]=e1e2;
      P[mrs][place_P+6]=e1e2;
      place_P+=7;
   break;
   case 6:
      P[mrr][place_P]  =e1e1;
      P[mrr][place_P+1]=e1e1e2e2;
      P[mss][place_P+2]=e2e2;
      P[mss][place_P+3]=e1e1e2e2;
      P[mrs][place_P+4]=e1e1;
      P[mrs][place_P+5]=e2e2;
      place_P+=6;
   break;
   }
/*----------------------------------------------------------------------
      DICKENRICHTUNG: E33 LINEAR (--> 7P - FORMULIERUNG)
  ----------------------------------------------------------------------*/
   switch(eas[2])
   {
   case 0:
   break;
   case 1:
      P[st][place_P] = 1.0;
      place_P+=1;
   break;
   case 3:
      P[st][place_P]   = 1.0;
      P[st][place_P+1] = e1;
      P[st][place_P+2] = e2;
      place_P+=3;
   break;
   case 4:
      P[st][place_P]   = 1.0;
      P[st][place_P+1] = e1;
      P[st][place_P+2] = e2;
      P[st][place_P+3] = e1e2;
      place_P+=4;
   break;
   case 6:
      P[st][place_P]   = 1.0;
      P[st][place_P+1] = e1;
      P[st][place_P+2] = e2;
      P[st][place_P+3] = e1e2;
      P[st][place_P+4] = e1e1;
      P[st][place_P+5] = e2e2;
      place_P+=6;
   break;
   case 8:
      P[st][place_P]   = 1.0;
      P[st][place_P+1] = e1;
      P[st][place_P+2] = e2;
      P[st][place_P+3] = e1e2;
      P[st][place_P+4] = e1e1;
      P[st][place_P+5] = e2e2;
      P[st][place_P+6] = e1e1e2;
      P[st][place_P+7] = e1e2e2;
      place_P+=8;
   break;
   case 9:
      P[st][place_P]   = 1.0;
      P[st][place_P+1] = e1;
      P[st][place_P+2] = e2;
      P[st][place_P+3] = e1e2;
      P[st][place_P+4] = 1.0-3.0*e1e1;
      P[st][place_P+5] = 1.0-3.0*e2e2;
      P[st][place_P+6] = e1e1e2;
      P[st][place_P+7] = e1e2e2;
      P[st][place_P+8] = 1.0-9.0*e1e1e2e2;
      place_P+=9;
   break;
   default:
      dserror("eas: DICKENRICHTUNG: E33 LINEAR other then 0,1,3,4,5,8,9");
   }
/*----------------------------------------------------------------------
      QUERSCHUB: E13,E23 KONSTANT
  ----------------------------------------------------------------------*/
   switch(eas[3])
   {
   case 0:
   break;
   case 2:
      P[qr][place_P]  =e1;
      P[qs][place_P+1]=e2;
      place_P+=2;
   break;
   case 4:
      P[qr][place_P]  =e1;
      P[qr][place_P+1]=e1e2;
      P[qs][place_P+2]=e2;
      P[qs][place_P+3]=e1e2;
      place_P+=4;
   break;
   }
/*----------------------------------------------------------------------
      QUERSCHUB: E13,E23 LINEAR
  ----------------------------------------------------------------------*/
   switch(eas[4])
   {
   case 0:
   break;
   case 2:
      P[sr][place_P]  =e1;
      P[ss][place_P+1]=e2;
      place_P+=2;
   break;
   case 4:
      P[sr][place_P]  =e1;
      P[sr][place_P+1]=e1e2;
      P[ss][place_P+2]=e2;
      P[ss][place_P+3]=e1e2;
      place_P+=4;
   break;
   }
}
/*------------------------------------------------------------- default */
else
{
   dserror("eas has 8,9 and 4 node elements only");
}
/*----------------------------------------------------------------------*/
if (place_P != nhyb) dserror("wrong parameter nhyb in EAS");
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8_eas */



/*----------------------------------------------------------------------*
 | transform the eas-strains from midpoint to gausspoint  m.gee 7/01    |
 *----------------------------------------------------------------------*/
void s8_transeas(DOUBLE      **P,
                    DOUBLE      **transP,
                    DOUBLE      **T,
                    DOUBLE      **akovr,
                    DOUBLE      **akonr0,
                    DOUBLE        detr,
                    DOUBLE        detr0,
                    INT           nhyb)
{
INT          i,j;
const DOUBLE two=2.0;
DOUBLE       t11,t12,t13,t21,t22,t23,t31,t32,t33;
DOUBLE       factor;
#ifdef DEBUG
dstrc_enter("s8_transeas");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------- components of the transformation matrix T */
t11=0.0;
t12=0.0;
t13=0.0;
t21=0.0;
t22=0.0;
t23=0.0;
t31=0.0;
t32=0.0;
t33=1.0;
for (i=0; i<3; i++)
{
   t11 += akovr[i][0]*akonr0[i][0];
   t12 += akovr[i][0]*akonr0[i][1];
/*   t13 += akovr[i][0]*akonr0[i][2]; */
   t21 += akovr[i][1]*akonr0[i][0];
   t22 += akovr[i][1]*akonr0[i][1];
/*   t23 += akovr[i][1]*akonr0[i][2]; */
/*   t31 += akovr[i][2]*akonr0[i][0]; */
/*   t32 += akovr[i][2]*akonr0[i][1]; */
/*   t33 += akovr[i][2]*akonr0[i][2]; */
}
factor = detr0/detr;

T[0][0] = factor*t11*t11;
T[1][0] = factor*two*t11*t21;
T[2][0] = factor*two*t11*t31;
T[3][0] = factor*t21*t21;
T[4][0] = factor*two*t21*t31;
T[5][0] = factor*t31*t31;

T[0][1] = factor*t11*t12;
T[1][1] = factor*(t11*t22+t21*t12);
T[2][1] = factor*(t11*t32+t31*t12);
T[3][1] = factor*t21*t22;
T[4][1] = factor*(t21*t32+t31*t22);
T[5][1] = factor*t31*t32;

T[0][2] = factor*t11*t13;
T[1][2] = factor*(t11*t23+t21*t13);
T[2][2] = factor*(t11*t33+t31*t13);
T[3][2] = factor*t21*t23;
T[4][2] = factor*(t21*t33+t31*t23);
T[5][2] = factor*t31*t33;

T[0][3] = factor*t12*t12;
T[1][3] = factor*two*t12*t22;
T[2][3] = factor*two*t12*t32;
T[3][3] = factor*t22*t22;
T[4][3] = factor*two*t22*t32;
T[5][3] = factor*t32*t32;

T[0][4] = factor*t12*t13          ;
T[1][4] = factor*(t12*t23+t22*t13);
T[2][4] = factor*(t12*t33+t32*t13);
T[3][4] = factor*t22*t23          ;
T[4][4] = factor*(t22*t33+t32*t23);
T[5][4] = factor*t32*t33          ;

T[0][5] = factor*t13*t13    ;
T[1][5] = factor*two*t13*t23;
T[2][5] = factor*two*t13*t33;
T[3][5] = factor*t23*t23    ;
T[4][5] = factor*two*t23*t33;
T[5][5] = factor*t33*t33    ;


T[6][6]  = factor*t11*t11     ;
T[7][6]  = factor*two*t11*t21 ;
T[8][6]  = factor*two*t11*t31 ;
T[9][6]  = factor*t21*t21    ;
T[10][6] = factor*two*t21*t31;
T[11][6] = factor*t31*t31    ;

T[6][7]  = factor*t11*t12           ;
T[7][7]  = factor*(t11*t22+t21*t12) ;
T[8][7]  = factor*(t11*t32+t31*t12) ;
T[9][7]  = factor*t21*t22          ;
T[10][7] = factor*(t21*t32+t31*t22);
T[11][7] = factor*t31*t32          ;

T[6][8]  = factor*t11*t13           ;
T[7][8]  = factor*(t11*t23+t21*t13) ;
T[8][8]  = factor*(t11*t33+t31*t13) ;
T[9][8]  = factor*t21*t23          ;
T[10][8] = factor*(t21*t33+t31*t23);
T[11][8] = factor*t31*t33          ;

T[6][9]  = factor*t12*t12     ;
T[7][9]  = factor*two*t12*t22 ;
T[8][9]  = factor*two*t12*t32 ;
T[9][9]  = factor*t22*t22    ;
T[10][9] = factor*two*t22*t32;
T[11][9] = factor*t32*t32    ;

T[6][10]  = factor*t12*t13           ;
T[7][10]  = factor*(t12*t23+t22*t13) ;
T[8][10]  = factor*(t12*t33+t32*t13) ;
T[9][10]  = factor*t22*t23          ;
T[10][10] = factor*(t22*t33+t32*t23);
T[11][10] = factor*t32*t33          ;

T[6][11]  = factor*t13*t13     ;
T[7][11]  = factor*two*t13*t23 ;
T[8][11]  = factor*two*t13*t33 ;
T[9][11]  = factor*t23*t23    ;
T[10][11] = factor*two*t23*t33;
T[11][11] = factor*t33*t33    ;

for (i=0; i<6; i++)
{
   for (j=0; j<6; j++)
   {
      T[i][j+6]=0.0;
      T[i+6][j]=0.0;
   }
}
/*--------------------------------------------------- multiply TP = T*P */
math_matmatdense(transP,T,P,12,12,nhyb,0,1.0);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8_transeas */
#endif

