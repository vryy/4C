/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1intg' which calculates the needed values
       for the integration points

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"

/*! 
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | integration points                                        al 6/01    |
 -----------------------------------------------------------------------|
 |     ngp[0] --> number of integration points r-direction              |
 |                for triangular elements :                             |
 |                number of over all integration points                 |
 |     ngp[1] --> number of integration points s-direction              |
 |                for triangular elements :                             |
 |                parameter for alternative integration                 |
 |                0 = standard    1 = gauss-radau                       |
 |     ntyp   --> 1 = triangular, 0 = rectangular  element              |
 |     xgrr   --> gauss sampling points        r-direction              |
 |     xgss   --> gauss sampling points        s-direction              |
 |     wgtr   --> weighting factors            r-direction              |
 |     wgts   --> weighting factors            s-direction              |
 *----------------------------------------------------------------------*/
void w1intg(ELEMENT   *ele,
            W1_DATA   *data,
            int        option)
{
int i, k;

double zero  = 0.0;
double one   = 1.0;
double two   = 2.0;
double three = 3.0;

double  q12, q13, q16, q23;
double  xgr[13][8],xgs[13][8],wgtt[13][8];
static double xg[6][6],wgt[6][6];
#ifdef DEBUG 
dstrc_enter("w1intg");
#endif
/*----------------------------------------------------------------------*/
if (option==0)
{                                                  /* initialize arrays */
  for (i=0; i<13; i++)
  {
    for (k=0; k<8; k++)
    {
       xgr[i][k] = 0.;
       xgs[i][k] = 0.;
      wgtt[i][k] = 0.;
    }
  }
  for (i=0; i<6; i++)
  {
    for (k=0; k<6; k++)
    {
       xg[i][k] = 0.;
       wgt[i][k] = 0.;
    }
  }
/*----------------------------------------------------------------------*  
 |     INTEGRATION PARAMETERS FOR    R E C T A N G U L A R   ELEMENTS   |        
 |     GAUSS SAMPLING POINTS  AT     R/S-COORDINATES     RESPECTIVELY   |      
 |                            AND    CORRESPONDING WEIGHTING  FACTORS   |      
 *----------------------------------------------------------------------*/       
      xg[0][1]  =  -0.5773502691896;                                            
      xg[1][1]  =  -xg[0][1]       ;
      xg[0][2]  =  -0.7745966692415;
      xg[2][2]  =  -xg[0][2]       ;
      xg[0][3]  =  -0.8611363115941;
      xg[1][3]  =  -0.3399810435849;
      xg[2][3]  =  -xg[1][3]       ;
      xg[3][3]  =  -xg[0][3]       ;
      xg[0][4]  =  -0.9061798459387;
      xg[1][4]  =  -0.5384693101057;
      xg[3][4]  =  -xg[1][4]       ;
      xg[4][4]  =  -xg[0][4]       ;
      xg[0][5]  =  -0.9324695142032;                                            
      xg[1][5]  =  -0.6612093864663;                                            
      xg[2][5]  =  -0.2386191860832;                                            
      xg[3][5]  =  -xg[2][5]       ;
      xg[4][5]  =  -xg[1][5]       ;
      xg[5][5]  =  -xg[0][5]       ;
                                 
      wgt[0][0] =  two             ;
      wgt[0][1] =  one             ;
      wgt[1][1] =  one             ;
      wgt[0][2] =  0.5555555555556 ;
      wgt[1][2] =  0.8888888888889 ;
      wgt[2][2] =  wgt[0][2]       ;
      wgt[0][3] =  0.3478548451375 ;
      wgt[1][3] =  0.6521451548625 ;
      wgt[2][3] =  wgt[1][3]       ;
      wgt[3][3] =  wgt[0][3]       ;
      wgt[0][4] =  0.2369268850562 ;
      wgt[1][4] =  0.4786286704994 ;
      wgt[2][4] =  0.5688888888889 ;
      wgt[3][4] =  wgt[1][4]       ;
      wgt[4][4] =  wgt[0][4]       ;
      wgt[0][5] =  0.1713244923792 ;
      wgt[1][5] =  0.3607615730481 ;
      wgt[2][5] =  0.4679139345727 ;
      wgt[3][5] =  wgt[2][5]       ;
      wgt[4][5] =  wgt[1][5]       ;
      wgt[5][5] =  wgt[0][5]       ;
/*----------------------------------------------------------------------*  
 |     INTEGRATION PARAMETERS FOR    T R I A N G U L A R     ELEMENTS   |
 |     GAUSS SAMPLING POINTS  AT     R/S-COORDINATES     RESPECTIVELY   |
 |                            AND    CORRESPONDING WEIGHTING  FACTORS   |
 *----------------------------------------------------------------------*/       
/*---------------------------------------------- form initial values ---*/
      q12 = one/two   ;
      q13 = one/three ;
      q16 = q13/two   ;
      q23 = q13*two   ;                                                         
/*------------------------------------------ initialize xgr,xgs,xgtt ---*/
      for(i=0;i<8 ;i++)
      {
        for(k=0;k<13;k++)
        {
         xgr[k][i]    =  zero ;
         xgs[k][i]    =  zero ;
         wgtt[k][i]   =  zero ;
         }
      }                                                       
/*----------------------------------------------------------------------*  
 |    GAUSS INTEGRATION         1 SAMPLING POINT, DEG.OF PRECISION 1    |        
 *----------------------------------------------------------------------*/       
      xgr[0][0]    =  q13 ;
      xgs[0][0]    =  q13 ;
      wgtt[0][0]   =  q12 ;
/*----------------------------------------------------------------------*  
 |    GAUSS INTEGRATION        3 SAMPLING POINTS, DEG.OF PRECISION 2    |        
 *----------------------------------------------------------------------*/       
       switch(ele->e.w1->nGP[1])/* direction s */
       {
       case 0:
          xgr[0][1]    =  q12  ;
          xgr[1][1]    =  q12  ;
          xgr[2][1]    =  zero ;
          xgs[0][1]    =  zero ;
          xgs[1][1]    =  q12  ;
          xgs[2][1]    =  q12  ;
          wgtt[0][1]   =  q16  ;
          wgtt[1][1]   =  q16  ;
          wgtt[2][1]   =  q16  ;
       break;
/*----------------------------------------------------------------------*  
 |    ALT.GAUSS INTEGRATION    3 SAMPLING POINTS, DEG.OF PRECISION 2    |        
 *----------------------------------------------------------------------*/       
       case 1:
          xgr[0][1]    =  q16  ;
          xgr[1][1]    =  q23  ;
          xgr[2][1]    =  q16  ;
          xgs[0][1]    =  q16  ;
          xgs[1][1]    =  q16  ;
          xgs[2][1]    =  q23  ;
          wgtt[0][1]   =  q16  ;
          wgtt[1][1]   =  q16  ;
          wgtt[2][1]   =  q16  ;
       break;
       }
/*----------------------------------------------------------------------*  
 |    GAUSS INTEGRATION        4 SAMPLING POINTS, DEG.OF PRECISION 3    |        
 *----------------------------------------------------------------------*/       
      xgr[0][2]    =  0.2                ;
      xgr[1][2]    =  0.6                ;
      xgr[2][2]    =  0.2                ;
      xgr[3][2]    =  q13                ;
      xgs[0][2]    =  0.2                ;
      xgs[1][2]    =  0.2                ;
      xgs[2][2]    =  0.6                ;
      xgs[3][2]    =  q13                ;
      wgtt[0][2]   =  0.2604166666667    ;
      wgtt[1][2]   =  wgtt[0][2]         ;
      wgtt[2][2]   =  wgtt[0][2]         ;
      wgtt[3][2]   = -0.28125            ;
/*----------------------------------------------------------------------*  
 |    GAUSS INTEGRATION        6 SAMPLING POINTS, DEG.OF PRECISION 4    |        
 *----------------------------------------------------------------------*/       
       switch(ele->e.w1->nGP[1])/* direction s */
       {
       case 0:
          xgr[0][3]    =  0.0915762135098   ;
          xgr[1][3]    =  0.8168475729805   ;
          xgr[2][3]    =  xgr[0][3]         ;
          xgr[3][3]    =  0.4459484909160   ;
          xgr[4][3]    =  xgr[3][3]         ;
          xgr[5][3]    =  0.1081030181681   ;
          xgs[0][3]    =  xgr[0][3]         ;
          xgs[1][3]    =  xgr[0][3]         ;
          xgs[2][3]    =  xgr[1][3]         ;
          xgs[3][3]    =  xgr[5][3]         ;
          xgs[4][3]    =  xgr[3][3]         ;
          xgs[5][3]    =  xgr[3][3]         ;
          wgtt[0][3]   =  0.0549758718277   ;
          wgtt[1][3]   =  wgtt[0][3]        ;
          wgtt[2][3]   =  wgtt[0][3]        ;
          wgtt[3][3]   =  0.1116907948390   ;
          wgtt[4][3]   =  wgtt[3][3]        ;
          wgtt[5][3]   =  wgtt[3][3]        ;                                      
       break;
/*----------------------------------------------------------------------*  
 |    ALT.GAUSS INTEGRATION    6 SAMPLING POINTS, DEG.OF PRECISION 3    |        
 *----------------------------------------------------------------------*/       
       case 1:
          xgr[0][3]    =  0.1090390090729   ;
          xgr[1][3]    =  0.2319333685530   ;
          xgr[2][3]    =  0.6590276223741   ;
          xgr[3][3]    =  xgr[2][3]         ;
          xgr[4][3]    =  xgr[1][3]         ;
          xgr[5][3]    =  xgr[0][3]         ;
          xgs[0][3]    =  xgr[1][3]         ;
          xgs[1][3]    =  xgr[0][3]         ;
          xgs[2][3]    =  xgr[0][3]         ;
          xgs[3][3]    =  xgr[1][3]         ;
          xgs[4][3]    =  xgr[2][3]         ;
          xgs[5][3]    =  xgr[2][3]         ;
          wgtt[0][3]   =  0.0833333333333   ;
          wgtt[1][3]   =  wgtt[0][3]        ;
          wgtt[2][3]   =  wgtt[0][3]        ;
          wgtt[3][3]   =  wgtt[0][3]        ;
          wgtt[4][3]   =  wgtt[0][3]        ;
          wgtt[5][3]   =  wgtt[0][3]        ;
       break;
       }
/*----------------------------------------------------------------------*  
 |    GAUSS INTEGRATION        7 SAMPLING POINTS, DEG.OF PRECISION 5    |        
 *----------------------------------------------------------------------*/       
       switch(ele->e.w1->nGP[1])/* direction s */
       {
       case 0:
          xgr[0][4]    =  0.1012865073235 ;
          xgr[1][4]    =  0.4701420641051 ;
          xgr[2][4]    =  0.7974269853531 ;
          xgr[3][4]    =  xgr[1][4]       ;
          xgr[4][4]    =  xgr[0][4]       ;
          xgr[5][4]    =  0.0597158717898 ;
          xgr[6][4]    =  q13             ;
          xgs[0][4]    =  xgr[0][4]       ;
          xgs[1][4]    =  xgr[5][4]       ;
          xgs[2][4]    =  xgr[0][4]       ;
          xgs[3][4]    =  xgr[1][4]       ;
          xgs[4][4]    =  xgr[2][4]       ;
          xgs[5][4]    =  xgr[1][4]       ;
          xgs[6][4]    =  q13             ;
          wgtt[0][4]   =  0.0629695902724 ;
          wgtt[1][4]   =  0.0661970763943 ;
          wgtt[2][4]   =  wgtt[0][4]      ;
          wgtt[3][4]   =  wgtt[1][4]      ;
          wgtt[4][4]   =  wgtt[0][4]      ;
          wgtt[5][4]   =  wgtt[1][4]      ;
          wgtt[6][4]   =  0.1125          ;
       break;                             
/*----------------------------------------------------------------------*  
 |    ALT.GAUSS INTEGRATION    7 SAMPLING POINTS, DEG.OF PRECISION 4    |        
 *----------------------------------------------------------------------*/       
       case 1:
          xgr[0][4]    =  0.2379323664724 ;
          xgr[1][4]    =  0.7367124989684 ;
          xgr[2][4]    =  xgr[1][4]       ;
          xgr[3][4]    =  xgr[0][4]       ;
          xgr[4][4]    =  0.0253551345591 ;
          xgr[5][4]    =  xgr[4][4]       ;
          xgr[6][4]    =  q13             ;
          xgs[0][4]    =  xgr[4][4]       ;
          xgs[1][4]    =  xgr[4][4]       ;
          xgs[2][4]    =  xgr[0][4]       ;
          xgs[3][4]    =  xgr[1][4]       ;
          xgs[4][4]    =  xgr[1][4]       ;
          xgs[5][4]    =  xgr[0][4]       ;
          xgs[6][4]    =  q13             ;
          wgtt[0][4]   =  0.0520833333333 ;
          wgtt[1][4]   =  wgtt[0][4]      ;
          wgtt[2][4]   =  wgtt[0][4]      ;
          wgtt[3][4]   =  wgtt[0][4]      ;
          wgtt[4][4]   =  wgtt[0][4]      ;
          wgtt[5][4]   =  wgtt[0][4]      ;
          wgtt[6][4]   =  0.1875          ;
       break;
       }
/*----------------------------------------------------------------------*  
 |    GAUSS INTEGRATION        9 SAMPLING POINTS, DEG.OF PRECISION 5    |        
 *----------------------------------------------------------------------*/       
      xgr[0][5]    =  0.1654099273898     ;
      xgr[1][5]    =  0.4375252483834     ;
      xgr[2][5]    =  0.7971126518601     ;
      xgr[3][5]    =  xgr[2][5]           ;
      xgr[4][5]    =  xgr[1][5]           ;
      xgr[5][5]    =  xgr[0][5]           ;
      xgr[6][5]    =  0.0374774207501     ;
      xgr[7][5]    =  0.1249495032332     ;
      xgr[8][5]    =  xgr[6][5]           ;
      xgs[0][5]    =  xgr[6][5]           ;
      xgs[1][5]    =  xgr[7][5]           ;
      xgs[2][5]    =  xgr[6][5]           ;
      xgs[3][5]    =  xgr[0][5]           ;
      xgs[4][5]    =  xgr[1][5]           ;
      xgs[5][5]    =  xgr[2][5]           ;
      xgs[6][5]    =  xgr[2][5]           ;
      xgs[7][5]    =  xgr[1][5]           ;
      xgs[8][5]    =  xgr[0][5]           ;
      wgtt[0][5]   =  0.0318457071431     ;
      wgtt[1][5]   =  0.1029752523804     ;
      wgtt[2][5]   =  wgtt[0][5]          ;
      wgtt[3][5]   =  wgtt[0][5]          ;
      wgtt[4][5]   =  wgtt[1][5]          ;
      wgtt[5][5]   =  wgtt[0][5]          ;
      wgtt[6][5]   =  wgtt[0][5]          ;
      wgtt[7][5]   =  wgtt[1][5]          ;
      wgtt[8][5]   =  wgtt[0][5]          ;                                       
/*----------------------------------------------------------------------*  
 |    GAUSS INTEGRATION       12 SAMPLING POINTS, DEG.OF PRECISION 6    |        
 *----------------------------------------------------------------------*/       
      xgr[0][6]    =  0.0630890144915     ;
      xgr[1][6]    =  0.3103524510338     ;
      xgr[2][6]    =  0.6365024991214     ;
      xgr[3][6]    =  0.8738219710170     ;
      xgr[4][6]    =  xgr[2][6]           ;
      xgr[5][6]    =  xgr[1][6]           ;
      xgr[6][6]    =  xgr[0][6]           ;
      xgr[7][6]    =  0.0531450498448     ;
      xgr[8][6]    =  xgr[7][6]           ;
      xgr[9][6]   =  0.2492867451709     ;
      xgr[10][6]   =  0.5014265096582     ;
      xgr[11][6]   =  xgr[9][6]          ;
      xgs[0][6]    =  xgr[0][6]           ;
      xgs[1][6]    =  xgr[7][6]           ;
      xgs[2][6]    =  xgr[7][6]           ;
      xgs[3][6]    =  xgr[0][6]           ;
      xgs[4][6]    =  xgr[1][6]           ;
      xgs[5][6]    =  xgr[2][6]           ;
      xgs[6][6]    =  xgr[3][6]           ;
      xgs[7][6]    =  xgr[2][6]           ;
      xgs[8][6]    =  xgr[1][6]           ;
      xgs[9][6]   =  xgr[9][6]          ;
      xgs[10][6]   =  xgr[9][6]          ;
      xgs[11][6]   =  xgr[10][6]          ;
      wgtt[0][6]   =  0.0254224531851     ;
      wgtt[1][6]   =  0.0414255378092     ;
      wgtt[2][6]   =  wgtt[1][6]          ;
      wgtt[3][6]   =  wgtt[0][6]          ;
      wgtt[4][6]   =  wgtt[1][6]          ;
      wgtt[5][6]   =  wgtt[1][6]          ;
      wgtt[6][6]   =  wgtt[0][6]          ;
      wgtt[7][6]   =  wgtt[1][6]          ;
      wgtt[8][6]   =  wgtt[1][6]          ;
      wgtt[9][6]  =  0.0583931378632     ;
      wgtt[10][6]  =  wgtt[9][6]         ;
      wgtt[11][6]  =  wgtt[9][6]         ;
/*----------------------------------------------------------------------*
 |    GAUSS INTEGRATION       13 SAMPLING POINTS, DEG.OF PRECISION 7    |
 *----------------------------------------------------------------------*/       
      xgr[0][7]    =  0.0651301029022     ;
      xgr[1][7]    =  0.3128654960049     ;
      xgr[2][7]    =  0.6384441885698     ;
      xgr[3][7]    =  0.8697397941956     ;
      xgr[4][7]    =  xgr[2][7]           ;
      xgr[5][7]    =  xgr[1][7]           ;
      xgr[6][7]    =  xgr[0][7]           ;
      xgr[7][7]    =  0.0486903154253     ;
      xgr[8][7]    =  xgr[7][7]           ;
      xgr[9][7]   =  0.2603459660790     ;
      xgr[10][7]   =  0.4793080678419     ;
      xgr[11][7]   =  xgr[9][7]          ;
      xgr[12][7]   =  q13                 ;
      xgs[0][7]    =  xgr[0][7]           ;
      xgs[1][7]    =  xgr[7][7]           ;
      xgs[2][7]    =  xgr[7][7]           ;
      xgs[3][7]    =  xgr[0][7]           ;
      xgs[4][7]    =  xgr[1][7]           ;
      xgs[5][7]    =  xgr[2][7]           ;
      xgs[6][7]    =  xgr[3][7]           ;
      xgs[7][7]    =  xgr[2][7]           ;
      xgs[8][7]    =  xgr[1][7]           ;
      xgs[9][7]   =  xgr[9][7]          ;
      xgs[10][7]   =  xgr[9][7]          ;
      xgs[11][7]   =  xgr[10][7]          ;
      xgs[12][7]   =  q13                 ;
      wgtt[0][7]   =  0.0266736178044     ;
      wgtt[1][7]   =  0.0385568804451     ;
      wgtt[2][7]   =  wgtt[1][7]          ;
      wgtt[3][7]   =  wgtt[0][7]          ;
      wgtt[4][7]   =  wgtt[1][7]          ;
      wgtt[5][7]   =  wgtt[1][7]          ;
      wgtt[6][7]   =  wgtt[0][7]          ;
      wgtt[7][7]   =  wgtt[1][7]          ;
      wgtt[8][7]   =  wgtt[1][7]          ;
      wgtt[9][7]  =  0.0878076287166     ;
      wgtt[10][7]  =  wgtt[9][7]         ;
      wgtt[11][7]  =  wgtt[9][7]         ;
      wgtt[12][7]  = -0.0747850222338     ;
/*----------------------------------------------------------------------*/
}                                                  /* initialize arrays */
else
{
/*-------------------------------------- rectangualar element values ---*/
      for(i=0;i<ele->e.w1->nGP[0] ;i++)
      {
        data->xgrr[i] = xg[ i][ele->e.w1->nGP[0]-1];
        data->wgtr[i] = wgt[i][ele->e.w1->nGP[0]-1];
      }                                                       

      for(i=0;i<ele->e.w1->nGP[1] ;i++)
      {
        data->xgss[i] = xg[ i][ele->e.w1->nGP[1]-1];
        data->wgts[i] = wgt[i][ele->e.w1->nGP[1]-1];
      }                                                       
/*----------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1intg */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
