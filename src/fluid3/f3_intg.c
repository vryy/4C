/*!----------------------------------------------------------------------
\file
\brief integration parameters for fluid3 element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID3 
#include "../headers/standardtypes.h"
#include "fluid3_prototypes.h"
#include "fluid3.h"
static DOUBLE Q12  = ONE/TWO;
static DOUBLE Q14  = ONE/FOUR;
static DOUBLE Q16  = ONE/SIX;
static DOUBLE Q124 = ONE/SIX/FOUR;
static DOUBLE Q430 = FOUR/FIVE/SIX;
static DOUBLE Q9120= NINE/FOUR/FIVE/SIX;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;   
/*----------------------------------------------------------------------*
 | integration points                                      genk 03/02   |
 | this routine is a try to organise the integration parameters         |
 | different. ALL paramters are saved in FLUID_DATA, so that this       |
 | routine has to be (hopefully) called only once!!!                    |
 *----------------------------------------------------------------------*/
void f3_intg(
             INT              option)
{
INT i, k;
DOUBLE  xgr[MAXTINTP][MAXTINTC],xgs[MAXTINTP][MAXTINTC],xgt[MAXTINTP][MAXTINTC];
DOUBLE  wgtt[MAXTINTP][MAXTINTC]; 
DOUBLE  xg[MAXQINTP][MAXQINTC],wgt[MAXQINTP][MAXQINTC];
DOUBLE  palpha, pbeta;
FLUID_DATA     *data; 
FLUID_DYNAMIC  *fdyn;

#ifdef DEBUG 
dstrc_enter("f3_intg");
#endif

fdyn   = alldyn[genprob.numff].fdyn;
data   = fdyn->data;

/*----------------------------------------------------------------------*/
if (option==0)
{                                                  /* initialize arrays */
  for (i=0; i<MAXTINTP; i++)
  {
    for (k=0; k<MAXTINTC; k++)
    {
       xgr[i][k] = ZERO;
       xgs[i][k] = ZERO;
       xgt[i][k] = ZERO;
      wgtt[i][k] = ZERO;
    }
  } 
  for (i=0; i<MAXQINTP; i++)
  {
    for (k=0; k<MAXQINTC; k++)
    {
       xg[i][k] = ZERO;
       wgt[i][k] = ZERO;
    }
  }
  palpha = (FIVE+THREE*sqrt(FIVE))/20.0;  
  pbeta  = (FIVE-sqrt(FIVE))/20.0;

   
/*----------------------------------------------------------------------*  
 |     INTEGRATION PARAMETERS FOR    H E X A H E D R A L    ELEMENTS    |        
 |     GAUSS SAMPLING POINTS  AT     R/S-COORDINATES     RESPECTIVELY   |      
 |                            AND    CORRESPONDING WEIGHTING  FACTORS   |      
 |     xg[i][j]                                                         |
 |    wgt[i][j]:  i+1 - actual number of gausspoint                     |
 |                j+1 - total number of gausspoints                     |
 *----------------------------------------------------------------------*/       
/* coordinates for two gauss points */
      xg[0][1]  =  -0.5773502691896;                                            
      xg[1][1]  =  -xg[0][1]       ;
/* coordinates for three gauss points */     
      xg[0][2]  =  -0.7745966692415;
      xg[2][2]  =  -xg[0][2]       ;
/* coordinates for four gauss points */      
      xg[0][3]  =  -0.8611363115941;
      xg[1][3]  =  -0.3399810435849;
      xg[2][3]  =  -xg[1][3]       ;
      xg[3][3]  =  -xg[0][3]       ;
/* coordinates for five gauss points */
      xg[0][4]  =  -0.9061798459387;
      xg[1][4]  =  -0.5384693101057;
      xg[3][4]  =  -xg[1][4]       ;
      xg[4][4]  =  -xg[0][4]       ;
/* coordinates for six gauss points */
      xg[0][5]  =  -0.9324695142032;                                            
      xg[1][5]  =  -0.6612093864663;                                            
      xg[2][5]  =  -0.2386191860832;                                            
      xg[3][5]  =  -xg[2][5]       ;
      xg[4][5]  =  -xg[1][5]       ;
      xg[5][5]  =  -xg[0][5]       ;

/* weights for one gauss points */                                 
      wgt[0][0] =  TWO             ;
/* weights for two gauss points */
      wgt[0][1] =  ONE             ;
      wgt[1][1] =  ONE             ;
/* weights for three gauss points */
      wgt[0][2] =  0.5555555555556 ;
      wgt[1][2] =  0.8888888888889 ;
      wgt[2][2] =  wgt[0][2]       ;
/* weights for four gauss points */
      wgt[0][3] =  0.3478548451375 ;
      wgt[1][3] =  0.6521451548625 ;
      wgt[2][3] =  wgt[1][3]       ;
      wgt[3][3] =  wgt[0][3]       ;
/* weights for five gauss points */
      wgt[0][4] =  0.2369268850562 ;
      wgt[1][4] =  0.4786286704994 ;
      wgt[2][4] =  0.5688888888889 ;
      wgt[3][4] =  wgt[1][4]       ;
      wgt[4][4] =  wgt[0][4]       ;
/* weights for six gauss points */
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
 |     xgr[i][j]                                                        |
 |    wgts[i][j]:  i+1 - actual number of gausspoint                    |
 |                 j+1 - number for integration case (from input)       |
 *----------------------------------------------------------------------*/      

/*----------------------------------------------------------------------*  
 |    GAUSS INTEGRATION         1 SAMPLING POINT, DEG.OF PRECISION 1    |        
 |                              CASE 0                                  |
 *----------------------------------------------------------------------*/       
      xgr[0][0]    =  Q14 ;
      xgs[0][0]    =  Q14 ;
      xgt[0][0]    =  Q14 ;
      wgtt[0][0]   =  Q16 ;       
/*----------------------------------------------------------------------*  
 |    GAUSS INTEGRATION        4 SAMPLING POINTS, DEG.OF PRECISION 2    |
 |                             CASE 1                                   |        
 *----------------------------------------------------------------------*/       
      xgr[0][1]    =    pbeta ;
      xgr[1][1]    =    palpha;
      xgr[2][1]    =    pbeta ;
      xgr[3][1]    =    pbeta ;
      xgs[0][1]    =    pbeta ;
      xgs[1][1]    =    pbeta ;
      xgs[2][1]    =    palpha;
      xgs[3][1]    =    pbeta ;
      xgt[0][1]    =    pbeta ;
      xgt[1][1]    =    pbeta ;
      xgt[2][1]    =    pbeta ;
      xgt[3][1]    =    palpha;
      wgtt[0][1]   =    Q124  ;
      wgtt[1][1]   =    Q124  ;
      wgtt[2][1]   =    Q124  ;
      wgtt[3][1]   =    Q124  ;
/*----------------------------------------------------------------------*  
 |    ALT.GAUSS INTEGRATION    4 SAMPLING POINTS, DEG.OF PRECISION 1    |
 |                             CASE 2                                   |        
 *----------------------------------------------------------------------*/       
      xgr[0][2]    =     ZERO;
      xgr[1][2]    =     ONE ;
      xgr[2][2]    =     ZERO;
      xgr[3][2]    =     ZERO;
      xgs[0][2]    =     ZERO;
      xgs[1][2]    =     ZERO;
      xgs[2][2]    =     ONE ;
      xgs[3][2]    =     ZERO;
      xgt[0][2]    =     ZERO;
      xgt[1][2]    =     ZERO;
      xgt[2][2]    =     ZERO;
      xgt[3][2]    =     ONE ;
      wgtt[0][2]   =     Q124;
      wgtt[1][2]   =     Q124;
      wgtt[2][2]   =     Q124;
      wgtt[3][2]   =     Q124;
/*----------------------------------------------------------------------*  
 |    GAUSS INTEGRATION        5 SAMPLING POINTS, DEG.OF PRECISION 3    |
 |                             CASE 3                                   |        
 *----------------------------------------------------------------------*/       
      xgr[0][3]    =     Q14  ;
      xgr[1][3]    =     Q12  ;
      xgr[2][3]    =     Q16  ;
      xgr[3][3]    =     Q16  ;
      xgr[4][3]    =     Q16  ;
      xgs[0][3]    =     Q14  ;
      xgs[1][3]    =     Q16  ;
      xgs[2][3]    =     Q16  ;
      xgs[3][3]    =     Q16  ;
      xgs[4][3]    =     Q12  ;
      xgt[0][3]    =     Q14  ;
      xgt[1][3]    =     Q16  ;
      xgt[2][3]    =     Q16  ;
      xgt[3][3]    =     Q12  ;
      xgt[4][3]    =     Q16  ;
      wgtt[0][3]   =    -Q430 ;
      wgtt[1][3]   =     Q9120;
      wgtt[2][3]   =     Q9120;
      wgtt[3][3]   =     Q9120;
      wgtt[4][3]   =     Q9120;
/*----------------------------------------------------------------------*/
/*-------------------------- save integration parameters in FLUID_DATA  */
/*----------------- HEXAEDER -------------------------------------------*/
   for (i=0;i<MAXQINTP;i++)
   {
     for (k=0;k<MAXQINTC;k++)
     {
        data->qxg[i][k]=xg[i][k];
        data->qwgt[i][k]=wgt[i][k];
     } /* end loop over k */
   } /* end loop over i */
/*----------------- TETRAEDER ------------------------------------------*/
   for (i=0;i<MAXTINTP;i++)
   {
     for (k=0;k<MAXTINTC;k++)
     {
        data->txgr[i][k]=xgr[i][k];
        data->txgs[i][k]=xgs[i][k];
        data->txgt[i][k]=xgt[i][k];
	data->twgt[i][k]=wgtt[i][k];
     } /* end loop over k */
   } /* end loop over i */
}
else
{
/* do nothing at the moment --------------------------------------------*/       
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f3_intg */


#endif

