/*!----------------------------------------------------------------------
\file
\brief integration parameters for fluid2 element

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID2 
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
static double Q12 = ONE/TWO;
static double Q13 = ONE/THREE;
static double Q16 = ONE/SIX;
static double Q23 = TWO/THREE;
/*!---------------------------------------------------------------------                                         
\brief integration parameters for fluid2 element

<pre>                                                         genk 06/02

this routine is a try to organise the integration parameters   
different. ALL paramters are stored in FLUID_DATA, so that this
routine has to be (hopefully) called only once!!!	       

</pre>
\param  *data 	  FLUID_DATA       (o)	   
\param   option	  int              (i)     flag (not used at the moment 
\return void                                                                       

------------------------------------------------------------------------*/
void f2_intg(FLUID_DATA         *data,
             int                option  
	    )
{
int     i, k;                          /* simply some counters          */
double  xgr[MAXTINTP][MAXTINTC];       /* coord. of integr. points QUAD */
double  xgs[MAXTINTP][MAXTINTC];       /* coord. of integr. points QUAD */
double  wgtt[MAXTINTP][MAXTINTC];      /* integr. weights QUAD          */
double  xg[MAXQINTP][MAXQINTC];        /* coord. of integr. points TRI  */
double  wgt[MAXQINTP][MAXQINTC];       /* integr. weights TRI    	*/

#ifdef DEBUG 
dstrc_enter("f2_intg");
#endif

/*----------------------------------------------------------------------*/
if (option==0)
{                                                  /* initialize arrays */
  for (i=0; i<MAXTINTP; i++)
  {
    for (k=0; k<MAXTINTC; k++)
    {
       xgr[i][k] = ZERO;
       xgs[i][k] = ZERO;
      wgtt[i][k] = ZERO;
    } /* end loop over k */
  } /* end loop over i */
  for (i=0; i<MAXQINTP; i++)
  {
    for (k=0; k<MAXQINTC; k++)
    {
       xg[i][k]  = ZERO;
       wgt[i][k] = ZERO;
    } /* end loop over k */
  } /* end loop over i */

/*----------------------------------------------------------------------*  
 |     INTEGRATION PARAMETERS FOR    R E C T A N G U L A R   ELEMENTS   |        
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
      xgr[0][0]    =  Q13 ;
      xgs[0][0]    =  Q13 ;
      wgtt[0][0]   =  Q12 ;
/*----------------------------------------------------------------------*  
 |    GAUSS INTEGRATION        3 SAMPLING POINTS, DEG.OF PRECISION 2    |
 |                             CASE 1                                   |        
 *----------------------------------------------------------------------*/       
      xgr[0][1]    =  Q12  ;
      xgr[1][1]    =  Q12  ;
      xgr[2][1]    =  ZERO ;
      xgs[0][1]    =  ZERO ;
      xgs[1][1]    =  Q12  ;
      xgs[2][1]    =  Q12  ;
      wgtt[0][1]   =  Q16  ;
      wgtt[1][1]   =  Q16  ;
      wgtt[2][1]   =  Q16  ;
/*----------------------------------------------------------------------*  
 |    ALT.GAUSS INTEGRATION    3 SAMPLING POINTS, DEG.OF PRECISION 2    |        
 |                             CASE 2                                   |    
 *----------------------------------------------------------------------*/       

      xgr[0][2]    =  Q16  ;
      xgr[1][2]    =  Q23  ;
      xgr[2][2]    =  Q16  ;
      xgs[0][2]    =  Q16  ;
      xgs[1][2]    =  Q16  ;
      xgs[2][2]    =  Q23  ;
      wgtt[0][2]   =  Q16  ;
      wgtt[1][2]   =  Q16  ;
      wgtt[2][2]   =  Q16  ;
/*----------------------------------------------------------------------*  
 |    GAUSS INTEGRATION        4 SAMPLING POINTS, DEG.OF PRECISION 3    |        
 |                             CASE 3                                   |
 *----------------------------------------------------------------------*/       
      xgr[0][3]    =  0.2                ;
      xgr[1][3]    =  0.6                ;
      xgr[2][3]    =  0.2                ;
      xgr[3][3]    =  Q13                ;
      xgs[0][3]    =  0.2                ;
      xgs[1][3]    =  0.2                ;
      xgs[2][3]    =  0.6                ;
      xgs[3][3]    =  Q13                ;
      wgtt[0][3]   =  0.2604166666667    ;
      wgtt[1][3]   =  wgtt[0][2]         ;
      wgtt[2][3]   =  wgtt[0][2]         ;
      wgtt[3][3]   = -0.28125            ;
/*----------------------------------------------------------------------*  
 |    GAUSS INTEGRATION        6 SAMPLING POINTS, DEG.OF PRECISION 4    |        
 |                             CASE 4                                   |
 *----------------------------------------------------------------------*/       
      xgr[0][4]    =  0.0915762135098	;
      xgr[1][4]    =  0.8168475729805	;
      xgr[2][4]    =  xgr[0][3] 	;
      xgr[3][4]    =  0.4459484909160	;
      xgr[4][4]    =  xgr[3][3] 	;
      xgr[5][4]    =  0.1081030181681	;
      xgs[0][4]    =  xgr[0][3] 	;
      xgs[1][4]    =  xgr[0][3] 	;
      xgs[2][4]    =  xgr[1][3] 	;
      xgs[3][4]    =  xgr[5][3] 	;
      xgs[4][4]    =  xgr[3][3] 	;
      xgs[5][4]    =  xgr[3][3] 	;
      wgtt[0][4]   =  0.0549758718277	;
      wgtt[1][4]   =  wgtt[0][3]	;
      wgtt[2][4]   =  wgtt[0][3]	;
      wgtt[3][4]   =  0.1116907948390	;
      wgtt[4][4]   =  wgtt[3][3]	;
      wgtt[5][4]   =  wgtt[3][3]	;				            
/*----------------------------------------------------------------------*  
 |    ALT.GAUSS INTEGRATION    6 SAMPLING POINTS, DEG.OF PRECISION 3    |        
 |                             CASE 5                                   |
 *----------------------------------------------------------------------*/       
      xgr[0][5]    =  0.1090390090729	;
      xgr[1][5]    =  0.2319333685530	;
      xgr[2][5]    =  0.6590276223741	;
      xgr[3][5]    =  xgr[2][3] 	;
      xgr[4][5]    =  xgr[1][3] 	;
      xgr[5][5]    =  xgr[0][3] 	;
      xgs[0][5]    =  xgr[1][3] 	;
      xgs[1][5]    =  xgr[0][3] 	;
      xgs[2][5]    =  xgr[0][3] 	;
      xgs[3][5]    =  xgr[1][3] 	;
      xgs[4][5]    =  xgr[2][3] 	;
      xgs[5][5]    =  xgr[2][3] 	;
      wgtt[0][5]   =  0.0833333333333	;
      wgtt[1][5]   =  wgtt[0][3]	;
      wgtt[2][5]   =  wgtt[0][3]	;
      wgtt[3][5]   =  wgtt[0][3]	;
      wgtt[4][5]   =  wgtt[0][3]	;
      wgtt[5][5]   =  wgtt[0][3]	;
/*----------------------------------------------------------------------*  
 |    GAUSS INTEGRATION        7 SAMPLING POINTS, DEG.OF PRECISION 5    |        
 |                             CASE 6                                   |
 *----------------------------------------------------------------------*/       
      xgr[0][6]    =  0.1012865073235 ;
      xgr[1][6]    =  0.4701420641051 ;
      xgr[2][6]    =  0.7974269853531 ;
      xgr[3][6]    =  xgr[1][4]       ;
      xgr[4][6]    =  xgr[0][4]       ;
      xgr[5][6]    =  0.0597158717898 ;
      xgr[6][6]    =  Q13	      ;
      xgs[0][6]    =  xgr[0][4]       ;
      xgs[1][6]    =  xgr[5][4]       ;
      xgs[2][6]    =  xgr[0][4]       ;
      xgs[3][6]    =  xgr[1][4]       ;
      xgs[4][6]    =  xgr[2][4]       ;
      xgs[5][6]    =  xgr[1][4]       ;
      xgs[6][6]    =  Q13	      ;
      wgtt[0][6]   =  0.0629695902724 ;
      wgtt[1][6]   =  0.0661970763943 ;
      wgtt[2][6]   =  wgtt[0][4]      ;
      wgtt[3][6]   =  wgtt[1][4]      ;
      wgtt[4][6]   =  wgtt[0][4]      ;
      wgtt[5][6]   =  wgtt[1][4]      ;
      wgtt[6][6]   =  0.1125	      ;                             
/*----------------------------------------------------------------------*  
 |    ALT.GAUSS INTEGRATION    7 SAMPLING POINTS, DEG.OF PRECISION 4    |        
 |                             CASE 7                                   |
 *----------------------------------------------------------------------*/       
      xgr[0][7]    =  0.2379323664724 ;
      xgr[1][7]    =  0.7367124989684 ;
      xgr[2][7]    =  xgr[1][4]       ;
      xgr[3][7]    =  xgr[0][4]       ;
      xgr[4][7]    =  0.0253551345591 ;
      xgr[5][7]    =  xgr[4][4]       ;
      xgr[6][7]    =  Q13	      ;
      xgs[0][7]    =  xgr[4][4]       ;
      xgs[1][7]    =  xgr[4][4]       ;
      xgs[2][7]    =  xgr[0][4]       ;
      xgs[3][7]    =  xgr[1][4]       ;
      xgs[4][7]    =  xgr[1][4]       ;
      xgs[5][7]    =  xgr[0][4]       ;
      xgs[6][7]    =  Q13	      ;
      wgtt[0][7]   =  0.0520833333333 ;
      wgtt[1][7]   =  wgtt[0][4]      ;
      wgtt[2][7]   =  wgtt[0][4]      ;
      wgtt[3][7]   =  wgtt[0][4]      ;
      wgtt[4][7]   =  wgtt[0][4]      ;
      wgtt[5][7]   =  wgtt[0][4]      ;
      wgtt[6][7]   =  0.1875	      ;
/*----------------------------------------------------------------------*  
 |    GAUSS INTEGRATION        9 SAMPLING POINTS, DEG.OF PRECISION 5    |        
 |                             CASE 8                                   |
 *----------------------------------------------------------------------*/       
      xgr[0][8]    =  0.1654099273898     ;
      xgr[1][8]    =  0.4375252483834     ;
      xgr[2][8]    =  0.7971126518601     ;
      xgr[3][8]    =  xgr[2][5]           ;
      xgr[4][8]    =  xgr[1][5]           ;
      xgr[5][8]    =  xgr[0][5]           ;
      xgr[6][8]    =  0.0374774207501     ;
      xgr[7][8]    =  0.1249495032332     ;
      xgr[8][8]    =  xgr[6][5]           ;
      xgs[0][8]    =  xgr[6][5]           ;
      xgs[1][8]    =  xgr[7][5]           ;
      xgs[2][8]    =  xgr[6][5]           ;
      xgs[3][8]    =  xgr[0][5]           ;
      xgs[4][8]    =  xgr[1][5]           ;
      xgs[5][8]    =  xgr[2][5]           ;
      xgs[6][8]    =  xgr[2][5]           ;
      xgs[7][8]    =  xgr[1][5]           ;
      xgs[8][8]    =  xgr[0][5]           ;
      wgtt[0][8]   =  0.0318457071431     ;
      wgtt[1][8]   =  0.1029752523804     ;
      wgtt[2][8]   =  wgtt[0][5]          ;
      wgtt[3][8]   =  wgtt[0][5]          ;
      wgtt[4][8]   =  wgtt[1][5]          ;
      wgtt[5][8]   =  wgtt[0][5]          ;
      wgtt[6][8]   =  wgtt[0][5]          ;
      wgtt[7][8]   =  wgtt[1][5]          ;
      wgtt[8][8]   =  wgtt[0][5]          ;                                       
/*----------------------------------------------------------------------*  
 |    GAUSS INTEGRATION       12 SAMPLING POINTS, DEG.OF PRECISION 6    |        
 |                            CASE 9                                    |
 *----------------------------------------------------------------------*/       
      xgr[0][9]    =  0.0630890144915     ;
      xgr[1][9]    =  0.3103524510338     ;
      xgr[2][9]    =  0.6365024991214     ;
      xgr[3][9]    =  0.8738219710170     ;
      xgr[4][9]    =  xgr[2][6]           ;
      xgr[5][9]    =  xgr[1][6]           ;
      xgr[6][9]    =  xgr[0][6]           ;
      xgr[7][9]    =  0.0531450498448     ;
      xgr[8][9]    =  xgr[7][6]           ;
      xgr[9][9]    =  0.2492867451709     ;
      xgr[10][9]   =  0.5014265096582     ;
      xgr[11][9]   =  xgr[9][6]           ;
      xgs[0][9]    =  xgr[0][6]           ;
      xgs[1][9]    =  xgr[7][6]           ;
      xgs[2][9]    =  xgr[7][6]           ;
      xgs[3][9]    =  xgr[0][6]           ;
      xgs[4][9]    =  xgr[1][6]           ;
      xgs[5][9]    =  xgr[2][6]           ;
      xgs[6][9]    =  xgr[3][6]           ;
      xgs[7][9]    =  xgr[2][6]           ;
      xgs[8][9]    =  xgr[1][6]           ;
      xgs[9][9]    =  xgr[9][6]           ;
      xgs[10][9]   =  xgr[9][6]           ;
      xgs[11][9]   =  xgr[10][6]          ;
      wgtt[0][9]   =  0.0254224531851     ;
      wgtt[1][9]   =  0.0414255378092     ;
      wgtt[2][9]   =  wgtt[1][6]          ;
      wgtt[3][9]   =  wgtt[0][6]          ;
      wgtt[4][9]   =  wgtt[1][6]          ;
      wgtt[5][9]   =  wgtt[1][6]          ;
      wgtt[6][9]   =  wgtt[0][6]          ;
      wgtt[7][9]   =  wgtt[1][6]          ;
      wgtt[8][9]   =  wgtt[1][6]          ;
      wgtt[9][9]   =  0.0583931378632     ;
      wgtt[10][9]  =  wgtt[9][6]          ;
      wgtt[11][6]  =  wgtt[9][6]          ;
/*----------------------------------------------------------------------*
 |    GAUSS INTEGRATION       13 SAMPLING POINTS, DEG.OF PRECISION 7    |
 |                            CASE 10                                   |
 *----------------------------------------------------------------------*/       
      xgr[0][10]    =  0.0651301029022     ;
      xgr[1][10]    =  0.3128654960049     ;
      xgr[2][10]    =  0.6384441885698     ;
      xgr[3][10]    =  0.8697397941956     ;
      xgr[4][10]    =  xgr[2][7]           ;
      xgr[5][10]    =  xgr[1][7]           ;
      xgr[6][10]    =  xgr[0][7]           ;
      xgr[7][10]    =  0.0486903154253     ;
      xgr[8][10]    =  xgr[7][7]           ;
      xgr[9][10]    =  0.2603459660790     ;
      xgr[10][10]   =  0.4793080678419     ;
      xgr[11][10]   =  xgr[9][7]           ;
      xgr[12][10]   =  Q13                 ;
      xgs[0][10]    =  xgr[0][7]           ;
      xgs[1][10]    =  xgr[7][7]           ;
      xgs[2][10]    =  xgr[7][7]           ;
      xgs[3][10]    =  xgr[0][7]           ;
      xgs[4][10]    =  xgr[1][7]           ;
      xgs[5][10]    =  xgr[2][7]           ;
      xgs[6][10]    =  xgr[3][7]           ;
      xgs[7][10]    =  xgr[2][7]           ;
      xgs[8][10]    =  xgr[1][7]           ;
      xgs[9][10]    =  xgr[9][7]           ;
      xgs[10][10]   =  xgr[9][7]           ;
      xgs[11][10]   =  xgr[10][7]          ;
      xgs[12][10]   =  Q13                 ;
      wgtt[0][10]   =  0.0266736178044     ;
      wgtt[1][10]   =  0.0385568804451     ;
      wgtt[2][10]   =  wgtt[1][7]          ;
      wgtt[3][10]   =  wgtt[0][7]          ;
      wgtt[4][10]   =  wgtt[1][7]          ;
      wgtt[5][10]   =  wgtt[1][7]          ;
      wgtt[6][10]   =  wgtt[0][7]          ;
      wgtt[7][10]   =  wgtt[1][7]          ;
      wgtt[8][10]   =  wgtt[1][7]          ;
      wgtt[9][10]   =  0.0878076287166     ;
      wgtt[10][10]  =  wgtt[9][7]          ;
      wgtt[11][10]  =  wgtt[9][7]          ;
      wgtt[12][10]  = -0.0747850222338     ;
/*----------------------------------------------------------------------*/
/*----------------------------- store integration parameters in F2_DATA */
/*----------------- RECTANGLES -----------------------------------------*/
   for (i=0;i<MAXQINTP;i++)
   {
     for (k=0;k<MAXQINTC;k++)
     {
        data->qxg[i][k]=xg[i][k];
        data->qwgt[i][k]=wgt[i][k];
     } /* end loop over k */
   } /* end loop over i */
/*----------------- TRIANGLES ------------------------------------------*/
   for (i=0;i<MAXTINTP;i++)
   {
     for (k=0;k<MAXTINTC;k++)
     {
        data->txgr[i][k]=xgr[i][k];
        data->txgs[i][k]=xgs[i][k];
	data->twgt[i][k]=wgtt[i][k];
     } /* end loop over k */
   } /* end loop over i */
} /* endif (option==0) */                                                
else
{
/* do nothing at the moment --------------------------------------------*/       
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of f2_intg */


/*!---------------------------------------------------------------------                                         
\brief local coordinates

<pre>                                                         genk 10/02

In this routine the local coordinates of nodes are determined
		     
</pre>
\param   node        int      (i)    number of actual node 
\param   irs         int      (i)    r/s identifier
\param   iel         int      (i)    number of nodes in actual element
\return  double                                                                      

------------------------------------------------------------------------*/
double f2_rsn(
	      int            node,     
	      int             irs,    
	      int             iel       
	    )
{

double c;        /*------ return value for the coord. of the gauss point */

/*-------array for a quad storing the local coord. of gauss points------*/
double quad[][2] = {
                   {ONE,ONE} ,{-ONE,ONE} ,{-ONE,-ONE},{ONE,-ONE},
		   {ZERO,ONE},{-ONE,ZERO},{ZERO,-ONE},{ONE,ZERO},
		   {ZERO,ZERO}
		   };
double tri[][2]  = {
                   {ZERO,ZERO},{ONE,ZERO},{ZERO,ONE},
		   {ONE/TWO,ZERO},{ONE/TWO,ONE/TWO},{ZERO,ONE/TWO}
		   };
/*----------------------------------------------------------------------*/

#ifdef DEBUG 
dstrc_enter("f2_rsn");
#endif


switch(iel)	/*-- switch to number of element nodes ---*/
{
case 6: case 3:
   c = tri[node][irs];
break;
case 4: case 8: case 9: 
   c = quad[node][irs]; 
break;
default:
   dserror("number of nodes 'iel' unknown");
}/* end of switch(iel) */

#ifdef DEBUG 
dstrc_exit();
#endif

return c;
} /* end of f2rsn */



#endif

/*! @} (documentation module close)*/
