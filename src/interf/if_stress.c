/*!-----------------------------------------------------------------------
\file
\brief contains the routine 'if_stress' which calculates the stresses of the
interface element in the lokal system
<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>

*-----------------------------------------------------------------------*/
#ifdef D_INTERF
#include "../headers/standardtypes.h"
#include "interf.h"
#include "interf_prototypes.h"

/*! 
\addtogroup INTERF
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculates stresses in the lokal system  

<pre>                                                              ah 05/03
This routine calculates the stresses in the lokal system
</pre>

\param   *ele          ELEMENT       (I)   actual element (macro)
\param   *data         INTERF_DATA   (I)   integration point data
\param   *mat          MATERIAL      (I)   actual material
\param    init         INT           (I)   flag

\warning There is nothing special to this routine
\return void                                               
\sa calling:  ---; 
    caled by: interf();

*----------------------------------------------------------------------*/


void if_stress(ELEMENT       *ele, 
               INTERF_DATA   *data, 
               MATERIAL      *mat,
               INT            init)
{
const DOUBLE    q12 = ONE/TWO;
INT             iel;         /* numnp to this element     */
INT             i,k,lr;      /* some loopers     */
DOUBLE          det;         /* determinant of jacobian matrix  */
DOUBLE         e1;
INT             nir;         /* number of gaussian points           */
DOUBLE          help;
DOUBLE          c_parabel=0.0;
DOUBLE          b_parabel=0.0;
DOUBLE          alpha;       /* angle between xi-axis and X-axis in [0,..2 pi] */
DOUBLE          gamma;       /* angle between X-axis and xi-axis in [0,..2 pi] */
DOUBLE          co,si;       /* cosinus alpha and sin alpha  */
INT             ip;

static ARRAY    xrefe_a;     /* coordinates of element nodes */     
static DOUBLE **xrefe;         
static ARRAY    funct_a;     /* shape functions */    
static DOUBLE  *funct;     

DOUBLE T[3];                 /* stress */
DOUBLE L[2];                 /* lengh of element edges */
DOUBLE sig[3];               /* global stresses */
DOUBLE x_mid[3];
DOUBLE y_mid[3];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("if_stress");
#endif
/*------------------------------------------------- some working arrays */
if (init==1)
{
  xrefe     = amdef("xrefe"  ,&xrefe_a, 2,8, "DA");
  funct     = amdef("funct"  ,&funct_a ,3,1, "DV");       
  goto end;
}
/*----------------------------------------------------------------------*/
/* uninit phase        (init=-1)                                        */
/*----------------------------------------------------------------------*/
else if (init==-1)
{
   amdel(&xrefe_a);   
   amdel(&funct_a);
   goto end;  
}
/*----------------------------------------------------------------------*/
/*----------- check orientation of element (which is my xi direction)---*/
iel     = ele->numnp;
for (k=0; k<iel; k++)
{
 xrefe[0][k] = ele->node[k]->x[0];          /* coordinates in x-direction */
 xrefe[1][k] = ele->node[k]->x[1];          /* coordinates in y-direction */              
}
L[0] = sqrt( (xrefe[0][1] - xrefe[0][0]) * (xrefe[0][1] - xrefe[0][0])
     +       (xrefe[1][1] - xrefe[1][0]) * (xrefe[1][1] - xrefe[1][0]));
L[1] = sqrt( (xrefe[0][2] - xrefe[0][1]) * (xrefe[0][2] - xrefe[0][1])
     +       (xrefe[1][2] - xrefe[1][1]) * (xrefe[1][2] - xrefe[1][1]));
/*--------------------------------------------- integration parameters---*/
ifintg(ele,data);

switch(ele->distyp)
{
case quad4:

     if(L[0]>L[1])
     {
      x_mid[0]   = q12*(xrefe[0][0] + xrefe[0][3]);
      y_mid[0]   = q12*(xrefe[1][0] + xrefe[1][3]);
      x_mid[1]   = q12*(xrefe[0][1] + xrefe[0][2]);
      y_mid[1]   = q12*(xrefe[1][1] + xrefe[1][2]);
     }
     else if (L[1]>L[0])
     {
      x_mid[0]   = q12*(xrefe[0][0] + xrefe[0][1]);
      y_mid[0]   = q12*(xrefe[1][0] + xrefe[1][1]);
      x_mid[1]   = q12*(xrefe[0][2] + xrefe[0][3]);
      y_mid[1]   = q12*(xrefe[1][2] + xrefe[1][3]);
     }
     
break;
/*-----------------------------------------------------------------------*/
case quad8:

     if(L[0]>L[1])
     {
      x_mid[0]   = q12*(xrefe[0][0] + xrefe[0][3]);
      y_mid[0]   = q12*(xrefe[1][0] + xrefe[1][3]);
      x_mid[1]   = q12*(xrefe[0][1] + xrefe[0][2]);
      y_mid[1]   = q12*(xrefe[1][1] + xrefe[1][2]);
      x_mid[2]   = q12*(xrefe[0][4] + xrefe[0][6]);
      y_mid[2]   = q12*(xrefe[1][4] + xrefe[1][6]);
     }
     else if (L[1]>L[0])
     {
      x_mid[0]   = q12*(xrefe[0][0] + xrefe[0][1]);
      y_mid[0]   = q12*(xrefe[1][0] + xrefe[1][1]);
      x_mid[1]   = q12*(xrefe[0][2] + xrefe[0][3]);
      y_mid[1]   = q12*(xrefe[1][2] + xrefe[1][3]);
      x_mid[2]   = q12*(xrefe[0][5] + xrefe[0][7]);
      y_mid[2]   = q12*(xrefe[1][5] + xrefe[1][7]);
     }
     help      = (x_mid[0]-x_mid[1])/(x_mid[0]-x_mid[2]);
     c_parabel = (y_mid[0]-y_mid[1]-(y_mid[0]-y_mid[2])*help)/
                 (x_mid[0]*x_mid[0]-x_mid[1]*x_mid[1]-
                  (x_mid[0]*x_mid[0]-x_mid[2]*x_mid[2])*help);
     b_parabel = (y_mid[0]-y_mid[1]-c_parabel*(x_mid[0]*x_mid[0]-x_mid[1]*x_mid[1]))/
                 (x_mid[0]-x_mid[1]);
     
break;
default:
   dserror("discretisation unknown for Interface");
break;
}
nir   = ele->e.interf->nGP;

/*========================================== loop over gauss points ===*/
ip = -1;
for (lr=0; lr<nir; lr++)
{   
   ip++;
   /*---------------------- get local stresses(calculated in update) ---*/
   T[0]   = ele->e.interf->elewa[0].ipwa[ip].Tt;
   T[1]   = ele->e.interf->elewa[0].ipwa[ip].Tn;
   /*-dummy->equal number of stress compon. for glob. and loc stresses -*/
   T[2]   = 0.0;
   /*--------------- check if local or global stresses are asked for ---*/
   switch(ele->e.interf->stresstyp)
   {
   /*---------------------------------- local stresses are asked for ---*/
     case if_tn:
       for (i=0; i<3; i++)
       {
         ele->e.interf->stress_GP.a.d3[0][i][ip]= T[i];
       }
     break;
   /*--------------------------- transformation into global stresses ---*/
     case if_xy:
     
      e1    = data->xgr[lr];
      if_funcderiv(e1,ele->distyp,x_mid,y_mid,b_parabel,c_parabel,funct,&co,&si,&det);
      alpha  = acos(co);
      gamma  = TWO*PI - alpha;
      sig[0] = (T[1]/TWO) * (ONE - cos(2*gamma)) + T[0] * sin(2*gamma); 
      sig[1] = (T[1]/TWO) * (ONE + cos(2*gamma)) - T[0] * sin(2*gamma); 
      sig[2] = (T[1]/TWO) * sin(2*gamma)         + T[0] * cos(2*gamma); 
      for (i=0; i<3; i++)
      {
        ele->e.interf->stress_GP.a.d3[0][i][ip]= sig[i];
      }
     break;
     default:
     break;
   }    
}/*================================================ end of loop over GP */
/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of if_stress */

/*! @} (documentation module close)*/

/*----------------------------------------------------------------------*/
#endif /*D_INTERF*/
