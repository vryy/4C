/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - s9_tvkg: which calculates the geometric stiffness matrix for a shell9 
            element

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief geometric striffness matrix                                      

<pre>                     m.gee 6/01              modified by    sh 10/02
This routine performs the geometric striffness matrix for a shell9 element
</pre>
\param  ARRAY     *estif     (i/o) element stiffness matrix to be modified
\param  double    *stress_r   (i)  stress resultants of actual step
\param  double    *funct      (i)  shape functions at GP
\param  double   **deriv      (i)  shape function derivatives at GP
\param  int        numdf      (i)  number of dofs to one node
\param  int        iel        (i)  number of nodes to this element
\param  double     weight     (i)  weight at GP
\param  int        num_klay   (i)  number of mat layers to this kin layer 
\param  int        klay       (i)  actual kin layer 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9static_keug()   [s9_static_keug.c]

*----------------------------------------------------------------------*/
void s9_tvkg(double **estif,
             double  *stress_r,
             double  *funct,
             double **deriv,
             int      numdf,
             int      iel,
             double   weight,
             int      klay,        /* actual kin layer */
             int      num_klay)    /* number of kin layers to this element */  
{
int     i,inode,jnode;
int     idof,jdof;
int     jlay;
double  fac1;

double  pi;
double  pj;
double  d11;
double  d12;
double  d21;
double  d22;

double  pd1ij;
double  pd1ji;
double  pd2ij;
double  pd2ji;

double  xn;
double  xm;
double  yu;
double  yo;
double  y1,y2;
double  y;
double  z;

double  sn11;
double  sn21;
double  sn31;
double  sn22;
double  sn32;
double  sn33;
double  sm11;
double  sm21;
double  sm31;
double  sm22;
double  sm32;

#ifdef DEBUG 
dstrc_enter("s9_tvkg");
#endif
/*----------------------------------------------------------------------*/
sn11 = stress_r[0];
sn21 = stress_r[1];
sn22 = stress_r[2];
sn31 = stress_r[3];
sn32 = stress_r[4];
sn33 = stress_r[5];
sm11 = stress_r[6];
sm21 = stress_r[7];
sm22 = stress_r[8];
sm31 = stress_r[9];
sm32 = stress_r[10];
/*----------------------------------------------------------------------*/
idof = 3 + (klay) * 3;
/*----------------------------------------------------------------------*/
/*---------- loop over all kinematic layers ---------------------------*/
for (jlay=0; jlay<num_klay; jlay++)
{
  /* is jlay on trajectory of reference layer to klay (=ilay) ? */
  fac1 = s9notr(num_klay,klay,jlay);
  if (fac1 == 0.0) continue;
  jdof = 3 + (jlay) * 3;

  for (inode=0; inode<iel; inode++)
  {
     for (jnode=0; jnode<=inode; jnode++)
     {
        pi = funct[inode];
        pj = funct[jnode];
        
        d11 = deriv[0][inode] * deriv[0][jnode];
        d12 = deriv[0][inode] * deriv[1][jnode];
        d21 = deriv[1][inode] * deriv[0][jnode];
        d22 = deriv[1][inode] * deriv[1][jnode];
        
        pd1ij = deriv[0][inode] * pj;
        pd1ji = deriv[0][jnode] * pi;
        pd2ij = deriv[1][inode] * pj;
        pd2ji = deriv[1][jnode] * pi;
        
        /*-------------------------------------E11,E12,E22(CONST)*/
        xn = (sn11*d11 + sn21*(d12+d21) + sn22*d22) * weight;
        /*---------------------------------------E11,E12,E22(LIN)*/
        xm = (sm11*d11 + sm21*(d12+d21) + sm22*d22) * weight;
        /*-----------------------------------------E13,E23(CONST)*/
        yu = (sn31*pd1ji + sn32*pd2ji) * weight;
        yo = (sn31*pd1ij + sn32*pd2ij) * weight;
        /*-------------------------------------------E13,E23(LIN)*/
        y1 = (sm31*pd1ij + sm32*pd2ij) * weight;
        y2 = (sm31*pd1ji + sm32*pd2ji) * weight;
        y  = y1 + y2;
        /*---------------------------------------------E33(CONST)*/
        z  = pi*pj*sn33*weight;
        
        if (klay == jlay)
        {
           estif[inode*numdf+0][jnode*numdf+0] += xn;
           estif[inode*numdf+1][jnode*numdf+1] += xn;
           estif[inode*numdf+2][jnode*numdf+2] += xn;
 
           estif[inode*numdf+idof+0][jnode*numdf+0] += (xm+yu);
           estif[inode*numdf+idof+1][jnode*numdf+1] += (xm+yu);
           estif[inode*numdf+idof+2][jnode*numdf+2] += (xm+yu);
 
           estif[inode*numdf+0][jnode*numdf+idof+0] += (xm+yo);
           estif[inode*numdf+1][jnode*numdf+idof+1] += (xm+yo);
           estif[inode*numdf+2][jnode*numdf+idof+2] += (xm+yo);
 
           estif[inode*numdf+idof+0][jnode*numdf+idof+0] += (y+z);
           estif[inode*numdf+idof+1][jnode*numdf+idof+1] += (y+z);
           estif[inode*numdf+idof+2][jnode*numdf+idof+2] += (y+z);
           
           /*make symmetric*/
           if (inode!=jnode)
           {
             estif[jnode*numdf+0][inode*numdf+0] += xn;
             estif[jnode*numdf+1][inode*numdf+1] += xn;
             estif[jnode*numdf+2][inode*numdf+2] += xn;
 
             estif[jnode*numdf+0][inode*numdf+idof+0] += (xm+yu);
             estif[jnode*numdf+1][inode*numdf+idof+1] += (xm+yu);
             estif[jnode*numdf+2][inode*numdf+idof+2] += (xm+yu);
 
             estif[jnode*numdf+idof+0][inode*numdf+0] += (xm+yo);
             estif[jnode*numdf+idof+1][inode*numdf+1] += (xm+yo);
             estif[jnode*numdf+idof+2][inode*numdf+2] += (xm+yo);
 
             estif[jnode*numdf+idof+0][inode*numdf+idof+0] += (y+z);
             estif[jnode*numdf+idof+1][inode*numdf+idof+1] += (y+z);
             estif[jnode*numdf+idof+2][inode*numdf+idof+2] += (y+z);
           }
        }
        else if (klay != jlay)
        {
           estif[jnode*numdf+jdof+0][inode*numdf+0] += xn;
           estif[jnode*numdf+jdof+1][inode*numdf+1] += xn;
           estif[jnode*numdf+jdof+2][inode*numdf+2] += xn;
        
           estif[jnode*numdf+0][inode*numdf+jdof+0] += xn;
           estif[jnode*numdf+1][inode*numdf+jdof+1] += xn;
           estif[jnode*numdf+2][inode*numdf+jdof+2] += xn;

           estif[jnode*numdf+jdof+0][inode*numdf+idof+0] += yo;
           estif[jnode*numdf+jdof+1][inode*numdf+idof+1] += yo;
           estif[jnode*numdf+jdof+2][inode*numdf+idof+2] += yo;

           estif[jnode*numdf+idof+0][inode*numdf+jdof+0] += yu;
           estif[jnode*numdf+idof+1][inode*numdf+jdof+1] += yu;
           estif[jnode*numdf+idof+2][inode*numdf+jdof+2] += yu;

           /*make symmetric*/
           if (inode!=jnode)
           {
             estif[inode*numdf+0][jnode*numdf+jdof+0] += xn;
             estif[inode*numdf+1][jnode*numdf+jdof+1] += xn;
             estif[inode*numdf+2][jnode*numdf+jdof+2] += xn;
        
             estif[inode*numdf+jdof+0][jnode*numdf+0] += xn;
             estif[inode*numdf+jdof+1][jnode*numdf+1] += xn;
             estif[inode*numdf+jdof+2][jnode*numdf+2] += xn;

             estif[inode*numdf+idof+0][jnode*numdf+jdof+0] += yo;
             estif[inode*numdf+idof+1][jnode*numdf+jdof+1] += yo;
             estif[inode*numdf+idof+2][jnode*numdf+jdof+2] += yo;

             estif[inode*numdf+jdof+0][jnode*numdf+idof+0] += yu;
             estif[inode*numdf+jdof+1][jnode*numdf+idof+1] += yu;
             estif[inode*numdf+jdof+2][jnode*numdf+idof+2] += yu;
           }
        }
        
     } /* end loop over jnode */
  } /* end loop over inode */

} /* end of loop over all kinematic layers*/

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_tvkg */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
