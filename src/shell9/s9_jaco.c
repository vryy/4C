/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - s9jaco: which calculates the jacobian matrix and the jacobi determinant. 
           This Routine is needed to calculate the shell shifter for 
           evaluating the loads applied on the shell (-> especially for 
           surface loads!)


<pre>
Maintainer: Stefan Hartmann
            hartmann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hartmann/
            0771 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculate jacobian matrix                                      

<pre>                     m.gee 10/01             modified by    sh 10/02
This routine calculates the jacobian matrix and the jacobi determinant,
which is needed to calculate the shell shifter for evaluating load vectors
if the load is applied on the surface of the shell.
</pre>
\param  DOUBLE   *funct   (i) shape functions at GP
\param  DOUBLE  **deriv   (i) shape function derivatives at GP 
\param  DOUBLE  **x       (i) coordinates of nodal points (global coordinate system)
\param  DOUBLE  **xjm     (o) jacobian matrix
\param  DOUBLE   *hte     (i) thickness of shell at nodal points
\param  DOUBLE ***a3r     (i) normed (not shared) director in ref/cur config. (->s9a3ref_extern)
\param  DOUBLE    e3      (i) thickness direction of shell (1.0=top; 0.0=mit; -1.0=bot)
\param  INT       iel     (i) number of nodes to this element
\param  DOUBLE   *deta    (o) Determinant of jacobian
\param  INT       init    (i) flag for initializing arrays
\param  INT       num_klay(i) number of kinematic layers to this element
\param  DOUBLE   *klayhgt (i) hight of kin layer in % of total thickness of shell 

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9eleload()   [s9_load1.c]
                             shell9()      [s9_main.c]

*----------------------------------------------------------------------*/
void s9jaco(DOUBLE    *funct,
            DOUBLE   **deriv,
            DOUBLE   **x,
            DOUBLE   **xjm,
            DOUBLE    *hte,
            DOUBLE  ***a3r,
            DOUBLE     e3,
            INT        iel,
            DOUBLE    *deta,
            INT        init,
            INT        num_klay,
            DOUBLE    *klayhgt)     /* hight of kin layer in % of total thickness of shell */
{
DOUBLE          x1r;
DOUBLE          x2r;
DOUBLE          x3r;
DOUBLE          x1s;
DOUBLE          x2s;
DOUBLE          x3s;

DOUBLE          hgt;                               /* element thickness */
DOUBLE          det_dummy;
INT             mlay,klay,num_mlay;
DOUBLE          mlayhgt_dummy[MAXKLAY_SHELL9];
INT             mod;

static ARRAY    gkov_a;  static DOUBLE **gkov;
static ARRAY    gkon_a;  static DOUBLE **gkon;
static ARRAY    gmkov_a; static DOUBLE **gmkov;
static ARRAY    gmkon_a; static DOUBLE **gmkon;

/* mid surface basis vectors and metric tensors -> help for s9_tvmr.c */
static ARRAY        akovh_a;     static DOUBLE **akovh;     
static ARRAY        akonh_a;     static DOUBLE **akonh;     
static ARRAY        amkovh_a;    static DOUBLE **amkovh;    
static ARRAY        amkonh_a;    static DOUBLE **amkonh;    

static ARRAY4D  akov_a;  static DOUBLE ***akov;    /* kovariant basis vectors at Int point ref.config. */
static ARRAY4D  akon_a;  static DOUBLE ***akon;    /* kontravar.--------------"----------- ref.config. */
static ARRAY4D  amkov_a; static DOUBLE ***amkov;   /* kovaraiant metric tensor at Int point ref.config. */
static ARRAY4D  amkon_a; static DOUBLE ***amkon;   /* kontravar.--------------"------------ ref.config. */
static ARRAY4D  a3kvp_a; static DOUBLE ***a3kvp;   /* partiel derivatives of normal vector ref.config. */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s9jaco");
#endif
/*----------------------------------------------------------------------*/
/*---------------------------------------------------------- init phase */
if (init==1)
{
   gkov     = amdef("gkov"  ,&gkov_a,3,3,"DA");         
   gkon     = amdef("gkon"  ,&gkon_a,3,3,"DA");         
   gmkov    = amdef("gmkov" ,&gmkov_a,3,3,"DA");         
   gmkon    = amdef("gmkon" ,&gmkon_a,3,3,"DA");         
   akov     = am4def("akov" ,&akov_a,3,3,MAXKLAY_SHELL9,0,"D3");      
   akon     = am4def("akon" ,&akon_a,3,3,MAXKLAY_SHELL9,0,"D3");       
   amkov    = am4def("amkov",&amkov_a,3,3,MAXKLAY_SHELL9,0,"D3");      
   amkon    = am4def("amkon",&amkon_a,3,3,MAXKLAY_SHELL9,0,"D3"); 
   a3kvp    = am4def("a3kvp",&a3kvp_a,3,2,MAXKLAY_SHELL9,0,"D3");              

   akovh     = amdef("akovh"  ,&akovh_a,3,3,"DA");      
   akonh     = amdef("akonh"  ,&akonh_a,3,3,"DA");      
   amkovh    = amdef("amkovh" ,&amkovh_a,3,3,"DA");     
   amkonh    = amdef("amkonh" ,&amkonh_a,3,3,"DA");     

   goto end;
}
/*-------------------------------------------------------- uninit phase */
else if (init==-1)
{
   amdel(&gkov_a);   
   amdel(&gkon_a);   
   amdel(&gmkov_a);   
   amdel(&gmkon_a);   
   am4del(&akov_a); 
   am4del(&akon_a); 
   am4del(&amkov_a);
   am4del(&amkon_a);
   am4del(&a3kvp_a);

   amdel(&akovh_a);  
   amdel(&akonh_a);  
   amdel(&amkovh_a); 
   amdel(&amkonh_a); 

   goto end;
}
/*------------------------------------------------ calculate gkov */
s9_tvmr(x,a3r,akov,akon,amkov,amkon,akovh,akonh,amkovh,amkonh,
        &det_dummy,funct,deriv,iel,a3kvp,num_klay);
/*----- make height at gaussian point -> variable thickness within the element is normally not allowed */
s9_xint(&hgt,hte,funct,iel);
/* check if top,bot or mid */
if (e3 == -1.0) /*bottom surface of shell body = bottom of first kinematic layer*/
{ 
   klay       = 0;
}
else if (e3 == 0.0) /* reference layer: NOTE: does not have to be the geometrical middle in thickness direction! */
{
   mod = num_klay % 2;    /*divide modulo 2*/
   if (mod == 0)
   {
      klay = (num_klay)/2; /* the one above the reference layer */
      e3   = -1.0;
   }
   else                    /* uneven number of kinematic layers */ 
   {
      klay = (num_klay - 1)/2; /*the middle kinematic layer*/
   }
}
else if (e3 == 1.0) /*top surface of shell body = top of last kinematic layer*/
{
   klay       = num_klay-1;
}
else dserror("wrong value for e3 in s9_jaco!");

mlayhgt_dummy[klay] = 100.0;
num_mlay   = 1;
mlay       = 0;

s9_tmtr(e3,gkov,gkon,gmkov,gmkon,&det_dummy,akov,a3kvp,hgt,
        klayhgt,mlayhgt_dummy,num_klay,klay,mlay,1.0);
/*--------------------------------------------------------------- */
xjm[0][0]=gkov[0][0];
xjm[0][1]=gkov[1][0];
xjm[0][2]=gkov[2][0];
xjm[1][0]=gkov[0][1];
xjm[1][1]=gkov[1][1];
xjm[1][2]=gkov[2][1];
xjm[2][0]=gkov[0][2];
xjm[2][1]=gkov[1][2];
xjm[2][2]=gkov[2][2];

x1r=xjm[0][0];
x2r=xjm[0][1];
x3r=xjm[0][2];
x1s=xjm[1][0];
x2s=xjm[1][1];
x3s=xjm[1][2];

*deta = DSQR(x1r*x2s - x2r*x1s) + DSQR(x3r*x1s - x3s*x1r) + DSQR(x2r*x3s - x3r*x2s);
*deta = sqrt(*deta);

if (*deta <= EPS14) dserror("Element Area equal 0.0 or negativ detected");
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s9jaco */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
