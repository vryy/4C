/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - s9init: which initializes a shell9 element


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
\brief vector of multilayer material law

<pre>                                                            sh 10/02
This structure struct _MULTIMAT  *multimat is defined in global_control.c
and the type is in materials.h                                                  
It holds all information about the layered section data
</pre>
*----------------------------------------------------------------------*/
extern struct _MULTIMAT  *multimat;


/*!----------------------------------------------------------------------
\brief initialize the element                                      

<pre>                     m.gee 6/01              modified by    sh 10/02
This routine initializes the multilayer shell9 element. The directors at
nodal poins are calculated and if necessary a shared director is made.
The directors are stored in 'actele->e.s9->a3ref'. Additionally some
memory for different arrays is allocated ("forces", "stresses"). 
</pre>
\param  FIELD *actfield   (i/o) structure which holds all the information of this field

\warning This routine is only called once for each element
\return void                                               
\sa calling: ---; called by: shell9()   [s9_main.c]

*----------------------------------------------------------------------*/
void s9init(FIELD *actfield)
{
INT          i,j,k;
INT          size_j;
ELEMENT     *actele;
NODE        *actnode;
S9_DATA      data;

INT          num_mlay,ml;     /*number of material layers to a kinematic layer*/
INT          num_klay,kl;     /*number of kineamtic layers to this element*/
INT          sumlay;          /*total number of layers to this element*/
INT          actlay;          /*actual layers */
MULTIMAT    *actmultimat;                            /* material of actual material layer */
INT          wa_flag;         /*flag if allocation of Working Array is necessary*/
INT          numa3;
DOUBLE       a3[3];
ARRAY        collaverdir_a;
DOUBLE     **collaverdir;
DOUBLE       h2;            /*thickness of actele -> evaluation of shared director of thickness of adjacent elements differ from each other*/

/*--- for s9_cdia -----------*/
ARRAY    funct_a;     DOUBLE  *funct;    /* shape functions */  
ARRAY    deriv_a;     DOUBLE **deriv;    /* derivatives of shape functions */ 
ARRAY4D  a3r_a;       DOUBLE ***a3r;     /* a3 in reference config -> for each kinematic layer */
ARRAY    x_a;         DOUBLE **x;        /* nodal coordinates */ 
ARRAY    xjm_a;       DOUBLE **xjm;      /* jacobian */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s9init");
#endif
/*----------------------------------------------------------------------*/
a3r     = am4def("a3r"   ,&a3r_a ,3  ,MAXNOD_SHELL9,1,0,"D3");         
x       = amdef("x"      ,&x_a   ,3  ,MAXNOD_SHELL9,"DA");
xjm     = amdef("xjm_a"  ,&xjm_a ,3  ,3            ,"DA");
funct   = amdef("funct"  ,&funct_a,MAXNOD_SHELL9,1,"DV");       
deriv   = amdef("deriv"  ,&deriv_a,2,MAXNOD_SHELL9,"DA");       
/*----------------------------------------------------------------------*/
for (i=0; i<actfield->dis[0].numele; i++)
{
   actele = &(actfield->dis[0].element[i]);
   if (actele->eltyp != el_shell9) continue;
   /*---------------------------------------- init directors of element */
   s9a3(actele);
   /*--------------- allocate the space for forces -> stress resultants */
/*   am4def("forces",&(actele->e.s9->forces),1,18,MAXGAUSS,0,"D3");*/
   /*------------------------- allocate the space for physical stresses */
   am4def("stresses",&(actele->e.s9->stresses),1,6,MAXNODESTRESS_SHELL9,0,"D3");
   /*amdef("energy",&(actele->e.s9->energy),MAXGAUSS,6,"DA");*/


  /*-------- init working array for each material layer if necessary ---*/
  num_klay = actele->e.s9->num_klay;
  /*check if there is a nonlinear material that needs a working array*/
  wa_flag = 0;
  sumlay = 0;
  for (kl=0; kl<num_klay; kl++) /*loop all kinematic layers*/
  {
     num_mlay = actele->e.s9->kinlay[kl].num_mlay;
     sumlay  += num_mlay; /*count total number of layers*/
     for (ml=0; ml<num_mlay; ml++) /*loop all material layers*/
     {
       actmultimat = &(multimat[actele->e.s9->kinlay[kl].mmatID[ml]-1]);
       if (actmultimat->mattyp == m_pl_mises  ||
           actmultimat->mattyp == m_pl_dp     ||
           actmultimat->mattyp == m_pl_hoff )     wa_flag = 1;
     }/*end loop material layers*/
  }/*end loop kinematic layers*/
  
  /*--- allocate memory for wa if wa_flag == 1 ; for each layer ----------*/
  if (wa_flag == 1)
  {
     actele->e.s9->elewa = (S9_ELE_WA*)CCACALLOC(sumlay,sizeof(S9_ELE_WA));

     actlay = 0;
     for (kl=0; kl<num_klay; kl++) /*loop all kinematic layers*/
     {
        num_mlay = actele->e.s9->kinlay[kl].num_mlay;
        for (ml=0; ml<num_mlay; ml++) /*loop all material layers*/
        {
          actmultimat = &(multimat[actele->e.s9->kinlay[kl].mmatID[ml]-1]);
          /*-------------- init working array for plasticity -------------*/
          if (actmultimat->mattyp == m_pl_mises ||
              actmultimat->mattyp == m_pl_dp    ||
              actmultimat->mattyp == m_pl_hoff )
          {
             size_j = actele->e.s9->nGP[0]*actele->e.s9->nGP[1]*actele->e.s9->nGP[2];
             actele->e.s9->elewa[actlay].ipwa = (S9_IP_WA*)CCACALLOC(size_j,sizeof(S9_IP_WA));

             /* mises or dp */
             if (actmultimat->mattyp == m_pl_mises ||
                 actmultimat->mattyp == m_pl_dp)
             {
                for (k=0; k<size_j; k++)/*initialize for every gausspoint*/
                {
                  actele->e.s9->elewa[actlay].ipwa[k].epstn = 0.;
                  actele->e.s9->elewa[actlay].ipwa[k].yip   = -1;
 
                  for (j=0; j<6; j++)
                  {
                    actele->e.s9->elewa[actlay].ipwa[k].sig[j] = 0.;
                    actele->e.s9->elewa[actlay].ipwa[k].eps[j] = 0.;
                    actele->e.s9->elewa[actlay].ipwa[k].qn[j]  = 0.;
                  }
                }
             }/*end mises or dp */

            /* hoff ... */
             if (actmultimat->mattyp == m_pl_hoff)
             {
                for (k=0; k<size_j; k++)/*initialize for every gausspoint*/
                {
                  actele->e.s9->elewa[actlay].ipwa[k].dhard = 0.;
                  actele->e.s9->elewa[actlay].ipwa[k].yip   = -1;
 
                  for (j=0; j<6; j++)
                  {
                    actele->e.s9->elewa[actlay].ipwa[k].sig[j]    = 0.;
                    actele->e.s9->elewa[actlay].ipwa[k].eps[j]    = 0.;
                    actele->e.s9->elewa[actlay].ipwa[k].dkappa[j] = 0.;
                    actele->e.s9->elewa[actlay].ipwa[k].gamma[j]  = 0.;
                  }
                  for (j=0; j<9; j++)
                  {
                    actele->e.s9->elewa[actlay].ipwa[k].rkappa[j] = 0.;
                  }
                }
             }/*end hoff ... */

          }/*end plasticity */

          actlay += 1;

        }/*end loop material layers*/
    }/*end loop kinematic layers*/

    /*--- calculate equivalent element length (dia) if a plastic material model to this element (wa_flag == 1) ----------*/
    s9intg(actele,&data,0);   
    s9_cdia(actele,&data,funct,deriv,xjm,x,a3r);

  }/*end wa_flag ==1 */

}
/* --------------------- now do modification of directors bischoff style */
 /*------------------------------ allocate space for directors at a node */
 collaverdir = amdef("averdir",&collaverdir_a,3,MAXELE,"DA");                                
 /*----------------------------------------------------------------------*/                  
 for (i=0; i<actfield->dis[0].numnp; i++)                                                    
 {                                                                                           
    actnode = &(actfield->dis[0].node[i]);                                                   
    numa3=0;                                                                                 
    for (j=0; j<actnode->numele; j++)                                                        
    {                                                                                        
       actele = actnode->element[j];                                                         
       h2 = actele->e.s9->thick;                                                             
       if (actele->eltyp != el_shell9) continue;                                             
       for (k=0; k<actele->numnp; k++)                                                       
       {                                                                                     
          if (actele->node[k] == actnode)                                                    
          {                                                                                  
             collaverdir[0][numa3] = actele->e.s9->a3ref.a.da[0][k] * h2;                    
             collaverdir[1][numa3] = actele->e.s9->a3ref.a.da[1][k] * h2;                    
             collaverdir[2][numa3] = actele->e.s9->a3ref.a.da[2][k] * h2;                    
             numa3++;                                                                        
             if (numa3 > MAXELE) dserror("Too many elements to a node, MAXELE too small");  
             break;                                                                          
          }                                                                                  
       }                                                                                     
    }                                                                                        
 /*----------- do nothing if there is only one shell9 element to actnode */                  
    if (numa3 <= 1) continue;                                                                
 /*------------------------------------------------ make shared director */                  
    s9averdir(collaverdir,numa3,a3,h2);                                                      
 /*---------------------------------put shared director back to elements */                  
    for (j=0; j<actnode->numele; j++)                                                        
    {                                                                                        
       actele = actnode->element[j];                                                         
       h2 = actele->e.s9->thick;                                                             
       if (actele->eltyp != el_shell9) continue;                                             
       for (k=0; k<actele->numnp; k++)                                                       
       {                                                                                     
          if (actele->node[k] == actnode)                                                    
          {                                                                                  
             actele->e.s9->a3ref.a.da[0][k] = a3[0]/h2;                                        
             actele->e.s9->a3ref.a.da[1][k] = a3[1]/h2;                                         
             actele->e.s9->a3ref.a.da[2][k] = a3[2]/h2;                                         
/*           make normation of director here instead of doing this in s9_a3.c !!!*/
/*             actele->e.s9->a3ref.a.da[0][k] = a3[0]; */                                       
/*             actele->e.s9->a3ref.a.da[1][k] = a3[1]; */                                        
/*             actele->e.s9->a3ref.a.da[2][k] = a3[2]; */                                        
             break;                                                                          
          }                                                                                  
       }                                                                                     
    }                                                                                        
 }/* end of loop over all nodes */                                                           
/*----------------------------------------------------------------------*/                   
amdel(&collaverdir_a);
amdel(&funct_a);       
amdel(&deriv_a);       
am4del(&a3r_a);
amdel(&x_a);
amdel(&xjm_a);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9init */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
