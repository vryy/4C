/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - s9init: which initializes a shell9 element

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/

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
int          i,j,k;
ELEMENT     *actele;
NODE        *actnode;
S9_DATA      data;
/*double     **a3ref;*/
int          numele;
int          numa3;
double       a3[3];
ARRAY        collaverdir_a;
double     **collaverdir;
double       h2;            /*thickness of actele -> evaluation of shared director of thickness of adjacent elements differ from each other*/

#ifdef DEBUG 
dstrc_enter("s9init");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<actfield->dis[0].numele; i++)
{
   actele = &(actfield->dis[0].element[i]);
   if (actele->eltyp != el_shell9) continue;
   /*---------------------------------------- init directors of element */
   s9a3(actele);
   /*--------------- allocate the space for forces -> stress resultants */
/*   am4def("forces",&(actele->e.s9->forces),1,18,MAXGAUSS,0,"D3");
   /*------------------------- allocate the space for physical stresses */
   am4def("stresses",&(actele->e.s9->stresses),1,6,MAXGAUSS,0,"D3");
   /*amdef("energy",&(actele->e.s9->energy),MAXGAUSS,6,"DA");*/
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
             if (numa3 == MAXELE) dserror("Too many elements to a node, MAXELE too small");  
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
       if (actele->eltyp != el_shell9) continue;                                             
       for (k=0; k<actele->numnp; k++)                                                       
       {                                                                                     
          if (actele->node[k] == actnode)                                                    
          {                                                                                  
             actele->e.s9->a3ref.a.da[0][k] = a3[0];                                         
             actele->e.s9->a3ref.a.da[1][k] = a3[1];                                         
             actele->e.s9->a3ref.a.da[2][k] = a3[2];                                         
             break;                                                                          
          }                                                                                  
       }                                                                                     
    }                                                                                        
 }/* end of loop over all nodes */                                                           
/*----------------------------------------------------------------------*/                   
amdel(&collaverdir_a);
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9init */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
