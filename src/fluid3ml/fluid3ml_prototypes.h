/*!----------------------------------------------------------------------
\file
\brief fluid3ml prototypes

------------------------------------------------------------------------*/
/* RULE HOW TO ADD NEW FILES AND FUNCTIONS: 
   1.) THE FILENAMES ARE IN ALPHABETICAL ORDER !!!
   2.) FUNCTIONS ARE IN THE SAME ORDER LIKE IN THE FILE!!!
*/
/************************************************************************
 | f3_mlbubmat.c                                                          |
 ************************************************************************/
/*!---------------------------------------------------------------------                                         
\brief evaluate bubble part of matrix Kvv for fluid3

<pre>                                                       gravem 07/03

In this routine, the bubble part of the matrix Kvv is calculated.

NOTE: there's only one elestif
      --> Kvv is stored in estif[0..(3*iel-1)][0..(3*iel-1)]
      
</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **estif     DOUBLE	   (i/o)  element stiffness matrix
\param  *velint    DOUBLE	   (i)    velocity at int point
\param **vderxy    DOUBLE	   (i)    global velocity derivatives
\param  *funct     DOUBLE	   (i)    natural shape functions
\param **derxy     DOUBLE	   (i)    global deriv. of shape fun.
\param  *vbubint   DOUBLE	   (i)    velocity bubble functions
\param **vbubderxy DOUBLE	   (i)    global deriv. of vel. bub. fun.
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   visc      DOUBLE	   (i)    fluid viscosity	     
\param   iel	   INT  	   (i)	  number of element nodes
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calbkvv(FLUID_DYN_CALC  *dynvar,
		DOUBLE         **estif,   
		DOUBLE          *velint, 
		DOUBLE         **vderxy, 
		DOUBLE          *funct,  
		DOUBLE         **derxy,  
		DOUBLE          *vbubint,  
		DOUBLE         **vbubderxy,  
		DOUBLE           fac,    
		DOUBLE           visc,   
		INT              iel);
/*!---------------------------------------------------------------------                                         
\brief evaluate bubble part of matrix Kvp for fluid3

<pre>                                                       gravem 07/03

In this routine, the bubble part of the matrix Kvp is calculated.

NOTE: there's only one elestif  				   
      --> Kvp is stored in estif[(0..(3*iel-1)][(3*iel)..(4*iel-1)]
      
</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **estif     DOUBLE	   (i/o)  element stiffness matrix
\param  *velint    DOUBLE	   (i)    velocity at int point
\param **vderxy    DOUBLE	   (i)    global velocity derivatives
\param  *funct     DOUBLE	   (i)    natural shape functions
\param **derxy     DOUBLE	   (i)    global deriv. of shape fun.
\param **pbubint   DOUBLE	   (i)    pressure bubble functions
\param***pbubderxy DOUBLE	   (i)    global deriv. of pre. bub. fun.
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   visc      DOUBLE	   (i)    fluid viscosity	     
\param   iel	   INT  	   (i)	  number of element nodes
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calbkvp(FLUID_DYN_CALC  *dynvar,
		DOUBLE         **estif,   
		DOUBLE          *velint, 
		DOUBLE         **vderxy, 
		DOUBLE          *funct,  
		DOUBLE         **derxy,  
		DOUBLE         **pbubint,  
		DOUBLE        ***pbubderxy,  
		DOUBLE           fac,    
		DOUBLE           visc,   
		INT              iel);
/*!---------------------------------------------------------------------                                         
\brief evaluate bubble part of matrix Kpv for fluid3

<pre>                                                       gravem 07/03

In this routine, the bubble part of the matrix Kpv is calculated.

NOTE: there's only one elestif  				   
      --> Kpv is stored in estif[(3*iel)..(4*iel-1)][0..(3*iel-1)]
      
</pre>
\param **estif     DOUBLE	   (i/o)  element stiffness matrix
\param  *funct     DOUBLE	   (i)    natural shape functions
\param **vbubderxy DOUBLE	   (i)    global deriv. of vel. bub. fun.
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   iel	   INT  	   (i)	  number of element nodes
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calbkpv(DOUBLE         **estif,   
		DOUBLE          *funct,  
		DOUBLE         **vbubderxy,  
		DOUBLE           fac,    
		INT              iel);
/*!---------------------------------------------------------------------                                         
\brief evaluate bubble part of matrix Kpp for fluid3

<pre>                                                       gravem 07/03

In this routine, the bubble part of the matrix Kpp is calculated.

NOTE: there's only one elestif  				   
      --> Kpp is stored in estif[(3*iel)..(4*iel-1)][(3*iel)..(4*iel-1)]
      
</pre>
\param **estif     DOUBLE	   (i/o)  element stiffness matrix
\param  *funct     DOUBLE	   (i)    natural shape functions
\param***pbubderxy DOUBLE	   (i)    global deriv. of pre. bub. fun.
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   iel	   INT  	   (i)	  number of element nodes
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calbkpp(DOUBLE         **estif,   
		DOUBLE          *funct,  
		DOUBLE        ***pbubderxy,  
		DOUBLE           fac,    
		INT              iel);
/*!---------------------------------------------------------------------                                         
\brief evaluate bubble part of matrix Mvv for fluid3

<pre>                                                       gravem 07/03

In this routine, the bubble part of the matrix Mvv is calculated.

NOTE: there's only one elemass
      --> Mvv is stored in emass[0..(3*iel-1)][0..(3*iel-1)]  
      
</pre>
\param **emass     DOUBLE	   (i/o)  element mass matrix
\param  *funct     DOUBLE	   (i)    natural shape functions
\param  *vbubint   DOUBLE	   (i)    velocity bubble functions
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   iel	   INT  	   (i)	  number of element nodes
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calbmvv(DOUBLE         **emass,  
	        DOUBLE          *funct, 
	        DOUBLE          *vbubint, 
	        DOUBLE           fac,   
	        INT              iel);
/*!---------------------------------------------------------------------                                         
\brief evaluate bubble part of matrix Mvp for fluid3

<pre>                                                       gravem 07/03

In this routine, the bubble part of the matrix Mvp is calculated.

NOTE: there's only one elemass
      --> Mvp is stored in emass[0..(3*iel-1)][(3*iel)..(4*iel-1)]  
      
</pre>
\param **emass     DOUBLE	   (i/o)  element mass matrix
\param  *funct     DOUBLE	   (i)    natural shape functions
\param **pbubint   DOUBLE	   (i)    pressure bubble functions
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   iel	   INT  	   (i)	  number of element nodes
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calbmvp(DOUBLE         **emass,  
	        DOUBLE          *funct, 
	        DOUBLE         **pbubint, 
	        DOUBLE           fac,   
	        INT              iel);
		
/************************************************************************
 | f3_mlele.c                                                          |
 ************************************************************************/
/*!---------------------------------------------------------------------                                         
\brief control routine for large-scale element integration of fluid3

<pre>                                                       gravem 07/03

This routine controls the element evaluation of the large-scale element:
-actual vel. and pres. variables are set
-for 3-level: control routine for dynamic subgrid viscosity is called
-the control routine for the small-scale solution is called
-the small-scale problem is solved
-the small-scale bubble functions are integrated for this l-s element
-for additional stabilization: stabilization parameters are calculated
-element integration is performed --> element stiffness matrix and 
                                  --> element load vectors
-stiffness matrix and load vectors are permuted for assembling
-element load vector due to dirichlet conditions is calculated				      
			     
</pre>
\param  *data	         FLUID_DATA     (i)
\param  *dynvar	         FLUID_DYN_CALC (i)
\param  *mlvar	         FLUID_DYN_ML   (i)
\param  *submesh	 FLUID_ML_SMESH (i)
\param  *ssmesh	         FLUID_ML_SMESH (i)
\param  *ele	         ELEMENT	(i)   actual element
\param  *estif_global    ARRAY	        (o)   ele stiffnes matrix
\param  *emass_global    ARRAY	        (o)   ele mass matrix
\param  *etforce_global  ARRAY	        (o)   element time force
\param  *eiforce_global  ARRAY	        (o)   ele iteration force
\param  *edforce_global  ARRAY	        (o)   ele dirichlet force
\param  *hasdirich       INT	        (o)   element flag
\param  *hasext          INT	        (o)   element flag
\param   init	         INT	        (i)   init flag
\return void                                               
                                 
------------------------------------------------------------------------*/
void f3_lsele(FLUID_DATA     *data, 
              FLUID_DYN_CALC *dynvar, 
              FLUID_DYN_ML   *mlvar, 
              FLUID_ML_SMESH *submesh, 
              FLUID_ML_SMESH *ssmesh, 
	      ELEMENT	     *ele,	       
              ARRAY	     *estif_global,   
              ARRAY	     *emass_global,   
	      ARRAY	     *etforce_global,	    
	      ARRAY	     *eiforce_global, 
	      ARRAY	     *edforce_global,	      
	      INT	     *hasdirich,      
              INT	     *hasext,
	      INT	      init);
/*!---------------------------------------------------------------------                                         
\brief control routine for submesh element integration of fluid3

<pre>                                                       gravem 07/03

This routine controls the element evaluation of the submesh element:
-element data is set
-stabilization parameter or subgrid viscosity is calculated if necessary
-element integration is performed --> element stiffness matrix and 
                                  --> element load vectors
			     
</pre>
\param  *data	         FLUID_DATA     (i)
\param  *dynvar	         FLUID_DYN_CALC (i)
\param  *mlvar	         FLUID_DYN_ML   (i)
\param  *submesh	 FLUID_ML_SMESH (i)
\param  *ele	         ELEMENT	(i)   actual large-scale element
\param   init	         INT	        (i)   init flag
\return void                                               
                                 
------------------------------------------------------------------------*/
void f3_smele(FLUID_DATA     *data, 
              FLUID_DYN_CALC *dynvar, 
              FLUID_DYN_ML   *mlvar, 
              FLUID_ML_SMESH *submesh, 
	      ELEMENT        *ele,             
              INT             init);
/*!---------------------------------------------------------------------                                         
\brief control routine for submesh bubble function integration of fluid3

<pre>                                                       gravem 07/03

This routine controls the element integration of bubble functions on
the submesh element:
-element data is set
-stabilization parameter or subgrid viscosity is calculated if necessary
-element integration is performed --> element stiffness matrix and 
                                  --> element load vectors
			     
</pre>
\param  *data	         FLUID_DATA     (i)
\param  *dynvar	         FLUID_DYN_CALC (i)
\param  *mlvar	         FLUID_DYN_ML   (i)
\param  *submesh	 FLUID_ML_SMESH (i)
\param  *ele	         ELEMENT	(i)   actual large-scale element
\return void                                               
                                 
------------------------------------------------------------------------*/
void f3_bubele(FLUID_DATA     *data, 
               FLUID_DYN_CALC *dynvar, 
               FLUID_DYN_ML   *mlvar, 
               FLUID_ML_SMESH *submesh, 
	       ELEMENT        *ele);
/*!---------------------------------------------------------------------                                         
\brief control routine for dynamic calc. of subgrid viscosity for fluid3

<pre>                                                       gravem 07/03

This routine controls the dynamic calculation of the subgrid viscosity
on the large-scale element:
-element data is set
-stabilization parameter or subgrid viscosity is calculated if necessary
-element integration is performed --> element stiffness matrix and 
                                  --> element load vectors
			     
</pre>
\param  *data	         FLUID_DATA     (i)
\param  *dynvar	         FLUID_DYN_CALC (i)
\param  *mlvar	         FLUID_DYN_ML   (i)
\param  *submesh	 FLUID_ML_SMESH (i)
\param  *ssmesh	         FLUID_ML_SMESH (i)
\param  *ele	         ELEMENT	(i)   actual large-scale element
\return void                                               
                                 
------------------------------------------------------------------------*/
void f3_dynsgv(FLUID_DATA     *data, 
               FLUID_DYN_CALC *dynvar, 
               FLUID_DYN_ML   *mlvar, 
               FLUID_ML_SMESH *submesh, 
               FLUID_ML_SMESH *ssmesh, 
	       ELEMENT	      *ele);
/*!---------------------------------------------------------------------                                         
\brief submesh global lhs and rhs element integration for fluid4

<pre>                                                       gravem 07/03

This routine controls the elementwise integration of the global lhs 
(only the diffusive part) and rhs on the submesh.
			     
</pre>
\param  *data	         FLUID_DATA     (i)
\param  *submesh	 FLUID_ML_SMESH (i)
\param  *ele	         ELEMENT	(i)  actual large-scale element
\param  *smidiff	 DOUBLE	        (o)  submesh global lhs integral
\param  *smirhs  	 DOUBLE	        (o)  submesh global rhs integral
\return void                                               
                                 
------------------------------------------------------------------------*/
void f3_smintele(FLUID_DATA     *data, 
                 FLUID_ML_SMESH *submesh, 
	         ELEMENT        *ele,
	         DOUBLE         *smidiff,
	         DOUBLE         *smirhs);
/*!---------------------------------------------------------------------                                         
\brief control routine for sub-submesh element integration for fluid3

<pre>                                                       gravem 07/03

This routine controls the element evaluation of the sub-submesh element:
-element data is set
-element integration is performed --> element stiffness matrix and 
                                  --> element load vectors
			     
</pre>
\param  *data	         FLUID_DATA     (i)
\param  *dynvar	         FLUID_DYN_CALC (i)
\param  *mlvar	         FLUID_DYN_ML   (i)
\param  *ssmesh	         FLUID_ML_SMESH (i)
\param  *ele	         ELEMENT	(i)   actual large-scale element
\param   init	         INT	        (i)   init flag
\return void                                               
                                 
------------------------------------------------------------------------*/
void f3_ssele(FLUID_DATA      *data, 
              FLUID_DYN_CALC  *dynvar, 
              FLUID_DYN_ML    *mlvar, 
              FLUID_ML_SMESH  *ssmesh, 
	      ELEMENT         *ele,
	      INT              init);
/*!---------------------------------------------------------------------                                         
\brief sub-submesh element integration of normalized bubble for fluid3

<pre>                                                       gravem 07/03

This routine controls the elementwise integration of the normalized 
bubble function on the sub-submesh.
			     
</pre>
\param  *data	         FLUID_DATA     (i)
\param  *ssmesh  	 FLUID_ML_SMESH (i)
\param  *ele	         ELEMENT	(i)  actual large-scale element
\param  *ssinbu    	 DOUBLE	        (o)  sub-submesh bubble integral
\return void                                               
                                 
------------------------------------------------------------------------*/
void f3_ssintele(FLUID_DATA     *data, 
                 FLUID_ML_SMESH *ssmesh, 
	         ELEMENT        *ele,
		 DOUBLE         *ssinbu);
/*!---------------------------------------------------------------------                                         
\brief control routine for l-s element-based error calculat. for fluid3

<pre>                                                       gravem 07/03

This routine controls the element-based error calculation of the 
large-scale element:
-final vel. and pres. variables are set
-for 3-level: control routine for dynamic subgrid viscosity is called
-the control routine for the small-scale solution is called
-the small-scale problem is solved
-the error calculation incorp. the small scales is done on the submesh
			     
</pre>
\param  *data	         FLUID_DATA     (i)
\param  *fdyn	         FLUID_DYNAMIC  (i)
\param  *dynvar	         FLUID_DYN_CALC (i)
\param  *mlvar	         FLUID_DYN_ML   (i)
\param  *submesh	 FLUID_ML_SMESH (i)
\param  *ssmesh	         FLUID_ML_SMESH (i)
\param  *ele	         ELEMENT	(i) actual large-scale element
\param   actpos          INT            (i) final position in sol. vector
\param  *dulinf  	 DOUBLE         (-) velocity error in Linf-norm 
\param  *dplinf  	 DOUBLE         (-) pressure error in Linf-norm 
\param  *dul2  	         DOUBLE         (-) velocity error in L2-norm 
\param  *dpl2  	         DOUBLE         (-) pressure error in L2-norm 
\param  *duh1  	         DOUBLE         (-) velocity error in H1-norm 
\param  *ul2  	         DOUBLE         (-) analytical velocity in L2-norm 
\param  *pl2  	         DOUBLE         (-) analytical pressure in L2-norm 
\param  *uh1  	         DOUBLE         (-) analytical velocity in H1-norm 
\return void                                               
                                 
------------------------------------------------------------------------*/
void f3_lselerror(FLUID_DATA     *data,
                  FLUID_DYNAMIC  *fdyn, 
                  FLUID_DYN_CALC *dynvar, 
                  FLUID_DYN_ML   *mlvar, 
                  FLUID_ML_SMESH *submesh, 
                  FLUID_ML_SMESH *ssmesh, 
                  ELEMENT        *ele,
                  INT		  actpos, 
                  DOUBLE	 *dulinf, 
                  DOUBLE	 *dplinf, 
                  DOUBLE	 *dul2,   
                  DOUBLE	 *dpl2,   
                  DOUBLE         *duh1,   
                  DOUBLE	 *ul2,	 
                  DOUBLE	 *pl2,	 
                  DOUBLE	 *uh1);	 
/*!---------------------------------------------------------------------                                         
\brief control routine for sm elem. integrat. for error cal. for fluid3

<pre>                                                       gravem 07/03

This routine controls the element evaluation of the submesh element:
-element data is set
-stabilization parameter or subgrid viscosity is calculated if necessary
-element integration is performed --> element stiffness matrix and 
                                  --> element load vectors
			     
</pre>
\param  *data	         FLUID_DATA     (i)
\param  *dynvar	         FLUID_DYN_CALC (i)
\param  *mlvar	         FLUID_DYN_ML   (i)
\param  *submesh	 FLUID_ML_SMESH (i)
\param  *ele	         ELEMENT	(i)   actual large-scale element
\return void                                               
                                 
------------------------------------------------------------------------*/
void f3_smelerror(FLUID_DATA     *data, 
                  FLUID_DYN_CALC *dynvar, 
                  FLUID_DYN_ML   *mlvar, 
                  FLUID_ML_SMESH *submesh, 
	          ELEMENT        *ele);
/*!---------------------------------------------------------------------                                         
\brief control routine for bubble-enhan. error calc. on subm. for fluid3

<pre>                                                       gravem 07/03

This routine controls the error calculation including the small-scale
bubble functions on the submesh element:
-element data is set
-stabilization parameter or subgrid viscosity is calculated if necessary
-element integration is performed --> element stiffness matrix and 
                                  --> element load vectors
			     
</pre>
\param  *data	         FLUID_DATA     (i)
\param  *fdyn	         FLUID_DYNAMIC  (i)
\param  *dynvar	         FLUID_DYN_CALC (i)
\param  *mlvar	         FLUID_DYN_ML   (i)
\param  *submesh	 FLUID_ML_SMESH (i)
\param  *ele	         ELEMENT	(i)   actual large-scale element
\param  *dulinf  	 DOUBLE         (-) velocity error in Linf-norm 
\param  *dplinf  	 DOUBLE         (-) pressure error in Linf-norm 
\param  *dul2  	         DOUBLE         (-) velocity error in L2-norm 
\param  *dpl2  	         DOUBLE         (-) pressure error in L2-norm 
\param  *duh1  	         DOUBLE         (-) velocity error in H1-norm 
\param  *ul2  	         DOUBLE         (-) analytical velocity in L2-norm 
\param  *pl2  	         DOUBLE         (-) analytical pressure in L2-norm 
\param  *uh1  	         DOUBLE         (-) analytical velocity in H1-norm 
\return void                                               
                                 
------------------------------------------------------------------------*/
void f3_bubelerror(FLUID_DATA     *data, 
                   FLUID_DYNAMIC  *fdyn, 
                   FLUID_DYN_CALC *dynvar, 
                   FLUID_DYN_ML   *mlvar, 
                   FLUID_ML_SMESH *submesh, 
	           ELEMENT        *ele,
                   DOUBLE         *dulinf,
                   DOUBLE         *dplinf,
                   DOUBLE         *dul2,  
                   DOUBLE         *dpl2,  
                   DOUBLE         *duh1,  
                   DOUBLE         *ul2,   
                   DOUBLE         *pl2,   
                   DOUBLE         *uh1);
		      
/************************************************************************
 | f3_mlelesize.c                                                          |
 ************************************************************************/
/*!---------------------------------------------------------------------
\brief routine to calculate submesh element size for fluid3

<pre>                                                       gravem 07/03

In this routine, the characteristic submesh element length is calculated
and the routine for the calculation of the stabilization parameter or 
the subgrid viscosity (depending on the respective flag) is called. 

</pre>
\param  *ele       ELEMENT	     (i)   actual element
\param  *data      FLUID_DATA	     (i)
\param  *dynvar    FLUID_DYN_CALC    (i/o)
\param  *mlvar     FLUID_DYN_ML      (i)
\param  *funct     DOUBLE 	     (-)   l-s shape functions
\param **deriv     DOUBLE 	     (-)   deriv. of l-s shape funcs
\param **deriv2    DOUBLE 	     (-)   2nd deriv. of l-s sh. funcs
\param  *smfunct   DOUBLE 	     (-)   submesh shape functions
\param **smderiv   DOUBLE 	     (-)   deriv. of submesh shape fun.
\param **smderiv2  DOUBLE 	     (-)   2nd deriv. of sm sh. fun.
\param **derxy     DOUBLE 	     (-)   global l-s sh. fun. deriv.
\param **xjm       DOUBLE 	     (-)   jacobian matrix
\param **evel      DOUBLE 	     (i)   l-s element velocities
\param  *velint    DOUBLE 	     (-)   l-s vel. at integr. point
\param **vderxy    DOUBLE 	     (-)   l-s vel. deriv. at int. point
\param **smxyze    DOUBLE 	     (i)   submesh element coordinates
\param **smxyzep   DOUBLE 	     (i)   sm ele. coord. on parent dom.
\param **wa1       DOUBLE 	     (-)   working array
\return void             

------------------------------------------------------------------------*/
void f3_smelesize(ELEMENT         *ele,    
		  FLUID_DATA	  *data, 
		  FLUID_DYN_CALC  *dynvar,
		  FLUID_DYN_ML    *mlvar,
	          DOUBLE	  *funct,  
	          DOUBLE	 **deriv,  
	          DOUBLE	 **deriv2,		
	          DOUBLE	  *smfunct,  
	          DOUBLE	 **smderiv,  
	          DOUBLE	 **smderiv2,		
	          DOUBLE	 **derxy,  
		  DOUBLE	 **xjm,    
		  DOUBLE	 **evel,		 
		  DOUBLE	  *velint, 
	          DOUBLE	 **vderxy,  
	          DOUBLE	 **smxyze,  
	          DOUBLE	 **smxyzep,  
		  DOUBLE	 **wa1);
		  
/************************************************************************
 | f3_mlgalmat.c                                                          |
 ************************************************************************/
/*!---------------------------------------------------------------------                                         
\brief evaluate Galerkin part of submesh stiffness matrix SMK for fluid3

<pre>                                                       gravem 07/03

In this routine, the Galerkin part of the submesh stiffness matrix SMK 
is calculated.

</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *mlvar     FLUID_DYN_ML    (i)
\param **smestif   DOUBLE	   (i/o)  submesh ele stiffness matrix
\param  *velint    DOUBLE	   (i)    velocity at int. point
\param **vderxy    DOUBLE	   (i)    global vel. deriv. at int. p.
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun. 
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   visc      DOUBLE	   (i)    fluid viscosity	     
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calsmk(FLUID_DYN_CALC  *dynvar,
	       FLUID_DYN_ML    *mlvar, 
	       DOUBLE         **smestif,   
	       DOUBLE          *velint, 
	       DOUBLE         **vderxy, 
	       DOUBLE          *smfunct,  
	       DOUBLE         **smderxy,  
	       DOUBLE           fac,    
	       DOUBLE           visc,   
	       INT              smiel);
/*!---------------------------------------------------------------------                                         
\brief evaluate Galerkin part of submesh mass matrix SMM for fluid3

<pre>                                                       gravem 07/03

In this routine, the Galerkin part of the submesh mass matrix SMM is 
calculated.

</pre>
\param **smemass   DOUBLE	   (i/o)  submesh element mass matrix
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calsmm(DOUBLE         **smemass,   
	       DOUBLE          *smfunct,  
	       DOUBLE           fac,    
	       INT              smiel);
/*!---------------------------------------------------------------------                                         
\brief evaluate diffusive part of submesh stiffn. matrix SMK for fluid3

<pre>                                                       gravem 07/03

In this routine, the diffusive Galerkin part of the submesh stiffness 
matrix SMK is calculated.

</pre>
\param **smiediff  DOUBLE	   (i/o)  sm ele stiffn. matrix (diff.)
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun. 
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calsmkd(DOUBLE         **smiediff,   
        	DOUBLE         **smderxy,  
	        DOUBLE           fac,    
	        INT              smiel);
		
/************************************************************************
 | f3_mlint.c                                                          |
 ************************************************************************/
/*!---------------------------------------------------------------------
\brief integration loop for submesh element for fluid3

<pre>                                                       gravem 07/03

In this routine, the element stiffness matrix, mass matrix, VMM-RHS and
Time-RHS for one submesh element is calculated.
      
</pre>
\param  *data        FLUID_DATA	   (i)    integration data
\param  *ele	     ELEMENT	   (i)    actual element
\param  *dynvar      FLUID_DYN_CALC(i)
\param  *mlvar       FLUID_DYN_ML  (i)
\param  *submesh     FLUID_ML_SMESH(i)   
\param **smestif      DOUBLE       (o)	sm element stiffness matrix
\param **smemass      DOUBLE       (o)	sm element mass matrix
\param  *smevfor      DOUBLE       (o)	sm element VMM force vector
\param  *smetfor      DOUBLE       (o)	sm element time force vector
\param  *smxyze       DOUBLE       (i)	submesh element coordinates
\param  *smxyzep      DOUBLE       (i)	sm ele. coord. on parent dom.
\param  *funct        DOUBLE       (-)	natural shape functions
\param **deriv        DOUBLE       (-)	deriv. of nat. shape funcs
\param **deriv2       DOUBLE       (-)	2nd deriv. of nat. shape f.
\param **xjm	      DOUBLE       (-)	jacobian matrix
\param **derxy        DOUBLE       (-)	global derivatives
\param **derxy2       DOUBLE       (-)	2nd global derivatives
\param  *smfunct      DOUBLE       (-)	sm natural shape functions
\param **smderiv      DOUBLE       (-)  sm deriv. of nat. shape funcs
\param **smderiv2     DOUBLE       (-)  sm 2nd deriv. of nat. shape f.
\param **smxjm	      DOUBLE       (-)  sm jacobian matrix
\param **smderxy      DOUBLE       (-)  sm global derivatives
\param **smderxy2     DOUBLE       (-)  sm 2nd global derivatives
\param **eveln        DOUBLE       (i)  ele vel. at time step n
\param **evel         DOUBLE       (i)  ele vel. at time step n+1
\param  *epren        DOUBLE       (i)  ele pres. at time step n
\param  *epre         DOUBLE       (i)  ele pres. at time step n+1
\param **evbub        DOUBLE       (i)  sm ele vel. bubble functions
\param **epbub        DOUBLE       (i)  sm ele pre. bubble functions
\param **efbub        DOUBLE       (i)  sm ele rhs bubble functions
\param **evbubn       DOUBLE       (i)  sm ele vel. bubble fun. at n
\param **epbubn       DOUBLE       (i)  sm ele pre. bubble fun. at n
\param **efbubn       DOUBLE       (i)  sm ele rhs bubble fun. at n
\param  *vbubint      DOUBLE       (-)  vel. bubble fun. at int. p.
\param  *vbubderxy    DOUBLE       (-)  vel. bub. fun. der. at int. p.
\param  *vbubderxy2   DOUBLE       (-)  2nd vel. bub. fun. der. at i.p.
\param  *pbubint      DOUBLE       (-)  pre. bubble fun. at int. p.
\param  *pbubderxy    DOUBLE       (-)  pre. bub. fun. der. at int. p.
\param  *pbubderxy2   DOUBLE       (-)  2nd pre. bub. fun. der. at i.p.
\param  *vbubintn     DOUBLE       (-)  vel. bubble fun. at int. p. at n
\param  *vbubderxyn   DOUBLE       (-)  vel. bub. fun. der. at int. p. at n
\param  *vbubderxy2n  DOUBLE       (-)  2nd vel. bub. fun. der. at i.p. at n
\param  *pbubintn     DOUBLE       (-)  pre. bubble fun. at int. p. at n
\param  *pbubderxyn   DOUBLE       (-)  pre. bub. fun. der. at int. p. at n
\param  *pbubderxy2n  DOUBLE       (-)  2nd pre. bub. fun. der. at i.p. at n
\param  *velint       DOUBLE       (-)  vel at integration point
\param  *velintn      DOUBLE       (-)  vel at integration point at n
\param  *velintnt     DOUBLE       (-)  'temporal' vel at int. p. at n
\param  *velintnc     DOUBLE       (-)  'convective' vel at int. p. at n
\param **vderxy       DOUBLE       (-)  global vel. derivatives
\param **vderxyn      DOUBLE       (-)  global vel. derivatives at n
\param **vderxync     DOUBLE       (-)  global 'convective' vel. der. at n
\param **vderxynv     DOUBLE       (-)  global 'viscous' vel. der. at n
\param **vderxy2      DOUBLE       (-)  2nd global vel. deriv.
\param **vderxy2n     DOUBLE       (-)  2nd global vel. derivatives at n
\param **vderxy2nv    DOUBLE       (-)  2nd global 'viscous' vel. der. at n
\param  *pderxyn      DOUBLE       (-)  global pres. derivatives at n
\param  *smvelint     DOUBLE       (-)  sm vel at integration point
\param **smvderxy     DOUBLE       (-)  sm global vel. derivatives
\param  *smpreint     DOUBLE       (-)  sm pre at integration point
\param **smpderxy     DOUBLE       (-)  sm global pre. derivatives
\param  *smvelintn    DOUBLE       (-)  sm vel at integration point at n
\param **smvderxyn    DOUBLE       (-)  sm global vel. derivatives at n
\param **smvderxy2n   DOUBLE       (-)  2nd sm global vel. derivatives at n
\param  *smpreintn    DOUBLE       (-)  sm pre at integration point at n
\param **smpderxyn    DOUBLE       (-)  sm global pre. derivatives at n
\param **smpderxy2n   DOUBLE       (-)  2nd sm global pre. derivatives at n
\param  *smfint       DOUBLE       (-)  sm rhs at integration point
\param **smfderxy     DOUBLE       (-)  sm global rhs. derivatives
\param  *smfintn      DOUBLE       (-)  sm rhs at integration point at n
\param **smfderxyn    DOUBLE       (-)  sm global rhs. derivatives at n
\param **smfderxy2n   DOUBLE       (-)  2nd sm global rhs. derivatives at n
\param **wa1	      DOUBLE       (-)  working array
\param **wa2	      DOUBLE       (-)  working array
\return void                                                   

------------------------------------------------------------------------*/
void f3_smint(FLUID_DATA      *data,     
	      ELEMENT	      *ele,	
	      FLUID_DYN_CALC  *dynvar, 
	      FLUID_DYN_ML    *mlvar, 
	      FLUID_ML_SMESH  *submesh, 
              DOUBLE	     **smestif,   
	      DOUBLE	     **smemass,   
	      DOUBLE	     **smevfor, 
	      DOUBLE	     **smetfor, 
	      DOUBLE	     **smxyze, 
	      DOUBLE	     **smxyzep, 
	      DOUBLE	      *funct,	
	      DOUBLE	     **deriv,	
	      DOUBLE	     **deriv2,  
	      DOUBLE	     **xjm,	
	      DOUBLE	     **derxy,	
	      DOUBLE	     **derxy2,  
	      DOUBLE	      *smfunct,   
	      DOUBLE	     **smderiv,   
	      DOUBLE	     **smderiv2,  
	      DOUBLE	     **smxjm,	  
	      DOUBLE	     **smderxy,   
	      DOUBLE	     **smderxy2,  
	      DOUBLE	     **eveln,	
	      DOUBLE	     **evel,  
	      DOUBLE	      *epren,	
	      DOUBLE	      *epre,
	      DOUBLE	     **evbub,	
	      DOUBLE	     **epbub,	
	      DOUBLE	     **efbub,	
	      DOUBLE	     **evbubn,   
	      DOUBLE	     **epbubn,   
	      DOUBLE	     **efbubn,   
              DOUBLE	      *vbubint,    
              DOUBLE	     **vbubderxy,  
              DOUBLE	     **vbubderxy2, 
              DOUBLE	     **pbubint,    
              DOUBLE	    ***pbubderxy,  
              DOUBLE	    ***pbubderxy2, 
              DOUBLE	      *vbubintn,   
              DOUBLE	     **vbubderxyn, 
              DOUBLE	     **vbubderxy2n,
              DOUBLE	     **pbubintn,   
              DOUBLE	    ***pbubderxyn, 
              DOUBLE	    ***pbubderxy2n,
	      DOUBLE	      *velint,  
              DOUBLE	      *velintn,   
              DOUBLE	      *velintnt,  
              DOUBLE	      *velintnc,  
	      DOUBLE	     **vderxy,  
              DOUBLE	     **vderxyn,   
              DOUBLE	     **vderxync,  
              DOUBLE	     **vderxynv,  
	      DOUBLE	     **vderxy2, 
              DOUBLE	     **vderxy2n,  
              DOUBLE	     **vderxy2nv,  
              DOUBLE	      *pderxyn,   
              DOUBLE	      *smvelint,  
              DOUBLE	     **smvderxy,  
              DOUBLE	      *smpreint,  
              DOUBLE	     **smpderxy,  
              DOUBLE	      *smvelintn, 
              DOUBLE	     **smvderxyn, 
              DOUBLE	     **smvderxy2n, 
              DOUBLE	      *smpreintn,  
              DOUBLE	     **smpderxyn,  
              DOUBLE	     **smpderxy2n, 
              DOUBLE	      *smfint,	 
              DOUBLE	     **smfderxy,  
              DOUBLE	      *smfintn,    
              DOUBLE	     **smfderxyn,  
              DOUBLE	     **smfderxy2n,
	      DOUBLE	     **wa1,	
	      DOUBLE	     **wa2);
/*!---------------------------------------------------------------------
\brief integration loop for bubble funct. on submesh element for fluid3

<pre>                                                       gravem 07/03

In this routine, the bubble function part of the large-scale element 
stiffness matrix, mass matrix, VMM-RHS and Time-RHS for one submesh 
element is calculated.
      
</pre>
\param  *data        FLUID_DATA	   (i)    integration data
\param  *ele	     ELEMENT	   (i)    actual element
\param  *dynvar      FLUID_DYN_CALC(i)
\param  *mlvar       FLUID_DYN_ML  (i)
\param  *submesh     FLUID_ML_SMESH(i)   
\param **estif        DOUBLE	   (o)  element stiffness matrix
\param **emass        DOUBLE	   (o)  element mass matrix
\param  *eiforce      DOUBLE	   (o)  element iteration force vector
\param  *smxyze	      DOUBLE	   (i)  submesh element coordinates
\param  *smxyzep      DOUBLE	   (i)  sm ele. coord. on parenyt dom.
\param  *funct        DOUBLE       (-)	natural shape functions
\param **deriv        DOUBLE       (-)	deriv. of nat. shape funcs
\param **deriv2       DOUBLE       (-)	2nd deriv. of nat. shape f.
\param **xjm	      DOUBLE       (-)	jacobian matrix
\param **derxy        DOUBLE       (-)	global derivatives
\param  *smfunct      DOUBLE       (-)	sm natural shape functions
\param **smderiv      DOUBLE       (-)  sm deriv. of nat. shape funcs
\param **smderiv2     DOUBLE       (-)  sm 2nd deriv. of nat. shape f.
\param **smxjm	      DOUBLE       (-)  sm jacobian matrix
\param **smderxy      DOUBLE       (-)  sm global derivatives
\param **evel         DOUBLE       (i)  ele vel. at time step n+1
\param  *epre         DOUBLE       (i)  ele pres. at time step n+1
\param **evbub        DOUBLE       (i)  sm ele vel. bubble functions
\param **epbub        DOUBLE       (i)  sm ele pre. bubble functions
\param **efbub        DOUBLE       (i)  sm ele rhs bubble functions
\param  *vbubint      DOUBLE       (-)  vel. bubble fun. at int. p.
\param  *vbubderxy    DOUBLE       (-)  vel. bub. fun. der. at int. p.
\param  *pbubint      DOUBLE       (-)  pre. bubble fun. at int. p.
\param  *pbubderxy    DOUBLE       (-)  pre. bub. fun. der. at int. p.
\param  *covint       DOUBLE       (-)  convective vel at int. point
\param  *velint       DOUBLE       (-)  vel at integration point
\param **vderxy       DOUBLE       (-)  global vel. derivatives
\param  *smvelint     DOUBLE       (-)  sm vel at integration point
\param **smvderxy     DOUBLE       (-)  sm global vel. derivatives
\param  *smpreint     DOUBLE       (-)  sm pre at integration point
\param **smpderxy     DOUBLE       (-)  sm global pre. derivatives
\param  *smfint       DOUBLE       (-)  sm rhs at integration point
\param **smfderxy     DOUBLE       (-)  sm global rhs. derivatives
\param **wa1	      DOUBLE       (-)  working array
\param **wa2	      DOUBLE       (-)  working array
\return void                                                   

------------------------------------------------------------------------*/
void f3_bubint(FLUID_DATA      *data,     
	       ELEMENT         *ele,     
	       FLUID_DYN_CALC  *dynvar, 
	       FLUID_DYN_ML    *mlvar, 
	       FLUID_ML_SMESH  *submesh, 
               DOUBLE	      **estif,	
	       DOUBLE	      **emass,	
	       DOUBLE	       *eiforce, 
	       DOUBLE         **smxyze, 
	       DOUBLE         **smxyzep, 
	       DOUBLE          *funct,   
	       DOUBLE         **deriv,   
	       DOUBLE         **deriv2,   
	       DOUBLE         **xjm,     
	       DOUBLE         **derxy,   
	       DOUBLE          *smfunct,   
	       DOUBLE         **smderiv,   
	       DOUBLE         **smderiv2,   
	       DOUBLE         **smxjm,     
	       DOUBLE         **smderxy,   
	       DOUBLE         **evel,  
	       DOUBLE          *epre,
	       DOUBLE         **evbub,   
	       DOUBLE         **epbub,   
	       DOUBLE         **efbub,   
               DOUBLE          *vbubint,    
               DOUBLE         **vbubderxy,  
               DOUBLE         **pbubint,    
               DOUBLE        ***pbubderxy,  
	       DOUBLE	       *covint,  
	       DOUBLE          *velint,  
	       DOUBLE         **vderxy,  
               DOUBLE          *smvelint,  
               DOUBLE         **smvderxy,  
               DOUBLE          *smpreint,  
               DOUBLE         **smpderxy,  
               DOUBLE          *smfint,    
               DOUBLE         **smfderxy,  
	       DOUBLE         **wa1,     
	       DOUBLE         **wa2);
/*!---------------------------------------------------------------------
\brief integration loop for one large-scale element for fluid3

<pre>                                                       gravem 07/03

In this routine the large-scale part of the large-scale element 
stiffness matrix, iteration-RHS, time-RHS and external-RHS for one 
large-scale element is calculated
      
</pre>
\param  *data      FLUID_DATA	   (i)	  integration data
\param  *ele	   ELEMENT	   (i)    actual element
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *mlvar     FLUID_DYN_ML    (i)
\param  *hasext    INT             (i)    element flag
\param **estif     DOUBLE	   (o)    element stiffness matrix
\param **emass     DOUBLE	   (o)    element mass matrix
\param  *eiforce   DOUBLE	   (o)    element iter force vector
\param  *etforce   DOUBLE	   (o)    element time force vector
\param  *funct     DOUBLE	   (-)    natural shape functions
\param **deriv     DOUBLE	   (-)	  deriv. of nat. shape funcs
\param **deriv2    DOUBLE	   (-)    2nd deriv. of nat. shape f.
\param **xjm	   DOUBLE	   (-)    jacobian matrix
\param **derxy     DOUBLE	   (-)	  global derivatives
\param **derxy2    DOUBLE	   (-)    2nd global derivatives
\param **evel      DOUBLE	   (i)    ele vel. at time n+1
\param **eveln     DOUBLE	   (i)    ele vel. at time n
\param  *epren     DOUBLE	   (-)    ele pres. at time n
\param  *edeadn    DOUBLE	   (-)    ele dead load at n 
\param  *edead     DOUBLE	   (-)    ele dead load at n+1
\param  *velint    DOUBLE	   (-)    vel at integration point
\param  *velintn   DOUBLE	   (-)    vel at integration point at n
\param  *covint    DOUBLE	   (-)    conv. vel. at integr. point
\param  *covintn   DOUBLE	   (-)    conv. vel. at integr. point at n
\param **vderxy    DOUBLE	   (-)    global vel. derivatives
\param **vderxyn   DOUBLE	   (-)    global vel. derivatives at n
\param  *pderxy    DOUBLE	   (-)    global pres. derivatives
\param **wa1	   DOUBLE	   (-)    working array
\param **wa2	   DOUBLE	   (-)    working array
\return void                                                   

------------------------------------------------------------------------*/
void f3_lsint(FLUID_DATA      *data,     
	      ELEMENT	      *ele,	
	      FLUID_DYN_CALC  *dynvar,
	      FLUID_DYN_ML    *mlvar, 
              INT             *hasext,
              DOUBLE	     **estif,	
	      DOUBLE	     **emass,	
	      DOUBLE	      *eiforce, 
	      DOUBLE	      *etforce, 
	      DOUBLE	      *funct,	
	      DOUBLE	     **deriv,	
	      DOUBLE	     **deriv2,  
	      DOUBLE	     **xjm,	
	      DOUBLE	     **derxy,	
	      DOUBLE	     **derxy2,  
	      DOUBLE	     **evel,  
	      DOUBLE	     **eveln,  
	      DOUBLE          *epren,   
	      DOUBLE          *edeadn,   
	      DOUBLE          *edead,   
	      DOUBLE	      *velint,  
              DOUBLE          *velintn,   
	      DOUBLE	      *covint,  
	      DOUBLE	      *covintn,  
	      DOUBLE	     **vderxy,  
              DOUBLE         **vderxyn,  
	      DOUBLE	     **wa1,	
	      DOUBLE	     **wa2);
/*!---------------------------------------------------------------------
\brief integration loop for dynamic procedure for sm ele. for fluid3

<pre>                                                       gravem 07/03

In this routine, the integral of the diffusive part of the lhs as well
as the normalized rhs are determined.
      
</pre>
\param  *data        FLUID_DATA	   (i)    integration data
\param  *ele	     ELEMENT	   (i)    actual element
\param  *submesh     FLUID_ML_SMESH(i)   
\param **smiediff     DOUBLE       (o)	sm ele. int. of diffusive lhs
\param  *smierhs      DOUBLE       (o)	sm ele. int. of rhs
\param  *smxyze       DOUBLE       (i)	submesh element coordinates
\param  *smfunct      DOUBLE       (-)	sm natural shape functions
\param **smderiv      DOUBLE       (-)  sm deriv. of nat. shape funcs
\param **smderiv2     DOUBLE       (-)  sm 2nd deriv. of nat. shape f.
\param **smxjm	      DOUBLE       (-)  sm jacobian matrix
\param **smderxy      DOUBLE       (-)  sm global derivatives
\return void                                                   

------------------------------------------------------------------------*/
void f3_smint2(FLUID_DATA      *data,     
	       ELEMENT	       *ele,	
	       FLUID_ML_SMESH  *submesh, 
               DOUBLE	      **smediff,   
	       DOUBLE	       *smerhs,   
	       DOUBLE	      **smxyze, 
	       DOUBLE	       *smfunct,   
	       DOUBLE	      **smderiv,   
	       DOUBLE	      **smderiv2,  
	       DOUBLE	      **smxjm,	  
	       DOUBLE	      **smderxy,
	       DOUBLE	      **wa1);
/*!---------------------------------------------------------------------
\brief integration loop for sub-submesh element for fluid2

<pre>                                                       gravem 07/03

In this routine, the element stiffness matrix and normalized RHS for 
one sub-submesh element is calculated.
      
</pre>
\param  *data        FLUID_DATA	   (i)    integration data
\param  *ele	     ELEMENT	   (i)    actual element
\param  *dynvar      FLUID_DYN_CALC(i)
\param  *mlvar       FLUID_DYN_ML  (i)
\param  *ssmesh      FLUID_ML_SMESH(i)   
\param **ssestif      DOUBLE       (o)	ssm element stiffness matrix
\param  *ssenfor      DOUBLE       (o)	ssm element norm. force vector
\param  *ssxyze       DOUBLE       (i)	sub-submesh element coordinates
\param  *ssxyzep      DOUBLE       (i)	ssm ele. coord. on parent dom.
\param  *funct        DOUBLE       (-)	natural shape functions
\param **deriv        DOUBLE       (-)	deriv. of nat. shape funcs
\param **deriv2       DOUBLE       (-)	2nd deriv. of nat. shape f.
\param **xjm	      DOUBLE       (-)	jacobian matrix
\param **derxy        DOUBLE       (-)	global derivatives
\param  *ssfunct      DOUBLE       (-)	sm natural shape functions
\param **ssderiv      DOUBLE       (-)  sm deriv. of nat. shape funcs
\param **ssderiv2     DOUBLE       (-)  sm 2nd deriv. of nat. shape f.
\param **ssxjm	      DOUBLE       (-)  sm jacobian matrix
\param **ssderxy      DOUBLE       (-)  sm global derivatives
\param **evel         DOUBLE       (i)  ele vel. at time step n+1
\param  *velint       DOUBLE       (-)  vel at integration point
\param **vderxy       DOUBLE       (-)  global vel. derivatives
\param **wa1          DOUBLE       (-)  working array
\return void                                                   

------------------------------------------------------------------------*/
void f3_ssint(FLUID_DATA      *data,     
	      ELEMENT	      *ele,	
	      FLUID_DYN_CALC  *dynvar, 
	      FLUID_DYN_ML    *mlvar, 
	      FLUID_ML_SMESH  *ssmesh, 
              DOUBLE	     **ssestif,   
	      DOUBLE	      *ssenfor,   
	      DOUBLE	     **ssxyze, 
	      DOUBLE	     **ssxyzep, 
	      DOUBLE	      *funct,	
	      DOUBLE	     **deriv,	
	      DOUBLE	     **deriv2,  
	      DOUBLE	     **xjm,	
	      DOUBLE	     **derxy,	
	      DOUBLE	      *ssfunct,   
	      DOUBLE	     **ssderiv,   
	      DOUBLE	     **ssderiv2,  
	      DOUBLE	     **ssxjm,	  
	      DOUBLE	     **ssderxy,   
	      DOUBLE	     **evel,  
	      DOUBLE	      *velint,  
	      DOUBLE	     **vderxy,
	      DOUBLE	     **wa1);
/*!---------------------------------------------------------------------
\brief integration loop for norm. bubble on sub-submesh ele. for fluid3

<pre>                                                       gravem 07/03

In this routine, the integral of the normalized bubble function for 
one sub-submesh element is calculated.
      
</pre>
\param  *data        FLUID_DATA	   (i)    integration data
\param  *ele	     ELEMENT	   (i)    actual element
\param  *ssmesh      FLUID_ML_SMESH(i)   
\param  *ssinbu       DOUBLE	   (o)  sub-submesh bubble integral
\param  *ebub         DOUBLE       (o)	ssm element bubble function
\param **ssxyze       DOUBLE       (i)	sub-submesh element coordinates
\param  *ssfunct      DOUBLE       (-)	sm natural shape functions
\param **ssderiv      DOUBLE       (-)  sm deriv. of nat. shape funcs
\param **ssderiv2     DOUBLE       (-)  sm 2nd deriv. of nat. shape f.
\param **ssxjm	      DOUBLE       (-)  sm jacobian matrix
\return void                                                   

------------------------------------------------------------------------*/
void f3_inbu(FLUID_DATA      *data,     
	     ELEMENT	     *ele,   
	     FLUID_ML_SMESH  *ssmesh, 
             DOUBLE	     *ssinbu,   
             DOUBLE	     *ebub,   
	     DOUBLE	    **ssxyze, 
	     DOUBLE	     *ssfunct,   
	     DOUBLE	    **ssderiv,   
	     DOUBLE	    **ssderiv2,  
	     DOUBLE	    **ssxjm);
/*!---------------------------------------------------------------------
\brief error calculation integration loop for submesh element for fluid3

<pre>                                                       gravem 07/03

In this routine, the element stiffness matrix and VMM-RHS for one 
submesh element is calculated for error calculation.
      
</pre>
\param  *data        FLUID_DATA	   (i)    integration data
\param  *ele	     ELEMENT	   (i)    actual element
\param  *dynvar      FLUID_DYN_CALC(i)
\param  *mlvar       FLUID_DYN_ML  (i)
\param  *submesh     FLUID_ML_SMESH(i)   
\param **smestif      DOUBLE       (o)	sm element stiffness matrix
\param  *smevfor      DOUBLE       (o)	sm element VMM force vector
\param  *smxyze       DOUBLE       (i)	submesh element coordinates
\param  *smxyzep      DOUBLE       (i)	sm ele. coord. on parent dom.
\param  *funct        DOUBLE       (-)	natural shape functions
\param **deriv        DOUBLE       (-)	deriv. of nat. shape funcs
\param **deriv2       DOUBLE       (-)	2nd deriv. of nat. shape f.
\param **xjm	      DOUBLE       (-)	jacobian matrix
\param **derxy        DOUBLE       (-)	global derivatives
\param **derxy2       DOUBLE       (-)	2nd global derivatives
\param  *smfunct      DOUBLE       (-)	sm natural shape functions
\param **smderiv      DOUBLE       (-)  sm deriv. of nat. shape funcs
\param **smderiv2     DOUBLE       (-)  sm 2nd deriv. of nat. shape f.
\param **smxjm	      DOUBLE       (-)  sm jacobian matrix
\param **smderxy      DOUBLE       (-)  sm global derivatives
\param **smderxy2     DOUBLE       (-)  sm 2nd global derivatives
\param **evel         DOUBLE       (i)  ele vel. at time step n+1
\param  *epre         DOUBLE       (i)  ele pres. at time step n+1
\param **evbub        DOUBLE       (i)  sm ele vel. bubble functions
\param **epbub        DOUBLE       (i)  sm ele pre. bubble functions
\param **efbub        DOUBLE       (i)  sm ele rhs bubble functions
\param  *vbubint      DOUBLE       (-)  vel. bubble fun. at int. p.
\param  *vbubderxy    DOUBLE       (-)  vel. bub. fun. der. at int. p.
\param  *vbubderxy2   DOUBLE       (-)  2nd vel. bub. fun. der. at i.p.
\param  *pbubint      DOUBLE       (-)  pre. bubble fun. at int. p.
\param  *pbubderxy    DOUBLE       (-)  pre. bub. fun. der. at int. p.
\param  *pbubderxy2   DOUBLE       (-)  2nd pre. bub. fun. der. at i.p.
\param  *velint       DOUBLE       (-)  vel at integration point
\param **vderxy       DOUBLE       (-)  global vel. derivatives
\param  *smvelint     DOUBLE       (-)  sm vel at integration point
\param **smvderxy     DOUBLE       (-)  sm global vel. derivatives
\param  *smpreint     DOUBLE       (-)  sm pre at integration point
\param **smpderxy     DOUBLE       (-)  sm global pre. derivatives
\param  *smfint       DOUBLE       (-)  sm rhs at integration point
\param **smfderxy     DOUBLE       (-)  sm global rhs. derivatives
\param **wa1	      DOUBLE       (-)  working array
\param **wa2	      DOUBLE       (-)  working array
\return void                                                   

------------------------------------------------------------------------*/
void f3_sminterr(FLUID_DATA      *data,     
	         ELEMENT	 *ele,	
	         FLUID_DYN_CALC  *dynvar, 
	         FLUID_DYN_ML    *mlvar, 
	         FLUID_ML_SMESH  *submesh, 
                 DOUBLE	        **smestif,   
	         DOUBLE	        **smevfor, 
	         DOUBLE	        **smxyze, 
	         DOUBLE	        **smxyzep, 
	         DOUBLE	         *funct,	
	         DOUBLE	        **deriv,	
	         DOUBLE	        **deriv2,  
	         DOUBLE	        **xjm,	
	         DOUBLE	        **derxy,	
	         DOUBLE	        **derxy2,  
	         DOUBLE	         *smfunct,   
	         DOUBLE	        **smderiv,   
	         DOUBLE	        **smderiv2,  
	         DOUBLE	        **smxjm,	  
	         DOUBLE	        **smderxy,   
	         DOUBLE	        **smderxy2,  
	         DOUBLE	        **evel,  
	         DOUBLE	         *epre,
	         DOUBLE	        **evbub,	
	         DOUBLE	        **epbub,	
	         DOUBLE	        **efbub,	
                 DOUBLE	         *vbubint,    
                 DOUBLE	        **vbubderxy,  
                 DOUBLE	        **vbubderxy2, 
                 DOUBLE	        **pbubint,    
                 DOUBLE	       ***pbubderxy,  
                 DOUBLE	       ***pbubderxy2, 
	         DOUBLE	         *velint,  
	         DOUBLE	        **vderxy,  
                 DOUBLE	         *smvelint,  
                 DOUBLE	        **smvderxy,  
                 DOUBLE	         *smpreint,  
                 DOUBLE	        **smpderxy,  
                 DOUBLE	         *smfint,	 
                 DOUBLE	        **smfderxy,  
	         DOUBLE	        **wa1,	
	         DOUBLE	        **wa2);
/*!---------------------------------------------------------------------
\brief error calc. integr. loop for bubble funct. on sm ele. for fluid2

<pre>                                                       gravem 07/03

In this routine, the bubble function part of the large-scale element 
stiffness matrix and VMM-RHS for one submesh element is calculated for
error calculation.
      
</pre>
\param  *data        FLUID_DATA	   (i)    integration data
\param  *fdyn        FLUID_DYNAMIC (i)
\param  *ele	     ELEMENT	   (i)    actual element
\param  *mlvar       FLUID_DYN_ML  (i)
\param  *submesh     FLUID_ML_SMESH(i)   
\param  *smxyze	      DOUBLE	   (i)  submesh element coordinates
\param  *smxyzep      DOUBLE	   (i)  sm ele. coord. on parent dom.
\param  *funct        DOUBLE       (-)	natural shape functions
\param **deriv        DOUBLE       (-)	deriv. of nat. shape funcs
\param **deriv2       DOUBLE       (-)	2nd deriv. of nat. shape f.
\param **xjm	      DOUBLE       (-)	jacobian matrix
\param **derxy        DOUBLE       (-)	global derivatives
\param  *smfunct      DOUBLE       (-)	sm natural shape functions
\param **smderiv      DOUBLE       (-)  sm deriv. of nat. shape funcs
\param **smderiv2     DOUBLE       (-)  sm 2nd deriv. of nat. shape f.
\param **smxjm	      DOUBLE       (-)  sm jacobian matrix
\param **smderxy      DOUBLE       (-)  sm global derivatives
\param **evel         DOUBLE       (i)  ele vel. at time step n+1
\param  *epre         DOUBLE       (i)  ele pres. at time step n+1
\param **evbub        DOUBLE       (i)  sm ele vel. bubble functions
\param **epbub        DOUBLE       (i)  sm ele pre. bubble functions
\param **efbub        DOUBLE       (i)  sm ele rhs bubble functions
\param  *vbubint      DOUBLE       (-)  vel. bubble fun. at int. p.
\param  *vbubderxy    DOUBLE       (-)  vel. bub. fun. der. at int. p.
\param  *pbubint      DOUBLE       (-)  pre. bubble fun. at int. p.
\param  *pbubderxy    DOUBLE       (-)  pre. bub. fun. der. at int. p.
\param  *velint       DOUBLE       (-)  vel at integration point
\param **vderxy       DOUBLE       (-)  global vel. derivatives
\param  *smvelint     DOUBLE       (-)  sm vel at integration point
\param **smvderxy     DOUBLE       (-)  sm global vel. derivatives
\param  *smpreint     DOUBLE       (-)  sm pre at integration point
\param **smpderxy     DOUBLE       (-)  sm global pre. derivatives
\param  *smfint       DOUBLE       (-)  sm rhs at integration point
\param **smfderxy     DOUBLE       (-)  sm global rhs. derivatives
\param **wa1	      DOUBLE       (-)  working array
\param **wa2	      DOUBLE       (-)  working array
\param  *dulinf       DOUBLE	   (o)  vel error in Linf-norm
\param  *dplinf       DOUBLE	   (o)  pre error in Linf-norm
\param  *dul2         DOUBLE	   (o)  vel error in L2-norm
\param  *dpl2         DOUBLE	   (o)  pre error in L2-norm
\param  *duh1         DOUBLE	   (o)  vel error in H1-norm
\param  *ul2          DOUBLE	   (o)  analytical vel in L2-norm
\param  *pl2          DOUBLE	   (o)  analytical pre in L2-norm
\param  *uh1          DOUBLE	   (o)  analytical vel in H1-norm
\return void                                                   

------------------------------------------------------------------------*/
void f3_bubinterr(FLUID_DATA      *data,     
                  FLUID_DYNAMIC   *fdyn, 
	          ELEMENT         *ele,     
	          FLUID_DYN_ML    *mlvar, 
	          FLUID_ML_SMESH  *submesh, 
   	          DOUBLE         **smxyze, 
	          DOUBLE         **smxyzep, 
	          DOUBLE          *funct,   
	          DOUBLE         **deriv,   
	          DOUBLE         **deriv2,   
	          DOUBLE         **xjm,     
	          DOUBLE         **derxy,   
	          DOUBLE          *smfunct,   
	          DOUBLE         **smderiv,   
	          DOUBLE         **smderiv2,   
	          DOUBLE         **smxjm,     
	          DOUBLE         **smderxy,   
	          DOUBLE         **evel,  
	          DOUBLE          *epre,
	          DOUBLE         **evbub,   
	          DOUBLE         **epbub,   
	          DOUBLE         **efbub,   
                  DOUBLE          *vbubint,    
                  DOUBLE         **vbubderxy,  
                  DOUBLE         **pbubint,    
                  DOUBLE        ***pbubderxy,  
	          DOUBLE          *velint,  
	          DOUBLE         **vderxy,  
                  DOUBLE          *smvelint,  
                  DOUBLE         **smvderxy,  
                  DOUBLE          *smpreint,  
                  DOUBLE         **smpderxy,  
                  DOUBLE          *smfint,    
                  DOUBLE         **smfderxy,  
	          DOUBLE         **wa1,     
	          DOUBLE         **wa2,
		  DOUBLE          *dulinf,  
		  DOUBLE          *dplinf,  
		  DOUBLE          *dul2,  
		  DOUBLE          *dpl2,  
		  DOUBLE          *duh1,  
		  DOUBLE          *ul2,  
		  DOUBLE          *pl2,  
		  DOUBLE          *uh1);
		    
/************************************************************************
 | f3_mlservice.c                                                            |
 ************************************************************************/
/*!--------------------------------------------------------------------- 
\brief set all arrays for large-scale element for multi-level fluid3

<pre>                                                       gravem 07/03

In this routine, the element velocities, pressure and external loads 
at time steps (n) and (n+1) are set. 

</pre>
\param   *dynvar   FLUID_DYN_CALC  (i)
\param   *ele      ELEMENT	   (i)    actual element
\param  **eveln    DOUBLE	   (o)    ele vels at time step n
\param  **evel     DOUBLE	   (o)    ele vels at time step n+1
\param   *epren    DOUBLE	   (o)    ele pres at time step n
\param   *epre     DOUBLE	   (o)    ele pres at time step n+1
\param   *edeadn   DOUBLE          (o)    ele dead load at time step n 
\param   *edead    DOUBLE          (o)    ele dead load at time step n+1 
\param   *hasext   INT             (o)    flag for external loads
\return void                                                                       

------------------------------------------------------------------------*/
void f3_lsset(FLUID_DYN_CALC  *dynvar, 
	      ELEMENT	      *ele,	
              DOUBLE	     **eveln,	 
	      DOUBLE	     **evel, 
	      DOUBLE	      *epren,
	      DOUBLE	      *epre,
	      DOUBLE	      *edeadn,
	      DOUBLE	      *edead,
	      INT	      *hasext);
/*!--------------------------------------------------------------------- 
\brief set all arrays for (sub-)submesh element for multi-level fluid3

<pre>                                                       gravem 07/03

In this routine, the (sub-)submesh element coordinates, topology array
and location matrix are set. 

</pre>
\param   *smesh    FLUID_ML_SMESH  (i)
\param   *ele      ELEMENT	   (i)    actual large-scale element
\param   *smlme    INT   	   (o)    (sub-)submesh location matrix
\param   *smitope  INT   	   (o)    (sub-)submesh topology array
\param   *smxyze   DOUBLE   	   (o)    (sub-)submesh coordinates
\param   *smxyzep  DOUBLE   	   (o)    (s)sm parent domain coordinates
\param    iele     INT             (i)    actual submesh element number
\param    flag     INT             (i)    flag: submesh or sub-submesh?
\return void                                                                       

------------------------------------------------------------------------*/
void f3_smset(FLUID_ML_SMESH  *smesh, 
	      ELEMENT	      *ele,	
              INT 	      *smlme,	 
	      INT	      *smitope, 
	      DOUBLE	     **smxyze,
	      DOUBLE	     **smxyzep,
	      INT	       iele,
	      INT	       flag);
/*!--------------------------------------------------------------------- 
\brief set all bubble functions for submesh element for fluid3

<pre>                                                       gravem 07/03

In this routine, the submesh element velocity, pressure and rhs bubble 
functions are set. 

</pre>
\param   *mlvar    FLUID_DYN_ML    (i)
\param   *submesh  FLUID_ML_SMESH  (i)
\param   *ele      ELEMENT	   (i)    actual large-scale element
\param   *smlme    INT   	   (i)    submesh location matrix
\param  **evbub    DOUBLE   	   (o)    submesh element velocity bubble
\param  **epbub    DOUBLE   	   (o)    submesh element pressure bubble
\param  **efbub    DOUBLE   	   (o)    submesh element rhs bubble
\param    flag     INT             (i)    flag: time step (n) or (n+1)?
\return void                                                                       

------------------------------------------------------------------------*/
void f3_bubset(FLUID_DYN_ML    *mlvar,  
               FLUID_ML_SMESH  *submesh, 
	       ELEMENT	       *ele,	
               INT	       *smlme,
	       DOUBLE         **evbub,
	       DOUBLE         **epbub,
	       DOUBLE         **efbub,
	       INT	        flag);
/*!--------------------------------------------------------------------- 
\brief set all arrays for integration of sub-submesh element for fluid3

<pre>                                                       gravem 07/03

In this routine, the sub-submesh element arrays for the elementwise
integration of a normalized bubble function are set. 

</pre>
\param   *ssmesh   FLUID_ML_SMESH  (i)
\param   *ele      ELEMENT	   (i)    actual large-scale element
\param   *sslme    INT   	   (o)    sub-submesh location matrix
\param   *ssitope  INT   	   (o)    sub-submesh topology array
\param   *ssxyze   DOUBLE   	   (o)    sub-submesh coordinates
\param   *ebub     DOUBLE   	   (o)    sub-submesh element bubble
\param    iele     INT             (i)    actual sub-submesh ele. number
\return void                                                                       

------------------------------------------------------------------------*/
void f3_ssset(FLUID_ML_SMESH  *ssmesh, 
	      ELEMENT	      *ele,	
              INT 	      *sslme,	 
	      INT	      *ssitope, 
	      DOUBLE	     **ssxyze,
	      DOUBLE          *ebub,
	      INT	       iele);
/*!--------------------------------------------------------------------- 
\brief copy submesh solution to element array for fluid3

<pre>                                                       gravem 07/03

In this routine, the submesh solution is copied to the respective
element array. 

</pre>
\param   *smrhs    DOUBLE          (i)    submesh solution array
\param   *ele      ELEMENT	   (o)    actual large-scale element
\param    numeq    INT             (i)    number of equations
\param    numrhs   INT             (i)    number of rhs
\return void                                                                       

------------------------------------------------------------------------*/
void f3_smcopy(DOUBLE  **smrhs,
               ELEMENT  *ele,
               INT       numeq,
               INT       numrhs);      
/*!--------------------------------------------------------------------- 
\brief copy submesh element solution at (n+1) to (n) for fluid3

<pre>                                                       gravem 07/03

In this routine, the submesh element solution array at time step (n+1)
is copied to the respective array at time step (n). 

</pre>
\param   *ele      ELEMENT	   (i/o)  actual large-scale element
\param    numeq    INT             (i)    number of equations
\param    numrhs   INT             (i)    number of rhs
\return void                                                                       

------------------------------------------------------------------------*/
void f3_smcopy2(ELEMENT  *ele,
                INT       numeq,
                INT       numrhs);      
/*!--------------------------------------------------------------------- 
\brief routine to calculate small-scale pressure at int. p. for fluid3

<pre>                                                       gravem 07/03
				      
</pre>
\param   *smpreint  DOUBLE     (o)   small-scale pressure at int. point
\param  **pbubint   DOUBLE     (i)   pressure bubble functions at int. p.
\param   *epre      DOUBLE     (i)   pressure at large-scale ele. nodes
\param    iel	    INT        (i)   number of large-scale ele. nodes 
\return void                                                                       

------------------------------------------------------------------------*/
void f3_smprei(DOUBLE  *smpreint,     
               DOUBLE **pbubint,    
	       DOUBLE  *epre,     
	       INT      iel); 
/*!--------------------------------------------------------------------- 
\brief routine to calculate s-s pressure deriv. at int. p. for fluid3

<pre>                                                       gravem 07/03
				      
In this routine, the derivatives of the small-scale pressure w.r.t x/y/z 
are calculated.
				      
</pre>
\param   *smpderxy  DOUBLE     (o)   s-s pressure deriv. at int. point
\param  **pbubderxy DOUBLE     (i)   pre. bubble fun. deriv. at int. p.
\param   *epre      DOUBLE     (i)   pressure at large-scale ele. nodes
\param    iel	    INT        (i)   number of large-scale ele. nodes 
\return void                                                                       

------------------------------------------------------------------------*/
void f3_smpder(DOUBLE  **smpderxy,     
               DOUBLE ***pbubderxy,    
	       DOUBLE   *epre,     
	       INT       iel); 
/*!--------------------------------------------------------------------- 
\brief routine to calc. 2nd s-s pressure deriv. at int. p. for fluid3

<pre>                                                       gravem 07/03
				      
In this routine, the 2nd derivatives of the small-scale pressure w.r.t 
x/y/z are calculated.
				      
</pre>
\param  **smpderxy2  DOUBLE     (o)   2nd s-s pre. deriv. at int. point
\param ***pbubderxy2 DOUBLE     (i)   2nd pre. bub. fun. deriv. at i. p.
\param   *epre       DOUBLE     (i)   pressure at large-scale ele. nodes
\param    iel	     INT        (i)   number of large-scale ele. nodes 
\return void                                                                       

------------------------------------------------------------------------*/
void f3_smpder2(DOUBLE  **smpderxy2,    
                DOUBLE ***pbubderxy2,    
	        DOUBLE   *epre,      
	        INT       iel);
		 
/************************************************************************
 | f3_mlstabmat.c                                                            |
 ************************************************************************/
/*!---------------------------------------------------------------------                                         
\brief evaluate stabilization part of sm stiffness matrix SMK for fluid3

<pre>                                                       gravem 07/03

In this routine, the stabilization part of the submesh stiffness matrix 
SMK is calculated.

</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *mlvar     FLUID_DYN_ML    (i)
\param **smestif   DOUBLE	   (i/o)  submesh ele stiffness matrix
\param  *velint    DOUBLE	   (i)    velocity at int. point
\param **vderxy    DOUBLE	   (i)    global vel. deriv. at int. p.
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun. 
\param **smderxy2  DOUBLE	   (i)    sm 2nd global deriv. of sh. fun. 
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   visc      DOUBLE	   (i)    fluid viscosity	     
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\param   ihoelsm   INT  	   (i)	  flag for higher-order elements
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabsmk(FLUID_DYN_CALC  *dynvar,
	           FLUID_DYN_ML    *mlvar, 
		   DOUBLE         **smestif,  
		   DOUBLE          *velint, 
		   DOUBLE         **vderxy, 
		   DOUBLE          *smfunct,  
		   DOUBLE         **smderxy,  
		   DOUBLE         **smderxy2, 
		   DOUBLE           fac,    
		   DOUBLE           visc,   
		   INT              smiel,    
                   INT              ihoelsm);
/*!---------------------------------------------------------------------                                         
\brief evaluate stabilization part of submesh mass matrix SMM for fluid3

<pre>                                                       gravem 07/03

In this routine, the stabilization part of the submesh mass matrix SMM 
is calculated.

</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *mlvar     FLUID_DYN_ML    (i)
\param **smemass   DOUBLE	   (i/o)  submesh element mass matrix
\param  *velint    DOUBLE	   (i)    velocity at int. point
\param **vderxy    DOUBLE	   (i)    global vel. deriv. at int. p.
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun. 
\param **smderxy2  DOUBLE	   (i)    sm 2nd global deriv. of sh. fun. 
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   visc      DOUBLE	   (i)    fluid viscosity	     
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\param   ihoelsm   INT  	   (i)	  flag for higher-order elements
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabsmm(FLUID_DYN_CALC  *dynvar,
	           FLUID_DYN_ML    *mlvar, 
		   DOUBLE         **smemass,  
		   DOUBLE          *velint, 
		   DOUBLE         **vderxy, 
		   DOUBLE          *smfunct,  
		   DOUBLE         **smderxy,  
		   DOUBLE         **smderxy2, 
		   DOUBLE           fac,    
		   DOUBLE           visc,   
		   INT              smiel,    
                   INT              ihoelsm);
		   
/************************************************************************
 | f3_mlstabpar.c                                                            |
 ************************************************************************/
/*!--------------------------------------------------------------------- 
\brief routine to calculate submesh stability parameter for fluid2                

<pre>                                                       gravem 07/03   
  									 
In this routine, the stability parameter for the actual submesh element
is calculated.

</pre>

\param   *ele         ELEMENT	      (i)    actual element
\param   *dynvar      FLUID_DYN_CALC  (i)
\param   *mlvar       FLUID_DYN_ML    (i/o)
\param   *velint      DOUBLE	      (i)    l-s velocity at ele. center
\param    visc        DOUBLE	      (i)    viscosity
\param    iel         INT	      (i)    number of nodes	     
\param	  ntyp        INT	      (i)    element type
\return void                                                                       

------------------------------------------------------------------------*/ 
void f3_smstabpar(ELEMENT         *ele,      
		  FLUID_DYN_CALC  *dynvar,
		  FLUID_DYN_ML    *mlvar,
		  DOUBLE	  *velint,  
		  DOUBLE	   visc,    
		  INT		   iel,     
		  INT		   ntyp);
/*!--------------------------------------------------------------------- 
\brief routine to calculate submesh subgrid viscosity for fluid3                

<pre>                                                       gravem 07/03   
  									 
In this routine, the subgrid viscosity for the actual submesh element
is calculated.

</pre>

\param   *ele         ELEMENT	      (i)    actual element
\param   *mlvar       FLUID_DYN_ML    (i/o)
\param   *velint      DOUBLE	      (i)    l-s velocity at ele. center
\param  **vderxy      DOUBLE	      (i)    l-s vel. der. at ele. center
\param    visc        DOUBLE	      (i)    viscosity
\param    iel         INT	      (i)    number of nodes	     
\param	  ntyp        INT	      (i)    element type
\return void                                                                       

------------------------------------------------------------------------*/ 
void f3_smsgvisc(ELEMENT         *ele,
                 FLUID_DYN_ML    *mlvar,
		 DOUBLE 	 *velint,  
		 DOUBLE 	**vderxy,  
		 DOUBLE 	  visc,    
		 INT		  iel,     
		 INT		  ntyp);
		 
/************************************************************************
 | f3_mlsubmesh.c                                                            |
 ************************************************************************/
/*!---------------------------------------------------------------------                                         
\brief creation of (sub-)submesh on parent domain for fluid3

<pre>                                                       gravem 07/03

In this routine, the creation of submesh or sub-submesh on the parent 
domain is performed, i.e. nodal coordinates on the parent domain, 
id-array and ien-array are established.
			     
</pre>   
\param  *smesh      FLUID_ML_SMESH   (i/o)  
\param   xele       INT              (i)    number of elements in x-dir.  
\param   yele       INT              (i)    number of elements in y-dir.  
\param   zele       INT              (i)    number of elements in z-dir.  
\param   order      INT              (i)    polyn. interpolation order  
\param   flag       INT              (i)    flag: submesh or sub-submesh?  
\return void 

------------------------------------------------------------------------*/
void f3_pdsubmesh(FLUID_ML_SMESH *smesh,
                  INT             xele,
                  INT             yele,
                  INT             zele,
                  INT             order,  
		  INT             flag);
/*!---------------------------------------------------------------------                                         
\brief creation of (sub-)submesh on individual element for fluid3

<pre>                                                       gravem 07/03

In this routine, the creation of submesh or sub-submesh on an individual 
element is performed, i.e. the nodal coordinates on this particular 
element are established.
			     
</pre>   
\param  *ele        ELEMENT          (i/o)  
\param  *smesh      FLUID_ML_SMESH   (i)  
\param   flag       INT              (i)    flag: submesh or sub-submesh?  
\return void 

------------------------------------------------------------------------*/
void f3_elesubmesh(ELEMENT        *ele,
                   FLUID_ML_SMESH *smesh,
		   INT             flag);
		   
/************************************************************************
 | f3_mltimerhs.c                                                       |
 ************************************************************************/
/*!---------------------------------------------------------------------                                         
\brief evaluate Galerkin part of submesh "Time" force vector for fluid3

<pre>                                                       gravem 07/03

In this routine, the Galerkin part of the submesh "Time" force vector
on the rhs is calculated.

</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *mlvar     FLUID_DYN_ML    (i)
\param **smetfor   DOUBLE	   (i/o)  submesh element time force vec.
\param  *velintn   DOUBLE          (i)    velocity at int. point at n
\param  *velintnt  DOUBLE          (i)    'temporal' vel. at int. p. at n
\param  *velintnc  DOUBLE          (i)    'convective' vel. at int. p. at n
\param **vderxyn   DOUBLE          (i)    global vel. derivatives at n
\param **vderxync  DOUBLE          (i)    global 'convective' vel. der. at n
\param **vderxynv  DOUBLE          (i)    global 'viscous' vel. der. at n
\param **vderxy2n  DOUBLE          (i)    2nd global vel. derivatives at n
\param  *pderxyn   DOUBLE          (i)    global pres. derivatives at n
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun.
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   visc 	   DOUBLE	   (i)    physical viscosity	      
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\param   iel	   INT  	   (i)	  number of nodes of l-s element
\param   ihoelsm   INT  	   (i)	  flag for higher-order sm ele.
\param   ihoel	   INT  	   (i)	  flag for higher-order l-s ele.
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calsmft(FLUID_DYN_CALC  *dynvar,
	        FLUID_DYN_ML	*mlvar, 
		DOUBLE         **smetfor,  
		DOUBLE  	*velintn, 
		DOUBLE  	*velintnt, 
		DOUBLE  	*velintnc, 
		DOUBLE         **vderxyn, 
		DOUBLE         **vderxync, 
		DOUBLE         **vderxynv, 
		DOUBLE         **vderxy2n, 
		DOUBLE          *pderxyn, 
		DOUBLE  	*smfunct,  
		DOUBLE         **smderxy,  
		DOUBLE  	 fac,	 
		DOUBLE  	 visc,   
		INT		 smiel,    
		INT		 iel,    
                INT		 ihoel);
/*!---------------------------------------------------------------------                                         
\brief evaluate stabilization part of sm "Time" force vector for fluid3

<pre>                                                       gravem 07/03

In this routine, the stabilization part of the submesh "Time" force 
vector on the rhs is calculated.

</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *mlvar     FLUID_DYN_ML    (i)
\param **smetfor   DOUBLE	   (i/o)  submesh element time force vec.
\param  *velintn   DOUBLE          (i)    velocity at int. point at n
\param  *velintnt  DOUBLE          (i)    'temporal' vel. at int. p. at n
\param  *velintnc  DOUBLE          (i)    'convective' vel. at int. p. at n
\param **vderxyn   DOUBLE          (i)    global vel. derivatives at n
\param **vderxync  DOUBLE          (i)    global 'convective' vel. der. at n
\param **vderxynv  DOUBLE          (i)    global 'viscous' vel. der. at n
\param **vderxy2n  DOUBLE          (i)    2nd global vel. derivatives at n
\param  *pderxyn   DOUBLE          (i)    global pres. derivatives at n
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun.
\param **smderxy2   DOUBLE	   (i)    sm 2nd global deriv. of shape fun.
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   visc 	   DOUBLE	   (i)    physical viscosity	      
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\param   iel	   INT  	   (i)	  number of nodes of l-s element
\param   ihoelsm   INT  	   (i)	  flag for higher-order sm ele.
\param   ihoel	   INT  	   (i)	  flag for higher-order l-s ele.
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabsmft(FLUID_DYN_CALC  *dynvar,
	            FLUID_DYN_ML    *mlvar, 
		    DOUBLE         **smetfor,  
		    DOUBLE  	    *velintn, 
		    DOUBLE  	    *velintnt, 
		    DOUBLE  	    *velintnc, 
		    DOUBLE         **vderxyn, 
		    DOUBLE         **vderxync, 
		    DOUBLE         **vderxynv, 
		    DOUBLE         **vderxy2n, 
		    DOUBLE          *pderxyn, 
	 	    DOUBLE  	    *smfunct,  
		    DOUBLE         **smderxy,  
		    DOUBLE         **smderxy2, 
		    DOUBLE  	     fac,	 
		    DOUBLE  	     visc,   
		    INT		     smiel,    
		    INT		     iel,    
                    INT		     ihoelsm,
                    INT		     ihoel);
		    
/************************************************************************
 | f3_mlvmmrhs.c                                                        |
 ************************************************************************/
/*!---------------------------------------------------------------------                                         
\brief evaluate Galerkin part of submesh "VMM" force vector for fluid3

<pre>                                                       gravem 07/03

In this routine, the Galerkin part of the submesh "VMM" force vector on
the rhs is calculated.

</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *mlvar     FLUID_DYN_ML    (i)
\param **smevfor   DOUBLE	   (i/o)  submesh element vmm force vec.
\param  *velint    DOUBLE	   (i)    velocity at int. point
\param **vderxy    DOUBLE	   (i)    global vel. deriv. at int. p.
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param  *funct     DOUBLE	   (i)    natural shape functions
\param **derxy     DOUBLE	   (i)    global deriv. of shape fun.
\param **derxy2    DOUBLE	   (i)    2nd global deriv. of shape fun.
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   visc      DOUBLE	   (i)    fluid viscosity	     
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\param   iel	   INT  	   (i)	  number of nodes of l-s element
\param   ihoel	   INT  	   (i)	  flag for higher-order l-s ele.
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calsmfv(FLUID_DYN_CALC  *dynvar,
	        FLUID_DYN_ML	*mlvar, 
		DOUBLE         **smevfor,  
		DOUBLE  	*velint, 
		DOUBLE         **vderxy, 
		DOUBLE  	*smfunct,  
		DOUBLE  	*funct,  
		DOUBLE         **derxy,  
		DOUBLE         **derxy2, 
		DOUBLE  	 fac,	 
		DOUBLE  	 visc,   
		INT		 smiel,    
		INT		 iel,    
                INT		 ihoel);
/*!---------------------------------------------------------------------                                         
\brief evaluate stabilization part of sm "VMM" force vector for fluid3

<pre>                                                       gravem 07/03

In this routine, the stabilization part of the submesh "VMM" force 
vector on the rhs is calculated.

</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *mlvar     FLUID_DYN_ML    (i)
\param **smevfor   DOUBLE	   (i/o)  submesh element vmm force vec.
\param  *velint    DOUBLE	   (i)    velocity at int. point
\param **vderxy    DOUBLE	   (i)    global vel. deriv. at int. p.
\param  *smfunct   DOUBLE	   (i)    sm natural shape functions
\param **smderxy   DOUBLE	   (i)    sm global deriv. of shape fun.
\param **smderxy2  DOUBLE	   (i)    sm 2nd glo. der. of shape fun.
\param  *funct     DOUBLE	   (i)    natural shape functions
\param **derxy     DOUBLE	   (i)    global deriv. of shape fun.
\param **derxy2    DOUBLE	   (i)    2nd global deriv. of shape fun.
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   visc      DOUBLE	   (i)    fluid viscosity	     
\param   smiel	   INT  	   (i)	  number of nodes of sm element
\param   iel	   INT  	   (i)	  number of nodes of l-s element
\param   ihoelsm   INT  	   (i)	  flag for higher-order sm ele.
\param   ihoel	   INT  	   (i)	  flag for higher-order l-s ele.
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calstabsmfv(FLUID_DYN_CALC  *dynvar,
	            FLUID_DYN_ML    *mlvar, 
		    DOUBLE         **smevfor,  
		    DOUBLE  	    *velint, 
		    DOUBLE         **vderxy, 
		    DOUBLE  	    *smfunct,  
		    DOUBLE         **smderxy,  
		    DOUBLE         **smderxy2, 
		    DOUBLE  	    *funct,  
		    DOUBLE         **derxy,  
		    DOUBLE         **derxy2, 
		    DOUBLE  	     fac,	 
		    DOUBLE  	     visc,   
		    INT		     smiel,    
		    INT		     iel,    
                    INT		     ihoelsm,
                    INT		     ihoel);
/*!---------------------------------------------------------------------                                         
\brief evaluate bubble part of large-scale rhs (vel. dofs) for fluid3

<pre>                                                       gravem 07/03

In this routine, the bubble part of the large-scale rhs for the
velocity dofs is calculated.

</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *mlvar     FLUID_DYN_ML    (i)
\param **eiforce   DOUBLE	   (i/o)  element iterative rhs
\param  *velint    DOUBLE	   (i)    velocity at int point
\param **vderxy    DOUBLE	   (i)    global velocity derivatives
\param  *funct     DOUBLE	   (i)    natural shape functions
\param **derxy     DOUBLE	   (i)    global deriv. of shape fun.
\param  *smfint    DOUBLE	   (i)    rhs bubble functions
\param **smfderxy  DOUBLE	   (i)    global deriv. of rhs bub. fun.
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   visc      DOUBLE	   (i)    fluid viscosity	     
\param   iel	   INT  	   (i)	  number of element nodes
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calbfv(FLUID_DYN_CALC  *dynvar,
	       FLUID_DYN_ML    *mlvar,
	       DOUBLE          *eiforce,   
	       DOUBLE          *velint, 
	       DOUBLE         **vderxy, 
	       DOUBLE          *funct,  
	       DOUBLE         **derxy,  
	       DOUBLE          *smfint,  
	       DOUBLE         **smfderxy,  
	       DOUBLE           fac,    
	       DOUBLE           visc,   
	       INT              iel);
/*!---------------------------------------------------------------------                                         
\brief evaluate bubble part of large-scale rhs (pre. dofs) for fluid2

<pre>                                                       gravem 07/03

In this routine, the bubble part of the large-scale rhs for the
pressure dofs is calculated.

</pre>
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **eiforce   DOUBLE	   (i/o)  element iterative rhs
\param  *funct     DOUBLE	   (i)    natural shape functions
\param **smfderxy  DOUBLE	   (i)    global deriv. of rhs bub. fun.
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   iel	   INT  	   (i)	  number of element nodes
\return void                                                                       

------------------------------------------------------------------------*/
void f3_calbfp(FLUID_DYN_CALC  *dynvar,
	       DOUBLE          *eiforce,   
	       DOUBLE          *funct,  
	       DOUBLE         **smfderxy,  
	       DOUBLE           fac,    
	       INT              iel);
