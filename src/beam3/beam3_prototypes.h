/*!----------------------------------------------------------------------
\file
\brief beam3 prototypes

*----------------------------------------------------------------------*/
#ifdef D_BEAM3

/*! 
\addtogroup BEAM3
*//*! @{ (documentation module open)*/

/*------------------------------------------------------------------------*/
/* RULE HOW TO ADD NEW FILES AND FUNCTIONS: 
   1.) THE FILENAMES ARE IN ALPHABETICAL ORDER !!!
   2.) FUNCTIONS ARE IN THE SAME ORDER LIKE IN THE FILE!!!
*/
/************************************************************************
 | b3_bop.c                                                             |
 ************************************************************************/
/*!----------------------------------------------------------------------
\brief calculates the B-Operator matrix for a 3d beam element

<pre>                                                              fh 02/03
This routine calculates the B-Operator matrix for a 3D-beam element with 3 
coordinates r,s,t. 

</pre>
\param **bop    DOUBLE  (o)   B-Operator
\param **deriv  DOUBLE  (i)   the derivatives of the shape functions
\param **func   DOUBLE  (i)   the shape functions
\param **xjm    DOUBLE  (i)   the Jacobian Matrix
\param **ijm    DOUBLE  (i)   the Inverse of the Jacobian Matrix
\param **V      DOUBLE  (i)   the vectors for local directions s,t
\param s        DOUBLE  (i)   s coordinate
\param t        DOUBLE  (i)   t coordinate
\param b        DOUBLE  (i)   width of beam cross section
\param a        DOUBLE  (i)   height of beam cross section
\param iel      INT     (i)   number of element nodes
\param init     INT     (i)   flag if initialization (init=1) or not
               

\warning This routine is not working yet, no call in other routines
\return void                                               
\sa calling: ---; called by: b3_cal_ele() 

*----------------------------------------------------------------------*/
void b3_boplin3D(DOUBLE    **bop,    
                 DOUBLE    **deriv,  
                 DOUBLE     *func,   
	         DOUBLE    **xjm,    
                 DOUBLE    **ijm,    
	         DOUBLE    **V,      
	         DOUBLE      s,      
	         DOUBLE      t,      
	         DOUBLE      b,      
	         DOUBLE      a,      
                 INT         iel,
	         INT         init);    
	    
/*!----------------------------------------------------------------------
\brief calculates the B-Operator matrix for a 1d spatial beam element

<pre>                                                              fh 01/03
This routine calculates the B-Operator matrix for a spatial 1D-Timoshenko
beam element w.r.t. coordinate r.

</pre>
\param **bop    DOUBLE  (o)   B-Operator
\param **deriv  DOUBLE  (i)   the derivatives of the shape functions
\param *func    DOUBLE  (i)   the shape functions
\param iel      INT     (i)   number of element nodes
\param l2       DOUBLE  (i)   dx/dr transformation x -> r               


\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; 
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_boplin(DOUBLE    **bop,    
               DOUBLE    **deriv,  
               DOUBLE     *func,   
               INT         iel,    
	       DOUBLE      l2);


/************************************************************************
 | b3_cal_force.c                                                       |
 ************************************************************************/
/*!----------------------------------------------------------------------
\brief calculates the internal forces of a Bernoulli beam element

<pre>                                                              fh 10/02
This routine calculates the internal force vector for a spatial 
1d-Bernoulli-beam element (direct stiffness method)

</pre>
\param *ele      ELEMENT  (i/o)  actual element
\param **estif   DOUBLE    (i)  local element stiffness matrix
\param *force    DOUBLE    (i)  local element load vector
\param init      INT       (i)  flag if init (1) or calculation	(2,3)
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: beam3() , b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_cal_force(ELEMENT  *ele,     
                  DOUBLE  **estif,     
	          DOUBLE   *force,     
	          INT	    init);
		   
/*!----------------------------------------------------------------------
\brief calculates the internal forces of a Timoshenko beam element

<pre>                                                              fh 01/03
This routine calculates the internal force vector for a spatial 
1d-Timoshenko-beam element (Finite Element method)

</pre>
\param *ele      ELEMENT  (i/o) actual element
\param *force    DOUBLE    (i)  stress vector at actual gauss point
\param lr        INT       (i)  number of actual gauss point
\param option    INT       (i)  =0 -> force_GP, =1 -> force_ND
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_cal_forcelin(ELEMENT  *ele,     
		     DOUBLE   *force,    
	             INT       lr,      
		     INT       option);
		      
/*!----------------------------------------------------------------------
\brief calculates the internal forces of a spatial 3D Timoshenko beam element

<pre>                                                              fh 02/03
This routine calculates the internal force vector for a spatial 
3D-Timoshenko-beam element (Finite Element method)

</pre>
\param *ele      ELEMENT  (i/o) actual element
\param *stress   DOUBLE    (i)  stress vector at actual integration point
\param facl      DOUBLE    (i)  weight for actual lobatto point
\param lr        INT       (i)  number of actual gauss point
\param s         DOUBLE    (i)  coordinate s of actual lobatto point
\param t         DOUBLE    (i)  coordinate t of actual lobatto point
\param option    INT       (i)  =0 -> force_GP, =1 -> force_ND
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_cal_forcelin3D(ELEMENT  *ele,     
		       DOUBLE   *stress,    
	               DOUBLE    facl,     
	               INT       lr,      
	               DOUBLE    s,
	       	       DOUBLE    t,
		       INT       option);

/*!----------------------------------------------------------------------
\brief calculates the strains at the actual integration point

<pre>                                                              fh 03/03
This routine calculates the strains at the actual integration point of
the actual beam element.

</pre>
\param *strain   DOUBLE    (o)  strain vector
\param pv        DOUBLE    (i)  possions ratio
\param *edisp    DOUBLE    (i)  vector of global element displacements
\param **bop     DOUBLE    (i)  B-Operator matrix
\param numeps    INT       (i)  number of strains per node
\param nedof     INT       (i)  number of dofs per element

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_cal_eps(DOUBLE    *strain,
                DOUBLE     pv,
		DOUBLE    *edisp,     
	        DOUBLE   **bop,     
		INT        numeps,
		INT        nedof);

/*!----------------------------------------------------------------------
\brief calculates the nodal element displacements

<pre>                                                              fh 03/03
This routine calculates the nodal element displacements.

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *edisp    DOUBLE    (i)  vector of global element displacements


\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_edisp(ELEMENT   *ele,
              DOUBLE    *edisp);

/*!----------------------------------------------------------------------
\brief calculates the nodal incremental element displacements

<pre>                                                              fh 05/03
This routine calculates the nodal incremental element displacements.

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *edisp    DOUBLE    (i)  vector of global element displacements


\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_dispinc(ELEMENT   *ele,
                DOUBLE    *edisp);

/*!----------------------------------------------------------------------
\brief extrapolates results from gauss points to the nodes

<pre>                                                              fh 09/02
This routine extrapolates the internal force values from the gauss point
to the nodes

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *data     B3_DATA   (o)  data for integration parameters
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele()

*-----------------------------------------------------------------------*/
void b3_exforce(ELEMENT   *ele,
                B3_DATA   *data);
		 

/************************************************************************
 | b3_cal_ele.c                                                       |
 ************************************************************************/
/*!----------------------------------------------------------------------
\brief performs all calculation for beam elements

<pre>                                                              fh 09/02
This routine performs all calculation for beam elements (stiffness matrix,
element load vector, internal forces)

</pre>
\param *ele            ELEMENT  (i/o) actual element
\param *data           B3_DATA   (o)  data set for gauss points
\param *mat            MATERIAL  (o)  actual material
\param **estif_global  ARRAY     (o)  global element stiffness matrix 
\param *force          DOUBLE    (o)  global force vector
\param *action         CALC_ACTION (i)  action to do
\param init            INT       (i)  initialization (1) or calculation (2,3)
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   b3_cal_sec() , b3_cal_lst() , b3_con_dof() , b3_cal_trn() ,
               b3_trans_stf() , b3_funct_deriv() , b3_jaco() , b3_boplin() ,
	       b3_boplin3D() , b3_call_mat() , b3_keku() , b3_load() ,
	       b3_loadlin() , b3_cal_stress() , b3_cal_stresslin(), b3_edisp(),
	       b3_exforce() , b3_cal_eps()
    called by: beam3()

*----------------------------------------------------------------------*/
void b3_cal_ele(ELEMENT     *ele, 
                B3_DATA     *data, 
                MATERIAL    *mat,
                ARRAY       *estif_global, 
                DOUBLE      *force,  
                CALC_ACTION *action,
		INT          init);

	       
/************************************************************************
 | b3_cal_fext.c                                                        |
 ************************************************************************/     
 /*!----------------------------------------------------------------------
\brief calculates the global element force vector for a spatial Bernoulli 
beam element

<pre>                                                              fh 09/02
This routine calculates the global element force vector for a spatial 
1d-Bernoulli-beam element (direct stiffness method)

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *mat      MATERIAL  (i)  actual material
\param **stiff   DOUBLE   (i/o) local element stiffness matrix
\param *loadvec  DOUBLE    (o)  global element load vector
\param *hinge    INT       (i)  Hinge Code for actual beam element
\param calcstep  INT       (i)  flag for calculation step
               

\warning done for n nodes per element
\return void                                               
\sa calling:   b3_fext() , b3_con_load() , b3_cal_trn() 
    called by: beam3() , b3_cal_ele() 

*----------------------------------------------------------------------*/
void b3_load(ELEMENT   *ele,     
             MATERIAL  *mat,     
	     DOUBLE   **stiff,
	     DOUBLE    *loadvec,
	     INT       *hinge,  
	     INT        calcstep);

/*!----------------------------------------------------------------------
\brief calculates the element force vector for a spatial Timoshenko beam 
element

<pre>                                                              fh 01/03
This routine calculates the global element force vector for a spatial 1d-
Timoshenko-beam element (Finite element method)

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *mat      MATERIAL  (i)  actual material
\param *func     DOUBLE    (i)  the shape functions
\param r         DOUBLE    (i)  actual gauss point
\param fac       DOUBLE    (i)  weight of actual gauss point
\param *loadvec  DOUBLE    (o)  global element load vector
\param calcstep  INT       (i)  flag for calculation step
               

\warning done for n nodes per element
\return void                                               
\sa calling:   b3_fextlin() 
    called by: beam3() , b3_cal_ele()

*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
|								       |
|    r---1-----------------------------------------------2---r  lin2   |
|								       |
|								       |
|    r---1-----------------------3-----------------------2---r  lin3   |
|								       |
|								       |
|    r---1---------------3-------------------4-----------2---r  lin4   |
|								       |
*----------------------------------------------------------------------*/ 
void b3_loadlin(ELEMENT  *ele,     
                MATERIAL *mat,     
        	DOUBLE   *func,    
		DOUBLE    r,       
		DOUBLE    fac,     
		DOUBLE	 *loadvec, 
		INT       calcstep);

/*!----------------------------------------------------------------------
\brief calculates the local element force vector for a spatial Bernoulli
beam element

<pre>                                                              fh 09/02
This routine calculates the local element force vector (line load)for a 
spatial 1d-Bernoulli-beam element (direct stiffness method)

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *mat      MATERIAL  (i)  actual material
\param *eload    DOUBLE    (o)  local  element load vector
               

\warning There is nothing special in this routine
\return void                                               
\sa calling: ---;
    called by: b3_load() 

*----------------------------------------------------------------------*/
void b3_fext(ELEMENT     *ele,
             MATERIAL    *mat,
             DOUBLE      *eload);
	     
/*!----------------------------------------------------------------------
\brief calculates the local element force vector for a spatial Timoshenko
beam element

<pre>                                                              fh 01/03
This routine calculates the local element force vector for a spatial 1d-
Timoshenko-beam element (Finite element method)

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *mat      MATERIAL  (i)  actual material
\param *func     DOUBLE    (i)  the shape functions
\param r         DOUBLE    (i)  actual gauss point
\param fac       DOUBLE    (i)  weight of actual gauss point
\param *eload    DOUBLE    (o)  global element load vector
               

\warning done for n nodes per element
\return void                                               
\sa calling:   ---; 
    called by: b3_loadlin() 

*----------------------------------------------------------------------*/
void b3_fextlin(ELEMENT     *ele,  
                MATERIAL    *mat,  
                DOUBLE      *func, 
		DOUBLE       r,    
		DOUBLE       fac,  
		DOUBLE      *eload);
		
/*!----------------------------------------------------------------------
\brief condensates the local element force vector if necessary

<pre>                                                              fh 11/02
This routine condensates the local spatial element force vector if there
are any hinge conditions at the nodes

</pre>
\param **stiff   DOUBLE   (i/o) local element stiffness matrix
\param *loadvec  DOUBLE    (o)  global element load vector
\param *hc       INT       (i)  Hinge Code for actual beam element
\param iel       INT       (i)  number of nodes per element
               
\warning This routine is adopted from CARAT (advanced to n nodes per ele)
\return void                                               
\sa calling:   ---; 
    called by: b3_load()

*----------------------------------------------------------------------*/
void b3_con_load(DOUBLE **estif,
		 DOUBLE  *loadvec, 
                 INT     *hc,
		 INT      iel);

/*!----------------------------------------------------------------------
\brief condensates the local element force vector of a Timoshenko beam 
if necessary

<pre>                                                              fh 02/03
This routine condensates the local spatial element force vector of a 
Timoshenko beam if there are any hinge conditions at the end nodes of
the element

</pre>
\param *loadvec  DOUBLE   (i/o) local element load vector
\param *hc       INT       (i)  Hinge Code for actual beam element
\param iel       INT       (i)  number of nodes per element
               
\warning This routine is adopted from CARAT (advanced to n nodes per ele)
\return void                                               
\sa calling:   ---; 
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_con_loadlin(DOUBLE  *loadvec, 
                    INT     *hc,
		    INT      iel);

		 
/************************************************************************
 | b3_cal_stiff.c                                                       |
 ************************************************************************/
/*!----------------------------------------------------------------------
\brief calculates the element stiffness matrix of 1d Bernoulli beam

<pre>                                                              fh 09/02
This routine calculates the local element stiffness matrix of a spatial
1d-Bernoulli beam w.r.t the main axes of the cross section

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *mat      MATERIAL  (i)  actual material
\param **estif   DOUBLE    (o)  local element stiffness matrix
\param *hinge    INT       (i)  Hinge Code for actual beam element
               

\warning done for n nodes per element
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele() 

*----------------------------------------------------------------------*/
void b3_cal_lst(ELEMENT  *ele, 
                MATERIAL *mat, 
		DOUBLE  **estif, 
		INT      *hinge);

/*!----------------------------------------------------------------------
\brief calculates the main axis of the cross section

<pre>                                                              fh 09/02
This routine calculates the main axis of the cross section of the actual
element and the main moments of inertia

</pre>
\param *ele      ELEMENT  (i/O)  actual element
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele() 

*----------------------------------------------------------------------*/
void b3_cal_sec(ELEMENT *ele);

/*!----------------------------------------------------------------------
\brief calculates the transformation matrix of the actual element

<pre>                                                              fh 09/02
This routine calculates the transformation matrix (local -> global) of the
actual element

</pre>
\param *ele      ELEMENT   (i)  actual element
\param **A       DOUBLE    (O)  transformation matrix

\warning done for n nodes per element
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele() 

*----------------------------------------------------------------------*/
void b3_cal_trn(ELEMENT *ele, 
                DOUBLE **A);

/*!----------------------------------------------------------------------
\brief condensates the local element stiffness matrix if necessary

<pre>                                                              fh 09/02
This routine calculates the condensed local element stiffness matrix if
there are any hinge conditions

</pre>
\param **estif   DOUBLE   (i/o) local element stiffness matrix 
\param *hc       INT       (i)  hinge code of actual element
\param nedof     INT       (i)  number of dofs per element

\warning done for n nodes per element
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele() 

*----------------------------------------------------------------------*/
void b3_con_dof(DOUBLE **estif, 
                INT     *hc,
		INT	 nedof);

/*!----------------------------------------------------------------------
\brief condensates the local element stiffness matrix for warping dofs

<pre>                                                              fh 06/03
This routine does the statical condensation of the additional warping dofs

</pre>
\param **estif   DOUBLE   (i/o) local element stiffness matrix 
\param iel       INT       (i)  number of nodes per element
\param nedof     INT       (i)  number of dofs per element

\warning done for n nodes per element
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele() 

*----------------------------------------------------------------------*/
void b3_con_warp(DOUBLE **estif, 
		 INT      iel,
		 INT	  nedof);
		   
/*!----------------------------------------------------------------------
\brief transforms local to global stiffness matrix

<pre>                                                              fh 09/02
This routine calculates the transformation of the local to the global 
element stiffness matrix

</pre>
\param **K       DOUBLE    (i)  local element stiffness matrix 
\param **TRN     DOUBLE    (i)  element transformation matrix 
\param **estif   DOUBLE    (o)  global element stiffness matrix
\param nedof     INT       (i)  number of dofs of the element
\param calcstep  INT       (i)  flag for calculation step


\warning done for n nodes per element
\return void                                               
\sa calling:   ---;
    called by: beam3() , b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_trans_stf(DOUBLE **K, 
                  DOUBLE **TRN,
		  DOUBLE **estif,
		  INT      nedof,
		  INT	   calcstep);

/*!----------------------------------------------------------------------
\brief calculates element stiffness matrix of Timoshenko beam element

<pre>                                                              fh 09/02
This routine calculates the local element stiffness matrix of a spatial
Timoshenko beam (Finite Element method)

</pre>
\param **S       DOUBLE    (o)  local element stiffness matrix 
\param **bs      DOUBLE    (i)  B-Operator matrix
\param **d       DOUBLE    (i)  Constitutive Matrix
\param fac       DOUBLE    (i)  integration factor
\param nd        INT       (i)  total number of degrees of freedom of element
\param neps      INT       (i)  actual number of strain components = 3

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_keku(DOUBLE  **s,    
             DOUBLE  **bs,   
             DOUBLE  **d,    
             DOUBLE    fac,  
             INT       nd,   
             INT       neps);


/************************************************************************
 | b3_call_mat.c                                                        |
 ************************************************************************/
/*!----------------------------------------------------------------------
\brief selects proper material law

<pre>                                                              fh 10/02
This routine selects the proper material law for a spatial 1d-beam-element

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *mat      MATERIAL  (o)  actual material
\param *eps      DOUBLE    (i)  strain vector
\param **bop     DOUBLE    (i)  B-Operator matrix
\param **D       DOUBLE    (o)  constitutive matrix
\param  *stress  DOUBLE    (o)  stress at actual gauss point
\param ip        INT       (i)  actual integration point
\param istore    INT       (i)  flag for storing of new stresses
\param newval    INT       (o)  flag for evaluating of new stresses

\warning There is nothing special in this routine
\return void                                               
\sa calling:   b3_mat_linel() , b3_mat_plast_mises() 
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_call_mat(ELEMENT   *ele,  
                 MATERIAL  *mat,  
                 DOUBLE    *eps,
		 DOUBLE   **bop,    
                 DOUBLE   **d,      
                 DOUBLE    *stress,
		 INT        ip,
		 INT        istore,
		 INT        newval);

		 
/************************************************************************
 | b3_dirich.c                                                          |
 ************************************************************************/
/*!----------------------------------------------------------------------
\brief puts all dirichlet values to the nodal solution

<pre>                                                              fh 12/02
This routine puts all values for dirichlet conditions a priori to the
nodal displacement solution vector

</pre>
\param *actfield FIELD    (i/o)  actual field
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: beam3() 

*----------------------------------------------------------------------*/
void b3_setdirich(FIELD     *actfield);


/************************************************************************
 | b3_funcderiv.c                                                      |
 ************************************************************************/
/*!----------------------------------------------------------------------
\brief calculates the values of the shape functions and the correspondent
derivatives at point r

<pre>                                                              fh 10/02
This routine calculates the values of the shape functions and the
correspondent derivatives if necessary at the actual gauss point r

</pre>
\param *funct    DOUBLE    (o)  value for shape functions at point r
\param **deriv   DOUBLE    (o)  value for derivatives at point r
\param r         DOUBLE    (i)  actual gauss point coordinate
\param typ       DIS_TYP   (i)  LIN2, LIN3 or LIN4
\param option    INT       (i)  flag for derivatives (1) or not (0)
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_funct_deriv(DOUBLE     *funct, 
                    DOUBLE    **deriv, 
                    DOUBLE      r,     
                    DIS_TYP     typ,     
                    INT         option);
		    

/************************************************************************
 | b3_init.c                                                            |
 ************************************************************************/
/*!----------------------------------------------------------------------
\brief initializes the element

<pre>                                                              fh 09/02
This routine initializes the actual beam element

</pre>
\param *actpart  PARTITION (i)  actual partition
\param *mat      MATERIAL  (i)  actual material
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   b3_intg()
    called by: beam3() 

*----------------------------------------------------------------------*/
void b3_init(PARTITION *actpart,
             MATERIAL  *mat );
	     

/************************************************************************
 | b3_intg.c                                                            |
 ************************************************************************/
/*!----------------------------------------------------------------------
\brief sets all integration parameters for actual beam element

<pre>                                                              fh 09/02
This routine sets the gauss-point coordinates and the corresponding
weighting factors for the actual Timoshenko-beam element

</pre>
\param *ele      ELEMENT   (i)  actual element
\param *data     B3_DATA   (o)  data for integration parameters
\param option    INT       (i)  initialization (0) or calculation (1)
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele() , b3_init()

*-----------------------------------------------------------------------*/
void b3_intg(ELEMENT   *ele,
            B3_DATA   *data,
            INT        option);
	    

/************************************************************************
 | b3_jaco.c                                                            |
 ************************************************************************/
/*!----------------------------------------------------------------------
\brief calculates the Jacobian Matrix at point r,s,t and the corresponding
determinant

<pre>                                                              fh 02/03
This routine calculates the Jacobian Matrix at point r,s,t and the 
corresponding determinant of the actual Timoshenko beam element

</pre>
\param *funct    DOUBLE    (i)  values of shape functions at actual GP
\param **deriv   DOUBLE    (i)  values of derivatives at actual GP
\param **xjm     DOUBLE    (o)  Jacobian Matrix
\param **ijm     DOUBLE    (o)  Inverse of Jacobian Matrix
\param **V       DOUBLE    (i)  Vs, Vt directions of local vectors s,t
\param *det      DOUBLE    (o)  determinant of Jacobian Matrix
\param s         DOUBLE    (i)  coordinate in local s-direction  
\param t         DOUBLE    (i)  coordinate in local t-direction  
\param *ele      ELEMENT   (i)  actual element
\param iel       INT       (i)  number of element nodes
               

\warning This routine is not used yet.
\return void                                               
\sa calling:   ---;
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_jaco(DOUBLE     *funct, 
             DOUBLE    **deriv, 
             DOUBLE    **xjm,   
             DOUBLE    **ijm,   
	     DOUBLE    **V,     
	     DOUBLE     *det,   
             DOUBLE      s,
	     DOUBLE      t,
	     ELEMENT    *ele,   
             INT         iel);   


/************************************************************************
 | b3_mat_linel.c                                                       |
 ************************************************************************/	      	    	     		    
/*!----------------------------------------------------------------------
\brief calculates the constitutive matrix for a beam

<pre>                                                              fh 10/02
This routine calculates the linear elastic constitutive matrix for a 
spatial 1d-Timoshenko-beam element 

</pre>
\param *ele      ELEMENT   (i)  actual element
\param ym        DOUBLE    (i)  Youngs Modulus
\param pv        DOUBLE    (i)  Poissons ratio 
\param **d       DOUBLE    (o)  Constitutive matrix


\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_call_mat()

*----------------------------------------------------------------------*/
void b3_mat_linel(ELEMENT *ele,    
                  DOUBLE   ym,     
                  DOUBLE   pv,     
		  DOUBLE **d);      

		 		  		     		  		      		   
/************************************************************************
 | b3_mat_plast_mises.c                                                 |
 ************************************************************************/	
/*!----------------------------------------------------------------------
\brief calculates the constitutive matrix - forces for von Mises material

<pre>                                                              fh 04/03
This routine calculates the constitutive matrix and the actual stresses
within a plastic calculation for von Mises material

</pre>
\param ym              DOUBLE      (i)  Youngs Modulus
\param pv              DOUBLE      (i)  Possions ratio
\param ALFAT           DOUBLE      (i)  alpha T
\param sigy            DOUBLE      (i)  1D yield stress
\param hard            DOUBLE      (i)  1D hardening modulus
\param gf              DOUBLE      (i)  Crack energy (Gf)
\param *ele            ELEMENT    (i/o) actual element
\param *strain         DOUBLE      (i)  strain vector of actual IP
\param ip              INT         (i)  actual integration point number
\param *stress         DOUBLE      (o)  stress vector of actual IP
\param **d             DOUBLE      (o)  constitutive matrix
\param istore          INT         (i)  flag to control storing of new val to WA
\param newval          INT         (i)  flog to control evaluation of new stresses
\param init            INT         (i)  initialization (1) or calculation (2,3)
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   b3_blowup_stress() , b3_compr_stress() , b3_blowup_strain() ,
               b3_compr_strain() , b3_compr_Cep() , b3_kronecker() , 
	       b3_cal_volCel() , b3_cal_devCel() , b3_cal_Cel() , b3_condense()
    called by: beam3() , b3_call_mat()

*----------------------------------------------------------------------*/
void b3_mat_plast_mises(DOUBLE    ym,
                        DOUBLE    pv,
                        DOUBLE    ALFAT,
                        DOUBLE    sigy,
                        DOUBLE    hard,
                        DOUBLE    gf,
                        ELEMENT  *ele,
                        DOUBLE   *strain,
			INT       ip,
                        DOUBLE   *stress,       
                        DOUBLE  **d,
                        INT       istore,
                        INT       newval,
			INT       init);

/*!----------------------------------------------------------------------
\brief calculates stress tensor (3*3) out of stress vector (6*1)

<pre>                                                              fh 03/03
This routine calculates the stress tensor (3*3) out of stress vector (6*1)

</pre>
\param *vector   DOUBLE    (i)  stress vector
\param **matrix  DOUBLE    (o)  stress tensor

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_blowup_stress(DOUBLE  *vector,
                      DOUBLE **matrix);
		      
/*!----------------------------------------------------------------------
\brief calculates stress vector (6*1) out of stress tensor (3*3)

<pre>                                                              fh 03/03
This routine calculates the stress vector (6*1) out of stress tensor (3*3)

</pre>
\param *vector   DOUBLE    (o)  stress vector
\param **matrix  DOUBLE    (i)  stress tensor

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_compr_stress(DOUBLE  *vector,
                     DOUBLE **matrix);
		     
/*!----------------------------------------------------------------------
\brief calculates strain tensor (3*3) out of strain vector (6*1)

<pre>                                                              fh 03/03
This routine calculates the strain tensor (3*3) out of strain vector (6*1)

</pre>
\param *vector   DOUBLE    (i)  strain vector
\param **matrix  DOUBLE    (o)  strain tensor

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_blowup_strain(DOUBLE  *vector,
                      DOUBLE **matrix);
		      
/*!----------------------------------------------------------------------
\brief calculates strain vector (6*1) out of strain tensor (3*3)

<pre>                                                              fh 03/03
This routine calculates the strain vector (6*1) out of strain tensor (3*3)

</pre>
\param *vector   DOUBLE    (o)  stress vector
\param **matrix  DOUBLE    (i)  stress tensor

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_compr_strain(DOUBLE  *vector,
                     DOUBLE **matrix);
		     
/*!----------------------------------------------------------------------
\brief calculates Cep matrix (6*6) out of Cep tensor (3*3*3*3)

<pre>                                                              fh 03/03
This routine calculates the elastoplastic tangent matrix (6*6) out of 
elastoplastic tangent tensor (3*3*3*3)

</pre>
\param **matrix    DOUBLE    (o)  Cep matrix
\param ****tensor  DOUBLE    (i)  Cep tensor

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_compr_Cep(DOUBLE   **matrix,
                  DOUBLE ****tensor);
		  
/*!----------------------------------------------------------------------
\brief calculates Kronecker-delta for (3*3) tensor

<pre>                                                              fh 03/03
This routine calculates the Kronecker-delta for (3*3) tensor

</pre>
\param **delta    DOUBLE    (i/o)  Kronecker-delta


\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_kronecker(DOUBLE **delta);

/*!----------------------------------------------------------------------
\brief calculates volumetric part of elasticity tensor

<pre>                                                              fh 03/03
This routine calculates the volumetric part of the elasticity tensor

</pre>
\param bulk       DOUBLE     (i)   bulk modulus
\param **delta    DOUBLE     (i)   Kronecker-delta
\param ****c_vol  DOUBLE     (o)   volumetric part of elasticity tensor


\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_cal_volCel(DOUBLE     bulk,
                   DOUBLE   **delta,
	           DOUBLE ****c_vol);
		   
/*!----------------------------------------------------------------------
\brief calculates volumetric part of elasticity tensor

<pre>                                                              fh 03/03
This routine calculates the volumetric part of the elasticity tensor

</pre>
\param bulk       DOUBLE     (i)   bulk modulus
\param **delta    DOUBLE     (i)   Kronecker-delta
\param ****c_dev  DOUBLE     (o)   volumetric part of elasticity tensor


\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_cal_devCel(DOUBLE     shear,
                   DOUBLE   **delta,
	           DOUBLE ****c_dev);
		   
/*!----------------------------------------------------------------------
\brief calculates elasticity tensor

<pre>                                                              fh 03/03
This routine calculates the elasticity tensor

</pre>
\param ym         DOUBLE     (i)   youngs modulus
\param pv         DOUBLE     (i)   possions value
\param **delta    DOUBLE     (i)   Kronecker-delta
\param ****c_el   DOUBLE     (o)   elasticity tensor


\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_cal_Cel(DOUBLE     ym,
                DOUBLE     pv,
                DOUBLE   **delta,
                DOUBLE ****c_el);
		
/*!----------------------------------------------------------------------
\brief condensates stresses and elastoplastic tangent

<pre>                                                              fh 03/03
This routine condensates the stresses and the elastoplastic tangent

</pre>
\param **D        DOUBLE    (i/o)   elastoplastic tangent
\param *sig       DOUBLE    (i/o)   stress vector


\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: b3_mat_plast_mises()

*----------------------------------------------------------------------*/
void b3_condense(DOUBLE **D,
                 DOUBLE  *sig);		


/************************************************************************
 | b3_restart.c 	                                                |
 ************************************************************************/
/*!----------------------------------------------------------------------
\brief writes the data needed to restart the beam element

<pre>                                                              fh 04/03
This routine writes all the data needed to restart the beam element

</pre>
\param *actele         ELEMENT    (i/o) actual element
\param nhandle         INT         (i)  number of handles
\param *handles        LONGINT     (i)  handles
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: beam3()

*----------------------------------------------------------------------*/
void b3_write_restart(ELEMENT *actele, INT nhandle, long int *handles);

/*!----------------------------------------------------------------------
\brief reads the data needed to restart the beam element

<pre>                                                              fh 04/03
This routine reads all the data needed to restart the beam element

</pre>
\param *actele         ELEMENT    (i/o) actual element
\param nhandle         INT         (i)  number of handles
\param *handles        LONGINT     (i)  handles
               

\warning There is nothing special in this routine
\return void                                               
\sa calling:   ---;
    called by: beam3()

*----------------------------------------------------------------------*/
void b3_read_restart(ELEMENT *actele, INT nhandle, long int *handles);

#endif
/*! @} (documentation module close)*/		   
		   		  		     		      		     		      			
