#ifdef D_XFEM
void xfem_f2_init(
  ELEMENT*
  );


void xfem_f2_array_init();


void xfem_f2_loc_con();


void xfem_f2_iand();


void xfem_f2_loc_ass_intforce(
  DOUBLE*,
  DOUBLE*
  );


void xfem_f2_loc_ass_tangent();


void xfem_polygon(
  XFEMPOLYFLAG,
  ELEMENT*
  );


void xfem_polygon_init(
  ELEMENT*
  );


void xfem_polygon_cons(
  ELEMENT*
  );


void xfem_polygon_GP(
  ELEMENT*
  );


void xfem_polygon_getsubp(
  INT,
  INT
  );


void xfem_polygon_compGP(
  DIS_TYP
  );


void xfem_polygon_funct();


void xfem_polygon_deriv();


void xfem_polygon_resNewton();


void xfem_polygon_tanNewton();


void xfem_polygon_write(
  ELEMENT*
  );


void xfem_polygon_open();


void xfem_polygon_close();


void xfem_polygon_area_rect();


DOUBLE xfem_polygon_area_subtri(
  INT
  );


DOUBLE xfem_polygon_area_tri();


void xfem_polygon_target_subtri(
  INT
  );


void xfem_polygon_target_tri();


void xfem_polygon_mat(
  ELEMENT*
  );


void xfem_f2_funct(
  DOUBLE*,     
  DOUBLE**,    
  DOUBLE**,   
  DOUBLE,        
  DOUBLE,        
  DIS_TYP,      
  DOUBLE*,
  INT,
  INT
  );


void xfem_f2_funct1(
  DOUBLE*,     
  DOUBLE,        
  DOUBLE,        
  DIS_TYP,      
  INT,
  DOUBLE*,
  INT
  );


void xfem_f2_derxy(
  DOUBLE**,     
  DOUBLE**,    
  DOUBLE**,      
  DOUBLE,      
  INT,
  DOUBLE*,
  DOUBLE*,     
  INT
  );


void xfem_f2_veli(
  DOUBLE*,
  DOUBLE*,    
  DOUBLE**,     
  INT
  );


void xfem_f2_vder(
  DOUBLE**,     
  DOUBLE**,    
  DOUBLE**,     
  INT
  );


void xfem_f2_vder2(
  DOUBLE**,    
  DOUBLE**,    
  DOUBLE**,      
  INT
  );


void xfem_f2_calkvv(
  ELEMENT*,
  DOUBLE**,   
  DOUBLE*,
  DOUBLE*,
  DOUBLE**, 
  DOUBLE*,  
  DOUBLE**,  
  DOUBLE,    
  DOUBLE,   
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_calkvp(
  DOUBLE**,   
  DOUBLE*,  
  DOUBLE**,  
  DOUBLE,    
  INT,
  INT*
  );


void xfem_f2_calkvg( 
  DOUBLE**,   
  DOUBLE**, 
  DOUBLE*,  
  DOUBLE,    
  INT,
  INT*
  );


void xfem_f2_calmvv(
  DOUBLE**,  
  DOUBLE*, 
  DOUBLE,   
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_derxy2(
  DOUBLE**,     
  DOUBLE**,      
  DOUBLE**,     
  DOUBLE**,  
  DOUBLE**,  
  DOUBLE**, 
  DOUBLE**, 
  INT,
  DOUBLE*,
  DOUBLE*,
  INT
  );


void xfem_f2_calstabkvv(			      
  ELEMENT*,    
  DOUBLE**,  
  DOUBLE*,
  DOUBLE*, 
  DOUBLE*,
  DOUBLE**, 
  DOUBLE*,  
  DOUBLE**,  
  DOUBLE**, 
  DOUBLE,    
  DOUBLE,   
  INT,    
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_calstabkvp(
  ELEMENT*,    
  DOUBLE**, 
  DOUBLE*,
  DOUBLE*, 
  DOUBLE**, 
  DOUBLE**,
  DOUBLE,   
  DOUBLE,  
  INT,   
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_calstabmvv(
  ELEMENT*,     
  DOUBLE**,  
  DOUBLE*, 
  DOUBLE*,  
  DOUBLE**,  
  DOUBLE**, 
  DOUBLE,    
  DOUBLE,   
  INT,    
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_calstabkpv(
  ELEMENT*,
  DOUBLE**,   
  DOUBLE*,
  DOUBLE*, 
  DOUBLE**, 
  DOUBLE*,  
  DOUBLE**,  
  DOUBLE**, 
  DOUBLE,    
  DOUBLE,   
  INT,    
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_calstabmpv(
  DOUBLE**,   
  DOUBLE*,  
  DOUBLE**,  
  DOUBLE,    
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_calstabkpp(
  DOUBLE**,   
  DOUBLE**,  
  DOUBLE,    
  INT
  );


void xfem_f2_calint(
  FLUID_DATA*,     
  ELEMENT*,     
  INT*,
  DOUBLE**,   
  DOUBLE**,   
  DOUBLE*, 
  DOUBLE*, 
  DOUBLE**,
  DOUBLE*,   
  DOUBLE**,   
  DOUBLE**,  
  DOUBLE**,     
  DOUBLE**,   
  DOUBLE**,  
  DOUBLE**,   
  DOUBLE**,  
  DOUBLE*,   
  DOUBLE*,
  DOUBLE*,	        	       
  DOUBLE*,  
  DOUBLE*, 
  DOUBLE*,  
  DOUBLE**,  
  DOUBLE*,  
  DOUBLE**, 
  DOUBLE**,     
  DOUBLE**
  );


void xfem_f2_calset( 
  ELEMENT*,     
  DOUBLE**,
  DOUBLE**,    
  DOUBLE**,   
  DOUBLE*,
  DOUBLE*,
  DOUBLE*,
  INT*
  );


void xfem_f2_calelesize(
  ELEMENT*,    
  FLUID_DATA*, 
  DOUBLE**,
  DOUBLE*,  
  DOUBLE**,  
  DOUBLE**,  		 
  DOUBLE**,    
  DOUBLE**, 
  DOUBLE**,
  DOUBLE**,    		  
  DOUBLE*, 
  DOUBLE**
  );


void xfem_f2_calgaltfv(
  DOUBLE*,    
  DOUBLE*,    
  DOUBLE*,   
  DOUBLE*,    
  DOUBLE**,    
  DOUBLE**,   
  DOUBLE,   
  DOUBLE,     
  DOUBLE,      
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_calstabtfv(
  ELEMENT*,      
  DOUBLE*,  
  DOUBLE*,  
  DOUBLE*, 
  DOUBLE*,  
  DOUBLE**,   
  DOUBLE**,  
  DOUBLE**,  
  DOUBLE**, 
  DOUBLE*,  
  DOUBLE,     
  DOUBLE,    
  INT,   
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_calstabtfp(
  DOUBLE*,    
  DOUBLE**,   
  DOUBLE**, 
  DOUBLE*,  
  DOUBLE*,  
  DOUBLE*,  
  DOUBLE,    
  DOUBLE,     
  INT,   
  INT,
  DOUBLE
  );


void xfem_f2_calele(
  FLUID_DATA*, 
  ELEMENT*,             
  ARRAY*,   
  ARRAY*,   
  ARRAY*,       
  ARRAY*, 
  ARRAY*,		
  INT*,      
  INT*,
  INT,
  INT,
  INT            
  );


void xfem_fluid2(
  PARTITION*,
  INTRA*,
  ELEMENT*,             
  ARRAY*,   
  ARRAY*,   
  ARRAY*, 
  ARRAY*, 
  ARRAY*, 
  CALC_ACTION*,
  INT*,
  INT*,
  CONTAINER*
  );


void xfem_f2_calgalexfv(
  DOUBLE*,     
  DOUBLE*,       
  DOUBLE*,
  DOUBLE*,
  DOUBLE,      
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_calstabexfv(
  ELEMENT*,  
  DOUBLE*,     
  DOUBLE**,
  DOUBLE**,      
  DOUBLE*,
  DOUBLE*,  
  DOUBLE,      
  DOUBLE,
  INT,
  INT,
  INT,
  INT*,
  DOUBLE
  ); 


void xfem_f2_calstabexfp(
  DOUBLE*,     
  DOUBLE**,       
  DOUBLE*,  
  DOUBLE,      
  INT,
  INT,
  DOUBLE
  );


void xfem_f2_intg(
  FLUID_DATA*
  );
#endif
