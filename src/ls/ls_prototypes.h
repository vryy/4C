#ifdef D_LS
void ls2inp(
  ELEMENT*
  );


void ls2_funct(
  DOUBLE*,     
  DOUBLE**,    
  DOUBLE,
  DOUBLE,
  DIS_TYP
  );


void ls2(
  PARTITION*,
  INTRA*,
  ELEMENT*,             
  ARRAY*,   
  ARRAY*,   
  ARRAY*,
  ARRAY*,
  CALC_ACTION*
  );


void ls2_calset( 
  ELEMENT*,
  DOUBLE **,
  DOUBLE*,
  DOUBLE*
  );


void ls2_calset1( 
  ELEMENT*,     
  INT,
  DOUBLE*
  );


void ls2_jaco(
  DOUBLE**,
  DOUBLE*,
  DOUBLE**,
  DOUBLE**,
  DOUBLE*,
  INT,
  ELEMENT*
  );


void ls2_gder(
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  DOUBLE,
  INT
  );


void ls_isi();


void ls_init(
  FIELD*,
  LS_DYNAMIC*,
  INT
  );


void ls_initdirich(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_setdirich(
  FIELD*,
  LS_DYNAMIC*,
  INT
  );


void ls_algout(
  LS_DYNAMIC*
  );


INT ls_convcheck(
  LS_DYNAMIC*,
  DOUBLE,
  INT,
  DOUBLE,
  DOUBLE
  );


INT ls_steadycheck(
  LS_DYNAMIC*,
  FIELD*,
  INT
  );


void ls_norm(
  LS_DYNAMIC*,
  FIELD*,
  INT,
  DOUBLE*
  );


void fluid_to_ls(
  const FIELD*,
  const FIELD*
  );


void ls_update(
  FRONTLSFLAG
  );


void ls_updatelement(
  ELEMENT*,
  DOUBLE*
  );


void is_tricut(
  DOUBLE*,
  INT*,
  INT*,
  INT*,
  INT*
  );


void ls_compintersection(
  DOUBLE*,
  DOUBLE*,
  INT,
  INT
  );


void ls_comppoint(
  DOUBLE*,
  DOUBLE*,
  INT,
  INT,
  INT
  );


void ls_updateneighbor(
  ELEMENT*, 
  INT,
  INT
  );


void ls_makepatch(
  ELEMENT*, 
  ELEMENT *[], 
  INT*
  );


void ls_construct(
  ELEMENT*
  );


void ls_reset();


void ls_resetintdata(
  LS_INT_DATA*
  );


void ls_resetpolydata(
  LS_POLY_DATA*
  );


void ls_initialize();


void ls_updt();


void ls_write();


void ls_finalize();


void ls_printelinfo(
  ELEMENT*
  );


void ls_printnodeinfo(
  NODE*
  );


void ls_setdata();


void ls2_init();


void ls2_calc(
  ELEMENT*
  );


void ls2_calc_reinitialized(
  ELEMENT*
  );


void ls2_calc_nonlocalized(
  ELEMENT*
  );


void ls2_calc_localized(
  ELEMENT*
  );


void ls_printelements();				 


void ls_printelement(
  ELEMENT*
  );


void ls_init_single_circle(
  FIELD*,
  LS_DYNAMIC*									
  );


void ls_init_double_circle(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_triple_circle(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_single_circle_sharp(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_double_circle_sharp(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_triple_circle_sharp(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_single_line(
  FIELD*,
  LS_DYNAMIC*									
  );


void ls_init_double_line(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_triple_line(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_single_line_sharp(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_double_line_sharp(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_triple_line_sharp(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_to_matlab();									  


void ls2_setfluidvel_byuser(
  INT,
  INT
  );


void ls2_setfluidvel_byfluid(
  ELEMENT*
  );


void ls2_setvelocity(
  ELEMENT*
  );


void ls_dyn();


void ls_levelset(
  INT
  );


void ls_levelset_init();


void ls_levelset_solv();


void ls_levelset_fina();


void ls_levelset_clea();


void ls_fluid(
  INT
  );


void ls_fluid_init();


void ls_fluid_solv();


void ls_fluid_fina();


void ls_fluid_clea();


void ls_activate();


void ls_localize();


void ls_polygon_con(
  INT*,
  INT,
  INT
  );

void ls_polygonize();


void ls_updt_material();


void ls_init_material();


void ls_main(
  LSFLAG
  );


void ls_inp_gendata(
  LS_GEN_DATA*
  );


void ls_find_compatible_ele(
  const ELEMENT*,
  const ELEMENT*,
  INT*
  );


void ls2_intg(
  LS2_INTG_DATA*
  );


void ls_check_profile();
#endif
