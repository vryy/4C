/*!---------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | PROBLEM TYPES                                          m.gee 7/01    |
 *----------------------------------------------------------------------*/
typedef enum _PROBLEM_TYP
{
                       prb_none,         /*  not a problem at all */
                       prb_fsi,          /*  fluid structure interaction problem */
                       prb_fsi_xfem,     /*  fluid structure interaction problem including XFEM interfaces*/
                       prb_fsi_lung,     /*  airway fsi problem with attached parenchyma balloon */
                       prb_gas_fsi,      /*  fsi with gas transport */
                       prb_tfsi_aero,    /*  thermo fluid structure interaction problem including the aero code */
                       prb_structure,    /*  structural problem */
                       prb_struct_ale,   /*  structural problem, ale formulation */
                       prb_fluid,        /*  fluid problem */
                       prb_fluid_xfem,   /*  fluid problem including XFEM interfaces */
                       prb_fluid_xfem2,  /*  fluid problem including XFEM interfaces */
                       prb_fluid_fluid_ale,/*fluid_fluid_ale problem */
                       prb_fluid_fluid,  /*  fluid_fluid problem */
                       prb_fluid_fluid_fsi,/*fluid_fluid fsi problem */
                       prb_fluid_ale,    /*  fluid on an ale mesh (no structure) */
                       prb_freesurf,     /*  free surface fluid */
                       prb_opt,          /*  structural optimization  problem */
                       prb_ale,          /*  pure ale problem */
                       prb_tsi,          /*  thermal structure interaction */
                       prb_thermo,       /*  thermal problem */
                       prb_fluid_pm,     /*  fluid with (any) projection method */
                       prb_scatra,       /*  scalar transport problem (e.g. convection-diffusion) */
                       prb_pfsi,         /*  projection fsi */
                       prb_loma,         /*  low-Mach-number flow problem */
                       prb_elch,         /*  electrochemical problem */
                       prb_combust,      /*  combustion problem */
                       prb_art_net,      /*  arterial network problem */ /*_1D_ARTERY_*/
                       prb_red_airways,  /*  reduced dimensional airways */
                       prb_biofilm_fsi   /*  biofilm growth problem */
} PROBLEM_TYP;
/* Mapping from problem type numbers to printable names. To be used to
 * initialize static variables. Keep in sync!
 * The trailing NULL is essential for the filters to read the problem
 * type! */
#define PROBLEMNAMES {"none","fsi","fsi_xfem","fsi_lung","gas_fsi","aero_tfsi","structure","structure_ale","fluid","fluid_xfem","fluid_xfem2","fluid_fluid_ale","fluid_fluid","fluid_fluid_fsi","fluid_ale","freesurf","opt","ale","tsi","thermo","fluid_pm","scatra","pfsi","loma","elch","combustion","art_net","red_airways","biofilm_fsi",NULL }
/*----------------------------------------------------------------------*
 | TIME TYPES                                             m.gee 7/01    |
 *----------------------------------------------------------------------*/
typedef enum _TIME_TYP
{
                       time_static,  /* time independent static analysis */
                       time_dynamic  /* time dependent analysis */
} TIME_TYP;


/*----------------------------------------------------------------------*
 | FIELD TYPES                                            m.gee 7/01    |
 *----------------------------------------------------------------------*/
typedef enum _FIELDTYP
{
                       none,        /* unknown type of mechanical field */
                       fluid,       /* fluid field */
                       xfluid,       /* xfluid field */
                       ale,         /* pseudo structural field */
                       structure,   /* structural field */
                       thermal,     /* thermal field */
                       pressure,    /* pure pressure field */
                       boundary,    /* boundary field */
                       scatra,      /* scalar transport field */
                       scatra1,      /* scalar transport field in case of multiple fields */
                       scatra2,      /* scalar transport field in case of multiple fields */
                       artery,      /* artery field*/
                       thermo,      /* thermal field */
                       fluidfluidboundary,  /*fluidfluidboundary field*/
                       red_airway  /* reduced dimensional airways */
} FIELDTYP;
/* Mapping from fieldtyp numbers to printable names. To be used to
 * initialize static variables. Keep in sync! */
#define FIELDNAMES {"none", "fluid", "xfluid", "ale", "structure", "thermal", "pressure", "boundary", "scatra", "scatra1", "scatra2", "artery", "thermo", "FluidFluidboundary", "red_airway", "inflow", NULL}


/*----------------------------------------------------------------------*
 | enum MATERIAL_TYP                                      m.gee 7/01    |
 | material laws                                                        |
 *----------------------------------------------------------------------*/
#if defined(CCADISCRET) && !defined(D_SHELL8)
#else
typedef enum _MATERIAL_TYP
{
                       m_stvenant,    /* St.Venant Kirchhoff material */
                       m_thermostvenant,    /* St.Venant Kirchhoff material with temperature */
                       m_pl_mises_3D, /* Stefans Mises*/
                       m_pl_mises,    /* von Mises material */
                       m_pl_hoff,     /* anisotropic plastic material based on hoffman criterion */
                       m_damage,      /* 3D damage matieral */
                       m_pl_foam,     /* foam material - large strains */
                       m_pl_mises_ls, /* von Mises material - large strains*/
                       m_pl_dp,       /* Drucker Prager material */
                       m_pl_epc,      /* elastoplastic concrete material */
                       m_pl_epc3D,    /* elastoplastic concrete material 3D formulation */
                       m_stvenpor,    /* porous St.Venant Kirchhoff material */
                       m_pl_por_mises,/* porous von Mises material */
                       m_neohooke,    /* Neo-Hooke material */
                       m_aaaneohooke, /* quasi Neo-Hooke material for aneurysmatic artery wall */
                       m_aaaraghavanvorp_damage, /* quasi Neo-Hooke material for aneurysmatic artery wall with damage*/
                       m_compogden,   /* compressible Ogden material (with shell8) */
                       m_viscohyper,  /* compressible viscous Ogden material (with shell8) */
                       m_fluid,       /* fluid */
                       m_sutherland_fluid,  /* fluid material with temperature dependence according to Sutherland law */
                       m_carreauyasuda,/* fluid with nonlinear viscosity according to Carreau-Yasuda*/
                       m_modpowerlaw,  /* fluid with nonlinear viscosity according to a modified power law*/
                       m_condif,      /* convection-diffusion */
                       m_sutherland_condif,  /* convection-diffusion material with temperature dependence according to Sutherland law */
                       m_pl_hash,     /* elpl. hashin delamination material */
                       m_el_orth,     /* elastic orthotropic material */
                       m_mfoc,        /* open cell metal foam */
                       m_mfcc,        /* closed cell metal foam */
                       m_nhmfcc,      /* foam, closed cell, based on modified Neo Hook */
                       m_multi_layer, /* multilayer material -> shell9*/
                       m_ifmat,        /* interface surface elasto-damage-plasto material*/
                       m_interf_therm, /* themodyn. based interface elasto-damage surface material*/
                       m_dam_mp,       /* isotropic damage model -> mazars/pijadier-cabot*/
                       m_damage_ge,    /* isotropic gradient enhanced damage model */
                       m_lung_penalty, /* lung tissue material with penalty function for incompressibility constraint*/
                       m_lung_ogden,   /* lung tissue material with compressible Ogden for volumetric part */
                       m_itskov,       /* hyperelastic polyconvex energy strain function following Itskov */
                       m_anisotropic_balzani,  /* anisotropic polyconvex material*/
                       m_mooneyrivlin,  /* Mooney-Rivlin material*/
                       m_yeoh,          /* Yeoh material*/
                       m_visconeohooke, /* Viscous NeoHookean Material */
                       m_viscoanisotropic, /* Viscous Anisotropic Fiber Material */
                       m_contchainnetw, /* Continuum Chain Network Material Law with remodeling */
                       m_artwallremod,  /* Arterial Wall Material Law (Holzapfel) with remodeling (Hariton) */
                       m_th_fourier_iso,  /* isotropic (linear) Fourier's law of heat conduction */
                       m_th_fourier_gen,  /* general (linear) Fourier's law of heat conduction */
                       m_vp_robinson,   /* Robinson's visco-plastic material */
                       m_struct_multiscale, /*  structural microscale approach */
                       m_matlist,       /* collection of single materials (used for scalar transport problems)*/
                       m_biocell,       /* biological cell model */
                       m_ion,           /* properties of an ion species in an electrolyte solution */
                       m_cnst_art,      /* 1D_Artery with constant material and geometric properties */
                       m_holzapfelcardiovascular, /* anisotropic fiber material for arteries */
                       m_humphreycardiovascular /* anisotropic material for arteries cf Humphrey */
} MATERIAL_TYP;
#endif

/*----------------------------------------------------------------------*/
/* The known time marching schemes for fluids. */
/*----------------------------------------------------------------------*/
typedef enum _FLUID_DYNTYPE
{
  dyntyp_nln_time_int=0,
  dyntyp_lin_time_int=1,
  dyntyp_pm_discont=2,
  dyntyp_pm_cont=3,
  dyntyp_pm_cont_laplace=4
} FLUID_DYNTYPE;

/*----------------------------------------------------------------------*/
/* The known solving strategies for fluids. */
/*----------------------------------------------------------------------*/
typedef enum _FLUID_SOLVINGSTRATEGIES
{
	fluid_solver_implicit=0,
	fluid_solver_pressurecorrection=1,
	fluid_solver_pressurecorrection_semiimplicit=2,
    fluid_solver_fluid_xfluid
} FLUID_SOLVINGSTRATEGIES;

/*----------------------------------------------------------------------*/
/* The coupling methods for FSI. */
/*----------------------------------------------------------------------*/
typedef enum _FSI_COUPLING
{
  fsi_coupling_freesurface=-1,
  fsi_coupling_undefined=0,
  fsi_basic_sequ_stagg=1,
  fsi_sequ_stagg_pred=2,
  fsi_sequ_stagg_shift=3,
  fsi_iter_stagg_fixed_rel_param=4,
  fsi_iter_stagg_AITKEN_rel_param=5,
  fsi_iter_stagg_steep_desc=6,
  fsi_iter_stagg_CHEB_rel_param=7,
  fsi_iter_stagg_AITKEN_rel_force=8,
  fsi_iter_stagg_steep_desc_force=9,
  fsi_iter_stagg_Newton_FD=10,
  fsi_iter_stagg_Newton_I=11,
  fsi_iter_nox=12,
  fsi_iter_monolithicfluidsplit=13,
  fsi_iter_monolithiclagrange,
  fsi_iter_monolithicstructuresplit,
  fsi_iter_lung_monolithicstructuresplit,
  fsi_iter_lung_monolithicfluidsplit,
  fsi_iter_mortar_monolithicstructuresplit,
  fsi_iter_mortar_monolithicfluidsplit,
  fsi_iter_constr_monolithicstructuresplit,
  fsi_iter_constr_monolithicfluidsplit,
  fsi_iter_monolithicxfem,
  fsi_iter_stagg_NLCG,
  fsi_iter_stagg_MFNK_FD,
  fsi_iter_stagg_MFNK_FSI,
  fsi_iter_stagg_MPE,
  fsi_iter_stagg_RRE,
  fsi_pseudo_structureale,
  fsi_iter_fluidfluid_monolithicstructuresplit,
} FSI_COUPLING;


/*!----------------------------------------------------------------------
\brief enum of arterial network dynamic types

<pre>                                                      ismail 11/08
This is the enumeration of all types of different integration schemes
</pre>

*-----------------------------------------------------------------------*/
typedef enum _ARTNET_DYNTYPE
{
  typ_tay_gal=0
} _ARTNET_DYNTYPE;


/*!----------------------------------------------------------------------
\brief enum of reduced dimensional airways dynamic types

<pre>                                                     ismail 01/10
This is the enumeration of all types of different integration schemes
</pre>

*-----------------------------------------------------------------------*/
typedef enum _RED_AIRWAYS_DYNTYPE
{
  typ_crank_nicolson,
  linear,
  nonlinear
} _RED_AIRWAYS_DYNTYPE;
