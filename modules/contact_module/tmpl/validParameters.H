
// Teuchos (Trilinos)
#include <Teuchos_Array.hpp>
#include <Teuchos_StrUtils.hpp>
#include <Teuchos_any.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

//-----------------------------------------------------------------------------
// function for definition of integer parameters
void IntParameter(std::string const &paramName,
                  int const value,
                  std::string const &docString,
                  Teuchos::ParameterList *paramList)
{
  Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes validator(false);
  validator.allowInt(true);
  Teuchos::setIntParameter(paramName,value,
                           docString,
                           paramList,validator);

  return;
}
//-----------------------------------------------------------------------------
// function for definition of double parameters
void DoubleParameter(std::string const &paramName,
                     double const &value,
                     std::string const &docString,
                     Teuchos::ParameterList *paramList)
{
  Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes validator(false);
  validator.allowDouble(true);
  validator.allowInt(true);
  Teuchos::setDoubleParameter(paramName,value,
                              docString,
                              paramList,validator);

  return;
}


//-----------------------------------------------------------------------------
// function for getting list with valid parameters
void getValidParameters(Teuchos::ParameterList& list)
{



  Teuchos::setStringToIntegralParameter<int>("LM_SHAPEFCN","Dual","Type of employed set of shape functions",
         Teuchos::tuple<std::string>("Dual", "dual", "Standard", "standard", "std", "PetrovGalerkin", "petrovgalerkin", "pg" ),
         Teuchos::tuple<int>(INPAR::MORTAR::shape_dual,INPAR::MORTAR::shape_dual,
        INPAR::MORTAR::shape_standard,INPAR::MORTAR::shape_standard,INPAR::MORTAR::shape_standard,
        INPAR::MORTAR::shape_petrovgalerkin,INPAR::MORTAR::shape_petrovgalerkin,INPAR::MORTAR::shape_petrovgalerkin),&list);

  Teuchos::setStringToIntegralParameter<int>("SEARCH_ALGORITHM","Binarytree","Type of contact search",
         Teuchos::tuple<std::string>("BruteForce", "bruteforce", "BruteForceEleBased", "bruteforceelebased", "BinaryTree", "Binarytree", "binarytree"),
         Teuchos::tuple<int>(INPAR::MORTAR::search_bfele,INPAR::MORTAR::search_bfele,
              INPAR::MORTAR::search_bfele,INPAR::MORTAR::search_bfele,
              INPAR::MORTAR::search_binarytree,INPAR::MORTAR::search_binarytree,
              INPAR::MORTAR::search_binarytree),&list);

  DoubleParameter("SEARCH_PARAM",0.3,"Radius / Bounding volume inflation for contact search",&list);

  Teuchos::setStringToIntegralParameter<int>("SEARCH_USE_AUX_POS","Yes","If chosen auxiliary position is used for computing dops",
         Teuchos::tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO"),
         Teuchos::tuple<int>(true,false,true,false,true,false),&list);

  Teuchos::setStringToIntegralParameter<int>("LM_QUAD","undefined","Type of LM interpolation for quadratic FE",
         Teuchos::tuple<std::string>("undefined", "quad", "quadratic", "pwlin", "piecewiselinear", "lin","linear"),
         Teuchos::tuple<int>(INPAR::MORTAR::lagmult_undefined, INPAR::MORTAR::lagmult_quad, INPAR::MORTAR::lagmult_quad, INPAR::MORTAR::lagmult_pwlin, INPAR::MORTAR::lagmult_pwlin, INPAR::MORTAR::lagmult_lin, INPAR::MORTAR::lagmult_lin),&list);

  Teuchos::setStringToIntegralParameter<int>("CROSSPOINTS","No","If chosen, multipliers are removed from crosspoints / edge nodes",
         Teuchos::tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO"),
         Teuchos::tuple<int>(true,false,true,false,true,false),&list);

  Teuchos::setStringToIntegralParameter<int>("LM_DUAL_CONSISTENT","Yes","If chosen consistent dual shape functions are calculated",
         Teuchos::tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO"),
         Teuchos::tuple<int>(true,false,true,false,true,false),&list);

  Teuchos::setStringToIntegralParameter<int>("LM_NODAL_SCALE","No","If chosen a nodal scaling factor is calculated for each LM",
         Teuchos::tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO"),
         Teuchos::tuple<int>(true,false,true,false,true,false),&list);

  Teuchos::setStringToIntegralParameter<int>("HERMITE_SMOOTHING","No","If chosen hermite interface smoothing is activated for line2 elements",
         Teuchos::tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO"),
         Teuchos::tuple<int>(true,false,true,false,true,false),&list);

  Teuchos::setStringToIntegralParameter<int>("MESH_RELOCATION","Initial","Type of mesh relocation",
         Teuchos::tuple<std::string>("Initial","initial", "Every_Timestep", "every_timestep", "No", "no"),
         Teuchos::tuple<int>(INPAR::MORTAR::relocation_initial,  INPAR::MORTAR::relocation_initial,
              INPAR::MORTAR::relocation_timestep, INPAR::MORTAR::relocation_timestep,
              INPAR::MORTAR::relocation_none,     INPAR::MORTAR::relocation_none),&list);

  Teuchos::setStringToIntegralParameter<int>("REDUNDANT_STORAGE","Master","Type of redundancy in interface storage",
         Teuchos::tuple<std::string>("All","all", "Master", "master", "None", "none"),
         Teuchos::tuple<int>(INPAR::MORTAR::redundant_all,    INPAR::MORTAR::redundant_all,
              INPAR::MORTAR::redundant_master, INPAR::MORTAR::redundant_master,
              INPAR::MORTAR::redundant_none,   INPAR::MORTAR::redundant_none),&list);

  Teuchos::setStringToIntegralParameter<int>("PARALLEL_STRATEGY","redundant_ghosting","Type of parallel interface evaluation",
         Teuchos::tuple<std::string>("rg","redundant_ghosting", "ghosting", "rre", "roundrobinevaluate","RoundRobinEvaluate", "rrg", "roundrobinghost", "RoundRobinGhost", "bs", "binningstrategy", "binning"),
         Teuchos::tuple<int>(INPAR::MORTAR::ghosting_redundant, INPAR::MORTAR::ghosting_redundant,INPAR::MORTAR::ghosting_redundant,
              INPAR::MORTAR::roundrobinevaluate, INPAR::MORTAR::roundrobinevaluate,INPAR::MORTAR::roundrobinevaluate,
              INPAR::MORTAR::roundrobinghost, INPAR::MORTAR::roundrobinghost,INPAR::MORTAR::roundrobinghost,
              INPAR::MORTAR::binningstrategy, INPAR::MORTAR::binningstrategy, INPAR::MORTAR::binningstrategy),&list);

  Teuchos::setStringToIntegralParameter<int>("PARALLEL_REDIST","Static","Type of redistribution algorithm",
         Teuchos::tuple<std::string>("None","none", "No", "no", "Static", "static", "Dynamic", "dynamic"),
         Teuchos::tuple<int>(INPAR::MORTAR::parredist_none, INPAR::MORTAR::parredist_none,
              INPAR::MORTAR::parredist_none, INPAR::MORTAR::parredist_none,
              INPAR::MORTAR::parredist_static, INPAR::MORTAR::parredist_static,
              INPAR::MORTAR::parredist_dynamic, INPAR::MORTAR::parredist_dynamic),&list);

  Teuchos::setStringToIntegralParameter<int>("ALGORITHM","Mortar","Type of meshtying/contact algorithm",
         Teuchos::tuple<std::string>("mortar","Mortar", "nts","NTS"),
         Teuchos::tuple<int>(INPAR::MORTAR::algorithm_mortar,INPAR::MORTAR::algorithm_mortar, INPAR::MORTAR::algorithm_nts,INPAR::MORTAR::algorithm_nts),&list);

  DoubleParameter("MAX_BALANCE",2.0,"Maximum value of load balance measure before parallel redistribution",&list);

  IntParameter("MIN_ELEPROC",0,"Minimum no. of elements per processor for parallel redistribution",&list);

  Teuchos::setStringToIntegralParameter<int>("INTTYPE","Segments","Type of numerical integration scheme",
         Teuchos::tuple<std::string>("Segments","segments", "Elements","elements", "Elements_BS", "elements_BS"),
         Teuchos::tuple<int>(INPAR::MORTAR::inttype_segments, INPAR::MORTAR::inttype_segments,
        INPAR::MORTAR::inttype_elements, INPAR::MORTAR::inttype_elements,
        INPAR::MORTAR::inttype_elements_BS, INPAR::MORTAR::inttype_elements_BS),&list);

  IntParameter("NUMGP_PER_DIM",0,"Number of employed integration points per dimension",&list);

  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for meshtying and contact. This parameter is not relevant for the standalone contact module.",&list);

  Teuchos::setStringToIntegralParameter<int>("RESTART_WITH_CONTACT","No","Must be chosen if a non-contact simulation is to be restarted with contact",
         Teuchos::tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO"),
         Teuchos::tuple<int>(true,false,true,false,true,false),&list);

  Teuchos::setStringToIntegralParameter<int>("ADHESION","None","Type of adhesion law",
         Teuchos::tuple<std::string>("None","none","bounded","b"),
         Teuchos::tuple<int>(INPAR::CONTACT::adhesion_none,INPAR::CONTACT::adhesion_none, INPAR::CONTACT::adhesion_bound,INPAR::CONTACT::adhesion_bound),&list);

  Teuchos::setStringToIntegralParameter<int>("STRATEGY","LagrangianMultipliers","Type of employed solving strategy",
         Teuchos::tuple<std::string>("LagrangianMultipliers","lagrange", "Lagrange",
                           "PenaltyMethod","penalty", "Penalty",
                           "UzawaAugementedLagrange","uzawa","Uzawa",
                           "AugmentedLagrange","augmented", "Augmented"),
         Teuchos::tuple<int>(INPAR::CONTACT::solution_lagmult, INPAR::CONTACT::solution_lagmult, INPAR::CONTACT::solution_lagmult,
                INPAR::CONTACT::solution_penalty, INPAR::CONTACT::solution_penalty, INPAR::CONTACT::solution_penalty,
                INPAR::CONTACT::solution_uzawa, INPAR::CONTACT::solution_uzawa, INPAR::CONTACT::solution_uzawa,
                INPAR::CONTACT::solution_augmented, INPAR::CONTACT::solution_augmented, INPAR::CONTACT::solution_augmented),&list);

  Teuchos::setStringToIntegralParameter<int>("SYSTEM","Condensed","Type of linear system setup / solution",
         Teuchos::tuple<std::string>("Condensed","condensed", "cond", "SaddlePoint","Saddlepoint","saddlepoint", "sp"),
         Teuchos::tuple<int>(INPAR::CONTACT::system_condensed, INPAR::CONTACT::system_condensed, INPAR::CONTACT::system_condensed,
                INPAR::CONTACT::system_saddlepoint, INPAR::CONTACT::system_saddlepoint,
                INPAR::CONTACT::system_saddlepoint, INPAR::CONTACT::system_saddlepoint),&list);

  DoubleParameter("PENALTYPARAM",0.0,"Penalty parameter for penalty / Uzawa augmented solution strategy",&list);

  DoubleParameter("PENALTYPARAMTAN",0.0,"Tangential penalty parameter for penalty / Uzawa augmented solution strategy",&list);

  IntParameter("UZAWAMAXSTEPS",10,"Maximum no. of Uzawa steps for Uzawa solution strategy",&list);

  DoubleParameter("UZAWACONSTRTOL",1.0e-8,"Tolerance of constraint norm for Uzawa solution strategy",&list);

  Teuchos::setStringToIntegralParameter<int>("SEMI_SMOOTH_NEWTON","Yes","If chosen semi-smooth Newton concept is applied",
         Teuchos::tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO"),
         Teuchos::tuple<int>(true,false,true,false,true,false),&list);

  DoubleParameter("SEMI_SMOOTH_CN",1.0,"Weighting factor cn for semi-smooth PDASS",&list);

  DoubleParameter("SEMI_SMOOTH_CT",1.0,"Weighting factor ct for semi-smooth PDASS",&list);

  Teuchos::setStringToIntegralParameter<int>("CONTACTFORCE_ENDTIME","No","If chosen, the contact force is not evaluated at the generalized midpoint, but at the end of the time step",
         Teuchos::tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO"),
         Teuchos::tuple<int>(true,false,true,false,true,false),&list);

  Teuchos::setStringToIntegralParameter<int>("VELOCITY_UPDATE","No","If chosen, velocity update method is applied",
         Teuchos::tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO"),
         Teuchos::tuple<int>(true,false,true,false,true,false),&list);

  Teuchos::setStringToIntegralParameter<int>("EMOUTPUT","None","Type of energy and momentum output",
         Teuchos::tuple<std::string>("None","none", "No", "no", "Screen", "screen", "File", "file", "Both", "both"),
         Teuchos::tuple<int>(INPAR::CONTACT::output_none, INPAR::CONTACT::output_none,
              INPAR::CONTACT::output_none, INPAR::CONTACT::output_none,
              INPAR::CONTACT::output_screen, INPAR::CONTACT::output_screen,
              INPAR::CONTACT::output_file, INPAR::CONTACT::output_file,
              INPAR::CONTACT::output_both, INPAR::CONTACT::output_both),&list);

  Teuchos::setStringToIntegralParameter<int>("ERROR_NORMS","None","Choice of analytical solution for error norm computation",
         Teuchos::tuple<std::string>("None","none", "No", "no", "Zero", "zero", "Bending", "bending", "Sphere", "sphere", "Thick", "thick"),
         Teuchos::tuple<int>(INPAR::CONTACT::errornorms_none, INPAR::CONTACT::errornorms_none,
              INPAR::CONTACT::errornorms_none, INPAR::CONTACT::errornorms_none,
              INPAR::CONTACT::errornorms_zero, INPAR::CONTACT::errornorms_zero,
              INPAR::CONTACT::errornorms_bending, INPAR::CONTACT::errornorms_bending,
              INPAR::CONTACT::errornorms_sphere, INPAR::CONTACT::errornorms_sphere,
              INPAR::CONTACT::errornorms_thicksphere, INPAR::CONTACT::errornorms_thicksphere),&list);

  Teuchos::setStringToIntegralParameter<int>("INITCONTACTBYGAP","No","Initialize init contact by weighted gap vector",
         Teuchos::tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO"),
         Teuchos::tuple<int>(true,false,true,false,true,false),&list);

  DoubleParameter("INITCONTACTGAPVALUE",0.0,"Value for initialization of init contact set with gap vector",&list);

  Teuchos::setStringToIntegralParameter<int>("MESH_ADAPTIVE_CN","No","use a scaling of cn with the local mesh size",
         Teuchos::tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO"),
         Teuchos::tuple<int>(true,false,true,false,true,false),&list);

  Teuchos::setStringToIntegralParameter<int>("MESH_ADAPTIVE_CT","No","use a scaling of ct with the local mesh size",
         Teuchos::tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO"),
         Teuchos::tuple<int>(true,false,true,false,true,false),&list);

  Teuchos::setStringToIntegralParameter<int>("CONSTRAINT_DIRECTIONS","ntt","formulation of constraints in normal/tangential or xyz-direction",
         Teuchos::tuple<std::string>("ntt", "xyz"),
         Teuchos::tuple<int>(INPAR::CONTACT::constr_ntt, INPAR::CONTACT::constr_xyz),&list);

  Teuchos::setStringToIntegralParameter<int>("WEARLAW","None","Type of wear law",
         Teuchos::tuple<std::string>("None","none", "Archard","archard"),
         Teuchos::tuple<int>(INPAR::WEAR::wear_none, INPAR::WEAR::wear_none, INPAR::WEAR::wear_archard, INPAR::WEAR::wear_archard),&list);

  Teuchos::setStringToIntegralParameter<int>("WEARCOEFF_CONF","material","configuration in which wcoeff is defined",
         Teuchos::tuple<std::string>("material","mat", "spatial","sp"),
         Teuchos::tuple<int>(INPAR::WEAR::wear_conf_mat, INPAR::WEAR::wear_conf_mat, INPAR::WEAR::wear_conf_sp, INPAR::WEAR::wear_conf_sp),&list);

  Teuchos::setStringToIntegralParameter<int>("WEAR_SHAPEFCN","std","Type of employed set of shape functions for wear",
         Teuchos::tuple<std::string>("Dual", "dual", "Standard", "standard", "std"),
         Teuchos::tuple<int>(INPAR::WEAR::wear_shape_dual, INPAR::WEAR::wear_shape_dual, INPAR::WEAR::wear_shape_standard, INPAR::WEAR::wear_shape_standard, INPAR::WEAR::wear_shape_standard),&list);

  DoubleParameter("WEARCOEFF",0.0,"Wear coefficient for slave surface",&list);

  DoubleParameter("WEARCOEFF_MASTER",0.0,"Wear coefficient for master surface",&list);

  DoubleParameter("WEAR_TIMERATIO",1.0,"Time step ratio between wear and spatial time scale",&list);

  DoubleParameter("SSSLIP",1.0,"Fixed slip for steady state wear",&list);

  Teuchos::setStringToIntegralParameter<int>("SSWEAR","No","flag for steady state wear",
         Teuchos::tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO"),
         Teuchos::tuple<int>(true,false,true,false,true,false),&list);

  Teuchos::setStringToIntegralParameter<int>("VOLMASS_OUTPUT","No","flag for output of mass/volume in ref,mat and cur. conf.",
         Teuchos::tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO"),
         Teuchos::tuple<int>(true,false,true,false,true,false),&list);

  Teuchos::setStringToIntegralParameter<int>("BOTH_SIDED_WEAR","No","Definition of wear side",
         Teuchos::tuple<std::string>("No","no", "none" , "Mapping","mapping", "map", "bothdiscr","bd", "sm"),
         Teuchos::tuple<int>(INPAR::WEAR::wear_slave, INPAR::WEAR::wear_slave, INPAR::WEAR::wear_slave,
                INPAR::WEAR::wear_both_map, INPAR::WEAR::wear_both_map,INPAR::WEAR::wear_both_map,
                INPAR::WEAR::wear_both_discr, INPAR::WEAR::wear_both_discr,INPAR::WEAR::wear_both_discr),&list);

  Teuchos::setStringToIntegralParameter<int>("WEARTYPE","internal_state_expl","Definition of wear algorithm",
         Teuchos::tuple<std::string>("intstate_impl","is_impl", "internal_state_impl", "intstate_expl","is_expl", "internal_state_expl", "primvar","pv", "primary_variable"),
         Teuchos::tuple<int>(INPAR::WEAR::wear_intstate_impl, INPAR::WEAR::wear_intstate_impl, INPAR::WEAR::wear_intstate_impl,
                INPAR::WEAR::wear_intstate_expl, INPAR::WEAR::wear_intstate_expl, INPAR::WEAR::wear_intstate_expl,
                INPAR::WEAR::wear_primvar, INPAR::WEAR::wear_primvar, INPAR::WEAR::wear_primvar),&list);

  Teuchos::setStringToIntegralParameter<int>("WEAR_COUPALGO","stagg","Definition of wear (ALE) coupling algorithm",
         Teuchos::tuple<std::string>("stagg","s", "iterstagg","is", "monolithic","mono"),
         Teuchos::tuple<int>(INPAR::WEAR::wear_stagg, INPAR::WEAR::wear_stagg, INPAR::WEAR::wear_iterstagg, INPAR::WEAR::wear_iterstagg, INPAR::WEAR::wear_monolithic, INPAR::WEAR::wear_monolithic),&list);

  Teuchos::setStringToIntegralParameter<int>("WEAR_TIMESCALE","equal","Definition wear time scale compares to std. time scale",
         Teuchos::tuple<std::string>("equal","e", "different","d"),
         Teuchos::tuple<int>(INPAR::WEAR::wear_time_equal, INPAR::WEAR::wear_time_equal, INPAR::WEAR::wear_time_different, INPAR::WEAR::wear_time_different),&list);

  Teuchos::setStringToIntegralParameter<int>("FRICTION","None","Type of friction law",
         Teuchos::tuple<std::string>("None","none","Stick","stick","Tresca","tresca","Coulomb","coulomb"),
         Teuchos::tuple<int>(INPAR::CONTACT::friction_none,INPAR::CONTACT::friction_none,
        INPAR::CONTACT::friction_stick,INPAR::CONTACT::friction_stick,
        INPAR::CONTACT::friction_tresca,INPAR::CONTACT::friction_tresca,
        INPAR::CONTACT::friction_coulomb,INPAR::CONTACT::friction_coulomb),&list);

  Teuchos::setStringToIntegralParameter<int>("FRLESS_FIRST","No","If chosen the first time step of a newly in contact slave node is regarded as frictionless",
         Teuchos::tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO"),
         Teuchos::tuple<int>(true,false,true,false,true,false),&list);

  Teuchos::setStringToIntegralParameter<int>("GP_SLIP_INCR","No","If chosen the slip increment is computed gp-wise which results to a non-objective quantity, but this would be consistent to wear and tsi calculations.",
         Teuchos::tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO"),
         Teuchos::tuple<int>(true,false,true,false,true,false),&list);

  DoubleParameter("FRCOEFF",0.0,"Friction coefficient for Coulomb friction",&list);

  return;
}

