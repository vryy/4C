/*----------------------------------------------------------------------*/
/*!
\file drt_validparameters.H

\brief Setup of the list of valid input parameters

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <iostream>

#include <Teuchos_Array.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <AztecOO.h>


#include "drt_validparameters.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
extern "C"
void PrintValidParameters()
{
  Teuchos::RCP<const Teuchos::ParameterList> list = DRT::ValidParameters();
  list->print(std::cout,
              Teuchos::ParameterList::PrintOptions()
              .showDoc(true)
              .showFlags(false)
              .indent(4)
              .showTypes(false));
  //std::cout << *list;
  //list->print(std::cout);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Teuchos::ParameterList> DRT::ValidParameters()
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::setIntParameter;
  using Teuchos::setDoubleParameter;
  using Teuchos::ParameterList;


  Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

  RCP<Teuchos::ParameterList> list = Teuchos::rcp(new Teuchos::ParameterList);

  /*----------------------------------------------------------------------*/
  ParameterList& discret = list->sublist("DISCRETISATION",false,"");

  setIntParameter("NUMFLUIDDIS",1,"Number of meshes in fluid field",&discret);
  setIntParameter("NUMSTRUCDIS",1,"Number of meshes in structural field",&discret);
  setIntParameter("NUMALEDIS",1,"Number of meshes in ale field",&discret);
  setIntParameter("NUMTHERMDIS",1,"Number of meshes in thermal field",&discret);

  /*----------------------------------------------------------------------*/
  ParameterList& size = list->sublist("PROBLEM SIZE",false,"");

  setIntParameter("ELEMENTS",0,"Total number of elements",&size);
  setIntParameter("NODES",0,"Total number of nodes",&size);
  setIntParameter("DIM",3,"2d or 3d problem",&size);
  setIntParameter("MATERIALS",0,"number of materials",&size);
  setIntParameter("NUMDF",3,"maximum number of degrees of freedom",&size);

  /*----------------------------------------------------------------------*/
  ParameterList& type = list->sublist("PROBLEM TYP",false,"");

  setStringToIntegralParameter("PROBLEMTYP","Fluid_Structure_Interaction","",
                               tuple<std::string>(
                                 "Structure",
                                 "Fluid",
                                 "Fluid_XFEM",
                                 "Convection_Diffusion",
                                 "Fluid_Structure_Interaction",
                                 "Ale",
                                 "Thermal_Structure_Interaction",
                                 "Structure_Multiscale"),
                                tuple<int>(
                                  prb_structure,
                                  prb_fluid,
                                  prb_fluid_xfem,
                                  prb_condif,
                                  prb_fsi,
                                  prb_ale,
                                  prb_tsi,
                                  prb_struct_multi),
                               &type);
  setIntParameter("NUMFIELD",1,"",&type);
  setStringToIntegralParameter("TIMETYP","Dynamic","",
                               tuple<std::string>("Static","Dynamic"),
                               tuple<int>(time_static,time_dynamic),
                               &type);
  //setIntParameter("GRADERW",0,"",&type);
  setIntParameter("MULTISC_STRUCT",0,"",&type);
  setIntParameter("RESTART",0,"",&type);
  setStringToIntegralParameter("ALGEBRA","Trilinos","outdated",
                               tuple<std::string>("Trilinos","ccarat"),
                               tuple<int>(1,0),
                               &type);

  /*----------------------------------------------------------------------*/
  ParameterList& io = list->sublist("IO",false,"");

  // are these needed?
  setStringToIntegralParameter("OUTPUT_OUT","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter("OUTPUT_GID","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter("OUTPUT_BIN","No","",yesnotuple,yesnovalue,&io);

  setStringToIntegralParameter("STRUCT_DISP","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter("STRUCT_STRESS","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter("STRUCT_STRESS_SMO","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter("STRUCT_SM_DISP","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter("STRUCT_SM_STRESS","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter("FLUID_SOL","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter("FLUID_STRESS","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter("FLUID_VIS","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter("ALE_DISP","No","",yesnotuple,yesnovalue,&io);

  setStringToIntegralParameter("THERM_TEMPERATURE","No","",yesnotuple,yesnovalue,&io);
  setStringToIntegralParameter("THERM_HEATFLUX","No","",yesnotuple,yesnovalue,&io);

  setIntParameter("FILESTEPS",1000,"",&io);

  //ParameterList& stat = list->sublist("STATIC",false,"");

  /*----------------------------------------------------------------------*/
  ParameterList& eigen = list->sublist("EIGENVALUE ANALYSIS",false,"");

  /*----------------------------------------------------------------------*/
  ParameterList& sdyn = list->sublist("STRUCTURAL DYNAMIC",false,"");

  setStringToIntegralParameter("DYNAMICTYP","Gen_Alfa",
                               "type of time integration control",
                               tuple<std::string>("Centr_Diff","Gen_EMM","Gen_Alfa"),
                               tuple<int>(STRUCT_DYNAMIC::centr_diff,STRUCT_DYNAMIC::Gen_EMM,STRUCT_DYNAMIC::gen_alfa),
                               &sdyn);

  setIntParameter("EIGEN",0,"EIGEN make eigenanalysis of the initial dynamic system",&sdyn);
  setIntParameter("RESEVRYDISP",1,"save displacements and contact forces every RESEVRYDISP steps",&sdyn);
  setIntParameter("RESEVRYSTRS",1,"save stresses every RESEVRYSTRS steps",&sdyn);
  setIntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&sdyn);
  setDoubleParameter("TIMESTEP",0.05,"time step size",&sdyn);
  setIntParameter("NUMSTEP",200,"maximum number of steps",&sdyn);
  setDoubleParameter("MAXTIME",5.0,"maximum time",&sdyn);
  setDoubleParameter("BETA",0.25,"generalized alpha factors, also used by explicit time integration",&sdyn);
  setDoubleParameter("GAMMA",0.5,"generalized alpha factors, also used by explicit time integration",&sdyn);
  setDoubleParameter("ALPHA_M",0.5,"generalized alpha factors",&sdyn);
  setDoubleParameter("ALPHA_F",0.5,"generalized alpha factors",&sdyn);
  setStringToIntegralParameter("DAMPING","No",
                               "build raleigh damping matrix and use it from M_DAMP x M + K_DAMP x K",
                               yesnotuple,yesnovalue,
                               &sdyn);
  setDoubleParameter("M_DAMP",0.5,"",&sdyn);
  setDoubleParameter("K_DAMP",0.5,"",&sdyn);

  setStringToIntegralParameter("ITERATION","full","unused",
                               tuple<std::string>("full","Full","FULL"),
                               tuple<int>(1,1,1),
                               &sdyn);

  setDoubleParameter("TOLDISP",1.0E-10,
                     "tolerance in the displacement norm for the newton iteration",
                     &sdyn);
  setIntParameter("MAXITER",50,
                  "maximum number of iterations allowed for newton iteration before failure",
                  &sdyn);
  setIntParameter("CONTACT",0,"contact algorithms",&sdyn);

  setStringToIntegralParameter("NLNSOL","fullnewton","",
                               tuple<std::string>(
                                 "fullnewton",
                                 "modnewton",
                                 "matfreenewton",
                                 "nlncg",
                                 "ptc"),
                               tuple<int>(
                                 STRUCT_DYNAMIC::fullnewton,
                                 STRUCT_DYNAMIC::modnewton,
                                 STRUCT_DYNAMIC::matfreenewton,
                                 STRUCT_DYNAMIC::nlncg,
                                 STRUCT_DYNAMIC::ptc),
                               &sdyn);

  setStringToIntegralParameter("PREDICT","ConstDis","",
                               tuple<std::string>(
                                 "Vague",
                                 "ConstDis",
                                 "ConstDisVelAcc"),
                               tuple<int>(
                                 STRUCT_DYNAMIC::pred_vague,
                                 STRUCT_DYNAMIC::pred_constdis,
                                 STRUCT_DYNAMIC::pred_constdisvelacc),
                               &sdyn);

  // time adaptivity (old style)
  setIntParameter("TIMEADAPT",0,"",&sdyn);
  setIntParameter("ITWANT",0,"",&sdyn);
  setDoubleParameter("MAXDT",0.0,"",&sdyn);
  setDoubleParameter("RESULTDT",0.0,"",&sdyn);

  SetValidTimeAdaptivityParameters(sdyn);

  /*----------------------------------------------------------------------*/
  ParameterList& fdyn = list->sublist("FLUID DYNAMIC",false,"");

  setStringToIntegralParameter("DYNAMICTYP","Nlin_Time_Int",
                               "Nonlinear Time Integraton Scheme",
                               tuple<std::string>("Nlin_Time_Int"),
                               tuple<int>(dyntyp_nln_time_int),
                               &fdyn);

  setStringToIntegralParameter("TIMEINTEGR","One_Step_Theta",
                               "Time Integration Scheme",
                               tuple<std::string>(
                                 "Stationary",
                                 "Gen_Alfa",
                                 "Gen_Alpha",
                                 "One_Step_Theta",
                                 "BDF2",
                                 "Inc_Acc_Gen_Alpha",
                                 "Theta_Adamsbashforth"
                                 ),
                               tuple<FLUID_TIMEINTTYPE>(
                                 timeint_stationary,
                                 timeint_gen_alpha,
                                 timeint_gen_alpha,
                                 timeint_one_step_theta,
                                 timeint_bdf2,
                                 timeint_inc_acc_gen_alpha,
                                 timeint_theta_adamsbashforth
                                 ),
                               &fdyn);
  setStringToIntegralParameter("STARTINGALGO","One_Step_Theta","",
                               tuple<std::string>(
                                 "One_Step_Theta"
                                 ),
                               tuple<int>(
                                 timeint_one_step_theta
                                 ),
                               &fdyn);
  setStringToIntegralParameter("NONLINITER","fixed_point_like",
                               "Nonlinear iteration scheme",
                               tuple<std::string>(
                                 "fixed_point_like",
                                 "Newton"
                                 ),
                               tuple<int>(1,2),
                               &fdyn);

  setStringToIntegralParameter("CONVCHECK","L_2_norm",
                               "norm for convergence check",
                               tuple<std::string>(
                                 "No",
                                 "L_infinity_norm",
                                 "L_1_norm",
                                 "L_2_norm"
                                 ),
                               tuple<int>(
                                 FLUID_DYNAMIC::fncc_no,
                                 FLUID_DYNAMIC::fncc_Linf,
                                 FLUID_DYNAMIC::fncc_L1,
                                 FLUID_DYNAMIC::fncc_L2
                                 ),
                               &fdyn);
  setStringToIntegralParameter("STEADYCHECK","L_2_norm",
                               "Norm of steady state check",
                               tuple<std::string>(
                                 "No",
                                 "L_infinity_norm",
                                 "L_1_norm",
                                 "L_2_norm"
                                 ),
                               tuple<int>(
                                 FLUID_DYNAMIC::fncc_no,
                                 FLUID_DYNAMIC::fncc_Linf,
                                 FLUID_DYNAMIC::fncc_L1,
                                 FLUID_DYNAMIC::fncc_L2
                                 ),
                               &fdyn);
  setStringToIntegralParameter("INITIALFIELD","zero_field",
                               "Initial Starting Field",
                               tuple<std::string>(
                                 "zero_field",
                                 "field_from_file",
                                 "field_by_function",
                                 "disturbed_field_from_function",
                                 "SOLWAVE",
                                 "WAVEBREAKING",
                                 "BELTRAMI-FLOW",
                                 "KIM-MOIN-FLOW",
                                 "BREAKING-DAM"),
                               tuple<int>(0,1,2,3,6,7,8,9,10),
                               &fdyn);

  setStringToIntegralParameter("VISCSTRESS","Yes",
                               "when calculating fluid stresses include viscose part (unused)",
                               yesnotuple,yesnovalue,&fdyn);

  setStringToIntegralParameter("FREESURFACE","No",
                               "Treatment of free surface",
                               tuple<std::string>(
                                 "No",
                                 "no",
                                 "loclag_exp",
                                 "loclag_imp",
                                 "hf_vert_sep",
                                 "hf_vert_imp",
                                 "genfs"
                                 ),
                               tuple<int>(0,0,1,2,3,5,6),
                               &fdyn);

  setStringToIntegralParameter("SURFTENSION","Yes",
                               "Include surfadce tension effects",
                               yesnotuple,yesnovalue,&fdyn);

  setStringToIntegralParameter("CHECKAREA","Yes",
                               "Monitor the size of the fluid domain",
                               yesnotuple,yesnovalue,&fdyn);

  setStringToIntegralParameter("LIFTDRAG","No",
                               "Calculate lift and drag forces along specified lines",
                               tuple<std::string>(
                                 "No",
                                 "no",
                                 "Yes",
                                 "yes",
                                 "Stress",
                                 "STRESS",
                                 "stress",
                                 "Nodeforce",
                                 "NODEFORCE",
                                 "nodeforce"
                                 ),
                               tuple<int>(
                                 FLUID_DYNAMIC::ld_none,
                                 FLUID_DYNAMIC::ld_none,
                                 FLUID_DYNAMIC::ld_stress,
                                 FLUID_DYNAMIC::ld_stress,
                                 FLUID_DYNAMIC::ld_stress,
                                 FLUID_DYNAMIC::ld_stress,
                                 FLUID_DYNAMIC::ld_stress,
                                 FLUID_DYNAMIC::ld_nodeforce,
                                 FLUID_DYNAMIC::ld_nodeforce,
                                 FLUID_DYNAMIC::ld_nodeforce
                                 ),
                               &fdyn);

  setStringToIntegralParameter("TURBULENCE","No","",
                               tuple<std::string>(
                                 "No",
                                 "algebraic",
                                 "kappa-eps",
                                 "kappa-omega",
                                 "VMM-LES"
                                 ),
                               tuple<int>(0,1,2,3,4),
                               &fdyn);

  setStringToIntegralParameter("HOMDIRECT","xy",
                               "homogeneous directions needed to evaluate in plane statistics for turbulent channel",
                               tuple<std::string>("xy","yz","xz"),
                               tuple<int>(2,0,1),
                               &fdyn);

  setStringToIntegralParameter("DISC_CAPT","No","",yesnotuple,yesnovalue,&fdyn);
  setStringToIntegralParameter("ADAPT_TIME","No",
                               "time stepping is adaptive",
                               yesnotuple,yesnovalue,&fdyn);

  setStringToIntegralParameter("CD_VELOCITY","Navier_Stokes","",
                               tuple<std::string>(
                                 "Navier_Stokes",
                                 "straight",
                                 "30_degree",
                                 "60_degree",
                                 "min60_degree"
                                 ),
                               tuple<int>(0,1,2,3,4),
                               &fdyn);

  setStringToIntegralParameter("CONVECTERM","convective","convective term",
                               tuple<std::string>(
                                 "convective",
                                 "divergence",
                                 "skew_symmetric"
                                 ),
                               tuple<int>(0,1,2),
                               &fdyn);
  setStringToIntegralParameter("VISCTERM","conventional","viscous term",
                               tuple<std::string>(
                                 "conventional",
                                 "stress_divergence"
                                 ),
                               tuple<int>(0,1),
                               &fdyn);

  setStringToIntegralParameter("SUBGRIDVISC","No","subgrid viscosity",
                               tuple<std::string>(
                                 "No",
                                 "artificial",
                                 "Smagorinsky"
                                 ),
                               tuple<int>(0,1,2),
                               &fdyn);

  setIntParameter("UPPSS",1,"Increment for visualisation (unused)",&fdyn);
  setIntParameter("UPOUT",1,"Increment for writing solution to output file",&fdyn);
  setIntParameter("UPRES",1,"Increment for writing solution",&fdyn);
  setIntParameter("RESSTEP",0,"Restart Step",&fdyn);
  setIntParameter("RESTARTEVRY",20,"Increment for writing restart",&fdyn);
  setIntParameter("NUMSTEP",1,"Total number of Timesteps",&fdyn);
  setIntParameter("STEADYSTEP",-1,"steady state check every step",&fdyn);
  setIntParameter("NUMSTASTEPS",0,"Number of Steps for Starting Scheme",&fdyn);
  setIntParameter("STARTFUNCNO",-1,"Function for Initial Starting Field",&fdyn);
  setIntParameter("ITEMAX",10,"max. number of nonlin. iterations",&fdyn);

  setDoubleParameter("TIMESTEP",0.01,"Time increment dt",&fdyn);
  setDoubleParameter("MAXTIME",1000.0,"Total simulation time",&fdyn);
  setDoubleParameter("ALPHA_M",1.0,"Time integration factor",&fdyn);
  setDoubleParameter("ALPHA_F",1.0,"Time integration factor",&fdyn);

  setDoubleParameter("THETA",0.66,"Time integration factor",&fdyn);

  setDoubleParameter("CONVTOL",1e-6,"Tolerance for convergence check",&fdyn);
  setDoubleParameter("STEADYTOL",1e-6,"Tolerance for steady state check",&fdyn);
  setDoubleParameter("START_THETA",1.0,"Time integraton factor for starting scheme",&fdyn);
  setDoubleParameter("INT_LENGHT",0.0,"",&fdyn);
  setDoubleParameter("ROUGHTNESS",0.0,"",&fdyn);
  setDoubleParameter("SC_COORD_X",0.0,"",&fdyn);
  setDoubleParameter("SC_COORD_Y",0.0,"",&fdyn);
  setDoubleParameter("MAX_DT",1.0,"Maximal Time increment dt_max in adaptive case",&fdyn);
  setDoubleParameter("MIN_DT",0.0,"Minimal Time increment dt_min in adaptive case",&fdyn);
  setDoubleParameter("LOC_TRUN_ERR",1e-3,"Local Truncation Error to rule adaptive time stepping",&fdyn);
  setDoubleParameter("SMAGCONST",0.0,"",&fdyn);

  /*----------------------------------------------------------------------*/
  ParameterList& adyn = list->sublist("ALE DYNAMIC",false,"");

  setDoubleParameter("TIMESTEP",0.1,"",&adyn);
  setIntParameter("NUMSTEP",41,"",&adyn);
  setDoubleParameter("MAXTIME",4.0,"",&adyn);
  setStringToIntegralParameter("ALE_TYPE","classic_lin","",
                               tuple<std::string>("classic_lin"),
                               tuple<int>(ALE_DYNAMIC::classic_lin),
                               &adyn);
  setIntParameter("NUM_INITSTEP",0,"",&adyn);
  setIntParameter("RESEVRYDISP",1,"",&adyn);

  setStringToIntegralParameter("QUALITY","none","unused",
                               tuple<std::string>("none","NONE"),
                               tuple<int>(
                                 ALE_DYNAMIC::no_quality,
                                 ALE_DYNAMIC::no_quality),
                               &adyn);

  /*----------------------------------------------------------------------*/
  ParameterList& fsidyn = list->sublist(
    "FSI DYNAMIC",false,
    "Fluid Structure Interaction\n"
    "Partitioned FSI solver with various coupling methods"
    );

  setStringToIntegralParameter("COUPALGO","iter_stagg_AITKEN_rel_param",
                               "Iteration Scheme over the fields",
                               tuple<std::string>(
                                 "basic_sequ_stagg",
                                 "sequ_stagg_pred",
                                 "sequ_stagg_shift",
                                 "iter_stagg_fixed_rel_param",
                                 "iter_stagg_AITKEN_rel_param",
                                 "iter_stagg_steep_desc",
                                 "iter_nox",
                                 "iter_monolithic"),
                               tuple<int>(
                                 fsi_basic_sequ_stagg,
                                 fsi_sequ_stagg_pred,
                                 fsi_sequ_stagg_shift,
                                 fsi_iter_stagg_fixed_rel_param,
                                 fsi_iter_stagg_AITKEN_rel_param,
                                 fsi_iter_stagg_steep_desc,
                                 fsi_iter_nox,
                                 fsi_iter_monolithic),
                               &fsidyn);

  setStringToIntegralParameter("PREDICTOR","d(n)",
                               "Predictor for interface displacements (unused)",
                               tuple<std::string>(
                                 "d(n)",
                                 "d(n)+dt*(1.5*v(n)-0.5*v(n-1))",
                                 "d(n)+dt*v(n)",
                                 "d(n)+dt*v(n)+0.5*dt^2*a(n)"
                                 ),
                               tuple<int>(1,2,3,4),
                               &fsidyn);

  setStringToIntegralParameter("CONVCRIT","||g(i)||:sqrt(neq)",
                               "Convergence criterium for iteration over fields (unused)",
                               tuple<std::string>(
                                 "||g(i)||:sqrt(neq)",
                                 "||g(i)||:||g(0)||"
                                 ),
                               tuple<int>(1,2),
                               &fsidyn);

  setStringToIntegralParameter("ENERGYCHECK","No",
                               "Energy check for iteration over fields",
                               yesnotuple,yesnovalue,&fsidyn);

  setStringToIntegralParameter("IALE","Pseudo_Structure",
                               "Treatment of ALE-field (outdated)",
                               tuple<std::string>(
                                 "Pseudo_Structure"
                                 ),
                               tuple<int>(1),
                               &fsidyn);

  setStringToIntegralParameter("COUPMETHOD","conforming",
                               "Coupling Method Mortar (mtr) or conforming nodes at interface (unused)",
                               tuple<std::string>(
                                 "MTR",
                                 "Mtr",
                                 "mtr",
                                 "conforming"
                                 ),
                               tuple<int>(0,0,0,1),
                               &fsidyn);

  setStringToIntegralParameter("COUPFORCE","nodeforce","",
                               tuple<std::string>(
                                 "none",
                                 "stress",
                                 "nodeforce"
                                 ),
                               tuple<int>(
                                 FSI_DYNAMIC::cf_none,
                                 FSI_DYNAMIC::cf_stress,
                                 FSI_DYNAMIC::cf_nodeforce),
                               &fsidyn);

  setIntParameter("ITECHAPP",1,"",&fsidyn);
  setIntParameter("ICHMAX",1,"",&fsidyn);
  setIntParameter("ISDMAX",1,"not used up to now",&fsidyn);
  setIntParameter("NUMSTEP",1,"Total number of Timesteps",&fsidyn);
  setIntParameter("ITEMAX",1,"Maximum number of iterations over fields",&fsidyn);
  setIntParameter("UPPSS",1,"Increment for visualisation",&fsidyn);
  setIntParameter("UPRES",1,"Increment for writing solution",&fsidyn);
  setIntParameter("RESTARTEVRY",1,"Increment for writing restart",&fsidyn);

  setDoubleParameter("TIMESTEP",0.1,"Time increment dt",&fsidyn);
  setDoubleParameter("MAXTIME",1000.0,"Total simulation time",&fsidyn);
  setDoubleParameter("TOLENCHECK",1e-6,"Tolerance for energy check",&fsidyn);
  setDoubleParameter("RELAX",1.0,"fixed relaxation parameter",&fsidyn);
  setDoubleParameter("CONVTOL",1e-6,"Tolerance for iteration over fields",&fsidyn);

  /*----------------------------------------------------------------------*/
  ParameterList& fluidsolver = list->sublist("FLUID SOLVER",false,"");
  SetValidSolverParameters(fluidsolver);

  /*----------------------------------------------------------------------*/
  ParameterList& structsolver = list->sublist("STRUCT SOLVER",false,"");
  SetValidSolverParameters(structsolver);

  /*----------------------------------------------------------------------*/
  ParameterList& alesolver = list->sublist("ALE SOLVER",false,"");
  SetValidSolverParameters(alesolver);

#ifdef D_TSI
  /*----------------------------------------------------------------------*/
  ParameterList& thermalsolver = list->sublist("THERMAL SOLVER",false,"");
  SetValidSolverParameters(thermalsolver);
#endif

  return list;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::SetValidSolverParameters(Teuchos::ParameterList& list)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::setIntParameter;
  using Teuchos::setDoubleParameter;
  using Teuchos::ParameterList;

  setStringToIntegralParameter("SOLVER","UMFPACK","",
                               tuple<std::string>(
                                 "Amesos_KLU_sym",
                                 "Amesos_KLU_nonsym",
                                 "Superlu",
                                 "vm3",
                                 "Aztec_MSR",
                                 "LAPACK_sym",
                                 "LAPACK_nonsym",
                                 "UMFPACK"
                                 ),
                               tuple<_SOLVER_TYP>(
                                 amesos_klu_sym,
                                 amesos_klu_nonsym,
                                 superlu,
                                 vm3,
                                 aztec_msr,
                                 lapack_sym,
                                 lapack_nonsym,
                                 umfpack),
                               &list);

  setStringToIntegralParameter(
    "AZSOLVE", "GMRES",
    "Type of linear solver algorithm to use.",
    tuple<std::string>("CG","GMRES","CGS","TFQMR","BiCGSTAB","LU"),
    tuple<_AZSOLVERTYP>(azsolv_CG,azsolv_GMRES,azsolv_CGS,azsolv_TFQMR,azsolv_BiCGSTAB,azsolv_LU),
    &list
    );
  setStringToIntegralParameter(
    "AZPREC", "ILU",
    "Type of internal preconditioner to use.\n"
    "Note! this preconditioner will only be used if the input operator\n"
    "supports the Epetra_RowMatrix interface and the client does not pass\n"
    "in an external preconditioner!",
    tuple<std::string>(
      "none",
      "ILU",
      "ILUT",
      "Jacobi",
      "SymmGaussSeidel",
      //"Polynomial",
      "Least_Squares",
      "Neumann",
      "ICC",
      "LU",
      "RILU").append("BILU").append("ML").append("MLFLUID").append("MLFLUID2").append("MLAPI"),
    //tuple<EAztecPreconditioner>(
    tuple<_AZPRECTYP>(
      azprec_none,
      azprec_ILU,
      azprec_ILUT,
      azprec_Jacobi,
      azprec_SymmGaussSeidel,
      //AZTEC_PREC_POLY,
      azprec_Least_Squares,
      azprec_Neumann,
      azprec_ICC,
      azprec_LU,
      azprec_RILU).append(azprec_BILU).append(azprec_ML).append(azprec_MLfluid).append(azprec_MLfluid2).append(azprec_MLAPI),
    &list
    );
  setIntParameter(
    "AZOVERLAP", 0,
    "The amount of overlap used for the internal \"ilu\" and \"ilut\" preconditioners.",
    &list
    );
  setIntParameter(
    "AZGFILL", 0,
    "The amount of fill allowed for the internal \"ilu\" preconditioner.",
    &list
    );
  setDoubleParameter(
    "AZDROP", 0.0,
    "The tolerance below which an entry from the factors of an internal \"ilut\"\n"
    "preconditioner will be dropped.",
    &list
    );
  setDoubleParameter(
    "AZFILL", 1.0,
    "The amount of fill allowed for an internal \"ilut\" preconditioner.",
    &list
    );
//   setIntParameter(
//     Steps_name, 3,
//     "Number of steps taken for the \"Jacobi\" or the \"Symmetric Gauss-Seidel\"\n"
//     "internal preconditioners for each preconditioner application.",
//     &list
//     );
  setIntParameter(
    "AZPOLY", 3,
    "The order for of the polynomials used for the \"Polynomial\" and\n"
    "\"Least-squares Polynomial\" internal preconditioners.",
    &list
    );
//   setStringToIntegralParameter(
//     RCMReordering_name, "Disabled",
//     "Determines if RCM reordering is used with the internal\n"
//     "\"ilu\" or \"ilut\" preconditioners.",
//     tuple<std::string>("Enabled","Disabled"),
//     tuple<int>(1,0),
//     &list
//     );
//   setStringToIntegralParameter(
//     Orthogonalization_name, "Classical",
//     "The type of orthogonalization to use with the \"GMRES\" solver.",
//     tuple<std::string>("Classical","Modified"),
//     tuple<int>(AZ_classic,AZ_modified),
//     &list
//     );
  setIntParameter(
    "AZSUB", 300,
    "The maximum size of the Krylov subspace used with \"GMRES\" before\n"
    "a restart is performed.",
    &list
    );
  setStringToIntegralParameter(
    "AZCONV", "AZ_r0", // Same as "rhs" when x=0
    "The convergence test to use for terminating the iterative solver.",
    tuple<std::string>(
      "AZ_r0",
      "AZ_rhs",
      "AZ_Anorm",
      "AZ_noscaled",
      "AZ_sol",
      "AZ_weighted",
      "AZ_expected_values",
      "AZTECOO_conv_test",
      "AZ_inf_noscaled"
      ),
    tuple<int>(
      AZ_r0,
      AZ_rhs,
      AZ_Anorm,
      AZ_noscaled,
      AZ_sol,
      AZ_weighted,
      AZ_expected_values,
      AZTECOO_conv_test,
      AZ_inf_noscaled
      ),
    &list
    );
//   setDoubleParameter(
//     IllConditioningThreshold_name, 1e+11,
//     "The threshold tolerance above which a system is considered\n"
//     "ill conditioned.",
//     &list
//     );
  setIntParameter(
    "AZOUTPUT", 0, // By default, no output from Aztec!
    "The number of iterations between each output of the solver's progress.",
    &list
    );

  setIntParameter("AZREUSE", 0, "how often to recompute some preconditioners", &list);
  setIntParameter("AZITER", 1000, "max iterations", &list);
  setIntParameter("AZGRAPH", 0, "unused", &list);
  setIntParameter("AZBDIAG", 0, "", &list);

  setDoubleParameter("AZTOL", 1e-8, "tolerance in (un)scaled residual", &list);
  setDoubleParameter("AZOMEGA", 0.0, "unused", &list);

  setStringToIntegralParameter(
    "AZSCAL","none","scaling of the system",
    tuple<std::string>("none","sym","infnorm"),
    tuple<int>(0,1,2),
    &list);

  // parameters of ML preconditioner

  setIntParameter("ML_PRINT",0,
                  "ML print-out level (0-10)",&list);
  setIntParameter("ML_MAXCOARSESIZE",5000,
                  "ML stop coarsening when coarse ndof smaller then this",&list);
  setIntParameter("ML_MAXLEVEL",5,
                  "ML max number of levels",&list);
  setIntParameter("ML_AGG_SIZE",27,
                  "objective size of an aggregate with METIS/VBMETIS, 2D: 9, 3D: 27",&list);

  setDoubleParameter("ML_DAMPFINE",1.,"damping fine grid",&list);
  setDoubleParameter("ML_DAMPMED",1.,"damping med grids",&list);
  setDoubleParameter("ML_DAMPCOARSE",1.,"damping coarse grid",&list);
  setDoubleParameter("ML_PROLONG_SMO",0.,"damping factor for prolongator smoother (usually 1.33 or 0.0)",&list);
  setDoubleParameter("ML_PROLONG_THRES",0.,"threshold for prolongator smoother/aggregation",&list);

  setNumericStringParameter("ML_SMOTIMES","1 1 1 1 1 1",
                            "no. smoothing steps or polynomial order on each level (at least ML_MAXLEVEL numbers)",&list);

  setStringToIntegralParameter(
    "ML_COARSEN","UC","",
    tuple<std::string>("UC","METIS","VBMETIS","MIS"),
    tuple<int>(0,1,2,3),
    &list);

  setStringToIntegralParameter(
    "ML_SMOOTHERFINE","ILU","",
    tuple<std::string>("SGS","Jacobi","Chebychev","MLS","ILU","KLU","Superlu"),
    tuple<int>(0,1,2,3,4,5,6),
    &list);

  setStringToIntegralParameter(
    "ML_SMOOTHERMED","ILU","",
    tuple<std::string>("SGS","Jacobi","Chebychev","MLS","ILU","KLU","Superlu"),
    tuple<int>(0,1,2,3,4,5,6),
    &list);

  setStringToIntegralParameter(
    "ML_SMOOTHERCOARSE","KLU","",
    tuple<std::string>("SGS","Jacobi","Chebychev","MLS","ILU","KLU","Superlu"),
    tuple<int>(0,1,2,3,4,5,6),
    &list);

  // unused
  setStringToIntegralParameter("PARTITION","Cut_Elements","unused",
                               tuple<std::string>("Cut_Elements"),
                               tuple<int>(0),
                               &list);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::SetValidTimeAdaptivityParameters(Teuchos::ParameterList& list)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::setIntParameter;
  using Teuchos::setDoubleParameter;
  using Teuchos::ParameterList;

  setStringToIntegralParameter(
    "TA_KIND","None","",
    tuple<std::string>("None","ZienkiewiczXie"),
    tuple<int>(
      TIMADA_DYNAMIC::timada_kind_none,
      TIMADA_DYNAMIC::timada_kind_zienxie),
    &list);

  setDoubleParameter("TA_STEPSIZEMAX", 0.0, "", &list);
  setDoubleParameter("TA_STEPSIZEMIN", 0.0, "", &list);
  setDoubleParameter("TA_SIZERATIOMAX", 0.0, "", &list);
  setDoubleParameter("TA_SIZERATIOMIN", 0.0, "", &list);
  setDoubleParameter("TA_SIZERATIOSCALE", 0.0, "", &list);

  setStringToIntegralParameter(
    "TA_ERRNORM", "Vague", "",
    tuple<std::string>("Vague",
                       "L1",
                       "L2",
                       "Rms",
                       "Inf"),
    tuple<int>(TIMADA_DYNAMIC::timada_err_norm_vague,
               TIMADA_DYNAMIC::timada_err_norm_l1,
               TIMADA_DYNAMIC::timada_err_norm_l2,
               TIMADA_DYNAMIC::timada_err_norm_rms,
               TIMADA_DYNAMIC::timada_err_norm_inf),
    &list);

  setDoubleParameter("TA_ERRTOL", 0.0, "", &list);
  setIntParameter("TA_ADAPTSTEPMAX", 0, "", &list);
}


#endif
