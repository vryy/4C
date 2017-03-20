/*----------------------------------------------------------------------*/
/*!
\file inpar_contact.cpp

\brief Input parameters for contact

\level 2

\maintainer Philipp Farah, Alexander Seitz

*/
/*----------------------------------------------------------------------*/
#include "drt_validparameters.H"
#include "inpar_contact.H"
#include "inpar_structure.H"



void INPAR::CONTACT::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

  /* parameters for structural meshtying and contact */
  Teuchos::ParameterList& scontact = list->sublist("CONTACT DYNAMIC",false,"");

  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for meshtying and contact",&scontact);

  setStringToIntegralParameter<int>("RESTART_WITH_CONTACT","No","Must be chosen if a non-contact simulation is to be restarted with contact",
                               yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<int>("ADHESION","None","Type of adhesion law",
      tuple<std::string>("None","none",
                         "bounded","b"),
      tuple<int>(
                 adhesion_none,adhesion_none,
                 adhesion_bound,adhesion_bound),
      &scontact);

  setStringToIntegralParameter<int>("DISCR_SMOOTHING","No","If chosen, interface smoothing with additional interface discr. is activated",
                               yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<int>("FRICTION","None","Type of friction law",
      tuple<std::string>("None","none",
                         "Stick","stick",
                         "Tresca","tresca",
                         "Coulomb","coulomb"),
      tuple<int>(
                 friction_none,friction_none,
                 friction_stick,friction_stick,
                 friction_tresca,friction_tresca,
                 friction_coulomb,friction_coulomb),
      &scontact);

  setStringToIntegralParameter<int>("FRLESS_FIRST","No",
      "If chosen the first time step of a newly in contact slave node is regarded as frictionless",
                               yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<int>("GP_SLIP_INCR","No",
      "If chosen the slip increment is computed gp-wise which results to a non-objective quantity, but this would be consistent to wear and tsi calculations.",
                               yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<int>("STRATEGY","LagrangianMultipliers","Type of employed solving strategy",
        tuple<std::string>("LagrangianMultipliers","lagrange", "Lagrange",
                           "PenaltyMethod","penalty", "Penalty",
                           "UzawaAugementedLagrange","uzawa","Uzawa",
                           "Augmented",
                           "SteepestAscent",
                           "XContact",
                           "Nitsche"),
        tuple<int>(
                solution_lagmult, solution_lagmult, solution_lagmult,
                solution_penalty, solution_penalty, solution_penalty,
                solution_uzawa, solution_uzawa, solution_uzawa,
                solution_augmented,
                solution_steepest_ascent,
                solution_xcontact,
                solution_nitsche),
        &scontact);

  setStringToIntegralParameter<int>("SYSTEM","Condensed","Type of linear system setup / solution",
        tuple<std::string>("Condensed","condensed", "cond",
                           "Condensedlagmult","condensedlagmult","condlm",
                           "SaddlePoint","Saddlepoint","saddlepoint", "sp","none"),
        tuple<int>(
                system_condensed, system_condensed, system_condensed,
                system_condensed_lagmult,system_condensed_lagmult,system_condensed_lagmult,
                system_saddlepoint, system_saddlepoint,
                system_saddlepoint, system_saddlepoint,
                system_none),
        &scontact);

  DoubleParameter("PENALTYPARAM",0.0,"Penalty parameter for penalty / Uzawa augmented solution strategy",&scontact);
  DoubleParameter("PENALTYPARAMTAN",0.0,"Tangential penalty parameter for penalty / Uzawa augmented solution strategy",&scontact);
  IntParameter("UZAWAMAXSTEPS",10,"Maximum no. of Uzawa steps for Uzawa solution strategy",&scontact);
  DoubleParameter("UZAWACONSTRTOL",1.0e-8,"Tolerance of constraint norm for Uzawa solution strategy",&scontact);

  setStringToIntegralParameter<int>("SEMI_SMOOTH_NEWTON","Yes","If chosen semi-smooth Newton concept is applied",
                               yesnotuple,yesnovalue,&scontact);

  DoubleParameter("SEMI_SMOOTH_CN",1.0,"Weighting factor cn for semi-smooth PDASS",&scontact);
  DoubleParameter("SEMI_SMOOTH_CT",1.0,"Weighting factor ct for semi-smooth PDASS",&scontact);

  setStringToIntegralParameter<int>(
      "CONTACTFORCE_ENDTIME",
      "No",
      "If chosen, the contact force is not evaluated at the generalized midpoint, but at the end of the time step",
      yesnotuple, yesnovalue, &scontact);

  setStringToIntegralParameter<int>("VELOCITY_UPDATE","No","If chosen, velocity update method is applied",
                               yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<int>("EMOUTPUT","None","Type of energy and momentum output",
      tuple<std::string>("None","none", "No", "no",
                         "Screen", "screen",
                         "File", "file",
                         "Both", "both"),
      tuple<int>(
              output_none, output_none,
              output_none, output_none,
              output_screen, output_screen,
              output_file, output_file,
              output_both, output_both),
      &scontact);

  setStringToIntegralParameter<int>("ERROR_NORMS","None","Choice of analytical solution for error norm computation",
      tuple<std::string>("None","none", "No", "no",
                         "Zero", "zero",
                         "Bending", "bending",
                         "Sphere", "sphere",
                         "Thick", "thick",
                         "Plate", "plate"),
      tuple<int>(
              errornorms_none, errornorms_none,
              errornorms_none, errornorms_none,
              errornorms_zero, errornorms_zero,
              errornorms_bending, errornorms_bending,
              errornorms_sphere, errornorms_sphere,
              errornorms_thicksphere, errornorms_thicksphere,
              errornorms_infiniteplate,errornorms_infiniteplate),
      &scontact);

  setStringToIntegralParameter<int>("INITCONTACTBYGAP","No","Initialize init contact by weighted gap vector",
                               yesnotuple,yesnovalue,&scontact);

  DoubleParameter("INITCONTACTGAPVALUE",0.0,"Value for initialization of init contact set with gap vector",&scontact);

  // solver convergence test parameters for contact/meshtying in saddlepoint formulation
  setStringToIntegralParameter<int>("NORMCOMBI_RESFCONTCONSTR","And",
    "binary operator to combine contact constraints and residual force values",
    tuple<std::string>(
      "And",
      "Or"),
    tuple<int>(
      INPAR::STR::bop_and,
      INPAR::STR::bop_or),
    &scontact
    );

  setStringToIntegralParameter<int>("NORMCOMBI_DISPLAGR","And",
      "binary operator to combine displacement increments and Lagrange multiplier increment values",
      tuple<std::string>(
        "And",
        "Or"),
      tuple<int>(
        INPAR::STR::bop_and,
        INPAR::STR::bop_or),
      &scontact
      );

  DoubleParameter("TOLCONTCONSTR",1.0E-6,
                  "tolerance in the contact constraint norm for the newton iteration (saddlepoint formulation only)",
                  &scontact);
  DoubleParameter("TOLLAGR",1.0E-6,
                  "tolerance in the LM norm for the newton iteration (saddlepoint formulation only)",
                  &scontact);

  setStringToIntegralParameter<int>("MESH_ADAPTIVE_CN","no",
                                     "use a scaling of cn with the local mesh size",
                                     yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<int>("MESH_ADAPTIVE_CT","no",
                                     "use a scaling of ct with the local mesh size",
                                     yesnotuple,yesnovalue,&scontact);

  setStringToIntegralParameter<int>("CONSTRAINT_DIRECTIONS","ntt",
      "formulation of constraints in normal/tangential or xyz-direction",
      tuple<std::string>(
        "ntt",
        "xyz"),
      tuple<int>(
        constr_ntt,
        constr_xyz),
      &scontact
      );

  setStringToIntegralParameter<int>("CONTACT_REGULARIZATION","no",
      "use regularized contact",
      tuple<std::string>(
        "no",
        "tanh"),
      tuple<int>(
        reg_none,
        reg_tanh),
      &scontact
      );

  setStringToIntegralParameter<int>("NONSMOOTH_GEOMETRIES","No",
      "If chosen the contact algorithm combines mortar and nts formulations.",
                               yesnotuple,yesnovalue,&scontact);

  DoubleParameter("HYBRID_ANGLE_MIN",-1.0,"Non-smooth contact: angle between cpp normal and element normal: begin transition (Mortar)",&scontact);
  DoubleParameter("HYBRID_ANGLE_MAX",-1.0,"Non-smooth contact: angle between cpp normal and element normal: end transition (NTS)",&scontact);

  setStringToIntegralParameter<int>("CPP_NORMALS","No",
      "If chosen the nodal normal field is created as averaged CPP normal field.",
                               yesnotuple,yesnovalue,&scontact);

  // --------------------------------------------------------------------------
  // sub-list "Augmented"
  Teuchos::ParameterList& augcontact=scontact.sublist("AUGMENTED");

  setStringToIntegralParameter<int>("PRINT_LINEAR_CONSERVATION","No",
      "Do and print the linear momentum conservation check.",
      yesnotuple,yesnovalue,&augcontact);

  setStringToIntegralParameter<int>("PRINT_ANGULAR_CONSERVATION","No",
      "Do and print the angular momentum conservation check.",
      yesnotuple,yesnovalue,&augcontact);

  // --------------------------------------------------------------------------
  // sub-sub-list "Augmented/SteepestAscent"
  Teuchos::ParameterList& sacontact=augcontact.sublist("STEEPESTASCENT");

  DoubleParameter("CN_UPPER_BOUND",std::numeric_limits<double>::max(),
      "Upper bound for the cn value. Used during the dynamic cn update.",&sacontact);

  // --------------------------------------------------------------------------
  // sub-list "eXtended contact formulation"
    Teuchos::ParameterList& xcontact=scontact.sublist("XCONTACT");
  // TODO
  setStringToIntegralParameter<int>("CONST_CPP_NORMAL", "No",
      "If chosen, closest point normal on master is assumed to be constant during"
      " variation and linearization.",
      yesnotuple, yesnovalue, &xcontact);

  // TODO
  setStringToIntegralParameter<int>("H1_DUALITY_PAIRING", "Yes",
      "If chosen, H1 duality pairing for contact potential is used.",
      yesnotuple, yesnovalue, &xcontact);

  // --------------------------------------------------------------------------
  DoubleParameter("NITSCHE_THETA",0.0,"+1: symmetric, 0: non-symmetric, -1: skew-symmetric",&scontact);

  setStringToIntegralParameter<int>("NITSCHE_WEIGHTING","harmonic",
      "how to weight consistency terms in Nitsche contact formulation",
      tuple<std::string>(
        "slave",
        "master",
        "harmonic"),
      tuple<int>(
        NitWgt_slave,
        NitWgt_master,
        NitWgt_harmonic),
      &scontact
      );

  setStringToIntegralParameter<int>("NITSCHE_PENALTY_ADAPTIVE","yes",
      "adapt penalty parameter after each converged time step",
      yesnotuple,yesnovalue,&scontact);

}


