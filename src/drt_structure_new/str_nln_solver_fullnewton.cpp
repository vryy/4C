/*
 * str_nln_solver_fullnewton.cpp
 *
 *  Created on: Sep 18, 2015
 *      Author: hiermeier
 */


#include "str_nln_solver_fullnewton.H"
#include "str_timint_base.H"

#include <Teuchos_XMLParameterListCoreHelpers.hpp>




/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::NLN::SOLVER::FullNewton::FullNewton()
    : Nox()
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::FullNewton::Setup()
{
  CheckInit();

  // setup the nox parameter list for a full Newton solution method
  SetFullNewtonParams();

  // Call the Setup() function of the base class
  // Note, that the issetup_ flag is also updated during this call.
  STR::NLN::SOLVER::Nox::Setup();

  if (not IsSetup())
    dserror("Should be \"true\" at this point!");

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::FullNewton::SetFullNewtonParams()
{
  CheckInit();

  // get the nox parameter list and set the necessary parameters for a
  // full Newton solution procedure
  Teuchos::ParameterList& p = DataSDyn().GetMutableNoxParams();

  // ---------------------------------------------------------------------------
  // Set-up the full Newton method
  // ---------------------------------------------------------------------------
  // Direction
  Teuchos::ParameterList& pdir = p.sublist("Direction",true);
  pdir.set("Method","Newton");

  // Line Search
  Teuchos::ParameterList& plinesearch = p.sublist("Line Search",true);
  plinesearch.set("Method","Full Step");

  // Line Search/Full Step
  Teuchos::ParameterList& pfullstep = plinesearch.sublist("Full Step",true);
  // check if the default value is set
  double fullstep = pfullstep.get<double>("Full Step");
  if (fullstep!=1.0)
  {
    std::string markerline;
    markerline.assign(40,'!');
    std::cout << markerline << std::endl
        << "WARNING: The Full Step length is " << fullstep
        << " (default=1.0)"<< std::endl
        << markerline << std::endl;
  }

  // ---------------------------------------------------------------------------
  // STATUS TEST
  // ---------------------------------------------------------------------------
  /*
   * This is only necessary for the special case, that you use no xml-file for
   * the definition of your convergence tests, but you use the dat-file instead.
   */
  const Teuchos::ParameterList& statustestparams =
      DataSDyn().GetNoxParams().sublist("Status Test");

  Teuchos::ParameterList xmlParams;
  const std::string& xmlfilename =
      statustestparams.get<const std::string>("XML File");

  // check the input: path to the "Status Test" xml-file
  if (xmlfilename != "none")
  {
    if (xmlfilename.length() and xmlfilename.rfind(".xml"))
      xmlParams = *(Teuchos::getParametersFromXmlFile(xmlfilename));
    else
      dserror("The file name \"%s\" is not a valid XML file name.",
          xmlfilename.c_str());
  }
  else
  {
    std::set<enum NOX::NLN::StatusTest::QuantityType> qtypes;
    CreateQuantityTypes(qtypes);
    SetStatusTestParams(qtypes);
  }

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::FullNewton::CreateQuantityTypes(
    std::set<enum NOX::NLN::StatusTest::QuantityType>& qtypes) const
{
  CheckInit();

  // ---------------------------------------------------------------------------
  // get the model types
  // ---------------------------------------------------------------------------
  const std::set<enum INPAR::STR::ModelType>& mtypes =
      DataSDyn().GetModelTypes();

  std::set<enum INPAR::STR::ModelType>::const_iterator miter;
  for (miter=mtypes.begin();miter!=mtypes.end();++miter)
  {
    const std::string name = INPAR::STR::ModelTypeString(*miter);
    enum NOX::NLN::StatusTest::QuantityType qtype =
            NOX::NLN::StatusTest::String2QuantityType(name);

    // check if the corresponding enum could be found.
    // Note: We do not throw an error, since there are model types,
    //       which have no representation in the quantity list
    //       (e.g. spring dashpot).
    if (qtype == NOX::NLN::StatusTest::quantity_unknown)
      continue;

    qtypes.insert(qtype);
  }
  // ---------------------------------------------------------------------------
  // get the element technologies
  // ---------------------------------------------------------------------------
  const std::set<enum INPAR::STR::EleTech>& eletechs =
      DataSDyn().GetElementTechnologies();

  std::set<enum INPAR::STR::EleTech>::const_iterator etiter;
  for (etiter=eletechs.begin();etiter!=eletechs.end();++etiter)
  {
    const std::string name = INPAR::STR::EleTechString(*etiter);
    const enum NOX::NLN::StatusTest::QuantityType qtype =
            NOX::NLN::StatusTest::String2QuantityType(name);

    // check if the corresponding enum could be found.
    // Note: We do not throw an error, since there are element technologies,
    //       which have no representation in the quantity list
    //       (e.g. fbar).
    if (qtype == NOX::NLN::StatusTest::quantity_unknown)
      continue;

    qtypes.insert(qtype);
  }

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::FullNewton::SetStatusTestParams(
    const std::set<enum NOX::NLN::StatusTest::QuantityType>& qt)
{
  CheckInit();

  // get the status test parameter list
  Teuchos::ParameterList& pstatus =
      DataSDyn().GetMutableNoxParams().sublist("Status Test",true);

  // ------ outer status test ------
  // Combo test: OR-combination:
  // Combines the maximum iteration and the remaining tests
  Teuchos::ParameterList& postatus = pstatus.sublist("Outer Status Test");
  postatus.set("Test Type","Combo");
  postatus.set("Combo Type","OR");

  // ---------------------------------------------------------------------------
  // lvl. 0: combo OR - Test 0:
  // Combo test: AND
  //             combination of increment AND/OR residual
  //    AND      constraints
  //    AND      windkessel
  //    AND      windkessel increment
  //    AND      contact    (active set)
  //    AND      plasticity (active set)
  // ---------------------------------------------------------------------------
  Teuchos::ParameterList& pcombo_incr_fres_constr =
      postatus.sublist("Test 0");
  pcombo_incr_fres_constr.set("Test Type","Combo");
  pcombo_incr_fres_constr.set("Combo Type","AND");

  // ---------------------------------------------------------------------------
  // | lvl. 1: combo AND - Test 0:
  // | Combo test: AND/OR
  // | Combines increment (i.e. displacement and others) AND/OR
  // | force residual tests
  // ---------------------------------------------------------------------------
  Teuchos::ParameterList& pcombo_incr_fres = pcombo_incr_fres_constr.sublist("Test 0");
  pcombo_incr_fres.set("Test Type","Combo");
  switch (DataSDyn().GetResIncrComboType(NOX::NLN::StatusTest::quantity_structure,
      NOX::NLN::StatusTest::quantity_structure))
  {
    case INPAR::STR::bop_and:
      pcombo_incr_fres.set("Combo Type","AND");
      break;
    case INPAR::STR::bop_or:
      pcombo_incr_fres.set("Combo Type","OR");
      break;
    default:
      dserror("Unknown structural combination type enum!");
      break;
  }

  // ---------------------------------------------------------------------------
  // || lvl. 2: combo AND/OR - Test 0:
  // || Tests the combination of the different solution increment norms
  // ---------------------------------------------------------------------------
  SetComboQuantityTestParams(pcombo_incr_fres,0,"NormIncr",qt);

  // ---------------------------------------------------------------------------
  // || lvl. 2: combo AND/OR - Test 1:
  // || Tests the combination of the different solution force norms
  // ---------------------------------------------------------------------------
  SetComboQuantityTestParams(pcombo_incr_fres,1,"NormF",qt);

  // ---------------------------------------------------------------------------
  // | lvl. 1: combo AND - Test 1:
  // | Enforces a minimal number of non-linear solution steps
  // ---------------------------------------------------------------------------
  Teuchos::ParameterList& pnstep = pcombo_incr_fres_constr.sublist("Test 1");
  pnstep.set("Test Type","NStep");
  pnstep.set("Number of Nonlinear Iterations",DataSDyn().GetIterMin());

  int opt_count = 2;
  // *** BEGIN: OPTIONAL STATUS TESTS ******************************************
  // ---------------------------------------------------------------------------
  // | lvl. 1: combo AND - Test OPTIONAL:
  // | Tests the constraint NormF
  // ---------------------------------------------------------------------------
  if (qt.find(NOX::NLN::StatusTest::quantity_lag_pen_constraint)!=qt.end())
  {
    SetQuantityTestParams(pcombo_incr_fres_constr,
        NOX::NLN::StatusTest::quantity_lag_pen_constraint,opt_count,"NormF");
    ++opt_count;
  }

  // ---------------------------------------------------------------------------
  // | lvl. 1: combo AND - Test OPTIONAL:
  // | Tests the windkessel NormF
  // ---------------------------------------------------------------------------
  if (qt.find(NOX::NLN::StatusTest::quantity_windkessel) != qt.end())
  {
    SetQuantityTestParams(pcombo_incr_fres_constr,
        NOX::NLN::StatusTest::quantity_windkessel,opt_count,"NormF");
    ++opt_count;

  // ---------------------------------------------------------------------------
  // | lvl. 1: combo AND - Test OPTIONAL:
  // | Tests the windkessel NormIncr
  // ---------------------------------------------------------------------------
    SetQuantityTestParams(pcombo_incr_fres_constr,
        NOX::NLN::StatusTest::quantity_windkessel,opt_count,"NormIncr");
    ++opt_count;
  }

  // ---------------------------------------------------------------------------
  // | lvl. 1: combo AND - Test OPTIONAL:
  // | Tests the semi-smooth contact active set
  // ---------------------------------------------------------------------------
  if (qt.find(NOX::NLN::StatusTest::quantity_contact) != qt.end())
  {
    SetQuantityTestParams(pcombo_incr_fres_constr,
        NOX::NLN::StatusTest::quantity_contact,opt_count,"ActiveSet");
    ++opt_count;
  }

  // ---------------------------------------------------------------------------
  // | lvl. 1: combo AND - Test OPTIONAL:
  // | Tests the semi-smooth plasticity active set
  // ---------------------------------------------------------------------------
  if (qt.find(NOX::NLN::StatusTest::quantity_plasticity) != qt.end())
  {
    SetQuantityTestParams(pcombo_incr_fres_constr,
        NOX::NLN::StatusTest::quantity_plasticity,opt_count,"ActiveSet");
    ++opt_count;
  }
  // *** END: OPTIONAL STATUS TESTS ********************************************

  // ---------------------------------------------------------------------------
  // lvl. 0: combo OR - Test 1:
  // Tests the maximal iteration counter
  // ---------------------------------------------------------------------------
  Teuchos::ParameterList& pmaxiters = postatus.sublist("Test 1");
  pmaxiters.set("Test Type","MaxIters");
  pmaxiters.set("Maximum Iterations",DataSDyn().GetIterMax());

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::FullNewton::SetComboQuantityTestParams(
    Teuchos::ParameterList& p,
    const std::size_t& count,
    const std::string& testname,
    const std::set<enum NOX::NLN::StatusTest::QuantityType>& qtypes)
{
  CheckInit();

  std::vector<enum NOX::NLN::StatusTest::QuantityType> combo_or(0);
  std::vector<enum NOX::NLN::StatusTest::QuantityType> combo_and(0);
  SplitAndOrCombo(combo_or,combo_and,testname,qtypes);
  std::ostringstream test_string;
  test_string << "Test " << count;
  Teuchos::ParameterList& ptest = p.sublist(test_string.str());

  std::vector<enum NOX::NLN::StatusTest::QuantityType>::const_iterator qtiter;
  // if there are any OR combinations
  std::size_t count_or = 0;
  std::size_t count_and = 0;
  if (combo_or.size()>0)
  {
    ptest.set("Test Type","Combo");
    ptest.set("Combo Type","OR");

    std::ostringstream test_string;
    for (qtiter=combo_or.begin();qtiter!=combo_or.end();++qtiter)
    {
      SetQuantityTestParams(ptest,*qtiter,count_or,testname);
      ++count_or;
    }
    // if there are only OR combinations
    if (combo_and.size()==0)
      SetQuantityTestParams(ptest,
          NOX::NLN::StatusTest::quantity_structure,count_or,testname);
  }
  // if there are any AND combinations
  if (combo_and.size()>0)
  {
    Teuchos::ParameterList& ptest_and = ptest;
    if (count_or>0)
    {
      test_string.str("");
      test_string << "Test " << count_or;
      ptest_and = ptest.sublist(test_string.str());
    }
    ptest_and.set("Test Type","Combo");
    ptest_and.set("Combo Type","AND");
    for (qtiter=combo_and.begin();qtiter!=combo_and.end();++qtiter)
    {
      SetQuantityTestParams(ptest_and,*qtiter,count_and,testname);
      ++count_and;
    }
  }
  // if there are neither AND nor OR combinations
  if (count_or == 0 and count_and == 0)
    SetQuantityTestParams(ptest,
        NOX::NLN::StatusTest::quantity_structure,testname);

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::FullNewton::SetQuantityTestParams(
    Teuchos::ParameterList& p,
    const enum NOX::NLN::StatusTest::QuantityType& qtype,
    const std::size_t& count,
    const std::string& testname
    )
{
  CheckInit();

  std::ostringstream test_string;
  test_string << "Test " << count;
  Teuchos::ParameterList& ptest = p.sublist(test_string.str());
  SetQuantityTestParams(ptest,qtype,testname);

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::FullNewton::SetQuantityTestParams(
    Teuchos::ParameterList& p,
    const enum NOX::NLN::StatusTest::QuantityType& qtype,
    const std::string& testname
    )
{
  CheckInit();


  if (testname == "NormIncr")
    SetNormIncrParams(p,qtype,
        DataSDyn().GetIncrToleranceType(qtype),
        DataSDyn().GetIncrTolerance(qtype),
        DataSDyn().GetNormType());
  else if (testname == "NormF")
    SetNormFParams(p,qtype,
        DataSDyn().GetResToleranceType(qtype),
        DataSDyn().GetResTolerance(qtype),
        DataSDyn().GetNormType());
  else if (testname == "ActiveSet")
    SetActiveSetParams(p,qtype);
  else
    dserror("Unknown/Unsupported status test name: %s",testname.c_str());

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::FullNewton::SplitAndOrCombo(
    std::vector<enum NOX::NLN::StatusTest::QuantityType>& combo_or,
    std::vector<enum NOX::NLN::StatusTest::QuantityType>& combo_and,
    const std::string& testname,
    const std::set<enum NOX::NLN::StatusTest::QuantityType>& qtypes
    )
{
  CheckInit();

  std::set<enum NOX::NLN::StatusTest::QuantityType>::const_iterator qtiter;

  for (qtiter=qtypes.begin();qtiter!=qtypes.end();++qtiter)
  {
    if (*qtiter == NOX::NLN::StatusTest::quantity_structure)
      continue;

    enum INPAR::STR::BinaryOp combotype =
        DataSDyn().GetIncrComboType(*qtiter);
    if (testname=="NormF")
      combotype = DataSDyn().GetResComboType(*qtiter);
    else if (testname!="NormIncr")
      dserror("The given test \"%s\" name is not supported!",testname);

    switch (combotype)
    {
      case INPAR::STR::bop_or:
        combo_or.push_back(*qtiter);
        break;
      case INPAR::STR::bop_and:
        combo_and.push_back(*qtiter);
        break;
      default:
        dserror("Unknown combination type. See list of valid "
            "\"INPAR::STR::BinaryOp\" enums for more information.");
        break;
    } // switch case
  } // loop over the model type vector

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::FullNewton::SetNormIncrParams(
    Teuchos::ParameterList& qlist,
    const enum NOX::NLN::StatusTest::QuantityType& qtype,
    const enum INPAR::STR::ConvNorm& toltype,
    const double& tol,
    const enum INPAR::STR::VectorNorm& normtype
    )
{
  if (!IsInit())
      dserror("You have to call Init() first.");

  SetNormIncrParams(qlist,qtype,1.0,0.5,toltype,tol,normtype,false);
  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::FullNewton::SetNormIncrParams(
    Teuchos::ParameterList& qlist,
    const enum NOX::NLN::StatusTest::QuantityType& qtype,
    const double& alpha,
    const double& beta,
    const enum INPAR::STR::ConvNorm& toltype,
    const double& tol,
    const enum INPAR::STR::VectorNorm& normtype,
    const bool& isscaled
    )
{
  CheckInit();

  /* Set the tolerance type
   * Be careful: This has to be done in first place because of the special
   *             treatment of the mixed tolerance type!
   */
  switch (toltype)
  {
    // ABSOLUTE TOLERANCE TYPE
    case INPAR::STR::convnorm_abs:
    {
      qlist.set("Test Type","NormIncr");
      qlist.set("Quantity Type",NOX::NLN::StatusTest::QuantityType2String(qtype).c_str());
      qlist.set("Tolerance Type","Absolute");
      break;
    }
    // RELATIVE TOLERANCE TYPE
    case INPAR::STR::convnorm_rel:
    {
      qlist.set("Test Type","NormIncr");
      qlist.set("Quantity Type",NOX::NLN::StatusTest::QuantityType2String(qtype).c_str());
      qlist.set("Tolerance Type","Relative");
      break;
    }
    // MIX TOLERANCE TYPE
    /* Actually, this should not be part of the tolerance type decision, but due
     * to historical reasons, we still support it. If you consider to implement a
     * related status test, please use a xml file in combination with the combo
     * status test instead.
     */
    case INPAR::STR::convnorm_mix:
    {
      qlist.set("Test Type","Combo");
      qlist.set("Combo Type","OR");
      Teuchos::ParameterList& qlist_abs = qlist.sublist("Test 0");
      // first recursive call
      SetNormIncrParams(qlist_abs,qtype,INPAR::STR::convnorm_abs,tol,normtype);
      Teuchos::ParameterList& qlist_rel = qlist.sublist("Test 1");
      // second recursive call
      SetNormIncrParams(qlist_rel,qtype,INPAR::STR::convnorm_rel,tol,normtype);
      break;
    }
    default:
      dserror("Unknown \"Tolerance Type\" defined!");
      break;
  } // switch case toltype

  // set tolerance
  qlist.set<double>("Tolerance",tol);

  // set norm type
  switch (normtype)
  {
    case INPAR::STR::norm_l2:
      qlist.set("Norm Type","Two Norm");
      break;
    case INPAR::STR::norm_l1:
      qlist.set("Norm Type","One Norm");
      break;
    case INPAR::STR::norm_inf:
      qlist.set("Norm Type","Max Norm");
      break;
    case INPAR::STR::norm_rms:
      dserror("The norm type \"Root Mean Square\" is no longer supported! "
          "Consider to use the \"NOX::NLN::StatusTest::NormWRMS\" test instead!");
      break;
    case INPAR::STR::norm_vague:
    default:
      dserror("Unknown vector norm type!");
      break;
  } // switch case normtype

  // set alpha
  qlist.set<double>("Alpha",alpha);

  // set beta
  qlist.set<double>("Beta",beta);

  // set scaling by the vector length
  if (isscaled)
    qlist.set("Scale Type","Scaled");

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::FullNewton::SetNormFParams(
    Teuchos::ParameterList& qlist,
    const enum NOX::NLN::StatusTest::QuantityType& qtype,
    const enum INPAR::STR::ConvNorm& toltype,
    const double& tol,
    const enum INPAR::STR::VectorNorm& normtype
    )
{
  CheckInit();

  SetNormFParams(qlist,qtype,toltype,tol,normtype,false);
  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::FullNewton::SetNormFParams(
    Teuchos::ParameterList& qlist,
    const enum NOX::NLN::StatusTest::QuantityType& qtype,
    const enum INPAR::STR::ConvNorm& toltype,
    const double& tol,
    const enum INPAR::STR::VectorNorm& normtype,
    const bool& isscaled
    )
{
  CheckInit();

  /* Set the tolerance type
   * Be careful: This has to be done in first place because of the special
   *             treatment of the mixed tolerance type!
   */
  switch (toltype)
  {
    // ABSOLUTE TOLERANCE TYPE
    case INPAR::STR::convnorm_abs:
    {
      qlist.set("Test Type","NormF");
      qlist.set("Quantity Type",NOX::NLN::StatusTest::QuantityType2String(qtype).c_str());
      qlist.set("Tolerance Type","Absolute");
      break;
    }
    // RELATIVE TOLERANCE TYPE
    case INPAR::STR::convnorm_rel:
    {
      qlist.set("Test Type","NormF");
      qlist.set("Quantity Type",NOX::NLN::StatusTest::QuantityType2String(qtype).c_str());
      qlist.set("Tolerance Type","Relative");
      break;
    }
    // MIX TOLERANCE TYPE
    /* Actually, this should not be part of the tolerance type decision, but due
     * to historical reasons, we still support it. If you consider to implement a
     * related status test, please use a xml file in combination with the combo
     * status test instead.
     */
    case INPAR::STR::convnorm_mix:
    {
      qlist.set("Test Type","Combo");
      qlist.set("Combo Type","OR");
      Teuchos::ParameterList& qlist_abs = qlist.sublist("Test 0");
      // first recursive call
      SetNormFParams(qlist_abs,qtype,INPAR::STR::convnorm_abs,tol,normtype);
      Teuchos::ParameterList& qlist_rel = qlist.sublist("Test 1");
      // second recursive call
      SetNormFParams(qlist_rel,qtype,INPAR::STR::convnorm_rel,tol,normtype);
      break;
    }
    default:
      dserror("Unknown \"Tolerance Type\" defined!");
      break;
  } // switch case toltype

  // set tolerance
  qlist.set<double>("Tolerance",tol);

  // set norm type
  switch (normtype)
  {
    case INPAR::STR::norm_l2:
      qlist.set("Norm Type","Two Norm");
      break;
    case INPAR::STR::norm_l1:
      qlist.set("Norm Type","One Norm");
      break;
    case INPAR::STR::norm_inf:
      qlist.set("Norm Type","Max Norm");
      break;
    case INPAR::STR::norm_rms:
      dserror("This norm type is no longer supported! "
          "Consider to use the \"NOX::NLN::StatusTest::NormWRMS\" test instead!");
      break;
    case INPAR::STR::norm_vague:
    default:
      dserror("Unknown vector norm type!");
      break;
  } // switch case normtype

  // set scaling by the vector length
  if (isscaled)
    qlist.set("Scale Type","Scaled");

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::FullNewton::SetActiveSetParams(
    Teuchos::ParameterList& qlist,
    const enum NOX::NLN::StatusTest::QuantityType& qtype)
{
  CheckInit();

  qlist.set("Test Type","ActiveSet");
  // set the quantity type
  qlist.set("Quantity Type",
      NOX::NLN::StatusTest::QuantityType2String(qtype).c_str());

  return;
}
