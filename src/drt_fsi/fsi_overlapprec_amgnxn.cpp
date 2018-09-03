/*----------------------------------------------------------------------*/
/*!
\file fsi_overlapprec_amgnxn.cpp

\brief Special version of block matrix that includes the AMGnxn preconditioner

<pre>
Maintainer: Francesc Verdugo
            verduo@lnm.mw.tum.de
            89 - 289-15262
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef HAVE_MueLu

#include <Teuchos_PtrDecl.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "fsi_overlapprec_amgnxn.H"
#include "../solver/solver_amgnxn_preconditioner.H"
#include "../linalg/linalg_precond.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::OverlappingBlockMatrixAMGnxn::OverlappingBlockMatrixAMGnxn(
    const LINALG::MultiMapExtractor& maps, ADAPTER::FSIStructureWrapper& structure,
    ADAPTER::Fluid& fluid, ADAPTER::AleFsiWrapper& ale, bool structuresplit, std::string xml_file,
    FILE* err, std::string ProblemType)
    : OverlappingBlockMatrix(Teuchos::null, maps, structure, fluid, ale, structuresplit, 0, 1.0, 1,
          1.0, 0, 1.0, 0, 1.0, 0,
          err),  // Many of this inputs are irrelevant for the AMGnxn preconditioner
      P_(Teuchos::null),
      A_(Teuchos::null),
      xml_file_(xml_file),
      ProblemType_(ProblemType)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixAMGnxn::SetupPreconditioner()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::OverlappingBlockMatrixAMGnxn::SetupPreconditioner");

  // Free old matrix and preconditioner
  A_ = Teuchos::null;
  P_ = Teuchos::null;

  // Create own copy of the system matrix in order to allow reusing the preconditioner
  A_ = Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(Teuchos::rcp(this, false));
  if (A_ == Teuchos::null) dserror("We expect here a block matrix!");
  A_ = A_->Clone(LINALG::Copy);
  A_->Complete();

  // Determine number of blocks
  int NumBlocks = A_->Rows();
  if (A_->Rows() != A_->Cols())
    dserror("The AMGnxn preconditioner works only for block square matrices");

  // choose the right interface
  Teuchos::RCP<AMGnxnInterfaceFSI> myInterface = Teuchos::null;
  if (ProblemType_ == "FSI")
  {
    if (NumBlocks != 3) dserror("We expect a 3x3 block matrix");
    myInterface =
        Teuchos::rcp(new AMGnxnInterfaceFSI(structuresolver_, fluidsolver_, alesolver_, xml_file_));
  }
  else if (ProblemType_ == "LungFSI")
  {
    if (NumBlocks != 4) dserror("We expect a 4x4 block matrix");
    myInterface = Teuchos::rcp(
        new AMGnxnInterfaceLungFSI(structuresolver_, fluidsolver_, alesolver_, xml_file_));
  }
  myInterface->ParseXML();
  myInterface->Check();


  // Create the Operator
  if (myInterface->GetPreconditionerType() == "AMG(BlockSmoother)")
  {
    P_ = Teuchos::rcp(new LINALG::SOLVER::AMGnxn_Operator(A_, myInterface->GetNumPdes(),
        myInterface->GetNullSpacesDim(), myInterface->GetNullSpacesData(),
        myInterface->GetPreconditionerParams(), myInterface->GetSmoothersParams(),
        myInterface->GetSmoothersParams()));
  }
  else if (myInterface->GetPreconditionerType() == "BlockSmoother(X)")
  {
    P_ = Teuchos::rcp(new LINALG::SOLVER::BlockSmoother_Operator(A_, myInterface->GetNumPdes(),
        myInterface->GetNullSpacesDim(), myInterface->GetNullSpacesData(),
        myInterface->GetPreconditionerParams(), myInterface->GetSmoothersParams()));
  }
  else
    dserror("Unknown preconditioner type");

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::AMGnxnInterfaceFSI::AMGnxnInterfaceFSI(Teuchos::RCP<LINALG::Preconditioner> structuresolver,
    Teuchos::RCP<LINALG::Preconditioner> fluidsolver,
    Teuchos::RCP<LINALG::Preconditioner> alesolver, std::string amgnxn_xml)
    : structuresolver_(structuresolver),
      fluidsolver_(fluidsolver),
      alesolver_(alesolver),
      amgnxn_xml_(amgnxn_xml)
{
  Setup();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::AMGnxnInterfaceFSI::Setup()
{
  // The blocks are ordered as follows
  //  structure dof < fluid dof < ale dof

  // Recover null space info
  // structuresolver_->Params() is the parameter list of the structure solver
  // structuresolver_->Params().sublist("MueLu Parameters")should contain the null space info
  //  if we have selected a MueLu solver for the structure
  //  Idem for the ale and fuid
  //  CONCLUSSION: We have to choose a Muelu preconditioner for each single field

  int NumBlocks = 3;
  num_pdes_.resize(NumBlocks);
  null_spaces_dim_.resize(NumBlocks);
  null_spaces_data_.resize(NumBlocks);

  // Structure
  {
    Teuchos::ParameterList mllist;
    if ((structuresolver_->Params().isSublist("MueLu Parameters")))
      mllist = structuresolver_->Params().sublist("MueLu Parameters");
    else if ((structuresolver_->Params().isSublist("ML Parameters")))
      mllist = structuresolver_->Params().sublist("ML Parameters");
    else
      dserror("MueLu Parameters nor ML Parameters not found in structure solver");
    int block = 0;
    num_pdes_[block] = mllist.get<int>("PDE equations", -1);
    null_spaces_dim_[block] = mllist.get<int>("null space: dimension", -1);
    null_spaces_data_[block] =
        mllist.get<Teuchos::RCP<std::vector<double>>>("nullspace", Teuchos::null);
  }

  // Fluid
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  //  At the moment we are not enriching the null space of the fluid with
  //  the structure rotations at the interface
  // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
  {
    Teuchos::ParameterList mllist;
    if ((fluidsolver_->Params().isSublist("MueLu Parameters")))
      mllist = fluidsolver_->Params().sublist("MueLu Parameters");
    else if ((fluidsolver_->Params().isSublist("ML Parameters")))
      mllist = fluidsolver_->Params().sublist("ML Parameters");
    else
      dserror("MueLu Parameters nor ML Parameters not found in fluid solver");
    int block = 1;
    num_pdes_[block] = mllist.get<int>("PDE equations", -1);
    null_spaces_dim_[block] = mllist.get<int>("null space: dimension", -1);
    null_spaces_data_[block] =
        mllist.get<Teuchos::RCP<std::vector<double>>>("nullspace", Teuchos::null);
  }

  // Ale
  {
    Teuchos::ParameterList mllist;
    if ((alesolver_->Params().isSublist("MueLu Parameters")))
      mllist = alesolver_->Params().sublist("MueLu Parameters");
    else if ((alesolver_->Params().isSublist("ML Parameters")))
      mllist = alesolver_->Params().sublist("ML Parameters");
    else
      dserror("MueLu Parameters nor ML Parameters not found in ale solver");
    int block = 2;
    num_pdes_[block] = mllist.get<int>("PDE equations", -1);
    null_spaces_dim_[block] = mllist.get<int>("null space: dimension", -1);
    null_spaces_data_[block] =
        mllist.get<Teuchos::RCP<std::vector<double>>>("nullspace", Teuchos::null);
  }


  return;
}

void FSI::AMGnxnInterfaceFSI::ParseXML()
{
  // Decide whether to choose a default type or parse a xml file
  if (amgnxn_xml_ == "AMG(BGS)")
    Default_AMG_BGS(smoo_params_);
  else
  {
    // Parse the whole file
    Teuchos::updateParametersFromXmlFile(
        amgnxn_xml_, Teuchos::Ptr<Teuchos::ParameterList>(&smoo_params_));

    // Find the path to the main xml file and include it in the param list to further usage
    std::string::size_type pos = amgnxn_xml_.rfind('/');
    if (pos != std::string::npos)
    {
      std::string tmp = amgnxn_xml_.substr(0, pos + 1);
      smoo_params_.set<std::string>("main xml path", tmp);
    }
    else
      smoo_params_.set<std::string>("main xml path", "");
  }

  return;
}

void FSI::AMGnxnInterfaceFSI::Check()
{
  // Find preconditioner type and parameters
  std::string myprec = smoo_params_.get<std::string>("Preconditioner", "none");
  if (myprec == "none") dserror("Not found \"Preconditioner\" parameter in your xml file.");
  if (!smoo_params_.isSublist(myprec))
    dserror("Not found your preconditioner list in your xml file.");
  Teuchos::ParameterList& myprec_list = smoo_params_.sublist(myprec);
  prec_type_ = myprec_list.get<std::string>("type", "none");
  if (!myprec_list.isSublist("parameters"))
    dserror("Not found the parameters list for your preconditioner. Fix your xml file.");
  prec_params_ = myprec_list.sublist("parameters");

  return;
}


void FSI::AMGnxnInterfaceFSI::Default_AMG_BGS(Teuchos::ParameterList& params)
{
  std::string file_contents = "";
  file_contents +=
      "<ParameterList name=\"dummy\">                                                              "
      "                                      ";
  file_contents +=
      "  <Parameter name=\"Preconditioner\" type=\"string\"  value=\"myAMGBGS\"/>                  "
      "                                      ";
  file_contents +=
      "  <!-- Preconditioners  -->                                                                 "
      "                                      ";
  file_contents +=
      "  <ParameterList name=\"myAMGBGS\">                                                         "
      "                                      ";
  file_contents +=
      "  <Parameter name=\"type\" type=\"string\"  value=\"AMG(BlockSmoother)\"/>                  "
      "                                      ";
  file_contents +=
      "    <ParameterList name=\"parameters\">                                                     "
      "                                      ";
  file_contents +=
      "      <Parameter name=\"verbosity\"                        type=\"string\"  value=\"on\"/>  "
      "                                      ";
  file_contents +=
      "      <Parameter name=\"number of levels\"                 type=\"int\"     value=\"3\"/>   "
      "                                      ";
  file_contents +=
      "      <Parameter name=\"smoother: all but coarsest level\" type=\"string\"  "
      "value=\"myFineSmoother\"/>                            ";
  file_contents +=
      "      <Parameter name=\"smoother: coarsest level\"         type=\"string\"  "
      "value=\"myFineSmoother\"/>                            ";
  file_contents +=
      "      <Parameter name=\"muelu parameters for block 0\"       type=\"string\"  "
      "value=\"SA-AMG-3D\"/>                               ";
  file_contents +=
      "      <Parameter name=\"muelu parameters for block 1\"       type=\"string\"  "
      "value=\"PG-AMG-3D\"/>                               ";
  file_contents +=
      "      <Parameter name=\"muelu parameters for block 2\"       type=\"string\"  "
      "value=\"SA-AMG-3D\"/>                               ";
  file_contents +=
      "    </ParameterList>                                                                        "
      "                                      ";
  file_contents +=
      "  </ParameterList>                                                                          "
      "                                      ";
  file_contents +=
      "  <!-- Smoothers  -->                                                                       "
      "                                      ";
  file_contents +=
      "  <ParameterList name=\"myFineSmoother\">                                                   "
      "                                      ";
  file_contents +=
      "    <Parameter name=\"type\"   type=\"string\"  value=\"BGS\"/>                             "
      "                                      ";
  file_contents +=
      "    <ParameterList name=\"parameters\">                                                     "
      "                                      ";
  file_contents +=
      "      <Parameter name=\"blocks\"      type=\"string\"  value=\"(0),(1),(2)\"/>              "
      "                                      ";
  file_contents +=
      "      <Parameter name=\"smoothers\"   type=\"string\"  "
      "value=\"REUSE_MUELU_SMOOTHER,REUSE_MUELU_SMOOTHER,REUSE_MUELU_SMOOTHER\"/> ";
  file_contents +=
      "      <Parameter name=\"sweeps\"      type=\"int\"     value=\"1\"/>                        "
      "                                      ";
  file_contents +=
      "      <Parameter name=\"omega\"       type=\"double\"  value=\"1.0\"/>                      "
      "                                      ";
  file_contents +=
      "    </ParameterList>                                                                        "
      "                                      ";
  file_contents +=
      "  </ParameterList>                                                                          "
      "                                      ";
  file_contents +=
      "  <!-- Muelu  -->                                                                           "
      "                                      ";
  file_contents +=
      "  <ParameterList name=\"SA-AMG-3D\">                                                        "
      "                                      ";
  file_contents +=
      "    <ParameterList name=\"Factories\">                                                      "
      "                                      ";
  file_contents +=
      "      <ParameterList name=\"myCoalesceDropFact\">                                           "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"CoalesceDropFactory\"/>                       ";
  file_contents +=
      "        <Parameter name=\"lightweight wrap\"                type=\"bool\"   "
      "value=\"false\"/>                                     ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                                      ";
  file_contents +=
      "      <ParameterList name=\"myUncoupledAggregationFact\">                                   "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"UncoupledAggregationFactory\"/>               ";
  file_contents +=
      "        <Parameter name=\"Graph\"                           type=\"string\" "
      "value=\"myCoalesceDropFact\"/>                        ";
  file_contents +=
      "        <Parameter name=\"DofsPerNode\"                     type=\"string\" "
      "value=\"myCoalesceDropFact\"/>                        ";
  file_contents +=
      "        <Parameter name=\"aggregation: max selected neighbors\"    type=\"int\"    "
      "value=\"0\"/>                                  ";
  file_contents +=
      "        <Parameter name=\"aggregation: min agg size\"              type=\"int\"    "
      "value=\"12\"/>                                 ";
  file_contents +=
      "        <Parameter name=\"aggregation: max agg size\"              type=\"int\"    "
      "value=\"27\"/>                                 ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                                      ";
  file_contents +=
      "      <ParameterList name=\"myTentativePFactory\">                                          "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"TentativePFactory\"/>                         ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                                      ";
  file_contents +=
      "      <ParameterList name=\"myProlongatorFact\">                                            "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"SaPFactory\"/>                                ";
  file_contents +=
      "        <Parameter name=\"P\"                               type=\"string\" "
      "value=\"myTentativePFactory\"/>                       ";
  file_contents +=
      "        <Parameter name=\"Damping factor\"                  type=\"double\" "
      "value=\"1.333\"/>                                     ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                                      ";
  file_contents +=
      "      <ParameterList name=\"myRestrictorFact\">                                             "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"TransPFactory\"/>                             ";
  file_contents +=
      "        <Parameter name=\"P\"                               type=\"string\" "
      "value=\"myProlongatorFact\"/>                         ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                                      ";
  file_contents +=
      "      <ParameterList name=\"myRAPFact\">                                                    "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"RAPFactory\"/>                                ";
  file_contents +=
      "        <Parameter name=\"P\"                               type=\"string\" "
      "value=\"myProlongatorFact\"/>                         ";
  file_contents +=
      "        <Parameter name=\"R\"                               type=\"string\" "
      "value=\"myRestrictorFact\"/>                          ";
  file_contents +=
      "        <Parameter name=\"RepairMainDiagonal\"              type=\"bool\"   "
      "value=\"true\"/>                                      ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                                      ";
  file_contents +=
      "      <ParameterList name=\"myForwardGaussSeidel\">                                         "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"TrilinosSmoother\"/>                          ";
  file_contents +=
      "        <Parameter name=\"type\"                            type=\"string\"  "
      "value=\"RELAXATION\"/>                               ";
  file_contents +=
      "        <ParameterList name=\"ParameterList\">                                              "
      "                                      ";
  file_contents +=
      "          <Parameter name=\"relaxation: type\"              type=\"string\"  "
      "value=\"Gauss-Seidel\"/>                             ";
  file_contents +=
      "          <Parameter name=\"relaxation: backward mode\"     type=\"bool\"    "
      "value=\"false\"/>                                    ";
  file_contents +=
      "          <Parameter name=\"relaxation: sweeps\"            type=\"int\"     value=\"3\"/>  "
      "                                      ";
  file_contents +=
      "          <Parameter name=\"relaxation: damping factor\"    type=\"double\"  value=\"1\"/>  "
      "                                      ";
  file_contents +=
      "        </ParameterList>                                                                    "
      "                                      ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                                      ";
  file_contents +=
      "      <ParameterList name=\"mySolverFact\">                                                 "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"DirectSolver\"/>                              ";
  file_contents +=
      "        <Parameter name=\"type\"                            type=\"string\" value=\"Klu\"/> "
      "                                      ";
  file_contents +=
      "        <ParameterList name=\"ParameterList\">                                              "
      "                                      ";
  file_contents +=
      "          <Parameter name=\"Reindex\"                       type=\"bool\" value=\"true\"/>  "
      "                                      ";
  file_contents +=
      "        </ParameterList>                                                                    "
      "                                      ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                                      ";
  file_contents +=
      "    </ParameterList>                                                                        "
      "                                      ";
  file_contents +=
      "    <ParameterList name=\"Hierarchy\">                                                      "
      "                                      ";
  file_contents +=
      "      <Parameter name=\"numDesiredLevel\"          type=\"int\"      value=\"3\"/>          "
      "                                      ";
  file_contents +=
      "      <Parameter name=\"maxCoarseSize\"            type=\"int\"      value=\"1000\"/>       "
      "                                      ";
  file_contents +=
      "      <Parameter name=\"verbosity\"                type=\"string\"   value=\"High\"/>       "
      "                                      ";
  file_contents +=
      "      <ParameterList name=\"All\">                                                          "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"startLevel\"             type=\"int\"      value=\"0\"/>          "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"Smoother\"               type=\"string\"   "
      "value=\"myForwardGaussSeidel\"/>                             ";
  file_contents +=
      "        <Parameter name=\"CoarseSolver\"           type=\"string\"   "
      "value=\"mySolverFact\"/>                                     ";
  file_contents +=
      "        <Parameter name=\"Aggregates\"             type=\"string\"   "
      "value=\"myUncoupledAggregationFact\"/>                       ";
  file_contents +=
      "        <Parameter name=\"Graph\"                  type=\"string\"   "
      "value=\"myCoalesceDropFact\"/>                               ";
  file_contents +=
      "        <Parameter name=\"A\"                      type=\"string\"   value=\"myRAPFact\"/>  "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"P\"                      type=\"string\"   "
      "value=\"myProlongatorFact\"/>                                ";
  file_contents +=
      "        <Parameter name=\"R\"                      type=\"string\"   "
      "value=\"myRestrictorFact\"/>                                 ";
  file_contents +=
      "        <Parameter name=\"Nullspace\"              type=\"string\"   "
      "value=\"myTentativePFactory\"/>                              ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                                      ";
  file_contents +=
      "    </ParameterList>                                                                        "
      "                                      ";
  file_contents +=
      "  </ParameterList>                                                                          "
      "                                      ";
  file_contents +=
      "  <ParameterList name=\"PG-AMG-3D\">                                                        "
      "                                      ";
  file_contents +=
      "    <ParameterList name=\"Factories\">                                                      "
      "                                      ";
  file_contents +=
      "      <ParameterList name=\"myCoalesceDropFact\">                                           "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"CoalesceDropFactory\"/>                       ";
  file_contents +=
      "        <Parameter name=\"lightweight wrap\"                type=\"bool\"   "
      "value=\"false\"/>                                     ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                                      ";
  file_contents +=
      "      <ParameterList name=\"myUncoupledAggregationFact\">                                   "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"UncoupledAggregationFactory\"/>               ";
  file_contents +=
      "        <Parameter name=\"Graph\"                           type=\"string\" "
      "value=\"myCoalesceDropFact\"/>                        ";
  file_contents +=
      "        <Parameter name=\"DofsPerNode\"                     type=\"string\" "
      "value=\"myCoalesceDropFact\"/>                        ";
  file_contents +=
      "        <Parameter name=\"aggregation: max selected neighbors\"    type=\"int\"    "
      "value=\"0\"/>                                  ";
  file_contents +=
      "        <Parameter name=\"aggregation: min agg size\"              type=\"int\"    "
      "value=\"12\"/>                                 ";
  file_contents +=
      "        <Parameter name=\"aggregation: max agg size\"              type=\"int\"    "
      "value=\"27\"/>                                 ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                                      ";
  file_contents +=
      "      <ParameterList name=\"myTentativePFactory\">                                          "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"TentativePFactory\"/>                         ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                                      ";
  file_contents +=
      "      <ParameterList name=\"myProlongatorFact\">                                            "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"PgPFactory\"/>                                ";
  file_contents +=
      "        <Parameter name=\"P\"                               type=\"string\" "
      "value=\"myTentativePFactory\"/>                       ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                                      ";
  file_contents +=
      "      <ParameterList name=\"myRestrictorFact\">                                             "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"GenericRFactory\"/>                           ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                                      ";
  file_contents +=
      "      <ParameterList name=\"myRAPFact\">                                                    "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"RAPFactory\"/>                                ";
  file_contents +=
      "        <Parameter name=\"P\"                               type=\"string\" "
      "value=\"myProlongatorFact\"/>                         ";
  file_contents +=
      "        <Parameter name=\"R\"                               type=\"string\" "
      "value=\"myRestrictorFact\"/>                          ";
  file_contents +=
      "        <Parameter name=\"RepairMainDiagonal\"              type=\"bool\"   "
      "value=\"true\"/>                                      ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                                      ";
  file_contents +=
      "      <ParameterList name=\"myForwardGaussSeidel\">                                         "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"TrilinosSmoother\"/>                          ";
  file_contents +=
      "        <Parameter name=\"type\"                            type=\"string\"  "
      "value=\"RELAXATION\"/>                               ";
  file_contents +=
      "        <ParameterList name=\"ParameterList\">                                              "
      "                                      ";
  file_contents +=
      "          <Parameter name=\"relaxation: type\"              type=\"string\"  "
      "value=\"Gauss-Seidel\"/>                             ";
  file_contents +=
      "          <Parameter name=\"relaxation: backward mode\"     type=\"bool\"    "
      "value=\"false\"/>                                    ";
  file_contents +=
      "          <Parameter name=\"relaxation: sweeps\"            type=\"int\"     value=\"3\"/>  "
      "                                      ";
  file_contents +=
      "          <Parameter name=\"relaxation: damping factor\"    type=\"double\"  value=\"1\"/>  "
      "                                      ";
  file_contents +=
      "        </ParameterList>                                                                    "
      "                                      ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                                      ";
  file_contents +=
      "      <ParameterList name=\"mySolverFact\">                                                 "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"DirectSolver\"/>                              ";
  file_contents +=
      "        <Parameter name=\"type\"                            type=\"string\" value=\"Klu\"/> "
      "                                      ";
  file_contents +=
      "        <ParameterList name=\"ParameterList\">                                              "
      "                                      ";
  file_contents +=
      "          <Parameter name=\"Reindex\"                       type=\"bool\" value=\"true\"/>  "
      "                                      ";
  file_contents +=
      "        </ParameterList>                                                                    "
      "                                      ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                                      ";
  file_contents +=
      "    </ParameterList>                                                                        "
      "                                      ";
  file_contents +=
      "    <ParameterList name=\"Hierarchy\">                                                      "
      "                                      ";
  file_contents +=
      "      <Parameter name=\"max levels\"               type=\"int\"      value=\"3\"/>          "
      "                                      ";
  file_contents +=
      "      <Parameter name=\"coarse: max size\"         type=\"int\"      value=\"1000\"/>       "
      "                                      ";
  file_contents +=
      "      <Parameter name=\"verbosity\"                type=\"string\"   value=\"High\"/>       "
      "                                      ";
  file_contents +=
      "      <ParameterList name=\"All\">                                                          "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"startLevel\"             type=\"int\"      value=\"0\"/>          "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"Smoother\"               type=\"string\"   "
      "value=\"myForwardGaussSeidel\"/>                             ";
  file_contents +=
      "        <Parameter name=\"CoarseSolver\"           type=\"string\"   "
      "value=\"mySolverFact\"/>                                     ";
  file_contents +=
      "        <Parameter name=\"Aggregates\"             type=\"string\"   "
      "value=\"myUncoupledAggregationFact\"/>                       ";
  file_contents +=
      "        <Parameter name=\"Graph\"                  type=\"string\"   "
      "value=\"myCoalesceDropFact\"/>                               ";
  file_contents +=
      "        <Parameter name=\"A\"                      type=\"string\"   value=\"myRAPFact\"/>  "
      "                                      ";
  file_contents +=
      "        <Parameter name=\"P\"                      type=\"string\"   "
      "value=\"myProlongatorFact\"/>                                ";
  file_contents +=
      "        <Parameter name=\"R\"                      type=\"string\"   "
      "value=\"myRestrictorFact\"/>                                 ";
  file_contents +=
      "        <Parameter name=\"Nullspace\"              type=\"string\"   "
      "value=\"myTentativePFactory\"/>                              ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                                      ";
  file_contents +=
      "    </ParameterList>                                                                        "
      "                                      ";
  file_contents +=
      "  </ParameterList>                                                                          "
      "                                      ";
  file_contents +=
      "</ParameterList>                                                                            "
      "                                      ";

  Teuchos::updateParametersFromXmlString(
      file_contents, Teuchos::Ptr<Teuchos::ParameterList>(&params));

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::AMGnxnInterfaceLungFSI::AMGnxnInterfaceLungFSI(
    Teuchos::RCP<LINALG::Preconditioner> structuresolver,
    Teuchos::RCP<LINALG::Preconditioner> fluidsolver,
    Teuchos::RCP<LINALG::Preconditioner> alesolver, std::string amgnxn_xml)
    : AMGnxnInterfaceFSI(structuresolver, fluidsolver, alesolver, amgnxn_xml)
{
  // The blocks are ordered as follows
  //  structure dof < fluid dof < ale dof < constrain dof

  // Here we append null space info for the constrain block
  // TODO Generate null space info for the second block
  num_pdes_.push_back(-1);
  null_spaces_dim_.push_back(-1);
  null_spaces_data_.push_back(Teuchos::null);

  // Some safety checks
  if (amgnxn_xml_ == "AMG(BGS)")
    dserror("A default AMG(BGS) preconditioner is not implemented for lung fsi simulations");
  if (prec_type_ == "AMG(BlockSmoother)")
    dserror("A AMG(BlockSmoother) preconditioner is not implemented yet for lung fsi simulations");
}


void FSI::AMGnxnInterfaceLungFSI::ParseXML()
{
  // Decide whether to choose a default type or parse a xml file
  if (amgnxn_xml_ == "SIMPLE(ILU,KLU)")
    Default_SCHUR_ILU_KLU(smoo_params_);
  else
  {
    // Parse the whole file
    Teuchos::updateParametersFromXmlFile(
        amgnxn_xml_, Teuchos::Ptr<Teuchos::ParameterList>(&smoo_params_));

    // Find the path to the main xml file and include it in the param list to further usage
    std::string::size_type pos = amgnxn_xml_.rfind('/');
    if (pos != std::string::npos)
    {
      std::string tmp = amgnxn_xml_.substr(0, pos + 1);
      smoo_params_.set<std::string>("main xml path", tmp);
    }
    else
      smoo_params_.set<std::string>("main xml path", "");
  }

  return;
}

void FSI::AMGnxnInterfaceLungFSI::Default_SCHUR_ILU_KLU(Teuchos::ParameterList& params)
{
  std::string file_contents = "";
  file_contents +=
      "<ParameterList name=\"dummy\">                                                              "
      "          ";
  file_contents +=
      "  <Parameter name=\"Preconditioner\" type=\"string\"  value=\"myBGS_SIMPLE\"/>              "
      "          ";
  file_contents +=
      "  <ParameterList name=\"myBGS_SIMPLE\">                                                     "
      "          ";
  file_contents +=
      "  <Parameter name=\"type\" type=\"string\"  value=\"BlockSmoother(X)\"/>                    "
      "          ";
  file_contents +=
      "    <ParameterList name=\"parameters\">                                                     "
      "          ";
  file_contents +=
      "      <Parameter name=\"smoother\"         type=\"string\"  value=\"mySIMPLE\"/>            "
      "          ";
  file_contents +=
      "      <Parameter name=\"verbosity\"        type=\"string\"  value=\"on\"/>                  "
      "          ";
  file_contents +=
      "    </ParameterList>                                                                        "
      "          ";
  file_contents +=
      "  </ParameterList>                                                                          "
      "          ";
  file_contents +=
      "  <ParameterList name=\"mySIMPLE\">                                                         "
      "          ";
  file_contents +=
      "    <Parameter name=\"type\"   type=\"string\"  value=\"SIMPLE\"/>                          "
      "          ";
  file_contents +=
      "    <ParameterList name=\"parameters\">                                                     "
      "          ";
  file_contents +=
      "      <Parameter name=\"predictor block\"     type=\"string\"  value=\"(0,1,2)\"/>          "
      "          ";
  file_contents +=
      "      <Parameter name=\"predictor smoother\"  type=\"string\"  value=\"myBGS\"/>            "
      "          ";
  file_contents +=
      "      <Parameter name=\"predictor inverse\"   type=\"string\"  value=\"row sums diagonal "
      "blocks\"/>   ";
  file_contents +=
      "      <Parameter name=\"schur block\"         type=\"string\"  value=\"(3)\"/>              "
      "          ";
  file_contents +=
      "      <Parameter name=\"schur smoother\"      type=\"string\"  value=\"DIRECT_SOLVER\"/>    "
      "          ";
  file_contents +=
      "      <Parameter name=\"correction\"          type=\"string\"  value=\"approximated "
      "inverse\"/>       ";
  file_contents +=
      "      <Parameter name=\"sweeps\"              type=\"int\"     value=\"2\"/>                "
      "          ";
  file_contents +=
      "      <Parameter name=\"alpha\"               type=\"double\"  value=\"0.8\"/>              "
      "          ";
  file_contents +=
      "    </ParameterList>                                                                        "
      "          ";
  file_contents +=
      "  </ParameterList>                                                                          "
      "          ";
  file_contents +=
      "  <ParameterList name=\"myBGS\">                                                            "
      "          ";
  file_contents +=
      "    <Parameter name=\"type\"   type=\"string\"  value=\"BGS\"/>                             "
      "          ";
  file_contents +=
      "    <ParameterList name=\"parameters\">                                                     "
      "          ";
  file_contents +=
      "      <Parameter name=\"blocks\"      type=\"string\"  value=\"(0),(1),(2)\"/>              "
      "          ";
  file_contents +=
      "      <Parameter name=\"smoothers\"   type=\"string\"  value=\"myILU,myILU,myILU\"/>        "
      "          ";
  file_contents +=
      "      <Parameter name=\"sweeps\"      type=\"int\"     value=\"2\"/>                        "
      "          ";
  file_contents +=
      "      <Parameter name=\"omega\"       type=\"double\"  value=\"1.0\"/>                      "
      "          ";
  file_contents +=
      "    </ParameterList>                                                                        "
      "          ";
  file_contents +=
      "  </ParameterList>                                                                          "
      "          ";
  file_contents +=
      "  <ParameterList name=\"myILU\">                                                            "
      "          ";
  file_contents +=
      "    <Parameter name=\"type\"   type=\"string\"  value=\"IFPACK\"/>                          "
      "          ";
  file_contents +=
      "    <ParameterList name=\"parameters\">                                                     "
      "          ";
  file_contents +=
      "      <Parameter name=\"type\"      type=\"string\"  value=\"ILU\"/>                        "
      "          ";
  file_contents +=
      "      <ParameterList name=\"ParameterList\">                                                "
      "          ";
  file_contents +=
      "        <Parameter name=\"fact: level-of-fill-in\"      type=\"int\"  value=\"0\"/>         "
      "          ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "          ";
  file_contents +=
      "    </ParameterList>                                                                        "
      "          ";
  file_contents +=
      "  </ParameterList>                                                                          "
      "          ";
  file_contents +=
      "</ParameterList>                                                                            "
      "          ";

  Teuchos::updateParametersFromXmlString(
      file_contents, Teuchos::Ptr<Teuchos::ParameterList>(&params));
  return;
}


#endif  // HAVE_MueLu
