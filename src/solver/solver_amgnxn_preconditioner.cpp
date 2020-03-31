/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15262
Created on: Feb 27, 2014
*----------------------------------------------------------------------*/

#ifdef HAVE_MueLu

#include <iostream>

#include <Teuchos_PtrDecl.hpp>
#include <Epetra_Time.h>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <MueLu_MLParameterListInterpreter_decl.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include "EpetraExt_RowMatrixOut.h"
#include "../drt_lib/drt_dserror.H"
#include "solver_amgnxn_preconditioner.H"
#include "solver_amgnxn_vcycle.H"

#include <Teuchos_TimeMonitor.hpp>

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

LINALG::SOLVER::AMGnxn_Preconditioner::AMGnxn_Preconditioner(
    FILE* outfile, Teuchos::ParameterList& params)
    : LINALG::SOLVER::PreconditionerType(outfile), params_(params)
{
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Epetra_Operator* LINALG::SOLVER::AMGnxn_Preconditioner::PrecOperator() const { return &*P_; }

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<Epetra_Operator> LINALG::SOLVER::AMGnxn_Preconditioner::PrecOperatorRCP() const
{
  return P_;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_Preconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  // Setup underlying linear system
  SetupLinearProblem(matrix, x, b);

  // Decide if the setup has to be done
  if (!create) return;

  // Check whether this is a block sparse matrix
  BlockSparseMatrixBase* A_bl = dynamic_cast<BlockSparseMatrixBase*>(matrix);
  if (A_bl == NULL)
    dserror("The AMGnxn preconditioner works only for BlockSparseMatrixBase or derived classes");

  // Do all the setup
  Setup(Teuchos::rcp(A_bl, false));

  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_Preconditioner::Setup(Teuchos::RCP<BlockSparseMatrixBase> A)
{
  TEUCHOS_FUNC_TIME_MONITOR("LINALG::SOLVER::AMGnxn_Preconditioner::Setup");

  Epetra_Time timer(A->Comm());
  timer.ResetStartTime();

  // Free old matrix and preconditioner
  A_ = Teuchos::null;
  P_ = Teuchos::null;

  // Create own copy of the system matrix in order to allow reusing the preconditioner
  A_ = A;
  A_ = A_->Clone(Copy);
  A_->Complete();

  // Determine number of blocks
  int NumBlocks = A_->Rows();
  if (A_->Rows() != A_->Cols())
    dserror("The AMGnxn preconditioner works only for block square matrices");

  // Pick-up the input parameters
  AMGnxn_Interface myInterface(params_, NumBlocks);

  // Create the Operator
  if (myInterface.GetPreconditionerType() == "AMG(BlockSmoother)")
  {
    P_ = Teuchos::rcp(
        new AMGnxn_Operator(A_, myInterface.GetNumPdes(), myInterface.GetNullSpacesDim(),
            myInterface.GetNullSpacesData(), myInterface.GetPreconditionerParams(),
            myInterface.GetSmoothersParams(), myInterface.GetSmoothersParams()));
  }
  else if (myInterface.GetPreconditionerType() == "BlockSmoother(X)")
  {
    P_ = Teuchos::rcp(new BlockSmoother_Operator(A_, myInterface.GetNumPdes(),
        myInterface.GetNullSpacesDim(), myInterface.GetNullSpacesData(),
        myInterface.GetPreconditionerParams(), myInterface.GetSmoothersParams()));
  }
  else if (myInterface.GetPreconditionerType() == "Merge_and_Ifpack")
  {
    P_ = Teuchos::rcp(new Merged_Operator(
        A_, myInterface.GetPreconditionerParams(), myInterface.GetSmoothersParams()));
  }
  else
    dserror("Unknown preconditioner type: %s", myInterface.GetPreconditionerType().c_str());

  double elaptime = timer.ElapsedTime();
  if (myInterface.GetPreconditionerParams().get<std::string>("verbosity", "off") == "on" and
      A->Comm().MyPID() == 0)
    std::cout << "       Calling LINALG::SOLVER::AMGnxn_Preconditioner::Setup takes "
              << std::setw(16) << std::setprecision(6) << elaptime << " s" << std::endl;

  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

LINALG::SOLVER::AMGnxn_Interface::AMGnxn_Interface(Teuchos::ParameterList& params, int NumBlocks)
{
  // Expected parameters in params
  //<ParameterList name="params">
  //
  //  <ParameterList name="AMGnxn Parameters">
  //    <Parameter name="AMGNXN_TYPE"          type="string"  value="AMG(BGS)"/>
  //    <!-- or -->
  //    <Parameter name="AMGNXN_TYPE"          type="string"  value="XML"/>
  //    <Parameter name="AMGNXN_XML_FILE"      type="string"  value=" ... .xml"/>
  //  </ParameterList>
  //
  //  <ParameterList name="Inverse1">
  //    <ParameterList name="MueLu Parameters">
  //    <Parameter name="PDE equations"          type="int"     value="..."/>
  //    <Parameter name="null space: dimension"  type="int"     value="..."/>
  //    <Parameter name="nullspace"              type="..."     value="..."/>
  //    </ParameterList>
  //  </ParameterList>
  //
  //</ParameterList>
  //
  //
  //
  // The xml file given in AMGNXN_XML_FILE should have the following format
  //<ParameterList name="dummy list which wraps everything">
  //
  //  <!-- Here we select which preconditioner we are going to use -->
  //  <Parameter name="Preconditioner"    type="string"  value="myPreconditioner"/>
  //
  //  <!-- Here we define our preconditioner -->
  //  <ParameterList name="myPreconditioner">
  //    <Parameter name="type"   type="string"  value="AMG(BlockSmoother)"/>
  //    <!-- or -- >
  //    <Parameter name="type"   type="string"  value="BlockSmoother(X)"/>
  //    <ParameterList name="parameters">
  //
  //     <!-- Fill this list with the parameters expected by your preconditioner -->
  //
  //    </ParameterList>
  //  </ParameterList>
  //
  //
  //  <!-- Here we put as many list as you need to define your smoothers -->
  //  <ParameterList name="mySmoother">
  //
  //     <!-- Fill this list with the parameters expected by your smoother -->
  //
  //  </ParameterList>
  //
  //</ParameterList>

  if (!params.isSublist("AMGnxn Parameters")) dserror("AMGnxn Parameters not found!");
  Teuchos::ParameterList& amglist = params.sublist("AMGnxn Parameters");


  // Decide whether to choose a default type or parse a xml file
  std::string amgnxn_type = amglist.get<std::string>("AMGNXN_TYPE", "AMG(BGS)");
  if (amgnxn_type == "TSI: AMG(BGS)")
    Params_TSI_AMG_BGS(smoo_params_);
  else if (amgnxn_type == "XML")
  {
    // Parse the whole file
    std::string amgnxn_xml = amglist.get<std::string>("AMGNXN_XML_FILE", "none");
    if (amgnxn_xml == "none") dserror("The input parameter AMGNXN_XML_FILE is empty.");
    if (not(amgnxn_xml == "none"))
    {
      Teuchos::updateParametersFromXmlFile(
          amgnxn_xml, Teuchos::Ptr<Teuchos::ParameterList>(&smoo_params_));

      // Find the path to the main xml file and include it in the param list to further usage
      std::string::size_type pos = amgnxn_xml.rfind('/');
      if (pos != std::string::npos)
      {
        std::string tmp = amgnxn_xml.substr(0, pos + 1);
        smoo_params_.set<std::string>("main xml path", tmp);
      }
      else
        smoo_params_.set<std::string>("main xml path", "");
    }
  }
  else
    dserror("\"%s\" is an invalid value for \"AMGNXN_TYPE\". Fix your .dat", amgnxn_type.c_str());



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

  // Find null spaces and relatives
  std::string Inverse_str = "Inverse";
  xml_files_.resize(NumBlocks);
  num_pdes_.resize(NumBlocks);
  null_spaces_dim_.resize(NumBlocks);
  null_spaces_data_.resize(NumBlocks);
  for (int block = 0; block < NumBlocks; block++)
  {
    if (!params.isSublist(Inverse_str + ConvertInt(block + 1)))
      dserror("Not found inverse list for block %d", block + 1);
    Teuchos::ParameterList& inverse_list = params.sublist(Inverse_str + ConvertInt(block + 1));

    if (!inverse_list.isSublist("MueLu Parameters")) dserror("MueLu Parameters not found");
    Teuchos::ParameterList& mllist = inverse_list.sublist("MueLu Parameters");

    xml_files_[block] = mllist.get<std::string>("xml file", "none");
    num_pdes_[block] = mllist.get<int>("PDE equations", -1);
    null_spaces_dim_[block] = mllist.get<int>("null space: dimension", -1);
    null_spaces_data_[block] =
        mllist.get<Teuchos::RCP<std::vector<double>>>("nullspace", Teuchos::null);

    // Some checks
    if (num_pdes_[block] < 1 or null_spaces_dim_[block] < 1)
      dserror("Error: PDE equations or null space dimension wrong.");
    if (null_spaces_data_[block] == Teuchos::null) dserror("Error: null space data is empty");
  }
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void LINALG::SOLVER::AMGnxn_Interface::Params_TSI_AMG_BGS(Teuchos::ParameterList& params)
{
  // This is a pre-defined AMG(BGS) preconditioner
  // TODO Add the possibility to the user of tuning some of the
  // more relevant parameters using the existing input params in the .dat file

  // This is the xml file which defines the preconditioner

  std::string file_contents = "";

  file_contents +=
      "<ParameterList name=\"dummy\">                                                              "
      "                       ";
  file_contents +=
      "  <Parameter name=\"Preconditioner\" type=\"string\"  value=\"myAMGBGS\"/>                  "
      "                       ";
  file_contents +=
      "  <ParameterList name=\"myAMGBGS\">                                                         "
      "                       ";
  file_contents +=
      "  <Parameter name=\"type\" type=\"string\"  value=\"AMG(BlockSmoother)\"/>                  "
      "                       ";
  file_contents +=
      "    <ParameterList name=\"parameters\">                                                     "
      "                       ";
  file_contents +=
      "      <Parameter name=\"verbosity\"                        type=\"string\"  value=\"on\"/>  "
      "                       ";
  file_contents +=
      "      <Parameter name=\"number of levels\"                 type=\"int\"     value=\"3\"/>   "
      "                       ";
  file_contents +=
      "      <Parameter name=\"smoother: all but coarsest level\" type=\"string\"  "
      "value=\"myFineSmoother\"/>             ";
  file_contents +=
      "      <Parameter name=\"smoother: coarsest level\"         type=\"string\"  "
      "value=\"myFineSmoother\"/>             ";
  file_contents +=
      "      <Parameter name=\"muelu parameters for block 0\"       type=\"string\"  "
      "value=\"SA-AMG-3D\"/>                ";
  file_contents +=
      "      <Parameter name=\"muelu parameters for block 1\"       type=\"string\"  "
      "value=\"SA-AMG-1D\"/>                ";
  file_contents +=
      "    </ParameterList>                                                                        "
      "                       ";
  file_contents +=
      "  </ParameterList>                                                                          "
      "                       ";
  file_contents +=
      "  <ParameterList name=\"myFineSmoother\">                                                   "
      "                       ";
  file_contents +=
      "    <Parameter name=\"type\"   type=\"string\"  value=\"BGS\"/>                             "
      "                       ";
  file_contents +=
      "    <ParameterList name=\"parameters\">                                                     "
      "                       ";
  file_contents +=
      "      <Parameter name=\"blocks\"      type=\"string\"  value=\"(0),(1)\"/>                  "
      "                       ";
  file_contents +=
      "      <Parameter name=\"smoothers\"   type=\"string\"  "
      "value=\"REUSE_MUELU_SMOOTHER,REUSE_MUELU_SMOOTHER\"/>       ";
  file_contents +=
      "      <Parameter name=\"sweeps\"      type=\"int\"     value=\"1\"/>                        "
      "                       ";
  file_contents +=
      "      <Parameter name=\"omega\"       type=\"double\"  value=\"1.0\"/>                      "
      "                       ";
  file_contents +=
      "    </ParameterList>                                                                        "
      "                       ";
  file_contents +=
      "  </ParameterList>                                                                          "
      "                       ";
  file_contents +=
      "  <ParameterList name=\"SA-AMG-3D\">                                                        "
      "                       ";
  file_contents +=
      "    <ParameterList name=\"Matrix\">                                                         "
      "                       ";
  file_contents +=
      "    <Parameter name=\"PDE equations\"                         type=\"int\" "
      "value=\"3\"/>        ";
  file_contents +=
      "    </ParameterList>                                                         "
      "                       ";
  file_contents +=
      "    <ParameterList name=\"Factories\">                                                      "
      "                       ";
  file_contents +=
      "      <ParameterList name=\"myCoalesceDropFact\">                                           "
      "                       ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"CoalesceDropFactory\"/>        ";
  file_contents +=
      "        <Parameter name=\"lightweight wrap\"                type=\"bool\"   "
      "value=\"false\"/>                      ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                       ";
  file_contents +=
      "      <ParameterList name=\"myUncoupledAggregationFact\">                                   "
      "                       ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"UncoupledAggregationFactory\"/>";
  file_contents +=
      "        <Parameter name=\"Graph\"                           type=\"string\" "
      "value=\"myCoalesceDropFact\"/>         ";
  file_contents +=
      "        <Parameter name=\"DofsPerNode\"                     type=\"string\" "
      "value=\"myCoalesceDropFact\"/>         ";
  file_contents +=
      "        <Parameter name=\"aggregation: max selected neighbors\"    type=\"int\"    "
      "value=\"0\"/>                   ";
  file_contents +=
      "        <Parameter name=\"aggregation: min agg size\"              type=\"int\"    "
      "value=\"12\"/>                  ";
  file_contents +=
      "        <Parameter name=\"aggregation: max agg size\"              type=\"int\"    "
      "value=\"27\"/>                  ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                       ";
  file_contents +=
      "      <ParameterList name=\"myTentativePFactory\">                                          "
      "                       ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"TentativePFactory\"/>          ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                       ";
  file_contents +=
      "      <ParameterList name=\"myProlongatorFact\">                                            "
      "                       ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"SaPFactory\"/>                 ";
  file_contents +=
      "        <Parameter name=\"P\"                               type=\"string\" "
      "value=\"myTentativePFactory\"/>        ";
  file_contents +=
      "        <Parameter name=\"Damping factor\"                  type=\"double\" "
      "value=\"1.333\"/>                      ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                       ";
  file_contents +=
      "      <ParameterList name=\"myRestrictorFact\">                                             "
      "                       ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"TransPFactory\"/>              ";
  file_contents +=
      "        <Parameter name=\"P\"                               type=\"string\" "
      "value=\"myProlongatorFact\"/>          ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                       ";
  file_contents +=
      "      <ParameterList name=\"myRAPFact\">                                                    "
      "                       ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"RAPFactory\"/>                 ";
  file_contents +=
      "        <Parameter name=\"P\"                               type=\"string\" "
      "value=\"myProlongatorFact\"/>          ";
  file_contents +=
      "        <Parameter name=\"R\"                               type=\"string\" "
      "value=\"myRestrictorFact\"/>           ";
  file_contents +=
      "        <Parameter name=\"RepairMainDiagonal\"              type=\"bool\"   "
      "value=\"true\"/>                       ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                       ";
  file_contents +=
      "      <ParameterList name=\"myForwardGaussSeidel\">                                         "
      "                       ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"TrilinosSmoother\"/>           ";
  file_contents +=
      "        <Parameter name=\"type\"                            type=\"string\"  "
      "value=\"RELAXATION\"/>                ";
  file_contents +=
      "        <ParameterList name=\"ParameterList\">                                              "
      "                       ";
  file_contents +=
      "          <Parameter name=\"relaxation: type\"              type=\"string\"  "
      "value=\"Gauss-Seidel\"/>              ";
  file_contents +=
      "          <Parameter name=\"relaxation: backward mode\"     type=\"bool\"    "
      "value=\"false\"/>                     ";
  file_contents +=
      "          <Parameter name=\"relaxation: sweeps\"            type=\"int\"     value=\"4\"/>  "
      "                       ";
  file_contents +=
      "          <Parameter name=\"relaxation: damping factor\"    type=\"double\"  "
      "value=\"0.79\"/>                      ";
  file_contents +=
      "        </ParameterList>                                                                    "
      "                       ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                       ";
  file_contents +=
      "      <ParameterList name=\"mySolverFact\">                                                 "
      "                       ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"DirectSolver\"/>               ";
  file_contents +=
      "        <Parameter name=\"type\"                            type=\"string\" value=\"Klu\"/> "
      "                       ";
  file_contents +=
      "        <ParameterList name=\"ParameterList\">                                              "
      "                       ";
  file_contents +=
      "          <Parameter name=\"Reindex\"                       type=\"bool\" value=\"true\"/>  "
      "                       ";
  file_contents +=
      "        </ParameterList>                                                                    "
      "                       ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                       ";
  file_contents +=
      "    </ParameterList>                                                                        "
      "                       ";
  file_contents +=
      "    <ParameterList name=\"Hierarchy\">                                                      "
      "                       ";
  file_contents +=
      "      <Parameter name=\"max levels\"               type=\"int\"      value=\"3\"/>          "
      "                       ";
  file_contents +=
      "      <Parameter name=\"coarse: max size\"         type=\"int\"      value=\"1\"/>          "
      "                       ";
  file_contents +=
      "      <Parameter name=\"verbosity\"                type=\"string\"   value=\"Low\"/>        "
      "                       ";
  file_contents +=
      "      <ParameterList name=\"All\">                                                          "
      "                       ";
  file_contents +=
      "        <Parameter name=\"startLevel\"             type=\"int\"      value=\"0\"/>          "
      "                       ";
  file_contents +=
      "        <Parameter name=\"Smoother\"               type=\"string\"   "
      "value=\"myForwardGaussSeidel\"/>              ";
  file_contents +=
      "        <Parameter name=\"CoarseSolver\"           type=\"string\"   "
      "value=\"mySolverFact\"/>                      ";
  file_contents +=
      "        <Parameter name=\"Aggregates\"             type=\"string\"   "
      "value=\"myUncoupledAggregationFact\"/>        ";
  file_contents +=
      "        <Parameter name=\"Graph\"                  type=\"string\"   "
      "value=\"myCoalesceDropFact\"/>                ";
  file_contents +=
      "        <Parameter name=\"A\"                      type=\"string\"   value=\"myRAPFact\"/>  "
      "                       ";
  file_contents +=
      "        <Parameter name=\"P\"                      type=\"string\"   "
      "value=\"myProlongatorFact\"/>                 ";
  file_contents +=
      "        <Parameter name=\"R\"                      type=\"string\"   "
      "value=\"myRestrictorFact\"/>                  ";
  file_contents +=
      "        <Parameter name=\"Nullspace\"              type=\"string\"   "
      "value=\"myTentativePFactory\"/>               ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                       ";
  file_contents +=
      "    </ParameterList>                                                                        "
      "                       ";
  file_contents +=
      "  </ParameterList>                                                                          "
      "                       ";
  file_contents +=
      "  <ParameterList name=\"SA-AMG-1D\">                                                        "
      "                       ";
  file_contents +=
      "    <ParameterList name=\"Matrix\">                                                         "
      "                       ";
  file_contents +=
      "    <Parameter name=\"PDE equations\"                         type=\"int\" "
      "value=\"1\"/>        ";
  file_contents +=
      "    </ParameterList>                                                         "
      "                       ";
  file_contents +=
      "    <ParameterList name=\"Factories\">                                                      "
      "                       ";
  file_contents +=
      "      <ParameterList name=\"myCoalesceDropFact\">                                           "
      "                       ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"CoalesceDropFactory\"/>        ";
  file_contents +=
      "        <Parameter name=\"lightweight wrap\"                type=\"bool\"   "
      "value=\"false\"/>                      ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                       ";
  file_contents +=
      "      <ParameterList name=\"myUncoupledAggregationFact\">                                   "
      "                       ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"UncoupledAggregationFactory\"/>";
  file_contents +=
      "        <Parameter name=\"Graph\"                           type=\"string\" "
      "value=\"myCoalesceDropFact\"/>         ";
  file_contents +=
      "        <Parameter name=\"DofsPerNode\"                     type=\"string\" "
      "value=\"myCoalesceDropFact\"/>         ";
  file_contents +=
      "        <Parameter name=\"aggregation: max selected neighbors\"    type=\"int\"    "
      "value=\"0\"/>                   ";
  file_contents +=
      "        <Parameter name=\"aggregation: min agg size\"              type=\"int\"    "
      "value=\"12\"/>                  ";
  file_contents +=
      "        <Parameter name=\"aggregation: max agg size\"              type=\"int\"    "
      "value=\"27\"/>                  ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                       ";
  file_contents +=
      "      <ParameterList name=\"myTentativePFactory\">                                          "
      "                       ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"TentativePFactory\"/>          ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                       ";
  file_contents +=
      "      <ParameterList name=\"myProlongatorFact\">                                            "
      "                       ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"SaPFactory\"/>                 ";
  file_contents +=
      "        <Parameter name=\"P\"                               type=\"string\" "
      "value=\"myTentativePFactory\"/>        ";
  file_contents +=
      "        <Parameter name=\"Damping factor\"                  type=\"double\" "
      "value=\"1.333\"/>                      ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                       ";
  file_contents +=
      "      <ParameterList name=\"myRestrictorFact\">                                             "
      "                       ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"TransPFactory\"/>              ";
  file_contents +=
      "        <Parameter name=\"P\"                               type=\"string\" "
      "value=\"myProlongatorFact\"/>          ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                       ";
  file_contents +=
      "      <ParameterList name=\"myRAPFact\">                                                    "
      "                       ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"RAPFactory\"/>                 ";
  file_contents +=
      "        <Parameter name=\"P\"                               type=\"string\" "
      "value=\"myProlongatorFact\"/>          ";
  file_contents +=
      "        <Parameter name=\"R\"                               type=\"string\" "
      "value=\"myRestrictorFact\"/>           ";
  file_contents +=
      "        <Parameter name=\"RepairMainDiagonal\"              type=\"bool\"   "
      "value=\"true\"/>                       ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                       ";
  file_contents +=
      "      <ParameterList name=\"myForwardGaussSeidel\">                                         "
      "                       ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"TrilinosSmoother\"/>           ";
  file_contents +=
      "        <Parameter name=\"type\"                            type=\"string\"  "
      "value=\"RELAXATION\"/>                ";
  file_contents +=
      "        <ParameterList name=\"ParameterList\">                                              "
      "                       ";
  file_contents +=
      "          <Parameter name=\"relaxation: type\"              type=\"string\"  "
      "value=\"Gauss-Seidel\"/>              ";
  file_contents +=
      "          <Parameter name=\"relaxation: backward mode\"     type=\"bool\"    "
      "value=\"false\"/>                     ";
  file_contents +=
      "          <Parameter name=\"relaxation: sweeps\"            type=\"int\"     value=\"4\"/>  "
      "                       ";
  file_contents +=
      "          <Parameter name=\"relaxation: damping factor\"    type=\"double\"  "
      "value=\"0.79\"/>                      ";
  file_contents +=
      "        </ParameterList>                                                                    "
      "                       ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                       ";
  file_contents +=
      "      <ParameterList name=\"mySolverFact\">                                                 "
      "                       ";
  file_contents +=
      "        <Parameter name=\"factory\"                         type=\"string\" "
      "value=\"DirectSolver\"/>               ";
  file_contents +=
      "        <Parameter name=\"type\"                            type=\"string\" value=\"Klu\"/> "
      "                       ";
  file_contents +=
      "        <ParameterList name=\"ParameterList\">                                              "
      "                       ";
  file_contents +=
      "          <Parameter name=\"Reindex\"                       type=\"bool\" value=\"true\"/>  "
      "                       ";
  file_contents +=
      "        </ParameterList>                                                                    "
      "                       ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                       ";
  file_contents +=
      "    </ParameterList>                                                                        "
      "                       ";
  file_contents +=
      "    <ParameterList name=\"Hierarchy\">                                                      "
      "                       ";
  file_contents +=
      "      <Parameter name=\"max levels\"               type=\"int\"      value=\"3\"/>          "
      "                       ";
  file_contents +=
      "      <Parameter name=\"coarse: max size\"         type=\"int\"      value=\"1\"/>          "
      "                       ";
  file_contents +=
      "      <Parameter name=\"verbosity\"                type=\"string\"   value=\"Low\"/>        "
      "                       ";
  file_contents +=
      "      <ParameterList name=\"All\">                                                          "
      "                       ";
  file_contents +=
      "        <Parameter name=\"startLevel\"             type=\"int\"      value=\"0\"/>          "
      "                       ";
  file_contents +=
      "        <Parameter name=\"Smoother\"               type=\"string\"   "
      "value=\"myForwardGaussSeidel\"/>              ";
  file_contents +=
      "        <Parameter name=\"CoarseSolver\"           type=\"string\"   "
      "value=\"mySolverFact\"/>                      ";
  file_contents +=
      "        <Parameter name=\"Aggregates\"             type=\"string\"   "
      "value=\"myUncoupledAggregationFact\"/>        ";
  file_contents +=
      "        <Parameter name=\"Graph\"                  type=\"string\"   "
      "value=\"myCoalesceDropFact\"/>                ";
  file_contents +=
      "        <Parameter name=\"A\"                      type=\"string\"   value=\"myRAPFact\"/>  "
      "                       ";
  file_contents +=
      "        <Parameter name=\"P\"                      type=\"string\"   "
      "value=\"myProlongatorFact\"/>                 ";
  file_contents +=
      "        <Parameter name=\"R\"                      type=\"string\"   "
      "value=\"myRestrictorFact\"/>                  ";
  file_contents +=
      "        <Parameter name=\"Nullspace\"              type=\"string\"   "
      "value=\"myTentativePFactory\"/>               ";
  file_contents +=
      "      </ParameterList>                                                                      "
      "                       ";
  file_contents +=
      "    </ParameterList>                                                                        "
      "                       ";
  file_contents +=
      "  </ParameterList>                                                                          "
      "                       ";
  file_contents +=
      "</ParameterList>                                                                            "
      "                       ";



  Teuchos::updateParametersFromXmlString(
      file_contents, Teuchos::Ptr<Teuchos::ParameterList>(&params));



  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
LINALG::SOLVER::AMGnxn_Operator::AMGnxn_Operator(Teuchos::RCP<BlockSparseMatrixBase> A,
    std::vector<int> num_pdes, std::vector<int> null_spaces_dim,
    std::vector<Teuchos::RCP<std::vector<double>>> null_spaces_data,
    const Teuchos::ParameterList& amgnxn_params, const Teuchos::ParameterList& smoothers_params,
    const Teuchos::ParameterList& muelu_params)
    : A_(A),
      num_pdes_(num_pdes),
      null_spaces_dim_(null_spaces_dim),
      null_spaces_data_(null_spaces_data),
      amgnxn_params_(amgnxn_params),
      smoothers_params_(smoothers_params),
      muelu_params_(muelu_params),
      is_setup_flag_(false)
{
  // Expected parameters in amgnxn_params (example)
  //<ParameterList name="amgnxn_params">
  //
  //  <Parameter name="number of levels"                 type="int"  value="..."/>
  //
  //  <Parameter name="smoother: all but coarsest level" type="string"  value="myFinestSmoother"/>
  //
  //  <Parameter name="smoother: coarsest level"         type="string"  value="myCoarsestSmoother"/>
  //
  //  <Parameter name="verbosity"                        type="string"  value="on"/>
  //
  //  <Parameter name="muelu parameters for block 0"       type="string"  value="myMuelu0"/>
  //
  //  <Parameter name="muelu parameters for block 1"       type="string"  value="myMuelu1"/>
  //
  //   ....
  //
  //  <Parameter name="muelu parameters for block N"       type="string"  value="myMueluN"/>
  //
  //</ParameterList>

  // Expected parameters in smoothers_params (example)
  //<ParameterList name="smoothers_params">
  //
  //  <ParameterList name="myFinestSmoother">
  //
  //   ...    ...    ...    ...    ...
  //
  //  </ParameterList>
  //
  //  <ParameterList name="myCoarsestSmoother">
  //
  //   ...    ...    ...    ...    ...
  //
  //  </ParameterList>
  //
  //</ParameterList>

  // Expected parameters in muelu_params (example)
  //<ParameterList name="muelu_params">
  //
  //   <ParameterList name="myMueluX">
  //     <Parameter name="xml file"      type="string"  value="myfile.xml"/>
  //   </ParameterList>
  //
  //   TODO or
  //
  //   <ParameterList name="myMueluX">
  //    ... ... list defining the muelue hierarcy (i.e.) the contents of the xml file
  //   </ParameterList>
  //
  //
  //</ParameterList>

  Setup();
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int LINALG::SOLVER::AMGnxn_Operator::ApplyInverse(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  TEUCHOS_FUNC_TIME_MONITOR("LINALG::SOLVER::AMGnxn_Operator::ApplyInverse");
  if (!is_setup_flag_)
    dserror("ApplyInverse cannot be called without a previous set up of the preconditioner");

  const MultiMapExtractor& range_ex = A_->RangeExtractor();
  const MultiMapExtractor& domain_ex = A_->DomainExtractor();

  int NumBlocks = A_->Rows();
  if (NumBlocks != A_->Cols()) dserror("The block matrix has to be square");

  AMGNXN::BlockedVector Xbl(NumBlocks);
  AMGNXN::BlockedVector Ybl(NumBlocks);

  int NV = X.NumVectors();
  for (int i = 0; i < NumBlocks; i++)
  {
    Teuchos::RCP<Epetra_MultiVector> Xi =
        Teuchos::rcp(new Epetra_MultiVector(*(range_ex.Map(i)), NV));

    Teuchos::RCP<Epetra_MultiVector> Yi =
        Teuchos::rcp(new Epetra_MultiVector(*(domain_ex.Map(i)), NV));

    range_ex.ExtractVector(X, i, *Xi);
    domain_ex.ExtractVector(Y, i, *Yi);
    Xbl.SetVector(Xi, i);
    Ybl.SetVector(Yi, i);
  }

  if (V_ == Teuchos::null) dserror("Null pointer. We cannot call the vcycle");

  V_->Solve(Xbl, Ybl, true);

  for (int i = 0; i < NumBlocks; i++) domain_ex.InsertVector(*(Ybl.GetVector(i)), i, Y);

  return 0;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_Operator::Setup()
{
  TEUCHOS_FUNC_TIME_MONITOR("LINALG::SOLVER::AMGnxn_Operator::Setup");


  int NumBlocks = A_->Rows();
  if (NumBlocks != A_->Cols()) dserror("We spect a square matrix here");

  // Extract the blockedMatrix
  Teuchos::RCP<AMGNXN::BlockedMatrix> Abl =
      Teuchos::rcp(new AMGNXN::BlockedMatrix(NumBlocks, NumBlocks));
  for (int i = 0; i < NumBlocks; i++)
  {
    for (int j = 0; j < NumBlocks; j++)
    {
      Teuchos::RCP<SparseMatrix> Aij = Teuchos::rcp(new SparseMatrix(A_->Matrix(i, j), View));
      Abl->SetMatrix(Aij, i, j);
    }
  }


  V_ = Teuchos::rcp(new AMGNXN::CoupledAmg(Abl, num_pdes_, null_spaces_dim_, null_spaces_data_,
      amgnxn_params_, smoothers_params_, muelu_params_));


  is_setup_flag_ = true;
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
LINALG::SOLVER::BlockSmoother_Operator::BlockSmoother_Operator(
    Teuchos::RCP<BlockSparseMatrixBase> A, std::vector<int> num_pdes,
    std::vector<int> null_spaces_dim,
    std::vector<Teuchos::RCP<std::vector<double>>> null_spaces_data,
    const Teuchos::ParameterList& amgnxn_params, const Teuchos::ParameterList& smoothers_params)
    : A_(A),
      num_pdes_(num_pdes),
      null_spaces_dim_(null_spaces_dim),
      null_spaces_data_(null_spaces_data),
      amgnxn_params_(amgnxn_params),
      smoothers_params_(smoothers_params),
      is_setup_flag_(false)
{
  // Expected parameters in amgnxn_params (example)
  //<ParameterList name="amgnxn_params">
  //
  //  <Parameter name="smoother"         type="string"  value="myBlockSmoother"/>
  //
  //  <Parameter name="verbosity"        type="string"  value="on"/>
  //
  //</ParameterList>

  // Expected parameters in smoothers_params (example)
  //<ParameterList name="smoothers_params">
  //
  //  <ParameterList name="myBlockSmoother">
  //
  //   ...    ...    ...    ...    ...
  //
  //  </ParameterList>
  //
  //
  //</ParameterList>

  Setup();
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int LINALG::SOLVER::BlockSmoother_Operator::ApplyInverse(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  TEUCHOS_FUNC_TIME_MONITOR("LINALG::SOLVER::BlockSmoother_Operator::ApplyInverse");

  if (!is_setup_flag_)
    dserror("ApplyInverse cannot be called without a previous set up of the preconditioner");


  const MultiMapExtractor& range_ex = A_->RangeExtractor();
  const MultiMapExtractor& domain_ex = A_->DomainExtractor();

  int NumBlocks = A_->Rows();
  if (NumBlocks != A_->Cols()) dserror("The block matrix has to be square");

  AMGNXN::BlockedVector Xbl(NumBlocks);
  AMGNXN::BlockedVector Ybl(NumBlocks);
  int NV = X.NumVectors();
  for (int i = 0; i < NumBlocks; i++)
  {
    Teuchos::RCP<Epetra_MultiVector> Xi =
        Teuchos::rcp(new Epetra_MultiVector(*(range_ex.Map(i)), NV));
    Teuchos::RCP<Epetra_MultiVector> Yi =
        Teuchos::rcp(new Epetra_MultiVector(*(domain_ex.Map(i)), NV));
    range_ex.ExtractVector(X, i, *Xi);
    domain_ex.ExtractVector(Y, i, *Yi);
    Xbl.SetVector(Xi, i);
    Ybl.SetVector(Yi, i);
  }


  if (S_ == Teuchos::null) dserror("Null pointer. We cannot call the smoother");

  // Ybl.Scale(0.0);
  S_->Solve(Xbl, Ybl, true);

  for (int i = 0; i < NumBlocks; i++) domain_ex.InsertVector(*(Ybl.GetVector(i)), i, Y);



  return 0;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::BlockSmoother_Operator::Setup()
{
  TEUCHOS_FUNC_TIME_MONITOR("LINALG::SOLVER::BlockSmoother_Operator::Setup");


  std::string verbosity = amgnxn_params_.get<std::string>("verbosity", "off");

  if (A_->Comm().MyPID() != 0) verbosity = "off";



  if (verbosity == "on")
  {
    std::cout << "===============================================" << std::endl;
    std::cout << "LINALG::SOLVER::BlockSmoother_Operator : debug info  (begin)" << std::endl;
    std::cout << std::endl;
  }


  // Prepare info
  int NumBlocks = A_->Rows();
  std::string smother_name = amgnxn_params_.get<std::string>("smoother", "BGS");
  std::vector<int> blocks(NumBlocks, 0);
  for (int i = 0; i < NumBlocks; i++) blocks[i] = i;
  std::vector<AMGNXN::NullSpaceInfo> null_space_blocks;
  for (int i = 0; i < NumBlocks; i++)
  {
    AMGNXN::NullSpaceInfo myNS(num_pdes_[i], null_spaces_dim_[i], null_spaces_data_[i]);
    null_space_blocks.push_back(myNS);
  }
  Teuchos::RCP<AMGNXN::BlockedMatrix> Abl =
      Teuchos::rcp(new AMGNXN::BlockedMatrix(NumBlocks, NumBlocks));
  for (int i = 0; i < NumBlocks; i++)
  {
    for (int j = 0; j < NumBlocks; j++)
    {
      Teuchos::RCP<SparseMatrix> Aij = Teuchos::rcp(new SparseMatrix(A_->Matrix(i, j), View));
      Abl->SetMatrix(Aij, i, j);
    }
  }

  // smoother factory
  AMGNXN::SmootherFactory mySmootherCreator;
  mySmootherCreator.SetOperator(Abl);
  mySmootherCreator.SetParamsSmoother(smoothers_params_);
  mySmootherCreator.SetLevel(0);
  mySmootherCreator.SetBlocks(blocks);
  mySmootherCreator.SetSmootherName(smother_name);
  mySmootherCreator.SetVerbosity(verbosity);
  mySmootherCreator.SetNullSpaceAllBlocks(null_space_blocks);

  // Create smoother
  Sbase_ = mySmootherCreator.Create();
  S_ = Teuchos::rcp_dynamic_cast<AMGNXN::BlockedSmoother>(Sbase_);
  if (S_ == Teuchos::null)
    dserror("We expect a blocked smoother. Fix the xml file defining the smoother");

  //// Print maps
  // for(int i=0;i<NumBlocks;i++)
  //{
  //  std::stringstream sstr;
  //  sstr << "Map_block" << i;
  //  PrintMap(A_->Matrix(i,i).RowMap(),sstr.str());
  //}


  if (verbosity == "on")
  {
    std::cout << std::endl;
    std::cout << "LINALG::SOLVER::BlockSmoother_Operator : debug info  (end)" << std::endl;
    std::cout << "===============================================" << std::endl;
  }

  is_setup_flag_ = true;
  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

LINALG::SOLVER::Merged_Operator::Merged_Operator(Teuchos::RCP<BlockSparseMatrixBase> A,
    const Teuchos::ParameterList& amgnxn_params, const Teuchos::ParameterList& smoothers_params)
    : A_(A),
      amgnxn_params_(amgnxn_params),
      smoothers_params_(smoothers_params),
      is_setup_flag_(false)
{
  // Expected parameters in amgnxn_params (example)
  //<ParameterList name="amgnxn_params">
  //  <Parameter name="verbosity"        type="string"  value="on"/>
  //  <Parameter name="smoother"         type="string"  value="mysmoother"/>
  //</ParameterList>
  //
  //
  // Expected parameters in smoothers_params (example)
  //<ParameterList name="smoothers_params">
  //
  //  <ParameterList name="mysmoother">
  //    <Parameter name="type"                           type="string"  value="point relaxation"/>
  //    <ParameterList name="ParameterList">
  //      <Parameter name="relaxation: type"             type="string"  value="Gauss-Seidel"/>
  //      <Parameter name="relaxation: backward mode"    type="bool"    value="false"/>
  //      <Parameter name="relaxation: sweeps"           type="int"     value="2"/>
  //      <Parameter name="relaxation: damping factor"   type="double"  value="1.0"/>
  //    </ParameterList>
  //  </ParameterList>
  //
  //</ParameterList>
  //
  //
  //
  //

  Setup();
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::Merged_Operator::Setup()
{
  TEUCHOS_FUNC_TIME_MONITOR("LINALG::SOLVER::Merged_Operator::Setup");

  // Read the parameter "verbosity" in amgnxn_params
  std::string verbosity = amgnxn_params_.get<std::string>("verbosity", "off");
  if (A_->Comm().MyPID() != 0) verbosity = "off";


  if (verbosity == "on")
  {
    std::cout << "===============================================" << std::endl;
    std::cout << "LINALG::SOLVER::Merged_Operator : debug info  (begin)" << std::endl;
    std::cout << std::endl;
  }

  // Merge the matrix
  Epetra_Time timer(A_->Comm());
  timer.ResetStartTime();

  Asp_ = A_->Merge();

  double elaptime = timer.ElapsedTime();
  if (A_->Comm().MyPID() == 0)
    std::cout << "   Merging the blocks takes " << std::setw(16) << std::setprecision(6) << elaptime
              << " s" << std::endl;

  // Safety check
  Teuchos::RCP<Epetra_Operator> Aop =
      Teuchos::rcp_dynamic_cast<Epetra_Operator>(Asp_->EpetraMatrix());
  Teuchos::RCP<Epetra_CrsMatrix> crsA = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(Aop);
  if (crsA == Teuchos::null) dserror("Houston, something went wrong in merging the matrix");


  // Read parameter called "smoother" in amgnxn_params
  std::string mysmoother = amgnxn_params_.get<std::string>("smoother", "none");
  if (mysmoother == "none") dserror("You have to set a parameter called smoother");

  // Read the parameters inside "smoothers_params"
  if (not smoothers_params_.isSublist(mysmoother))
    dserror("Not found a sublist with name %s", mysmoother.c_str());
  Teuchos::ParameterList myparams = smoothers_params_.sublist(mysmoother);

  // Some output
  if (verbosity == "on")
  {
    std::cout << std::endl;
    std::cout << "Creating an IFPACK smoother for a merged fine-level matrix" << std::endl;
    std::cout << "The Ifpack type is: " << myparams.get<std::string>("type") << std::endl;
    int overlap = myparams.get<int>("overlap", 0);
    std::cout << "The overlap is: " << overlap << std::endl;
    std::cout << "The parameters are: " << std::endl;
    std::cout << myparams.sublist("ParameterList");
  }

  // Create the smoother
  S_ = Teuchos::rcp(new AMGNXN::IfpackWrapper(Asp_, myparams));

  is_setup_flag_ = true;

  if (verbosity == "on")
  {
    std::cout << std::endl;
    std::cout << "LINALG::SOLVER::Merged_Operator : debug info  (end)" << std::endl;
    std::cout << "===============================================" << std::endl;
  }


  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int LINALG::SOLVER::Merged_Operator::ApplyInverse(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  TEUCHOS_FUNC_TIME_MONITOR("LINALG::SOLVER::Merged_Operator::ApplyInverse");

  if (!is_setup_flag_)
    dserror("ApplyInverse cannot be called without a previous set up of the preconditioner");

  if (S_ == Teuchos::null) dserror("Null pointer. We cannot call the smoother");

  Epetra_MultiVector Ybis(Y.Map(), Y.NumVectors(), false);
  S_->Apply(X, Ybis, true);
  Y.Update(1., Ybis, 0.);

  return 0;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::PrintMap(const Epetra_Map& Map, std::string prefix)
{
  const int pid = Map.Comm().MyPID();
  std::stringstream strr;
  strr << prefix << "_pid" << pid << ".txt";
  std::ofstream ofile(strr.str().c_str());
  int NumLID = Map.NumMyElements();
  int NumGID = Map.NumGlobalElements();
  ofile << "NumLID " << NumLID << " NumGID " << NumGID << std::endl;
  for (int LID = 0; LID < NumLID; LID++) ofile << LID << " " << Map.GID(LID) << std::endl;
  return;
}


#endif  // HAVE_MueLu
