/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_linear_solver_amgnxn_preconditioner.hpp"

#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linear_solver_amgnxn_vcycle.hpp"
#include "4C_utils_exceptions.hpp"

#include <EpetraExt_RowMatrixOut.h>
#include <MueLu_MLParameterListInterpreter_decl.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <Teuchos_PtrDecl.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include <iostream>

FOUR_C_NAMESPACE_OPEN

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Core::LinearSolver::AmGnxnPreconditioner::AmGnxnPreconditioner(Teuchos::ParameterList& params)
    : params_(params)
{
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<Epetra_Operator> Core::LinearSolver::AmGnxnPreconditioner::PrecOperator() const
{
  return p_;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void Core::LinearSolver::AmGnxnPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  // Decide if the setup has to be done
  if (!create) return;

  // Check whether this is a block sparse matrix
  Core::LinAlg::BlockSparseMatrixBase* A_bl =
      dynamic_cast<Core::LinAlg::BlockSparseMatrixBase*>(matrix);
  if (A_bl == nullptr)
    FOUR_C_THROW(
        "The AMGnxn preconditioner works only for BlockSparseMatrixBase or derived classes");

  // Do all the setup
  Setup(Teuchos::rcp(A_bl, false));

  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void Core::LinearSolver::AmGnxnPreconditioner::Setup(
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> A)
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::LinAlg::SOLVER::AMGnxn_Preconditioner::Setup");

  Teuchos::Time timer("", true);
  timer.reset();

  // Free old matrix and preconditioner
  a_ = Teuchos::null;
  p_ = Teuchos::null;

  // Create own copy of the system matrix in order to allow reusing the preconditioner
  a_ = A;
  a_ = a_->Clone(Core::LinAlg::Copy);
  a_->Complete();

  // Determine number of blocks
  int NumBlocks = a_->Rows();
  if (a_->Rows() != a_->Cols())
    FOUR_C_THROW("The AMGnxn preconditioner works only for block square matrices");

  // Pick-up the input parameters
  AmGnxnInterface myInterface(params_, NumBlocks);

  // Create the Operator
  if (myInterface.get_preconditioner_type() == "AMG(BlockSmoother)")
  {
    p_ = Teuchos::rcp(
        new AmGnxnOperator(a_, myInterface.GetNumPdes(), myInterface.GetNullSpacesDim(),
            myInterface.GetNullSpacesData(), myInterface.get_preconditioner_params(),
            myInterface.GetSmoothersParams(), myInterface.GetSmoothersParams()));
  }
  else if (myInterface.get_preconditioner_type() == "BlockSmoother(X)")
  {
    p_ = Teuchos::rcp(new BlockSmootherOperator(a_, myInterface.GetNumPdes(),
        myInterface.GetNullSpacesDim(), myInterface.GetNullSpacesData(),
        myInterface.get_preconditioner_params(), myInterface.GetSmoothersParams()));
  }
  else if (myInterface.get_preconditioner_type() == "Merge_and_Ifpack")
  {
    p_ = Teuchos::rcp(new MergedOperator(
        a_, myInterface.get_preconditioner_params(), myInterface.GetSmoothersParams()));
  }
  else
    FOUR_C_THROW("Unknown preconditioner type: %s", myInterface.get_preconditioner_type().c_str());

  double elaptime = timer.totalElapsedTime(true);
  if (myInterface.get_preconditioner_params().get<std::string>("verbosity", "off") == "on" and
      A->Comm().MyPID() == 0)
    std::cout << "       Calling Core::LinAlg::SOLVER::AMGnxn_Preconditioner::Setup takes "
              << std::setw(16) << std::setprecision(6) << elaptime << " s" << std::endl;

  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Core::LinearSolver::AmGnxnInterface::AmGnxnInterface(Teuchos::ParameterList& params, int NumBlocks)
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

  if (!params.isSublist("AMGnxn Parameters")) FOUR_C_THROW("AMGnxn Parameters not found!");
  Teuchos::ParameterList& amglist = params.sublist("AMGnxn Parameters");


  // Decide whether to choose a default type or parse a xml file
  std::string amgnxn_type = amglist.get<std::string>("AMGNXN_TYPE", "AMG(BGS)");
  if (amgnxn_type == "XML")
  {
    // Parse the whole file
    std::string amgnxn_xml = amglist.get<std::string>("AMGNXN_XML_FILE", "none");
    if (amgnxn_xml == "none") FOUR_C_THROW("The input parameter AMGNXN_XML_FILE is empty.");
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
    FOUR_C_THROW(
        "\"%s\" is an invalid value for \"AMGNXN_TYPE\". Fix your .dat", amgnxn_type.c_str());



  // Find preconditioner type and parameters
  std::string myprec = smoo_params_.get<std::string>("Preconditioner", "none");
  if (myprec == "none") FOUR_C_THROW("Not found \"Preconditioner\" parameter in your xml file.");
  if (!smoo_params_.isSublist(myprec))
    FOUR_C_THROW("Not found your preconditioner list in your xml file.");
  Teuchos::ParameterList& myprec_list = smoo_params_.sublist(myprec);
  prec_type_ = myprec_list.get<std::string>("type", "none");
  if (!myprec_list.isSublist("parameters"))
    FOUR_C_THROW("Not found the parameters list for your preconditioner. Fix your xml file.");
  prec_params_ = myprec_list.sublist("parameters");

  // Find null spaces and relatives
  std::string Inverse_str = "Inverse";
  xml_files_.resize(NumBlocks);
  num_pdes_.resize(NumBlocks);
  null_spaces_dim_.resize(NumBlocks);
  null_spaces_data_.resize(NumBlocks);
  for (int block = 0; block < NumBlocks; block++)
  {
    if (!params.isSublist(Inverse_str + convert_int(block + 1)))
      FOUR_C_THROW("Not found inverse list for block %d", block + 1);
    Teuchos::ParameterList& inverse_list = params.sublist(Inverse_str + convert_int(block + 1));

    if (!inverse_list.isSublist("MueLu Parameters")) FOUR_C_THROW("MueLu Parameters not found");
    Teuchos::ParameterList& mllist = inverse_list.sublist("MueLu Parameters");

    xml_files_[block] = mllist.get<std::string>("xml file", "none");
    num_pdes_[block] = mllist.get<int>("PDE equations", -1);
    null_spaces_dim_[block] = mllist.get<int>("null space: dimension", -1);

    Teuchos::RCP<Epetra_MultiVector> nullspace =
        mllist.get<Teuchos::RCP<Epetra_MultiVector>>("nullspace", Teuchos::null);
    if (nullspace == Teuchos::null) FOUR_C_THROW("Nullspace vector is null!");

    Teuchos::RCP<std::vector<double>> ns =
        Teuchos::rcp(new std::vector<double>(nullspace->MyLength() * nullspace->NumVectors()));

    Core::LinAlg::EpetraMultiVectorToStdVector(nullspace, *ns, null_spaces_dim_[block]);
    null_spaces_data_[block] = ns;

    // Some checks
    if (num_pdes_[block] < 1 or null_spaces_dim_[block] < 1)
      FOUR_C_THROW("Error: PDE equations or null space dimension wrong.");
    if (null_spaces_data_[block] == Teuchos::null) FOUR_C_THROW("Error: null space data is empty");
  }
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
Core::LinearSolver::AmGnxnOperator::AmGnxnOperator(
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> A, std::vector<int> num_pdes,
    std::vector<int> null_spaces_dim,
    std::vector<Teuchos::RCP<std::vector<double>>> null_spaces_data,
    const Teuchos::ParameterList& amgnxn_params, const Teuchos::ParameterList& smoothers_params,
    const Teuchos::ParameterList& muelu_params)
    : a_(A),
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

int Core::LinearSolver::AmGnxnOperator::ApplyInverse(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::LinAlg::SOLVER::AMGnxn_Operator::ApplyInverse");
  if (!is_setup_flag_)
    FOUR_C_THROW("ApplyInverse cannot be called without a previous set up of the preconditioner");

  const Core::LinAlg::MultiMapExtractor& range_ex = a_->RangeExtractor();
  const Core::LinAlg::MultiMapExtractor& domain_ex = a_->DomainExtractor();

  int NumBlocks = a_->Rows();
  if (NumBlocks != a_->Cols()) FOUR_C_THROW("The block matrix has to be square");

  AMGNxN::BlockedVector Xbl(NumBlocks);
  AMGNxN::BlockedVector Ybl(NumBlocks);

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

  if (v_ == Teuchos::null) FOUR_C_THROW("Null pointer. We cannot call the vcycle");

  v_->Solve(Xbl, Ybl, true);

  for (int i = 0; i < NumBlocks; i++) domain_ex.InsertVector(*(Ybl.GetVector(i)), i, Y);

  return 0;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void Core::LinearSolver::AmGnxnOperator::Setup()
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::LinAlg::SOLVER::AMGnxn_Operator::Setup");


  int NumBlocks = a_->Rows();
  if (NumBlocks != a_->Cols()) FOUR_C_THROW("We spect a square matrix here");

  // Extract the blockedMatrix
  Teuchos::RCP<AMGNxN::BlockedMatrix> Abl =
      Teuchos::rcp(new AMGNxN::BlockedMatrix(NumBlocks, NumBlocks));
  for (int i = 0; i < NumBlocks; i++)
  {
    for (int j = 0; j < NumBlocks; j++)
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> Aij =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(a_->Matrix(i, j), Core::LinAlg::View));
      Abl->SetMatrix(Aij, i, j);
    }
  }


  v_ = Teuchos::rcp(new AMGNxN::CoupledAmg(Abl, num_pdes_, null_spaces_dim_, null_spaces_data_,
      amgnxn_params_, smoothers_params_, muelu_params_));


  is_setup_flag_ = true;
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
Core::LinearSolver::BlockSmootherOperator::BlockSmootherOperator(
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> A, std::vector<int> num_pdes,
    std::vector<int> null_spaces_dim,
    std::vector<Teuchos::RCP<std::vector<double>>> null_spaces_data,
    const Teuchos::ParameterList& amgnxn_params, const Teuchos::ParameterList& smoothers_params)
    : a_(A),
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

int Core::LinearSolver::BlockSmootherOperator::ApplyInverse(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::LinAlg::SOLVER::BlockSmoother_Operator::ApplyInverse");

  if (!is_setup_flag_)
    FOUR_C_THROW("ApplyInverse cannot be called without a previous set up of the preconditioner");


  const Core::LinAlg::MultiMapExtractor& range_ex = a_->RangeExtractor();
  const Core::LinAlg::MultiMapExtractor& domain_ex = a_->DomainExtractor();

  int NumBlocks = a_->Rows();
  if (NumBlocks != a_->Cols()) FOUR_C_THROW("The block matrix has to be square");

  AMGNxN::BlockedVector Xbl(NumBlocks);
  AMGNxN::BlockedVector Ybl(NumBlocks);
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


  if (s_ == Teuchos::null) FOUR_C_THROW("Null pointer. We cannot call the smoother");

  s_->Solve(Xbl, Ybl, true);

  for (int i = 0; i < NumBlocks; i++) domain_ex.InsertVector(*(Ybl.GetVector(i)), i, Y);



  return 0;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void Core::LinearSolver::BlockSmootherOperator::Setup()
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::LinAlg::SOLVER::BlockSmoother_Operator::Setup");


  std::string verbosity = amgnxn_params_.get<std::string>("verbosity", "off");

  if (a_->Comm().MyPID() != 0) verbosity = "off";



  if (verbosity == "on")
  {
    std::cout << "===============================================" << std::endl;
    std::cout << "Core::LinAlg::SOLVER::BlockSmoother_Operator : debug info  (begin)" << std::endl;
    std::cout << std::endl;
  }


  // Prepare info
  int NumBlocks = a_->Rows();
  std::string smother_name = amgnxn_params_.get<std::string>("smoother", "BGS");
  std::vector<int> blocks(NumBlocks, 0);
  for (int i = 0; i < NumBlocks; i++) blocks[i] = i;
  std::vector<AMGNxN::NullSpaceInfo> null_space_blocks;
  for (int i = 0; i < NumBlocks; i++)
  {
    AMGNxN::NullSpaceInfo myNS(num_pdes_[i], null_spaces_dim_[i], null_spaces_data_[i]);
    null_space_blocks.push_back(myNS);
  }
  Teuchos::RCP<AMGNxN::BlockedMatrix> Abl =
      Teuchos::rcp(new AMGNxN::BlockedMatrix(NumBlocks, NumBlocks));
  for (int i = 0; i < NumBlocks; i++)
  {
    for (int j = 0; j < NumBlocks; j++)
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> Aij =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(a_->Matrix(i, j), Core::LinAlg::View));
      Abl->SetMatrix(Aij, i, j);
    }
  }

  // smoother factory
  AMGNxN::SmootherFactory mySmootherCreator;
  mySmootherCreator.SetOperator(Abl);
  mySmootherCreator.SetParamsSmoother(smoothers_params_);
  mySmootherCreator.SetLevel(0);
  mySmootherCreator.SetBlocks(blocks);
  mySmootherCreator.SetSmootherName(smother_name);
  mySmootherCreator.SetVerbosity(verbosity);
  mySmootherCreator.set_null_space_all_blocks(null_space_blocks);

  // Create smoother
  sbase_ = mySmootherCreator.Create();
  s_ = Teuchos::rcp_dynamic_cast<AMGNxN::BlockedSmoother>(sbase_);
  if (s_ == Teuchos::null)
    FOUR_C_THROW("We expect a blocked smoother. Fix the xml file defining the smoother");

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
    std::cout << "Core::LinAlg::SOLVER::BlockSmoother_Operator : debug info  (end)" << std::endl;
    std::cout << "===============================================" << std::endl;
  }

  is_setup_flag_ = true;
  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Core::LinearSolver::MergedOperator::MergedOperator(
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> A,
    const Teuchos::ParameterList& amgnxn_params, const Teuchos::ParameterList& smoothers_params)
    : a_(A),
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

void Core::LinearSolver::MergedOperator::Setup()
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::LinAlg::SOLVER::Merged_Operator::Setup");

  // Read the parameter "verbosity" in amgnxn_params
  std::string verbosity = amgnxn_params_.get<std::string>("verbosity", "off");
  if (a_->Comm().MyPID() != 0) verbosity = "off";


  if (verbosity == "on")
  {
    std::cout << "===============================================" << std::endl;
    std::cout << "Core::LinAlg::SOLVER::Merged_Operator : debug info  (begin)" << std::endl;
    std::cout << std::endl;
  }

  // Merge the matrix
  Teuchos::Time timer("", true);
  timer.reset();

  asp_ = a_->Merge();

  double elaptime = timer.totalElapsedTime(true);
  if (a_->Comm().MyPID() == 0)
    std::cout << "   Merging the blocks takes " << std::setw(16) << std::setprecision(6) << elaptime
              << " s" << std::endl;

  // Safety check
  Teuchos::RCP<Epetra_Operator> Aop =
      Teuchos::rcp_dynamic_cast<Epetra_Operator>(asp_->EpetraMatrix());
  Teuchos::RCP<Epetra_CrsMatrix> crsA = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(Aop);
  if (crsA == Teuchos::null) FOUR_C_THROW("Houston, something went wrong in merging the matrix");


  // Read parameter called "smoother" in amgnxn_params
  std::string mysmoother = amgnxn_params_.get<std::string>("smoother", "none");
  if (mysmoother == "none") FOUR_C_THROW("You have to set a parameter called smoother");

  // Read the parameters inside "smoothers_params"
  if (not smoothers_params_.isSublist(mysmoother))
    FOUR_C_THROW("Not found a sublist with name %s", mysmoother.c_str());
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
  s_ = Teuchos::rcp(new AMGNxN::IfpackWrapper(asp_, myparams));

  is_setup_flag_ = true;

  if (verbosity == "on")
  {
    std::cout << std::endl;
    std::cout << "Core::LinAlg::SOLVER::Merged_Operator : debug info  (end)" << std::endl;
    std::cout << "===============================================" << std::endl;
  }


  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int Core::LinearSolver::MergedOperator::ApplyInverse(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::LinAlg::SOLVER::Merged_Operator::ApplyInverse");

  if (!is_setup_flag_)
    FOUR_C_THROW("ApplyInverse cannot be called without a previous set up of the preconditioner");

  if (s_ == Teuchos::null) FOUR_C_THROW("Null pointer. We cannot call the smoother");

  Epetra_MultiVector Ybis(Y.Map(), Y.NumVectors(), false);
  s_->Apply(X, Ybis, true);
  Y.Update(1., Ybis, 0.);

  return 0;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void Core::LinearSolver::PrintMap(const Epetra_Map& Map, std::string prefix)
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

FOUR_C_NAMESPACE_CLOSE
