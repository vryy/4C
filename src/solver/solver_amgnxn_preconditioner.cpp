/*!----------------------------------------------------------------------
\file solver_amgnxn_preconditioner.cpp

<pre>
Maintainer: Francesc Verdugo
            verdugo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
Created on: Feb 27, 2014
</pre>
*----------------------------------------------------------------------*/

#ifdef HAVE_MueLu
#ifdef HAVE_Trilinos_Q1_2014

#include <iostream>

//#include "../drt_lib/drt_globalproblem.H"
#include <Teuchos_PtrDecl.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <MueLu_MLParameterListInterpreter_decl.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include "EpetraExt_RowMatrixOut.h"
#include "../drt_lib/drt_dserror.H"
#include "solver_amgnxn_preconditioner.H"


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

LINALG::SOLVER::AMGnxn_Preconditioner::AMGnxn_Preconditioner
(
  FILE * outfile,
  Teuchos::ParameterList & params
) :
  LINALG::SOLVER::PreconditionerType( outfile ),
  params_(params)
{}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Epetra_Operator * LINALG::SOLVER::AMGnxn_Preconditioner::PrecOperator() const
{
  return &*P_;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_Preconditioner::Setup
(
  bool create,
  Epetra_Operator * matrix,
  Epetra_MultiVector * x,
  Epetra_MultiVector * b
)
{

  // Setup underlying linear system
  SetupLinearProblem( matrix, x, b );

  // Decide if the setup has to be done
  if(!create)
    return;

  // Check whether this is a block sparse matrix
  BlockSparseMatrixBase* A_bl = dynamic_cast<BlockSparseMatrixBase*>(matrix);
  if(A_bl==NULL)
    dserror("The AMGnxn preconditioner works only for BlockSparseMatrixBase or derived classes");

  // Do all the setup
  Setup(Teuchos::rcp(A_bl,false));

  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_Preconditioner::Setup(Teuchos::RCP<BlockSparseMatrixBase> A)
{

  // Free old matrix and preconditioner
  A_ = Teuchos::null;
  P_ = Teuchos::null;

  // Create own copy of the system matrix in order to allow reusing the preconditioner
  A_ = A;
  A_ = A_->Clone(Copy);
  A_->Complete();

  // Determine number of blocks
  int NumBlocks = A_->Rows();
  if(A_->Rows() != A_->Cols())
    dserror("The AMGnxn preconditioner works only for block square matrices");

  // Pick-up the input parameters
  AMGnxn_Interface myInterface(params_,NumBlocks);

  // Create the Operator
  if (myInterface.GetPreconditionerType() == "AMG(BlockSmoother)")
  {
    P_ = Teuchos::rcp(new AMGnxn_Operator(
          A_,
          myInterface.GetNumPdes(),
          myInterface.GetNullSpacesDim(),
          myInterface.GetNullSpacesData(),
          myInterface.GetPreconditionerParams(),
          myInterface.GetSmoothersParams()));
  }
  else if (myInterface.GetPreconditionerType() == "BlockSmoother(X)")
  {
    P_ = Teuchos::rcp(new BlockSmoother_Operator(
          A_,
          myInterface.GetNumPdes(),
          myInterface.GetNullSpacesDim(),
          myInterface.GetNullSpacesData(),
          myInterface.GetPreconditionerParams(),
          myInterface.GetSmoothersParams()));
  }
  else
    dserror("Unknown preconditioner type");

  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

LINALG::SOLVER::AMGnxn_Interface::AMGnxn_Interface(Teuchos::ParameterList& params,int NumBlocks)
{

  // Expected parameters in params
  //<ParameterList name="params">
  //
  //  <ParameterList name="AMGnxn Parameters">
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
  //  <!-- Here we which preconditioner we are going to use -->
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

  if(!params.isSublist("AMGnxn Parameters"))
    dserror("AMGnxn Parameters not found!");
  Teuchos::ParameterList& amglist = params.sublist("AMGnxn Parameters");

  // Parse the whole file
  std::string amgnxn_xml =amglist.get<std::string>("AMGNXN_XML_FILE","none");
  if(amgnxn_xml=="none")
    dserror("The input parameter AMGNXN_XML_FILE is empty.");
  if(not(amgnxn_xml=="none"))
    Teuchos::updateParametersFromXmlFile(
        amgnxn_xml,Teuchos::Ptr<Teuchos::ParameterList>(&smoo_params_));

  // Find preconditioner type and parameters
  std::string myprec =  smoo_params_.get<std::string>("Preconditioner","none");
    if(myprec == "none")
      dserror("Not found \"Preconditioner\" parameter in your xml file.");
  if(!smoo_params_.isSublist(myprec))
    dserror("Not found your preconditioner list in your xml file.");
  Teuchos::ParameterList& myprec_list = smoo_params_.sublist(myprec);
  prec_type_ =  myprec_list.get<std::string>("type","none");
  if(!myprec_list.isSublist("parameters"))
    dserror("Not found the parameters list for your preconditioner. Fix your xml file.");
  prec_params_ = myprec_list.sublist("parameters");

  // Find null spaces and relatives
  std::string Inverse_str = "Inverse";
  xml_files_.resize(NumBlocks);
  num_pdes_.resize(NumBlocks);
  null_spaces_dim_.resize(NumBlocks);
  null_spaces_data_.resize(NumBlocks);
  for(int block=0;block<NumBlocks;block++)
  {

    if(!params.isSublist(Inverse_str + ConvertInt(block+1)))
      dserror("Not found inverse list for block %d", block+1);
    Teuchos::ParameterList& inverse_list = params.sublist(Inverse_str + ConvertInt(block+1));

    if(!inverse_list.isSublist("MueLu Parameters"))
      dserror("MueLu Parameters not found");
    Teuchos::ParameterList& mllist = inverse_list.sublist("MueLu Parameters");

    xml_files_[block] = mllist.get<std::string>("xml file","none");
    num_pdes_[block] = mllist.get<int>("PDE equations",-1);
    null_spaces_dim_[block] = mllist.get<int>("null space: dimension",-1);
    null_spaces_data_[block] =
      mllist.get<Teuchos::RCP<std::vector<double> > >("nullspace",Teuchos::null);

    //Some cheks
    if(num_pdes_[block]<1 or num_pdes_[block]<1)
      dserror("Error: PDE equations or null space dimension wrong.");
    if(null_spaces_data_[block]==Teuchos::null)
      dserror("Error: null space data is empty");
  }


}



/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
LINALG::SOLVER::AMGnxn_Operator::AMGnxn_Operator(
    Teuchos::RCP<BlockSparseMatrixBase> A,
    std::vector<int> num_pdes,
    std::vector<int> null_spaces_dim,
    std::vector<Teuchos::RCP<std::vector<double> > > null_spaces_data,
    const Teuchos::ParameterList& amgnxn_params,
    const Teuchos::ParameterList& smoothers_params):
    A_                   (A               ),
    num_pdes_            (num_pdes        ),
    null_spaces_dim_     (null_spaces_dim ),
    null_spaces_data_    (null_spaces_data),
    amgnxn_params_       (amgnxn_params   ),
    smoothers_params_    (smoothers_params),
    is_setup_flag_       (false           )
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
  //  <Parameter name="muelu xml file for block 0"       type="string"  value="file0.xml"/>
  //
  //  <Parameter name="muelu xml file for block 1"       type="string"  value="file1.xml"/>
  //
  //   ....
  //
  //  <Parameter name="muelu xml file for block N"       type="string"  value="fileN.xml"/>
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

  Setup();
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int  LINALG::SOLVER::AMGnxn_Operator::ApplyInverse(
    const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  if(!is_setup_flag_)
    dserror("ApplyInverse cannot be called without a previous set up of the preconditioner");
  V_->Apply(X,Y);
  return 0;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void  LINALG::SOLVER::AMGnxn_Operator::Setup()
{

  std::string verbosity = amgnxn_params_.get<std::string>("verbosity","off");
  if (verbosity=="on")
  {
    std::cout << "===============================================" << std::endl;
    std::cout << "LINALG::SOLVER::AMGnxn_Operator : debug info  (begin)" << std::endl;
    std::cout << std::endl;
    std::cout << ">>>>>> Creating MueLu AMG Hierarchies for each one of the blocks" << std::endl;
    std::cout << std::endl;
  }

  // Parse the name of the xml files
  int NumBlocks = A_->Rows();
  for(int i=0;i<NumBlocks;i++)
  {
    std::stringstream ss;
    ss << i;
    std::string param_name = "muelu xml file for block " + ss.str();
    std::string xml_name =amgnxn_params_.get<std::string>(param_name,"none");
    if (xml_name == "none")
      dserror("You must specify an xml file for creating the AMG on block %d",i);
    xml_files_.push_back(xml_name);
  }


  int num_levels_amg = amgnxn_params_.get<int>("number of levels",20);
  //if(num_levels_amg==-1)
  //  dserror("Missing \"number of levels\" in your xml file");
  H_ = Teuchos::rcp(new  AMGnxn_Hierarchies(
        A_,
        xml_files_,
        num_pdes_,
        null_spaces_dim_,
        null_spaces_data_,
        num_levels_amg,
        verbosity));

  if (verbosity=="on")
  {
    std::cout << std::endl;
    std::cout << ">>>>>> Creating the monolithic hierarchy" << std::endl;
    std::cout << std::endl;
  }

  M_ = Teuchos::rcp( new AMGnxn_MonolithicHierarchy( H_,amgnxn_params_,smoothers_params_));
  V_ = M_->BuildVCycle();

  if (verbosity=="on")
  {
    std::cout << std::endl;
    std::cout << "LINALG::SOLVER::AMGnxn_Operator : debug info  (end)" << std::endl;
    std::cout << "===============================================" << std::endl;
  }

  is_setup_flag_=true;
  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
LINALG::SOLVER::BlockSmoother_Operator::BlockSmoother_Operator(
    Teuchos::RCP<BlockSparseMatrixBase> A,
    std::vector<int> num_pdes,
    std::vector<int> null_spaces_dim,
    std::vector<Teuchos::RCP<std::vector<double> > > null_spaces_data,
    const Teuchos::ParameterList& amgnxn_params,
    const Teuchos::ParameterList& smoothers_params):
    A_                   (A               ),
    num_pdes_            (num_pdes        ),
    null_spaces_dim_     (null_spaces_dim ),
    null_spaces_data_    (null_spaces_data),
    amgnxn_params_       (amgnxn_params   ),
    smoothers_params_    (smoothers_params),
    is_setup_flag_       (false           )
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

int  LINALG::SOLVER::BlockSmoother_Operator::ApplyInverse(
    const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  if(!is_setup_flag_)
    dserror("ApplyInverse cannot be called without a previous set up of the preconditioner");
  if(Sbase_==Teuchos::null) dserror("Something wrong");
  Sbase_->Apply(X,Y);
  // TODO This does not work, Why??
  //if(S_==Teuchos::null) dserror("Something wrong");
  //S_->Apply(X,Y);
  return 0;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void  LINALG::SOLVER::BlockSmoother_Operator::Setup()
{

  std::string verbosity = amgnxn_params_.get<std::string>("verbosity","off");
  if (verbosity=="on")
  {
    std::cout << "===============================================" << std::endl;
    std::cout << "LINALG::SOLVER::BlockSmoother_Operator : debug info  (begin)" << std::endl;
    std::cout << std::endl;
  }


  // Prepare info
  int NumBlocks = A_->Rows();
  std::string smother_name = amgnxn_params_.get<std::string>("smoother","BGS");
  std::vector<int> blocks(NumBlocks,0);
  for(int i=0;i<NumBlocks;i++)
    blocks[i]=i;
  std::vector<NullSpaceInfo> null_space_blocks;
  for(int i=0;i<NumBlocks;i++)
  {
    NullSpaceInfo myNS(num_pdes_[i],null_spaces_dim_[i],null_spaces_data_[i]);
    null_space_blocks.push_back(myNS);
  }

  // smoother factory
  AMGnxn_SmootherFactory mySmootherCreator;
  mySmootherCreator.SetOperator(A_);
  mySmootherCreator.SetParamsSmoother(smoothers_params_);
  mySmootherCreator.SetLevel(0);
  mySmootherCreator.SetBlocks(blocks);
  mySmootherCreator.SetSmootherName(smother_name);
  mySmootherCreator.SetVerbosity(verbosity);
  mySmootherCreator.SetNullSpaceAllBlocks(null_space_blocks);

  // Create smoother
  Sbase_ = mySmootherCreator.Create();
  Teuchos::RCP<AMGnxn_BlockSmootherBase>
    S_ = Teuchos::rcp_dynamic_cast<AMGnxn_BlockSmootherBase>(Sbase_);
  if(S_ == Teuchos::null)
    dserror("Something wrong. Fix the xml file defining the smoother");


  if (verbosity=="on")
  {
    std::cout << std::endl;
    std::cout << "LINALG::SOLVER::BlockSmoother_Operator : debug info  (end)" << std::endl;
    std::cout << "===============================================" << std::endl;
  }

  is_setup_flag_=true;
  return;
}




#endif // HAVE_MueLu
#endif // HAVE_Trilinos_Q1_2014

