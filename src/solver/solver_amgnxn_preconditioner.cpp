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

#include <iostream>

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
          myInterface.GetSmoothersParams(),
          myInterface.GetSmoothersParams()
          ));
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


  // Decide whether to choose a default type or parse a xml file
  std::string amgnxn_type =amglist.get<std::string>("AMGNXN_TYPE","AMG(BGS)");
  if(amgnxn_type=="AMG(BGS)") Params_AMG_BGS(smoo_params_);
  else if (amgnxn_type=="XML")
  {
    // Parse the whole file
    std::string amgnxn_xml =amglist.get<std::string>("AMGNXN_XML_FILE","none");
    if(amgnxn_xml=="none")
      dserror("The input parameter AMGNXN_XML_FILE is empty.");
    if(not(amgnxn_xml=="none"))
      Teuchos::updateParametersFromXmlFile(
          amgnxn_xml,Teuchos::Ptr<Teuchos::ParameterList>(&smoo_params_));
  }
  else dserror("\"%s\" is an invalid value for \"AMGNXN_TYPE\". Fix your .dat",amgnxn_type.c_str());



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
void LINALG::SOLVER::AMGnxn_Interface::Params_AMG_BGS(Teuchos::ParameterList& params)
{

// This is a pre-defined AMG(BGS) preconditioner
// TODO Add the possibility to the user of tuning some of the
// more relevant parameters using the existing input params in the .dat file

// This is the xml file which defines the preconditioner

//<ParameterList name="dummy">
//
//  <Parameter name="Preconditioner" type="string"  value="default_AMG(BGS)"/>
//
//  <!-- Global setup of the preconditioner  -->
//
//  <ParameterList name="default_AMG(BGS)">
//  <Parameter name="type" type="string"  value="AMG(BlockSmoother)"/>
//    <ParameterList name="parameters">
//      <Parameter name="number of levels"                 type="int"     value="5"                 />
//      <Parameter name="smoother: all but coarsest level" type="string"  value="BGS"               />
//      <Parameter name="smoother: coarsest level"         type="string"  value="MERGE_AND_SOLVE"   />
//      <Parameter name="verbosity"                        type="string"  value="off"               />
//      <Parameter name="muelu parameters for block 0"     type="string"  value="myMueluLISTDefault"/>
//      <Parameter name="muelu parameters for block 1"     type="string"  value="myMueluLISTDefault"/>
//      <Parameter name="muelu parameters for block 2"     type="string"  value="myMueluLISTDefault"/>
//      <Parameter name="muelu parameters for block 3"     type="string"  value="myMueluLISTDefault"/>
//      <Parameter name="muelu parameters for block 5"     type="string"  value="myMueluLISTDefault"/>
//      <Parameter name="muelu parameters for block 6"     type="string"  value="myMueluLISTDefault"/>
//      <Parameter name="muelu parameters for block 7"     type="string"  value="myMueluLISTDefault"/>
//      <Parameter name="muelu parameters for block 8"     type="string"  value="myMueluLISTDefault"/>
//      <Parameter name="muelu parameters for block 9"     type="string"  value="myMueluLISTDefault"/>
//    </ParameterList>
//  </ParameterList>
//
//  <!-- Muelu Hierarachies  -->
//
//  <ParameterList name="myMueluLISTDefault">
//    <ParameterList name="Factories">
//      <ParameterList name="myForwardGaussSeidel">
//        <Parameter name="factory"                      type="string"  value="TrilinosSmoother"/>
//        <Parameter name="type"                         type="string"  value="RELAXATION"      />
//        <ParameterList name="ParameterList">
//          <Parameter name="relaxation: type"           type="string"  value="Gauss-Seidel"/>
//        </ParameterList>
//      </ParameterList>
//    </ParameterList>
//    <ParameterList name="Hierarchy">
//      <Parameter name="numDesiredLevel" type="int"      value="3"   />
//      <Parameter name="maxCoarseSize"   type="int"      value="700" />
//      <Parameter name="verbosity"       type="string"   value="None"/>
//      <ParameterList name="All">
//      <Parameter name="startLevel"   type="int"      value="0"                   />
//      <Parameter name="Smoother"     type="string"   value="myForwardGaussSeidel"/>
//      <Parameter name="CoarseSolver" type="string"   value="myForwardGaussSeidel"/>
//      </ParameterList>
//    </ParameterList>
//  </ParameterList>
//
//
//</ParameterList>

  // This is the parameter list

  params.set<std::string>("Preconditioner","default_AMG(BGS)");

  Teuchos::ParameterList& default_AMG_BGS = params.sublist("default_AMG(BGS)");
  default_AMG_BGS.set<std::string>("type","AMG(BlockSmoother)");
  Teuchos::ParameterList& default_AMG_BGS_params = default_AMG_BGS.sublist("parameters");
  default_AMG_BGS_params.set             ("number of levels"                ,5                   );
  default_AMG_BGS_params.set<std::string>("smoother: all but coarsest level","BGS"               );
  default_AMG_BGS_params.set<std::string>("smoother: coarsest level"        ,"MERGE_AND_SOLVE"   );
  default_AMG_BGS_params.set<std::string>("verbosity"                       ,"off"               );
  default_AMG_BGS_params.set<std::string>("muelu parameters for block 0"    ,"myMueluLISTDefault");
  default_AMG_BGS_params.set<std::string>("muelu parameters for block 1"    ,"myMueluLISTDefault");
  default_AMG_BGS_params.set<std::string>("muelu parameters for block 2"    ,"myMueluLISTDefault");
  default_AMG_BGS_params.set<std::string>("muelu parameters for block 3"    ,"myMueluLISTDefault");
  default_AMG_BGS_params.set<std::string>("muelu parameters for block 5"    ,"myMueluLISTDefault");
  default_AMG_BGS_params.set<std::string>("muelu parameters for block 6"    ,"myMueluLISTDefault");
  default_AMG_BGS_params.set<std::string>("muelu parameters for block 7"    ,"myMueluLISTDefault");
  default_AMG_BGS_params.set<std::string>("muelu parameters for block 8"    ,"myMueluLISTDefault");
  default_AMG_BGS_params.set<std::string>("muelu parameters for block 9"    ,"myMueluLISTDefault");

  Teuchos::ParameterList& muelu = params.sublist("myMueluLISTDefault");
  Teuchos::ParameterList& muelu_fact = muelu.sublist("Factories");
  Teuchos::ParameterList& muelu_fact_gs = muelu_fact.sublist("myForwardGaussSeidel");
  muelu_fact_gs.set<std::string>("factory","TrilinosSmoother");
  muelu_fact_gs.set<std::string>("type"   ,"RELAXATION"      );
  muelu_fact_gs.sublist("ParameterList").set<std::string>("relaxation: type","Gauss-Seidel");

  Teuchos::ParameterList& muelu_hier = muelu.sublist("Hierarchy");
  muelu_hier.set             ("numDesiredLevel",3     );
  muelu_hier.set             ("maxCoarseSize"  ,700   );
  muelu_hier.set<std::string>("verbosity"      ,"None");
  muelu_hier.sublist("All").set             ("startLevel"  ,0                     );
  muelu_hier.sublist("All").set<std::string>("Smoother"    ,"myForwardGaussSeidel");
  muelu_hier.sublist("All").set<std::string>("CoarseSolver","myForwardGaussSeidel");




  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
LINALG::SOLVER::AMGnxn_Operator::AMGnxn_Operator(
    Teuchos::RCP<BlockSparseMatrixBase> A,
    std::vector<int> num_pdes,
    std::vector<int> null_spaces_dim,
    std::vector<Teuchos::RCP<std::vector<double> > > null_spaces_data,
    const Teuchos::ParameterList& amgnxn_params,
    const Teuchos::ParameterList& smoothers_params,
    const Teuchos::ParameterList& muelu_params):
    A_                   (A               ),
    num_pdes_            (num_pdes        ),
    null_spaces_dim_     (null_spaces_dim ),
    null_spaces_data_    (null_spaces_data),
    amgnxn_params_       (amgnxn_params   ),
    smoothers_params_    (smoothers_params),
    muelu_params_        (muelu_params),
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
  //  <Parameter name="muelu parameters for block 0"       type="string"  value="myMuelu0"/>
  //
  //  <Parameter name="muelu parameters for block 1"       type="string"  value="myMuelu1"/>
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

  //recover the muelu params
  int NumBlocks = A_->Rows();
  for(int i=0;i<NumBlocks;i++)
  {

    // recover the name of the list
    std::stringstream ss;
    ss << i;
    std::string param_name = "muelu parameters for block " + ss.str();
    std::string list_name =amgnxn_params_.get<std::string>(param_name,"none");
    if (list_name == "none")
      dserror("You must specify the parameters for creating the AMG on block %d",i);

    // Parse contents of the list
    Teuchos::ParameterList muelu_list_this_block;
    if (not muelu_params_.isSublist(list_name)) dserror("list %s not found",list_name.c_str());
    std::string xml_file = muelu_params_.sublist(list_name).get<std::string>("xml file","none");
    if (xml_file!="none")
      Teuchos::updateParametersFromXmlFile(
          xml_file,Teuchos::Ptr<Teuchos::ParameterList>(&muelu_list_this_block));
    else
      muelu_list_this_block=muelu_params_.sublist(list_name);

    muelu_lists_.push_back(muelu_list_this_block);
  }


  int num_levels_amg = amgnxn_params_.get<int>("number of levels",20);
  //if(num_levels_amg==-1)
  //  dserror("Missing \"number of levels\" in your xml file");
  H_ = Teuchos::rcp(new  AMGnxn_Hierarchies(
        A_,
        muelu_lists_,
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

