/*----------------------------------------------------------------------*/
/*!
\file pre_exodus.cpp

\brief preprocessor for exodusII format 

<pre>
Maintainer: Moritz & Georg
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089 - 289-15240
</pre>

Pre_exodus contains classes to open and preprocess exodusII files into the
drt of Baci. It uses the "valid-parameters"-list defined in Baci for preparing
a up-to-date Baci header and another file specifying element and boundary
specifications based on "valid-conditions". As result either a preliminary
input file set is suggestioned, or the well-known .dat file is created.
Addionally, specify an already existing BACI input file in order to validate
its parameters and conditions.

*/
/*----------------------------------------------------------------------*/
#ifdef D_EXODUS
#ifdef CCADISCRET
#include "pre_exodus.H"
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include "../drt_lib/drt_validparameters.H"
#include "../drt_lib/drt_validconditions.H"
#include "../drt_lib/drt_conditiondefinition.H"
#include "pre_exodus_reader.H"
#include "pre_exodus_soshextrusion.H"
#include "pre_exodus_writedat.H"
#include "pre_exodus_readbc.H"
#include "pre_exodus_validate.H"


using namespace std;
using namespace Teuchos;

int main(
        int argc,
        char** argv)
{

#ifdef PARALLEL
  MPI_Init(&argc,&argv);
#endif

    string exofile;
    string bcfile;
    string headfile;
    string datfile;
    
    // related to solid shell extrusion
    double soshthickness = 0.0;
    int soshnumlayer = 1;
    int soshseedid = 0;
    int soshgmsh = -1;

    Teuchos::CommandLineProcessor My_CLP;
    My_CLP.setDocString("Preprocessor Exodus2Baci \n");
    My_CLP.throwExceptions(false);
    My_CLP.setOption("exo",&exofile,"exodus file to open");
    My_CLP.setOption("bc",&bcfile,"bc's and ele's file to open");
    My_CLP.setOption("head",&headfile,"baci header file to open");
    My_CLP.setOption("dat",&datfile,"output .dat file name [defaults to exodus file name]");
    // here options related to solid shell extrusion are defined
    My_CLP.setOption("gensosh",&soshthickness,"generate solid-shell body with given thickness");
    My_CLP.setOption("numlayer",&soshnumlayer,"number of layers of generated solid-shell body");
    My_CLP.setOption("seedid",&soshseedid,"id where to start extrusion, default is first");
    My_CLP.setOption("gmsh",&soshgmsh,"gmsh output of xxx elements, default off, 0 all eles");
    
    CommandLineProcessor::EParseCommandLineReturn
      parseReturn = My_CLP.parse(argc,argv);

    if (parseReturn == CommandLineProcessor::PARSE_HELP_PRINTED)
    {
      exit(0);
    }
    if (parseReturn != CommandLineProcessor::PARSE_SUCCESSFUL)
    {
      dserror("CommandLineProcessor reported an error");
    }

    /**************************************************************************
     * Start with the preprocessing
     **************************************************************************/
    if (exofile=="")
    {
      if (datfile!="")
      {
        // just validate a given BACI input file
        EXODUS::ValidateInputFile(datfile);
        return 0;
      }
      else
      {
        cout<<"No Exodus II file was found"<<endl;
        My_CLP.printHelpMessage(argv[0],cout);
        exit(1);
      }
    }

    // create mesh object based on given exodus II file
    EXODUS::Mesh mymesh(exofile.c_str());
    // print infos to cout
    mymesh.Print(cout);

    /**************************************************************************
     * Edit a existing Mesh, e.g. extrusion of surface
     **************************************************************************/

    // generate solid shell extrusion based on exodus file
    if (soshthickness!=0.0){
      if (exofile=="") dserror("no exofile specified for extrusion");
      EXODUS::Mesh mysosh = EXODUS::SolidShellExtrusion(mymesh, soshthickness, soshnumlayer, soshseedid, soshgmsh);
      mysosh.WriteMesh("extr_" + exofile);
      exit(0);
    }

    /**************************************************************************
     * Read ControlFile for Boundary and Element descriptions
     **************************************************************************/

    // declare empty vectors for holding "boundary" conditions
    vector<EXODUS::elem_def> eledefs;
    vector<EXODUS::cond_def> condefs;

    if (bcfile=="")
    {
      int error = EXODUS::CreateDefaultBCFile(mymesh);
      if (error!=0) dserror("Creation of default bc-file not successful.");
    }
    else
    {
      // read provided bc-file
      EXODUS::ReadBCFile(bcfile,eledefs,condefs);
      
      int sum = mymesh.GetNumElementBlocks() + mymesh.GetNumNodeSets() + mymesh.GetNumSideSets();
      int test = eledefs.size() + condefs.size();
      if (test != sum) cout << "Your " << test << " definitions do not match the " << sum << " entities in your mesh!" <<endl; 
    }

    /**************************************************************************
     * Read HeaderFile for 'header' parameters, e.g. solver, dynamic, material
     **************************************************************************/
    if (headfile=="")
    {
      string defaultheadfilename = "default.head";
      cout << "found no header file           --> creating "<<defaultheadfilename<< endl;

      // open default header file
      ofstream defaulthead(defaultheadfilename.c_str());
      if (!defaulthead) dserror("failed to open file: %s", defaultheadfilename.c_str());

      // get valid input parameters
      Teuchos::RCP<const Teuchos::ParameterList> list = DRT::INPUT::ValidParameters();

      // write default .dat header into file 
      DRT::INPUT::PrintDatHeader(defaulthead,*list);

      // close default header file
      if (defaulthead.is_open()) defaulthead.close();
    }

    /**************************************************************************
     * Finally, create and validate the BACI input file
     **************************************************************************/
    if ((headfile!="") && (bcfile!="") && (exofile!=""))
    {
      // set default dat-file name if needed
      if (datfile=="")
      {
        string exofilebasename = exofile.substr(0,exofile.find_first_of("."));
        datfile=exofilebasename+".dat";
      }

      // write the BACI input file
      cout << "creating BACI input file       --> " << datfile << endl;
      EXODUS::WriteDatFile(datfile, mymesh, headfile, eledefs, condefs);

      //validate the generated BACI input file
      EXODUS::ValidateInputFile(datfile);
    }

#ifdef PARALLEL
  MPI_Finalize();
#endif
  
  return 0;

} //main.cpp


// create default bc file
int EXODUS::CreateDefaultBCFile(EXODUS::Mesh& mymesh)
{
  string defaultbcfilename = "default.bc";
  cout << "found no BC specification file --> creating " <<defaultbcfilename<< endl;

  // open default bc specification file
  ofstream defaultbc(defaultbcfilename.c_str());
  if (!defaultbc)
    dserror("failed to open file: %s", defaultbcfilename.c_str());

  // write mesh verbosely
  defaultbc<<"----------- Mesh contents -----------"<<endl<<endl;
  mymesh.Print(defaultbc, false);
  
  // give examples for element and boundary condition syntax
  defaultbc<<"---------- Syntax examples ----------"<<endl<<endl<<
  "Element Block, named: "<<endl<<
  "of Shape: TET4"<<endl<<
  "has 9417816 Elements"<<endl<<
  "\"*eb0=\"ELEMENT\""<<endl<<
  "sectionname=\"FLUID\""<<endl<<
  "description=\"MAT 1 NA Euler GP 2 2 2\""<<endl<<
  "elementname=\"FLUID3\" \n"<<endl<<
  "Node Set, named:"<<endl<<
  "Property Name: INFLOW"<<endl<<
  "has 45107 Nodes"<<endl<<
  "\"*ns0=\"CONDITION\""<<endl<<
  "sectionname=\"DESIGN SURF DIRICH CONDITIONS\""<<endl<<
  "description=\"E 1 - 1 1 1 0 0 0 2.0 0.0 0.0 0.0 0.0 0.0  1 none none none none none  1 0 0 0 0 0\""
  <<endl<<endl;

  defaultbc << "MIND that you can specify a condition also on an ElementBlock, just replace 'ELEMENT' with 'CONDITION'"<<endl;
  defaultbc<<"------------------------------------------------BCSPECS"<<endl<<endl;

  // write ElementBlocks with specification proposal
  RCP<const map<int,EXODUS::ElementBlock> > myblocks = mymesh.GetElementBlocks();
  map<int,EXODUS::ElementBlock>::const_iterator it;
  for (it = myblocks->begin(); it != myblocks->end(); ++it){
    it->second.Print(defaultbc);
    defaultbc<<"*eb"<< it->first << "=\"ELEMENT\""<<endl
    <<"sectionname=\"\""<<endl
    <<"description=\"\""<<endl
    <<"elementname=\"\""<<endl
    //<<"elementshape=\""
    //<< DRT::DistypeToString(PreShapeToDrt(it->second.GetShape()))<<"\""<<endl
    <<endl;
  }

  // write NodeSets with specification proposal
  const map<int,EXODUS::NodeSet> mynodesets = mymesh.GetNodeSets();
  map<int,EXODUS::NodeSet>::const_iterator ins;
  for (ins =mynodesets.begin(); ins != mynodesets.end(); ++ins){
    ins->second.Print(defaultbc);
    defaultbc<<"*ns"<< ins->first << "=\"CONDITION\""<<endl
    <<"sectionname=\"\""<<endl
    <<"description=\"\""<<endl
    <<endl;
  }

  // write SideSets with specification proposal
  const map<int,EXODUS::SideSet> mysidesets = mymesh.GetSideSets();
  map<int,EXODUS::SideSet>::const_iterator iss;
  for (iss = mysidesets.begin(); iss!=mysidesets.end(); ++iss){
    iss->second.Print(defaultbc);
    defaultbc<<"*ss"<< iss->first << "=\"CONDITION\""<<endl
    <<"sectionname=\"\""<<endl
    <<"description=\"\""<<endl
    <<endl;
  }

  // print validconditions as proposal
  defaultbc << "-----------------------------------------VALIDCONDITIONS"<< endl;
  Teuchos::RCP<std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> > > condlist = DRT::INPUT::ValidConditions();
  DRT::INPUT::PrintEmptyConditionDefinitions(defaultbc,*condlist,false);

  // close default bc specification file
  if (defaultbc.is_open()) 
    defaultbc.close();

  return 0;
}


#endif
#endif
