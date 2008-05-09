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
specifications. As result either a preliminary input file set is suggestioned,
or the well-known .dat file is created.

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

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#endif
#include <Epetra_SerialComm.h>
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_inputreader.H"


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

    if (exofile=="")
    {
      cout<<"No Exodus II file was found"<<endl;
      My_CLP.printHelpMessage(argv[0],cout);
      exit(1);
    }
    
    
    // create mesh object based on given exodus II file
    EXODUS::Mesh mymesh(exofile.c_str());
    // print infos to cout
    mymesh.Print(cout);
    //mymesh.CloseExo();

    // declare empty vectors for holding "boundary" conditions
    vector<EXODUS::elem_def> eledefs;
    vector<EXODUS::cond_def> condefs;
    
    // generate solid shell extrusion based on exodus file
    if (soshthickness!=0.0){
      if (exofile=="") dserror("no exofile specified for extrusion");
      //if (soshthickness < 0.0) dserror("thickness specified for solid-shell extrusion is negative");
      EXODUS::Mesh mysosh = EXODUS::SolidShellExtrusion(mymesh, soshthickness, soshnumlayer, soshseedid, soshgmsh);
      string extrudefile;
      extrudefile = "extr_" + exofile;
      //Mesh mysosh(exofile.c_str());
      //mysosh.Print(cout,true);
      mysosh.WriteMesh(extrudefile);
      exit(1);
    }

    mymesh.CloseExo();
    if (bcfile=="")
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

      /*
      // show all necessary bc and element specifications (suggestion for user)     
      defaultbc<<"-----------Please specify:-----------"<<endl<<endl;
	    defaultbc<<"\\\\SYNTAX:"<<endl;
      defaultbc<<"\\\\"<<"matrID"<<"=\"<DESIGN TYPE or ELEMENT>\"\"<DESIGN-number or ELEMENT-TYPE>\""<<endl
      <<"\\\\"<<"boundr_cond=\"<BC TYPE>\"\"<the BC itself>\""<<endl
      <<"\\\\"<<"type=\"<ELEMENT-type>\"\"<ELEMENT-type2>\"\"<ELEMENT-type3>\""<<endl
      <<"\\\\"<<"prop=\"<element flags>\""<<endl<<endl;
      
      defaultbc<<"\\\\ For Example:"<<endl
      <<"\\\\matr1 =\"ELEMENT\"\"STRUCTURE\""<<endl
      <<"\\\\boundr_cond=\" \"\" \" "<<endl
      <<"\\\\type=\"STRUCTURE\"\"SOLIDSH8\"\"HEX8\" "<<endl
      <<"\\\\prop=\"MAT       1 EAS sosh8 THICKDIR auto\" "<<endl<<endl;
      defaultbc<<"\\\\matr1 =\"ELEMENT\"\"FLUID\" "<<endl
      <<"\\\\boundr_cond=\" \"\" \" "<<endl
      <<"\\\\type=\"FLUID\"\"FLUID3\"\"HEX8\" " <<endl
      <<"\\\\prop=\"MAT       1 NA ALE GP 2 2 2\" "<<endl<<endl;
      defaultbc<<"\\\\matr2 =\"DSURF\"\"1\" "<<endl
      << "\\\\boundr_cond=\"NEUMANN\"\"E 1 - 1  0 0 1 0 0 0 0.0 0.0 1.0 0.0 0.0 0.0 Live Mid\" "<<endl<<endl;
      defaultbc<<"\\\\matr6 =\"DSURF\"\"4\" "<<endl
      <<"\\\\boundr_cond=\"DIRICH\"\"E 4 - 1 1 1 0 0 0 0.0 0.0 0.0 0.0 0.0 0.0 none none none none none none 0 0 0 0 0 0\" "<<endl<<endl;
      */
      
      defaultbc << "MIND that you can specify a condittion also on an ElementBlock, just replace 'ELEMENT' with 'CONDITION'"<<endl;
      defaultbc<<"------------------------------------------------BCSPECS"<<endl<<endl;
      
      
      // write ElementBlocks with specification proposal
      const map<int,EXODUS::ElementBlock> myblocks = mymesh.GetElementBlocks();
      map<int,EXODUS::ElementBlock>::const_iterator it;
      for (it = myblocks.begin(); it != myblocks.end(); ++it){
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
      
      //(*condlist)[0]->Print(cout,NULL,true);


      
      /*
      int numEntities =  + mymesh.GetNumNodeSets() + mymesh.GetNumSideSets();
      for (int i = 0; i < numEntities; ++i) 
      {    
          defaultbc<<"*matr"<<i+1<<"=\"\"\"\""<<endl
          <<"boundr_cond=\"\"\"\""<<endl
          <<"type=\"\"\"\"\"\""<<endl
          <<"prop=\"\""<<endl<<endl;
      }
      */

      // close default bc specification file
      if (defaultbc.is_open()) 
        defaultbc.close();
    }
    else
    {
      // read provided bc-file
      EXODUS::ReadBCFile(bcfile,eledefs,condefs);
      
//      if ((signed) conditions[0].size() != mymesh.GetNumElementBlocks()){
//        cout << "You specified only " << conditions[0].size() << " ElementBlocks of " << mymesh.GetNumElementBlocks() << " in the provided mesh!" << endl;
//      } else if ((signed) conditions[1].size() != mymesh.GetNumNodeSets()){
//        cout << "You specified only " << conditions[1].size() << " NodeSets of " << mymesh.GetNumNodeSets() << " in the provided mesh!" << endl;
//      } else if ((signed) conditions[2].size() != mymesh.GetNumSideSets()){
//        cout << "You specified only " << conditions[2].size() << " SideSets of " << mymesh.GetNumSideSets() << " in the provided mesh!" << endl;
//      }
    }

    if (headfile=="")
    {
      string defaultheadfilename = "default.head";
      cout << "found no header file           --> creating "<<defaultheadfilename<< endl;
      // open default header file
      ofstream defaulthead(defaultheadfilename.c_str());
      if (!defaulthead)
        dserror("failed to open file: %s", defaultheadfilename.c_str());

      // get valid input parameters
      Teuchos::RCP<const Teuchos::ParameterList> list = DRT::INPUT::ValidParameters();

      // write default .dat header into file 
      DRT::INPUT::PrintDatHeader(defaulthead,*list);

      // add additional line at the end of the default header file
      // (needed in order to tell the DatFileReader that the last section has finished)
      defaulthead<<"-----------------------------------------------------EOF HEADERFILE"<<endl;

      // close default header file
      if (defaulthead.is_open())
         defaulthead.close();
    }
    else
    {
      // read and check the provided header file
      cout << "checking given header file "<<headfile<< endl;

#ifdef PARALLEL
      int myrank = 0;
      int nproc  = 1;
      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
      MPI_Comm_size(MPI_COMM_WORLD, &nproc);
      Epetra_MpiComm* com = new Epetra_MpiComm(MPI_COMM_WORLD);
      RefCountPtr<Epetra_Comm> comm = rcp(com);
#else
        Epetra_SerialComm* com = new Epetra_SerialComm();
        RefCountPtr<Epetra_Comm> comm = rcp(com);
#endif

      Teuchos::RCP<DRT::Problem> problem = DRT::Problem::Instance();
      DRT::INPUT::DatFileReader reader(headfile.c_str(), comm, 0, false);

      // do reading AND validation!
      problem->ReadParameter(reader);

      //the header file is valid --> go on

      // clean up
      problem->Done();
    }

#ifdef PARALLEL
  MPI_Finalize();
#endif

  if (datfile=="")
  {
    // default dat file, later
  }
  else
  {
    cout << "creating datfile " << datfile << endl;
    EXODUS::WriteDatFile(datfile, mymesh, headfile, eledefs, condefs);
  }

  
  return 0;

}

#endif
#endif
