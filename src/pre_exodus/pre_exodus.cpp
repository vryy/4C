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
#include "pre_exodus_reader.H"
#include "pre_exodus_soshextrusion.H"

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
    
    // generate solid shell extrusion based on exodus file
    if (soshthickness!=0.0){
      if (exofile=="") dserror("no exofile specified for extrusion");
      if (soshthickness < 0.0) dserror("thickness specified for solid-shell extrusion is negative");
      //Soshextrusion mysosh(exofile.c_str(),soshthickness,soshnumlayer);
      string extrudefile;
      extrudefile = "extr_" + exofile;
      Mesh mysosh(exofile.c_str());
      mysosh.Print(cout,true);
      //mysosh.WriteMesh(extrudefile);
      exit(1);
    }
    
    // create mesh object based on given exodus II file
    Mesh mymesh(exofile.c_str());
    // print infos to cout
    mymesh.Print(cout);
    //mymesh.CloseExo();

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
      mymesh.Print(defaultbc, true);

      // show all necessary bc and element specifications (suggestion for user)     
      defaultbc<<"-----------Please specify:-----------"<<endl<<endl;
	    defaultbc<<"\\\\SYNTAX:"<<endl;
      defaultbc<<"\\\\"<<"*matrID"<<"=\"<DESIGN TYPE or ELEMENT>\"\"<DESIGN-number or ELEMENT-TYPE>\""<<endl
      <<"\\\\"<<"boundr_cond=\"<BC TYPE>\"\"<the BC itself>\""<<endl
      <<"\\\\"<<"type=\"<ELEMENT-type>\"\"<ELEMENT-type2>\"\"<ELEMENT-type3>\""<<endl
      <<"\\\\"<<"prop=\"<element flags>\""<<endl<<endl;
      
      defaultbc<<"\\\\ For Example:"<<endl
      <<"\\\\*matr1 =\"ELEMENT\"\"STRUCTURE\""<<endl
      <<"\\\\boundr_cond=\" \"\" \" "<<endl
      <<"\\\\type=\"STRUCTURE\"\"SOLIDSH8\"\"HEX8\" "<<endl
      <<"\\\\prop=\"MAT       1 GP 2 2 2 STRESS None KINEM Totlag  EAS sosh8 THICKDIR AUTO\" "<<endl<<endl;
      defaultbc<<"\\\\*matr1 =\"ELEMENT\"\"FLUID\" "<<endl
      <<"\\\\boundr_cond=\" \"\" \" "<<endl
      <<"\\\\type=\"FLUID\"\"FLUID3\"\"HEX8\" " <<endl
      <<"\\\\prop=\"MAT       1 NA ALE GP 2 2 2\" "<<endl<<endl;
      defaultbc<<"\\\\*matr2 =\"DSURF\"\"1\" "<<endl
      << "\\\\boundr_cond=\"NEUMANN\"\"E 1 - 1  0 0 1 0 0 0 0.0 0.0 1.0 0.0 0.0 0.0 Live Mid\" "<<endl<<endl;
      defaultbc<<"\\\\*matr6 =\"DSURF\"\"4\" "<<endl
      <<"\\\\boundr_cond=\"DIRICH\"\"E 4 - 1 1 1 0 0 0 0.0 0.0 0.0 0.0 0.0 0.0 none none none none none none 0 0 0 0 0 0\" "<<endl<<endl;

      defaultbc<<"------------------------------------------------"<<endl<<endl;
      
      int numEntities = mymesh.GetNumElementBlocks() + mymesh.GetNumNodeSets() + mymesh.GetNumSideSets();
      for (int i = 0; i < numEntities; ++i) 
      {    
    	  defaultbc<<"*matr"<<i+1<<"=\"\"\"\""<<endl
    	  <<"boundr_cond=\"\"\"\""<<endl
    	  <<"type=\"\"\"\"\"\""<<endl
    	  <<"prop=\"\""<<endl<<endl;
      }

      // close default bc specification file
      if (defaultbc.is_open()) 
    	  defaultbc.close();
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
      Teuchos::RCP<const Teuchos::ParameterList> list = DRT::ValidParameters();

      // write default .dat header into file 
      DRT::PrintDatHeader(defaulthead,*list);

      // close default header file    
      if (defaulthead.is_open()) 
    	  defaulthead.close();   
    }
    
#ifdef PARALLEL
  MPI_Finalize();
#endif
}

#endif
#endif
