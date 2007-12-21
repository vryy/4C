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
#ifdef EXODUS
#ifdef CCADISCRET
#include "pre_exodus.H"
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include "../drt_lib/drt_validparameters.H"
#include "pre_exodus_reader.H"

using namespace std;
using namespace Teuchos;

int main(
        int argc,
        char** argv)
{
    string exofile;
    string bcfile;
    string headfile;
    string datfile;

    Teuchos::CommandLineProcessor My_CLP;
    My_CLP.setDocString("Preprocessor Exodus2Baci \n");
    My_CLP.throwExceptions(false);
    My_CLP.setOption("exo",&exofile,"exodus file to open");
    My_CLP.setOption("bc",&bcfile,"bc's and ele's file to open");
    My_CLP.setOption("head",&headfile,"baci header file to open");
    My_CLP.setOption("dat",&datfile,"output .dat file name [defaults to exodus file name]");
    
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
    Mesh mymesh(exofile.c_str());
    // print infos to cout
    mymesh.Print(cout);

    if (bcfile=="")
    {
      string defaultbcfilename = "default.bc";
      cout << "found no BC specification file --> creating " <<defaultbcfilename<< endl;

      // open default bc specification file
      ofstream defaultbc(defaultbcfilename.c_str());
      if (!defaultbc)
    	  dserror("failed to open file: %s", defaultbcfilename.c_str());

      // write general mesh info
      mymesh.Print(defaultbc);
      // write infos of each mesh entity
      for (int i = 0; i < mymesh.GetNumberEntities(); ++i) 
      {
        RCP<Entity> actEntity = mymesh.GetEntity(i);
        actEntity->Print(defaultbc);
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

      exit(0);    
    }
}

#endif
#endif
