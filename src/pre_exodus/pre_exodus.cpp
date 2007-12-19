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
      exit(1);
    }

    
    if (exofile=="")
    {
      My_CLP.printHelpMessage(argv[0],cout);
      exit(1);
    }
    else
    {
      mesh(exofile.c_str());
    }
      
   
    if (headfile=="")
    {
      cout << "Headfile leer, jetzt validparameters" << endl;
      string defaultheadfilename = "default.head";
      ofstream defaulthead(defaultheadfilename.c_str());
      if (!defaulthead);
          //dserror("failed to open file: %s", defaultheadfilename.c_str());
      
      Teuchos::RCP<const Teuchos::ParameterList> list = DRT::ValidParameters();
      //cout << *list << endl;
      //PrintDefaultDatHeader();
      DRT::PrintDatHeader(defaulthead,*list);
      if (defaulthead.is_open()) defaulthead.close();

      exit(0);
      
    }
    
    

    
    cout << "Hello World" << endl;
    
}



#endif
#endif