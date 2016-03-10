/*!----------------------------------------------------------------------
\file statmech_output.cpp
\brief output methods for statistical mechanics

<pre>
Maintainer: Kei Müller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>

*----------------------------------------------------------------------*/

#include "statmech_manager.H"

#include <Teuchos_Time.hpp>

#include "../drt_inpar/inpar_statmech.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_beamcontact/beam3contact_manager.H"
#include "../drt_beamcontact/beam3contact_octtree.H"

#include "../drt_beam3/beam3.H"
#include "../drt_beam3r/beam3r.H"
#include "../drt_beam3eb/beam3eb.H"
#include "../drt_truss3/truss3.H"
#include "../drt_truss3cl/truss3cl.H"
#include "../drt_beam3cl/beam3cl.H"
#include "../drt_torsion3/torsion3.H"

//MEASURETIME activates measurement of computation time for certain parts of the code
//#define MEASURETIME

void STATMECH::StatMechManager::BuildStatMechRootPath()
{
    // create file name and check existence of the required output folder
    outputrootpath_ = DRT::Problem::Instance()->OutputControlFile()->FileName();
    size_t pos = outputrootpath_.rfind('/');
    outputrootpath_ = outputrootpath_.substr(0,pos);
    // replace last folder by new pattern
    std::string::iterator it = outputrootpath_.end();
    while(it!=outputrootpath_.begin())
    {
      if(*it=='/')
        break;
      it--;
    }
    if(it==outputrootpath_.begin())
      outputrootpath_.replace(it,outputrootpath_.end(),".");
    else
      outputrootpath_.replace(it,outputrootpath_.end(),"");

    // Check for existence of the folder StatMechOutput
    std::ostringstream statmechfilepath;
    statmechfilepath << outputrootpath_ << "/StatMechOutput/";
    struct stat st;
    if(stat(statmechfilepath.str().c_str(), &st) !=0 && DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")!=INPAR::STATMECH::statout_none)
      dserror("The folder %s was not found but is required for statistical mechanics output!", statmechfilepath.str().c_str());
    std::ostringstream gmshfilepath;
    gmshfilepath << outputrootpath_ <<"/GmshOutput/";
    if(stat(gmshfilepath.str().c_str(), &st) !=0 && DRT::INPUT::IntegralValue<int>(statmechparams_, "GMSHOUTPUT"))
      dserror("The folder %s was not found but is required for Gmsh output!", gmshfilepath.str().c_str());
  return;
}

/*----------------------------------------------------------------------*
 | initialize special output for statistical mechanics(public)cyron 12/08|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::InitOutput(const int& ndim,
                                           const Epetra_Vector& dis,
                                           const int& istep,
                                           const double& dt)
{
  //initializing special output for statistical mechanics by looking for a suitable name of the outputfile and setting up an empty file with this name
  switch (DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT"))
  {
    case INPAR::STATMECH::statout_endtoendlog:
    {

      //output is written on proc 0 only
      if(!discret_->Comm().MyPID())
      {

        FILE* fp = NULL; //file pointer for statistical output file

        //defining name of output file
        std::ostringstream outputfilename;

        outputfilenumber_ = 0;

        //file pointer for operating with numbering file
        FILE* fpnumbering = NULL;
        std::ostringstream numberingfilename;

        //look for a numbering file where number of already existing output files is stored:
        numberingfilename.str("NumberOfRealizationsLog");
        fpnumbering = fopen(numberingfilename.str().c_str(), "r");

        //if there is no such numbering file: look for a not yet existing output file name (numbering upwards)
        if (fpnumbering == NULL)
        {
          do
          {
            outputfilenumber_++;
            outputfilename.str("");
            outputfilename << outputrootpath_ << "/StatMechOutput/EndToEnd" << outputfilenumber_ << ".dat";
            fp = fopen(outputfilename.str().c_str(), "r");
          }
          while (fp != NULL);

          //set up new file with name "outputfilename" without writing anything into this file
          fp = fopen(outputfilename.str().c_str(), "w");
          fclose(fp);
        }
        //if there already exists a numbering file
        else
        {
          fclose(fpnumbering);

          //read the number of the next realization out of the file into the variable testnumber
          std::fstream f(numberingfilename.str().c_str());
          while (f)
          {
            std::string tok;
            f >> tok;
            if (tok == "Next")
            {
              f >> tok;
              if (tok == "Number:")
                f >> outputfilenumber_;
            }
          } //while(f)

          //defining outputfilename by means of new testnumber
          outputfilename.str("");
          outputfilename << outputrootpath_ << "/StatMechOutput/EndToEnd" << outputfilenumber_ << ".dat";
        }

        //increasing the number in the numbering file by one
        fpnumbering = fopen(numberingfilename.str().c_str(), "w");
        std::stringstream filecontent;
        filecontent << "Next Number: " << (outputfilenumber_ + 1);
        fputs(filecontent.str().c_str(), fpnumbering);
        fclose(fpnumbering);
      }

    }
      break;
    case INPAR::STATMECH::statout_endtoendconst:
    {
      //output is written on proc 0 only
      if(!discret_->Comm().MyPID())
      {

        FILE* fp = NULL; //file pointer for statistical output file

        //defining name of output file
        std::ostringstream outputfilename;
        outputfilename.str("");
        outputfilenumber_ = 0;

        //file pointer for operating with numbering file
        FILE* fpnumbering = NULL;
        std::ostringstream numberingfilename;

        //look for a numbering file where number of already existing output files is stored:
        numberingfilename.str("NumberOfRealizationsConst");
        fpnumbering = fopen(numberingfilename.str().c_str(), "r");

        /*in the following we assume that there is only a pulling force point Neumann condition of equal absolute
         *value on either filament length; we get the absolute value of the first one of these two conditions */
        double neumannforce;
        std::vector<DRT::Condition*> pointneumannconditions(0);
        discret_->GetCondition("PointNeumann", pointneumannconditions);
        if (pointneumannconditions.size() > 0)
        {
          const std::vector<double>* val = pointneumannconditions[0]->Get<std::vector<double> > ("val");
          neumannforce = fabs((*val)[0]);
        }
        else
          neumannforce = 0;

        //if there is no such numbering file: look for a not yet existing output file name (numbering upwards)
        if (fpnumbering == NULL)
        {
          do
          {
            //defining name of output file
            outputfilenumber_++;
            outputfilename.str("");
            outputfilename << outputrootpath_ << "/StatMechOutput/E2E_" << discret_->NumMyRowElements() << '_' << dt<< '_' << neumannforce << '_' << outputfilenumber_ << ".dat";
            fp = fopen(outputfilename.str().c_str(), "r");
          }
          while (fp != NULL);

          //set up new file with name "outputfilename" without writing anything into this file
          fp = fopen(outputfilename.str().c_str(), "w");
          fclose(fp);
        }
        //if there already exists a numbering file
        else
        {
          fclose(fpnumbering);

          //read the number of the next realization out of the file into the variable testnumber
          std::fstream f(numberingfilename.str().c_str());
          while (f)
          {
            std::string tok;
            f >> tok;
            if (tok == "Next")
            {
              f >> tok;
              if (tok == "Number:")
                f >> outputfilenumber_;
            }
          } //while(f)

          //defining outputfilename by means of new testnumber
          outputfilename.str("");
          outputfilename << outputrootpath_ << "/StatMechOutput/E2E_" << discret_->NumMyRowElements() << '_' << dt << '_' << neumannforce << '_' << outputfilenumber_ << ".dat";

          //set up new file with name "outputfilename" without writing anything into this file
          fp = fopen(outputfilename.str().c_str(), "w");
          fclose(fp);
        }

        //increasing the number in the numbering file by one
        fpnumbering = fopen(numberingfilename.str().c_str(), "w");
        std::stringstream filecontent;
        filecontent << "Next Number: " << (outputfilenumber_ + 1);
        //write fileconent into file!
        fputs(filecontent.str().c_str(), fpnumbering);
        //close file
        fclose(fpnumbering);
      }

    }
    break;
    /*computing and writing into file data about correlation of orientation of different elements
     *as it is considered in the context of the persistence length*/
    case INPAR::STATMECH::statout_orientationcorrelation:
    {
      //output is written on proc 0 only
      if(!discret_->Comm().MyPID())
      {

        FILE* fp = NULL; //file pointer for statistical output file

        //defining name of output file
        std::ostringstream outputfilename;

        outputfilenumber_ = 0;

        //file pointer for operating with numbering file
        FILE* fpnumbering = NULL;
        std::ostringstream numberingfilename;

        //look for a numbering file where number of already existing output files is stored:
        numberingfilename.str("NumberOfRealizationsOrientCorr");
        fpnumbering = fopen(numberingfilename.str().c_str(), "r");

        //if there is no such numbering file: look for a not yet existing output file name (numbering upwards)
        if (fpnumbering == NULL)
        {
          do
          {
            outputfilenumber_++;
            outputfilename.str("");
            outputfilename << outputrootpath_ << "/StatMechOutput/OrientationCorrelation" << outputfilenumber_<< ".dat";
            fp = fopen(outputfilename.str().c_str(), "r");
          }
          while (fp != NULL);

          //set up new file with name "outputfilename" without writing anything into this file
          fp = fopen(outputfilename.str().c_str(), "w");
          fclose(fp);
        }
        //if there already exists a numbering file
        else
        {
          fclose(fpnumbering);

          //read the number of the next realization out of the file into the variable testnumber
          std::fstream f(numberingfilename.str().c_str());
          while (f)
          {
            std::string tok;
            f >> tok;
            if (tok == "Next")
            {
              f >> tok;
              if (tok == "Number:")
                f >> outputfilenumber_;
            }
          } //while(f)

          //defining outputfilename by means of new testnumber
          outputfilename.str("");
          outputfilename << outputrootpath_ << "/StatMechOutput/OrientationCorrelation" << outputfilenumber_<< ".dat";
        }

        //set up new file with name "outputfilename" without writing anything into this file
        fp = fopen(outputfilename.str().c_str(), "w");
        fclose(fp);

        //increasing the number in the numbering file by one
        fpnumbering = fopen(numberingfilename.str().c_str(), "w");
        std::stringstream filecontent;
        filecontent << "Next Number: " << (outputfilenumber_ + 1);
        //write fileconent into file!
        fputs(filecontent.str().c_str(), fpnumbering);
        //close file
        fclose(fpnumbering);
      }
    }
    break;
    //simulating diffusion coefficient for anisotropic friction
    case INPAR::STATMECH::statout_anisotropic:
    {
      #ifndef CHOSENOUTPUTNODE
        //output is written on proc 0 only
        if(!discret_->Comm().MyPID())
        {

          FILE* fp = NULL; //file pointer for statistical output file

          //defining name of output file
          std::ostringstream outputfilename;

          outputfilenumber_ = 0;

          //file pointer for operating with numbering file
          FILE* fpnumbering = NULL;
          std::ostringstream numberingfilename;

          //look for a numbering file where number of already existing output files is stored:
          numberingfilename.str("NumberOfRealizationsAniso");
          fpnumbering = fopen(numberingfilename.str().c_str(), "r");

          //if there is no such numbering file: look for a not yet existing output file name (numbering upwards)
          if (fpnumbering == NULL)
          {
            do
            {
              outputfilenumber_++;
              outputfilename.str("");
              outputfilename << outputrootpath_ << "/StatMechOutput/AnisotropicDiffusion" << outputfilenumber_ << ".dat";
              fp = fopen(outputfilename.str().c_str(), "r");
            }
            while (fp != NULL);

            //set up new file with name "outputfilename" without writing anything into this file
            fp = fopen(outputfilename.str().c_str(), "w");
            fclose(fp);
          }
          //if there already exists a numbering file
          else
          {
            fclose(fpnumbering);

            //read the number of the next realization out of the file into the variable testnumber
            std::fstream f(numberingfilename.str().c_str());
            while (f)
            {
              std::string tok;
              f >> tok;
              if (tok == "Next")
              {
                f >> tok;
                if (tok == "Number:")
                  f >> outputfilenumber_;
              }
            } //while(f)

            //defining outputfilename by means of new testnumber
            outputfilename.str("");
            outputfilename << outputrootpath_ << "/StatMechOutput/AnisotropicDiffusion" << outputfilenumber_ << ".dat";
          }

          //set up new file with name "outputfilename" without writing anything into this file
          fp = fopen(outputfilename.str().c_str(), "w");
          fclose(fp);

          //increasing the number in the numbering file by one
          fpnumbering = fopen(numberingfilename.str().c_str(), "w");
          std::stringstream filecontent;
          filecontent << "Next Number: " << (outputfilenumber_ + 1);
          //write fileconent into file!
          fputs(filecontent.str().c_str(), fpnumbering);
          //close file
          fclose(fpnumbering);

          //initializing variables for positions of first and last node at the beginning
          beginold_.PutScalar(0);
          endold_.PutScalar(0);

          for (int i=0; i<ndim; i++)
          {
            beginold_(i) = (discret_->gNode(0))->X()[i];
            endold_(i) = (discret_->gNode(discret_->NumMyRowNodes() - 1))->X()[i];
          }

          for (int i=0; i<3; i++)
            sumdispmiddle_(i, 0) = 0.0;

          sumsquareincpar_ = 0.0;
          sumsquareincort_ = 0.0;
          sumrotmiddle_ = 0.0;
          sumsquareincmid_ = 0.0;
          sumsquareincrot_ = 0.0;
        }
      #else
        //output is written on proc 0 only
        if(!discret_->Comm().MyPID())
        {

          FILE* fp = NULL; //file pointer for statistical output file

          //defining name of output file
          std::ostringstream outputfilename;

          outputfilenumber_ = 0;

          //file pointer for operating with numbering file
          FILE* fpnumbering = NULL;
          std::ostringstream numberingfilename;

          //look for a numbering file where number of already existing output files is stored:
          numberingfilename.str("NumberOfRealizationsAniso");
          fpnumbering = fopen(numberingfilename.str().c_str(), "r");

          //if there is no such numbering file: look for a not yet existing output file name (numbering upwards)
          if (fpnumbering == NULL)
          {
            do
            {
              outputfilenumber_++;
              outputfilename.str("");
              outputfilename << outputrootpath_ << "/StatMechOutput/NodalAverageDisp" << outputfilenumber_ << ".dat";
              fp = fopen(outputfilename.str().c_str(), "r");
            }
            while (fp != NULL);

            //set up new file with name "outputfilename" without writing anything into this file
            fp = fopen(outputfilename.str().c_str(), "w");
            fclose(fp);
          }
          //if there already exists a numbering file
          else
          {
            fclose(fpnumbering);

            //read the number of the next realization out of the file into the variable testnumber
            std::fstream f(numberingfilename.str().c_str());
            while (f)
            {
              std::string tok;
              f >> tok;
              if (tok == "Next")
              {
                f >> tok;
                if (tok == "Number:")
                  f >> outputfilenumber_;
              }
            } //while(f)

            //defining outputfilename by means of new testnumber
            outputfilename.str("");
            outputfilename << outputrootpath_ << "/StatMechOutput/NodalAverageDisp" << outputfilenumber_ << ".dat";
          }

          //set up new file with name "outputfilename" without writing anything into this file
          fp = fopen(outputfilename.str().c_str(), "w");
          fclose(fp);

          //increasing the number in the numbering file by one
          fpnumbering = fopen(numberingfilename.str().c_str(), "w");
          std::stringstream filecontent;
          filecontent << "Next Number: " << (outputfilenumber_ + 1);
          //write fileconent into file!
          fputs(filecontent.str().c_str(), fpnumbering);
          //close file
          fclose(fpnumbering);

          //initializing variables for positions of first and last node at the beginning
          chosennodeold_.PutScalar(0);
          sumsquarechosennode_=0.0;

          transdispold_ = Teuchos::rcp( new Epetra_Vector(*(discret_->DofRowMap())) );
          transdispold_->PutScalar(0.0);
          sumsquareallnodes_=0.0;
        }
      #endif
    }
    break;
    case INPAR::STATMECH::statout_viscoelasticity:
    {
      if(!discret_->Comm().MyPID())
      {
        //pointer to file into which each processor writes the output related with the dof of which it is the row map owner
        FILE* fp = NULL;

        //content to be written into the output file
        std::stringstream filecontent;

        //defining name of output file related to processor Id
        std::ostringstream outputfilename;
        outputfilename.str("");
        outputfilename << outputrootpath_ << "/StatMechOutput/ViscoElOutputProc.dat";

        fp = fopen(outputfilename.str().c_str(), "w");

        //filecontent << "Output for measurement of viscoelastic properties written by processor "<< discret_->Comm().MyPID() << std::endl;

        // move temporary std::stringstream to file and close it
        fputs(filecontent.str().c_str(), fp);
        fclose(fp);
      }
      linkernodepairs_ = Teuchos::null;
    }
      break;
    case INPAR::STATMECH::statout_networkrelax:
    {
      if(!discret_->Comm().MyPID())
      {
        //pointer to file into which each processor writes the output related with the dof of which it is the row map owner
        FILE* fp = NULL;

        //content to be written into the output file
        std::stringstream filecontent;

        //defining name of output file related to processor Id
        std::ostringstream outputfilename;
        outputfilename.str("");
        outputfilename << outputrootpath_ << "/StatMechOutput/CreepForces.dat";

        fp = fopen(outputfilename.str().c_str(), "w");

        //filecontent << "Output for measurement of viscoelastic properties written by processor "<< discret_->Comm().MyPID() << std::endl;

        // move temporary std::stringstream to file and close it
        fputs(filecontent.str().c_str(), fp);
        fclose(fp);

        std::ostringstream outpufilename2;
        outpufilename2 << outputrootpath_ << "/StatMechOutput/StructCOGInertia.dat";
        fp = fopen(outpufilename2.str().c_str(), "w");
        fputs(filecontent.str().c_str(), fp);
        fclose(fp);
      }
    }
    break;
    case INPAR::STATMECH::statout_none:
    default:
      break;
  }

  if(DRT::INPUT::IntegralValue<int>(statmechparams_,"GMSHOUTPUT"))
  {
    std::ostringstream filename;
    filename << outputrootpath_<<"/GmshOutput/network000000.pos";
    GmshOutput(dis,filename,istep,0.0);
  }

  return;
} // STATMECH::StatMechManager::InitOutput()

/*----------------------------------------------------------------------*
 | write special output for statistical mechanics (public)    cyron 09/08|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::Output(const int                            ndim,
                                       const double&                        time,
                                       const int&                           istep,
                                       const double&                        dt,
                                       const Epetra_Vector&                 dis,
                                       const Epetra_Vector&                 fint,
                                       Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager,
                                       bool                                 printscreen)
{
  /*in general simulations in statistical mechanics run over so many time steps that the amount of data stored in the error file
   * may exceed the capacity even of a server hard disk; thus, we rewind the error file in each time step so that the amount of data
   * does not increase after the first time step any longer*/
  Teuchos::ParameterList params = DRT::Problem::Instance()->StructuralDynamicParams();
  int numstep = params.get<int>("NUMSTEP", -1);
  bool printerr = params.get<bool> ("print to err", false);
  FILE* errfile = params.get<FILE*> ("err file", NULL);
  if (printerr)
    rewind(errfile);

  //the following variable makes sense in case of serial computing only; its use is not allowed for parallel computing!
  int num_dof = dis.GlobalLength();

  double starttime = statmechparams_.get<double>("STARTTIMEOUT",0.0);

  switch (DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT"))
  {
    case INPAR::STATMECH::statout_endtoendlog:
    {

      double endtoend = 0; //end to end length at a certain time step in single filament dynamics
      double DeltaR2 = 0;

      //as soon as system is equilibrated (after time STARTTIMEOUT) a new file for storing output is generated
      if ( (time > starttime && fabs(time-starttime)>dt/1e4) && (starttimeoutput_ == -1.0))
      {
        endtoendref_ = std::pow(pow((dis)[num_dof - 3] + 10 - (dis)[0], 2) + pow(
            (dis)[num_dof - 2] - (dis)[1], 2), 0.5);
        starttimeoutput_ = time;
        istart_ = istep;
      }
      if (time > starttimeoutput_ && starttimeoutput_ > -1.0)
      {
        endtoend = std::pow(pow((dis)[num_dof - 3] + 10 - (dis)[0], 2) + pow(
            (dis)[num_dof - 2] - (dis)[1], 2), 0.5);

        //applying in the following a well conditioned substraction formula according to Crisfield, Vol. 1, equ. (7.53)
        DeltaR2 = std::pow((endtoend * endtoend - endtoendref_ * endtoendref_) / (endtoend + endtoendref_), 2);

        //writing output: writing Delta(R^2) according to PhD thesis Hallatschek, eq. (4.60), where t=0 corresponds to starttimeoutput_
        if ((istep - istart_) % int(ceil(pow(10, floor(log10((time - starttimeoutput_) / (10* dt )))))) == 0)
        {

          //proc 0 write complete output into file, all other proc inactive
          if(!discret_->Comm().MyPID())
          {
            FILE* fp = NULL; //file pointer for statistical output file

            //name of output file
            std::ostringstream outputfilename;
            outputfilename << outputrootpath_ << "/StatMechOutput/EndToEnd" << outputfilenumber_ << ".dat";

            // open file and append new data line
            fp = fopen(outputfilename.str().c_str(), "a");

            //defining temporary std::stringstream variable
            std::stringstream filecontent;
            filecontent << std::scientific << std::setprecision(15) << time - starttimeoutput_ << "  " << DeltaR2 << std::endl;

            // move temporary std::stringstream to file and close it
            fputs(filecontent.str().c_str(), fp);
            fclose(fp);
          }
        }
      }

    }
    break;
    case INPAR::STATMECH::statout_endtoendconst:
    {
      /*in the following we assume that there is only a pulling force point Neumann condition of equal absolute
       *value on either filament length; we get the absolute value of the first one of these two conditions */
      double neumannforce;
      std::vector<DRT::Condition*> pointneumannconditions(0);
      discret_->GetCondition("PointNeumann", pointneumannconditions);
      if (pointneumannconditions.size() > 0)
      {
        const std::vector<double>* val = pointneumannconditions[0]->Get<std::vector<double> > ("val");
        neumannforce = fabs((*val)[0]);
      }
      else
        neumannforce = 0;

      double endtoend = 0.0; //end to end length at a certain time step in single filament dynamics

      //as soon as system is equilibrated (after time STARTTIMEOUT) a new file for storing output is generated
      if ((time > starttime && fabs(time-starttime)>dt/1e4) && (starttimeoutput_ == -1.0))
      {
        starttimeoutput_ = time;
        istart_ = istep;
      }

      if (time > starttimeoutput_ && starttimeoutput_ > -1.0)
      {
        //end to end vector
        LINALG::Matrix<3, 1> endtoendvector(true);
        for (int i = 0; i < ndim; i++)
        {
          endtoendvector(i) -= (discret_->gNode(0))->X()[i] + dis[i];
          endtoendvector(i) += (discret_->gNode(discret_->NumMyRowNodes() - 1))->X()[i] + dis[num_dof - discret_->NumDof(discret_->gNode(discret_->NumMyRowNodes() - 1)) + i];
        }

        endtoend = endtoendvector.Norm2();

        //writing output: current time and end to end distance are stored at each 100th time step
        if ((istep - istart_) % statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0)
        {

          //proc 0 write complete output into file, all other proc inactive
          if(discret_->Comm().MyPID()==0)
          {
            FILE* fp = NULL; //file pointer for statistical output file

            //name of output file
            std::ostringstream outputfilename;
            outputfilename << outputrootpath_ << "/StatMechOutput/E2E_" << discret_->NumMyRowElements() << '_' << dt<< '_' << neumannforce << '_' << outputfilenumber_ << ".dat";

            // open file and append new data line
            fp = fopen(outputfilename.str().c_str(), "a");
            //defining temporary std::stringstream variable
            std::stringstream filecontent;
            filecontent << std::scientific << std::setprecision(15) << time << "  "
                        << endtoend << " " << fint[num_dof - discret_->NumDof(
                           discret_->gNode(discret_->NumMyRowNodes() - 1))] << std::endl;
            // move temporary std::stringstream to file and close it
            fputs(filecontent.str().c_str(), fp);
            fclose(fp);
          }
        }
      }
    }
    break;
    /*computing and writing into file data about correlation of orientation of different elements
     *as it is considered in the context of the persistence length*/
    case INPAR::STATMECH::statout_orientationcorrelation:
    {

      //we need displacements also of ghost nodes and hence export displacment vector to column map format
      Epetra_Vector discol(*(discret_->DofColMap()), true);
      LINALG::Export(dis, discol);

      std::vector<double> arclength(discret_->NumMyColNodes() - 1, 0);
      std::vector<double> cosdiffer(discret_->NumMyColNodes() - 1, 0);

      //after initilization time write output cosdiffer in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps,
      //when discret_->NumMyRowNodes()-1 = 0,cosdiffer is always equil to 1!!
      if ((time > starttime && fabs(time-starttime)>dt/1e4) && (istep% statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0))
      {
        Epetra_SerialDenseMatrix coord;
        coord.Shape(discret_->NumMyColNodes(), ndim);

        for (int id=0; id<discret_->NumMyColNodes(); id++)
          for (int j=0; j<ndim; j++)
            coord(id, j) = (discret_->lColNode(id))->X()[j] + (discol)[(id) * (ndim - 1) * 3 + j];

        for (int id=0; id < discret_->NumMyColNodes() - 1; id++)
        {

          //calculate the deformed length of every element
          for (int j=0; j<ndim; j++)
          {
            arclength[id] += std::pow((coord(id + 1, j) - coord(id, j)), 2);
            cosdiffer[id] += (coord(id + 1, j) - coord(id, j)) * (coord(0 + 1, j) - coord(0, j));
          }

          //calculate the cosine difference referring to the first element
          //Dot product of the (id+1)th element with the 1st element and devided by the length of the (id+1)th element and the 1st element
          arclength[id] = std::pow(arclength[id], 0.5);
          cosdiffer[id] = cosdiffer[id] / (arclength[id] * arclength[0]);
        }

        //proc 0 write complete output into file, all other proc inactive
        if(!discret_->Comm().MyPID())
        {
          FILE* fp = NULL; //file pointer for statistical output file

          //name of output file
          std::ostringstream outputfilename;
          outputfilename.str("");
          outputfilename << outputrootpath_ << "/StatMechOutput/OrientationCorrelation" << outputfilenumber_ << ".dat";

          fp = fopen(outputfilename.str().c_str(), "a");
          std::stringstream filecontent;
          filecontent << istep;
          filecontent << std::scientific << std::setprecision(10);

          for (int id = 0; id < discret_->NumMyColNodes() - 1; id++)
            filecontent << " " << cosdiffer[id];

          filecontent << std::endl;
          fputs(filecontent.str().c_str(), fp);
          fclose(fp);
        }

      }

    }
    break;
    //the following output allows for anisotropic diffusion simulation of a quasi stiff polymer
    case INPAR::STATMECH::statout_anisotropic:
    {
      #ifndef CHOSENOUTPUTNODE
        //output in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps
        if ((time > starttime && fabs(time-starttime)>dt/1e4) && (istep % statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0))
        {
          //positions of first and last node in current time step (note: always 3D variables used; in case of 2D third component just constantly zero)
          LINALG::Matrix<3, 1> beginnew;
          LINALG::Matrix<3, 1> endnew;

          beginnew.PutScalar(0);
          endnew.PutScalar(0);

          if(printscreen)
            std::cout << "ndim: " << ndim << "\n";

          for (int i = 0; i < ndim; i++)
          {
            beginnew(i) = (discret_->gNode(0))->X()[i] + dis[i];
            endnew(i) = (discret_->gNode(discret_->NumMyRowNodes() - 1))->X()[i] + dis[num_dof - discret_->NumDof(discret_->gNode(discret_->NumMyRowNodes() - 1)) + i];
          }

          //unit direction vector for filament axis in last time step
          LINALG::Matrix<3, 1> axisold;
          axisold = endold_;
          axisold -= beginold_;
          axisold.Scale(1 / axisold.Norm2());

          //unit direction vector for filament axis in current time step
          LINALG::Matrix<3, 1> axisnew;
          axisnew = endnew;
          axisnew -= beginnew;
          axisnew.Scale(1 / axisnew.Norm2());

          //displacement of first and last node between last time step and current time step
          LINALG::Matrix<3, 1> dispbegin;
          LINALG::Matrix<3, 1> dispend;
          dispbegin = beginnew;
          dispbegin -= beginold_;
          dispend = endnew;
          dispend -= endold_;

          //displacement of middle point
          LINALG::Matrix<3, 1> dispmiddle;
          dispmiddle = dispbegin;
          dispmiddle += dispend;
          dispmiddle.Scale(0.5);
          sumdispmiddle_ += dispmiddle;

          //update sum of square displacement increments of middle point
          double incdispmiddle = dispmiddle.Norm2() * dispmiddle.Norm2();
          sumsquareincmid_ += incdispmiddle;

          //update sum of square displacement increments of middle point parallel to new filament axis (computed by scalar product)
          double disppar_square = std::pow(axisnew(0) * dispmiddle(0) + axisnew(1) * dispmiddle(1) + axisnew(2) * dispmiddle(2), 2);
          sumsquareincpar_ += disppar_square;

          //update sum of square displacement increments of middle point orthogonal to new filament axis (via crossproduct)
          LINALG::Matrix<3, 1> aux;
          aux(0) = dispmiddle(1) * axisnew(2) - dispmiddle(2) * axisnew(1);
          aux(1) = dispmiddle(2) * axisnew(0) - dispmiddle(0) * axisnew(2);
          aux(2) = dispmiddle(0) * axisnew(1) - dispmiddle(1) * axisnew(0);
          double disport_square = aux.Norm2() * aux.Norm2();
          sumsquareincort_ += disport_square;

          //total displacement of rotational angle (in 2D only)
          double incangle = 0;
          if (ndim == 2)
          {
            //angle of old axis relative to x-axis
            double phiold = acos(axisold(0) / axisold.Norm2());
            if (axisold(1) < 0)
              phiold *= -1;

            //angle of new axis relative to x-axis
            double phinew = acos(axisnew(0) / axisnew.Norm2());
            if (axisnew(1) < 0)
              phinew *= -1;

            //angle increment
            incangle = phinew - phiold;
            if (incangle > M_PI)
            {
              incangle -= 2*M_PI;
              incangle *= -1;
            }
            if (incangle < -M_PI)
            {
              incangle += 2*M_PI;
              incangle *= -1;
            }

            //update absolute rotational displacement compared to reference configuration
            sumsquareincrot_ += incangle * incangle;
            sumrotmiddle_ += incangle;
          }

          //proc 0 write complete output into file, all other proc inactive
          if(!discret_->Comm().MyPID())
          {
            FILE* fp = NULL; //file pointer for statistical output file

            //name of output file
            std::ostringstream outputfilename;
            outputfilename << outputrootpath_ << "/StatMechOutput/AnisotropicDiffusion" << outputfilenumber_ << ".dat";

            // open file and append new data line
            fp = fopen(outputfilename.str().c_str(), "a");

            //defining temporary std::stringstream variable
            std::stringstream filecontent;
            filecontent << std::scientific << std::setprecision(15) << dt << " "<< sumsquareincmid_ << " "<< sumdispmiddle_.Norm2() * sumdispmiddle_.Norm2() <<std::endl;
  //              << sumsquareincmid_ << " " << sumsquareincpar_ << " "
  //              << sumsquareincort_ << " " << sumsquareincrot_ << " "
  //              << sumdispmiddle_.Norm2() * sumdispmiddle_.Norm2() << " "
  //              << sumrotmiddle_ * sumrotmiddle_ << std::endl;

            // move temporary std::stringstream to file and close it
            fputs(filecontent.str().c_str(), fp);
            fclose(fp);

            //write position
            LINALG::Matrix<3,1> midpoint(beginnew);
            midpoint += endnew;
            midpoint.Scale(0.5);
            std::ostringstream outputfilename2;
            outputfilename2 << outputrootpath_ << "/StatMechOutput/AnisotropicMidPosition" << outputfilenumber_ << ".dat";

            // open file and append new data line
            fp = fopen(outputfilename2.str().c_str(), "a");

            //defining temporary std::stringstream variable
            std::stringstream filecontent2;
            filecontent2 << std::scientific << std::setprecision(15) <<midpoint(0)<<" "<<midpoint(1)<<" "<<midpoint(2)<< std::endl;

            // move temporary std::stringstream to file and close it
            fputs(filecontent2.str().c_str(), fp);
            fclose(fp);
          }

          //new positions in this time step become old positions in last time step
          beginold_ = beginnew;
          endold_ = endnew;
        }
      #else
        if(discret_->Comm().MyPID()!=0)
          dserror("The special output type nodalaveragedisp is only implemented in serial so far!");

        const double outputtimeinterval = 0.0005;
        //output in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps
        if (time >= numoutputstep_*outputtimeinterval -1.0e-8)
        {
          numoutputstep_++;
          std::cout << "Output of nodalaveragedisp written!" << std::endl;
          LINALG::Matrix<3, 1> chosennodenew;
          const int nodeid = CHOSENOUTPUTNODE;

          Epetra_Vector trans_dis_new(dis);

          int numnode = discret_->NodeRowMap()->NumMyElements();

          //Set rotational/tangential displacements to zero
          for(int i=0;i<numnode;i++)
          {
            // get node pointer
            DRT::Node* node = discret_->lRowNode(i);

            // get GIDs of this node's degrees of freedom
            std::vector<int> dofnode = discret_->Dof(node);

            if(dofnode.size()!=6)
              dserror("This methode only works for beam elements with 6 DoFs per node so far!");

            //Set displacements in rotational/tangential DoFs to zero
            trans_dis_new[discret_->DofRowMap()->LID(dofnode[3])]=0.0;
            trans_dis_new[discret_->DofRowMap()->LID(dofnode[4])]=0.0;
            trans_dis_new[discret_->DofRowMap()->LID(dofnode[5])]=0.0;
          }

          Epetra_Vector delta_dis(trans_dis_new);
          delta_dis.Update(-1.0,*transdispold_,1.0);

          double norm_delta_dis=0.0;
          delta_dis.Norm2(&norm_delta_dis);

          sumsquareallnodes_+=norm_delta_dis*norm_delta_dis;

          chosennodenew.PutScalar(0);

          for (int i = 0; i < ndim; i++)
          {
            chosennodenew(i) = dis[discret_->DofColMap()->LID((discret_->Dof(discret_->gNode(nodeid)))[i])];
          }

          LINALG::Matrix<3, 1> chosennodedisp(chosennodenew);
          chosennodedisp-=chosennodeold_;
          sumsquarechosennode_+=chosennodedisp.Norm2()*chosennodedisp.Norm2();

          //proc 0 write complete output into file, all other proc inactive
          if(!discret_->Comm().MyPID())
          {
            FILE* fp = NULL; //file pointer for statistical output file

            std::ostringstream outputfilename;
            outputfilename << outputrootpath_ << "/StatMechOutput/NodalAverageDisp" << outputfilenumber_ << ".dat";


            // open file and append new data line
            fp = fopen(outputfilename.str().c_str(), "a");

            //defining temporary std::stringstream variable
            std::stringstream filecontent;
            filecontent << std::scientific << std::setprecision(15) << sumsquarechosennode_ << " " << sumsquareallnodes_ << std::endl;

            // move temporary std::stringstream to file and close it
            fputs(filecontent.str().c_str(), fp);
            fclose(fp);
          }

          //new positions in this time step become old positions in last time step
          chosennodeold_ = chosennodenew;
          transdispold_->Update(1.0,trans_dis_new,0.0);
        }
      #endif
    }
    break;
    case INPAR::STATMECH::statout_viscoelasticity:
    {
      //output in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps (or for the very last step)
      if ((time>=starttime && istep<numstep && (istep-istart_) % statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0) || fabs(time-starttime)<1e-10)
      {
        // name of file into which output is written
        std::ostringstream filename;
        filename << outputrootpath_ << "/StatMechOutput/ViscoElOutputProc.dat";
        ViscoelasticityOutput(time, dis, fint, filename);

        // additional Density-Density-Correlation Output
        if((istep-istart_) % (10*statmechparams_.get<int> ("OUTPUTINTERVALS", 1)) == 0)
        {
        //  std::ostringstream ddcorrfilename;
         // ddcorrfilename << outputrootpath_ << "/StatMechOutput/DensityDensityCorrFunction_"<<std::setw(6) << std::setfill('0') << istep <<".dat";
          //DDCorrOutput(dis, ddcorrfilename, istep, dt);
        }

        std::ostringstream linkerfilename;
        linkerfilename << outputrootpath_ << "/StatMechOutput/LinkerPosition_"<<std::setw(6) << std::setfill('0') << istep <<".dat";
        OutputLinkerPositions(linkerfilename);

        std::ostringstream nodepairfilename;
        nodepairfilename << outputrootpath_ << "/StatMechOutput/NodePairPosition_"<<std::setw(6) << std::setfill('0') << istep <<".dat";
        OutputSlidingMotion(dis, nodepairfilename);

        std::ostringstream freefillengthname;
        freefillengthname << outputrootpath_ << "/StatMechOutput/FreeFilLength_"<<std::setw(6) << std::setfill('0') << istep <<".dat";
        OutputFreeFilamentLength(dis,freefillengthname);

        std::ostringstream matforcefilename;
        matforcefilename << outputrootpath_ << "/StatMechOutput/IntMatForces_"<<std::setw(6) << std::setfill('0') << istep <<".dat";
        OutputElementMaterialInternalForces(dis,matforcefilename);

        if(statmechparams_.get<double>("DELTABELLSEQ", 0.0)!=0.0 && (linkermodel_==statmech_linker_bellseq || linkermodel_==statmech_linker_bellseqintpol))
        {
          std::ostringstream forcedepfilename;
          forcedepfilename << outputrootpath_ << "/StatMechOutput/UnbindingProbability_"<<std::setw(6) << std::setfill('0') << istep <<".dat";
          BellsEquationOutput(dis, forcedepfilename, dt);
        }
      }
    }
    break;
    case INPAR::STATMECH::statout_networkrelax:
    {
      //output in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps (or for the very last step)
      if ((time>=starttime && istep<numstep && (istep-istart_) % statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0) || fabs(time-starttime)<1e-8)
      {
        // name of file into which output is written
        std::ostringstream filename;
        filename << outputrootpath_ << "/StatMechOutput/CreepForces.dat";
        ViscoelasticityOutput(time, dis, fint, filename);

        std::ostringstream filename2;
        filename2 << outputrootpath_ << "/StatMechOutput/StructCOGInertia.dat";
        StructureCOGInertiaTensorOutput(istep, time, dis, filename2);
      }
    }
    break;
    case INPAR::STATMECH::statout_networkcreep:
    {
      if ((time>=starttime && istep<numstep && (istep-istart_) % statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0) || fabs(time-starttime)<1e-10)
      {
        if(!dbcnodesets_.size())
          dserror("No Dirichlet node sets have been found! Check DBCTYPE and SPECIAL_OUTPUT in your input file.");

        Epetra_Vector discol(*(discret_->DofColMap()),true);
        LINALG::Export(dis,discol);

        if(!discret_->Comm().MyPID())
        {
          std::ostringstream outputfilename;
          outputfilename << outputrootpath_ << "/StatMechOutput/CreepDisplacement_"<<std::setw(6) << std::setfill('0') << istep <<".dat";
          CreepDisplacementOutput(discol, outputfilename);
        }
      }
    }
    break;
    case INPAR::STATMECH::statout_networkdispfield:
    {
      if ((time>=starttime && istep<numstep && (istep-istart_) % statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0) || fabs(time-starttime)<1e-10)
      {
          std::ostringstream dispfilename;
          std::ostringstream nodeposfilename;
          dispfilename << outputrootpath_ << "/StatMechOutput/Displacements_"<<std::setw(6) << std::setfill('0') << istep <<".dat";
          nodeposfilename << outputrootpath_ << "/StatMechOutput/NodePositions_"<<std::setw(6) << std::setfill('0') << istep <<".dat";
          OutputNodalDisplacements(dis, dispfilename);
          OutputNodalPositions(dis,nodeposfilename);

          INPAR::STATMECH::NBCType nbctype = DRT::INPUT::IntegralValue<INPAR::STATMECH::NBCType>(statmechparams_,"NBCTYPE");
          if(nbctype == INPAR::STATMECH::nbctype_randompointforce)
          {
            std::ostringstream forcefilename;
            forcefilename << outputrootpath_ <<"/StatMechOutput/PointForces.dat";
            OutputNeumannPointForce(time,forcefilename);
          }
      }
    }
    break;
    case INPAR::STATMECH::statout_structanaly:
    {
      //output in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps
      if( istep % statmechparams_.get<int>("OUTPUTINTERVALS",1) == 0 )
      {
        if(periodlength_->at(0) == periodlength_->at(1) && periodlength_->at(0) == periodlength_->at(2))
        {
          std::ostringstream filename;
          filename << outputrootpath_ << "/StatMechOutput/DensityDensityCorrFunction_"<<std::setw(6) << std::setfill('0') << istep <<".dat";
          DDCorrOutput(dis, filename, istep, dt);
        }
        else
          dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));
      }

      std::ostringstream filename02;
      filename02 << outputrootpath_ <<"/StatMechOutput/LinkerUnbindingTimes.dat";
      LinkerUnbindingTimes(filename02);
    }
    break;
    case INPAR::STATMECH::statout_octree:
    {
      if(beamcmanager!=Teuchos::null && beamcmanager->OcTree()==Teuchos::null)
        dserror("Beam Contact withtout Octree search detected! Either set BEAMS_OCTREE to something other than -None- or disable this output!");

      //output in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps
      if( ((time>=starttime && (istep-istart_) % statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0) || fabs(time-starttime)<1e-8) && beamcmanager!=Teuchos::null)
      {
        std::map<int, LINALG::Matrix<3, 1> > currentpositions;
        std::map<int, LINALG::Matrix<3, 1> > currentrotations;
        Epetra_Vector discol(*discret_->DofColMap(), true);
        LINALG::Export(dis, discol);
        GetNodePositionsFromDisVec(discol, currentpositions, currentrotations, true);
        beamcmanager->OcTree()->OctTreeSearch(currentpositions, istep);
      }
    }
    break;
    case INPAR::STATMECH::statout_loom:
    {
      if ((time>=starttime && (istep-istart_) % statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0) || fabs(time-starttime)<1e-8)
      {
        std::ostringstream filename1;
        std::ostringstream filename2;
        filename1 << outputrootpath_ << "/StatMechOutput/DoubleBondDistances.dat";
        filename2 << outputrootpath_ << "/StatMechOutput/LoomLinkerpositions.dat";
        LoomOutput(dis,filename1,filename2);

        std::ostringstream filename3;
        filename3 << outputrootpath_ << "/StatMechOutput/ForceMeasurement.dat";
        LoomOutputAttraction(dis, filename3, istep);

        std::ostringstream filename4;
        filename4 << outputrootpath_ << "/StatMechOutput/CrosslinkerCoverage.dat";
        CrosslinkCoverageOutput(dis, filename4);
        // node coords output for modal analysis
      }
    }
    break;
    case INPAR::STATMECH::statout_loomelnrg:
    {
      if ((time>=starttime && (istep-istart_) % statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0) || fabs(time-starttime)<1e-8)
      {
        std::ostringstream filename;
        filename << outputrootpath_ << "/StatMechOutput/LoomElasticEnergy.dat";
        LoomOutputElasticEnergy(dis,dt,filename);

        std::ostringstream filename2;
        filename2 << outputrootpath_ << "/StatMechOutput/NodeDisp_"<<std::setw(6) << std::setfill('0') << istep <<".dat";
        OutputNodalDisplacements(dis,filename2);

        std::ostringstream filename3;
        filename3 << outputrootpath_ << "/StatMechOutput/InternalForces_"<<std::setw(6) << std::setfill('0') << istep <<".dat";
        OutputElementSpatialInternalForces(filename3);
      }
    }
    break;
    case INPAR::STATMECH::statout_motassay:
    {
      if ((time>=starttime && (istep-istart_) % statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0) || fabs(time-starttime)<1e-8)
      {
        std::ostringstream filename;
        filename  << outputrootpath_ << "/StatMechOutput/filcoords.dat";
        std::ostringstream filename2;
        filename2  << outputrootpath_ << "/StatMechOutput/linkerforces.dat";
        MotilityAssayOutput(dis,istep,time,dt,filename,filename2);
      }
    }
    break;
    case INPAR::STATMECH::statout_linkerlength:
    {
      if ((time>=starttime && (istep-istart_) % statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0) || fabs(time-starttime)<1e-8)
      {
        std::ostringstream linkerlengthfilename;
        linkerlengthfilename << outputrootpath_ << "/StatMechOutput/LinkerLength_"<<std::setw(6) << std::setfill('0') << istep <<".dat";
        OutputLinkerLength(linkerlengthfilename);
      }
    }
    break;
    case INPAR::STATMECH::statout_deltatheta:
    {
      if ((time>=starttime && (istep-istart_) % statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0) || fabs(time-starttime)<1e-8)
      {
        std::ostringstream deltathetafilename;
        deltathetafilename << outputrootpath_ << "/StatMechOutput/DeltaTheta_"<<std::setw(6) << std::setfill('0') << istep <<".dat";
        OutputDeltaTheta(deltathetafilename);
      }
    }
    break;
    case INPAR::STATMECH::statout_none:
    default:
    break;
  }

  // handling gmsh output seperately
  if(DRT::INPUT::IntegralValue<int>(statmechparams_,"GMSHOUTPUT") && (time>=starttime && (istep-istart_) % statmechparams_.get<int> ("GMSHOUTINTERVALS", 100) == 0) )
  {
    /*construct unique filename for gmsh output with two indices: the first one marking the time step number
     * and the second one marking the newton iteration number, where numbers are written with zeros in the front
     * e.g. number one is written as 000001, number fourteen as 000014 and so on;*/

    // first index = time step index
    std::ostringstream filename;

    //creating complete file name dependent on step number with 6 digits and leading zeros
    if (istep<1000000)
      filename << outputrootpath_ << "/GmshOutput/network"<< std::setw(6) << std::setfill('0') << istep <<".pos";
    else
      dserror("Gmsh output implemented for a maximum of 999999 steps");

    //calling method for writing Gmsh output
    if(beamcmanager!=Teuchos::null)
      GmshOutput(dis,filename,istep,time,beamcmanager);
    else
      GmshOutput(dis,filename,istep,time);
  }

  return;
} // STATMECH::StatMechManager::Output()

/*------------------------------------------------------------------------------*
 | Output of distances between doubly bound linkers of the horizontal filament  |
 | in a loom type network                                                       |
 |                                                       (private) mueller 2/12 |
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::LoomOutput(const Epetra_Vector& disrow,
                                           const std::ostringstream& nearestneighborfilename,
                                           const std::ostringstream& linkerpositionsfilename)
{
  if(networktype_ != statmech_network_casimir)
      dserror("For Casimir force related output, set NETWORKTYPE to 'casimir' in your input file!");

  //detect loom type
  enum LoomType
  {
    loom_none,
    loom_std, // standard loom
    loom_singlefil,  // one single filament (with pinning linkers)
    loom_twovfil_free, // 2 vertical freely filaments
    loom_twovfil_spring, // 2 vertival filaments connected by a linear spring (truss)
  };

  LoomType loomtype = loom_none;
  int filnum = (int)(*filamentnumber_)[filamentnumber_->MyLength()-1];
  switch(filnum)
  {
    case 0:
      loomtype = loom_singlefil;
    break;
    case 2:
    {
      int tmp = 0;
      for(int i=0; i<discret_->NumMyRowElements(); i++)
      {
        DRT::Element* element = discret_->lRowElement(i);
        const DRT::ElementType & eot = element->ElementType();
        // assumption: the Truss3 element has one explicit owner and therefore occurs only once in the row map
        if(eot == DRT::ELEMENTS::Truss3Type::Instance())
        {
          tmp = 1;
          break;
        }
      }
      // Check for greates value and return to all processors in order to correctly determine the loom type
      int max = -1;
      discret_->Comm().MaxAll(&tmp,&max,1);
      if(max==1)
        loomtype = loom_twovfil_spring;
      else
        loomtype = loom_twovfil_free;
    }
    break;
    default: loomtype = loom_std;
  }

  Epetra_Vector discol(*discret_->DofColMap(),true);
  LINALG::Export(disrow,discol);
  std::map<int, LINALG::Matrix<3, 1> > currentpositions;
  std::map<int, LINALG::Matrix<3, 1> > currentrotations;
  GetNodePositionsFromDisVec(discol, currentpositions, currentrotations, true);

  // We assume in the following that the horizontal filament more or less stays horizontal without tilting too much.
  // Also, we assume that the vertical filaments are discretized in the same manner, i.e. same element lengths
  if(discret_->Comm().MyPID()==0)
  {
    FILE* fp = NULL;
    fp = fopen(nearestneighborfilename.str().c_str(), "a");
    FILE* fp2 = NULL;
    fp2 = fopen(linkerpositionsfilename.str().c_str(), "a");
    std::stringstream distances;
    std::stringstream linkerpositions;

    // step 1: determine the first double bond between vertical and horizontal filaments
    int firstvfil = -1;
    // node IDs at first double bond (smallest coordinate value) between horizontal and a vertical filament
    std::vector<int> nodeIDs;
    std::vector<int> evalnodeIDs;
    nodeIDs.clear();
    switch(loomtype)
    {
      case loom_singlefil:
      {
        for(int i=0; i<filamentnumber_->MyLength(); i++)
        {
          if((int)(*bspotstatus_)[i]>-1)
            if((*numbond_)[(int)(*bspotstatus_)[i]]>0.9)
              evalnodeIDs.push_back(bspotcolmap_->GID(i));
        }
      }
      break;
      case loom_twovfil_free:
      {
        // choose the GIDs of the two center nodes of the vertical filaments for  distance measurements
        int nodeID1 = (int)(floor((double)(filamentnumber_->MyLength())/4.0));
        // assumption that both filaments are of the same discretization length
        int nodeID2 = nodeID1 + filamentnumber_->MyLength()/2;
        evalnodeIDs.push_back(nodeID1);
        evalnodeIDs.push_back(nodeID2);

      }
      break;
      default:
      {
        for(int i=0; i<filamentnumber_->MyLength(); i++)
        {
          if((int)(*filamentnumber_)[i]==0) // horizontal filament
          {
            if(nodeIDs.empty() && (*numbond_)[(int)(*bspotstatus_)[i]]>1.9)
            {
              nodeIDs.push_back(bspotcolmap_->GID(i));
              int currlink = (int)(*bspotstatus_)[i];
              for(int j=0; j<crosslinkerbond_->NumVectors(); j++)
                if((int)(*crosslinkerbond_)[j][currlink]!=nodeIDs[0])
                {
                  nodeIDs.push_back((int)(*crosslinkerbond_)[j][currlink]);
                  firstvfil = (int)(*filamentnumber_)[discret_->NodeColMap()->LID(nodeIDs[1])];
                  break;
                }
            }
          }
          else if((int)(*filamentnumber_)[i] == (int)evalnodeIDs.size()+1)
            evalnodeIDs.push_back(bspotcolmap_->GID(i));
        }
      }
    }

    // evaluate distances between binding spots / vertical filaments
    if(!nodeIDs.empty() || loomtype == loom_singlefil)
    {
      // starting at at least two occupied binding spots
      if((int)evalnodeIDs.size()>1)
      {
        // single filament simulations, i.e. no real crosslinker elements -> nodeIDs are IDs on horizontal filament
        if(loomtype == loom_singlefil)
        {
          nodeIDs = evalnodeIDs;
        }
        else // standard case with real crosslinker elements
        {
          int vnodeoffset = nodeIDs[1] - evalnodeIDs[firstvfil-1];
          for(int i=0; i<(int)evalnodeIDs.size(); i++)
            evalnodeIDs.at(i) += vnodeoffset;

          // sort nodes by ascending x-coordinate
          for(int i=0; i<(int)evalnodeIDs.size()-1; i++)
          {
            std::map< int,LINALG::Matrix<3,1> >::const_iterator posi = currentpositions.find(evalnodeIDs[i]);
            for(int j=i+1; j<(int)evalnodeIDs.size(); j++)
            {
              std::map< int,LINALG::Matrix<3,1> >::const_iterator posj = currentpositions.find(evalnodeIDs[j]);
              if((posj->second)(0)<(posi->second)(0))
              {
                int tmp = evalnodeIDs[i];
                evalnodeIDs.at(i) = evalnodeIDs[j];
                evalnodeIDs.at(j) = tmp;
              }
            }
          }
        }

        // gather positions in one vector
        std::vector<LINALG::Matrix<3,1> > evalnodepositions;
        evalnodepositions.clear();
        for(int i=0; i<(int)evalnodeIDs.size(); i++)
        {
          std::map< int,LINALG::Matrix<3,1> >::const_iterator pos = currentpositions.find(evalnodeIDs[i]);
          evalnodepositions.push_back(pos->second);
        }

        // pre-evaluation of maximal number of nonequal nearest neighbors (max 3)
        int maxnumevalneighbors = 4;

        if((int)evalnodepositions.size()>1)
        {
          int maxnearneighbors = (int)(ceil((double(evalnodepositions.size())/2.0)))-1;
          if(maxnearneighbors>maxnumevalneighbors)
            maxnearneighbors = maxnumevalneighbors;

          // calculate nearest neighbor to 4th nearest neighbor distances , no shift ac
          double globalarclength = 0.0;
          for(int i=0; i<(int)evalnodepositions.size(); i++)
          {
            std::map< int,LINALG::Matrix<3,1> >::const_iterator pos0 = currentpositions.find(evalnodeIDs[i]);
            // write linker positions to stream
            if(loomtype == loom_singlefil)
            {
              if(i==0) // first binding spot
              {
                for(int k=0; k<evalnodeIDs[i]; k++)
                {
                  std::map< int,LINALG::Matrix<3,1> >::const_iterator posk = currentpositions.find(k);
                  std::map< int,LINALG::Matrix<3,1> >::const_iterator poskp = currentpositions.find(k+1);
                  LINALG::Matrix<3,1> diffk = (poskp->second);
                  diffk -= (posk->second);
                  globalarclength += diffk.Norm2();
                }
              }
              else // subsequent binding spots
              {
                for(int k=evalnodeIDs[i-1]; k<evalnodeIDs[i]; k++)
                {
                  std::map< int,LINALG::Matrix<3,1> >::const_iterator posk = currentpositions.find(k);
                  std::map< int,LINALG::Matrix<3,1> >::const_iterator poskp = currentpositions.find(k+1);
                  LINALG::Matrix<3,1> diffk = (poskp->second);
                  diffk -= (posk->second);
                  globalarclength += diffk.Norm2();
                }
              }
            }
            else
            {
              globalarclength = 0.0;
              for(int k=0; k<evalnodeIDs[i]; k++)
              {
                std::map< int,LINALG::Matrix<3,1> >::const_iterator posk = currentpositions.find(k);
                std::map< int,LINALG::Matrix<3,1> >::const_iterator poskp = currentpositions.find(k+1);
                LINALG::Matrix<3,1> diffk = (poskp->second);
                diffk -= (posk->second);
                globalarclength += diffk.Norm2();
              }
            }
            linkerpositions<<std::scientific<<std::setprecision(6)<<(pos0->second)(0)<<"\t"<<(pos0->second)(1)<<"\t"<<(pos0->second)(2)<<"\t"<<globalarclength<<std::endl;

            // number of emtpy blocks in a line of the output file
            int emptyblocks = maxnumevalneighbors;
            for(int j=i+1; j<i+maxnearneighbors+1; j++)
            {
              // reached the end of the vector
              if(j>=(int)evalnodepositions.size())
                break;

              emptyblocks--;

              std::map< int,LINALG::Matrix<3,1> >::const_iterator pos1 = currentpositions.find(evalnodeIDs.at(j));

              // spatial distance vector
              LINALG::Matrix<3,1> diff = (pos1->second);
              diff -= (pos0->second);

              // arc length
              if(nodeIDs[j]<=nodeIDs[i])
                dserror("Node IDs are not sorted in an ascending manner!");

              double diffarclength = 0.0;
              if(loomtype == loom_singlefil)
              {
                for(int k=nodeIDs[i] ; k<nodeIDs[j]; k++)
                {
                  std::map< int,LINALG::Matrix<3,1> >::const_iterator posk = currentpositions.find(k);
                  std::map< int,LINALG::Matrix<3,1> >::const_iterator poskp = currentpositions.find(k+1);
                  LINALG::Matrix<3,1> diffk = (poskp->second);
                  diffk -= (posk->second);
                  diffarclength += diffk.Norm2();
                }
                distances<<std::scientific<<std::setprecision(6)<<diffarclength<<"\t";
              }
              else
                distances<<std::scientific<<std::setprecision(6)<<diff(0)<<"\t"<<diff.Norm2()<<"\t";
            }
            if(emptyblocks>0)
            {
              // determine dummy entries for this line
              if(loomtype==loom_singlefil)
                for(int i=0; i<emptyblocks; i++)
                  distances<<std::scientific<<std::setprecision(6)<<-1e9<<"\t";
              else
                for(int i=0; i<emptyblocks; i++)
                  distances<<std::scientific<<std::setprecision(6)<<-1e9<<"\t"<<-1e9<<"\t";
            }
            distances<<std::endl;
          }
          // time step divider line
          linkerpositions<<std::scientific<<std::setprecision(6)<<-1e9<<"\t"<<-1e9<<"\t"<<-1e9<<"\t"<<-1e9<<std::endl;
        }
      }
    }
    fputs(distances.str().c_str(), fp);
    fclose(fp);
    fputs(linkerpositions.str().c_str(), fp2);
    fclose(fp2);
  }
  return;
}

/*------------------------------------------------------------------------------*                                                 |
 | Measure the force between two vertical filaments by using a truss element    |
 |                                                        (public) mueller 5/12 |
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::LoomOutputAttraction(const Epetra_Vector& disrow, const std::ostringstream& filename, const int& step)
{
  if(networktype_ != statmech_network_casimir)
    dserror("For Casimir force related output, set NETWORKTYPE to 'casimir' in your input file!");

  for(int i=0; i<discret_->NumMyRowElements(); i++)
  {
    DRT::Element* element = discret_->lRowElement(i);
    const DRT::ElementType & eot = element->ElementType();

    // assumption: the Truss3 element has one explicit owner and therefore occurs only once in the row map
    if(eot == DRT::ELEMENTS::Truss3Type::Instance())
    {
      // Get information on Boundary Conditions, i.e. the number of spatial dimensions (2D,3D)
      std::vector<DRT::Condition*> dirichlet;
      discret_->GetCondition("Dirichlet",dirichlet);
      if (!dirichlet.size())
        dserror("No Dirichlet boundary conditions in discretization");
      // choose the first node which is not the (first) end node of the horizontal filament
      const std::vector<int>* onoff = NULL;
      for(int j=0; j<(int)dirichlet.size(); j++)
        if(dirichlet.at(j)->Type() == DRT::Condition::PointDirichlet && dirichlet.at(j)->ContainsNode(0))
        {
          onoff = dirichlet.at(j)->Get<std::vector<int> >("onoff");
          break;
        }

      // retrieve the internal force vector of the element
      Teuchos::RCP<Epetra_SerialDenseVector> force = dynamic_cast<DRT::ELEMENTS::Truss3*>(element)->InternalForceVector();
      LINALG::Matrix<3,1> fint;
      for(int j=0; j<(int)fint.M(); j++)
        fint(j) = (double)onoff->at(j)*(*force)[j];
      // retrieve the current positions of the element's nodes
      //get pointer at a node
      const DRT::Node* node0 = element->Nodes()[0];
      const DRT::Node* node1 = element->Nodes()[1];

      // calculate the current length of the truss element (should be close to the x-difference in the 2D-case)
      std::vector<int> dofnode0 = discret_->Dof(node0);
      std::vector<int> dofnode1 = discret_->Dof(node1);

      // nodal coordinates
      LINALG::SerialDenseMatrix coord(3,2,true);
      // indicator for cuts due to periodic BCs
      LINALG::SerialDenseMatrix cut(3,1,true);

      for(int j=0; j<coord.M(); j++)
      {
        coord(j,0) = node0->X()[j] + disrow[discret_->DofRowMap()->LID(dofnode0[j])];
        coord(j,1) = node1->X()[j] + disrow[discret_->DofRowMap()->LID(dofnode1[j])];
      }

      // check if truss element is cut by  periodic boundary conditions (assumption: l<PeriodLength/2) and unshift
      for (int dof=0; dof<coord.M(); dof++)
      {
        if (fabs(coord(dof,1)-periodlength_->at(dof)-coord(dof,0)) < fabs(coord(dof,1) - coord(dof,0)))
          coord(dof,1) -= periodlength_->at(dof);
        if (fabs(coord(dof,1)+periodlength_->at(dof) - coord(dof,0)) < fabs(coord(dof,1)-coord(dof,0)))
          coord(dof,1) += periodlength_->at(dof);
      }

      // get the distance vector
      LINALG::Matrix<3,1> distance;
      for(int j=0; j<(int)distance.M(); j++)
        distance(j) = coord(j,1)-coord(j,0);

      // retrieve reference length of the truss
      double l0 = dynamic_cast<DRT::ELEMENTS::Truss3*>(element)->L0();
      if(distance.Norm2()<l0)
      {
        if(fint(0)>0.0)
          fint(0) *= -1.0;
      }
      else
      {
        if(fint(0)<0.0)
          fint(0) *= -1.0;
      }

      FILE* fp = NULL;
      fp = fopen(filename.str().c_str(), "a");
      std::stringstream internalforce;
      internalforce << std::scientific << std::setprecision(15) << distance.Norm2()<< "  " << (distance.Norm2()-l0)/l0<<"  "<< fint(0) <<std::endl;

      fputs(internalforce.str().c_str(), fp);
      fclose(fp);

//      // every 20000 steps (in the 2 vertical filament case only), change the stiffness of the truss
//      if(step % 20000 == 0 && (*filamentnumber_)[filamentnumber_->MyLength()-1]==2)
//      {
//        double csecnp = dynamic_cast<DRT::ELEMENTS::Truss3*>(element)->CSec()/1.1;
//        dynamic_cast<DRT::ELEMENTS::Truss3*>(element)->SetCrossSec(csecnp);
//      }
    }
  }
  return;
}
/*------------------------------------------------------------------------------*
 | output of elastic internal energy (special version for Loom setups           |
 |                                                        (public) mueller 4/11 |
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::LoomOutputElasticEnergy(const Epetra_Vector&                 disrow,
                                                        const double&                        dt,
                                                        const std::ostringstream&            filename,
                                                        Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager)
{
  Epetra_Vector discol(*discret_->DofColMap(), true);
  LINALG::Export(disrow, discol);

  // Compute internal energy
  std::vector<double> internalenergy;
  internalenergy.clear();
  const Teuchos::RCP<Epetra_Vector> disp = Teuchos::rcp(new Epetra_Vector(disrow));
  ComputeInternalEnergy(disp, internalenergy,dt, filename, false, false);

  // retrieve distance between fixed end (at x=0) and hoop position (x=x_hoop)
  if(!discret_->Comm().MyPID())
  {
    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "a");
    std::stringstream elasticenergy;

    if(beamcmanager!=Teuchos::null)
    {
      for(int i=0; i<filamentnumber_->MyLength(); i++)
      {
        if((*filamentnumber_)[i]>0)
        {
          DRT::Node* ringnode = discret_->lColNode(i);
          std::vector<int> dofnode = discret_->Dof(ringnode);
          double xring = discret_->lColNode(i)->X()[0] + discol[discret_->DofColMap()->LID(dofnode[0])];
          elasticenergy << std::scientific << std::setprecision(15) << xring << " ";
          break;
        }
      }
    }

    for(int j=0; j<(int)internalenergy.size(); j++)
      elasticenergy << std::scientific << std::setprecision(15)<< internalenergy[j] <<" ";
    elasticenergy<<std::endl;
    fputs(elasticenergy.str().c_str(), fp);
    fclose(fp);
  }
  discret_->Comm().Barrier();
  return;
}

/*------------------------------------------------------------------------------*
| motility assay output                                  (private) mueller 10/13|
*-------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::MotilityAssayOutput(const Epetra_Vector&      disrow,
                                                    const int&                istep,
                                                    const double&             timen,
                                                    const double&             dt,
                                                    const std::ostringstream& nodeposfilename,
                                                    const std::ostringstream& forcefilename)
{
  double kt = statmechparams_.get<double>("KT",0.00404531);
  double delta = statmechparams_.get<double>("DELTABELLSEQ", 0.0);
  double koff0 = 0.0;
  if (timen <= actiontime_->at(1) || (timen>actiontime_->at(1) && fabs(actiontime_->at(1))<dt/1e4))
    koff0 = statmechparams_.get<double> ("K_OFF_start", 0.0);
  else
    koff0 = statmechparams_.get<double> ("K_OFF_end", 0.0);
  int subfil = statmechparams_.get<int>("NUMSUBSTRATEFIL",0);

  // get filament number conditions
  std::vector<DRT::Condition*> filaments(0);
  discret_->GetCondition("FilamentNumber", filaments);

  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(disrow, discol);

  // Output of filament node positions
  if(!discret_->Comm().MyPID())
  {
    FILE* fp = NULL;
    // open file and append new data line
    fp = fopen(nodeposfilename.str().c_str(), "a");
    //defining temporary stringstream variable
    std::stringstream nodecoords;

    nodecoords<<istep<<"\t"<<-1e9<<"\t"<<-1e9<<"\t"<<-1e9<<std::endl;

    for(int i=subfil; i<(int)filaments.size(); i++)
    {
      for(int j=0; j<(int)filaments[i]->Nodes()->size(); j++)
      {
        int nodelid = discret_->NodeColMap()->LID((int)filaments[i]->Nodes()->at(j));
        DRT::Node* node = discret_->lColNode(nodelid);
        LINALG::Matrix<3, 1> coord;
        for(int dof=0; dof<3; dof++)
        {
          int dofgid = discret_->Dof(node)[dof];
          coord(dof) = node->X()[dof]+discol[discret_->DofColMap()->LID(dofgid)];
        }
        nodecoords<<i<<"\t"<<std::scientific<<std::setprecision(8)<<coord(0)<<"\t"<<coord(1)<<"\t"<<coord(2)<<std::endl;
      }
    }
    fputs(nodecoords.str().c_str(), fp);
    fclose(fp);
  }

  discret_->Comm().Barrier();

  //output of forces within crosslinkers
  FILE* fp02 = NULL;
  fp02 = fopen(forcefilename.str().c_str(), "a");
  //defining temporary stringstream variable
  std::stringstream linkerforces;

  Epetra_SerialDenseVector force;
  double eps = 0.0;
  LINALG::Matrix<3,1> matforces;

  for(int pid=0; pid<discret_->Comm().NumProc(); pid++)
  {
    if(pid==discret_->Comm().MyPID())
    {
      if(!pid && istep==statmechparams_.get<int>("OUTPUTINTERVALS",1))
      {
        linkerforces<<kt<<"\t"<<koff0<<"\t"<<delta<<"\t"<<dt<<"\t"<<statmechparams_.get<int>("OUTPUTINTERVALS",1)<<"\t"<<subfil<<"\t"<<filaments[0]->Nodes()->size()<<std::endl;

        linkerforces<<periodlength_->at(0)<<"\t"<<periodlength_->at(1)<<"\t"<<periodlength_->at(2)<<"\t"<<riseperbspot_->at(0)<<"\t";
        if((int)riseperbspot_->size()>1)
          linkerforces<<riseperbspot_->at(1)<<"\t";
        else
          linkerforces<<-1e9<<"\t";
        linkerforces<<statmechparams_.get<double>("R_LINK", 0.1)<<"\t"<<statmechparams_.get<double>("STROKEDISTANCE",0.005)<<std::endl;
      }
      for(int i=0; i<crosslink2element_->MyLength(); i++)
      {
        int rowlid = discret_->ElementRowMap()->LID((int)(*crosslink2element_)[i]);
        if((*crosslink2element_)[i]>-0.9 && rowlid!=-1)
        {
          DRT::Element* crosslinker = discret_->lRowElement(rowlid);
          int bspotgid0 = (int)(*crosslinkerbond_)[0][i];
          int bspotgid1 = (int)(*crosslinkerbond_)[1][i];

          const DRT::ElementType &eot = discret_->lRowElement(rowlid)->ElementType();
          if(eot == DRT::ELEMENTS::Beam3Type::Instance())
          {
            force.Resize(crosslinker->NumNode()*6);
            force = (dynamic_cast<DRT::ELEMENTS::Beam3*>(crosslinker))->InternalForceVector();
            matforces = ((dynamic_cast<DRT::ELEMENTS::Beam3*>(crosslinker))->MatForceGp());
            eps = (dynamic_cast<DRT::ELEMENTS::Beam3*>(crosslinker))->EpsilonSgn();
          }
          else if(eot == DRT::ELEMENTS::BeamCLType::Instance())
          {
            if(crosslinker->NumNode()!=4)
              dserror("Currently only implemented for BEAM3CL with four nodes.");
            force.Resize(crosslinker->NumNode()/2*6);
            force = (dynamic_cast<DRT::ELEMENTS::BeamCL*>(crosslinker))->InternalForceVector();
            matforces = ((dynamic_cast<DRT::ELEMENTS::BeamCL*>(crosslinker))->MatForceGp());
            eps = (dynamic_cast<DRT::ELEMENTS::BeamCL*>(crosslinker))->EpsilonSgn();
          }
          else
            dserror("Unknown crosslinker beam element!");

          LINALG::Matrix<3,1> f0;
          for(int j=0; j<(int)f0.M(); j++)
            f0(j) = force[j];

          double Fbspot0 = f0.Norm2();

          linkerforces<<bspotgid0<<"\t"<<bspotgid1<<"\t"<<std::scientific<<std::setprecision(8)<<Fbspot0<<"\t"<<matforces(0)<<"\t"<<matforces(1)<<"\t"<<matforces(2)<<"\t"<<eps<<std::endl;
        }
      }
    }
    discret_->Comm().Barrier();
  }

  fputs(linkerforces.str().c_str(), fp02);
  fclose(fp02);

  return;
}

/*------------------------------------------------------------------------------*
 | Output of Creep displacement                           (public) mueller 05/13|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::CreepDisplacementOutput(Epetra_Vector&      discol,
                                                        std::ostringstream& filename)
{
  FILE* fp = NULL;
  std::stringstream filecontent;
  fp = fopen(filename.str().c_str(), "w");

  int freedir = statmechparams_.get<int>("DBCDISPDIR", -1)-1;

  // write only the reference position and the displacement of the unconstrained direction
  for(int i=0; i<(int)dbcnodesets_[0].size(); i++)
  {
    DRT::Node* dbcnode = discret_->gNode(dbcnodesets_[0].at(i));
    std::vector<int> dofnode = discret_->Dof(dbcnode);
    int freedof = dofnode.at(freedir);
    filecontent << dbcnode->Id() <<"\t"<<std::scientific<<std::setprecision(8)<<dbcnode->X()[freedir]<<"\t"<<discol[freedof]<<std::endl;
  }

  fputs(filecontent.str().c_str(), fp);
  fclose(fp);

  return;
}

/*------------------------------------------------------------------------------*
| distribution of spherical coordinates                  (public) mueller 11/12|
*------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::BellsEquationOutput(const Epetra_Vector&      disrow,
                                                    const std::ostringstream& filename,
                                                    const double&             dt)
{
  if(!discret_->Comm().MyPID() && unbindingprobability_!=Teuchos::null)
  {
    FILE* fp = NULL;
    std::stringstream filecontent;
    fp = fopen(filename.str().c_str(), "a");

    double tol = 1e-9;
    for(int i=0; i<unbindingprobability_->MyLength(); i++)
      if((*unbindingprobability_)[0][i]>tol) // there are either two non-zero entries or none
        filecontent<<std::scientific<<std::setprecision(6)<<(*unbindingprobability_)[0][i]<<" "<<(*unbindingprobability_)[1][i]<<" "<<(*crosslinkerpositions_)[0][i]<<" "<<(*crosslinkerpositions_)[1][i]<<" "<<(*crosslinkerpositions_)[2][i]<<std::endl;

    fputs(filecontent.str().c_str(), fp);
    fclose(fp);
  }

  return;
}

/*------------------------------------------------------------------------------*                                                 |
 | Output of linker unbinding times, i.e. time that a linker stays bound        |
 |                                                      (private) mueller 01/13 |
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::LinkerUnbindingTimes(const std::ostringstream& filename)
{
  if(!discret_->Comm().MyPID())
  {
    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "a");

    std::stringstream unbindingtimes;

    for(int i=0; i<crosslinkunbindingtimes_->MyLength(); i++)
    {
      if((*crosslinkunbindingtimes_)[0][i]<-0.9 && (*crosslinkunbindingtimes_)[1][i]>=0.0)
        unbindingtimes << (*crosslinkunbindingtimes_)[1][i] <<std::endl;
    }
    fputs(unbindingtimes.str().c_str(), fp);
    fclose(fp);
  }
  return;
}
