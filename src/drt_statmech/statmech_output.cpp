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

#include <Teuchos_Time.hpp>

#include "statmech_manager.H"
#include "../drt_inpar/inpar_statmech.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_element.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_beamcontact/beam3contact_manager.H"
#include "../drt_beamcontact/beam3contact_octtree.H"

#include "../drt_beam3/beam3.H"
#include "../drt_beam3ii/beam3ii.H"
#include "../drt_truss3/truss3.H"
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
    string::iterator it = outputrootpath_.end();
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
        fprintf(fpnumbering, filecontent.str().c_str());
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
        vector<DRT::Condition*> pointneumannconditions(0);
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
        fprintf(fpnumbering, filecontent.str().c_str());
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
        fprintf(fpnumbering, filecontent.str().c_str());
        //close file
        fclose(fpnumbering);
      }
    }
    break;
    //simulating diffusion coefficient for anisotropic friction
    case INPAR::STATMECH::statout_anisotropic:
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
        fprintf(fpnumbering, filecontent.str().c_str());
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

        //filecontent << "Output for measurement of viscoelastic properties written by processor "<< discret_->Comm().MyPID() << endl;

        // move temporary stringstream to file and close it
        fprintf(fp, filecontent.str().c_str());
        fclose(fp);
      }
    }
      break;
    case INPAR::STATMECH::statout_networkcreep:
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

        //filecontent << "Output for measurement of viscoelastic properties written by processor "<< discret_->Comm().MyPID() << endl;

        // move temporary stringstream to file and close it
        fprintf(fp, filecontent.str().c_str());
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
    GmshOutput(dis,filename,istep);
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

            //defining temporary stringstream variable
            std::stringstream filecontent;
            filecontent << std::scientific << std::setprecision(15) << time - starttimeoutput_ << "  " << DeltaR2 << endl;

            // move temporary stringstream to file and close it
            fprintf(fp, filecontent.str().c_str());
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
      vector<DRT::Condition*> pointneumannconditions(0);
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
            //defining temporary stringstream variable
            std::stringstream filecontent;
            filecontent << std::scientific << std::setprecision(15) << time << "  "
                        << endtoend << " " << fint[num_dof - discret_->NumDof(
                           discret_->gNode(discret_->NumMyRowNodes() - 1))] << endl;
            // move temporary stringstream to file and close it
            fprintf(fp, filecontent.str().c_str());
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

          filecontent << endl;
          fprintf(fp, filecontent.str().c_str());
          fclose(fp);
        }

      }

    }
    break;
    //the following output allows for anisotropic diffusion simulation of a quasi stiff polymer
    case INPAR::STATMECH::statout_anisotropic:
    {
      //output in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps
      if ((time > starttime && fabs(time-starttime)>dt/1e4) && (istep % statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0))
      {


        //positions of first and last node in current time step (note: always 3D variables used; in case of 2D third compoenent just constantly zero)
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

          //defining temporary stringstream variable
          std::stringstream filecontent;
          filecontent << std::scientific << std::setprecision(15) << dt << " "<< sumsquareincmid_ << " "<< sumdispmiddle_.Norm2() * sumdispmiddle_.Norm2() <<endl;
//              << sumsquareincmid_ << " " << sumsquareincpar_ << " "
//              << sumsquareincort_ << " " << sumsquareincrot_ << " "
//              << sumdispmiddle_.Norm2() * sumdispmiddle_.Norm2() << " "
//              << sumrotmiddle_ * sumrotmiddle_ << endl;

          // move temporary stringstream to file and close it
          fprintf(fp, filecontent.str().c_str());
          fclose(fp);

          //write position
          LINALG::Matrix<3,1> midpoint(beginnew);
          midpoint += endnew;
          midpoint.Scale(0.5);
          std::ostringstream outputfilename2;
          outputfilename2 << outputrootpath_ << "/StatMechOutput/AnisotropicMidPosition" << outputfilenumber_ << ".dat";

          // open file and append new data line
          fp = fopen(outputfilename2.str().c_str(), "a");

          //defining temporary stringstream variable
          std::stringstream filecontent2;
          filecontent2 << std::scientific << std::setprecision(15) <<midpoint(0)<<" "<<midpoint(1)<<" "<<midpoint(2)<< endl;

          // move temporary stringstream to file and close it
          fprintf(fp, filecontent2.str().c_str());
          fclose(fp);
        }

        //new positions in this time step become old positions in last time step
        beginold_ = beginnew;
        endold_ = endnew;
      }
    }
    break;
    case INPAR::STATMECH::statout_viscoelasticity:
    {
      Teuchos::ParameterList sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
      int numstep = sdyn.get<int>("NUMSTEP", -1);
      //output in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps (or for the very last step)
      if ((time>=starttime && istep<numstep && (istep-istart_) % statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0) || fabs(time-starttime)<1e-8)
      {
        // name of file into which output is written
        std::ostringstream filename;
        filename << outputrootpath_ << "/StatMechOutput/ViscoElOutputProc.dat";
        ViscoelasticityOutput(time, dis, fint, filename);

        // additional Density-Density-Correlation Output
        if((istep-istart_) % (10*statmechparams_.get<int> ("OUTPUTINTERVALS", 1)) == 0)
        {
          std::ostringstream ddcorrfilename;
          ddcorrfilename << outputrootpath_ << "/StatMechOutput/DensityDensityCorrFunction_"<<std::setw(6) << std::setfill('0') << istep <<".dat";
          DDCorrOutput(dis, ddcorrfilename, istep, dt);
        }

        if(DRT::INPUT::IntegralValue<int>(statmechparams_,"FORCEDEPUNLINKING"))
        {
          std::ostringstream forcedepfilename;
          forcedepfilename << outputrootpath_ << "/StatMechOutput/UnbindingProbability_"<<std::setw(6) << std::setfill('0') << istep <<".dat";
          BellsEquationOutput(dis, forcedepfilename, dt);
        }
      }
    }
    break;
    case INPAR::STATMECH::statout_networkcreep:
    {
      Teuchos::ParameterList sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
      int numstep = sdyn.get<int>("NUMSTEP", -1);
      //output in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps (or for the very last step)
      if ((time>=starttime && istep<numstep && (istep-istart_) % statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0) || fabs(time-starttime)<1e-8)
      {
        // name of file into which output is written
        std::ostringstream filename;
        filename << outputrootpath_ << "/StatMechOutput/CreepForces.dat";
        ViscoelasticityOutput(time, dis, fint, filename);
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
      //output in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps
      if( ((time>=starttime && (istep-istart_) % statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0) || fabs(time-starttime)<1e-8) && beamcmanager!=Teuchos::null)
      {
        std::map<int, LINALG::Matrix<3, 1> > currentpositions;
        std::map<int, LINALG::Matrix<3, 1> > currentrotations;
        Epetra_Vector discol(*discret_->DofColMap(), true);
        LINALG::Export(dis, discol);
        GetNodePositionsFromDisVec(discol, currentpositions, currentrotations, true);
        beamcmanager->OcTree()->OctTreeSearch(currentpositions, istep);
        beamcmanager->ResetPairs();
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
      GmshOutput(dis,filename,istep,beamcmanager);
    else
      GmshOutput(dis,filename,istep);
  }

  return;
} // STATMECH::StatMechManager::Output()

/*----------------------------------------------------------------------*
 | writing Gmsh data for current step                 public)cyron 01/09|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GmshOutput(const Epetra_Vector& disrow,const std::ostringstream& filename, const int& step, RCP<CONTACT::Beam3cmanager> beamcmanager)
{
  /*the following method writes output data for Gmsh into file with name "filename"; all line elements are written;
   * the nodal displacements are handed over in the variable "dis"; note: in case of parallel computing only
   * processor 0 writes; it is assumed to have a fully overlapping column map and hence all the information about
   * all the nodal position; parallel output is now possible with the restriction that the nodes(processors) in question
   * are of the same machine*/
  GmshPrepareVisualization(disrow);

  //we need displacements also of ghost nodes and hence export displacement vector to column map format
  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(disrow, discol);

  // do output to file in c-style
  FILE* fp = NULL;

  //number of solid elements by which a round line is depicted
  const int nline = 16;

  // first processor starts by opening the file and writing the header, other processors have to wait
  if (!discret_->Comm().MyPID())
  {
    //open file to write output data into
    fp = fopen(filename.str().c_str(), "w");

    // write output to temporary stringstream;
    std::stringstream gmshfileheader;
    /*the beginning of the stream is "View \"" to indicate Gmsh that the following data is in order to create an image and
     * this command is followed by the name of that view displayed during it's shown in the video; in the following example
     * this name is for the 100th time step: Step00100; then the data to be presented within this view is written within { ... };
     * in the following example this data consists of scalar lines defined by the coordinates of their end points*/
    gmshfileheader <<"General.BackgroundGradient = 0;\n";
    gmshfileheader <<"View.LineType = 1;\n";
    gmshfileheader <<"View.LineWidth = 1.4;\n";
    gmshfileheader <<"View.PointType = 1;\n";
    gmshfileheader <<"View.PointSize = 3;\n";
    gmshfileheader <<"General.ColorScheme = 1;\n";
    gmshfileheader <<"General.Color.Background = {255,255,255};\n";
    gmshfileheader <<"General.Color.Foreground = {255,255,255};\n";
    //gmshfileheader <<"General.Color.BackgroundGradient = {128,147,255};\n";
    gmshfileheader <<"General.Color.Foreground = {85,85,85};\n";
    gmshfileheader <<"General.Color.Text = {0,0,0};\n";
    gmshfileheader <<"General.Color.Axes = {0,0,0};\n";
    gmshfileheader <<"General.Color.SmallAxes = {0,0,0};\n";
    gmshfileheader <<"General.Color.AmbientLight = {25,25,25};\n";
    gmshfileheader <<"General.Color.DiffuseLight = {255,255,255};\n";
    gmshfileheader <<"General.Color.SpecularLight = {255,255,255};\n";
    gmshfileheader <<"View.ColormapAlpha = 1;\n";
    gmshfileheader <<"View.ColormapAlphaPower = 0;\n";
    gmshfileheader <<"View.ColormapBeta = 0;\n";
    gmshfileheader <<"View.ColormapBias = 0;\n";
    gmshfileheader <<"View.ColormapCurvature = 0;\n";
    gmshfileheader <<"View.ColormapInvert = 0;\n";
    gmshfileheader <<"View.ColormapNumber = 2;\n";
    gmshfileheader <<"View.ColormapRotation = 0;\n";
    gmshfileheader <<"View.ColormapSwap = 0;\n";
    gmshfileheader <<"View.ColorTable = {Black,Yellow,Blue,Orange,Red,Cyan,Purple,Brown,Green};\n";
    gmshfileheader << "View \" Step " << step << " \" {" << endl;

    //write content into file and close it
    fprintf(fp, gmshfileheader.str().c_str());
    fclose(fp);
  }

  // wait for all processors to arrive at this point
  discret_->Comm().Barrier();

  // loop over the participating processors each of which appends its part of the output to one output file
  for (int proc = 0; proc < discret_->Comm().NumProc(); proc++)
  {
    if (discret_->Comm().MyPID() == proc)
    {
      //open file again to append ("a") output data into
      fp = fopen(filename.str().c_str(), "a");
      // write output to temporary stringstream;
      std::stringstream gmshfilecontent;

      //looping through all elements on the processor
      for (int i=0; i<discret_->NumMyColElements(); ++i)
      {
        // coordinates of nodes or binding spots (depending on element)
        LINALG::SerialDenseMatrix coord(3, 2, true);
        //getting pointer to current element
        DRT::Element* element = discret_->lColElement(i);
        const DRT::ElementType & eot = element->ElementType();

        // interpolated beam element
        if(eot==DRT::ELEMENTS::BeamCLType::Instance())
        {
          DRT::ELEMENTS::BeamCL* currele = NULL;
          currele = dynamic_cast<DRT::ELEMENTS::BeamCL*> (discret_->gElement(discret_->ElementColMap()->GID(i)));
          // Alternative way to get nodal positions
          DRT::Element* filelement = discret_->gElement(discret_->ElementColMap()->GID(i));
          DRT::Node* node0 = NULL;
          DRT::Node* node1 = NULL;
          std::vector<double> position(6);
          std::vector<int> dofnode0;
          std::vector<int> dofnode1;

          std::vector<LINALG::Matrix<1,2> > Ibp(2);
          for(int filament=0; filament<2; filament++)
             DRT::UTILS::shape_function_1D(Ibp[filament],currele->MyBindingPosition()[filament],currele->Shape());

          // determine positions of the interpolated nodes
          for(int ifil=0; ifil<coord.N(); ifil++)
          {
            node0 = discret_->gNode(filelement->NodeIds()[0+2*ifil]);
            node1 = discret_->gNode(filelement->NodeIds()[1+2*ifil]);

            dofnode0 = discret_->Dof(node0);
            dofnode1 = discret_->Dof(node1);
            position[0] = node0->X()[0] + discol[discret_->DofColMap()->LID(dofnode0[0])];
            position[1] = node0->X()[1] + discol[discret_->DofColMap()->LID(dofnode0[1])];
            position[2] = node0->X()[2] + discol[discret_->DofColMap()->LID(dofnode0[2])];
            position[3] = node1->X()[0] + discol[discret_->DofColMap()->LID(dofnode1[0])];
            position[4] = node1->X()[1] + discol[discret_->DofColMap()->LID(dofnode1[1])];
            position[5] = node1->X()[2] + discol[discret_->DofColMap()->LID(dofnode1[2])];

            //shift nodes back in case they have been shifted due to PBC
            UnshiftPositions(position);
            for(int id=0;id<3;id++)
             coord(id, ifil) = Ibp[ifil](0)*position[id]+Ibp[ifil](1)*position[id+3];
            //shift Crosslinker intperpolated nodal positions into Periodic Boundary Domain
            for(int fil=0;fil<2;fil++)
              for(int j=0;j<3;j++)
              {  if(coord(j,fil) > periodlength_->at(j))
                  coord(j,fil) -= periodlength_->at(j)*floor(coord(j,fil)/periodlength_->at(j));
                 if(coord(j,fil) < 0.0)
                   coord(j,fil) -= periodlength_->at(j)*floor(coord(j,fil)/periodlength_->at(j));
              }
          }
        }
        else  // standard beam or truss element
        {
          for (int id=0; id<3; id++)
            for (int jd=0; jd<element->NumNode(); jd++)
            {
              double referenceposition = ((element->Nodes())[jd])->X()[id];
              vector<int> dofnode = discret_->Dof((element->Nodes())[jd]);
              double displacement = discol[discret_->DofColMap()->LID(dofnode[id])];
              coord(id, jd) = referenceposition + displacement;
            }
        }

        //declaring variable for color of elements
        double color;

        //apply different colors for elements representing filaments and those representing dynamics crosslinkers
        if (element->Id() < basisnodes_)
          color = 1.0;
        else
          color = 0.5;

        // highlight contacting elements
        if(beamcmanager!=Teuchos::null)
          for(int j=0; j<(int)(beamcmanager->Pairs()).size(); j++)
            if(beamcmanager->Pairs()[j]->GetContactFlag() && (element->Id()==(beamcmanager->Pairs())[j]->Element1()->Id() || element->Id()==(beamcmanager->Pairs())[j]->Element2()->Id()))
              color = 1.0; //0.375;

        //if no periodic boundary conditions are to be applied, we just plot the current element
        if (periodlength_->at(0) == 0.0)
        {
          // check whether the kinked visualization is to be applied
          bool kinked = CheckForKinkedVisual(element->Id());
          if (eot == DRT::ELEMENTS::Beam3Type::Instance() ||
              eot==DRT::ELEMENTS::Beam3iiType::Instance() ||
              eot==DRT::ELEMENTS::BeamCLType::Instance() )
          {
            if (!kinked)
            {
              int numnode = element->NumNode();
              if(eot==DRT::ELEMENTS::BeamCLType::Instance())
                numnode = 2;
              for (int j=0; j<numnode - 1; j++)
              {
                //define output coordinates
                LINALG::SerialDenseMatrix coordout(3,2);
                for(int m=0; m<coordout.M(); m++)
                  for(int n=0; n<coordout.N(); n++)
                    coordout(m,n)=coord(m,j+n);

                 GmshWedge(nline,coordout,element,gmshfilecontent,color);
              }
            }
            else
              GmshKinkedVisual(coord, 0.875, element->Id(), gmshfilecontent);
          }
          else if (eot == DRT::ELEMENTS::Truss3Type::Instance())
          {
            if (!kinked)
            {
              for (int j=0; j<element->NumNode()-1; j++)
              {
                gmshfilecontent << "SL(" << std::scientific;
                gmshfilecontent << coord(0, j) << "," << coord(1, j) << ","<< coord(2, j) << ","
                                << coord(0, j + 1) << "," << coord(1,j + 1) << "," << coord(2, j + 1);
                gmshfilecontent << ")" << "{" << std::scientific << color << ","<< color << "};" << endl;
              }
            }
            else
              GmshKinkedVisual(coord, 0.875, element->Id(), gmshfilecontent);
          }
          else if (eot == DRT::ELEMENTS::Torsion3Type::Instance())
          {
            double beadcolor = 0.75;
            for (int j=0; j<element->NumNode(); j++)
            {
              gmshfilecontent << "SP(" << std::scientific;
              gmshfilecontent << coord(0, j) << "," << coord(1, j) << ","<< coord(2, j);
              gmshfilecontent << ")" << "{" << std::scientific << beadcolor << ","<< beadcolor << "};" << endl;
            }
          }
          else
          {
            //nothing!
          }
        }
        //in case of periodic boundary conditions we have to take care to plot correctly an element broken at some boundary plane
        else
          GmshOutputPeriodicBoundary(coord, color, gmshfilecontent,element->Id(),false);
      }
      //write content into file and close it (this way we make sure that the output is written serially)
      fprintf(fp, gmshfilecontent.str().c_str());
      fclose(fp);
    }
    discret_->Comm().Barrier();
  }
  // plot the periodic boundary box
  LINALG::Matrix<3,1> center;
  for(int i=0; i<(int)center.M(); i++)
    center(i) = periodlength_->at(i)/2.0;
  GmshOutputBox(0.0, &center, *periodlength_, &filename);
  // plot crosslink molecule diffusion and (partial) bonding
  GmshOutputCrosslinkDiffusion(0.125, &filename, disrow);

  // finish data section of this view by closing curly brackets
  if (discret_->Comm().MyPID() == 0)
  {
    fp = fopen(filename.str().c_str(), "a");
    std::stringstream gmshfileend;

    if(DRT::INPUT::IntegralValue<int>(statmechparams_,"GMSHNETSTRUCT") &&
       DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_structanaly)
    {
      // plot the rotated material triad (with axis length 0.5)
      for(int i=0; i<trafo_->M(); i++)
      {
        LINALG::SerialDenseMatrix coord(3,2);
        for(int j=0;j<coord.M(); j++)
        {
          coord(j,0) = cog_(j);
          coord(j,1) = cog_(j)+0.5*(*trafo_)(i,j);
        }
        GmshWedge(1,coord,discret_->lRowElement(0),gmshfileend,0.0,true,true);
      }
      // plot the cog
      std::vector<double> dimension(3,0.05);
      GmshOutputBox(0.75, &cog_, dimension, &filename, false);
      gmshfileend << "SP(" << std::scientific;
      gmshfileend << cog_(0)<<","<<cog_(1)<<","<<cog_(2)<<"){" << std::scientific << 0.75 << ","<< 0.75 <<"};"<<endl;

      // gmsh output of detected network structure volume
      int nline = 16;
      if(step>0)
        GmshNetworkStructVolume(nline, gmshfileend, 0.875);
    }
    // add black dot for correct coloring...
    if(periodlength_->at(0)==0.0)
      gmshfileend << "SP("<<periodlength_->at(0)<<",0,"<<periodlength_->at(2)<<"){0,0};"<<endl;
    gmshfileend << "};" << endl;
    fprintf(fp, gmshfileend.str().c_str());
    fclose(fp);
  }

  // return simultaneously (not sure if really needed)
  discret_->Comm().Barrier();

  return;
} // STATMECH::StatMechManager::GmshOutput()


/*----------------------------------------------------------------------*
 | gmsh output data in case of periodic boundary conditions             |
 |                                                    public)cyron 02/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GmshOutputPeriodicBoundary(const LINALG::SerialDenseMatrix& coord, const double& color, std::stringstream& gmshfilecontent, int eleid, bool ignoreeleid)
{
  //number of solid elements by which a round line is depicted
  const int nline = 16;

  //number of spatial dimensions
  const int ndim = 3;
  // get Element Type of the first Element to determine the graphics output
  DRT::Element* element = discret_->gElement(eleid);

  bool dotline = false;
  bool kinked = false;

  if (ignoreeleid)
    dotline = true;
  else
  {
    // draw colored lines between two nodes of a beam or a truss element (meant for filaments/crosslinks/springs)
    const DRT::ElementType & eot = element->ElementType();
    if(element->ElementType().Name()=="Beam3iiType")
      dotline = eot==DRT::ELEMENTS::Beam3iiType::Instance();
    if (element->ElementType().Name() == "Beam3Type")
      dotline = eot == DRT::ELEMENTS::Beam3Type::Instance();
    else if(element->ElementType().Name()=="BeamCLType")
      dotline = eot==DRT::ELEMENTS::BeamCLType::Instance();
    else if (element->ElementType().Name() == "Truss3Type")
      dotline = dotline or eot == DRT::ELEMENTS::Truss3Type::Instance();
    // draw spheres at node positions ("beads" of the bead spring model)
    else if (eot == DRT::ELEMENTS::Torsion3Type::Instance())
    {
      double beadcolor = 0.75;
      for (int i=0; i<element->NumNode(); i++)
      {
        //writing element by nodal coordinates as a sphere
        gmshfilecontent << "SP(" << std::scientific;
        gmshfilecontent << coord(0,i) << "," << coord(1, i) << "," << coord(2,i);
        gmshfilecontent << ")" << "{" << std::scientific << beadcolor << ","<< beadcolor << "};" << endl;
      }
    }
    /* determine whether crosslink connects two filaments or occupies two binding spots on the same filament;
     * variable triggers different visualizations.*/
    kinked = CheckForKinkedVisual(element->Id());
  }

  if (dotline)
  {
    /*detect and save in vector "cut", at which boundaries the element is broken due to periodic boundary conditions;
     * the entries of cut have the following meaning: 0: element not broken in respective coordinate direction, 1:
     * element broken in respective coordinate direction (node 0 close to zero boundary and node 1 close to boundary
     * at periodlength_);  2: element broken in respective coordinate direction (node 1 close to zero boundary and node
     * 0 close to boundary at periodlength_);*/
    LINALG::SerialDenseMatrix cut= LINALG::SerialDenseMatrix(3, coord.N()-1, true);

    /* "coord" currently holds the shifted set of coordinates.
     * In order to determine the correct vector "dir" of the visualization at the boundaries,
     * a copy of "coord" with adjustments in the proper places is introduced*/
    for (int i=0; i<cut.N(); i++)
    {
      LINALG::SerialDenseMatrix unshift(3,2);

      int numshifts = 0;
      int shiftdof = -1;
      for (int dof=0; dof<ndim; dof++)
      {
        // initialize unshift with coord values
        unshift(dof,0) = coord(dof,i);
        unshift(dof,1) = coord(dof,i+1);
        if (fabs(coord(dof,i+1)-periodlength_->at(dof)-coord(dof,i)) < fabs(coord(dof,i+1) - coord(dof,i)))
        {
          cut(dof, i) = 1.0;
          shiftdof = dof;
          unshift(dof,1) -= periodlength_->at(dof);
          numshifts++;
        }
        if (fabs(coord(dof,1)+periodlength_->at(dof) - coord(dof,i)) < fabs(coord(dof,i+1)-coord(dof,i)))
        {
          cut(dof,i) = 2.0;
          shiftdof = dof;
          unshift(dof,1) += periodlength_->at(dof);
          numshifts++;
        }
      }

      // write special output for broken elements
      if (numshifts > 0)
      {
        // directional vector
        LINALG::Matrix<3, 1> dir;
        for (int dof = 0; dof < ndim; dof++)
          dir(dof) = unshift(dof,1) - unshift(dof,0);
        dir.Scale(1.0/dir.Norm2());

        /* determine the intersection points of the line through unshift(:,0) and direction dir with the faces of the boundary cube
         * and sort them by distance. Thus, we obtain an order by which we have to shift the element back into the cube so that all
         * segments that arise by multiple shifts remain within the volume (see my notes on 6/12/2011).*/
        LINALG::Matrix<3,2> lambdaorder;
        lambdaorder.PutScalar(1e6);
        // collect lambdas
        for(int dof=0; dof<(int)lambdaorder.M(); dof++)
        {
          switch((int)cut(dof,i))
          {
            case 1:
            {
              lambdaorder(dof,0) = -coord(dof, i) / dir(dof);
              lambdaorder(dof,1) = dof;
            }
            break;
            case 2:
            {
              lambdaorder(dof,0) = (periodlength_->at(dof) - coord(dof,i)) / dir(dof);
              lambdaorder(dof,1) = dof;
            }
            break;
            default:
            {
              lambdaorder(dof,1) = dof;
            }
            break;
          }
        }
        // sort the lambdas (ascending values) and indices accordingly
        // in case of multiple shifts
        if(numshifts>1)
        {
          for(int j=0; j<(int)lambdaorder.M()-1; j++)
            for(int k=j+1; k<(int)lambdaorder.M(); k++)
              if(lambdaorder(k,0)<lambdaorder(j,0))
              {
                double temp = lambdaorder(j,0);
                int tempindex = (int)lambdaorder(j,1);
                lambdaorder(j,0) = lambdaorder(k,0);
                lambdaorder(j,1) = lambdaorder(k,1);
                lambdaorder(k,0) = temp;
                lambdaorder(k,1) = tempindex;
              }
        }
        else  // for a single shift (the majority of broken elements), just put the index and the lambda of the broken dof in front
          for(int n=0; n<(int)lambdaorder.N(); n++)
          {
            double tmp = lambdaorder(shiftdof,n);
            lambdaorder(0,n) = tmp;
          }

        // calculate segment lambdas
        for(int dof=numshifts-1; dof>0; dof--)
          lambdaorder(dof,0) -= lambdaorder(dof-1,0);

        // the idea is to gradually shift the matrix "unshift" back into the volume and, while doing so,
        // calculate the segments except for the last one
        // determine closest boundary component-wise
        for(int shift=0; shift<numshifts; shift++)
        {
          //second point
          for(int j=0 ;j<unshift.M(); j++)
            unshift(j,1) = unshift(j,0) + lambdaorder(shift,0)*dir(j);
          //GmshOutput
          if(!ignoreeleid)
            GmshWedge(nline,unshift,element,gmshfilecontent,color);
          else
            GmshWedge(nline,unshift,element,gmshfilecontent,color,true);

          int currshift = (int)lambdaorder(shift,1);
          // shift the coordinates of the second point
          if(cut(currshift,i)==1.0)
            unshift(currshift,1) += periodlength_->at(currshift);
          else if(cut(currshift,i)==2.0)
            unshift(currshift,1) -= periodlength_->at(currshift);
          // make second point the first and calculate new second point in the next iteration!
          for(int j=0; j<unshift.M(); j++)
            unshift(j,0) = unshift(j,1);
        }
        //last segment
        // the last segment
        for(int dof=0; dof<unshift.M(); dof++)
          unshift(dof,1) = coord(dof,i+1);
        if(!ignoreeleid)
          GmshWedge(nline,unshift,element,gmshfilecontent,color);
        else
          GmshWedge(nline,unshift,element,gmshfilecontent,color,true);
      }
      else // output for continuous elements
      {
        if (!kinked)
        {
          if(!ignoreeleid)
            GmshWedge(nline,coord,element,gmshfilecontent,color);
          else
            GmshWedge(nline,coord,element,gmshfilecontent,color,true);
        }
        else
          GmshKinkedVisual(coord, 0.875, element->Id(), gmshfilecontent);
      }
    }
  }
  return;
} // STATMECH::StatMechManager::GmshOutputPeriodicBoundary()

/*----------------------------------------------------------------------*
 | plot the periodic boundary box                  (public) mueller 7/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GmshOutputBox(double boundarycolor, LINALG::Matrix<3,1>* boxcenter,
                                              std::vector<double>& dimension,
                                              const std::ostringstream *filename,
                                              bool barrier)
{
  // plot the periodic box in case of periodic boundary conditions (first processor)
  if (periodlength_->at(0) > 0.0 && discret_->Comm().MyPID() == 0)
  {
    FILE *fp = fopen(filename->str().c_str(), "a");
    std::stringstream gmshfilefooter;
    // get current period length

    double xmin = (*boxcenter)(0)-dimension[0]/2.0;
    double xmax = (*boxcenter)(0)+dimension[0]/2.0;
    double ymin = (*boxcenter)(1)-dimension[1]/2.0;
    double ymax = (*boxcenter)(1)+dimension[1]/2.0;
    double zmin = (*boxcenter)(2)-dimension[2]/2.0;
    double zmax = (*boxcenter)(2)+dimension[2]/2.0;

    // define boundary lines
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmin << "," << ymin << "," << zmin << "," << xmax << ","<< ymin << "," << zmin;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
    // line 2
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmax << "," << ymin << "," << zmin << "," << xmax << "," << ymax << "," << zmin;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
    // line 3
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmax << "," << ymax << "," << zmin << "," << xmax << "," << ymax << "," << zmax;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
    // line 4
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmax << "," << ymax << "," << zmax << "," << xmin << "," << ymax << "," << zmax;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
    // line 5
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmin << "," << ymax << "," << zmax << "," << xmin << "," << ymin << "," << zmax;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
    // line 6
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmin << "," << ymin << "," << zmax << "," << xmin << ","<< ymin << "," << zmin;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
    // line 7
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmin << "," << ymin << "," << zmin << "," << xmin << ","<< ymax << "," << zmin;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
    // line 8
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmin << "," << ymax << "," << zmin << "," << xmax << "," << ymax << "," << zmin;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
    // line 9
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmin << "," << ymax << "," << zmin << "," << xmin << "," << ymax << "," << zmax;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
    // line 10
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmax << "," << ymin << "," << zmin << "," << xmax << "," << ymin << "," << zmax;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
    // line 11
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmax << "," << ymin << "," << zmax << "," << xmax << "," << ymax<< "," << zmax;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
    // line 12
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmax << "," << ymin << "," << zmax << "," << xmin << "," << ymin<< "," << zmax;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << endl;

    fprintf(fp, gmshfilefooter.str().c_str());
    fclose(fp);
  }
  // wait for Proc 0 to catch up to the others
  if(barrier)
    discret_->Comm().Barrier();
}// STATMECH::StatMechManager::GmshOutputBoundaryBox

/*----------------------------------------------------------------------*
 | Gmsh output for crosslink molecule diffusion    (public) mueller 7/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GmshOutputCrosslinkDiffusion(double color, const std::ostringstream *filename, const Epetra_Vector& disrow)
{
  // export row displacement to column map format
  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(disrow, discol);

  //In case of the 4-noded Crosslinker Element we need the following vectors
  Teuchos::RCP<Epetra_MultiVector> bspotpositions = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));
  Teuchos::RCP<Epetra_MultiVector> bspotrotations = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));


  GetBindingSpotPositions(discol, bspotpositions, bspotrotations);


  if (discret_->Comm().MyPID() == 0)
  {
    FILE *fp = fopen(filename->str().c_str(), "a");
    //special visualization for crosslink molecules with one/two bond(s); going through the Procs

    //fp = fopen(filename->str().c_str(), "a");
    std::stringstream gmshfilebonds;

    // first, just update positions: redundant information on all procs
    for (int i=0; i<numbond_->MyLength(); i++)
    {
      switch ((int) (*numbond_)[i])
      {
        // free linkers
        case 0:
        {
//          double beadcolor = 5*color;
//          //writing element by nodal coordinates as a sphere
//          gmshfilebonds << "SP(" << scientific;
//          gmshfilebonds<< (*visualizepositions_)[0][i]<< "," << (*visualizepositions_)[1][i] << "," << (*visualizepositions_)[2][i];
//          gmshfilebonds << ")" << "{" << scientific << beadcolor << "," << beadcolor << "};" << endl;
        }
        break;
        // crosslink molecule with one bond
        case 1:
        {
          // determine position of bpspotLID entry
          int bspotLID = bspotcolmap_->LID((int)(*crosslinkerbond_)[0][i]);
          if((*crosslinkerbond_)[0][i]<-0.9)
            bspotLID = bspotcolmap_->LID((int)(*crosslinkerbond_)[1][i]);

          LINALG::SerialDenseMatrix coord(3, 2, true);
          double length = 0.0;
          for (int j=0; j<coord.M(); j++)
          {
            coord(j, 0) = (*bspotpositions)[j][bspotLID];
            coord(j, 1) = (*visualizepositions_)[j][i];
            length += (coord(j,1)-coord(j,0))*(coord(j,1)-coord(j,0));
          }
          double beadcolor = 2*color; //blue
          // in case of periodic boundary conditions
          if (periodlength_->at(0) > 0.0)
          {
            // get arbitrary element (we just need it to properly visualize)
            DRT::Element* tmpelement=discret_->lRowElement(0);
            GmshOutputPeriodicBoundary(coord, 2*color, gmshfilebonds, tmpelement->Id(), true);
            // visualization of "real" crosslink molecule positions (attention: shifted by one step since StatMechUpdate is called before the Newton scheme)
            //beadcolor = 0.0; //black
            //gmshfilebonds << "SP(" << scientific;
            //gmshfilebonds << (*crosslinkerpositions_)[0][i] << ","<< (*crosslinkerpositions_)[1][i] << ","<< (*crosslinkerpositions_)[2][i];
            //gmshfilebonds << ")" << "{" << scientific << beadcolor << ","<< beadcolor << "};" << endl;
          }
          else
          {
            gmshfilebonds << "SL(" << scientific;
            gmshfilebonds << coord(0, 0) << "," << coord(1, 0) << ","<< coord(2, 0) << "," << coord(0, 1) << "," << coord(1, 1)<< "," << coord(2, 1);
            gmshfilebonds << ")" << "{" << scientific << 2*color << ","<< 2*color << "};" << endl;
            gmshfilebonds << "SP(" << scientific;
            gmshfilebonds << coord(0, 1) << "," << coord(1, 1) << ","<< coord(2, 1);
            gmshfilebonds << ")" << "{" << scientific << beadcolor << ","<< beadcolor << "};" << endl;
          }
        }
        break;
        // crosslinker element: crosslink molecule (representation) position (Proc 0 only)
        case 2:
        {
          // actual crosslinker element connecting two filaments (self-binding kinked crosslinkers are visualized in GmshKinkedVisual())
          if((*searchforneighbours_)[i] > 0.9)
          {
            if((*crosslinkonsamefilament_)[i] < 0.1)
            {
              double beadcolor = 4* color ; //red
              gmshfilebonds << "SP(" << scientific;
              gmshfilebonds << (*visualizepositions_)[0][i] << ","<< (*visualizepositions_)[1][i] << ","<< (*visualizepositions_)[2][i];
              gmshfilebonds << ")" << "{" << scientific << beadcolor << ","<< beadcolor << "};" << endl;
            }
          }
          else  // passive crosslink molecule
          {
            if(DRT::INPUT::IntegralValue<int>(statmechparams_, "INTERNODALBSPOTS") && statmechparams_.get<double>("K_ON_SELF",0.0)>0.0)
              dserror("Passified links not implemented for BeamCL elements!");
            // determine position of nodeGID entry
            int bspotLID = bspotcolmap_->LID((int)(*crosslinkerbond_)[0][i]);
            if((*crosslinkerbond_)[0][i]<-0.9)
              bspotLID = bspotcolmap_->LID((int)(*crosslinkerbond_)[1][i]);

            DRT::Node *node = discret_->lColNode(bspotLID);
            LINALG::SerialDenseMatrix coord(3, 2, true);
            for (int j=0; j<coord.M(); j++)
            {
              int dofgid = discret_->Dof(node)[j];
              coord(j, 0) = node->X()[j] + discol[dofgid];
              coord(j, 1) = (*visualizepositions_)[j][i];
            }

            double beadcolor = 3*color;
            // in case of periodic boundary conditions
            if (periodlength_->at(0) > 0.0)
            {
              // get arbitrary element (we just need it to properly visualize)
              DRT::Element* tmpelement=discret_->lRowElement(0);
              GmshOutputPeriodicBoundary(coord, 3*color, gmshfilebonds, tmpelement->Id(), true);
              // visualization of "real" crosslink molecule positions (attention: shifted by one step since StatMechUpdate is called before the Newton scheme)
              //beadcolor = 0.0; //black
              //gmshfilebonds << "SP(" << scientific;
              //gmshfilebonds << (*crosslinkerpositions_)[0][i] << ","<< (*crosslinkerpositions_)[1][i] << ","<< (*crosslinkerpositions_)[2][i];
              //gmshfilebonds << ")" << "{" << scientific << beadcolor << ","<< beadcolor << "};" << endl;
            }
            else
            {
              gmshfilebonds << "SL(" << scientific;
              gmshfilebonds << coord(0, 0) << "," << coord(1, 0) << ","<< coord(2, 0) << "," << coord(0, 1) << "," << coord(1, 1)<< "," << coord(2, 1);
              gmshfilebonds << ")" << "{" << scientific << 3*color << ","<< 3*color << "};" << endl;
              gmshfilebonds << "SP(" << scientific;
              gmshfilebonds << coord(0, 1) << "," << coord(1, 1) << ","<< coord(2, 1);
              gmshfilebonds << ")" << "{" << scientific << beadcolor << ","<< beadcolor << "};" << endl;
            }
          }
        }
          break;
        default:
          continue;
      }
    }
    fprintf(fp, gmshfilebonds.str().c_str());
    fclose(fp);
  }

  discret_->Comm().Barrier();
}// GmshOutputCrosslinkDiffusion

/*----------------------------------------------------------------------*
 | Special Gmsh output for crosslinkers occupying two binding spots on  |
 | the same filament                              (public) mueller 07/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GmshKinkedVisual(const LINALG::SerialDenseMatrix& coord, const double& color, int eleid, std::stringstream& gmshfilecontent)
{
  /* We need a third point in order to visualize the crosslinker.
   * It marks the location of the kink.
   */
  std::vector<double> thirdpoint(3,0.0);
  // get the element
  DRT::Element* element = discret_->gElement(eleid);

  // calculate unit tangent
  LINALG::Matrix<3,1> t;
  for (int j=0; j<coord.M(); j++)
    t(j) = coord(j, coord.N() - 1) - coord(j, 0);
  t.Scale(1.0/t.Norm2());

  // calculate normal via cross product: [0 0 1]x[tx ty tz]
  LINALG::Matrix<3,1> n;
  n(0) = -t(1);
  n(1) = t(0);
  // norm it since the cross product does not keep the length
  n.Scale(1.0/n.Norm2());

  // by modulo operation involving the node IDs
  double alpha = fmod((double) (element->Nodes()[element->NumNode() - 1]->Id() + element->Nodes()[0]->Id()), 2*M_PI);

  // rotate the normal by alpha
  RotationAroundFixedAxis(t,n,alpha);

  // calculation of the third point lying in the direction of the rotated normal
  // height of the third point above the filament
  double h = 0.33 * (statmechparams_.get<double> ("R_LINK", 0.0) + statmechparams_.get<double> ("DeltaR_LINK", 0.0));
  for (int j=0; j<3; j++)
    thirdpoint.at(j) = (coord(j, 0) + coord(j, element->NumNode() - 1)) / 2.0 + h * n(j);

  gmshfilecontent << "SL(" << std::scientific << coord(0,0) << ","<< coord(1,0) << "," << coord(2,0) << ","
                  << thirdpoint.at(0) << ","<< thirdpoint.at(1) << "," << thirdpoint.at(2) << ")"
                  << "{" << std::scientific<< color << "," << color << "};" << endl;
  gmshfilecontent << "SL(" << std::scientific << thirdpoint.at(0) << "," << thirdpoint.at(1)<< "," << thirdpoint.at(2) << ","
                  << coord(0, 1) << "," << coord(1,1) << "," << coord(2,1) << ")"
                  << "{" << std::scientific<< color << "," << color << "};" << endl;
  gmshfilecontent << "SP(" << std::scientific
                  << thirdpoint.at(0) << "," << thirdpoint.at(1) << ","<< thirdpoint.at(2)
                  << ")" << "{" << std::scientific << color << ","<< color << "};" << endl;
  return;
}// STATMECH::StatMechManager::GmshKinkedVisual

/*----------------------------------------------------------------------*
 | prepare visualization vector for Gmsh Output  (private) mueller 08/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GmshPrepareVisualization(const Epetra_Vector& dis)
{
  double ronebond = statmechparams_.get<double> ("R_LINK", 0.0) / 2.0;

  // column map displacement and displacement increment

  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(dis, discol);

  //In case of the 4-noded Crosslinker Element we need the following vectors
  Teuchos::RCP<Epetra_MultiVector> bspotpositions = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));
  Teuchos::RCP<Epetra_MultiVector> bspotrotations = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));
  GetBindingSpotPositions(discol, bspotpositions, bspotrotations);

  // get binding spot triads
  Teuchos::RCP<Epetra_MultiVector> bspottriadscol = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,4,true));
  GetBindingSpotTriads(bspotrotations, bspottriadscol);

  if (discret_->Comm().MyPID() == 0)
  {
    for (int i=0; i<numbond_->MyLength(); i++)
    {
      switch((int)(*numbond_)[i])
      {
        // diffusion
        case 0:
        {
          for (int j=0; j<visualizepositions_->NumVectors(); j++)
            (*visualizepositions_)[j][i] = (*crosslinkerpositions_)[j][i];
        }
        break;
        // one bonded crosslink molecule
        case 1:
        {
          // determine position of bspotLID entry
          int bspotLID = bspotcolmap_->LID((int)(*crosslinkerbond_)[0][i]);
          if((*crosslinkerbond_)[0][i]<-0.9)
            bspotLID = bspotcolmap_->LID((int)(*crosslinkerbond_)[1][i]);

          //calculate unit tangent
          LINALG::Matrix<3,1> bspotpos0;
          for (int j=0; j<3; j++)
            bspotpos0(j,0) = (*bspotpositions)[j][bspotLID];

          // first and second vector of the nodal triad
          LINALG::Matrix<3, 1> tangent;
          LINALG::Matrix<3, 1> normal;
          // rotation angle
          double alpha = 0.0;

          LINALG::Matrix<3,3> bspottriad;
          // auxiliary variable for storing a triad in quaternion form
          LINALG::Matrix<4, 1> qnode;
          // triad of node on first filament which is affected by the new crosslinker
          for (int j=0; j<4; j++)
            qnode(j) = (*bspottriadscol)[j][bspotLID];
          LARGEROTATIONS::quaterniontotriad(qnode, bspottriad);

          for(int j=0; j<(int)tangent.M(); j++)
          {
            tangent(j) = bspottriad(j,0);
            normal(j) = bspottriad(j,1);
          }
        // rotation angle
        if(DRT::INPUT::IntegralValue<int>(statmechparams_, "HELICALBINDINGSTRUCT"))
         alpha = (*bspotorientations_)[bspotLID];
        else
          alpha = fmod((double) i, 2*M_PI);

          // rotate the normal by alpha and store the new direction into the same Matrix ("normal")
          RotationAroundFixedAxis(tangent,normal,alpha);

          // calculation of the visualized point lying in the direction of the rotated normal
          for (int j=0; j<visualizepositions_->NumVectors(); j++)
            (*visualizepositions_)[j][i] = bspotpos0(j,0) + ronebond*normal(j,0);

        }

        break;
        case 2:
        {
          // actual crosslinker element (not kinked)
          if((*searchforneighbours_)[i]>0.9)
          {
            // loop over filament node components
            for (int j=0; j<visualizepositions_->NumVectors(); j++)
            {
              std::vector<double> dofbspotpositions(crosslinkerbond_->NumVectors(), 0.0);
              // loop over filament node GIDs
              for (int k=0; k<crosslinkerbond_->NumVectors(); k++)
              {
                int bspotGID = (int) (*crosslinkerbond_)[k][i];
                dofbspotpositions.at(k) = (*bspotpositions)[j][bspotGID];
              }
              /* Check if the crosslinker element is broken/discontinuous; if so, reposition the second nodal value.
               * It does not matter which value is shifted as long the shift occurs in a consistent way.*/
              if (periodlength_->at(j) > 0.0)
              {
                (*visualizepositions_)[j][i] = dofbspotpositions.at(0);
                for (int k=0; k<1; k++)
                {
                  // shift position if it is found to be outside the boundary box
                  if (fabs(dofbspotpositions.at(k+1) - periodlength_->at(j) - dofbspotpositions.at(k))< fabs(dofbspotpositions.at(k+1) - dofbspotpositions.at(k)))
                    dofbspotpositions.at(k+1) -= periodlength_->at(j);
                  if (fabs(dofbspotpositions.at(k+1) + periodlength_->at(j) - dofbspotpositions.at(k))< fabs(dofbspotpositions.at(k+1) - dofbspotpositions.at(k)))
                    dofbspotpositions.at(k+1) += periodlength_->at(j);

                  (*visualizepositions_)[j][i] += dofbspotpositions.at(k+1);
                }
                // new crosslink molecule position
                (*visualizepositions_)[j][i] /= 2.0;
              }
              else
                (*visualizepositions_)[j][i] /= 2.0;
            }
          }
          else  // passive crosslink molecule
          {
            if(DRT::INPUT::IntegralValue<int>(statmechparams_, "INTERNODALBSPOTS") && statmechparams_.get<double>("K_ON_SELF",0.0)>0.0)
              dserror("Passified links not implemented for BeamCL elements!");
            // determine position of bspotLID entry
            int bspotLID = bspotcolmap_->LID((int)(*crosslinkerbond_)[0][i]);
            if ((*crosslinkerbond_)[0][i] < -0.9)
              bspotLID = bspotcolmap_->LID((int)(*crosslinkerbond_)[1][i]);

            const DRT::Node *node0 = discret_->lColNode(bspotLID);
            //calculate unit tangent
            LINALG::Matrix<3,1> nodepos0;

            for (int j=0; j<3; j++)
            {
              int dofgid0 = discret_->Dof(node0)[j];
              nodepos0(j,0) = node0->X()[j] + discol[discret_->DofColMap()->LID(dofgid0)];
            }

            // first and second vector of the nodal triad
            LINALG::Matrix<3, 1> tangent;
            LINALG::Matrix<3, 1> normal;
            // rotation angle
            double alpha = 0.0;

            // determine tangent and normal direction when helical binding spot geometry is applied
            if(DRT::INPUT::IntegralValue<int>(statmechparams_, "HELICALBINDINGSTRUCT"))
            {
              LINALG::Matrix<3,3> bspottriad;
              // auxiliary variable for storing a triad in quaternion form
              LINALG::Matrix<4, 1> qnode;
              // triad of node on first filament which is affected by the new crosslinker
              for (int j=0; j<4; j++)
                qnode(j) = (*bspottriadscol)[j][bspotLID];
              LARGEROTATIONS::quaterniontotriad(qnode, bspottriad);

              for(int j=0; j<(int)tangent.M(); j++)
              {
                tangent(j) = bspottriad(j,0);
                normal(j) = bspottriad(j,1);
              }
              // rotation angle
              alpha = (*bspotorientations_)[bspotLID];
            }
            else  // conventional case
            {
              // choose a second (neighbour) node
              const DRT::Node *node1 = NULL;
              int currfilament = (int)(*filamentnumber_)[bspotLID];
              if(bspotLID < basisnodes_-1)
              {
                if((*filamentnumber_)[bspotLID+1]==currfilament)
                  node1 = discret_->lColNode(bspotLID+1);
                else
                  node1 = discret_->lColNode(bspotLID-1);
              }
              if(bspotLID == basisnodes_-1)
                if((*filamentnumber_)[bspotLID-1]==currfilament)
                  node1 = discret_->lColNode(bspotLID-1);

              //calculate unit tangent
              for (int j=0; j<3; j++)
              {
                int dofgid1 = discret_->Dof(node1)[j];
                double nodeposj1 = node1->X()[j] + discol[discret_->DofColMap()->LID(dofgid1)];
                tangent(j) = nodeposj1 - nodepos0(j,0);
              }
              tangent.Scale(1 / tangent.Norm2());
              // calculate normal via cross product: [0 0 1]x[tx ty tz]
              normal.Clear();
              normal(0) = -tangent(1);
              normal(1) = tangent(0);
              // norm it since the cross product does not keep the length
              normal.Scale(1 / normal.Norm2());
              // random angle
              // by modulo operation involving the crosslink molecule number
              alpha = fmod((double) i, 2*M_PI);
            }

            // rotate the normal by alpha and store the new direction into the same Matrix ("normal")
            RotationAroundFixedAxis(tangent,normal,alpha);

            // calculation of the visualized point lying in the direction of the rotated normal
            for (int j=0; j<visualizepositions_->NumVectors(); j++)
              (*visualizepositions_)[j][i] = nodepos0(j,0) + ronebond*normal(j,0);
          }
        }
        break;
      }
    }
    if (periodlength_->at(0) > 0.0)
      CrosslinkerPeriodicBoundaryShift(visualizepositions_);
  }
  // synchronize results
  Teuchos::RCP<Epetra_MultiVector> visualizepositionstrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_, 3, true));
  CommunicateMultiVector(visualizepositionstrans, visualizepositions_);
  // debug couts
  //cout<<*visualizepositions_<<endl;
}//GmshPrepareVisualization

/*----------------------------------------------------------------------*
 | wedge output for two-noded beams                        cyron   11/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GmshWedge(const int& n,
                                  const Epetra_SerialDenseMatrix& coord,
                                  DRT::Element* thisele,
                                  std::stringstream& gmshfilecontent,
                                  const double color,
                                  bool ignoreeleid,
                                  bool drawsphere)
{
  //if this element is a line element capable of providing its radius get that radius
  double radius = 0.0;
  if(!ignoreeleid)
  {
    const DRT::ElementType & eot = thisele->ElementType();
    if(eot == DRT::ELEMENTS::Beam3Type::Instance())
      radius = sqrt(sqrt(4 * ((dynamic_cast<DRT::ELEMENTS::Beam3*>(thisele))->Izz()) / M_PI));
    else if(eot == DRT::ELEMENTS::Beam3iiType::Instance())
      radius = sqrt(sqrt(4 * ((dynamic_cast<DRT::ELEMENTS::Beam3ii*>(thisele))->Izz()) / M_PI));
    else if(eot == DRT::ELEMENTS::BeamCLType::Instance())
      radius = sqrt(sqrt(4 * ((dynamic_cast<DRT::ELEMENTS::BeamCL*>(thisele))->Izz()) / M_PI));
    else if(eot == DRT::ELEMENTS::Truss3Type::Instance())
      radius = sqrt((dynamic_cast<DRT::ELEMENTS::Truss3*>(thisele))->CSec() / M_PI);
    else
      dserror("thisele is not a line element providing its radius. Check your input file and your defines flags!");
    // case: crosslinker
    if(thisele->Id()>basisnodes_)
      radius = sqrt(statmechparams_.get<double>("ALINK",4.75166e-06) / M_PI); //default value according to diss. Tharmann
  }
  //line elements are plotted by a factor PlotFactorThick thicker than they are actually to allow for better visibility in gmsh pictures
  radius *= statmechparams_.get<double>("PlotFactorThick", 1.0);

  //solid lines plotted if PlotFactorThick != 0; otherwise output utilizes gmsh line visualization (lines width with fixed number of pixels, not a physical diameter)
  if(radius > 0.0)
  {

    // some local variables
    LINALG::Matrix<3,6> prism;
    LINALG::Matrix<3,1> axis;
    LINALG::Matrix<3,1> radiusvec1;
    LINALG::Matrix<3,1> radiusvec2;
    LINALG::Matrix<3,1> auxvec;
    LINALG::Matrix<3,1> theta;
    LINALG::Matrix<3,3> R;


    // compute three dimensional angle theta
    for (int j=0;j<3;++j)
      axis(j) = coord(j,1) - coord(j,0);
    double norm_axis = axis.Norm2();
    for (int j=0;j<3;++j)
      theta(j) = axis(j) / norm_axis * 2 * M_PI / n;

    // Compute rotation matirx R from rotation angle theta
    LARGEROTATIONS::angletotriad(theta,R);

    // Now the first prism will be computed via two radiusvectors, that point from each of
    // the nodes to two points on the beam surface. Further prisms will be computed via a
    // for-loop, where the second node of the previous prism is used as the first node of the
    // next prism, whereas the central points (=nodes) stay  identic for each prism. The
    // second node will be computed by a rotation matrix and a radiusvector.

    // compute radius vector for first surface node of first prims
    for (int j=0;j<3;++j) auxvec(j) = coord(j,0) + norm_axis;

    // radiusvector for point on surface
    radiusvec1(0) = auxvec(1)*axis(2) - auxvec(2)*axis(1);
    radiusvec1(1) = auxvec(2)*axis(0) - auxvec(0)*axis(2);
    radiusvec1(2) = auxvec(0)*axis(1) - auxvec(1)*axis(0);

    // initialize all prism points to nodes
    for (int j=0;j<3;++j)
    {
      prism(j,0) = coord(j,0);
      prism(j,1) = coord(j,0);
      prism(j,2) = coord(j,0);
      prism(j,3) = coord(j,1);
      prism(j,4) = coord(j,1);
      prism(j,5) = coord(j,1);
    }

    // get first point on surface for node1 and node2
    for (int j=0;j<3;++j)
    {
      prism(j,1) += radiusvec1(j) / radiusvec1.Norm2() * radius;
      prism(j,4) += radiusvec1(j) / radiusvec1.Norm2() * radius;
    }

    // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
    radiusvec2.Multiply(R,radiusvec1);

    // get second point on surface for node1 and node2
    for(int j=0;j<3;j++)
    {
      prism(j,2) += radiusvec2(j) / radiusvec2.Norm2() * radius;
      prism(j,5) += radiusvec2(j) / radiusvec2.Norm2() * radius;
    }

    // now first prism is built -> put coordinates into filecontent-stream
    // Syntax for gmsh input file
    // SI(x,y--,z,  x+.5,y,z,    x,y+.5,z,   x,y,z+.5, x+.5,y,z+.5, x,y+.5,z+.5){1,2,3,3,2,1};
    // SI( coordinates of the six corners ){colors}
    gmshfilecontent << "SI("<<scientific;
    gmshfilecontent << prism(0,0) << "," << prism(1,0) << "," << prism(2,0) << ",";
    gmshfilecontent << prism(0,1) << "," << prism(1,1) << "," << prism(2,1) << ",";
    gmshfilecontent << prism(0,2) << "," << prism(1,2) << "," << prism(2,2) << ",";
    gmshfilecontent << prism(0,3) << "," << prism(1,3) << "," << prism(2,3) << ",";
    gmshfilecontent << prism(0,4) << "," << prism(1,4) << "," << prism(2,4) << ",";
    gmshfilecontent << prism(0,5) << "," << prism(1,5) << "," << prism(2,5);
    gmshfilecontent << "){" << std::scientific;
    gmshfilecontent << color << "," << color << "," << color << "," << color << "," << color << "," << color << "};" << endl;

    // now the other prisms will be computed
    for (int sector=0;sector<n-1;++sector)
    {
      // initialize for next prism
      // some nodes of last prims can be taken also for the new next prism
      for (int j=0;j<3;++j)
      {
        prism(j,1)=prism(j,2);
        prism(j,4)=prism(j,5);
        prism(j,2)=prism(j,0);
        prism(j,5)=prism(j,3);
      }

      // old radiusvec2 is now radiusvec1; radiusvec2 is set to zero
      for (int j=0;j<3;++j)
      {
        radiusvec1(j) = radiusvec2(j);
        radiusvec2(j) = 0.0;
      }

      // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
      radiusvec2.Multiply(R,radiusvec1);


      // get second point on surface for node1 and node2
      for (int j=0;j<3;++j)
      {
        prism(j,2) += radiusvec2(j) / radiusvec2.Norm2() * radius;
        prism(j,5) += radiusvec2(j) / radiusvec2.Norm2() * radius;
      }

      // put coordinates into filecontent-stream
      // Syntax for gmsh input file
      // SI(x,y--,z,  x+.5,y,z,    x,y+.5,z,   x,y,z+.5, x+.5,y,z+.5, x,y+.5,z+.5){1,2,3,3,2,1};
      // SI( coordinates of the six corners ){colors}
      gmshfilecontent << "SI("<<scientific;
      gmshfilecontent << prism(0,0) << "," << prism(1,0) << "," << prism(2,0) << ",";
      gmshfilecontent << prism(0,1) << "," << prism(1,1) << "," << prism(2,1) << ",";
      gmshfilecontent << prism(0,2) << "," << prism(1,2) << "," << prism(2,2) << ",";
      gmshfilecontent << prism(0,3) << "," << prism(1,3) << "," << prism(2,3) << ",";
      gmshfilecontent << prism(0,4) << "," << prism(1,4) << "," << prism(2,4) << ",";
      gmshfilecontent << prism(0,5) << "," << prism(1,5) << "," << prism(2,5);
      gmshfilecontent << "){" << std::scientific;
      gmshfilecontent << color << "," << color << "," << color << "," << color << "," << color << "," << color << "};" << endl;
    }
  }
  //no thickness >0 specified; plot elements as line segements without physical volume
  else
  {
    gmshfilecontent << "SL(" << std::scientific;
    gmshfilecontent << coord(0,0) << "," << coord(1,0) << "," << coord(2,0) << ","
                    << coord(0,1) << "," << coord(1,1) << "," << coord(2,1);
    gmshfilecontent << ")" << "{" << std::scientific << color << ","<< color << "};" << endl;
  }
  // crosslink molecules are marked with an additional small ball if they are plotted as volumeless lines
  if(ignoreeleid && drawsphere)
  {
    double beadcolor = color;
    gmshfilecontent << "SP(" << std::scientific;
    gmshfilecontent << coord(0,1) << "," << coord(1,1) << ","<< coord(2,1);
    gmshfilecontent << ")" << "{" << std::scientific << beadcolor << ","<< beadcolor << "};" << endl;
  }
  return;
}//GmshWedge

/*----------------------------------------------------------------------*
 | Gmsh Output of detected network structure volume        mueller 12/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GmshNetworkStructVolume(const int& n, std::stringstream& gmshfilecontent, const double color)
{
  if(DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_structanaly)
  {
    cout<<"Visualizing test volume: ";
    switch(structuretype_)
    {
      // either cluster or homogeneous network
      case 0:
      {
        cout<<"Cluster"<<endl;
        // draw three octagons/hexadecagon lying in the base planes
        for(int i=0; i<3; i++)//spatial comp i
          for(int j=0; j<3; j++)//spatial comp j
            if(j<i)
              for(int k=0; k<3; k++)//spatial comp k perp. to the others
                if(k!=i && k!=j)
                {
                  double radius = characlength_/2.0;
                  // some local variables
                  LINALG::SerialDenseMatrix edge(3,2,true);
                  LINALG::Matrix<3,1> radiusvec1;
                  LINALG::Matrix<3,1> radiusvec2;
                  LINALG::Matrix<3,1> auxvec;
                  LINALG::Matrix<3,1> theta;
                  LINALG::Matrix<3,3> R;

                  // compute three dimensional angle theta
                  // unit vector perpendicular to jk-plane
                  LINALG::Matrix<3,1> axis;
                  for (int l=0;l<3;++l)
                    if(l==k)
                      axis(l) = 1.0;
                    else
                      axis(l) = 0.0;
                  for (int l=0;l<3;++l)
                    theta(l) = axis(l) * 2 * M_PI / n;

                  // Compute rotation matrix R from rotation angle theta
                  LARGEROTATIONS::angletotriad(theta,R);

                  // compute radius vector for first surface node of first edges
                  auxvec.Clear();
                  for (int l=0;l<3;++l)
                    if(l==j)
                      auxvec(l) = 1.0;

                  // radiusvector for point on surface
                  radiusvec1(0) = auxvec(1)*axis(2) - auxvec(2)*axis(1);
                  radiusvec1(1) = auxvec(2)*axis(0) - auxvec(0)*axis(2);
                  radiusvec1(2) = auxvec(0)*axis(1) - auxvec(1)*axis(0);

                  // initialize all edge points to nodes
                  for (int l=0;l<3;++l)
                  {
                    edge(l,0) = cog_(l);
                    edge(l,1) = cog_(l);
                  }

                  // get first point on surface for node1 and node2
                  for (int l=0;l<3;++l)
                    edge(l,0) += radiusvec1(l) / radiusvec1.Norm2() * radius;


                  // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
                  radiusvec2.Multiply(R,radiusvec1);

                  // get second point on surface for node1 and node2
                  for(int l=0;l<3;l++)
                    edge(l,1) += radiusvec2(l) / radiusvec2.Norm2() * radius;

                  // shift the coordinates according to periodic boundary conditions
                  LINALG::SerialDenseMatrix edgeshift = edge;
                  for(int l=0; l<(int)edge.M(); l++)
                    for(int m=0; m<(int)edge.N(); m++)
                    {
                      if(edge(l,m)>periodlength_->at(l))
                        edgeshift(l,m) -= periodlength_->at(l);
                      else if(edge(l,m)<0.0)
                        edgeshift(l,m) += periodlength_->at(l);
                    }
                  // write only the edge of the triangle (i.e. the line connecting two corners of the octagon/hexadecagon)
                  GmshOutputPeriodicBoundary(edgeshift, color, gmshfilecontent, discret_->lRowElement(0)->Id(),true);

                  // now the other edges will be computed
                  for (int sector=0;sector<n-1;++sector)
                  {
                    // initialize for next edge
                    for (int l=0;l<3;++l)
                    {
                      edge(l,0)=edge(l,1);
                      edge(l,1)=cog_(l);
                    }

                    // old radiusvec2 is now radiusvec1; radiusvec2 is set to zero
                    radiusvec1 = radiusvec2;
                    radiusvec2.Clear();

                    // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
                    radiusvec2.Multiply(R,radiusvec1);

                    // get second point on surface for node1 and node2
                    for (int l=0;l<3;++l)
                      edge(l,1) += radiusvec2(l) / radiusvec2.Norm2() * radius;
                    edgeshift = edge;
                    for(int l=0; l<(int)edge.M(); l++)
                      for(int m=0; m<(int)edge.N(); m++)
                      {
                        if(edge(l,m)>periodlength_->at(l))
                          edgeshift(l,m) -= periodlength_->at(l);
                        else if(edge(l,m)<0.0)
                          edgeshift(l,m) += periodlength_->at(l);
                      }

                    GmshOutputPeriodicBoundary(edgeshift, color, gmshfilecontent, discret_->lRowElement(0)->Id(),true);
                  }
                }
      }
      break;
      // bundle network
      case 1:
      {
        cout<<"Bundle"<<endl;

        double radius = characlength_;
        LINALG::Matrix<3,2> coord;
        for(int i=0; i<(int)coord.M(); i++)
          for(int j=0; j<(int)coord.N(); j++)
          coord(i,j) = testvolumepos_[j](i);

        // some local variables
        LINALG::Matrix<3,4> edges;
        LINALG::Matrix<3,1> axis;
        LINALG::Matrix<3,1> radiusvec1;
        LINALG::Matrix<3,1> radiusvec2;
        LINALG::Matrix<3,1> auxvec;
        LINALG::Matrix<3,1> theta;
        LINALG::Matrix<3,3> R;

        // compute three dimensional angle theta
        for (int j=0;j<3;++j)
          axis(j) = coord(j,1) - coord(j,0);
        double norm_axis = axis.Norm2();
        for (int j=0;j<3;++j)
          theta(j) = axis(j) / norm_axis * 2 * M_PI / n;

        // Compute rotation matrix R from rotation angle theta
        LARGEROTATIONS::angletotriad(theta,R);

        // compute radius vector for first surface node of first edges
        for (int j=0;j<3;++j) auxvec(j) = coord(j,0) + norm_axis;

        // radiusvector for point on surface
        radiusvec1(0) = auxvec(1)*axis(2) - auxvec(2)*axis(1);
        radiusvec1(1) = auxvec(2)*axis(0) - auxvec(0)*axis(2);
        radiusvec1(2) = auxvec(0)*axis(1) - auxvec(1)*axis(0);

        // initialize all edge points to nodes
        for (int j=0;j<3;++j)
        {
          edges(j,0) = coord(j,0);
          edges(j,1) = coord(j,1);
          edges(j,2) = coord(j,1);
          edges(j,3) = coord(j,0);
        }

        // get first point on surface for node1 and node2
        for (int j=0;j<3;++j)
        {
          edges(j,0) += radiusvec1(j) / radiusvec1.Norm2() * radius;
          edges(j,1) += radiusvec1(j) / radiusvec1.Norm2() * radius;
        }

        // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
        radiusvec2.Multiply(R,radiusvec1);

        // get second point on surface for node1 and node2
        for(int j=0;j<3;j++)
        {
          edges(j,2) += radiusvec2(j) / radiusvec2.Norm2() * radius;
          edges(j,3) += radiusvec2(j) / radiusvec2.Norm2() * radius;
        }
        // write only the edge of the triangle (i.e. the line connecting two corners of the octagon/hexadecagon)
        // and the connecting lines between the prism bases -> rectangle
        for(int j=1; j<(int)edges.N()+1; j++)
        {
          int jlow = (j-1)%(int)edges.N();
          int jhigh= j%(int)edges.N();
          int numsections = 10;
          // original (unshifted) edges
          LINALG::SerialDenseMatrix curredge(3,2);
          for(int k=0; k<curredge.M(); k++)
          {
            curredge(k,0) = edges(k, jlow);
            curredge(k,1) = edges(k, jhigh);
          }

          GmshNetworkStructVolumePeriodic(curredge, numsections, gmshfilecontent,color-0.5);
          /*gmshfilecontent << "SL(" << std::scientific;
          gmshfilecontent << edges(0,ilow) << "," << edges(1,ilow) << "," << edges(2,ilow) << ",";
          gmshfilecontent << edges(0,ihigh) << "," << edges(1,ihigh) << "," << edges(2,ihigh);
          gmshfilecontent << ")" << "{" << std::scientific << color-0.125 << ","<< color-0.125 << "};" << endl;*/
        }

        // now the other edges will be computed
        for (int sector=0;sector<n-1;++sector)
        {
          // initialize for next edge
          for (int j=0;j<3;++j)
          {
            edges(j,0)=edges(j,3);
            edges(j,1)=edges(j,2);
            edges(j,2)=coord(j,1);
            edges(j,3)=coord(j,0);
          }

          // old radiusvec2 is now radiusvec1; radiusvec2 is set to zero
          radiusvec1 = radiusvec2;
          radiusvec2.Clear();

          // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
          radiusvec2.Multiply(R,radiusvec1);

          // get second point on surface for node1 and node2
          for (int j=0;j<3;++j)
          {
            edges(j,2) += radiusvec2(j) / radiusvec2.Norm2() * radius;
            edges(j,3) += radiusvec2(j) / radiusvec2.Norm2() * radius;
          }

          // put coordinates into filecontent-stream
          for(int j=1; j<(int)edges.N()+1; j++)
          {
            int jlow = (j-1)%(int)edges.N();
            int jhigh= j%(int)edges.N();
            int numsections = 10;
            LINALG::SerialDenseMatrix curredge(3,2);
            for(int k=0; k<curredge.M(); k++)
            {
              curredge(k,0) = edges(k,jlow);
              curredge(k,1) = edges(k,jhigh);
            }

            GmshNetworkStructVolumePeriodic(curredge, numsections, gmshfilecontent,color-0.5);
            /*gmshfilecontent << "SL(" << std::scientific;
            gmshfilecontent << edges(0,ilow) << "," << edges(1,ilow) << "," << edges(2,ilow) << ",";
            gmshfilecontent << edges(0,ihigh) << "," << edges(1,ihigh) << "," << edges(2,ihigh);
            gmshfilecontent << ")" << "{" << std::scientific << color-0.125 << ","<< color-0.125 << "};" << endl;*/
          }
        }
      }
      break;
      // layer
      case 2:
      {
        cout<<"Layer ("<<testvolumepos_.size()<<"-noded)"<<endl;
        double halfthick = characlength_/2.0;
        // compute normal
        // cross product v_1 x v_2, plane normal
        LINALG::Matrix<3,1> normal;
        LINALG::Matrix<3,1> firstdir = testvolumepos_[1];
        LINALG::Matrix<3,1> secdir = testvolumepos_[1];
        firstdir -= testvolumepos_[0];
        secdir -= testvolumepos_[2];
        normal(0) = firstdir(1)*secdir(2) - firstdir(2)*secdir(1);
        normal(1) = firstdir(2)*secdir(0) - firstdir(0)*secdir(2);
        normal(2) = firstdir(0)*secdir(1) - firstdir(1)*secdir(0);
        normal.Scale(1.0/normal.Norm2());
        // upper and lower bound
        // componentwise comparison between two points to determine whether or not their connection is an edge
        // of the volume to visualize
        for(int i=0; i<(int)testvolumepos_.size(); i++)
          for(int j=0; j<(int)testvolumepos_.size(); j++)
            if(j>i) // consider entries above diagonal only
              for(int k=0; k<3; k++)
                if(fabs((testvolumepos_[i])(k)-(testvolumepos_[j])(k))<1e-7)
                {
                  int numsections = 10;

                  LINALG::SerialDenseMatrix loweredge(3,2);
                  LINALG::SerialDenseMatrix upperedge(3,2);
                  LINALG::SerialDenseMatrix connection0(3,2);
                  LINALG::SerialDenseMatrix connection1(3,2);
                  for(int l=0; l<loweredge.M(); l++)
                  {
                    loweredge(l,0) = testvolumepos_[i](l)-halfthick*normal(l);
                    loweredge(l,1) = testvolumepos_[j](l)-halfthick*normal(l);
                    upperedge(l,0) = testvolumepos_[i](l)+halfthick*normal(l);
                    upperedge(l,1) = testvolumepos_[j](l)+halfthick*normal(l);
                    connection0(l,0) = loweredge(l,0);
                    connection0(l,1) = upperedge(l,0);
                    connection1(l,0) = loweredge(l,1);
                    connection1(l,1) = upperedge(l,1);
                  }
                  // (long) edges
                  GmshNetworkStructVolumePeriodic(loweredge,numsections,gmshfilecontent,color-0.5);
                  GmshNetworkStructVolumePeriodic(upperedge,numsections,gmshfilecontent,color-0.5);
                  // connecting edges
                  GmshNetworkStructVolumePeriodic(connection0,4,gmshfilecontent,color-0.5);
                  GmshNetworkStructVolumePeriodic(connection1,4,gmshfilecontent,color-0.5);
                  // draw continous box also (for now)
                  // upper edge
                  gmshfilecontent << "SL(" << std::scientific;
                  gmshfilecontent << testvolumepos_[i](0)+halfthick*normal(0) << "," << testvolumepos_[i](1)+halfthick*normal(1) << "," << testvolumepos_[i](2)+halfthick*normal(2) << ",";
                  gmshfilecontent << testvolumepos_[j](0)+halfthick*normal(0) << "," << testvolumepos_[j](1)+halfthick*normal(1) << "," << testvolumepos_[j](2)+halfthick*normal(2);
                  gmshfilecontent << ")" << "{" << std::scientific << color-0.25 << ","<< color-0.25 << "};" << endl;
                  // lower edge
                  gmshfilecontent << "SL(" << std::scientific;
                  gmshfilecontent << testvolumepos_[i](0)-halfthick*normal(0) << "," << testvolumepos_[i](1)-halfthick*normal(1) << "," << testvolumepos_[i](2)-halfthick*normal(2) << ",";
                  gmshfilecontent << testvolumepos_[j](0)-halfthick*normal(0) << "," << testvolumepos_[j](1)-halfthick*normal(1) << "," << testvolumepos_[j](2)-halfthick*normal(2);
                  gmshfilecontent << ")" << "{" << std::scientific << color-0.25 << ","<< color-0.25 << "};" << endl;
                  // connections
                  gmshfilecontent << "SL(" << std::scientific;
                  gmshfilecontent << testvolumepos_[i](0)+halfthick*normal(0) << "," << testvolumepos_[i](1)+halfthick*normal(1) << "," << testvolumepos_[i](2)+halfthick*normal(2) << ",";
                  gmshfilecontent << testvolumepos_[i](0)-halfthick*normal(0) << "," << testvolumepos_[i](1)-halfthick*normal(1) << "," << testvolumepos_[i](2)-halfthick*normal(2);
                  gmshfilecontent << ")" << "{" << std::scientific << color-0.25 << ","<< color-0.25 << "};" << endl;
                  gmshfilecontent << "SL(" << std::scientific;
                  gmshfilecontent << testvolumepos_[j](0)+halfthick*normal(0) << "," << testvolumepos_[j](1)+halfthick*normal(1) << "," << testvolumepos_[j](2)+halfthick*normal(2) << ",";
                  gmshfilecontent << testvolumepos_[j](0)-halfthick*normal(0) << "," << testvolumepos_[j](1)-halfthick*normal(1) << "," << testvolumepos_[j](2)-halfthick*normal(2);
                  gmshfilecontent << ")" << "{" << std::scientific << color-0.25 << ","<< color-0.25 << "};" << endl;
                }
      }
      break;
      case 3:
        cout<<"Homogeneous network"<<endl;
      break;
    }
  }
}//GmshNetworkStructVolume()

/*----------------------------------------------------------------------*
 | Gmsh periodic bound. output for test volumes  (protected)mueller 1/11|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GmshNetworkStructVolumePeriodic(const Epetra_SerialDenseMatrix& coord, const int numsections, std::stringstream& gmshfilecontent,const double color)
{
  // direction of edge

  LINALG::Matrix<3, 1> dir;
  for(int k=0; k<(int)dir.M(); k++)
    dir(k) = coord(k,1)-coord(k,0);
  dir.Scale(1.0/dir.Norm2());
  // divide given edge into intervals for periodic boundary shift and visualization

  for(int k=0; k<numsections; k++)
  {
    // k-th intersection point
    LINALG::SerialDenseMatrix linepart(3,2);
    for(int l=0; l<(int)linepart.M(); l++)
    {
      linepart(l,0) = coord(l,0)+(coord(l,1)-coord(l,0))/(double)numsections*(double)k;
      linepart(l,1) = coord(l,0)+(coord(l,1)-coord(l,0))/(double)numsections*(double)(k+1);
    }

    // determine how the line is broken
    // Periodic boundary shift
    for(int l=0; l<linepart.M(); l++)
      for(int m=0; m<linepart.N(); m++)
      {
        if(linepart(l,m)>periodlength_->at(l))
          while(linepart(l,m)>periodlength_->at(l))
            linepart(l,m) -= periodlength_->at(l);
        if(linepart(l,m)<0.0)
          while(linepart(l,m)<0.0)
            linepart(l,m) += periodlength_->at(l);
      }

    LINALG::Matrix<3,1> cut;
    cut.Clear();
    for (int dof=0; dof<linepart.M(); dof++)
    {
      // point 0 close to zero boundary
      if (fabs(linepart(dof,1)-periodlength_->at(dof)-linepart(dof,0)) < fabs(linepart(dof,1) - linepart(dof,0)))
        cut(dof) = 1.0;
      // point 1 close to zero boundary
      if (fabs(linepart(dof,1)+periodlength_->at(dof) - linepart(dof,0)) < fabs(linepart(dof,1)-linepart(dof,0)))
        cut(dof) = 2.0;
    }

    // write special output for broken elements
    if (cut(0) + cut(1) + cut(2) > 0.0)
    {
      //from node 0 to nearest boundary where element is broken you get by vector X + lambda0*dir
      double lambda0 = 1e4;
      for (int dof = 0; dof < linepart.M(); dof++)
      {
        if (cut(dof) == 1.0)
        {
          if (fabs(-linepart(dof, 0) / dir(dof)) < fabs(lambda0))
            lambda0 = -linepart(dof, 0) / dir(dof);
        }
        else if (cut(dof) == 2.0)
        {
          if (fabs((periodlength_->at(dof) - linepart(dof, 0)) / dir(dof)) < fabs(lambda0))
            lambda0 = (periodlength_->at(dof) - linepart(dof, 0)) / dir(dof);
        }
      }
      //from node 1 to nearest boundary where element is broken you get by vector X + lambda1*dir
      double lambda1 = 1e4;
      for (int dof = 0; dof < linepart.M(); dof++)
      {
        if (cut(dof) == 2.0)
        {
          if (fabs(-linepart(dof, 1) / dir(dof)) < fabs(lambda1))
            lambda1 = -linepart(dof, 1) / dir(dof);
        }
        else if (cut(dof) == 1.0)
        {
          if (fabs((periodlength_->at(dof) - linepart(dof, 1)) / dir(dof)) < fabs(lambda1))
            lambda1 = (periodlength_->at(dof) - linepart(dof, 1)) / dir(dof);
        }
      }
      //define output coordinates for broken elements, first segment
      LINALG::SerialDenseMatrix coordout=linepart;
      for(int l=0 ;l<coordout.M(); l++)
        coordout(l,1) = linepart(l,0) + lambda0*dir(l);
      GmshWedge(1,coordout,discret_->lRowElement(0),gmshfilecontent,color,true);

      //define output coordinates for broken elements, second segment
      for(int l=0; l<coordout.M(); l++)
      {
        coordout(l,0) = linepart(l,1);
        coordout(l,1) = linepart(l,1)+lambda1*dir(l);
      }
      GmshWedge(1,coordout,discret_->lRowElement(0),gmshfilecontent,color,true);
    }
    else // output for continuous elements
      GmshWedge(1,linepart,discret_->lRowElement(0),gmshfilecontent,color,true,false);
  }
}//GmshNetworkStructVolumePeriodic()

/*----------------------------------------------------------------------*
 | output for density-density-correlation-function(public) mueller 07/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::DDCorrOutput(const Epetra_Vector&      disrow,
                                             const std::ostringstream& filename,
                                             const int&                istep,
                                             const double&             dt)
{
  /*Output:
   * (1) structure number and characteristic length (radius, thickness)
   * (2) internal energy
   * (3) histogram of inter-crosslinker distance of crosslinker elements -> Density-Density-Correlation
   *     in three columns (length of one column: numbins) due to periodic continuation of the boundary volume
   * (4) histograms of spherical coordinates (azimuth angle phi, polar angle theta/ cos(theta)
   * (5) radial density distribution
   */
  if(!discret_->Comm().MyPID())
    cout<<"\n\n====================== Analysis of structural polymorphism ======================"<<endl;

  if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2))
    dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

  int numbins = statmechparams_.get<int>("HISTOGRAMBINS", 1);
  // storage vector for shifted crosslinker LIDs(crosslinkermap)
  LINALG::Matrix<3,1> boxcenter;
  LINALG::Matrix<3,1> centershift;
  std::vector<int> crosslinkerentries;

#ifdef MEASURETIME
  double t0 = Teuchos::Time::wallTime();
#endif

  // determine the new center point of the periodic volume
  DDCorrShift(&boxcenter, &centershift, &crosslinkerentries);

#ifdef MEASURETIME
  double t1 = Teuchos::Time::wallTime();
#endif

  if(!discret_->Comm().MyPID())
    cout<<crosslinkerentries.size()<<" crosslinker elements found...\n"<<endl;

  // Determine current network structure
  LINALG::Matrix<3,1> cog;
  DDCorrCurrentStructure(disrow, &cog, &centershift, &crosslinkerentries, istep, filename);
  // store center of gravity for later use in GmshOutput()
  cog_ = cog;

#ifdef MEASURETIME
  double t2 = Teuchos::Time::wallTime();
#endif

  // Compute internal energy
  std::vector<double> internalenergy;
  internalenergy.clear();
  const RCP<Epetra_Vector> disp = Teuchos::rcp(new Epetra_Vector(disrow));
  ComputeInternalEnergy(disp, internalenergy,dt, filename);

#ifdef MEASURETIME
  double t3 = Teuchos::Time::wallTime();
#endif

  //calculcate the distance of crosslinker elements to all crosslinker elements (crosslinkermap)
  // testwise, base centershift upon calculated cog_ (rise in adequacy?)
  LINALG::Matrix<3,1> newcentershift;
  for(int i=0; i<(int)cog.M(); i++)
    newcentershift(i) = cog(i) - periodlength_->at(i)/2.0;

  // write the numbers of free, one-bonded, and two-bonded crosslink molecules
  CrosslinkCount(filename);

#ifdef MEASURETIME
  double t4 = Teuchos::Time::wallTime();
#endif

  // write the numbers of free, one-bonded, and two-bonded crosslink molecules
  OrientationCorrelation(disrow, istep);

#ifdef MEASURETIME
  double t5 = Teuchos::Time::wallTime();
#endif

  // Compute the the network mesh size in dependency to the radial distance to a given COG
  ComputeLocalMeshSize(disrow, newcentershift, istep);

#ifdef MEASURETIME
  double t6 = Teuchos::Time::wallTime();
#endif

  // MultiVector because result vector will be of length 3*ddcorrrowmap_->MyLength()
  Epetra_MultiVector crosslinksperbinrow(*ddcorrrowmap_,9 , true);
  Epetra_MultiVector crosslinksperbinrotrow(*ddcorrrowmap_,3 , true);
  DDCorrFunction(crosslinksperbinrow, crosslinksperbinrotrow, &newcentershift);

#ifdef MEASURETIME
  double t7 = Teuchos::Time::wallTime();
#endif

  // calculation of filament element orientation in spherical coordinates, sorted into histogram
  Epetra_Vector phibinsrow(*ddcorrrowmap_, true);
  Epetra_Vector thetabinsrow(*ddcorrrowmap_, true);
  Epetra_Vector costhetabinsrow(*ddcorrrowmap_, true);
  SphericalCoordsDistribution(disrow, phibinsrow, thetabinsrow, costhetabinsrow);

#ifdef MEASURETIME
  double t8 = Teuchos::Time::wallTime();
#endif

  Epetra_Vector radialdistancesrow(*ddcorrrowmap_, true);
  RadialDensityDistribution(radialdistancesrow, centershift);

#ifdef MEASURETIME
  double t9 = Teuchos::Time::wallTime();
#endif

  // Import
  Epetra_MultiVector crosslinksperbincol(*ddcorrcolmap_,crosslinksperbinrow.NumVectors() , true);
  Epetra_MultiVector crosslinksperbinrotcol(*ddcorrcolmap_,crosslinksperbinrotrow.NumVectors() , true);
  Epetra_Vector phibinscol(*ddcorrcolmap_, true);
  Epetra_Vector thetabinscol(*ddcorrcolmap_, true);
  Epetra_Vector costhetabinscol(*ddcorrcolmap_, true);
  Epetra_Vector radialdistancescol(*ddcorrcolmap_, true);
  Epetra_Import importer(*ddcorrcolmap_, *ddcorrrowmap_);
  crosslinksperbincol.Import(crosslinksperbinrow, importer, Insert);
  crosslinksperbinrotcol.Import(crosslinksperbinrotrow, importer, Insert);
  phibinscol.Import(phibinsrow, importer, Insert);
  thetabinscol.Import(thetabinsrow, importer, Insert);
  costhetabinscol.Import(costhetabinsrow, importer, Insert);
  radialdistancescol.Import(radialdistancesrow, importer, Insert);

  // Add the processor-specific data up
  std::vector<std::vector<int> > crosslinksperbin(numbins, std::vector<int>(crosslinksperbincol.NumVectors(),0));
  std::vector<std::vector<int> > crosslinksperbinrot(numbins, std::vector<int>(crosslinksperbinrotcol.NumVectors(),0));
  std::vector<int> phibins(numbins, 0);
  std::vector<int> thetabins(numbins, 0);
  std::vector<int> costhetabins(numbins, 0);
  std::vector<int> radialdistbins(numbins, 0);
  //int total = 0;
  for(int i=0; i<numbins; i++)
    for(int pid=0; pid<discret_->Comm().NumProc(); pid++)
    {
      for(int col=0; col<crosslinksperbincol.NumVectors(); col++)
        crosslinksperbin[i][col] += (int)crosslinksperbincol[col][pid*numbins+i];
      for(int col=0; col<crosslinksperbinrotcol.NumVectors(); col++)
        crosslinksperbinrot[i][col] += (int)crosslinksperbinrotcol[col][pid*numbins+i];
      phibins[i] += (int)phibinscol[pid*numbins+i];
      thetabins[i] += (int)thetabinscol[pid*numbins+i];
      costhetabins[i] += (int)costhetabinscol[pid*numbins+i];
      radialdistbins[i] += (int)radialdistancescol[pid*numbins+i];
    }

  // write data to file
  if(discret_->Comm().MyPID()==0)
  {
    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "a");
    std::stringstream histogram;
    // first part of output
    for(int i=0; i<numbins; i++)
    {
      histogram<<i+1<<"    ";
      for(int j=0; j<(int)crosslinksperbin[i].size(); j++)
        histogram<<crosslinksperbin[i][j]<<"    ";
      for(int j=0; j<(int)crosslinksperbinrot[i].size(); j++)
        histogram<<crosslinksperbinrot[i][j]<<"    ";
      histogram<<"  "<<phibins[i]<<"    "<<thetabins[i]<<"    "<<costhetabins[i]<<"    "<<radialdistbins[i]<<endl;
      //cout<<"  "<<phibins[i]<<"    "<<thetabins[i]<<"    "<<costhetabins[i]<<"    "<<radialdistbins[i]<<endl;
    }

    fprintf(fp, histogram.str().c_str());
    fclose(fp);
  }
  if(!discret_->Comm().MyPID())
  {
#ifdef MEASURETIME
    cout<<"\n=================Time  Measurement================"<<endl;
    cout<<"StatMechOutput::DDCorrOutput"<<endl;
    cout<<"DDCorrShift                 :\t"<<std::setprecision(4)<<t1-t0<<"\ts"<<endl;
    cout<<"DDCorrCurrentStructure      :\t"<<std::setprecision(4)<<t2-t1<<"\ts"<<endl;
    cout<<"ComputeInternalEnergy       :\t"<<std::setprecision(4)<<t3-t2<<"\ts"<<endl;
    cout<<"CrosslinkCount              :\t"<<std::setprecision(4)<<t4-t3<<"\ts"<<endl;
    cout<<"OrientationCorrelation      :\t"<<std::setprecision(4)<<t5-t4<<"\ts"<<endl;
    cout<<"ComputeLocalMeshSize        :\t"<<std::setprecision(4)<<t6-t5<<"\ts"<<endl;
    cout<<"DDCorrFunction              :\t"<<std::setprecision(4)<<t7-t6<<"\ts"<<endl;
    cout<<"SphericalCoordsDistribution :\t"<<std::setprecision(4)<<t8-t7<<"\ts"<<endl;
    cout<<"RadialDensityDistribution   :\t"<<std::setprecision(4)<<t9-t8<<"\ts"<<endl;
    cout<<"Communication               :\t"<<std::setprecision(4)<<Teuchos::Time::wallTime()-t9<<"\ts"<<endl;
    cout<<"=================================================="<<endl;
    cout<<"total time                  :\t"<<std::setprecision(4)<<Teuchos::Time::wallTime()-t0<<"\ts"<<endl;
#endif
    cout<<"================================================================================="<<endl;
  }
}//STATMECH::StatMechManager::DDCorrOutput()

/*------------------------------------------------------------------------------*
 | Selects raster point with the smallest average distance to all crosslinker   |
 | elements, makes it the new center of the boundary box and shifts crosslinker |
 | positions.                                                                   |
 |                                                        (public) mueller 11/10|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::DDCorrShift(LINALG::Matrix<3,1>* boxcenter, LINALG::Matrix<3,1>* centershift, std::vector<int>* crosslinkerentries)
{
  if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2))
    dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

  double periodlength = periodlength_->at(0);
  int numrasterpoints = statmechparams_.get<int>("NUMRASTERPOINTS", 3);
  // smallest average distance among all the average distances between the rasterpoints and all crosslinker elements
  // (init with 2*pl,so that it is definitely overwritten by the first "real" value)
  double  smallestdistance = 2.0*periodlength_->at(0);

  //store crosslinker element position within crosslinkerbond to crosslinkerentries
  for(int i=0; i<crosslinkerbond_->MyLength(); i++)
    if((*crosslinkerbond_)[0][i]>-0.9 && (*crosslinkerbond_)[1][i]>-0.9)
      crosslinkerentries->push_back(i);

  int numcrossele = (int)crosslinkerentries->size();

  // determine the new center of the box
  if(numcrossele>0)
  {
    for(int i=0; i<numrasterpoints; i++)
      for(int j=0; j<numrasterpoints; j++)
        for(int k=0; k<numrasterpoints; k++)
        {
          double averagedistance = 0.0;
          LINALG::Matrix<3,1> currentrasterpoint;
          LINALG::Matrix<3,1> currentcentershift;

          // calculate current raster point
          currentrasterpoint(0) = i*periodlength/(numrasterpoints-1);
          currentrasterpoint(1) = j*periodlength/(numrasterpoints-1);
          currentrasterpoint(2) = k*periodlength/(numrasterpoints-1);

          // calculate the center shift (difference vector between regular center and new center of the boundary box)
          for(int l=0; l<(int)currentrasterpoint.M(); l++)
            currentcentershift(l) = currentrasterpoint(l)-periodlength/2.0;

          // calculate average distance of crosslinker elements to raster point
          for(int l=0; l<numcrossele; l++)
          {
            // get the crosslinker position in question and shift it according to new boundary box center
            LINALG::Matrix<3,1> distance;
            for(int m=0; m<(int)distance.M(); m++)
            {
              distance(m) = (*crosslinkerpositions_)[m][(*crosslinkerentries)[l]];
              if (distance(m) > periodlength+currentcentershift(m))
                distance(m) -= periodlength;
              if (distance(m) < 0.0+currentcentershift(m))
                distance(m) += periodlength;
            }
            distance -= currentrasterpoint;
            averagedistance += distance.Norm2();
          }
          averagedistance /= (double)crosslinkerentries->size();

          if(averagedistance<smallestdistance)
          {
            smallestdistance = averagedistance;
            *boxcenter = currentrasterpoint;
            *centershift = currentcentershift;
          }
        }
  }
  else
  {
    boxcenter->PutScalar(periodlength/2.0);
    centershift->PutScalar(0.0);
  }
  //if(!discret_->Comm().MyPID())
    //cout<<"Box Center(2): "<<(*boxcenter)[0]<<", "<<(*boxcenter)[1]<<", "<<(*boxcenter)[2]<<endl;
  return;
}//STATMECH::StatMechManager::DDCorrShift()

/*------------------------------------------------------------------------------*
 | Determine current network structure and output network type as single        |
 | characteristic number. Also, output filament orientations.                   |
 |                                                        (public) mueller 11/10|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::DDCorrCurrentStructure(const Epetra_Vector& disrow,
                                             LINALG::Matrix<3,1>* cog,
                                             LINALG::Matrix<3,1>* centershift,
                                             std::vector<int>* crosslinkerentries,
                                             const int& istep,
                                             const std::ostringstream& filename,
                                             bool filorientoutput)
{
  if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2))
    dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

  double periodlength = periodlength_->at(0);
  // number of crosslinker elements
  int numcrossele = (int)crosslinkerentries->size();
  std::vector<int> crosslinksinvolume(3,0);

  // get column map displacements
  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(disrow, discol);

  /// calculate center of gravity with respect to shiftedpositions for bound crosslinkers
  std::vector<LINALG::Matrix<3,1> > shiftedpositions;
  cog->Clear();
  // shift positions according to new center point (raster point)
  for(int i=0; i<numcrossele; i++)
  {
    LINALG::Matrix<3,1> currposition;
    for(int j=0; j<(int)currposition.M(); j++)
    {
      currposition(j) = (*crosslinkerpositions_)[j][(*crosslinkerentries)[i]];
      if (currposition(j) > periodlength+(*centershift)(j))
        currposition(j) -= periodlength;
      if (currposition(j) < 0.0+(*centershift)(j))
        currposition(j) += periodlength;
    }
    (*cog) += currposition;
    if(discret_->Comm().MyPID()==0)
      shiftedpositions.push_back(currposition);
  }
  if(numcrossele != 0)
    cog->Scale(1.0/(double)numcrossele);

  // zero out for new run
  trafo_->Zero();

  // calculations done by Proc 0 only
  if(discret_->Comm().MyPID()==0)
  {
    // number indicating structure type
    int structurenumber = 0;
    // indices for layer plane vectors
    int dir1 = -1, dir2 = -1;
    double rlink = statmechparams_.get<double>("R_LINK", 1.0);
    // number of nested intervals
    int maxexponent = (int)ceil(log(periodlength/rlink)/log(2.0))*2;
    // vector for test volumes (V[0]-sphere, V[1]-cylinder, V[2]-layer/homogenous network)
    std::vector<double> volumes(3,pow(10*periodlength, 3.0));
    // characteristic lengths of "volumes"
    std::vector<double> characlength(3,10*periodlength);
    // crosslinker fraction included in test volume
    std::vector<double> crossfraction(3,0.0);
    // number of iterations until crosslinker fraction lies within the given threshold fraction +/- tolerance
    std::vector<int> niter(3,0);

    // iterated vectors
    // bundle
    LINALG::Matrix<3,1> cylvec;
    // layer
    std::vector<LINALG::Matrix<3,1> > layervectors;
    // cluster-layer vectors (layer vectors determined when actually, a cluster phase is detected)
    std::vector<LINALG::Matrix<3,1> > clusterlayervecs;

    // get the intersections of the axis of a cylinder with the (two) cube faces
    std::vector<LINALG::Matrix<3,1> > intersections;
    // coordinates of intersection points of layer-type volume with the cube edges
    std::vector<LINALG::Matrix<3,1> > interseccoords;

    if(numcrossele>0)
    {
  /// calculate normed vectors and output filament element orientations
      // normed vectors for structural analysis (projections of base vectors e1, e2, e3 onto structure)
      std::vector<LINALG::Matrix<3,1> > normedvectors;
      for(int i=0; i<3; i++)
      {
        LINALG::Matrix<3,1> normedi;
        normedi.Clear();
        normedvectors.push_back(normedi);
      }

  /// determine normed vectors as well as output filament element vectors
      std::ostringstream orientfilename;
      orientfilename << "./FilamentOrientations_"<<std::setw(6) << std::setfill('0') << istep <<".dat";
      FilamentOrientations(discol, &normedvectors, orientfilename, filorientoutput);

      // select the vector best fitting the axis of the bundle cylinder: vector with the greatest length -> smallest
      // enclosed angle with the cylinder axis
      int startiterindex = 0;
      double veclength = 0.0;
      cout<<"Normed vectors:"<<endl;
      for(int i=0; i<(int)normedvectors.size(); i++)
      {
        if(normedvectors[i].Norm2() > veclength)
        {
          veclength = normedvectors[i].Norm2();
          startiterindex = i;
        }
        // Scale normed vectors to unit length
        normedvectors[i].Scale(1.0/normedvectors[i].Norm2());
        cout<<i<<": "<<normedvectors[i](0)<<" "<<normedvectors[i](1)<<" "<<normedvectors[i](2)<<endl;
      }

  /// determine network structure
      // threshold fraction of crosslinkers
      double pthresh = 0.9;

      // calculate smallest possible test volumes that fulfill the requirement /numcrossele >= pthresh
      for(int i=0; i<(int)volumes.size(); i++)
      {
        switch(i)
        {
          // spherical volume
          case 0:
          {
            bool leaveloop = false;
            // tolerance
            double lowertol = 0.05;
            double uppertol = 0.01;
            // initial search radius
            double radius = periodlength/2.0;
            // fraction of crosslinks within test volume
            double pr = 0.0;
            int exponent = 1;

            // determine the two normed vectors with the largest cross product 2-Norm (closest to being perpendicular)
            // this is done in order to capture the smooth transformation from cluster to layer. More precisely, we
            // calculate the layer triad. The information we gain from this procedure is useful once we observe a
            // filament orientation distribution which is no longer homogeneous.
            double crossprodnorm = 0.0;
            for(int j=0; j<3; j++)
              for(int k=0; k<3; k++)
                if(k>j)
                {
                  // cross product
                  LINALG::Matrix<3,1> crossvec;
                  crossvec(0) = normedvectors[j](1)*normedvectors[k](2) - normedvectors[j](2)*normedvectors[k](1);
                  crossvec(1) = normedvectors[j](2)*normedvectors[k](0) - normedvectors[j](0)*normedvectors[k](2);
                  crossvec(2) = normedvectors[j](0)*normedvectors[k](1) - normedvectors[j](1)*normedvectors[k](0);
                  if(crossvec.Norm2()>crossprodnorm)
                  {
                    crossprodnorm = crossvec.Norm2();
                    dir1=j;
                    dir2=k;
                  }
                }

            // store initial vectors to be iterated
            clusterlayervecs.push_back(normedvectors[dir1]);
            clusterlayervecs.push_back(normedvectors[dir2]);

            // iterate both vectors in order to obtain projections into the layer plane
            const int maxiterations = 25;
            cout<<"\nVector iteration:"<<endl;
            cout<<"Cluster1: ";
            DDCorrIterateVector(discol, &clusterlayervecs[0], maxiterations);
            cout<<"Cluster2: ";
            DDCorrIterateVector(discol, &clusterlayervecs[1], maxiterations);

            // cross product n_1 x n_2, plane normal
            LINALG::Matrix<3,1> normal;
            normal(0) = clusterlayervecs[0](1)*clusterlayervecs[1](2) - clusterlayervecs[0](2)*clusterlayervecs[1](1);
            normal(1) = clusterlayervecs[0](2)*clusterlayervecs[1](0) - clusterlayervecs[0](0)*clusterlayervecs[1](2);
            normal(2) = clusterlayervecs[0](0)*clusterlayervecs[1](1) - clusterlayervecs[0](1)*clusterlayervecs[1](0);
            normal.Scale(1.0/normal.Norm2());
            clusterlayervecs.push_back(normal);

            // loop as long as pr has not yet reached pthresh
            // flag indicating convergence (set to false, if max no. of iterations is reached)
            bool converged = true;
            while(!leaveloop)
            {
              int rcount = 0;
              // loop over crosslinker elements
              for(int j=0; j<numcrossele; j++)
              {
                // get distance of crosslinker element to center of gravity
                LINALG::Matrix<3,1> dist = shiftedpositions[j];
                dist -= *cog;
                if(dist.Norm2()<=radius)
                  rcount++;
              }
              pr = double(rcount)/double(numcrossele);

              exponent++;
              niter[0]++;
              // new radius
              if(exponent<=maxexponent)
              {
                if((pr<pthresh-lowertol || pr>pthresh+uppertol))
                {
                  // determine "growth direction"
                  double sign;
                  if(pr<pthresh)
                    sign = 1.0;
                  else
                    sign = -1.0;
                  radius += sign*periodlength/pow(2.0,(double)exponent);
                }
                else
                {
                  crosslinksinvolume[i] = rcount;
                  crossfraction[i] = pr;
                  leaveloop = true;
                }
              }
              else
              {
                converged = false;
                crosslinksinvolume[i] = rcount;
                crossfraction[i] = pr;
                leaveloop = true;
              }
            }
            // store characteristic length and test sphere volume
            if(converged)
            {
              // store cluster diameter
              characlength[i] = 2.0*radius;
              volumes[i] = 4/3 * M_PI * pow(radius, 3.0);
            }
          }
          break;
          // cylindrical volume
          case 1:
          {
            bool leaveloop = false;
            // tolerance
            double lowertol = 0.02;
            double uppertol = 0.02;
            double radius = periodlength/2.0;
            double cyllength = 0.0;
            double pr = 0.0;
            int exponent = 1;

            // iterate to obtain fitting normed direction
            LINALG::Matrix<3,1> normj = normedvectors[startiterindex];
            const int maxiterations = 25;
            cout<<"Bundle:   ";
            DDCorrIterateVector(discol, &normj, maxiterations);
            // for output
            cylvec = normj;

            // cube face boundaries of jk-surface of cubical volume
            LINALG::Matrix<3,2> surfaceboundaries;
            for(int j=0; j<(int)surfaceboundaries.M(); j++)
            {
              surfaceboundaries(j,0) = (*centershift)(j);
              surfaceboundaries(j,1) = (*centershift)(j)+periodlength;
            }
            for(int j=0; j<(int)surfaceboundaries.N(); j++)
              for(int k=0; k<3; k++)
                for(int l=0; l<3; l++)
                  if(l>k)
                    for(int m=0; m<3; m++)
                      if(m!=k && m!=l)
                      {
                        LINALG::Matrix<3,1> currentintersection;
                        // known intersection component
                        currentintersection(m) = surfaceboundaries(m,j);
                        // get line parameter
                        double lambdaline = (currentintersection(m)-(*cog)(m))/normj(m);
                        currentintersection(k) = (*cog)(k)+lambdaline*normj(k);
                        currentintersection(l) = (*cog)(l)+lambdaline*normj(l);
                        // check if intersection lies on volume boundary
                        if(currentintersection(k)<=surfaceboundaries(k,1) && currentintersection(k)>=surfaceboundaries(k,0) &&
                           currentintersection(l)<=surfaceboundaries(l,1) && currentintersection(l)>=surfaceboundaries(l,0))
                          intersections.push_back(currentintersection);
                      }
            LINALG::Matrix<3,1> deltaisecs = intersections[1];
            deltaisecs -= intersections[0];
            cyllength = deltaisecs.Norm2();

            // calculate fraction of crosslinkers within cylinder
            // flag indicating convergence (set to false, if max no. of iterations is reached)
            bool converged = true;
            while(!leaveloop)
            {
              int rcount = 0;
              // loop over crosslinker elements
              for(int j=0; j<numcrossele; j++)
              {
                // get distance of crosslinker element to normed1 through center of gravity
                // intersection line-plane
                LINALG::Matrix<3,1> crosstocog = shiftedpositions[j];
                crosstocog -= *cog;

                double numerator = crosstocog.Dot(normj);
                double denominator = (normj).Dot(normj);
                double lambdaisec = numerator/denominator;
                // intersection and distance of crosslinker to intersection
                LINALG::Matrix<3,1> isecpt = normj;
                isecpt.Scale(lambdaisec);
                isecpt += (*cog);
                LINALG::Matrix<3,1> deltaiseccog = shiftedpositions[j];
                deltaiseccog -= isecpt;
                double distance = deltaiseccog.Norm2();

                if(distance<=radius)
                  rcount++;
              }
              pr = double(rcount)/double(numcrossele);

              exponent++;
              niter[1]++;
              // new radius
              if(exponent<=maxexponent)
              {
                if((pr<pthresh-lowertol || pr>pthresh+uppertol))
                {
                  // determine "growth direction"
                  double sign;
                  if(pr<pthresh)
                    sign = 1.0;
                  else
                    sign = -1.0;
                  radius += sign*periodlength/pow(2.0,(double)exponent);
                }
                else
                {
                  crosslinksinvolume[i] = rcount;
                  crossfraction[i] = pr;
                  leaveloop = true;
                }
              }
              else
              {
                converged = false;
                crosslinksinvolume[i] = rcount;
                crossfraction[i] = pr;
                leaveloop = true;
              }
            }
            if(converged)
            {
              // store bundle diameter
              characlength[i] = 2.0*radius;
              volumes[i] = M_PI*radius*radius*cyllength;
            }
          }
          break;
          // cuboid layer volume
          case 2:
          {
            bool leaveloop = false;
            // tolerance
            double lowertol = 0.02;
            double uppertol = 0.02;
            double thickness = periodlength/2.0;
            double pr = 0.0;
            int exponent = 1;

            //determine the two normed vectors with the largest cross product 2-Norm (closest to being perpendicular)
            double crossprodnorm = 0.0;
            for(int j=0; j<3; j++)
              for(int k=0; k<3; k++)
                if(k>j)
                {
                  // cross product
                  LINALG::Matrix<3,1> crossvec;
                  crossvec(0) = normedvectors[j](1)*normedvectors[k](2) - normedvectors[j](2)*normedvectors[k](1);
                  crossvec(1) = normedvectors[j](2)*normedvectors[k](0) - normedvectors[j](0)*normedvectors[k](2);
                  crossvec(2) = normedvectors[j](0)*normedvectors[k](1) - normedvectors[j](1)*normedvectors[k](0);
                  if(crossvec.Norm2()>crossprodnorm)
                  {
                    crossprodnorm = crossvec.Norm2();
                    dir1=j;
                    dir2=k;
                  }
                }

            // store initial vectors to be iterated
            layervectors.push_back(normedvectors[dir1]);
            layervectors.push_back(normedvectors[dir2]);

            // iterate both vectors in order to obtain projections into the layer plane
            const int maxiterations = 25;
            cout<<"Layer1:  ";
            DDCorrIterateVector(discol, &layervectors[0], maxiterations);
            cout<<"Layer2:  ";
            DDCorrIterateVector(discol, &layervectors[1], maxiterations);

            // cross product n_1 x n_2, plane normal
            LINALG::Matrix<3,1> normal;
            normal(0) = layervectors[0](1)*layervectors[1](2) - layervectors[0](2)*layervectors[1](1);
            normal(1) = layervectors[0](2)*layervectors[1](0) - layervectors[0](0)*layervectors[1](2);
            normal(2) = layervectors[0](0)*layervectors[1](1) - layervectors[0](1)*layervectors[1](0);
            normal.Scale(1.0/normal.Norm2());
            layervectors.push_back(normal);

            // crosslinkerpositions which are considered to be within the cuboid
            std::vector<LINALG::Matrix<3,1> > crosslinkswithinvolume;

            // flag indicating convergence (set to false, if max no. of iterations is reached)
            bool converged = true;
            while(!leaveloop)
            {
              for(int j=0; j<numcrossele; j++)
              {
                // given, that cog E plane with normal vector "normal"
                // constant in Hessian normal form
                double d = normal.Dot((*cog));
                // distance of crosslinker element to plane
                double pn = normal.Dot(shiftedpositions[j]);
                double disttoplane = fabs(pn-d);

                if(disttoplane <= thickness)
                  crosslinkswithinvolume.push_back(shiftedpositions[j]);
              }
              pr = double(crosslinkswithinvolume.size())/double(numcrossele);

              exponent++;
              niter[2]++;

              if(exponent<=maxexponent)
              {
                if((pr<pthresh-lowertol || pr>pthresh+uppertol))
                {
                  double sign;
                  if(pr<pthresh)
                    sign = 1.0;
                  else
                    sign = -1.0;
                  thickness += sign*periodlength/pow(2.0,(double)exponent);
                  crosslinkswithinvolume.clear();
                }
                else
                {
                  crosslinksinvolume[i] = (int)crosslinkswithinvolume.size();
                  crossfraction[i] = pr;
                  leaveloop = true;
                }
              }
              else
              {
                converged = false;
                crosslinksinvolume[i] = (int)crosslinkswithinvolume.size();
                crossfraction[i] = pr;
                leaveloop = true;
              }
            }
            // calculation of the volume
            // according to number of intersection points (only if crosslinker fraction calculation converged)
            if(converged)
            {
              // cube face boundaries of jk-surface of cubical volume
              LINALG::Matrix<3,2> surfaceboundaries;
              for(int j=0; j<(int)surfaceboundaries.M(); j++)
              {
                surfaceboundaries(j,0) = (*centershift)(j);
                surfaceboundaries(j,1) = (*centershift)(j)+periodlength;
              }

              // first step: find the intersection points layer x volume egdes.
              // At first, we assume a layer delimited by the volume boundaries)
              for(int surf=0; surf<(int)surfaceboundaries.N(); surf++) // (two) planes perpendicular to l-direction
              {
                for(int j=0; j<3; j++) // spatial component j
                  for(int k=0; k<3; k++) // spatial component k
                    if(k>j)// above diagonal
                      for(int l=0; l<3; l++) // spatial component l
                        if(l!=j && l!=k)
                        {
                          LINALG::Matrix<3,1> coords;
                          coords(l) = surfaceboundaries(l,surf);
                          for(int edge=0; edge<(int)surfaceboundaries.N(); edge++)
                            for(int m=0; m<(int)surfaceboundaries.M(); m++)
                              if(l!=m)
                              {
                                coords(m) = surfaceboundaries(m,edge);
                                for(int n=0; n<(int)surfaceboundaries.M(); n++)
                                  if(n!=m && n!=l)
                                  {
                                    coords(n) = (((*cog)(l)-coords(l))*normal(l) - (coords(m)-(*cog)(m))*normal(m))/normal(n) + (*cog)(n);
                                    double lowerbound = surfaceboundaries(n,0);
                                    double upperbound = surfaceboundaries(n,1);
                                    if((coords(n)>= lowerbound || fabs(coords(n)-lowerbound)<1e-6) && (coords(n)<= upperbound || fabs(coords(n)-upperbound)<1e-6))
                                    {
                                      bool redundant = false;
                                      // check for redundant entries
                                      if((int)interseccoords.size()>0)
                                      {
                                        LINALG::Matrix<3,1> check = coords;
                                        for(int p=0; p<(int)interseccoords.size(); p++)
                                        {
                                          check -= interseccoords[p];
                                          if(check.Norm2()<1e-4)
                                            redundant = true;
                                          check = coords;
                                        }
                                      }
                                      if(!redundant)
                                        interseccoords.push_back(coords);
                                    }
                                  }
                              }
                        }
              }

              //second step: Approximate the base area of the volume/prism
              // calculate and store unit direction vectors corners-COG
              std::vector<LINALG::Matrix<3,1> > itercoords = interseccoords;
              std::vector<LINALG::Matrix<3,1> > dirvecs;
              std::vector<double> initdistances;
              //cout<<"initdist = ";
              for(int j=0; j<(int)interseccoords.size(); j++)
              {
                // directional vectors pointing from center of gravity towards the corners
                LINALG::Matrix<3,1> dirvec = interseccoords[j];
                dirvec -= *cog;
                initdistances.push_back(dirvec.Norm2());
                //cout<<dirvec.Norm2()<<" ";
                dirvec.Scale(1.0/dirvec.Norm2());
                dirvecs.push_back(dirvec);
              }
              //cout<<endl;

              // calculate crosslinker projection positions in layer plane
              std::vector<LINALG::Matrix<3,1> > crossprojections;
              for(int j=0; j<(int)crosslinkswithinvolume.size(); j++)
              {
                // calculate projection of crosslinker onto layer plane
                LINALG::Matrix<3,1> projection = crosslinkswithinvolume[j];
                LINALG::Matrix<3,1> diffvec = normal;
                // line paramater
                double mu = (interseccoords[0].Dot(normal) - crosslinkswithinvolume[j].Dot(normal)) / (normal.Dot(normal));
                diffvec.Scale(mu);
                projection += diffvec;
                crossprojections.push_back(projection);
              }

              // calculate edge directions of the layer poygon
              std::vector<LINALG::Matrix<3,1> > edgedirs;
              std::vector<int> neworder;
              for(int j=0; j<(int)interseccoords.size(); j++)
                for(int k=0; k<(int)interseccoords.size(); k++)
                  if(k>j)
                    for(int l=0; l<(int)interseccoords[j].M(); l++)
                      if(fabs((interseccoords[j])(l)-(interseccoords[k])(l))<1e-7)
                      {
                        LINALG::Matrix<3,1> jdir = interseccoords[k];
                        jdir -= interseccoords[j];
                        jdir.Scale(1.0/jdir.Norm2());
                        edgedirs.push_back(jdir);
                        neworder.push_back(j);
                      }

              // iterate as long as layer volume contains 90-95% of the crosslinks
              leaveloop = false;
              int numiter = 30;
              int iter = 0;
              int rcount = 0;
              double uplim = 0.01;
              double lowlim = 0.03;
              double threshold = 0.95;
              // start with 1.0 (100% of crosslinkers within layer volume)
              double crossfrac = 1.0;
              while(!leaveloop)
              {
                if(iter<numiter)
                {
                  if(crossfrac<threshold-lowlim || crossfrac>threshold+uplim)
                  {
                    //calculate new set of layer corner coordinates
                    for(int j=0; j<(int)itercoords.size(); j++)
                    {
                      LINALG::Matrix<3,1> deltaj = dirvecs[j];
                      deltaj.Scale(initdistances[j]/pow(1.4,(double)(iter+1)));
                      // depending on crossfrac, choose in which direction the triangle grows:
                      if(crossfrac>threshold)
                        itercoords[j] -= deltaj;
                      if(crossfrac<threshold)
                        itercoords[j] += deltaj;
                    }

                    // check if crosslinkers lie within the volume with the new base area
                    // number of crosslinkers
                    rcount = 0;
                    for(int j=0; j<(int)crosslinkswithinvolume.size(); j++)
                    {
                      // distance of crosslinker projection to center of gravity
                      LINALG::Matrix<3,1> crosstocog = crossprojections[j];
                      crosstocog -= *cog;
                      double dcrosstocog = crosstocog.Norm2();
                      crosstocog.Scale(1.0/dcrosstocog);

                      // calculate shortest positive distance for intersections of the line (cog->crosslinker) with the layer polygon from the center of gravity
                      double dclosestisec = 10*periodlength;
                      for(int k=0; k<(int)edgedirs.size(); k++)
                      {
                        // line paramete (=length)
                        double numerator = edgedirs[k](1)*(itercoords[neworder[k]](0)-(*cog)(0)) - edgedirs[k](0)*(itercoords[neworder[k]](1)-(*cog)(1));
                        double denominator = crosstocog(0)*edgedirs[k](1) - crosstocog(1)*edgedirs[k](0);
                        // lambda either positive (i.e. intersection lies in the direction of the crosslinker) or negative
                        double lambda = numerator/denominator;
                        // in case of parallel directional vectors, lambda takes big value (does not matter anyway as this case is not important for the coming calculations)
                        if(fabs(denominator) < 1e-8)
                          lambda = 1e8;
                        // save smallest positive lambda. It marks the closest intersection of the (cog->crosslinker)-line with an edge
                        if(lambda>0.0 && lambda<dclosestisec)
                          dclosestisec = lambda;
                      }
                      // if the distance to the closest intersection is bigger than the distance between crosslinker and center of gravity,
                      // the linker in question lies within the iterated volume
                      if(dclosestisec>=dcrosstocog)
                        rcount++;
                    }
                    //new crosslinker fraction in volume
                    crossfrac = double(rcount)/double(crosslinkswithinvolume.size());
                    //cout<<"i="<<iter<<", crossfrac = "<<crossfrac<<", rcount = "<<rcount<<"/"<<crosslinkswithinvolume.size()<<endl;
                    iter++;
                  }
                  else
                  {
                    // set itercoords as new interseccoords
                    interseccoords = itercoords;
                    crosslinksinvolume[i] = rcount;
                    crossfraction[i] *= pr;
                    cout<<"-> adjusted volume after "<<iter<<" iterations with "<<rcount<<"/"<<crosslinkswithinvolume.size()<<"( p="<<crossfrac<<" )"<<endl;
                    leaveloop = true;
                  }
                }
                else
                {
                  // set itercoords as new interseccoords
                  interseccoords = itercoords;
                  crosslinksinvolume[i] = rcount;
                  crossfraction[i] *= pr;
                  cout<<"-> adjusted volume after maxiter = "<<iter<<" iterations with "<<rcount<<"/"<<crosslinkswithinvolume.size()<<"( p="<<crossfrac<<" )"<<endl;
                  leaveloop = true;
                }
              }

              // calculation of layer volume
              switch((int)interseccoords.size())
              {
                // triangle
                case 3:
                {
                  // area of the triangle via cosine rule
                  LINALG::Matrix<3,1> c = interseccoords[0];
                  c -= interseccoords[1];
                  LINALG::Matrix<3,1> a = interseccoords[1];
                  a -= interseccoords[2];
                  LINALG::Matrix<3,1> b = interseccoords[0];
                  b -= interseccoords[2];
                  double cl = c.Norm2();
                  double al = a.Norm2();
                  double bl = b.Norm2();
                  double alpha = acos((cl*cl+bl*bl-al*al)/(2.0*bl*cl));
                  double h = bl * sin(alpha);

                  // volume of layer (factor 0.5 missing, since "real" thickness is thickness*2.0)
                  volumes[i] = cl*h*thickness;
                }
                break;
                // square/rectangle/trapezoid
                case 4:
                {
                  // edges
                  LINALG::Matrix<3,1> a = interseccoords[1];
                  a -= interseccoords[0];
                  LINALG::Matrix<3,1> c = interseccoords[3];
                  c -= interseccoords[2];
                  LINALG::Matrix<3,1> d = interseccoords[2];
                  d -= interseccoords[0];
                  // diagonal
                  LINALG::Matrix<3,1> f = interseccoords[2];
                  f -= interseccoords[1];
                  double al = a.Norm2();
                  double cl = c.Norm2();
                  double dl = d.Norm2();
                  double fl = f.Norm2();
                  double alpha = acos((al*al+dl*dl-fl*fl)/(2.0*al*dl));
                  double h = dl * sin(alpha);

                  volumes[i] = (al+cl) * h * thickness;
                  //cout<<"layer volume_rec = "<<volumes[i]<<endl;
                }
                break;
                // hexahedron
                case 6:
                {
                  double hexvolume = 0.0;
                  // indices mapping correct order of hexagon edges
                  for(int j=0; j<(int)interseccoords.size(); j++)
                    for(int k=0; k<(int)interseccoords.size(); k++)
                      if(k<j)
                        for(int l=0; l<3; l++)
                          if(fabs(interseccoords[i](l)-interseccoords[j](l))<1e-7)// components identical
                          {
                            // get edge of j-th triangle within hexagon
                            LINALG::Matrix<3,1> a = interseccoords[k];
                            a -= *cog;
                            LINALG::Matrix<3,1> b = interseccoords[j];
                            b -= *cog;
                            LINALG::Matrix<3,1> c = interseccoords[k];
                            c -= interseccoords[j];

                            double al = a.Norm2();
                            double bl = b.Norm2();
                            double cl = c.Norm2();
                            double alpha = acos((cl*cl+bl*bl-al*al)/(2.0*bl*cl));
                            double h = bl * sin(alpha);

                            hexvolume += cl*h*thickness;
                          }
                  volumes[i] = hexvolume;
                }
                break;
              }
              characlength[i] = 2.0*thickness;
            }
          }
          break;
        }
      }
    }
    // smallest volume
    double minimalvol = 9e99;
    int minimum = 0;
    for(int i=0; i<(int)volumes.size(); i++)
      if(volumes[i]<minimalvol)
      {
        minimalvol = volumes[i];
        structurenumber = i;
        minimum = i;
      }
    if(structurenumber==0 && characlength[0]>=periodlength)
      structurenumber = 3;

    cout<<"\nVolumes: "<<endl;
    for(int i=0; i<(int)volumes.size(); i++)
      cout<<fixed<<setprecision(6)<<"V("<<i<<"): "<<volumes[i]<<"  l_c: "<<characlength[i]<<"  p_cross: "<<crossfraction[i]<<" crosslinks: "<<crosslinksinvolume[i]<<"/"<<numcrossele<<"  niter: "<<niter[i]<<endl;

  /// cout and return network structure
    // write to output files
    // append structure number to DDCorr output
    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "w");
    std::stringstream structuretype;
    structuretype<<structurenumber<<"    "<<characlength[minimum];
    for(int j=0; j<15; j++)
      structuretype<<"    "<<0.0;
    structuretype<<endl;
    fprintf(fp, structuretype.str().c_str());
    fclose(fp);

    // clear the vector first
    testvolumepos_.clear();
    if(numcrossele>0)
    {
      switch(structurenumber)
      {
        // cluster
        case 0:
        {
          cout<<"\nNetwork structure: Cluster"<<endl;
          characlength_ = characlength[structurenumber];
          structuretype_ = structurenumber;

          // calculate the trafo_ matrix (as if we had a layer)

          // adjust second plane direction so that we get an orthonormal basis
          clusterlayervecs[1](0) = clusterlayervecs[0](1)*clusterlayervecs[2](2) - clusterlayervecs[0](2)*clusterlayervecs[2](1);
          clusterlayervecs[1](1) = clusterlayervecs[0](2)*clusterlayervecs[2](0) - clusterlayervecs[0](0)*clusterlayervecs[2](2);
          clusterlayervecs[1](2) = clusterlayervecs[0](0)*clusterlayervecs[2](1) - clusterlayervecs[0](1)*clusterlayervecs[2](0);
          clusterlayervecs[1].Scale(1/clusterlayervecs[1].Norm2()); // hm, not necessary

          // build the base
          for(int i=0; i<trafo_->M(); i++)
            for(int j=0; j<trafo_->N(); j++)
              (*trafo_)(i,j) = clusterlayervecs[i](j);
        }
        break;
        // bundle
        case 1:
        {
          cout<<"\nNetwork structure: Bundle"<<endl;
          cout<<"axis vector: "<<cylvec(0)<<" "<<cylvec(1)<<" "<<cylvec(2)<<endl;
          structuretype_ = structurenumber;
          characlength_ = characlength[structurenumber];
          for(int i=0; i<(int)intersections.size(); i++)
            testvolumepos_.push_back(intersections[i]);

          // save trafo matrix for later use in DDCorrFunction()
          // second direction
          LINALG::Matrix<3,1> raddir1;
          raddir1(0) = 0.0;
          raddir1(1) = 1.0;
          raddir1(2) = cylvec(1)/cylvec(2);
          raddir1.Scale(1.0/raddir1.Norm2());
          // third direction
          LINALG::Matrix<3,1> raddir2;
          raddir2(0) = raddir1(1)*cylvec(2) - raddir1(2)*cylvec(1);
          raddir2(1) = raddir1(2)*cylvec(0) - raddir1(0)*cylvec(2);
          raddir2(2) = raddir1(0)*cylvec(1) - raddir1(1)*cylvec(0);
          raddir2.Scale(1.0/raddir2.Norm2());

          // build the base
          for(int j=0; j<trafo_->N(); j++)
          {
            (*trafo_)(0,j) = cylvec(j);
            (*trafo_)(1,j) = raddir1(j);
            (*trafo_)(2,j)= raddir2(j);
          }
        }
        break;
        // layer
        case 2:
        {
          cout<<"\nNetwork structure: Layer ( ";
          switch((int)(interseccoords.size()))
          {
            case 3: cout<<"triangular shape )"; break;
            case 4: cout<<"rectangular shape )"; break;
            case 6: cout<<"haxagonal shape )"; break;
          }
          cout<<endl;
          for(int i=0; i<(int)layervectors.size(); i++)
            cout<<"layer vector "<<i+1<<": "<<layervectors[i](0)<<" "<<layervectors[i](1)<<" "<<layervectors[i](2)<<endl;
          structuretype_ = structurenumber;
          characlength_ = characlength[structurenumber];
          for(int i=0; i<(int)interseccoords.size(); i++)
            testvolumepos_.push_back(interseccoords[i]);

          // calculate second plane vector (which is now exactly orthogonal to vec1)
          layervectors[1](0) = layervectors[0](1)*layervectors[2](2) - layervectors[0](2)*layervectors[2](1);
          layervectors[1](1) = layervectors[0](2)*layervectors[2](0) - layervectors[0](0)*layervectors[2](2);
          layervectors[1](2) = layervectors[0](0)*layervectors[2](1) - layervectors[0](1)*layervectors[2](0);
          layervectors[1].Scale(1/layervectors[1].Norm2()); // hm, not necessary

          // build the base
          for(int i=0; i<trafo_->M(); i++)
            for(int j=0; j<trafo_->N(); j++)
              (*trafo_)(i,j) = layervectors[i](j);
        }
        break;
        // homogeneous
        case 3:
        {
          cout<<"\nNetwork structure: Homogeneous network"<<endl;
          structuretype_ = structurenumber;
          characlength_ = characlength[0];

          // save trafo matrix for later use in DDCorrFunction()
          for(int i=0; i<trafo_->M(); i++)
            (*trafo_)(i,i) = 1.0;
        }
        break;
      }
    }
    else
    {
      cout<<"\nNetwork structure: Homogeneous network"<<endl;
      structuretype_ = structurenumber;
      characlength_ = characlength[0];

      // save trafo matrix for later use in DDCorrFunction()
      for(int i=0; i<trafo_->M(); i++)
        (*trafo_)(i,i) = 1.0;
    }
  }

  // Communicate trafo_ to other procs
  std::vector<double> localtrafo(9,0.0);
  std::vector<double> globaltrafo(9,0.0);
  for(int i=0; i<trafo_->M(); i++)
    for(int j=0; j<trafo_->N(); j++)
    {
      localtrafo.at(3*i+j) = (*trafo_)(i,j);
      discret_->Comm().SumAll(&localtrafo[3*i+j], &globaltrafo[3*i+j], 1);
      (*trafo_)(i,j) = globaltrafo.at(3*i+j);
    }
  //cout<<*trafo_<<endl;
}//DDCorrCurrentStructure()

/*------------------------------------------------------------------------------*                                                 |
 | density-density correlation function                   (private) mueller 01/11|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::DDCorrIterateVector(const Epetra_Vector& discol, LINALG::Matrix<3,1>* vectorj, const int& maxiterations)
{
  if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2))
    dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

  // get filament number conditions
  vector<DRT::Condition*> filaments(0);
  discret_->GetCondition("FilamentNumber", filaments);
  bool vectorconverged = false;
  int iteration = 0;
  double periodlength = periodlength_->at(0);
  double tolangle = M_PI/180.0; // 1°

  while(!vectorconverged)
  {
    if(iteration<maxiterations)
    {
      LINALG::Matrix<3,1> vectorjp;
      vectorjp.Clear();
      for(int fil=0; fil<(int)filaments.size(); fil++)
      {
        // get next filament
        DRT::Condition* currfilament = filaments[fil];
        for(int node=1; node<(int)currfilament->Nodes()->size(); node++)
        {
          // obtain column map LIDs
          int gid0 = currfilament->Nodes()->at(node-1);
          int gid1 = currfilament->Nodes()->at(node);
          int nodelid0 = discret_->NodeColMap()->LID(gid0);
          int nodelid1 = discret_->NodeColMap()->LID(gid1);
          DRT::Node* node0 = discret_->lColNode(nodelid0);
          DRT::Node* node1 = discret_->lColNode(nodelid1);

          // calculate directional vector between nodes
          LINALG::Matrix<3, 1> dirvec;
          for(int dof=0; dof<3; dof++)
          {
            int dofgid0 = discret_->Dof(node0)[dof];
            int dofgid1 = discret_->Dof(node1)[dof];
            double poscomponent0 = node0->X()[dof]+discol[discret_->DofColMap()->LID(dofgid0)];
            double poscomponent1 = node1->X()[dof]+discol[discret_->DofColMap()->LID(dofgid1)];
            // check for periodic boundary shift and correct accordingly
            if (fabs(poscomponent1 - periodlength - poscomponent0) < fabs(poscomponent1 - poscomponent0))
              poscomponent1 -= periodlength;
            else if (fabs(poscomponent1 + periodlength - poscomponent0) < fabs(poscomponent1 - poscomponent0))
              poscomponent1 += periodlength;

            dirvec(dof) = poscomponent1-poscomponent0;
          }
          // normed vector
          dirvec.Scale(1.0/dirvec.Norm2());

          if(acos(dirvec.Dot((*vectorj)))>(M_PI/2.0))
            dirvec.Scale(-1.0);
          vectorjp += dirvec;
        }
      }
      vectorjp.Scale(1.0/vectorjp.Norm2());
      // check angle between old n_j and n_(j+1)
      double vecvecangle = acos(vectorjp.Dot((*vectorj)));
      if(vecvecangle < tolangle)
      {
        cout<<" vector converged after "<<iteration+1<<" iteration(s) with angle "<<vecvecangle/M_PI*180.0<<" deg"<<endl;
        vectorconverged = true;
      }
      *vectorj = vectorjp;
    }
    else
    {
      cout<<"...Vector did not converge after "<<maxiterations<<" iterations. Continuing..."<<endl;
      vectorconverged = true;
    }
    iteration++;
  }
}//DDCorrIterateVector()
/*------------------------------------------------------------------------------*                                                 |
 | density-density correlation function                  (private) mueller 12/10|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::DDCorrFunction(Epetra_MultiVector& crosslinksperbinrow, Epetra_MultiVector& crosslinksperbinrotrow, LINALG::Matrix<3,1>* centershift)
{
  if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2))
    dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

  int numbins = statmechparams_.get<int>("HISTOGRAMBINS", 1);
  double periodlength = periodlength_->at(0);


/// preliminary operations to set workframe
  // exporter and importer
  Epetra_Export crosslinkexporter(*crosslinkermap_, *transfermap_);
  Epetra_Import crosslinkimporter(*crosslinkermap_, *transfermap_);
  // tranfer map format of crosslinkerbond_ and crosslinkerpositions_
  Epetra_MultiVector crosslinkerbondtrans(*transfermap_,2,true);
  Epetra_MultiVector crosslinkerpositionstrans(*transfermap_,3,true);
  if(discret_->Comm().MyPID()!=0)
  {
    crosslinkerbond_->PutScalar(0.0);
    crosslinkerpositions_->PutScalar(0.0);
  }
  // distribute information to processors
  crosslinkerbondtrans.Export(*crosslinkerbond_, crosslinkexporter, Add);
  crosslinkerpositionstrans.Export(*crosslinkerpositions_, crosslinkexporter, Add);

  // reproduce redundancy prior to Export (not needed in this method, but in the course of DDCorroutput())
  crosslinkerbond_->Import(crosslinkerbondtrans, crosslinkimporter, Insert);
  crosslinkerpositions_->Import(crosslinkerpositionstrans, crosslinkimporter, Insert);

/// preparations for calculation of inter-crosslink distances between central box and surrounding box crosslinker elements
  /* rotate crosslinker position into new (material) coordinate system, then construct
   * new surrounding box crosslinker position with periodic continuations according to the new
   * coordinate directions. This will in some cases lead to a cutoff of crosslinkers which come to lie
   * outside of the rotated volume. These crosslinkers are neglected.
   * The reasons why we are allowed do this, lie with the structure of the network as well as the object of observation.
   * When considering the worst case, we lose 19.6% of the volume (rotation around two axis, 45° each). Still, this does
   * not matter in our case. Since we want to study stationary phases like the cluster or the layer phase,
   * we only consider time intervals where the phase in question does not change its shape anymore. In addition,
   * we try to center the structure as well as possible within the boundary volume by shifting the center of the
   * periodic box. Clusters are not affected by the rotation and cutoff since they constitute a structural singularity,
   * so to speak, and therefore have no contace to the volume boundaries. They are not affected by periodic BCs and are
   * strongly localized density peaks. Layers tend to span over at least one direction of the periodic boundary volume.
   * If this stage is reached, the cutoff will not affect the layer, i.e. its crosslinkers at all, since there is no material
   * in the corners of the boundary volume anymore.
   * Ergo, we wait until this state is reached and then start to evaluate (done via Matlab, not directly implemented).
   * For homogeneous networks, of course, things are different. However, we do not need to worry about it since we have
   * equally distributed filament directions and hence no need to rotate neither the coordinate system nor the entire volume.
   */
/// calculate  shifted as well as shifted and rotated crosslinker positions of the center box
  Epetra_MultiVector centerboxpostrans(*transfermap_,3);
  Epetra_MultiVector centerboxrotpostrans(*transfermap_,3);
  for(int i=0; i<transfermap_->NumMyElements(); i++)
  {
    LINALG::SerialDenseMatrix crosspos(3,1);
    for(int j=0; j<crosspos.M(); j++)
    {
      // 1. original position
      crosspos(j,0) = crosslinkerpositionstrans[j][i];
      // 2. shift, so that crosspos comes to lie within not yet rotated new center box
      if (crosspos(j,0) > periodlength+(*centershift)(j))
        crosspos(j,0) -= periodlength;
      if (crosspos(j,0) < 0.0+(*centershift)(j))
        crosspos(j,0) += periodlength;
      // 2.5 store standard shifted center box position
      centerboxpostrans[j][i] = crosspos(j,0);
      // 3. translate origin
      crosspos(j,0) -= cog_(j);
    }
    // 4. transform the crosslink molecule position into material coordinates
    LINALG::SerialDenseMatrix crossposrot(3,1);
    trafo_->Multiply(false, crosspos, crossposrot);
    for(int j=0; j<centerboxrotpostrans.NumVectors(); j++)
      centerboxrotpostrans[j][i] = crossposrot(j,0);
  }

/// calculate shifted and rotated crosslinker positions for the surrounding boxes (including the center box)
  Epetra_MultiVector positionstrans(*transfermap_, 3*27, true);
  Epetra_MultiVector rotatedpositionstrans(*transfermap_, 3*27, true);
  std::vector<std::vector<int> > boxindices;
  int boxnumber=0;
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
      {
        // box indices
        std::vector<int> ijk;
        ijk.push_back(i);
        ijk.push_back(j);
        ijk.push_back(k);
        boxindices.push_back(ijk);
        // obtain standard crosslinker position and shift it according to current box with respect to box (i=1;j=1;k=1)
        for(int l=0; l<transfermap_->NumMyElements(); l++)
          for(int m=0; m<centerboxpostrans.NumVectors(); m++)
          {
            // shift standard crosslinker positions to box with indices (i,j,k)
            positionstrans[boxnumber*3+m][l] = centerboxpostrans[m][l] + (ijk[m]-1)*periodlength;
            // shift the rotated crosslinker positions to box with indices (i,j,k)
            rotatedpositionstrans[boxnumber*3+m][l] = centerboxrotpostrans[m][l] + (ijk[m]-1)*periodlength;
          }
        boxnumber++;
      }
  // communicate positions for redundancy (crosslinkermap format)
  Epetra_MultiVector positions(*crosslinkermap_, 3*27, true);
  Epetra_MultiVector rotatedpositions(*crosslinkermap_, 3*27, true);
  positions.Import(positionstrans, crosslinkimporter, Insert);
  rotatedpositions.Import(rotatedpositionstrans, crosslinkimporter, Insert);

  // retrieve center of gravity (global coordinates)
  Epetra_SerialDenseMatrix cog(3,1);
  for(int i=0; i<cog.M(); i++)
    cog(i,0) = cog_(i);

/// sort inter-crosslink distances into respective bins
  //parallel brute force from here on
  /* 1. loop over boxes
   * 2. loop over the crosslinkermap format crosslinkers of the boxes
   * 3. loop over the transfermap format crosslinkers of the actual (center box)
   * 4. if(): only calculate deltaxij if the following if both indices i and j stand for crosslinker elements or periodic
   *          images of the those elements
   */
  for(int boxnum=0; boxnum<boxnumber; boxnum++)
    for(int j=0; j<crosslinkermap_->NumMyElements(); j++)
      for(int i=0; i<transfermap_->NumMyElements(); i++)
        if(crosslinkerbondtrans[0][i]>-0.9 && crosslinkerbondtrans[1][i]>-0.9 && (*crosslinkerbond_)[0][j]>-0.9 && (*crosslinkerbond_)[1][j]>-0.9)
        {
          // sort standard (GLOBAL coordinates) inter-crosslink distances into bins
          for(int m=0; m<centerboxpostrans.NumVectors(); m++)
          {
            // centerbox crosslinker component m
            double cboxposm = centerboxpostrans[m][i];
            // surr. box crosslinker component m
            double surrboxposm = positions[boxnum*3+m][j];

            double distm = fabs(surrboxposm-cboxposm);
            int absbin = (int)floor(distm/periodlength*(double)numbins);
            if(absbin==3*numbins)
              absbin--;
            int thebin = absbin%numbins;
            int thecol = (int)floor((double)absbin/(double)numbins)+3*m;
            crosslinksperbinrow[thecol][thebin] += 1.0;
          }

          // sort inter-crosslink distances of ROTATED system into bins
          for(int m=0; m<centerboxrotpostrans.NumVectors(); m++)
          {
            //centerbox crosslinker component m
            double cboxposrotm = centerboxrotpostrans[m][i];
            // surrounding box crosslinker position component m
            double surrboxposrotm = rotatedpositions[boxnum*3+m][j];
            int surrboxindex = boxindices[boxnum][m];
            // check if both crosslinkers lie within the shifted and rotated boundary volume
            double deltacbox = fabs(cboxposrotm-(*centershift)(m));
            double deltasbox = fabs(surrboxposrotm-((*centershift)(m)+(surrboxindex-1)*periodlength));
            if(deltacbox<=periodlength/2.0 && deltasbox<=periodlength/2.0) // inside the rotated volume
            {
              // determine bin for rotated fixed system coordinates
              double distm = fabs(surrboxposrotm-cboxposrotm);
              int thebin = (int)floor(distm/periodlength*(double)numbins);
              int thecol = m;
              // only distances [0;H[
              if(thebin<numbins)
                crosslinksperbinrotrow[thecol][thebin] += 1.0;
            }
          }
        }
  return;
}//STATMECH::StatMechManager::DDCorrFunction()

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
        unbindingtimes << (*crosslinkunbindingtimes_)[1][i] <<endl;
    }
    fprintf(fp, unbindingtimes.str().c_str());
    fclose(fp);
  }
  return;
}

/*------------------------------------------------------------------------------*                                                 |
 | Output of distances between doubly bound linkers of the horizontal filament  |
 | in a loom type network                                                       |
 |                                                       (private) mueller 2/12 |
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::LoomOutput(const Epetra_Vector& disrow,
                                           const std::ostringstream& nearestneighborfilename,
                                           const std::ostringstream& linkerpositionsfilename)
{
  if(!DRT::INPUT::IntegralValue<int>(statmechparams_, "LOOMSETUP"))
      dserror("For Loom related output, activate LOOMSETUP in your input file!");

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
            linkerpositions<<std::scientific<<std::setprecision(6)<<(pos0->second)(0)<<"\t"<<(pos0->second)(1)<<"\t"<<(pos0->second)(2)<<"\t"<<globalarclength<<endl;

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
            distances<<endl;
          }
          // time step divider line
          linkerpositions<<std::scientific<<std::setprecision(6)<<-1e9<<"\t"<<-1e9<<"\t"<<-1e9<<"\t"<<-1e9<<endl;
        }
      }
    }
    fprintf(fp, distances.str().c_str());
    fclose(fp);
    fprintf(fp2, linkerpositions.str().c_str());
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
  if(!DRT::INPUT::IntegralValue<int>(statmechparams_, "LOOMSETUP"))
    dserror("For Loom related output, activate LOOMSETUP in your input file!");

  for(int i=0; i<discret_->NumMyRowElements(); i++)
  {
    DRT::Element* element = discret_->lRowElement(i);
    const DRT::ElementType & eot = element->ElementType();

    // assumption: the Truss3 element has one explicit owner and therefore occurs only once in the row map
    if(eot == DRT::ELEMENTS::Truss3Type::Instance())
    {
      // Get information on Boundary Conditions, i.e. the number of spatial dimensions (2D,3D)
      vector<DRT::Condition*> dirichlet;
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
      Teuchos::RCP<Epetra_SerialDenseVector> force = dynamic_cast<DRT::ELEMENTS::Truss3*>(element)->InternalForces();
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
      internalforce << std::scientific << std::setprecision(15) << distance.Norm2()<< "  " << (distance.Norm2()-l0)/l0<<"  "<< fint(0) <<endl;

      fprintf(fp, internalforce.str().c_str());
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
/*------------------------------------------------------------------------------*                                                 |
 | simply counts the number of free, one-bonded, and two-bonded crosslinkers    |
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
  const RCP<Epetra_Vector> disp = Teuchos::rcp(new Epetra_Vector(disrow));
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
    elasticenergy<<endl;
    fprintf(fp, elasticenergy.str().c_str());
    fclose(fp);
  }
  discret_->Comm().Barrier();
  return;
}

/*------------------------------------------------------------------------------*                                                 |
 | output nodal displacements                                      mueller 5/12 |
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::OutputNodalDisplacements(const Epetra_Vector&                 disrow,
                                                         const std::ostringstream&            filename)
{
  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(disrow, discol);

  if(!discret_->Comm().MyPID())
  {
    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "a");

    std::stringstream dispnode;
    // retrieve translational node displacements
    for(int i=0; i<discret_->DofColMap()->NumMyElements(); i=i+6)
        dispnode << discol[i]<<" "<< discol[i+1]<<" "<< discol[i+2]<<endl;
    fprintf(fp, dispnode.str().c_str());
    fclose(fp);
  }
  return;
}

/*------------------------------------------------------------------------------*                                                 |
 | output the coverage of crosslinker binding sites (nodes) and the distribution|
 |of bound crosslinkers                                  (private) mueller 5/12 |
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::CrosslinkCoverageOutput(const Epetra_Vector& disrow, const std::ostringstream& filename, bool coverageonly)
{
  if(coverageonly)
  {
    RCP<Epetra_Vector> bspotstatustrans = Teuchos::rcp(new Epetra_Vector(*bspotrowmap_));
    CommunicateVector(bspotstatustrans, bspotstatus_,true,false,false,true);

    // check for occupied binding spots
    int partialsum = 0;
    for(int i=0; i<bspotstatustrans->MyLength(); i++)
    {
      if((*filamentnumber_)[bspotcolmap_->LID(bspotrowmap_->GID(i))]==0 && (*bspotstatustrans)[i]>-0.1)
        partialsum++;
      else if((int)(*filamentnumber_)[bspotcolmap_->LID(bspotrowmap_->GID(i))]>0)
        break;
    }
    // sum up processor-specific values and communicate
    int sum = 0;
    discret_->Comm().SumAll(&partialsum,&sum,1);

    if(!discret_->Comm().MyPID())
    {
      // file pointer
      FILE* fp = NULL;
      fp = fopen(filename.str().c_str(), "a");
      std::stringstream coverage;

      coverage << sum << " "<< std::setprecision(5) <<(double)(sum)/(double)(numbspots_) <<endl;
      // print to file and close
      fprintf(fp, coverage.str().c_str());
      fclose(fp);
    }
  }
  else
  {
    /* Consider the horizontal filament of a loom setup. Go along this filament and count the
     * number of occupied binding spots. Also, store the spatial distribution of occupied spots
     * for analysis of cluster size, etc.*/
    Epetra_Vector discol(*discret_->DofColMap(),true);
    LINALG::Export(disrow,discol);
    std::map<int, LINALG::Matrix<3, 1> > currentpositions;
    std::map<int, LINALG::Matrix<3, 1> > currentrotations;
    GetNodePositionsFromDisVec(discol, currentpositions, currentrotations, true);

    if(!discret_->Comm().MyPID())
    {
      // file pointer
      FILE* fp = NULL;
      fp = fopen(filename.str().c_str(), "a");
      std::stringstream coverage;

      for(int i=0; i<bspotstatus_->MyLength(); i++)
      {
        // consider only filament 0 (horizontal filament)
        if((int)(*filamentnumber_)[i]==0 && (*bspotstatus_)[i]>-0.1)
          // store both singly and doubly bound linkers
          coverage<<bspotcolmap_->GID(i)<<"  "<<(int)(*numbond_)[(int)(*bspotstatus_)[i]]<<"  "<< std::scientific << std::setprecision(15) <<(currentpositions.find(i)->second)(0)<<endl;
        else if((int)(*filamentnumber_)[i]>0)
          break;
      }
      // print to file and close
      fprintf(fp, coverage.str().c_str());
      fclose(fp);
    }
  }
  return;
}

/*------------------------------------------------------------------------------*                                                 |
 | simply counts the number of free, one-bonded, and two-bonded crosslinkers    |
 |                                                        (public) mueller 4/11 |
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::CrosslinkCount(const std::ostringstream& filename)
{
  if(discret_->Comm().MyPID()==0)
  {
    int free = 0;
    int onebond = 0;
    int twobond = 0;

    for(int i=0; i<crosslinkerbond_->MyLength(); i++)
    {
      // count number of node IDs in crosslinkerbond_ entry i
      int numofnodes = 0;
      for(int j=0; j<crosslinkerbond_->NumVectors(); j++)
        if((*crosslinkerbond_)[j][i]>-0.9)
          numofnodes++;
      // increment according to numofnodes
      switch(numofnodes)
      {
        // free
        case 0:
          free++;
        break;
        // crosslink with one bond
        case 1:
          onebond++;
        break;
        // crosslink with two bonds
        case 2:
          twobond++;
        break;
      }
    }

    // write to file
    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "a");
    std::stringstream ccount;
    ccount<<free<<"    "<<onebond<<"    "<<twobond<<"    ";
    for(int i=0; i<13; i++)
      ccount<<"    0";
    ccount<<endl;
    fprintf(fp, ccount.str().c_str());
    fclose(fp);
  }

  return;
}

/*------------------------------------------------------------------------------*                                                 |
 | linker spot counter and check of interfilament orient. (public) mueller 12/10|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::OrientationCorrelation(const Epetra_Vector& disrow, const int &istep)
{
  /* what this does:
   * 1) orientation correlation function for all proximal binding spot pairs (regardless of orientation)
   * 2) count the number of overall possible binding spots for crosslinkers (proximal and correctly oriented)
   * 3) order parameter
   */
  if (DRT::INPUT::IntegralValue<int>(statmechparams_, "CHECKORIENT"))
  {
    if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2))
      dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

    Epetra_Vector discol(*discret_->DofColMap(), true);
    LINALG::Export(disrow, discol);

    Teuchos::RCP<Epetra_MultiVector> bspotpositions = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));
    Teuchos::RCP<Epetra_MultiVector> bspotrotations = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));
    GetBindingSpotPositions(discol,bspotpositions, bspotrotations);

    // NODAL TRIAD UPDATE
    Teuchos::RCP<Epetra_MultiVector> bspottriadscol = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,4,true));
    GetBindingSpotTriads(bspotrotations, bspottriadscol);

    // distance and orientation checks
    int numbins = statmechparams_.get<int>("HISTOGRAMBINS", 1);
    // max. distance between two binding spots
    double maxdist = periodlength_->at(0)*sqrt(3.0);
    // max. angle
    double maxangle = M_PI/2.0;
    // minimal and maximal linker search radii
    double rmin = statmechparams_.get<double>("R_LINK",0.0)-statmechparams_.get<double>("DeltaR_LINK",0.0);
    double rmax = statmechparams_.get<double>("R_LINK",0.0)+statmechparams_.get<double>("DeltaR_LINK",0.0);
    // number of overall independent combinations
    int numnodes = discret_->NodeColMap()->NumMyElements();
    int numcombinations = (numnodes*numnodes-numnodes)/2;
    // combinations on each processor
    int combinationsperproc = (int)floor((double)numcombinations/(double)discret_->Comm().NumProc());
    int remainder = numcombinations%combinationsperproc;

    int bindingspots = 0;
    int combicount = 0;
    Teuchos::RCP<Epetra_Vector> anglesrow = Teuchos::rcp(new Epetra_Vector(*ddcorrrowmap_, true));
    // including correlation on the same filament
    Teuchos::RCP<Epetra_MultiVector> orderparameterbinsincrow = Teuchos::rcp(new Epetra_MultiVector(*ddcorrrowmap_,2, true));
    // excluding correlation on identical filament
    Teuchos::RCP<Epetra_MultiVector> orderparameterbinsexcrow = Teuchos::rcp(new Epetra_MultiVector(*ddcorrrowmap_,2, true));
    // loop over crosslinkermap_ (column map, same for all procs: maps all crosslink molecules)
    for(int mypid=0; mypid<discret_->Comm().NumProc(); mypid++)
    {
      bool quitloop = false;
      if(mypid==discret_->Comm().MyPID())
      {
        bool continueloop = false;
        int appendix = 0;
        if(mypid==discret_->Comm().NumProc()-1)
          appendix = remainder;

        for(int i=0; i<bspotcolmap_->NumMyElements(); i++)
        {
          for(int j=0; j<bspotcolmap_->NumMyElements(); j++)
          {
            // start adding crosslink from here
            if(i==(*startindex_)[2*mypid] && j==(*startindex_)[2*mypid+1])
              continueloop = true;
            // only entries above main diagonal and within limits of designated number of crosslink molecules per processor
            if(j>i && continueloop)
            {
              if(combicount<combinationsperproc+appendix)
              {
                combicount++;

                Epetra_SerialDenseMatrix LID(2,1);
                LID(0,0) = i;
                LID(1,0) = j;

                LINALG::Matrix<3,1> distance;
                for(int k=0; k<(int)distance.M(); k++)
                  distance(k) = (*bspotpositions)[k][i]-(*bspotpositions)[k][j];

                // current distance bin
                int currdistbin = (int)floor(distance.Norm2()/maxdist*numbins);
                // reduce bin if distance == maxdist
                if(currdistbin==numbins)
                  currdistbin--;

                // direction between currently considered two nodes
                LINALG::Matrix<3,1> direction(distance);
                direction.Scale(1.0/direction.Norm2());
                Teuchos::RCP<double> phifil = Teuchos::rcp(new double(0.0));
                bool orientation = CheckOrientation(direction,*bspottriadscol,LID,phifil);

                // increment count for that bin
                (*orderparameterbinsincrow)[0][currdistbin] += 1.0;
                // order parameter
                (*orderparameterbinsincrow)[1][currdistbin] += (3.0*cos((*phifil))*cos((*phifil))-1)/2.0;
                if((*filamentnumber_)[i]!=(*filamentnumber_)[j])
                {
                  (*orderparameterbinsexcrow)[0][currdistbin] += 1.0;
                  (*orderparameterbinsexcrow)[1][currdistbin] += (3.0*cos((*phifil))*cos((*phifil))-1)/2.0;
                }
                // proximity check
                if(distance.Norm2()>rmin && distance.Norm2()<rmax)
                {
                  // if angular constraints are met, increase binding spot count
                  if(orientation)
                    bindingspots++;

                  // determine the bin
                  int curranglebin = (int)floor((*phifil)/maxangle*numbins);
                  // in case the distance is exactly periodlength*sqrt(3)
                  if(curranglebin==numbins)
                    curranglebin--;
                  (*anglesrow)[curranglebin] += 1.0;
                }
              }
              else
              {
                quitloop = true;
                break;
              }
            }
          }
          if(quitloop)
            break;
        }
        if(quitloop)
          break;
      }
      else
        continue;
    }
    // Export
    // add up binding spots
    int bspotsglob = 0;
    discret_->Comm().SumAll(&bindingspots,&bspotsglob,1);

    Teuchos::RCP<Epetra_Vector> anglescol = Teuchos::rcp(new Epetra_Vector(*ddcorrcolmap_, true));
    Teuchos::RCP<Epetra_MultiVector> orderparameterbinsinccol = Teuchos::rcp(new Epetra_MultiVector(*ddcorrcolmap_,2, true));
    Teuchos::RCP<Epetra_MultiVector> orderparameterbinsexccol = Teuchos::rcp(new Epetra_MultiVector(*ddcorrcolmap_,2, true));
    CommunicateVector(anglesrow,anglescol, false, true,false);
    CommunicateMultiVector(orderparameterbinsincrow, orderparameterbinsinccol,false,true,false);
    CommunicateMultiVector(orderparameterbinsexcrow, orderparameterbinsexccol,false,true,false);

    // Add the processor-specific data up
    std::vector<int> angles(numbins, 0);
    std::vector<std::vector<double> > orderparameterinc(numbins, std::vector<double>(2,0.0));
    std::vector<std::vector<double> > orderparameterexc(numbins, std::vector<double>(2,0.0));
    for(int i=0; i<numbins; i++)
      for(int pid=0; pid<discret_->Comm().NumProc(); pid++)
      {
        angles[i] += (int)(*anglescol)[pid*numbins+i];
        orderparameterinc[i][0] += (*orderparameterbinsinccol)[0][pid*numbins+i];
        orderparameterinc[i][1] += (*orderparameterbinsinccol)[1][pid*numbins+i];
        orderparameterexc[i][0] += (*orderparameterbinsexccol)[0][pid*numbins+i];
        orderparameterexc[i][1] += (*orderparameterbinsexccol)[1][pid*numbins+i];
      }

    // average values
    for(int i=0; i<numbins; i++)
      if(orderparameterinc[i][0]>0.0) // i.e. >0
      {
        orderparameterinc[i][1] /= orderparameterinc[i][0];
        orderparameterexc[i][1] /= orderparameterexc[i][0];
      }
      else
      {
        orderparameterinc[i][1] = -99.0;
        orderparameterexc[i][1] = -99.0;
      }

    // write data to file
    if(!discret_->Comm().MyPID())
    {
      std::ostringstream orientfilename;
      orientfilename << outputrootpath_ << "/StatMechOutput/LinkerSpotsOrCorr_"<<std::setw(6) << std::setfill('0') << istep <<".dat";

      FILE* fp = NULL;
      fp = fopen(orientfilename.str().c_str(), "w");
      std::stringstream histogram;

      histogram<<bspotsglob<<"    "<<-99<<"    "<<-99<<endl;
      for(int i=0; i<numbins; i++)
        histogram<<i+1<<"    "<<angles[i]<<"    "<<std::setprecision(12)<<orderparameterinc[i][1]<<"    "<<orderparameterinc[i][0]<<"    "<<orderparameterexc[i][1]<<"    "<<orderparameterexc[i][0]<<endl;
      //write content into file and close it
      fprintf(fp, histogram.str().c_str());
      fclose(fp);
    }
  }
}//NumLinkerSpotsAndOrientation()

/*------------------------------------------------------------------------------*                                                 |
 | computes the mesh size of the network dep. on distance to cog                |
 |                                                        (public) mueller 12/10|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::ComputeLocalMeshSize(const Epetra_Vector& disrow, LINALG::Matrix<3,1>& centershift, const int &istep)
{
  if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2))
    dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

  double periodlength = periodlength_->at(0);
  double maxdist = periodlength/2.0*sqrt(3.0);
  int numbins = statmechparams_.get<int>("HISTOGRAMBINS",1);
  // center of gravity
  LINALG::Matrix<3,1> cog;
  for(int i=0; i<(int)cog.M(); i++)
    cog(i) = centershift(i) + periodlength/2.0;

  // 1. set up of a row map vector containing the global node positions
  std::vector<LINALG::Matrix<3,1> > xrelnodes;
  xrelnodes.clear();
  for (int i=0; i<discret_->NumMyRowNodes(); i++)
  {
    //get pointer at a node
    const DRT::Node* node = discret_->lRowNode(i);
    //get GIDs of this node's degrees of freedom
    std::vector<int> dofnode = discret_->Dof(node);
    // global position and shift according to new centershift
    LINALG::Matrix<3, 1> xglob;

    // shift according to given center (cog)
    for(int j=0; j<(int)xglob.M(); j++)
    {
      // get the global node position
      xglob(j) = node->X()[j] + disrow[discret_->DofRowMap()->LID(dofnode[j])];
      // shift the j-th component according to new center
      if (xglob(j) > periodlength+centershift(j))
        xglob(j) -= periodlength;
      if (xglob(j) < centershift(j))
        xglob(j) += periodlength;
    }
    xglob -= cog;
    xrelnodes.push_back(xglob);
  }

  Epetra_Vector fillengthrow(*ddcorrrowmap_, true);
  // calculate sqrt(DV/DL)
  for(int i=1; i<discret_->NumMyRowNodes(); i++)
  {
    // column map LID of the two nodes in question
    int gid0 = discret_->NodeRowMap()->GID(i-1);
    int gid1 = discret_->NodeRowMap()->GID(i);
    int collid0 = discret_->NodeColMap()->LID(gid0);
    int collid1 = discret_->NodeColMap()->LID(gid1);

    // make sure both nodes lie on the same filament
    if((*filamentnumber_)[collid0] == (*filamentnumber_)[collid1])
    {
      // calculate node bins
      int index0 = i-1;
      int index1 = i;
      int bin0 = (int)floor(xrelnodes[i-1].Norm2()/maxdist*numbins);
      int bin1 = (int)floor(xrelnodes[i].Norm2()/maxdist*numbins);
      if(bin0 == numbins)
        bin0--;
      if(bin1 == numbins)
        bin1--;

      // case: the two nodes lie within different bins: add length segments binwise
      if(bin0 != bin1)
      {
        // switch in order to follow convention
        if(bin1<bin0)
        {
          index0 = i;
          index1 = i-1;
          bin0 = (int)floor(xrelnodes[i].Norm2()/maxdist*numbins);
          bin1 = (int)floor(xrelnodes[i-1].Norm2()/maxdist*numbins);
        }

        // check, if element is broken and unshift in order to obtain correct directional vector
        // calculate directional vector between nodes
        LINALG::Matrix<3,1> unshift = xrelnodes[index1];

        for(int j=0; j<(int)xrelnodes[index0].M(); j++)
        {
          // check for periodic boundary shift and correct accordingly
          if (fabs(xrelnodes[index1](j) - periodlength - xrelnodes[index0](j)) < fabs(xrelnodes[index1](j) - xrelnodes[index0](j)))
            unshift(j) -= periodlength;
          else if (fabs(xrelnodes[index1](j) + periodlength - xrelnodes[index0](j)) < fabs(xrelnodes[index1](j) - xrelnodes[index0](j)))
            unshift(j) += periodlength;
        }


        LINALG::Matrix<3,1> dirvec = unshift;
        dirvec -= xrelnodes[index0];
        // directional unit vector
        dirvec.Scale(1.0/dirvec.Norm2());

        // number of bin intersections
        int numisecs = abs(bin1-bin0);
        // starting point
        LINALG::Matrix<3,1> xstartj = xrelnodes[index0];

        /*-----------------------------------
         * short example:
         *    bin:   16    17    18    19
         *         |  o--|-----|-----|--o  |
         * abslim: 16    17    18    19    20
         * binlim:       0     1     2
          -----------------------------------*/
        //cout<<"numisecs = "<<numisecs<<endl;
        for(int j=0; j<numisecs; j++)
        {
          // all segments except last segment
          if(j<=numisecs-1)
          {
            int currbin = bin0+j;
            // j-th bin limit
            double nextbinlimit = (double)(currbin+1)/(double)numbins*maxdist;

            // line parameter mu for intersection of line with spherical shell:
            // note: we take the solution that is (may be?) >=0 since we chose the directional vector accordingly

            //cout<<"  j="<<j<<", xstartj: "<<xstartj(0)<<" "<<xstartj(1)<<" "<<xstartj(2)<<endl;
            double a = xstartj.Dot(dirvec);
            double b = sqrt(a*a-xstartj.Norm2()*xstartj.Norm2()+nextbinlimit*nextbinlimit);
            double mu = -a + b;

            //cout<<"  mu("<<j<<") = "<<mu<<endl;

            // intersection of the line with the spherical shell
            LINALG::Matrix<3,1> isection = dirvec;
            isection.Scale(mu);
            isection += xstartj;
            // calculate element (filament) segment added to the respective bin
            LINALG::Matrix<3,1> segment = isection;
            segment -= xstartj;
            fillengthrow[currbin] += segment.Norm2();
            // store current intersection in case there is more than one
            xstartj = isection;
            //if(numisecs>0)
            //  cout<<"    j="<<j+1<<", previs: "<<xstartj(0)<<" "<<xstartj(1)<<" "<<xstartj(2)<<endl;
          }
          // last segment
          if(j==numisecs-1)
          {
            //cout<<"  xfinale: "<<xstartj(0)<<" "<<xstartj(1)<<" "<<xstartj(2)<<endl;
            LINALG::Matrix<3,1> segment = xrelnodes[index1];
            segment -= xstartj;
            fillengthrow[bin1] += segment.Norm2();
            xstartj.Clear();
          }
        }
      }
      else // just add the entire element length
      {
        LINALG::Matrix<3,1> elength = xrelnodes[i];
        elength -= xrelnodes[i-1];
        fillengthrow[bin0] += elength.Norm2();
      }
    }
  }

  Epetra_Vector fillengthcol(*ddcorrcolmap_,true);
  Epetra_Import importer(*ddcorrcolmap_, *ddcorrrowmap_);
  fillengthcol.Import(fillengthrow,importer,Insert);
  // Add the processor-specific data up
  std::vector<double> fillength(numbins, 0.0);
  for(int i=0; i<numbins; i++)
    for(int pid=0; pid<discret_->Comm().NumProc(); pid++)
      fillength[i] += fillengthcol[pid*numbins+i];

  // Proc 0 section
  if(discret_->Comm().MyPID()==0)
  {
    std::ostringstream filename;
    filename << outputrootpath_ << "/StatMechOutput/LocalMeshSize_"<<std::setw(6) << std::setfill('0') << istep <<".dat";

    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "w");
    std::stringstream histogram;

    for(int i=0; i<numbins; i++)
      histogram<<i+1<<"    "<<std::setprecision(12)<<fillength[i]<<endl;

    //write content into file and close it
    fprintf(fp, histogram.str().c_str());
    fclose(fp);
  }
}

/*------------------------------------------------------------------------------*                                                 |
| Viscoelasticity ouput                                  (public) mueller 11/12|
*------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::ViscoelasticityOutput(const double& time, const Epetra_Vector& dis, const Epetra_Vector& fint, std::ostringstream& filename)
{
#ifdef DEBUG
  if (forcesensor_ == Teuchos::null)
    dserror("forcesensor_ is NULL pointer; possible reason: dynamic crosslinkers not activated and forcesensor applicable in this case only");
#endif  // #ifdef DEBUG
  double f = 0;//mean value of force
  double d = 0;//Displacement

  for(int i=0; i<forcesensor_->MyLength(); i++)// forcesensor_ is unique on each Proc (see UpdateForceSensors() !)
    if((*forcesensor_)[i]>0.9)
    {
      // translate i to DofRowMap LID
      int dofgid = discret_->DofColMap()->GID(i);
      int rowid = discret_->DofRowMap()->LID(dofgid);
      f += fint[rowid];
      d = dis[rowid];
    }

  //f is the sum of all forces at the top on this processor; compute the sum fglob on all processors all together
  double fglob = 0;
  discret_->Comm().SumAll(&f,&fglob,1);

  if(!discret_->Comm().MyPID())
  {
    //pointer to file into which each processor writes the output related with the dof of which it is the row map owner
    FILE* fp = NULL;

    //content to be written into the output file
    std::stringstream filecontent;

    fp = fopen(filename.str().c_str(), "a");

    /*the output to be written consists of internal forces at exactly those degrees of freedom
     * marked in *forcesensor_ by a one entry*/

    filecontent << std::scientific << std::setprecision(10) << time;//changed

    //Putting time, displacement, meanforce  in Filestream
    filecontent << "   "<< d << "   " << fglob << "   " << discret_->NumGlobalElements() << endl; //changed
    //writing filecontent into output file and closing it
    fprintf(fp,filecontent.str().c_str());
    fclose(fp);
  }
}

/*------------------------------------------------------------------------------*                                                 |
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
        filecontent<<std::scientific<<std::setprecision(6)<<(*unbindingprobability_)[0][i]<<" "<<(*unbindingprobability_)[1][i]<<" "<<(*crosslinkerpositions_)[0][i]<<" "<<(*crosslinkerpositions_)[1][i]<<" "<<(*crosslinkerpositions_)[2][i]<<endl;

    fprintf(fp,filecontent.str().c_str());
    fclose(fp);
  }

  return;
}

/*------------------------------------------------------------------------------*                                                 |
| distribution of spherical coordinates                  (public) mueller 12/10|
*------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::SphericalCoordsDistribution(const Epetra_Vector& disrow,
                                                  Epetra_Vector& phibinsrow,
                                                  Epetra_Vector& thetabinsrow,
                                                  Epetra_Vector& costhetabinsrow)
{
  if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2))
    dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

  int numbins = statmechparams_.get<int>("HISTOGRAMBINS", 1);
  double periodlength = periodlength_->at(0);

  for(int i=0; i<discret_->NumMyRowElements(); i++)
  {
    DRT::Element* element = discret_->lRowElement(i);
    // consider filament elements only
    if(element->Id()<basisnodes_)
    {
      int gid0 = element->Nodes()[0]->Id();
      int gid1 = element->Nodes()[1]->Id();
      int lid0 = discret_->NodeRowMap()->LID(gid0);
      int lid1 = discret_->NodeRowMap()->LID(gid1);
      DRT::Node* node0 = discret_->lRowNode(lid0);
      DRT::Node* node1 = discret_->lRowNode(lid1);

      // calculate directional vector between nodes
      double dirlength = 0.0;
      Epetra_SerialDenseMatrix dirvec(3,1);
      for(int dof=0; dof<dirvec.M(); dof++)
      {
        int dofgid0 = discret_->Dof(node0)[dof];
        int dofgid1 = discret_->Dof(node1)[dof];
        double poscomponent0 = node0->X()[dof]+disrow[discret_->DofRowMap()->LID(dofgid0)];
        double poscomponent1 = node1->X()[dof]+disrow[discret_->DofRowMap()->LID(dofgid1)];
        // check for periodic boundary shift and correct accordingly
        if (fabs(poscomponent1 - periodlength - poscomponent0) < fabs(poscomponent1 - poscomponent0))
          poscomponent1 -= periodlength;
        else if (fabs(poscomponent1 + periodlength - poscomponent0) < fabs(poscomponent1 - poscomponent0))
          poscomponent1 += periodlength;

        dirvec(dof,0) = poscomponent1-poscomponent0;
        dirlength += dirvec(dof,0)*dirvec(dof,0);
      }
      dirlength = sqrt(dirlength);
      // normed directional vector
      dirvec.Scale(1.0/dirlength);

      Epetra_SerialDenseMatrix dirvecrot(3,1,true);
      trafo_->Multiply(false,dirvec,dirvecrot);

      // transform into spherical coordinates (phi E [-pi;pi], theta E [0; pi]) and sort into appropriate bin
      double phi = atan2(dirvecrot(1,0),dirvecrot(0,0)) + M_PI;
      double theta = acos(dirvecrot(2,0));

      int phibin = (int)floor(phi/(2*M_PI)*numbins);
      int thetabin = (int)floor(theta/M_PI*numbins);
      int costhetabin = (int)floor((cos(theta)+1.0)/2.0*numbins);
      if(phibin == numbins)
        phibin--;
      if(thetabin == numbins)
        thetabin--;
      if(costhetabin == numbins)
        costhetabin--;
      if(phibin<0 || thetabin<0)
        dserror("bin smaller zero");
      phibinsrow[phibin] += 1.0;
      thetabinsrow[thetabin] += 1.0;
      costhetabinsrow[costhetabin] += 1.0;
    }
  }
  return;
}//SphericalCoordsDistribution()

/*------------------------------------------------------------------------------*
 | radial crosslinker density distribution                (public) mueller 12/10|
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::RadialDensityDistribution(Epetra_Vector& radialdistancesrow, LINALG::Matrix<3,1>& centershift)
{
  if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2))
    dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

  // simpler version taking into account only the original boundary volume
  double periodlength = periodlength_->at(0);
  double maxdistance = periodlength/2.0*sqrt(3.0);
  int numbins = statmechparams_.get<int>("HISTOGRAMBINS", 1);

  LINALG::Matrix<3,1> cog;
  for(int i=0; i<(int)cog.M(); i++)
    cog(i) = centershift(i) + periodlength/2.0;

  // Export to transfer map format
  Epetra_MultiVector crosslinkerbondtrans(*transfermap_, 2, true);
  Epetra_MultiVector crosslinkerpositionstrans(*transfermap_, 3, true);
  Epetra_Export crosslinkexporter(*crosslinkermap_, *transfermap_);
  Epetra_Export crosslinkimporter(*crosslinkermap_, *transfermap_);
  // note: it seems to be necessary to clear all vectors other than on Proc 0.
  // otherwise, i.e. using method "Insert" on Export, from time to time (no clear pattern has emerged so far)
  // incorrect data is written to transfer vectors. Odd!
  if(discret_->Comm().MyPID()!=0)
  {
    crosslinkerbond_->PutScalar(0.0);
    crosslinkerpositions_->PutScalar(0.0);
  }
  //Export
  crosslinkerbondtrans.Export(*crosslinkerbond_, crosslinkexporter, Add);
  crosslinkerpositionstrans.Export(*crosslinkerpositions_, crosslinkexporter, Add);
  // Reimport to create pre-Export status on all procs
  crosslinkerbond_->Import(crosslinkerbondtrans, crosslinkimporter, Insert);
  crosslinkerpositions_->Import(crosslinkerpositionstrans, crosslinkimporter, Insert);

  // calculate distance of crosslinkers to center of gravity
  for(int i=0; i<crosslinkerbondtrans.MyLength(); i++)
    if(crosslinkerbondtrans[0][i]>-0.9 && crosslinkerbondtrans[1][i]>-0.9)
    {
      LINALG::Matrix<3,1> distance;
      // shift coordinates according to new boundary box center
      for(int j=0; j<crosslinkerpositionstrans.NumVectors(); j++)
      {
        distance(j) = crosslinkerpositionstrans[j][i];
        if (distance(j) > periodlength+centershift(j))
          distance(j) -= periodlength;
        if (distance(j) < centershift(j))
          distance(j) += periodlength;
      }
      distance -= cog;
      // determine histogram bin for current crosslinker element
      int currbin = (int)floor(distance.Norm2()/maxdistance * numbins);
      if(currbin==numbins)
        currbin--;
      radialdistancesrow[currbin] += 1.0;
    }
  return;
}//RadialDensityDistribution()

/*----------------------------------------------------------------------*
 | filament orientations and output               (public) mueller 07/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::FilamentOrientations(const Epetra_Vector& discol, std::vector<LINALG::Matrix<3,1> >* normedvectors, const std::ostringstream& filename, bool fileoutput)
{
  /* Output of filament element orientations (Proc 0 only):
   * format: filamentnumber    d_x  d_y  d_z
   */

  if(discret_->Comm().MyPID()==0)
  {
    if(periodlength_->at(0) != periodlength_->at(1) || periodlength_->at(0) != periodlength_->at(2))
      dserror("For this analysis, we require a cubic periodic box! In your input file, PERIODLENGTH = [ %4.2f, %4.2f, %4.2f]", periodlength_->at(0), periodlength_->at(1), periodlength_->at(2));

    double periodlength = periodlength_->at(0);

    FILE* fp = NULL;
    if(fileoutput)
      fp = fopen(filename.str().c_str(), "w");
    std::stringstream fileleorientation;

    // get filament number conditions
    vector<DRT::Condition*> filaments(0);
    discret_->GetCondition("FilamentNumber", filaments);

    for(int fil=0; fil<(int)filaments.size(); fil++)
    {
      // get next filament
      DRT::Condition* currfilament = filaments[fil];
      for(int node=1; node<(int)currfilament->Nodes()->size(); node++)
      {
        // obtain column map LIDs
        int gid0 = currfilament->Nodes()->at(node-1);
        int gid1 = currfilament->Nodes()->at(node);
        int nodelid0 = discret_->NodeColMap()->LID(gid0);
        int nodelid1 = discret_->NodeColMap()->LID(gid1);
        DRT::Node* node0 = discret_->lColNode(nodelid0);
        DRT::Node* node1 = discret_->lColNode(nodelid1);

        // calculate directional vector between nodes
        LINALG::Matrix<3, 1> dirvec;
        for(int dof=0; dof<3; dof++)
        {
          int dofgid0 = discret_->Dof(node0)[dof];
          int dofgid1 = discret_->Dof(node1)[dof];
          double poscomponent0 = node0->X()[dof]+discol[discret_->DofColMap()->LID(dofgid0)];
          double poscomponent1 = node1->X()[dof]+discol[discret_->DofColMap()->LID(dofgid1)];
          // check for periodic boundary shift and correct accordingly
          if (fabs(poscomponent1 - periodlength - poscomponent0) < fabs(poscomponent1 - poscomponent0))
            poscomponent1 -= periodlength;
          else if (fabs(poscomponent1 + periodlength - poscomponent0) < fabs(poscomponent1 - poscomponent0))
            poscomponent1 += periodlength;

          dirvec(dof) = poscomponent1-poscomponent0;
        }
        // normed vector
        dirvec.Scale(1/dirvec.Norm2());

        // add element directional vectors up
        for(int i=0; i<(int)normedvectors->size(); i++)
        {
          // base vector
          LINALG::Matrix<3,1> ei;
          LINALG::Matrix<3,1> vi = dirvec;

          ei.Clear();
          ei(i) = 1.0;

          if(acos(vi.Dot(ei))>(M_PI/2.0))
            vi.Scale(-1.0);
          (*normedvectors)[i] += vi;
        }

        // write normed directional vector to stream
        if(fileoutput)
          fileleorientation<<fil<<"    "<<std::setprecision(12)<<dirvec(0)<<" "<<dirvec(1)<<" "<<dirvec(2)<<endl;
      }
    }
    if(fileoutput)
    {
      fprintf(fp, fileleorientation.str().c_str());
      fclose(fp);
    }
  }
}// StatMechManager:FilamentOrientations()

/*----------------------------------------------------------------------*
 | check the visualization mode of a crosslinker depending on the       |
 | filament(s) it is bound to.                    (public) mueller 07/10|
 *----------------------------------------------------------------------*/
bool STATMECH::StatMechManager::CheckForKinkedVisual(int eleid)
{
  bool kinked = true;
  // if element is a crosslinker
  if(eleid>basisnodes_)
  {
    DRT::Element *element = discret_->gElement(eleid);
    int lid = discret_->NodeColMap()->LID(element->NodeIds()[0]);
    for(int i=0; i<element->NumNode(); i++)
      if((*filamentnumber_)[lid]!=(*filamentnumber_)[element->NodeIds()[i]])
        kinked = false;
    return kinked;
  }
  else
    kinked = false;
  return kinked;
}//STATMECH::StatMechManager::CheckForKinkedVisual

