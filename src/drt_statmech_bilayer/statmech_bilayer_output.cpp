/*!----------------------------------------------------------------------
\file statmech_bilayer_output.cpp
\brief management and auxiliary functions for statistical mechanics of lipid
       bilayer

<pre>
Maintainer: Dhrubajyoti Mukherjee
            mukherjee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>

*----------------------------------------------------------------------*/
#include "statmech_manager_bilayer.H"

#include <Teuchos_Time.hpp>

#include "../drt_inpar/inpar_statmech.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_discsh3/discsh3.H"


/*--------------------------------------------------------------------------*
 | write special output for statistical mechanics (public)   mukherjee 09/15|
 *--------------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::Output(const int                     ndim,
                                       const double&                        time,
                                       const int&                           istep,
                                       const double&                        dt,
                                       const Epetra_Vector&                 dis,
                                       const Epetra_Vector&                 fint,
                                       bool                                 printscreen)
{
  /*in general simulations in statistical mechanics run over so many time steps that the amount of data stored in the error file
   * may exceed the capacity even of a server hard disk; thus, we rewind the error file in each time step so that the amount of data
   * does not increase after the first time step any longer*/
  Teuchos::ParameterList params = DRT::Problem::Instance()->StructuralDynamicParams();
//  int numstep = params.get<int>("NUMSTEP", -1);
  bool printerr = params.get<bool> ("print to err", false);
  FILE* errfile = params.get<FILE*> ("err file", NULL);
  if (printerr)
    rewind(errfile);

  //the following variable makes sense in case of serial computing only; its use is not allowed for parallel computing!
//  int num_dof = dis.GlobalLength();

  double starttime = statmechBilayerparams_.get<double>("STARTTIMEOUT",0.0);

    // handling gmsh output seperately
  if(DRT::INPUT::IntegralValue<int>(statmechBilayerparams_,"GMSHOUTPUT") && (time>=starttime && (istep-istart_) % statmechBilayerparams_.get<int> ("GMSHOUTINTERVALS", 100) == 0) )
  {
    /*construct unique filename for gmsh output with two indices: the first one marking the time step number
     * and the second one marking the newton iteration number, where numbers are written with zeros in the front
     * e.g. number one is written as 000001, number fourteen as 000014 and so on;*/

    // first index = time step index
    std::ostringstream filename;

    //creating complete file name dependent on step number with 6 digits and leading zeros
    if (istep<1000000)
      filename << outputrootpath_ << "/GmshOutput/bilayer"<< std::setw(6) << std::setfill('0') << istep <<".pos";
    else
      dserror("Gmsh output implemented for a maximum of 999999 steps");

    //calling method for writing Gmsh output
      GmshOutput(dis,filename,istep,time);
  }

  // Different tests on Memrane/vesicles
  switch (DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechBilayerparams_, "SPECIAL_OUTPUT"))
  {
    case INPAR::STATMECH::statout_vesceqshapes:
    {
      if(dis.Comm().MyPID()==0)
      {
        std::cout<<"----------------------------------------------------------------------"<<std::endl;
        std::cout<<"Reduced volume at the end of Equilibrium Shapes Simulation:"<<(1-time)<<std::endl;
        std::cout<<"---------------------------------------------------------------------"<<std::endl;
      }

      PrintTotalArea(dis);

      PrintTotalVol(dis);

    }
    break;
    default:
    break;
  }
  return;
} // STATMECH::StatMechManagerBilayer::Output()

/*--------------------------------------------------------------------------------------*
 | print total area of lipid bilayer (effect of area penalty) (public)   mukherjee 09/15|
 *--------------------------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::PrintTotalArea(const Epetra_Vector&  dis)
{
  //we need displacements also of ghost nodes and hence export displacment vector to column map format
  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(dis, discol);

  double area_ele=0;
  double area_total=0;
  for (int i=0; i < discret_->NumMyRowElements() - 1; i++)
  {
    //getting pointer to current element
    DRT::Element* element = discret_->lRowElement(i);
    const DRT::ElementType & eot = element->ElementType();

    if(eot == DRT::ELEMENTS::DiscSh3Type::Instance())    // discrete shell element
    {
      DRT::ELEMENTS::DiscSh3* ele = dynamic_cast<DRT::ELEMENTS::DiscSh3*>(element);
      double area=ele->CalcSurfArea(*discret_,discol,false);
      area_ele += area;
    }
  }
    discret_->Comm().SumAll(&area_ele, &area_total, 1);

  //proc 0 write complete output into file, all other proc inactive
  if(!discret_->Comm().MyPID())
  {

    FILE* fp = NULL; //file pointer for statistical output file

    //name of output file
    std::ostringstream outputfilename;
    outputfilename.str("");
    outputfilename << outputrootpath_<< "/StatMechOutput/SurfaceArea"<<outputfilenumber_+1 << ".dat";

    fp = fopen(outputfilename.str().c_str(), "a");
    std::stringstream filecontent;
    //          filecontent << iter_-1;
    filecontent << std::scientific << std::setprecision(10);
    filecontent << " " << area_total;

    filecontent << std::endl;
    fputs(filecontent.str().c_str(), fp);
    fclose(fp);
  }

  return;

}

/*---------------------------------------------------------------------------------------*
 | print enclosed vol of lipid bilayer (effect of vol penalty) (public)   mukherjee 09/15|
 *---------------------------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::PrintTotalVol(const Epetra_Vector&  dis)
{
  //we need displacements also of ghost nodes and hence export displacment vector to column map format
  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(dis, discol);

  double vol_ele=0;
  double vol_total=0;
  for (int i=0; i < discret_->NumMyRowElements() - 1; i++)
  {
    //getting pointer to current element
    DRT::Element* element = discret_->lRowElement(i);
    const DRT::ElementType & eot = element->ElementType();

    if(eot == DRT::ELEMENTS::DiscSh3Type::Instance())    // discrete shell element
    {
      DRT::ELEMENTS::DiscSh3* ele = dynamic_cast<DRT::ELEMENTS::DiscSh3*>(element);
      double vol=ele->CalcVolume(*discret_,discol,false);
      vol_ele +=vol;
    }
  }
    discret_->Comm().SumAll(&vol_ele, &vol_total, 1);

  //proc 0 write complete output into file, all other proc inactive
  if(!discret_->Comm().MyPID())
  {

    FILE* fp = NULL; //file pointer for statistical output file

    //name of output file
    std::ostringstream outputfilename;
    outputfilename.str("");
    outputfilename << outputrootpath_ << "/StatMechOutput/Volume"<<outputfilenumber_+1 << ".dat";

    fp = fopen(outputfilename.str().c_str(), "a");
    std::stringstream filecontent;
    //          filecontent << iter_-1;
    filecontent << std::scientific << std::setprecision(10);
    filecontent << " " << vol_total;

    filecontent << std::endl;
    fputs(filecontent.str().c_str(), fp);
    fclose(fp);
  }

  return;

}

/*----------------------------------------------------------------------------*
 | initialize special output for statistical mechanics(public) mukherjee 09/15|
 *----------------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::InitOutput(const int& ndim,
                                           const Epetra_Vector& dis,
                                           const int& istep,
                                           const double& dt)
{
//  dis.Print(std::cout);
  if(DRT::INPUT::IntegralValue<int>(statmechBilayerparams_,"GMSHOUTPUT"))
  {
    std::ostringstream filename;
    filename << outputrootpath_<<"/GmshOutput/bilayer000000.pos";
    GmshOutput(dis,filename,istep,0.0);
  }

  return;
} // STATMECH::StatMechManagerBilayer::InitOutput()


/*----------------------------------------------------------------------------*
 | Build output path for statistical mechanics(public)         mukherjee 09/15|
 *----------------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::BuildStatMechRootPath()
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
    if(stat(statmechfilepath.str().c_str(), &st) !=0 && DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechBilayerparams_, "SPECIAL_OUTPUT")!=INPAR::STATMECH::statout_none)
      dserror("The folder %s was not found but is required for statistical mechanics output!", statmechfilepath.str().c_str());
    std::ostringstream gmshfilepath;
    gmshfilepath << outputrootpath_ <<"/GmshOutput/";
    if(stat(gmshfilepath.str().c_str(), &st) !=0 && DRT::INPUT::IntegralValue<int>(statmechBilayerparams_, "GMSHOUTPUT"))
      dserror("The folder %s was not found but is required for Gmsh output!", gmshfilepath.str().c_str());
  return;
}


