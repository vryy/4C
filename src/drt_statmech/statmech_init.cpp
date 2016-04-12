/*!----------------------------------------------------------------------
\file statmech_init.cpp

\maintainer Kei Müller

\brief intialize statistical mechanics problems
*----------------------------------------------------------------------*/

#include "statmech_manager.H"

#include <Teuchos_Time.hpp>

#include "../drt_inpar/inpar_statmech.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../drt_beamcontact/beam3contact_manager.H"
#include "../drt_beamcontact/beam3contact_octtree.H"


/*----------------------------------------------------------------------*
 | Set Period Length and Search Resolution       mueller (public)  04/12|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::InitializeStatMechValues()
{
  if(!discret_->Comm().MyPID())
    std::cout<<"=================== StatMechManager Init ======================="<<std::endl;

  // determine network type / simulation type
  switch(DRT::INPUT::IntegralValue<INPAR::STATMECH::NetworkType>(statmechparams_, "NETWORKTYPE"))
  {
    case INPAR::STATMECH::networktype_std:
    {
      networktype_ = statmech_network_std;
      if(!discret_->Comm().MyPID())
        std::cout<<"-- standard random network simulation"<<std::endl;
    }
    break;
    case INPAR::STATMECH::networktype_casimir:
    {
      networktype_ = statmech_network_casimir;
      if(!discret_->Comm().MyPID())
        std::cout<<"-- casimir-force simulation setup..."<<std::endl;
    }
    break;
    case INPAR::STATMECH::networktype_motassay:
    {
      networktype_ = statmech_network_motassay;
      if(!discret_->Comm().MyPID())
        std::cout<<"-- motility assay simulation setup..."<<std::endl;
    }
    break;
    default:
      dserror("Unknown network type %i", DRT::INPUT::IntegralValue<INPAR::STATMECH::FilamentModel>(statmechparams_, "NETWORKTYPE"));
  }

  // determine filament model
  switch(DRT::INPUT::IntegralValue<INPAR::STATMECH::FilamentModel>(statmechparams_, "FILAMENTMODEL"))
  {
    case INPAR::STATMECH::filamentmodel_std:
    {
      filamentmodel_ = statmech_filament_std;
      if(!discret_->Comm().MyPID())
        std::cout<<"  -- standard filaments"<<std::endl;
    }
    break;
    case INPAR::STATMECH::filamentmodel_helical:
    {
      filamentmodel_ = statmech_filament_helical;
      if(!discret_->Comm().MyPID())
        std::cout<<"  -- filaments with helical binding site orientation..."<<std::endl;
    }
    break;
    default:
      dserror("Unknown filament model %i", DRT::INPUT::IntegralValue<INPAR::STATMECH::FilamentModel>(statmechparams_, "FILAMENTMODEL"));
  }

  //determine linker model
  switch(DRT::INPUT::IntegralValue<INPAR::STATMECH::LinkerModel>(statmechparams_, "LINKERMODEL"))
  {
    case INPAR::STATMECH::linkermodel_none:
      linkermodel_ = statmech_linker_none;
      if(!discret_->Comm().MyPID())
        std::cout<<"  -- no linkers"<<std::endl;
      break;
    case INPAR::STATMECH::linkermodel_std:
    {
      linkermodel_ = statmech_linker_std;
      if(statmechparams_.get<double>("ILINK",0.0) == 0.0 && statmechparams_.get<double>("ALINK",0.0) != 0.0)
      {
        if(!discret_->Comm().MyPID())
          std::cout<<"  -- standard Truss3 linkers"<<std::endl;
      }
      else if(statmechparams_.get<double>("ILINK",0.0) == 0.0 && statmechparams_.get<double>("ALINK",0.0) == 0.0)
      {
        if(!discret_->Comm().MyPID())
          std::cout<<"  -- standard Spring3 linkers"<<std::endl;
      }
      else
      {
        linkermodel_ = statmech_linker_std;
        if(!discret_->Comm().MyPID())
          std::cout<<"  -- standard Beam3 linkers"<<std::endl;
      }
    }
      break;
    case INPAR::STATMECH::linkermodel_stdintpol:
    {
      linkermodel_ = statmech_linker_stdintpol;
      if(statmechparams_.get<double>("ILINK",0.0) == 0.0)
      {
        if(!discret_->Comm().MyPID())
          std::cout<<"  -- standard Truss3 linkers with interpolated binding site positions for BEAM3 elements"<<std::endl;
      }
      else
      {
        linkermodel_ = statmech_linker_stdintpol;
        if(!discret_->Comm().MyPID())
          std::cout<<"  -- standard Beam3 linkers with interpolated binding site positions"<<std::endl;
      }
    }
      break;
    case INPAR::STATMECH::linkermodel_bellseq:
      linkermodel_ = statmech_linker_bellseq;
      if(!discret_->Comm().MyPID())
        std::cout<<"  -- Beam3 linkers accounting for Bell's equation"<<std::endl;
      break;
    case INPAR::STATMECH::linkermodel_bellseqintpol:
      linkermodel_ = statmech_linker_bellseqintpol;
      if(!discret_->Comm().MyPID())
        std::cout<<"  -- Beam3 linkers accounting for Bell's equation (interpolated binding site positions)"<<std::endl;
      break;
    case INPAR::STATMECH::linkermodel_active:
      linkermodel_ = statmech_linker_active;
      if(!discret_->Comm().MyPID())
        std::cout<<"  -- active Beam3 linkers"<<std::endl;
      break;
    case INPAR::STATMECH::linkermodel_activeintpol:
      linkermodel_ = statmech_linker_activeintpol;
      if(!discret_->Comm().MyPID())
        std::cout<<"  -- active Beam3 linkers (interpolated binding site positions)"<<std::endl;
      break;
    case INPAR::STATMECH::linkermodel_myosinthick:
      linkermodel_ = statmech_linker_myosinthick;
      if(!discret_->Comm().MyPID())
      {
        std::cout<<"  -- active Beam3 myosin thick filament (interpolated binding site positions)"<<std::endl;
        std::cout<<"  ===WARNING! Implementatation incomplete!!!==="<<std::endl;
      }
      break;
    default:
      dserror("Unknown linker model %i", DRT::INPUT::IntegralValue<INPAR::STATMECH::LinkerModel>(statmechparams_, "LINKERMODEL"));
  }

  // check linker fraction
  riseperbspot_ = Teuchos::rcp(new std::vector<double>());
  riseperbspot_->clear();
  {
    std::istringstream RPB(Teuchos::getNumericStringParameter(statmechparams_,"RISEPERBSPOT"));
    std::string word;
    char* input;
    while (RPB >> word)
      riseperbspot_->push_back(std::strtod(word.c_str(), &input));
  }
  if((int)riseperbspot_->size()>2)
    dserror("Currently only a maximum of two values is supported!");

  if(statmechparams_.get<double>("ACTIVELINKERFRACTION",0.0)>1.0 || statmechparams_.get<double>("ACTIVELINKERFRACTION",0.0)<0.0)
    dserror("The entered value %d for ACTIVELINKERFRACTION lies outside the valid interval [0;1]!");
  if(statmechparams_.get<double>("ACTIVELINKERFRACTION",0.0)>0.0 && linkermodel_!=statmech_linker_active && linkermodel_!=statmech_linker_activeintpol && linkermodel_!=statmech_linker_myosinthick)
    dserror("The parameter LINKERMODEL has to be set to active or activeintpol in order to work with ACTIVELINKERFRACTION>0.0!");
  if(statmechparams_.get<double>("ACTIVELINKERFRACTION",0.0)==0.0 && (linkermodel_==statmech_linker_active || linkermodel_==statmech_linker_activeintpol || linkermodel_==statmech_linker_myosinthick))
    dserror("Set ACTIVELINKERFRACTION to a value >0.0! Otherwise, the choice of LINKERMODEL is meaningless!");
  if((linkermodel_==statmech_linker_active || linkermodel_==statmech_linker_activeintpol || linkermodel_==statmech_linker_myosinthick) && statmechparams_.get<double>("DELTABELLSEQ",0.0)==0.0)
    dserror("Active linkers need a non-zero value DELTABELLSEQ for a force-dependent off-rate");
  if((linkermodel_==statmech_linker_std || linkermodel_==statmech_linker_stdintpol) && statmechparams_.get<double>("DELTABELLSEQ",0.0)!=0.0)
    dserror("This linker model requires DELTABELLSEQ==0.0! Check your input file!");
  if((linkermodel_==statmech_linker_stdintpol || linkermodel_ == statmech_linker_bellseqintpol || linkermodel_ == statmech_linker_activeintpol) && riseperbspot_->at(0)<=0.0)
    dserror("The input parameter RISEPERBSPOT has an invalid value %f", riseperbspot_->at(0));
  if(statmechparams_.get<double>("MAXRANDFORCE",-1.0)< 0.0 && statmechparams_.get<double>("MAXRANDFORCE",-1.0)!= -1.0)
    dserror("Choose a positive value for MAXRANDFORCE! Default value: -1.0 (no threshold for random forces!)");

  switch(DRT::INPUT::IntegralValue<INPAR::STATMECH::FrictionModel>(statmechparams_, "FRICTION_MODEL"))
  {
    case INPAR::STATMECH::frictionmodel_isotropiclumped:
      if(!discret_->Comm().MyPID())
        std::cout<<"  -- friction model: isotropic lumped"<<std::endl;
      break;
    case INPAR::STATMECH::frictionmodel_isotropicconsistent:
      if(!discret_->Comm().MyPID())
        std::cout<<"  -- friction model: isotropic consistent"<<std::endl;
      break;
    case INPAR::STATMECH::frictionmodel_anisotropicconsistent:
      if(!discret_->Comm().MyPID())
        std::cout<<"  -- friction model: anisotropic consistent"<<std::endl;
      break;
    default: dserror("No friction model (i.e. FRICTION_MODEL == none) declared in your input file!");
  }

  periodlength_ = Teuchos::rcp(new std::vector<double>);
  periodlength_->clear();
  {
    std::istringstream PL(Teuchos::getNumericStringParameter(statmechparams_,"PERIODLENGTH"));
    std::string word;
    char* input;
    while (PL >> word)
      periodlength_->push_back(std::strtod(word.c_str(), &input));
  }

  if((int)periodlength_->size()<3 && (int)periodlength_->size()!=1)
    dserror("You only gave %d values for PERIODLENGTH! Check your input file.", (int)periodlength_->size());
  else if((int)periodlength_->size()==1)
    for(int i=0; i<2; i++)
      periodlength_->push_back(periodlength_->at(0));
  for(int i=0; i<(int)periodlength_->size(); i++)
    if(periodlength_->at(i)<0.0)
      dserror("PERIODLENGTH( %d ) = %4.2f < 0.0 does not make any sense! Check your input file.", i, periodlength_->at(i));

  // set spatial resolution for search algorithm binding spots x crosslinkers
  if(statmechparams_.get<int>("SEARCHRES",1)<1)
    dserror("Please give a plausible value (>0) for SEARCHRES!");

  // determine search resolution in each spatial direction according to the periodlength_ vector
  searchres_ = Teuchos::rcp(new std::vector<int>(3,statmechparams_.get<int>("SEARCHRES",1)));
  // in case of a non-cubic periodic volume
  if(fabs(pow(periodlength_->at(0)*periodlength_->at(1)*periodlength_->at(2), 1.0/3.0)-periodlength_->at(0))>1e-4)
  {
    double Hmax = std::max(periodlength_->at(0), periodlength_->at(1));
    Hmax = std::max(Hmax, periodlength_->at(2));
    for(int i=0; i<(int)searchres_->size(); i++)
      searchres_->at(i) = (int)(floor((periodlength_->at(i)/Hmax) * (double)(statmechparams_.get<int>("SEARCHRES",1))));
  }

  if(!discret_->Comm().MyPID())
  {
    std::cout<<"search resolution: [x y z] = [";
    for(int i=0; i<(int)searchres_->size(); i++)
      std::cout<<searchres_->at(i)<<" ";
    std::cout<<"]"<<std::endl;
  }

  // read times for actions and the corresponding step sizes from input file
  actiontime_ = Teuchos::rcp(new std::vector<double>);
  actiontime_->clear();
  {
    std::istringstream TIME(Teuchos::getNumericStringParameter(statmechparams_,"ACTIONTIME"));
    std::string word;
    char* input;
    while (TIME >> word)
      actiontime_->push_back(std::strtod(word.c_str(), &input));
  }
  timestepsizes_ = Teuchos::rcp(new std::vector<double>);
  timestepsizes_->clear();
  {
    std::istringstream DT(Teuchos::getNumericStringParameter(statmechparams_,"ACTIONDT"));
    std::string word;
    char* input;
    while (DT >> word)
      timestepsizes_->push_back(std::strtod(word.c_str(), &input));
  }

  // sanity checks for time points and time step sizes
  Teuchos::ParameterList sdynparams = DRT::Problem::Instance()->StructuralDynamicParams();
  double iodt = sdynparams.get<double>("TIMESTEP", 1.0e9);
  double maxtime = sdynparams.get<double>("MAXTIME",-1.0);

  for(int i=0; i<(int)actiontime_->size()-1; i++)
    if(actiontime_->at(i)>actiontime_->at(i+1))
      dserror("ACTIONTIME values must be monotonously increasing!");
  if((int)actiontime_->size()!=(int)timestepsizes_->size())
    dserror("ACTIONTIME and ACTIONDT have to be equal in number in the input file!");
  if((int)actiontime_->size()<2)
  {
    if(timestepsizes_->at(0)>0.0) // unwanted stuff
      dserror("Give at least two values for ACTIONTIME AND ACTIONDT");
    else if(timestepsizes_->at(0)==0.0) //non-sensical stuff
      dserror("Given ACTIONDT is 0.0!");
    else // handled as default case from drt_validparameters.cpp and the rest of the non-sensical stuff
    {
      timestepsizes_->at(0) = iodt;
      timestepsizes_->push_back(iodt);
      actiontime_->at(0) = maxtime;
      actiontime_->push_back(maxtime);
    }
  }

  // set bctimeindex_ (position in actiontime_ where DBCs start being applied
  bctimeindex_ = statmechparams_.get<int>("BCTIMEINDEX", -1);

  //Element sanity checks
  ElementSanityCheck();

  // BC sanity checks
  BCSanityCheck();

  if(!discret_->Comm().MyPID())
  {
    std::cout<<"t_equilib  = t(0) = "<<std::setprecision(10)<<actiontime_->at(0)<<" @ dt = "<<timestepsizes_->at(0)<<std::endl;
    std::cout<<"t_ktswitch = t(1) = "<<std::setprecision(10)<<actiontime_->at(1)<<" @ dt = "<<timestepsizes_->at(1)<<std::endl;
    if(bctimeindex_>-1)
      std::cout<<"t_dbc     = t("<<bctimeindex_<<") ="<<std::setprecision(10)<<actiontime_->at(bctimeindex_)<<" @ dt = "<<timestepsizes_->at(bctimeindex_)<<std::endl;
    if(actiontime_->size()>3)
      std::cout<<"other action times..."<<std::endl;
    for(int i=0; i<(int)actiontime_->size(); i++)
      if(i>1 && i!=bctimeindex_)
        std::cout<<"t("<<i<<")        = "<<std::setprecision(10)<<actiontime_->at(i)<<" @ dt = "<<timestepsizes_->at(i)<<std::endl;
    std::cout<<"================================================================"<<std::endl;
  }
  return;
}

/*--------------------------------------------------------------------------------------------------*
 | Sanity checks for elementypes for filament and linkers                (private)  mukherjee 10/14 |
 *--------------------------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::ElementSanityCheck()
{
  // Not all Beam element type for filaments are compatible with all type of element types currently used
  // for Linkers. This method makes sure that the correct element type for filament and  linker is chosen.
  // An element used to browse through local Row Elements

  DRT::Element * element = discret_->gElement(discret_->lRowElement(0)->Id());
//  if (element->ElementType().Name()=="Beam3rType" && linkermodel_ != statmech_linker_none && statmechparams_.get<double>("ILINK",0.0) != 0.0)
//    dserror("Truss linkers are not currently configured to bind with filaments discretized with Beam3r elements.\nPlease choose Beam3Type element for linkers.");
  if (element->ElementType().Name()=="Beam3ebType" && linkermodel_ != statmech_linker_none && statmechparams_.get<double>("ILINK",0.0) != 0.0)
    dserror("Currently only Truss and spring linkers are configured to bind with filaments discretized with Beam3eb elements.\nPlease set ILINK=0 to activate truss linkers, ILINK=ALINK=0 activate spring linkers.");
  return;
}

/*----------------------------------------------------------------------*
 | Some simple sanity checks                   (private)  mueller 01/13 |
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::BCSanityCheck()
{
  // Get the time index signalling the start of DBC application ("-1", alternatively bctimeindex_--, due to the C++ counting convention)
  if(bctimeindex_<0) // default
    bctimeindex_ = (int)actiontime_->size()-1;
  else if(bctimeindex_==0)
    dserror("Given index BCTIMEINDEX = %i ! Start counting at 1!", bctimeindex_);
  else if(bctimeindex_>(int)actiontime_->size())
    dserror("Given index BCTIMEINDEX = %i lies outside the ACTIONTIME vector! Check your input file!", bctimeindex_);
  else
    bctimeindex_--;

  INPAR::STATMECH::DBCType dbctype = DRT::INPUT::IntegralValue<INPAR::STATMECH::DBCType>(statmechparams_,"DBCTYPE");
  switch(dbctype)
  {
    case INPAR::STATMECH::dbctype_std:
      break;
    case INPAR::STATMECH::dbctype_pinnodes:
      break;
    case INPAR::STATMECH::dbctype_none:
      break;
    // default: everything involving periodic boundary conditions and rheology
    default:
    {
      // "-1" because of counting convention
      int curvenumber = statmechparams_.get<int>("CURVENUMBER",0)-1;
      int displacementdir = statmechparams_.get<int>("DBCDISPDIR",-1)-1;

      if(periodlength_->at(0)<=0.0)
        dserror("PERIODLENGTH = %f! This method works only in conjunction with periodic boundary conditions!", periodlength_->at(0));
      if(periodlength_->at(0)!=periodlength_->at(1) || periodlength_->at(0)!=periodlength_->at(2) || periodlength_->at(1)!=periodlength_->at(2))
        dserror("Check PERIODLENGTH in your input file! This method only works for a cubic boundary volume");
      if(displacementdir>2 || displacementdir<0 || (curvenumber<0 && dbctype != INPAR::STATMECH::dbctype_movablesupport1d))
        dserror("In case of imposed DBC values, please define the StatMech parameters DBCDISPDIR={1,2,3} and/or CURVENUMBER correctly");
      if(dbctype==INPAR::STATMECH::dbctype_affineshear)
      {
        if((int)actiontime_->size()<3)
          dserror("For affine deformation, give three time values for ACTIONTIME! 1. equilibration time 2. dbc application time 3. release of DBC nodes other than nodes of elements shifted in z-direction!");
        else if(bctimeindex_!=1)
          dserror("BCTIMEINDEX must be set to 2 for affine shear displacement!");
      }
    }
  }
  INPAR::STATMECH::NBCType nbctype = DRT::INPUT::IntegralValue<INPAR::STATMECH::NBCType>(statmechparams_, "NBCTYPE");
  switch(nbctype)
  {
    case INPAR::STATMECH::nbctype_std:
      break;
    case INPAR::STATMECH::nbctype_constcreep:
    {
      if(dbctype != INPAR::STATMECH::dbctype_movablesupport1d)
        dserror("If you intend to run a creep simulation, set NBCTYPE to constcreep and DBCTYPE to movablesupport1d!");
      int displacementdir = statmechparams_.get<int>("DBCDISPDIR",-1)-1;
      if(displacementdir>2 || displacementdir<0)
        dserror("In case of of a Dirichlet-based movable support, please define DBCDISPDIR={1,2,3} as the direction which remains ->unconstrained<-!");
      int nbccurvenumber = statmechparams_.get<int>("NBCCURVENUMBER",0)-1;
      if(nbccurvenumber<0)
        dserror("Please give a Neumann BC curve number >0");
    }
    break;
    case INPAR::STATMECH::nbctype_randompointforce:
    {
      // for now: hard-coded number of (ten) consecutive nodes that are probed
      int N = statmechparams_.get<int>("NUMNBCNODES", 0);
      if(N<=0)
        dserror("Provide the number of Neumann nodes by means of the input parameter NUMNBCNODES (>0)!");
      if(N>discret_->NumMyColNodes())
        dserror("You have provided NUMNBCNODES greater than the number of nodes in the discretization!");
      if((int)actiontime_->size()-bctimeindex_<N)
        dserror("Starting from BCTIMEINDEX %i, only %i values in ACTIONTIME are given! Ten values are required!");
    }
    break;
    default:
      break;
  }
  return;
}

/*----------------------------------------------------------------------*
 | Set start step                                (private) mueller 08/13|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::SetStartStep(Teuchos::ParameterList& parameters)
{
  // set istart_, the number of steps after which we start writing output
  double starttimeout = statmechparams_.get<double>("STARTTIMEOUT", 0.0);
  // intermediate time intervals
  for(int i=0; i<(int)actiontime_->size(); i++)
  {
    // first time interval
    if(i==0)
    {
      if(starttimeout>=actiontime_->front())
        istart_ = (int)(round(actiontime_->front()/parameters.get<double>("delta time",0.01)));
      else
      {
        istart_ = (int)(round((starttimeout)/parameters.get<double>("delta time",0.01)));
        break;
      }
    }
    // last time interval
    else if(i==(int)actiontime_->size()-1)
    {
      if(starttimeout>=actiontime_->at(i))
        istart_ += (int)(round((actiontime_->back()-actiontime_->at(i-1))/timestepsizes_->at(i-1)) + (int)((starttimeout-actiontime_->back())/timestepsizes_->back()));
      else
        istart_ += (int)(round((starttimeout-actiontime_->at(i-1)) / timestepsizes_->at(i-1)));
    }
    // intermediate time intervals
    else
    {
      if(starttimeout>=actiontime_->at(i))
        istart_ += (int)(round((actiontime_->at(i)-actiontime_->at(i-1)) / timestepsizes_->at(i-1)));
      else
      {
        istart_ += (int)(round((starttimeout-actiontime_->at(i-1)) / timestepsizes_->at(i-1)));
        break;
      }
    }
  }
  return;
}

/*---------------------------------------------------------------------------------*
 | Initialize crosslinker positions                                                |
 |(point-like particles and not doubly-bound)                (public) mueller 07/10|
 *---------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::CrosslinkerMoleculeInit()
{
  int ncrosslink = statmechparams_.get<int> ("N_crosslink", 0);
  int numbins = statmechparams_.get<int>("HISTOGRAMBINS", 1);

  // create crosslinker maps
  std::vector<int> gids;
  for (int i=0; i<ncrosslink; i++)
    gids.push_back(i);
  // crosslinker column and row map
  crosslinkermap_ = Teuchos::rcp(new Epetra_Map(-1, ncrosslink, &gids[0], 0, discret_->Comm()));
  transfermap_    = Teuchos::rcp(new Epetra_Map(ncrosslink, 0, discret_->Comm()));
  startindex_ = Teuchos::rcp(new std::vector<double>);

  // create maps for binding spots in case of helical binding spot geometry of the actin filament
  if(filamentmodel_ == statmech_filament_helical)
  {
    // rise and rotation per binding spot (default values: 2.77nm and -166.15deg (left-handed, one-start helix))
    double rotperbspot = statmechparams_.get<double>("ROTPERBSPOT", -2.8999);

    //getting a vector consisting of pointers to all filament number conditions set
    std::vector<DRT::Condition*> filaments(0);
    discret_->GetCondition("FilamentNumber",filaments);

    // apply new beam element with intermediate binding spot positions
    if(linkermodel_ == statmech_linker_stdintpol  ||
        linkermodel_ == statmech_linker_activeintpol ||
        linkermodel_ == statmech_linker_bellseqintpol ||
        linkermodel_ == statmech_linker_myosinthick)
    {
      // temporary vectors
      // binding spot orientations relative to second triad vector
      std::vector<double> bspotorientations;
      // value of binding spot in parameter space
      std::vector<double> bspotxi;
      // maps binding spot to element
      std::vector<int>    bspot2element;
      // vector holding all bspot gids
      std::vector<int>    bspotgids;
      // vector signaling the existence of a binding spot on this proc
      std::vector<int>     bspotonproc;
      // CHECK vector
      std::vector<int>   filnodes;


      //loop over filaments to create redundant column map
      for (int i=0; i<(int)filaments.size(); i++)
      {
        double riseperbspot = riseperbspot_->at(0);
        if(riseperbspot_->size()>1 && statmechparams_.get<int>("NUMSUBSTRATEFIL",0)>0 && i<statmechparams_.get<int>("NUMSUBSTRATEFIL",0))
          riseperbspot = riseperbspot_->at(1);

        //get a pointer to nodal cloud covered by the current condition
        const std::vector<int>* nodeids = filaments[i]->Nodes();
        // retrieve filament length using first and last node of the current filament
        DRT::Node* node0 = discret_->lColNode(discret_->NodeColMap()->LID((*nodeids)[0]));
        double lfil = 0.0;
        std::vector<double> pos(6);
        // loop through all Nodes of the filaments, consider breakage by Periodic Boundary and calculate length of filament i
        for(int j=1;j<(int)nodeids->size();j++)
        {
          DRT::Node* tmpnode = discret_->lColNode(discret_->NodeColMap()->LID((*nodeids)[j]));
          for(int l=0;l<3;l++)
            pos[l+3]= tmpnode->X()[l];

          if(j==1)
          { for(int l=0;l<3;l++)
             pos[l]= node0->X()[l];
          }

          //shift position in case of breakage due to Periodic Boundary
          for(int k=0;k<3;k++)
          {
            if( fabs( pos[k+3] + periodlength_->at(k) - pos[k] ) < fabs( pos[k+3] - pos[k] ) )
            {
              pos[k+3] += periodlength_->at(k);
            }
            if( fabs( pos[k+3] - periodlength_->at(k) - pos[k] ) < fabs( pos[k+3] - pos[k] ) )
            {
              pos[k+3] -= periodlength_->at(k);
            }
          }

          lfil += sqrt((pos[3]-pos[0])*(pos[3]-pos[0])+
                       (pos[4]-pos[1])*(pos[4]-pos[1])+
                       (pos[5]-pos[2])*(pos[5]-pos[2]));

          for(int l=0;l<3;l++)
            pos[l] = pos[l+3];
        }
        // element length (assuming equal node spacing)
        double elelength = lfil/(double)(filaments[i]->Nodes()->size()-1);

        // add as many binding spots to the i-th filament as possible
        int bspot = 0;
        while((double)bspot*riseperbspot <= lfil)
        {
          // put binding spot id into the vector in order to create a fully overlapping binding spot map afterwards
          bspotgids.push_back((int)bspotgids.size());
          //std::cout<< " BSPOTGIDS "<< (int)bspotgids << std::endl;


          // calculate the orientation of the binding spot as well as the respective curve parameter
          /* determine element current element GID:
           * As set up by the Inputfile-Generator, the filament element GID can be retrieved by knowing the filament number and the lesser node number.
           * Attention: If that property is changed, one should of course refer to node->Elements(). In the meanwhile, this method seems faster*/
          //int elegid = node0->Id() - i + (int)(floor((double)bspot*riseperbspot/lfil));
          int ival;
          // we need to make sure that bspot belongs to existing Filaments. In case a bspot is located on a node we always choose Element with the lower GID
          ival = (int)(floor((double)bspot*riseperbspot/elelength));
          // There might be a more elegant way
          if((int)(ceil((double)bspot*riseperbspot/elelength))==ival && bspot != 0)
            ival-=1;
          int elegid = node0->Id() - i + ival;

          // add element GID if this element is a row map element on this processor
          if(discret_->ElementRowMap()->LID(elegid)>-1)
          {
            // mark binding spot as being on this proc
            bspotonproc.push_back(1);
            bspot2element.push_back(elegid);

            // orientation [0;2pi] and parameter xi_bs of the binding spot
            double phibs = 0.0;
            double xibs = -1;
            if(bspot==0)
            {
              bspotorientations.push_back(phibs);
              bspotxi.push_back(xibs);
            }
            else
            {
              phibs = (((double)(bspot)*rotperbspot)/(2.0*M_PI)-floor(((double)(bspot-1)*rotperbspot)/(2.0*M_PI)))*2.0*M_PI;
              ival = (int)floor((double)bspot*riseperbspot/elelength);
              if((int)ceil((double)bspot*riseperbspot/elelength)==ival)
                  ival-=1;
              xibs = 2*((double)bspot*riseperbspot/elelength - ival - 0.5);
              bspotorientations.push_back(phibs);
              bspotxi.push_back(xibs);
            }

          }
          else // only bspotonproc stores the information about off-proc elements since the structure can be easily rebuilt for the other vectors
            bspotonproc.push_back(0);
          // next binding spot on filament
          bspot++;
        }
      }

      // maps
      // create redundant binding spot column map based upon the bspotids
      bspotcolmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)bspotgids.size(), &bspotgids[0], 0, discret_->Comm()));
      // create processor-specific row maps
      std::vector<int> bspotrowgids;
      for(int i=0; i<(int)bspotonproc.size(); i++)
        if(bspotonproc[i]==1)
          bspotrowgids.push_back(i);  // note: since column map is fully overlapping: i=col. LID = GID
      bspotrowmap_ = Teuchos::rcp(new Epetra_Map((int)bspotgids.size(), (int)bspotrowgids.size(), &bspotrowgids[0], 0, discret_->Comm()));

      // vectors
      // initialize class vectors
      bspotstatus_ = Teuchos::rcp(new Epetra_Vector(*bspotcolmap_));
      bspotstatus_->PutScalar(-1.0);
      bspotorientations_ = Teuchos::rcp(new Epetra_Vector(*bspotcolmap_,true));
      bspotxi_ = Teuchos::rcp(new Epetra_Vector(*bspotcolmap_,true));
      bspot2element_ = Teuchos::rcp(new Epetra_Vector(*bspotrowmap_, true));
      DRT::Element *ele0 = discret_->gElement(discret_->ElementRowMap()->GID(0));
      bspot2nodes_= Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,ele0->NumNode(),true));

      Epetra_Vector bspotorientationsrow(*bspotrowmap_);
      Epetra_Vector bspotxirow(*bspotrowmap_);
      Epetra_MultiVector bspot2nodesrow(*bspotrowmap_,bspot2nodes_->NumVectors());
      for(int i=0; i<bspotrowmap_->NumMyElements(); i++)
      {
        bspotorientationsrow[i] = bspotorientations[i];
        bspotxirow[i] = bspotxi[i];
        (*bspot2element_)[i] = bspot2element[i];
        DRT::Element *ele = discret_->gElement(bspot2element[i]);
        for(int l=0;l<ele->NumNode();l++)
          bspot2nodesrow[l][i]=ele->NodeIds()[l];
      }
      // make information on orientations and curve parameters redundant on all procs
      Epetra_Import bspotimporter(*bspotcolmap_, *bspotrowmap_);
      bspotorientations_->Import(bspotorientationsrow,bspotimporter,Insert);
      bspotxi_->Import(bspotxirow,bspotimporter,Insert);
      bspot2nodes_->Import(bspot2nodesrow,bspotimporter,Insert);
    }
    else // beam3/beam3r element: new binding spot maps are equivalent to node maps
    {
      // maps
      bspotcolmap_ = Teuchos::rcp(new Epetra_Map(*(discret_->NodeColMap())));
      bspotrowmap_ = Teuchos::rcp(new Epetra_Map(*(discret_->NodeRowMap())));
      // initialize class vectors (we do not need nspotxi_ here since the nodes are the binding spots)
      bspotstatus_ = Teuchos::rcp(new Epetra_Vector(*bspotcolmap_));
      bspotstatus_->PutScalar(-1.0);
      bspotorientations_ = Teuchos::rcp(new Epetra_Vector(*bspotcolmap_));
      bspot2element_ = Teuchos::rcp(new Epetra_Vector(*bspotrowmap_));

      // orientations and vectors
      int bspotgid = 0;
      for(int i=0; i<(int)filaments.size(); i++)
      {
        //get a pointer to nodal cloud covered by the current condition
        const std::vector<int>* nodeids = filaments[i]->Nodes();
        // retrieve filament length using first and last node of the current filament
        DRT::Node* node0 = discret_->lColNode(discret_->NodeColMap()->LID((*nodeids)[0]));
        //DRT::Node* node1 = discret_->lColNode(discret_->NodeColMap()->LID((*nodeids)[((int)nodeids->size())-1]));
       double lfil= 0.0;
       std::vector<double> pos(6);
        // loop through all Nodes of the filaments, consider breakage by Periodic Boundary and calculate length of filament i
        for(int j=1;j<(int)nodeids->size();j++)
        {
          DRT::Node* tmpnode = discret_->lColNode(discret_->NodeColMap()->LID((*nodeids)[j]));
          for(int l=0;l<3;l++)
            pos[l+3]= tmpnode->X()[l];

          if(j==1)
          {
            for(int l=0;l<3;l++)
             pos[l]= node0->X()[l];
          }
          //shift position in case of breakage due to Periodic Boundary
          for(int k=0;k<3;k++)
          {
            if( fabs( pos[k+3] + periodlength_->at(k) - pos[k] ) < fabs( pos[k+3] - pos[k] ) )
            {
              pos[k+3] += periodlength_->at(k);
            }
            if( fabs( pos[k+3] - periodlength_->at(k) - pos[k] ) < fabs( pos[k+3] - pos[k] ) )
            {
              pos[k+3] -= periodlength_->at(k);
            }
          }

          lfil += sqrt(pow(pos[3]-pos[0],2)+pow(pos[4]-pos[1],2)+pow(pos[5]-pos[2],2));

          for(int l=0;l<3;l++)
           pos[l] = pos[l+3];
        }

        // element length (assuming equal node spacing)
        double elelength = lfil/(double)(filaments[i]->Nodes()->size()-1);

        for(int j=0; j<(int)filaments[i]->Nodes()->size(); j++)
        {
          double phibs = 0.0;
          // first node
          if(j==0)
            (*bspotorientations_)[bspotcolmap_->LID(bspotgid)] = phibs;
          else
          {
            // determine orientation (relative to second triad vector)
            phibs = (((double)(j)*(elelength/riseperbspot_->at(0))*rotperbspot)/(2.0*M_PI)-floor(((double)(j-1)*(elelength/riseperbspot_->at(0))*rotperbspot)/(2.0*M_PI)))*2.0*M_PI;
            (*bspotorientations_)[bspotcolmap_->LID(bspotgid)] = phibs;
          }
          // assumption:
          if(bspotrowmap_->LID(bspotgid)>-1)
            (*bspot2element_)[bspotrowmap_->LID(bspotgid)] = discret_->lRowNode(discret_->NodeRowMap()->LID(bspotgid))->Elements()[0]->Id();
          bspotgid++;
        }
      }
    }
  }
  else // without helical orientation
  {
    if(linkermodel_ == statmech_linker_stdintpol ||
        linkermodel_ == statmech_linker_activeintpol ||
        linkermodel_ == statmech_linker_bellseqintpol ||
        linkermodel_ == statmech_linker_myosinthick)
    {
      // rise and rotation per binding spot (default values: 2.77nm and -166.15deg (left-handed, one-start helix))
      //getting a vector consisting of pointers to all filament number conditions set
      std::vector<DRT::Condition*> filaments(0);
      discret_->GetCondition("FilamentNumber",filaments);

      // value of binding spot in parameter space
      std::vector<double> bspotxi;
      // maps binding spot to element
      std::vector<int>    bspot2element;
      // vector holding all bspot gids
      std::vector<int>    bspotgids;
      // vector signaling the existence of a binding spot on this proc
      std::vector<int>    bspotonproc;
      // CHECK vector
      std::vector<int>   filnodes;

      //loop over filaments to create redundant column map
      for (int i=0; i<(int)filaments.size(); i++)
      {
        double riseperbspot = riseperbspot_->at(0);
        if(riseperbspot_->size()>1 && statmechparams_.get<int>("NUMSUBSTRATEFIL",0)>0 && i<statmechparams_.get<int>("NUMSUBSTRATEFIL",0))
          riseperbspot = riseperbspot_->at(1);
        // filament length
        double lfil = 0;
        //get a pointer to nodal cloud covered by the current condition
        const std::vector<int>* nodeids = filaments[i]->Nodes();
        // retrieve filament length using first and last node of the current filament
        DRT::Node* node0 = discret_->lColNode(discret_->NodeColMap()->LID((*nodeids)[0]));
        //DRT::Node* node1 = discret_->lColNode(discret_->NodeColMap()->LID((*nodeids)[((int)nodeids->size())-1]));

        std::vector<double> pos(6);
        // loop through all Nodes of the filaments, consider breakage by Periodic Boundary and calculate length of filament i
        for(int j=1;j<(int)nodeids->size();j++)
        {
          DRT::Node* tmpnode = discret_->lColNode(discret_->NodeColMap()->LID((*nodeids)[j]));
          for(int l=0;l<3;l++)
            pos[l+3]= tmpnode->X()[l];

          if(j==1)
          { for(int l=0;l<3;l++)
             pos[l]= node0->X()[l];
          }
          //shift position in case of breakage due to Periodic Boundary
          for(int k=0;k<3;k++)
          {
            if( fabs( pos[k+3] + periodlength_->at(k) - pos[k] ) < fabs( pos[k+3] - pos[k] ) )
            {
              pos[k+3] += periodlength_->at(k);
            }
            if( fabs( pos[k+3] - periodlength_->at(k) - pos[k] ) < fabs( pos[k+3] - pos[k] ) )
            {
              pos[k+3] -= periodlength_->at(k);
            }
          }

          lfil += sqrt((pos[3]-pos[0])*(pos[3]-pos[0])+
                       (pos[4]-pos[1])*(pos[4]-pos[1])+
                      (pos[5]-pos[2])*(pos[5]-pos[2]));

          for(int l=0;l<3;l++)
            pos[l] = pos[l+3];
        }
        // element length (assuming equal node spacing)
        double elelength = lfil/(double)(filaments[i]->Nodes()->size()-1);

        // add as many binding spots to the i-th filament as possible
        int bspot = 0;
        //more elegant way for tolerance?
        double ltol = 1e-7;
        while((double)bspot*riseperbspot < lfil+ltol)
        {
          // put binding spot id into the vector in order to create a fully overlapping binding spot map afterwards
          bspotgids.push_back((int)bspotgids.size());
          // calculate the orientation of the binding spot as well as the respective curve parameter
          /* determine element current element GID:
           * As set up by the Inputfile-Generator, the filament element GID can be retrieved by knowing the filament number and the lesser node number.
           * Attention: If that property is changed, one should of course refer to node->Elements(). In the meanwhile, this method seems faster*/
          //int elegid = node0->Id() - i + (int)(floor((double)bspot*riseperbspot/lfil));
          int ival;
          // we need to make sure that bspot belongs to existing Filaments. In case a bspot is located on a node we always choose Element with the lower GID
          ival = (int)(floor((double)bspot*riseperbspot/elelength));
          // There might be a more elegant way
          if((int)(ceil((double)bspot*riseperbspot/elelength))==ival && bspot != 0)
            ival-=1;
          int elegid = node0->Id() - i + ival;

          // add element GID if this element is a row map element on this processor
          if(discret_->ElementRowMap()->LID(elegid)>-1)
          {
            // mark binding spot as being on this proc
            bspotonproc.push_back(1);
            bspot2element.push_back(elegid);


            // orientation [0;2pi] and parameter xi_bs of the binding spot
            double xibs = -1; // 0
            if(bspot==0)
            {
              bspotxi.push_back(xibs);
            }
            else
            {
              ival = (int)floor((double)bspot*riseperbspot/elelength);
              if((int)ceil((double)bspot*riseperbspot/elelength)==ival)
                  ival-=1;
              xibs = 2*((double)bspot*riseperbspot/elelength - ival - 0.5);
              bspotxi.push_back(xibs);
            }

          }
          else // only bspotonproc stores the information about off-proc elements since the structure can be easily rebuilt for the other vectors
            bspotonproc.push_back(0);
          // next binding spot on filament
          bspot++;
        }
      }
      // maps
      // create redundant binding spot column map based upon the bspotids
      bspotcolmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)bspotgids.size(), &bspotgids[0], 0, discret_->Comm()));
      // create processor-specific row maps
      std::vector<int> bspotrowgids;
      for(int i=0; i<(int)bspotonproc.size(); i++)
        if(bspotonproc[i]==1)
          bspotrowgids.push_back(i);  // note: since column map is fully overlapping: i=col. LID = GID

      bspotrowmap_ = Teuchos::rcp(new Epetra_Map((int)bspotgids.size(), (int)bspotrowgids.size(), &bspotrowgids[0], 0, discret_->Comm()));

      // vectors
      // initialize class vectors
      bspotstatus_ = Teuchos::rcp(new Epetra_Vector(*bspotcolmap_, true));
      bspotstatus_->PutScalar(-1.0);
      bspotxi_ = Teuchos::rcp(new Epetra_Vector(*bspotcolmap_,true));
      bspot2element_ = Teuchos::rcp(new Epetra_Vector(*bspotrowmap_, true));
      bspot2nodes_= Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,2,true));

      Epetra_Vector bspotxirow(*bspotrowmap_);
      Epetra_MultiVector bspot2nodesrow(*bspotrowmap_,2);
      for(int i=0; i<bspotrowmap_->NumMyElements(); i++)
      {
        bspotxirow[i] = bspotxi[i];
        (*bspot2element_)[i] = bspot2element[i];
        DRT::Element *ele = discret_->gElement(bspot2element[i]);
        for(int l=0;l<2;l++)
          bspot2nodesrow[l][i]=ele->NodeIds()[l];
      }
      // make information on orientations and curve parameters redundant on all procs
      Epetra_Import bspotimporter(*bspotcolmap_, *bspotrowmap_);
      bspotxi_->Import(bspotxirow,bspotimporter,Insert);
      bspot2nodes_->Import(bspot2nodesrow,bspotimporter,Insert);
    }
    else
    {
      bspotcolmap_ = Teuchos::rcp(new Epetra_Map(*(discret_->NodeColMap())));
      bspotrowmap_ = Teuchos::rcp(new Epetra_Map(*(discret_->NodeRowMap())));
      bspotstatus_ = Teuchos::rcp(new Epetra_Vector(*bspotcolmap_));
      bspotstatus_->PutScalar(-1.0);
    }
  }

  // create density-density-correlation-function map with
  if(DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_structanaly ||
     DRT::INPUT::IntegralValue<int>(statmechparams_, "GMSHOUTPUT"))
  {
    std::vector<int> bins;
    for(int i=0; i<discret_->Comm().NumProc()*numbins; i++)
      bins.push_back(i);
    ddcorrcolmap_ = Teuchos::rcp(new Epetra_Map(-1, discret_->Comm().NumProc()*numbins, &bins[0], 0, discret_->Comm()));
    // create processor-specific density-density-correlation-function map
    ddcorrrowmap_ = Teuchos::rcp(new Epetra_Map(discret_->Comm().NumProc()*numbins, 0, discret_->Comm()));
    // create new trafo matrix (for later use in DDCorr Function where we evaluate in layer directions), initialize with identity matrix
    trafo_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(3,3,true));
    for(int i=0; i<trafo_->M(); i++)
      (*trafo_)(i,i) = 1.0;
    for(int i=0; i<(int)cog_.M(); i++)
      cog_(i) = periodlength_->at(i)/2.0;

    // start indices for parallel handling of orientation correlation etc in NumLinkerSpotsAndOrientation()
    // calculation of start indices for each processor
    // number of overall independent combinations
    int numnodes = discret_->NodeColMap()->NumMyElements();
    int numcombinations = (numnodes*numnodes-numnodes)/2;
    // combinations on each processor
    int combinationsperproc = (int)floor((double)numcombinations/(double)discret_->Comm().NumProc());
    int remainder = numcombinations%combinationsperproc;

    // get starting index tuples for later use
    startindex_->assign(2*discret_->Comm().NumProc(), 0.0);

    for(int mypid=0; mypid<discret_->Comm().NumProc()-1; mypid++)
    {
      std::vector<int> start(2,0);
      bool continueloop = false;
      bool quitloop = false;
      int counter = 0;
      int appendix = 0;
      if(mypid==discret_->Comm().NumProc()-1)
        appendix = remainder;

      // loop over crosslinker pairs
      for(int i=0; i<discret_->NodeColMap()->NumMyElements(); i++)
      {
        for(int j=0; j<discret_->NodeColMap()->NumMyElements(); j++)
        {
          if(i==(*startindex_)[2*mypid] && j==(*startindex_)[2*mypid+1])
            continueloop = true;
          if(j>i && continueloop)
          {
            if(counter<combinationsperproc+appendix)
              counter++;
            else
            {
              // new start index j
              if(j==discret_->NodeColMap()->NumMyElements()-1)
                start[1] = 0;
              else
                start[1] = j;
              quitloop = true;
              break;
            }
          }
        }
        if(quitloop)
        {
          // new start index i
          if(start[1]==0)
            start[0] = i+1;
          else
            start[0] = i;
          // new start tuple
          (*startindex_)[2*(mypid+1)] = (double)(start[0]);
          (*startindex_)[2*(mypid+1)+1] = (double)(start[1]);
          break;
        }
      }
    }

    if(DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_structanaly)
    {
      // initialize unbinding times vector to "-1"
      crosslinkunbindingtimes_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_,2));
      crosslinkunbindingtimes_->PutScalar(-1.0);
    }
  }

  // in case the internal forces of the crosslinker affect the off-rate
  if(linkermodel_ == statmech_linker_bellseq ||
     linkermodel_ == statmech_linker_bellseqintpol ||
     linkermodel_ == statmech_linker_active ||
     linkermodel_ == statmech_linker_activeintpol ||
     linkermodel_ == statmech_linker_myosinthick)
  {
    // note: needs to be in the column map format because we change element-specific properties which need to be communicated to ghosted elements also...
    element2crosslink_ = Teuchos::rcp(new Epetra_Vector(*(discret_->ElementColMap())));
    element2crosslink_->PutScalar(-1.0);
    unbindingprobability_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_,2,true));
  }
  else
  {
    element2crosslink_ = Teuchos::null;
    unbindingprobability_ = Teuchos::null;
  }

  std::vector<double> upperbound = *periodlength_;
  // handling both cases: with and without periodic boundary conditions
  if (periodlength_->at(0) == 0.0)
    for(int i=0; i<(int)upperbound.size(); i++)
      upperbound.at(i) = statmechparams_.get<double> ("MaxRandValue", 0.0);

  crosslinkerpositions_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_, 3, true));

  Teuchos::RCP<Epetra_MultiVector> crosslinkerpositionstrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_,3,true));

  if(DRT::INPUT::IntegralValue<int>(statmechparams_, "PLANELINKERMOTION"))
  {
    for (int i=0; i<crosslinkerpositionstrans->MyLength(); i++)
    {
      for (int j=0; j<crosslinkerpositionstrans->NumVectors()-1; j++)
        (*crosslinkerpositionstrans)[j][i] = upperbound.at(j) * (*uniformgen_)();
      (*crosslinkerpositionstrans)[2][i] = upperbound.at(2)/2.0;
    }
  }
  else
  {
    for (int i=0; i<crosslinkerpositionstrans->MyLength(); i++)
      for (int j=0; j<crosslinkerpositionstrans->NumVectors(); j++)
        (*crosslinkerpositionstrans)[j][i] = upperbound.at(j) * (*uniformgen_)();
  }
  CommunicateMultiVector(crosslinkerpositionstrans, crosslinkerpositions_,false,true);

  // initial bonding status is set (no bonds)
  crosslinkerbond_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_, 2));
  crosslinkerbond_->PutScalar(-1.0);
  // initial bond counter is set (no bonds)
  numbond_ = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_, true));

  // crosslinker element IDs of the crosslink molecules
  crosslink2element_ = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_));
  crosslink2element_->PutScalar(-1.0);
  if(linkermodel_==statmech_linker_myosinthick)
  {
    additionalcross2ele_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_,2));
    additionalcross2ele_->PutScalar(-1.0);
  }

  // initialize the beautiful visuals vector (aka beevee-vector)
  visualizepositions_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkermap_, 3, false));
  for (int i=0; i<visualizepositions_->MyLength(); i++)
    for (int j=0; j<visualizepositions_->NumVectors(); j++)
      (*visualizepositions_)[j][i] = (*crosslinkerpositions_)[j][i];

  if(linkermodel_ == statmech_linker_active ||
      linkermodel_ == statmech_linker_activeintpol ||
      linkermodel_ == statmech_linker_myosinthick ||
      statmechparams_.get<double>("ACTIVELINKERFRACTION",0.0)>0.0)
  {
    // set management vector to store the actual length of an active crosslinker: 0=long, 1=short
    // initial set long (=0)
    crosslinkeractlength_ = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_, true));
    //set management vector to store the time since the last time the crosslinker was linked
    crosslinkeractcycletime_ = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_));
    crosslinkeractcycletime_->PutScalar(statmechparams_.get<double>("ACTIVELINKERCYCLE",0.04));

    // signals linker type
    crosslinkertype_ = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_, true));
    if(!discret_->Comm().MyPID())
    {
      std::vector<int> activeorder = Permutation(statmechparams_.get<int>("N_crosslink",0));

      for(int i=0; i<(int)(floor(statmechparams_.get<double>("ACTIVELINKERFRACTION",0.0)*statmechparams_.get<int>("N_crosslink",0))); i++)
        (*crosslinkertype_)[activeorder[i]] = 1.0;
    }
    Teuchos::RCP<Epetra_Vector> crosslinkertypetrans = Teuchos::rcp(new Epetra_Vector(*transfermap_,true));
    CommunicateVector(crosslinkertypetrans,crosslinkertype_);
  }

  // calculate  number of total binding spot depending on BSPOTINTERVAL
  numbspots_ = 0;
  int bspotinterval = statmechparams_.get<int>("BSPOTINTERVAL",1);
  if(bspotinterval>1)
  {
    for(int i=0; i<bspotstatus_->MyLength(); i++)
      if(i%bspotinterval==0)
        numbspots_++;
  }
  else if(bspotinterval==1)
    numbspots_ = bspotstatus_->MyLength();
  else
    dserror("Check your input file parameter BSPOTINTERVAL! It's %i !", bspotinterval);

#ifdef DEBUGCOUT
  std::cout<<"Proc "<<discret_->Comm().MyPID()<<": total number of binding spots = "<<numbspots_<<std::endl;
#endif

  return;
}//StatMechManager::CrosslinkerMoleculeInit

/*----------------------------------------------------------------------*
| Set crosslinkers wherever possible before the first time step         |
|                                                (private) mueller 11/11|
*----------------------------------------------------------------------*/
void STATMECH::StatMechManager::SetInitialCrosslinkers(Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager)
{
  if(linkermodel_ != statmech_linker_none)
  {
    if(beamcmanager!=Teuchos::null && (linkermodel_ == statmech_linker_stdintpol ||
                                       linkermodel_ == statmech_linker_activeintpol ||
                                       linkermodel_ == statmech_linker_bellseqintpol ||
                                       linkermodel_ == statmech_linker_myosinthick))
      dserror("Beam contact not implemented yet for the interpolated crosslinker element");

    const Epetra_Map noderowmap = *discret_->NodeRowMap();
    const Epetra_Map nodecolmap = *discret_->NodeColMap();
    // TODO
//==NEW
    // new node positions and rotations
    Teuchos::RCP<Epetra_MultiVector> bspotpositions = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));
    // In case of Truss C/L, bspotrotations indicate the tangent of nodal displacement vector
    Teuchos::RCP<Epetra_MultiVector> bspotrotations = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));;
    Epetra_Vector disrow(*discret_->DofRowMap(), true);
    Epetra_Vector discol(*discret_->DofColMap(), true);
    GetBindingSpotPositions(discol, bspotpositions, bspotrotations);

    Teuchos::RCP<Epetra_MultiVector> bspottriadscol = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,4,true));;
    GetBindingSpotTriads(bspotrotations, bspottriadscol);
//==ENDNEW
//==OLD
    //node positions and rotations (binding spot positions and rotations when applying 4-noded beam element)
    std::map<int, LINALG::Matrix<3, 1> > currentpositions;

    {
      std::map<int, LINALG::Matrix<3, 1> > currentrotations;
      // Vectors hold zero->ok, since no displacement at this stage of the simulation
      GetNodePositionsFromDisVec(discol, currentpositions, currentrotations, true);
      // do initial octree build
      if(beamcmanager!=Teuchos::null)
        beamcmanager->OcTree()->OctTreeSearch(currentpositions);
    }
//==ENDOLD

    // generate a random order of binding spots and crosslinkers (only on Proc 0)
    std::vector<int> randbspot;
    std::vector<int> randlink;

    // consider only every BSPOTINTERVAL-th binding spot for binding
    int bspotinterval = statmechparams_.get<int>("BSPOTINTERVAL",1);

    //1. set singly bound linkers
    int numsinglybound = 0;
    if(discret_->Comm().MyPID()==0)
    {
      randbspot = Permutation(bspotcolmap_->NumMyElements());
      randlink = Permutation(statmechparams_.get<int>("N_crosslink", 0));
      int numbspots = statmechparams_.get<int>("INITOCCUPIEDBSPOTS",0);

      if(numbspots != 0 && statmechparams_.get<int>("N_crosslink", 0)==0)
        dserror("You can't have finite number of initial binding spots and zero crosslinkers! Please check your input file.");

      if(linkermodel_ != statmech_linker_none && numbspots>bspotcolmap_->NumMyElements())
        dserror("Given number of initially occupied binding spots (%i) exceeds the total binding spot count (%i)! Check your input file!",numbspots, bspotcolmap_->NumMyElements());

      // attach linkers with one binding domain to filaments
      if(networktype_ == statmech_network_motassay)
      {
        // attach linkers to substrate filaments
        int subfil = statmechparams_.get<int>("NUMSUBSTRATEFIL",0);
        int crosslinknum=0;
        for(int i=0; i<numbspots; i++)
        {
          // depending on whether binding spots are real nodes or virtual binding sites along a filament, we have to map the binding spot LID to its corresponding node LID
          int nodelid = i;
          if(linkermodel_ == statmech_linker_stdintpol ||
              linkermodel_ == statmech_linker_activeintpol ||
              linkermodel_ == statmech_linker_bellseqintpol ||
              linkermodel_ == statmech_linker_myosinthick)
            nodelid = (*bspot2nodes_)[0][i];
          if (((int)(*filamentnumber_)[nodelid])<subfil) //vector starts with zero element
          {
            // if this binding spot is still unoccupied
            if((*bspotstatus_)[i]<-0.9)
            {
              // attach it to the first binding spot (i.e. update of relevant class vectors)
              (*bspotstatus_)[i] = crosslinknum;
              (*numbond_)[crosslinknum] = 1.0;
              for(int j=0; j<crosslinkerbond_->NumVectors(); j++)
                if((*crosslinkerbond_)[j][crosslinknum]<-0.1)
                {
                  (*crosslinkerbond_)[j][crosslinknum] = bspotcolmap_->GID(i);
                  numsinglybound++;
                  break;
                }

              //update crosslinker position
              Epetra_SerialDenseMatrix LID(1,1);
              LID(0,0) = i;
              CrosslinkerIntermediateUpdate(*bspotpositions, LID, crosslinknum);
              // update visualization
              for (int j=0; j<visualizepositions_->NumVectors(); j++)
                (*visualizepositions_)[j][crosslinknum] = (*crosslinkerpositions_)[j][crosslinknum];
              crosslinknum++;
            }
          }
          else
            break;
        }
      }
      else
      {
        // first, establish specified number of singly bound crosslinkers
        for(int i=0; i<numbspots; i++)
        {
          int firstbspot = randbspot[i];
          // if this binding spot is still unoccupied
          if((*bspotstatus_)[firstbspot]<-0.9 && bspotcolmap_->GID(firstbspot)%bspotinterval==0)
          {
            // get the ilink-th random crosslinker
            int currlink = randlink[i];
            // attach it to the first binding spot (i.e. update of relevant class vectors)
            (*bspotstatus_)[firstbspot] = currlink;
            (*numbond_)[currlink] = 1.0;
            for(int j=0; j<crosslinkerbond_->NumVectors(); j++)
              if((*crosslinkerbond_)[j][currlink]<0.1)
              {
                (*crosslinkerbond_)[j][currlink] = bspotcolmap_->GID(firstbspot);
                numsinglybound++;
                break;
              }

            //update crosslinker position
            Epetra_SerialDenseMatrix LID(1,1);
            LID(0,0) = firstbspot;
            CrosslinkerIntermediateUpdate(*bspotpositions, LID, currlink);
            // update visualization
            for (int j=0; j<visualizepositions_->NumVectors(); j++)
              (*visualizepositions_)[j][currlink] = (*crosslinkerpositions_)[j][currlink];
          }
        }
      }
    }

    // Communication, first stage (not sure if we need to communicate all of these vectors here (or later)
    // transfer vectors
    Teuchos::RCP<Epetra_Vector> bspotstatusrow = Teuchos::rcp(new Epetra_Vector(*discret_->NodeRowMap(),true));
    Teuchos::RCP<Epetra_Vector> numbondtrans = Teuchos::rcp(new Epetra_Vector(*transfermap_, true));
    Teuchos::RCP<Epetra_MultiVector> crosslinkerpositionstrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_,3,true));
    Teuchos::RCP<Epetra_MultiVector> visualizepositionstrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_,3,true));
    Teuchos::RCP<Epetra_MultiVector> crosslinkerbondtrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_,2,true));

    CommunicateVector(bspotstatusrow, bspotstatus_);
    CommunicateVector(numbondtrans, numbond_);
    CommunicateMultiVector(crosslinkerpositionstrans, crosslinkerpositions_);
    CommunicateMultiVector(visualizepositionstrans, visualizepositions_);
    CommunicateMultiVector(crosslinkerbondtrans, crosslinkerbond_);

    // this section creates crosslinker finite elements
    // Commented because of performance issues  when adding a large number of elements at once.
    // 2. Now, parallely search for neighbour nodes
    Teuchos::RCP<Epetra_MultiVector> neighbourslid = Teuchos::null;
    if(statmechparams_.get<int>("SEARCHRES",1)>0)
      PartitioningAndSearch(bspotpositions,bspottriadscol, neighbourslid);

    // 3. create double bonds
    // a vector indicating the crosslink molecule which is going to constitute a crosslinker element
    Teuchos::RCP<Epetra_Vector> addcrosselement = Teuchos::rcp(new Epetra_Vector(*crosslinkermap_, true));
    int numsetelements = 0;

    if(discret_->Comm().MyPID()==0 && networktype_ != statmech_network_motassay)
    {
      int numbspots = statmechparams_.get<int>("INITOCCUPIEDBSPOTS",0);
      for(int i=0; i<numbspots; i++)
      {
        // get current linker
        int currlink = randlink[i];
        // obtain a random order of neighboursLID indices
        std::vector<int> neighbourorder = Permutation(neighbourslid->NumVectors());
        // if there is a potential second binding spot in the vicinity of the crosslinker (sufficient to check the first entry)
        for(int j=0; j<neighbourslid->NumVectors(); j++)
        {
          int currneighbour = neighbourorder[j];
          int secondbspot = (int)(*neighbourslid)[currneighbour][currlink];

          // if second binding binding spot exists and spot is unoccupied
          if((*neighbourslid)[currneighbour][currlink]>-0.1 && (*bspotstatus_)[secondbspot] < -0.9)// && bspotcolmap_->GID(secondbspot)%bspotinterval==0)
          {
            int free = -1;
            Epetra_SerialDenseMatrix LID(2,1);
            for(int k=0; k<crosslinkerbond_->NumVectors(); k++)
            {
              if((*crosslinkerbond_)[k][currlink]<-0.9)
              {
                free = k;
                LID(k,0) = secondbspot;
              }
              else
                LID(k,0) = (*crosslinkerbond_)[k][currlink];
            }

            // deal with crosslinkers which are about to occupy two binding sites on the same filament
            if(statmechparams_.get<double>("K_ON_SELF",0.0)==0.0)
            {
              int nodelid0 = -1;
              int nodelid1 = -1;
              bool linkonsamefilament = false;
              switch(linkermodel_)
              {
                case statmech_linker_stdintpol:
                case statmech_linker_activeintpol:
                case statmech_linker_bellseqintpol:
                case statmech_linker_myosinthick:
                {
                  nodelid0 = discret_->NodeColMap()->LID((int)(*bspot2nodes_)[0][(int)LID(0,0)]);
                  nodelid1 = discret_->NodeColMap()->LID((int)(*bspot2nodes_)[0][(int)LID(1,0)]);
                  if((*filamentnumber_)[nodelid0] == (*filamentnumber_)[nodelid1])
                    linkonsamefilament = true;
                  break;
                }
                default:
                {
                  if((*filamentnumber_)[(int)LID(0,0)]==(*filamentnumber_)[(int)LID(1,0)])
                    linkonsamefilament = true;
                }
              }
              // leave neighourslid-loop
              if(linkonsamefilament)
                break;
            }

            // check for intersection in case of beam contact (CURRENTLY ONLY FOR CONVENTIONAL BEAM3 ELEMENT)
            bool intersection = false;
            if(beamcmanager!=Teuchos::null)
            {
              Epetra_SerialDenseMatrix nodecoords(3,2);
              for(int k=0; k<nodecoords.M(); k++)
              {
                nodecoords(k,0) = (*bspotpositions)[k][(int)LID(0,0)];
                nodecoords(k,1) = (*bspotpositions)[k][(int)LID(1,0)];
              }
              for(int k=0; k<(int)LID.M(); k++)
              {
                intersection = beamcmanager->OcTree()->IntersectBBoxesWith(nodecoords, LID);
                if(intersection)
                  break;
              }
            }

            // orientation check
            //unit direction vector between currently considered two nodes
            LINALG::Matrix<3,1> direction;
            for(int j=0;j<3;j++)
              direction(j) = (*bspotpositions)[j][(int)LID(0,0)]-(*bspotpositions)[j][(int)LID(1,0)];
            direction.Scale(1.0/direction.Norm2());

            if(CheckOrientation(direction,discol,bspottriadscol,LID) && !intersection)
            {
              numsetelements++;
              (*addcrosselement)[currlink] = 1.0;
              // establish double bond to the first given neighbour
              // attach it to the second binding spot
              (*bspotstatus_)[secondbspot] = currlink;
              (*crosslinkerbond_)[free][currlink] = bspotcolmap_->GID(secondbspot);
              (*numbond_)[currlink] = 2.0;

              if(DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_structanaly)
                for(int k=0; k<crosslinkunbindingtimes_->NumVectors(); k++)
                (*crosslinkunbindingtimes_)[k][currlink] = 0.0;

              //update crosslinker position
              CrosslinkerIntermediateUpdate(*bspotpositions, LID, currlink);
              // update visualization
              for (int k=0; k<visualizepositions_->NumVectors(); k++)
                (*visualizepositions_)[k][currlink] = (*crosslinkerpositions_)[k][currlink];
              break;
            }
          }
        }
      }
    }

    // Communication, second stage
    Teuchos::RCP<Epetra_Vector> addcrosselementtrans = Teuchos::rcp(new Epetra_Vector(*transfermap_,true));
    bspotstatusrow->PutScalar(0.0);
    crosslinkerpositionstrans->PutScalar(0.0);
    visualizepositionstrans->PutScalar(0.0);
    crosslinkerbondtrans->PutScalar(0.0);
    numbondtrans->PutScalar(0.0);

    // export and reimport
    CommunicateVector(addcrosselementtrans, addcrosselement);
    CommunicateVector(bspotstatusrow, bspotstatus_);
    CommunicateVector(numbondtrans, numbond_);
    CommunicateMultiVector(crosslinkerpositionstrans, crosslinkerpositions_);
    CommunicateMultiVector(visualizepositionstrans, visualizepositions_);
    CommunicateMultiVector(crosslinkerbondtrans, crosslinkerbond_);

    int numsetelementsall = -1;
    discret_->Comm().MaxAll(&numsetelements, &numsetelementsall, 1);

    if(numsetelementsall>0)
    {
      Teuchos::RCP<Epetra_MultiVector> nodalpositions = Teuchos::null;
      // rotations from displacement vector (not necessary for Truss3CL crosslinkers)
      Teuchos::RCP<Epetra_MultiVector> nodalrotations = Teuchos::null;
      // NODAL quaternions from filament elements
      Teuchos::RCP<Epetra_MultiVector> nodalquaternions = Teuchos::null;
      if((linkermodel_ == statmech_linker_stdintpol  ||
          linkermodel_ == statmech_linker_activeintpol ||
          linkermodel_ == statmech_linker_bellseqintpol ||
          linkermodel_ == statmech_linker_myosinthick) && statmechparams_.get<double>("ILINK",0.0)>0.0)
      {
        if(linkermodel_==statmech_linker_myosinthick)
          nodalpositions = Teuchos::rcp(new Epetra_MultiVector(*(discret_->NodeColMap()),3));
        nodalrotations = Teuchos::rcp(new Epetra_MultiVector(*(discret_->NodeColMap()),3));
        GetNodalBindingSpotPositionsFromDisVec(discol, Teuchos::null, nodalrotations);
        nodalquaternions = Teuchos::rcp(new Epetra_MultiVector(*(discret_->NodeColMap()),4));
        GetElementBindingSpotTriads(nodalquaternions);
      }
      else if(linkermodel_ == statmech_linker_std  && statmechparams_.get<double>("ILINK",0.0)==0.0 && statmechparams_.get<double>("ALINK",0.0)==0.0)
      {
        nodalquaternions = Teuchos::rcp(new Epetra_MultiVector(*(discret_->NodeColMap()),4));
        GetElementBindingSpotTriads(nodalquaternions);
        std::cout<<"nodalquaternions="<<*nodalquaternions<<std::endl;
      }

      // add elements to problem discretization (processor specific)
      for(int i=0; i<addcrosselement->MyLength(); i++)
      {
        if((*addcrosselement)[i]>0.9)
        {
          /*there is the problem of how to assign to a crosslinker element a GID which is certainly not used by any
           * other element; we know that for a network at the beginning (without crosslinkers) each node has a
           * connectivity of 1 or 2 so that the number of all elements used to discretize the filaments is smaller
           * than basisnodes_. Thus we choose crosslinker GIDs >= basisnodes_ for the additionally added crosslinker
           * elements. To make sure that two crosslinkers never have the same GID we add to basisnodes_ the value
           * basisnodes_*GID1, where GID1 is the GID of the first node of the crosslinker element. Then we add GID2
           * (note: GID2 < basisnodes_), which is the GID of the second node of the crosslinker element.
           * Hence basisnodes_ + GID1*basisnodes_ + GID2 always gives a GID which cannot be used by any other element;
           * the first node of the crosslinker element is always assumed to be the one with the greater GID; the owner
           * of the first node is assumed to be the owner of the crosslinker element*/

          // copied binding spot map (in order to handle the addition of elements based upon an invariant map just in case)
          Epetra_Map bspotcolmap(*bspotcolmap_);

          // obtain binding spot GID
          std::vector<int> bspotgid(2);
          // determine smaller and larger of the GIDs
          bspotgid.at(1) = std::min((int)(*crosslinkerbond_)[0][i],(int)(*crosslinkerbond_)[1][i]);
          bspotgid.at(0) = std::max((int)(*crosslinkerbond_)[0][i],(int)(*crosslinkerbond_)[1][i]);

          // different sizes due to different linker elements
          Teuchos::RCP<std::vector<int> > globalnodeids;
          if(linkermodel_ == statmech_linker_stdintpol ||
             linkermodel_ == statmech_linker_activeintpol ||
             linkermodel_ == statmech_linker_bellseqintpol ||
             linkermodel_ == statmech_linker_myosinthick)
          {
            globalnodeids = Teuchos::rcp(new std::vector<int>(4,0));
            for(int l=0;l<2;l++)
            {
              (*globalnodeids)[l]=(int)(*bspot2nodes_)[l][bspotgid[0]];
              (*globalnodeids)[l+2]=(int)(*bspot2nodes_)[l][bspotgid[1]];
            }
          }
          else
          {
            globalnodeids = Teuchos::rcp(new std::vector<int>(2,0));
            (*globalnodeids)[0] = bspotgid[0];
            (*globalnodeids)[1] = bspotgid[1];
          }

          // calculate element GID
          int newcrosslinkerGID = GenerateNewLinkerGID(bspotgid);
          /* Create mapping from crosslink molecule to crosslinker element GID
           * Note: No need for the usual procedure of exporting and reimporting to make things redundant
           * because info IS already redundant by design here.*/
          (*crosslink2element_)[i] = newcrosslinkerGID;

          //save positions of nodes between which a crosslinker has to be established in variables xrefe and rotrefe:
          std::vector<double> rotrefe(6);
          std::vector<double> xrefe(6);
          // resize in case of interpolated crosslinker element (not necessary for Truss3CL crosslinkers)
          if (linkermodel_ == statmech_linker_stdintpol ||
              linkermodel_ == statmech_linker_activeintpol ||
              linkermodel_ == statmech_linker_bellseqintpol ||
              linkermodel_ == statmech_linker_myosinthick)
            rotrefe.resize(12);

          for(int k=0; k<3; k++)
          {
            xrefe[k ] = (*bspotpositions)[k][bspotcolmap.LID(bspotgid.at(0))];
            xrefe[k+3] = (*bspotpositions)[k][bspotcolmap.LID(bspotgid.at(1))];

            //set nodal rotations (not true ones, only those given in the displacement vector)
            // In case of Truss C/L, bspotrotations indicate the tangent of nodal displacement vector
            if(statmechparams_.get<double>("ILINK",0.0)>0.0 || statmechparams_.get<double>("ILINK",0.0)==0.0)
            {
              if (linkermodel_ == statmech_linker_stdintpol  ||
                  linkermodel_ == statmech_linker_activeintpol ||
                  linkermodel_ == statmech_linker_bellseqintpol ||
                  linkermodel_ == statmech_linker_myosinthick)
              {
                rotrefe[k]   = (*nodalrotations)[k][bspotcolmap.LID((*globalnodeids)[0])];
                rotrefe[k+3] = (*nodalrotations)[k][bspotcolmap.LID((*globalnodeids)[1])];
                rotrefe[k+6] = (*nodalrotations)[k][bspotcolmap.LID((*globalnodeids)[2])];
                rotrefe[k+9] = (*nodalrotations)[k][bspotcolmap.LID((*globalnodeids)[3])];
              }
              else
              {
                rotrefe[k] = (*bspotrotations)[k][bspotcolmap.LID((*globalnodeids)[0])];
                rotrefe[k+3] = (*bspotrotations)[k][bspotcolmap.LID((*globalnodeids)[1])];
              }
            }
          }

          /*a crosslinker is added on each processor which is row node owner of at least one of its nodes;
           *the processor which is row map owner of the node with the larger GID will be the owner of the new crosslinker element; on the other processors the new
           *crosslinker element will be present as a ghost element only*/
          bool hasrownode = false;
          for(int k=0; k<(int)globalnodeids->size(); k++)
            if(noderowmap.LID(globalnodeids->at(k)) > -1)
            {
              hasrownode = true;
              break;
            }

          if(hasrownode)
            AddNewCrosslinkerElement(newcrosslinkerGID, globalnodeids,bspotgid, xrefe,rotrefe,*discret_,nodalquaternions);

          // add all new elements to contact discretization on all Procs
          if(beamcmanager!=Teuchos::null)
            AddNewCrosslinkerElement(newcrosslinkerGID, globalnodeids,bspotgid, xrefe,rotrefe,beamcmanager->BTSolDiscret(),nodalquaternions);
        }
      }
    }

    // synchronization for problem discretization
    discret_->CheckFilledGlobally();
    discret_->FillComplete(true, false, false);

    if(beamcmanager!=Teuchos::null)
    {
      beamcmanager->BTSolDiscret().CheckFilledGlobally();
      beamcmanager->BTSolDiscret().FillComplete(true, false, false);
    }

    //Gmsh output
    if(DRT::INPUT::IntegralValue<int>(statmechparams_,"GMSHOUTPUT"))
    {
      std::ostringstream filename;
      filename << StatMechRootPath() <<"/GmshOutput/InitLinks.pos";
      Epetra_Vector disrow(*discret_->DofRowMap(), true);
      GmshOutput(disrow,filename,0,0.0);
    }
    if(beamcmanager!=Teuchos::null && DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_octree)
      // "-2" for initial octree output
      beamcmanager->OcTree()->OctTreeSearch(currentpositions,-2);
    //std::couts
    if(!discret_->Comm().MyPID())
    {
      std::cout<<"====setting initial crosslinkers===="<<std::endl;
      std::cout<<"singly bound: "<<numsinglybound-numsetelements<<std::endl;
      std::cout<<"doubly bound: "<<numsetelements<<std::endl;
      std::cout<<"===================================="<<std::endl;
    }
  }
  return;
} // SetInitialCrosslinkers()

/*----------------------------------------------------------------------*
 | intialize filament numbers in ctor            (private) mueller 08/13|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::InitializeFilamentNumbers()
{
  filamentnumber_ = Teuchos::rcp( new Epetra_Vector(*(discret_->NodeColMap())) );
  filamentnumber_->PutScalar(-1.0);

  //getting a vector consisting of pointers to all filament number conditions set
  std::vector<DRT::Condition*> filamentnumberconditions(0);
  discret_->GetCondition("FilamentNumber",filamentnumberconditions);

  //next all the pointers to all the different conditions are looped
  for (int i=0; i<(int)filamentnumberconditions.size(); ++i)
  {
    //get filament number described by the current condition
    int filamentnumber = filamentnumberconditions[i]->GetInt("Filament Number");

    //get a pointer to nodal cloud covered by the current condition
    const std::vector<int>* nodeids = filamentnumberconditions[i]->Nodes();

    //loop through all the nodes of the nodal cloud
    for(int j=0; j<(int)nodeids->size() ; j++)
    {
      //get the node's global id
      int nodenumber = (*nodeids)[j];

      //turning global id into local one
      nodenumber = discret_->NodeColMap()->LID(nodenumber);

      //if the node does not belong to current processor Id is set by LID() to -1 and the node is ignored
      if(nodenumber > -1)
        (*filamentnumber_)[nodenumber] = filamentnumber;
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | initialize force sensors in ctor              (private) mueller 08/13|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::InitializeForceSensors()
{
  forcesensor_ = Teuchos::rcp( new Epetra_Vector(*(discret_->DofColMap())) );
  forcesensor_->PutScalar(-1.0);

  //Getting a vector consisting of pointers to all filament number conditions set
  std::vector<DRT::Condition*> forcesensorconditions(0);
  discret_->GetCondition("ForceSensor",forcesensorconditions);

  //next all the pointers to all the different conditions are looped
  for (int i=0; i<(int)forcesensorconditions.size(); ++i)
  {
    //get number of nodal dof with respect to which force is to be measured; note: numbering starts with zero
    int nodedofnumber = forcesensorconditions[i]->GetInt("DOF Number") ;

    //get a pointer to nodal cloud covered by the current condition
    const std::vector<int>* nodeids = forcesensorconditions[i]->Nodes();

    //loop through all the nodes of the nodal cloud
    for (int j=0; j<(int)nodeids->size(); j++)
    {
      //get the node's global id
      int nodenumber = (*nodeids)[j];

      //testing whether current nodedofnumber makes sense for current node
      if (nodedofnumber < 0 || nodedofnumber >= discret_->NumDof(discret_->gNode(nodenumber)))
        dserror("ForceSensor condition applied with improper local dof number");

        //global id of degree of freedom at which force is to be measured
        int dofnumber = discret_->Dof( discret_->gNode(nodenumber), nodedofnumber-1 );

        /*if the node does not belong to current processor Id is set by LID() to -1 and the node is ignored; otherwise the degrees of
         * freedom affected by this condition are marked in the vector *forcesensor_ by a one entry*/
        if(nodenumber > -1)
          (*forcesensor_)[dofnumber] = 1.0;
    }
  }
  return;
}
