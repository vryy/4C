/*!
\file str_timint_databiopolynetdyn.cpp

\brief Biopolymer dynamics data container for the structural (time)
       integration

\maintainer Jonas Eichinger

\date Jun 22, 2016

\level 3

*/
/*---------------------------------------------------------------------*/


#include "str_timint_databiopolynetdyn.H"

#include "str_model_evaluator_data.H"
#include "../drt_lib/drt_globalproblem.H"
#include "str_timint_base.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::TIMINT::DataSMDyn::DataSMDyn()
: isinit_(false),
  issetup_(false),
  smdynparams_(Teuchos::null),
  statmechprob_(false),
  dynloadbalanceevery_(0),
  dbctype_(INPAR::STATMECH::dbctype_none),
  nbctype_(INPAR::STATMECH::nbctype_std),
  networktype_(INPAR::STATMECH::networktype_std),
  searchres_(1),
  checkorient_(false),
  bindingsitesearch_(INPAR::STATMECH::bsstype_volpart),
  maxrandforce_(-1.0),
  gmshoutput_(false),
  gmshnumintpt_(10),
  gmshnetstruct_(false),
  plotfactorthick_(0.0),
  specialoutput_(INPAR::STATMECH::statout_none),
  gmshoutinterval_(100),
  outputinterval_(1),
  starttimeout_(0.0),
  actiontime_(Teuchos::null),
  actiondt_(Teuchos::null),
  bctimeindex_(-1),
  timeintconstrandnumb_(-1.0),
  periodlength_(Teuchos::null),
  histogrambins_(1),
  numrasterpoints_(3),
  shearamplitude_(0.0),
  eta_(0.0),
  kt_(0.0),
  ktact_(0.0),
  k_on_start_(0.0),
  k_on_end_(0.0),
  k_on_self_(0.0),
  k_off_start_(0.0),
  k_off_end_(0.0),
  dbcdispdir_(0),
  curvenumber_(0),
  nbccurvenumber_(0),
  nbcforceamp_(0.0),
  numnbcnodes_(0),
  filamentmodel_(INPAR::STATMECH::filamentmodel_std),
  filamentpolarity_(false),
  num_eval_elements_(-1),
  riseperbinspot_(Teuchos::null),
  rotperbinspot_(-2.8999),
  phibinspot_(0.524),
  binspotinterval_(1),
  linkermodel_(INPAR::STATMECH::linkermodel_none),
  planelinkermotion_(false),
  activelinkerfraction_(0.0),
  numcrosslink_(0),
  initoccupiedbinspots_(0),
  numsubstratefil_(0),
  r_link_(0.0),
  deltar_link_(0.0),
  alink_(0.0),
  ilink_(0.0),
  iplink_(0.0),
  elink_(0.0),
  linkerscalefactor_(0.0),
  strokedistance_(0.005),
  deltabellseq_(0.0),
  corient_(0.0),
  phizero_(0.0),
  phizerodev_(6.28),
  activelinkercycle_(0.04),
  activerecoveryfraction_(0.95),
  k_active_short_start_(0.0),
  k_active_short_end_(0.0),
  k_active_long_start_(0.0),
  k_active_long_end_(0.0),
  reducecrosslinksby_(0),
  crossbridgemodel_(false)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 |  Init statmech data                         (public)  eichinger 06/16|
 *----------------------------------------------------------------------*/
void STR::TIMINT::DataSMDyn::Init(
    const Teuchos::RCP<const STR::TIMINT::BaseDataSDyn>& data_sdyn_ptr)
{
  issetup_ = false;
  //---------------------------------------------------------------------------
  // get structure data container for global state and struct dyn params
  data_sdyn_ptr_ = data_sdyn_ptr;
  //---------------------------------------------------------------------------
  // get ID of actual processor in parallel
  //---------------------------------------------------------------------------
  const int& myrank = DRT::Problem::Instance()->GetDis("structure")->Comm().MyPID();
  //---------------------------------------------------------------------------
  // screen data output
  //---------------------------------------------------------------------------
  if(!myrank)
    std::cout<<"=================== StatMechData Init ======================="<<std::endl;
  //---------------------------------------------------------------------------
  /// parameter list of the statistical dynamics (read only)
  //---------------------------------------------------------------------------
  smdynparams_ = Teuchos::rcpFromRef<const Teuchos::ParameterList>
      (DRT::Problem::Instance()->StatisticalMechanicsParams());

  //---------------------------------------------------------------------------
  // flag for statistical mechanics problem
  //---------------------------------------------------------------------------
  statmechprob_ = DRT::INPUT::IntegralValue<int>(*smdynparams_, "STATMECHPROB");

  // set spatial resolution for search algorithm binding spots x crosslinkers
  dynloadbalanceevery_ = smdynparams_->get<int>("DYNLOADBALANCEEVERY");
  if(dynloadbalanceevery_ < 0)
    dserror("Please give a plausible value (>=0) for DYNLOADBALANCEEVERY!");

  //---------------------------------------------------------------------------
  // general statmech simulation flags
  //---------------------------------------------------------------------------

  dbctype_ =
      DRT::INPUT::IntegralValue<INPAR::STATMECH::DBCType>(*smdynparams_,"DBCTYPE");
  nbctype_ =
      DRT::INPUT::IntegralValue<INPAR::STATMECH::NBCType>(*smdynparams_,"NBCTYPE");

  switch(DRT::INPUT::IntegralValue<INPAR::STATMECH::NetworkType>(*smdynparams_, "NETWORKTYPE"))
  {
    case INPAR::STATMECH::networktype_std:
    {
      networktype_ = INPAR::STATMECH::networktype_std;
      if(!myrank)
        std::cout<<"-- standard random network simulation"<<std::endl;
    }
    break;
    case INPAR::STATMECH::networktype_casimir:
    {
      networktype_ = INPAR::STATMECH::networktype_casimir;
      if(!myrank)
        std::cout<<"-- casimir-force simulation setup..."<<std::endl;
    }
    break;
    case INPAR::STATMECH::networktype_motassay:
    {
      networktype_ = INPAR::STATMECH::networktype_motassay;
      if(!myrank)
        std::cout<<"-- motility assay simulation setup..."<<std::endl;
    }
    break;
    default:
      dserror("Unknown network type. Fix your Input file ...");
  }

  periodlength_ = Teuchos::rcp(new std::vector<double>);
  periodlength_->clear();
  {
    std::istringstream PL(Teuchos::getNumericStringParameter(*smdynparams_,"PERIODLENGTH"));
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
  searchres_ = smdynparams_->get<int>("SEARCHRES");
  if(searchres_ < 1)
    dserror("Please give a plausible value (>0) for SEARCHRES!");

  // determine search resolution in each spatial direction according to the periodlength_ vector
  searchresdir_ = Teuchos::rcp(new std::vector<int>(3,searchres_));
  // in case of a non-cubic periodic volume
  if(fabs(pow(periodlength_->at(0)*periodlength_->at(1)*periodlength_->at(2), 1.0/3.0)-periodlength_->at(0))>1e-4)
  {
    double Hmax = std::max(periodlength_->at(0), periodlength_->at(1));
    Hmax = std::max(Hmax, periodlength_->at(2));
    for(int i=0; i<(int)searchresdir_->size(); i++)
      searchresdir_->at(i) = (int)(floor((periodlength_->at(i)/Hmax) * (double)(searchres_)));
  }

  // write to screen
  if(!myrank)
  {
    std::cout<<"search resolution: [x y z] = [";
    for(int i=0; i<(int)searchresdir_->size(); i++)
      std::cout<<searchresdir_->at(i)<<" ";
    std::cout<<"]"<<std::endl;
  }

  checkorient_ =
      (DRT::INPUT::IntegralValue<int>(*smdynparams_, "CHECKORIENT") == 1);
  bindingsitesearch_ =
      DRT::INPUT::IntegralValue<INPAR::STATMECH::BSSearchType>(*smdynparams_,"BINDINGSITESEARCH");
  maxrandforce_ =
      smdynparams_->get<double>("MAXRANDFORCE");

  //---------------------------------------------------------------------------
  // output set up
  //---------------------------------------------------------------------------
  gmshoutput_ =
      (DRT::INPUT::IntegralValue<int>(*smdynparams_, "GMSHOUTPUT") == 1);
  gmshnumintpt_ =
      smdynparams_->get<int>("GMSHNUMINTPT");
  gmshnetstruct_ =
      (DRT::INPUT::IntegralValue<int>(*smdynparams_,"GMSHNETSTRUCT") == 1);
  plotfactorthick_=
      smdynparams_->get<double>("PlotFactorThick");
  specialoutput_ =
      DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(*smdynparams_, "SPECIALOUTPUT");
  gmshoutinterval_ =
      smdynparams_->get<int> ("GMSHOUTINTERVAL");
  outputinterval_ =
      smdynparams_->get<int> ("OUTPUTINTERVAL");

  //---------------------------------------------------------------------------
  // event triggers
  //---------------------------------------------------------------------------
  starttimeout_ =
      smdynparams_->get<double>("STARTTIMEOUT");

  // read times for actions and the corresponding step sizes from input file
  actiontime_ = Teuchos::rcp(new std::vector<double>);
  actiontime_->clear();
  {
    std::istringstream TIME(Teuchos::getNumericStringParameter(*smdynparams_,"ACTIONTIME"));
    std::string word;
    char* input;
    while (TIME >> word)
      actiontime_->push_back(std::strtod(word.c_str(), &input));
  }
  actiondt_ = Teuchos::rcp(new std::vector<double>);
  actiondt_->clear();
  {
    std::istringstream DT(Teuchos::getNumericStringParameter(*smdynparams_,"ACTIONDT"));
    std::string word;
    char* input;
    while (DT >> word)
      actiondt_->push_back(std::strtod(word.c_str(), &input));
  }

  double iodt = DRT::Problem::Instance()->StructuralDynamicParams().get<double>("TIMESTEP");
  double maxtime = DRT::Problem::Instance()->StructuralDynamicParams().get<double>("MAXTIME");

  if(actiontime_->at(0)==0.0 and not (actiondt_->at(0)==iodt))
    dserror("TIMESTEP and ACTIONDT(0) have to be equal for ACTIONTIME(0)=0.0. Fix your input file");

  for(int i=0; i<(int)actiontime_->size()-1; i++)
    if(actiontime_->at(i)>actiontime_->at(i+1))
      dserror("ACTIONTIME values must be monotonously increasing!");
  if((int)actiontime_->size()!=(int)actiondt_->size())
    dserror("ACTIONTIME and ACTIONDT have to be equal in number in the input file!");
  if((int)actiontime_->size()<2)
  {
    if(actiondt_->at(0)>0.0) // unwanted stuff
      dserror("Give at least two values for ACTIONTIME AND ACTIONDT");
    else if(actiondt_->at(0)==0.0) //non-sensical stuff
      dserror("Given ACTIONDT is 0.0!");
    else // handled as default case from drt_validparameters.cpp and the rest of the non-sensical stuff
    {
      actiondt_->at(0) = iodt;
      actiondt_->push_back(iodt);
      actiontime_->at(0) = maxtime;
      actiontime_->push_back(maxtime);
    }
  }

  // set bctimeindex_ (position in actiontime_ where DBCs start being applied
  bctimeindex_ = smdynparams_->get<int>("BCTIMEINDEX");

  // Get the time index signalling the start of DBC application ("-1", alternatively bctimeindex_--, due to the C++ counting convention)
  if(bctimeindex_<0) // default
    bctimeindex_ = (int)actiontime_->size()-1;
  else if(bctimeindex_==0)
    dserror("Given index BCTIMEINDEX = %i ! Start counting at 1!", bctimeindex_);
  else if(bctimeindex_>(int)actiontime_->size())
    dserror("Given index BCTIMEINDEX = %i lies outside the ACTIONTIME vector! Check your input file!", bctimeindex_);
  else
    bctimeindex_--;

  if(!myrank)
  {
    std::cout<<"t_equilib  = t(0) = "<<std::setprecision(10)<<actiontime_->at(0)<<" @ dt = "<<actiondt_->at(0)<<std::endl;
    std::cout<<"t_ktswitch = t(1) = "<<std::setprecision(10)<<actiontime_->at(1)<<" @ dt = "<<actiondt_->at(1)<<std::endl;
    if(bctimeindex_>-1)
      std::cout<<"t_dbc     = t("<<bctimeindex_<<") ="<<std::setprecision(10)<<actiontime_->at(bctimeindex_)<<" @ dt = "<<actiondt_->at(bctimeindex_)<<std::endl;
    if(actiontime_->size()>3)
      std::cout<<"other action times..."<<std::endl;
    for(int i=0; i<(int)actiontime_->size(); i++)
      if(i>1 && i!=bctimeindex_)
        std::cout<<"t("<<i<<")        = "<<std::setprecision(10)<<actiontime_->at(i)<<" @ dt = "<<actiondt_->at(i)<<std::endl;
  }

  timeintconstrandnumb_ =
      smdynparams_->get<double>("TIMEINTCONSTRANDNUMB");

  //---------------------------------------------------------------------------
  // simulation specifications
  //---------------------------------------------------------------------------
  histogrambins_ =
      smdynparams_->get<int>("HISTOGRAMBINS");
  numrasterpoints_ =
      smdynparams_->get<int>("NUMRASTERPOINTS");
  shearamplitude_ =
      smdynparams_->get<double>("SHEARAMPLITUDE");
  eta_ =
      smdynparams_->get<double> ("ETA");
  kt_ =
      smdynparams_->get<double> ("KT");
  ktact_ =
      smdynparams_->get<double>("KTACT");
  k_on_start_ =
      smdynparams_->get<double>("K_ON_start");
  k_on_end_ =
      smdynparams_->get<double>("K_ON_end");
  k_on_self_ =
      smdynparams_->get<double>("K_ON_SELF");
  k_off_start_ =
      smdynparams_->get<double> ("K_OFF_start");
  k_off_end_ =
      smdynparams_->get<double> ("K_OFF_end");
  dbcdispdir_ =
      smdynparams_->get<int>("DBCDISPDIR");
  curvenumber_ =
      smdynparams_->get<int>("CURVENUMBER");
  nbccurvenumber_ =
      smdynparams_->get<int>("NBCCURVENUMBER");
  nbcforceamp_ =
      smdynparams_->get<double>("NBCFORCEAMP");
  numnbcnodes_ =
      smdynparams_->get<int>("NUMNBCNODES");

  // nbc sanity check
  switch(nbctype_)
  {
    case INPAR::STATMECH::nbctype_std:
      break;
    case INPAR::STATMECH::nbctype_constcreep:
    {
      if(dbctype_ != INPAR::STATMECH::dbctype_movablesupport1d)
        dserror("If you intend to run a creep simulation, set NBCTYPE to constcreep and DBCTYPE to movablesupport1d!");
      int displacementdir = dbcdispdir_ - 1;
      if(displacementdir>2 || displacementdir<0)
        dserror("In case of of a Dirichlet-based movable support, please define DBCDISPDIR={1,2,3} as the direction which remains ->unconstrained<-!");
      int nbccurvenumber = nbccurvenumber_ - 1;
      if(nbccurvenumber<0)
        dserror("Please give a Neumann BC curve number >0");
    }
    break;
    case INPAR::STATMECH::nbctype_randompointforce:
    {
      // for now: hard-coded number of (ten) consecutive nodes that are probed
      int N = numnbcnodes_;
      if(N<=0)
        dserror("Provide the number of Neumann nodes by means of the input parameter NUMNBCNODES (>0)!");
      if(N> DRT::Problem::Instance()->GetDis("structure")->NumMyColNodes())
        dserror("You have provided NUMNBCNODES greater than the number of nodes in the discretization!");
      if((int)actiontime_->size() - bctimeindex_ < N)
        dserror("Starting from BCTIMEINDEX %i, only %i values in ACTIONTIME are given! Ten values are required!");
    }
    break;
    default:
      break;
  }

  //---------------------------------------------------------------------------
  // filament model
  //---------------------------------------------------------------------------
  switch(DRT::INPUT::IntegralValue<INPAR::STATMECH::FilamentModel>(*smdynparams_, "FILAMENTMODEL"))
  {
    case INPAR::STATMECH::filamentmodel_std:
    {
      filamentmodel_ = INPAR::STATMECH::filamentmodel_std;
      if(!myrank)
        std::cout<<"  -- standard filaments"<<std::endl;
      break;
    }
    case INPAR::STATMECH::filamentmodel_helical:
    {
      filamentmodel_ = INPAR::STATMECH::filamentmodel_helical;
      if(!myrank)
        std::cout<<"  -- filaments with helical binding site orientation..."<<std::endl;
      break;
    }
    default:
      dserror("Unknown filament model. Fix your Input file ...");
  }

  filamentpolarity_ =
      (DRT::INPUT::IntegralValue<int>(*smdynparams_,"FILAMENTPOLARITY") == 1);

  // todo: this is dumb and needs to go (only used in dbc management, read from elemap)
  num_eval_elements_ =
      smdynparams_->get<int>("NUM_EVAL_ELEMENTS");

  // check linker fraction
  riseperbinspot_ = Teuchos::rcp(new std::vector<double>());
  riseperbinspot_->clear();
  {
    std::istringstream RPB(Teuchos::getNumericStringParameter(*smdynparams_,"RISEPERBINSPOT"));
    std::string word;
    char* input;
    while (RPB >> word)
      riseperbinspot_->push_back(std::strtod(word.c_str(), &input));
  }
  if((int)riseperbinspot_->size()>2)
    dserror("Currently only a maximum of two values is supported!");

  rotperbinspot_ =
      smdynparams_->get<double>("ROTPERBINSPOT");
  phibinspot_ =
      smdynparams_->get<double>("PHIBINSPOT");
  binspotinterval_ =
      smdynparams_->get<int>("BINSPOTINTERVAL");

  //---------------------------------------------------------------------------
  // crosslinker model
  //---------------------------------------------------------------------------
  switch(DRT::INPUT::IntegralValue<INPAR::STATMECH::LinkerModel>(*smdynparams_, "LINKERMODEL"))
  {
    case INPAR::STATMECH::linkermodel_none:
    {
      linkermodel_ = INPAR::STATMECH::linkermodel_none;
      if(!myrank)
        std::cout<<"  -- no linkers"<<std::endl;
      break;
    }
    case INPAR::STATMECH::linkermodel_std:
    {
      linkermodel_ = INPAR::STATMECH::linkermodel_std;
      if(ilink_ == 0.0 && alink_ != 0.0)
      {
        if(!myrank)
          std::cout<<"  -- standard Truss3 linkers"<<std::endl;
      }
      else if(ilink_ == 0.0 && alink_ == 0.0)
      {
        if(!myrank)
          std::cout<<"  -- standard Spring3 linkers"<<std::endl;
      }
      else
      {
        if(!myrank)
          std::cout<<"  -- standard Beam3 linkers"<<std::endl;
      }
      break;
    }
    case INPAR::STATMECH::linkermodel_stdintpol:
    {
      linkermodel_ = INPAR::STATMECH::linkermodel_stdintpol;
      if(ilink_ == 0.0)
      {
        if(!myrank)
          std::cout<<"  -- standard Truss3 linkers with interpolated binding site positions for BEAM3 elements"<<std::endl;
      }
      else
      {
        if(!myrank)
          std::cout<<"  -- standard Beam3 linkers with interpolated binding site positions"<<std::endl;
      }
      break;
    }
    case INPAR::STATMECH::linkermodel_bellseq:
    {
      linkermodel_ = INPAR::STATMECH::linkermodel_bellseq;
      if(!myrank)
        std::cout<<"  -- Beam3 linkers accounting for Bell's equation"<<std::endl;
      break;
    }
    case INPAR::STATMECH::linkermodel_bellseqintpol:
    {
      linkermodel_ = INPAR::STATMECH::linkermodel_bellseqintpol;
      if(!myrank)
        std::cout<<"  -- Beam3 linkers accounting for Bell's equation (interpolated binding site positions)"<<std::endl;
      break;
    }
    case INPAR::STATMECH::linkermodel_active:
    {
      linkermodel_ = INPAR::STATMECH::linkermodel_active;
      if(!myrank)
        std::cout<<"  -- active Beam3 linkers"<<std::endl;
      break;
    }
    case INPAR::STATMECH::linkermodel_activeintpol:
    {
      linkermodel_ = INPAR::STATMECH::linkermodel_activeintpol;
      if(!myrank)
        std::cout<<"  -- active Beam3 linkers (interpolated binding site positions)"<<std::endl;
      break;
    }
    case INPAR::STATMECH::linkermodel_myosinthick:
    {
      linkermodel_ = INPAR::STATMECH::linkermodel_myosinthick;
      if(!myrank)
      {
        std::cout<<"  -- active Beam3 myosin thick filament (interpolated binding site positions)"<<std::endl;
        std::cout<<"  ===WARNING! Implementatation incomplete!!!==="<<std::endl;
      }
      break;
    }
    default:
      dserror("Unknown linker model. Fix your Input file ...");
  }

  planelinkermotion_ =
      (DRT::INPUT::IntegralValue<int>(*smdynparams_, "PLANELINKERMOTION") == 1);
  activelinkerfraction_ =
      smdynparams_->get<double>("ACTIVELINKERFRACTION");
  numcrosslink_ =
      smdynparams_->get<int> ("NUMCROSSLINK");

  if (numcrosslink_>0 and linkermodel_ == INPAR::STATMECH::linkermodel_none)
    dserror("You have Crosslinker in your Volume but not choosen a linkermodel. Fix your Input file ...");

  initoccupiedbinspots_ =
      smdynparams_->get<int>("INITOCCUPIEDBINSPOTS");
  numsubstratefil_ =
      smdynparams_->get<int>("NUMSUBSTRATEFIL");
  r_link_ =
      smdynparams_->get<double>("R_LINK");
  deltar_link_ =
      smdynparams_->get<double>("DeltaR_LINK");
  alink_ =
      smdynparams_->get<double>("ALINK");
  ilink_ =
      smdynparams_->get<double>("ILINK");
  iplink_ =
      smdynparams_->get<double>("IPLINK");
  linkerscalefactor_ =
      smdynparams_->get<double>("LINKERSCALEFACTOR");
  strokedistance_ =
      smdynparams_->get<double>("STROKEDISTANCE");
  deltabellseq_ =
      smdynparams_->get<double>("DELTABELLSEQ");
  corient_ =
      smdynparams_->get<double> ("CORIENT");
  phizero_ =
      smdynparams_->get<double> ("PHIZERO");
  phizerodev_ =
      smdynparams_->get<double>("PHIZERODEV");
  activelinkercycle_ =
      smdynparams_->get<double>("ACTIVELINKERCYCLE");
  activerecoveryfraction_ =
      smdynparams_->get<double>("ACTIVERECOVERYFRACTION");
  reducecrosslinksby_ =
      smdynparams_->get<int>("REDUCECROSSLINKSBY");
  crossbridgemodel_ =
      (DRT::INPUT::IntegralValue<int>(*smdynparams_,"CROSSBRIDGEMODEL") == 1);

  // some safety checks
  if(activelinkerfraction_ > 1.0 ||
     activelinkerfraction_ < 0.0)
    dserror("The entered value %d for ACTIVELINKERFRACTION lies outside the valid interval [0;1]!");

  if(activelinkerfraction_ > 0.0 &&
     linkermodel_!=INPAR::STATMECH::linkermodel_active &&
     linkermodel_!=INPAR::STATMECH::linkermodel_activeintpol &&
     linkermodel_!=INPAR::STATMECH::linkermodel_myosinthick)
    dserror("The parameter LINKERMODEL has to be set to active or activeintpol in order to work with ACTIVELINKERFRACTION>0.0!");

  if(activelinkerfraction_ == 0.0 &&
     (linkermodel_ == INPAR::STATMECH::linkermodel_active ||
      linkermodel_ == INPAR::STATMECH::linkermodel_activeintpol ||
      linkermodel_ == INPAR::STATMECH::linkermodel_myosinthick))
    dserror("Set ACTIVELINKERFRACTION to a value >0.0! Otherwise, the choice of LINKERMODEL is meaningless!");

  if((linkermodel_ == INPAR::STATMECH::linkermodel_active ||
      linkermodel_ == INPAR::STATMECH::linkermodel_activeintpol ||
      linkermodel_ == INPAR::STATMECH::linkermodel_myosinthick) &&
      deltabellseq_ == 0.0)
    dserror("Active linkers need a non-zero value DELTABELLSEQ for a force-dependent off-rate");

  if((linkermodel_ == INPAR::STATMECH::linkermodel_std ||
      linkermodel_ == INPAR::STATMECH::linkermodel_stdintpol) &&
      deltabellseq_!=0.0)
    dserror("This linker model requires DELTABELLSEQ==0.0! Check your input file!");

  if((linkermodel_ == INPAR::STATMECH::linkermodel_stdintpol ||
      linkermodel_ == INPAR::STATMECH::linkermodel_bellseqintpol ||
      linkermodel_ == INPAR::STATMECH::linkermodel_activeintpol) &&
      riseperbinspot_->at(0)<=0.0)
    dserror("The input parameter RISEPERBINSPOT has an invalid value %f", riseperbinspot_->at(0));

  if(maxrandforce_ < 0.0 &&
     maxrandforce_ != -1.0)
    dserror("Choose a positive value for MAXRANDFORCE! Default value: -1.0 (no threshold for random forces!)");

  if(!myrank)
    std::cout<<"================================================================"<<std::endl;

  // set flag
  isinit_ = true;

  return;

} // Init()

/*----------------------------------------------------------------------*
 |  Setup                                      (public)  eichinger 06/16|
 *----------------------------------------------------------------------*/
void STR::TIMINT::DataSMDyn::Setup()
{
  CheckInit();

  // set flag
  issetup_ = true;

  return;
} // Setup()

