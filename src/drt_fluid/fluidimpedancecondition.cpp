/*!----------------------------------------------------------------------
\file fluidimpedancecondition.cpp

\brief Method to deal with windkessel and other vascular BCs

\level 3

<pre>
\maintainer  Moritz Thon
             thon@mhpc.mw.tum.de
</pre>

*----------------------------------------------------------------------*/

#include <stdio.h>

#include "fluidimpedancecondition.H"
#include "fluid_utils.H"

#include "../linalg/linalg_ana.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_sparsematrixbase.H"

#include "../drt_fluid_ele/fluid_ele_action.H"



/*----------------------------------------------------------------------*
 |  Constructor (public)                                     chfoe 06/08|
 *----------------------------------------------------------------------*/
FLD::UTILS::FluidImpedanceWrapper::FluidImpedanceWrapper(const Teuchos::RCP<DRT::Discretization> actdis,
                                                         const double dta)
{
  std::vector<DRT::Condition*> impedancecond;
  actdis->GetCondition("ImpedanceCond",impedancecond);

  // the number of lines of impedance boundary conditions found in the input
  // note that some of these lines could belong to the same physical condition
  // which is then marked by the same 'ConditionID'
  // The corresponding resorting of result values has to be done later
  int numcondlines = impedancecond.size();


  if (numcondlines < 1)
    dserror("this function should just be called if there is a least one impedance condition.");

  // -------------------------------------------------------------------
  // get time period length of first condition, this should always be
  // the same!
  // -------------------------------------------------------------------
  double period = (impedancecond[0])->GetDouble("timeperiod");

  // -------------------------------------------------------------------
  // now care for the fact that there could be more than one input line
  // belonging to the same impedance boundary condition
  // -------------------------------------------------------------------
  for (int i=0; i<numcondlines; i++)
  {
    int condid = (impedancecond[i])->GetInt("ConditionID");

    double thisperiod = (impedancecond[i])->GetDouble("timeperiod");
    if (thisperiod != period)
      dserror("all periods of impedance conditions in one problem have to be the same!!!");

    // -------------------------------------------------------------------
    // test that we have an integer number of time steps per cycle
    // -------------------------------------------------------------------
    // something more intelligent could possibly be found one day ...
    double doublestepnum = period/dta;

    int cyclesteps = (int)(doublestepnum+0.5);
    double diff = doublestepnum - (double)cyclesteps;

    if ( abs(diff) > 1.0E-5 )
      dserror("Make sure that the cycle can be calculated within an integer number of steps!!!");


    // -------------------------------------------------------------------
    // allocate the impedance bc class members for every case
    // -------------------------------------------------------------------
    Teuchos::RCP<FluidImpedanceBc> impedancebc = Teuchos::rcp(new FluidImpedanceBc(actdis, dta, condid, impedancecond[i]) );

    // -----------------------------------------------------------------
    // sort impedance bc's in map and test, if one condition ID appears
    // more than once. Currently this case is forbidden.
    // -----------------------------------------------------------------
    bool inserted = impmap_.insert( std::make_pair( condid, impedancebc ) ).second;
    if ( !inserted )
      dserror("There are more than one impedance condition lines with the same ID. This can not yet be handled.");
  } // end loop over condition lines from input
  return;
} // end FluidImpedanceWrapper

/*----------------------------------------------------------------------*
 |  Split linearization matrix to a BlockSparseMatrixBase    Thon 07/15 |
 *----------------------------------------------------------------------*/

void FLD::UTILS::FluidImpedanceWrapper::UseBlockMatrix(Teuchos::RCP<std::set<int> >     condelements,
                                                       const LINALG::MultiMapExtractor& domainmaps,
                                                       const LINALG::MultiMapExtractor& rangemaps,
                                                       bool                             splitmatrix)
{
  std::map<const int, Teuchos::RCP<class FluidImpedanceBc> >::iterator mapiter;

  for (mapiter = impmap_.begin(); mapiter != impmap_.end(); mapiter++ )
  {
    mapiter->second->FluidImpedanceBc::UseBlockMatrix(condelements,domainmaps,rangemaps,splitmatrix);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Return relative vector of relative pressure errors      Thon 07/150 |
 *----------------------------------------------------------------------*/
std::vector<double> FLD::UTILS::FluidImpedanceWrapper::getWKrelerrors( )
{
  std::vector<double> wk_rel_error;

  // get an iterator to my map
    std::map<const int, Teuchos::RCP<class FluidImpedanceBc> >::iterator mapiter;

  for (mapiter = impmap_.begin(); mapiter != impmap_.end(); mapiter++ )
  {
    wk_rel_error.push_back( mapiter->second->FluidImpedanceBc::getWKrelerror() );
  }

 return wk_rel_error;
}

/*----------------------------------------------------------------------*
 |  Wrap for time step prepare of impedance conditions       Thon 07/15 |
 *----------------------------------------------------------------------*/
void FLD::UTILS::FluidImpedanceWrapper::PrepareTimeStepImpedances(const double time)
{
  std::map<const int, Teuchos::RCP<class FluidImpedanceBc> >::iterator mapiter;

  for (mapiter = impmap_.begin(); mapiter != impmap_.end(); mapiter++ )
  {
    mapiter->second->FluidImpedanceBc::PrepareTimeStepImpedance(time);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Wrap update of residual                                 chfoe 06/08 |
 *----------------------------------------------------------------------*/
void FLD::UTILS::FluidImpedanceWrapper::AddImpedanceBCToResidualAndSysmat(const double time, const double dta, Teuchos::RCP<Epetra_Vector>& residual, Teuchos::RCP<LINALG::SparseOperator>& sysmat )
{
  std::map<const int, Teuchos::RCP<class FluidImpedanceBc> >::iterator mapiter;

  for (mapiter = impmap_.begin(); mapiter != impmap_.end(); mapiter++ )
  {
    //calc flux
    mapiter->second->FluidImpedanceBc::FlowRateCalculation(time,mapiter->first);
    //calc pressure and traction vector
    mapiter->second->FluidImpedanceBc::OutflowBoundary(dta,mapiter->first);
    //add traction vector and linearisation to fluid residual and sysmat
    mapiter->second->FluidImpedanceBc::UpdateResidualAndSysmat(residual,sysmat);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Wrap for time update of impedance conditions             Thon 07/15 |
 *----------------------------------------------------------------------*/
void FLD::UTILS::FluidImpedanceWrapper::TimeUpdateImpedances(const double time, const double dta)
{
  std::map<const int, Teuchos::RCP<class FluidImpedanceBc> >::iterator mapiter;

  for (mapiter = impmap_.begin(); mapiter != impmap_.end(); mapiter++ )
  {
    //update time step
    mapiter->second->FluidImpedanceBc::TimeUpdateImpedance(time);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Wrap restart writing                                    chfoe 06/08 |
 *----------------------------------------------------------------------*/
void FLD::UTILS::FluidImpedanceWrapper::WriteRestart( IO::DiscretizationWriter& output )
{
  std::map<const int, Teuchos::RCP<class FluidImpedanceBc> >::iterator mapiter;

  for (mapiter = impmap_.begin(); mapiter != impmap_.end(); mapiter++ )
  {
    mapiter->second->FluidImpedanceBc::WriteRestart(output,mapiter->first);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Wrap restart reading                                    chfoe 06/08 |
 *----------------------------------------------------------------------*/
void FLD::UTILS::FluidImpedanceWrapper::ReadRestart( IO::DiscretizationReader& reader)
{
  std::map<const int, Teuchos::RCP<class FluidImpedanceBc> >::iterator mapiter;

  for (mapiter = impmap_.begin(); mapiter != impmap_.end(); mapiter++ )
    mapiter->second->FluidImpedanceBc::ReadRestart(reader,mapiter->first);

  return;
}



/*----------------------------------------------------------------------*
 |  Constructor (public)                                     chfoe 04/08|
 *----------------------------------------------------------------------*/
FLD::UTILS::FluidImpedanceBc::FluidImpedanceBc(const Teuchos::RCP<DRT::Discretization> actdis,
                                               const double dta,
                                               const int condid,
                                               DRT::Condition* impedancecond
):  discret_(actdis),
    myrank_(discret_->Comm().MyPID()),
    period_(impedancecond->GetDouble("timeperiod")),
    termradius_(impedancecond->GetDouble("termradius")),
    treetype_(*(impedancecond->Get<std::string>("tree"))),
    R1_(impedancecond->GetDouble("R1")),
    R2_(impedancecond->GetDouble("R2")),
    C_(impedancecond->GetDouble("C")),
    E_(impedancecond->GetDouble("stiffness")),
    H1_(impedancecond->GetDouble("H1")),
    H2_(impedancecond->GetDouble("H2")),
    H3_(impedancecond->GetDouble("H3")),
    P_np_(0.0),
    P_n_(0.0),
    Pc_np_(0.0),
    Pc_n_(0.0),
    Q_np_(0.0),
    Q_n_(0.0)
{
  // ---------------------------------------------------------------------
  // get time period length, steps per cycle and initialise flowratespos
  // ---------------------------------------------------------------------

  cyclesteps_   = (int)(period_/dta+0.5);
  flowratespos_ = 0;
  pressurespos_ = 0;
  dP_           = 100.0;
  P_0_           = 1.0;
  dta_          = dta;

  // ---------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // ---------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  impedancetbc_ = LINALG::CreateVector(*dofrowmap,true);
  impedancetbcsysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,108,false,true));
  //NOTE: do not call impedancetbcsysmat_->Complete() before it is filled, since
  // this is our check if it has already been initialized


  if(treetype_ == "windkessel_freq_indp")
  {
    // ---------------------------------------------------------------------
    // read in all possible impedance calibrating conditions
    // ---------------------------------------------------------------------
    std::vector<DRT::Condition*> impedance_calb_cond;
    discret_->GetCondition("ImpedanceCalbCond",impedance_calb_cond);

    for(unsigned int i = 0; i<impedance_calb_cond.size(); i++)
    {
      if(impedance_calb_cond[i]->GetInt("ConditionID") == condid)
      {
        P_n_  = impedance_calb_cond[i]->GetDouble("Pin_n");
        P_np_ = impedance_calb_cond[i]->GetDouble("Pin_np");
        Pc_n_   = impedance_calb_cond[i]->GetDouble("Pc_n");
        Pc_np_  = impedance_calb_cond[i]->GetDouble("Pc_np");
      }
    }

    Q_np_ = (P_np_ - Pc_np_)/R1_;
    // This par might look little messy but could be fixed in the future
    if (myrank_ == 0)
    {
      printf(" Pin initially is: %f  --  %f",P_n_,P_np_);
      printf("Frequency independent windkessel condition(%d) with:\n",condid);
      printf("          R1 = %f\n",R1_);
      printf("          R2 = %f\n",R2_);
      printf("          C  = %f\n",C_);
      printf("          Pc(initialt)= %f:\n",termradius_);
    }
  }
  else if ( treetype_ == "lung" or treetype_ == "artery")
  {
    flowrates_    = Teuchos::rcp(new std::vector<double>);
    flowrates_->push_back(0.0);

    // -------------------------------------------------------------------
    // determine area of actual outlet and get material data
    // -------------------------------------------------------------------
    double density=0.0, viscosity=0.0;
    double area = Area(density,viscosity,condid);

    // -------------------------------------------------------------------
    // calculate impedance values and fill vector 'impvalues_'
    // -------------------------------------------------------------------
    impvalues_.resize(cyclesteps_);
    Impedances(area,density,viscosity); //initial calculations happen in here
  }
  else if ( treetype_ == "windkessel")
  {
    flowrates_    = Teuchos::rcp(new std::vector<double>);
    flowrates_->push_back(0.0);

    // -------------------------------------------------------------------
    // determine area of actual outlet and get material data
    // -------------------------------------------------------------------
    double density=0.0, viscosity=0.0;
    double area = Area(density,viscosity,condid);

    // -------------------------------------------------------------------
    // calculate impedance values and fill vector 'impvalues_'
    // -------------------------------------------------------------------
    impvalues_.resize(cyclesteps_);
    Impedances(area,density,viscosity); //initial calculations happen in here
  }
  else
    dserror("Treetype %s not supported!",treetype_.c_str());

  // ---------------------------------------------------------------------
  // initialize the pressures vectors
  // ---------------------------------------------------------------------
  pressures_ = Teuchos::rcp(new  std::vector<double>);
  pressures_->push_back(0.0);


  return;
}


/*----------------------------------------------------------------------*
 |  Restart writing                                         chfoe 05/08 |
 *----------------------------------------------------------------------*/
void FLD::UTILS::FluidImpedanceBc::WriteRestart( IO::DiscretizationWriter&  output, const int condnum )
{
  // condnum contains the number of the present condition
  // condition Id numbers must not change at restart!!!!

  std::stringstream stream1, stream4, stream5, dpstream;


  if(treetype_ == "windkessel_freq_indp")
  {
    stream1 << "P_n" << condnum;
    // write the input pressure at time step n
    output.WriteDouble(stream1.str(), P_n_);
    stream1 << "Pc_n" << condnum;
    // write the capacitor pressure at time step n
    output.WriteDouble(stream1.str(), Pc_n_);
    stream1 << "P_np" << condnum;
    // write the input pressure at time step n
    output.WriteDouble(stream1.str(), P_np_);
    stream1 << "Pc_np" << condnum;
    // write the capacitor pressure at time step n
    output.WriteDouble(stream1.str(), Pc_np_);

    dpstream<<"dP"<<condnum;
    output.WriteDouble(dpstream.str(), dP_);
  }
  else if ( treetype_ == "windkessel")
  {
    stream1 << "P_n" << condnum;
    // write the input pressure at time step n
    output.WriteDouble(stream1.str(), P_n_);
    stream1 << "Q_n" << condnum;
    // write the flux pressure at time step n
    output.WriteDouble(stream1.str(), Q_n_);
    stream1 << "P_np" << condnum;
    // write the input pressure at time step n
    output.WriteDouble(stream1.str(), P_np_);
    stream1 << "Q_np" << condnum;
    // write the flux pressure at time step n
    output.WriteDouble(stream1.str(), Q_np_);

    dpstream<<"dP"<<condnum;
    output.WriteDouble(dpstream.str(), dP_);

    dpstream <<"P_0"<<condnum;
    output.WriteDouble(dpstream.str(), P_0_);

    stream1 << "flowratesId" << condnum;
    // write the flowrates of the previous period
    output.WriteRedundantDoubleVector(stream1.str(),flowrates_);

    if ( flowratespos_ != pressurespos_)
      dserror("Positions do not match!");

    // write the pressures
    stream4<< "pressuresId"<<condnum;
    output.WriteRedundantDoubleVector(stream4.str(), pressures_);

    // also write pressuresposition of this outlet
    stream5 << "pressuresposId" << condnum;
    output.WriteInt(stream5.str(), pressurespos_);
  }
  else
  {
    stream1 << "flowratesId" << condnum;
    // write the flowrates of the previous period
    output.WriteRedundantDoubleVector(stream1.str(),flowrates_);

    if ( flowratespos_ != pressurespos_)
      dserror("Positions do not match!");

    // write the pressures
    stream4<< "pressuresId"<<condnum;
    output.WriteRedundantDoubleVector(stream4.str(), pressures_);

    // also write pressuresposition of this outlet
    stream5 << "pressuresposId" << condnum;
    output.WriteInt(stream5.str(), pressurespos_);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Restart reading                                         chfoe 05/08 |
 *----------------------------------------------------------------------*/
void FLD::UTILS::FluidImpedanceBc::ReadRestart( IO::DiscretizationReader& reader, const int condnum  )
{
  // condnum contains the number of the present condition
  // condition Id numbers must not change at restart!!!!
  std::stringstream stream1, stream2, stream3, stream4, stream5, dpstream;

  // old time step size
  double odta = reader.ReadDouble("timestep");

  // get time step of the current problems
  double ndta = dta_;


  // time of restart
  double time = reader.ReadDouble("time");

  // -------------------------------------------------------------------
  // Read in the pressure values and the pressure difference
  // -------------------------------------------------------------------
  stream4 << "pressuresId"<<condnum;
  stream5 << "pressuresposId" << condnum;

  // read in pressures
  reader.ReadRedundantDoubleVector(pressures_ ,stream4.str());

  // read in the pressures' position
  pressurespos_ = reader.ReadInt(stream5.str());
  flowratespos_ = pressurespos_;

  // Get old Pressure Vector size
  int oPSize = (int)pressures_->size();

  // Calculate new Pressure Vector size
  int nPSize = (int)(double(oPSize)*odta/ndta);


  // evaluate the new pressure vector
  int np_pos = 0;
  Teuchos::RCP<std::vector<double> > np = Teuchos::rcp(new std::vector<double>(nPSize,0.0));
  this->interpolate(pressures_,np,pressurespos_,np_pos,time);

  // store new values in class
  pressurespos_ = np_pos;
  pressures_    = np;

  if(treetype_ == "windkessel_freq_indp")
  {
    stream1 << "Pin_n" << condnum;
    // read the input pressure at time step n
    P_n_ = reader.ReadDouble(stream1.str());
    stream1 << "Pc_n" << condnum;
    // read the capacitor pressure at time step n
    Pc_n_  = reader.ReadDouble(stream1.str());
    stream1 << "Pin_np" << condnum;
    // read the input pressure at time step n
    P_np_ = reader.ReadDouble(stream1.str());
    stream1 << "Pc_np" << condnum;
    // read the capacitor pressure at time step n
    Pc_np_  = reader.ReadDouble(stream1.str());
    // read in pressure difference
    dpstream <<"dP"<<condnum;
    dP_ = reader.ReadDouble(dpstream.str());

    // new time step size
    //    double ndta = dta_;
    if (dta_ == odta)
    {
      if (!myrank_)
      printf("Impedance restart old and new time step are the same - life is good\n");
      return;
    }
    // Get the new Pc_n and Pin_n through linear interpolation mapping
    double t = odta - dta_;
    Pc_n_  = (odta - t)/odta*Pc_n_  + t/odta*Pc_np_;
    P_n_ = (odta - t)/odta*P_n_ + t/odta*P_np_;
  }
  else if(treetype_ == "windkessel")
  {
    stream1 << "P_n" << condnum;
    // read the input pressure at time step n
    P_n_ = reader.ReadDouble(stream1.str());
    stream1 << "Q_n" << condnum;
    // read the flux pressure at time step n
    Q_n_ = reader.ReadDouble(stream1.str());
    stream1 << "P_np" << condnum;
    // read the input pressure at time step n
    P_np_ = reader.ReadDouble(stream1.str());
    stream1 << "Q_np" << condnum;
    // read the flux pressure at time step n
    Q_np_ = reader.ReadDouble(stream1.str());

    //get pressure difference and pressure of last period
    if ( abs(dP_ - 100.0) < 1e-14) //if we just initialized the class (in context of AC-FS3I this is not guaranteed!)
    {
      // read in pressure difference
      dpstream <<"dP"<<condnum;
      dP_ = reader.ReadDouble(dpstream.str());

      // read in pressure difference
      dpstream <<"P_0"<<condnum;
      P_0_ = reader.ReadDouble(dpstream.str());
    }
    else
    {
      dP_ = 0.0;
      P_0_ = 1.0;
    }

    stream1 << "flowratesId" << condnum;
    // read in flow rates
    reader.ReadRedundantDoubleVector(flowrates_ ,stream1.str());

    // Get old flowrates Vector size
    int oQSize = (int)flowrates_->size();

    // Calculate new flowrates Vector size
    int nQSize = (int)(double(oQSize)*odta/ndta);


    // check if vector of flowrates is not empty
    if (flowrates_->size() == 0)
    dserror("could not re-read vector of flowrates");

    // evaluate the new flow rate vector
    int nfr_pos = 0;
    Teuchos::RCP<std::vector<double> > nfr = Teuchos::rcp(new std::vector<double>(nQSize,0.0));
    this->interpolate(flowrates_,nfr,flowratespos_,nfr_pos,time);
    // store new values in class
    flowratespos_ = nfr_pos;
    flowrates_    = nfr;

    // finally, recompute the outflow boundary condition from last step
    // this way the vector need not to be stored
    OutflowBoundary(ndta,condnum);
  }
  else
  {
    stream1 << "flowratesId" << condnum;
    // read in flow rates
    reader.ReadRedundantDoubleVector(flowrates_ ,stream1.str());

    // Get old flowrates Vector size
    int oQSize = (int)flowrates_->size();

    // Calculate new flowrates Vector size
    int nQSize = (int)(double(oQSize)*odta/ndta);


    // check if vector of flowrates is not empty
    if (flowrates_->size() == 0)
    dserror("could not re-read vector of flowrates");

    // evaluate the new flow rate vector
    int nfr_pos = 0;
    Teuchos::RCP<std::vector<double> > nfr = Teuchos::rcp(new std::vector<double>(nQSize,0.0));
    this->interpolate(flowrates_,nfr,flowratespos_,nfr_pos,time);
    // store new values in class
    flowratespos_ = nfr_pos;
    flowrates_    = nfr;

    // finally, recompute the outflow boundary condition from last step
    // this way the vector need not to be stored
    OutflowBoundary(ndta,condnum);
  }

  return;
}


/*----------------------------------------------------------------------*
 | Area calculation                                         chfoe 05/08 |
 *----------------------------------------------------------------------*/
double FLD::UTILS::FluidImpedanceBc::Area( double& density, double& viscosity, const int condid )
{
  // fill in parameter list for subsequent element evaluation
  // there's no assembly required here
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",FLD::calc_area);
  eleparams.set<double>("area",0.0);
  eleparams.set<double>("viscosity",0.0);
  eleparams.set<double>("density",0.0);

  const std::string condstring("ImpedanceCond");

  discret_->EvaluateCondition(eleparams,condstring,condid);

  double actarea = eleparams.get<double>("area");
  density = eleparams.get<double>("density");
  viscosity = eleparams.get<double>("viscosity");

  // find the lowest proc number that knows the material data
  int numproc = discret_->Comm().NumProc();
  int theproc = -1;   // the lowest proc that has the desired information
  std::vector<double> alldens(numproc);

  discret_->Comm().GatherAll( &density,&(alldens[0]),1 );
  for(int i=0; i<numproc; i++)
    if( alldens[i] > 0.0 )
    {
      theproc = i;
      break;
    }
  if(theproc < 0)
    dserror("Something parallel went terribly wrong!");

  // do the actual communication of density ...
  discret_->Comm().Broadcast(&density,1,theproc);
  // ... and viscosity
  discret_->Comm().Broadcast(&viscosity,1,theproc);

  // get total area in parallel case
  double pararea = 0.0;
  discret_->Comm().SumAll(&actarea,&pararea,1);

  if (myrank_ == 0)
  {
    std::cout << "Impedance condition Id: " << condid << " area = " << pararea << std::endl;
  }
  return pararea;
}//FluidImplicitTimeInt::Area



/*----------------------------------------------------------------------*
 |  Calculate Impedance depending on history           ac | chfoe 05/08 |
 *----------------------------------------------------------------------*/
/*

  This routine contains major parts of the following paper:

  Olufsen et al.: "Numerical Simulation and Experimental Validation of
  Blood Flow in Arteries with Structured-Tree Outflow Conditions",
  Annals of Biomedical Eingineering, Vol. 28, pp. 1281--1299, 2000.

  Basic Idea:
  (1) Calculate the root impedance of the cut-off vascular tree from an
      analytical solution within the tree for different frequencies
  (2) Transfer the frequency-domain solution back to the time domain via
      inverse Fourier transform

*/
void FLD::UTILS::FluidImpedanceBc::Impedances( const double area, const double density, const double viscosity )
{
  // setup variables
  std::vector<std::complex<double> > frequencydomain;  // impedances in the frequency domain
  frequencydomain.resize(cyclesteps_,0);
  std::vector<std::complex<double> > timedomain;       // impedances in the time domain
  timedomain.resize(cyclesteps_,0);

  // size: number of generations
  std::complex<double> zparent (0,0);  // only as argument

  // to store already calculated impedances (to be developed)
  std::map<const double, std::complex<double> > zstored;     // as argument

  // set up some geometry data
  double radius = sqrt(area/PI);  // the radius to which the artery tree is connected

  // calculate DC (direct current) component depending on type of tree
  if ( treetype_ == "lung" )
    frequencydomain[0] = DCLungImpedance(0,radius,termradius_,density,viscosity,zstored);

  if ( treetype_ == "artery" )
    frequencydomain[0] = DCArteryImpedance(0,radius,termradius_,density,viscosity,zstored);

  if ( treetype_ == "windkessel" )
    frequencydomain[0] = WindkesselImpedance(0);

  // erase all DC entities from the stored data
  zstored.clear();

  //Loop over some frequencies making a call to Impedance, w=2*pi*k/T
  // this results in a field like
  // frequencydomain =
  // [ Z(w_0)  Z(w_1)  Z(w_2)  Z(w_3)  ...  Z(w_cyclesteps) ]
  //
  // note that we also need the complex conjugated ones, i.e.
  // Z(-w_1) = conj(Z(w_1)  to Z(-w_cyclesteps) = conj(Z(w_cyclesteps)
   for (int k=1; k<cyclesteps_; k++)
   {
     int generation=0;
     if ( treetype_ == "lung" )
       frequencydomain[k] = LungImpedance(k,generation,radius,termradius_,density,viscosity,zstored);

     if ( treetype_ == "artery" )
       frequencydomain[k] = ArteryImpedance(k,generation,radius,termradius_,density,viscosity,zstored);

     if ( treetype_ == "windkessel" )
       frequencydomain[k] = WindkesselImpedance(k);
   }


  // --------------------------------------------------------------------------
  // inverse Fourier transform
  // idea:
  //
  //                    cyclesteps                i omega_k t
  //      imp(x,t) = sum          IMP(x,omega_k) e
  //                   -cyclesteps
  //
  // with omega_k = 2 PI k / T
  // --------------------------------------------------------------------------

  double constexp = 2*PI/cyclesteps_;
  double realpart = cos(constexp);
  double imagpart = sin(constexp);
  std::complex<double> eiwt (realpart,imagpart);  // e^{i w_k / T}

  // now for all discrete times get the influence of all frequencies that have
  // been pre-computed
  for (int timefrac = 0; timefrac < cyclesteps_; timefrac++)
  {
    for (int k=0; k<cyclesteps_/2; k++)
      timedomain[timefrac] += std::pow(eiwt,timefrac*k) * frequencydomain[k];

    // and now the conjugated ones are added
    for (int k=1; k<cyclesteps_/2; k++)
      timedomain[timefrac] += std::pow(eiwt,-timefrac*k) * conj(frequencydomain[k]);
  }

  //Get Real component of TimeDomain, imaginary part should be anyway zero.
  for (int i=0; i<cyclesteps_; i++)
  {
    impvalues_[i]=real(timedomain[i]);
    if (abs(imag(timedomain[i])) > 1E-4 )
      std::cout << "error in imaginary part is = " << timedomain[i] << std::endl;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Flow rate calculation                                      ac 03/08 |
 *----------------------------------------------------------------------*/
/*!
  modified by chfoe 04/08

  Calculate the flow rate across an impedance boundary surface

  Flow rates are
  (1) calculated for single element surfaces
  (2) added up over the elements of the single procs
  (3) communicated and added over the procs
  (4) and finally stored within the vector 'flowrates_'

  The vector of the flowrates holds the flow rate history of the
  very last cycle!

*/
void FLD::UTILS::FluidImpedanceBc::FlowRateCalculation(const double time, const int condid)
{
  // fill in parameter list for subsequent element evaluation
  // there's no assembly required here
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",FLD::calc_flowrate);

  // get a vector layout from the discretization to construct matching
  // vectors and matrices local <-> global dof numbering
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // create vector (+ initialization with zeros)
  Teuchos::RCP<Epetra_Vector> flowrates = LINALG::CreateVector(*dofrowmap,true);

  const std::string condstring("ImpedanceCond");
  discret_->EvaluateCondition(eleparams,flowrates,condstring,condid);

  double local_flowrate = 0.0;
  for (int i=0; i < dofrowmap->NumMyElements(); i++)
  {
    local_flowrate +=((*flowrates)[i]);
  }

  double flowrate = 0.0;
  dofrowmap->Comm().SumAll(&local_flowrate,&flowrate,1);

  if(treetype_ == "windkessel_freq_indp" || treetype_ == "resistive")
  {
    Q_np_ = flowrate;
  }
  else
  {
    Q_np_ = flowrate;

    // fill vector of flowrates calculated within the last cycle
    if (time < period_) // we are within the very first cycle
    {
      // we are now in the initial fill-in phase
      flowrates_->at(flowratespos_) = flowrate;
    }
    else
    {
      // we are now in the post-initial phase
      // replace the element that was computed exactly a cycle ago
      int pos = flowratespos_ % cyclesteps_;
      flowrates_->at(pos) = flowrate;
    }
  }

//  if (myrank_ == 0)
//    printf("Impedance condition Id: %d Flowrate = %f \t time: %f \n",condid,flowrate, time);

  return;
}//FluidImplicitTimeInt::FlowRateCalculation



/*----------------------------------------------------------------------*
 |  Apply Impedance to outflow boundary                      thon 10/15 |
 *----------------------------------------------------------------------*/
/*

  This routine contains major parts of the following paper:

  Olufsen et al.: "Numerical Simulation and Experimental Validation of
  Blood Flow in Arteries with Structured-Tree Outflow Conditions",
  Annals of Biomedical Eingineering, Vol. 28, pp. 1281--1299, 2000.

  Basic Idea:
  (1) Evaluate convolution integral over one cycle to obtain outflow
      pressure
  (2) Apply this pressure as a Neumann-load type at the outflow boundary

*/
double FLD::UTILS::FluidImpedanceBc::OutflowBoundary(const double dta,const int condid)
{
  double pressure=0.0;

  if(treetype_ == "windkessel_freq_indp")
  {
    double Qc_n   = (P_n_ - Pc_n_)/R1_ - Pc_n_/R2_;
    double Pceq_n =  Pc_n_ + dta/(2.0*C_)*Qc_n;
    Pc_np_  = (Q_np_ + Pceq_n*2.0*C_/dta)/(2.0*C_/dta + 1.0/R2_);
    P_np_ = Pc_n_ + Q_np_*R1_;

    pressure = P_np_;
  }
  else if(treetype_ == "resistive")
  {
    double R = R1_;
    pressure = R*Q_np_;
  }
  else if(treetype_ == "lung" or treetype_ == "artery")
  {
    // evaluate convolution integral

    // the convolution integral
    for (int j=0; j<cyclesteps_; j++)
    {
      int qindex = ( flowratespos_+j ) % cyclesteps_;

      // flowrate is zero if not yet a full cycle is calculated
      double actflowrate = 0.0;
      if (qindex > (int)flowrates_->size()-1)
      actflowrate = 0.0;
      else
      actflowrate = (*flowrates_)[qindex];

      int zindex = -1-j+cyclesteps_;
      pressure += impvalues_[zindex] * actflowrate * dta; // units: pressure x time
    }

    pressure = pressure/period_; // this cures the dimension; missing in Olufsen paper

    P_np_=pressure;
  }
  else if(treetype_ == "windkessel")
  {
    // evaluate convolution integral

    // the convolution integral
    for (int j=0; j<cyclesteps_; j++)
    {
      int qindex = ( flowratespos_+j ) % cyclesteps_;

      // flowrate is zero if not yet a full cycle is calculated
      double actflowrate = 0.0;
      if (qindex > (int)flowrates_->size()-1)
      actflowrate = 0.0;
      else
      actflowrate = (*flowrates_)[qindex];

      int zindex = -1-j+cyclesteps_;
      pressure += impvalues_[zindex] * actflowrate * dta; // units: pressure x time
    }

    pressure = pressure/period_; // this cures the dimension; missing in Olufsen paper

    P_np_=pressure;
  }
  // call the element to apply the pressure
  Teuchos::ParameterList eleparams;
  // action for elements
  eleparams.set<int>("action",FLD::Outletimpedance);

  eleparams.set("WindkesselPressure",pressure);

//  if (myrank_ == 0)
//    printf("Impedance condition Id: %d Pressure from convolution = %f\t time = %f\n",condid,pressure, time);


  impedancetbc_->PutScalar(0.0);

  discret_->EvaluateCondition(eleparams,impedancetbc_,"ImpedanceCond",condid);


  // ---------------------------------------------------------------------
  // initialize the linearization matrix (iff not done already)
  // ---------------------------------------------------------------------
  if ( not impedancetbcsysmat_->Filled() )
  {
    if(treetype_ == "windkessel")
    {
      //calculate dQ/du = ( \phi o n )_Gamma
      const Epetra_Map* dofrowmap = discret_->DofRowMap();
      Teuchos::RCP<Epetra_Vector> dQdu = LINALG::CreateVector(*dofrowmap,true);

      Teuchos::ParameterList eleparams2;
      // action for elements
      eleparams2.set<int>("action",FLD::dQdu);

      discret_->EvaluateCondition(eleparams2,dQdu,"ImpedanceCond",condid);

      //now move dQdu to one proc
      Teuchos::RCP<Epetra_Map> dofrowmapred =  LINALG::AllreduceEMap(*dofrowmap);
      Teuchos::RCP<Epetra_Vector> dQdu_full = Teuchos::rcp(new Epetra_Vector(*dofrowmapred,true));

      LINALG::Export(*dQdu, *dQdu_full); //!!! add off proc components


      //calculate d wk/du = d/du ( (v,n)_gamma n,phi)_Gamma were (d wk/du)_i,j= timefacs* (phi_i,n)_Gamma * (phi_j,n)_Gamma
      // Note: this derivative cannot be build on element level, hence we have to do it here!
      impedancetbcsysmat_->Zero();

      const double tfaclhs = eleparams2.get<double>("tfaclhs",0.0);
      const double tfaclhs_wkfac = tfaclhs * dta / period_ * impvalues_[cyclesteps_-1];

      for (int lid = 0; lid <dQdu->MyLength();lid++)
      {
        const int gid    =  dofrowmap->GID(lid);
        const double val = (*dQdu)[lid];
        if (abs(val)>1e-15)
        {
          for (int lid2 = 0; lid2 <dQdu_full->MyLength();lid2++)
          {
            const int gid2    =  dofrowmapred->GID(lid2);
            const double val2 = (*dQdu_full)[lid2];
            const double actmatentry = tfaclhs_wkfac*val*val2;
            if (abs(actmatentry)>1e-15)
              impedancetbcsysmat_->Assemble(actmatentry,gid,gid2);
          }
        }
      }

      impedancetbcsysmat_->Complete();
      //  std::cout<<__FILE__<<__LINE__<<*((Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(impedancetbcsysmat_))->EpetraMatrix())<<std::endl;
    }
    else //still needs to be implemented
    {
      impedancetbcsysmat_->Zero();
      impedancetbcsysmat_->Complete();
    }
  }

  return pressure;
} //FluidImplicitTimeInt::OutflowBoundary

/*----------------------------------------------------------------------*
 |  Update flowrate and pressure vector                       Thon 07/15 |
 *----------------------------------------------------------------------*/
void FLD::UTILS::FluidImpedanceBc::TimeUpdateImpedance(const double time)
{
  const double actpressure = P_np_;
  // -------------------------------------------------------------------
  // fill the pressure vector
  // -------------------------------------------------------------------

  // we are now in the post-initial phase

  // set the begining of the flowrate vector as the end
  // this is due to the periodicity reason

  // fill vector of flowrates calculated within the last cycle
  if (time < period_) // we are within the very first cycle
  {
    // we are now in the initial fill-in phase
    pressures_->at(pressurespos_) = actpressure;
  }
  else
  {
    // we are now in the post-initial phase
    // replace the element that was computed exactly a cycle ago
    int pos = pressurespos_ % (cyclesteps_);

    if ( (fmod(time+1e-14,period_)-1e-14) < 1e-14*time) //iff we are at the beginning of a new period
    {
      dP_ = actpressure - P_0_;
      P_0_ = actpressure;
    }

    pressures_->at(pos) = actpressure;
  }

  P_n_ = P_np_;
  Pc_n_  = Pc_np_;
  Q_n_  = Q_np_;

  return;
} //FluidImplicitTimeInt::OutflowBoundary


/*----------------------------------------------------------------------*
|  prepare time step of impedance conditions                  Thon 07/15 |
*----------------------------------------------------------------------*/
void FLD::UTILS::FluidImpedanceBc::PrepareTimeStepImpedance(const double time)
{
  pressurespos_++;
  flowratespos_++;

  if (time < period_) // we are within the very first cycle
  {
    // we are now in the initial fill-in phase

    // new data is appended to our pressure vector
    pressures_->push_back(0.0);
    if ( (unsigned)pressurespos_ != (pressures_->size()-1) )
      dserror("size of pressures_ vector does not fit!");

    // new data is appended to our flowrates vector
    flowrates_->push_back(0.0);
    if ( (unsigned)flowratespos_ != (flowrates_->size()-1) )
      dserror("size of flowrates_ vector does not fit!");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Split linearization matrix to a BlockSparseMatrixBase   Thon 07/15  |
 *----------------------------------------------------------------------*/
/*!
*/
void FLD::UTILS::FluidImpedanceBc::UseBlockMatrix(Teuchos::RCP<std::set<int> >     condelements,
                                                  const LINALG::MultiMapExtractor& domainmaps,
                                                  const LINALG::MultiMapExtractor& rangemaps,
                                                  bool                             splitmatrix)
{
  Teuchos::RCP<LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy> > mat;

  if (splitmatrix)
  {
    // (re)allocate system matrix
    mat = Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>(domainmaps,rangemaps,108,false,true));
    mat->SetCondElements(condelements);
    impedancetbcsysmat_ = mat;
  }

}

/*----------------------------------------------------------------------*
 |  Update of residual vector and its linearization         Thon 07/15  |
 *----------------------------------------------------------------------*/
/*!
*/
void FLD::UTILS::FluidImpedanceBc::UpdateResidualAndSysmat(Teuchos::RCP<Epetra_Vector>& residual, Teuchos::RCP<LINALG::SparseOperator>& sysmat)
{
  residual->Update(1.0,*impedancetbc_,1.0);

  sysmat->Add(*impedancetbcsysmat_,true,1.0,1.0);
}



/*----------------------------------------------------------------------*
 |  windkessel impedance for wave number k                  chfoe 06/08 |
 *----------------------------------------------------------------------*/
/*!
determine impedance for every wave number from simple windkessel model

this method is an alternative to the impedance calculation from arterial trees
 */
std::complex<double> FLD::UTILS::FluidImpedanceBc::WindkesselImpedance(const double k)
{
  double R1 = R1_;  // proximal resistance
  double R2 = R2_;  // distal resistance
  double C = C_;  // capacitance
  std::complex<double> imag(0,1), imp;

  double omega = 2.0*PI*k/period_; // circular frequency

  imp = ( R1+R2 + imag*omega*C*R1*R2 ) / ( 1.0 + imag*omega*C*R2 );

  return imp;
}



/*----------------------------------------------------------------------*
 |  arterial tree impedance for wave number k               chfoe 05/08 |
 *----------------------------------------------------------------------*/
/*!
build up artery tree and calculate root impedance recursively for a given
frequency determined by the wave number k

Data taken from
  Olufsen et al.: "Numerical Simulation and Experimental Validation of
  Blood Flow in Arteries with Structured-Tree Outflow Conditions",
  Annals of Biomedical Eingineering, Vol. 28, pp. 1281--1299, 2000.

Further also needed
  Olufsen, M.S.: "Structured tree outflow condition for blood flow
  in larger systematic arteries", American Physiological Society, 1999.

parameter:

int     k           (i)   wavenumber
int     generation  (i)   generation of actual vessel (root is generation 1)
double  radius      (i)   radius of present vessel
double  termradius  (i)   termination radius (minimal radius)
double  density     (i)   the fluid's density
double  viscosity   (i)   the fluid's viscosity

returns impedance of present vessel, after recursive call: impedance of
root vessel for given frequency

*/
std::complex<double> FLD::UTILS::FluidImpedanceBc::ArteryImpedance(int k,
                   int generation,
                   double radius,
                   double termradius,
                   double density,
                   double viscosity,
                   std::map<const double,std::complex<double> > zstored)
{
  // general data
  double lscale = 50.0; // length to radius ratio
  double alpha = 0.50;   // right daughter vessel ratio
  double beta = 0.85;    // left daughter vessel ratio

  // some auxiliary stuff
  std::complex<double> koeff, imag(0,1), cwave;

  // terminal resistance is assumed zero
  std::complex<double> zterminal (0,0);

  double omega = 2.0*PI*k/period_;

  // build up geometry of present generation
  double area = radius*radius*PI;
  double length = lscale * radius;


  // get impedances of downward vessels ...
  //*****************************************
  generation++;  // this is the next generation

  // left hand side:
  double leftradius  = alpha*radius;
  double rightradius = beta*radius;
  std::complex<double> zleft;
  std::complex<double> zright;
  bool terminated = false;

  // only if both vessels are smaller than the limit truncate
  if (leftradius < termradius && rightradius < termradius)
    terminated = true;
  else
  {
    std::map<const double,std::complex<double> >::iterator iter = zstored.find(leftradius);
    if(iter != zstored.end()) // impedance of this left radius was already computed, is in map
      zleft = iter->second;
    else                      // left hand side impedance not yet stored
    {
      zleft  = ArteryImpedance(k,generation,leftradius,termradius,density,viscosity,zstored);
      zstored.insert( make_pair( leftradius, zleft ) );
    }

    iter = zstored.find(rightradius);
    if(iter != zstored.end()) // impedance of this right radius was already computed, is in map
      zright = iter->second;
    else                      // right hand side impedance not yet stored
    {
      zright = ArteryImpedance(k,generation,rightradius,termradius,density,viscosity,zstored);
      zstored.insert( make_pair( rightradius, zright ) );
    }
  }


  // ... combine this to the impedance at my downstream end ...
  //*************************************************************
  // note, we truncate both terminal vessels at once!
  std::complex<double> zdown;
  if (terminated)
    zdown = zterminal;
  else
    zdown = 1.0 / (1.0/zleft + 1.0/zright);


  // ... and compute impedance at my upstream end!
  //*************************************************************
  double compliance = 1.5*area / ( R1_ * exp(R2_*radius) + C_ );

  double sqrdwo = radius*radius*omega/(viscosity/density);  // square of Womersley number
  double wonu = sqrt(sqrdwo);                     // Womersley number itself

  if (wonu > 4.0)
    koeff = 1.0 - (2.0/wonu)/sqrt(imag);
  else
    koeff = 1.0 / ( 1.333333333333333333 - 8.0*imag/ (wonu*wonu) );

  // wave speed of this frequency in present vessel
  cwave=sqrt( area*koeff / (density*compliance) );

  //Convenience coefficient
  std::complex<double> gcoeff = compliance * cwave;

  // calculate impedance of this, the present vessel
  std::complex<double> argument = omega*length/cwave;
  std::complex<double> zparent  = (imag/gcoeff * sin(argument) + zdown*cos(argument) ) /
                             ( cos(argument) + imag*gcoeff*zdown*sin(argument) );

  return zparent;
}




/*----------------------------------------------------------------------*
 |  impedance w.r.t. constant flow                          chfoe 05/08 |
 *----------------------------------------------------------------------*/
/*!
  The special case of omega=0 is covered here, i.e. what is the impedance
  of the tree for constant flow that is a very special harmonic function
  with zero frequency.

  This case is also called direct current component in the classical
  analogy.
*/
std::complex<double> FLD::UTILS::FluidImpedanceBc::DCArteryImpedance(int generation,
               double radius,
               double termradius,
               double density,
               double viscosity,
               std::map<const double,std::complex<double> > zstored)
{
  // general data
  double lscale = 50.0; // length to radius ratio
  double alpha = 0.50;   // right daughter vessel ratio
  double beta = 0.85;    // left daughter vessel ratio
  double mu = viscosity; // dynamic (physical) viscosity

  // terminal resistance is assumed zero
  std::complex<double> zterminal (0,0);


  // get impedances of downward vessels ...
  //*****************************************
  generation++;  // this is the next generation

  double leftradius  = alpha*radius;
  double rightradius = beta*radius;
  std::complex<double> zleft;
  std::complex<double> zright;
  bool terminated = false;

  // only if both vessels are smaller than the limit truncate
  if (leftradius < termradius && rightradius < termradius)
    terminated = true;
  else
  {
    std::map<const double,std::complex<double> >::iterator iter = zstored.find(leftradius);
    if(iter != zstored.end()) // impedance of this left radius was already computed, is in map
      zleft = iter->second;
    else                      // left hand side impedance not yet stored
    {
      zleft  = DCArteryImpedance(generation,leftradius,termradius,density,viscosity,zstored);
      zstored.insert( make_pair( leftradius, zleft ) );
    }

    iter = zstored.find(rightradius);
    if(iter != zstored.end()) // impedance of this right radius was already computed, is in map
      zright = iter->second;
    else                      // right hand side impedance not yet stored
    {
      zright = DCArteryImpedance(generation,rightradius,termradius,density,viscosity,zstored);
      zstored.insert( make_pair( rightradius, zright ) );
    }
  }


  // ... combine this to the impedance at my downstream end ...
  //*************************************************************
  // note, we truncate both terminal vessels at once!
  std::complex<double> zdown;
  if (terminated)
    zdown = zterminal;
  else
    zdown = 1.0 / (1.0/zleft + 1.0/zright);


  // ... and compute dc impedance at my upstream end!
  //*************************************************************

  // calculate dc impedance of this, the present, vessel
  std::complex<double> zparentdc = 8.0 * mu * lscale / ( M_PI*radius*radius*radius ) + zdown;
  // DEBUG output
  //std::cout << "generation: " << generation << std::endl;
  return zparentdc;
}


/*----------------------------------------------------------------------*
 |  ?????????????????????????                                  ac 06/08 |
 *----------------------------------------------------------------------*/
/*!
What is supposed to happen within these lines?
*/
std::complex<double> FLD::UTILS::FluidImpedanceBc::LungImpedance(int k,
                 int generation,
                 double radius,
                 double termradius,
                 double density,
                 double viscosity,
                 std::map<const double,std::complex<double> > zstored)
{
   // general data
  double lscale = 5.8; // length to radius ratio
  double alpha = 0.876;   // right daughter vessel ratio
  double beta = 0.686;    // left daughter vessel ratio

  // some auxiliary stuff
  std::complex<double> imag(0,1), Z1, Z2, Z3, ZW;
  std::complex<double> koeff, cwave;
  // terminal resistance is assumed zero
  std::complex<double> zterminal (0,0);

  double omega = 2.0*PI*k/period_;

  // Alveolar impedance
  //double G=0.96e-2;
  //double H=8e-2;
  //double athing=(2.0/PI)*atan(H/G);
  //complex<double> zterminal=(G-imag*H)/(pow(omega,athing));

  // build up geometry of present generation
  double area = radius*radius*PI;
  double length = lscale * radius;

#if 0
  double h=-0.0057*radius*radius+0.2096*radius+0.0904;
  double E= 0.0033;
  double compliance=2.0*PI*radius*radius*radius/(3.0*h*E);
#else
  double h= H1_*radius*radius+H2_*radius+H3_; //calculate wallthickness as quadratic polynomial
  double E= E_;
  double compliance=2.0*PI*radius*radius*radius/(3.0*h*E);
#endif

  // get impedances of downward vessels ...
  //*****************************************
  generation++;  // this is the next generation
  //std::cout<<"generation "<<generation<<std::endl;
  // left hand side:
    double leftradius  = alpha*radius;
    double rightradius = beta*radius;
    std::complex<double> zleft;
    std::complex<double> zright;
    bool terminated = false;

    if (leftradius < termradius && rightradius < termradius)
            terminated = true;
          else
          {
            std::map<const double,std::complex<double> >::iterator iter = zstored.find(leftradius);
            if(iter != zstored.end()) // impedance of this left radius was already computed, is in map
              zleft = iter->second;
            else                      // left hand side impedance not yet stored
            {
              zleft  = LungImpedance(k,generation,leftradius,termradius,density,viscosity,zstored);
              zstored.insert( make_pair( leftradius, zleft ) );
            }

            iter = zstored.find(rightradius);
            if(iter != zstored.end()) // impedance of this right radius was already computed, is in map
              zright = iter->second;
            else                      // right hand side impedance not yet stored
            {
              zright = LungImpedance(k,generation,rightradius,termradius,density,viscosity,zstored);
              zstored.insert( make_pair( rightradius, zright ) );
            }
          }
    // only if both vessels are smaller than the limit truncate

    // ... combine this to the impedance at my downstream end ...
    //*************************************************************
    // note, we truncate both terminal vessels at once!
    std::complex<double> zdown;
    if (terminated)
      zdown = zterminal;
    else
      zdown = 1.0 / (1.0/zleft + 1.0/zright);

  // ... and compute impedance at my upstream end!
  //*************************************************************

    double sqrdwo = radius*radius*omega/(viscosity/density);  // square of Womersley number
    double wonu = sqrt(sqrdwo);                     // Womersley number itself

    if (wonu > 4.0)
    koeff = 1.0 - (2.0/wonu)/sqrt(imag);
    else
    koeff = 1.0 / ( 1.333333333333333333 - 8.0*imag/ (wonu*wonu) );

    // wave speed of this frequency in present vessel
    cwave=sqrt( area*koeff / (density*compliance) );

    //Convenience coefficient
    std::complex<double> gcoeff = compliance * cwave;

    // calculate impedance of this, the present vessel
    std::complex<double> argument = omega*length/cwave;
    std::complex<double> zparent  = (imag/gcoeff * sin(argument) + zdown*cos(argument) ) /
                                     ( cos(argument) + imag*gcoeff*zdown*sin(argument) );

  // calculate impedance of this, the present vessel

  /*ZW=rw+1.0/(imag*omega*cw)+imag*omega*lw;
  Z1=1.0/(imag*omega*C+1.0/ZW);
  Z2=(R+imag*omega*L)/2.0+zdown;
  Z3=1.0/(1.0/Z1+1.0/Z2);
  zparent=(R+imag*omega*L)/2.0+Z3;*/

  return zparent;
} //BioFluidImplicitTimeInt::LungImpedance





/*----------------------------------------------------------------------*
 |  ?????????????????????????                                  ac 06/08 |
 *----------------------------------------------------------------------*/
/*!
What is supposed to happen within these lines?
*/
std::complex<double> FLD::UTILS::FluidImpedanceBc::DCLungImpedance(int generation,
                 double radius,
                 double termradius,
                 double density,
                 double viscosity,
                 std::map<const double,std::complex<double> > zstored)
{
  //general data
  double lscale = 5.8; // length to radius ratio
  double alpha = 0.876;   // right daughter vessel ratio
  double beta = 0.686;    // left daughter vessel ratio

  double mu = viscosity; // dynamic (physical) viscosity
  generation++;
  // terminal resistance is assumed zero
  std::complex<double> zterminal (0,0);

    double leftradius  = alpha*radius;
    double rightradius = beta*radius;
    std::complex<double> zleft;
    std::complex<double> zright;
    bool terminated = false;

    if (leftradius < termradius && rightradius < termradius)
           terminated = true;
         else
         {
           std::map<const double,std::complex<double> >::iterator iter = zstored.find(leftradius);
           if(iter != zstored.end()) // impedance of this left radius was already computed, is in map
             zleft = iter->second;
           else                      // left hand side impedance not yet stored
           {
             zleft  = DCLungImpedance(generation,leftradius,termradius,density,viscosity,zstored);
             zstored.insert( make_pair( leftradius, zleft ) );
           }

           iter = zstored.find(rightradius);
           if(iter != zstored.end()) // impedance of this right radius was already computed, is in map
             zright = iter->second;
           else                      // right hand side impedance not yet stored
           {
             zright = DCLungImpedance(generation,rightradius,termradius,density,viscosity,zstored);
             zstored.insert( make_pair( rightradius, zright ) );
           }
         }

    // ... combine this to the impedance at my downstream end ...
    //*************************************************************
    // note, we truncate both terminal vessels at once!
    std::complex<double> zdown;
    if (terminated)
      zdown = zterminal;
    else
      zdown = 1.0 / (1.0/zleft + 1.0/zright);


    // ... and compute dc impedance at my upstream end!
    //*************************************************************

    // calculate dc impedance of this, the present vessel
    std::complex<double>  zparentdc = 8.0 * mu * lscale / ( M_PI*radius*radius*radius ) + zdown;
    return zparentdc;
}//FluidImplicitTimeInt::DCLungImpedance


/*----------------------------------------------------------------------*
 |  Set windkessel parameters                              ismail 02/10 |
 *----------------------------------------------------------------------*/
void FLD::UTILS::FluidImpedanceBc::SetWindkesselParams(Teuchos::ParameterList & params)
{

  R1_ = params.get<double> ("R1");
  R2_ = params.get<double> ("R2");
  C_ = params.get<double> ("C");
}// FluidImpedanceWrapper::SetWindkesselParams


/*----------------------------------------------------------------------*
 |  Get windkessel parameters                              ismail 02/10 |
 *----------------------------------------------------------------------*/
void FLD::UTILS::FluidImpedanceBc::GetWindkesselParams(Teuchos::ParameterList & params)
{

  params.set<double> ("R1",R1_);
  params.set<double> ("R2",R2_);
  params.set<double> ("C",C_);

}// FluidImpedanceWrapper::SetWindkesselParams


/*----------------------------------------------------------------------*
 |  Returns results of one cardiac period                  ismail 02/10|
 *----------------------------------------------------------------------*/
void FLD::UTILS::FluidImpedanceBc::getResultsOfAPeriod(
  Teuchos::ParameterList & params,
  const int                condid)
{

  std::stringstream pstream, qstream, dpstream, endCystream;

  pstream<<"pressures"<<condid;
  qstream<<"flowrates"<<condid;
  dpstream<<"dP"<<condid;
  endCystream<<"EndOfCycle"<<condid;

  params.set<Teuchos::RCP<std::vector<double> > >(pstream.str(),pressures_ );
  params.set<Teuchos::RCP<std::vector<double> > >(qstream.str(),flowrates_);
  params.set<double> (dpstream.str(),dP_);
  bool EndOfCycle = ((pressurespos_ % (cyclesteps_)) == 0);
  params.set<bool> (endCystream.str(),EndOfCycle);

}

/*----------------------------------------------------------------------*
 |  Interpolate values of Vector1 to fit in Vector2         ismail 04/10|
 |                                                                      |
 | V ^                                                                  |
 |   |                              +                                   |
 |   |                            , .                                   |
 |   |                          ,   .                                   |
 |   |                        ,     .                                   |
 |   |                      ,+      .                                   |
 |   |                    ,  .      .                                   |
 |   |                  //   .      .                                   |
 |   |                ,      .      .                                   |
 |   |              ,+       .      .                                   |
 |   |            ,  .       .      .                                   |
 |   |          ,    .       .      .                                   |
 |   |        +      .       .      .                                   |
 |   |        .      .       .      .                                   |
 |   |        .      .       .      .                                   |
 |   |        .      .       .      .                                   |
 |   |        .      .       .      .                                   |
 |   |        .      .       .      .                                   |
 |   |        .      .       .      .                                   |
 |   |        .      .       .      .                                   |
 |   +--------+------o---//--o------+----------->                       |
 |            T(i)   .              T(i+1)      t                       |
 |            .      .              .                                   |
 |            t(m)   t(m+j)         t(m+k)                              |
 |                                                                      |
 |                                                                      |
 |  (T) is the time step of the original vector                         |
 |  (t) is the time step of the new vector                              |
 |  1 - Loop over all intervals (T(i) and T(i+1))                       |
 |  2 - Check if V2 has any time steps between T(i) and T(i+1)          |
 |  3 - Using linear interpolation check get the value of V2 at t(m+j)  |
 |                                                                      |
 | *The advantage of this method is that is works for finer and coarser |
 |  interpolations.                                                     |
 |                                                                      |
 *----------------------------------------------------------------------*/
void FLD::UTILS::FluidImpedanceBc::interpolate(const Teuchos::RCP<const std::vector<double> > V1,
                                               Teuchos::RCP<std::vector<double> > V2,
                                               const int index1,
                                               int & index2,
                                               const double time)
{
  // Get size of V1 and V2
  int n1 = V1->size();
  int n2 = V2->size();

  double TotalTime;
  if (time<period_)
    TotalTime = time;
  else
    TotalTime = period_;

  // Get time step size of V1 and V2
  double dt1 = (TotalTime)/double(n1 - 1);
  double dt2 = (TotalTime)/double(n2 - 1);

  // defining some necessary variables
  double t1_1, t1_2;
  double v1_1, v1_2;

  // define t (time step of V2) and k (index of V2)
  double t = 0.0;
  int k    = 0;
  for (int i = 0; i< n1-1; i++)
  {
    // -----------------------------------------------------------------
    // Get V1 values at T(i) and T(i+1)
    // -----------------------------------------------------------------
    v1_1 = (*V1)[i];
    v1_2 = (*V1)[i+1];

    // -----------------------------------------------------------------
    // Calculate T(i) and T(i+1)
    // -----------------------------------------------------------------
    t1_1 = double(i)*dt1;
    t1_2 = double(i+1)*dt1;

    // -----------------------------------------------------------------
    // Evaluate V2 values between T(i) and  T(i+1)
    // -----------------------------------------------------------------
    while (t < t1_2)
    {
      // Evaluate value of V2 using Interpolation
      (*V2)[k] = (t1_2 - t)/dt1*v1_1 + (t - t1_1)/dt1*v1_2;
      // Increment k
      k++;
      // Increment t
      t += dt2;
    }
  }

  // -------------------------------------------------------------------
  // Finally resolve the last step where V2(n2) = V1(n1)
  // -------------------------------------------------------------------
  (*V2)[V2->size()-1] = (*V1)[V1->size()-1];


  // -------------------------------------------------------------------
  // Get the index of V2
  // where t = t     => dt1*index1 = dt2*index2
  //                              dt1
  //                 => index2 = ----- index1
  //                              dt2
  // -------------------------------------------------------------------
  index2 = int(double(index1)*(dt1/dt2)+1.0e-14);

}//FLD::UTILS::FluidImpedanceBc::interpolate

