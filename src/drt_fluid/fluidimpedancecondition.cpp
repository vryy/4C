/*!----------------------------------------------------------------------
\file fluidimpedancecondition.cpp
\brief evaluation of impedance vascular bc

<pre>
Maintainer: Christiane Förster
            foerster@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*----------------------------------------------------------------------*/

#include <stdio.h>

#include "fluidimpedancecondition.H"

#include "../linalg/linalg_ana.H"
#include "../drt_fluid_ele/fluid_ele_action.H"



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                     chfoe 06/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::UTILS::FluidImpedanceWrapper::FluidImpedanceWrapper(Teuchos::RCP<DRT::Discretization> actdis,
                                                         IO::DiscretizationWriter& output,
                                                         double dta) :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  output_ (output)
{
  std::vector<DRT::Condition*> impedancecond;
  discret_->GetCondition("ImpedanceCond",impedancecond);

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
    Teuchos::RCP<FluidImpedanceBc> impedancebc = Teuchos::rcp(new FluidImpedanceBc(discret_, output_, dta, condid, i) );

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


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Destructor dtor (public)                                chfoe 06/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::UTILS::FluidImpedanceWrapper::~FluidImpedanceWrapper()
{
  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Return a pointer to the pressures of condition         ismail 02/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
Teuchos::RCP<std::vector<double> > FLD::UTILS::FluidImpedanceWrapper::getPressures(int condid)
{
  return (impmap_[condid])->FluidImpedanceBc::getPressures();
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Returns results of one cardiac period                   ismail 02/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidImpedanceWrapper::getResultsOfAPeriod(
  Teuchos::ParameterList & params)
{
  // get an iterator to my map
  std::map<const int, Teuchos::RCP<class FluidImpedanceBc> >::iterator mapiter;

  for (mapiter = impmap_.begin(); mapiter != impmap_.end(); mapiter++ )
  {
    mapiter->second->FluidImpedanceBc::getResultsOfAPeriod(params,mapiter->first);
  }
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Wrap flow rate calculation                              chfoe 06/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidImpedanceWrapper::FlowRateCalculation(double time, double dta)
{
  // get an iterator to my map
  std::map<const int, Teuchos::RCP<class FluidImpedanceBc> >::iterator mapiter;

  for (mapiter = impmap_.begin(); mapiter != impmap_.end(); mapiter++ )
  {
    mapiter->second->FluidImpedanceBc::FlowRateCalculation(time,dta,mapiter->first);
  }
  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Wrap outflow boundary pressure application              chfoe 06/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidImpedanceWrapper::OutflowBoundary(double time, double dta, double theta)
{
  std::map<const int, Teuchos::RCP<class FluidImpedanceBc> >::iterator mapiter;

  for (mapiter = impmap_.begin(); mapiter != impmap_.end(); mapiter++ )
  {
    mapiter->second->FluidImpedanceBc::OutflowBoundary(time,dta,theta,mapiter->first);
  }
  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Wrap impedances calculation                            ismail 02/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidImpedanceWrapper::Impedances()
{
  // get an iterator to my map
  std::map<const int, Teuchos::RCP<class FluidImpedanceBc> >::iterator mapiter;

  for (mapiter = impmap_.begin(); mapiter != impmap_.end(); mapiter++ )
  {
    double density=0.0, viscosity=0.0;
    int    condid = mapiter->first;
    double area = mapiter->second->FluidImpedanceBc::Area(density,viscosity,condid);
    mapiter->second->FluidImpedanceBc::Impedances(area,density,viscosity);
  }
  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Set windkessel parameters                              ismail 02/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidImpedanceWrapper::SetWindkesselParams(
  Teuchos::ParameterList  & params,
  int              condid)
{
  // -------------------------------------------------------------------
  // set the windkessel params associated with the coressponding
  // condition ID
  // -------------------------------------------------------------------
  impmap_[condid]->SetWindkesselParams(params);

  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Get windkessel parameters                              ismail 02/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidImpedanceWrapper::GetWindkesselParams(
  Teuchos::ParameterList  & params,
  int                       condid)
{
  // -------------------------------------------------------------------
  // set the windkessel params associated with the coressponding
  // condition ID
  // -------------------------------------------------------------------
  impmap_[condid]->GetWindkesselParams(params);

  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Wrap update of residual                                 chfoe 06/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidImpedanceWrapper::UpdateResidual(Teuchos::RCP<Epetra_Vector>  residual )
{
  std::map<const int, Teuchos::RCP<class FluidImpedanceBc> >::iterator mapiter;

  for (mapiter = impmap_.begin(); mapiter != impmap_.end(); mapiter++ )
  {
    mapiter->second->FluidImpedanceBc::UpdateResidual(residual);
  }
  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Wrap restart writing                                    chfoe 06/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidImpedanceWrapper::WriteRestart( IO::DiscretizationWriter&  output )
{
  std::map<const int, Teuchos::RCP<class FluidImpedanceBc> >::iterator mapiter;

  for (mapiter = impmap_.begin(); mapiter != impmap_.end(); mapiter++ )
  {
    mapiter->second->FluidImpedanceBc::WriteRestart(output,mapiter->first);
  }
  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Wrap restart reading                                    chfoe 06/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidImpedanceWrapper::ReadRestart( IO::DiscretizationReader& reader)
{
  std::map<const int, Teuchos::RCP<class FluidImpedanceBc> >::iterator mapiter;

  for (mapiter = impmap_.begin(); mapiter != impmap_.end(); mapiter++ )
    mapiter->second->FluidImpedanceBc::ReadRestart(reader,mapiter->first);

  return;
}



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                     chfoe 04/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::UTILS::FluidImpedanceBc::FluidImpedanceBc(Teuchos::RCP<DRT::Discretization> actdis,
                                               IO::DiscretizationWriter& output,
                                               double dta,
                                               int condid,
                                               int numcond) :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  output_ (output)
{
  // ---------------------------------------------------------------------
  // read in all impedance conditions
  // ---------------------------------------------------------------------
  std::vector<DRT::Condition*> impedancecond;
  discret_->GetCondition("ImpedanceCond",impedancecond);

  // ---------------------------------------------------------------------
  // read in all possible impedance calibrating conditions
  // ---------------------------------------------------------------------
  std::vector<DRT::Condition*> impedance_calb_cond;
  discret_->GetCondition("ImpedanceCalbCond",impedance_calb_cond);
  IsPrecalibrated_ = false;

  // ---------------------------------------------------------------------
  // get time period length, steps per cycle and initialise flowratespos
  // ---------------------------------------------------------------------
  period_       = (impedancecond[numcond])->GetDouble("timeperiod");
  cyclesteps_   = (int)(period_/dta+0.5);
  flowratespos_ = 0;
  pressurespos_ = 0;
  dP_           = 0.0;
  dta_          = dta;

  // ---------------------------------------------------------------------
  // get relevant data from impedance condition
  // ---------------------------------------------------------------------
  treetype_ = *((impedancecond[numcond])->Get<std::string>("tree"));
  termradius_ = (impedancecond[numcond])->GetDouble("termradius");

  // 'material' parameters required for artery tree and for 3 element windkessel
  R1_ = (impedancecond[numcond])->GetDouble("R1");
  R2_ = (impedancecond[numcond])->GetDouble("R2");
  C_ = (impedancecond[numcond])->GetDouble("C");

  E_  = (impedancecond[numcond])->GetDouble("stiffness");
  H1_ = (impedancecond[numcond])->GetDouble("H1");
  H2_ = (impedancecond[numcond])->GetDouble("H2");
  H3_ = (impedancecond[numcond])->GetDouble("H3");
  // ---------------------------------------------------------------------
  // get the processor ID from the communicator
  // ---------------------------------------------------------------------
  myrank_  = discret_->Comm().MyPID();

  // ---------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // ---------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  impedancetbc_ = LINALG::CreateVector(*dofrowmap,true);

  // initialize all of the variables
  Pin_n_  = 0.0;
  Pin_np_ = 0.0;
  Pc_n_   = 0.0;
  Pc_np_  = 0.0;
  Qin_np_ = 0.0;

  if(treetype_ == "windkessel_freq_indp")
  {
    // initialize all of the variables
    Pin_n_  = 0.0;
    Pin_np_ = 0.0;
    Pc_n_   = 0.0;
    Pc_np_  = 0.0;

    for(unsigned int i = 0; i<impedance_calb_cond.size(); i++)
    {
      if(impedance_calb_cond[i]->GetInt("ConditionID") == condid)
      {
        IsPrecalibrated_ = true;
        Pin_n_  = impedance_calb_cond[i]->GetDouble("Pin_n");
        Pin_np_ = impedance_calb_cond[i]->GetDouble("Pin_np");
        Pc_n_   = impedance_calb_cond[i]->GetDouble("Pc_n");
        Pc_np_  = impedance_calb_cond[i]->GetDouble("Pc_np");
      }
    }

    Qin_np_ = (Pin_np_ - Pc_np_)/R1_;
    // This par might look little messy but could be fixed in the future
    if (myrank_ == 0)
    {
      printf(" Pin initially is: %f  --  %f",Pin_n_,Pin_np_);
      printf("Frequency independent windkessel condition(%d) with:\n",condid);
      printf("          R1 = %f\n",R1_);
      printf("          R2 = %f\n",R2_);
      printf("          C  = %f\n",C_);
      printf("          Pc(initialt)= %f:\n",termradius_);
    }
  }
  else
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
    Impedances(area,density,viscosity);
  }

  // ---------------------------------------------------------------------
  // initialize the pressures vecotrs
  // ---------------------------------------------------------------------
  pressures_ = Teuchos::rcp(new  std::vector<double>);
  pressures_->push_back(0.0);


  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Restart writing                                         chfoe 05/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidImpedanceBc::WriteRestart( IO::DiscretizationWriter&  output, int condnum )
{
  // condnum contains the number of the present condition
  // condition Id numbers must not change at restart!!!!

  std::stringstream stream1, stream2, stream3, stream4, stream5, dpstream;


  if(treetype_ == "windkessel_freq_indp")
  {
    stream1 << "Pin_n" << condnum;
    // write the input pressure at time step n
    output.WriteDouble(stream1.str(), Pin_n_);
    stream1 << "Pc_n" << condnum;
    // write the capacitor pressure at time step n
    output.WriteDouble(stream1.str(), Pc_n_);
    stream1 << "Pin_np" << condnum;
    // write the input pressure at time step n
    output.WriteDouble(stream1.str(), Pin_np_);
    stream1 << "Pc_np" << condnum;
    // write the capacitor pressure at time step n
    output.WriteDouble(stream1.str(), Pc_np_);

  }
  else
  {
    stream1 << "flowratesId" << condnum;
    // write the flowrates of the previous period
    output.WriteRedundantDoubleVector(stream1.str(),flowrates_);

    // also write flowratesposition of this outlet
    stream2 << "flowratesposId" << condnum;
    output.WriteInt(stream2.str(), flowratespos_);

    // write the pressures
    stream4<< "pressuresId"<<condnum;
    output.WriteRedundantDoubleVector(stream4.str(), pressures_);

    // also write pressuresposition of this outlet
    stream5 << "pressuresposId" << condnum;
    output.WriteInt(stream5.str(), pressurespos_);

    // write cyclesteps_
    output.WriteInt("ImpedanceBC_cyclesteps", cyclesteps_);
  }
  // also write vector impedancetbc_ (previously missing, gee)
  stream3 << "impedancetbc" << condnum;
  output.WriteVector(stream3.str(), impedancetbc_);
  // write time step size dta_
  output.WriteDouble("ImpedanceBC_dta", dta_);

  dpstream<<"dP"<<condnum;
  output.WriteDouble(dpstream.str(), dP_);


  /*
  double norm = 0.0;
  impedancetbc_->Norm2(&norm);
  if (!myrank_)
  {
    printf("impedancetbc_NORM: %10.5e\n",norm);
  }
  */
  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Restart reading                                         chfoe 05/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidImpedanceBc::ReadRestart( IO::DiscretizationReader& reader, int condnum  )
{
  // condnum contains the number of the present condition
  // condition Id numbers must not change at restart!!!!
  std::stringstream stream1, stream2, stream3, stream4, stream5, dpstream;

  // also read vector impedancetbc_ (previously missing, gee)
  stream3 << "impedancetbc" << condnum;
  reader.ReadVector(impedancetbc_,stream3.str());
  /*
  double norm = 0.0;
  impedancetbc_->Norm2(&norm);
  if (!myrank_)
  {
    printf("impedancetbc_NORM: %10.5e\n",norm);
  }
  */
  // old time step size
  double odta = reader.ReadDouble("ImpedanceBC_dta");

  // get time step of the current problems
  double ndta = dta_;


  // time of restart
  double t = reader.ReadDouble("time");

  // -------------------------------------------------------------------
  // Read in the pressure values and the pressure difference
  // -------------------------------------------------------------------
  stream4 << "pressuresId"<<condnum;
  stream5 << "pressuresposId" << condnum;

  // read in pressure difference
  dpstream <<"dP"<<condnum;
  dP_ = reader.ReadDouble(dpstream.str());

  // read in pressures
  reader.ReadRedundantDoubleVector(pressures_ ,stream4.str());

  // read in the pressures' position
  pressurespos_ = reader.ReadInt(stream5.str());

  // Get old Pressure Vector size
  int oPSize = (int)pressures_->size();

  // Calculate new Pressure Vector size
  int nPSize = (int)(double(oPSize)*odta/ndta);


  // evaluate the new pressure vector
  int np_pos = 0;
  Teuchos::RCP<std::vector<double> > np = Teuchos::rcp(new std::vector<double>(nPSize,0.0));
  this->interpolate(pressures_,np,pressurespos_,np_pos,t);

  // store new values in class
  pressurespos_ = np_pos;
  pressures_    = np;

  if(treetype_ == "windkessel_freq_indp")
  {
    if (IsPrecalibrated_)
      return;
    stream1 << "Pin_n" << condnum;
    // read the input pressure at time step n
    Pin_n_ = reader.ReadDouble(stream1.str());
    stream1 << "Pc_n" << condnum;
    // read the capacitor pressure at time step n
    Pc_n_  = reader.ReadDouble(stream1.str());
    stream1 << "Pin_np" << condnum;
    // read the input pressure at time step n
    Pin_np_ = reader.ReadDouble(stream1.str());
    stream1 << "Pc_np" << condnum;
    // read the capacitor pressure at time step n
    Pc_np_  = reader.ReadDouble(stream1.str());

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
    Pin_n_ = (odta - t)/odta*Pin_n_ + t/odta*Pin_np_;
  }
  else
  {
    stream1 << "flowratesId" << condnum;
    stream2 << "flowratesposId" << condnum;

    // read in flow rates
    reader.ReadRedundantDoubleVector(flowrates_ ,stream1.str());
    flowratespos_ = reader.ReadInt(stream2.str());

    // Get old flowrates Vector size
    int oQSize = (int)flowrates_->size();

    // Calculate new flowrates Vector size
    int nQSize = (int)(double(oQSize)*odta/ndta);


    // check if vector of flowrates is not empty
    if (flowrates_->size() == 0)
    dserror("could not re-read vector of flowrates");

    // old number of flowrates in vector

#if 1
    if (dta_ == odta)
    {
      if (!myrank_)
        printf("Impedance restart old and new time step are the same - life is good\n");

      OutflowBoundary(t,ndta,0.66,condnum);
      return;
    }
#endif

    // evaluate the new flow rate vector
    int nfr_pos = 0;
    Teuchos::RCP<std::vector<double> > nfr = Teuchos::rcp(new std::vector<double>(nQSize,0.0));
    this->interpolate(flowrates_,nfr,flowratespos_,nfr_pos,t);
    // store new values in class
    flowratespos_ = nfr_pos;
    flowrates_    = nfr;

    // finally, recompute the outflow boundary condition from last step
    // this way the vector need not to be stored
    OutflowBoundary(t,ndta,0.66,condnum);
  }

  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Area calculation                                         chfoe 05/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*!

*/
double FLD::UTILS::FluidImpedanceBc::Area( double& density, double& viscosity, int condid )
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



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Calculate Impedance depending on history           ac | chfoe 05/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
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
void FLD::UTILS::FluidImpedanceBc::Impedances( double area, double density, double viscosity )
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

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Flow rate calculation                                      ac 03/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
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
void FLD::UTILS::FluidImpedanceBc::FlowRateCalculation(double time, double dta, int condid)
{
  // fill in parameter list for subsequent element evaluation
  // there's no assembly required here
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",FLD::calc_flowrate);
  eleparams.set("total time",time);

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
    Qin_np_ = flowrate;
  }
  else
  {
    // fill vector of flowrates calculated within the last cycle
    if (time < period_) // we are within the very first cycle
    {
      // we are now in the initial fill-in phase
      // new data is appended to our flowrates vector
      flowratespos_++;
      flowrates_->push_back(flowrate);
    }
    else
    {
      // we are now in the post-initial phase
      // replace the element that was computed exactly a cycle ago
      flowratespos_++;
      int pos = flowratespos_ % cyclesteps_;
      (*flowrates_)[pos] = flowrate;
    }
  }
  if (myrank_ == 0)
  {
    printf("Impedance condition Id: %d Flowrate = %f \t time: %f \n",condid,flowrate, time);
  }

#if 0 // This is kept for some minor debugging purposes
  eleparams.set<int>("action",FLD::calc_pressure_bou_int);
  eleparams.set<double>("pressure boundary integral",0.0);
  eleparams.set("total time",time);


  // get elemental flowrates ...
  Teuchos::RCP<Epetra_Vector> myStoredPressures=Teuchos::rcp(new Epetra_Vector(*dofrowmap,100));

  discret_->EvaluateCondition(eleparams,myStoredPressures,condstring,condid);

  // ... as well as actual total flowrate on this proc
  double actpressure = eleparams.get<double>("pressure boundary integral");

  // get total flowrate in parallel case
  double parpressure = 0.0;
  discret_->Comm().SumAll(&actpressure,&parpressure,1);

  double density=0.0, viscosity=0.0;
  double area = Area(density,viscosity,condid);
  if (myrank_ == 0)
  {
    printf("Pressure calculation is:  %10.5e | %10.5e\n",parpressure, parpressure/area);
  }

#endif

  return;
}//FluidImplicitTimeInt::FlowRateCalculation



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Apply Impedance to outflow boundary                ac | chfoe 05/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
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
void FLD::UTILS::FluidImpedanceBc::OutflowBoundary(double time, double dta, double theta,int condid)
{
  double pressure=0.0;
  if(treetype_ == "windkessel_freq_indp")
  {
    double R1 = R1_;
    double R2 = R2_;
    double C  = C_;
    double Qc_n   = (Pin_n_ - Pc_n_)/R1 - Pc_n_/R2;
    double Pceq_n =  Pc_n_ + dta/(2.0*C)*Qc_n;
    Pc_np_  = (Qin_np_ + Pceq_n*2.0*C/dta)/(2.0*C/dta + 1.0/R2);
    Pin_np_ = Pc_n_ + Qin_np_*R1;

    pressure = Pin_np_;
  }
  else if(treetype_ == "resistive")
  {
    double R = R1_;
    pressure = R*Qin_np_;
  }
  else
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
  }
  // call the element to apply the pressure
  Teuchos::ParameterList eleparams;
  // action for elements
  eleparams.set<int>("action",FLD::Outletimpedance);

  eleparams.set("total time",time);
  eleparams.set("delta time",dta);
  eleparams.set("thsl",theta*dta);
  eleparams.set("ConvolutedPressure",pressure);

  if (myrank_ == 0)
  printf("Impedance condition Id: %d Pressure from convolution = %f\t time = %f\n",condid,pressure, time);


  impedancetbc_->PutScalar(0.0);

  const std::string condstring("ImpedanceCond");
  discret_->EvaluateCondition(eleparams,impedancetbc_,condstring,condid);

  /*
  norm = 0.0;
  impedancetbc_->Norm2(&norm);
  if (!myrank_)
  {
    std::cout<<"time: "<<time<<std::endl;
    std::cout<<"dta: "<<dta<<std::endl;
    std::cout<<"theta*dta: "<<theta*dta<<std::endl;
    std::cout<<"pressure: "<<pressure<<std::endl;
    printf("IMPEDANCE_NORM: %10.5e\n",norm);
    printf("PRESSURE:  %10.5e\n",pressure);
    //    exit(0);
  }
  */

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
    // new data is appended to our flowrates vector
    pressurespos_++;
    pressures_->push_back(pressure);
  }
  else
  {
    // we are now in the post-initial phase
    // replace the element that was computed exactly a cycle ago
    pressurespos_++;
    int pos = pressurespos_ % (cyclesteps_);
    if (pos == 0)
    {
      dP_ = pressure - (*pressures_)[pos];
      endOfCycle_ = true;
    }
    else
    {
      endOfCycle_ = false;
    }

    (*pressures_)[pos] = pressure;
  }

  discret_->ClearState();

  return;
} //FluidImplicitTimeInt::OutflowBoundary



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Update of residual vector                               chfoe 03/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*!
*/
void FLD::UTILS::FluidImpedanceBc::UpdateResidual(Teuchos::RCP<Epetra_Vector>  residual )
{
  residual->Update(1.0,*impedancetbc_,1.0);
  Pin_n_ = Pin_np_;
  Pc_n_  = Pc_np_;

  /*
  double norm = 0.0;
  impedancetbc_->Norm2(&norm);
  if (!myrank_)
  {
    printf("RES--_NORM: %10.5e\n",norm);
  }
  */

}



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  windkessel impedance for wave number k                  chfoe 06/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*!
determine impedance for every wave number from simple windkessel model

this method is an alternative to the impedance calculation from arterial trees
 */
std::complex<double> FLD::UTILS::FluidImpedanceBc::WindkesselImpedance(double k)
{
  double pr = R1_;  // proximal resistance
  double dr = R2_;  // distal resistance
  double ct = C_;  // capacitance

  std::complex<double> imag(0,1), imp;

  double omega = 2.0*PI*k/period_; // circular frequency

  imp = ( pr+dr + imag*omega*ct*pr*dr ) / ( 1.0 + imag*omega*ct*dr );

  return imp;
}



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  arterial tree impedance for wave number k               chfoe 05/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
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




//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  impedance w.r.t. constant flow                          chfoe 05/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
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


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  ?????????????????????????                                  ac 06/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
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





//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  ?????????????????????????                                  ac 06/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
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


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Set windkessel parameters                              ismail 02/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidImpedanceBc::SetWindkesselParams(Teuchos::ParameterList & params)
{

  R1_ = params.get<double> ("R1");
  R2_ = params.get<double> ("R2");
  C_ = params.get<double> ("C");
}// FluidImpedanceWrapper::SetWindkesselParams


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Get windkessel parameters                              ismail 02/10 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidImpedanceBc::GetWindkesselParams(Teuchos::ParameterList & params)
{

  params.set<double> ("R1",R1_);
  params.set<double> ("R2",R2_);
  params.set<double> ("C",C_);

}// FluidImpedanceWrapper::SetWindkesselParams


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Returns results of one cardiac period                  ismail 02/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidImpedanceBc::getResultsOfAPeriod(
  Teuchos::ParameterList & params,
  int             condid)
{

  std::stringstream pstream, qstream, dpstream, endCystream;

  pstream<<"pressures"<<condid;
  qstream<<"flowrates"<<condid;
  dpstream<<"dP"<<condid;
  endCystream<<"EndOfCycle"<<condid;

  params.set<Teuchos::RCP<std::vector<double> > >(pstream.str(),pressures_ );
  params.set<Teuchos::RCP<std::vector<double> > >(qstream.str(),flowrates_);
  params.set<double> (dpstream.str(),dP_);
  params.set<bool> (endCystream.str(),endOfCycle_);

}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
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
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidImpedanceBc::interpolate(Teuchos::RCP<std::vector<double> > V1,
                                               Teuchos::RCP<std::vector<double> > V2,
                                               int index1,
                                               int & index2,
                                               double time)
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
  index2 = int(double(index1)*(dt1/dt2));

}//FLD::UTILS::FluidImpedanceBc::interpolate

