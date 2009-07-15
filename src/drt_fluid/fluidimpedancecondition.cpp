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
#ifdef CCADISCRET

#include <stdio.h>

#include "fluidimpedancecondition.H"

#include "../drt_lib/linalg_ana.H"



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                     chfoe 06/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::UTILS::FluidImpedanceWrapper::FluidImpedanceWrapper(RefCountPtr<DRT::Discretization> actdis,
                                                         IO::DiscretizationWriter& output,
                                                         double dta) :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  output_ (output)
{
  vector<DRT::Condition*> impedancecond;
  discret_->GetCondition("ImpedanceCond",impedancecond);

  // the number of lines of impedance boundary conditions found in the input
  // note that some of these lines could belong to the same physical condition
  // which is then marked by the same 'ConditionID'
  // The corresponding resorting of result values has to be done later
  int numcondlines = impedancecond.size();

  if (numcondlines > 0) // if there is at least one impedance condition
  {
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
      RCP<FluidImpedanceBc> impedancebc = rcp(new FluidImpedanceBc(discret_, output_, dta, condid, i) );

      // -----------------------------------------------------------------
      // sort impedance bc's in map and test, if one condition ID appears
      // more than once. Currently this case is forbidden.
      // -----------------------------------------------------------------
      bool inserted = impmap_.insert( make_pair( condid, impedancebc ) ).second;
      if ( !inserted )
	dserror("There are more than one impedance condition lines with the same ID. This can not yet be handled.");
    } // end loop over condition lines from input
  } // end if there were conditions
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
 |  Wrap flow rate calculation                              chfoe 06/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidImpedanceWrapper::FlowRateCalculation(double time, double dta)
{
  // get an iterator to my map
  map<const int, RCP<class FluidImpedanceBc> >::iterator mapiter;

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
  map<const int, RCP<class FluidImpedanceBc> >::iterator mapiter;

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
 |  Wrap update of residual                                 chfoe 06/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidImpedanceWrapper::UpdateResidual(RCP<Epetra_Vector>  residual )
{
  map<const int, RCP<class FluidImpedanceBc> >::iterator mapiter;

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
  map<const int, RCP<class FluidImpedanceBc> >::iterator mapiter;

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
  map<const int, RCP<class FluidImpedanceBc> >::iterator mapiter;

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
FLD::UTILS::FluidImpedanceBc::FluidImpedanceBc(RefCountPtr<DRT::Discretization> actdis,
				   IO::DiscretizationWriter& output,
				   double dta,
				   int condid,
				   int numcond) :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  output_ (output)
{
  vector<DRT::Condition*> impedancecond;
  discret_->GetCondition("ImpedanceCond",impedancecond);

  // ---------------------------------------------------------------------
  // get time period length, steps per cycle and initialise flowratespos
  // ---------------------------------------------------------------------
  period_ = (impedancecond[numcond])->GetDouble("timeperiod");
  cyclesteps_ = (int)(period_/dta+0.5);
  flowratespos_ = 0;
  dta_ = dta;

  // ---------------------------------------------------------------------
  // get relevant data from impedance condition
  // ---------------------------------------------------------------------
  treetype_ = *((impedancecond[numcond])->Get<string>("tree"));
  termradius_ = (impedancecond[numcond])->GetDouble("termradius");
  // 'material' parameters required for artery tree
  k1_ = (impedancecond[numcond])->GetDouble("k1");
  k2_ = (impedancecond[numcond])->GetDouble("k2");
  k3_ = (impedancecond[numcond])->GetDouble("k3");

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

  flowrates_    = rcp(new vector<double>);
  impedancetbc_ = LINALG::CreateVector(*dofrowmap,true);

  // ---------------------------------------------------------------------
  // determine area of actual outlet and get material data
  // ---------------------------------------------------------------------
  double density=0.0, viscosity=0.0;
  double area = Area(density,viscosity,condid);

  // ---------------------------------------------------------------------
  // calculate impedance values and fill vector 'impvalues_'
  // ---------------------------------------------------------------------
  impvalues_.resize(cyclesteps_);
  Impedances(area,density,viscosity);

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

  std::stringstream stream1, stream2, stream3;

  stream1 << "flowratesId" << condnum;

  // write the flowrates of the previous period
  output.WriteRedundantDoubleVector(stream1.str(),flowrates_);

  // also write flowratesposition of this outlet
  stream2 << "flowratesposId" << condnum;
  output.WriteInt(stream2.str(), flowratespos_);

  // also write vector impedancetbc_ (previously missing, gee)
  stream3 << "impedancetbc" << condnum;
  output.WriteVector(stream3.str(), impedancetbc_);

  // write time step size dta_ and cyclesteps_
  output.WriteInt("ImpedanceBC_cyclesteps", cyclesteps_);
  output.WriteDouble("ImpedanceBC_dta", dta_);

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

  std::stringstream stream1, stream2, stream3;
  stream1 << "flowratesId" << condnum;
  stream2 << "flowratesposId" << condnum;

  flowratespos_ = reader.ReadInt(stream2.str());

  reader.ReadRedundantDoubleVector(flowrates_ ,stream1.str());

  // check if vector of flowrates is not empty
  if (flowrates_->size() == 0)
    dserror("could not re-read vector of flowrates");

  // also read vector impedancetbc_ (previously missing, gee)
  stream3 << "impedancetbc" << condnum;
  reader.ReadVector(impedancetbc_,stream3.str());

  // old time step size
  double odta = reader.ReadDouble("ImpedanceBC_dta");
  // new time step size
  double ndta = dta_;
  // old cyclesteps_ (no steps in a cycle)
  int oncycle = reader.ReadInt("ImpedanceBC_cyclesteps");
  // new cyclesteps_ (no steps in a cycle)
  int nncycle = cyclesteps_;


  // if the time step size has changed, need to interpolate flowrates
  // accordingly
#if 1
  if (dta_ == odta)
  {
    if (!myrank_)
      printf("Impedance restart old and new time step are the same - life is good\n");
    return;
  }
#endif

  // time of restart
  double t = reader.ReadDouble("time");
  // old number of flowrates in vector
  int onfr = (int)flowrates_->size();
  // new number of flowrates in vector
  int nnfr = (int)(onfr*odta/ndta);
  // old current position (valid for first and subsequent periods)
  //int opos = flowratespos_ % oncycle;
  int opos = flowratespos_;
  // current time within period is opos * odta
  double tinperiod = opos * odta;
  // new flowratesposition = (opos*odta)/ndta
  //int npos = (int)(tinperiod / ndta);
  int npos = (int)(t / ndta + 0.5);

  if (!myrank_)
  {
    printf("Impedance restart with time step change:\n");
    printf("time               %10.5e \n",t);
    printf("time in period     %10.5e \n",tinperiod);
    printf("old time step      %10.5e \n",odta);
    printf("new time step      %10.5e \n",ndta);
    printf("old steps in cycle        %d \n",oncycle);
    printf("new steps in cycle        %d \n",nncycle);
    printf("old length of flowrates   %d \n",onfr);
    printf("new length of flowrates   %d \n",nnfr);
    printf("old position in flowrates %d \n",opos);
    printf("new position in flowrates %d \n",npos);
    fflush(stdout);
  }

  RCP<std::vector<double> > rcpfr = rcp(new vector<double>(nnfr,0.0));
  std::vector<double>& fr = *rcpfr;

  // loop through the new time intervals
  for (int i=0; i<nnfr; ++i)
  {
    const double acttime = (i+1)*ndta; // time we interpolate to
    int j=0;
    while ((j+1)*odta<=acttime) ++j;
    const double x1 = j*odta;
    const double x2 = (j+1)*odta;
    double y1;
    if (j-1<0)
    {
      if (t <= period_) y1 = 0.0;                    // within first period starting value
      else              y1 = (*flowrates_)[onfr-1];  // subsequent periods starting value
    }
    else
      y1 = (*flowrates_)[j-1];
    double y2;
    if (j>=onfr) y2 = 0.0;
    else         y2 = (*flowrates_)[j];
    const double a = (y1-y2)/(x1-x2);
    const double b = y2 - a*x2;
    fr[i] = a * acttime + b;
    //if (myrank_==0) printf("j %4d      x1 %10.5e y1 %15.10e\n",j,x1,y1);
    //if (myrank_==0) printf("i %4d acttime %10.5e fr %15.10e\n",i,acttime,fr[i]);
    //if (myrank_==0) printf("j %4d      x2 %10.5e y2 %15.10e\n\n",j,x2,y2);
  }

  // store new values in class
  flowratespos_ = npos;
  flowrates_    = rcpfr;

  // finally, recompute the outflow boundary condition from last step
  // this way the vector need not to be stored
  OutflowBoundary(t,ndta,0.66,condnum);

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
  ParameterList eleparams;
  eleparams.set("action","area calculation");
  eleparams.set<double>("Area calculation", 0.0);
  eleparams.set<double>("viscosity", 0.0);
  eleparams.set<double>("density", 0.0);

  const string condstring("ImpedanceCond");

  discret_->EvaluateCondition(eleparams,condstring,condid);

  double actarea = eleparams.get<double>("Area calculation");
  density = eleparams.get<double>("density");
  viscosity = eleparams.get<double>("viscosity");

  // find the lowest proc number that knows the material data
  int numproc = discret_->Comm().NumProc();
  int theproc = -1;   // the lowest proc that has the desired information
  vector<double> alldens(numproc);

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
    cout << "Impedance condition Id: " << condid << " area = " << pararea << endl;
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
  vector<complex<double> > frequencydomain;  // impedances in the frequency domain
  frequencydomain.resize(cyclesteps_,0);
  vector<complex<double> > timedomain;       // impedances in the time domain
  timedomain.resize(cyclesteps_,0);

  // size: number of generations
  std::complex<double> zparent (0,0);  // only as argument

  // to store already calculated impedances (to be developed)
  map<const double, std::complex<double> > zstored;     // as argument

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
      timedomain[timefrac] += pow(eiwt,timefrac*k) * frequencydomain[k];

    // and now the conjugated ones are added
    for (int k=1; k<cyclesteps_/2; k++)
      timedomain[timefrac] += pow(eiwt,-timefrac*k) * conj(frequencydomain[k]);
  }

  //Get Real component of TimeDomain, imaginary part should be anyway zero.
  for (int i=0; i<cyclesteps_; i++)
  {
    impvalues_[i]=real(timedomain[i]);
    if (abs(imag(timedomain[i])) > 1E-4 )
      cout << "error in imaginary part is = " << timedomain[i] << endl;
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

  Calculate the flow rate accros an impedance boundary surface

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
  ParameterList eleparams;
  eleparams.set("action","flowrate calculation");
  eleparams.set<double>("Outlet flowrate", 0.0);
  eleparams.set("total time",time);

  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // get elemental flowrates ...
  RCP<Epetra_Vector> myStoredFlowrates=rcp(new Epetra_Vector(*dofrowmap,100));
  const string condstring("ImpedanceCond");
  discret_->EvaluateCondition(eleparams,myStoredFlowrates,condstring,condid);

  // ... as well as actual total flowrate on this proc
  double actflowrate = eleparams.get<double>("Outlet flowrate");

  // get total flowrate in parallel case
  double parflowrate = 0.0;
  discret_->Comm().SumAll(&actflowrate,&parflowrate,1);

  // fill vector of flowrates calculated within the last cycle
  if (time <= period_) // we are within the very first cycle
  {
    // we are now in the initial fill-in phase
    // new data is appended to our flowrates vector
    flowrates_->push_back(parflowrate);
    flowratespos_++;
  }
  else
  {
    // we are now in the post-initial phase
    // replace the element that was computed exactly a cycle ago
    int pos = flowratespos_ % cyclesteps_;
    (*flowrates_)[pos] = parflowrate;
    flowratespos_++;
  }

  if (myrank_ == 0)
  {
    cout << "Impedance condition Id: " << condid << " current Flowrate = " << parflowrate << endl;
  }

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
  // evaluate convolution integral
  double pressure=0.0;

  // the convolution integral
  for (int j=0; j<cyclesteps_; j++)
  {
    int qindex = ( flowratespos_+j ) % cyclesteps_;

    // flowrate is zero if not yet a full cycle is calculated
    double actflowrate;
    if (qindex > (int)flowrates_->size()-1)
      actflowrate = 0.0;
    else
      actflowrate = (*flowrates_)[qindex];

    int zindex = -1-j+cyclesteps_;
    pressure += impvalues_[zindex] * actflowrate * dta; // units: pressure x time
  }

  pressure = pressure/period_; // this cures the dimension; missing in Olufsen paper

  // call the element to apply the pressure
  ParameterList eleparams;
  // action for elements
  eleparams.set("action","Outlet impedance");

  eleparams.set("total time",time);
  eleparams.set("delta time",dta);
  eleparams.set("thsl",theta*dta);
  eleparams.set("ConvolutedPressure",pressure);

  if (myrank_ == 0)
    printf("Impedance condition Id: %d Pressure from convolution = %f\n",condid,pressure);


  impedancetbc_->PutScalar(0.0);
  const string condstring("ImpedanceCond");
  discret_->EvaluateCondition(eleparams,impedancetbc_,condstring,condid);
/*
  double norm = 0.0;
  impedancetbc_->Norm2(&norm);
  if (!myrank_)
  {
    printf("FLD::UTILS::FluidImpedanceBc::OutflowBoundary Norm of tbc: %10.5e\n",norm);
    fflush(stdout);
  }
*/

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
void FLD::UTILS::FluidImpedanceBc::UpdateResidual(RCP<Epetra_Vector>  residual )
{
  residual->Update(1.0,*impedancetbc_,1.0);
/*
  double norm = 0.0;
  impedancetbc_->Norm2(&norm);
  if (!myrank_)
  {
    printf("FLD::UTILS::FluidImpedanceBc::UpdateResidual Norm of tbc: %10.5e\n",norm);
    fflush(stdout);
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
  double pr = k1_;  // proximal resistance
  double dr = k2_;  // distal resistance
  double ct = k3_;  // capacitance

  complex<double> imag(0,1), imp;

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
						       map<const double,complex<double> > zstored)
{
  // general data
  double lscale = 50.0; // length to radius ratio
  double alpha = 0.50;   // right daughter vessel ratio
  double beta = 0.85;    // left daughter vessel ratio

  // some auxiliary stuff
  complex<double> koeff, imag(0,1), cwave;

  // terminal resistance is assumed zero
  complex<double> zterminal (0,0);

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
  complex<double> zleft;
  complex<double> zright;
  bool terminated = false;

  // only if both vessels are smaller than the limit truncate
  if (leftradius < termradius && rightradius < termradius)
    terminated = true;
  else
  {
    map<const double,complex<double> >::iterator iter = zstored.find(leftradius);
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
  complex<double> zdown;
  if (terminated)
    zdown = zterminal;
  else
    zdown = 1.0 / (1.0/zleft + 1.0/zright);


  // ... and compute impedance at my upstream end!
  //*************************************************************
  double compliance = 1.5*area / ( k1_ * exp(k2_*radius) + k3_ );

  double sqrdwo = radius*radius*omega/viscosity;  // square of Womersley number
  double wonu = sqrt(sqrdwo);                     // Womersley number itself

  if (wonu > 4.0)
    koeff = 1.0 - (2.0/wonu)/sqrt(imag);
  else
    koeff = 1.0 / ( 1.333333333333333333 - 8.0*imag/ (wonu*wonu) );

  // wave speed of this frequency in present vessel
  cwave=sqrt( area*koeff / (density*compliance) );

  //Convenience coefficient
  complex<double> gcoeff = compliance * cwave;

  // calculate impedance of this, the present vessel
  complex<double> argument = omega*length/cwave;
  complex<double> zparent  = (imag/gcoeff * sin(argument) + zdown*cos(argument) ) /
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
							 map<const double,complex<double> > zstored)
{
  // general data
  double lscale = 50.0; // length to radius ratio
  double alpha = 0.50;   // right daughter vessel ratio
  double beta = 0.85;    // left daughter vessel ratio
  double mu = viscosity * density; // dynamic (physical) viscosity

  // terminal resistance is assumed zero
  complex<double> zterminal (0,0);


  // get impedances of downward vessels ...
  //*****************************************
  generation++;  // this is the next generation

  double leftradius  = alpha*radius;
  double rightradius = beta*radius;
  complex<double> zleft;
  complex<double> zright;
  bool terminated = false;

  // only if both vessels are smaller than the limit truncate
  if (leftradius < termradius && rightradius < termradius)
    terminated = true;
  else
  {
    map<const double,complex<double> >::iterator iter = zstored.find(leftradius);
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
  complex<double> zdown;
  if (terminated)
    zdown = zterminal;
  else
    zdown = 1.0 / (1.0/zleft + 1.0/zright);


  // ... and compute dc impedance at my upstream end!
  //*************************************************************

  // calculate dc impedance of this, the present, vessel
  complex<double> zparentdc = 8.0 * mu * lscale / ( PI*radius*radius*radius ) + zdown;
  // DEBUG output
  //cout << "generation: " << generation << endl;
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
						     map<const double,complex<double> > zstored)
{
   // general data
  double lscale = 5.8; // length to radius ratio
  double alpha = 0.876;   // right daughter vessel ratio
  double beta = 0.686;    // left daughter vessel ratio

  // some auxiliary stuff
  complex<double> imag(0,1), Z1, Z2, Z3, ZW;
  complex<double> koeff, cwave;
  // terminal resistance is assumed zero
  complex<double> zterminal (0,0);

  double omega = 2.0*PI*k/period_;

  // Alveolar impedance
  //double G=0.96e-2;
  //double H=8e-2;
  //double athing=(2.0/PI)*atan(H/G);
  //complex<double> zterminal=(G-imag*H)/(pow(omega,athing));

  // build up geometry of present generation
  double area = radius*radius*PI;
  double length = lscale * radius;

  double h=-0.0057*radius*radius+0.2096*radius+0.0904;
  double E=0.0033;
  double compliance=2.0*PI*radius*radius*radius/(3.0*h*E);

  // get impedances of downward vessels ...
  //*****************************************
  generation++;  // this is the next generation
  //cout<<"generation "<<generation<<endl;
  // left hand side:
    double leftradius  = alpha*radius;
    double rightradius = beta*radius;
    complex<double> zleft;
    complex<double> zright;
    bool terminated = false;

    if (leftradius < termradius && rightradius < termradius)
            terminated = true;
          else
          {
            map<const double,complex<double> >::iterator iter = zstored.find(leftradius);
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
    complex<double> zdown;
    if (terminated)
      zdown = zterminal;
    else
      zdown = 1.0 / (1.0/zleft + 1.0/zright);

  // ... and compute impedance at my upstream end!
  //*************************************************************

    double sqrdwo = radius*radius*omega/viscosity;  // square of Womersley number
    double wonu = sqrt(sqrdwo);                     // Womersley number itself

    if (wonu > 4.0)
    koeff = 1.0 - (2.0/wonu)/sqrt(imag);
    else
    koeff = 1.0 / ( 1.333333333333333333 - 8.0*imag/ (wonu*wonu) );

    // wave speed of this frequency in present vessel
    cwave=sqrt( area*koeff / (density*compliance) );

    //Convenience coefficient
    complex<double> gcoeff = compliance * cwave;

    // calculate impedance of this, the present vessel
    complex<double> argument = omega*length/cwave;
    complex<double> zparent  = (imag/gcoeff * sin(argument) + zdown*cos(argument) ) /
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
								 map<const double,complex<double> > zstored)
{
  //general data
  double lscale = 5.8; // length to radius ratio
  double alpha = 0.876;   // right daughter vessel ratio
  double beta = 0.686;    // left daughter vessel ratio

  double mu = viscosity * density; // dynamic (physical) viscosity
  generation++;
  // terminal resistance is assumed zero
  complex<double> zterminal (0,0);

    double leftradius  = alpha*radius;
    double rightradius = beta*radius;
    complex<double> zleft;
    complex<double> zright;
    bool terminated = false;

    if (leftradius < termradius && rightradius < termradius)
           terminated = true;
         else
         {
           map<const double,complex<double> >::iterator iter = zstored.find(leftradius);
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
    complex<double> zdown;
    if (terminated)
      zdown = zterminal;
    else
      zdown = 1.0 / (1.0/zleft + 1.0/zright);


    // ... and compute dc impedance at my upstream end!
    //*************************************************************

    // calculate dc impedance of this, the present vessel
    complex<double>  zparentdc = 8.0 * mu * lscale / ( PI*radius*radius*radius ) + zdown;
    return zparentdc;
}//FluidImplicitTimeInt::DCLungImpedance





#endif /* CCADISCRET       */
