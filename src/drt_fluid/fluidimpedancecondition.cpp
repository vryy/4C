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
 |  Constructor (public)                                     chfoe 04/08|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FluidImpedanceBc::FluidImpedanceBc(RefCountPtr<DRT::Discretization> actdis,
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
  numcondlines_ = impedancecond.size();

  // debugging output
  //  cout << *impedancecond[0] << endl;

  if (numcondlines_ > 0) // if there is at least one impedance condition
  {
    // -------------------------------------------------------------------
    // now care for the fact that there could be more than one input line
    // belonging to the same impedance boundary condition
    // -------------------------------------------------------------------
    // still to come ... will cause restructuring of the remainder ...

    // -------------------------------------------------------------------
    // get relevant data from impedance condition
    // -------------------------------------------------------------------
    condid_ = (impedancecond[0])->Getint("ConditionID");
    period_ = (impedancecond[0])->GetDouble("timeperiod");
    treetype_ = *((impedancecond[0])->Get<string>("tree"));
    // 'material' parameters required for artery tree
    k1_ = (impedancecond[0])->GetDouble("k1");
    k2_ = (impedancecond[0])->GetDouble("k2");
    k3_ = (impedancecond[0])->GetDouble("k3");

    // -------------------------------------------------------------------
    // test that we have an integer number of time steps per cycle
    // -------------------------------------------------------------------
    // something more intelligent could possibly be found one day ...
    double doublestepnum = period_/dta;

    cyclesteps_ = (int)(doublestepnum+0.5);
    double diff = doublestepnum - (double)cyclesteps_;

    if ( abs(diff) > 1.0E-5 )
      dserror("Make sure that the cycle can be calculated within an integer number of steps!!!");

    flowratespos_ = 0;

    // -------------------------------------------------------------------
    // get the processor ID from the communicator
    // -------------------------------------------------------------------
    myrank_  = discret_->Comm().MyPID();

    // -------------------------------------------------------------------
    // get a vector layout from the discretization to construct matching
    // vectors and matrices
    //                 local <-> global dof numbering
    // -------------------------------------------------------------------
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    // history data
    hist_         = LINALG::CreateVector(*dofrowmap,true);

    flowrates_    = rcp(new vector<double>);
    impedancetbc_ = LINALG::CreateVector(*dofrowmap,true);

    // -------------------------------------------------------------------
    // determine area of actual outlet and get material data
    // -------------------------------------------------------------------
    double density=0.0, viscosity=0.0;
    double area = Area(density,viscosity);

    // -------------------------------------------------------------------
    // calculate impedance values and fill vector 'impvalues_'
    // -------------------------------------------------------------------
    impvalues_.resize(cyclesteps_);
    Impedances(area,density,viscosity);
  } // end if numcondlines_ > 0
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Destructor dtor (public)                                chfoe 04/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FluidImpedanceBc::~FluidImpedanceBc()
{
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
void FluidImpedanceBc::WriteRestart( IO::DiscretizationWriter&  output )
{
  if ( numcondlines_ > 0 )
  {
    // write the flowrates of the previous period
    const string outstring1 = "flowrates";
    output.WriteRedundantDoubleVector(outstring1,flowrates_);

    // also write
    const string outstring2 = "flowratespos";
    output.WriteInt(outstring2, flowratespos_);
  }
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
void FluidImpedanceBc::ReadRestart( IO::DiscretizationReader& reader )
{
  if ( numcondlines_ > 0 )
  {
    flowratespos_ = reader.ReadInt("flowratespos");

    reader.ReadRedundantDoubleVector(flowrates_ ,"flowrates");

    if (flowrates_->size() == 0)
      dserror("could not re-read vector of flowrates");

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
double FluidImpedanceBc::Area( double& density, double& viscosity )
{
  // fill in parameter list for subsequent element evaluation
  // there's no assembly required here
  ParameterList eleparams;
  eleparams.set("action","area calculation");
  eleparams.set<double>("Area calculation", 0.0);
  eleparams.set<double>("viscosity", 0.0);
  eleparams.set<double>("density", 0.0);
  eleparams.set("assemble matrix 1",false);
  eleparams.set("assemble matrix 2",false);
  eleparams.set("assemble vector 1",false);
  eleparams.set("assemble vector 2",false);
  eleparams.set("assemble vector 3",false);

  const string condstring("ImpedanceCond");
  discret_->EvaluateCondition(eleparams,condstring);

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
    cout << "This Area: " << pararea << endl;
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
void FluidImpedanceBc::Impedances( double area, double density, double viscosity )
{
  // setup variables
  vector<complex<double> > frequencydomain;  // impedances in the frequency domain
  frequencydomain.resize(cyclesteps_,0);
  vector<complex<double> > timedomain;       // impedances in the time domain
  timedomain.resize(cyclesteps_,0);

  // size: number of generations
  std::complex<double> zparent (0,0);  // only as argument

  // to store already calculated impedances (to be developed)
  map<const double*, std::complex<double> > zstored;

  //Loop over some frequencies making a call to Impedance, w=2*pi*k/T
  //where k=-N/2,....,N/2,

  double radius = sqrt(area/PI);   // the radius to which the artery tree is connected
  //double radius = 3.0;
  //  double termradius = radius/5000; // the radius at which the artery tree is terminated
  double termradius = 0.5; // the radius at which the artery tree is terminated

  // calculate DC (direct current) component depending on type of tree
  if ( treetype_ == "lung" )
    frequencydomain[0] = DCLungImpedance(0,radius,termradius,density,viscosity,zparent);

  if ( treetype_ == "artery" )
    frequencydomain[0] = DCArteryImpedance(0,radius,termradius,density,viscosity,zparent);

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
      frequencydomain[k] = LungImpedance(k,generation,radius,termradius,density,viscosity,zparent);

    if ( treetype_ == "artery" )
      frequencydomain[k] = ArteryImpedance(k,generation,radius,termradius,density,viscosity,zparent);
  }

  for (int k=0; k<cyclesteps_; k++)
    cout << "wavenumber k = " << k << "  |Z| = " << abs(frequencydomain[k]) << endl;

  cout << endl << endl;

  for (int k=0; k<cyclesteps_; k++)
    cout << k*2.0*PI << "   " << abs(frequencydomain[k]) << endl;

//   ArteryImpedance(3,0,radius,termradius,density,viscosity,zparent);

//   dserror("erst einmal Schluß");

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
    for (int k=0; k<cyclesteps_; k++)
      timedomain[timefrac] += pow(eiwt,timefrac*k) * frequencydomain[k];

    // and now the conjugated ones are added
    for (int k=1; k<cyclesteps_; k++)
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
  modifyed by chfoe 04/08

  Calculate the flow rate accros an impedance boundary surface

  Flow rates are
  (1) calculated for single element surfaces
  (2) added up over the elements of the single procs
  (3) communicated and added over the procs
  (4) and finally stored within the vector 'flowrates_'

  The vector of the flowrates holds the flow rate history of the
  very last cycle!

*/
void FluidImpedanceBc::FlowRateCalculation(double time, double dta)
{
  // act only, if there is some impedance condition
  if ( numcondlines_ > 0 )
  {
    // fill in parameter list for subsequent element evaluation
    // there's no assembly required here
    ParameterList eleparams;
    eleparams.set("action","flowrate calculation");
    eleparams.set<double>("Outlet flowrate", 0.0);
    eleparams.set("total time",time);
    eleparams.set("assemble matrix 1",false);
    eleparams.set("assemble matrix 2",false);
    eleparams.set("assemble vector 1",false);
    eleparams.set("assemble vector 2",false);
    eleparams.set("assemble vector 3",false);

    // get a vector layout from the discretization to construct matching
    // vectors and matrices
    //                 local <-> global dof numbering
    const Epetra_Map* dofrowmap = discret_->DofRowMap();

    // get elemental flowrates ...
    RCP<Epetra_Vector> myStoredFlowrates=rcp(new Epetra_Vector(*dofrowmap,100));
    const string condstring("ImpedanceCond");
    discret_->EvaluateCondition(eleparams,myStoredFlowrates,condstring);

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
      cout << "Current Flowrate: " << parflowrate << endl;
    }
  } // end if ( numcondlines_ > 0 )
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
void FluidImpedanceBc::OutflowBoundary(double time, double dta, double theta)
{
  // act only if we have impedance conditions
  if ( numcondlines_ > 0 )
  {
    // calculate outflow boundary only for cycles past the first one
    if ( time > period_ )
    {
      // evaluate convolution integral
      double pressure=0.0;

      // the convolution integral
      for (int j=0; j<cyclesteps_; j++)
      {
	int qindex = ( flowratespos_+j ) % cyclesteps_;
	int zindex = -1-j+cyclesteps_;

	pressure += impvalues_[zindex] * (*flowrates_)[qindex] * dta; // units: pressure x time
      }

      pressure = pressure/period_; // this cures the dimension; missing in Olufsen paper

      // call the element to apply the pressure
      ParameterList eleparams;
      // action for elements
      eleparams.set("action","Outlet impedance");

      //Only assemble a single vector
      eleparams.set("assemble matrix 1",false);
      eleparams.set("assemble matrix 2",false);
      eleparams.set("assemble vector 1",true);
      eleparams.set("assemble vector 2",false);
      eleparams.set("assemble vector 3",false);

      eleparams.set("total time",time);
      eleparams.set("delta time",dta);
      eleparams.set("thsl",theta*dta);
      eleparams.set("ConvolutedPressure",pressure);

      if (myrank_ == 0)
      {
	printf("Pressure from convolution: %f\n",pressure);
      }

      impedancetbc_->PutScalar(0.0); // ??
      const string condstring("ImpedanceCond");
      discret_->EvaluateCondition(eleparams,impedancetbc_,condstring);
      discret_->ClearState();

    } // end if ( time too small )
  } // end if ( numcondlines_ > 0 )
  return;
}//FluidImplicitTimeInt::OutflowBoundary


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
void FluidImpedanceBc::UpdateResidual(RCP<Epetra_Vector>  residual )
{
  if ( numcondlines_ > 0 )
  {
    residual->Update(1.0,*impedancetbc_,1.0);
  }
  return;
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
build up artery tree and calculate root impedance for a given frequency

Data taken from
  Olufsen et al.: "Numerical Simulation and Experimental Validation of
  Blood Flow in Arteries with Structured-Tree Outflow Conditions",
  Annals of Biomedical Eingineering, Vol. 28, pp. 1281--1299, 2000.

Further also needed
  Olufsen, M.S.: "Structured tree outflow condition for blood flow
  in larger systematic arteries", American Physiological Society, 1999.

*/
std::complex<double> FluidImpedanceBc::ArteryImpedance(int k,
						       int generation,
						       double radius,
						       double termradius,
						       double density,
						       double viscosity,
						       std::complex<double> zparent)
{
  // general data
  double lscale = 50.0; // length to radius ratio
  double alpha = 0.9;   // right daughter vessel ratio
  double beta = 0.6;    // left daughter vessel ratio

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
    zleft  = ArteryImpedance(k,generation,leftradius,termradius,density,viscosity,zparent);
    zright = ArteryImpedance(k,generation,rightradius,termradius,density,viscosity,zparent);
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
  zparent = (imag/gcoeff * sin(argument) + zdown*cos(argument) ) /
            ( cos(argument) + imag*gcoeff*zdown*sin(argument) );

//   cout << "Generation: " << generation << endl;
//   cout << " argument = " << argument << "   g = " << gcoeff << endl;
//   cout << " Zup = " << abs(zparent) << "   Zdown = " << abs(zdown) << "   Zleft = " << abs(zleft) << "   Zright = " << abs(zright) << endl;

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
std::complex<double> FluidImpedanceBc::DCArteryImpedance(int generation,
							 double radius,  
							 double termradius,
							 double density,
							 double viscosity,
						         std::complex<double> zparentdc)
{
  // general data
  double lscale = 50.0; // length to radius ratio
  double alpha = 0.9;   // right daughter vessel ratio
  double beta = 0.6;    // left daughter vessel ratio
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
    zleft  = DCArteryImpedance(generation,leftradius,termradius,density,viscosity,zparentdc);
    zright = DCArteryImpedance(generation,rightradius,termradius,density,viscosity,zparentdc);
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
  zparentdc = 8.0 * mu * lscale / ( PI*radius*radius*radius ) + zdown;
  cout << "generation: " << generation << endl;
  return zparentdc;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  ?????????????????????????                               ????? 04/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*! 
What is supposed to happen within these lines? 
*/
std::complex<double> FluidImpedanceBc::LungImpedance(int k, 
						       int generation,
						       double radius,
						       double termradius,
							   double density,
							   double viscosity,
						       std::complex<double> zparent)
{
   // general data
  double lscale = 5.8; // length to radius ratio
  double alpha = 0.876;   // right daughter vessel ratio
  double beta = 0.686;    // left daughter vessel ratio
  double h=0.2;
  // some auxiliary stuff
  complex<double> koeff, imag(0,1), cwave;
  // terminal resistance is assumed zero ff
  complex<double> zterminal (0,0);
  
  double omega = 2.0*PI*k/period_;

  // this has to be moved!!
  double E=0.0033;

  // build up geometry of present generation
  double area = radius*radius*PI;
  double length = lscale * radius;


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

    // only if both vessels are smaller than the limit truncate
    if (leftradius < termradius && rightradius < termradius)
      terminated = true;
    else
    {
      zleft  = LungImpedance(k,generation,leftradius,termradius,density,viscosity,zparent);
      zright = LungImpedance(k,generation,rightradius,termradius,density,viscosity,zparent);
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
  double compliance = (4.0*E*h)/(3*radius);
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
  zparent = (imag/gcoeff * sin(argument) + zdown*cos(argument) ) /
            ( cos(argument) + imag*gcoeff*zdown*sin(argument) );

  return zparent;
} //BioFluidImplicitTimeInt::LungImpedance





//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  ?????????????????????????                               ????? 05/08 |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*! 
What is supposed to happen within these lines? 
*/
std::complex<double> FluidImpedanceBc::DCLungImpedance(int generation,
								 double radius,
								 double termradius,
								 double density,
								 double viscosity,
						         std::complex<double> zparentdc)
{
// general data
  double lscale = 5.8; // length to radius ratio
  double alpha = 0.876;   // right daughter vessel ratio
  double beta = 0.686;    // left daughter vessel ratio

  double mu = viscosity * density; // dynamic (physical) viscosity

  // terminal resistance is assumed zero
  complex<double> zterminal (0,0);

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
      zleft  = DCLungImpedance(generation,leftradius,termradius,density,viscosity,zparentdc);
      zright = DCLungImpedance(generation,rightradius,termradius,density,viscosity,zparentdc);
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
    zparentdc = 8.0 * mu * lscale / ( PI*radius*radius*radius ) + zdown;
    return zparentdc;
}//FluidImplicitTimeInt::DCLungImpedance





#endif /* CCADISCRET       */
