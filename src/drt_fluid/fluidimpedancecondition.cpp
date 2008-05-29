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
  std::complex<double> storage[35]={0}; // only as argument of other methods
  std::complex<double> zleft, zparent (0,0);  // only as argument
  std::complex<double> entry;

  // **************** apparently this does not work as it should ... **********

//       //Loop over all frequencies making a call to LungImpedance, w=2*pi*k/T
//       //where k=-N/2,....,N/2, however k is self adjointed.
//       for (int k=1; k<=cyclesteps_/2; k++)
//       {
// 	int generation=0;
// 	entry = LungImpedance(k,generation,zparent,zleft,storage);
// 	frequencydomain[k]             = conj(entry);
// 	frequencydomain[cyclesteps_-k] = entry;
//       }

//       // calculate DC (direct current) component ...
//       entry              = DCLungImpedance(0,0,zparent,storage);
//       frequencydomain[0] = entry;

//       // --------------------------------------------------------------------------
//       // inverse Fourier transform
//       // idea:
//       //              
//       //                    cyclesteps                i omega_k t
//       //      imp(x,t) = sum          IMP(x,omega_k) e
//       //                   -cyclesteps
//       //
//       // with omega_k = 2 PI k / T

//       double trigonometrystuff = 2*PI/cyclesteps_;
//       double realpart = cos(trigonometrystuff);
//       double imagpart = sin(trigonometrystuff);
//       std::complex<double> impcomplex (realpart,imagpart);

//       for (int PeriodFraction=0; PeriodFraction<cyclesteps_; PeriodFraction++)
//       {
// 	for (int k=0; k<cyclesteps_; k++)
// 	{
// 	  timedomain[PeriodFraction] += pow(impcomplex,-PeriodFraction*k) * frequencydomain[k];
// 	  timedomain[PeriodFraction] += pow(impcomplex,PeriodFraction*k) * frequencydomain[k];
// 	}
//       }

  // **************** ... therefore an alternative approach ****************

  //Loop over some frequencies making a call to Impedance, w=2*pi*k/T
  //where k=-N/2,....,N/2, 

  // calculate DC (direct current) component ...
  frequencydomain[0] = DCLungImpedance(0,0,zparent,storage);
  //frequencydomain[0] = DCArteryImpedance(0,0,area,density,viscosity,zparent,storage);

  // this results in a field like
  // frequencydomain = 
  // [ Z(w_0)  Z(w_1)  Z(w_2)  Z(w_3)  ...  Z(w_cyclesteps) ]
  //
  // note that we also need the complex conjugated ones, i.e.
  // Z(-w_1) = conj(Z(w_1)  to Z(-w_cyclesteps) = conj(Z(w_cyclesteps)
  for (int k=1; k<cyclesteps_; k++)
  {
    int generation=0;
    frequencydomain[k] = LungImpedance(k,generation,zparent,zleft,storage);
    //frequencydomain[k] = ArteryImpedance(k,generation,area,density,viscosity,zparent,zleft,storage);
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
  //  vector<double> values;   // real impedance values
  //  values.resize(cyclesteps_);

  for (int i=0; i<cyclesteps_; i++)
  {
    //	values[i]=real(timedomain[i])/cyclesteps_;  // WHY division by cyclesteps?
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

  //     cout << "flowrates:" << endl;
//       for (int i=0; i<flowrates_->size(); i++)
//       {
// 	cout << "[" << i <<"]" << (*flowrates_)[i] << endl;
//       }
    }
    else
    {
      // we are now in the post-initial phase
      // replace the element that was computed exactly a cycle ago
    	flowrates_->erase(flowrates_->begin());
    	flowrates_->push_back(parflowrate);

      //int pos = flowratespos_ % cyclesteps_;
      //(*flowrates_)[pos] = parflowrate;

      //flowratespos_++;
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
	int zindex = ( flowratespos_-1-j+cyclesteps_ ) % cyclesteps_;
	int qindex = ( flowratespos_+j ) % cyclesteps_;

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
 |  ?????????????????????????                               chfoe 05/08 |
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
						       double area,
						       double density,
						       double viscosity,
						       std::complex<double> zparent,
						       std::complex<double> zleft,
						       std::complex<double> storage[])
{
  // general data
  double lscale = 50.0; // length to radius ratio
  double alpha = 0.9;   // right daughter vessel ratio
  double beta = 0.6;    // left daughter vessel ratio

  double omega = 2.0*PI*k/period_;

  double k1 = 2.0;    // g/ms/ms/mm
  double k2 = -2.253; // per mm
  double k3 = 0.0865; // g/ms/ms/mm

  // build up geometry of root generation
  double rootrad = sqrt(area/PI); // an estimate for the root radius of the tree
  double rootlength = lscale* rootrad;

  // just as a test: root impedance only

  double compliance = 1.5*area / ( k1 * exp(k2*rootrad) + k3 );

  complex<double> koeff, imag(0,1), cwave;

  double sqrdwo = rootrad*rootrad*omega/viscosity;
//   complex<double> aux (0,-imgpart);
//   complex<double> wonu = sqrt(aux);

  double wonu = sqrt(sqrdwo);

  if (wonu > 4.0)
  {
    complex<double> realone(1,0);
    koeff = realone - (2.0/wonu)*sqrt(imag);
  }
  else
  {
    complex<double> number(1.333333333333333333,0);
    koeff = 1.0 / ( number - 8.0*imag/ (wonu*wonu) );
  }

  // wave speed of this frequency in present vessel
  cwave=sqrt( area*koeff / (density*compliance) );

  //Convenience coefficient
  complex<double> gcoeff = compliance * cwave;

  // terminal resistance is assumed zero
  complex<double> zterminal (0,0);

  // calculate impedance of this, the parent, vessel
  complex<double> argument = omega*rootlength/cwave;
  zparent = (imag/gcoeff * sin(argument) + zterminal*cos(argument) ) /
            ( cos(argument) + imag*gcoeff*zterminal*sin(argument) );

  //  cout << "Impedance : " << zparent << endl;

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
std::complex<double> FluidImpedanceBc::DCArteryImpedance(int ImpedanceCounter, 
						         int generation,  
							 double area,
							 double density,
							 double viscosity,
						         std::complex<double> zparentdc,
						         std::complex<double> storage[])
{
  // general data
  double lscale = 50.0; // length to radius ratio
  double mu = viscosity * density; // dynamic (physical) viscosity

  // build up geometry of root generation
  double rootrad = sqrt(area/PI); // an estimate for the root radius of the tree

  zparentdc = 8.0 * mu * lscale / ( PI*rootrad*rootrad*rootrad );

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
std::complex<double> FluidImpedanceBc::LungImpedance(int ImpedanceCounter, 
						     int generation,
						     std::complex<double> zparent,
						     std::complex<double> zleft,
						     std::complex<double> storage[])
{
  std::complex<double> imag (0,1), realone (1,0), koeff, wave_c, StoredTmpVar, g, zright, zend, StorageEntry=0;
  double area, omega, E=1e6, rho=1.206, nue=1.50912106e-5, compliance;
  int generationLimit=33;
  storage[generation]=zleft;

  //Lung morphological data
  int delta[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 1};
  double diameter[]={0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0002, 0.0003, 0.0003, 0.0004, 0.0005, 0.0006, 0.0008, 0.001, 0.0012, 0.0014, 0.0015, 0.0017, 0.0019, 0.0020, 0.0022, 0.0023, 0.0024, 0.0026, 0.0030, 0.0030, 0.0038, 0.0049, 0.0054, 0.0054, 0.0069, 0.0077, 0.011, 0.0121, 0.0167};
  double length[]={0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0007, 0.0009, 0.0012, 0.0015, 0.0012, 0.0028, 0.0035, 0.0041, 0.0047, 0.0054, 0.0058, 0.0071, 0.0072, 0.0087, 0.0092, 0.0093, 0.0104, 0.0090, 0.0112, 0.0097, 0.0107, 0.0122, 0.0110, 0.0128, 0.0128, 0.0119, 0.0124, 0.0249, 0.0565, 0.1130};
  double h[]={0.00014, 0.00014, 0.00014, 0.00014, 0.00014, 0.00014, 0.00010, 0.00010, 0.00011, 0.00012, 0.00013, 0.00015, 0.00017, 0.00020, 0.00022, 0.00024, 0.00026, 0.00028, 0.00030, 0.00030, 0.00032, 0.00033, 0.00034, 0.00026, 0.00040, 0.00040, 0.00047, 0.00057, 0.00062, 0.00062, 0.00074, 0.00080, 0.00105, 0.00113, 0.00146};

  //Area of vessel
  area=(PI*diameter[generation])/4;
  //Frequency
  omega=2*PI*ImpedanceCounter/period_;
  //Womersley number
  double sqrdwo=(diameter[generation]/2)*(diameter[generation]/2)*omega/nue;
  //K depends on the womersley number
  
  double wonu = sqrt(sqrdwo);

    if (wonu > 4.0)
    {
      complex<double> realone(1,0);
      koeff = realone - (2.0/(wonu*sqrt(imag)));
    }
    else
    {
      complex<double> number(1.333333333333333333,0);
      koeff = 1.0 / ( number - 8.0*imag/ (wonu*wonu) );
    }

  //Compliance
  compliance=(3*area*diameter[generation]/2)/(2*E*h[generation]);

  //Wave speed
  wave_c=sqrt(area*koeff/(rho*compliance));

  //Convenience coefficient
  g=sqrt(compliance*area*koeff/rho);

  
  
  if (abs(zleft) == 0)
  {
  	zleft = (imag)*((realone/g)*sin(omega*length[generation]/wave_c))/(cos(omega*length[generation]/wave_c));
  	storage[generation]=zleft;
  }

  //Right side is always pre calculated!
  if ( abs(zleft) != 0 && generation-delta[generation] == 0)
  {
  	zright = (imag)*((realone/g)*sin(omega*length[generation]/wave_c))/(cos(omega*length[generation]/wave_c));
  }
  else if (abs(storage[generation-delta[generation]]) == 0)
  {
  	zright = (imag)*(realone/g*sin(omega*length[generation]/wave_c))/(cos(omega*length[generation]/wave_c));
  }
  else
  {
  	zright = storage[generation-delta[generation]];
  }

  //Bifurcation condition
  zend = (zright*zleft)/(zleft+zright);
  //Impedance at the parent level
  zparent = (imag)*(realone/g*sin(omega*length[generation]/wave_c)+zend*cos(omega*length[generation]/wave_c))/(cos(omega*length[generation]/wave_c)+imag*g*zend*sin(omega*length[generation]/wave_c));
  zleft=zparent;
  //Right side is always pre calculated!
  if (generation < generationLimit)
  {
  	generation++;
  	LungImpedance(ImpedanceCounter, generation, zparent, zleft, storage);
  }
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
std::complex<double> FluidImpedanceBc::DCLungImpedance(int ImpedanceCounter, 
						       int generation, 
						       std::complex<double> zparentdc,
						       std::complex<double> storage[])
{
  //DC impedance
  std::complex<double> zleft, zright, StorageEntry, zend;
  //double zend;
  int generationLimit=33;

  zleft=zparentdc;
  storage[generation]=zleft;
  int delta[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 1};
  double diameter[]={0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0002, 0.0003, 0.0003, 0.0004, 0.0005, 0.0006, 0.0008, 0.001, 0.0012, 0.0014, 0.0015, 0.0017, 0.0019, 0.0020, 0.0022, 0.0023, 0.0024, 0.0026, 0.0030, 0.0030, 0.0038, 0.0049, 0.0054, 0.0054, 0.0069, 0.0077, 0.011, 0.0121, 0.0167};
  double length[]={0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0007, 0.0009, 0.0012, 0.0015, 0.0012, 0.0028, 0.0035, 0.0041, 0.0047, 0.0054, 0.0058, 0.0071, 0.0072, 0.0087, 0.0092, 0.0093, 0.0104, 0.0090, 0.0112, 0.0097, 0.0107, 0.0122, 0.0110, 0.0128, 0.0128, 0.0119, 0.0124, 0.0249, 0.0565, 0.1130};

  if ( abs(zleft) == 0)
  {
	zleft = 8*length[generation]/(PI*diameter[generation]/2);
  	storage[generation]=zleft;
  }

  if ( abs(zleft) != 0 && generation-delta[generation] == 0)
  {
  	zright = 8*length[generation-delta[generation]]/(PI*diameter[generation-delta[generation]]/2);
  }
  else if ( abs(storage[generation-delta[generation]]) == 0)     //Right side is always pre calculated!
  {
  	zright = 8*length[generation-delta[generation]]/(PI*diameter[generation-delta[generation]]/2);
  }
  else
  {
  	zright = storage[generation-delta[generation]];
  }

  //Bifurcation condition
  zend = (zright*zleft)/(zleft+zright);
  //Impedance at the parent level
  zparentdc = zend;

  if (generation < generationLimit) //Number of Generations to model, better for small vessels ie. not the trachea
  {
  	generation++;
  	DCLungImpedance(ImpedanceCounter, generation, zparentdc, storage);
  }
  //zparentdc= 8*length[generation]/(PI*diameter[generation]/2)+zparentdc;
  return StorageEntry=zparentdc;
}//FluidImplicitTimeInt::DCLungImpedance


#endif /* CCADISCRET       */
