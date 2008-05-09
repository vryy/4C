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
  }
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
    eleparams.set<double>("Area calculation", 0.0); // required for diameter calculation
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

 //      if (myrank_ == 0) {
//       cout << "actual flowrates: " << endl;
//       for( int i = 0; i < flowrates_->size(); i++ ) 
//       {
// 	cout << "Element " << i << " is " << (*flowrates_)[i] << endl;
//       } }
    }
    else
    {
      // we are now in the post-initial phase
//       cout << "else block met\n";
//       vector<double> ::iterator vecbegin;
//       vecbegin = flowrates_->begin();

//       flowrates_->erase(vecbegin);
//       cout << "erase done\n";
//       flowrates_->push_back(parflowrate);
//       cout << "push back done\n";

      // replace the element that was computed exactly a cycle ago
      int pos = flowratespos_ % cyclesteps_;
      (*flowrates_)[pos] = parflowrate;

      flowratespos_++;

   //    if (myrank_ == 0) {
//       cout << "actual flowrates: " << endl;;
//       for( int i = 0; i < flowrates_->size(); i++ ) 
//       {
// 	cout << "Element " << i << " is " << (*flowrates_)[i] << endl;
//       } }

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
 |  ??????????????????????????                                     ???? |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/* 

  do we have any chance of getting additional information here?        

*/
void FluidImpedanceBc::OutflowBoundary(double time, double dta, double theta)
{
  // act only if we have impedance conditions
  if ( numcondlines_ > 0 )
  {
    // calculate outflow boundary only for cycles past the first one
    if ( time > period_ )
    {

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

      // setup variables
      vector<double> values;   // impedance values, what ever this means
      values.resize(cyclesteps_);

      vector<complex<double> > frequencydomain;
      frequencydomain.resize(cyclesteps_,0);
      vector<complex<double> > timedomain;
      timedomain.resize(cyclesteps_,0);

      // size: number of generations
      std::complex<double> storage[35]={0}; // only as argument of other methods

      std::complex<double> zleft, zparent (0,0);  // only as arguments

      std::complex<double> entry;
      double PressureFromConv=0; 
      double pressuretmp;

      //Loop over all frequencies making a call to LungImpedance, w=2*pi*k/T
      //where k=-N/2,....,N/2, however k is self adjointed.
      for (int i=1; i<=cyclesteps_/2; i++)
      {
	int generation=0;
	entry = LungImpedance(i,generation,zparent,zleft,storage);
	frequencydomain[i]             = conj(entry);
	frequencydomain[cyclesteps_-i] = entry;
      }

      // calculate DC (direct current) component ...
      entry              = DCLungImpedance(0,0,zparent,storage);
      frequencydomain[0] = entry;

      // --------------------------------------------------------------------------
      // inverse Fourier transform
      // idea:
      //              
      //                    cyclesteps                i omega_k t
      //      imp(x,t) = sum          IMP(x,omega_k) e
      //                   -cyclesteps
      //
      // with omega_k = 2 PI k / T

      double trigonometrystuff = 2*PI/cyclesteps_;
      double realpart = cos(trigonometrystuff);
      double imagpart = sin(trigonometrystuff);
      std::complex<double> impcomplex (realpart,imagpart);

      for (int PeriodFraction=0; PeriodFraction<cyclesteps_; PeriodFraction++)
      {
	for (int k=0; k<cyclesteps_; k++)
	{
	  timedomain[PeriodFraction] += pow(impcomplex,-PeriodFraction*k) * frequencydomain[k];
	  timedomain[PeriodFraction] += pow(impcomplex,PeriodFraction*k) * frequencydomain[k];
	}
      }

      //Get Real component of TimeDomain, imaginary part should be anyway zero.
      for (int i=0; i<cyclesteps_; i++)
      {
	values[i]=real(timedomain[i])/cyclesteps_;
// 	if (abs(imag(timedomain[i])) > 1E-4 )
// 	  cout << "error in imaginary part is = " << timedomain[i] << endl;
      }

      // this is still somewhat mysterious
      pressuretmp=0;
      // the integral over the last period
      for (int j=0; j<=cyclesteps_; j++)
      {
	pressuretmp += values[j]*(*flowrates_)[cyclesteps_-1-j]*dta; // units: pressure x time
      }
      PressureFromConv=pressuretmp;


      eleparams.set("ConvolutedPressure",PressureFromConv);
      if (myrank_ == 0)
      {
	printf("Pressure from convolution: %f\n",PressureFromConv);
      }

      impedancetbc_->PutScalar(0.0); // ??
      const string condstring("ImpedanceCond");
      discret_->EvaluateCondition(eleparams,impedancetbc_,condstring);
      discret_->ClearState();
    } // end if ( time too small )
  } // end if ( numcondlines_ > 0 )

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
  std::complex<double> Imag (0,1), RealNumberOne (1,0), RealThreeQuarters (0.75,0), RealNumberZero (0,0), K, wave_c, StoredTmpVar, g, zright, zend, StorageEntry=0;
  double area, omega, WomersleyNumber, E=1e6, rho=1.206, nue=1.50912106e-5, compliance;
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
  omega=2*PI*ImpedanceCounter/0.8;
  //Womersley number
  WomersleyNumber=(diameter[generation]/2)*sqrt(omega/nue);
  //K depends on the womersley number
  if (WomersleyNumber > 4)
  {
  	StoredTmpVar=(2/WomersleyNumber)*Imag;
  	K=RealNumberOne+StoredTmpVar;
  }
  else
  {
  	StoredTmpVar=8/WomersleyNumber*Imag;
  	K=RealThreeQuarters+StoredTmpVar;
  }

  //Compliance
  compliance=(3*area*diameter[generation]/2)/(2*E*h[generation]);

  //Wave speed
  wave_c=sqrt(area*K/(rho*compliance));

  //Convenience coefficient
  g=sqrt(compliance*area*K/rho);

  if (real(zleft) == 0)
  {
  	zleft = (Imag)*((RealNumberOne/g)*sin(omega*length[generation]/wave_c))/(cos(omega*length[generation]/wave_c));
  	storage[generation]=zleft;
  }

  //Right side is always pre calculated!
  if (real(zleft) != 0 && generation-delta[generation] == 0)
  {
  	zright = (Imag)*((RealNumberOne/g)*sin(omega*length[generation]/wave_c))/(cos(omega*length[generation]/wave_c));
  }
  else if (real(storage[generation-delta[generation]]) == 0)
  {
  	zright = (Imag)*(RealNumberOne/g*sin(omega*length[generation]/wave_c))/(cos(omega*length[generation]/wave_c));
  }
  else
  {
  	zright = storage[generation-delta[generation]];
  }

  //Bifurcation condition
  zend = (zright*zleft)/(zleft+zright);
  //Impedance at the parent level
  zparent = (Imag)*(RealNumberOne/g*sin(omega*length[generation]/wave_c)+zend*cos(omega*length[generation]/wave_c))/(cos(omega*length[generation]/wave_c)+Imag*g*zend*sin(omega*length[generation]/wave_c));
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
  std::complex<double> Imag (0,1), RealNumberOne (1,0), RealThreeQuarters (0.75,0), RealNumberZero (0,0), K, wave_c, StoredTmpVar, g, zleft, zright, StorageEntry, zend;
  //double zend;
  int generationLimit=33;

  zleft=zparentdc;
  storage[generation]=zleft;
  int delta[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 1};
  double diameter[]={0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0002, 0.0003, 0.0003, 0.0004, 0.0005, 0.0006, 0.0008, 0.001, 0.0012, 0.0014, 0.0015, 0.0017, 0.0019, 0.0020, 0.0022, 0.0023, 0.0024, 0.0026, 0.0030, 0.0030, 0.0038, 0.0049, 0.0054, 0.0054, 0.0069, 0.0077, 0.011, 0.0121, 0.0167};
  double length[]={0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0007, 0.0009, 0.0012, 0.0015, 0.0012, 0.0028, 0.0035, 0.0041, 0.0047, 0.0054, 0.0058, 0.0071, 0.0072, 0.0087, 0.0092, 0.0093, 0.0104, 0.0090, 0.0112, 0.0097, 0.0107, 0.0122, 0.0110, 0.0128, 0.0128, 0.0119, 0.0124, 0.0249, 0.0565, 0.1130};

  if (real(zleft) == 0)
  {
	zleft = 8*length[generation]/(PI*diameter[generation]/2);
  	storage[generation]=zleft;
  }

  if (real(zleft) != 0 && generation-delta[generation] == 0)
  {
  	zright = 8*length[generation-delta[generation]]/(PI*diameter[generation-delta[generation]]/2);
  }
  else if (real(storage[generation-delta[generation]]) == 0)     //Right side is always pre calculated!
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
