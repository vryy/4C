
/*!
\file timeadaptivity.cpp

\brief Time step adaptivity

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#ifdef CCADISCRET

#include "timeadaptivity.H"
#include "iostream"


/*----------------------------------------------------------------------*/
/*!
\brief Constructor (public)
\author bborn
\date 10/07
*/
TimeAdaptivity::TimeAdaptivity
(
   double timeinitial,
   double timefinal,
   int timestepinitial,
   int timestepfinal,
   double stepsizeinitial,
   //
   double stepsizemax,
   double stepsizemin,
   double sizeratiomax,
   double sizeratiomin,
   double sizeratioscale,
   TAErrNorm errnorm,
   double errtol,
   int errorder,
   int adaptstepmax,
   //
   DRT::Discretization& discret,
   LINALG::Solver& solver,
   IO::DiscretizationWriter& output
)
:  timeinitial_(timeinitial),
   timefinal_(timefinal),
   timestepinitial_(timestepinitial),
   timestepfinal_(timestepfinal),
   stepsizeinitial_(stepsizeinitial),
   //
   stepsizemax_(stepsizemax),
   stepsizemin_(stepsizemin),
   sizeratiomax_(sizeratiomax),
   sizeratiomin_(sizeratiomin),
   sizeratioscale_(sizeratioscale),
   errnorm_(errnorm),
   errtol_(errtol),
   errorder_(errorder),
   adaptstepmax_(adaptstepmax),
   //
   discret_(discret),
   solver_(solver),
   output_(output),
   mypid_(discret_.Comm().MyPID())
{
   // initialise variables
   time_ = timeinitial_;
   stepsizepre_ = stepsizeinitial_;
   stepsize_ = stepsizeinitial_;
   adaptstep_ = 0;

   // get a vector layout from the discretization to construct matching
   // vectors and matrices
   if (!discret_.Filled())
   {
      discret_.FillComplete();
   }
   const Epetra_Map* dofrowmap = discret_.DofRowMap();

   // local discretisation error
   locdiserrn_ = LINALG::CreateVector(*dofrowmap, true);

   // leave for good
   return;
}  // TimeAdaptivity::TimeAdaptivity(...)


/*----------------------------------------------------------------------*/
/*!
\brief Destructor
\author bborn
\date 10/07
*/
TimeAdaptivity::~TimeAdaptivity()
{
   return;
}  // TimeAdaptivity::~TimeAdaptivity


/*----------------------------------------------------------------------*/
/*!
\brief Set parameter list
\author bborn
\date 10/07
*/
//void TimeAdaptivity::SetParaList(ParameterList& params)
//{
//   *params_ = params;
//   return;
//}

/*----------------------------------------------------------------------*/
/*!
\brief Indicate error and determine new step size
\author bborn
\date 10/07
*/
void TimeAdaptivity::Indicate
(
   bool& accepted,
   double& stpsiznew
)
{
   // norm of local discretisation error vector
   double norm;
   switch (errnorm_)
   {
   case norm_l1:
      locdiserrn_->Norm1(&norm);
      break;
   case norm_l2:
      locdiserrn_->Norm2(&norm);
      break;
   case norm_rms:
      locdiserrn_->Norm2(&norm);
      norm /= sqrt((double) (*locdiserrn_).GlobalLength());
      break;
   case norm_inf:
      locdiserrn_->NormInf(&norm);
      break;
   default:
      dserror("Cannot handle requested error norm");
   }

   // check if acceptable
   if (norm < errtol_)
   {
      accepted = true;
   }
   else
   {
      accepted = false;
   }

// debug
   cout << "ErrNorm " << norm << ", ErrTol " << errtol_ << ", Accept " << accepted << endl;

   // optimal size ration with respect to given tolerance
   double sizrat = pow(errtol_/norm, 1.0/(errorder_+1.0));
//debug
   printf("sizrat %g, stepsize %g, stepsizepre %g\n", sizrat, stepsize_, stepsizepre_);
   // scaled by safety parameter
   sizrat *= sizeratioscale_;
   // optimal new step size
   stpsiznew = sizrat * stepsize_;
   // redefine sizrat to be dt*_{n}/dt_{n-1}, ie true optimal ratio
   sizrat = stpsiznew/stepsizepre_;
   // limit sizrat by maximum and minimum
   if (sizrat > sizeratiomax_)
   {
     stpsiznew = sizeratiomax_ * stepsizepre_;
   }
   else if (sizrat < sizeratiomin_)
   {
     stpsiznew = sizeratiomin_ * stepsizepre_;
   }
   // new step size subject to safety measurements 
   if (stpsiznew > stepsizemax_)
   {
      stpsiznew = stepsizemax_;
   }
   else if (stpsiznew < stepsizemin_)
   {
      stpsiznew = stepsizemin_;
   }

   // get away from here
   return;
}


/*----------------------------------------------------------------------*/
/*!
\brief Print string for error norm
\author bborn
\date 10/07
*/
string TimeAdaptivity::PrintErrNorm() const
{
   string str;

   switch (errnorm_)
   {
   case norm_vague:
      str = "norm_vague";
      break;
   case norm_l1:
      str = "norm_l1";
      break;
   case norm_l2:
      str = "norm_l2";
      break;
   case norm_rms:
      str = "norm_rms";
      break;
   case norm_inf:
      str = "norm_inf";
      break;
   default:
      str = "norm is undefined";
      break;
   }

   return str;
}


/*----------------------------------------------------------------------*/
/*!
\brief Print constants
\author bborn
\date 10/07
*/
void TimeAdaptivity::PrintConstants(std::ostream& str) const
{
  str << "TimeAdaptivity:  Constants" << endl
      << "   Initial time = " << timeinitial_ << endl
      << "   Final time = " << timefinal_ << endl
      << "   Initial Step = " << timestepinitial_ << endl
      << "   Final Step = " << timestepfinal_ << endl
      << "   Initial step size = " << stepsizeinitial_ << endl
      << "   Max step size = " << stepsizemax_ << endl
      << "   Min step size = " << stepsizemin_ << endl
      << "   Max size ratio = " << sizeratiomax_ << endl
      << "   Min size ratio = " << sizeratiomin_ << endl
      << "   Size ratio scale = " << sizeratioscale_ << endl
      << "   Error norm = " << PrintErrNorm() << endl
      << "   Error order = " << errorder_ << endl
      << "   Error tolerance = " << errtol_ << endl
      << "   Max adaptive step = " << adaptstepmax_ << endl;
  return;
}

/*----------------------------------------------------------------------*/
/*!
\brief Print constants
\author bborn
\date 10/07
*/
void TimeAdaptivity::PrintVariables(std::ostream& str) const
{
  str << "TimeAdaptivity:  Variables" << endl
      << "   Current time = " << time_ << endl
      << "   Previous step size = " << stepsizepre_ << endl
      << "   Current step size = " << stepsize_ << endl
      << "   Current adaptive step = " << adaptstep_ << endl;
  return;
}


/*----------------------------------------------------------------------*/
/*!
\brief Print (public)
\author bborn
\date 10/07
*/
void TimeAdaptivity::Print(std::ostream& str) const
{
  str << "TimeAdaptivity" << endl;
  PrintConstants(str);
  PrintVariables(str);
  return;
}


/*======================================================================*/
/*!
\brief Out stream
\author
\date 10/07
*/
std::ostream& operator<<
(
   std::ostream& str, 
   const TimeAdaptivity::TimeAdaptivity& ta
)
{
   ta.Print(str);
   return str;
}


#endif  // #ifdef CCADISCRET
