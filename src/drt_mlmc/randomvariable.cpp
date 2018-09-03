/*----------------------------------------------------------------------------*/
/*!
\file randomvariable.cpp
\brief class for generating random variables

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
*/
/*-----------------------------------------------------------------------------*
 *

 *----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

#include "randomvariable.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_mlmc.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_io/io_pstream.H"


#ifdef HAVE_FFTW

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
UQ::RandomVariable::RandomVariable(
    const Teuchos::ParameterList& rfp, int random_variable_id, int seed)
    : pdf_(none),
      param_1_(0),
      param_2_(0),
      is_bounded_(false),
      lower_bound_(-10000),
      upper_bound_(100000),
      randomnumbergen_(),
      normalgen_(Teuchos::null),
      value_(0.0),
      det_value_paracont_(0.0),
      Id_(random_variable_id),
      myrank_(-1)
{
  // initialize parameters
  param_1_ = rfp.get<double>("PARAM_1");
  param_2_ = rfp.get<double>("PARAM_2");

  is_bounded_ = DRT::INPUT::IntegralValue<int>(rfp, "BOUNDED");

  lower_bound_ = rfp.get<double>("LOWERBOUND");

  upper_bound_ = rfp.get<double>("UPPERBOUND");

  det_value_paracont_ = rfp.get<double>("CONTBLENDVALUE");

  // init pdf
  INPAR::MLMC::MarginalPdf mpdf = DRT::INPUT::IntegralValue<INPAR::MLMC::MarginalPdf>(rfp, "PDF");
  switch (mpdf)
  {
    case INPAR::MLMC::pdf_gaussian:
      pdf_ = normal;
      break;
    case INPAR::MLMC::pdf_beta:
      dserror("Beta distribution not implemented yet");
      break;
    case INPAR::MLMC::pdf_lognormal:
      pdf_ = lognormal;
      break;
    default:
      dserror("Unknown pdf");
      break;
  }
  IO::cout << "Contructor of random variable" << IO::endl;
  IO::cout << "pdf_ " << pdf_ << IO::endl;

  randomnumbergen_.seed((unsigned int)seed);

  boost::normal_distribution<> normaldist(0.0, 1.0);
  normalgen_ = Teuchos::rcp(new boost::variate_generator<randnumgen&, boost::normal_distribution<>>(
      randomnumbergen_, normaldist));

  // draw random value
  value_ = (*normalgen_)();

  // get myrank without discretization
  myrank_ = DRT::Problem::Instance()->GetNPGroup()->LocalComm()->MyPID();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void UQ::RandomVariable::CreateNewSample(unsigned int seed)
{
  randomnumbergen_.seed(seed);
  // draw random number
  value_ = (*normalgen_)();
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double UQ::RandomVariable::EvalVariable(
    const double paracont_parameter, const bool writetofile, const bool output)
{
  double temp_rv_val = 0.0;
  if (pdf_ == normal)
    temp_rv_val = det_value_paracont_ * (1.0 - paracont_parameter) +
                  (paracont_parameter) * (param_1_ + param_2_ * value_);
  else if (pdf_ == lognormal)
  {
    temp_rv_val = det_value_paracont_ * (1.0 - paracont_parameter) +
                  (paracont_parameter)*exp((param_1_ + param_2_ * value_));
  }
  else
  {
    dserror("unknown marginal density ");
  }
  // quick check if the value if outside of specified bounds
  if (is_bounded_)
  {
    if (temp_rv_val > upper_bound_) temp_rv_val = upper_bound_;
    if (temp_rv_val < lower_bound_) temp_rv_val = lower_bound_;
  }
  if (writetofile)
  {
    if (myrank_ == 0)
    {
      std::stringstream filename_helper;
      filename_helper << "RandomVariable_" << Id_ << ".txt";
      std::ofstream File;
      File.open(filename_helper.str().c_str(), std::ios::app);
      // use at() to get an error massage just in case
      File << std::setprecision(9) << temp_rv_val << std::endl;
      File.close();
    }
  }
  if (output) dserror("Not implemented yet");

  physical_value_ = temp_rv_val;
  return temp_rv_val;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void UQ::RandomVariable::WriteRandomVariablesToFile(std::string filename, int numrun)
{
  const char* c = filename.c_str();
  if (myrank_ == 0)
  {
    std::ofstream File;
    File.open(c, std::ios::app);
    File << "run: " << numrun << " val: " << std::setprecision(9) << physical_value_ << std::endl;
    File.close();
  }
}

#endif
