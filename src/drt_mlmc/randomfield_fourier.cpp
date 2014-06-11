/*!----------------------------------------------------------------------
\file randomfield_fourier.cpp
Created on: 15 May, 2014
\brief
Class for generating samples of gaussian and lognormal random fields with
a squared exponential covariance structure based on a fourier series expansion.
This particular type of expansion is taken from the PhD Thesis from Tamellini
  "Polynomial approximation of PDEs with stochastic coefficients"

 <pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
 *!----------------------------------------------------------------------*/
#ifdef HAVE_FFTW

#include "../drt_inpar/inpar_mlmc.H"
#include "randomfield.H"
#include "randomfield_fourier.H"
#include "mlmc.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_io/io_pstream.H"
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <cmath>

#include <boost/random.hpp>



#include <fstream>


/*----------------------------------------------------------------------*/
/* standard constructor */
STR::UQ::RandomFieldFourier::RandomFieldFourier(unsigned int seed,Teuchos::RCP<DRT::Discretization> discret, const Teuchos::ParameterList& rfp)
{

  pi_=M_PI;
  seed_ = seed;
  myrank_ = discret->Comm().MyPID();

  // spatial dimension of random field  only 2 adn 3 supported
  dim_ = rfp.get<int>("RANDOM_FIELD_DIMENSION");

  if(dim_!=3)
      dserror("Dimension of random field for fourier series representation must be 3, fix your input file");

  // do we want to perform spectral matching with PSD
  bool temp = DRT::INPUT::IntegralValue<int>(rfp ,"SPECTRAL_MATCHING");
  if(temp)
    dserror("Spectral matching works only with spectral representation method");

  m_= rfp.get<int>("NUM_COS_TERMS");
  d_ = rfp.get<double>("CORRLENGTH");


  mean_ =rfp.get<double>("MEAN");
  sigma_0_= rfp.get<double>("SIGMA");

  ComputeBoundingBox(discret);
  periodicity_=1.1*largestlength_;

  // quick check whether correlation length of the field is < 0.35 periodicity
  // otherwise the fourier expansion cannot be used
  if(d_>0.35*periodicity_)
    dserror("correlation length is to large, choose a smaller correlation length or another representatoin method ");

  k_trunk_threshold_= rfp.get<int>("FOURIER_TRUNCATION_THRESHOLD");

  det_value_paracont_= rfp.get<double>("CONTBLENDVALUE");

  // Get correlation structure here can only use gaussian
  INPAR::MLMC::CorrStruct cstruct = DRT::INPUT::IntegralValue<INPAR::MLMC::CorrStruct>(rfp,"CORRSTRUCT");
  switch(cstruct)
  {
    case INPAR::MLMC::corr_gaussian:
    // nothing to do here
    break;
    default:
      dserror("Unknown Correlation structure");
    break;
  }



  // compute parameters for non gaussian pdf
  INPAR::MLMC::MarginalPdf mpdf = DRT::INPUT::IntegralValue<INPAR::MLMC::MarginalPdf>(rfp,"MARGINALPDF");
  switch(mpdf){
    case INPAR::MLMC::pdf_gaussian:
      marginal_pdf_=normal;
      break;
    case INPAR::MLMC::pdf_beta:
      dserror("Beta distribution not implemented for fourier representation");
     break;
    case INPAR::MLMC::pdf_lognormal:
      marginal_pdf_=lognormal;
       break;
    default:
      dserror("Unknown Marginal pdf");
      break;
  }

  // setup the truncation
  num_k_=0;
  double sum_ck=0.0;
  // determine size of truncation
  for (int k1=0;k1<=m_;k1++)
  {
    for (int k2=0;k2<=m_;k2++)
    {
      for (int k3=0;k3<=m_;k3++)
      {
        if( ((k1*k1)+(k2*k2)+(k3*k3)) < k_trunk_threshold_)
        {
          sum_ck=sum_ck+(ComputeLambda3D(k1)*ComputeLambda3D(k2)*ComputeLambda3D(k3));
          num_k_=num_k_+1;
        }
      }
    }
  }

  // now initialize vector containing k'i with known size
  kb_=Teuchos::rcp(new std::vector < std::vector<int> > (num_k_, std::vector<int>(3,0)));
  int index=0;
    for (int k1=0;k1<m_;k1++)
    {
      for (int k2=0;k2<m_;k2++)
      {
        for (int k3=0;k3<m_;k3++)
        {
          if( ((k1*k1)+(k2*k2)+(k3*k3)) < k_trunk_threshold_)
          {
            kb_->at(index).at(0)=k1;
            kb_->at(index).at(1)=k2;
            kb_->at(index).at(2)=k3;
            index++;
          }
        }
      }
    }

  // init vectors with random amplitudes
  xi_=Teuchos::rcp(new std::vector < std::vector<double> > (num_k_, std::vector<double>(8,0.0)));

  randomnumbergen_.seed((unsigned int)seed);

  boost::normal_distribution<>  normaldist(0.0,1.0);
  normalgen_ = Teuchos::rcp(new boost::variate_generator<randnumgen&,boost::normal_distribution<> >(randomnumbergen_,normaldist));

  // draw random amplitudes
  for (int i=0; i<num_k_; i++)
  {
    for (int j=0; j<8; j++)
    {
      xi_->at(i).at(j) = (*normalgen_)();
    }
  }
  // display some info
  IO::cout << "fraction of variance " << sum_ck << IO::endl;
  IO::cout << "Stochastic dimension is " << num_k_*8 << IO::endl;

}

void STR::UQ::RandomFieldFourier::CreateNewSample(unsigned int seed)
{
  CreateNewPhaseAngles(seed);
}

void STR::UQ::RandomFieldFourier::CreateNewPhaseAngles(unsigned int seed)
{
  randomnumbergen_.seed(seed);
  // draw random amplitudes
  for (int i=0; i<num_k_; i++)
    for (int j=0; j<8; j++)
      xi_->at(i).at(j) = (*normalgen_)();
}

double STR::UQ::RandomFieldFourier::EvalFieldAtLocation(std::vector<double> location,double paracont_parameter, bool writetofile, bool output)
{
  double tempgp=0.0;
  for (unsigned int i=0;i<kb_->size();i++)
  {
    double wk1=(*kb_)[i][0]*pi_/periodicity_;
    double wk2=(*kb_)[i][1]*pi_/periodicity_;
    double wk3=(*kb_)[i][2]*pi_/periodicity_;
    double lambda1=ComputeLambda3D(((*kb_)[i][0]));
    double lambda2=ComputeLambda3D(((*kb_)[i][1]));
    double lambda3=ComputeLambda3D(((*kb_)[i][2]));
   tempgp=tempgp+sqrt(lambda1*lambda2*lambda3)*(
        (*xi_)[i][0]*cos(wk1*location[0])*cos(wk2*location[1])*cos(wk3*location[2])
        +(*xi_)[i][1]*sin(wk1*location[0])*sin(wk2*location[1])*cos(wk3*location[2])
        +(*xi_)[i][2]*cos(wk1*location[0])*sin(wk2*location[1])*cos(wk3*location[2])
        +(*xi_)[i][3]*sin(wk1*location[0])*cos(wk2*location[1])*cos(wk3*location[2])
        +(*xi_)[i][4]*cos(wk1*location[0])*cos(wk2*location[1])*sin(wk3*location[2])
        +(*xi_)[i][5]*sin(wk1*location[0])*sin(wk2*location[1])*sin(wk3*location[2])
        +(*xi_)[i][6]*cos(wk1*location[0])*sin(wk2*location[1])*sin(wk3*location[2])
        +(*xi_)[i][7]*sin(wk1*location[0])*cos(wk2*location[1])*sin(wk3*location[2]) );

  }
  if(marginal_pdf_==normal)
    return det_value_paracont_*(1-paracont_parameter)+(paracont_parameter)*(mean_+sigma_0_*tempgp);
  else if (marginal_pdf_ == lognormal)
   return det_value_paracont_*(1.0-paracont_parameter)+(paracont_parameter)*exp((mean_+sigma_0_*tempgp));
  else
  {
    dserror("unknown marginal density ");
    return -1;
  }

}



double STR::UQ::RandomFieldFourier::ComputeLambda3D(int k)
{
  double lambda;
  if (k==0)
  {
    lambda=d_*sqrt(pi_)/(2*periodicity_);
  }
  else
  {
    lambda=d_*sqrt(pi_)/(periodicity_)*exp(-(pow(k*pi_*d_,2)/(4*pow(periodicity_,2.0))));
    // lambda=periodicity_*sqrt(pi_)/(periodicity_)*exp(-((k*pi_*d_)^2/(4*periodicity_)));
  }
  return lambda;
}


#endif
