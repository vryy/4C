/*!----------------------------------------------------------------------
\file randomfield.cpp
Created on: 15 May, 2014
\brief Class for generating samples of gaussian and non gaussian random fields based on spectral representation


 <pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
*/
#ifdef HAVE_FFTW

#include "randomfield.H"

//#include "mlmc.H"
//#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_io/io_pstream.H"
//#include <ctime>
//#include <cstdlib>
//#include <iostream>
//#include <complex>
//#include <cmath>
// For colored std::couts
//#include "../drt_lib/drt_colors.H"
//#include <boost/random.hpp>

// include fftw++ stuff for multidimensional FFT
#include"fftw3.h"

//#include <fstream>


#include "../drt_inpar/inpar_mlmc.H"


void STR::UQ::RandomField::ComputeBoundingBox(Teuchos::RCP<DRT::Discretization> discret)
{
  // root bounding Box
  std::vector<double> maxrbb;

  maxrbb.push_back(-10.0e19);
  maxrbb.push_back(-10.0e19);
  maxrbb.push_back(-10.0e19);

  std::vector<double> minrbb;
  minrbb.push_back(10.0e19);
  minrbb.push_back(10.0e19);
  minrbb.push_back(10.0e19);

  bb_max_.push_back(-10.0e19);
  bb_max_.push_back(-10.0e19);
  bb_max_.push_back(-10.0e19);

  bb_min_.push_back(10.0e19);
  bb_min_.push_back(10.0e19);
  bb_min_.push_back(10.0e19);
  {
    for (int lid = 0; lid <discret->NumMyColNodes(); ++lid)
    {
      const DRT::Node* node = discret->lColNode(lid);
      // check if greater than maxrbb
      if (maxrbb[0]<node->X()[0])
        maxrbb[0]=node->X()[0];
      if (maxrbb[1]<node->X()[1])
        maxrbb[1]=node->X()[1];
      if (maxrbb[2]<node->X()[2])
        maxrbb[2]=node->X()[2];
      // check if smaller than minrbb
      if (minrbb[0]>node->X()[0])
        minrbb[0]=node->X()[0];
      if (minrbb[1]>node->X()[1])
        minrbb[1]=node->X()[1];
      if (minrbb[2]>node->X()[2])
        minrbb[2]=node->X()[2];
    }
  }

  discret->Comm().MaxAll(&maxrbb[0],&bb_max_[0],3);
  discret->Comm().MinAll(&minrbb[0],&bb_min_[0],3);


  discret->Comm().Barrier();

  //compute largest dimension of the problem
  largestlength_=-10e19;

  for(int i=0;i<3;i++)
  {
    if(bb_max_[i]-bb_min_[i]>largestlength_)
      largestlength_=bb_max_[i]-bb_min_[i];
  }

  if (myrank_ == 0)
  {
    IO::cout<< "min " << bb_min_[0] << " "<< bb_min_[1]  << " "<< bb_min_[2] << IO::endl;
    IO::cout<< "max " << bb_max_[0] << " "<< bb_max_[1]  << " "<< bb_max_[2] << IO::endl;
    IO::cout<< "largest length " << largestlength_ << IO::endl;
  }






}


#endif // HAVE_FFTw
