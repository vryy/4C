/*!----------------------------------------------------------------------
\file randomfield.h
Created on: Apr 12, 2011
\brief Class for generating samples of random fields based on spectral representation

 <pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
 *!----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef RANDOMFIELD_H_
#define RANDOMFIELD_H_
//testing
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_LAPACK.h"
#include "Epetra_Vector.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_RefCountPtr.hpp>
#include <Epetra_SerialDenseVector.h>
#include <math.h>
#include "../drt_adapter/adapter_structure.H"

#include <cmath>
#include <ctime>

//including random number library of blitz for generation of random phase angles
#include <random/uniform.h>




class RandomField {

private:
  // some variables

  double pi_;// = M_PI;
  double d_; //correlation length
  // set variance to 0.25 by default
  double sigma_0_; //= 0.25;
  int N_;
  double dkappa_;
  unsigned int  seed_;
  vector<double>  Phi_ ;

  //! \brief generates uniformly distributed random numbers in interval [0;1], taken from statmech
   ranlib::UniformClosed<double> uniformclosedgen_;



public:
  // Standard constructor
  RandomField(unsigned int seed,double sigma_0, double corr_length);

  // Evaluate RandomField at specific location
  double EvalRandomField( double x, double y, double z);

  // Calculate Power Spectral Density
  double  PowerSpectralDensity( double kappa_x, double kappa_y, double kappa_z );

}; // class RandomField

#endif /* RANDOMFIELD_H_ */
#endif  // #ifdef CCADISCRET
