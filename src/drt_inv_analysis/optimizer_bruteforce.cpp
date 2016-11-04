/*----------------------------------------------------------------------*/
/*!
 * \file optimizer_bruteforce.cpp
 * \brief Brute force sampling of the parameter space
 *
<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
#include "optimizer_bruteforce.H"

#include "Epetra_Vector.h"

#include "invana_base.H"
#include "matpar_manager.H"
#include "../drt_lib/drt_dserror.H"

#include <fstream>


/*----------------------------------------------------------------------*/
INVANA::OptimizerBruteForce::OptimizerBruteForce(const Teuchos::ParameterList& invp) :
  OptimizerBase(invp)
{}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerBruteForce::Setup()
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerBruteForce::Integrate()
{
  // make sure parameter space dimension is 1
  if (not (GetSolution()->GlobalLength() == 1) )
    dserror("BruteForce only works for a parameter space of dim 1");

  // only the first vector has optimization parameters
  Epetra_Vector params(*(OptProb()->Matman()->InitialParams()(0)));

  int nsamples=1000;

  double rlower = -2.0;
  double rupper = 1.0;

  double dr = (rupper-rlower)/nsamples;

  std::vector<double> func(nsamples);

  double* p;
  params.ExtractView(&p);

  std::cout << "sampling in the range " << p[0]+rlower << " - " << p[0]+rupper << std::endl;

  double p0 = p[0]+rlower;
  p[0] = p0;

  for (int i=0; i<nsamples; i++)
  {
    OptProb()->Evaluate(params,&func[i],Teuchos::null);
    p[0] += dr;
  }

  std::ofstream myfile ("sampling.txt");
  if (myfile.is_open())
  {
    for (int i=0; i<nsamples; i++)
    {
      myfile << p0 <<  " " << func[i] << "\n";
      p0 += dr;
    }
    myfile.close();
  }
  else std::cout << "Unable to open file";

  return;
}


/*----------------------------------------------------------------------*/
void INVANA::OptimizerBruteForce::ReadRestart(int run)
{
  dserror("no restart for the brute force sampler");
  return;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerBruteForce::WriteRestart()
{
  // nothing written
  return;
}
