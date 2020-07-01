/*----------------------------------------------------------------------*/
/*! \file
\brief Brute force sampling of the parameter space

\level 3

*/
/*----------------------------------------------------------------------*/
#include "optimizer_bruteforce.H"

#include "Epetra_Vector.h"

#include "invana_base.H"
#include "matpar_manager.H"
#include "objective_funct.H"
#include "../drt_lib/drt_dserror.H"

#include <fstream>


/*----------------------------------------------------------------------*/
INVANA::OptimizerBruteForce::OptimizerBruteForce(const Teuchos::ParameterList& invp)
    : OptimizerBase(invp)
{
}

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
  if (not(GetSolution()->GlobalLength() == 1))
    dserror("BruteForce only works for a parameter space of dim 1");

  // only the first vector has optimization parameters
  Epetra_Vector params(*(OptProb()->Matman()->InitialParams()(0)));

  // create linspace
  double rlower = 0.0;
  double rupper = 1.0;
  int nsamples = 100;
  std::istringstream bflinspace(Teuchos::getNumericStringParameter(Inpar(), "BFLINSPACE"));
  std::string word;
  int wordi = 0;
  while (bflinspace >> word)
  {
    if (wordi == 0)
      rlower = std::atof(word.c_str());
    else if (wordi == 1)
      rupper = std::atof(word.c_str());
    else if (wordi == 2)
      nsamples = std::atoi(word.c_str());

    wordi += 1;
  }
  double dr = (rupper - rlower) / nsamples;

  std::vector<double> func(nsamples);

  double* p;
  params.ExtractView(&p);

  if (OptProb()->Comm().MyPID() == 0)
  {
    printf("Sampling a %d steps in the interval [%.2f,%.2f]\n", nsamples, rlower, rupper);
    fflush(stdout);
  }

  double p0 = rlower;
  if (OptProb()->Comm().MyPID() == 0) p[0] = p0;

  double prob_scale_fac = OptProb()->ObjectiveFunct()->GetScaleFac();
  for (int i = 0; i < nsamples; i++)
  {
    OptProb()->Evaluate(params, &func[i], Teuchos::null);
    func[i] /= prob_scale_fac;
    if (OptProb()->Comm().MyPID() == 0) p[0] += dr;
  }

  if (OptProb()->Comm().MyPID() == 0)
  {
    std::ofstream myfile("sampling.txt");
    if (myfile.is_open())
    {
      for (int i = 0; i < nsamples; i++)
      {
        myfile << p0 << " " << func[i] << "\n";
        p0 += dr;
      }
      myfile.close();
    }
    else
      std::cout << "Unable to open file";
  }

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
