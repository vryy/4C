/*----------------------------------------------------------------------*/
/*!
 * \file pat_imagereconstruction.H

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/staff/svenja-schoeder/
            089 - 289-15271
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef PAT_IMAGERECONSTRUCTION_H
#define PAT_IMAGERECONSTRUCTION_H

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseVector.h"
#include "../drt_inpar/inpar_acou.H"

//{ forward declarations
namespace DRT
{
  class Discretization;
  class DiscretizationHDG;
  class ResultTest;
}
namespace ACOU
{
  class AcouTimeInt;
  class PatMatParManager;
}
namespace SCATRA
{
  class TimIntStationary;
}
namespace LINALG
{
  class Solver;
  class MapExtractor;
}
namespace IO
{
  class DiscretizationWriter;
}
namespace INVANA
{
  class MatParManager;
  class RegularizationBase;
}
//} forward declarations

namespace ACOU
{

class PATSearchDirection;
class PATLineSearch;

class PatImageReconstruction
{
  friend class PATLineSearch;

public:

  //! standard constructor
  PatImageReconstruction(
    Teuchos::RCP<DRT::Discretization>      scatradis,
    Teuchos::RCP<DRT::DiscretizationHDG>   acoudis,
    Teuchos::RCP<Teuchos::ParameterList>   scatrapara,
    Teuchos::RCP<Teuchos::ParameterList>   acoupara,
    Teuchos::RCP<LINALG::Solver>           scatrasolv,
    Teuchos::RCP<LINALG::Solver>           acousolv,
    Teuchos::RCP<IO::DiscretizationWriter> scatraout,
    Teuchos::RCP<IO::DiscretizationWriter> acouout);

  //! standard destructor
  virtual ~PatImageReconstruction() {}

  //! iterative optimization loop
  void Optimize();

  //! first solve with given material parameters
  void InitialRun();

  //! read restart from given iteration
  virtual void ReadRestart(int restartiteration) = 0;

  //! create test field
  Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  //! return the scatra discretization (for result test frame work)
  const Teuchos::RCP<DRT::Discretization> ScatraDiscretization() { return scatra_discret_; }

  //! output of element based absorption coefficient vector (for result test)
  const Teuchos::RCP<Epetra_MultiVector> ElementMatVec();

  //! replace parameters
  virtual void ReplaceParams(Teuchos::RCP<Epetra_MultiVector> params) = 0;

protected:

  //@{ everything we need for the time integrations (scatra and acoustic)
  //! scatra discretization to work with
  Teuchos::RCP<DRT::Discretization> scatra_discret_;

  //! acoustic discretization to work with
  Teuchos::RCP<DRT::DiscretizationHDG> acou_discret_;

  //! scatra parameters
  Teuchos::RCP<Teuchos::ParameterList> scatraparams_;

  //! scatra additional parameters
  Teuchos::RCP<Teuchos::ParameterList> scatraextraparams_;

  //! acoustic parameters
  Teuchos::RCP<Teuchos::ParameterList> acouparams_;

  //! scatra solver
  Teuchos::RCP<LINALG::Solver> scatrasolver_;

  //! acoustic solver
  Teuchos::RCP<LINALG::Solver> acousolver_;

  //! scatra output
  Teuchos::RCP<IO::DiscretizationWriter> scatraoutput_;

  //! acoustic output
  Teuchos::RCP<IO::DiscretizationWriter> acououtput_;

  //! scatra base algorithm to work with
  Teuchos::RCP<SCATRA::TimIntStationary> scatraalgo_;

  //! acoustic base algorithms to work with
  Teuchos::RCP<ACOU::AcouTimeInt> acoualgo_;

  //! time integration scheme
  INPAR::ACOU::DynamicType dyna_;

  //! physical type
  INPAR::ACOU::PhysicalType phys_;
  //@}

  //@{ general
  //! output name given in command line
  std::string name_;

  //! time measurement
  double tstart_;

  //! prescribed objective function tolerance
  double tol_;

  //! actual iteration number
  int    iter_;

  //! maximal iteration number for inverse analysis
  int    maxiter_;

  //! processor id
  int myrank_;

  //! auxiliary variable for nicely named output files
  int output_count_;

  //! last forward acoustic output count
  int last_acou_fw_output_count_;

  //! is acoustical mesh conform with optical
  bool meshconform_;

  //! the objective function value
  double J_;

  //! the initial objective function value
  double J_start_;

  //! the error value
  double error_;

  //! the initial error value
  double error_start_;

  //! list of observed node gids
  std::vector<int> nodes_;

  //! simulated "measurement" values of acoustic forward problem
  Teuchos::RCP<Epetra_MultiVector> acou_rhs_;

  //! measurement values from monitor file
  Teuchos::RCP<Epetra_MultiVector> acou_rhsm_;

  //! map belonging to acou_rhs_ and acou_rhsm_ holding the nodes of the measurement boundary
  Teuchos::RCP<Epetra_Map> abcnodes_map_;

  //! node based vector holding absorption coefficient
  Teuchos::RCP<Epetra_Vector> node_reac_;
  //@}

  //@{ solution vectors
  //! solution of the optical forward problem
  Teuchos::RCP<Epetra_Vector> phi_;

  //! solution of the optical adjoint problem
  Teuchos::RCP<Epetra_Vector> adjoint_phi_;

  //! solution of the acoustical adjoint problem at initial step
  Teuchos::RCP<Epetra_Vector> adjoint_psi_;
  //@}

  //@{ helper classes
  //! class to calculate directions
  Teuchos::RCP<PATSearchDirection> searchdirection_;

  //! line search class
  Teuchos::RCP<PATLineSearch> linesearch_;

  //! material handling for scatra parameters
  Teuchos::RCP<INVANA::MatParManager> scatra_matman_;
  //@}

  //@{ read of monitored values and helper functions
  //! read in measured values from the monitor file
  void ReadMonitor(std::string monitorfilename, double dtacou);

  //! calculate the distance between mics and nods in x,y,z
  double ReadMonitorDelta(double coord_M_x,double coord_M_y, double coord_M_z, double coord_N_x,double coord_N_y, double coord_N_z);

  //! write interpolated data
  void ReadMonitorInterpol(const double nod_coords[3],std::vector<std::vector<double> > mic_coords, unsigned int mic_1, unsigned int mic_2, int nmic, unsigned int nod, int timesteps, Epetra_SerialDenseVector& curve, Epetra_SerialDenseVector& inter_curve);

  //! calculate epsilon
  double ReadMonitorGetEpsilon(int nnodes);
  //@}

  //! the underlying subproblems
  //@{
  //! optical forward problem
  void SolveStandardScatra();

  //! acoustical forward problem
  void SolveStandardAcou();

  //! acoustical adjoint problem
  void SolveAdjointAcou();

  //! optical adjoint problem
  void SolveAdjointScatra();
  //@}

  //! calculate the value of the objective function
  virtual double EvalulateObjectiveFunction() = 0;

  //! calculate the value of the error in the objective function
  void EvaluateError();

  //! build the gradient
  virtual void EvaluateGradient() = 0;

  //! calculate the norm of the dot product of the gradient and a given direction
  virtual void CalculateGradDirNorm(const Epetra_MultiVector& bvector, const Epetra_Map& uniquemap, double* result) = 0;

  //! perform one iteration
  virtual bool PerformIteration() = 0;

  //! output of optical material parameters
  virtual void OutputReactionAndDiffusion();

  //! compute node based reaction vector
  void ComputeNodeBasedReactionCoefficient();

  //! calculate right hand side vector for adjoint optical Problem
  Teuchos::RCP<Epetra_Vector> CalculateAdjointOptiRhsvec(Teuchos::RCP<Epetra_Vector> acounodevec);

  //! output statistics
  void OutputStats();
};

class PatImageReconstructionOpti : public PatImageReconstruction
{
public:

  //! standard constructor
  PatImageReconstructionOpti(
    Teuchos::RCP<DRT::Discretization>      scatradis,
    Teuchos::RCP<DRT::DiscretizationHDG>   acoudis,
    Teuchos::RCP<Teuchos::ParameterList>   scatrapara,
    Teuchos::RCP<Teuchos::ParameterList>   acoupara,
    Teuchos::RCP<LINALG::Solver>           scatrasolv,
    Teuchos::RCP<LINALG::Solver>           acousolv,
    Teuchos::RCP<IO::DiscretizationWriter> scatraout,
    Teuchos::RCP<IO::DiscretizationWriter> acouout);

  //! standard destructor
  virtual ~PatImageReconstructionOpti() {}

  //! read restart from given iteration
  virtual void ReadRestart(int restartiteration) {return;}

  //! replace parameters
  virtual void ReplaceParams(Teuchos::RCP<Epetra_MultiVector> params);

protected:

  //! regularization for scatra parameters
  Teuchos::RCP<INVANA::RegularizationBase> scatra_regman_;

  //! gradient of the objective function with respect to scatra parameters
  Teuchos::RCP<Epetra_MultiVector> scatra_objgrad_;

  //! build the gradient
  virtual void EvaluateGradient();

private:

  //! calculate the value of the objective function
  virtual double EvalulateObjectiveFunction();

  //! perform one iteration
  virtual bool PerformIteration();

  //! calculate the norm of the dot product of the gradient and a given direction
  virtual void CalculateGradDirNorm(const Epetra_MultiVector& bvector, const Epetra_Map& uniquemap, double* result);
};

class PatImageReconstructionOptiAcou : public PatImageReconstructionOpti
{
public:

  //! standard constructor
  PatImageReconstructionOptiAcou(
    Teuchos::RCP<DRT::Discretization>      scatradis,
    Teuchos::RCP<DRT::DiscretizationHDG>   acoudis,
    Teuchos::RCP<Teuchos::ParameterList>   scatrapara,
    Teuchos::RCP<Teuchos::ParameterList>   acoupara,
    Teuchos::RCP<LINALG::Solver>           scatrasolv,
    Teuchos::RCP<LINALG::Solver>           acousolv,
    Teuchos::RCP<IO::DiscretizationWriter> scatraout,
    Teuchos::RCP<IO::DiscretizationWriter> acouout);

  //! standard destructor
  virtual ~PatImageReconstructionOptiAcou() {}

  //! read restart from given iteration
  virtual void ReadRestart(int restartiteration) {return;}

  //! replace parameters
  virtual void ReplaceParams(Teuchos::RCP<Epetra_MultiVector> params);

private:

  //! material handling for acoustic parameters
  Teuchos::RCP<INVANA::MatParManager> acou_matman_;

  //! regularization for acoustic parameters
  Teuchos::RCP<INVANA::RegularizationBase> acou_regman_;

  //! gradient of the objective function with respect to acoustic parameters
  Teuchos::RCP<Epetra_MultiVector> acou_objgrad_;

  //! class to calculate directions (needs an additional one for the different parameter layouts)
  Teuchos::RCP<PATSearchDirection> acou_searchdirection_;

  //! indicator if acoustic or optic material is optimized (true for acou)
  bool optoracou_;

  //! calculate the value of the objective function
  virtual double EvalulateObjectiveFunction();

  //! build the gradient
  virtual void EvaluateGradient();

  //! perform one iteration
  virtual bool PerformIteration();

  //! calculate the norm of the dot product of the gradient and a given direction
  virtual void CalculateGradDirNorm(const Epetra_MultiVector& bvector, const Epetra_Map& uniquemap, double* result);
};

class PatImageReconstructionSegmentation : public PatImageReconstruction
{
public:

  //! standard constructor
  PatImageReconstructionSegmentation(
    Teuchos::RCP<DRT::Discretization>      scatradis,
    Teuchos::RCP<DRT::DiscretizationHDG>   acoudis,
    Teuchos::RCP<Teuchos::ParameterList>   scatrapara,
    Teuchos::RCP<Teuchos::ParameterList>   acoupara,
    Teuchos::RCP<LINALG::Solver>           scatrasolv,
    Teuchos::RCP<LINALG::Solver>           acousolv,
    Teuchos::RCP<IO::DiscretizationWriter> scatraout,
    Teuchos::RCP<IO::DiscretizationWriter> acouout);

  //! standard destructor
  virtual ~PatImageReconstructionSegmentation() {}

  //! read restart from given iteration
  virtual void ReadRestart(int restartiteration) {return;}

  //! replace parameters
  virtual void ReplaceParams(Teuchos::RCP<Epetra_MultiVector> params);

private:

  //! number of materials
  unsigned int nummats_;

  // table holding the material values
  std::vector<std::vector<double> > materialtable_;

  //! penalty parameter
  double penaltyparam_;

  //! regularization for scatra parameters
  Teuchos::RCP<INVANA::RegularizationBase> scatra_regman_;

  //! material handling for acoustic parameters
  Teuchos::RCP<INVANA::MatParManager> acou_matman_;

  //! regularization for acoustic parameters
  Teuchos::RCP<INVANA::RegularizationBase> acou_regman_;

  //! vector holding the concentration values
  Teuchos::RCP<Epetra_MultiVector> k_;

  //! gradient of the objective function with respect to the concentrations
  Teuchos::RCP<Epetra_MultiVector> k_objgrad_;

  //! calculate the value of the objective function
  virtual double EvalulateObjectiveFunction();

  //! build the gradient
  virtual void EvaluateGradient();

  //! read of predefined materials
  void ReadMaterials(std::string materialfilename);

  //! initialization of concentrations
  void InitConcentrations();

  //! perform one iteration
  virtual bool PerformIteration();

  //! calculate the norm of the dot product of the gradient and a given direction
  virtual void CalculateGradDirNorm(const Epetra_MultiVector& bvector, const Epetra_Map& uniquemap, double* result);
};

} // namespace ACOU


#endif /* PAT_IMAGERECONSTRUCTION_H */
