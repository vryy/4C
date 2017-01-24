/*----------------------------------------------------------------------*/
/*!
\file pat_imagereconstruction.H

\brief image reconstruction

\level 3

<pre>
\level 3

\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
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
#include "../drt_lib/drt_dserror.H"


//{ forward declarations
namespace DRT
{
  class Discretization;
  class DiscretizationHDG;
  class ResultTest;
  class Element;
}
namespace ACOU
{
  class AcouTimeInt;
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
//} forward declarations

namespace ACOU
{

class PATSearchDirection;
class PATLineSearch;
class PATRegula;

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
  virtual void Optimize();

  //! sample objective function
  virtual void SampleObjectiveFunction() { dserror("not implemented for base class"); }

  //! first solve with given material parameters
  virtual void InitialRun();

  //! read restart from given iteration
  void ReadRestart(int restartoutputcount);

  //! create test field
  Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  //! return the scatra discretization (for result test frame work)
  const Teuchos::RCP<DRT::Discretization> ScatraDiscretization() { return scatra_discret_; }

  //! output of element based absorption coefficient vector (for result test)
  const Teuchos::RCP<Epetra_MultiVector> ElementMatVec();

  //! replace parameters
  virtual void ReplaceParams(Teuchos::RCP<Epetra_Vector> params);

  //! finite difference gradient check
  virtual void FDCheck();

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
  //! vector holding the mat ids
  std::vector<int> opti_matids_;

  //! vector indicating if element is optimized or not
  Teuchos::RCP<Epetra_Vector> opti_opt_ind_;

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

  //! time reversal as initial guess
  bool timereversal_;

  //! patch type
  INPAR::ACOU::PatchType patchtype_;

  //! fdcheck
  bool fdcheck_;

  //! the objective function value
  double J_;

  //! the initial objective function value
  double J_start_;

  //! the error value
  double error_;

  //! the initial error value
  double error_start_;

  //! consider impulse response or not
  bool conv_imp_resp_;

  //! overwrite output (recommended for big simulations)
  bool overwrite_output_;

  //! the impulse respone already calculated for dtacou_
  Epetra_SerialDenseVector imp_resp_;

  //! list of observed node gids
  std::vector<int> nodes_;

  //! simulated "measurement" values of acoustic forward problem
  Teuchos::RCP<Epetra_MultiVector> acou_rhs_;

  //! measurement values from monitor file
  Teuchos::RCP<Epetra_MultiVector> acou_rhsm_;

  //! map belonging to acou_rhs_ and acou_rhsm_ holding the nodes of the measurement boundary
  Teuchos::RCP<Epetra_Map> abcnodes_map_;
  Teuchos::RCP<Epetra_Map> abcnodes_colmap_;

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
  //! line search
  Teuchos::RCP<PATLineSearch> linesearch_;

  //! material handling for acoustic parameters
  Teuchos::RCP<Epetra_Vector> reac_vals_;

  //! gradient of the objective function with respect to acoustic parameters
  Teuchos::RCP<Epetra_Vector> reac_objgrad_;

  //! class to calculate directions (needs an additional one for the different parameter layouts)
  Teuchos::RCP<PATSearchDirection> reac_searchdirection_;

  //! regularization for reaction
  Teuchos::RCP<PATRegula> reac_regula_;

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

  //! initial estimation with time reversal
  void TimeReversalEstimate();

  //! helper funciton for time reversal
  void UpdateAbsorptionCoefficientFromTimeReversal();

  //! calculate the value of the objective function
  virtual double EvalulateObjectiveFunction();

  //! calculate the value of the error in the objective function
  void EvaluateError();

  //! build the gradient
  virtual void EvaluateGradient();

  //! calculate the norm of the dot product of the gradient and a given direction
  virtual void CalculateGradDirNorm(const Epetra_Vector& bvector, const Epetra_Map& uniquemap, double* result);

  //! perform one iteration
  virtual bool PerformIteration();

  //! output of optical material parameters
  virtual void OutputReactionAndDiffusion();

  //! compute node based reaction vector
  void ComputeNodeBasedReactionCoefficient();

  //! calculate right hand side vector for adjoint optical Problem
  Teuchos::RCP<Epetra_Vector> CalculateAdjointOptiRhsvec(Teuchos::RCP<Epetra_Vector> acounodevec);

  //! output statistics
  void OutputStats();

  //! evaluate absorption coefficient
  void EvaluateReacGrad();

  //! build basis patches
  void CheckNeighborsReacGrad(DRT::Element* actele, int owner, Teuchos::RCP<Epetra_Vector> setsids, double set, double reacval, double interval, Teuchos::RCP<Epetra_Vector> auxvals);

  //! calculation of error in parameter field for analytic examples
  virtual void ComputeParameterError() { std::cout<<"here could be a parameter error calculation"<<std::endl; return; }

  //! convert gradient to different basis
  void ConvertGradient(Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<Epetra_Vector> gradient);

  //! bring the restart parameters to the vectors and elements (each derived class implements this, since derived classes have different parameter fields)
  virtual void SetRestartParameters(Teuchos::RCP<Epetra_Vector> reacs, Teuchos::RCP<Epetra_Vector> diffs, Teuchos::RCP<Epetra_Vector> cs, Teuchos::RCP<Epetra_Vector> rhos);

};

class PatImageReconstructionReduction: public PatImageReconstruction
{
public:
  //! standard constructor
  PatImageReconstructionReduction(
    Teuchos::RCP<DRT::Discretization>      scatradis,
    Teuchos::RCP<DRT::DiscretizationHDG>   acoudis,
    Teuchos::RCP<Teuchos::ParameterList>   scatrapara,
    Teuchos::RCP<Teuchos::ParameterList>   acoupara,
    Teuchos::RCP<LINALG::Solver>           scatrasolv,
    Teuchos::RCP<LINALG::Solver>           acousolv,
    Teuchos::RCP<IO::DiscretizationWriter> scatraout,
    Teuchos::RCP<IO::DiscretizationWriter> acouout);

  //! iterative optimization loop
  virtual void Optimize();

private:
  double reductioncuttime_;
};

class PatImageReconstructionOptiSplit : public PatImageReconstruction
{
public:

  //! standard constructor
  PatImageReconstructionOptiSplit(
    Teuchos::RCP<DRT::Discretization>      scatradis,
    Teuchos::RCP<DRT::DiscretizationHDG>   acoudis,
    Teuchos::RCP<Teuchos::ParameterList>   scatrapara,
    Teuchos::RCP<Teuchos::ParameterList>   acoupara,
    Teuchos::RCP<LINALG::Solver>           scatrasolv,
    Teuchos::RCP<LINALG::Solver>           acousolv,
    Teuchos::RCP<IO::DiscretizationWriter> scatraout,
    Teuchos::RCP<IO::DiscretizationWriter> acouout);

  //! standard destructor
  virtual ~PatImageReconstructionOptiSplit() {}

  //! replace parameters
  virtual void ReplaceParams(Teuchos::RCP<Epetra_Vector> params);

  //! finite difference gradient check
  virtual void FDCheck();

protected:

  //! number of iterations for sequential optimization
  int sequenzeiter_;

  //! material handling for acoustic parameters
  Teuchos::RCP<Epetra_Vector> diff_vals_;

  //! gradient of the objective function with respect to acoustic parameters
  Teuchos::RCP<Epetra_Vector> diff_objgrad_;

  //! class to calculate directions (needs an additional one for the different parameter layouts)
  Teuchos::RCP<PATSearchDirection> diff_searchdirection_;

  //! regularization for diffusion
  Teuchos::RCP<PATRegula> diff_regula_;

  //! build the gradient
  virtual void EvaluateGradient();

  //! calculate the value of the objective function
  virtual double EvalulateObjectiveFunction();

  //! calculate the norm of the dot product of the gradient and a given direction
  virtual void CalculateGradDirNorm(const Epetra_Vector& bvector, const Epetra_Map& uniquemap, double* result);

  //! perform one iteration
  virtual bool PerformIteration();

  //! evaluate diffusion coefficient gradient
  void EvaluateDiffGrad();

  //! build basis patches
  void CheckNeighborsDiffGrad(DRT::Element* actele, int owner, Teuchos::RCP<Epetra_Vector> setsids, double set, double diffval, double interval, Teuchos::RCP<Epetra_Vector> auxvals);

  //! indicator if reaction or diffusion or speed of sound or density is optimized
  int reacordifforcorrho_;

  //! calculation of error in parameter field for analytic examples
  virtual void ComputeParameterError() { std::cout<<"here could be a parameter error evaluation"<<std::endl; return; }

  //! bring restart parameters to vectors and elements
  virtual void SetRestartParameters(Teuchos::RCP<Epetra_Vector> reacs, Teuchos::RCP<Epetra_Vector> diffs, Teuchos::RCP<Epetra_Vector> cs, Teuchos::RCP<Epetra_Vector> rhos);

private:

};

class PatImageReconstructionOptiSplitAcouSplit : public PatImageReconstructionOptiSplit
{
public:

  //! standard constructor
  PatImageReconstructionOptiSplitAcouSplit(
    Teuchos::RCP<DRT::Discretization>      scatradis,
    Teuchos::RCP<DRT::DiscretizationHDG>   acoudis,
    Teuchos::RCP<Teuchos::ParameterList>   scatrapara,
    Teuchos::RCP<Teuchos::ParameterList>   acoupara,
    Teuchos::RCP<LINALG::Solver>           scatrasolv,
    Teuchos::RCP<LINALG::Solver>           acousolv,
    Teuchos::RCP<IO::DiscretizationWriter> scatraout,
    Teuchos::RCP<IO::DiscretizationWriter> acouout);

  //! standard destructor
  virtual ~PatImageReconstructionOptiSplitAcouSplit() {}

  //! replace parameters
  virtual void ReplaceParams(Teuchos::RCP<Epetra_Vector> params);

  //! sample objective function
  virtual void SampleObjectiveFunction();

  //! finite difference gradient check
  virtual void FDCheck();

protected:

  //! vector indicating if element is optimized or not
  Teuchos::RCP<Epetra_Vector> acou_opt_ind_;

  //! material handling for acoustic parameters
  Teuchos::RCP<Epetra_Vector> c_vals_;

  //! gradient of the objective function with respect to acoustic parameters
  Teuchos::RCP<Epetra_Vector> c_objgrad_;

  //! class to calculate directions (needs an additional one for the different parameter layouts)
  Teuchos::RCP<PATSearchDirection> c_searchdirection_;

  //! regularization for c
  Teuchos::RCP<PATRegula> c_regula_;

  //! material handling for acoustic parameters
  Teuchos::RCP<Epetra_Vector> rho_vals_;

  //! gradient of the objective function with respect to acoustic parameters
  Teuchos::RCP<Epetra_Vector> rho_objgrad_;

  //! class to calculate directions (needs an additional one for the different parameter layouts)
  Teuchos::RCP<PATSearchDirection> rho_searchdirection_;

  //! regularization for density
  Teuchos::RCP<PATRegula> rho_regula_;

  //! vector holding the mat ids
  std::vector<int> acou_matids_;

  //! build the gradient
  virtual void EvaluateGradient();

  //! reduce basis
  void ReduceBasis();

  //! recalculate gradients
  void ReducedGradientsCalculation(Teuchos::RCP<Epetra_Vector> setids, bool acouonly);

  //! patch self
  void ReducePatchSelf(Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<Epetra_Vector> patchvec, Teuchos::RCP<Epetra_Vector> opt_ind);

  //! patch reaction values
  void ReducePatchReacVals();

  //! check neighbors patch self
  void CheckNeighborsPatchSelf(Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<Epetra_Vector> patchvec, DRT::Element* actele, int owner, Teuchos::RCP<Epetra_Vector> setsids, double set, double val, double interval, Teuchos::RCP<Epetra_Vector> auxvals);

  //! check neighbors patch reaction values
  void CheckNeighborsPatchReacVals(DRT::Element* actele, int owner, Teuchos::RCP<Epetra_Vector> setsids, double set, double reacval, double interval, Teuchos::RCP<Epetra_Vector> auxvals);

  //! calculate the value of the objective function
  virtual double EvalulateObjectiveFunction();

  //! calculate the norm of the dot product of the gradient and a given direction
  virtual void CalculateGradDirNorm(const Epetra_Vector& bvector, const Epetra_Map& uniquemap, double* result);

  //! calculation of error in parameter field for analytic examples
  virtual void ComputeParameterError();

  //! bring restart parameters to vectors and elements
  virtual void SetRestartParameters(Teuchos::RCP<Epetra_Vector> reacs, Teuchos::RCP<Epetra_Vector> diffs, Teuchos::RCP<Epetra_Vector> cs, Teuchos::RCP<Epetra_Vector> rhos);

private:

  //! perform one iteration
  virtual bool PerformIteration();

  //! evaluate c gradient and rho gradient
  void EvaluateCGrad();
  void EvaluateRhoGrad();

};

class PatImageReconstructionOptiSplitAcouIdent : public PatImageReconstructionOptiSplit
{
public:

  //! standard constructor
  PatImageReconstructionOptiSplitAcouIdent(
    Teuchos::RCP<DRT::Discretization>      scatradis,
    Teuchos::RCP<DRT::DiscretizationHDG>   acoudis,
    Teuchos::RCP<Teuchos::ParameterList>   scatrapara,
    Teuchos::RCP<Teuchos::ParameterList>   acoupara,
    Teuchos::RCP<LINALG::Solver>           scatrasolv,
    Teuchos::RCP<LINALG::Solver>           acousolv,
    Teuchos::RCP<IO::DiscretizationWriter> scatraout,
    Teuchos::RCP<IO::DiscretizationWriter> acouout);

  //! standard destructor
  virtual ~PatImageReconstructionOptiSplitAcouIdent() {}

  //! replace parameters
  virtual void ReplaceParams(Teuchos::RCP<Epetra_Vector> params);

private:

  //! number of materials
  unsigned int nummats_;

  // table holding the material values
  std::vector<std::vector<double> > materialtable_;

  //! number of iterations for sequential optimization
  int sequenzeiter_;

  //! material handling for acoustic parameters
  Teuchos::RCP<Epetra_Vector> c_vals_;

  //! material handling for acoustic parameters
  Teuchos::RCP<Epetra_Vector> rho_vals_;

  //! vector holding the mat ids
  std::vector<int> acou_matids_;

  //! bool determines how acoustical materials are identified
  bool acouident_avg_;

  //! perform one iteration
  virtual bool PerformIteration();

  //! update acoustical parameters according to the optical ones
  void UpdateAcousticalParameters();

  //! read of predefined materials
  void ReadMaterials(std::string materialfilename);

  //! calculation of error in parameter field for analytic examples
  virtual void ComputeParameterError();

  //! bring restart parameters to vectors and elements
  virtual void SetRestartParameters(Teuchos::RCP<Epetra_Vector> reacs, Teuchos::RCP<Epetra_Vector> diffs, Teuchos::RCP<Epetra_Vector> cs, Teuchos::RCP<Epetra_Vector> rhos);

};

} // namespace ACOU


#endif /* PAT_IMAGERECONSTRUCTION_H */
