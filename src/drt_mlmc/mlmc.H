/*!----------------------------------------------------------------------
\file  mlmc.H
\brief Class for performing Multi Level Monte Carlo (MLMC)analysis of structure


 <pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
 *!----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef HAVE_FFTW
#ifndef MLMC_H
#define MLMC_H

#include "Teuchos_RefCountPtr.hpp"
#include "../drt_adapter/adapter_structure.H"

namespace STR
{

/// Multi Level Monte Carlo Simulation of Structures with
/// some sort of stochastic property


class MLMC {

private:

  /// my structure discretisation to work at
  Teuchos::RCP<DRT::Discretization> discret_;

  /// linear algebraic solver
  Teuchos::RCP<LINALG::Solver> solver_;

  /// context for output and restart
  Teuchos::RCP<IO::DiscretizationWriter> output_;

  // another writer
  Teuchos::RCP<IO::OutputControl> output_control_fine_;
  Teuchos::RCP<IO::OutputControl> output_control_coarse_;
  // init writer
  Teuchos::RCP<IO::DiscretizationWriter> output_fine_ ;

  ///  time integrator for structural dynamics
  Teuchos::RCP<STR::TimInt> sti_;

  // prolongator which maps all nodal results to finest grid

  Teuchos::RCP<Epetra_MultiVector> prolongator_disp_;
  // for stresses
  Teuchos::RCP<Epetra_MultiVector> prolongator_stress_;

  // prolongators based on crs
   RCP<Epetra_FECrsMatrix> prolongator_disp_crs_;
   RCP<Epetra_FECrsMatrix> prolongator_stress_crs_;

  // Vetors to store mean stresses and disp
  Teuchos::RCP<Epetra_MultiVector> mean_stress_;
  Teuchos::RCP<Epetra_MultiVector> mean_strain_;
  Teuchos::RCP<Epetra_MultiVector> mean_disp_;
  // Vectors to store standard deviation
  Teuchos::RCP<Epetra_MultiVector> var_stress_;
  Teuchos::RCP<Epetra_MultiVector> var_strain_;
  Teuchos::RCP<Epetra_MultiVector> var_disp_;
  // Vectors needed onthe way to calc standard dev
  Teuchos::RCP<Epetra_MultiVector> delta_stress_;
  Teuchos::RCP<Epetra_MultiVector> delta_strain_;
  Teuchos::RCP<Epetra_MultiVector> delta_disp_;
  Teuchos::RCP<Epetra_MultiVector> m2_var_stress_;
  Teuchos::RCP<Epetra_MultiVector> m2_var_strain_;
  Teuchos::RCP<Epetra_MultiVector> m2_var_disp_;
  Teuchos::RCP<Epetra_MultiVector> m2_helper_var_disp_;
  Teuchos::RCP<Epetra_MultiVector> m2_helper_var_stress_;
  Teuchos::RCP<Epetra_MultiVector> m2_helper_var_strain_;
  // Same Vectors for Difference between levels again
  // Vetors to store mean stresses and disp
  Teuchos::RCP<Epetra_MultiVector> diff_mean_stress_;
  Teuchos::RCP<Epetra_MultiVector> diff_mean_strain_;
  Teuchos::RCP<Epetra_MultiVector> diff_mean_disp_;
  // Vectors to store standard deviation
  Teuchos::RCP<Epetra_MultiVector> diff_var_stress_;
  Teuchos::RCP<Epetra_MultiVector> diff_var_strain_;
  Teuchos::RCP<Epetra_MultiVector> diff_var_disp_;
  // Vectors needed onthe way to calc standard dev
  Teuchos::RCP<Epetra_MultiVector> diff_delta_stress_;
  Teuchos::RCP<Epetra_MultiVector> diff_delta_strain_;
  Teuchos::RCP<Epetra_MultiVector> diff_delta_disp_;
  Teuchos::RCP<Epetra_MultiVector> diff_m2_var_stress_;
  Teuchos::RCP<Epetra_MultiVector> diff_m2_var_strain_;
  Teuchos::RCP<Epetra_MultiVector> diff_m2_var_disp_;
  Teuchos::RCP<Epetra_MultiVector> diff_m2_helper_var_disp_;
  Teuchos::RCP<Epetra_MultiVector> diff_m2_helper_var_stress_;
  Teuchos::RCP<Epetra_MultiVector> diff_m2_helper_var_strain_;
  // vectors to store data from lower level
  Teuchos::RCP<Epetra_MultiVector> stress_lower_level_;
  Teuchos::RCP<Epetra_MultiVector> strain_lower_level_;
  Teuchos::RCP<Epetra_MultiVector> disp_lower_level_;



  // filename of standard control file
  string filename_ ;
  // meshfilename name which is written to controlfiel every run
  // this should be the name of meshfile of first run
  string meshfilename_;

  // Parameters from input file
  // Get number of Newton iterations
   int num_newton_it_;
   // Get convergence tolerance
   double convtol_;

   double InEleRange_ ;

   // number of timesteps
   int tsteps_;

  ///  number of measurment points the
  double tstep_;                    // timestep

   // number of run in MLMC
  int numb_run_;
  int start_run_;
  // starting random seed
  int start_random_seed_;
  // output filename of lower level
  string filename_lower_level_;

  // number of current level
  int num_level_;

  int write_stats_;

  int reset_out_count_;

  bool calc_diff_;
  bool prolongate_res_;
  // the two discretizations
  Teuchos::RCP<DRT::Discretization> actdis_fine_;
  Teuchos::RCP<DRT::Discretization> actdis_coarse_;
  // output
 // Teuchos::RCP<IO::DiscretizationWriter> output_fine_ ;

  ;


public:

  /// standard constructor
  MLMC(
    Teuchos::RCP<DRT::Discretization> dis,  ///< the discretisation
    Teuchos::RCP<LINALG::Solver> solver,  ///< algebraic solver
    Teuchos::RCP<IO::DiscretizationWriter> output  ///< outputer
    );

  // Setup prolongator to transfer results to finest grid
  void SetupProlongator();
  void SetupProlongatorParallel();
  // Map results from one Discretization to another using the prolongator
  void ProlongateResults();
  ///
  void Integrate();

  // Find backgroundelement of a node
  int FindBackgroundElement(DRT::Node node, Teuchos::RCP<DRT::Discretization> background_dis, int* bg_ele_id, double* xsi);

  // check if a node lies within an element
  // returns double based on minimum value of local coordinates abscissa
  double CheckIfNodeInElement(DRT::Node& node ,DRT::Element& ele, double* xsi);

  bool EvaluateF(double* f, DRT::Node& node, DRT::Element& ele,const double* xsi);

  bool EvaluateGradF(LINALG::Matrix<3,3>& fgrad,DRT::Node& node, DRT::Element& ele,const double* xsi);

  void ReadResultsFromLowerLevel();

  void CalcDifferenceToLowerLevel(RCP< Epetra_MultiVector> stress,RCP< Epetra_MultiVector> strain, RCP<Epetra_MultiVector> disp);
  // setup stoch mat law
  void SetupStochMat(unsigned int random_seed);

  void CalcStatStressDisp(RCP<Epetra_MultiVector> curr_stress,RCP<Epetra_MultiVector> curr_strain,RCP<Epetra_MultiVector> curr_disp);

  void WriteStatOutput();

  // Reset prestress F f present in discret
  void ResetPrestress();

  void HelperForDebuggin();

  // Funtion that writes the disp and sigma values of specific DOF to a seperate outputfile
  void HelperFunctionOutput(RCP< Epetra_MultiVector> stress,RCP< Epetra_MultiVector> strain, RCP<Epetra_MultiVector> disp);
  void HelperFunctionOutputTube(RCP< Epetra_MultiVector> stress,RCP< Epetra_MultiVector> strain, RCP<Epetra_MultiVector> disp);

};  // class MLMC

}  // namespace STR

#endif /*MLMC_H_*/
#endif // HAVE_FFTw
#endif  // #ifdef CCADISCRET
