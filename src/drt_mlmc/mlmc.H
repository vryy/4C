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
#ifdef HAVE_FFTW
#ifndef MLMC_H
#define MLMC_H
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include "Epetra_FECrsMatrix.h"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_structure.H"


using namespace std;
namespace Teuchos
{
  class ParameterList;
}

namespace DRT
{
  class Discretization;
  class Node;
  class Element;
}

namespace LINALG
{
  class Solver;
  class SparseMatrix;
  class SparseOperator;
}

namespace IO

{
  class DiscretizationWriter;
  class OutputControl;
  class ErrorFileControl;
}
class GenRandomField;

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


  // prolongator which maps all nodal results to finest grid
  Teuchos::RCP<Epetra_MultiVector> prolongator_disp_;
  // for stresses
  Teuchos::RCP<Epetra_MultiVector> prolongator_stress_;

  // prolongators based on crs
   Teuchos::RCP<Epetra_FECrsMatrix> prolongator_disp_crs_;
   Teuchos::RCP<Epetra_FECrsMatrix> prolongator_stress_crs_;

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
  // runs per nested par group
  int local_numruns_;
  int numruns_pergroup_;
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

  // need to store randomfield in class variable so that we can reuse it in every run
  Teuchos::RCP<GenRandomField> random_field_;
  vector <int> AllMyOutputEleIds_;
  // stores gids of elements we want the outpuut of
  vector <int> my_output_elements_;
  // Map for communicating the element data to proc 0
  Teuchos::RCP<const Epetra_Map> OutputMap_;
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

  void CalcDifferenceToLowerLevel(Teuchos::RCP< Epetra_MultiVector> stress,Teuchos::RCP< Epetra_MultiVector> strain, Teuchos::RCP<Epetra_MultiVector> disp);
  // setup stoch mat law
  void SetupStochMat(unsigned int random_seed);

  void CalcStatStressDisp(Teuchos::RCP<Epetra_MultiVector> curr_stress,Teuchos::RCP<Epetra_MultiVector> curr_strain,Teuchos::RCP<Epetra_MultiVector> curr_disp);

  void WriteStatOutput();

  // Reset prestress F f present in discret
  void ResetPrestress();

  void HelperForDebuggin();

  // Funtion that writes the disp and sigma values of specific DOF to a seperate outputfile
  void HelperFunctionOutput(Teuchos::RCP< Epetra_MultiVector> stress,Teuchos::RCP< Epetra_MultiVector> strain, Teuchos::RCP<Epetra_MultiVector> disp);
  void HelperFunctionOutputTube(Teuchos::RCP< Epetra_MultiVector> stress,Teuchos::RCP< Epetra_MultiVector> strain, Teuchos::RCP<Epetra_MultiVector> disp);
  // Function to Evaluate the Discretization at ceratin nodes which are for now hard coded in this function
  void EvalDisAtNodes(Teuchos::RCP<const Epetra_Vector> displacement_of_discretization );
  void EvalDisAtEleCenters(Teuchos::RCP<const Epetra_Vector> disp,INPAR::STR::StressType iostress, INPAR::STR::StrainType iostrain);

  void SetupEvalDisAtEleCenters(vector <int> AllOutputEleIds);



 void Evaluate2(
                         Teuchos::ParameterList&              params,
                         Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                         Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                         Teuchos::RCP<Epetra_Vector>          systemvector1,
                         Teuchos::RCP<Epetra_Vector>          systemvector2,
                         Teuchos::RCP<Epetra_Vector>          systemvector3);
};  // class MLMC

}  // namespace STR

#endif /*MLMC_H_*/
#endif // HAVE_FFTw
