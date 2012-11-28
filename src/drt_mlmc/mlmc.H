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
  // to be sorted and/or renamed/ deleted

  //! \brief my structure discretisation to work at
  Teuchos::RCP<DRT::Discretization> discret_;

  //! \brief the fine discretization
  Teuchos::RCP<DRT::Discretization> actdis_fine_;

  //! \brief the coarse discretization
  Teuchos::RCP<DRT::Discretization> actdis_coarse_;

  //! \brief linear algebraic solver
  //Teuchos::RCP<LINALG::Solver> solver_;

  //! \brief context for output and restart
  Teuchos::RCP<IO::DiscretizationWriter> output_;

  //! \brief context for output and restart fine discretization
  Teuchos::RCP<IO::DiscretizationWriter> output_fine_ ;

  //! \brief another writer
  Teuchos::RCP<IO::OutputControl> output_control_fine_;

  //! \brief another writer
  Teuchos::RCP<IO::OutputControl> output_control_coarse_;

  // #########################################################################
  // ##     variables needed for plain  Monte Carlo                         ##
  // #########################################################################

  //! \brief filename of standard control file
  string filename_ ;

  //! \brief meshfilename name which is written to controlfile every run this should be the name of meshfile of first run
  string meshfilename_;

  // Parameters from input file

  //! brief number of Newton iterations
  int num_newton_it_;

  //! brief convergence tolerance for loop in CheckIfNodeInElement routine
  double convtol_;

  //! brief maximum value of xi_i which is still considered to be "within" the element
  double InEleRange_ ;

  //! brief number of monte carlo runs
  int numb_run_;

  //! brief run to start with
  int start_run_;

  //! brief number of mc runs per npar group
  int numruns_pergroup_;

  //! brief starting random seed
  int start_random_seed_;

  //! brief use a global deterministic value as material parameter instead of random field
  bool use_det_value_;

  //! brief global deterministic value that is used instead of random field
  double det_value_;

  //! brief flag to reduce output to files as well as stdout
  bool reduced_output_;

  //! brief randomfield stored in class variable so that we can reuse it in every run
  Teuchos::RCP<GenRandomField> random_field_;

  //! brief vector which contains all element numbers we need the output from
  std::vector <int> AllMyOutputEleIds_;

  //! brief vector stores gids of elements we want the output of specific to proc
  std::vector <int> my_output_elements_;

  //! brief map for communicating the element data to proc 0
  Teuchos::RCP<const Epetra_Map> OutputMap_;



  // #########################################################################
  // ##     variables needed for Multilevel Monte Carlo only                ##
  // #########################################################################

  //! brief output filename of lower level
  string filename_lower_level_;

  //! brief number of current level
  int num_level_;

  //! brief write statistics every write_stats step
  int write_stats_;

  //! brief do we want to calculate difference between the to levels
  bool calc_diff_;

  //! brief do we want to prolongate the results to the lower levels
  bool prolongate_res_;

  //! brief number of timesteps
  int tsteps_;


  //! \brief prolongator which maps displacements to finest grid
  Teuchos::RCP<Epetra_MultiVector> prolongator_disp_;

  //! \brief prolongator which maps stresses to finest grid
  Teuchos::RCP<Epetra_MultiVector> prolongator_stress_;

  //! \brief prolongators based on crs which maps displacements to finest grid
  Teuchos::RCP<Epetra_FECrsMatrix> prolongator_disp_crs_;

  //! \brief prolongators based on crs which maps stresses to finest grid
  Teuchos::RCP<Epetra_FECrsMatrix> prolongator_stress_crs_;

   // vectors needed to calc this and that

  // means
  //! \brief Vector to store mean stresses
  Teuchos::RCP<Epetra_MultiVector> mean_stress_;

  //! \brief Vector to store mean strains
  Teuchos::RCP<Epetra_MultiVector> mean_strain_;

  //! \brief Vector to store mean disp
  Teuchos::RCP<Epetra_MultiVector> mean_disp_;


  // std
  //! \brief Vector to store standard deviation of stress
  Teuchos::RCP<Epetra_MultiVector> var_stress_;

  //! \brief Vector to store standard deviation of strain
  Teuchos::RCP<Epetra_MultiVector> var_strain_;

  //! \brief Vector to store standard deviation of disp
  Teuchos::RCP<Epetra_MultiVector> var_disp_;


  // Vectors needed on the way to calculate std
  //! \brief some helper vector
  Teuchos::RCP<Epetra_MultiVector> delta_stress_;

  //! \brief some helper vector
  Teuchos::RCP<Epetra_MultiVector> delta_strain_;

  //! \brief some helper vector
  Teuchos::RCP<Epetra_MultiVector> delta_disp_;

  //! \brief some helper vector
  Teuchos::RCP<Epetra_MultiVector> m2_var_stress_;

  //! \brief some helper vector
  Teuchos::RCP<Epetra_MultiVector> m2_var_strain_;

  //! \brief some helper vector
  Teuchos::RCP<Epetra_MultiVector> m2_var_disp_;

  //! \brief some helper vector
  Teuchos::RCP<Epetra_MultiVector> m2_helper_var_disp_;

  //! \brief some helper vector
  Teuchos::RCP<Epetra_MultiVector> m2_helper_var_stress_;

  //! \brief some helper vector
  Teuchos::RCP<Epetra_MultiVector> m2_helper_var_strain_;


  // We need all those vectors for the difference between levels again


  // vectors to store mean

  //! \brief Vector to store diff of mean stresses
  Teuchos::RCP<Epetra_MultiVector> diff_mean_stress_;

  //! \brief Vector to store diff of mean strains
  Teuchos::RCP<Epetra_MultiVector> diff_mean_strain_;

  //! \brief Vector to store diff of mean disp
  Teuchos::RCP<Epetra_MultiVector> diff_mean_disp_;


  // vectors to store standard deviation
  //! \brief Vector to store standard deviation diff of stress
  Teuchos::RCP<Epetra_MultiVector> diff_var_stress_;

  //! \brief Vector to store standard deviation diff of strain
  Teuchos::RCP<Epetra_MultiVector> diff_var_strain_;

  //! \brief Vector to store standard deviation diff of disp
  Teuchos::RCP<Epetra_MultiVector> diff_var_disp_;


  // Vectors needed on the way to calculate std

  //! \brief some helper vector
  Teuchos::RCP<Epetra_MultiVector> diff_delta_stress_;

  //! \brief some helper vector
  Teuchos::RCP<Epetra_MultiVector> diff_delta_strain_;

  //! \brief some helper vector
  Teuchos::RCP<Epetra_MultiVector> diff_delta_disp_;

  //! \brief some helper vector
  Teuchos::RCP<Epetra_MultiVector> diff_m2_var_stress_;

  //! \brief some helper vector
  Teuchos::RCP<Epetra_MultiVector> diff_m2_var_strain_;

  //! \brief some helper vector
  Teuchos::RCP<Epetra_MultiVector> diff_m2_var_disp_;

  //! \brief some helper vector
  Teuchos::RCP<Epetra_MultiVector> diff_m2_helper_var_disp_;

  //! \brief some helper vector
  Teuchos::RCP<Epetra_MultiVector> diff_m2_helper_var_stress_;

  //! \brief some helper vector
  Teuchos::RCP<Epetra_MultiVector> diff_m2_helper_var_strain_;

  // vectors to store data from lower level
  //! \brief vector to store stresses from lower level
  Teuchos::RCP<Epetra_MultiVector> stress_lower_level_;

  //! \brief vector to store strains from lower level
  Teuchos::RCP<Epetra_MultiVector> strain_lower_level_;

  //! \brief vector to store disp from lower level
  Teuchos::RCP<Epetra_MultiVector> disp_lower_level_;
  ;

public:

  // #########################################################################
  // ##     functions needed for plain and simple Monte Carlo               ##
  // #########################################################################

  //! \brief standard constructor
  MLMC( Teuchos::RCP<DRT::Discretization> dis,  ///< the discretisation
    Teuchos::RCP<IO::DiscretizationWriter> output  ///< outputer
    );

  //! \brief perform time integration and eval (production version)
  void Integrate();

  //! \brief perform time integration and eval (testing version to try out new stuff, to be merged some day)
  void IntegrateNoReset();

  //! \brief routine to assign material parameters to each element based on random field
  void SetupStochMat(unsigned int random_seed);

  //! \brief reset prestress deformation gradient F present in discretisation
  void ResetPrestress();

  //! \brief  helper function that writes the displacement and stress values of specific DOFs to a separate outputfile
  void HelperFunctionOutput(Teuchos::RCP< Epetra_MultiVector> stress,
      Teuchos::RCP< Epetra_MultiVector> strain,
      Teuchos::RCP<Epetra_MultiVector> disp);

  //! \brief  helper function that writes the displacement and stress values of specific DOFs to a separate outputfile
  void HelperFunctionOutputTube(Teuchos::RCP< Epetra_MultiVector> stress,
      Teuchos::RCP< Epetra_MultiVector> strain,
      Teuchos::RCP<Epetra_MultiVector> disp);

  //! function to evaluate the discretization at certain nodes which are for now hard coded in this function
  void EvalDisAtNodes(Teuchos::RCP<const Epetra_Vector> displacement_of_discretization );

  //! \brief evaluate the discretization at specific elements and write to txt. files vector of elements passed from input file
  void EvalDisAtEleCenters(Teuchos::RCP<const Epetra_Vector> disp,
      INPAR::STR::StressType iostress,
      INPAR::STR::StrainType iostrain);

  //! \brief determine which of the output elemenst are on which proc
  void SetupEvalDisAtEleCenters(std::vector <int> AllOutputEleIds);

  //! \brief more or less a copy of the evaluate function of discretization class, allwos us to evaluate specifi elements only
  void Evaluate2(
                         Teuchos::ParameterList&              params,
                         Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                         Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                         Teuchos::RCP<Epetra_Vector>          systemvector1,
                         Teuchos::RCP<Epetra_Vector>          systemvector2,
                         Teuchos::RCP<Epetra_Vector>          systemvector3);

  // #########################################################################
  // ##     functions needed for Multilevel Monte Carlo only                ##
  // #########################################################################

  // someday these should move to another class

  //! \brief find backgroundelement of a node in the background discretization
  int FindBackgroundElement(DRT::Node node,
      Teuchos::RCP<DRT::Discretization> background_dis,
      int* bg_ele_id, double* xsi);

  //! \brief check if a node lies within a specific element, returns double based on minimum value of local coordinates abscissa
  double CheckIfNodeInElement(DRT::Node& node ,DRT::Element& ele, double* xsi);

  //! \brief helper function
  bool EvaluateF(double* f, DRT::Node& node, DRT::Element& ele,const double* xsi);

  //! \brief helper function
  bool EvaluateGradF(LINALG::Matrix<3,3>& fgrad,DRT::Node& node, DRT::Element& ele,const double* xsi);

  //! \brief calculate some statistics
  void CalcStatStressDisp(Teuchos::RCP<Epetra_MultiVector> curr_stress,
      Teuchos::RCP<Epetra_MultiVector> curr_strain,Teuchos::RCP<Epetra_MultiVector> curr_disp);

  //! \brief Setup prolongator to transfer results to finest grid
  void SetupProlongator();

  //! \brief Setup prolongator to transfer results to finest grid
  void SetupProlongatorParallel();

  //! brief Map results from one Discretization to another using the prolongator
  void ProlongateResults();

  //! brief read in results of lower level from file system
  void ReadResultsFromLowerLevel();

  //! brief calculate difference to lower level discretization
  void CalcDifferenceToLowerLevel(Teuchos::RCP< Epetra_MultiVector> stress,
      Teuchos::RCP< Epetra_MultiVector> strain, Teuchos::RCP<Epetra_MultiVector> disp);

  //! brief write calculated statistics to file
  void WriteStatOutput();


};  // class MLMC

}  // namespace STR

#endif /*MLMC_H_*/
#endif // HAVE_FFTw
