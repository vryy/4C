/*----------------------------------------------------------------------*/
/*! \file
\brief Class for performing Multi Level Monte Carlo (MLMC)analysis of structure

\maintainer Jonas Nitzler

\level 2
*/
/*----------------------------------------------------------------------*/

#ifdef HAVE_FFTW
#ifndef MLMC_H
#define MLMC_H
#include <list>
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include "Epetra_FECrsMatrix.h"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_structure.H"

// forward declarations
namespace ADAPTER
{
  class Structure;
}

namespace Teuchos
{
  class ParameterList;
}

namespace DRT
{
  class Discretization;
  class Node;
  class Element;
}  // namespace DRT

namespace LINALG
{
  class SparseMatrix;
  class SparseOperator;
}  // namespace LINALG

namespace IO

{
  class DiscretizationWriter;
  class OutputControl;
  class ErrorFileControl;
}  // namespace IO

namespace UQ
{
  class RandomField;
  class MCMatParManager;
  class MCVarThicknessManager;

  /// Multi Level Monte Carlo Simulation of Structures with
  /// some sort of stochastic property
  class MLMC
  {
   private:
    // to be sorted and/or renamed/ deleted

    //! \brief my structure discretisation to work at
    Teuchos::RCP<DRT::Discretization> discret_;

    //! \brief context for output and restart
    Teuchos::RCP<IO::DiscretizationWriter> output_;

    // #########################################################################
    // ##     variables needed for plain  Monte Carlo                         ##
    // #########################################################################

    //! \brief filename of standard control file
    std::string filename_;

    //! \brief meshfilename name which is written to controlfile every run this should be the name
    //! of meshfile of first run
    std::string meshfilename_;

    // Parameters from input file

    //! brief maximum value of xi_i which is still considered to be "within" the element
    double InEleRange_;

    //! brief number of monte carlo runs
    int numb_run_;

    //! brief run to start with
    int start_run_;

    //! brief number of mc runs per npar group
    int numruns_pergroup_;

    //! brief starting random seed
    int start_random_seed_;

    //! brief how many times do half the stepsize in parameter continuation before we give up
    int cont_num_maxtrials_;

    //! brief flag to reduce output to files as well as stdout
    bool reduced_output_;

    //! do we have stochastic wall thickness
    bool stoch_wall_thickness_;

    //! brief vector which contains all element numbers we need the output from
    std::vector<int> AllMyOutputEleIds_;

    //! brief vector stores gids of elements we want the output of specific to proc
    std::vector<int> my_output_elements_;

    //! brief vector stores gids of elements in the AAA wall
    std::vector<int> my_wall_elements_;

    //! brief vector store location of elements in the AAA wall
    std::vector<int> my_wall_elements_ele_c_loc_;

    //! brief total number of elements in the AAA wall
    int tot_num_wall_elements_;

    //! brief map for communicating the element data to proc 0
    Teuchos::RCP<const Epetra_Map> OutputMap_;

    //! brief compares the second entry of a pair<int,double>. To be passed to the sorting algo
    static bool MyComparePairs(
        const std::pair<int, double>& first, const std::pair<int, double>& second)
    {
      if (first.second > second.second)
        return true;
      else
        return false;
    }
    //! brief do we want to write rv to some text file
    bool write_rv_to_file_;

    //! store some time
    double t0_;
    double t1_;
    double t2_;
    double t3_;
    double t4_;
    double t5_;
    double t6_;

    // #########################################################################
    // ##     variables needed for parameter continuation               ##
    // #########################################################################

    // init variables to store restart information for parameter continuation
    Teuchos::RCP<int> cont_step_;
    Teuchos::RCP<double> cont_time_;
    Teuchos::RCP<Epetra_Vector> cont_disn_;
    Teuchos::RCP<Epetra_Vector> cont_disn_init_;
    Teuchos::RCP<Epetra_Vector> cont_veln_;
    Teuchos::RCP<Epetra_Vector> cont_accn_;
    // this should also have row map layout
    // we need to of those because we need to store and keep eleementdata of initial run with
    // constant beta
    Teuchos::RCP<std::vector<char>> cont_elementdata_init_;
    Teuchos::RCP<std::vector<char>> cont_elementdata_;
    Teuchos::RCP<std::vector<char>> cont_nodedata_init_;
    Teuchos::RCP<std::vector<char>> cont_nodedata_;

    // ###########################################################################
    // ##     variables regarding setup and management of stochastic mat params  ##
    // ###########################################################################

    //! \brief object to manage uncertain material parameters
    Teuchos::RCP<MCMatParManager> my_matpar_manager_;

    //! \brief object to manage uncertain wall thickness
    Teuchos::RCP<MCVarThicknessManager> my_thickness_manager_;

   public:
    // #########################################################################
    // ##     functions needed for plain and simple Monte Carlo               ##
    // #########################################################################

    //! \brief standard constructor
    MLMC(Teuchos::RCP<DRT::Discretization> dis);  ///< the discretisation

    //! \brief perform time integration and eval (production version)
    void Integrate();

    //! \brief perform time integration and eval using parameter continuation
    void IntegrateNoReset();

    //! \brief do not integrate only scale von Mises stress by thickness
    void IntegrateScaleByThickness();


    //! brief compute solution of nonlinear problem using numerical continuation scheme
    int ParameterContinuation(unsigned int num_cont_steps, unsigned int random_seed, bool re_use_rf,
        Teuchos::RCP<ADAPTER::Structure> structadaptor);

    //! \brief routine to assign material parameters to each element based on random field
    void SetupStochMat(unsigned int random_seed);

    // Setup Material Parameters in each element based on deterministic value
    void SetupStochMatDet(double value);

    //! \brief routine to assign material parameters to each element based on random field and
    void BlendStochMat(
        unsigned int random_seed, bool reuse_rf_values, double beta_old, double gamma);

    //! \brief reset prestress deformation gradient F present in discretisation
    void ResetPrestress();

    //! \brief evaluate the discretization at specific elements  vector of elements passed from
    //! input file
    void EvalDisAtEleCenters(Teuchos::RCP<const Epetra_Vector> disp,
        INPAR::STR::StressType iostress, INPAR::STR::StrainType iostrain,
        std::vector<int>* output_elements,
        Teuchos::RCP<std::vector<std::vector<double>>> my_output_elements_c_disp,
        // std::vector <std::vector<double > >   *my_output_elements_c_stresses,
        Teuchos::RCP<std::vector<std::vector<double>>> my_output_elements_c_stresses,
        Teuchos::RCP<std::vector<std::vector<double>>> my_output_elements_c_strains,
        Teuchos::RCP<std::vector<std::vector<double>>> my_output_elements_mat_params);

    //! \brief write max values and 99% quantile values to file
    void ExportPeakStressDataAndWriteToFile(
        Teuchos::RCP<std::list<std::pair<int, double>>> peak_values,
        Teuchos::RCP<std::list<std::pair<int, double>>> quantile99_values);

    //! \brief data of elementsto txt. files vector of elements passed from input file
    void ExportEleDataAndWriteToFile(Teuchos::RCP<const Epetra_Map> OutputMap,
        INPAR::STR::StressType iostress, INPAR::STR::StrainType iostrain,
        std::vector<int>* output_elements,
        Teuchos::RCP<std::vector<std::vector<double>>> my_output_elements_c_disp,
        Teuchos::RCP<std::vector<std::vector<double>>> my_output_elements_c_stresses,
        Teuchos::RCP<std::vector<std::vector<double>>> my_output_elements_c_strains,
        Teuchos::RCP<std::vector<std::vector<double>>> my_output_elements_mat_params);

    //! \brief compute von mises stress on the fly and write to file
    void ExportElementVonMisesStressAndWriteToFile(Teuchos::RCP<const Epetra_Map> OutputMap,
        std::vector<int>* output_elements,
        Teuchos::RCP<std::vector<std::vector<double>>> my_output_elements_von_mises_stresses);

    //! \brief write some info about parametercontinuation scheme to file
    void WriteParamContInfoToFile(Teuchos::RCP<std::vector<int>> paramcont_info);

    //! \brief determine which of the output elemenst are on which proc
    void SetupEvalDisAtEleCenters(std::vector<int> AllOutputEleIds);

    //! \brief more or less a copy of the evaluate function of discretization class, allwos us to
    //! evaluate specific elements only
    void Evaluate2(Teuchos::ParameterList& params, std::vector<int>* eval_elements,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
        Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
        Teuchos::RCP<Epetra_Vector> systemvector3);


    //! brief setup map to communicate element stresses and strains in the wall to proc 0
    void SetupEvalPeakWallStress();

    //! brief compute peak wall stress
    void EvalPeakWallStress(Teuchos::RCP<const Epetra_Vector> disp, INPAR::STR::StressType iostress,
        INPAR::STR::StrainType iostrain);

    //! brief compute von Mises stresses or strains
    void CalcVonMises(Teuchos::RCP<std::vector<std::vector<double>>> input_components,
        Teuchos::RCP<std::vector<double>> output_vM);

    //! brief compute maximum value as well as some quantile of some quantity (over all procs)
    void ComputePeakAndQuantile(Teuchos::RCP<std::vector<double>> values,
        std::vector<int>* corresponding_ele_ids, Teuchos::RCP<std::pair<int, double>> peak,
        Teuchos::RCP<std::pair<int, double>> quantile99);


    // ###########################################################################
    // ##     general helper functions  ##
    // ###########################################################################

    void WriteStdVectorToFile(std::vector<double>, std::string FileNameWithPath);

  };  // class MLMC

}  // namespace UQ


#endif  /*MLMC_H_*/
#endif  // HAVE_FFTw
