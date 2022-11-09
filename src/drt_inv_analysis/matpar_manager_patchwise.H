/*----------------------------------------------------------------------*/
/*! \file
\brief Creating patches from an elementwise layout

\level 3

*/


/*----------------------------------------------------------------------*/
/* definitions */
#ifndef MATPAR_MANAGER_PERPATCH_H
#define MATPAR_MANAGER_PERPATCH_H

/*----------------------------------------------------------------------*/
/* headers */
#include "matpar_manager_elementwise.H"


// forward declarations
namespace DRT
{
  class Discretization;
}

namespace INVANA
{
  /*! \class MatParManagerPerPatch
   *  \brief Class to have patchwise constant parameters
   *
   *  \author kehl \date 08/2016
   */
  class MatParManagerPerPatch : public MatParManagerPerElement
  {
   public:
    //! constructor
    MatParManagerPerPatch(Teuchos::RCP<DRT::Discretization> discret);

    //! destructor
    virtual ~MatParManagerPerPatch(){};

    //! Setup specific parametrization layout
    virtual void Setup();

    //! Account for distributed chain-rule application
    virtual void Finalize(
        Teuchos::RCP<Epetra_MultiVector> source, Teuchos::RCP<Epetra_MultiVector> target);

    //! Return initial covariance guess
    virtual Teuchos::RCP<Epetra_CrsMatrix> InitialCovariance();

   protected:
    //! Get current material parameters
    virtual void FillParameters(Teuchos::RCP<Epetra_MultiVector> params);

    //! apply projection to INVANA::DcsMatrix and get diagonal
    virtual void ApplyParametrization(
        DcsMatrix& matrix, Teuchos::RCP<Epetra_MultiVector> diagonals);

    //! Initialize optimization parameters
    virtual void InitParameters();

    //! Chain rule application
    virtual void ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint, double val, int elepos,
        int paraposglobal, int paraposlocal);

   private:
    //! @name Setup the parametrization
    //@{

    //! Create projection to patch-wise atoms
    void CreateProjection();
    //@}

    //! @name Helper to setup the parametrization
    //@{

    //! Read elementwise parameter distribution from input
    void ReadMAPApproximation();

    //! Create a dictionary of patch atoms for a specific number of patches
    void CreateLevelDictionary(int patchlevel  // [in] patchlevel to create atoms for
    );

    //! check the approximation power of this dictionary
    double CheckApproximation();

    //! find connectivities within levels
    void FindLevelConnectivity(
        std::map<int, std::vector<int>>& neighbourhood,  // [in] neighbourinformation
        std::map<int, std::vector<int>>& patches         // [out] patches within this level
    );

    //! Make a histogram of the elementwise solution
    void MakeHistogram();

    //! Find bin with maximum number of entries
    void FindHistogramMax(std::map<int, int>& hisbins, int& index);

    //! quality threshold for the approximation
    double qthresh_;

    //! histogram of the elementwise parameters
    std::vector<double> histvalues_;

    //! bins of the histogram
    std::map<int, int> histbins_;

    // file to read the map approximation from
    std::string map_restart_file_;

    // step from which to read the map approximation
    int map_restart_step_;
    //@}

    //! @name Parameters of the initial elementwise distribution
    //@{

    //! elementwise set of optimization parameters
    Teuchos::RCP<Epetra_MultiVector> optparams_elewise_;

    //! elementwise layout
    Teuchos::RCP<Epetra_Map> elewise_map_;

    //! Covariance approximation
    Teuchos::RCP<DcsMatrix> fullcovariance_;
    //@}

    //! @name Projection between sparse appoximation and full basis
    //@{

    //! maps of the atoms of the final dictionary
    Teuchos::RCP<LINALG::MultiMapExtractor> patchmap_;
    //@}

    //! global optimization parameter id to patch id
    std::map<int, int> pidtopatch_;

    //! number of patch levels
    int max_num_levels_;
  };
}  // namespace INVANA

#endif
