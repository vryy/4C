/*----------------------------------------------------------------------*/
/*!
\file mc_mat_par_manager.H
\brief Class to handle calls to material parameters from a uq routine

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            089 - 289-15276
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef HAVE_FFTW
#ifndef MC_MATPAR_MANAGER_H
#define MC_MATPAR_MANAGER_H

#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_dserror.H"
// forward declarations
namespace DRT
{
  class Discretization;
}

namespace LINALG
{
  class MultiMapExtractor;
}


namespace STR
{

namespace UQ
{

class RandomField;
/*----------------------------------------------------------------------*/
class MCMatParManager
{
public:

  //! Constructor
  MCMatParManager(Teuchos::RCP<DRT::Discretization> discret);

  //! Destructor
  ~MCMatParManager(){};

  //! evaluate derivative of the residual force vector w.r.t to the material parameters using duals
  //void Evaluate(double time, Teuchos::RCP<Epetra_MultiVector> dfint, bool consolidate=false);

  //! brief setup all random field needed for MC
  void SetupRandomFields(unsigned int myseed);

  //brief create new realization of all random fields and push parameters down to material
  void SetUpStochMats(unsigned int myseed, double para_cont_parameter, bool reuse_rf);


private:

  //! brief initialze map that contains <matid> <vector of stochastic params>
  void InitStochParaMap();

  void SetParams(double para_cont_parameter);

  //brief compute new set of random Phase angles/ new realizations of all random fields
  void CreateNewRealizationOfRandomFields(unsigned int myseed);

  void ComputeMatParamsMultivectorFromRandomFields(Teuchos::RCP<Epetra_MultiVector> params, double para_cont_parameter);

  //! brief create new random field based on input data
  Teuchos::RCP<RandomField> CreateRandomField(int random_field_id, unsigned int myseed);

public:


private:

  //! discretization for convenience
  Teuchos::RCP<DRT::Discretization> discret_;

  // make those temporarily public only
  //! map from stochparaid (int) to the random field which models this parameter (RCP<RandomField>)
  std::map<int, Teuchos::RCP<RandomField> > randomfields_;

  //! parametrization layout of the optimization parameters
  //Teuchos::RCP<Epetra_Map> paramlayoutmap_;

  //! a uniqueGIDS version to be used for e.g.: writing restart data
  //Teuchos::RCP<Epetra_Map> paramlayoutmapunique_;

  //! number of stochastic material parameters
  int numstochparams_;

  //! map of names for each material parameter
  //! key: matid, value: value-enum
  //! value is of enum type and can be used to set the respective material parameter directly
  //! e.g. < 1 ; ( youngs, beta)>
  //!      < 2 ; ( k1 , k2 )>
  std::map<int, std::vector<int> > stochparamap_;

  //! map of matids to the id ofthe  stochastic parameter
  //! e.g.  < 1 ; (0 , 1)>
  //!       < 2 ; (2)>
  //! thereby if say matid 1 represents AAAneohook with two stochasitc parameters therby the first is modelled by random field number
  // 0 and the second by random field number 1
  // the second material has only one stochastic parameter which is modelled by field number 2
  std::map<int,std::vector<int> > stochparaid_;

  //! actual set of material parameters in elementwise colmap layout
  Teuchos::RCP<Epetra_MultiVector> params_;


  //! Map Extractor to extract Blocks of parameters of the parameter vector according to the value in paramap_
  //Teuchos::RCP<LINALG::MultiMapExtractor> paramapextractor_;

  //! the two following 'maps' are processor specific: each proc has its own
  //! map keeping correspondence between params in the parameter vector (parameter LID in optparams_) and
  //! eleids (i.e gids of the discretizations ele layout)
  //std::vector<int> paramsLIDtoeleGID_;

  //! map keeping correspondence between eleid (GID) + parameter block id and parameter LID
  //std::map< int, std::vector<int> > eleGIDtoparamsLID_;


};
} // namespace MLMC
} // namespace STR

#endif /* PAT_MATPAR_MANAGER_H */
#endif // FFTW
