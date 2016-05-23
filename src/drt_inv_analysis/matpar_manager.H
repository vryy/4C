/*----------------------------------------------------------------------*/
/*!
\file matpar_manager.H
\brief Class to handle calls to material parameters from an optimization routine

<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifndef MATPAR_MANAGER_H
#define MATPAR_MANAGER_H

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "../drt_inpar/inpar_statinvanalysis.H"

// forward declarations
namespace DRT
{
  class Discretization;
  class Element;
}

namespace LINALG
{
  class MultiMapExtractor;
}

namespace INVANA
{
  /*!
  \brief a container for the connectivity among optimization parameters.

  This class is thought to be constructed upon request by the MatParManager, e.g. to be used by
  the regularization. The main purpose of it is to encapsulate the connectivity information among
  the optimzation parameters.

  \author kehl
   */
class ConnectivityData
{

public:
   // constructor
  ConnectivityData(Teuchos::RCP<LINALG::MultiMapExtractor> mapextractor, Teuchos::RCP<Epetra_CrsMatrix> graph)
  {
    mapextractor_=mapextractor;
    graph_=graph;
  }

  // destructor
  virtual ~ConnectivityData() {}

  Teuchos::RCP<LINALG::MultiMapExtractor> MapExtractor() { return mapextractor_;}

  Teuchos::RCP<Epetra_CrsMatrix> AdjacencyMatrix() { return graph_;}

private:

  // maps of different blocks of parameters
  Teuchos::RCP<LINALG::MultiMapExtractor> mapextractor_;

  //! graph to hold the connectivity of parameters inside a vector of parameters
  Teuchos::RCP<Epetra_CrsMatrix> graph_;

};

  /*! \brief Meta parametrization of material parameters
   *
   * To manage non-negativity/bound constraints at a central place.
   *
   * \author kehl \date 05/2016
 */
class MetaParametrization
{
public:

  //! Constructor
  MetaParametrization()
  {
    metatype_= INPAR::INVANA::stat_inv_meta_none;
  }

  //! Destructor
  virtual ~MetaParametrization(){}

  void SetParametrization(INPAR::INVANA::StatInvMetaParamsType metatype)
  {
    metatype_= metatype;
  }

  //! Material parameter value -> optimization parameter value
  double Material2Meta(double matval);

  //! Optimization parameter value -> material parameter value
  void Meta2Material(Teuchos::RCP<Epetra_MultiVector> meta,
      Teuchos::RCP<Epetra_MultiVector> material);

  //! derivative of material parameter wrt optimization parameter
  double DMaterialDMeta(double metaval);


private:

  INPAR::INVANA::StatInvMetaParamsType metatype_;
};


  /*! \class MatParManager
  \brief A class to handle material parameters for optimization

  The MatParManager handles the mapping between the elemenwise material parameters and
  the vector of material parameters held in the optimization process. Specific choices of
  parametrizations are to be derived from this class (e.g. elementwise distribution,
  uniform distribution, kernel based parametrization, ...). Basic functionality common to
  all parametrizations is given here.

  <b>Material Parameters</b>(MatParManager::params_) are "stored" as Epetra_MultiVector where
  each vector refers to a physical material property (e.g. youngs modulus for St. Venant
  Kirchhoff, alpha "stiffness" and beta "stiffness" in the AAA material) corresponding to
  one input file material. The layout of each vector is given by the underlying discretization's
  elecolmap. It is not stored but only filled with the current material parameters upon request.

  <b>Optimization parameters</b>(MatParManager::optparams_) are stored as Epetra_MultiVector
  for historical reasons too but only the first vector is used. The respective parametrizations
  need to take care about the mapping between optparams_ and params_ The row-wise
  layout (MatParManager::paramlayoutmap_) must be provided  by the respective parametrizations.
  Parameterizations should only use the first vector in the multivector. Multiple-vector support might go away.

  Core functionality is the evaluation of the derivative of some residual vector wrt to the
  optimization parameters \f$ \frac{\partial \underline{R}}{\underline{p}_o}\f$ and the
  assembly into the gradient \f$ \frac{\partial \pi}{\underline{p}_o}\f$ of some objective function
  \f$ \pi \f$ via directly evaluating \f$ \lambda^T \cdot \frac{\partial \pi}{\underline{p}_o}\f$.
  With \f$ \lambda \f$ being the lagrange multiplier given by the solution of an adjoint problem.
  Parametrizations are respected via application of the chain rule. MatParManager::ContractGradient()
  therefore needs to be implemented in derived classes. To allow for parallel processing, the
  different procs results are accounted for finally in MatParManager::Finalize().

  \author kehl \date 2012

 */
class MatParManager
{

public:

  //! @name Construction/destruction
  //@{

  //! Constructor
  MatParManager(Teuchos::RCP<DRT::Discretization> discret);

  //!Destructor
  virtual ~MatParManager() {}

  //@}

  //! @name Initialization
  //@{

  //! Init function for base class
  void Init(const Teuchos::ParameterList& invp);

  //! Setup specific parametrization layout for derived classes
  virtual void Setup() = 0;

  //@}

  //! @name Evaluation of gradients
  //@{

  /*! \brief Evaluation of gradients
   *
   * Evaluate derivative of the residual force vector
   * wrt. to the material parameters using dual variables
   *
   *  \author kehl \date 01/2014
   */
  virtual void AddEvaluate(double time, Teuchos::RCP<Epetra_MultiVector> dfint);

  /*! \brief Evaluation of gradients (FD)
   *
   * Finite difference version of INVANA::AddEvaluate
   *
   *  \author kehl \date 01/2014
   */
  void AddEvaluateFD(double time, Teuchos::RCP<Epetra_MultiVector> dfint);

  //@}


  //! @name Set/Get parameters
  //@{

  //! Update optimization parameters
  void ReplaceParams(const Epetra_MultiVector& toreplace);

  //! Return initial optimization parameters
  const Epetra_MultiVector& InitialParams(){return *optparams_initial_;}

  //@}


  //! return the unique map of the optimization parameter layout
  Teuchos::RCP<Epetra_Map> ParamLayoutMapUnique(){return paramlayoutmapunique_;}

  //! return connectivity object with mapextractor and adjacency matrix
  Teuchos::RCP<ConnectivityData> GetConnectivityData();

  //! @name Queries used for result testing
  //@{

  //! Return material parameters
  Teuchos::RCP<Epetra_MultiVector> GetMatParams();

  //! Return vector index of a given material parameter <string> into loc
  int GetParameterLocation(int eleid, std::string name);

  //@}


protected:

  //! @name Mapping between parameter vectors
  //@{

  /*! \brief Initialize optimization parameters
   *
   * Initialize optimization parameters from the input file
   *
   * \author kehl
   *
   */
  virtual void InitParameters(int parapos, double val) = 0;

  /*! \brief Get current material parameters
   *
   * Let the vector of material parameters be filled according
   * to the current set of optimization parameters
   */
  virtual void FillParameters(Teuchos::RCP<Epetra_MultiVector> params) = 0;

  //! Call INVANA::FillParameters
  void SetParams();

  //! Call INVANA::InitParameters
  void InitParams();

  //@}

  //! @name Help to evaluate the gradient
  //@{

  /*! \brief Chain rule application
   *
   * Apply chain rule for the derivation of the residual
   * wrt the optimization parameters
   *
   * \author kehl
   */
  virtual void ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint, /// [in/out] to be filled
      double val, /// [in] value to be assembled
      int elepos, /// [in] eleid
      int paraposglobal, /// [in] gid
      int paraposlocal /// [in] lid
      ) = 0;

  /*! \brief Finalize gradient evaluation
   *
   * account for partial chain rule application per
   * proc only and sum up contributions accross procs
   *
   * \author kehl
   */
  virtual void Finalize(Teuchos::RCP<Epetra_MultiVector> source, /// [in] distributed source
      Teuchos::RCP<Epetra_MultiVector> target /// [out] finalized on unique layout
      ) = 0;

  //@}

  //! @name Access to maps
  //@{

  //! Get the underlying discretization
  Teuchos::RCP<DRT::Discretization> Discret(){return discret_;}

  //! Get map of input file parameters to be optimized for
  const std::map<int, std::vector<int> >& ParaMap(){return paramap_;}

  //! Get map of the location input file parameters in storage
  const std::map<int,std::vector<int> >& ParaPos(){return parapos_;}

  //! Return number of physically different parameters
  int NumParams(){ return numparams_; };

  //! Return the material to be optimized for this element, -1 otherwise
  int ElementOptMat(DRT::Element* ele);

  //@}

  //! actual set of material parameters in some parametrization layout
  Teuchos::RCP<Epetra_MultiVector> optparams_;

  //! initial guess of material parameters in some parametrization layout
  Teuchos::RCP<Epetra_MultiVector> optparams_initial_;

  //! parametrization layout of the optimization parameters
  Teuchos::RCP<Epetra_Map> paramlayoutmap_;

  //! a uniqueGIDS version to be used for e.g.: writing restart data
  Teuchos::RCP<Epetra_Map> paramlayoutmapunique_;

  //! Map Extractor to extract Blocks of material parameters from the optimization parameter vector
  Teuchos::RCP<LINALG::MultiMapExtractor> paramapextractor_;

  //! inidicates meta parametrization of the material parameters
  INVANA::MetaParametrization metaparams_;

  //! @name To go away
  //@{

  //! Return const pointer to parameters without to call export etc
  const Teuchos::RCP<Epetra_MultiVector> GetMatParamsVec(){return params_;}

  //! Return const pointer to parameters without to call export etc
  Teuchos::RCP<Epetra_MultiVector> WriteParamsVec(){return params_;}

  //@}

private:

  //! @name Setup
  //@{

  //! setup paramap_, parapos_ from input and check for validity
  void SetupMatOptMap(const Teuchos::ParameterList& invp);

  /*! \brief graph of optimization parameters
   *
   * A graph having the connectivity of the
   * optimization parameters. Used for regularization.
  */
  virtual void FillAdjacencyMatrix(const Epetra_Map& paramrowmap, /// [in] a map from which to build the graph
      Teuchos::RCP<Epetra_CrsMatrix> graph /// [out] filled graph
      );

  //@}

  //! @name Maps keeping track of the layout of MatParManager::params_
  //@{

  /*! \brief push params_ to the elements
   *
   * set the material parameters stored in params_ to the MAT::Par classes
   * used by the element's materials.
   *
   * \author kehl
   * \date 07/15
   */
  void PushParamsToElements();

  /*! \brief Material ids to parameter ids
   *
   *  map of matids to parameter-ids for the material parameter vector<BR>
   *  e.g.  < 1 ; (0 , 1)> <BR>
   *        < 2 ; (2 , 3)> <BR>
   *  wrt the above example this means youngs is stored in
   *  vector 0, beta in vector 1, k1 in vector 3, etc.
   *
   */
  std::map<int,std::vector<int> > parapos_;

  /*! \brief Materials to be optimized
   * map of names for each material parameter
   * key: matid, value: value-enum <BR>
   * value is of enum type and can be used to set the respective material parameter directly<BR>
   * e.g. < 1 ; ( youngs, beta)> <BR>
   *      < 2 ; ( k1 , k2 )>
   *
   */
  std::map<int, std::vector<int> > paramap_;

  /*! \brief Map between materials and parameters
   * map of parameter block ids to local parameter block ids, i.e. ids per material <BR>
   * for  the above examle this would be\n
   * < (0 , 0); (1 , 1); (2 , 0);(3 , 1) >
   * wrt the above example this means youngs is stored in vector 0, beta in vector 1, k1 in vector 3, etc.
   *
   */
  std::map<int,int> paraposGIDtoLID_;

  //! number of physically different material parameters
  int numparams_;

  //@}

  //! actual set of material parameters in elecolmap layout
  Teuchos::RCP<Epetra_MultiVector> params_;

  //! discretization for convenience
  Teuchos::RCP<DRT::Discretization> discret_;

};
}  // namespace INVANA

#endif
