/*----------------------------------------------------------------------*/
/*!
\file matpar_manager.H
\brief Class to handle calls to material parameters from an optimization routine

<pre>
Maintainer: Sebastian Kehl
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


  /*!
  \brief A class to handle material parameters for optimization

  The MatParManager handles the mapping between the elemenwise material parameters and
  the vector of material parameters held in the optimization process. Specific choices of
  parametrizations are to be derived from this class (e.g. elementwise distribution,
  uniform distribution, kernel based parametrization, ...). Basic functionality common to
  all parametrizations is given here.

  <b>Material Parameters</b>(MatParManager::params_) are "stored" as Epetra_MultiVector where each vector refers to a
  physical material property (e.g. youngs modulus for St. Venant Kirchhoff, alpha "stiffness"
  and beta "stiffness" in the AAA material) corresponding to one input file material.
  The layout of each vector is given by the underlying discretization's elecolmap. It is not stored but only
  filled with the current material parameters upon request.

  <b>Optimization parameters</b>(MatParManager::optparams_) are stored as Epetra_MultiVector for historical reasons too.
  The respective parametrizations need to take care about the mapping between optparams_ and params_
  The row-wise layout (MatParManager::paramlayoutmap_) must be provided  by the respective parametrizations.
  Parameterizations should only use the first vector in the multivector. Multiple-vector support might go away.

  Core functionality is the evaluation of the derivative of some residual vector wrt to the
  optimization parameters \f$ \frac{\partial \underline{R}}{\underline{p}_o}\f$ and the assembly into the gradient
  \f$ \frac{\partial \pi}{\underline{p}_o}\f$ of some objective function \f$ \pi \f$ via directly evaluating
  \f$ \lambda^T \cdot \frac{\partial \pi}{\underline{p}_o}\f$. With \f$ \lambda \f$ being the lagrange multiplier
  given by the solution of an adjoint problem. Parametrizations are respected via application of
  the chain rule. MatParManager::ContractGradient() therefore needs to be implemented in derived classes. To allow
  for parallel processing, the different procs results are accounted for finally in MatParManager::Consolidate(),
  which is also derived classes business.

 */
class MatParManager
{

public:

  //! constructor
  MatParManager(Teuchos::RCP<DRT::Discretization> discret);

  //!destructor
  virtual ~MatParManager() {}

  //! return number of vectors of the optimization parameter vector layout
  virtual int NumVectors() {return 1;}

  //! init function for base class
  void Init(const Teuchos::ParameterList& invp);

  //! Setup specific parametrization layout for derived classes
  virtual void Setup() = 0;

  //! evaluate derivative of the residual force vector w.r.t to the material parameters using duals
  virtual void AddEvaluate(double time, Teuchos::RCP<Epetra_MultiVector> dfint);

  //! account for partial chain rule application per proc only
  //! and sum up contributions accross procs
  virtual void Finalize(Teuchos::RCP<Epetra_MultiVector> dfint) = 0;

  //! evaluate derivative of the residual force vector w.r.t to the material parameters using FD
  void AddEvaluateFD(double time, Teuchos::RCP<Epetra_MultiVector> dfint);

  //! update optimization parameters
  void ReplaceParams(const Epetra_MultiVector& toreplace);

  //! return material parameters
  Teuchos::RCP<Epetra_MultiVector> GetMatParams();

  //! return optimization parameters for initial guess
  const Epetra_MultiVector& InitialParams(){return *optparams_initial_;}

  //! return the map of the optimization parameter layout
  Teuchos::RCP<Epetra_Map> ParamLayoutMap(){return paramlayoutmap_;}

  //! return the unique map of the optimization parameter layout
  Teuchos::RCP<Epetra_Map> ParamLayoutMapUnique(){return paramlayoutmapunique_;}

  //! return connectivity object with mapextractor and adjacency matrix
  Teuchos::RCP<ConnectivityData> GetConnectivityData();

  //! put vector index of a given material parameter <string> into loc (used for result testing)
  int GetParameterLocation(int eleid, std::string name);

  //! Return optimization parameters
  const Teuchos::RCP<Epetra_MultiVector> GetParams(){return optparams_;}

protected:

  /*----------------------------------------------------------------------*/
  //! functions to be provided by the derived classes

  //! initialize parameters based on parametrization ("local routine")
  virtual void InitParameters(int parapos, double val) = 0;

  //! let the vector of material parameters be filled from optimization parameters
  virtual void FillParameters(Teuchos::RCP<Epetra_MultiVector> params) = 0;

  //! apply chain rule for derivation wrt the optimization parameters
  virtual void ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,
      double val,
      int elepos,
      int paraposglobal,
      int paraposlocal) = 0;

  Teuchos::RCP<DRT::Discretization> Discret(){return discret_;}

  const std::map<int, std::vector<int> >& ParaMap(){return paramap_;}

  const std::map<int,std::vector<int> >& ParaPos(){return parapos_;}

  //! return the material to be optimized for this element, -1 otherwise
  int ElementOptMat(DRT::Element* ele);

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

  /*----------------------------------------------------------------------*/
  //! functions used by the derived classes

  //! Return const pointer to parameters without to call export etc
  const Teuchos::RCP<Epetra_MultiVector> GetMatParamsVec(){return params_;}

  //! Return number of physically different parameters
  int NumParams(){ return numparams_; };

  //! set elementwise material parameters from current set of optimization parameters
  void SetParams();

  //! initialze vector of optimization parameters ("global routine")
  void InitParams();

  //! Return const pointer to parameters without to call export etc
  Teuchos::RCP<Epetra_MultiVector> WriteParamsVec(){return params_;}

  //! inidicates meta parametrization of the material parameters. Hard coded: \f$ p_m = 0.5*p_o^2+0.1\f$
  INPAR::INVANA::StatInvMetaParamsType metaparams_;

private:

  //! setup paramap_, parapos_ from input and check for validity
  void SetupMatOptMap(const Teuchos::ParameterList& invp);

  //! build a graph of elements in the optimization parameter vector
  /*!
    \param In
    paramrowmap - a map from which to build the graph.
    \param In
    blockid - id of the parameter block.
    \param Out
    graph - a reference to itself to be passed to the optimization algorithm.
  */
  virtual void FillAdjacencyMatrix(const Epetra_Map& paramrowmap, Teuchos::RCP<Epetra_CrsMatrix> graph);

  /*! \brief push params_ to the elements
   *
   * set the material parameters stored in params_ to the MAT::Par classes
   * used by the element's materials.
   *
   * \author kehl
   * \date 07/15
   */
  void PushParamsToElements();

  //! map of matids to parameter-ids for the material parameter vector<BR>
  //! e.g.  < 1 ; (0 , 1)> <BR>
  //!       < 2 ; (2 , 3)> <BR>
  //! wrt the above example this means youngs is stored in vector 0, beta in vector 1, k1 in vector 3, etc.
  std::map<int,std::vector<int> > parapos_;

  //! actual set of material parameters in elecolmap layout
  Teuchos::RCP<Epetra_MultiVector> params_;

  //! number of physically different material parameters
  int numparams_;

  //! discretization for convenience
  Teuchos::RCP<DRT::Discretization> discret_;

  //! map of names for each material parameter
  //! key: matid, value: value-enum <BR>
  //! value is of enum type and can be used to set the respective material parameter directly<BR>
  //! e.g. < 1 ; ( youngs, beta)> <BR>
  //!      < 2 ; ( k1 , k2 )>
  std::map<int, std::vector<int> > paramap_;

  //! map of parameter block ids to local parameter block ids, i.e. ids per material <BR>
  //! for  the above examle this would be\n
  //!       < (0 , 0); (1 , 1); (2 , 0);(3 , 1) >
  //! wrt the above example this means youngs is stored in vector 0, beta in vector 1, k1 in vector 3, etc.
  std::map<int,int> paraposGIDtoLID_;
};
}  // namespace INVANA

/*----------------------------------------------------------------------*/
#endif
