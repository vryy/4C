/*----------------------------------------------------------------------*/
/*! \file
\brief Class to handle calls to material parameters from an optimization routine

\level 3

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
#include "DcsMatrix.H"
#include "objective_funct.H"

// forward declarations
namespace DRT
{
  class Discretization;
  class Element;
}  // namespace DRT

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
    ConnectivityData(Teuchos::RCP<Epetra_CrsMatrix> graph, Teuchos::RCP<Epetra_CrsMatrix> projector)
    {
      graph_ = graph;
      projector_ = projector;
    }

    // destructor
    virtual ~ConnectivityData() {}

    Teuchos::RCP<Epetra_CrsMatrix> AdjacencyMatrix() { return graph_; }

    Teuchos::RCP<Epetra_CrsMatrix> Projector() { return projector_; }

   private:
    //! graph to hold the connectivity of parameters inside a vector of parameters
    Teuchos::RCP<Epetra_CrsMatrix> graph_;

    /*! \projector
     * a projection operator from some optimization parameter layout
     */
    Teuchos::RCP<Epetra_CrsMatrix> projector_;
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
    MetaParametrization() { metatype_ = INPAR::INVANA::stat_inv_meta_none; }

    //! Destructor
    virtual ~MetaParametrization() {}

    void SetParametrization(INPAR::INVANA::StatInvMetaParamsType metatype) { metatype_ = metatype; }

    //! Material parameter value -> optimization parameter value
    double Material2Meta(double matval);

    //! Optimization parameter value -> material parameter value
    void Meta2Material(
        Teuchos::RCP<Epetra_MultiVector> meta, Teuchos::RCP<Epetra_MultiVector> material);

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
  Parameterizations should only use the first vector in the multivector. Multiple-vector support
  might go away.

  Core functionality is the evaluation of the derivative of some residual vector wrt to the
  optimization parameters \f$ \frac{\partial \underline{R}}{\underline{p}_o}\f$ and the
  assembly into the gradient \f$ \frac{\partial \pi}{\underline{p}_o}\f$ of some objective function
  \f$ \pi \f$ via directly evaluating \f$ \lambda^T \cdot \frac{\partial \pi}{\underline{p}_o}\f$.
  With \f$ \lambda \f$ being the lagrange multiplier given by the solution of an adjoint problem.
  Parametrizations are respected via application of the chain rule.
  MatParManager::ContractGradient() therefore needs to be implemented in derived classes. To allow
  for parallel processing, the different procs results are accounted for finally in
  MatParManager::Finalize().

  \author kehl \date 2012

 */
  class MatParManager
  {
   public:
    //! @name Construction/destruction
    //@{

    //! Constructor
    MatParManager(Teuchos::RCP<DRT::Discretization> discret);

    //! Destructor
    virtual ~MatParManager() {}

    //@}

    //! @name Initialization
    //@{

    //! Init function for base class
    void Init(const Teuchos::ParameterList& invp,     /// [in] inverse analysis input
        Teuchos::RCP<INVANA::ObjectiveFunct> objfunc  /// [in] objective function pointer
    );

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

    //! Update optimization parameters shallow (don't push to discretization)
    void ReplaceParamsShallow(const Epetra_MultiVector& toreplace);

    //! Return initial optimization parameters
    const Epetra_MultiVector& InitialParams() { return *optparams_initial_; }

    //! Return initial guess covariance
    virtual Teuchos::RCP<Epetra_CrsMatrix> InitialCovariance() = 0;

    //@}


    //! return the unique map of the optimization parameters layout
    Teuchos::RCP<Epetra_Map> ParamLayoutMapUnique() { return paramlayoutmapunique_; }

    //! return connectivity object with mapextractor and adjacency matrix
    Teuchos::RCP<ConnectivityData> GetConnectivityData();

    //! @name Queries for result testing, output, etc.
    //@{

    //! Get material parameters from the parametrization
    Teuchos::RCP<Epetra_MultiVector> GetMatParams();

    //! Get 'raw' optimization parameters mapped to the elements
    Teuchos::RCP<Epetra_MultiVector> GetRawParams();

    //! Get projected matrix diagonal from an INVANA::DcsMatrix
    Teuchos::RCP<Epetra_MultiVector> GetMatrixDiagonal(DcsMatrix& matrix);

    //! Return vector index of a given material parameter <string> into loc
    int GetParameterLocation(int eleid, std::string name);
    //@}


   protected:
    //! @name Mapping between parameter vectors
    //@{

    /*! \brief Initialize optimization parameters
     *
     * Initialize optimization parameters as the actual
     * parametrization suggests.
     *
     * \author kehl
     */
    virtual void InitParameters(int parapos, double val) = 0;

    /*! \brief Get current material parameters
     *
     * Let a vector of element-wise row-layout be filled
     * according to the projection of the of the actual
     * optimization parameters
     */
    virtual void FillParameters(Teuchos::RCP<Epetra_MultiVector> params) = 0;

    /*! \brief Get current material parameters
     *
     * Let a vector of element-wise row-layout be filled
     * according to the projection of the actual parametrization
     * to an INVANA::DcsMatrix. And return the diagonal.
     */
    virtual void ApplyParametrization(
        DcsMatrix& matrix, Teuchos::RCP<Epetra_MultiVector> diagonals) = 0;

    //! Set material parameters to the discretization
    void SetParams();

    //! Initialize optimization parameters
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
    virtual void ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,  /// [in/out] to be filled
        double val,         /// [in] value to be assembled
        int elepos,         /// [in] eleid
        int paraposglobal,  /// [in] gid
        int paraposlocal    /// [in] lid
        ) = 0;

    /*! \brief Finalize gradient evaluation
     *
     * account for partial chain rule application per
     * proc only and sum up contributions accross procs
     *
     * \author kehl
     */
    virtual void Finalize(Teuchos::RCP<Epetra_MultiVector> source,  /// [in] distributed source
        Teuchos::RCP<Epetra_MultiVector> target  /// [out] finalized on unique layout
        ) = 0;

    //@}

    //! @name Access to maps
    //@{

    //! Get the underlying discretization
    Teuchos::RCP<DRT::Discretization> Discret() { return discret_; }

    //! Get map of input file parameters to be optimized for
    const std::map<int, std::vector<int>>& ParaMap() { return paramap_; }

    //! Get map of the location input file parameters in storage
    const std::map<int, std::vector<int>>& ParaPos() { return parapos_; }

    //! Return number of physically different parameters
    int NumParams() { return numparams_; };

    //! Return the material to be optimized for this element, -1 otherwise
    int ElementOptMat(DRT::Element* ele);

    //! Elementwise connectivity graph
    virtual Teuchos::RCP<Epetra_CrsMatrix> Graph() = 0;

    //! Projection from full to reduced
    virtual Teuchos::RCP<Epetra_CrsMatrix> Projector() = 0;
    //@}

    //! @name Access to inverse analysis stuff
    //@{

    //! Access to the input file parameters
    const Teuchos::ParameterList& Inpar() { return inpar_invana_; }

    INVANA::ObjectiveFunct& Objfunct() { return *objfunct_; }
    //@}

    //! actual set of material parameters in some parametrization layout
    Teuchos::RCP<Epetra_MultiVector> optparams_;

    //! initial guess of material parameters in some parametrization layout
    Teuchos::RCP<Epetra_MultiVector> optparams_initial_;

    //! parametrization layout of the optimization parameters
    Teuchos::RCP<Epetra_Map> paramlayoutmap_;

    //! a uniqueGIDS version to be used for e.g.: writing restart data
    Teuchos::RCP<Epetra_Map> paramlayoutmapunique_;

    //! Map Extractor to extract Blocks of material parameters from the optimization parameter
    //! vector
    Teuchos::RCP<LINALG::MultiMapExtractor> paramapextractor_;

    //! inidicates meta parametrization of the material parameters
    INVANA::MetaParametrization metaparams_;

    //!@name  Communication utilities
    //@{

    //! The local problems communicator
    const Epetra_Comm& Comm();

    //! The overall global Communicator (e.g. in case of nested parallelity)
    const Epetra_Comm& GComm();
    //@}


    //! @name To go away
    //@{

    //! Return const pointer to parameters without to call export etc
    const Teuchos::RCP<Epetra_MultiVector> GetMatParamsVec() { return params_; }

    //! Return const pointer to parameters without to call export etc
    Teuchos::RCP<Epetra_MultiVector> WriteParamsVec() { return params_; }

    //@}

   private:
    //! @name Setup
    //@{

    //! setup paramap_, parapos_ from input and check for validity
    void SetupMatOptMap(const Teuchos::ParameterList& invp);
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
    std::map<int, std::vector<int>> parapos_;

    /*! \brief Materials to be optimized
     * map of names for each material parameter
     * key: matid, value: value-enum <BR>
     * value is of enum type and can be used to set the respective material parameter directly<BR>
     * e.g. < 1 ; ( youngs, beta)> <BR>
     *      < 2 ; ( k1 , k2 )>
     *
     */
    std::map<int, std::vector<int>> paramap_;

    /*! \brief Map between materials and parameters
     * map of parameter block ids to local parameter block ids, i.e. ids per material <BR>
     * for  the above examle this would be\n
     * < (0 , 0); (1 , 1); (2 , 0);(3 , 1) >
     * wrt the above example this means youngs is stored in vector 0, beta in vector 1, k1 in vector
     * 3, etc.
     *
     */
    std::map<int, int> paraposGIDtoLID_;

    //! number of physically different material parameters
    int numparams_;

    //@}

    //! actual set of material parameters in elecolmap layout
    Teuchos::RCP<Epetra_MultiVector> params_;

    //! discretization for convenience
    Teuchos::RCP<DRT::Discretization> discret_;

    //! Input file parameters
    Teuchos::ParameterList inpar_invana_;

    //! pointer to the objective function
    Teuchos::RCP<INVANA::ObjectiveFunct> objfunct_;
  };
}  // namespace INVANA

#endif
