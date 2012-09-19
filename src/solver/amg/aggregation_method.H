/*
 * aggregation_method.H
 *
 *  Created on: May 20, 2010
 *      Author: wiesner
 */

/*!
  \file aggregation_method.H
  \brief basic uncoupled aggregation method
 */


#ifndef AGGREGATION_METHOD_H_
#define AGGREGATION_METHOD_H_

#include "../linalg/linalg_utils.H"

using namespace std;
using namespace Teuchos;

namespace LINALG
{
  /*! \brief implementation of uncoupled coarsening
   *
   *  implementation of uncoupled coarsening for a given Epetra_CrsMatrix
   */
  class AggregationMethod
  {
      static const int AGGR_READY = -11;
      static const int AGGR_NOTSEL = -12;
      static const int AGGR_SELECTED = -13;
      static const int AGGR_BDRY = -15;

      class Aggregate
      {
        public: int id;
          std::vector<int> list;
      };

    public:
      AggregationMethod(FILE* outfile = NULL);
      virtual ~AggregationMethod() {};


      /*!
       *  \brief: performs aggregation process
       *  performs aggregation process for given Epetra_CrsMatrix A and parameters params
       *  using an Uncoupled aggregation method. We're supposing a constant block size over all nodes.
       *
       *  \param const RCP<Epetra_CrsMatrix> A: input matrix
       *  \param ParameterList params: parameter list with parameters for aggregation
       *  \param RCP<Epetra_IntVector>& aggrinfo: returns vector with aggregation info. lives on row map of input matrix A.
       *  \param int& naggregates_local: returns number of built aggregates on current processor
       *  \param const RCP<Epetra_MultiVector>& ThisNS: (input) approximation for nullspace of this level (default: null, only needed for ML)
       *  \return number of built aggregates (over all processors)
       *
       *  The following parameters are used for the aggregation process
       *  - "PDE equations": number of real underlaying PDE equations of the non-amalgamated matrix, e.g. 3 for 2D Navier-Stokes examples (2 vel dofs and 1 pressure dof per node).
       *    This is used for the Amalgamation_Matrix routine.
       *  - "Unamalgamated BlockSize": number of degrees of freedom per "node" in the amalgamated map (constant block size). In many cases this should be the same as "PDE equations". However for the AMG(Braess-Sarazin)
       *    coarsening this must be set to "PDE equations" - 1.
       *  - "phase 1: min nodes per aggregate": minimal allowed number of nodes that can build an aggregate (default 9)
       *  - "phase 1: max neighbour nodes": maximal allowed number of neighbour nodes, that already has been added to (other) aggregates. The idea of this parameter is
       *    that such nodes should be handled in phase 2 of aggregation process, where they are attached to existing aggregates. (default 2)
       *  - "phase 2: node attachement scheme": attachement scheme for non-aggregated nodes in phase 2. Can be either "MinRank" or "MaxLink".
       *    With "MinRank" the node is attached to the smallest neighbour aggregate. With "MaxLink" the node is added to the aggregate with the maximum number of links to it.
       */
      virtual int GetGlobalAggregates(const RCP<Epetra_CrsMatrix>& A, ParameterList& params, RCP<Epetra_IntVector>& aggrinfo, int& naggregates_local,const RCP<Epetra_MultiVector>& ThisNS = null);

      /*!
       *  \brief exports Epetra_IntVector in matlab format to file
       *  basic debug routine.
       */
      virtual void PrintIntVectorInMatlabFormat(std::string fname, const Epetra_IntVector& V, const bool newfile = true);

      /*!
       *  \brief return number of dirichlet blocks
       *  returns number of dirichlet blocks, that is detected during amalgamation process. (global number over all processors)
       *  a dirichlet block corresponds to one row in the amalgamated matrix.
       */
      virtual int getNumGlobalDirichletBlocks() {return nGlobalDirichletBlocks; };

    protected:

      int nGlobalDirichletBlocks;         ///< number of Dirichlet boundary blocks
      int nVerbose_;                      ///< verbosity level
  };

  class AggregationMethodFactory
  {
    public:
      AggregationMethodFactory() {};
      ~AggregationMethodFactory() {};

      static RCP<AggregationMethod> Create(const string AggregationMethodType, FILE* outfile = NULL);
  };
}

#endif /* AGGREGATION_METHOD_H_ */
