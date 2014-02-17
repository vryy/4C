/*
 * aggregation_method_ml.H
 *
 *  Created on: May 28, 2010
 *      Author: wiesner
 */

#ifndef AGGREGATION_METHOD_ML_H_
#define AGGREGATION_METHOD_ML_H_

#include "aggregation_method.H"

#include "../../drt_lib/drt_dserror.H"

#include "MLAPI_Error.h"
#include "MLAPI_CompObject.h"
#include "MLAPI_TimeObject.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Operator_Utils.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_InverseOperator.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_BaseOperator.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_Aggregation.h"
#include "MLAPI_Eig.h"

using namespace Teuchos;

namespace LINALG
{

  /*! \brief implementation of uncoupled coarsening (based on ML)
   *
   *  implementation of uncoupled coarsening for a given Epetra_CrsMatrix.
   *  uses the ML Uncoupled Aggregation
   */
  class AggregationMethod_ML : public AggregationMethod
  {
    public:
      AggregationMethod_ML(FILE* outfile = NULL);
      virtual ~AggregationMethod_ML() {};


      /*!
       *  \brief: performs aggregation process
       *  performs aggregation process for given Epetra_CrsMatrix A and parameters params
       *  using an Uncoupled aggregation method. We're supposing a constant block size over all nodes.
       *
       *  \param const Teuchos::RCP<Epetra_CrsMatrix> A: input matrix
       *  \param ParameterList params: parameter list with parameters for aggregation
       *  \param Teuchos::RCP<Epetra_IntVector>& aggrinfo: returns vector with aggregation info. lives on row map of input matrix A.
       *  \param int& naggregates_local: returns number of built aggregates on current processor
       *  \param const Teuchos::RCP<Epetra_MultiVector>& ThisNS: null space of current level, needed for ML
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
      virtual int GetGlobalAggregates(const Teuchos::RCP<Epetra_CrsMatrix>& A, ParameterList& params, Teuchos::RCP<Epetra_IntVector>& aggrinfo, int& naggregates_local, const Teuchos::RCP<Epetra_MultiVector>& ThisNS = null);

      /*!
       *  \brief return number of dirichlet blocks
       *  returns number of dirichlet blocks, that is detected during amalgamation process. (global number over all processors)
       *  a dirichlet block corresponds to one row in the amalgamated matrix.
       */
      virtual int getNumGlobalDirichletBlocks() {return nGlobalDirichletBlocks; };

    private:
      /*!
      \brief Call ML aggregation on A according to parameters supplied in List. Return
             aggregates in aggrinfo

             On input, map of aggrinfo has to map row map of A. On output, aggrinfo[i]
             contains number of aggregate the row belongs to, where aggregates are
             numbered starting from 0.
             Return value is the processor-local number of aggregates build.
             If aggrinfo[i] < 0 the dofs are not considered by the aggregation process.
             This can happen since ML Uncoupled Coarsening uses some MIS routines, that
             allow non-aggregated nodes. These are put together in the aggregate with number -1

      \param A (in): Epetra_CrsMatrix to be aggregated on
      \param List (in): ParameterList containing ML options
      \param ThisNS (in): MultiVector with null space of current level
      \param aggrinfo(out): vector containing aggregation information

      \note Map of aggrinfo has to match rowmap of A on input.

      \return returns processor-local number of aggregates

      \author Michael Gee (gee@lnm.mw.tum.de)
      */
      int GetAggregates(const Teuchos::RCP<Epetra_CrsMatrix>& A, ParameterList& List, const Epetra_MultiVector& ThisNS, Teuchos::RCP<Epetra_IntVector>& aggrinfo);


    protected:
      int nGlobalDirichletBlocks;         ///< number of Dirichlet boundary blocks
  };
}

#endif /* AGGREGATION_METHOD_ML_H_ */
