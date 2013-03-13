/*
 * aggregation_method_uncoupled.H
 *
 *  Created on: May 28, 2010
 *      Author: wiesner
 */


#ifndef AGGREGATION_METHOD_UNCOUPLED_H_
#define AGGREGATION_METHOD_UNCOUPLED_H_

#include "aggregation_method.H"

using namespace Teuchos;

namespace LINALG
{
  /*! \brief implementation of uncoupled coarsening
   *
   *  implementation of uncoupled coarsening for a given Epetra_CrsMatrix
   */
  class AggregationMethod_Uncoupled : public AggregationMethod
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
      AggregationMethod_Uncoupled(FILE* outfile = NULL);
      virtual ~AggregationMethod_Uncoupled() {};


      /*!
       *  \brief: performs aggregation process
       *  performs aggregation process for given Epetra_CrsMatrix A and parameters params
       *  using an Uncoupled aggregation method. We're supposing a constant block size over all nodes.
       *
       *  \param const RCP<Epetra_CrsMatrix> A: input matrix
       *  \param ParameterList params: parameter list with parameters for aggregation
       *  \param RCP<Epetra_IntVector>& aggrinfo: returns vector with aggregation info. lives on row map of input matrix A.
       *  \param int& naggregates_local: returns number of built aggregates on current processor
       *  \param const RCP<Epetra_MultiVector>& ThisNS: null space of current level -> NOT NEEDED -> set to null
       *  \return number of built aggregates (over all processors)
       *
       *  The following parameters are used for the aggregation process
       *  - "Unamalgamated BlockSize": number of unamalgamated original block size of Input Matrix.
       *    should be the number of real underlaying PDE equations of the non-amalgamated matrix, e.g. 3 for 2D Navier-Stokes examples (2 vel dofs and 1 pressure dof per node).
       *    This is used for the Amalgamation_Matrix routine.
       *  - "PDE equations": number of real underlaying PDE equations of the non-amalgamated matrix, e.g. 3 for 2D Navier-Stokes examples (2 vel dofs and 1 pressure dof per node).
       *    This is used for the Amalgamation_Matrix routine. Represents the number of degrees of freedom per "node" in the amalgamated map (constant block size).
       *    In many cases this should be the same as "Unamalgamated BlockSize". However for the AMG(Braess-Sarazin)
       *    coarsening this must be set to "Unamalgamated BlockSize" - 1.
       *  - "phase 1: min nodes per aggregate": minimal allowed number of nodes that can build an aggregate (default 9)
       *  - "phase 1: max neighbour nodes": maximal allowed number of neighbour nodes, that already has been added to (other) aggregates. The idea of this parameter is
       *    that such nodes should be handled in phase 2 of aggregation process, where they are attached to existing aggregates. (default 2)
       *  - "phase 2: node attachement scheme": attachement scheme for non-aggregated nodes in phase 2. Can be either "MinRank" or "MaxLink".
       *    With "MinRank" the node is attached to the smallest neighbour aggregate. With "MaxLink" the node is added to the aggregate with the maximum number of links to it.
       *  - "aggregation method": controls the pruning algorithm. "isotropic aggregation": eliminate all zero entries in matrix (plain pruning), "anisotropic aggregation": uncoupled anisotropic aggregation, "anisotropic aggregation (coupled)": expensive variant of anisotropic aggregation that tries to do anisotropic aggregation over processor boundaries.
       *  - "anisotropic aggregation: epsilon": threshold parameter for anisotropic aggregation. used for "anisotropic aggregation" and "anisotropic aggregation (coupled)".
       *  - "Current Level": number of current level, should be provided of multigrid algorithm. only used for anisotropic aggregation process now. 0 = finest level.
       */
      virtual int GetGlobalAggregates(const RCP<Epetra_CrsMatrix>& A, ParameterList& params, RCP<Epetra_IntVector>& aggrinfo, int& naggregates_local, const RCP<Epetra_MultiVector>& ThisNS = null);

      /*!
       *  \brief return number of dirichlet blocks
       *  returns number of dirichlet blocks, that is detected during amalgamation process. (global number over all processors)
       *  a dirichlet block corresponds to one row in the amalgamated matrix.
       */
      virtual int getNumGlobalDirichletBlocks() {return nGlobalDirichletBlocks; };
    private:
      /*! \brief Amalgamation routine for given Epetra_CrsMatrix A
       *
       * \param A (in): Epetra_CrsMatrix for amalgamation.
       * \param params (in): parameter list
       * \param amalA (out): Graph of amalgamated matrix. should be null on input
       * \param bdry_array (out): Epetra_IntVector with boundary block information (same row map as amalA)
       * \return 0
       *
       * The input matrix must be a quadratic Epetra_CrsMatrix. You have to provide the "PDE equations" parameter and the
       * "BlockSize" parameter in the parameter list params. "PDE equations" is the number of real PDE equations, whereas
       * "BlockSize" denotes the local number of dofs per node. In most cases there is "PDE equatons" equal to "BlockSize".
       * However for the SaddlePoint preconditioner coarsening, we have to choose "PDE equations" to be the real number of
       * PDE equations in A (eg. 3 for 2D) and nPDE-1 for "BlockSize" on the finest level. For the other grid levels you
       * have to choose "PDE equations" = "BlockSize" = 2 (in 2D), since there are no holes for the pressure dofs in the row map
       * of A11 on the coarser levels.
       */
      int AmalgamateMatrix(const RCP<Epetra_CrsMatrix>& A, ParameterList& params, RCP<Epetra_CrsGraph>& amalA, RCP<Epetra_IntVector>& bdry_array);

      /*!
       * \brief Simple Matrix pruning method
       *
       * All matrix graph entries with a zero value are eliminated
       *
       * \param A (in): Epetra_CrsMatrix for pruning
       * \param prunedA (out): Epetra_CrsMatrix with pruned matrix
       * \return 0
       */
      int PruneMatrix(const RCP<Epetra_CrsMatrix>& A, RCP<Epetra_CrsMatrix>& prunedA);

      /*!
       * \brief Matrix pruning method with anisotropic aggregation algorithm
       *
       * All matrix graph entries with a zero value are eliminated
       * Furthermore all entries of matrix A with a(i,j) < epsilon \sqrt{a(i,i)*a(j,j)} are eliminated.
       * The diagonal entries of A (a(i,i)) are stored global on all processors, that is we can do anisotropic aggregation for processor overlapping aggregates, too. Is somewhat expensive as additional inter processor communication is necessary.
       *
       * \param A (in): Epetra_CrsMatrix for pruning
       * \param prunedA (out): Epetra_CrsMatrix with pruned matrix
       * \param epsilon (in): threshold parameter for anisotropic aggregation
       * \param curlevel (in): level number of current multigrid level (0=finest level). used for adaption of epsilon -> epsilon = epsilon * pow(0.5,curlevel)
       * \return 0
       */
      int PruneMatrixAnisotropic(const RCP<Epetra_CrsMatrix>& A, RCP<Epetra_CrsMatrix>& prunedA, double epsilon, int curlevel = 0);

      /*!
       * \brief Matrix pruning method with anisotropic aggregation algorithm
       *
       * All matrix graph entries with a zero value are eliminated
       * Furthermore all entries of matrix A with a(i,j) < epsilon \sqrt{a(i,i)*a(j,j)} are eliminated.
       * If a(j,j) doesn't belong to same processor as a(i,i) the anisotropic aggregation is just ignored. Then simple pruning is applied.
       *
       * \param A (in): Epetra_CrsMatrix for pruning
       * \param prunedA (out): Epetra_CrsMatrix with pruned matrix
       * \param epsilon (in): threshold parameter for anisotropic aggregation
       * \param curlevel (in): level number of current multigrid level (0=finest level). used for adaption of epsilon -> epsilon = epsilon * pow(0.5,curlevel)
       * \return 0
       */
      int PruneMatrixAnisotropicSimple(const RCP<Epetra_CrsMatrix>& A, RCP<Epetra_CrsMatrix>& prunedA, double epsilon, int curlevel = 0);

      /*! \brief unamalgamtes vector
       *  build unamalgamated vector
       *
       *  \param amal_Aggregates: amalgamated vector with aggregate information (lives on amalgamated map amal_map_)
       *  \return Aggregates: unamalgamated vector with aggregate information (lives on map_)
       */
      RCP<Epetra_IntVector> UnamalgamateVector(const RCP<Epetra_IntVector>& amal_Aggregates);

      /*!
       *  \brief this routine controls the coarsening process
       *
       *  coarsening process
       *
       *  \param const RCP<Epetra_CrsGraph>& amalA: the graph of the amalgamated matrix that is to be coarsened
       *  \param const RCP<Epetra_IntVector>& bdry_array: vector array with information whether node is a dirichlet boundary node or not
       *  \param ParameterList& params: parameter list (see GetGlobalAggregates for more information about the accepted parameters)
       *  \param RCP<Epetra_IntVector>& amal_Aggregates: returns the local aggregation information (the local aggregate ids for the nodes)
       *  \param int& nLocalAggregates: returns the local number of nodes, that have been found on this processor
       *  \return 0, if successful
       *
       *
       */
      int Coarsen(const RCP<Epetra_CrsGraph>& amalA, const RCP<Epetra_IntVector>& bdry_array, ParameterList& params, RCP<Epetra_IntVector>& amal_Aggregates, int& nLocalAggregates);

      /*! \brief phase 1 of uncoupled coarsening
       * Phase 1 of uncoupled coarsening. build as many as possible tentative aggregates using the "phase 1" parameters
       *
       * \param const RCP<Epetra_CrsGraph>& amalA: amalgamated graph
       * \param ParameterList& params: parameter list. see GetGlobalAggregates for more information about allowed parameters
       * \param RCP<Epetra_IntVector>& amal_Aggregates: returns vector with aggregation information (local aggid for each node on current proc. -1, if node is not aggregated yet)
       * \param RCP<Epetra_IntVector>& aggr_stat: stores internal information about the aggregation status of each node
       * \param int& nLocalAggregates: number of local aggregates that has been built on current processor so far.
       * \return number of nodes, that are not aggregated (and also no dirichlet boundary nodes!)
       */
      int Phase1(const RCP<Epetra_CrsGraph>& amalA, ParameterList& params, RCP<Epetra_IntVector>& amal_Aggregates, RCP<Epetra_IntVector>& aggr_stat, int& nLocalAggregates);
      /*! \brief phase 2 of uncoupled coarsening
       * Phase 2 of uncoupled coarsening. use MinRank node attachement method
       *
       * \param const RCP<Epetra_CrsGraph>& amalA: amalgamated graph
       * \param ParameterList& params: parameter list. see GetGlobalAggregates for more information about allowed parameters
       * \param RCP<Epetra_IntVector>& amal_Aggregates: returns vector with aggregation information (local aggid for each node on current proc. -1, if node is not aggregated yet)
       * \param RCP<Epetra_IntVector>& aggr_stat: stores internal information about the aggregation status of each node
       * \param int& nLocalAggregates: number of local aggregates that has been built on current processor so far.
       * \return number of nodes, that are not aggregated (and also no dirichlet boundary nodes!)
       */
      int Phase2_minrank(const RCP<Epetra_CrsGraph>& amalA, ParameterList& params, RCP<Epetra_IntVector>& amal_Aggregates, RCP<Epetra_IntVector>& aggr_stat, int& nLocalAggregates);
      /*! \brief phase 2 of uncoupled coarsening
       * Phase 2 of uncoupled coarsening. use MaxLink as node attachment method
       *
       * \param const RCP<Epetra_CrsGraph>& amalA: amalgamated graph
       * \param ParameterList& params: parameter list. see GetGlobalAggregates for more information about allowed parameters
       * \param RCP<Epetra_IntVector>& amal_Aggregates: returns vector with aggregation information (local aggid for each node on current proc. -1, if node is not aggregated yet)
       * \param RCP<Epetra_IntVector>& aggr_stat: stores internal information about the aggregation status of each node
       * \param int& nLocalAggregates: number of local aggregates that has been built on current processor so far.
       * \return number of nodes, that are not aggregated (and also no dirichlet boundary nodes!)
       */
      int Phase2_maxlink(const RCP<Epetra_CrsGraph>& amalA, ParameterList& params, RCP<Epetra_IntVector>& amal_Aggregates, RCP<Epetra_IntVector>& aggr_stat, int& nLocalAggregates);
      /*! \brief phase 3 of uncoupled coarsening
       * Phase 3 of uncoupled coarsening. put all left non-Dirichlet boundary nodes into (one-node) aggregates.
       * The method should return 0, as all nodes have to be aggregated now (unless Dirichlet boundary nodes)
       *
       * \param const RCP<Epetra_CrsGraph>& amalA: amalgamated graph
       * \param ParameterList& params: parameter list. see GetGlobalAggregates for more information about allowed parameters
       * \param RCP<Epetra_IntVector>& amal_Aggregates: returns vector with aggregation information (local aggid for each node on current proc. -1, if node is not aggregated yet)
       * \param RCP<Epetra_IntVector>& aggr_stat: stores internal information about the aggregation status of each node
       * \param int& nLocalAggregates: number of local aggregates that has been built on current processor so far.
       * \return number of nodes, that are not aggregated (and also no dirichlet boundary nodes!)
       */
      int Phase3(const RCP<Epetra_CrsGraph>& amalA, ParameterList& params, RCP<Epetra_IntVector>& amal_Aggregates, RCP<Epetra_IntVector>& aggr_stat, int& nLocalAggregates);
      /*! \brief phase 4 of uncoupled coarsening
       * Phase 4 of uncoupled coarsening. Put all Dirichlet boundary nodes into one-node aggregates.
       *
       * \param const RCP<Epetra_CrsGraph>& amalA: amalgamated graph
       * \param ParameterList& params: parameter list. see GetGlobalAggregates for more information about allowed parameters
       * \param RCP<Epetra_IntVector>& amal_Aggregates: returns vector with aggregation information (local aggid for each node on current proc. -1, if node is not aggregated yet)
       * \param RCP<Epetra_IntVector>& aggr_stat: stores internal information about the aggregation status of each node
       * \param int& nLocalAggregates: number of local aggregates that has been built on current processor so far.
       * \return 0)
       */
      int Phase4(const RCP<Epetra_CrsGraph>& amalA, ParameterList& params, RCP<Epetra_IntVector>& amal_Aggregates, RCP<Epetra_IntVector>& aggr_stat, int& nLocalAggregates);

    protected:
      int nGlobalDirichletBlocks;         ///< number of Dirichlet boundary blocks

    private:
      std::map<int,int> myrowid2globalamalblockid_;               ///< map: local row id of unamalgamated matrix on current proc -> global block id
      std::map<int,std::vector<int> > globalamalblockid2myrowid_;  ///< map: global block id -> myrowid of unamalgamated matrix (only for global block ids of current proc)
      std::vector<int>  globalamalblockids_;                      ///< vector with global block ids for amalgamated matrix on current proc

      RCP<Epetra_Map> map_;               ///< row map of input matrix
      RCP<Epetra_Map> amal_map_;          ///< row map of amalgamated matrix
      int nVerbose_;                      ///< verbosity level
  };
}

#endif /* AGGREGATION_METHOD_UNCOUPLED_H_ */
