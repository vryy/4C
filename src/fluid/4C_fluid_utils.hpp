/*-----------------------------------------------------------*/
/*! \file

\brief utility functions for fluid problems


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_UTILS_HPP
#define FOUR_C_FLUID_UTILS_HPP

#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_linalg_blocksparsematrix.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::DOFSets
{
  class DofSet;
  class DofSetInterface;
}  // namespace Core::DOFSets

namespace Core::LinAlg
{
  class MultiMapExtractor;
  class MapExtractor;
}  // namespace Core::LinAlg

namespace FLD
{
  namespace UTILS
  {
    /// velocity-pressure block matrix split strategy
    class VelPressSplitStrategy
    {
     public:
      /// construct with a block matrix base
      explicit VelPressSplitStrategy(Core::LinAlg::BlockSparseMatrixBase& mat)
          : mat_(mat),
            matrix00_(mat_.matrix(0, 0)),
            matrix01_(mat_.matrix(0, 1)),
            matrix10_(mat_.matrix(1, 0)),
            matrix11_(mat_.matrix(1, 1)),
            numdim_(-1),
            numdofpernode_(-1)
      {
      }

      /// find row block to a given row gid
      int row_block(int lrow, int rgid)
      {
        if ((lrow % numdofpernode_) < numdim_) return 0;
        return 1;
      }

      /// find column block to a given column gid
      int col_block(int rblock, int lcol, int cgid)
      {
        if ((lcol % numdofpernode_) < numdim_) return 0;
        return 1;
      }

      /// assemble into the given block
      void assemble(int eid, int myrank, const std::vector<int>& lmstride,
          const Core::LinAlg::SerialDenseMatrix& Aele, const std::vector<int>& lmrow,
          const std::vector<int>& lmrowowner, const std::vector<int>& lmcol)
      {
        const int lrowdim = (int)lmrow.size();
        const int lcoldim = (int)lmcol.size();

        if (mat_.filled())
        {
          // We use the maps of the matrix to gain fast access to the LID's.
          // Assembling with SumIntoMyValues based on LID's is two times faster
          // than inserting single values based on the global row and column ids.

          // There is the case of nodes without dofs (XFEM).
          // If no row dofs are present on this proc, their is nothing to assemble.
          // However, the subsequent check for coldofs (in DEBUG mode) would incorrectly fail.
          bool doit = false;
          for (int lrow = 0; lrow < lrowdim; ++lrow)
            if (lmrowowner[lrow] == myrank)
            {
              doit = true;  // This proc owns at least one row of interest
              break;
            }
          if (!doit) return;

          // get the maps
          const Epetra_Map& colmap00 = mat_.matrix(0, 0).col_map();
          const Epetra_Map& colmap01 = mat_.matrix(0, 1).col_map();
          const Epetra_Map& colmap10 = mat_.matrix(1, 0).col_map();
          const Epetra_Map& colmap11 = mat_.matrix(1, 1).col_map();
          const Epetra_Map& rowmap00 = mat_.matrix(0, 0).row_map();
          const Epetra_Map& rowmap01 = mat_.matrix(0, 1).row_map();
          const Epetra_Map& rowmap10 = mat_.matrix(1, 0).row_map();
          const Epetra_Map& rowmap11 = mat_.matrix(1, 1).row_map();

          // prepare vectors for holding column local ids and the values to be assembled
          const int nnode = lcoldim / numdofpernode_;
          std::vector<double> values0(numdim_ * nnode);
          std::vector<double> values1(nnode);
          std::vector<int> localcol00(numdim_ * nnode);
          std::vector<int> localcol01(nnode);
          std::vector<int> localcol10(numdim_ * nnode);
          std::vector<int> localcol11(nnode);

          // fill vectors with the LID's
          int nodespassed = 0;
          for (int lcol = 0; lcol < lcoldim; ++lcol)
          {
            const int cgid = lmcol[lcol];
            const int rest = (lcol % numdofpernode_);
            if (rest < numdim_)
            {
              const int pos = nodespassed * numdim_ + rest;
              localcol00[pos] = (colmap00.LID(cgid));
              localcol10[pos] = (colmap10.LID(cgid));
            }
            else
            {
              const int pos = nodespassed;
              localcol01[pos] = (colmap01.LID(cgid));
              localcol11[pos] = (colmap11.LID(cgid));
              nodespassed++;
            }
          }

          // loop rows of local matrix
          for (int lrow = 0; lrow < lrowdim; ++lrow)
          {
            // check ownership of row
            if (lmrowowner[lrow] != myrank) continue;

            const int rgid = lmrow[lrow];
            int rlid0;
            int rlid1;
            int rowblock = row_block(lrow, rgid);
            if (rowblock == 0)
            {
              rlid0 = rowmap00.LID(rgid);
              rlid1 = rowmap01.LID(rgid);
            }
            else
            {
              rlid0 = rowmap10.LID(rgid);
              rlid1 = rowmap11.LID(rgid);
            }
#ifdef FOUR_C_ENABLE_ASSERTIONS
            if (rlid0 < 0) FOUR_C_THROW("Sparse matrix A does not have global row %d", rgid);
            if (rlid1 < 0) FOUR_C_THROW("Sparse matrix A does not have global row %d", rgid);
#endif
            int errone = 0;
            // separate the values of the current row
            nodespassed = 0;
            for (int lcol = 0; lcol < lcoldim; ++lcol)
            {
              double val = Aele(lrow, lcol);
              const int rest = lcol % numdofpernode_;
              if (rest < numdim_)
              {
                int pos = nodespassed * numdim_ + rest;
                values0[pos] = val;
              }
              else
              {
                values1[nodespassed] = val;
                nodespassed++;
              }
            }

            // now assemble
            if (rowblock == 0)
            {  // rowblock 0
              errone = matrix00_.epetra_matrix()->SumIntoMyValues(
                  rlid0, nnode * numdim_, values0.data(), localcol00.data());
              if (errone)
                FOUR_C_THROW("Epetra_CrsMatrix::SumIntoMyValues returned error code %d", errone);
              errone = matrix01_.epetra_matrix()->SumIntoMyValues(
                  rlid1, nnode, values1.data(), localcol01.data());
              if (errone)
                FOUR_C_THROW("Epetra_CrsMatrix::SumIntoMyValues returned error code %d", errone);
            }
            else
            {  // rowblock 1
              errone = matrix10_.epetra_matrix()->SumIntoMyValues(
                  rlid0, nnode * numdim_, values0.data(), localcol10.data());
              if (errone)
                FOUR_C_THROW("Epetra_CrsMatrix::SumIntoMyValues returned error code %d", errone);
              errone = matrix11_.epetra_matrix()->SumIntoMyValues(
                  rlid1, nnode, values1.data(), localcol11.data());
              if (errone)
                FOUR_C_THROW("Epetra_CrsMatrix::SumIntoMyValues returned error code %d", errone);
            }
          }  // for (int lrow=0; lrow<ldim; ++lrow)
        }
        else
        {
          // the maps of the blockmatrix are not yet available; use global ids instead

          // loop rows of local matrix
          for (int lrow = 0; lrow < lrowdim; ++lrow)
          {
            // check ownership of row
            if (lmrowowner[lrow] != myrank) continue;

            int rgid = lmrow[lrow];
            int rblock = row_block(lrow, rgid);

            for (int lcol = 0; lcol < lcoldim; ++lcol)
            {
              double val = Aele(lrow, lcol);
              int cgid = lmcol[lcol];
              int cblock = col_block(rblock, lcol, cgid);

              Core::LinAlg::SparseMatrix& matrix = mat_.matrix(rblock, cblock);
              matrix.assemble(val, rgid, cgid);
            }
          }
        }
      }

      /// assemble into the given block
      void assemble(double val, int rgid, int cgid)
      {
        int rblock = row_block(0, rgid);
        int cblock = col_block(rblock, 0, cgid);
        Core::LinAlg::SparseMatrix& matrix = mat_.matrix(rblock, cblock);
        matrix.assemble(val, rgid, cgid);
      }

      /// assemble the remaining ghost entries
      void complete() {}

      /// set number of velocity dofs
      void set_numdim(int numdim)
      {
        numdim_ = numdim;
        numdofpernode_ = numdim + 1;
      }

     private:
      /// my block matrix base
      Core::LinAlg::BlockSparseMatrixBase& mat_;

      // the four sub-matrices of the whole matrix
      Core::LinAlg::SparseMatrix& matrix00_;
      Core::LinAlg::SparseMatrix& matrix01_;
      Core::LinAlg::SparseMatrix& matrix10_;
      Core::LinAlg::SparseMatrix& matrix11_;

      /// number of velocity dofs
      int numdim_;

      /// number of dofs per node (= numdim_ +1)
      int numdofpernode_;
    };


    /// (FSI) interface block matrix split strategy
    class InterfaceSplitStrategy : public Core::LinAlg::DefaultBlockMatrixStrategy
    {
     public:
      explicit InterfaceSplitStrategy(Core::LinAlg::BlockSparseMatrixBase& mat)
          : Core::LinAlg::DefaultBlockMatrixStrategy(mat)
      {
      }

      /// assemble into the given block
      void assemble(int eid, int myrank, const std::vector<int>& lmstride,
          const Core::LinAlg::SerialDenseMatrix& Aele, const std::vector<int>& lmrow,
          const std::vector<int>& lmrowowner, const std::vector<int>& lmcol)
      {
        if (condelements_->find(eid) != condelements_->end())
        {
          // if we have an element with conditioned nodes, we have to do the
          // default assembling
          Core::LinAlg::DefaultBlockMatrixStrategy::assemble(
              eid, myrank, lmstride, Aele, lmrow, lmrowowner, lmcol);
        }
        else
        {
          // if there are no conditioned nodes we can simply assemble to the
          // internal matrix
          Core::LinAlg::SparseMatrix& matrix = mat().matrix(0, 0);
          matrix.assemble(eid, lmstride, Aele, lmrow, lmrowowner, lmcol);
        }
      }

      void assemble(double val, int rgid, int cgid)
      {
        // forward single value assembling
        Core::LinAlg::DefaultBlockMatrixStrategy::assemble(val, rgid, cgid);
      }

      void set_cond_elements(Teuchos::RCP<std::set<int>> condelements)
      {
        condelements_ = condelements;
      }

     private:
      Teuchos::RCP<std::set<int>> condelements_;
    };


    /// Stress manager manages everything to do with stresses and wallshearstresses
    class StressManager
    {
     public:
      /// constructor
      StressManager(Teuchos::RCP<Core::FE::Discretization> discret,
          Teuchos::RCP<Epetra_Vector> dispnp, const bool alefluid, const int numdim);

      /// initialize smoothing of stresses
      void init_aggr(Teuchos::RCP<Core::LinAlg::SparseOperator> sysmat);

      /// update and return WSS vector
      Teuchos::RCP<Epetra_Vector> get_wall_shear_stresses(
          Teuchos::RCP<const Epetra_Vector> trueresidual, const double dt);

      /// return WSS vector (without updating the mean wss vector)
      Teuchos::RCP<Epetra_Vector> get_pre_calc_wall_shear_stresses(
          Teuchos::RCP<const Epetra_Vector> trueresidual);

      /// return WSS vector always without aggregation, even if scale separation matrix exists
      Teuchos::RCP<Epetra_Vector> get_wall_shear_stresses_wo_agg(
          Teuchos::RCP<const Epetra_Vector> trueresidual);

      /// update and return stress vector
      Teuchos::RCP<Epetra_Vector> get_stresses(
          Teuchos::RCP<const Epetra_Vector> trueresidual, const double dt);

      /// return stress vector (without updating the mean stress vector)
      Teuchos::RCP<Epetra_Vector> get_pre_calc_stresses(
          Teuchos::RCP<const Epetra_Vector> trueresidual);

      /// return stress vector always without aggregation, even if scale separation matrix exists
      Teuchos::RCP<Epetra_Vector> get_stresses_wo_agg(
          Teuchos::RCP<const Epetra_Vector> trueresidual);

      /// return flag if StressManager has already been initialized
      bool is_init() { return isinit_; };

     private:
      /// return stress vector
      Teuchos::RCP<Epetra_Vector> calc_stresses(Teuchos::RCP<const Epetra_Vector> trueresidual);

      /// integrate shape functions at nodes marked by condition
      Teuchos::RCP<Epetra_Vector> integrate_interface_shape(std::string condname);

      /// calculate WSS based on residual
      Teuchos::RCP<Epetra_Vector> calc_wall_shear_stresses(Teuchos::RCP<Epetra_Vector> stresses);

      /// smooth stress/wss via ML-aggregation
      Teuchos::RCP<Epetra_Vector> aggreagte_stresses(Teuchos::RCP<Epetra_Vector> wss);

      /// time average stresses
      Teuchos::RCP<Epetra_Vector> time_average_stresses(
          Teuchos::RCP<const Epetra_Vector> stresses, double dt);

      /// time average wss
      Teuchos::RCP<Epetra_Vector> time_average_wss(
          Teuchos::RCP<const Epetra_Vector> wss, double dt);

      /// Calculate Aggregation Matrix
      void calc_sep_enr(Teuchos::RCP<Core::LinAlg::SparseOperator> sysmat);

      /// fluid discretization
      const Teuchos::RCP<Core::FE::Discretization> discret_;

      /// displacement at time \f$t^{n+1}\f$
      const Teuchos::RCP<Epetra_Vector> dispnp_;

      /// do we move the fluid mesh and calculate the fluid on this moving mesh?
      const bool alefluid_;

      /// number of spatial dimensions
      const int numdim_;

      /// filtering matrix for wall shear stress
      Teuchos::RCP<Core::LinAlg::SparseMatrix> sep_enr_;

      /// wss calculation type
      const Inpar::FLUID::WSSType wss_type_;

      /// weighted sum of all prior stresses
      Teuchos::RCP<Epetra_Vector> sum_stresses_;

      /// weighted sum of all prior wss
      Teuchos::RCP<Epetra_Vector> sum_wss_;

      /// time the stresses are averaged for
      double sum_dt_stresses_;

      /// time the wss are averaged for
      double sum_dt_wss_;

      bool isinit_;
    };


    void SetupFluidFluidVelPresSplit(const Core::FE::Discretization& fluiddis, int ndim,
        const Core::FE::Discretization& alefluiddis, Core::LinAlg::MapExtractor& extractor,
        Teuchos::RCP<Epetra_Map> fullmap);

    /*!

    \brief calculate lift&drag forces and angular momenta

    Lift and drag forces are based upon the right hand side true-residual entities
    of the corresponding nodes. The contribution of the end node of a line is entirely
    added to a present L&D force.

    Idea of this routine:

    create

    map< label, std::set<Core::Nodes::Node*> >

    which is a set of nodes to each L&D Id
    nodal forces of all the nodes within one set are added to one L&D force

    Notice: Angular moments obtained from lift&drag forces currently refere to the
            initial configuration, i.e. are built with the coordinates X of a particular
            node irrespective of its current position.

          \author chfoe
      \date 11/07
     */
    void LiftDrag(const Teuchos::RCP<const Core::FE::Discretization>
                      dis,  //! fluid discretization (node distribution, conditions)
        const Teuchos::RCP<const Epetra_Vector> trueresidual,  //! vector of nodalforces
        const Teuchos::RCP<const Epetra_Vector>
            dispnp,      ///< solution vector with velocities (and pressure)
        const int ndim,  //! number of dimensions
        Teuchos::RCP<std::map<int, std::vector<double>>>&
            liftdragvals,  //! the computed values for lift and drag in an array
        bool alefluid      //! flag for ALE case
    );

    /*!
     * \brief proc 0 writes transient liftdrag values to files (1 file per label)
     *
     * \author axel
     * \date 02/09
     */
    void WriteLiftDragToFile(const double time,  ///< current real time
        const int step,                          ///< time step
        const std::map<int, std::vector<double>>&
            liftdragvals  ///< the computed values for lift and drag in an array
    );

    /*!
      \brief integrate mass flow over surfaces

      for each condition Id, compute the flow through the surface


      \return a map, where for each condition Id, we get the flowrate.
              positive and negative signs indicate net inflow and outflow


      \author Axel Gerstenberger
      \date 10/08
     */
    // std::map<int,double> ComputeSurfaceFlowRates(
    //    Core::FE::Discretization&       dis  ,      ///< the discretisation (node
    //    distribution, conditions) const Teuchos::RCP<Epetra_Vector>   velnp       ///< solution
    //    vector with velocities and pressure
    ///< (only velocities are used)
    //    );


    /*!
      \brief integrate mass flow over surfaces

      for each condition Id, compute the flow through a boundary condition

      \return a map, where for each condition Id, we get the flowrate.
              positive and negative signs indicate net inflow and outflow

      \author Ursula Mayer
      \date 10/08
     */
    std::map<int, double> ComputeFlowRates(
        Core::FE::Discretization& dis,  ///< the discretisation (node distribution, conditions)
        const Teuchos::RCP<Epetra_Vector>&
            velnp,                      ///< solution vector with velocities (and pressure)
        const std::string& condstring,  ///< name of the condition (LineFlowRate or SurfaceFlowRate)
        const Inpar::FLUID::PhysicalType physicaltype  ///< physical type of flow
    );

    std::map<int, double> ComputeFlowRates(
        Core::FE::Discretization& dis,  ///< the discretisation (node distribution, conditions)
        const Teuchos::RCP<Epetra_Vector>&
            velnp,  ///< solution vector with velocities (and pressure)
        const Teuchos::RCP<Epetra_Vector>& gridvel,  ///< solution vector with grid velocities (ALE)
        const Teuchos::RCP<Epetra_Vector>&
            dispnp,                     ///< solution vector with mesh displacements (ALE)
        const std::string& condstring,  ///< name of the condition (LineFlowRate or SurfaceFlowRate)
        const Inpar::FLUID::PhysicalType physicaltype  ///< physical type of flow
    );

    std::map<int, double> compute_volume(
        Core::FE::Discretization& dis,  ///< the discretisation (node distribution, conditions)
        const Teuchos::RCP<Epetra_Vector>&
            velnp,  ///< solution vector with velocities (and pressure)
        const Teuchos::RCP<Epetra_Vector>& gridvel,  ///< solution vector with grid velocities (ALE)
        const Teuchos::RCP<Epetra_Vector>&
            dispnp,  ///< solution vector with mesh displacements (ALE)
        const Inpar::FLUID::PhysicalType physicaltype  ///< physical type of flow
    );

    /*!
     * \brief proc 0 writes the flow rate values for each condition ID to a file
     *
     * \author mayer
     * \date 01/10
     */
    void WriteDoublesToFile(const double time, const int step, const std::map<int, double>& data,
        const std::string& name);


    void WriteVolumeToFile(
        const double time, const int step, const std::map<int, double>& flowrates);

    /*!
    \brief Project gradient and store vector in param list

    */
    void ProjectGradientAndSetParam(Teuchos::RCP<Core::FE::Discretization> discret,
        Teuchos::ParameterList& eleparams, Teuchos::RCP<const Epetra_Vector> vel,
        const std::string paraname, bool alefluid);

    /*!
    \brief Project velocity gradient, depends on time integrator used

    */
    Teuchos::RCP<Epetra_MultiVector> ProjectGradient(Teuchos::RCP<Core::FE::Discretization> discret,
        Teuchos::RCP<const Epetra_Vector> vel, bool alefluid);

    /*!
      \brief integrate impulse rate over surfaces

      for each condition Id, compute the impulse rate vector through the surface

      integral over surface: rho * u_i * u_j * n_j dx

      \return a map, where for each condition Id, we get the impulse rate.
              positive and negative signs indicate net inflow and outflow


      \author Axel Gerstenberger
      \date 10/08
     */
    std::map<int, Core::LinAlg::Matrix<3, 1>> ComputeSurfaceImpulsRates(
        Core::FE::Discretization& dis,  ///< the discretisation (node distribution, conditions)
        const Teuchos::RCP<Epetra_Vector> velnp  ///< solution vector with velocities and pressure
                                                 ///< (only velocities are used)
    );


  }  // namespace UTILS
}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
