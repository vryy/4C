/*-----------------------------------------------------------*/
/*!

\brief utility functions for fluid problems

\maintainer Martin Kronbichler

\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FLUID_UTILS_H
#define FLUID_UTILS_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "../drt_inpar/inpar_fluid.H"
#include <Teuchos_TimeMonitor.hpp>

// forward declarations
namespace DRT
{
  class DofSet;
  class DofSetInterface;
}  // namespace DRT

namespace LINALG
{
  class MultiMapExtractor;
  class MapExtractor;
}  // namespace LINALG

namespace FLD
{
  namespace UTILS
  {
    /// velocity-pressure block matrix split strategy
    class VelPressSplitStrategy
    {
     public:
      /// construct with a block matrix base
      explicit VelPressSplitStrategy(LINALG::BlockSparseMatrixBase& mat)
          : mat_(mat),
            matrix00_(mat_.Matrix(0, 0)),
            matrix01_(mat_.Matrix(0, 1)),
            matrix10_(mat_.Matrix(1, 0)),
            matrix11_(mat_.Matrix(1, 1)),
            numdim_(-1),
            numdofpernode_(-1)
      {
      }

      /// find row block to a given row gid
      int RowBlock(int lrow, int rgid)
      {
        if ((lrow % numdofpernode_) < numdim_) return 0;
        return 1;
      }

      /// find column block to a given column gid
      int ColBlock(int rblock, int lcol, int cgid)
      {
        if ((lcol % numdofpernode_) < numdim_) return 0;
        return 1;
      }

      /// assemble into the given block
      void Assemble(int eid, int myrank, const std::vector<int>& lmstride,
          const Epetra_SerialDenseMatrix& Aele, const std::vector<int>& lmrow,
          const std::vector<int>& lmrowowner, const std::vector<int>& lmcol)
      {
        const int lrowdim = (int)lmrow.size();
        const int lcoldim = (int)lmcol.size();

        if (mat_.Filled())
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
          const Epetra_Map& colmap00 = mat_.Matrix(0, 0).ColMap();
          const Epetra_Map& colmap01 = mat_.Matrix(0, 1).ColMap();
          const Epetra_Map& colmap10 = mat_.Matrix(1, 0).ColMap();
          const Epetra_Map& colmap11 = mat_.Matrix(1, 1).ColMap();
          const Epetra_Map& rowmap00 = mat_.Matrix(0, 0).RowMap();
          const Epetra_Map& rowmap01 = mat_.Matrix(0, 1).RowMap();
          const Epetra_Map& rowmap10 = mat_.Matrix(1, 0).RowMap();
          const Epetra_Map& rowmap11 = mat_.Matrix(1, 1).RowMap();

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
            int rowblock = RowBlock(lrow, rgid);
            ;
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
#ifdef DEBUG
            if (rlid0 < 0) dserror("Sparse matrix A does not have global row %d", rgid);
            if (rlid1 < 0) dserror("Sparse matrix A does not have global row %d", rgid);
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
              errone = matrix00_.EpetraMatrix()->SumIntoMyValues(
                  rlid0, nnode * numdim_, &values0[0], &localcol00[0]);
              if (errone)
                dserror("Epetra_CrsMatrix::SumIntoMyValues returned error code %d", errone);
              errone = matrix01_.EpetraMatrix()->SumIntoMyValues(
                  rlid1, nnode, &values1[0], &localcol01[0]);
              if (errone)
                dserror("Epetra_CrsMatrix::SumIntoMyValues returned error code %d", errone);
            }
            else
            {  // rowblock 1
              errone = matrix10_.EpetraMatrix()->SumIntoMyValues(
                  rlid0, nnode * numdim_, &values0[0], &localcol10[0]);
              if (errone)
                dserror("Epetra_CrsMatrix::SumIntoMyValues returned error code %d", errone);
              errone = matrix11_.EpetraMatrix()->SumIntoMyValues(
                  rlid1, nnode, &values1[0], &localcol11[0]);
              if (errone)
                dserror("Epetra_CrsMatrix::SumIntoMyValues returned error code %d", errone);
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
            int rblock = RowBlock(lrow, rgid);

            for (int lcol = 0; lcol < lcoldim; ++lcol)
            {
              double val = Aele(lrow, lcol);
              int cgid = lmcol[lcol];
              int cblock = ColBlock(rblock, lcol, cgid);

              LINALG::SparseMatrix& matrix = mat_.Matrix(rblock, cblock);
              matrix.Assemble(val, rgid, cgid);
            }
          }
        }
      }

      /// assemble into the given block
      void Assemble(double val, int rgid, int cgid)
      {
        int rblock = RowBlock(0, rgid);
        int cblock = ColBlock(rblock, 0, cgid);
        LINALG::SparseMatrix& matrix = mat_.Matrix(rblock, cblock);
        matrix.Assemble(val, rgid, cgid);
      }

      /// assemble the remaining ghost entries
      void Complete() {}

      /// set number of velocity dofs
      void SetNumdim(int numdim)
      {
        numdim_ = numdim;
        numdofpernode_ = numdim + 1;
      }

     private:
      /// my block matrix base
      LINALG::BlockSparseMatrixBase& mat_;

      // the four sub-matrices of the whole matrix
      LINALG::SparseMatrix& matrix00_;
      LINALG::SparseMatrix& matrix01_;
      LINALG::SparseMatrix& matrix10_;
      LINALG::SparseMatrix& matrix11_;

      /// number of velocity dofs
      int numdim_;

      /// number of dofs per node (= numdim_ +1)
      int numdofpernode_;
    };


    /// (FSI) interface block matrix split strategy
    class InterfaceSplitStrategy : public LINALG::DefaultBlockMatrixStrategy
    {
     public:
      explicit InterfaceSplitStrategy(LINALG::BlockSparseMatrixBase& mat)
          : LINALG::DefaultBlockMatrixStrategy(mat)
      {
      }

      /// assemble into the given block
      void Assemble(int eid, int myrank, const std::vector<int>& lmstride,
          const Epetra_SerialDenseMatrix& Aele, const std::vector<int>& lmrow,
          const std::vector<int>& lmrowowner, const std::vector<int>& lmcol)
      {
        if (condelements_->find(eid) != condelements_->end())
        {
          // if we have an element with conditioned nodes, we have to do the
          // default assembling
          LINALG::DefaultBlockMatrixStrategy::Assemble(
              eid, myrank, lmstride, Aele, lmrow, lmrowowner, lmcol);
        }
        else
        {
          // if there are no conditioned nodes we can simply assemble to the
          // internal matrix
          LINALG::SparseMatrix& matrix = Mat().Matrix(0, 0);
          matrix.Assemble(eid, lmstride, Aele, lmrow, lmrowowner, lmcol);
        }
      }

      void Assemble(double val, int rgid, int cgid)
      {
        // forward single value assembling
        LINALG::DefaultBlockMatrixStrategy::Assemble(val, rgid, cgid);
      }

      void SetCondElements(Teuchos::RCP<std::set<int>> condelements)
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
      StressManager(Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<Epetra_Vector> dispnp,
          const bool alefluid, const int numdim);

      /// initialize smoothing of stresses
      void InitAggr(Teuchos::RCP<LINALG::SparseOperator> sysmat);

      /// update and return WSS vector
      Teuchos::RCP<Epetra_Vector> GetWallShearStresses(
          Teuchos::RCP<const Epetra_Vector> trueresidual, const double dt);

      /// return WSS vector (without updating the mean wss vector)
      Teuchos::RCP<Epetra_Vector> GetPreCalcWallShearStresses(
          Teuchos::RCP<const Epetra_Vector> trueresidual);

      /// return WSS vector always without aggregation, even if scale separation matrix exists
      Teuchos::RCP<Epetra_Vector> GetWallShearStressesWOAgg(
          Teuchos::RCP<const Epetra_Vector> trueresidual);

      /// update and return stress vector
      Teuchos::RCP<Epetra_Vector> GetStresses(
          Teuchos::RCP<const Epetra_Vector> trueresidual, const double dt);

      /// return stress vector (without updating the mean stress vector)
      Teuchos::RCP<Epetra_Vector> GetPreCalcStresses(
          Teuchos::RCP<const Epetra_Vector> trueresidual);

      /// return stress vector always without aggregation, even if scale separation matrix exists
      Teuchos::RCP<Epetra_Vector> GetStressesWOAgg(Teuchos::RCP<const Epetra_Vector> trueresidual);

      /// return flag if StressManager has already been initialized
      bool IsInit() { return isinit_; };

     private:
      /// return stress vector
      Teuchos::RCP<Epetra_Vector> CalcStresses(Teuchos::RCP<const Epetra_Vector> trueresidual);

      /// integrate shape functions at nodes marked by condition
      Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape(std::string condname);

      /// calculate WSS based on residual
      Teuchos::RCP<Epetra_Vector> CalcWallShearStresses(Teuchos::RCP<Epetra_Vector> stresses);

      /// smooth stress/wss via ML-aggregation
      Teuchos::RCP<Epetra_Vector> AggreagteStresses(Teuchos::RCP<Epetra_Vector> wss);

      /// time average stresses
      Teuchos::RCP<Epetra_Vector> TimeAverageStresses(
          Teuchos::RCP<const Epetra_Vector> stresses, double dt);

      /// time average wss
      Teuchos::RCP<Epetra_Vector> TimeAverageWss(
          Teuchos::RCP<const Epetra_Vector> stresses, double dt);

      /// Calculate Aggregation Matrix
      void CalcSepEnr(Teuchos::RCP<LINALG::SparseOperator> sysmat);

      /// fluid discretization
      const Teuchos::RCP<DRT::Discretization> discret_;

      /// displacement at time \f$t^{n+1}\f$
      const Teuchos::RCP<Epetra_Vector> dispnp_;

      /// do we move the fluid mesh and calculate the fluid on this moving mesh?
      const bool alefluid_;

      /// number of spatial dimensions
      const int numdim_;

      /// filtering matrix for wall shear stress
      Teuchos::RCP<LINALG::SparseMatrix> SepEnr_;

      /// wss calculation type
      const INPAR::FLUID::WSSType WssType_;

      /// weighted sum of all prior stresses
      Teuchos::RCP<Epetra_Vector> SumStresses_;

      /// weighted sum of all prior wss
      Teuchos::RCP<Epetra_Vector> SumWss_;

      /// time the stresses are averaged for
      double SumDtStresses_;

      /// time the wss are averaged for
      double SumDtWss_;

      bool isinit_;
    };


    /// Split dof row map of discretization into velocity and pressure maps
    /*!

      We assume that each node possesses ndim velocity dofs and (optionally)
      one pressure dof. The dof row map is thus split in two.
      If one node has more than ndim dofs, only the ndimth dof is included in the pressure part

      The velocity dofs are assigned to map 0 (the other map) and the pressure
      dofs are assigned to map 1 (the condition map).

      \param dis : (in) Discretization
      \param ndim : (in) dimensions of problem (2 or 3)
      \param extractor : (out) ready made map splitter

      \author u.kue
      \date 02/08
     */
    void SetupFluidSplit(
        const DRT::DiscretizationInterface& dis, int ndim, LINALG::MultiMapExtractor& extractor);

    /// Split dof row map of discretization into velocity and pressure maps
    /*!

      We assume that each node possesses ndim velocity dofs and (optionally)
      one pressure dof. The dof row map is thus split in two.

      The velocity dofs are assigned to map 0 (the other map) and the pressure
      dofs are assigned to map 1 (the condition map).

      \param dis : (in) Discretization
      \param ndim : (in) dimensions of problem (2 or 3)
      \param extractor : (out) ready made map splitter

      \author u.kue
      \date 02/08
     */
    void SetupFluidSplit(const DRT::Discretization& dis, const DRT::DofSetInterface& dofset,
        int ndim, LINALG::MapExtractor& extractor);


    /// Split dof row map of discretization into velocity and pressure maps
    /*!
      This versions is adopted for XFluid and XFluidFluid where more dofsets are used for one node.
      To get a unique split for velocity and pressure we have to define velocity ndof and pressure
      ndof. The pres_ndof defines the offset between velocity dofs for two nodes. This splitter
      works also for cases when nodes do not carry dofs (nodes in void domain!).

      The velocity dofs are assigned to map 0 (the other map) and the pressure
      dofs are assigned to map 1 (the condition map).

      \param dis : (in) Discretization
      \param fluid_ndof : (in) ndofs for fluid vector (2 or 3)
      \param pres_ndof : (in) ndofs for pressure (usually 1)
      \param extractor : (out) ready made map splitter

      \author schott
      \date 12/11
     */
    void SetupFluidSplit(const DRT::Discretization& dis, int fluid_ndof, int pres_ndof,
        LINALG::MultiMapExtractor& extractor);



    void SetupFluidFluidVelPresSplit(const DRT::Discretization& fluiddis, int ndim,
        const DRT::Discretization& alefluiddis, LINALG::MapExtractor& extractor,
        Teuchos::RCP<Epetra_Map> fullmap);

    /*!

    \brief calculate lift&drag forces and angular momenta

    Lift and drag forces are based upon the right hand side true-residual entities
    of the corresponding nodes. The contribution of the end node of a line is entirely
    added to a present L&D force.

    Idea of this routine:

    create

    map< label, std::set<DRT::Node*> >

    which is a set of nodes to each L&D Id
    nodal forces of all the nodes within one set are added to one L&D force

    Notice: Angular moments obtained from lift&drag forces currently refere to the
            initial configuration, i.e. are built with the coordinates X of a particular
            node irrespective of its current position.

        Maintainer: Ursula Rasthofer
                rasthofer@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15236
          \author chfoe
      \date 11/07
     */
    void LiftDrag(const Teuchos::RCP<const DRT::Discretization>
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
    //    DRT::Discretization&       dis  ,      ///< the discretisation (node distribution,
    //    conditions) const Teuchos::RCP<Epetra_Vector>   velnp       ///< solution vector with
    //    velocities and pressure
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
        DRT::Discretization& dis,  ///< the discretisation (node distribution, conditions)
        const Teuchos::RCP<Epetra_Vector>&
            velnp,                      ///< solution vector with velocities (and pressure)
        const std::string& condstring,  ///< name of the condition (LineFlowRate or SurfaceFlowRate)
        const INPAR::FLUID::PhysicalType physicaltype  ///< physical type of flow
    );

    std::map<int, double> ComputeFlowRates(
        DRT::Discretization& dis,  ///< the discretisation (node distribution, conditions)
        const Teuchos::RCP<Epetra_Vector>&
            velnp,  ///< solution vector with velocities (and pressure)
        const Teuchos::RCP<Epetra_Vector>& gridvel,  ///< solution vector with grid velocities (ALE)
        const Teuchos::RCP<Epetra_Vector>&
            dispnp,                     ///< solution vector with mesh displacements (ALE)
        const std::string& condstring,  ///< name of the condition (LineFlowRate or SurfaceFlowRate)
        const INPAR::FLUID::PhysicalType physicaltype  ///< physical type of flow
    );

    std::map<int, double> ComputeVolume(
        DRT::Discretization& dis,  ///< the discretisation (node distribution, conditions)
        const Teuchos::RCP<Epetra_Vector>&
            velnp,  ///< solution vector with velocities (and pressure)
        const Teuchos::RCP<Epetra_Vector>& gridvel,  ///< solution vector with grid velocities (ALE)
        const Teuchos::RCP<Epetra_Vector>&
            dispnp,  ///< solution vector with mesh displacements (ALE)
        const INPAR::FLUID::PhysicalType physicaltype  ///< physical type of flow
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
    void ProjectGradientAndSetParam(Teuchos::RCP<DRT::Discretization> discret,
        Teuchos::ParameterList& eleparams, Teuchos::RCP<const Epetra_Vector> vel,
        const std::string paraname, bool alefluid);

    /*!
    \brief Project velocity gradient, depends on time integrator used

    */
    Teuchos::RCP<Epetra_MultiVector> ProjectGradient(Teuchos::RCP<DRT::Discretization> discret,
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
    std::map<int, LINALG::Matrix<3, 1>> ComputeSurfaceImpulsRates(
        DRT::Discretization& dis,  ///< the discretisation (node distribution, conditions)
        const Teuchos::RCP<Epetra_Vector> velnp  ///< solution vector with velocities and pressure
                                                 ///< (only velocities are used)
    );


  }  // namespace UTILS
}  // namespace FLD

#endif  // #ifndef FLUID_UTILS_H
