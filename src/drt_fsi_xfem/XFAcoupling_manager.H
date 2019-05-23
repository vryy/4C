/*!----------------------------------------------------------------------
\file XFAcoupling_manager.H
\brief Coupling Manager for eXtended Fluid Ale Coupling

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249

*----------------------------------------------------------------------*/

#ifndef XFACOUPLING_MANAGER_H
#define XFACOUPLING_MANAGER_H

#include "coupling_manager.H"
#include "coupling_comm_manager.H"

namespace FLD
{
  class XFluid;
}

namespace ADAPTER
{
  class AleFpsiWrapper;
  class Structure;
}  // namespace ADAPTER

namespace XFEM
{
  class XFACoupling_Manager : public Coupling_Manager, public Coupling_Comm_Manager
  {
   public:
    /// constructor
    // in idx ... idx[0] structureal discretization index , idx[1] fluid discretization index in the
    // blockmatrix
    explicit XFACoupling_Manager(Teuchos::RCP<FLD::XFluid> xfluid,
        Teuchos::RCP<ADAPTER::AleFpsiWrapper> ale, std::vector<int> idx,
        Teuchos::RCP<ADAPTER::Structure> structure = Teuchos::null);
    //! @name Destruction
    //@{

    /// virtual to get polymorph destruction
    virtual ~XFACoupling_Manager() {}

    //! predict states in the coupling object
    virtual void PredictCouplingStates();

    //! Set required displacement & velocity states in the coupling object
    virtual void SetCouplingStates();

    //! Initializes the couplings (done at the beginning of the algorithm after fields have their
    //! state for timestep n) -- not yet done here
    virtual void InitCouplingStates() { return; }

    //! Add the coupling matrixes to the global systemmatrix
    // in ... scaling between xfluid evaluated coupling matrixes and coupled systemmatrix
    virtual void AddCouplingMatrix(LINALG::BlockSparseMatrixBase& systemmatrix, double scaling);

    //! Add the coupling rhs

    // in scaling ... scaling between xfluid evaluated coupling rhs and coupled rhs
    // in me ... global map extractor of coupled problem (same index used as for idx)
    virtual void AddCouplingRHS(
        Teuchos::RCP<Epetra_Vector> rhs, const LINALG::MultiMapExtractor& me, double scaling);

    //! Update (Perform after Each Timestep) -- nothing to do here
    virtual void Update(double scaling) { return; }

    //! Write Output -- nothing to do here
    virtual void Output(IO::DiscretizationWriter& writer) { return; }

    //! Read Restart -- nothing to do here
    virtual void ReadRestart(IO::DiscretizationReader& reader) { return; }

   private:
    //! Ale Object
    Teuchos::RCP<ADAPTER::AleFpsiWrapper> ale_;
    //! eXtendedFluid
    Teuchos::RCP<FLD::XFluid> xfluid_;

    // Global Index in the blockmatrix of the coupled sytem [0] = fluid-, [1] = ale- block,
    // [2]struct- block
    std::vector<int> idx_;

    //! Structural Object (just set if Ale is coupled to a structure)
    Teuchos::RCP<ADAPTER::Structure> structure_;
  };
}  // namespace XFEM
#endif
