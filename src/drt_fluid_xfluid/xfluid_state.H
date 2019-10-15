/*----------------------------------------------------------------------*/
/*! \file

\brief State class for (in)stationary XFEM fluid problems

\level 0

\maintainer  Martin Kronbichler
             kronbichler@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
 */
/*----------------------------------------------------------------------*/

#ifndef XFLUID_STATE_H_
#define XFLUID_STATE_H_

#include <Epetra_Vector.h>
#include <Epetra_Map.h>

#include <Teuchos_RCP.hpp>

#include "../drt_lib/drt_dserror.H"

#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_cut.H"

// forward declarations
namespace DRT
{
  class Discretization;
  class DiscretizationXFEM;
}  // namespace DRT

namespace GEO
{
  class CutWizard;
}

namespace LINALG
{
  class SparseMatrix;
  class MultiMapExtractor;
  class MapExtractor;
}  // namespace LINALG

namespace IO
{
  class DiscretizationWriter;
}

namespace XFEM
{
  class ConditionManager;
  class XFEMDofSet;
}  // namespace XFEM

namespace FLD
{
  namespace UTILS
  {
    class KSPMapExtractor;
  }


  /**
   * Container class for the state vectors and maps of the intersected background
   * fluid - tied to a specific intersection state (interface position,
   * independent from the fact, in which form the interface is given (boundary mesh or level-set
   * field)).
   */
  class XFluidState
  {
   public:
    /*!
     * \brief Container class for coupling state matrices and vectors for a certain coupling object
     */
    class CouplingState
    {
     public:
      //! ctor
      CouplingState()
          : C_sx_(Teuchos::null),
            C_xs_(Teuchos::null),
            C_ss_(Teuchos::null),
            rhC_s_(Teuchos::null),
            rhC_s_col_(Teuchos::null),
            is_active_(false)
      {
      }

      //! ctor  Initialize coupling matrices with fluid sysmat and fluid rhs vector
      CouplingState(const Teuchos::RCP<LINALG::SparseMatrix>& C_sx,
          const Teuchos::RCP<LINALG::SparseMatrix>& C_xs,
          const Teuchos::RCP<LINALG::SparseMatrix>& C_ss, const Teuchos::RCP<Epetra_Vector>& rhC_s,
          const Teuchos::RCP<Epetra_Vector>& rhC_s_col)
          : C_sx_(C_sx),  // this rcp constructor using const references copies also the ownership
                          // and increases the specific weak/strong reference counter
            C_xs_(C_xs),
            C_ss_(C_ss),
            rhC_s_(rhC_s),
            rhC_s_col_(rhC_s_col),
            is_active_(true)
      {
      }

      //! ctor  Initialize coupling matrices
      CouplingState(const Teuchos::RCP<const Epetra_Map>& xfluiddofrowmap,
          const Teuchos::RCP<DRT::Discretization>& slavediscret_mat,
          const Teuchos::RCP<DRT::Discretization>& slavediscret_rhs);

      //! zero coupling matrices and rhs vectors
      void ZeroCouplingMatricesAndRhs();

      //! complete coupling matrices and rhs vectors
      void CompleteCouplingMatricesAndRhs(
          const Epetra_Map& xfluiddofrowmap, const Epetra_Map& slavedofrowmap);

      //! destroy the coupling objects and it's content
      void Destroy(bool throw_exception = true);

      //! @name coupling matrices x: xfluid, s: coupling slave (structure, ALE-fluid, xfluid-element
      //! with other active dofset, etc.)
      //@{
      Teuchos::RCP<LINALG::SparseMatrix> C_sx_;  ///< slave - xfluid coupling block
      Teuchos::RCP<LINALG::SparseMatrix> C_xs_;  ///< xfluid - slave coupling block
      Teuchos::RCP<LINALG::SparseMatrix> C_ss_;  ///< slave - slave coupling block
      Teuchos::RCP<Epetra_Vector> rhC_s_;        ///< slave rhs block
      Teuchos::RCP<Epetra_Vector> rhC_s_col_;    ///< slave rhs block
      //@}

      bool is_active_;
    };

    /*!
     ctor for one-sided problems
     @param xfluiddofrowmap dof-rowmap of intersected fluid
     */
    explicit XFluidState(const Teuchos::RCP<XFEM::ConditionManager>& condition_manager,
        const Teuchos::RCP<GEO::CutWizard>& wizard, const Teuchos::RCP<XFEM::XFEMDofSet>& dofset,
        const Teuchos::RCP<const Epetra_Map>& xfluiddofrowmap,
        const Teuchos::RCP<const Epetra_Map>& xfluiddofcolmap);

    /// dtor
    virtual ~XFluidState() {}

    /// setup map extractors for dirichlet maps & velocity/pressure maps
    void SetupMapExtractors(
        const Teuchos::RCP<DRT::Discretization>& xfluiddiscret, const double& time);

    /// zero system matrix and related rhs vectors
    virtual void ZeroSystemMatrixAndRhs();

    /// zero all coupling matrices and rhs vectors for all coupling objects
    virtual void ZeroCouplingMatricesAndRhs();

    /// Complete coupling matrices and rhs vectors
    virtual void CompleteCouplingMatricesAndRhs();

    /// destroy the stored objects
    virtual bool Destroy();

    /// update the coordinates of the cut boundary cells
    void UpdateBoundaryCellCoords();


    //! @name Accessors
    //@{

    /// access to the cut wizard
    Teuchos::RCP<GEO::CutWizard> Wizard() const
    {
      if (wizard_ == Teuchos::null) dserror("Cut wizard is uninitialized!");
      return wizard_;
    }


    //! access to the xfem-dofset
    Teuchos::RCP<XFEM::XFEMDofSet> DofSet() { return dofset_; }

    virtual Teuchos::RCP<LINALG::MapExtractor> DBCMapExtractor() { return dbcmaps_; }

    virtual Teuchos::RCP<LINALG::MapExtractor> VelPresSplitter() { return velpressplitter_; }

    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return sysmat_; }
    virtual Teuchos::RCP<Epetra_Vector>& Residual() { return residual_; }
    virtual Teuchos::RCP<Epetra_Vector>& Zeros() { return zeros_; }
    virtual Teuchos::RCP<Epetra_Vector>& IncVel() { return incvel_; }
    virtual Teuchos::RCP<Epetra_Vector>& Velnp() { return velnp_; }
    virtual Teuchos::RCP<Epetra_Vector>& Veln() { return veln_; }
    virtual Teuchos::RCP<Epetra_Vector>& Accnp() { return accnp_; }
    virtual Teuchos::RCP<Epetra_Vector>& Accn() { return accn_; }
    //@}

    /// dof-rowmap of intersected fluid
    Teuchos::RCP<const Epetra_Map> xfluiddofrowmap_;

    /// dof-colmap of intersected fluid
    Teuchos::RCP<const Epetra_Map> xfluiddofcolmap_;

    /// system matrix (internally EpetraFECrs)
    Teuchos::RCP<LINALG::SparseMatrix> sysmat_;

    /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
    Teuchos::RCP<Epetra_Vector> zeros_;

    /// the vector containing body and surface forces
    Teuchos::RCP<Epetra_Vector> neumann_loads_;

    /// (standard) residual vector (rhs for the incremental form),
    Teuchos::RCP<Epetra_Vector> residual_;

    /// (standard) residual vector (rhs for the incremental form) wrt col map,
    Teuchos::RCP<Epetra_Vector> residual_col_;

    /// true (rescaled) residual vector without zeros at dirichlet positions
    Teuchos::RCP<Epetra_Vector> trueresidual_;

    /// nonlinear iteration increment vector
    Teuchos::RCP<Epetra_Vector> incvel_;

    //! @name acceleration/(scalar time derivative) at time n+1, n and n+alpha_M/(n+alpha_M/n)
    //@{
    Teuchos::RCP<Epetra_Vector> accnp_;
    Teuchos::RCP<Epetra_Vector> accn_;
    Teuchos::RCP<Epetra_Vector> accam_;
    //@}

    //! @name velocity and pressure at time n+1, n, n-1 and n+alpha_F
    //@{
    Teuchos::RCP<Epetra_Vector> velnp_;
    Teuchos::RCP<Epetra_Vector> veln_;
    Teuchos::RCP<Epetra_Vector> velnm_;
    Teuchos::RCP<Epetra_Vector> velaf_;
    //@}

    //! @name scalar at time n+alpha_F/n+1 and n+alpha_M/n
    //@{
    Teuchos::RCP<Epetra_Vector> scaaf_;
    Teuchos::RCP<Epetra_Vector> scaam_;
    //@}

    //! @name displacemets at time n+1, n and n-1 (if we have an XFEM-ALE-fluid)
    //@{
    Teuchos::RCP<Epetra_Vector> dispnp_;
    //  Teuchos::RCP<Epetra_Vector>         dispn_;
    //  Teuchos::RCP<Epetra_Vector>         dispnm_;
    //@}

    /// grid velocity (if we have an XFEM-ALE-fluid) (set from the adapter!)
    Teuchos::RCP<Epetra_Vector> gridvnp_;

    /// histvector --- a linear combination of velnm, veln (BDF)
    ///                or veln, accn (One-Step-Theta)
    Teuchos::RCP<Epetra_Vector> hist_;

    //! @name map extractors
    //@{
    /// maps for extracting Dirichlet and free DOF sets
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

    /// velocity/pressure map extractor, used for convergence check
    Teuchos::RCP<LINALG::MapExtractor> velpressplitter_;

    //@}


    //! @name coupling matrices for each coupling object (=key); x: xfluid, s: coupling slave
    //! (structure, ALE-fluid, xfluid-element with other active dofset, etc.)
    //@{
    std::map<int, Teuchos::RCP<CouplingState>> coup_state_;
    //@}

   protected:
    /// initialize all state members based on the xfluid dof-rowmap
    void InitStateVectors();

    /// initialize the system matrix of the intersected fluid
    void InitSystemMatrix();

    /// initialize coupling matrices and rhs vectors for all coupling objects
    void InitCouplingMatricesAndRhs();

    /// Complete coupling matrices and rhs vectors
    void CompleteCouplingMatricesAndRhs(const Teuchos::RCP<const Epetra_Map>& fluiddofrowmap);


   public:
    /*!
     \brief initialize ALE state vectors
     @param dispnp and grivnp vectors w.r.t initial full dofrowmap
     */
    void InitALEStateVectors(const Teuchos::RCP<DRT::DiscretizationXFEM>& xdiscret,
        Teuchos::RCP<const Epetra_Vector> dispnp_initmap,
        Teuchos::RCP<const Epetra_Vector> gridvnp_initmap);

    /// XFEM dofset
    Teuchos::RCP<XFEM::XFEMDofSet> dofset_;

    /// cut wizard
    Teuchos::RCP<GEO::CutWizard> wizard_;

    /// condition manager
    Teuchos::RCP<XFEM::ConditionManager> condition_manager_;
  };

}  // namespace FLD

#endif
