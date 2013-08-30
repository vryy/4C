/*----------------------------------------------------------------------*/
/*!
\file adapter_coupling_mortar.H

\brief A class providing coupling capabilities based on mortar methods

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FSI_COUPLING_MORTAR_H
#define FSI_COUPLING_MORTAR_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

// forward declarations

namespace LINALG
{
  class SparseMatrix;
}

namespace DRT
{
  class Discretization;
}

namespace MORTAR
{
  class MortarInterface;
}

namespace ADAPTER
{
  /// Couple non-matching interface meshes using mortar method
  /*!
    This is a generic class used to couple any non-matching meshes
    (or more general: discretizations) at interfaces. The current
    applications in BACI emcompass FSI coupling algorithms (i.e. to
    interpolate between fluid and structure fields at the interface)
    and fluid mesh tying algorithms (i.e. to couple non-matching
    Eulerian fluid meshes). All the hard work is actually done by
    the MORTAR::MortarInterface class (thus we use the mortar method).

    The major part of this code is the Setup() method that gets the
    non-matching interface meshes on input, initializes the mortar
    interface and computes the so-called coupling matrices D and M.

    The actual coupling methods MasterToSlave() and SlaveToMaster()
    just evaluate one simple equation each, i.e. primal variables
    are projected from master to slave side via \f$D^{-1} M\f$ when
    calling MasterToSlaver(), and dual variables are projected from
    slave to master side via \f$M^T D^{-T}\f$ when calling SlaveToMaster().

    Whenever you want to add a new problem class, check whether you
    can re-use one of the already existing Setup() methods. If not,
    feel free to write your own tailored Setup() method.
   */
  class CouplingMortar
  {
  public:

    /// empty constructor
    CouplingMortar();

    //! setup the machinery (generalized version)
    /// master and slave discretization are identical for sliding ale and fluid/scatra meshtying
    /// ale discretization is Techos::null for  sliding ale and fluid/scatra meshtying
    void Setup( Teuchos::RCP<DRT::Discretization>   masterdis,      ///> master discretization
                Teuchos::RCP<DRT::Discretization>   slavedis,       ///> slave discretization
                Teuchos::RCP<DRT::Discretization>   aledis,         ///> ale discretization
                std::vector<int>                    coupleddof,     ///> vector defining coupled degree's of freedom
                const std::string&                       couplingcond,   ///> string for coupling condition
                const Epetra_Comm&                  comm,           ///> communicate
                bool                                structslave=false,  ///> flag defining if structure is slave side
                bool                                slidingale=false);  //> flag defining sliding ale case

    /// perform mesh initialization
    void MeshInit( Teuchos::RCP<DRT::Discretization> slavedis,
                   Teuchos::RCP<DRT::Discretization> aledis,
                   Teuchos::RCP<Epetra_Map> masterdofrowmap,
                   Teuchos::RCP<Epetra_Map> slavedofrowmap,
                   const Epetra_Comm& comm,
                   bool   structslave);

    /// perform mesh initialization
    void ReEvaluate(const Epetra_Comm& comm);

    void Evaluate(Teuchos::RCP<Epetra_Vector> idisp);

    /// filling matrices (case of transferring same dofs on two different meshes)
    void Evaluate(Teuchos::RCP<Epetra_Vector> idispma, Teuchos::RCP<Epetra_Vector> idispsl);

    /// transfer a dof vector from master to slave
    Teuchos::RCP<Epetra_Vector> MasterToSlave
    (
        Teuchos::RCP<Epetra_Vector> mv ///< e.g. for FSI: structural interface displacements
    ) const;

    /// const version
    Teuchos::RCP<Epetra_Vector> MasterToSlave
    (
        Teuchos::RCP<const Epetra_Vector> mv ///< e.g. for FSI: structural interface displacements
    ) const;

    /// transfer a dof vector from slave to master
    Teuchos::RCP<Epetra_Vector> SlaveToMaster
    (
        Teuchos::RCP<Epetra_Vector> sv   ///< e.g. for FSI: fluid interface loads
    ) const;

    /// const version
    Teuchos::RCP<Epetra_Vector> SlaveToMaster
    (
        Teuchos::RCP<const Epetra_Vector> sv   ///< e.g. for FSI: fluid interface loads
    ) const;

    /// the mortar interface itself
    Teuchos::RCP<MORTAR::MortarInterface> Interface() const { return interface_;}

    /// the interface dof map of the master side
    Teuchos::RCP<Epetra_Map> MasterDofRowMap() const { return masterdofrowmap_; }

    /// the interface dof map of the slave side
    Teuchos::RCP<Epetra_Map> SlaveDofRowMap()  const { return slavedofrowmap_; }

    /// return mortar D matrix
    Teuchos::RCP<LINALG::SparseMatrix> GetDMatrix() const { return D_;}

    /// return inverted mortar matrix \f$D^{-1}\f$
    Teuchos::RCP<LINALG::SparseMatrix> GetDinvMatrix() const { return Dinv_;}

    /// return mortar M matrix
    Teuchos::RCP<LINALG::SparseMatrix> GetMMatrix() const { return M_;}

    /// return the mortar projection operator \f$P=D^{-1} M\f$
    Teuchos::RCP<LINALG::SparseMatrix> GetMortarTrafo() const { return DinvM_;}


  private:

    ///filling matrices by using mortar interface
    void Evaluate();

    Teuchos::RCP<MORTAR::MortarInterface> interface_;       ///< interface
    Teuchos::RCP<Epetra_Map>              masterdofrowmap_; ///< map of master row dofs
    Teuchos::RCP<Epetra_Map>              slavedofrowmap_;  ///< map of slave row dofs


    Teuchos::RCP<LINALG::SparseMatrix>    D_;               ///< mortar matrix
    Teuchos::RCP<LINALG::SparseMatrix>    Dinv_;            ///< inverse of mortar matrix D_
    Teuchos::RCP<LINALG::SparseMatrix>    M_;               ///< mortar matrix
    Teuchos::RCP<LINALG::SparseMatrix>    DinvM_;           ///< product of Dinv_ \times M_
  };
}

#endif // FSI_COUPLING_MORTAR_CONTACT_H
