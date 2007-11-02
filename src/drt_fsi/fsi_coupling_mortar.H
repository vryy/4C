
#ifdef CCADISCRET

#ifndef FSI_COUPLING_MORTAR_H
#define FSI_COUPLING_MORTAR_H

#include <Teuchos_RefCountPtr.hpp>
#include <Epetra_Map.h>

#include "../drt_lib/drt_discret.H"


namespace FSI
{
  /// Couple nonmatching FSI interfaces using mortar method (trilinos::Moertel)
  /*!
    This is a special class used by the Dirichlet-Neumann coupling algorithm
    to interpolate between fluid and structure fields at the interface. The
    hard work is done by the Moertel package from trilinos. (Thus we use the
    mortar method here.)

    The major part of this code is the Setup() method that feeds the FSI
    interface into Moertel and extracts the coupling matrices D and M in
    return. A lot of the work involved is due to the different node, element
    and dof numbering in baci and moertel. Whereas baci uses independant
    numberings for both fields (fluid and structure) moertel thinks in terms
    of one big system of equations, so all participating nodes, elements, dofs
    have to have an unique number. The returned matrices D and M also come
    within the global mortar map (the system of three blocks of equations that
    couples the two sides and contains the Lagrangian multiplier.) So we have
    to extract D and M to get interface matrices.

    The actual coupling methods MasterToSlave() and SlaveToMaster() just
    evaluate one simple equation each.

    \note This class is specific to FSI coupling with structure being the
    master and fluid the slave.

    \node The force coupling in SlaveToMaster() expacts interface loads on the
    fluid side and returns consistant interface forces on the structural side.

    \author h.kue
    \date 09/07
   */
  class CouplingMortar
  {
  public:

    /// empty constructor
    CouplingMortar();

    /// setup the machinery
    void Setup( const DRT::Discretization& masterdis,
                const DRT::Discretization& slavedis,
                Epetra_Comm& comm );

    /// transfer a dof vector from master to slave
    /*!
      Expects the structural interface displacements and returns interface
      displacements on the fluid side.
     */
    RefCountPtr<Epetra_Vector> MasterToSlave( RefCountPtr<Epetra_Vector> mv );

    /// transfer a dof vector from slave to master
    /*!
      Expects the fluid interface loads and returns consistant interface
      forces at the structural side. So here we have a subtle change of the
      meaning these numbers.
     */
    RefCountPtr<Epetra_Vector> SlaveToMaster( RefCountPtr<Epetra_Vector> sv );

    /// the interface dof map of the master side
    RefCountPtr<Epetra_Map> MasterDofMap() const { return masterdofmap_; }

    /// the interface dof map of the slave side
    RefCountPtr<Epetra_Map> SlaveDofMap()  const { return slavedofmap_; }

  private:
    /// find all local nodes and elements in the interface
    void FindInterfaceObjects(
      const DRT::Discretization& dis,
      map<int, DRT::Node*>& nodes,
      map< int, RefCountPtr<DRT::Element> >& elements );

    RefCountPtr<Epetra_Map> masterdofmap_;
    RefCountPtr<Epetra_Map> slavedofmap_;
    RefCountPtr<Epetra_Map> slavemortardofmap_;

    RefCountPtr<Epetra_CrsMatrix> D_;
    RefCountPtr<Epetra_CrsMatrix> Dinv_;
    RefCountPtr<Epetra_CrsMatrix> M_;
  };
}

#endif // FSI_COUPLING_MORTAR_H
#endif // CCADISCRET
