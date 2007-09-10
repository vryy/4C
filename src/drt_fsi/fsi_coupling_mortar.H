
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#ifndef FSI_COUPLING_MORTAR_H
#define FSI_COUPLING_MORTAR_H

#include <string>
#include <vector>
#include <set>

#include <Teuchos_RefCountPtr.hpp>
#include <Epetra_Export.h>
#include <Epetra_Map.h>

#include "../drt_lib/drt_discret.H"


namespace FSI
{

  /// Management of coupling between two (matching) sets of nodes from two discretizations
  /*!
    Interface coupled problems with matching meshes need to transfer
    results between their interfaces. That is values which belong to
    the dofs of one side have to be accessed by the other side as
    well. In a parallel setting, of course, these sides do not in
    general reside on the same processors.

    This class handles the transfer of dof values at the interface
    (Epetra_Vectors build on the interface dof map of either side) to
    the other side. To distinguish both sides lets speak of master and
    slave, even though no side really dominates the other. On the
    contrary we provide perfect symmetry once the setup is done.

    The idea is simple: We have a master dof map that describes the
    distribution of the master's interface dofs. And we have a slave
    dof map that describes the distribution of the slave's interface
    dofs. Both maps, however, live in there respective
    communicator. We cannot transfer data between them. Furthermore
    the dof gids are totally independent of each other.

    So we build a mapping between the nodes from both meshes during
    the setup phase. This establishes the connection. From that
    connection we deduce (in a somewhat painful way) a permutated
    master dof map and a permutated slave dof map. These permutated
    maps are bound to have the same layout as the normal maps from the
    other side. So we can exchange dof values between fields by simply
    copying from a normal Epetra_Vector to the permutated
    Epetra_Vector from the other side without actually looking at the
    respective maps. Afterwards the communication happens within one
    field in the usual fashion. So the transfer functions
    MasterToSlave() and SlaveToMaster() are quite simple. The hard
    work happens (once) during setup.

    \author u.kue
    \date 06/07
   */
  class CouplingMortar
  {
  public:

    /// empty constructor
    CouplingMortar();

    /// setup coupling of nodes marked with condition
    /*!
      Setup the whole thing. Find matching nodes via octtree and build
      appropiate dof maps.

      \param masterdis  (i) master side mesh
      \param slavedis   (i) slave side mesh
      \param coupname   (i) name of condition on coupling nodes
     */
    void Setup( const DRT::Discretization& masterdis, const DRT::Discretization& slavedis );

    /// transfer a dof vector from master to slave
    Teuchos::RefCountPtr<Epetra_Vector> MasterToSlave(Teuchos::RefCountPtr<Epetra_Vector> mv);

    /// transfer a dof vector from slave to master
    Teuchos::RefCountPtr<Epetra_Vector> SlaveToMaster(Teuchos::RefCountPtr<Epetra_Vector> sv);

    /// the interface dof map of the master side
    RefCountPtr<Epetra_Map> MasterDofMap() const { return masterdofmap_; }

      /// the interface dof map of the slave side
      RefCountPtr<Epetra_Map> SlaveDofMap()  const { return slavedofmap_; }

  private:
      /// find all local nodes in the interface
      void FindInterfaceNodes( const DRT::Discretization& dis, std::set<int>& nodes );

  };

}

#endif // FSI_COUPLING_MORTAR_H
#endif // TRILINOS_PACKAGE
#endif // CCADISCRET
