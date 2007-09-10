
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
  class CouplingMortar
  {
  public:

      /// empty constructor
      CouplingMortar();

      void Setup( const DRT::Discretization& masterdis,
                  const DRT::Discretization& slavedis,
                  Epetra_Comm& comm );

      /// transfer a dof vector from master to slave
      RefCountPtr<Epetra_Vector> MasterToSlave( RefCountPtr<Epetra_Vector> mv );

      /// transfer a dof vector from slave to master
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
#endif // TRILINOS_PACKAGE
#endif // CCADISCRET
