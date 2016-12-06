/*!
\file xdofmapcreation_parallel_utils.cpp

\brief defines unknowns based on the intersection pattern from the xfem intersection

this is related to the physics of the fluid problem and therefore should not be part of the standard xfem routines

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
 *------------------------------------------------------------------------------------------------*/


#include "xdofmapcreation_parallel_utils.H"

#include "../drt_combust/combust_interface.H"
#include "dofkey.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_io/io_pstream.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::fillNodalDofKeySet(
    const COMBUST::InterfaceHandleCombust& ih,
    const std::map<int, std::set<XFEM::FieldEnr> >&  nodalDofSet,
    std::set<XFEM::DofKey >&      nodaldofkeyset
)
{
  nodaldofkeyset.clear();
  // loop all (non-overlapping = Row)-Nodes and store the DOF information w.t.h. of DofKeys
  for (int i=0; i<ih.FluidDis()->NumMyColNodes(); ++i)
  {
    const DRT::Node* actnode = ih.FluidDis()->lColNode(i);
    const int gid = actnode->Id();
    std::map<int, std::set<XFEM::FieldEnr> >::const_iterator entry = nodalDofSet.find(gid);
    if (entry == nodalDofSet.end())
    {
      // no dofs for this node... must be a hole or somethin'
      continue;
    }
    const std::set<XFEM::FieldEnr> & dofset = entry->second;

    std::set<XFEM::FieldEnr>::const_iterator fieldenr;
    for(fieldenr = dofset.begin(); fieldenr != dofset.end(); ++fieldenr )
    {
      nodaldofkeyset.insert(XFEM::DofKey(gid, *fieldenr));
    }
  };
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::updateNodalDofMap(
    const COMBUST::InterfaceHandleCombust& ih,
    std::map<int, std::set<XFEM::FieldEnr> >&  nodalDofSet,
    const std::set<XFEM::DofKey >&      nodaldofkeyset
)
{
  // pack data on all processors
  for(std::set<XFEM::DofKey >::const_iterator dofkey=nodaldofkeyset.begin();
      dofkey != nodaldofkeyset.end();
      ++dofkey)
  {
    const int nodegid = dofkey->getGid();
    const XFEM::FieldEnr fieldenr = dofkey->getFieldEnr();
    // is node a row node?
//    if (ih.xfemdis()->HaveGlobalNode(nodegid))
      nodalDofSet[nodegid].insert(fieldenr);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::packDofKeys(
    const std::set<XFEM::DofKey >&     dofkeyset,
    DRT::PackBuffer&                            dataSend )
{
  // pack data on all processors
  for(std::set<XFEM::DofKey >::const_iterator dofkey=dofkeyset.begin(); dofkey != dofkeyset.end(); ++dofkey)
  {
    DRT::ParObject::AddtoPack(dataSend,*dofkey);
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::unpackDofKeys(
    const std::vector<char>&                     dataRecv,
    std::set<XFEM::DofKey >&       dofkeyset )
{
  std::vector<char>::size_type index = 0;
  while (index < dataRecv.size())
  {
    std::vector<char> data;
    DRT::ParObject::ExtractfromPack(index, dataRecv, data);
    dofkeyset.insert(XFEM::DofKey(data));
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::syncNodalDofs(
    const COMBUST::InterfaceHandleCombust& ih,
    std::map<int, std::set<XFEM::FieldEnr> >&  nodalDofSet)
{
  const Epetra_Comm& comm = ih.FluidDis()->Comm();
  const int myrank = comm.MyPID();
  const int numproc = comm.NumProc();

  int size_one = 1;

  DRT::Exporter exporter(comm);

  int dest = myrank+1;
  if(myrank == (numproc-1))
    dest = 0;

  int source = myrank-1;
  if(myrank == 0)
    source = numproc-1;

  std::set<XFEM::DofKey >  original_dofkeyset;
  XFEM::fillNodalDofKeySet(ih, nodalDofSet, original_dofkeyset);

  std::set<XFEM::DofKey >  new_dofkeyset;
  //  new_dofkeyset = original_dofkeyset;
  for (std::set<XFEM::DofKey >::const_iterator dofkey = original_dofkeyset.begin(); dofkey != original_dofkeyset.end(); ++dofkey)
  {
    new_dofkeyset.insert(*dofkey);
  }

  // pack date for initial send
  DRT::PackBuffer data;
  XFEM::packDofKeys(original_dofkeyset, data);
  data.StartPacking();
  XFEM::packDofKeys(original_dofkeyset, data);

  std::vector<char> dataSend;
  swap( dataSend, data() );

#ifdef DEBUG
  IO::cout << "proc " << myrank << ": sending "<< original_dofkeyset.size() << " dofkeys to proc " << dest << IO::endl;
#endif

  // send data in a circle
  for(int num = 0; num < numproc-1; num++)
  {
    std::vector<int> lengthSend(1,0);
    lengthSend[0] = dataSend.size();

#ifdef DEBUG
    IO::cout << "proc " << myrank << ": sending "<< lengthSend[0] << " bytes to proc " << dest << IO::endl;
#endif

    // send length of the data to be received ...
    MPI_Request req_length_data;
    int length_tag = 0;
    exporter.ISend(myrank, dest, &(lengthSend[0]) , size_one, length_tag, req_length_data);

    // ... and receive length
    std::vector<int> lengthRecv(1,0);
    exporter.Receive(source, length_tag, lengthRecv, size_one);
    exporter.Wait(req_length_data);

    //    if(lengthRecv[0] > 0) ??
    //    {
    // send actual data ...
    int data_tag = 4;
    MPI_Request req_data;
    exporter.ISend(myrank, dest, &(dataSend[0]), lengthSend[0], data_tag, req_data);
    // ... and receive date
    std::vector<char> dataRecv(lengthRecv[0]);
    exporter.ReceiveAny(source, data_tag, dataRecv, lengthRecv[0]);
    exporter.Wait(req_data);


    // unpack dofkeys from char array
    std::set<XFEM::DofKey >       dofkeyset;
    XFEM::unpackDofKeys(dataRecv, dofkeyset);
#ifdef DEBUG
    IO::cout << "proc " << myrank << ": receiving "<< lengthRecv[0] << " bytes from proc " << source << IO::endl;
    IO::cout << "proc " << myrank << ": receiving "<< dofkeyset.size() << " dofkeys from proc " << source << IO::endl;
#endif
    // get all dofkeys whose nodegid is on this proc in the coloumnmap
    for (std::set<XFEM::DofKey >::const_iterator dofkey = dofkeyset.begin(); dofkey != dofkeyset.end(); ++dofkey)
    {
      const int nodegid = dofkey->getGid();
      if (ih.FluidDis()->HaveGlobalNode(nodegid))
      {
        new_dofkeyset.insert(*dofkey);
      }
    }
    // make received data the new 'to be sent' data
    dataSend = dataRecv;
    //    }
    //    else
    //    {
    //      dataSend.clear();
    //    }
    comm.Barrier();
  }   // loop over procs

#ifdef DEBUG
  IO::cout << "sync nodal dofs on proc " << myrank << ": before/after -> " << original_dofkeyset.size()<< "/" << new_dofkeyset.size() << IO::endl;
#endif

  XFEM::updateNodalDofMap(ih, nodalDofSet, new_dofkeyset);

}
