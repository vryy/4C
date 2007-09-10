/*!----------------------------------------------------------------------
\file drt_utils.cpp
\brief A collection of helper methods for namespace DRT

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#ifdef PARALLEL
#ifdef LINUX_MUENCH
typedef int idxtype;
extern "C"
{
  void METIS_PartGraphRecursive(int *, idxtype *, idxtype *,
				idxtype *, idxtype *, int *, int *, int *,
				int *, int *, idxtype *);
  void METIS_PartGraphKway(int *, idxtype *, idxtype *, idxtype *,
			   idxtype *, int *, int *, int *, int *, int *,
			   idxtype *);
}
#endif
#include <mpi.h>
#endif // PARALLEL

#include <algorithm>
#include <numeric>
#include <vector>

#include "drt_utils.H"
#include "drt_node.H"
#include "drt_dofset.H"
#include "drt_discret.H"
#include "../drt_s8/shell8.H"
#include "../drt_f2/fluid2.H"
#include "../drt_f2/condif2.H"
#include "../drt_f3/fluid3.H"
#include "../drt_ale2/ale2.H"
#include "../drt_ale3/ale3.H"
#include "../drt_w1/wall1.H"
#include "../drt_so3/so_hex8.H"
#include "../drt_so3/so_sh8.H"
#include "../drt_so3/so_tet10.H"
#include "../drt_so3/so_weg6.H"
#include "../drt_so3/so_disp.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/micromaterial.H"
#include "../drt_mat/neohooke.H"
#include "../drt_mat/convecdiffus.H"
#include "drt_dserror.H"


/*----------------------------------------------------------------------*
 |  allocate an instance of a specific impl. of ParObject (public) mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ParObject* DRT::Utils::Factory(const vector<char>& data)
{
  // mv ptr behind the size record
  const int* ptr = (const int*)(&data[0]);
  // get the type
  const int type = *ptr;
  switch(type)
  {
    case ParObject_Container:
    {
      DRT::Container* object = new DRT::Container();
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Condition:
    {
      DRT::Condition* object = new DRT::Condition();
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Node:
    {
      double dummycoord[3] = {999.,999.,999.};
      DRT::Node* object = new DRT::Node(-1,dummycoord,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Element:
    {
      dserror("DRT::Element is pure virtual, cannot create instance");
    }
    break;
#ifdef D_SHELL8
    case ParObject_Shell8:
    {
      DRT::Elements::Shell8* object = new DRT::Elements::Shell8(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Shell8Register:
    {
      DRT::Elements::Shell8Register* object =
                      new DRT::Elements::Shell8Register(DRT::Element::element_shell8);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_WALL1
    case ParObject_Wall1:
    {
      DRT::Elements::Wall1* object = new DRT::Elements::Wall1(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Wall1Register:
    {
      DRT::Elements::Wall1Register* object =
                      new DRT::Elements::Wall1Register(DRT::Element::element_wall1);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_FLUID2
    case ParObject_Fluid2:
    {
      DRT::Elements::Fluid2* object = new DRT::Elements::Fluid2(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Fluid2Register:
    {
      DRT::Elements::Fluid2Register* object =
                      new DRT::Elements::Fluid2Register(DRT::Element::element_fluid2);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Condif2:
    {
      DRT::Elements::Condif2* object = new DRT::Elements::Condif2(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Condif2Register:
    {
      DRT::Elements::Condif2Register* object =
                      new DRT::Elements::Condif2Register(DRT::Element::element_condif2);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_FLUID3
    case ParObject_Fluid3:
    {
      DRT::Elements::Fluid3* object = new DRT::Elements::Fluid3(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Fluid3Register:
    {
      DRT::Elements::Fluid3Register* object =
                      new DRT::Elements::Fluid3Register(DRT::Element::element_fluid3);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_ALE
    case ParObject_Ale3:
    {
      DRT::Elements::Ale3* object = new DRT::Elements::Ale3(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Ale3Register:
    {
      DRT::Elements::Ale3Register* object =
                      new DRT::Elements::Ale3Register(DRT::Element::element_ale3);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_ALE
    case ParObject_Ale2:
    {
      DRT::Elements::Ale2* object = new DRT::Elements::Ale2(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Ale2Register:
    {
      DRT::Elements::Ale2Register* object =
                      new DRT::Elements::Ale2Register(DRT::Element::element_ale2);
      object->Unpack(data);
      return object;
    }
    break;
#endif
#ifdef D_SOH8
    case ParObject_So_hex8:
    {
      DRT::Elements::So_hex8* object = new DRT::Elements::So_hex8(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Soh8Register:
    {
      DRT::Elements::Soh8Register* object =
                new DRT::Elements::Soh8Register(DRT::Element::element_so_hex8);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_So_sh8:
    {
      DRT::Elements::So_sh8* object = new DRT::Elements::So_sh8(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Sosh8Register:
    {
      DRT::Elements::Sosh8Register* object =
                new DRT::Elements::Sosh8Register(DRT::Element::element_sosh8);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_So_weg6:
    {
      DRT::Elements::So_weg6* object =
                new DRT::Elements::So_weg6(-1,-1);
      object->Unpack(data);
      return object;
    }
    case ParObject_Sow6Register:
    {
      DRT::Elements::Sow6Register* object =
                new DRT::Elements::Sow6Register(DRT::Element::element_so_weg6);
      object->Unpack(data);
      return object;
    }
    case ParObject_SoDisp:
    {
      DRT::Elements::SoDisp* object =
                new DRT::Elements::SoDisp(-1,-1);
      object->Unpack(data);
      return object;
    }
    case ParObject_SoDispRegister:
    {
      DRT::Elements::SoDispRegister* object =
                new DRT::Elements::SoDispRegister(DRT::Element::element_sodisp);
      object->Unpack(data);
      return object;
    }
#endif
#ifdef D_SOTET10
    case ParObject_So_tet10:
    {
      DRT::Elements::So_tet10* object = new DRT::Elements::So_tet10(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Sotet10Register:
    {
      DRT::Elements::Sotet10Register* object =
                new DRT::Elements::Sotet10Register(DRT::Element::element_so_tet10);
      object->Unpack(data);
      return object;
    }
    break;
#endif
    case ParObject_ElementRegister:
    {
      dserror("DRT::ElementRegister is pure virtual, cannot create instance");
    }
    break;
    case ParObject_NewtonianFluid:
    {
      MAT::NewtonianFluid* fluid = new MAT::NewtonianFluid();
      fluid->Unpack(data);
      return fluid;
    }
    case ParObject_StVenantKirchhoff:
    {
      MAT::StVenantKirchhoff* stvenantk = new MAT::StVenantKirchhoff();
      stvenantk->Unpack(data);
      return stvenantk;
    }
    case ParObject_MicroMaterial:
    {
      MAT::MicroMaterial* micro = new MAT::MicroMaterial();
      micro->Unpack(data);
      return micro;
    }
    case ParObject_NeoHooke:
    {
	MAT::NeoHooke* neo = new MAT::NeoHooke();
	neo->Unpack(data);
	return neo;
    }
    case ParObject_ConvecDiffus:
    {
      MAT::ConvecDiffus* condif = new MAT::ConvecDiffus();
      condif->Unpack(data);
      return condif;
    }
    default:
      dserror("Unknown type of ParObject instance: %d",type);
    break;
  }

  return NULL;
}

/*----------------------------------------------------------------------*
 |  allocate an element of a specific type (public)          mwgee 03|07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::Element> DRT::Utils::Factory(const string eletype,
                                              const int id,
                                              const int owner)
{
  enum TypeofElement
  {
    none,
    shell8,
    wall1,
    fluid2,
    condif2,
    fluid3,
    ale2,
    ale3,
    so_hex8,
    so_sh8,
    so_tet10,
    so_weg6,
    sodisp,
  };

  TypeofElement type = none;
  if (eletype=="none"); // dummy
  else if (eletype=="SHELL8") type = shell8;
  else if (eletype=="WALL")  type = wall1;
  else if (eletype=="FLUID2") type = fluid2;
  else if (eletype=="CONDIF2") type = condif2;
  else if (eletype=="FLUID3") type = fluid3;
  else if (eletype=="ALE2") type = ale2;
  else if (eletype=="ALE3") type = ale3;
  else if (eletype=="SOLIDH8") type = so_hex8;
  else if (eletype=="SOLIDSH8") type = so_sh8;
  else if (eletype=="SOLIDTET10") type = so_tet10;
  else if (eletype=="SOLIDW6") type = so_weg6;
  else if (eletype=="SOLID3") type = sodisp;
  // continue to add elements here....
  else dserror("Unknown type of finite element");


  switch (type)
  {
#ifdef D_SHELL8
    case shell8:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::Elements::Shell8(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_FLUID2
    case fluid2:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::Elements::Fluid2(id,owner));
      return ele;
    }
    break;
    case condif2:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::Elements::Condif2(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_FLUID3
    case fluid3:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::Elements::Fluid3(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_ALE
    case ale3:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::Elements::Ale3(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_ALE
    case ale2:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::Elements::Ale2(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_WALL1
    case wall1:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::Elements::Wall1(id,owner));
      return ele;
    }
    break;
#endif
#ifdef D_SOH8
    case so_hex8:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::Elements::So_hex8(id,owner));
      return ele;
    }
    break;
    case so_sh8:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::Elements::So_sh8(id,owner));
      return ele;
    }
    break;
    case so_weg6:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::Elements::So_weg6(id,owner));
      return ele;
    }
    break;
    case sodisp:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::Elements::SoDisp(id,owner));
      return ele;
    }
    break;
#ifdef D_SOTET10
    case so_tet10:
    {
      RefCountPtr<DRT::Element> ele = rcp(new DRT::Elements::So_tet10(id,owner));
      return ele;
    }
    break;
#endif
#endif
    // continue to add types of elements here....
    default:
      dserror("Unknown type '%s' of finite element", eletype.c_str());
    break;
  }

  return null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int DRT::Utils::FindMyPos(int nummyelements, const Epetra_Comm& comm)
{
  const int myrank  = comm.MyPID();
  const int numproc = comm.NumProc();

  vector<int> snum(numproc);
  vector<int> rnum(numproc);
  fill(snum.begin(), snum.end(), 0);
  snum[myrank] = nummyelements;

  comm.SumAll(&snum[0],&rnum[0],numproc);

  return std::accumulate(&rnum[0], &rnum[myrank], 0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Utils::AllreduceEMap(vector<int>& rredundant, const Epetra_Map& emap)
{
  int mynodepos = FindMyPos(emap.NumMyElements(), emap.Comm());

  vector<int> sredundant(emap.NumGlobalElements());
  fill(sredundant.begin(), sredundant.end(), 0);

  int* gids = emap.MyGlobalElements();
  copy(gids, gids+emap.NumMyElements(), &sredundant[mynodepos]);

  rredundant.resize(emap.NumGlobalElements());
  emap.Comm().SumAll(&sredundant[0], &rredundant[0], emap.NumGlobalElements());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::Utils::AllreduceEMap(map<int,int>& idxmap, const Epetra_Map& emap)
{
  idxmap.clear();

  vector<int> rredundant(emap.NumGlobalElements());
  AllreduceEMap(rredundant, emap);

  for (unsigned i=0; i<rredundant.size(); ++i)
  {
    idxmap[rredundant[i]] = i;
  }
}


/*----------------------------------------------------------------------*
 |  partition a graph using metis  (public)                  mwgee 11/06|
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_CrsGraph> DRT::Utils::PartGraphUsingMetis(
                                             const Epetra_CrsGraph& graph,
                                             const Epetra_Vector& weights)
{
  const int myrank   = graph.Comm().MyPID();
  const int numproc  = graph.Comm().NumProc();

  if (numproc==1)
  {
    RefCountPtr<Epetra_CrsGraph> outgraph = rcp(new Epetra_CrsGraph(graph));
    return outgraph;
  }

  // proc that will do the serial partitioning
  // the graph is collapsed to this proc
  // Normally this would be proc 0 but 0 always has so much to do.... ;-)
  int workrank = 0;
  if (graph.Comm().NumProc()>1) workrank=1;

  // get rowmap of the graph
  const Epetra_BlockMap& tmp = graph.RowMap();
  Epetra_Map rowmap(tmp.NumGlobalElements(),tmp.NumMyElements(),
                    tmp.MyGlobalElements(),0,graph.Comm());

  // build a target map that stores everything on proc workrank
  // We have arbirtary gids here and we do not tell metis about
  // them. So we have to keep rowrecv until the redistributed map is
  // build.
  vector<int> rowrecv(rowmap.NumGlobalElements());
  AllreduceEMap(rowrecv, rowmap);
  Epetra_Map tmap(rowmap.NumGlobalElements(),
                  (myrank == workrank) ? (int)rowrecv.size() : 0,
                  &rowrecv[0],
                  0,
                  rowmap.Comm());

  // export the graph to tmap
  Epetra_CrsGraph tgraph(Copy,tmap,108,false);
  Epetra_Export exporter(rowmap,tmap);
  int err = tgraph.Export(graph,exporter,Add);
  if (err<0) dserror("Graph export returned err=%d",err);
  tgraph.FillComplete();
  tgraph.OptimizeStorage();

  // export the weights to tmap
  Epetra_Vector tweights(tmap,false);
  err = tweights.Export(weights,exporter,Insert);
  if (err<0) dserror("Vector export returned err=%d",err);

  // do partitioning using metis on workrank
  vector<int> part(tmap.NumMyElements());
  if (myrank==workrank)
  {
    // metis requests indexes. So we need a reverse lookup from gids
    // to indexes.
    map<int,int> idxmap;
    for (unsigned i=0; i<rowrecv.size(); ++i)
    {
      idxmap[rowrecv[i]] = i;
    }

    vector<int> xadj(tmap.NumMyElements()+1);
    vector<int> adjncy(tgraph.NumGlobalNonzeros()); // the size is an upper bound
    vector<int> vwgt(tweights.MyLength());
    for (int i=0; i<tweights.MyLength(); ++i) vwgt[i] = (int)tweights[i];

    int count=0;
    xadj[0] = 0;
    for (int row=0; row<tgraph.NumMyRows(); ++row)
    {
      //cout << "xadj[" << row << "] = " << xadj[row] << endl;
      int grid = tgraph.RowMap().GID(row);
      int numindices;
      int* lindices;
      int err = tgraph.ExtractMyRowView(row,numindices,lindices);
      if (err) dserror("Epetra_CrsGraph::ExtractMyRowView returned err=%d",err);
      //cout << "adjncy: ";
      for (int col=0; col<numindices; ++col)
      {
        int gcid = tgraph.ColMap().GID(lindices[col]);
        if (gcid==grid) continue;
        adjncy[count] = idxmap[gcid];
        //cout << adjncy[count] << " ";
        ++count;
      }
      //cout << endl;
      xadj[row+1] = count;
    }
    //cout << "xadj[" << xadj.size()-1 << "] = " << xadj[xadj.size()-1] << endl;
    //cout << "tgraph.NumGlobalNonzeros() " << tgraph.NumGlobalNonzeros() << endl
    //     << "tmap.NumMyElements()       " << tmap.NumMyElements() << endl
    //     << "count                      " << count << endl;

    idxmap.clear();

    if (numproc<8) // better for smaller no. of partitions
    {
#ifdef PARALLEL
      int wgtflag=2;
      int numflag=0;
      int npart=numproc;
      int options[5] = { 0,3,1,1,0 };
      int edgecut=0;
      int nummyele = tmap.NumMyElements();
      METIS_PartGraphRecursive(&nummyele,
                               &xadj[0],
                               &adjncy[0],
                               &vwgt[0],
                               NULL,
                               &wgtflag,
                               &numflag,
                               &npart,
                               options,
                               &edgecut,
                               &part[0]);
#endif
    }
    else
    {
#ifdef PARALLEL
      int wgtflag=2;
      int numflag=0;
      int npart=numproc;
      int options[5] = { 0,3,1,1,0 };
      int edgecut=0;
      int nummyele = tmap.NumMyElements();
      METIS_PartGraphKway(&nummyele,
                          &xadj[0],
                          &adjncy[0],
                          &vwgt[0],
                          NULL,
                          &wgtflag,
                          &numflag,
                          &npart,
                          options,
                          &edgecut,
                          &part[0]);
#endif
    }
  } // if (myrank==workrank)

  // broadcast partitioning result
  int size = tmap.NumMyElements();
  tmap.Comm().Broadcast(&size,1,workrank);
  part.resize(size);
  tmap.Comm().Broadcast(&part[0],size,workrank);

  // loop part and count no. of nodes belonging to me
  // (we reuse part to save on memory)
  int count=0;
  for (int i=0; i<size; ++i)
    if (part[i]==myrank)
    {
      part[count] = rowrecv[i];
      ++count;
    }

  rowrecv.clear();

  // create map with new layout
  Epetra_Map newmap(size,count,&part[0],0,graph.Comm());

  // create the output graph and export to it
  RefCountPtr<Epetra_CrsGraph> outgraph =
                           rcp(new Epetra_CrsGraph(Copy,newmap,108,false));
  Epetra_Export exporter2(graph.RowMap(),newmap);
  err = outgraph->Export(graph,exporter2,Add);
  if (err<0) dserror("Graph export returned err=%d",err);
  outgraph->FillComplete();
  outgraph->OptimizeStorage();

  return outgraph;
}



/*----------------------------------------------------------------------*
 |  locallly extract a subset of values  (public)            mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Utils::ExtractMyValues(const Epetra_Vector& global,
                                 vector<double>& local,
                                 const vector<int> lm)
{
  const int ldim = (int)lm.size();
  local.resize(ldim);
  for (int i=0; i<ldim; ++i)
  {
    const int lid = global.Map().LID(lm[i]);
    if (lid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",global.Comm().MyPID(),lm[i]);
    local[i] = global[lid];
  }
  return;
}

/*
 * Send and receive lists of ints.  (heiner 09/07)
 */
void DRT::Utils::AllToAllCommunication( const Epetra_Comm& comm,
                                        const vector< vector<int> >& send,
                                        vector< vector<int> >& recv )
{
#ifndef PARALLEL

  dsassert(send.size()==1, "there has to be just one entry for sending");

  // make a copy
  recv.clear();
  recv.push_back(send[0]);

#else
    const Epetra_MpiComm& mpicomm = dynamic_cast<const Epetra_MpiComm&>(comm);

    vector<int> sendbuf;
    vector<int> sendcounts;
    sendcounts.reserve( comm.NumProc() );
    vector<int> sdispls;
    sdispls.reserve( comm.NumProc() );

    int displacement = 0;
    sdispls.push_back( 0 );
    for ( vector< vector<int> >::const_iterator iter = send.begin();
          iter != send.end(); ++iter )
    {
        sendbuf.insert( sendbuf.end(), iter->begin(), iter->end() );
        sendcounts.push_back( iter->size() );
        displacement += iter->size();
        sdispls.push_back( displacement );
    }

    vector<int> recvcounts( comm.NumProc() );

    // initial communication: Request. Send and receive the number of
    // ints we communicate with each process.

    int status = MPI_Alltoall( &sendcounts[0], 1, MPI_INT,
                               &recvcounts[0], 1, MPI_INT, mpicomm.GetMpiComm() );

    if ( status != MPI_SUCCESS )
        dserror( "MPI_Alltoall returned status=%d", status );

    vector<int> rdispls;
    rdispls.reserve( comm.NumProc() );

    displacement = 0;
    rdispls.push_back( 0 );
    for ( vector<int>::const_iterator iter = recvcounts.begin();
          iter != recvcounts.end(); ++iter )
    {
        displacement += *iter;
        rdispls.push_back( displacement );
    }

    vector<int> recvbuf( rdispls.back() );

    // transmit communication: Send and get the data.

    status = MPI_Alltoallv ( &sendbuf[0], &sendcounts[0], &sdispls[0], MPI_INT,
                             &recvbuf[0], &recvcounts[0], &rdispls[0], MPI_INT,
                             mpicomm.GetMpiComm() );
    if ( status != MPI_SUCCESS )
        dserror( "MPI_Alltoallv returned status=%d", status );

    recv.clear();
    for ( int proc = 0; proc < comm.NumProc(); ++proc )
    {
        recv.push_back( vector<int>( &recvbuf[rdispls[proc]], &recvbuf[rdispls[proc+1]] ) );
    }
#endif // PARALLEL
}

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
