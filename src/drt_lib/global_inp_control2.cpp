/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <list>
#include <string>
#include <algorithm>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "Epetra_SerialDenseMatrix.h"
#include "global_inp_control2.H"
#include "Epetra_SerialComm.h"


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

#ifdef DEBUG
/*!----------------------------------------------------------------------
  \brief the tracing variable

  <pre>                                                         m.gee 8/00
  defined in pss_ds.c, declared in tracing.h
  </pre>
 *----------------------------------------------------------------------*/
extern struct _CCA_TRACE         trace;
#endif

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;

/*!----------------------------------------------------------------------
  \brief file pointers

  <pre>                                                         m.gee 8/00
  This structure struct _FILES allfiles is defined in input_control_global.c
  and the type is in standardtypes.h
  It holds all file pointers and some variables needed for the FRSYSTEM
  </pre>
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;

/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;

/*----------------------------------------------------------------------*/
/*!
  \brief Section start positions of excluded section in input file

  Another global variable!

  \author u.kue
  \date 03/07
 */
/*----------------------------------------------------------------------*/
extern map<string,ifstream::pos_type> ExcludedSectionPositions;


/*----------------------------------------------------------------------*/
/*!
  \brief helper class to read the elements of a discretization

  Together with NodeReader this class constitutes a (almost) parallel
  and efficient reading mechanism for discretizations from dat files.

  We face the following problem:

  - There are elements and nodes. One set of elements per
    discretization. One set of nodes with the nodes from all
    discretizations.

  - Elements and nodes have ids. These are unique but otherwise
    arbitrary.

  - We cannot afford to read all elements or all nodes on one
    processor.

  - Only processor 0 can actually read the (ascii) input file

  - We do not want to setup (that is read) elements more than once.

  The idea is to read blocks of elements and nodes and distribute them
  to different processors at first. Afterwards a reasonable
  distribution can be calculated and the discretizations can be
  redistributed. How this work is done is a mere
  technicality. However, we need to be able to use the discretization
  in a partially constructed state. In particular we need to read, add
  and distribute elements even if nodes not yet are available.

  \author u.kue
  \date 05/07
 */
/*----------------------------------------------------------------------*/
class ElementReader
{
  friend class NodeReader;

public:

  /// construct element reader for a given field that reads a given section
  /*!
    Create empty discretization and append it to given field.
   */
  ElementReader(string name, int pos, RefCountPtr<Epetra_Comm> comm, string sectionname);

  /// give the discretization this reader fills
  RefCountPtr<DRT::Discretization> Discretization() const { return dis_; }

private:

  /// read elements, partition node graph
  void Partition();

  /// allreduce node row map to simplify node reading
  void CollapseRowMap();

  /// finialize reading. FillComplete.
  void Complete();

  /// discretization name
  string name_;

  /// my comm
  RefCountPtr<Epetra_Comm> comm_;

  /// my section to read
  string sectionname_;

  /// my discretization
  RefCountPtr<DRT::Discretization> dis_;

  /// node row map
  RefCountPtr<Epetra_Map> rownodes_;

  /// node col map
  RefCountPtr<Epetra_Map> colnodes_;

  /// allreduced node map
  RefCountPtr<Epetra_Map> collapsedrows_;

  /// element row map
  RefCountPtr<Epetra_Map> roweles_;

  /// element col map
  RefCountPtr<Epetra_Map> coleles_;
};


/*----------------------------------------------------------------------*/
/*!
  \brief helper class to read all nodes shared by many discretizations

  The interface class to the reading subsystem.
  \author u.kue
  \date 05/07
 */
/*----------------------------------------------------------------------*/
class NodeReader
{
public:

  /// construct a reader that reads a given section
  NodeReader(RefCountPtr<Epetra_Comm> comm, string sectionname);

  /// add an element reader for each discretization
  void AddElementReader(RefCountPtr<ElementReader> er) { ereader_.push_back(er); }

  /// do the actual reading
  void Read();

private:

  /// find the discretization that contains the given node
  RefCountPtr<DRT::Discretization> FindDisNode(int nodeid);

  /// my comm
  RefCountPtr<Epetra_Comm> comm_;

  /// my section to read
  string sectionname_;

  /// my element readers
  vector<RefCountPtr<ElementReader> > ereader_;
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ElementReader::ElementReader(string name,
                             int pos,
                             RefCountPtr<Epetra_Comm> comm,
                             string sectionname)
  : name_(name), comm_(comm), sectionname_(sectionname),
    dis_(rcp(new DRT::Discretization(name,comm)))
{
  vector<RefCountPtr<DRT::Discretization> >* discretization;

  // We can test for NULL because the field variable was allocated by CCACALLOC.
  if (field[pos].ccadis==NULL)
  {
    discretization = new vector<RefCountPtr<DRT::Discretization> >();
    field[pos].ccadis = (void*)discretization;
  }
  else
  {
    discretization = (vector<RefCountPtr<DRT::Discretization> >*)field[pos].ccadis;
  }
  discretization->push_back(dis_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ElementReader::Partition()
{
  const int myrank  = comm_->MyPID();
  const int numproc = comm_->NumProc();

  // - read global ids of elements of this discretization
  // - determine a preliminary element distribution
  // - read elements of this discretization and distribute
  // - construct nodal row map with everything on proc 0
  // - construct and partition graph
  // - build final nodal row map

  Epetra_Time time(*comm_);

  vector<int> eids;             // global element ids
  int numele = 0;

  // We read all elements at proc 0 and keep them. This way we do not
  // need to read (initialize) them again.
  //
  // If we happen to consume too much memory here, we have to
  // distribute at this point already and use parmetis.
  if (myrank==0)
  {
    cout << "Entering jumbo reading mode for " << name_ << " discretization ...\n"
         << "Read, create and partition elements      in....";
    fflush(stdout);

    // open input file at the right position
    ifstream file(allfiles.inputfile_name);
    file.seekg(ExcludedSectionPositions[sectionname_]);

    // loop all element lines
    // Comments in the element section are not supported!

    // Ok. We do this twice here! The first time we just gather the
    // element numbers. With those we construct a preliminary element
    // row map.
    string line;
    for (int i=0; getline(file, line); ++i)
    {
      if (line.find("--")==0)
      {
        break;
      }
      else
      {
        istringstream t;
        t.str(line);
        int elenumber;
        t >> elenumber;
        elenumber -= 1;
        eids.push_back(elenumber);
      }
    }
    numele = static_cast<int>(eids.size());
  }

  // Simply allreduce the element ids
  comm_->Broadcast(&numele,1,0);
  eids.resize(numele);
  comm_->Broadcast(&eids[0],numele,0);

  // number of element junks to split the reading process in
  // approximate block size (just a guess!)
  const int nblock = numproc;
  int bsize = max(static_cast<int>(eids.size())/nblock, 1);

  // create a simple (pseudo linear) map
  int mysize = bsize;
  if (myrank==numproc-1)
    mysize = eids.size()-(numproc-1)*bsize;

  roweles_ = rcp(new Epetra_Map(-1,mysize,&eids[myrank*bsize],0,*comm_));
  eids.clear();

  // We need to remember all node ids on processor 0 in order to
  // create the preliminary nodal row map.
  set<int> nodes;

  // For simplicity we remember all node ids of all elements on
  // processor 0. This way we can create the graph on processor 0 and
  // use serial metis. If this turns out to be too memory consuming,
  // we have to use parmetis and use the distributed elements from the
  // discretization.
  list<vector<int> > elementnodes;

  // open input file at correct position,
  // valid on proc 0 only!
  ifstream file;
  if (0==myrank)
  {
    file.open(allfiles.inputfile_name);
    file.seekg(ExcludedSectionPositions[sectionname_]);
  }
  string line;
  int filecount=0;

  // note that the last block is special....
  for (int block=0; block<nblock; ++block)
  {
    if (0==myrank)
    {
      int bcount=0;
      for (;getline(file, line); ++filecount)
      {
        if (line.find("--")==0)
          break;
        else
        {
          istringstream t;
          t.str(line);
          int elenumber;
          string eletype;
          // read element id and type
          t >> elenumber >> eletype;
          elenumber -= 1;

          // Set the current row to the empty slot after the file rows
          // and store the current line. This way the elements can use
          // the normal fr* functions to read the line.
          // Of course this is a hack.
          allfiles.actrow = allfiles.numrows;
          allfiles.actplace = allfiles.input_file[allfiles.actrow] = const_cast<char*>(line.c_str());
          // let the factory create a matching empty element
          RefCountPtr<DRT::Element> ele = DRT::Utils::Factory(eletype,elenumber,0);
          // let this element read its input line
          ele->ReadElement();
          // add element to discretization
          dis_->AddElement(ele);

          //elements_[elenumber] = ele;

          // get the node ids of this element
          const int  numnode = ele->NumNode();
          const int* nodeids = ele->NodeIds();

          copy(nodeids, nodeids+numnode, inserter(nodes, nodes.begin()));
          elementnodes.push_back(vector<int>(nodeids, nodeids+numnode));

          ++bcount;
          if (block != nblock-1) // last block is different....
            if (bcount==bsize)
            {
              filecount++;
              break;
            }
        }
      } // for (;getline(file, line); ++filecount)
    } // if (0==myrank)
    // export block of elements to other processors as reflected in the linear
    // map roweles

    // export junk of elements to other processors
    dis_->ExportRowElements(*roweles_);
  }

  vector<int> nids;             // global node ids

  if (myrank==0)
  {
    // Reset fr* functions. Still required.
    frrewind();

    cout << time.ElapsedTime() << " secs\n";
    time.ResetStartTime();
    cout << "Read, create and partition problem graph in....";
    fflush(stdout);

    copy(nodes.begin(), nodes.end(), back_inserter(nids));
    nodes.clear();
  }

  // create preliminary node row map
  rownodes_ = rcp(new Epetra_Map(-1,nids.size(),&nids[0],0,*comm_));
  nids.clear();

  // construct graph
  RefCountPtr<Epetra_CrsGraph> graph = rcp(new Epetra_CrsGraph(Copy,*rownodes_,81,false));
  if (myrank==0)
  {
    for (list<vector<int> >::iterator i=elementnodes.begin();
         i!=elementnodes.end();
         ++i)
    {
      // get the node ids of this element
      int  numnode = static_cast<int>(i->size());
      int* nodeids = &(*i)[0];

      // loop nodes and add this topology to the row in the graph of every node
      for (int i=0; i<numnode; ++i)
      {
        int err = graph->InsertGlobalIndices(nodeids[i],numnode,nodeids);
        if (err<0) dserror("graph->InsertGlobalIndices returned %d",err);
      }
    }
  }

  // finalize construction of this graph
  int err = graph->FillComplete(*rownodes_,*rownodes_);
  if (err) dserror("graph->FillComplete returned %d",err);

  // partition graph using metis
  Epetra_Vector weights(graph->RowMap(),false);
  weights.PutScalar(1.0);
  RefCountPtr<Epetra_CrsGraph> newgraph = DRT::Utils::PartGraphUsingMetis(*graph,weights);
  graph = newgraph;
  newgraph = null;

  // replace rownodes, colnodes with row and column maps from the graph
  // do stupid conversion from Epetra_BlockMap to Epetra_Map
  const Epetra_BlockMap& brow = graph->RowMap();
  const Epetra_BlockMap& bcol = graph->ColMap();
  rownodes_ = rcp(new Epetra_Map(-1,
                                 brow.NumMyElements(),
                                 brow.MyGlobalElements(),
                                 0,
                                 *comm_));
  colnodes_ = rcp(new Epetra_Map(-1,
                                 bcol.NumMyElements(),
                                 bcol.MyGlobalElements(),
                                 0,
                                 *comm_));

  graph = null;

  // now we have all elements in a linear map roweles
  // build resonable maps for elements from the
  // already valid and final node maps
  // note that nothing is actually redistributed in here
  dis_->BuildElementRowColumn(*rownodes_,*colnodes_,roweles_,coleles_);

  // we can now export elements to resonable row element distribution
  dis_->ExportRowElements(*roweles_);

  // export to the column map / create ghosting of elements
  dis_->ExportColumnElements(*coleles_);

  if (comm_->MyPID()==0)
  {
    cout << time.ElapsedTime() << " secs\n";
    fflush(stdout);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ElementReader::CollapseRowMap()
{
  // collapse rownodes on proc 0 so it knows what to read and what to ignore
  const int* mygids = rownodes_->MyGlobalElements();
  vector<int> sbuff(mygids, mygids+rownodes_->NumMyElements());
  vector<int> rbuff;

  // gather all global ids on proc 0
  int target = 0;
  LINALG::Gather(sbuff,rbuff,1,&target,*comm_);
  collapsedrows_ = rcp(new Epetra_Map(-1,(int)rbuff.size(),&rbuff[0],0,*comm_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ElementReader::Complete()
{
  const int myrank  = comm_->MyPID();

  Epetra_Time time(*comm_);

  if (!myrank)
  {
    cout << "Complete discretization                  in....";
    fflush(stdout);
  }

  int err = dis_->FillComplete();
  if (err)
    dserror("aledis->FillComplete() returned %d",err);

  if (!myrank)
  {
    cout << time.ElapsedTime() << " secs\n";
    fflush(stdout);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NodeReader::NodeReader(RefCountPtr<Epetra_Comm> comm, string sectionname)
  : comm_(comm), sectionname_(sectionname)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
RefCountPtr<DRT::Discretization> NodeReader::FindDisNode(int nodeid)
{
  for (unsigned i=0; i<ereader_.size(); ++i)
  {
    if (ereader_[i]->collapsedrows_->MyGID(nodeid))
    {
      return ereader_[i]->dis_;
    }
  }
  return null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void NodeReader::Read()
{
  const int myrank  = comm_->MyPID();
  const int numproc = comm_->NumProc();

  int numnodes = 0;
  for (unsigned i=0; i<ereader_.size(); ++i)
  {
    ereader_[i]->Partition();
    ereader_[i]->CollapseRowMap();
    numnodes += ereader_[i]->rownodes_->NumGlobalElements();
  }

  Epetra_Time time(*comm_);

  if (!myrank)
  {
    cout << "Read, create and partition nodes         in....";
    fflush(stdout);
  }

  // We will read the nodes block wise. we will use one block per processor
  // so the number of blocks is numproc
  // determine a rough blocksize
  const int nblock = numproc;
  int bsize = max(numnodes/nblock, 1);

  // open input file at the right position
  // note that stream is valid on proc 0 only!
  ifstream file;
  if (myrank==0)
  {
    file.open(allfiles.inputfile_name);
    file.seekg(ExcludedSectionPositions[sectionname_]);
  }
  string tmp;

  // note that the last block is special....
  int filecount=0;
  for (int block=0; block<nblock; ++block)
  {
    if (0==myrank)
    {
      int bcount=0;
      for (; file; ++filecount)
      {
        file >> tmp;

        if (tmp=="NODE")
        {
          double coords[3];
          int nodeid;
          file >> nodeid >> tmp >> coords[0] >> coords[1] >> coords[2];
          nodeid--;
          RefCountPtr<DRT::Discretization> dis = FindDisNode(nodeid);
          if (dis==null)
            continue;
          if (nodeid != filecount)
            dserror("Reading of nodes failed: Nodes must be numbered consecutive!!");
          if (tmp!="COORD")
            dserror("failed to read node %d",nodeid);
          // create node and add to discretization
          RefCountPtr<DRT::Node> node = rcp(new DRT::Node(nodeid,coords,myrank));
          dis->AddNode(node);
          ++bcount;
          if (block != nblock-1) // last block takes all the rest
            if (bcount==bsize)   // block is full
            {
              ++filecount;
              break;
            }
        }
        else if (tmp.find("--")==0)
          break;
        else
          dserror("unexpected word '%s'",tmp.c_str());
      } // for (filecount; file; ++filecount)
    } // if (0==myrank)

    // export block of nodes to other processors as reflected in rownodes,
    // changes ownership of nodes
    for (unsigned i=0; i<ereader_.size(); ++i)
    {
      ereader_[i]->dis_->ExportRowNodes(*ereader_[i]->rownodes_);
    }

  } // for (int block=0; block<nblock; ++block)

  // last thing to do here is to produce nodal ghosting/overlap
  for (unsigned i=0; i<ereader_.size(); ++i)
  {
    ereader_[i]->dis_->ExportColumnNodes(*ereader_[i]->colnodes_);
  }

  if (!myrank)
  {
    cout << time.ElapsedTime() << " secs\n";
    fflush(stdout);
  }

  for (unsigned i=0; i<ereader_.size(); ++i)
  {
    ereader_[i]->Complete();
  }
}


/*----------------------------------------------------------------------*
  | input of control, element and load information         m.gee 10/06  |
  | This version of the routine uses the new discretization subsystem   |
  | ccadiscret                                                          |
 *----------------------------------------------------------------------*/
void ntainp_ccadiscret()
{
  /* the input of the tracing option has not been done yet, so
     we have to make the dstrc_enter 'by hand'
     */
#ifdef DEBUG
  trace.actroutine = trace.actroutine->next;
  trace.actroutine->name = "ntainp";
  trace.actroutine->dsroutcontrol=TRACEROUT::dsin;
  trace.deepness++;
#endif

  /* input of not mesh or time based problem data  */
  inpctr();

  /* input of materials */
  inp_material();
  /* input of multilayer materials -> shell9  (sh 10/02) */
  inp_multimat();

  /* input of fields */
  // depending on the size of the problem we enter 2 different input routines
  // here
  inpfield_ccadiscret_jumbo(); // jumbo mode input for large problems

  // read dynamic control data
  if (genprob.timetyp==time_dynamic) inpctrdyn();

  // read static control data
  else inpctrstat();

  // read input of eigensolution control data
  inpctreig();

  // read all types of geometry related conditions (e.g. boundary conditions)
  // Also read time and space functions and local coord systems
  input_conditions();

  /*-------------------------------------------- input of monitoring data */
  inp_monitor();

#ifdef RESULTTEST
  /*---------------------------------------- input of result descriptions */
  inp_resultdescr();
#endif

  // all reading is done at this point!

  return;
} // end of ntainp_ccadiscret()




/*----------------------------------------------------------------------*
  | input of fields                                        m.gee 03/07  |
  | This version of the routine uses the new discretization subsystem   |
  | ccadiscret                                                          |
 *----------------------------------------------------------------------*/
void inpfield_ccadiscret_jumbo()
{
  DSTraceHelper dst("inpfield_ccadiscret_jumbo");

  int myrank = 0;
  int nproc  = 1;
  fflush(stdout);

#ifdef PARALLEL
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  Epetra_MpiComm* com = new Epetra_MpiComm(MPI_COMM_WORLD);
  RefCountPtr<Epetra_Comm> comm = rcp(com);
#else
  Epetra_SerialComm* com = new Epetra_SerialComm();
  RefCountPtr<Epetra_Comm> comm = rcp(com);
#endif

  genprob.create_dis = 0;
  genprob.create_ale = 0;
  genprob.maxnode    = 0;
  genprob.nodeshift  = genprob.nnode;

  // read elements the first time to create graph object
  // row distribution of nodes
  // column distribution of nodes
  // graph of problem
  RefCountPtr<Epetra_Map> rownodes   = null;
  RefCountPtr<Epetra_Map> colnodes   = null;
  RefCountPtr<Epetra_Map> roweles    = null;
  RefCountPtr<Epetra_Map> coleles    = null;
  RefCountPtr<Epetra_CrsGraph> graph = null;

  RefCountPtr<DRT::Discretization> structdis = null;
  RefCountPtr<DRT::Discretization> fluiddis  = null;
  RefCountPtr<DRT::Discretization> aledis    = null;

#ifdef STRUCT_MULTI
  RefCountPtr<DRT::Discretization> structdis_macro = null;
  RefCountPtr<DRT::Discretization> structdis_micro = null;
  RefCountPtr<DRT::Discretization> structdis_micro_serial = null;
#endif

  if (genprob.probtyp == prb_fsi)
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=3) dserror("numfld != 3 for fsi problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numsf].fieldtyp = structure;
    inpdis(&(field[genprob.numsf]));
    field[genprob.numff].fieldtyp = fluid;
    inpdis(&(field[genprob.numff]));
    field[genprob.numaf].fieldtyp = ale;
    inpdis(&(field[genprob.numaf]));

    NodeReader nodereader(comm, "--NODE COORDS");

    nodereader.AddElementReader(rcp(new ElementReader("Structure", genprob.numsf, comm, "--STRUCTURE ELEMENTS")));
    nodereader.AddElementReader(rcp(new ElementReader("Fluid", genprob.numff, comm, "--FLUID ELEMENTS")));
    nodereader.AddElementReader(rcp(new ElementReader("Ale", genprob.numaf, comm, "--ALE ELEMENTS")));

    nodereader.Read();
  }

  else if (genprob.probtyp==prb_ale)
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=1) dserror("numfld != 1 for ale problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numaf].fieldtyp = ale;
    inpdis(&(field[genprob.numaf]));

    NodeReader nodereader(comm, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new ElementReader("Ale", genprob.numaf, comm, "--ALE ELEMENTS")));
    nodereader.Read();
  }

  else if (genprob.probtyp==prb_fluid)
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=1) dserror("numfld != 1 for fluid problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numff].fieldtyp = fluid;
    inpdis(&(field[genprob.numff]));

    NodeReader nodereader(comm, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new ElementReader("Fluid", genprob.numff, comm, "--FLUID ELEMENTS")));
    nodereader.Read();
  }

  else if (genprob.probtyp==prb_fluid_pm)
    dserror("prb_fluid_pm not yet impl.");

  else if (genprob.probtyp == prb_tsi)
    dserror("prb_tsi not yet impl.");

  else if (genprob.probtyp==prb_structure)
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=1) dserror("numfld != 1 for structural problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numsf].fieldtyp = structure;
    inpdis(&(field[genprob.numsf]));

    NodeReader nodereader(comm, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new ElementReader("Structure", genprob.numsf, comm, "--STRUCTURE ELEMENTS")));
    nodereader.Read();
  } // end of else if (genprob.probtyp==prb_structure)


#ifdef STRUCT_MULTI
  else if (genprob.probtyp==prb_struct_multi)
  {
    // allocate and input general old stuff....
    if (genprob.numfld!=1) dserror("numfld != 1 for structural multi-scale problem");
    field = (FIELD*)CCACALLOC(genprob.numfld,sizeof(FIELD));
    field[genprob.numsf].fieldtyp = structure;
    inpdis(&(field[genprob.numsf]));

    NodeReader nodereader(comm, "--NODE COORDS");
    nodereader.AddElementReader(rcp(new ElementReader("Structure", genprob.numsf, comm, "--STRUCTURE ELEMENTS")));
    nodereader.Read();

    NodeReader micronodereader(comm, "--MICROSTRUCTURE NODE COORDS");
    micronodereader.AddElementReader(rcp(new ElementReader("Micro Structure", genprob.numsf, comm, "--MICROSTRUCTURE ELEMENTS")));
    micronodereader.Read();

    // microscale discretization is distributed over processors but it
    // is needed on every processor redundantly

    vector<RefCountPtr<DRT::Discretization> >* discretization;
    discretization = (vector<RefCountPtr<DRT::Discretization> >*)field[genprob.numsf].ccadis;
    RefCountPtr<Discretization> structdis_micro = (*discretization)[1];

    RefCountPtr<Epetra_SerialComm> serialcomm = rcp(new Epetra_SerialComm());
    structdis_micro_serial = rcp(new DRT::Discretization("Micro Structure Serial",serialcomm));
    (*discretization)[1] = structdis_micro_serial;

    RefCountPtr<Epetra_Map> parallel_rownodes = rcp(new Epetra_Map(*structdis_micro->NodeRowMap()));
    RefCountPtr<Epetra_Map> parallel_roweles  = rcp(new Epetra_Map(*structdis_micro->ElementRowMap()));

    // build redundant colnodes
    vector<int> mygid(parallel_rownodes->NumMyElements());
    for (int i=0; i<parallel_rownodes->NumMyElements(); ++i) mygid[i] = parallel_rownodes->MyGlobalElements()[i];
    vector<int> rmygid(0);
    vector<int> targetprocs(structdis_micro->Comm().NumProc());
    for (int i=0; i<structdis_micro->Comm().NumProc(); ++i) targetprocs[i] = i;
    LINALG::Gather(mygid,rmygid,structdis_micro->Comm().NumProc(),&targetprocs[0], structdis_micro->Comm());
    RefCountPtr<Epetra_Map> redundantmap = rcp(new Epetra_Map(-1,(int)rmygid.size(),&rmygid[0],0,structdis_micro->Comm()));

    // build redundant coleles
    vector<int> mygid_ele(parallel_roweles->NumMyElements());
    for (int i=0; i<parallel_roweles->NumMyElements(); ++i) mygid_ele[i] = parallel_roweles->MyGlobalElements()[i];
    vector<int> rmygid_ele(0);
    LINALG::Gather(mygid_ele,rmygid_ele,structdis_micro->Comm().NumProc(),&targetprocs[0], structdis_micro->Comm());
    RefCountPtr<Epetra_Map> redundantmap_ele = rcp(new Epetra_Map(-1,(int)rmygid_ele.size(),&rmygid_ele[0],0,structdis_micro->Comm()));

    structdis_micro->ExportColumnNodes(*redundantmap);
    structdis_micro->ExportColumnElements(*redundantmap_ele);
    err = structdis_micro->FillComplete();
    if (err) dserror("structdis_micro->FillComplete() returned %d",err);

    for (int i=0; i<structdis_micro->NumMyColElements(); ++i)
    {
      DRT::Element* actele = structdis_micro->lColElement(i);
      RefCountPtr<DRT::Element> newele = rcp(actele->Clone());
      newele->SetOwner(0);
      structdis_micro_serial->AddElement(newele);
    }
    for (int i=0; i<structdis_micro->NumMyColNodes(); ++i)
    {
      DRT::Node* actnode = structdis_micro->lColNode(i);
      RefCountPtr<DRT::Node> newnode = rcp(actnode->Clone());
      newnode->SetOwner(0);
      structdis_micro_serial->AddNode(newnode);
    }
    if (!myrank) cout << "Complete microscale discretization in serial        in...."; fflush(stdout);
    err = structdis_micro_serial->FillComplete();

    if (0)
    {
    for (int i=0; i<structdis_micro->Comm().NumProc(); ++i)
    {
      cout << "\n";
      if (structdis_micro->Comm().MyPID()==i)
      {
        cout << "Proc " << i << "meine serielle Diskretisierung:\n";
        cout << *structdis_micro_serial;
      }
      structdis_micro->Comm().Barrier();
    }
    cout << "\n";
    }

    if (!myrank) cout << time.ElapsedTime() << " secs\n"; fflush(stdout);
  } // end of else if (genprob.probtyp==prb_struct_multi)
#endif


  else dserror("Type of problem unknown");

  return;
} // void inpfield_ccadiscret_jumbo()



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
