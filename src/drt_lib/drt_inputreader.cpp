/*----------------------------------------------------------------------*/
/*!
\file drt_inputreader.H

\brief Internal classes to read elements and nodes

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE


#include "drt_inputreader.H"
#include "drt_utils.H"
#include "linalg_utils.H"

#include <Epetra_Time.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

/*!----------------------------------------------------------------------
  \brief file pointers

  <pre>                                                         m.gee 8/00
  This structure struct _FILES allfiles is defined in input_control_global.c
  and the type is in standardtypes.h
  It holds all file pointers and some variables needed for the FRSYSTEM
  </pre>
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*/
/*!
  \brief Section start positions of excluded section in input file

  Another global variable!

  \author u.kue
  \date 03/07
 */
/*----------------------------------------------------------------------*/
extern map<string,ifstream::pos_type> ExcludedSectionPositions;


namespace DRT
{

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

}

#endif
#endif
