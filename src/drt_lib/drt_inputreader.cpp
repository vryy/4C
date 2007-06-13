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
#include <iterator>

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
  \brief remove all leading and trailing whitespaces from a string.

  Note: consecutive whitespaces inside the string will be reduced to a
  single space.

  \author u.kue
  \date 03/07
*/
/*----------------------------------------------------------------------*/
static std::string trim(const std::string& line)
{
  std::istringstream t;
  std::string s;
  std::string newline;
  t.str(line);
  while (t >> s)
  {
    newline.append(s);
    newline.append(" ");
  }
  if (newline.size()>0)
    newline.resize(newline.size()-1);
  return newline;
}


namespace DRT
{


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DatFileReader::DatFileReader(string filename, RefCountPtr<Epetra_Comm> comm, int outflag)
  : filename_(filename), comm_(comm), outflag_(outflag)
{
  ReadDat();
  DumpInput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
string DatFileReader::MyInputfileName() const
{
  return filename_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int DatFileReader::MyOutputFlag() const
{
  return outflag_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DatFileReader::Activate()
{
  // Publish (some) internal data to the old ccarat reading system
  // Note that these links remain intact even when the reader goes
  // away.

  allfiles.numrows = numrows_;
  allfiles.input_file_hook = &inputfile_[0];
  allfiles.input_file = &lines_[0];
  //allfiles.inputfile_name = const_cast<char*>(filename_.c_str());

  // set fr-system to begin of input_file
  allfiles.actrow = 0;
  allfiles.actplace = &(allfiles.input_file[0][0]);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DatFileReader::ReadDat()
{
  vector<string> exclude;

  exclude.push_back("--NODE COORDS");
  exclude.push_back("--STRUCTURE ELEMENTS");
  exclude.push_back("--FLUID ELEMENTS");
  exclude.push_back("--ALE ELEMENTS");

  int arraysize = 0;

  if (comm_->MyPID()==0)
  {
    ifstream file(filename_.c_str());
    if (not file)
      dserror("unable to open file: %s", filename_.c_str());

    list<string> content;
    bool ignoreline = false;
    string line;

    // loop all input lines
    while (getline(file, line))
    {
      // remove comments
      string::size_type loc = line.find("//");
      if (loc != string::npos)
      {
        line = line.substr(0,loc);
      }

      // remove trailing and leading whitespaces
      // compact internal whitespaces
      line = trim(line);

      // exclude all special sections
      // this includes the section header and all lines in that section
      if (ignoreline)
      {
        if (line.find("--")==0)
        {
          ignoreline = false;
        }
      }

      // Two sections to be ignored can follow each other. We need
      // independent tests.
      if (!ignoreline)
      {
        for (vector<int>::size_type i=0; i<exclude.size(); ++i)
        {
          if (line.find(exclude[i]) != string::npos)
          {
            positions_[exclude[i]] = file.tellg();
            ignoreline = true;
            break;
          }
        }
      }

      // remember line
      if (!ignoreline && line.length() > 0)
      {
        content.push_back(line);
        arraysize += line.length()+1;
      }
    }

    // setup global variables
    numrows_ = static_cast<int>(content.size());

    // allocate space for copy of file
    inputfile_.resize(arraysize);

    // CAUTION: We allocate one more row pointer that necessary. This
    // pointer will be used to point to temporary lines when the
    // excluded section are read (on proc 0). Of course that's just
    // another EVIL HACK. Don't tell anybody.
    lines_.reserve(numrows_+1);

    // fill file buffer
    char* ptr = &inputfile_[0];
    for (list<string>::iterator i=content.begin(); i!=content.end(); ++i)
    {
      strcpy(ptr,i->c_str());
      lines_.push_back(ptr);
      ptr += i->length()+1;
    }

    if (ptr - &inputfile_[0] != static_cast<int>(arraysize))
      dserror("internal error in file read");

    // add the slot for the temporary line...
    lines_.push_back(NULL);
  }

  // Now lets do all the parallel setup. Afterwards all processors
  // have to be the same.

  if (comm_->NumProc()>1)
  {
    /* Now that we use a variable number of bytes per line we have to
     * communicate the buffer size as well. */
    comm_->Broadcast(&arraysize,1,0);
    comm_->Broadcast(&numrows_,1,0);

    if (comm_->MyPID()>0)
    {
      /*--------------------------------------allocate space for copy of file */
      inputfile_.resize(arraysize);
      lines_.reserve(numrows_+1);
    }

#ifdef PARALLEL

    // There are no char based functions available! Do it by hand!
    //comm_->Broadcast(&inputfile_[0],arraysize,0);

    const Epetra_MpiComm& mpicomm = dynamic_cast<const Epetra_MpiComm&>(*comm_);

    MPI_Bcast(&inputfile_[0], arraysize, MPI_CHAR, 0, mpicomm.GetMpiComm());
#else
    dserror("How did you get here? Go away!");
#endif

    /* We have not yet set the row pointers on procs > 0. So do it now. */
    if (comm_->MyPID()>0)
    {
      lines_.push_back(&inputfile_[0]);
      for (int i=0; i<arraysize; ++i)
      {
        if (inputfile_[i]=='\0')
        {
          lines_.push_back(&inputfile_[i+1]);
        }
      }

      if (static_cast<int>(lines_.size()) != numrows_+1)
        dserror("line count mismatch: %d lines expected but %d lines received",
                numrows_+1, lines_.size());
    }

#ifdef PARALLEL
    // distribute excluded section positions
    for (vector<int>::size_type i=0; i<exclude.size(); ++i)
      //comm_->Broadcast(&positions_[exclude[i]],1,0);
      MPI_Bcast(&positions_[exclude[i]],1,MPI_INT,0,mpicomm.GetMpiComm());
#endif
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DatFileReader::DumpInput()
{
#if defined(DEBUG) || defined(OUTPUT_INPUT)
  if (comm_->MyPID()==0)
  {
    fprintf(allfiles.out_err,
            "============================================================================\n"
            "broadcasted copy of input file:\n"
            "============================================================================\n"
      );
    for (unsigned i=0; i<lines_.size()-1; ++i)
    {
      fprintf(allfiles.out_err,"%s\n", lines_[i]);
    }
    fprintf(allfiles.out_err,
            "============================================================================\n"
            "end of broadcasted copy of input file\n"
            "============================================================================\n"
      );
    fflush(allfiles.out_err);
  }
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ElementReader::ElementReader(RefCountPtr<Discretization> dis,
                             const DRT::DatFileReader& reader,
                             string sectionname)
  : name_(dis->Name()),
    reader_(reader),
    comm_(reader.Comm()),
    sectionname_(sectionname),
    dis_(dis)
{
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
  string inputfile_name = reader_.MyInputfileName();

  // We read all elements at proc 0 and keep them. This way we do not
  // need to read (initialize) them again.
  //
  // If we happen to consume too much memory here, we have to
  // distribute at this point already and use parmetis.
  if (myrank==0)
  {
    if (!reader_.MyOutputFlag())
    {
      cout << "Entering jumbo reading mode for " << name_ << " discretization ...\n"
           << "Read, create and partition elements      in....";
      fflush(stdout);
    }

    // open input file at the right position
    ifstream file(inputfile_name.c_str());
    file.seekg(reader_.ExcludedSectionPosition(sectionname_));

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
    file.open(inputfile_name.c_str());
    file.seekg(reader_.ExcludedSectionPosition(sectionname_));
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

          copy(nodeids, nodeids+numnode, inserter(nodes_, nodes_.begin()));
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

    if (!reader_.MyOutputFlag())
      cout << time.ElapsedTime() << " secs\n";
    time.ResetStartTime();
    if (!reader_.MyOutputFlag())
    {
      cout << "Read, create and partition problem graph in....";
      fflush(stdout);
    }

    copy(nodes_.begin(), nodes_.end(), back_inserter(nids));
  }

  // create preliminary node row map
  rownodes_ = rcp(new Epetra_Map(-1,nids.size(),&nids[0],0,*comm_));
  nids.clear();

#if 0
  if (myrank==0)
  {
    cout << "\n\nelementnodes: size=" << elementnodes.size() << endl;
    for (list<vector<int> >::iterator i=elementnodes.begin();
         i!=elementnodes.end();
         ++i)
    {
      copy(i->begin(), i->end(), ostream_iterator<int>(cout, " "));
      cout << endl;
    }
  }
#endif

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

  elementnodes.clear();

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
  rownodes_ = rcp(new Epetra_Map(brow.NumGlobalElements(),
                                 brow.NumMyElements(),
                                 brow.MyGlobalElements(),
                                 0,
                                 *comm_));
  colnodes_ = rcp(new Epetra_Map(bcol.NumGlobalElements(),
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

  if (comm_->MyPID()==0 && reader_.MyOutputFlag() == 0)
  {
    cout << time.ElapsedTime() << " secs\n";
    fflush(stdout);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ElementReader::Complete()
{
  const int myrank  = comm_->MyPID();

  Epetra_Time time(*comm_);

  if (!myrank && !reader_.MyOutputFlag())
  {
    cout << "Complete discretization                  in....";
    fflush(stdout);
  }

  int err = dis_->FillComplete();
  if (err)
    dserror("aledis->FillComplete() returned %d",err);

  if (!myrank && !reader_.MyOutputFlag())
  {
    cout << time.ElapsedTime() << " secs\n";
    fflush(stdout);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NodeReader::NodeReader(const DRT::DatFileReader& reader, string sectionname)
  : reader_(reader), comm_(reader.Comm()), sectionname_(sectionname)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
RefCountPtr<DRT::Discretization> NodeReader::FindDisNode(int nodeid)
{
  for (unsigned i=0; i<ereader_.size(); ++i)
  {
    if (ereader_[i]->HasNode(nodeid))
    {
      return ereader_[i]->MyDis();
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
  string inputfile_name = reader_.MyInputfileName();

  int numnodes = 0;
  for (unsigned i=0; i<ereader_.size(); ++i)
  {
    ereader_[i]->Partition();
    numnodes += ereader_[i]->rownodes_->NumGlobalElements();
  }

  Epetra_Time time(*comm_);

  if (!myrank && !reader_.MyOutputFlag())
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
    file.open(inputfile_name.c_str());
    file.seekg(reader_.ExcludedSectionPosition(sectionname_));
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

  if (!myrank && !reader_.MyOutputFlag())
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
