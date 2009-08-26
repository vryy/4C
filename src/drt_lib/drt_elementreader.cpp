
#ifdef CCADISCRET

#include "drt_elementreader.H"
#include "standardtypes_cpp.H"
#include "drt_utils.H"
#include "drt_elementdefinition.H"

#include <Epetra_Time.h>

/*!----------------------------------------------------------------------
  \brief file pointers

  <pre>                                                         m.gee 8/00
  This structure struct _FILES allfiles is defined in input_control_global.c
  and the type is in standardtypes.h
  It holds all file pointers and some variables needed for the FRSYSTEM
  </pre>
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

namespace DRT
{
namespace INPUT
{

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ElementReader::ElementReader(Teuchos::RCP<Discretization> dis,
                             const DRT::INPUT::DatFileReader& reader,
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

  Epetra_Time time(*comm_);

  // - read global ids of elements of this discretization
  //   (this is one fully redundant vector for elements)
  // - determine a preliminary element distribution. The fully redundant
  //   vector is trashed after the construction.
  // - define blocksizes for blocks of elements we read (not necessarily
  //   the same as it was used to construct the map --- we may have a
  //   smaller blocksize here).
  // - read elements of this discretization and distribute according
  //   to a linear map. While reading, remember node gids an assemble
  //   them into a second fully redundant vector (mapping vertex id->gid).
  //   In addition, we keep them in a fully redundant set (required by
  //   node reader). Construct reverse lookup from gids to vertex ids.
  //   Again, this is a global, fully redundant map!
  // - define preliminary linear distributed nodal row map
  // - determine adjacency array (i.e. the infos for the node graph)
  //   using the nodal row distribution and a round robin communication
  //   of element connectivity information.
  //   Use adjacency array to build an initial Crsgraph on the linear map.
  // - do partitioning using parmetis
  //   Results are distributed to other procs using two global vectors!
  // - build final nodal row map, export graph to the new map
  //

  // --------------------------------------------------
  // - read global ids of elements of this discretization

  // vector of all global element ids
  vector<int> eids;
  int numele = 0;
  string inputfile_name = reader_.MyInputfileName();

  // all reading is done on proc 0
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

  // --------------------------------------------------
  // - determine a preliminary element distribution

  // number of element junks to split the reading process in
  // approximate block size (just a guess!)
  int nblock = numproc;
  int bsize = static_cast<int>(numele)/nblock;

  // create a simple (pseudo linear) map
  int mysize = bsize;
  if (myrank==numproc-1)
    mysize = numele-(numproc-1)*bsize;

  // construct the map
  roweles_ = rcp(new Epetra_Map(-1,mysize,&eids[myrank*bsize],0,*comm_));

  // throw away redundant vector of elements
  eids.clear();

  // --------------------------------------------------
  // - define blocksizes for blocks of elements we read

  // for block sizes larger than about 250000 elements (empirical value !)
  // the code sometimes hangs during ExportRowElements call for the
  // second block (block 1).
  // Therefore an upper limit of 200000 for bsize is ensured below.
  int maxblocksize = 200000;

  if (bsize > maxblocksize)
  {
    // without an additional increase of nblock by 1 the last block size
    // could reach a maximum value of (2*maxblocksize)-1, potentially
    // violating the intended upper limit!
    nblock = 1+ static_cast<int>(numele/maxblocksize);
    bsize = maxblocksize;
  }

#ifndef PARMETIS
  // For simplicity we remember all node ids of all elements on
  // processor 0. This way we can create the graph on processor 0 and
  // use serial metis. If this turns out to be too memory consuming,
  // we have to use parmetis and use the distributed elements from the
  // discretization.
  list<vector<int> > elementnodes;
#endif

  // --------------------------------------------------
  // - read elements of this discretization and distribute according
  //   to a linear map. While reading, remember node gids an assemble
  //   them into a second fully redundant vector.

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
  bool endofsection = false;

  DRT::INPUT::ElementDefinition ed;
  ed.SetupValidElementLines();

  // note that the last block is special....
  for (int block=0; block<nblock; ++block)
  {
    if (not endofsection and 0==myrank)
    {
      int bcount=0;
      for (;getline(file, line); ++filecount)
      {
        if (line.find("--")==0)
        {
          // If we have an empty element section (or fewer elements
          // that processors) we cannot read on. But we cannot exit
          // the block loop beforehand either because the other
          // processors need to syncronize nblock times.
          endofsection = true;
          break;
        }
        else
        {
          istringstream t;
          t.str(line);
          int elenumber;
          string eletype;
          string distype;
          // read element id type and distype
          t >> elenumber >> eletype >> distype;

#if 0
          DRT::INPUT::LineDefinition* linedef = ed.ElementLines(eletype,distype);
          if (linedef!=NULL)
          {
            if (not linedef->Read(t))
              dserror("failed to read element %d %s %s",elenumber,eletype.c_str(),distype.c_str());
            linedef->Print(std::cout);
            std::cout << "\n";
          }
#endif

          elenumber -= 1;

          // Set the current row to the empty slot after the file rows
          // and store the current line. This way the elements can use
          // the normal fr* functions to read the line.
          // Of course this is a hack.
          allfiles.actrow = allfiles.numrows;
          allfiles.actplace = allfiles.input_file[allfiles.actrow] = const_cast<char*>(line.c_str());
          // let the factory create a matching empty element
          Teuchos::RCP<DRT::Element> ele = DRT::UTILS::Factory(eletype,distype,elenumber,0);
          // let this element read its input line
          ele->ReadElement();
          // add element to discretization
          dis_->AddElement(ele);

          // get the node ids of this element
          const int  numnode = ele->NumNode();
          const int* nodeids = ele->NodeIds();

          // all node gids of this element are inserted into a set of
          // node ids --- it will be used later during reading of nodes
          // to add the node to one or more discretisations
          copy(nodeids, nodeids+numnode, inserter(nodes_, nodes_.begin()));
#ifndef PARMETIS
          elementnodes.push_back(vector<int>(nodeids, nodeids+numnode));
#endif

          ++bcount;
          if (block != nblock-1) // last block is different....
          {
            if (bcount==bsize)
            {
	      filecount++;
              break;
            }
          }
        }
      } // for (;getline(file, line); ++filecount)
    } // if (0==myrank)

    // export junk of elements to other processors as reflected in the linear
    // map roweles
    dis_->ExportRowElements(*roweles_);

  } // end loop blocks

  // global node ids --- this will be a fully redundant vector!
  int numnodes=0;
  vector<int> nids;

  if (myrank==0)
  {
    // Reset fr* functions. Still required.
    frrewind();

    numnodes = (int)nodes_.size();

    // copy set content into nids vector
    nids.reserve(numnodes);
    copy(nodes_.begin(), nodes_.end(), back_inserter(nids));
  }

  // create preliminary node row map
  rownodes_ = rcp(new Epetra_Map(-1,nids.size(),&nids[0],0,*comm_));

  // We want to be able to read empty fields. If we have such a beast
  // just skip the partitioning.
  if (rownodes_->NumGlobalElements()>0)
  {
#ifdef PARMETIS

    // Simply allreduce the node ids --- the vector is ordered according
    // to the < operator from the set which was used to constuct it
    comm_->Broadcast(&numnodes,1,0);
    nids.resize(numnodes);
    comm_->Broadcast(&nids[0],numnodes,0);

    DRT::UTILS::PartUsingParMetis(dis_,roweles_,rownodes_,colnodes_,nids,nblock,
                                  comm_,time,not reader_.MyOutputFlag());

#else
    nids.clear();
    DRT::UTILS::PartUsingMetis(rownodes_,colnodes_,elementnodes,comm_);
#endif
  }
  else
  {
    // We are empty. Just a proper initialization.
    colnodes_ = rownodes_;
  }

  // now we have all elements in a linear map roweles
  // build reasonable maps for elements from the
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

  return;
} // end Partition()


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ElementReader::Complete()
{
  const int myrank  = comm_->MyPID();

  Epetra_Time time(*comm_);

  if (!myrank && !reader_.MyOutputFlag())
  {
    cout << "Complete discretization ";
    printf("%-16s",name_.c_str());
    cout << " in...." << flush;
  }

  int err = dis_->FillComplete(false,false,false);
  if (err)
    dserror("dis_->FillComplete() returned %d",err);

  if (!myrank && !reader_.MyOutputFlag())
  {
    cout << time.ElapsedTime() << " secs" << endl;
  }

  DRT::UTILS::PrintParallelDistribution(*dis_);
}

}
}

#endif
