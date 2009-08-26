
#ifdef CCADISCRET

#include "drt_nodereader.H"

#include <Epetra_Time.h>

namespace DRT
{
namespace INPUT
{

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NodeReader::NodeReader(const DRT::INPUT::DatFileReader& reader, string sectionname)
  : reader_(reader), comm_(reader.Comm()), sectionname_(sectionname)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Discretization> > NodeReader::FindDisNode(int nodeid)
{
  std::vector<Teuchos::RCP<DRT::Discretization> > v;
  for (unsigned i=0; i<ereader_.size(); ++i)
  {
    if (ereader_[i]->HasNode(nodeid))
    {
      v.push_back(ereader_[i]->MyDis());
    }
  }
  return v;
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
    cout << "Read, create and partition nodes " << flush;
#ifdef PARMETIS
    cout << "block " << flush;
#else
    cout << "        " << flush;
#endif
  }

  // We will read the nodes block wise. we will use one block per processor
  // so the number of blocks is numproc
  // determine a rough blocksize
  int nblock = numproc;
  int bsize = max(numnodes/nblock, 1);

  // for block sizes larger than about 50000 elements (empirical value !)
  // the code sometimes hangs during ExportRowElements
  // Therefore an upper limit for bsize is ensured below.
  int maxblocksize = 50000;

  if (bsize > maxblocksize)
  {
    // without an additional increase of nblock by 1 the last block size
    // could reach a maximum value of (2*maxblocksize)-1, potentially
    // violating the intended upper limit!
    nblock = 1+ static_cast<int>(numnodes/maxblocksize);
    bsize = maxblocksize;
  }

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
#ifdef PARMETIS
      if (!reader_.MyOutputFlag())
      {
        printf("%d",block);
        if(block != nblock-1)
        {
          printf(",");
        }
        else
        {
          printf(" ");
        }
        fflush(stdout);
      }
#endif

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
          if (nodeid != filecount)
            dserror("Reading of nodes failed: Nodes must be numbered consecutive!!");
          if (tmp!="COORD")
            dserror("failed to read node %d",nodeid);
          std::vector<Teuchos::RCP<DRT::Discretization> > diss = FindDisNode(nodeid);
          for (unsigned i=0; i<diss.size(); ++i)
          {
            Teuchos::RCP<DRT::Discretization> dis = diss[i];
            // create node and add to discretization
            Teuchos::RCP<DRT::Node> node = rcp(new DRT::Node(nodeid,coords,myrank));
            dis->AddNode(node);
          }
          ++bcount;
          if (block != nblock-1) // last block takes all the rest
            if (bcount==bsize)   // block is full
            {
              ++filecount;
              break;
            }
        }
        else if (tmp=="CP")
        {
          // read control points for isogeometric analysis
          double coords[3];
          double weight;

          int cpid;
          file >> cpid >> tmp >> coords[0] >> coords[1] >> coords[2] >> weight;
          cpid--;
          if (cpid != filecount)
            dserror("Reading of control points failed: They must be numbered consecutive!!");
          if (tmp!="COORD")
            dserror("failed to read control point %d",cpid);
          std::vector<Teuchos::RCP<DRT::Discretization> > diss = FindDisNode(cpid);
          for (unsigned i=0; i<diss.size(); ++i)
          {
            Teuchos::RCP<DRT::Discretization> dis = diss[i];
            // create node/control point and add to discretization
            Teuchos::RCP<DRT::NURBS::ControlPoint> node = rcp(new DRT::NURBS::ControlPoint(cpid,coords,weight,myrank));
            dis->AddNode(node);
          }
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
    cout << "in...." << time.ElapsedTime() << " secs\n" << endl;
  }

  for (unsigned i=0; i<ereader_.size(); ++i)
  {
    ereader_[i]->Complete();
  }
} // NodeReader::Read

}
}

#endif
