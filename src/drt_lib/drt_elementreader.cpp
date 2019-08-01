/*----------------------------------------------------------------------*/
/*!

\brief Read element sections of dat files.

\level 0

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/


#include "drt_elementreader.H"
#include "standardtypes_cpp.H"
#include "drt_elementdefinition.H"
#include "drt_globalproblem.H"
#include "drt_utils_parmetis.H"
#include "drt_utils_factory.H"
#include "drt_utils_parallel.H"

#include <Epetra_Time.h>


namespace DRT
{
  namespace INPUT
  {
    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    ElementReader::ElementReader(Teuchos::RCP<Discretization> dis,
        const DRT::INPUT::DatFileReader& reader, std::string sectionname)
        : name_(dis->Name()),
          reader_(reader),
          comm_(reader.Comm()),
          sectionname_(sectionname),
          dis_(dis)
    {
    }


    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    ElementReader::ElementReader(Teuchos::RCP<Discretization> dis,
        const DRT::INPUT::DatFileReader& reader, std::string sectionname, std::string elementtype)
        : name_(dis->Name()),
          reader_(reader),
          comm_(reader.Comm()),
          sectionname_(sectionname),
          dis_(dis)
    {
      elementtypes_.insert(elementtype);
    }


    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    ElementReader::ElementReader(Teuchos::RCP<Discretization> dis,
        const DRT::INPUT::DatFileReader& reader, std::string sectionname,
        const std::set<std::string>& elementtypes)
        : name_(dis->Name()),
          reader_(reader),
          comm_(reader.Comm()),
          sectionname_(sectionname),
          dis_(dis)
    {
      std::copy(elementtypes.begin(), elementtypes.end(),
          std::inserter(elementtypes_, elementtypes_.begin()));
    }


    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    void ElementReader::Partition()
    {
      const int myrank = comm_->MyPID();
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
      std::vector<int> eids;
      int numele = 0;
      std::string inputfile_name = reader_.MyInputfileName();

      // all reading is done on proc 0
      if (myrank == 0)
      {
        if (!reader_.MyOutputFlag())
        {
          std::cout << "Entering jumbo reading mode for " << name_ << " discretization ...\n"
                    << "Read, create and partition elements      in....\n";
          fflush(stdout);
        }

        // open input file at the right position
        std::ifstream file(inputfile_name.c_str());
        std::ifstream::pos_type pos = reader_.ExcludedSectionPosition(sectionname_);
        if (pos != std::ifstream::pos_type(-1))
        {
          file.seekg(pos);

          // loop all element lines
          // Comments in the element section are not supported!

          // Ok. We do this twice here! The first time we just gather the
          // element numbers. With those we construct a preliminary element
          // row map.
          std::string line;
          for (int i = 0; getline(file, line); ++i)
          {
            if (line.find("--") == 0)
            {
              break;
            }
            else
            {
              std::istringstream t;
              t.str(line);
              int elenumber;
              std::string eletype;
              t >> elenumber >> eletype;
              elenumber -= 1;

              // only read registered element types or all elements if nothing is
              // registered
              if (elementtypes_.size() == 0 or elementtypes_.count(eletype) > 0)
              {
                eids.push_back(elenumber);
              }
            }
          }
          numele = static_cast<int>(eids.size());
        }
        else
        {
          // No such section. Leave as soon as possible.
        }
      }

      // Simply allreduce the element ids
      comm_->Broadcast(&numele, 1, 0);

      if (numele == 0)
      {
        // This is it. Build an empty reader and leave.
        coleles_ = roweles_ = colnodes_ = rownodes_ =
            Teuchos::rcp(new Epetra_Map(-1, 0, NULL, 0, *comm_));
        if (comm_->MyPID() == 0 && reader_.MyOutputFlag() == 0)
        {
          std::cout << time.ElapsedTime() << " secs\n";
          fflush(stdout);
        }
        return;
      }

      eids.resize(numele);
      comm_->Broadcast(&eids[0], numele, 0);

      // --------------------------------------------------
      // - determine a preliminary element distribution

      // number of element chunks to split the reading process in
      // approximate block size (just a guess!)
      int nblock = numproc;
      int bsize = static_cast<int>(numele) / nblock;

      // create a simple (pseudo linear) map
      int mysize = bsize;
      if (myrank == numproc - 1) mysize = numele - (numproc - 1) * bsize;

      // construct the map
      roweles_ = Teuchos::rcp(new Epetra_Map(-1, mysize, &eids[myrank * bsize], 0, *comm_));

      // throw away redundant vector of elements
      eids.clear();

      // --------------------------------------------------
      // - define blocksizes for blocks of elements we read

      // for block sizes larger than about 250000 elements (empirical value !)
      // the code sometimes hangs during ExportRowElements call for the
      // second block (block 1).
      // Therefore an upper limit of 200000 for bsize is ensured below.
      int maxblocksize = 100000;  // 200000;

      if (bsize > maxblocksize)
      {
        // without an additional increase of nblock by 1 the last block size
        // could reach a maximum value of (2*maxblocksize)-1, potentially
        // violating the intended upper limit!
        nblock = 1 + static_cast<int>(numele / maxblocksize);
        bsize = maxblocksize;
      }

#if !defined(PARALLEL) || !defined(HAVE_PARMETIS)
      // For simplicity we remember all node ids of all elements on
      // processor 0. This way we can create the graph on processor 0 and
      // use serial metis. If this turns out to be too memory consuming,
      // we have to use parmetis and use the distributed elements from the
      // discretization.
      std::list<std::vector<int>> elementnodes;
#endif

      // --------------------------------------------------
      // - read elements of this discretization and distribute according
      //   to a linear map. While reading, remember node gids an assemble
      //   them into a second fully redundant vector.

      // open input file at correct position,
      // valid on proc 0 only!
      std::ifstream file;
      if (0 == myrank)
      {
        file.open(inputfile_name.c_str());
        file.seekg(reader_.ExcludedSectionPosition(sectionname_));
      }
      std::string line;
      int filecount = 0;
      bool endofsection = false;

      DRT::INPUT::ElementDefinition ed;
      ed.SetupValidElementLines();

      if (!myrank)
      {
        printf("numele %d nblock %d bsize %d\n", numele, nblock, bsize);
        fflush(stdout);
      }
      Epetra_Time timer(*comm_);

      // note that the last block is special....
      for (int block = 0; block < nblock; ++block)
      {
        double t1 = timer.ElapsedTime();
        std::vector<int> gidlist;
        if (!endofsection && 0 == myrank)
        {
          gidlist.reserve(bsize);
          int bcount = 0;
          for (; getline(file, line); ++filecount)
          {
            if (line.find("--") == 0)
            {
              // If we have an empty element section (or fewer elements
              // than processors) we cannot read on. But we cannot exit
              // the block loop beforehand either because the other
              // processors need to syncronize nblock times.
              endofsection = true;
              break;
            }
            else
            {
              std::istringstream t;
              t.str(line);
              int elenumber;
              std::string eletype;
              std::string distype;
              // read element id type and distype
              t >> elenumber >> eletype >> distype;
              elenumber -= 1;
              gidlist.push_back(elenumber);

              // only read registered element types or all elements if nothing is
              // registered
              if (elementtypes_.size() == 0 or elementtypes_.count(eletype) > 0)
              {
                // let the factory create a matching empty element
                Teuchos::RCP<DRT::Element> ele =
                    DRT::UTILS::Factory(eletype, distype, elenumber, 0);
                if (ele.is_null()) dserror("element creation failed");

                // For the time being we support old and new input facilities. To
                // smooth transition.

                DRT::INPUT::LineDefinition* linedef = ed.ElementLines(eletype, distype);
                if (linedef != NULL)
                {
                  if (not linedef->Read(t))
                  {
                    std::cout << "\n" << elenumber << " " << eletype << " " << distype << " ";
                    linedef->Print(std::cout);
                    std::cout << "\n";
                    std::cout << line << "\n";
                    dserror("failed to read element %d %s %s", elenumber, eletype.c_str(),
                        distype.c_str());
                  }

                  ele->SetNodeIds(distype, linedef);
                  ele->ReadElement(eletype, distype, linedef);
                }
                else
                {
                  dserror("a matching line definition is needed for %s %s", eletype.c_str(),
                      distype.c_str());
                }

                // add element to discretization
                dis_->AddElement(ele);

                // get the node ids of this element
                const int numnode = ele->NumNode();
                const int* nodeids = ele->NodeIds();

                // all node gids of this element are inserted into a set of
                // node ids --- it will be used later during reading of nodes
                // to add the node to one or more discretisations
                std::copy(nodeids, nodeids + numnode, std::inserter(nodes_, nodes_.begin()));
#if !defined(PARALLEL) || !defined(HAVE_PARMETIS)
                elementnodes.push_back(std::vector<int>(nodeids, nodeids + numnode));
#endif

                ++bcount;
                if (block != nblock - 1)  // last block is different....
                {
                  if (bcount == bsize)
                  {
                    filecount++;
                    break;
                  }
                }
              }
            }
          }  // for (;getline(file, line); ++filecount)
        }    // if (0==myrank)

        double t2 = timer.ElapsedTime();
        if (!myrank) printf("ele block %d reading %10.5e secs / ", block, t2 - t1);

        // export junk of elements to other processors as reflected in the linear
        // map roweles
        dis_->ProcZeroDistributeElementsToAll(*roweles_, gidlist);
        // this also works but is slower
        // dis_->ExportRowElements(*roweles_);
        double t3 = timer.ElapsedTime();
        if (!myrank)
        {
          printf("distrib time %10.5e secs\n", t3 - t2);
          fflush(stdout);
        }

      }  // end loop blocks

      // global node ids --- this will be a fully redundant vector!
      int numnodes = 0;
      std::vector<int> nids;
      numnodes = (int)nodes_.size();
      comm_->Broadcast(&numnodes, 1, 0);

      const double imbalance_tol =
          DRT::Problem::Instance()->MeshPartitioningParams().get<double>("IMBALANCE_TOL");

      // We want to be able to read empty fields. If we have such a beast
      // just skip the partitioning.
      if (numnodes)
      {
#if defined(PARALLEL) && defined(HAVE_PARMETIS)

        rownodes_ = Teuchos::null;
        colnodes_ = Teuchos::null;
        nids.clear();
        DRT::UTILS::RedistributeGraphOfDiscretization(dis_, roweles_, rownodes_, colnodes_, comm_,
            !reader_.MyOutputFlag(), comm_->NumProc(), imbalance_tol);

#else
        nids.clear();
        DRT::UTILS::PartUsingMetis(rownodes_, colnodes_, elementnodes, comm_);
#endif
      }
      else
      {
        // We are empty. Just a proper initialization.
        int zero = 0;
        colnodes_ = Teuchos::rcp(new Epetra_Map(-1, zero, &zero, 0, *comm_));
        rownodes_ = colnodes_;
      }

      // now we have all elements in a linear map roweles
      // build reasonable maps for elements from the
      // already valid and final node maps
      // note that nothing is actually redistributed in here
      dis_->BuildElementRowColumn(*rownodes_, *colnodes_, roweles_, coleles_);

      // we can now export elements to resonable row element distribution
      dis_->ExportRowElements(*roweles_);

      // export to the column map / create ghosting of elements
      dis_->ExportColumnElements(*coleles_);

      if (!myrank && reader_.MyOutputFlag() == 0)
      {
        printf("............................................... %10.5e secs\n", time.ElapsedTime());
        // cout << time.ElapsedTime() << " secs\n";
        fflush(stdout);
      }

      return;
    }  // end Partition()


    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    void ElementReader::Complete()
    {
      const int myrank = comm_->MyPID();

      Epetra_Time time(*comm_);

      if (!myrank && !reader_.MyOutputFlag())
      {
        std::cout << "Complete discretization ";
        printf("%-16s", name_.c_str());
        std::cout << " in...." << std::flush;
      }

      int err = dis_->FillComplete(false, false, false);
      if (err) dserror("dis_->FillComplete() returned %d", err);

      if (!myrank && !reader_.MyOutputFlag())
      {
        std::cout << time.ElapsedTime() << " secs" << std::endl;
      }

      DRT::UTILS::PrintParallelDistribution(*dis_);
    }


    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    OldParticleReader::OldParticleReader(
        Teuchos::RCP<Discretization> dis, const DRT::INPUT::DatFileReader& reader)
        : ElementReader(dis, reader, "--DUMMY ELEMENTS")
    {
    }


    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    void OldParticleReader::Partition()
    {
      const int myrank = comm_->MyPID();
      const int numproc = comm_->NumProc();

      // --------------------------------------------------
      // - read global ids of particles of this discretization

      // vector of all global particle ids
      std::vector<int> particleids;
      int numparticles = 0;
      std::string inputfile_name = reader_.MyInputfileName();

      // all reading is done on proc 0
      if (myrank == 0)
      {
        if (!reader_.MyOutputFlag())
        {
          std::cout << "Preliminary node maps are created for " << name_ << " discretization ...\n";
          fflush(stdout);
        }

        // open input file at the node section
        std::ifstream file(inputfile_name.c_str());
        std::ifstream::pos_type pos = reader_.ExcludedSectionPosition("--NODE COORDS");
        if (pos != std::ifstream::pos_type(-1))
        {
          file.seekg(pos);

          // loop all node lines
          // Comments in the node section are not supported!

          // Here we construct a preliminary particle row map.
          std::string line;
          for (int i = 0; getline(file, line); ++i)
          {
            if (line.find("--") == 0)
            {
              break;
            }
            else
            {
              std::istringstream t;
              t.str(line);
              int nodeid;
              std::string nodetype;
              t >> nodetype >> nodeid;
              nodeid -= 1;

              // this node is a particle
              if (nodetype == "PARTICLE" or nodetype == "RPARTICLE" or nodetype == "ELLIPSOID")
              {
                particleids.push_back(nodeid);
              }
            }
          }
          numparticles = static_cast<int>(particleids.size());
        }
        else
        {
          dserror("--NODE COORDS section is missing in input file");
        }
      }

      // Simply allreduce the particle ids
      comm_->Broadcast(&numparticles, 1, 0);

      particleids.resize(numparticles);
      comm_->Broadcast(&particleids[0], numparticles, 0);

      // --------------------------------------------------
      // - determine a preliminary particle distribution

      // number of element junks to split the reading process in
      // approximate block size (just a guess!)
      int nblock = numproc;
      int bsize = static_cast<int>(numparticles) / nblock;

      // create a simple (pseudo linear) map
      int mysize = bsize;
      if (myrank == numproc - 1) mysize = numparticles - (numproc - 1) * bsize;

      // construct the map
      rownodes_ = Teuchos::rcp(new Epetra_Map(-1, mysize, &particleids[myrank * bsize], 0, *comm_));

      // so far there is no information about ghosting available
      colnodes_ = Teuchos::rcp(new Epetra_Map(*rownodes_));

      // throw away redundant vector of particles
      particleids.clear();

      return;
    }

  }  // namespace INPUT
}  // namespace DRT
