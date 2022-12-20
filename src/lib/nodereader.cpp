/*---------------------------------------------------------------------*/
/*! \file

\brief Functionality for reading nodes

\level 0


*/
/*---------------------------------------------------------------------*/


#include "discret.H"
#include "nodereader.H"
#include "globalproblem.H"
#include "control_point.H"
#include "immersed_node.H"
#include "fiber_node.H"
#include "io_pstream.H"

#include <Epetra_Time.h>
#include <istream>
#include <string>

namespace DRT
{
  namespace INPUT
  {
    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    NodeReader::NodeReader(const DRT::INPUT::DatFileReader& reader, std::string sectionname)
        : reader_(reader), comm_(reader.Comm()), sectionname_(sectionname)
    {
    }


    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    void NodeReader::AddAdvancedReader(Teuchos::RCP<Discretization> dis,
        const DRT::INPUT::DatFileReader& reader, const std::string& sectionname,
        const std::set<std::string>& elementtypes, const INPAR::GeometryType geometrysource,
        const std::string* geofilepath)
    {
      switch (geometrysource)
      {
        case INPAR::geometry_full:
        {
          std::string fullsectionname("--" + sectionname + " ELEMENTS");
          Teuchos::RCP<ElementReader> er = Teuchos::rcp(
              new DRT::INPUT::ElementReader(dis, reader, fullsectionname, elementtypes));
          ereader_.push_back(er);
          break;
        }
        case INPAR::geometry_box:
        {
          std::string fullsectionname("--" + sectionname + " DOMAIN");
          Teuchos::RCP<DomainReader> dr = Teuchos::rcp(
              new DRT::INPUT::DomainReader(dis, reader, fullsectionname, elementtypes));
          dreader_.push_back(dr);
          break;
        }
        case INPAR::geometry_file:
        {
          dserror("Unfortunately not yet implemented, but feel free ...");
          break;
        }
        default:
          dserror("Unknown geometry source");
          break;
      }
      return;
    }


    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    void NodeReader::AddAdvancedReader(Teuchos::RCP<Discretization> dis,
        const DRT::INPUT::DatFileReader& reader, const std::string& sectionname,
        const INPAR::GeometryType geometrysource, const std::string* geofilepath)
    {
      std::set<std::string> dummy;
      AddAdvancedReader(dis, reader, sectionname, dummy, geometrysource, geofilepath);
      return;
    }


    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    std::vector<Teuchos::RCP<DRT::Discretization>> NodeReader::FindDisNode(int nodeid)
    {
      std::vector<Teuchos::RCP<DRT::Discretization>> v;
      for (unsigned i = 0; i < ereader_.size(); ++i)
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
      const int myrank = comm_->MyPID();
      const int numproc = comm_->NumProc();
      std::string inputfile_name = reader_.MyInputfileName();

      // we need this as the offset for the domain reader
      int maxnodeid = 0;

      // First we let the element reader do their thing
      int numnodes = reader_.ExcludedSectionLength(sectionname_);

      /**************************************************************************
       * first we process all domains that are read from the .dat-file as usual *
       **************************************************************************/

      if (numnodes > 0)  // skip it altogether, if no nodes are in the .dat-file
      {
        for (unsigned i = 0; i < ereader_.size(); ++i)
        {
          ereader_[i]->Partition();
        }

        Epetra_Time time(*comm_);

        if (!myrank && !reader_.MyOutputFlag())
          IO::cout << "Read, create and partition nodes\n" << IO::flush;

        // We will read the nodes block wise. we will use one block per processor
        // so the number of blocks is numproc
        // OR number of blocks is numnodes if less nodes than procs are read in
        // determine a rough blocksize
        int nblock = std::min(numproc, numnodes);
        int bsize = std::max(numnodes / nblock, 1);

        // an upper limit for bsize
        int maxblocksize = 200000;

        if (bsize > maxblocksize)
        {
          // without an additional increase of nblock by 1 the last block size
          // could reach a maximum value of (2*maxblocksize)-1, potentially
          // violating the intended upper limit!
          nblock = 1 + static_cast<int>(numnodes / maxblocksize);
          bsize = maxblocksize;
        }

        // open input file at the right position
        // note that stream is valid on proc 0 only!
        std::ifstream file;
        if (myrank == 0)
        {
          file.open(inputfile_name.c_str());
          file.seekg(reader_.ExcludedSectionPosition(sectionname_));
        }
        std::string tmp;
        std::string tmp2;

        if (!myrank && !reader_.MyOutputFlag())
        {
          printf("numnode %d nblock %d bsize %d\n", numnodes, nblock, bsize);
          fflush(stdout);
        }


        // note that the last block is special....
        int filecount = 0;
        for (int block = 0; block < nblock; ++block)
        {
          double t1 = time.ElapsedTime();
          if (0 == myrank)
          {
            if (!reader_.MyOutputFlag()) printf("block %d ", block);

            int bcount = 0;
            for (; file; ++filecount)
            {
              file >> tmp;

              if (tmp == "NODE")
              {
                double coords[6];
                int nodeid;
                // read in the node coordinates
                file >> nodeid >> tmp >> coords[0] >> coords[1] >> coords[2];
                // store current position of file reader
                int length = file.tellg();
                file >> tmp2;
                nodeid--;
                maxnodeid = std::max(maxnodeid, nodeid) + 1;
                std::vector<Teuchos::RCP<DRT::Discretization>> diss = FindDisNode(nodeid);
                if (tmp2 != "ROTANGLE")  // Common (Boltzmann) Nodes with 3 DoFs
                {
                  // go back with file reader in order to make the expression in tmp2 available to
                  // tmp in the next iteration step
                  file.seekg(length);
                  for (unsigned i = 0; i < diss.size(); ++i)
                  {
                    // create node and add to discretization
                    Teuchos::RCP<DRT::Node> node =
                        Teuchos::rcp(new DRT::Node(nodeid, coords, myrank));
                    diss[i]->AddNode(node);
                  }
                }
                else  // Cosserat Nodes with 6 DoFs
                {
                  // read in the node ritations in case of cosserat nodes
                  file >> coords[3] >> coords[4] >> coords[5];
                  for (unsigned i = 0; i < diss.size(); ++i)
                  {
                    // create node and add to discretization
                    Teuchos::RCP<DRT::Node> node =
                        Teuchos::rcp(new DRT::Node(nodeid, coords, myrank, true));
                    diss[i]->AddNode(node);
                  }
                }
                ++bcount;
                if (block != nblock - 1)  // last block takes all the rest
                  if (bcount == bsize)    // block is full
                  {
                    ++filecount;
                    break;
                  }
              }
              // this is a specialized node for immersed problems
              else if (tmp == "INODE")
              {
                double coords[6];
                int nodeid;
                // read in the node coordinates
                file >> nodeid >> tmp >> coords[0] >> coords[1] >> coords[2];
                // store current position of file reader
                int length = file.tellg();
                file >> tmp2;
                nodeid--;
                maxnodeid = std::max(maxnodeid, nodeid) + 1;
                std::vector<Teuchos::RCP<DRT::Discretization>> diss = FindDisNode(nodeid);

                if (tmp2 != "ROTANGLE")  // Common (Boltzmann) Nodes with 3 DoFs
                {
                  // go back with file reader in order to make the expression in tmp2 available to
                  // tmp in the next iteration step
                  file.seekg(length);
                  for (unsigned i = 0; i < diss.size(); ++i)
                  {
                    // create node and add to discretization
                    Teuchos::RCP<DRT::Node> node =
                        Teuchos::rcp(new IMMERSED::ImmersedNode(nodeid, coords, myrank));
                    diss[i]->AddNode(node);
                  }
                }
                else  // Cosserat Nodes with 6 DoFs
                {
                  dserror("no valid immersed node definition");
                }
                ++bcount;
                if (block != nblock - 1)  // last block takes all the rest
                  if (bcount == bsize)    // block is full
                  {
                    ++filecount;
                    break;
                  }
              }
              // this node is a Nurbs control point
              else if (tmp == "CP")
              {
                // read control points for isogeometric analysis (Nurbs)
                double coords[3];
                double weight;

                int cpid;
                file >> cpid >> tmp >> coords[0] >> coords[1] >> coords[2] >> weight;
                cpid--;
                maxnodeid = std::max(maxnodeid, cpid) + 1;
                if (cpid != filecount)
                  dserror("Reading of control points failed: They must be numbered consecutive!!");
                if (tmp != "COORD") dserror("failed to read control point %d", cpid);
                std::vector<Teuchos::RCP<DRT::Discretization>> diss = FindDisNode(cpid);
                for (unsigned i = 0; i < diss.size(); ++i)
                {
                  Teuchos::RCP<DRT::Discretization> dis = diss[i];
                  // create node/control point and add to discretization
                  Teuchos::RCP<DRT::NURBS::ControlPoint> node =
                      Teuchos::rcp(new DRT::NURBS::ControlPoint(cpid, coords, weight, myrank));
                  dis->AddNode(node);
                }
                ++bcount;
                if (block != nblock - 1)  // last block takes all the rest
                  if (bcount == bsize)    // block is full
                  {
                    ++filecount;
                    break;
                  }
              }
              // this is a special node with additional fiber information
              else if (tmp == "FNODE")
              {
                enum class FiberType
                {
                  Unknown,
                  Angle,
                  Fiber,
                  CosyDirection
                };

                // read fiber node
                std::array<double, 3> coords = {0.0, 0.0, 0.0};
                std::map<FIBER::CoordinateSystemDirection, std::array<double, 3>> cosyDirections;
                std::vector<std::array<double, 3>> fibers;
                std::map<FIBER::AngleType, double> angles;

                int nodeid;
                // read in the node coordinates and fiber direction
                file >> nodeid >> tmp >> coords[0] >> coords[1] >> coords[2];
                nodeid--;
                maxnodeid = std::max(maxnodeid, nodeid) + 1;

                while (true)
                {
                  // store current position of file reader
                  std::ifstream::pos_type length = file.tellg();
                  // try to read new fiber direction or coordinate system
                  file >> tmp2;

                  DRT::FIBER::CoordinateSystemDirection coordinateSystemDirection;
                  DRT::FIBER::AngleType angleType;
                  FiberType type = FiberType::Unknown;

                  if (tmp2 == "FIBER" + std::to_string(1 + fibers.size()))
                  {
                    type = FiberType::Fiber;
                  }
                  else if (tmp2 == "CIR")
                  {
                    coordinateSystemDirection = DRT::FIBER::CoordinateSystemDirection::Circular;
                    type = FiberType::CosyDirection;
                  }
                  else if (tmp2 == "TAN")
                  {
                    coordinateSystemDirection = DRT::FIBER::CoordinateSystemDirection::Tangential;
                    type = FiberType::CosyDirection;
                  }
                  else if (tmp2 == "RAD")
                  {
                    coordinateSystemDirection = DRT::FIBER::CoordinateSystemDirection::Radial;
                    type = FiberType::CosyDirection;
                  }
                  else if (tmp2 == "HELIX")
                  {
                    angleType = DRT::FIBER::AngleType::Helix;
                    type = FiberType::Angle;
                  }
                  else if (tmp2 == "TRANS")
                  {
                    angleType = DRT::FIBER::AngleType::Transverse;
                    type = FiberType::Angle;
                  }
                  else
                  {
                    // No more fiber information. Jump to last position.
                    file.seekg(length);
                    break;
                  }

                  // add fiber / angle to the map
                  switch (type)
                  {
                    case FiberType::Unknown:
                    {
                      dserror(
                          "Unknown fiber node attribute. Numbered fibers must be in order, i.e. "
                          "FIBER1, FIBER2, ...");
                    }
                    case FiberType::Angle:
                    {
                      file >> angles[angleType];
                      break;
                    }
                    case FiberType::Fiber:
                    {
                      std::array<double, 3> fiber_components;
                      file >> fiber_components[0] >> fiber_components[1] >> fiber_components[2];
                      fibers.emplace_back(fiber_components);
                      break;
                    }
                    case FiberType::CosyDirection:
                    {
                      file >> cosyDirections[coordinateSystemDirection][0] >>
                          cosyDirections[coordinateSystemDirection][1] >>
                          cosyDirections[coordinateSystemDirection][2];
                      break;
                    }
                    default:
                      dserror("Unknown number of components");
                  }
                }

                // add fiber information to node
                std::vector<Teuchos::RCP<DRT::Discretization>> discretizations =
                    FindDisNode(nodeid);
                for (auto& dis : discretizations)
                {
                  auto node = Teuchos::rcp(new DRT::FIBER::FiberNode(
                      nodeid, coords, cosyDirections, fibers, angles, myrank));
                  dis->AddNode(node);
                }

                ++bcount;
                if (block != nblock - 1)  // last block takes all the rest
                {
                  if (bcount == bsize)  // block is full
                  {
                    ++filecount;
                    break;
                  }
                }
              }
              else if (tmp.find("--") == 0)
                break;
              else
                dserror("unexpected word '%s'", tmp.c_str());
            }  // for (filecount; file; ++filecount)
          }    // if (0==myrank)

          double t2 = time.ElapsedTime();
          if (!myrank && !reader_.MyOutputFlag()) printf("reading %10.5e secs", t2 - t1);

          // export block of nodes to other processors as reflected in rownodes,
          // changes ownership of nodes
          for (unsigned i = 0; i < ereader_.size(); ++i)
          {
            ereader_[i]->dis_->ProcZeroDistributeNodesToAll(*ereader_[i]->rownodes_);
            // this does the same job but slower
            // ereader_[i]->dis_->ExportRowNodes(*ereader_[i]->rownodes_);
          }
          double t3 = time.ElapsedTime();
          if (!myrank && !reader_.MyOutputFlag())
          {
            printf(" / distrib %10.5e secs\n", t3 - t2);
            fflush(stdout);
          }

        }  // for (int block=0; block<nblock; ++block)

        // last thing to do here is to produce nodal ghosting/overlap
        for (unsigned i = 0; i < ereader_.size(); ++i)
        {
          ereader_[i]->dis_->ExportColumnNodes(*ereader_[i]->colnodes_);
        }

        if (!myrank && !reader_.MyOutputFlag())
          printf(
              "in............................................. %10.5e secs\n", time.ElapsedTime());

        for (unsigned i = 0; i < ereader_.size(); ++i)
        {
          ereader_[i]->Complete();
        }
      }

      /********************************************************
       * now we process all domains that use the domainreader *
       ********************************************************/
      for (size_t i = 0; i < dreader_.size(); ++i)
      {
        // communicate node offset to all procs
        int lmaxnodeid = maxnodeid;
        comm_->MaxAll(&lmaxnodeid, &maxnodeid, 1);
        dreader_[i]->Partition(&maxnodeid);
        maxnodeid++;
      }

      for (size_t i = 0; i < dreader_.size(); ++i)
      {
        dreader_[i]->Complete();
      }

      int lmaxnodeid = maxnodeid;
      comm_->MaxAll(&lmaxnodeid, &maxnodeid, 1);
      if ((maxnodeid < numproc) && (numnodes != 0))
        dserror("Bad idea: Simulation with %d procs for problem with %d nodes", numproc, maxnodeid);

    }  // NodeReader::Read

  }  // namespace INPUT
}  // namespace DRT
