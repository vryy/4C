/*----------------------------------------------------------------------*/
/*! \file

\brief Read domain sections of dat files.

\level 0


*/
/*----------------------------------------------------------------------*/

#include "lib_domainreader.H"
#include "lib_gridgenerator.H"
#include "lib_elementdefinition.H"
#include "lib_utils_parallel.H"
#include "lib_utils_reader.H"
#include "lib_discret.H"
#include "lib_parobject.H"
#include "io_pstream.H"
#include "rebalance_utils.H"

#include <Teuchos_Time.hpp>
#include <algorithm>

namespace DRT
{
  namespace INPUT
  {
    /// forward declarations
    void BroadcastInputDataToAllProcs(
        Teuchos::RCP<Epetra_Comm> comm, DRT::GRIDGENERATOR::RectangularCuboidInputs& inputData);

    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    DomainReader::DomainReader(Teuchos::RCP<Discretization> dis,
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
    void DomainReader::CreatePartitionedMesh(int nodeGIdOfFirstNewNode) const
    {
      const int myrank = comm_->MyPID();

      Teuchos::Time time("", true);

      if (!reader_.MyOutputFlag() && myrank == 0)
        IO::cout << "Entering domain generation mode for " << name_
                 << " discretization ...\nCreate and partition elements      in...." << IO::endl;

      DRT::GRIDGENERATOR::RectangularCuboidInputs inputData = ReadRectangularCuboidInputData();
      inputData.node_gid_of_first_new_node_ = nodeGIdOfFirstNewNode;

      DRT::GRIDGENERATOR::CreateRectangularCuboidDiscretization(
          *dis_, inputData, static_cast<bool>(reader_.MyOutputFlag()));

      if (!myrank && reader_.MyOutputFlag() == 0)
        IO::cout << "............................................... " << std::setw(10)
                 << std::setprecision(5) << std::scientific << time.totalElapsedTime(true)
                 << " secs" << IO::endl;

      return;
    }

    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    DRT::GRIDGENERATOR::RectangularCuboidInputs DomainReader::ReadRectangularCuboidInputData() const
    {
      DRT::GRIDGENERATOR::RectangularCuboidInputs inputData;
      // all reading is done on proc 0
      if (comm_->MyPID() == 0)
      {
        // open input file at the right position
        std::string inputFileName = reader_.MyInputfileName();
        std::ifstream file(inputFileName.c_str());
        std::ifstream::pos_type pos = reader_.ExcludedSectionPosition(sectionname_);
        if (pos != std::ifstream::pos_type(-1))
        {
          file.seekg(pos);

          // read domain info
          std::string line;
          for (int i = 0; getline(file, line); ++i)
          {
            // remove comments, trailing and leading whitespaces
            // compact internal whitespaces
            line = DRT::UTILS::StripComment(line);

            // line is now empty
            if (line.size() == 0) continue;

            if (line.find("--") == 0)
            {
              break;
            }
            else
            {
              std::istringstream t;
              t.str(line);
              std::string key;
              t >> key;
              if (key == "LOWER_BOUND")
                t >> inputData.bottom_corner_point_[0] >> inputData.bottom_corner_point_[1] >>
                    inputData.bottom_corner_point_[2];
              else if (key == "UPPER_BOUND")
                t >> inputData.top_corner_point_[0] >> inputData.top_corner_point_[1] >>
                    inputData.top_corner_point_[2];
              else if (key == "INTERVALS")
                t >> inputData.interval_[0] >> inputData.interval_[1] >> inputData.interval_[2];
              else if (key == "ROTATION")
                t >> inputData.rotation_angle_[0] >> inputData.rotation_angle_[1] >>
                    inputData.rotation_angle_[2];
              else if (key == "ELEMENTS")
              {
                t >> inputData.elementtype_ >> inputData.distype_;
                getline(t, inputData.elearguments_);
              }
              else if (key == "PARTITION")
              {
                std::string tmp;
                t >> tmp;
                std::transform(tmp.begin(), tmp.end(), tmp.begin(),
                    [](unsigned char c) { return std::tolower(c); });
                if (tmp == "auto")
                  inputData.autopartition_ = true;
                else if (tmp == "structured")
                  inputData.autopartition_ = false;
                else
                  dserror(
                      "Invalid argument for PARTITION in DOMAIN reader. Valid options are \"auto\" "
                      "and \"structured\".");
              }
              else
                dserror("Unknown Key in DOMAIN section");
            }
          }
        }
        else
        {
          dserror("No DOMAIN specified but box geometry selected!");
        }
      }

      // broadcast if necessary
      if (comm_->NumProc() > 1)
      {
        BroadcastInputDataToAllProcs(comm_, inputData);
      }

      return inputData;
    }

    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    void BroadcastInputDataToAllProcs(
        Teuchos::RCP<Epetra_Comm> comm, DRT::GRIDGENERATOR::RectangularCuboidInputs& inputData)
    {
      const int myrank = comm->MyPID();

      std::vector<char> data;
      if (myrank == 0)
      {
        DRT::PackBuffer buffer;
        DRT::ParObject::AddtoPack<double, 3>(buffer, inputData.bottom_corner_point_);
        DRT::ParObject::AddtoPack<double, 3>(buffer, inputData.top_corner_point_);
        DRT::ParObject::AddtoPack<int, 3>(buffer, inputData.interval_);
        DRT::ParObject::AddtoPack<double, 3>(buffer, inputData.rotation_angle_);
        DRT::ParObject::AddtoPack(buffer, static_cast<int>(inputData.autopartition_));
        DRT::ParObject::AddtoPack(buffer, inputData.elementtype_);
        DRT::ParObject::AddtoPack(buffer, inputData.distype_);
        DRT::ParObject::AddtoPack(buffer, inputData.elearguments_);
        buffer.StartPacking();
        DRT::ParObject::AddtoPack<double, 3>(buffer, inputData.bottom_corner_point_);
        DRT::ParObject::AddtoPack<double, 3>(buffer, inputData.top_corner_point_);
        DRT::ParObject::AddtoPack<int, 3>(buffer, inputData.interval_);
        DRT::ParObject::AddtoPack<double, 3>(buffer, inputData.rotation_angle_);
        DRT::ParObject::AddtoPack(buffer, static_cast<int>(inputData.autopartition_));
        DRT::ParObject::AddtoPack(buffer, inputData.elementtype_);
        DRT::ParObject::AddtoPack(buffer, inputData.distype_);
        DRT::ParObject::AddtoPack(buffer, inputData.elearguments_);
        std::swap(data, buffer());
      }

      ssize_t data_size = data.size();
      comm->Broadcast(&data_size, 1, 0);
      if (myrank != 0) data.resize(data_size, 0);
      comm->Broadcast(data.data(), data.size(), 0);

      if (myrank != 0)
      {
        size_t pos = 0;
        DRT::ParObject::ExtractfromPack<double, 3>(pos, data, inputData.bottom_corner_point_);
        DRT::ParObject::ExtractfromPack<double, 3>(pos, data, inputData.top_corner_point_);
        DRT::ParObject::ExtractfromPack<int, 3>(pos, data, inputData.interval_);
        DRT::ParObject::ExtractfromPack<double, 3>(pos, data, inputData.rotation_angle_);
        int autopartitionInteger;
        DRT::ParObject::ExtractfromPack(pos, data, autopartitionInteger);
        inputData.autopartition_ = autopartitionInteger;
        DRT::ParObject::ExtractfromPack(pos, data, inputData.elementtype_);
        DRT::ParObject::ExtractfromPack(pos, data, inputData.distype_);
        DRT::ParObject::ExtractfromPack(pos, data, inputData.elearguments_);
      }
    }

    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    void DomainReader::Complete() const
    {
      const int myrank = comm_->MyPID();

      Teuchos::Time time("", true);

      if (!myrank && !reader_.MyOutputFlag())
        IO::cout << "Complete discretization " << std::left << std::setw(16) << name_ << " in...."
                 << IO::flush;

      int err = dis_->FillComplete(false, false, false);
      if (err) dserror("dis_->FillComplete() returned %d", err);

      if (!myrank && !reader_.MyOutputFlag())
        IO::cout << time.totalElapsedTime(true) << " secs" << IO::endl;

      REBALANCE::UTILS::PrintParallelDistribution(*dis_);
    }

  }  // namespace INPUT
}  // namespace DRT
