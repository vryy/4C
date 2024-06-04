/*----------------------------------------------------------------------*/
/*! \file

\brief Read domain sections of dat files.

\level 0


*/
/*----------------------------------------------------------------------*/

#include "4C_io_domainreader.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_discretization_fem_general_element_definition.hpp"
#include "4C_io_gridgenerator.hpp"
#include "4C_io_pstream.hpp"
#include "4C_io_utils_reader.hpp"
#include "4C_lib_discret.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_rebalance_print.hpp"

#include <Teuchos_Time.hpp>

#include <algorithm>

FOUR_C_NAMESPACE_OPEN

namespace CORE::IO
{
  /// forward declarations
  void BroadcastInputDataToAllProcs(
      Teuchos::RCP<Epetra_Comm> comm, CORE::IO::GRIDGENERATOR::RectangularCuboidInputs& inputData);

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  DomainReader::DomainReader(Teuchos::RCP<DRT::Discretization> dis,
      const CORE::IO::DatFileReader& reader, std::string sectionname)
      : name_(dis->Name()),
        reader_(reader),
        comm_(reader.Comm()),
        sectionname_(sectionname),
        dis_(dis)
  {
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void DomainReader::create_partitioned_mesh(int nodeGIdOfFirstNewNode) const
  {
    const int myrank = comm_->MyPID();

    Teuchos::Time time("", true);

    if (!reader_.MyOutputFlag() && myrank == 0)
      CORE::IO::cout << "Entering domain generation mode for " << name_
                     << " discretization ...\nCreate and partition elements      in...."
                     << CORE::IO::endl;

    GRIDGENERATOR::RectangularCuboidInputs inputData =
        DomainReader::read_rectangular_cuboid_input_data();
    inputData.node_gid_of_first_new_node_ = nodeGIdOfFirstNewNode;

    CORE::IO::GRIDGENERATOR::CreateRectangularCuboidDiscretization(
        *dis_, inputData, static_cast<bool>(reader_.MyOutputFlag()));

    if (!myrank && reader_.MyOutputFlag() == 0)
      CORE::IO::cout << "............................................... " << std::setw(10)
                     << std::setprecision(5) << std::scientific << time.totalElapsedTime(true)
                     << " secs" << CORE::IO::endl;

    return;
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  GRIDGENERATOR::RectangularCuboidInputs DomainReader::read_rectangular_cuboid_input_data() const
  {
    CORE::IO::GRIDGENERATOR::RectangularCuboidInputs inputData;
    // all reading is done on proc 0
    if (comm_->MyPID() == 0)
    {
      // open input file at the right position
      std::string inputFileName = reader_.MyInputfileName();
      std::ifstream file(inputFileName.c_str());
      std::ifstream::pos_type pos = reader_.excluded_section_position(sectionname_);
      if (pos != std::ifstream::pos_type(-1))
      {
        file.seekg(pos);

        // read domain info
        std::string line;
        while (getline(file, line))
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
                FOUR_C_THROW(
                    "Invalid argument for PARTITION in DOMAIN reader. Valid options are \"auto\" "
                    "and \"structured\".");
            }
            else
              FOUR_C_THROW("Unknown Key in DOMAIN section");
          }
        }
      }
      else
      {
        FOUR_C_THROW("No DOMAIN specified but box geometry selected!");
      }
    }

    // broadcast if necessary
    if (comm_->NumProc() > 1)
    {
      IO::BroadcastInputDataToAllProcs(comm_, inputData);
    }

    return inputData;
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void BroadcastInputDataToAllProcs(
      Teuchos::RCP<Epetra_Comm> comm, CORE::IO::GRIDGENERATOR::RectangularCuboidInputs& inputData)
  {
    const int myrank = comm->MyPID();

    std::vector<char> data;
    if (myrank == 0)
    {
      CORE::COMM::PackBuffer buffer;
      CORE::COMM::ParObject::AddtoPack<double, 3>(buffer, inputData.bottom_corner_point_);
      CORE::COMM::ParObject::AddtoPack<double, 3>(buffer, inputData.top_corner_point_);
      CORE::COMM::ParObject::AddtoPack<int, 3>(buffer, inputData.interval_);
      CORE::COMM::ParObject::AddtoPack<double, 3>(buffer, inputData.rotation_angle_);
      CORE::COMM::ParObject::AddtoPack(buffer, static_cast<int>(inputData.autopartition_));
      CORE::COMM::ParObject::AddtoPack(buffer, inputData.elementtype_);
      CORE::COMM::ParObject::AddtoPack(buffer, inputData.distype_);
      CORE::COMM::ParObject::AddtoPack(buffer, inputData.elearguments_);
      buffer.StartPacking();
      CORE::COMM::ParObject::AddtoPack<double, 3>(buffer, inputData.bottom_corner_point_);
      CORE::COMM::ParObject::AddtoPack<double, 3>(buffer, inputData.top_corner_point_);
      CORE::COMM::ParObject::AddtoPack<int, 3>(buffer, inputData.interval_);
      CORE::COMM::ParObject::AddtoPack<double, 3>(buffer, inputData.rotation_angle_);
      CORE::COMM::ParObject::AddtoPack(buffer, static_cast<int>(inputData.autopartition_));
      CORE::COMM::ParObject::AddtoPack(buffer, inputData.elementtype_);
      CORE::COMM::ParObject::AddtoPack(buffer, inputData.distype_);
      CORE::COMM::ParObject::AddtoPack(buffer, inputData.elearguments_);
      std::swap(data, buffer());
    }

    ssize_t data_size = data.size();
    comm->Broadcast(&data_size, 1, 0);
    if (myrank != 0) data.resize(data_size, 0);
    comm->Broadcast(data.data(), data.size(), 0);

    if (myrank != 0)
    {
      size_t pos = 0;
      CORE::COMM::ParObject::ExtractfromPack<double, 3>(pos, data, inputData.bottom_corner_point_);
      CORE::COMM::ParObject::ExtractfromPack<double, 3>(pos, data, inputData.top_corner_point_);
      CORE::COMM::ParObject::ExtractfromPack<int, 3>(pos, data, inputData.interval_);
      CORE::COMM::ParObject::ExtractfromPack<double, 3>(pos, data, inputData.rotation_angle_);
      int autopartitionInteger;
      CORE::COMM::ParObject::ExtractfromPack(pos, data, autopartitionInteger);
      inputData.autopartition_ = autopartitionInteger;
      CORE::COMM::ParObject::ExtractfromPack(pos, data, inputData.elementtype_);
      CORE::COMM::ParObject::ExtractfromPack(pos, data, inputData.distype_);
      CORE::COMM::ParObject::ExtractfromPack(pos, data, inputData.elearguments_);
    }
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void DomainReader::complete() const
  {
    const int myrank = comm_->MyPID();

    Teuchos::Time time("", true);

    if (!myrank && !reader_.MyOutputFlag())
      CORE::IO::cout << "Complete discretization " << std::left << std::setw(16) << name_
                     << " in...." << CORE::IO::flush;

    int err = dis_->fill_complete(false, false, false);
    if (err) FOUR_C_THROW("dis_->fill_complete() returned %d", err);

    if (!myrank && !reader_.MyOutputFlag())
      CORE::IO::cout << time.totalElapsedTime(true) << " secs" << CORE::IO::endl;

    CORE::REBALANCE::UTILS::print_parallel_distribution(*dis_);
  }

}  // namespace CORE::IO

FOUR_C_NAMESPACE_CLOSE
