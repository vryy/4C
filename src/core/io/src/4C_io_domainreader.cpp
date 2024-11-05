// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_io_domainreader.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element_definition.hpp"
#include "4C_io_gridgenerator.hpp"
#include "4C_io_input_file.hpp"
#include "4C_io_pstream.hpp"
#include "4C_rebalance_print.hpp"

#include <Teuchos_Time.hpp>

#include <algorithm>

FOUR_C_NAMESPACE_OPEN

namespace
{
  void broadcast_input_data_to_all_procs(
      const Epetra_Comm& comm, Core::IO::GridGenerator::RectangularCuboidInputs& inputData)
  {
    const int myrank = comm.MyPID();

    std::vector<char> data;
    if (myrank == 0)
    {
      Core::Communication::PackBuffer buffer;
      add_to_pack(buffer, inputData.bottom_corner_point_);
      add_to_pack(buffer, inputData.top_corner_point_);
      add_to_pack(buffer, inputData.interval_);
      add_to_pack(buffer, inputData.rotation_angle_);
      add_to_pack(buffer, static_cast<int>(inputData.autopartition_));
      add_to_pack(buffer, inputData.elementtype_);
      add_to_pack(buffer, inputData.distype_);
      add_to_pack(buffer, inputData.elearguments_);
      std::swap(data, buffer());
    }

    ssize_t data_size = data.size();
    comm.Broadcast(&data_size, 1, 0);
    if (myrank != 0) data.resize(data_size, 0);
    comm.Broadcast(data.data(), data.size(), 0);

    Core::Communication::UnpackBuffer buffer(data);
    if (myrank != 0)
    {
      extract_from_pack(buffer, inputData.bottom_corner_point_);
      extract_from_pack(buffer, inputData.top_corner_point_);
      extract_from_pack(buffer, inputData.interval_);
      extract_from_pack(buffer, inputData.rotation_angle_);
      int autopartitionInteger;
      extract_from_pack(buffer, autopartitionInteger);
      inputData.autopartition_ = autopartitionInteger;
      extract_from_pack(buffer, inputData.elementtype_);
      extract_from_pack(buffer, inputData.distype_);
      extract_from_pack(buffer, inputData.elearguments_);
    }
  }
}  // namespace

namespace Core::IO
{

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  DomainReader::DomainReader(Teuchos::RCP<Core::FE::Discretization> dis, Core::IO::InputFile& input,
      std::string sectionname)
      : name_(dis->name()),
        input_(input),
        comm_(dis->get_comm()),
        sectionname_(sectionname),
        dis_(dis)
  {
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void DomainReader::create_partitioned_mesh(int nodeGIdOfFirstNewNode) const
  {
    const int myrank = comm_.MyPID();

    Teuchos::Time time("", true);

    if (!input_.my_output_flag() && myrank == 0)
      Core::IO::cout << "Entering domain generation mode for " << name_
                     << " discretization ...\nCreate and partition elements      in...."
                     << Core::IO::endl;

    GridGenerator::RectangularCuboidInputs inputData =
        DomainReader::read_rectangular_cuboid_input_data();
    inputData.node_gid_of_first_new_node_ = nodeGIdOfFirstNewNode;

    Core::IO::GridGenerator::create_rectangular_cuboid_discretization(
        *dis_, inputData, static_cast<bool>(input_.my_output_flag()));

    if (!myrank && input_.my_output_flag() == 0)
      Core::IO::cout << "............................................... " << std::setw(10)
                     << std::setprecision(5) << std::scientific << time.totalElapsedTime(true)
                     << " secs" << Core::IO::endl;

    return;
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  GridGenerator::RectangularCuboidInputs DomainReader::read_rectangular_cuboid_input_data() const
  {
    Core::IO::GridGenerator::RectangularCuboidInputs inputData;
    // all reading is done on proc 0
    if (comm_.MyPID() == 0)
    {
      bool any_lines_read = false;
      // read domain info
      for (const auto& line : input_.lines_in_section(sectionname_))
      {
        any_lines_read = true;
        std::istringstream t{std::string{line}};
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
          std::transform(
              tmp.begin(), tmp.end(), tmp.begin(), [](unsigned char c) { return std::tolower(c); });
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

      if (!any_lines_read)
      {
        FOUR_C_THROW("No DOMAIN specified but box geometry selected!");
      }
    }

    // broadcast if necessary
    if (comm_.NumProc() > 1)
    {
      broadcast_input_data_to_all_procs(comm_, inputData);
    }

    return inputData;
  }


  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void DomainReader::complete() const
  {
    const int myrank = comm_.MyPID();

    Teuchos::Time time("", true);

    if (!myrank && !input_.my_output_flag())
      Core::IO::cout << "Complete discretization " << std::left << std::setw(16) << name_
                     << " in...." << Core::IO::flush;

    int err = dis_->fill_complete(false, false, false);
    if (err) FOUR_C_THROW("dis_->fill_complete() returned %d", err);

    if (!myrank && !input_.my_output_flag())
      Core::IO::cout << time.totalElapsedTime(true) << " secs" << Core::IO::endl;

    Core::Rebalance::Utils::print_parallel_distribution(*dis_);
  }

}  // namespace Core::IO

FOUR_C_NAMESPACE_CLOSE
