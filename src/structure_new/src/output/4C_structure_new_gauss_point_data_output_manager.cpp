/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of a Gauss point data output handler

\level 3
*/
/*----------------------------------------------------------------------*/
#include "4C_structure_new_gauss_point_data_output_manager.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_comm_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_structure_new_model_evaluator.hpp"

#include <Epetra_IntVector.h>
#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCPDecl.hpp>

FOUR_C_NAMESPACE_OPEN

STR::MODELEVALUATOR::GaussPointDataOutputManager::GaussPointDataOutputManager(
    Inpar::STR::GaussPointDataOutputType output_type)
    : output_type_(output_type),
      max_num_gp_(0),
      data_nodes_({}),
      data_nodes_count_({}),
      data_element_center_({}),
      data_gauss_point_({}),
      quantities_({})
{
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::add_quantity_if_not_existant(
    const std::string& name, int size)
{
  const auto item = quantities_.find(name);
  if (item != quantities_.end())
  {
    if (item->second != size)
    {
      FOUR_C_THROW(
          "The quantity %s is already registered, but with a different size (%d vs. %d). This is "
          "fatal!",
          name.c_str(), size, item->second);
    }
  }
  else
  {
    if (name.find(MPI_DELIMITER) != std::string::npos)
    {
      FOUR_C_THROW(
          "The quantity name %s for Gauss Point VTK runtime output contains the delimiter %s that "
          "is used for MPI communication. This is not allowed.",
          name.c_str(), MPI_DELIMITER);
    }
    quantities_.insert({name, size});
  }
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::merge_quantities(
    const std::unordered_map<std::string, int>& quantities)
{
  for (const auto& name_and_size : quantities)
  {
    const std::string& name = name_and_size.first;
    const int size = name_and_size.second;

    add_quantity_if_not_existant(name, size);
  }
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::add_element_number_of_gauss_points(
    const int numgp)
{
  if (numgp > max_num_gp_)
  {
    max_num_gp_ = numgp;
  }
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::prepare_data(
    const Epetra_Map& node_col_map, const Epetra_Map& element_row_map)
{
  switch (output_type_)
  {
    case Inpar::STR::GaussPointDataOutputType::nodes:
      prepare_nodal_data_vectors(node_col_map);
      break;
    case Inpar::STR::GaussPointDataOutputType::element_center:
      prepare_element_center_data_vectors(element_row_map);
      break;
    case Inpar::STR::GaussPointDataOutputType::gauss_points:
      prepare_gauss_point_data_vectors(element_row_map);
      break;
    case Inpar::STR::GaussPointDataOutputType::none:
      FOUR_C_THROW("Your Gauss point data output type is none, so you don't need to prepare data!");
    default:
      FOUR_C_THROW("Unknown Gauss point data output type");
  }
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::prepare_nodal_data_vectors(
    const Epetra_Map& node_col_map)
{
  for (const auto& name_and_size : quantities_)
  {
    const std::string& name = name_and_size.first;
    const int size = name_and_size.second;

    data_nodes_[name] = Teuchos::rcp(new Epetra_MultiVector(node_col_map, size, true));
    data_nodes_count_[name] = Teuchos::rcp(new Epetra_IntVector(node_col_map, true));
  }
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::prepare_element_center_data_vectors(
    const Epetra_Map& element_col_map)
{
  for (const auto& name_and_size : quantities_)
  {
    const std::string& name = name_and_size.first;
    const int size = name_and_size.second;

    data_element_center_[name] = Teuchos::rcp(new Epetra_MultiVector(element_col_map, size, true));
  }
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::prepare_gauss_point_data_vectors(
    const Epetra_Map& element_col_map)
{
  for (const auto& name_and_size : quantities_)
  {
    const std::string& name = name_and_size.first;
    const int size = name_and_size.second;

    data_gauss_point_.emplace(name, std::vector<Teuchos::RCP<Epetra_MultiVector>>());
    std::vector<Teuchos::RCP<Epetra_MultiVector>>& gp_data = data_gauss_point_[name];

    gp_data.resize(max_num_gp_);
    for (auto& data_i : gp_data)
    {
      data_i = Teuchos::rcp(new Epetra_MultiVector(element_col_map, size, true));
    }
  }
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::post_evaluate()
{
  if (output_type_ == Inpar::STR::GaussPointDataOutputType::nodes)
  {
    // divide nodal quantities by the nodal count
    for (const auto& name_and_size : quantities_)
    {
      const std::string& name = name_and_size.first;

      Epetra_MultiVector& nodal_data = *data_nodes_[name];
      const Epetra_IntVector& nodal_count = *data_nodes_count_[name];

      for (int col = 0; col < nodal_data.NumVectors(); ++col)
      {
        Epetra_Vector& data_item = *nodal_data(col);

        for (int i = 0; i < data_item.MyLength(); ++i)
        {
          if (nodal_count[i] != 0)
          {
            data_item[i] /= nodal_count[i];
          }
        }
      }
    }
  }
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::distribute_quantities(
    const Epetra_Comm& comm)
{
  const Core::Communication::Exporter exporter(comm);

  int max_quantities = quantities_.size();
  comm.MaxAll(&max_quantities, &max_quantities, 1);

  if (max_quantities == 0)
  {
    // Nothing to distribute
    return;
  }

  // Communicate max number of Gauss points
  comm.MaxAll(&max_num_gp_, &max_num_gp_, 1);

  // Collect all quantities on proc 0
  if (comm.MyPID() == 0)
  {
    // receive everything from all other procs
    for (int i = 1; i < comm.NumProc(); ++i)
    {
      std::unique_ptr<std::unordered_map<std::string, int>> received_quantities =
          receive_quantities_from_proc(exporter, i);

      merge_quantities(*received_quantities);
    }
  }
  else
  {
    send_my_quantities_to_proc(exporter, 0);
  }

  // Broadcast merged quantities to every proc
  broadcast_my_quantitites(exporter);
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::send_my_quantities_to_proc(
    const Core::Communication::Exporter& exporter, int to_proc) const
{
  // Pack quantities
  std::vector<char> sdata(0);
  pack_my_quantities(sdata);

  MPI_Request request;
  exporter.i_send(exporter.Comm().MyPID(), 0, sdata.data(), sdata.size(), MPI_TAG, request);
  exporter.Wait(request);
}

std::unique_ptr<std::unordered_map<std::string, int>>
STR::MODELEVALUATOR::GaussPointDataOutputManager::receive_quantities_from_proc(
    const Core::Communication::Exporter& exporter, int from_proc) const
{
  std::vector<char> rdata(0);
  int size;
  exporter.Receive(from_proc, MPI_TAG, rdata, size);

  auto quantities = std::unique_ptr<std::unordered_map<std::string, int>>(
      new std::unordered_map<std::string, int>());

  std::size_t pos = 0;
  unpack_quantities(pos, rdata, *quantities);

  return quantities;
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::broadcast_my_quantitites(
    const Core::Communication::Exporter& exporter)
{
  std::vector<char> data(0);
  if (exporter.Comm().MyPID() == 0)
  {
    pack_my_quantities(data);
  }

  exporter.Broadcast(0, data, MPI_TAG);

  if (exporter.Comm().MyPID() != 0)
  {
    std::size_t pos = 0;
    std::unordered_map<std::string, int> received_quantities{};
    unpack_quantities(pos, data, received_quantities);

    merge_quantities(received_quantities);
  }
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::pack_my_quantities(
    std::vector<char>& data) const
{
  Core::Communication::PackBuffer packBuffer;
  Core::Communication::ParObject::add_to_pack(packBuffer, quantities_);
  std::swap(data, packBuffer());
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::unpack_quantities(std::size_t pos,
    const std::vector<char>& data, std::unordered_map<std::string, int>& quantities) const
{
  Core::Communication::ParObject::extract_from_pack(pos, data, quantities);
}

FOUR_C_NAMESPACE_CLOSE
