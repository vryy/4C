/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of a Gauss point data output handler

\level 3
*/
/*----------------------------------------------------------------------*/
#include "gauss_point_data_output_manager.H"
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_IntVector.h>
#include <Epetra_Map.h>
#include <mpi.h>
#include <Teuchos_RCPDecl.hpp>
#include "str_model_evaluator.H"
#include "drt_globalproblem.H"
#include "comm_utils.H"
#include "drt_exporter.H"

STR::MODELEVALUATOR::GaussPointDataOutputManager::GaussPointDataOutputManager(
    INPAR::STR::GaussPointDataOutputType output_type)
    : output_type_(output_type),
      max_num_gp_(0),
      data_nodes_({}),
      data_nodes_count_({}),
      data_element_center_({}),
      data_gauss_point_({}),
      quantities_({})
{
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::AddQuantityIfNotExistant(
    const std::string& name, int size)
{
  const auto item = quantities_.find(name);
  if (item != quantities_.end())
  {
    if (item->second != size)
    {
      dserror(
          "The quantity %s is already registered, but with a different size (%d vs. %d). This is "
          "fatal!",
          name.c_str(), size, item->second);
    }
  }
  else
  {
    if (name.find(MPI_DELIMITER) != std::string::npos)
    {
      dserror(
          "The quantity name %s for Gauss Point VTK runtime output contains the delimiter %s that "
          "is used for MPI communication. This is not allowed.",
          name.c_str(), MPI_DELIMITER);
    }
    quantities_.insert({name, size});
  }
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::MergeQuantities(
    const std::unordered_map<std::string, int>& quantities)
{
  for (const auto& name_and_size : quantities)
  {
    const std::string& name = name_and_size.first;
    const int size = name_and_size.second;

    AddQuantityIfNotExistant(name, size);
  }
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::AddElementNumberOfGaussPoints(
    const int numgp)
{
  if (numgp > max_num_gp_)
  {
    max_num_gp_ = numgp;
  }
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::PrepareData(
    const Epetra_Map& node_col_map, const Epetra_Map& element_row_map)
{
  switch (output_type_)
  {
    case INPAR::STR::GaussPointDataOutputType::nodes:
      PrepareNodalDataVectors(node_col_map);
      break;
    case INPAR::STR::GaussPointDataOutputType::element_center:
      PrepareElementCenterDataVectors(element_row_map);
      break;
    case INPAR::STR::GaussPointDataOutputType::gauss_points:
      PrepareGaussPointDataVectors(element_row_map);
      break;
    case INPAR::STR::GaussPointDataOutputType::none:
      dserror("Your Gauss point data output type is none, so you don't need to prepare data!");
    default:
      dserror("Unknown Gauss point data output type");
  }
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::PrepareNodalDataVectors(
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

void STR::MODELEVALUATOR::GaussPointDataOutputManager::PrepareElementCenterDataVectors(
    const Epetra_Map& element_col_map)
{
  for (const auto& name_and_size : quantities_)
  {
    const std::string& name = name_and_size.first;
    const int size = name_and_size.second;

    data_element_center_[name] = Teuchos::rcp(new Epetra_MultiVector(element_col_map, size, true));
  }
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::PrepareGaussPointDataVectors(
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

void STR::MODELEVALUATOR::GaussPointDataOutputManager::PostEvaluate()
{
  if (output_type_ == INPAR::STR::GaussPointDataOutputType::nodes)
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

void STR::MODELEVALUATOR::GaussPointDataOutputManager::DistributeQuantities(const Epetra_Comm& comm)
{
  const DRT::Exporter exporter(comm);

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
          ReceiveQuantitiesFromProc(exporter, i);

      MergeQuantities(*received_quantities);
    }
  }
  else
  {
    SendMyQuantitiesToProc(exporter, 0);
  }

  // Broadcast merged quantities to every proc
  BroadcastMyQuantitites(exporter);
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::SendMyQuantitiesToProc(
    const DRT::Exporter& exporter, int to_proc) const
{
  // Pack quantities
  std::vector<char> sdata(0);
  PackMyQuantities(sdata);

  MPI_Request request;
  exporter.ISend(exporter.Comm().MyPID(), 0, &sdata[0], sdata.size(), MPI_TAG, request);
  exporter.Wait(request);
}

std::unique_ptr<std::unordered_map<std::string, int>>
STR::MODELEVALUATOR::GaussPointDataOutputManager::ReceiveQuantitiesFromProc(
    const DRT::Exporter& exporter, int from_proc) const
{
  std::vector<char> rdata(0);
  int size;
  exporter.Receive(from_proc, MPI_TAG, rdata, size);

  auto quantities = std::unique_ptr<std::unordered_map<std::string, int>>(
      new std::unordered_map<std::string, int>());

  std::size_t pos = 0;
  UnpackQuantities(pos, rdata, *quantities);

  return quantities;
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::BroadcastMyQuantitites(
    const DRT::Exporter& exporter)
{
  std::vector<char> data(0);
  if (exporter.Comm().MyPID() == 0)
  {
    PackMyQuantities(data);
  }

  exporter.Broadcast(0, data, MPI_TAG);

  if (exporter.Comm().MyPID() != 0)
  {
    std::size_t pos = 0;
    std::unordered_map<std::string, int> received_quantities{};
    UnpackQuantities(pos, data, received_quantities);

    MergeQuantities(received_quantities);
  }
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::PackMyQuantities(
    std::vector<char>& data) const
{
  DRT::PackBuffer packBuffer;
  packBuffer.StartPacking();
  DRT::ParObject::AddtoPack(packBuffer, quantities_);
  std::swap(data, packBuffer());
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::UnpackQuantities(std::size_t pos,
    const std::vector<char>& data, std::unordered_map<std::string, int>& quantities) const
{
  DRT::ParObject::ExtractfromPack(pos, data, quantities);
}