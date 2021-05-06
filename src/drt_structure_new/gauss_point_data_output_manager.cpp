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
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"

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

void STR::MODELEVALUATOR::GaussPointDataOutputManager::DistributeQuantities()
{
  const Epetra_Comm& comm = *DRT::Problem::Instance()->GetNPGroup()->GlobalComm();
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
    // receive averything from all other procs
    for (int i = 1; i < comm.NumProc(); ++i)
    {
      std::unique_ptr<std::unordered_map<std::string, int>> received_quantities =
          ReceiveQuantitiesFromProc(i);

      MergeQuantities(*received_quantities);
    }
  }
  else
  {
    SendMyQuantitiesToProc(0);
  }

  // Broadcast merged quantities to every proc
  BroadcastMyQuantitites(comm);
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::SendMyQuantitiesToProc(int to_proc) const
{
  // Send num quantities
  int num_quantitites = quantities_.size();
  MPI_Send(&num_quantitites, 1, MPI_INT, to_proc, MPI_TAG, MPI_COMM_WORLD);

  if (num_quantitites == 0) return;

  std::vector<int> quantities_size(0);
  std::vector<char> quantities_names_raw(0);
  quantities_size.reserve(num_quantitites);

  PackMyQuantityNames(quantities_names_raw);
  PackMyQuantitySizes(quantities_size);
  int length = quantities_names_raw.size();

  // Send total length of quantities names
  MPI_Send(&length, 1, MPI_INT, to_proc, MPI_TAG, MPI_COMM_WORLD);

  // Send quantities names
  MPI_Send(&quantities_names_raw[0], length, MPI_CHAR, to_proc, MPI_TAG, MPI_COMM_WORLD);

  // send size of quantities
  MPI_Send(&quantities_size[0], num_quantitites, MPI_INT, to_proc, MPI_TAG, MPI_COMM_WORLD);
}

std::unique_ptr<std::unordered_map<std::string, int>>
STR::MODELEVALUATOR::GaussPointDataOutputManager::ReceiveQuantitiesFromProc(int from_proc) const
{
  auto quantities = std::unique_ptr<std::unordered_map<std::string, int>>(
      new std::unordered_map<std::string, int>());
  MPI_Status status;

  // Receive num quantities
  int num_quantitites = 0;
  MPI_Recv(&num_quantitites, 1, MPI_INT, from_proc, MPI_TAG, MPI_COMM_WORLD, &status);

  if (num_quantitites == 0) return quantities;

  // Receive total length of quantity names
  int length = 0;
  MPI_Recv(&length, 1, MPI_INT, from_proc, MPI_TAG, MPI_COMM_WORLD, &status);

  // Receive quantities names
  std::vector<char> received_names_raw(length);
  MPI_Recv(&received_names_raw[0], length, MPI_CHAR, from_proc, MPI_TAG, MPI_COMM_WORLD, &status);

  // Receive sizes of quantities
  std::vector<int> quantities_size(num_quantitites, -1);
  MPI_Recv(
      &quantities_size[0], num_quantitites, MPI_INT, from_proc, MPI_TAG, MPI_COMM_WORLD, &status);


  UnpackQuantities(received_names_raw, quantities_size, *quantities);

  return quantities;
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::BroadcastMyQuantitites(
    const Epetra_Comm& comm)
{
  // Receive num quantities
  int num_quantitites = quantities_.size();
  MPI_Bcast(&num_quantitites, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // setup character representation of items to send
  std::vector<int> quantities_size(0);
  quantities_size.reserve(num_quantitites);
  std::vector<char> quantities_names_raw(0);
  int length;
  if (comm.MyPID() == 0)
  {
    PackMyQuantityNames(quantities_names_raw);
    PackMyQuantitySizes(quantities_size);
    length = quantities_names_raw.size();
  }
  else
  {
    quantities_size.resize(num_quantitites);
  }
  MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (comm.MyPID() != 0)
  {
    quantities_names_raw.resize(length);
  }

  MPI_Bcast(&quantities_names_raw[0], length, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&quantities_size[0], num_quantitites, MPI_INT, 0, MPI_COMM_WORLD);

  if (comm.MyPID() != 0)
  {
    std::unordered_map<std::string, int> received_quantities{};
    UnpackQuantities(quantities_names_raw, quantities_size, received_quantities);

    MergeQuantities(received_quantities);
  }

  // Make check, whether number of the quantities are as expected
  if (num_quantitites != static_cast<int>(quantities_.size()))
  {
    dserror("The number of quantities do not match after communication (%d vs. %d)",
        num_quantitites, quantities_.size());
  }
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::PackMyQuantityNames(
    std::vector<char>& names) const
{
  bool is_first = true;
  for (const auto& name_and_size : quantities_)
  {
    if (!is_first)
    {
      char delimiter = MPI_DELIMITER;
      names.emplace_back(delimiter);
    }
    std::copy(name_and_size.first.begin(), name_and_size.first.end(), std::back_inserter(names));
    is_first = false;
  }
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::PackMyQuantitySizes(
    std::vector<int>& sizes) const
{
  for (const auto& name_and_size : quantities_)
  {
    sizes.emplace_back(name_and_size.second);
  }
}

void STR::MODELEVALUATOR::GaussPointDataOutputManager::UnpackQuantities(
    std::vector<char>& names_raw, std::vector<int>& sizes,
    std::unordered_map<std::string, int>& quantities) const
{
  std::size_t i = 0;
  size_t pos = 0;
  std::string received_names(names_raw.begin(), names_raw.end());
  while ((pos = received_names.find(MPI_DELIMITER)) != std::string::npos)
  {
    std::string name = received_names.substr(0, pos);

    quantities.insert({name, sizes[i]});

    received_names.erase(0, pos + 1);
    ++i;
  }
}