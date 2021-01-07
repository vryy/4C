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
#include <Teuchos_RCPDecl.hpp>
#include "str_model_evaluator.H"

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
      const int size = name_and_size.second;

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