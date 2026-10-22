// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_graph.hpp"

#include "4C_linalg_utils_exceptions.hpp"

#include <algorithm>

FOUR_C_NAMESPACE_OPEN

Core::LinAlg::Graph::Graph(const Epetra_CrsGraph& Source)
    : graphtype_(CRS_GRAPH), graph_(std::make_unique<Epetra_CrsGraph>(Source))
{
}

Core::LinAlg::Graph::Graph(const Epetra_FECrsGraph& Source)
    : graphtype_(CRS_GRAPH), graph_(std::make_unique<Epetra_CrsGraph>(Source))
{
}

Core::LinAlg::Graph::Graph(const Map& RowMap, const int* NumIndicesPerRow, GraphType graphtype)
    : graphtype_(graphtype)
{
  if (graphtype_ == CRS_GRAPH)
  {
    graph_ =
        std::make_unique<Epetra_CrsGraph>(::Copy, RowMap.get_epetra_block_map(), NumIndicesPerRow);
  }
  else if (graphtype_ == FE_GRAPH)
  {
    graph_ = std::make_unique<Epetra_FECrsGraph>(
        ::Copy, RowMap.get_epetra_block_map(), const_cast<int*>(NumIndicesPerRow));
  }
  else
  {
    FOUR_C_THROW("Construction of a graph with an unknown graph type.");
  }
}

Core::LinAlg::Graph::Graph(const Map& RowMap, int NumIndicesPerRow, GraphType graphtype)
    : graphtype_(graphtype)
{
  if (graphtype_ == CRS_GRAPH)
  {
    graph_ =
        std::make_unique<Epetra_CrsGraph>(::Copy, RowMap.get_epetra_block_map(), NumIndicesPerRow);
  }
  else if (graphtype_ == FE_GRAPH)
  {
    graph_ = std::make_unique<Epetra_FECrsGraph>(
        ::Copy, RowMap.get_epetra_block_map(), NumIndicesPerRow);
  }
  else
  {
    FOUR_C_THROW("Construction of a graph with an unknown graph type.");
  }
}

Core::LinAlg::Graph::Graph(const Graph& other)
    : graphtype_(CRS_GRAPH), graph_(std::make_unique<Epetra_CrsGraph>(other.get_epetra_crs_graph()))
{
}

Core::LinAlg::Graph& Core::LinAlg::Graph::operator=(const Graph& other)
{
  *graph_ = other.get_epetra_crs_graph();
  return *this;
}

void Core::LinAlg::Graph::fill_complete()
{
  if (graphtype_ == CRS_GRAPH)
  {
    ASSERT_EPETRA_CALL(graph_->FillComplete());
  }
  else if (graphtype_ == FE_GRAPH)
  {
    ASSERT_EPETRA_CALL(static_cast<Epetra_FECrsGraph*>(graph_.get())->GlobalAssemble());
  }
}

void Core::LinAlg::Graph::fill_complete(const Map& domain_map, const Map& range_map)
{
  if (graphtype_ == CRS_GRAPH)
  {
    ASSERT_EPETRA_CALL(
        graph_->FillComplete(domain_map.get_epetra_block_map(), range_map.get_epetra_block_map()));
  }
  else if (graphtype_ == FE_GRAPH)
  {
    ASSERT_EPETRA_CALL(static_cast<Epetra_FECrsGraph*>(graph_.get())
            ->GlobalAssemble(domain_map.get_epetra_map(), range_map.get_epetra_map()));
  }
}

void Core::LinAlg::Graph::optimize_storage() { ASSERT_EPETRA_CALL(graph_->OptimizeStorage()); }

void Core::LinAlg::Graph::export_to(const Core::LinAlg::Graph& A,
    const Core::LinAlg::Export& Exporter, Core::LinAlg::CombineMode CombineMode)
{
  ASSERT_EPETRA_CALL(graph_->Export(A.get_epetra_crs_graph(), Exporter.get_epetra_export(),
      Core::LinAlg::to_epetra_combine_mode(CombineMode)));
}

void Core::LinAlg::Graph::import_from(const Core::LinAlg::Graph& A,
    const Core::LinAlg::Import& Importer, Core::LinAlg::CombineMode CombineMode)
{
  ASSERT_EPETRA_CALL(graph_->Import(A.get_epetra_crs_graph(), Importer.get_epetra_import(),
      Core::LinAlg::to_epetra_combine_mode(CombineMode)));
}

void Core::LinAlg::Graph::insert_global_indices(int GlobalRow, std::span<int>& Indices)
{
  if (graphtype_ == CRS_GRAPH)
  {
    ASSERT_EPETRA_CALL(graph_->InsertGlobalIndices(GlobalRow, Indices.size(), Indices.data()));
  }
  else if (graphtype_ == FE_GRAPH)
  {
    ASSERT_EPETRA_CALL(static_cast<Epetra_FECrsGraph*>(graph_.get())
            ->InsertGlobalIndices(1, &GlobalRow, Indices.size(), Indices.data()));
  }
}

void Core::LinAlg::Graph::extract_local_row_view(int LocalRow, std::span<int>& Indices) const
{
  int num_indices;
  int* indices;
  ASSERT_EPETRA_CALL(graph_->ExtractMyRowView(LocalRow, num_indices, indices));
  Indices = std::span(indices, num_indices);
}

void Core::LinAlg::Graph::extract_global_row_view(int GlobalRow, std::span<int>& Indices) const
{
  int num_indices;
  int* indices;
  ASSERT_EPETRA_CALL(graph_->ExtractGlobalRowView(GlobalRow, num_indices, indices));
  Indices = std::span(indices, num_indices);
}

bool Core::LinAlg::Graph::is_identical_to(const Graph& other) const
{
  FOUR_C_ASSERT(filled() && other.filled(), "Only filled graphs can be compared.");

  // Comparing maps can be collective, so none of the comparisons is skipped.
  const bool same_row_map = graph_->RowMap().SameAs(other.graph_->RowMap());
  const bool same_col_map = graph_->ColMap().SameAs(other.graph_->ColMap());
  const bool same_domain_map = graph_->DomainMap().SameAs(other.graph_->DomainMap());
  const bool same_range_map = graph_->RangeMap().SameAs(other.graph_->RangeMap());
  if (!(same_row_map && same_col_map && same_domain_map && same_range_map)) return false;

  const int num_rows = num_local_rows();
  const int num_nonzeros = num_local_nonzeros();
  if (num_rows != other.num_local_rows() || num_nonzeros != other.num_local_nonzeros())
    return false;

  if (graph_->StorageOptimized() && other.graph_->StorageOptimized())
  {
    // The expert accessors are not const but only read here.
    const int* offsets = graph_->ExpertExtractIndexOffset().Values();
    const int* other_offsets = other.graph_->ExpertExtractIndexOffset().Values();
    const int* indices = graph_->ExpertExtractIndices().Values();
    const int* other_indices = other.graph_->ExpertExtractIndices().Values();

    return std::equal(offsets, offsets + num_rows + 1, other_offsets) &&
           std::equal(indices, indices + num_nonzeros, other_indices);
  }

  for (int row = 0; row < num_rows; ++row)
  {
    std::span<int> indices;
    std::span<int> other_indices;
    extract_local_row_view(row, indices);
    other.extract_local_row_view(row, other_indices);

    if (!std::ranges::equal(indices, other_indices)) return false;
  }
  return true;
}

FOUR_C_NAMESPACE_CLOSE
