/*----------------------------------------------------------------------*/
/*! \file
\brief Computes a deal.II Triangulation from DRT::Discretization

\level 2

\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
*----------------------------------------------------------------------*/

#include "deal_tria_wrapper.H"

#ifdef HAVE_DEAL_II

#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>

#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_exporter.H"


namespace dealii
{
  template <int dim>
  DistributedTriangulation<dim>::DistributedTriangulation(
      Teuchos::RCP<DRT::DiscretizationHDG> discret)
      : discret(discret)
  {
    // 1. copy everything to processor zero:

    // the element node ids
    Teuchos::RCP<Epetra_Map> proc0elemap = LINALG::AllreduceEMap(*discret->ElementRowMap(), 0);
    Epetra_Import proc0eledataimporter(*proc0elemap, *discret->ElementRowMap());
    Teuchos::RCP<Epetra_MultiVector> elenodeids = Teuchos::rcp(
        new Epetra_MultiVector(*discret->ElementRowMap(), GeometryInfo<dim>::vertices_per_cell));
    for (int i = 0; i < discret->NumMyRowElements(); ++i)
      for (unsigned int j = 0; j < GeometryInfo<dim>::vertices_per_cell; ++j)
        elenodeids->operator()(j)->operator[](i) = discret->lRowElement(i)->NodeIds()[j];
    Teuchos::RCP<Epetra_MultiVector> proc0elenodeids =
        Teuchos::rcp(new Epetra_MultiVector(*proc0elemap, GeometryInfo<dim>::vertices_per_cell));
    int err = proc0elenodeids->Import(*elenodeids, proc0eledataimporter, Insert);
    if (err) dserror("Importing everything to proc 0 went wrong. Import returns %d", err);

    // and the node coordinates
    Teuchos::RCP<Epetra_Map> proc0nodemap = LINALG::AllreduceEMap(*discret->NodeRowMap(), 0);
    Epetra_Import proc0nodedataimporter(*proc0nodemap, *discret->NodeRowMap());
    Teuchos::RCP<Epetra_MultiVector> nodecoords =
        Teuchos::rcp(new Epetra_MultiVector(*discret->NodeRowMap(), dim));
    for (int i = 0; i < discret->NumMyRowNodes(); ++i)
      for (unsigned int j = 0; j < dim; ++j)
        nodecoords->operator()(j)->operator[](i) = discret->lRowNode(i)->X()[j];
    Teuchos::RCP<Epetra_MultiVector> proc0nodecoords =
        Teuchos::rcp(new Epetra_MultiVector(*proc0nodemap, dim));
    err = proc0nodecoords->Import(*nodecoords, proc0nodedataimporter, Insert);
    if (err) dserror("Importing everything to proc 0 went wrong. Import returns %d", err);

    std::vector<CellData<dim>> cells;
    SubCellData subcelldata;
    std::vector<Point<dim>> vertices(discret->NumGlobalNodes());

    int minallelegid = discret->ElementRowMap()->MinAllGID();
    int minallnodegid = discret->NodeRowMap()->MinAllGID();
    if (discret->Comm().MyPID() == 0)
    {
      // create the cells
      cells.resize(discret->NumGlobalElements());
      for (int cell = 0; cell < discret->NumGlobalElements(); ++cell)
      {
        for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        {
          cells[proc0elemap->GID(cell) - minallelegid].vertices[i] =
              proc0elenodeids->operator()(i)->operator[](cell) - minallnodegid;
        }
      }
      // create the vertices
      for (int node = 0; node < discret->NumGlobalNodes(); ++node)
      {
        for (unsigned int d = 0; d < dim; ++d)
        {
          vertices[proc0nodemap->GID(node) - minallnodegid][d] =
              proc0nodecoords->operator()(d)->operator[](node);
        }
      }
      Assert(subcelldata.check_consistency(dim), ExcInternalError());
    }
    // communicate the vertices
    for (int node = 0; node < discret->NumGlobalNodes(); ++node)
      discret->Comm().Broadcast(&vertices[node][0], dim, 0);

    // 2. invert_all_cells_of_negative_grid
    GridReordering<dim>::invert_all_cells_of_negative_grid(vertices, cells);

    // 3. reorder_cells
    if (discret->Comm().MyPID() ==
        0)  // need to enter this with only one processor, otherwise debug mode throws error
      GridReordering<dim>::reorder_cells(cells);

    // 4. send the col elements to the owner
    // the cells write back into the baci vector which is exported and then each processor can fill
    // its list of col elements
    if (discret->Comm().MyPID() == 0)
      for (int cell = 0; cell < discret->NumGlobalElements(); ++cell)
        for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
          proc0elenodeids->operator()(i)->operator[](cell) =
              cells[proc0elemap->GID(cell) - minallelegid].vertices[i];

    Teuchos::RCP<Epetra_MultiVector> elecolnodeids = Teuchos::rcp(
        new Epetra_MultiVector(*discret->ElementColMap(), GeometryInfo<dim>::vertices_per_cell));
    LINALG::Export(*proc0elenodeids, *elecolnodeids);

    std::vector<CellData<dim>> col_cells(discret->NumMyColElements());
    for (int cell = 0; cell < discret->NumMyColElements(); ++cell)
      for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
        col_cells[cell].vertices[i] = elecolnodeids->operator()(i)->operator[](cell);

    // 5. delete unused vertices (necessary since all vertices are stored on every processor) but
    // only col elements
    GridTools::delete_unused_vertices(vertices, col_cells, subcelldata);

    // 6. eventually create the triangulation
    this->create_triangulation_compatibility(vertices, col_cells, subcelldata);

    // 7. set the element owner
    for (int cell = 0; cell < discret->NumMyColElements(); ++cell)
    {
      cell_iterator dcell(this, 0, cell);
      dcell->set_subdomain_id(discret->lColElement(cell)->Owner());
    }
  }

  // This overrides the respective method in deal.II because we need to access
  // otherwise inaccessible data.
  template <int dim, int spacedim>
  __attribute__((always_inline)) inline void DoFHandler<dim, spacedim>::distribute_dofs(
      const FiniteElement<dim, spacedim> &ff)
  {
    std::cout << "Now I'm there!" << std::endl;
    selected_fe = &ff;
    clear_space();

    // this is internal::DoFHandler::Implementation::reserve_space (*this);
    vertex_dofs.resize(tria->n_vertices() * selected_fe->dofs_per_vertex,
        DoFHandler<dim, spacedim>::invalid_dof_index);

    for (unsigned int i = 0; i < tria->n_levels(); ++i)
    {
      levels.push_back(new internal::DoFHandler::DoFLevel<dim>);

      levels.back()->dof_object.dofs.resize(
          tria->n_raw_cells(i) * dim == 2
              ? selected_fe->dofs_per_quad
              : (dim == 3 ? selected_fe->dofs_per_hex : selected_fe->dofs_per_line),
          DoFHandler<dim, spacedim>::invalid_dof_index);

      levels.back()->cell_dof_indices_cache.resize(
          tria->n_raw_cells(i) * selected_fe->dofs_per_cell,
          DoFHandler<dim, spacedim>::invalid_dof_index);
    }

    faces = new internal::DoFHandler::DoFFaces<dim>;
    // avoid access to n_raw_lines when there are no cells
    if (tria->n_cells() > 0)
    {
      if (dim > 1)
        faces->lines.dofs.resize(tria->n_raw_lines() * selected_fe->dofs_per_line,
            DoFHandler<dim, spacedim>::invalid_dof_index);
      if (dim == 3)
        reinterpret_cast<dealii::internal::DoFHandler::DoFFaces<3> *>(faces)->quads.dofs.resize(
            tria->n_raw_quads() * selected_fe->dofs_per_quad,
            DoFHandler<dim, spacedim>::invalid_dof_index);
    }

    const DistributedTriangulation<dim> *tria =
        dynamic_cast<const DistributedTriangulation<dim> *>(&get_triangulation());
    if (tria == 0) dserror("Only works for baci wrapper into DistributedTriangulation<dim>!");

    // then go through the cells again an manually set the degrees of freedom
    // deal.II is flexible enough to allow doing this... :-)
    for (typename DoFHandler<dim>::active_cell_iterator cell = begin_active(); cell != end();
         ++cell)
    {
      const types::global_dof_index index_base =
          selected_fe->dofs_per_cell *
          types::global_dof_index(tria->discret->lColElement(cell->index())->Id());
      for (unsigned int d = 0; d < selected_fe->dofs_per_cell; ++d)
        dealii::internal::DoFAccessor::Implementation::set_dof_index(*this, cell->level(),
            cell->index(), 0, d, dealii::internal::int2type<dim>(), index_base + d);
    }

    // finally update the cache for dof indices
    for (typename DoFHandler<dim>::active_cell_iterator cell = begin_active(); cell != end();
         ++cell)
      cell->update_cell_dof_indices_cache();

    const Epetra_Map *rowmap = tria->discret->ElementRowMap();
    number_cache.n_global_dofs =
        types::global_dof_index(rowmap->NumGlobalElements()) * selected_fe->dofs_per_cell;
    number_cache.n_locally_owned_dofs =
        types::global_dof_index(rowmap->NumMyElements()) * selected_fe->dofs_per_cell;

    if (!rowmap->LinearMap()) dserror("Only linear maps are supported");

    number_cache.locally_owned_dofs = IndexSet(number_cache.n_global_dofs);
    number_cache.locally_owned_dofs.add_range(
        types::global_dof_index(rowmap->MinMyGID()) * selected_fe->dofs_per_cell,
        types::global_dof_index(rowmap->MaxMyGID()) * selected_fe->dofs_per_cell + 1);
    number_cache.locally_owned_dofs.compress();

    // TODO: fill in correct info here!
    number_cache.n_locally_owned_dofs_per_processor =
        std::vector<types::global_dof_index>(rowmap->Comm().NumProc(), -1);

    number_cache.locally_owned_dofs_per_processor =
        std::vector<IndexSet>(rowmap->Comm().NumProc(), IndexSet(number_cache.n_global_dofs));
  }

  template <int dim>
  void assign_dg_dofs(const FiniteElement<dim> &fe, DoFHandler<dim> &dof_handler)
  {
    const DistributedTriangulation<dim> *tria =
        dynamic_cast<const DistributedTriangulation<dim> *>(&dof_handler.get_triangulation());
    if (tria == 0) dserror("Only works for baci wrapper into DistributedTriangulation<dim>!");

    if ((dim == 2 && fe.dofs_per_quad != fe.dofs_per_cell) ||
        (dim == 3 && fe.dofs_per_hex != fe.dofs_per_cell))
      dserror("This code only works for DG elements with all DoFs sitting inside the element.");

    // first let deal.II distribute the degrees of freedom
    dof_handler.distribute_dofs(fe);

    // then go through the cells again and manually set the degrees of freedom
    // deal.II is flexible enough to allow doing this... :-)
    //
    // for the degrees of freedom, we want to assign contiguous ranges of numbers
    // even though the elements in DRT::Discretization might have arbitrary numberings
    // in the element ids. Therefore, we create an index vector of 'renumbered'
    // element ids.
    //
    // step 1: count the number of row elements in each processor and let all processors know
    const Epetra_Map *rowmap = tria->discret->ElementRowMap();
    int nids = rowmap->NumMyElements();
    std::vector<int> ids_per_proc(rowmap->Comm().NumProc());
    rowmap->Comm().GatherAll(&nids, &ids_per_proc[0], 1);
    int shift = 0;
    for (int i = 0; i < rowmap->Comm().MyPID(); ++i) shift += ids_per_proc[i];

    // step 2: assign contiguous number to row elements of current processor. Since only
    // double vectors support communication, fill the data there. Luckily, 64-bit doubles
    // can exactly hold all 32-bit integers. These numbers will then be set to the deal.II
    // elements
    Epetra_Vector eleids(*rowmap);
    for (int i = 0; i < rowmap->NumMyElements(); ++i) eleids[i] = shift + i;
    Epetra_Vector ghosted(*tria->discret->ElementColMap());
    LINALG::Export(eleids, ghosted);

    for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
         cell != dof_handler.end(); ++cell)
    {
      const types::global_dof_index index_base =
          fe.dofs_per_cell *
          types::global_dof_index(
              ghosted[ghosted.Map().LID(tria->discret->lColElement(cell->index())->Id())]);
      for (unsigned int d = 0; d < fe.dofs_per_cell; ++d)
        dealii::internal::DoFAccessor::Implementation::set_dof_index(dof_handler, cell->level(),
            cell->index(), 0, d, dealii::internal::int2type<dim>(), index_base + d);
    }

    // finally update the cache for dof indices
    for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
         cell != dof_handler.end(); ++cell)
      cell->update_cell_dof_indices_cache();

    dof_handler.number_cache.n_global_dofs =
        types::global_dof_index(rowmap->NumGlobalElements()) * fe.dofs_per_cell;
    dof_handler.number_cache.n_locally_owned_dofs =
        types::global_dof_index(rowmap->NumMyElements()) * fe.dofs_per_cell;

    dof_handler.number_cache.locally_owned_dofs = IndexSet(dof_handler.number_cache.n_global_dofs);
    dof_handler.number_cache.locally_owned_dofs.add_range(
        types::global_dof_index(shift) * fe.dofs_per_cell,
        types::global_dof_index(shift + nids) * fe.dofs_per_cell);
    dof_handler.number_cache.locally_owned_dofs.compress();

    dof_handler.number_cache.n_locally_owned_dofs_per_processor =
        std::vector<types::global_dof_index>(rowmap->Comm().NumProc(), -1);
    for (int i = 0; i < rowmap->Comm().NumProc(); ++i)
      dof_handler.number_cache.n_locally_owned_dofs_per_processor[i] =
          types::global_dof_index(ids_per_proc[i]) * fe.dofs_per_cell;

    dof_handler.number_cache.locally_owned_dofs_per_processor = std::vector<IndexSet>(
        rowmap->Comm().NumProc(), IndexSet(dof_handler.number_cache.n_global_dofs));
    types::global_dof_index last_size = 0;
    for (int i = 0; i < rowmap->Comm().NumProc(); ++i)
    {
      dof_handler.number_cache.locally_owned_dofs_per_processor[i].add_range(
          last_size, last_size + dof_handler.number_cache.n_locally_owned_dofs_per_processor[i]);
      dof_handler.number_cache.locally_owned_dofs_per_processor[i].compress();
      last_size = last_size + dof_handler.number_cache.n_locally_owned_dofs_per_processor[i];
    }
    dsassert(last_size == types::global_dof_index(rowmap->NumGlobalElements()) * fe.dofs_per_cell,
        "Unexpected result in accumulation");
  }


  template class DistributedTriangulation<2>;
  template class DistributedTriangulation<3>;
  template void assign_dg_dofs(const FiniteElement<2> &, DoFHandler<2> &);
  template void assign_dg_dofs(const FiniteElement<3> &, DoFHandler<3> &);
}  // namespace dealii

#endif
