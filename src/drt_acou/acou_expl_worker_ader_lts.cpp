/*!----------------------------------------------------------------------
\file acou_expl_worker_ader_lts.cpp
\brief Control routine for acoustic explicit time integration with ADER LTS

<pre>
\level 3

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
*/
/*----------------------------------------------------------------------*/

#include "acou_expl_worker.H"

#ifdef HAVE_DEAL_II

#include <deal.II/lac/parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>
//#include <deal.II/matrix_free/fe_evaluation.h>
#include "fe_evaluation.h"
#include <deal.II/matrix_free/operators.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/base/timer.h>

#include <Epetra_MpiComm.h>

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/acoustic.H"
#include "../drt_mat/acoustic_sol.H"
#include "acou_ele.H"
#include "acou_sol_ele.H"


namespace ACOU
{
  template <typename Number>
  template <typename Operator>
  void ClusterManager<Number>::setup_mf_index_to_cell_index(const Operator &op)
  {
    mf_index.resize(op.get_matrix_free().get_dof_handler(0).get_triangulation().n_active_cells());

    for (unsigned int i = 0; i < n_cells_with_ghosts; ++i)
      for (unsigned int v = 0; v < op.get_matrix_free().n_components_filled(i); ++v)
        mf_index[op.get_matrix_free().get_cell_iterator(i, v)->active_cell_index()] = i;
  }


  template <typename Number>
  template <typename Operator>
  void ClusterManager<Number>::setup(const Operator &op, const unsigned int max_clusters,
      const unsigned int max_diff, Number cfl_number_in)
  {
    // the input value for dt is the smallest allowed time step, all other time
    // steps now must be set to multiples of this dt
    fastest_time_step = op.time_step;
    std::cout << cfl_number_in << std::endl;
    cfl_number = cfl_number_in;

    n_cells_with_ghosts =
        op.get_matrix_free().n_macro_cells() + op.get_matrix_free().n_macro_ghost_cells();
    n_cells = op.get_matrix_free().n_macro_cells();

    setup_mf_index_to_cell_index(op);

    const unsigned dofs_per_cell = op.get_matrix_free().get_dof_handler().get_fe().dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    {
      distr_cell_index.resize(n_cells_with_ghosts);
      for (unsigned int i = 0; i < n_cells_with_ghosts; ++i)
        distr_cell_index[i].resize(op.get_matrix_free().n_components_filled(i));
    }

    IndexSet elerowset(op.flux_memory[0].size() / dofs_per_cell);
    IndexSet elecolset(op.flux_memory[0].size() / dofs_per_cell);

    for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      for (unsigned int v = 0; v < op.get_matrix_free().n_components_filled(cell); ++v)
      {
        op.get_matrix_free().get_cell_iterator(cell, v)->get_dof_indices(local_dof_indices);
        elerowset.add_index(local_dof_indices[0] / dofs_per_cell);
        distr_cell_index[cell][v] = local_dof_indices[0] / dofs_per_cell;
      }
    }
    for (unsigned int cell = n_cells; cell < n_cells_with_ghosts; ++cell)
    {
      for (unsigned int v = 0; v < op.get_matrix_free().n_components_filled(cell); ++v)
      {
        op.get_matrix_free().get_cell_iterator(cell, v)->get_dof_indices(local_dof_indices);
        elecolset.add_index(local_dof_indices[0] / dofs_per_cell);
        distr_cell_index[cell][v] = local_dof_indices[0] / dofs_per_cell;
      }
    }

    elerowset.compress();
    elecolset.compress();

    cell_cluster_ids.resize(n_cells_with_ghosts);
    cell_have_faster_neighbor.resize(n_cells_with_ghosts);
    cell_have_slower_neighbor.resize(n_cells_with_ghosts);
    n_clusters = 0;

    unsigned int smallest_cluster_id = std::numeric_limits<unsigned int>::max();
    for (unsigned int i = 0; i < n_cells_with_ghosts; ++i)
    {
      // all cells vectorized together get the same time step (the smallest of all)
      double vectcelltimestep = std::numeric_limits<double>::max();
      for (unsigned int v = 0; v < op.get_matrix_free().n_components_filled(i); ++v)
      {
        // std::cout<<"p "<<Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)<<" i "<<i<<" v "<<v<<"
        // center "<<op.get_matrix_free().get_cell_iterator(i,v)->center()<<std::endl;
        double h = op.get_matrix_free().get_cell_iterator(i, v)->minimum_vertex_distance();
        double c = op.speed_of_sound(i, v);
        double acttimestep = cfl_number * h / c;

        if (acttimestep < vectcelltimestep) vectcelltimestep = acttimestep;
      }
      // now cast the optimal time step to a multiple of the smallest time step
      unsigned int cluster_id = int(vectcelltimestep / fastest_time_step);
      if (cluster_id > n_clusters) n_clusters = cluster_id;
      if (cluster_id < smallest_cluster_id) smallest_cluster_id = cluster_id;
      cell_cluster_ids[i] = cluster_id;
    }

    // communicate the number of clusters and the smallest cluster
    n_clusters = Utilities::MPI::max(n_clusters, MPI_COMM_WORLD);
    smallest_cluster_id = Utilities::MPI::min(smallest_cluster_id, MPI_COMM_WORLD);

    n_clusters = n_clusters - smallest_cluster_id + 1;
    // adapt for missmatch in cluster ids
    if (smallest_cluster_id != 0)
      for (unsigned int i = 0; i < n_cells; ++i)
        cell_cluster_ids[i] = cell_cluster_ids[i] - smallest_cluster_id;

    // more than allowed clusters for diff=1
    unsigned int cluster_diff = 1;
    if (n_clusters > max_clusters)
    {
      // the differences between neighboring clusters should be diff
      cluster_diff = n_clusters / max_clusters + 1;
      cluster_diff =
          std::min(cluster_diff, max_diff);  // but this difference is also restricted by user input

      // the cluster ids are reduced by this factor
      for (unsigned int i = 0; i < n_cells_with_ghosts; ++i)
        cell_cluster_ids[i] = cell_cluster_ids[i] / cluster_diff;
    }

    // bring information into a distributed vector
    distr_cluster_ids.reinit(elerowset, elecolset, MPI_COMM_WORLD);
    distr_cluster_ids = 0.;
    for (unsigned int i = 0; i < n_cells_with_ghosts; ++i)
      for (unsigned int v = 0; v < op.get_matrix_free().n_components_filled(i); ++v)
        distr_cluster_ids[distr_cell_index[i][v]] = cell_cluster_ids[i];

    // communicate the information, important: use VectorOperation::min
    // distr_cluster_ids.compress(VectorOperation::min);
    distr_cluster_ids.compress(
        VectorOperation::insert);  // TODO: update deal in lnm folder to allow for min evaluation
    distr_cluster_ids.update_ghost_values();

    // bring the communicated values back to the cell_cluster_ids
    for (unsigned int i = 0; i < n_cells_with_ghosts; ++i)
    {
      Number minval = std::numeric_limits<Number>::max();
      for (unsigned int v = 0; v < op.get_matrix_free().n_components_filled(i); ++v)
      {
        if (distr_cluster_ids[distr_cell_index[i][v]] < minval)
          minval = distr_cluster_ids[distr_cell_index[i][v]];
      }
      cell_cluster_ids[i] = minval;
    }

    // a difference of one minimum dt and not more (especially in parallel) AND
    // they are only allowed to have slower XOR faster neighbor not both
    bool faster_and_slower_neighbor = true;
    unsigned int count = 0;
    while (faster_and_slower_neighbor && count < 100)
    {
      count++;
      cell_have_faster_neighbor.clear();
      cell_have_slower_neighbor.clear();
      for (unsigned int c = 0; c < n_clusters; ++c)
      {
        for (unsigned int i = 0; i < n_cells; ++i)
        {
          if (cell_cluster_ids[i] == c)
          {
            bool is_anyone_faster = false;
            for (unsigned int v = 0; v < op.get_matrix_free().n_components_filled(i); ++v)
              for (unsigned int n = 0; n < GeometryInfo<Operator::dimension>::faces_per_cell; ++n)
                if (op.get_matrix_free().get_cell_iterator(i, v)->neighbor_index(n) >= 0)
                {
                  if (op.get_matrix_free().get_cell_iterator(i, v)->neighbor(n)->has_children())
                  {
                    for (unsigned int subfaces = 0;
                         subfaces < GeometryInfo<Operator::dimension>::max_children_per_face;
                         ++subfaces)
                    {
                      if (cell_cluster_ids[mf_index[op.get_matrix_free()
                                                        .get_cell_iterator(i, v)
                                                        ->neighbor_child_on_subface(n, subfaces)
                                                        ->active_cell_index()]] < c)
                        is_anyone_faster = true;
                    }
                  }
                  else if (cell_cluster_ids[mf_index[op.get_matrix_free()
                                                         .get_cell_iterator(i, v)
                                                         ->neighbor(n)
                                                         ->active_cell_index()]] < c)
                    is_anyone_faster = true;
                }
            for (unsigned int v = 0; v < op.get_matrix_free().n_components_filled(i); ++v)
            {
              for (unsigned int n = 0; n < GeometryInfo<Operator::dimension>::faces_per_cell; ++n)
              {
                if (op.get_matrix_free().get_cell_iterator(i, v)->neighbor_index(n) >= 0)
                {
                  if (op.get_matrix_free().get_cell_iterator(i, v)->neighbor(n)->has_children())
                  {
                    for (unsigned int subfaces = 0;
                         subfaces < GeometryInfo<Operator::dimension>::max_children_per_face;
                         ++subfaces)
                    {
                      if (cell_cluster_ids[mf_index[op.get_matrix_free()
                                                        .get_cell_iterator(i, v)
                                                        ->neighbor_child_on_subface(n, subfaces)
                                                        ->active_cell_index()]] > c)
                      {
                        if (is_anyone_faster)
                          cell_cluster_ids[mf_index[op.get_matrix_free()
                                                        .get_cell_iterator(i, v)
                                                        ->neighbor_child_on_subface(n, subfaces)
                                                        ->active_cell_index()]] = c;
                        else
                          cell_cluster_ids[mf_index[op.get_matrix_free()
                                                        .get_cell_iterator(i, v)
                                                        ->neighbor_child_on_subface(n, subfaces)
                                                        ->active_cell_index()]] = c + 1;
                      }
                    }
                  }
                  else
                  {
                    if (cell_cluster_ids[mf_index[op.get_matrix_free()
                                                      .get_cell_iterator(i, v)
                                                      ->neighbor(n)
                                                      ->active_cell_index()]] > c)
                    {
                      if (is_anyone_faster)
                        cell_cluster_ids[mf_index[op.get_matrix_free()
                                                      .get_cell_iterator(i, v)
                                                      ->neighbor(n)
                                                      ->active_cell_index()]] = c;
                      else
                        cell_cluster_ids[mf_index[op.get_matrix_free()
                                                      .get_cell_iterator(i, v)
                                                      ->neighbor(n)
                                                      ->active_cell_index()]] = c + 1;
                    }
                  }
                }
              }
            }
          }
        }
      }

      int commcount = 0;
      std::vector<unsigned int> cell_cluster_ids_old(n_cells_with_ghosts, false);
      do
      {
        cell_cluster_ids_old = cell_cluster_ids;
        communicate_cluster_ids();
        commcount++;
      } while (Utilities::MPI::max(int(cell_cluster_ids_old != cell_cluster_ids), MPI_COMM_WORLD) &&
               commcount < 100);
      if (commcount == 100) Assert(false, ExcMessage("clustering communication is too bad"));

      // determine if a cell has a faster neighbor
      cell_have_faster_neighbor.resize(n_cells, false);
      cell_have_slower_neighbor.resize(n_cells, false);

      faster_and_slower_neighbor = false;
      for (unsigned int i = 0; i < n_cells; ++i)
      {
        unsigned int current_c = cell_cluster_ids[i];
        bool is_anyone_faster = false;
        bool is_anyone_slower = false;
        for (unsigned int v = 0; v < op.get_matrix_free().n_components_filled(i); ++v)
        {
          for (unsigned int n = 0; n < GeometryInfo<Operator::dimension>::faces_per_cell; ++n)
          {
            if (op.get_matrix_free().get_cell_iterator(i, v)->neighbor_index(n) >= 0)
            {
              if (op.get_matrix_free().get_cell_iterator(i, v)->neighbor(n)->has_children())
              {
                for (unsigned int subfaces = 0;
                     subfaces < GeometryInfo<Operator::dimension>::max_children_per_face;
                     ++subfaces)
                {
                  if (cell_cluster_ids[mf_index[op.get_matrix_free()
                                                    .get_cell_iterator(i, v)
                                                    ->neighbor_child_on_subface(n, subfaces)
                                                    ->active_cell_index()]] < current_c)
                  {
                    is_anyone_faster = true;
                    cell_have_faster_neighbor[i] = true;
                  }
                  else if (cell_cluster_ids[mf_index[op.get_matrix_free()
                                                         .get_cell_iterator(i, v)
                                                         ->neighbor_child_on_subface(n, subfaces)
                                                         ->active_cell_index()]] > current_c)
                  {
                    is_anyone_slower = true;
                    cell_have_slower_neighbor[i] = true;
                  }
                  if (int(cell_cluster_ids[mf_index[op.get_matrix_free()
                                                        .get_cell_iterator(i, v)
                                                        ->neighbor_child_on_subface(n, subfaces)
                                                        ->active_cell_index()]]) <
                          int(current_c - 1) ||
                      cell_cluster_ids[mf_index[op.get_matrix_free()
                                                    .get_cell_iterator(i, v)
                                                    ->neighbor_child_on_subface(n, subfaces)
                                                    ->active_cell_index()]] > current_c + 1)
                  {
                    faster_and_slower_neighbor = true;
                  }
                }
              }
              else
              {
                if (cell_cluster_ids[mf_index[op.get_matrix_free()
                                                  .get_cell_iterator(i, v)
                                                  ->neighbor(n)
                                                  ->active_cell_index()]] < current_c)
                {
                  is_anyone_faster = true;
                  cell_have_faster_neighbor[i] = true;
                }
                else if (cell_cluster_ids[mf_index[op.get_matrix_free()
                                                       .get_cell_iterator(i, v)
                                                       ->neighbor(n)
                                                       ->active_cell_index()]] > current_c)
                {
                  is_anyone_slower = true;
                  cell_have_slower_neighbor[i] = true;
                }
                if (int(cell_cluster_ids[mf_index[op.get_matrix_free()
                                                      .get_cell_iterator(i, v)
                                                      ->neighbor(n)
                                                      ->active_cell_index()]]) <
                        int(current_c - 1) ||
                    cell_cluster_ids[mf_index[op.get_matrix_free()
                                                  .get_cell_iterator(i, v)
                                                  ->neighbor(n)
                                                  ->active_cell_index()]] > current_c + 1)
                {
                  faster_and_slower_neighbor = true;
                }
              }
            }
          }
        }
        if (is_anyone_faster && is_anyone_slower)
        {
          faster_and_slower_neighbor = true;
        }
      }
      faster_and_slower_neighbor =
          Utilities::MPI::max(double(faster_and_slower_neighbor), MPI_COMM_WORLD);
    }
    if (count == 100)
      Assert(false, ExcMessage("could not find suited clustering in 100 iterations"));
    std::cout << "needed to iterate " << count << " times to get valid clustering" << std::endl;

    // communicate the cell_have_faster_neighbor and cell_have_slower_neighbor vector to a vector
    // with ghost entries
    //{
    distr_cluster_ids = 0.;
    for (unsigned int i = 0; i < n_cells; ++i)
      for (unsigned int v = 0; v < op.get_matrix_free().n_components_filled(i); ++v)
        distr_cluster_ids[distr_cell_index[i][v]] = cell_have_faster_neighbor[i];
    distr_cluster_ids.compress(VectorOperation::add);
    distr_cluster_ids.update_ghost_values();
    cell_have_faster_neighbor.clear();
    cell_have_faster_neighbor.resize(n_cells_with_ghosts, false);
    for (unsigned int i = 0; i < n_cells_with_ghosts; ++i)
    {
      Number max = 0;
      for (unsigned int v = 0; v < op.get_matrix_free().n_components_filled(i); ++v)
        if (distr_cluster_ids[distr_cell_index[i][v]] > max)
          max = distr_cluster_ids[distr_cell_index[i][v]];
      cell_have_faster_neighbor[i] = max;
    }
    distr_cluster_ids = 0.;
    for (unsigned int i = 0; i < n_cells; ++i)
      for (unsigned int v = 0; v < op.get_matrix_free().n_components_filled(i); ++v)
        distr_cluster_ids[distr_cell_index[i][v]] = cell_have_slower_neighbor[i];
    distr_cluster_ids.compress(VectorOperation::add);
    distr_cluster_ids.update_ghost_values();
    cell_have_slower_neighbor.clear();
    cell_have_slower_neighbor.resize(n_cells_with_ghosts, false);
    for (unsigned int i = 0; i < n_cells_with_ghosts; ++i)
    {
      Number max = 0;
      for (unsigned int v = 0; v < op.get_matrix_free().n_components_filled(i); ++v)
        if (distr_cluster_ids[distr_cell_index[i][v]] > max)
          max = distr_cluster_ids[distr_cell_index[i][v]];
      cell_have_slower_neighbor[i] = max;
    }
    //}

    // some statistics
    std::vector<unsigned int> num_cell_cluster(n_clusters, 0);
    for (unsigned int i = 0; i < n_cells; ++i) num_cell_cluster[cell_cluster_ids[i]]++;

    // communicate this thing
    Utilities::MPI::sum(num_cell_cluster, MPI_COMM_WORLD, num_cell_cluster);

    int prev = n_clusters;
    for (unsigned int c = prev - 1; c > 0; c--)
      if (num_cell_cluster[c] == 0) n_clusters = c;

    // next stuff to init
    cluster_timelevels.clear();
    cell_timelevels.clear();
    evaluate_cell.clear();
    update_cell.clear();
    evaluate_face.clear();
    phi_to_dst.clear();
    phi_neighbor_to_dst.clear();
    phi_to_fluxmemory.clear();
    phi_neighbor_to_fluxmemory.clear();
    cluster_timelevels.resize(n_clusters, op.time);
    cell_timelevels.resize(n_cells_with_ghosts, op.time);
    evaluate_cell.resize(n_cells_with_ghosts, false);
    update_cell.resize(n_cells_with_ghosts, false);
    evaluate_face.resize(op.get_matrix_free().faces.size(), false);
    phi_to_dst.resize(op.get_matrix_free().faces.size());
    phi_neighbor_to_dst.resize(op.get_matrix_free().faces.size());
    phi_to_fluxmemory.resize(op.get_matrix_free().faces.size());
    phi_neighbor_to_fluxmemory.resize(op.get_matrix_free().faces.size());


    // set the cluster time steps
    cluster_timestepmultiples.resize(n_clusters);
    for (unsigned c = 0; c < n_clusters; ++c)
      cluster_timestepmultiples[c] = cluster_diff * int(c) + 1;

    // for (unsigned c=0; c<n_clusters; ++c)
    //  std::cout<<"cluster  "<<c<<" contains "<<num_cell_cluster[c]<<" cells and works with time
    //  step "<<cluster_timestepmultiples[c]*fastest_time_step<<std::endl;

    /*std::cout<<"cluster_timestepmultiples"<<std::endl;
    for (unsigned c=0; c<n_clusters; ++c)
      std::cout<<cluster_timestepmultiples[c]<<std::endl;
    std::cout<<"n_clusters "<<n_clusters<<" clusterdiff "<<cluster_diff<<" levelmax
    "<<cluster_diff*(n_clusters-1)+1<<std::endl;*/

    // determine how one updates from one global time step to the next
    // for a general clustering, we have to check the update criteria for each cluster
    // restriction is: neighboring clusters are only one apart
    std::vector<unsigned int> templevels(n_clusters);
    std::vector<Number> t_a_b(2, 0.0);

    // need this clear for adaptivity
    cluster_update_times.clear();
    cluster_update_order.clear();
    unsigned int max_level = cluster_diff * (n_clusters - 1) + 1;
    while (*std::min_element(templevels.begin(), templevels.end()) < max_level)
    {
      for (unsigned c = 0; c < n_clusters; ++c)
      {
        if (templevels[c] < cluster_diff * (n_clusters - 1) + 1)
        {
          if (c == 0)
          {
            if (std::min(max_level, templevels[c] + cluster_timestepmultiples[c]) <=
                std::min(max_level, templevels[c + 1] + cluster_timestepmultiples[c + 1]))
            {
              t_a_b[0] = templevels[c] * fastest_time_step;
              templevels[c] += cluster_timestepmultiples[c];
              cluster_update_order.push_back(c);
              t_a_b[1] = templevels[c] * fastest_time_step;
              cluster_update_times.push_back(t_a_b);
            }
          }
          else if (c == n_clusters - 1)
          {
            if (std::min(max_level, templevels[c] + cluster_timestepmultiples[c]) <=
                std::min(max_level, templevels[c - 1] + cluster_timestepmultiples[c - 1]))
            {
              t_a_b[0] = templevels[c] * fastest_time_step;
              templevels[c] += cluster_timestepmultiples[c];
              cluster_update_order.push_back(c);
              t_a_b[1] = templevels[c] * fastest_time_step;
              cluster_update_times.push_back(t_a_b);
            }
          }
          else
          {
            if (std::min(max_level, templevels[c] + cluster_timestepmultiples[c]) <=
                    std::min(max_level, templevels[c + 1] + cluster_timestepmultiples[c + 1]) &&
                std::min(max_level, templevels[c] + cluster_timestepmultiples[c]) <=
                    std::min(max_level, templevels[c - 1] + cluster_timestepmultiples[c - 1]))
            {
              t_a_b[0] = templevels[c] * fastest_time_step;
              templevels[c] += cluster_timestepmultiples[c];
              cluster_update_order.push_back(c);
              t_a_b[1] = templevels[c] * fastest_time_step;
              cluster_update_times.push_back(t_a_b);
            }
          }
        }
      }
    }

    // correct for too long updates
    for (unsigned c = 0; c < n_clusters; ++c)
      if (templevels[c] > cluster_diff * (n_clusters - 1) + 1)
      {
        templevels[c] = cluster_diff * (n_clusters - 1) + 1;
      }

    // tell the wave equation problem with which maximum time step we update
    op.set_time_step_size(fastest_time_step * (cluster_diff * (n_clusters - 1) + 1));

    n_updates = cluster_update_order.size();

    for (unsigned c = 0; c < cluster_update_times.size(); ++c)
      if (cluster_update_times[c][1] > (cluster_diff * (n_clusters - 1) + 1) * fastest_time_step)
        cluster_update_times[c][1] = (cluster_diff * (n_clusters - 1) + 1) * fastest_time_step;

    for (unsigned c = 0; c < cluster_update_times.size(); ++c)
    {
      cluster_update_times[c][0] += op.time;
      cluster_update_times[c][1] += op.time;
    }
    // for(unsigned c=0; c<cluster_update_order.size(); ++c)
    //  std::cout<<"c "<<c<<" updatecluster "<<cluster_update_order[c]<<std::endl;
  }



  template <typename Number>
  template <typename Operator>
  void ClusterManager<Number>::perform_time_step(const Operator &op,
      const std::vector<parallel::distributed::Vector<Number>> &src,
      std::vector<parallel::distributed::Vector<Number>> &dst) const
  {
    for (unsigned int d = 0; d < src.size(); ++d)
    {
      dst[d] = 0.;
      state[d] = src[d];
    }

    if (op.use_ader_post)
    {
      // standard reconstruction:
      op.apply(state, improvedgraddiv, op.time, op.time_step);
    }

    // update cell times
    for (unsigned int e = 0; e < n_cells_with_ghosts; ++e)
      cell_timelevels[e] = cluster_timelevels[cell_cluster_ids[e]];

    for (unsigned int cycle = 0; cycle < n_updates; ++cycle)
    {
      // we want to update the elements of cluster cluster_update_order[cycle]
      unsigned int actual_cluster = cluster_update_order[cycle];
      double actual_cluster_time = cluster_timelevels[actual_cluster];

      // security check
      if (actual_cluster_time != cluster_update_times[cycle][0])
      {
        std::cout << "actual cluster time " << actual_cluster_time << " cluster update times "
                  << cluster_update_times[cycle][0] << std::endl;
        std::cout << "cycle " << cycle << " actual cluster " << actual_cluster << std::endl;
        Assert(false, ExcMessage("cluster time missmatch"));
      }
      is_fluxmemory_considered = true;

      for (unsigned int e = 0; e < n_cells_with_ghosts; ++e)
      {
        if (cell_cluster_ids[e] == actual_cluster)
          update_cell[e] = true;
        else
          update_cell[e] = false;
      }

      double faster_cluster_time = 0.0;
      if (actual_cluster > 0)
      {
        faster_cluster_time = cluster_timelevels[actual_cluster - 1];
        t1fa = std::max(faster_cluster_time, actual_cluster_time);
        t2fa = std::min(
            faster_cluster_time + cluster_timestepmultiples[actual_cluster - 1] * fastest_time_step,
            actual_cluster_time + cluster_timestepmultiples[actual_cluster] * fastest_time_step);
        t2fa = std::min(t2fa, op.time);
      }
      else
      {
        t1fa = 0.0;
        t2fa = 0.0;
      }

      if (actual_cluster < n_clusters - 1)
      {
        double slower_cluster_time = cluster_timelevels[actual_cluster + 1];
        t1sl = std::max(slower_cluster_time, actual_cluster_time);
        t2sl = std::min(
            slower_cluster_time + cluster_timestepmultiples[actual_cluster + 1] * fastest_time_step,
            actual_cluster_time + cluster_timestepmultiples[actual_cluster] * fastest_time_step);
        t2sl = std::min(t2sl, op.time);
      }
      else
      {
        t1sl = 0.0;
        t2sl = 0.0;
      }
      t1sa = actual_cluster_time;
      t2sa = std::min(Number(actual_cluster_time +
                             cluster_timestepmultiples[actual_cluster] * fastest_time_step),
          Number(op.time));

      dt = std::min(cluster_timestepmultiples[actual_cluster] * fastest_time_step,
          op.time - cluster_timelevels[actual_cluster]);

      update_elements(op, dst, state, actual_cluster);

      // update time level of cluster
      cluster_timelevels[actual_cluster] = cluster_update_times[cycle][1];

      // update cell times
      for (unsigned int e = 0; e < n_cells_with_ghosts; ++e)
        cell_timelevels[e] = cluster_timelevels[cell_cluster_ids[e]];

      // in case we want superconvergence, we need to do the reconstruction step as explained in the
      // paper
      if (op.use_ader_post)
      {
        // init
        is_fluxmemory_considered = false;
        temporary_recon_state = state;

        // current cluster advanced in time
        actual_cluster_time = cluster_timelevels[actual_cluster];

        // which cells do we have to update to be able to perform reconstruction
        // all cells who are neighbor of actual cluster and  have lower or higher time level
        // start with faster cluster:
        std::vector<bool> temp_faster(n_cells);
        std::vector<bool> temp_slower(n_cells);
        if (actual_cluster > 0)
          if (std::abs(actual_cluster_time - cluster_timelevels[actual_cluster - 1]) >
              relative_tolerance * fastest_time_step)
          {
            for (unsigned int e = 0; e < n_cells; ++e)
            {
              if (cell_cluster_ids[e] == actual_cluster - 1 && cell_have_slower_neighbor[e])
                temp_faster[e] = true;
              else
                temp_faster[e] = false;
            }
            // expand this selection by the neighbors of those already set (later we need neighbor
            // of neighbor info!)
            for (unsigned int e = 0; e < n_cells; ++e) update_cell[e] = false;
            for (unsigned int e = 0; e < n_cells; ++e)
            {
              for (unsigned int v = 0; v < op.get_matrix_free().n_components_filled(e); ++v)
                for (unsigned int n = 0; n < GeometryInfo<Operator::dimension>::faces_per_cell; ++n)
                  if (op.get_matrix_free().get_cell_iterator(e, v)->neighbor_index(n) >= 0)
                  {
                    if (temp_faster[op.get_matrix_free().get_cell_iterator(e, v)->neighbor_index(
                                        n) /
                                    Operator::n_vect] &&
                        cell_cluster_ids[e] < actual_cluster)
                      update_cell[e] = true;
                    else
                      update_cell[e] = false;
                  }
            }
            // combine both
            for (unsigned int e = 0; e < n_cells; ++e)
              update_cell[e] = update_cell[e] || temp_faster[e];
            temp_faster = update_cell;

            t1fa = cluster_timelevels[actual_cluster - 1];
            t2fa = actual_cluster_time;

            t1sl = t1fa;
            t2sl = t2fa;

            t1sa = t1fa;
            t2sa = t2fa;

            dt = t2fa - t1fa;

            update_elements(op, dst, temporary_recon_state, actual_cluster - 1, false);
          }
        // slower cluster
        if (actual_cluster < n_clusters - 1)
          if (std::abs(actual_cluster_time - cluster_timelevels[actual_cluster + 1]) >
              relative_tolerance * fastest_time_step)
          {
            for (unsigned int e = 0; e < n_cells; ++e)
            {
              if (cell_cluster_ids[e] == actual_cluster + 1 && cell_have_faster_neighbor[e])
                temp_slower[e] = true;
              else
                temp_slower[e] = false;
            }
            // expand this selection by the neighbors of those already set (later we need neighbor
            // of neighbor info!)
            for (unsigned int e = 0; e < n_cells; ++e) update_cell[e] = false;
            for (unsigned int e = 0; e < n_cells; ++e)
            {
              for (unsigned int v = 0; v < op.get_matrix_free().n_components_filled(e); ++v)
                for (unsigned int n = 0; n < GeometryInfo<Operator::dimension>::faces_per_cell; ++n)
                  if (op.get_matrix_free().get_cell_iterator(e, v)->neighbor_index(n) >= 0)
                    if (temp_slower[op.get_matrix_free().get_cell_iterator(e, v)->neighbor_index(
                                        n) /
                                    Operator::n_vect] &&
                        cell_cluster_ids[e] > actual_cluster)
                    {
                      update_cell[e] = true;
                    }
            }
            // combine both
            for (unsigned int e = 0; e < n_cells; ++e)
              update_cell[e] = update_cell[e] || temp_slower[e];
            temp_slower = update_cell;

            t1fa = cluster_timelevels[actual_cluster + 1];
            t2fa = actual_cluster_time;

            t1sl = t1fa;
            t2sl = t2fa;

            t1sa = t1fa;
            t2sa = t2fa;

            dt = t2fa - t1fa;

            update_elements(op, dst, temporary_recon_state, actual_cluster + 1, false);
          }
        // now, the neighboring elements are at the required time level

        // reconstruction:
        for (unsigned int e = 0; e < n_cells; ++e)
        {
          if (cell_cluster_ids[e] == actual_cluster)
            evaluate_cell[e] = true;
          else
            evaluate_cell[e] = false;
        }

        for (unsigned int f = 0; f < op.get_matrix_free().faces.size(); ++f)
        {
          evaluate_face[f] = false;
          phi_to_dst[f] = std::bitset<Operator::n_vect>(false);
          phi_neighbor_to_dst[f] = std::bitset<Operator::n_vect>(false);
          if (f < op.get_matrix_free().n_macro_inner_faces())
          {
            for (unsigned int v = 0; v < Operator::n_vect; ++v)
            {
              if (op.get_matrix_free().faces[f].left_cell[v] != numbers::invalid_unsigned_int &&
                  op.get_matrix_free().faces[f].right_cell[v] != numbers::invalid_unsigned_int)
                if (evaluate_cell[op.get_matrix_free().faces[f].left_cell[v] / Operator::n_vect] ||
                    evaluate_cell[op.get_matrix_free().faces[f].right_cell[v] / Operator::n_vect])
                {
                  evaluate_face[f] = true;
                  // set masks
                  if (evaluate_cell[op.get_matrix_free().faces[f].left_cell[v] / Operator::n_vect])
                    phi_to_dst[f][v] = true;
                  if (evaluate_cell[op.get_matrix_free().faces[f].right_cell[v] / Operator::n_vect])
                    phi_neighbor_to_dst[f][v] = true;
                }
            }
          }
          else
          {
            for (unsigned int v = 0; v < Operator::n_vect; ++v)
            {
              if (op.get_matrix_free().faces[f].left_cell[v] != numbers::invalid_unsigned_int)
                if (evaluate_cell[op.get_matrix_free().faces[f].left_cell[v] / Operator::n_vect])
                {
                  evaluate_face[f] = true;
                  phi_to_dst[f][v] = true;
                }
            }
          }
        }
        op.reconstruct_div_grad(temporary_recon_state, improvedgraddiv);
      }

    }  // for(unsigned int cycle = 0; cycle<n_updates; ++cycle)

    // update the update times
    for (unsigned c = 0; c < cluster_update_times.size(); ++c)
      for (unsigned int n = 0; n < 2; ++n) cluster_update_times[c][n] += op.time_step;

    // tell the time integrator about the new state
    for (unsigned int d = 0; d < src.size(); ++d) dst[d] = state[d];
  }



  template <typename Number>
  template <typename Operator>
  void ClusterManager<Number>::update_elements(const Operator &op,
      std::vector<parallel::distributed::Vector<Number>> &dst,
      std::vector<parallel::distributed::Vector<Number>> &local_state,
      const unsigned int actual_cluster, const bool write_to_fluxmemory) const
  {
    double actual_cluster_time = cluster_timelevels[actual_cluster];

    // security  checks
    for (unsigned int e = 0; e < n_cells_with_ghosts; ++e)
    {
      if (update_cell[e] == true && cell_cluster_ids[e] != actual_cluster)
        Assert(false,
            ExcMessage("it is not allowed to call update elements on cells of differing cluster!"));
      if (update_cell[e] == true)
        if (cell_timelevels[e] != actual_cluster_time)
          Assert(false,
              ExcMessage(
                  "it is not allowed to call update elements on cells of differing time level!"));
    }

    // fill neighbor vector
    std::vector<bool> is_neighbor_of_update_cell(n_cells_with_ghosts, false);
    for (unsigned int e = 0; e < n_cells_with_ghosts; ++e)
    {
      if (update_cell[e] == false)
      {
        for (unsigned int v = 0; v < op.get_matrix_free().n_components_filled(e); ++v)
          for (unsigned int n = 0; n < GeometryInfo<Operator::dimension>::faces_per_cell; ++n)
            if (op.get_matrix_free().get_cell_iterator(e, v)->neighbor_index(n) >= 0)
            {
              if (op.get_matrix_free().get_cell_iterator(e, v)->neighbor(n)->has_children())
              {
                for (unsigned int subfaces = 0;
                     subfaces < GeometryInfo<Operator::dimension>::max_children_per_face;
                     ++subfaces)
                {
                  if (update_cell[mf_index[op.get_matrix_free()
                                               .get_cell_iterator(e, v)
                                               ->neighbor_child_on_subface(n, subfaces)
                                               ->active_cell_index()]])
                    is_neighbor_of_update_cell[e] = true;
                }
              }
              else
              {
                if (update_cell[mf_index[op.get_matrix_free()
                                             .get_cell_iterator(e, v)
                                             ->neighbor(n)
                                             ->active_cell_index()]])
                  is_neighbor_of_update_cell[e] = true;
              }
            }
      }
    }

    // contribution from the faster cluster
    if (actual_cluster > 0)
    {
      // prepare everything for face evaluation between actual_cluster and actual_cluster-1
      t1 = t1fa;
      t2 = t2fa;

      if (std::abs(t2 - t1) > relative_tolerance * fastest_time_step)  // only do this if necessary
      {
        // setup the element list that need evaluation (all of update_cell who have faster neighbor
        // and all of actual_cluster-1 who are neighbor of update_cell)
        for (unsigned int e = 0; e < n_cells_with_ghosts; ++e)
        {
          if (cell_cluster_ids[e] == actual_cluster && cell_have_faster_neighbor[e] &&
              update_cell[e])
            evaluate_cell[e] = true;
          else if (cell_cluster_ids[e] == actual_cluster - 1 && cell_have_slower_neighbor[e] &&
                   is_neighbor_of_update_cell[e])
            evaluate_cell[e] = true;
          else
            evaluate_cell[e] = false;
        }

        // setup the face list that need evaluation (if face is in between update cell and faster
        // cluster) also, setup the masks for the write functions
        for (unsigned int f = 0; f < op.get_matrix_free().faces.size(); ++f)
        {
          evaluate_face[f] = false;
          phi_to_dst[f] = std::bitset<Operator::n_vect>(false);
          phi_neighbor_to_dst[f] = std::bitset<Operator::n_vect>(false);
          phi_to_fluxmemory[f] = std::bitset<Operator::n_vect>(false);
          phi_neighbor_to_fluxmemory[f] = std::bitset<Operator::n_vect>(false);
          // inner faces
          if (f < op.get_matrix_free().n_macro_inner_faces())
          {
            for (unsigned int v = 0; v < Operator::n_vect; ++v)
              if (op.get_matrix_free().faces[f].left_cell[v] != numbers::invalid_unsigned_int &&
                  op.get_matrix_free().faces[f].right_cell[v] != numbers::invalid_unsigned_int)
                if (evaluate_cell[op.get_matrix_free().faces[f].left_cell[v] / Operator::n_vect] &&
                    evaluate_cell[op.get_matrix_free().faces[f].right_cell[v] / Operator::n_vect])
                {
                  // set face evaluation
                  if ((cell_cluster_ids[op.get_matrix_free().faces[f].left_cell[v] /
                                        Operator::n_vect] == actual_cluster &&
                          cell_cluster_ids[op.get_matrix_free().faces[f].right_cell[v] /
                                           Operator::n_vect] != actual_cluster &&
                          is_neighbor_of_update_cell[op.get_matrix_free().faces[f].right_cell[v] /
                                                     Operator::n_vect]) ||
                      (cell_cluster_ids[op.get_matrix_free().faces[f].left_cell[v] /
                                        Operator::n_vect] != actual_cluster &&
                          cell_cluster_ids[op.get_matrix_free().faces[f].right_cell[v] /
                                           Operator::n_vect] == actual_cluster &&
                          is_neighbor_of_update_cell[op.get_matrix_free().faces[f].left_cell[v] /
                                                     Operator::n_vect]))
                  {
                    evaluate_face[f] = true;
                    // set masks
                    if (cell_cluster_ids[op.get_matrix_free().faces[f].left_cell[v] /
                                         Operator::n_vect] == actual_cluster)
                    {
                      phi_to_dst[f][v] = true;
                      if (write_to_fluxmemory) phi_neighbor_to_fluxmemory[f][v] = true;
                    }
                    else
                    {
                      if (write_to_fluxmemory) phi_to_fluxmemory[f][v] = true;
                      phi_neighbor_to_dst[f][v] = true;
                    }
                  }
                }
          }
          else  // boundary faces
          {
            // no contribution from boundary faces in this stage
          }
        }
        // do first ader
        op.evaluate_cells_and_faces_first_ader(local_state, dst);
      }
    }

    // contribution from the slower cluster
    if (actual_cluster < n_clusters - 1)
    {
      // prepare everything for face evaluation between actual_cluster and actual_cluster+1
      t1 = t1sl;
      t2 = t2sl;
      if (std::abs(t2 - t1) > relative_tolerance * fastest_time_step)  // only do this if necessary
      {
        // setup the element list that need evaluation (all of actual_cluster who have slower
        // neighbor and all of actual_cluster+1 who have faster neighbor)
        for (unsigned int e = 0; e < n_cells_with_ghosts; ++e)
        {
          if (cell_cluster_ids[e] == actual_cluster && cell_have_slower_neighbor[e] &&
              update_cell[e])
            evaluate_cell[e] = true;
          else if (cell_cluster_ids[e] == actual_cluster + 1 && cell_have_faster_neighbor[e] &&
                   is_neighbor_of_update_cell[e])
            evaluate_cell[e] = true;
          else
            evaluate_cell[e] = false;
        }

        // setup the face list that need evaluation (if face is in between actual cluster and slower
        // cluster) also, setup the masks for the write functions
        for (unsigned int f = 0; f < op.get_matrix_free().faces.size(); ++f)
        {
          evaluate_face[f] = false;
          phi_to_dst[f] = std::bitset<Operator::n_vect>(false);
          phi_neighbor_to_dst[f] = std::bitset<Operator::n_vect>(false);
          phi_to_fluxmemory[f] = std::bitset<Operator::n_vect>(false);
          phi_neighbor_to_fluxmemory[f] = std::bitset<Operator::n_vect>(false);
          // inner faces
          if (f < op.get_matrix_free().n_macro_inner_faces())
          {
            for (unsigned int v = 0; v < Operator::n_vect; ++v)
              if (op.get_matrix_free().faces[f].left_cell[v] != numbers::invalid_unsigned_int &&
                  op.get_matrix_free().faces[f].right_cell[v] != numbers::invalid_unsigned_int)
              {
                if (evaluate_cell[op.get_matrix_free().faces[f].left_cell[v] / Operator::n_vect] &&
                    evaluate_cell[op.get_matrix_free().faces[f].right_cell[v] / Operator::n_vect])
                {
                  // set face evaluation
                  if ((cell_cluster_ids[op.get_matrix_free().faces[f].left_cell[v] /
                                        Operator::n_vect] == actual_cluster &&
                          cell_cluster_ids[op.get_matrix_free().faces[f].right_cell[v] /
                                           Operator::n_vect] != actual_cluster &&
                          is_neighbor_of_update_cell[op.get_matrix_free().faces[f].right_cell[v] /
                                                     Operator::n_vect]) ||
                      (cell_cluster_ids[op.get_matrix_free().faces[f].left_cell[v] /
                                        Operator::n_vect] != actual_cluster &&
                          cell_cluster_ids[op.get_matrix_free().faces[f].right_cell[v] /
                                           Operator::n_vect] == actual_cluster &&
                          is_neighbor_of_update_cell[op.get_matrix_free().faces[f].left_cell[v] /
                                                     Operator::n_vect]))
                  {
                    evaluate_face[f] = true;
                    // set masks
                    if (cell_cluster_ids[op.get_matrix_free().faces[f].left_cell[v] /
                                         Operator::n_vect] == actual_cluster)
                    {
                      phi_to_dst[f][v] = true;
                      if (write_to_fluxmemory) phi_neighbor_to_fluxmemory[f][v] = true;
                    }
                    else
                    {
                      if (write_to_fluxmemory) phi_to_fluxmemory[f][v] = true;
                      phi_neighbor_to_dst[f][v] = true;
                    }
                  }
                }
              }
          }
          else  // boundary faces
          {
            // no contribution from boundary faces in this stage
          }
        }
        // do first ader
        op.evaluate_cells_and_faces_first_ader(local_state, dst);
      }
    }  // check the slower cluster

    // contribution from the same cluster
    {
      t1 = t1sa;
      t2 = t2sa;

      // setup the element list that need evaluation (all of actual_cluster)
      for (unsigned int e = 0; e < n_cells_with_ghosts; ++e)
      {
        if (cell_cluster_ids[e] == actual_cluster && update_cell[e])
          evaluate_cell[e] = true;
        else if (cell_cluster_ids[e] == actual_cluster && is_neighbor_of_update_cell[e])
          evaluate_cell[e] = true;
        else
          evaluate_cell[e] = false;
      }

      // setup the face list that need evaluation (if face is in between two cells with actual
      // cluster) also, setup the masks for the write functions
      for (unsigned int f = 0; f < op.get_matrix_free().faces.size(); ++f)
      {
        evaluate_face[f] = false;
        phi_to_dst[f] = std::bitset<Operator::n_vect>(false);
        phi_neighbor_to_dst[f] = std::bitset<Operator::n_vect>(false);
        phi_to_fluxmemory[f] = std::bitset<Operator::n_vect>(false);
        phi_neighbor_to_fluxmemory[f] = std::bitset<Operator::n_vect>(false);
        // inner faces
        if (f < op.get_matrix_free().n_macro_inner_faces())
        {
          for (unsigned int v = 0; v < Operator::n_vect; ++v)
            if (op.get_matrix_free().faces[f].left_cell[v] != numbers::invalid_unsigned_int &&
                op.get_matrix_free().faces[f].right_cell[v] != numbers::invalid_unsigned_int)
            {
              if ((update_cell[op.get_matrix_free().faces[f].left_cell[v] / Operator::n_vect] &&
                      update_cell[op.get_matrix_free().faces[f].right_cell[v] /
                                  Operator::n_vect]) ||
                  (update_cell[op.get_matrix_free().faces[f].left_cell[v] / Operator::n_vect] &&
                      evaluate_cell[op.get_matrix_free().faces[f].right_cell[v] /
                                    Operator::n_vect] &&
                      cell_cluster_ids[op.get_matrix_free().faces[f].right_cell[v] /
                                       Operator::n_vect] == actual_cluster) ||
                  (evaluate_cell[op.get_matrix_free().faces[f].left_cell[v] / Operator::n_vect] &&
                      update_cell[op.get_matrix_free().faces[f].right_cell[v] / Operator::n_vect] &&
                      cell_cluster_ids[op.get_matrix_free().faces[f].left_cell[v] /
                                       Operator::n_vect] == actual_cluster))
              {
                evaluate_face[f] = true;
                // set masks
                phi_to_dst[f][v] = true;
                phi_neighbor_to_dst[f][v] = true;
              }
            }
        }
        else  // boundary faces
        {
          for (unsigned int v = 0; v < Operator::n_vect; ++v)
          {
            if (op.get_matrix_free().faces[f].left_cell[v] != numbers::invalid_unsigned_int)
              if (evaluate_cell[op.get_matrix_free().faces[f].left_cell[v] / Operator::n_vect])
              {
                evaluate_face[f] = true;
                phi_to_dst[f][v] = true;
              }
          }
        }
      }
      // do first ader
      op.evaluate_cells_and_faces_first_ader(local_state, dst);
    }

    // sum the flux memory contribution from all processors
    {
      op.communicate_flux_memory();
    }

    // do the update
    {
      // the evaluate cell vector can be reused
      // only the time step size is needed
      op.evaluate_cells_second_ader(local_state, dst);
    }

    for (unsigned int d = 0; d < dst.size(); ++d)
    {
      local_state[d] -= dst[d];
      dst[d] = 0.;
    }
  }

  template <int dim, int fe_degree, typename Number>
  WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree,
      Number>::WaveEquationOperationAcousticWaveADERLTS(const std::vector<const DoFHandler<dim> *>
                                                            &dof_handlers,
      Teuchos::RCP<DRT::DiscretizationHDG> &discret,
      Teuchos::RCP<Function<dim>> boundary_conditions, Teuchos::RCP<Function<dim>> source_term,
      value_type time_step_in, value_type cfl_number, int sourceno,
      Teuchos::RCP<PATMonitorManager> monitormanagerin)
      : WaveEquationOperationAcousticWaveADER<dim, fe_degree, Number>(dof_handlers, discret,
            boundary_conditions, source_term, time_step_in, sourceno, monitormanagerin)
  {
    const Teuchos::ParameterList &acouparams = DRT::Problem::Instance()->AcousticParams();
    int max_n_clusters = acouparams.get<int>("MAX_NUM_CLUSTERS");
    int max_diff_clusters = acouparams.get<int>("MAX_DIFF_CLUSTERS");

    // initialize flux memory
    flux_memory.resize(dim + 1);
    this->data.initialize_dof_vector(flux_memory[0]);
    for (unsigned int d = 1; d < flux_memory.size(); ++d) flux_memory[d] = flux_memory[0];
    cluster_manager.state = flux_memory;
    if (this->use_ader_post) cluster_manager.improvedgraddiv = flux_memory;

    cluster_manager.setup(*this, max_n_clusters, max_diff_clusters, cfl_number);
  }

  template <int dim, int fe_degree, typename Number>
  void
  WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree, Number>::local_apply_firstader_domain(
      const MatrixFree<dim, value_type> &,
      std::vector<parallel::distributed::Vector<value_type>> &dst,
      const std::vector<parallel::distributed::Vector<value_type>> &src,
      const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, value_type> &phi_eval =
        this->mass_matrix_data->phi[0];
    FEEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, value_type> phi_spectral(this->data, 2);

    // cell loop
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      if (cluster_manager.is_evaluate_cell(cell))
      {
        value_type t1 = cluster_manager.get_t1();
        value_type t2 = cluster_manager.get_t2();
        value_type te = cluster_manager.get_te(cell);

        this->integrate_taylor_cauchykovalewski(
            cell, phi_eval, phi_spectral, src, t2, t1, te, cluster_manager.improvedgraddiv);
        phi_eval.set_dof_values(dst, 0);
      }
    }  // for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
  }



  template <int dim, int fe_degree, typename Number>
  void
  WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree, Number>::local_apply_secondader_domain(
      const MatrixFree<dim, value_type> &,
      std::vector<parallel::distributed::Vector<value_type>> &dst,
      const std::vector<parallel::distributed::Vector<value_type>> &src,
      const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    // for calculation of higher spatial derivatives
    //{
    const unsigned int n_q_points =
        FEEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, value_type>::n_q_points;
    FEEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, value_type> &phi_eval =
        this->mass_matrix_data->phi[0];
    FEEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, value_type> help_eval(
        this->data);  // for memory variable and update of src
    FEEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, value_type> phi_spectral(this->data, 2);
    //}

    // cell loop
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      if (cluster_manager.is_update_cell(cell))
      {
        value_type dt = cluster_manager.get_dt();
        this->integrate_taylor_cauchykovalewski(
            cell, phi_eval, phi_spectral, src, dt, 0., 0., cluster_manager.improvedgraddiv);

        // the standard business analog to local_apply_firstaderlts is done
        // now comes the update!
        {
          phi_eval.evaluate(true, true, false);

          const VectorizedArray<value_type> rho = this->densities[cell];
          const VectorizedArray<value_type> rho_inv = 1. / this->densities[cell];
          const VectorizedArray<value_type> c_sq = this->speeds[cell] * this->speeds[cell];

          for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const Tensor<1, dim + 1, Tensor<1, dim, VectorizedArray<value_type>>> v_and_p_grad =
                phi_eval.get_gradient(q);
            const Tensor<1, dim + 1, VectorizedArray<value_type>> v_and_p = phi_eval.get_value(q);

            Tensor<1, dim + 1, VectorizedArray<value_type>> temp_value;
            for (unsigned int d = 0; d < dim; ++d) temp_value[d] = rho_inv * v_and_p_grad[dim][d];
            phi_eval.submit_value(temp_value, q);

            Tensor<1, dim + 1, Tensor<1, dim, VectorizedArray<value_type>>> temp_gradient;
            for (unsigned int d = 0; d < dim; ++d) temp_gradient[dim][d] = -rho * c_sq * v_and_p[d];

            phi_eval.submit_gradient(temp_gradient, q);
          }

          phi_eval.integrate(true, true);

          // add memory variable
          unsigned int dofs_per_cell = phi_eval.dofs_per_cell;
          help_eval.reinit(cell);
          if (cluster_manager.is_fluxmemory_considered)
          {
            help_eval.read_dof_values(flux_memory, 0);
            for (unsigned j = 0; j < dofs_per_cell; ++j)
              for (unsigned int d = 0; d < dim + 1; ++d)
              {
                phi_eval.begin_dof_values()[d * dofs_per_cell + j] +=
                    help_eval.begin_dof_values()[d * dofs_per_cell + j];
                help_eval.begin_dof_values()[d * dofs_per_cell + j] = 0.;
              }
            help_eval.set_dof_values(flux_memory,
                0);  // tell the flux_memory variable, that some of its values are reset
          }

          // add face contribution (stored in dst) and reset dst to zero
          help_eval.read_dof_values(dst, 0);
          for (unsigned j = 0; j < dofs_per_cell; ++j)
            for (unsigned int d = 0; d < dim + 1; ++d)
            {
              phi_eval.begin_dof_values()[d * dofs_per_cell + j] +=
                  help_eval.begin_dof_values()[d * dofs_per_cell + j];
              help_eval.begin_dof_values()[d * dofs_per_cell + j] = 0.;
            }
          help_eval.set_dof_values(dst, 0);

          // apply inverse mass matrix
          this->mass_matrix_data->inverse.apply(this->mass_matrix_data->coefficients, dim + 1,
              phi_eval.begin_dof_values(), phi_eval.begin_dof_values());
          //}
          phi_eval.set_dof_values(dst, 0);
        }
      }
    }  // for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
  }



  // does nothing
  template <int dim, int fe_degree, typename Number>
  void WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree, Number>::local_apply_dummy_domain(
      const MatrixFree<dim, value_type> &, std::vector<parallel::distributed::Vector<value_type>> &,
      const std::vector<parallel::distributed::Vector<value_type>> &,
      const std::pair<unsigned int, unsigned int> &) const
  {
  }



  template <int dim, int fe_degree, typename Number>
  void WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree, Number>::local_apply_ader_face(
      const MatrixFree<dim, value_type> &,
      std::vector<parallel::distributed::Vector<value_type>> &dst,
      const std::vector<parallel::distributed::Vector<value_type>> &src,
      const std::pair<unsigned int, unsigned int> &face_range) const
  {
    // There is some overhead in the methods in FEEvaluation, so it is faster
    // to combine pressure and velocity in the same object and just combine
    // them at the level of quadrature points
    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, value_type> phi(
        this->data, true, 0, 0, 0, true);
    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, value_type> phi_neighbor(
        this->data, false, 0, 0, 0, true);

    for (unsigned int face = face_range.first; face < face_range.second; face++)
    {
      if (cluster_manager.is_evaluate_face(face))
      {
        this->evaluate_inner_face(phi, phi_neighbor, src, face, -1.0);

        // get bitsets from cluster manager
        phi.distribute_local_to_global(dst, 0, cluster_manager.get_phi_to_dst(face));
        phi.distribute_local_to_global(flux_memory, 0, cluster_manager.get_phi_to_fluxmemory(face));
        phi_neighbor.distribute_local_to_global(
            dst, 0, cluster_manager.get_phi_neighbor_to_dst(face));
        phi_neighbor.distribute_local_to_global(
            flux_memory, 0, cluster_manager.get_phi_neighbor_to_fluxmemory(face));
      }
    }
  }

  template <int dim, int fe_degree, typename Number>
  void
  WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree, Number>::local_apply_ader_boundary_face(
      const MatrixFree<dim, value_type> &,
      std::vector<parallel::distributed::Vector<value_type>> &dst,
      const std::vector<parallel::distributed::Vector<value_type>> &src,
      const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, value_type> phi(
        this->data, true, 0, 0, 0, true);

    for (unsigned int face = face_range.first; face < face_range.second; face++)
    {
      if (cluster_manager.is_evaluate_face(face))
      {
        this->evaluate_boundary_face(phi, src, face, -1.0);

        // write only to the cell who is active
        phi.distribute_local_to_global(dst, 0, cluster_manager.get_phi_to_dst(face));
      }
    }
  }

  template <int dim, int fe_degree, typename Number>
  void WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree, Number>::
      evaluate_cells_and_faces_first_ader(
          const std::vector<parallel::distributed::Vector<value_type>> &src,
          std::vector<parallel::distributed::Vector<value_type>> &dst) const
  {
    for (unsigned int d = 0; d < dim + 1; ++d) this->tempsrc[d] = 0.;

    this->data.cell_loop(&WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree,
                             Number>::local_apply_firstader_domain,
        this, this->tempsrc, src);

    // evaluate faces for these cells
    this->data.loop(
        &WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree, Number>::local_apply_dummy_domain,
        &WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree, Number>::local_apply_ader_face,
        &WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree,
            Number>::local_apply_ader_boundary_face,
        this, dst, this->tempsrc);
  }


  template <int dim, int fe_degree, typename Number>
  void WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree, Number>::evaluate_cells_second_ader(
      const std::vector<parallel::distributed::Vector<value_type>> &src,
      std::vector<parallel::distributed::Vector<value_type>> &dst) const
  {
    this->data.cell_loop(&WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree,
                             Number>::local_apply_secondader_domain,
        this, dst, src);
  }

  template <int dim, int fe_degree, typename Number>
  void WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree, Number>::apply_ader(
      const std::vector<parallel::distributed::Vector<value_type>> &src,
      std::vector<parallel::distributed::Vector<value_type>> &dst, const double &cur_time,
      const double &dt) const
  {
    this->time = cur_time + this->time_step;
    cluster_manager.perform_time_step(*this, src, dst);
  }

  template <int dim, int fe_degree, typename Number>
  void WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree, Number>::reconstruct_div_grad(
      const std::vector<parallel::distributed::Vector<value_type>> &src,
      std::vector<parallel::distributed::Vector<value_type>> &dst) const
  {
    this->data.loop(&WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree,
                        Number>::local_apply_postprocessing_domain,
        &WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree,
            Number>::local_apply_postprocessing_face,
        &WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree,
            Number>::local_apply_postprocessing_boundary_face,
        this, dst, src);

    this->data.cell_loop(&WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree,
                             Number>::local_apply_postprocessing_mass_matrix,
        this, dst, dst);
  }

  template <int dim, int fe_degree, typename Number>
  void WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree,
      Number>::local_apply_postprocessing_domain(const MatrixFree<dim, value_type> &,
      std::vector<parallel::distributed::Vector<value_type>> &dst,
      const std::vector<parallel::distributed::Vector<value_type>> &src,
      const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, dim, value_type> velocity(this->data);
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, value_type> pressure(this->data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      if (cluster_manager.is_evaluate_cell(cell))
      {
        this->evaluate_cell(velocity, pressure, src, cell);
        velocity.set_dof_values(dst, 0);
        pressure.set_dof_values(dst, dim);
      }
    }
  }

  template <int dim, int fe_degree, typename Number>
  void
  WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree, Number>::local_apply_postprocessing_face(
      const MatrixFree<dim, value_type> &,
      std::vector<parallel::distributed::Vector<value_type>> &dst,
      const std::vector<parallel::distributed::Vector<value_type>> &src,
      const std::pair<unsigned int, unsigned int> &face_range) const
  {
    // There is some overhead in the methods in FEEvaluation, so it is faster
    // to combine pressure and velocity in the same object and just combine
    // them at the level of quadrature points
    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, value_type> phi(
        this->data, true, 0, 0, 0, true);
    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, value_type> phi_neighbor(
        this->data, false, 0, 0, 0, true);

    for (unsigned int face = face_range.first; face < face_range.second; face++)
    {
      if (cluster_manager.is_evaluate_face(face))
      {
        this->evaluate_inner_face(phi, phi_neighbor, src, face, 1.0);
        phi.distribute_local_to_global(dst, 0, cluster_manager.get_phi_to_dst(face));
        phi_neighbor.distribute_local_to_global(
            dst, 0, cluster_manager.get_phi_neighbor_to_dst(face));
      }
    }
  }

  template <int dim, int fe_degree, typename Number>
  void WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree,
      Number>::local_apply_postprocessing_boundary_face(const MatrixFree<dim, value_type> &,
      std::vector<parallel::distributed::Vector<value_type>> &dst,
      const std::vector<parallel::distributed::Vector<value_type>> &src,
      const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, value_type> phi(
        this->data, true, 0, 0, 0, true);
    for (unsigned int face = face_range.first; face < face_range.second; face++)
    {
      if (cluster_manager.is_evaluate_face(face))
      {
        this->evaluate_boundary_face(phi, src, face, 1.0);
        phi.distribute_local_to_global(dst, 0, cluster_manager.get_phi_to_dst(face));
      }
    }
  }

  template <int dim, int fe_degree, typename Number>
  void WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree,
      Number>::local_apply_postprocessing_mass_matrix(const MatrixFree<dim, value_type> &,
      std::vector<parallel::distributed::Vector<value_type>> &dst,
      const std::vector<parallel::distributed::Vector<value_type>> &src,
      const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      if (cluster_manager.is_evaluate_cell(cell))
      {
        this->mass_matrix_data->phi[0].reinit(cell);
        this->mass_matrix_data->phi[0].read_dof_values(src, 0);

        this->mass_matrix_data->inverse.fill_inverse_JxW_values(
            this->mass_matrix_data->coefficients);
        this->mass_matrix_data->inverse.apply(this->mass_matrix_data->coefficients, dim + 1,
            this->mass_matrix_data->phi[0].begin_dof_values(),
            this->mass_matrix_data->phi[0].begin_dof_values());

        this->mass_matrix_data->phi[0].set_dof_values(dst, 0);
      }
    }
  }


  // explicit instantiations
  template class WaveEquationOperationAcousticWaveADERLTS<2, 1, double>;
  template class WaveEquationOperationAcousticWaveADERLTS<2, 2, double>;
  template class WaveEquationOperationAcousticWaveADERLTS<2, 3, double>;
  template class WaveEquationOperationAcousticWaveADERLTS<2, 4, double>;
  template class WaveEquationOperationAcousticWaveADERLTS<2, 5, double>;
  template class WaveEquationOperationAcousticWaveADERLTS<2, 6, double>;
  template class WaveEquationOperationAcousticWaveADERLTS<3, 1, double>;
  template class WaveEquationOperationAcousticWaveADERLTS<3, 2, double>;
  template class WaveEquationOperationAcousticWaveADERLTS<3, 3, double>;
  template class WaveEquationOperationAcousticWaveADERLTS<3, 4, double>;
  template class WaveEquationOperationAcousticWaveADERLTS<3, 5, double>;
  template class WaveEquationOperationAcousticWaveADERLTS<3, 6, double>;
  template class WaveEquationOperationAcousticWaveADERLTS<2, 1, float>;
  template class WaveEquationOperationAcousticWaveADERLTS<2, 2, float>;
  template class WaveEquationOperationAcousticWaveADERLTS<2, 3, float>;
  template class WaveEquationOperationAcousticWaveADERLTS<2, 4, float>;
  template class WaveEquationOperationAcousticWaveADERLTS<2, 5, float>;
  template class WaveEquationOperationAcousticWaveADERLTS<2, 6, float>;
  template class WaveEquationOperationAcousticWaveADERLTS<3, 1, float>;
  template class WaveEquationOperationAcousticWaveADERLTS<3, 2, float>;
  template class WaveEquationOperationAcousticWaveADERLTS<3, 3, float>;
  template class WaveEquationOperationAcousticWaveADERLTS<3, 4, float>;
  template class WaveEquationOperationAcousticWaveADERLTS<3, 5, float>;
  template class WaveEquationOperationAcousticWaveADERLTS<3, 6, float>;
}  // namespace ACOU


#endif  // HAVE_DEAL_II
