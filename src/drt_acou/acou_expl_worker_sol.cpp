/*!----------------------------------------------------------------------
\file acou_expl_worker_sol.cpp
\brief Control routine for acoustic explicit time integration for solids
       or aka elastodynamics
\level 2

<pre>
\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
*----------------------------------------------------------------------*/

#include "acou_expl_worker.H"

#ifdef HAVE_DEAL_II

#include <deal.II/lac/parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>
//#include "fe_evaluation.h"
#include <deal.II/matrix_free/operators.h>
#include <deal.II/base/conditional_ostream.h>
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
  template <int dim, int fe_degree, typename Number>
  WaveEquationOperationElasticWave<dim, fe_degree, Number>::WaveEquationOperationElasticWave(
      const std::vector<const DoFHandler<dim> *> &dof_handlers,
      Teuchos::RCP<DRT::DiscretizationHDG> &discret,
      Teuchos::RCP<Function<dim>> boundary_conditions, Teuchos::RCP<Function<dim>> source_term,
      value_type time_step_in, int sourceno, Teuchos::RCP<PATMonitorManager> monitormanagerin)
      : WaveEquationOperation<dim, fe_degree, Number>(dof_handlers, discret, boundary_conditions,
            source_term, time_step_in, sourceno, monitormanagerin)
  {
    this->viscs.resize(this->data.n_macro_cells() + this->data.n_macro_ghost_cells());

    for (unsigned int i = 0; i < this->data.n_macro_cells() + this->data.n_macro_ghost_cells(); ++i)
    {
      this->densities[i] = make_vectorized_array<value_type>(1.);
      this->speeds[i] = make_vectorized_array<value_type>(1.);
      this->viscs[i] = make_vectorized_array<value_type>(1.);
      for (unsigned int v = 0; v < this->data.n_components_filled(i); ++v)
      {
        if (this->data.get_cell_iterator(i, v)->level() != 0)
          dserror("Refined meshes currently not implemented!");

        const int element_index = this->data.get_cell_iterator(i, v)->index();
        Teuchos::RCP<MAT::Material> mat = discret->lColElement(element_index)->Material();

        MAT::AcousticSolMat *actmat = static_cast<MAT::AcousticSolMat *>(mat.get());
        this->densities[i][v] = actmat->Density();
        this->speeds[i][v] = actmat->SpeedofSound();
        this->viscs[i][v] = actmat->Viscosity();
      }
    }
  }

  template <int dim, int fe_degree, typename Number>
  void WaveEquationOperationElasticWave<dim, fe_degree, Number>::local_apply_domain(
      const MatrixFree<dim, value_type> &data,
      std::vector<parallel::distributed::Vector<value_type>> &dst,
      const std::vector<parallel::distributed::Vector<value_type>> &src,
      const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, dim, value_type> v(data);
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, value_type> p(data);
    FEEvaluation<dim, fe_degree, fe_degree + 1, dim * dim, value_type> H(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // It is faster to evaluate values of the vector-valued velocity and
      // gradients of the scalar pressure than divergence of velocity and
      // values of pressure
      v.reinit(cell);
      v.read_dof_values(src, 0);
      v.evaluate(true, true, false);

      p.reinit(cell);
      p.read_dof_values(src, dim);
      p.evaluate(true, false, false);

      H.reinit(cell);
      H.read_dof_values(src, dim + 1);
      H.evaluate(true, false, false);

      const VectorizedArray<value_type> rho = this->densities[cell];
      const VectorizedArray<value_type> rho_inv = 1. / this->densities[cell];
      const VectorizedArray<value_type> c_sq = this->speeds[cell] * this->speeds[cell];
      const VectorizedArray<value_type> visc = this->viscs[cell];

      for (unsigned int q = 0; q < v.n_q_points; ++q)
      {
        const VectorizedArray<value_type> p_value = p.get_value(q);
        const Tensor<1, dim, VectorizedArray<value_type>> v_value = v.get_value(q);
        const Tensor<2, dim, VectorizedArray<value_type>> v_gradient = v.get_gradient(q);
        const Tensor<1, dim * dim, VectorizedArray<value_type>> H_value = H.get_value(q);

        Point<dim, VectorizedArray<value_type>> q_points = v.quadrature_point(q);
        Tensor<1, dim, VectorizedArray<value_type>> rhs;
        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int n = 0; n < rhs[d].n_array_elements; ++n)
          {
            Point<dim> q_point;
            for (unsigned int d = 0; d < dim; ++d) q_point[d] = q_points[d][n];
            rhs[d][n] = this->source_term->value(q_point, d);
          }
        v.submit_value(rho_inv * rhs, q);

        Tensor<2, dim, VectorizedArray<value_type>> muHminp;
        for (unsigned int m = 0; m < dim; ++m)
          for (unsigned int n = 0; n < dim; ++n)
            if (m == n)
              muHminp[n][n] = visc * H_value[n + dim * m] + p_value;
            else
              muHminp[n][m] = visc * H_value[m + dim * n];

        // contribution to H
        Tensor<1, dim * dim, VectorizedArray<value_type>> help;
        for (unsigned int m = 0; m < dim; ++m)
          for (unsigned int n = 0; n < dim; ++n) help[m + dim * n] = v_gradient[n][m];

        /*Tensor<1,dim*dim,Tensor<1,dim,VectorizedArray<value_type> > > contribH;
        contribH[0][0]=v_value[0];
        contribH[1][1]=v_value[0];
        contribH[2][0]=v_value[1];
        contribH[3][1]=v_value[1];*/

        // note that, v_gradient is sorted as follows
        //
        //  [ v_grad[0][0]   v_grad[0][1] ]
        //  [ v_grad[1][0]   v_grad[1][1] ]
        //
        //       =
        //
        //  [ v_x,x   v_x,y ]
        //  [ v_y,x   v_y,y ]


        if (this->adjoint_eval == false)
        {
          // H.submit_gradient(-contribH,q);
          H.submit_value(help, q);
          v.submit_gradient(-rho_inv * muHminp, q);
          p.submit_gradient(-rho * c_sq * v_value, q);
        }
        else
        {
          dserror("think about it");
        }
      }

      v.integrate(true, true);
      v.distribute_local_to_global(dst, 0);

      p.integrate(false, true);
      p.distribute_local_to_global(dst, dim);

      H.integrate(true, false);
      H.distribute_local_to_global(dst, dim + 1);
    }
  }



  template <int dim, int fe_degree, typename Number>
  void WaveEquationOperationElasticWave<dim, fe_degree, Number>::local_apply_face(
      const MatrixFree<dim, value_type> &,
      std::vector<parallel::distributed::Vector<value_type>> &dst,
      const std::vector<parallel::distributed::Vector<value_type>> &src,
      const std::pair<unsigned int, unsigned int> &face_range) const
  {
    // There is some overhead in the methods in FEEvaluation, so it is faster
    // to combine pressure and velocity in the same object and just combine
    // them at the level of quadrature points
    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, dim * dim + dim + 1, value_type> phi(
        this->data, true, 0, 0, 0, true);
    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, dim * dim + dim + 1, value_type> phi_neighbor(
        this->data, false, 0, 0, 0, true);

    for (unsigned int face = face_range.first; face < face_range.second; face++)
    {
      phi.reinit(face);
      phi.read_dof_values(src, 0);
      phi.evaluate(true, false);
      const VectorizedArray<value_type> rho_plus = phi.read_cell_data(this->densities);
      const VectorizedArray<value_type> rho_inv_plus = 1. / rho_plus;
      const VectorizedArray<value_type> c_plus = phi.read_cell_data(this->speeds);
      const VectorizedArray<value_type> c_sq_plus = c_plus * c_plus;
      const VectorizedArray<value_type> tau_plus = make_vectorized_array<value_type>(1.0);
      const VectorizedArray<value_type> visc_plus = phi.read_cell_data(viscs);

      phi_neighbor.reinit(face);
      phi_neighbor.read_dof_values(src, 0);
      phi_neighbor.evaluate(true, false);
      const VectorizedArray<value_type> rho_minus = phi_neighbor.read_cell_data(this->densities);
      const VectorizedArray<value_type> rho_inv_minus = 1. / rho_minus;
      const VectorizedArray<value_type> c_minus = phi_neighbor.read_cell_data(this->speeds);
      const VectorizedArray<value_type> c_sq_minus = c_minus * c_minus;
      const VectorizedArray<value_type> tau_minus = make_vectorized_array<value_type>(1.0);
      const VectorizedArray<value_type> visc_minus = phi_neighbor.read_cell_data(this->viscs);

      const VectorizedArray<value_type> tau_inv = 1. / (tau_plus + tau_minus);

      AssertDimension(phi.n_q_points, this->data.get_n_q_points_face(0));

      for (unsigned int q = 0; q < phi.n_q_points; ++q)
      {
        Tensor<1, dim * dim + dim + 1, VectorizedArray<value_type>> val_plus = phi.get_value(q);
        Tensor<1, dim * dim + dim + 1, VectorizedArray<value_type>> val_minus =
            phi_neighbor.get_value(q);
        Tensor<1, dim, VectorizedArray<value_type>> normal = phi.get_normal_vector(q);

        Tensor<1, dim, VectorizedArray<value_type>> vhat;
        for (unsigned int d = 0; d < dim; ++d)
        {
          vhat[d] = 1. / 2. * (val_plus[d] + val_minus[d]) -
                    tau_inv * (val_plus[dim] - val_minus[dim]) * normal[d];
          for (unsigned int e = 0; e < dim; ++e)
            vhat[d] -= tau_inv *
                       (visc_plus * val_plus[dim + 1 + d * dim + e] -
                           visc_minus * val_minus[dim + 1 + d * dim + e]) *
                       normal[e];
        }

        // vhat*normal+
        VectorizedArray<value_type> normal_vhat = vhat[0] * normal[0];
        for (unsigned int d = 1; d < dim; ++d) normal_vhat += vhat[d] * normal[d];

        Tensor<1, dim * dim + dim + 1, VectorizedArray<value_type>> submit_plus;
        Tensor<1, dim * dim + dim + 1, VectorizedArray<value_type>> submit_minus;

        // 1. contribution to pressure
        submit_plus[dim] = c_sq_plus * rho_plus * normal_vhat;
        submit_minus[dim] = -c_sq_minus * rho_minus * normal_vhat;

        // 2. contribution to velocity
        for (unsigned int d = 0; d < dim; ++d)
        {
          submit_plus[d] = -rho_inv_plus * tau_plus * (val_plus[d] - vhat[d]);
          submit_minus[d] = -rho_inv_minus * tau_minus * (val_minus[d] - vhat[d]);
          for (unsigned int e = 0; e < dim; ++e)
            if (e == d)
            {
              submit_plus[d] += rho_inv_plus *
                                (visc_plus * val_plus[dim + 1 + e * dim + d] + val_plus[dim]) *
                                normal[e];
              submit_minus[d] -= rho_inv_minus *
                                 (visc_minus * val_minus[dim + 1 + e * dim + d] + val_minus[dim]) *
                                 normal[e];
            }
            else
            {
              submit_plus[d] +=
                  rho_inv_plus * (visc_plus * val_plus[dim + 1 + d * dim + e]) * normal[e];
              submit_minus[d] -=
                  rho_inv_minus * (visc_minus * val_minus[dim + 1 + d * dim + e]) * normal[e];
            }
        }

        // 3. contribution to H
        for (unsigned int d = 0; d < dim; ++d)
          for (unsigned int e = 0; e < dim; ++e)
          {
            submit_plus[dim + 1 + d * dim + e] -= (val_plus[d] - vhat[d]) * normal[e];
            submit_minus[dim + 1 + d * dim + e] += (val_minus[d] - vhat[d]) * normal[e];
            // submit_plus[dim+1+d*dim+e] += vhat[d]*normal[e];
            // submit_minus[dim+1+d*dim+e] -= vhat[d]*normal[e];
          }


        if (this->adjoint_eval) dserror("think about adjoint problem");

        phi.submit_value(submit_plus, q);
        phi_neighbor.submit_value(submit_minus, q);
      }

      phi.integrate(true, false);
      phi.distribute_local_to_global(dst, 0);

      phi_neighbor.integrate(true, false);
      phi_neighbor.distribute_local_to_global(dst, 0);
    }
  }

  template <int dim, int fe_degree, typename Number>
  void WaveEquationOperationElasticWave<dim, fe_degree, Number>::local_apply_boundary_face(
      const MatrixFree<dim, value_type> &,
      std::vector<parallel::distributed::Vector<value_type>> &dst,
      const std::vector<parallel::distributed::Vector<value_type>> &src,
      const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, dim * dim + dim + 1, value_type> phi(
        this->data, true, 0, 0, 0, true);

    // quantities we need in the loop
    Point<dim> point;
    std::vector<value_type> node_values;
    std::vector<std::vector<value_type>> node_coords;
    node_coords.resize(GeometryInfo<dim>::vertices_per_face);
    node_values.resize(GeometryInfo<dim>::vertices_per_face);
    for (unsigned int n = 0; n < GeometryInfo<dim>::vertices_per_face; ++n)
      node_coords[n].resize(dim);

    for (unsigned int face = face_range.first; face < face_range.second; face++)
    {
      phi.reinit(face);
      phi.read_dof_values(src, 0);
      phi.evaluate(true, false);
      const VectorizedArray<value_type> rho = phi.read_cell_data(this->densities);
      const VectorizedArray<value_type> rho_inv = 1. / rho;
      const VectorizedArray<value_type> c_sq =
          phi.read_cell_data(this->speeds) * phi.read_cell_data(this->speeds);
      const VectorizedArray<value_type> c = phi.read_cell_data(this->speeds);
      const VectorizedArray<value_type> tau = make_vectorized_array<value_type>(1.0);
      const VectorizedArray<value_type> visc = phi.read_cell_data(this->viscs);

      const types::boundary_id boundary_index = this->data.get_boundary_indicator(face);
      const int int_boundary_id = int(boundary_index);

      for (unsigned int q = 0; q < phi.n_q_points; ++q)
      {
        Tensor<1, dim, VectorizedArray<value_type>> normal = phi.get_normal_vector(q);
        Tensor<1, dim * dim + dim + 1, VectorizedArray<value_type>> val_plus = phi.get_value(q);
        Point<dim, VectorizedArray<value_type>> q_point = phi.quadrature_point(q);

        // calculation of vhat dependent on boundary type
        Tensor<1, dim, VectorizedArray<value_type>> vhat;
        if (int_boundary_id == 0)  // absorbing boundary
        {
          for (unsigned int d = 0; d < dim; ++d)
          {
            vhat[d] = tau / (rho * c + tau) * val_plus[d] -
                      1. / (tau + rho * c) * val_plus[dim] * normal[d];
            for (unsigned int e = 0; e < dim; ++e)
              vhat[d] -= 1. / (tau + rho * c) * visc * val_plus[dim + 1 + d * dim + e] * normal[e];
          }
        }
        else if (int_boundary_id == 1)  // monitored
        {
          if (this->adjoint_eval == false)
            for (unsigned int d = 0; d < dim; ++d)
            {
              vhat[d] = val_plus[d] - 1. / tau * val_plus[dim] * normal[d];
              for (unsigned int e = 0; e < dim; ++e)
                vhat[d] -= 1. / tau * visc * val_plus[dim + 1 + d * dim + e] * normal[e];
            }
          else
            dserror("todo");
        }
        else if (int_boundary_id == 2)  // monitored and absorbing
        {
          if (this->adjoint_eval == false)
            for (unsigned int d = 0; d < dim; ++d)
            {
              vhat[d] = tau / (rho * c + tau) * val_plus[d] -
                        1. / (tau + rho * c) * val_plus[dim] * normal[d];
              for (unsigned int e = 0; e < dim; ++e)
                vhat[d] -=
                    1. / (tau + rho * c) * visc * val_plus[dim + 1 + d * dim + e] * normal[e];
            }
          else
            dserror("todo");
        }
        else if (int_boundary_id == 3)  // free boundary
          for (unsigned int d = 0; d < dim; ++d)
          {
            vhat[d] = val_plus[d] - 1. / tau * val_plus[dim] * normal[d];
            for (unsigned int e = 0; e < dim; ++e)
              vhat[d] -= 1. / tau * visc * val_plus[dim + 1 + d * dim + e] * normal[e];
          }
        else if (int_boundary_id == 4)  // dbc from time reversal
        {
          dserror("todo");
        }
        else if (int_boundary_id >= 5)  // dbcs
        {
          if (this->adjoint_eval == false)
            for (unsigned int v = 0; v < VectorizedArray<value_type>::n_array_elements; ++v)
            {
              Point<dim> point;
              for (unsigned int d = 0; d < dim; ++d) point[d] = q_point[d][v];
              for (unsigned int d = 0; d < dim; ++d)
                vhat[d][v] = this->dirichlet_boundary_conditions->value(
                    point, (int_boundary_id - 5) * dim + d);
            }
          else
            for (unsigned int d = 0; d < dim; ++d) vhat[d] = VectorizedArray<value_type>();
        }

        // now the actual boundary term evaluation
        Tensor<1, dim * dim + dim + 1, VectorizedArray<value_type>> submit_plus;
        if (this->adjoint_eval == false)
        {
          // vhat*normal
          VectorizedArray<value_type> normal_vhat = vhat[0] * normal[0];
          for (unsigned int d = 1; d < dim; ++d) normal_vhat += vhat[d] * normal[d];

          // 1. contribution to pressure
          submit_plus[dim] = c_sq * rho * normal_vhat;

          // 2. contribution to velocity
          for (unsigned int d = 0; d < dim; ++d)
          {
            submit_plus[d] = -rho_inv * tau * (val_plus[d] - vhat[d]);
            for (unsigned int e = 0; e < dim; ++e)
              if (e == d)
                submit_plus[d] +=
                    rho_inv * (visc * val_plus[dim + 1 + e * dim + d] + val_plus[dim]) * normal[e];
              else
                submit_plus[d] += rho_inv * (visc * val_plus[dim + 1 + d * dim + e]) * normal[e];
          }

          // 3. contribution to H
          for (unsigned int d = 0; d < dim; ++d)
            for (unsigned int e = 0; e < dim; ++e)
            {
              // submit_plus[dim+1+d*dim+e] += vhat[d]*normal[e];
              submit_plus[dim + 1 + d * dim + e] -= (val_plus[d] - vhat[d]) * normal[e];
            }
        }
        else
        {
          dserror("todo");
        }

        phi.submit_value(submit_plus, q);
      }
      phi.integrate(true, false);
      phi.distribute_local_to_global(dst, 0);
    }
  }

  template <int dim, int fe_degree, typename Number>
  void WaveEquationOperationElasticWave<dim, fe_degree, Number>::local_apply_mass_matrix(
      const MatrixFree<dim, value_type> &data,
      std::vector<parallel::distributed::Vector<value_type>> &dst,
      const std::vector<parallel::distributed::Vector<value_type>> &src,
      const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      this->mass_matrix_data_solid->phi[0].reinit(cell);
      this->mass_matrix_data_solid->phi[0].read_dof_values(src, 0);

      this->mass_matrix_data_solid->inverse.fill_inverse_JxW_values(
          this->mass_matrix_data_solid->coefficients);
      this->mass_matrix_data_solid->inverse.apply(this->mass_matrix_data_solid->coefficients,
          dim + 1, this->mass_matrix_data_solid->phi[0].begin_dof_values(),
          this->mass_matrix_data_solid->phi[0].begin_dof_values());

      this->mass_matrix_data_solid->phi[0].set_dof_values(dst, 0);
    }
  }

  template <int dim, int fe_degree, typename Number>
  void WaveEquationOperationElasticWave<dim, fe_degree, Number>::write_deal_cell_values(
      Teuchos::RCP<DRT::DiscretizationHDG> &discret,
      const std::vector<parallel::distributed::Vector<value_type>> &src) const
  {
    const unsigned dofs_per_cell = this->data.get_dof_handler().get_fe().dofs_per_cell;
    std::vector<types::global_dof_index> indices, local_dof_indices(dofs_per_cell);
    for (typename DoFHandler<dim>::active_cell_iterator cell =
             this->data.get_dof_handler().begin_active();
         cell != this->data.get_dof_handler().end(); ++cell)
    {
      cell->get_dof_indices(local_dof_indices);
      indices.insert(indices.end(), local_dof_indices.begin(), local_dof_indices.end());
    }

    IndexSet relevant_dofs(src[0].size());
    relevant_dofs.add_indices(indices.begin(), indices.end());
    relevant_dofs.compress();
    std::vector<parallel::distributed::Vector<value_type>> ghosted_vector(src.size());
    for (unsigned int i = 0; i < src.size(); ++i)
    {
      ghosted_vector[i].reinit(this->data.get_dof_handler().locally_owned_dofs(), relevant_dofs,
          src[0].get_mpi_communicator());
      ghosted_vector[i] = src[i];
      ghosted_vector[i].update_ghost_values();
    }

    unsigned int ndofs1d;
    if (dim == 2)
      ndofs1d = std::sqrt(dofs_per_cell);
    else if (dim == 3)
      ndofs1d = int(std::pow(dofs_per_cell, 1.0 / 3.0));
    unsigned int ndofs2d = ndofs1d * ndofs1d;

    unsigned int nodes_per_cell = GeometryInfo<dim>::vertices_per_cell;
    std::vector<Point<dim>> baci_vals_loc(nodes_per_cell);
    std::vector<Point<dim>> deal_vals_loc(nodes_per_cell);

    Vector<value_type> local_values(dofs_per_cell);
    for (int i = 0; i < discret->NumMyColElements(); ++i)
    {
      typename DoFHandler<dim>::active_cell_iterator cell(
          &this->data.get_dof_handler().get_triangulation(), 0, i, &this->data.get_dof_handler());
      DRT::ELEMENTS::AcouSol *acouele =
          dynamic_cast<DRT::ELEMENTS::AcouSol *>(discret->lColElement(i));

      for (unsigned int n = 0; n < nodes_per_cell; ++n)
      {
        for (int d = 0; d < dim; ++d)
        {
          deal_vals_loc[n](d) = cell->vertex(n)(d);
          baci_vals_loc[n](d) = acouele->Nodes()[n]->X()[d];
        }
      }

      // perform permutation: step 2: swap it
      switch (acouele->Shape())
      {
        case DRT::Element::quad4:
        {
          if ((deal_vals_loc[0].distance(baci_vals_loc[0]) < 1e-10 &&
                  deal_vals_loc[1].distance(baci_vals_loc[1]) < 1e-10 &&
                  deal_vals_loc[2].distance(baci_vals_loc[3]) < 1e-10 &&
                  deal_vals_loc[3].distance(baci_vals_loc[2]) < 1e-10))
          {
            // everything is alright
            for (unsigned int d = 0; d < dim; ++d)
            {
              cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                acouele->eleinteriorVelnp_(d * dofs_per_cell + j) = local_values[j];
            }
            cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              acouele->eleinteriorPressnp_(j) = local_values[j];
            // stresses
            for (unsigned d = 0; d < dim; ++d)
              for (unsigned e = 0; e < dim; ++e)
              {
                cell->get_interpolated_dof_values(
                    ghosted_vector[dim + 1 + d * dim + e], local_values);
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  acouele->eleinteriorGradVelnp_((d * dim + e) * dofs_per_cell + j) =
                      local_values[j];
              }
          }
          else if (deal_vals_loc[0].distance(baci_vals_loc[3]) < 1e-10 &&
                   deal_vals_loc[1].distance(baci_vals_loc[0]) < 1e-10 &&
                   deal_vals_loc[2].distance(baci_vals_loc[2]) < 1e-10 &&
                   deal_vals_loc[3].distance(baci_vals_loc[1]) < 1e-10)
          {
            for (unsigned int d = 0; d < dim; ++d)
            {
              cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                const int ax = j % ndofs1d;
                const int ay = j / ndofs1d;
                int permute = (ndofs1d - 1 - ax) * ndofs1d + ay;
                acouele->eleinteriorVelnp_(d * dofs_per_cell + permute) = local_values[j];
              }
            }
            cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              const int ax = j % ndofs1d;
              const int ay = j / ndofs1d;
              int permute = (ndofs1d - 1 - ax) * ndofs1d + ay;
              acouele->eleinteriorPressnp_(permute) = local_values[j];
            }

            for (unsigned int d = 0; d < dim; ++d)
              for (unsigned e = 0; e < dim; ++e)
              {
                cell->get_interpolated_dof_values(
                    ghosted_vector[dim + 1 + d * dim + e], local_values);
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  const int ax = j % ndofs1d;
                  const int ay = j / ndofs1d;
                  int permute = (ndofs1d - 1 - ax) * ndofs1d + ay;
                  acouele->eleinteriorGradVelnp_((d * dim + e) * dofs_per_cell + permute) =
                      local_values[j];
                }
              }
          }
          else if (deal_vals_loc[0].distance(baci_vals_loc[2]) < 1e-10 &&
                   deal_vals_loc[1].distance(baci_vals_loc[3]) < 1e-10 &&
                   deal_vals_loc[2].distance(baci_vals_loc[1]) < 1e-10 &&
                   deal_vals_loc[3].distance(baci_vals_loc[0]) < 1e-10)
          {
            for (unsigned int d = 0; d < dim; ++d)
            {
              cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                const int ax = j % ndofs1d;
                const int ay = j / ndofs1d;
                int permute = (ndofs1d - 1 - ax) + (ndofs1d - 1 - ay) * ndofs1d;
                acouele->eleinteriorVelnp_(d * dofs_per_cell + permute) = local_values[j];
              }
            }
            cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              const int ax = j % ndofs1d;
              const int ay = j / ndofs1d;
              int permute = (ndofs1d - 1 - ax) + (ndofs1d - 1 - ay) * ndofs1d;
              acouele->eleinteriorPressnp_(permute) = local_values[j];
            }

            for (unsigned int d = 0; d < dim; ++d)
              for (unsigned e = 0; e < dim; ++e)
              {
                cell->get_interpolated_dof_values(
                    ghosted_vector[dim + 1 + d * dim + e], local_values);
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  const int ax = j % ndofs1d;
                  const int ay = j / ndofs1d;
                  int permute = (ndofs1d - 1 - ax) + (ndofs1d - 1 - ay) * ndofs1d;
                  acouele->eleinteriorGradVelnp_((d * dim + e) * dofs_per_cell + permute) =
                      local_values[j];
                }
              }
          }
          else if (deal_vals_loc[0].distance(baci_vals_loc[1]) < 1e-10 &&
                   deal_vals_loc[1].distance(baci_vals_loc[2]) < 1e-10 &&
                   deal_vals_loc[2].distance(baci_vals_loc[0]) < 1e-10 &&
                   deal_vals_loc[3].distance(baci_vals_loc[3]) < 1e-10)
          {
            for (unsigned int d = 0; d < dim; ++d)
            {
              cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                const int ax = j % ndofs1d;
                const int ay = j / ndofs1d;
                int permute = (ax)*ndofs1d + (ndofs1d - 1 - ay);
                acouele->eleinteriorVelnp_(d * dofs_per_cell + permute) = local_values[j];
              }
            }
            cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              const int ax = j % ndofs1d;
              const int ay = j / ndofs1d;
              int permute = (ax)*ndofs1d + (ndofs1d - 1 - ay);
              acouele->eleinteriorPressnp_(permute) = local_values[j];
            }

            for (unsigned int d = 0; d < dim; ++d)
              for (unsigned e = 0; e < dim; ++e)
              {
                cell->get_interpolated_dof_values(
                    ghosted_vector[dim + 1 + d * dim + e], local_values);
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                  const int ax = j % ndofs1d;
                  const int ay = j / ndofs1d;
                  int permute = (ax)*ndofs1d + (ndofs1d - 1 - ay);
                  acouele->eleinteriorGradVelnp_((d * dim + e) * dofs_per_cell + permute) =
                      local_values[j];
                }
              }
          }
          else
            dserror("unknown permutation");
          break;
        }
        case DRT::Element::hex8:
        {
          if (deal_vals_loc[0].distance(baci_vals_loc[0]) < 1e-10 &&
              deal_vals_loc[1].distance(baci_vals_loc[1]) < 1e-10 &&
              deal_vals_loc[2].distance(baci_vals_loc[3]) < 1e-10 &&
              deal_vals_loc[3].distance(baci_vals_loc[2]) < 1e-10 &&
              deal_vals_loc[4].distance(baci_vals_loc[4]) < 1e-10 &&
              deal_vals_loc[5].distance(baci_vals_loc[5]) < 1e-10 &&
              deal_vals_loc[6].distance(baci_vals_loc[7]) < 1e-10 &&
              deal_vals_loc[7].distance(baci_vals_loc[6]) < 1e-10)
          {
            // everything is alright
            for (unsigned int d = 0; d < dim; ++d)
            {
              cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                acouele->eleinteriorVelnp_(d * dofs_per_cell + j) = local_values[j];
            }
            cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              acouele->eleinteriorPressnp_(j) = local_values[j];
          }
          else if (deal_vals_loc[0].distance(baci_vals_loc[4]) < 1e-10 &&
                   deal_vals_loc[1].distance(baci_vals_loc[5]) < 1e-10 &&
                   deal_vals_loc[2].distance(baci_vals_loc[0]) < 1e-10 &&
                   deal_vals_loc[3].distance(baci_vals_loc[1]) < 1e-10 &&
                   deal_vals_loc[4].distance(baci_vals_loc[7]) < 1e-10 &&
                   deal_vals_loc[5].distance(baci_vals_loc[6]) < 1e-10 &&
                   deal_vals_loc[6].distance(baci_vals_loc[3]) < 1e-10 &&
                   deal_vals_loc[7].distance(baci_vals_loc[2]) < 1e-10)
          {
            for (unsigned int d = 0; d < dim; ++d)
            {
              cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                const int ax = j % ndofs1d;
                const int ay = int(j / ndofs1d) % ndofs1d;
                const int az = j / ndofs2d;
                int permute = ax + (ndofs1d - 1 - ay) * ndofs2d + az * ndofs1d;
                acouele->eleinteriorVelnp_(d * dofs_per_cell + permute) = local_values[j];
              }
            }
            cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              const int ax = j % ndofs1d;
              const int ay = int(j / ndofs1d) % ndofs1d;
              const int az = j / ndofs2d;
              int permute = ax + (ndofs1d - 1 - ay) * ndofs2d + az * ndofs1d;
              acouele->eleinteriorPressnp_(permute) = local_values[j];
            }
          }
          else if (deal_vals_loc[0].distance(baci_vals_loc[5]) < 1e-10 &&
                   deal_vals_loc[1].distance(baci_vals_loc[6]) < 1e-10 &&
                   deal_vals_loc[2].distance(baci_vals_loc[1]) < 1e-10 &&
                   deal_vals_loc[3].distance(baci_vals_loc[2]) < 1e-10 &&
                   deal_vals_loc[4].distance(baci_vals_loc[4]) < 1e-10 &&
                   deal_vals_loc[5].distance(baci_vals_loc[7]) < 1e-10 &&
                   deal_vals_loc[6].distance(baci_vals_loc[0]) < 1e-10 &&
                   deal_vals_loc[7].distance(baci_vals_loc[3]) < 1e-10)
          {
            for (unsigned int d = 0; d < dim; ++d)
            {
              cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                const int ax = j % ndofs1d;
                const int ay = int(j / ndofs1d) % ndofs1d;
                const int az = j / ndofs2d;
                int permute = ax * ndofs1d + (ndofs1d - 1 - ay) * ndofs2d + (ndofs1d - 1 - az);
                acouele->eleinteriorVelnp_(d * dofs_per_cell + permute) = local_values[j];
              }
            }
            cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              const int ax = j % ndofs1d;
              const int ay = int(j / ndofs1d) % ndofs1d;
              const int az = j / ndofs2d;
              int permute = ax * ndofs1d + (ndofs1d - 1 - ay) * ndofs2d + (ndofs1d - 1 - az);
              acouele->eleinteriorPressnp_(permute) = local_values[j];
            }
          }
          else if (deal_vals_loc[0].distance(baci_vals_loc[7]) < 1e-10 &&
                   deal_vals_loc[1].distance(baci_vals_loc[4]) < 1e-10 &&
                   deal_vals_loc[2].distance(baci_vals_loc[3]) < 1e-10 &&
                   deal_vals_loc[3].distance(baci_vals_loc[0]) < 1e-10 &&
                   deal_vals_loc[4].distance(baci_vals_loc[6]) < 1e-10 &&
                   deal_vals_loc[5].distance(baci_vals_loc[5]) < 1e-10 &&
                   deal_vals_loc[6].distance(baci_vals_loc[2]) < 1e-10 &&
                   deal_vals_loc[7].distance(baci_vals_loc[1]) < 1e-10)
          {
            for (unsigned int d = 0; d < dim; ++d)
            {
              cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                const int ax = j % ndofs1d;
                const int ay = int(j / ndofs1d) % ndofs1d;
                const int az = j / ndofs2d;
                int permute = (ndofs1d - 1 - ax) * ndofs1d + (ndofs1d - 1 - ay) * ndofs2d + az;
                acouele->eleinteriorVelnp_(d * dofs_per_cell + permute) = local_values[j];
              }
            }
            cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              const int ax = j % ndofs1d;
              const int ay = int(j / ndofs1d) % ndofs1d;
              const int az = j / ndofs2d;
              int permute = (ndofs1d - 1 - ax) * ndofs1d + (ndofs1d - 1 - ay) * ndofs2d + az;
              acouele->eleinteriorPressnp_(permute) = local_values[j];
            }
          }
          else if (deal_vals_loc[0].distance(baci_vals_loc[6]) < 1e-10 &&
                   deal_vals_loc[1].distance(baci_vals_loc[7]) < 1e-10 &&
                   deal_vals_loc[2].distance(baci_vals_loc[2]) < 1e-10 &&
                   deal_vals_loc[3].distance(baci_vals_loc[3]) < 1e-10 &&
                   deal_vals_loc[4].distance(baci_vals_loc[5]) < 1e-10 &&
                   deal_vals_loc[5].distance(baci_vals_loc[4]) < 1e-10 &&
                   deal_vals_loc[6].distance(baci_vals_loc[1]) < 1e-10 &&
                   deal_vals_loc[7].distance(baci_vals_loc[0]) < 1e-10)
          {
            for (unsigned int d = 0; d < dim; ++d)
            {
              cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {
                const int ax = j % ndofs1d;
                const int ay = int(j / ndofs1d) % ndofs1d;
                const int az = j / ndofs2d;
                int permute = (ndofs1d - 1 - ax) + (ndofs1d - 1 - ay) * ndofs2d +
                              (ndofs1d - 1 - az) * ndofs1d;
                acouele->eleinteriorVelnp_(d * dofs_per_cell + permute) = local_values[j];
              }
            }
            cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              const int ax = j % ndofs1d;
              const int ay = int(j / ndofs1d) % ndofs1d;
              const int az = j / ndofs2d;
              int permute =
                  (ndofs1d - 1 - ax) + (ndofs1d - 1 - ay) * ndofs2d + (ndofs1d - 1 - az) * ndofs1d;
              acouele->eleinteriorPressnp_(permute) = local_values[j];
            }
          }
          else
            dserror("unknown permutation");
          break;
        }
        default:
          dserror("other distypes not yet implemented!");
          break;
      }
    }
    return;
  }

  template <int dim, int fe_degree, typename Number>
  void WaveEquationOperationElasticWave<dim, fe_degree, Number>::read_initial_conditions(
      Teuchos::RCP<DRT::DiscretizationHDG> &discret,
      std::vector<parallel::distributed::Vector<value_type>> &dst)
  {
    FEEvaluation<dim, fe_degree, fe_degree + 1, dim + 1 + dim * dim, value_type> phi(this->data);
    for (unsigned int j = 0; j < phi.dofs_per_cell; ++j)
      phi.submit_dof_value(Tensor<1, dim + 1 + dim * dim, VectorizedArray<value_type>>(), j);

    unsigned int dofs_per_cell = phi.dofs_per_cell;  // i assume, that it is the same for all cells
    unsigned int ndofs1d;
    if (dim == 2)
      ndofs1d = std::sqrt(dofs_per_cell);
    else if (dim == 3)
      ndofs1d = int(std::pow(dofs_per_cell, 1.0 / 3.0));
    // unsigned int ndofs2d = ndofs1d * ndofs1d;

    unsigned int nodes_per_cell = GeometryInfo<dim>::vertices_per_cell;
    std::vector<Point<dim>> baci_vals_loc(nodes_per_cell);
    std::vector<Point<dim>> deal_vals_loc(nodes_per_cell);

    for (unsigned int i = 0; i < this->data.n_macro_cells(); ++i)
    {
      phi.reinit(i);
      for (unsigned int v = 0; v < this->data.n_components_filled(i); ++v)
      {
        const int element_index = this->data.get_cell_iterator(i, v)->index();
        DRT::ELEMENTS::AcouSol *acouele =
            dynamic_cast<DRT::ELEMENTS::AcouSol *>(discret->lColElement(element_index));
        if (acouele == NULL) dserror("No acoustic solid element given!");

        // perform permutation: step 1: get the node coordinates
        for (unsigned int n = 0; n < nodes_per_cell; ++n)
        {
          for (int d = 0; d < dim; ++d)
          {
            deal_vals_loc[n](d) = this->data.get_cell_iterator(i, v)->vertex(n)(d);
            baci_vals_loc[n](d) = acouele->Nodes()[n]->X()[d];
          }
        }

        // perform permutation: step 2: swap it
        switch (acouele->Shape())
        {
          case DRT::Element::quad4:
          {
            if (deal_vals_loc[0].distance(baci_vals_loc[0]) < 1e-10 &&
                deal_vals_loc[1].distance(baci_vals_loc[1]) < 1e-10 &&
                deal_vals_loc[2].distance(baci_vals_loc[3]) < 1e-10 &&
                deal_vals_loc[3].distance(baci_vals_loc[2]) < 1e-10)
            {
              // everything is alright
              for (unsigned j = 0; j < dofs_per_cell; ++j)
              {
                for (unsigned int d = 0; d < dim; ++d)
                  phi.begin_dof_values()[d * dofs_per_cell + j][v] =
                      acouele->eleinteriorVelnp_(d * dofs_per_cell + j);
                for (unsigned int d = 0; d < dim; ++d)
                  for (unsigned int e = 0; e < dim; ++e)
                    phi.begin_dof_values()[(d * dim + e + dim + 1) * dofs_per_cell + j][v] =
                        acouele->eleinteriorGradVelnp_((d * dim + e) * dofs_per_cell + j);
                phi.begin_dof_values()[dim * dofs_per_cell + j][v] =
                    acouele->eleinteriorPressnp_(j);
              }
            }
            else if (deal_vals_loc[0].distance(baci_vals_loc[3]) < 1e-10 &&
                     deal_vals_loc[1].distance(baci_vals_loc[0]) < 1e-10 &&
                     deal_vals_loc[2].distance(baci_vals_loc[2]) < 1e-10 &&
                     deal_vals_loc[3].distance(baci_vals_loc[1]) < 1e-10)
            {
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                const int ax = i % ndofs1d;
                const int ay = i / ndofs1d;
                int permute = (ndofs1d - 1 - ax) * ndofs1d + ay;
                phi.begin_dof_values()[dim * dofs_per_cell + i][v] =
                    acouele->eleinteriorPressnp_(permute);
                for (unsigned int d = 0; d < dim; ++d)
                  phi.begin_dof_values()[d * dofs_per_cell + i][v] =
                      acouele->eleinteriorVelnp_(d * dofs_per_cell + permute);
                for (unsigned int d = 0; d < dim; ++d)
                  for (unsigned int e = 0; e < dim; ++e)
                    phi.begin_dof_values()[(d * dim + e + dim + 1) * dofs_per_cell + i][v] =
                        acouele->eleinteriorGradVelnp_((d * dim + e) * dofs_per_cell + permute);
              }
            }
            else if (deal_vals_loc[0].distance(baci_vals_loc[2]) < 1e-10 &&
                     deal_vals_loc[1].distance(baci_vals_loc[3]) < 1e-10 &&
                     deal_vals_loc[2].distance(baci_vals_loc[1]) < 1e-10 &&
                     deal_vals_loc[3].distance(baci_vals_loc[0]) < 1e-10)
            {
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                const int ax = i % ndofs1d;
                const int ay = i / ndofs1d;
                int permute = (ndofs1d - 1 - ax) + (ndofs1d - 1 - ay) * ndofs1d;
                phi.begin_dof_values()[dim * dofs_per_cell + i][v] =
                    acouele->eleinteriorPressnp_(permute);
                for (unsigned int d = 0; d < dim; ++d)
                  phi.begin_dof_values()[d * dofs_per_cell + i][v] =
                      acouele->eleinteriorVelnp_(d * dofs_per_cell + permute);
                for (unsigned int d = 0; d < dim; ++d)
                  for (unsigned int e = 0; e < dim; ++e)
                    phi.begin_dof_values()[(d * dim + e + dim + 1) * dofs_per_cell + i][v] =
                        acouele->eleinteriorGradVelnp_((d * dim + e) * dofs_per_cell + permute);
              }
            }
            else if (deal_vals_loc[0].distance(baci_vals_loc[1]) < 1e-10 &&
                     deal_vals_loc[1].distance(baci_vals_loc[2]) < 1e-10 &&
                     deal_vals_loc[2].distance(baci_vals_loc[0]) < 1e-10 &&
                     deal_vals_loc[3].distance(baci_vals_loc[3]) < 1e-10)
            {
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                const int ax = i % ndofs1d;
                const int ay = i / ndofs1d;
                int permute = (ax)*ndofs1d + (ndofs1d - 1 - ay);
                phi.begin_dof_values()[dim * dofs_per_cell + i][v] =
                    acouele->eleinteriorPressnp_(permute);
                for (unsigned int d = 0; d < dim; ++d)
                  phi.begin_dof_values()[d * dofs_per_cell + i][v] =
                      acouele->eleinteriorVelnp_(d * dofs_per_cell + permute);
                for (unsigned int d = 0; d < dim; ++d)
                  for (unsigned int e = 0; e < dim; ++e)
                    phi.begin_dof_values()[(d * dim + e + dim + 1) * dofs_per_cell + i][v] =
                        acouele->eleinteriorGradVelnp_((d * dim + e) * dofs_per_cell + permute);
              }
            }
            else
            {
              std::cout << "d " << deal_vals_loc[0](0) << " " << deal_vals_loc[0](1) << " b "
                        << baci_vals_loc[0](0) << " " << baci_vals_loc[0](1) << " " << std::endl;
              std::cout << "d " << deal_vals_loc[1](0) << " " << deal_vals_loc[1](1) << " b "
                        << baci_vals_loc[1](0) << " " << baci_vals_loc[1](1) << " " << std::endl;
              std::cout << "d " << deal_vals_loc[2](0) << " " << deal_vals_loc[2](1) << " b "
                        << baci_vals_loc[2](0) << " " << baci_vals_loc[2](1) << " " << std::endl;
              std::cout << "d " << deal_vals_loc[3](0) << " " << deal_vals_loc[3](1) << " b "
                        << baci_vals_loc[3](0) << " " << baci_vals_loc[3](1) << " " << std::endl;
              dserror("unknown permutation");
            }
            break;
          }
          default:
            dserror("other distypes not yet implemented!");
            break;
        }
      }
      phi.set_dof_values(dst);
    }

    for (unsigned int i = 0; i < dim + 1 + dim * dim; ++i) dst[i].update_ghost_values();
  }

  template <int dim, int fe_degree, typename Number>
  void WaveEquationOperationElasticWave<dim, fe_degree, Number>::apply(
      const std::vector<parallel::distributed::Vector<value_type>> &src,
      std::vector<parallel::distributed::Vector<value_type>> &dst, const double &cur_time,
      const double &dt) const
  {
    Timer timer;
    this->time = cur_time;
    this->dirichlet_boundary_conditions->set_time(this->time);
    this->source_term->set_time(this->time);

    this->data.loop(&WaveEquationOperationElasticWave<dim, fe_degree, Number>::local_apply_domain,
        &WaveEquationOperationElasticWave<dim, fe_degree, Number>::local_apply_face,
        &WaveEquationOperationElasticWave<dim, fe_degree, Number>::local_apply_boundary_face, this,
        dst, src);

    this->computing_times[0] += timer.wall_time();
    timer.restart();

    this->data.cell_loop(
        &WaveEquationOperationElasticWave<dim, fe_degree, Number>::local_apply_mass_matrix, this,
        dst, dst);

    this->computing_times[1] += timer.wall_time();
    this->computing_times[2] += 1.;
  }

  // explicit instantiations
  template class WaveEquationOperationElasticWave<2, 1, double>;
  template class WaveEquationOperationElasticWave<2, 2, double>;
  template class WaveEquationOperationElasticWave<2, 3, double>;
  template class WaveEquationOperationElasticWave<2, 4, double>;
  template class WaveEquationOperationElasticWave<2, 5, double>;
  template class WaveEquationOperationElasticWave<2, 6, double>;
  template class WaveEquationOperationElasticWave<3, 1, double>;
  template class WaveEquationOperationElasticWave<3, 2, double>;
  template class WaveEquationOperationElasticWave<3, 3, double>;
  template class WaveEquationOperationElasticWave<3, 4, double>;
  template class WaveEquationOperationElasticWave<3, 5, double>;
  template class WaveEquationOperationElasticWave<3, 6, double>;
  template class WaveEquationOperationElasticWave<2, 1, float>;
  template class WaveEquationOperationElasticWave<2, 2, float>;
  template class WaveEquationOperationElasticWave<2, 3, float>;
  template class WaveEquationOperationElasticWave<2, 4, float>;
  template class WaveEquationOperationElasticWave<2, 5, float>;
  template class WaveEquationOperationElasticWave<2, 6, float>;
  template class WaveEquationOperationElasticWave<3, 1, float>;
  template class WaveEquationOperationElasticWave<3, 2, float>;
  template class WaveEquationOperationElasticWave<3, 3, float>;
  template class WaveEquationOperationElasticWave<3, 4, float>;
  template class WaveEquationOperationElasticWave<3, 5, float>;
  template class WaveEquationOperationElasticWave<3, 6, float>;
}  // namespace ACOU


#endif  // HAVE_DEAL_II
