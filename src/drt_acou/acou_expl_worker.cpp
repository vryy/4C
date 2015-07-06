/*!----------------------------------------------------------------------
\file acou_expl_worker.cpp
\brief Control routine for acoustic explicit time integration.

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*----------------------------------------------------------------------*/

#include "acou_expl_worker.H"

#ifdef HAVE_DEAL_II

#include <deal.II/lac/parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/base/timer.h>

#include <Epetra_MpiComm.h>

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_mat/acoustic.H"
#include "acou_ele.H"


namespace ACOU
{
  namespace internal
  {
    template <int dim, typename Number>
    MatrixFree<dim,Number>
    create_matrix_free(const DoFHandler<dim> &dof_handler,
                       const unsigned int     fe_degree,
                       const Epetra_Comm     &comm)
    {
      if (fe_degree != dof_handler.get_fe().degree)
        dserror("Internal error in element degree detection");

      QGauss<1> quadrature (fe_degree+1);
      typename MatrixFree<dim,Number>::AdditionalData additional_data;

      const Epetra_MpiComm* mpi_comm = dynamic_cast<const Epetra_MpiComm*>(&comm);
      if (mpi_comm == 0)
        dserror("The Epetra MPI communicator is not derived from Epetra_MpiComm. Fatal error.");

      additional_data.mpi_communicator = mpi_comm->Comm();
      additional_data.tasks_parallel_scheme = MatrixFree<dim,Number>::AdditionalData::partition_partition;
      additional_data.build_face_info = true;
      additional_data.mapping_update_flags = (update_gradients | update_JxW_values |
                                              update_quadrature_points | update_normal_vectors |
                                              update_values);
      ConstraintMatrix dummy;
      dummy.close();
      MatrixFree<dim,Number> data;
      data.reinit (dof_handler, dummy, quadrature, additional_data);

      return data;
    }
  }



// TODO: also need to have a Mapping for representing curved boundaries
template<int dim, int fe_degree>
WaveEquationOperation<dim,fe_degree>::
WaveEquationOperation(const DoFHandler<dim> &dof_handler,
                      Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                      Teuchos::RCP<Function<dim> > boundary_conditions,
                      Teuchos::RCP<Function<dim> > source_term)
  :
  data(internal::create_matrix_free<dim,value_type>(dof_handler, fe_degree,
                                                    discret->Comm())),
  time(0.),
  computing_times(3),
  dirichlet_boundary_conditions(boundary_conditions),
  source_term(source_term),
  mass_matrix_data(data)
{
  densities.resize(data.n_macro_cells()+data.n_macro_ghost_cells());
  speeds.resize(data.n_macro_cells()+data.n_macro_ghost_cells());
  for (unsigned int i=0; i<data.n_macro_cells()+data.n_macro_ghost_cells(); ++i)
  {
    densities[i] = make_vectorized_array<value_type>(1.);
    speeds[i] = make_vectorized_array<value_type>(1.);
    for (unsigned int v=0; v<data.n_components_filled(i); ++v)
      {
        if (data.get_cell_iterator(i,v)->level() != 0)
          dserror("Refined meshes currently not implemented!");

        const int element_index = data.get_cell_iterator(i,v)->index();
        Teuchos::RCP<MAT::Material> mat = discret->lColElement(element_index)->Material();
        MAT::AcousticMat* actmat = static_cast<MAT::AcousticMat*>(mat.get());

        densities[i][v] = actmat->Density();
        speeds[i][v] = actmat->SpeedofSound();
      }
  }

  ConditionalOStream pcout(std::cout, discret->Comm().MyPID() == 0);
  data.print_memory_consumption(pcout);

}



template <int dim, int fe_degree>
WaveEquationOperation<dim,fe_degree>::~WaveEquationOperation()
{
  if (computing_times[2] > 0)
    std::cout << "Computing " << (std::size_t)computing_times[2]
              << " times: evaluate "
              << computing_times[0] << "s, inv mass: " << computing_times[1]
              << "s" << std::endl;
}



template<int dim, int fe_degree>
void
WaveEquationOperation<dim,fe_degree>::
read_initial_conditions(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                        std::vector<parallel::distributed::Vector<value_type> > &dst) const
{
  FEEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> phi(data);
  for (unsigned int j=0; j<phi.dofs_per_cell; ++j)
    phi.submit_dof_value(Tensor<1,dim+1,VectorizedArray<value_type> >(), j);

  unsigned int dofs_per_cell = phi.dofs_per_cell; // i assume, that it is the same for all cells
  unsigned int ndofs1d;
  if(dim==2)
    ndofs1d = std::sqrt(dofs_per_cell);
  else if(dim==3)
    ndofs1d = int(std::pow(dofs_per_cell,1.0/3.0));

  unsigned int nodes_per_cell = GeometryInfo< dim >::vertices_per_cell;
  std::vector<Point<dim> > baci_vals_loc(nodes_per_cell);
  std::vector<Point<dim> > deal_vals_loc(nodes_per_cell);

  for (unsigned int i=0; i<data.n_macro_cells(); ++i)
  {
    phi.reinit(i);
    for (unsigned int v=0; v<data.n_components_filled(i); ++v)
    {

      const int element_index = data.get_cell_iterator(i,v)->index();
      DRT::ELEMENTS::Acou * acouele = dynamic_cast<DRT::ELEMENTS::Acou*>(discret->lColElement(element_index));
      if (acouele == NULL)
        dserror("No acoustic element given!");

      // perform permutation: step 1: get the node coordinates
      for (unsigned int n=0; n<nodes_per_cell; ++n)
      {
        for(int d=0; d<dim; ++d)
        {
          deal_vals_loc[n](d) = data.get_cell_iterator(i,v)->vertex(n)(d);
          baci_vals_loc[n](d) = acouele->Nodes()[n]->X()[d];
        }
      }

      // perform permutation: step 2: swap it
      switch(acouele->Shape())
      {
      case DRT::Element::quad4:
      {
        if(deal_vals_loc[0].distance(baci_vals_loc[0])<1e-10 &&
           deal_vals_loc[1].distance(baci_vals_loc[1])<1e-10 &&
           deal_vals_loc[2].distance(baci_vals_loc[3])<1e-10 &&
           deal_vals_loc[3].distance(baci_vals_loc[2])<1e-10)
        {
          // everything is alright
          for (unsigned j=0; j<dofs_per_cell; ++j)
          {
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[d*dofs_per_cell+j][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+j);
            phi.begin_dof_values()[dim*dofs_per_cell+j][v] = acouele->eleinteriorPressnp_(j);
          }
        }
        else if(deal_vals_loc[0].distance(baci_vals_loc[3])<1e-10 &&
                deal_vals_loc[1].distance(baci_vals_loc[0])<1e-10 &&
                deal_vals_loc[2].distance(baci_vals_loc[2])<1e-10 &&
                deal_vals_loc[3].distance(baci_vals_loc[1])<1e-10)
        {
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const int ax = i%ndofs1d;
            const int ay = i/ndofs1d;
            int permute = (ndofs1d-1-ax)*ndofs1d + ay;
            phi.begin_dof_values()[dim*dofs_per_cell+i][v] = acouele->eleinteriorPressnp_(permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[d*dofs_per_cell+i][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+permute);
          }
        }
        else if(deal_vals_loc[0].distance(baci_vals_loc[2])<1e-10 &&
                deal_vals_loc[1].distance(baci_vals_loc[3])<1e-10 &&
                deal_vals_loc[2].distance(baci_vals_loc[1])<1e-10 &&
                deal_vals_loc[3].distance(baci_vals_loc[0])<1e-10)
        {
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const int ax = i%ndofs1d;
            const int ay = i/ndofs1d;
            int permute = (ndofs1d-1-ax) + (ndofs1d-1-ay) * ndofs1d;
            phi.begin_dof_values()[dim*dofs_per_cell+i][v] = acouele->eleinteriorPressnp_(permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[d*dofs_per_cell+i][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+permute);
          }
        }
        else if(deal_vals_loc[0].distance(baci_vals_loc[1])<1e-10 &&
                deal_vals_loc[1].distance(baci_vals_loc[2])<1e-10 &&
                deal_vals_loc[2].distance(baci_vals_loc[0])<1e-10 &&
                deal_vals_loc[3].distance(baci_vals_loc[3])<1e-10)
        {
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const int ax = i%ndofs1d;
            const int ay = i/ndofs1d;
            int permute = (ax) * ndofs1d + (ndofs1d-1-ay);
            phi.begin_dof_values()[dim*dofs_per_cell+i][v] = acouele->eleinteriorPressnp_(permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[d*dofs_per_cell+i][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+permute);
          }
        }
        else
        {
          std::cout<<"d "<<deal_vals_loc[0](0)<<" "<<deal_vals_loc[0](1)<<" b "<<baci_vals_loc[0](0)<<" "<<baci_vals_loc[0](1)<<" "<<std::endl;
          std::cout<<"d "<<deal_vals_loc[1](0)<<" "<<deal_vals_loc[1](1)<<" b "<<baci_vals_loc[1](0)<<" "<<baci_vals_loc[1](1)<<" "<<std::endl;
          std::cout<<"d "<<deal_vals_loc[2](0)<<" "<<deal_vals_loc[2](1)<<" b "<<baci_vals_loc[2](0)<<" "<<baci_vals_loc[2](1)<<" "<<std::endl;
          std::cout<<"d "<<deal_vals_loc[3](0)<<" "<<deal_vals_loc[3](1)<<" b "<<baci_vals_loc[3](0)<<" "<<baci_vals_loc[3](1)<<" "<<std::endl;
          dserror("unkown permutation");
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

  for (unsigned int i=0; i<dim+1; ++i)
  {
    dst[i].update_ghost_values();
  }
}

template<int dim, int fe_degree>
void
WaveEquationOperation<dim,fe_degree>::write_deal_cell_values(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
    const std::vector<parallel::distributed::Vector<value_type> >   &src) const
{
  const unsigned dofs_per_cell = data.get_dof_handler().get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> indices, local_dof_indices (dofs_per_cell);
  for (typename DoFHandler<dim>::active_cell_iterator cell = data.get_dof_handler().begin_active();
       cell != data.get_dof_handler().end(); ++cell)
  {
    cell->get_dof_indices(local_dof_indices);
    indices.insert(indices.end(), local_dof_indices.begin(), local_dof_indices.end());
  }
  IndexSet relevant_dofs(src[0].size());
  relevant_dofs.add_indices(indices.begin(), indices.end());
  relevant_dofs.compress();
  std::vector<parallel::distributed::Vector<value_type> > ghosted_vector(src.size());
  for (unsigned int i=0; i<src.size(); ++i)
  {
    ghosted_vector[i].reinit(data.get_dof_handler().locally_owned_dofs(),
                             relevant_dofs, src[0].get_mpi_communicator());
    ghosted_vector[i] = src[i];
    ghosted_vector[i].update_ghost_values();
  }

  unsigned int ndofs1d;
  if(dim==2)
    ndofs1d = std::sqrt(dofs_per_cell);
  else if(dim==3)
    ndofs1d = int(std::pow(dofs_per_cell,1.0/3.0));

  unsigned int nodes_per_cell = GeometryInfo< dim >::vertices_per_cell;
  std::vector<Point<dim> > baci_vals_loc(nodes_per_cell);
  std::vector<Point<dim> > deal_vals_loc(nodes_per_cell);

  Vector<value_type> local_values(dofs_per_cell);
  for (int i=0; i<discret->NumMyColElements(); ++i)
  {
    typename DoFHandler<dim>::active_cell_iterator cell(&data.get_dof_handler().get_tria(), 0, i, &data.get_dof_handler());
    DRT::ELEMENTS::Acou * acouele = dynamic_cast<DRT::ELEMENTS::Acou*>(discret->lColElement(i));

    for (unsigned int n=0; n<nodes_per_cell; ++n)
    {
      for(int d=0; d<dim; ++d)
      {
        deal_vals_loc[n](d) = cell->vertex(n)(d);
        baci_vals_loc[n](d) = acouele->Nodes()[n]->X()[d];
      }
    }

    // perform permutation: step 2: swap it
    switch(acouele->Shape())
    {
    case DRT::Element::quad4:
    {
      if((deal_vals_loc[0].distance(baci_vals_loc[0])<1e-10 &&
         deal_vals_loc[1].distance(baci_vals_loc[1])<1e-10 &&
         deal_vals_loc[2].distance(baci_vals_loc[3])<1e-10 &&
         deal_vals_loc[3].distance(baci_vals_loc[2])<1e-10))
      {
        // everything is alright
        for (unsigned int d=0; d<dim; ++d)
        {
          cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            acouele->eleinteriorVelnp_(d*dofs_per_cell+j) = local_values[j];
        }
        cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          acouele->eleinteriorPressnp_(j) = local_values[j];
      }
      else if(deal_vals_loc[0].distance(baci_vals_loc[3])<1e-10 &&
              deal_vals_loc[1].distance(baci_vals_loc[0])<1e-10 &&
              deal_vals_loc[2].distance(baci_vals_loc[2])<1e-10 &&
              deal_vals_loc[3].distance(baci_vals_loc[1])<1e-10)
      {
        for (unsigned int d=0; d<dim; ++d)
        {
          cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = j/ndofs1d;
            int permute = (ndofs1d-1-ax)*ndofs1d + ay;
            acouele->eleinteriorVelnp_(d*dofs_per_cell+permute) = local_values[j];
          }
        }
        cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
        for (unsigned int j=0; j<dofs_per_cell; ++j)
        {
          const int ax = j%ndofs1d;
          const int ay = j/ndofs1d;
          int permute = (ndofs1d-1-ax)*ndofs1d + ay;
          acouele->eleinteriorPressnp_(permute) = local_values[j];
        }
      }
      else if(deal_vals_loc[0].distance(baci_vals_loc[2])<1e-10 &&
              deal_vals_loc[1].distance(baci_vals_loc[3])<1e-10 &&
              deal_vals_loc[2].distance(baci_vals_loc[1])<1e-10 &&
              deal_vals_loc[3].distance(baci_vals_loc[0])<1e-10)
      {
        for (unsigned int d=0; d<dim; ++d)
        {
          cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = j/ndofs1d;
            int permute = (ndofs1d-1-ax) + (ndofs1d-1-ay) * ndofs1d;
            acouele->eleinteriorVelnp_(d*dofs_per_cell+permute) = local_values[j];
          }
        }
        cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
        for (unsigned int j=0; j<dofs_per_cell; ++j)
        {
          const int ax = j%ndofs1d;
          const int ay = j/ndofs1d;
          int permute = (ndofs1d-1-ax) + (ndofs1d-1-ay) * ndofs1d;
          acouele->eleinteriorPressnp_(permute) = local_values[j];
        }
      }
      else if(deal_vals_loc[0].distance(baci_vals_loc[1])<1e-10 &&
              deal_vals_loc[1].distance(baci_vals_loc[2])<1e-10 &&
              deal_vals_loc[2].distance(baci_vals_loc[0])<1e-10 &&
              deal_vals_loc[3].distance(baci_vals_loc[3])<1e-10)
      {
        for (unsigned int d=0; d<dim; ++d)
        {
          cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = j/ndofs1d;
            int permute = (ax) * ndofs1d + (ndofs1d-1-ay);
            acouele->eleinteriorVelnp_(d*dofs_per_cell+permute) = local_values[j];
          }
        }
        cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
        for (unsigned int j=0; j<dofs_per_cell; ++j)
        {
          const int ax = j%ndofs1d;
          const int ay = j/ndofs1d;
          int permute = (ax) * ndofs1d + (ndofs1d-1-ay);
          acouele->eleinteriorPressnp_(permute) = local_values[j];
        }
      }
      else
        dserror("unkown permutation");
      break;
    }
    default:
      dserror("other distypes not yet implemented!");
      break;
    }
  }

  return;
}

template<int dim, int fe_degree>
void WaveEquationOperation<dim, fe_degree>::
local_apply_domain(const MatrixFree<dim,value_type>                               &data,
                   std::vector<parallel::distributed::Vector<value_type> >         &dst,
                   const std::vector<parallel::distributed::Vector<value_type> >   &src,
                   const std::pair<unsigned int,unsigned int>                 &cell_range) const
{
  FEEvaluation<dim,fe_degree,fe_degree+1,dim,value_type> velocity(data);
  FEEvaluation<dim,fe_degree,fe_degree+1,1,value_type> pressure(data);

  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
  {
    // It is faster to evaluate values of the vector-valued velocity and
    // gradients of the scalar pressure than divergence of velocity and
    // values of pressure
    velocity.reinit(cell);
    velocity.read_dof_values(src, 0);
    velocity.evaluate (true, false, false);

    pressure.reinit(cell);
    pressure.read_dof_values(src, dim);
    pressure.evaluate(false, true, false);

    const VectorizedArray<value_type> rho = densities[cell];
    const VectorizedArray<value_type> rho_inv = 1./densities[cell];
    const VectorizedArray<value_type> c_sq = speeds[cell]*speeds[cell];

    for (unsigned int q=0; q<velocity.n_q_points; ++q)
    {
      const Tensor<1,dim,VectorizedArray<value_type> >
      pressure_gradient = pressure.get_gradient(q);
      const Tensor<1,dim,VectorizedArray<value_type> >
      velocity_value = velocity.get_value(q);

      Point<dim,VectorizedArray<value_type> > q_points = velocity.quadrature_point(q);
      VectorizedArray<value_type> rhs;
      for (unsigned int n=0; n<rhs.n_array_elements; ++n)
      {
        Point<dim> q_point;
        for (unsigned int d=0; d<dim; ++d)
          q_point[d] = q_points[d][n];
        rhs[n] = source_term->value(q_point);
      }
      velocity.submit_value(-rho_inv*pressure_gradient,q);

      pressure.submit_value(c_sq*rhs,q);
      pressure.submit_gradient(rho*c_sq*velocity_value,q);
    }

    velocity.integrate (true, false);
    velocity.distribute_local_to_global (dst, 0);

    pressure.integrate(true, true);
    pressure.distribute_local_to_global (dst,dim);
  }
}



template <int dim, int fe_degree>
void
WaveEquationOperation<dim,fe_degree>::
local_apply_face (const MatrixFree<dim,value_type> &,
  std::vector<parallel::distributed::Vector<value_type> >        &dst,
  const std::vector<parallel::distributed::Vector<value_type> >  &src,
  const std::pair<unsigned int,unsigned int>                &face_range) const
{
  // There is some overhead in the methods in FEEvaluation, so it is faster
  // to combine pressure and velocity in the same object and just combine
  // them at the level of quadrature points
  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> phi(this->data, true, 0, 0, true);
  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> phi_neighbor(this->data, false, 0, 0, true);

  for (unsigned int face=face_range.first; face<face_range.second; face++)
  {
    phi.reinit(face);
    phi.read_dof_values(src, 0);
    phi.evaluate(true,false);
    const VectorizedArray<value_type> rho_plus = phi.read_cell_data(densities);
    const VectorizedArray<value_type> rho_inv_plus = 1./rho_plus;
    const VectorizedArray<value_type> c_plus = phi.read_cell_data(speeds);
    const VectorizedArray<value_type> c_sq_plus = c_plus * c_plus;
    const VectorizedArray<value_type> tau_plus = 1./c_plus;

    phi_neighbor.reinit(face);
    phi_neighbor.read_dof_values(src, 0);
    phi_neighbor.evaluate(true,false);
    const VectorizedArray<value_type> rho_minus = phi_neighbor.read_cell_data(densities);
    const VectorizedArray<value_type> rho_inv_minus = 1./rho_minus;
    const VectorizedArray<value_type> c_minus = phi_neighbor.read_cell_data(speeds);
    const VectorizedArray<value_type> c_sq_minus = c_minus * c_minus;
    const VectorizedArray<value_type> tau_minus = 1./c_minus;

    const VectorizedArray<value_type> tau_inv = 1./(tau_plus + tau_minus);

    AssertDimension(phi.n_q_points, data.get_n_q_points_face(0));

    for (unsigned int q=0; q<phi.n_q_points; ++q)
    {
      Tensor<1,dim+1,VectorizedArray<value_type> > val_plus = phi.get_value(q);
      Tensor<1,dim+1,VectorizedArray<value_type> > val_minus = phi_neighbor.get_value(q);
      Tensor<1,dim,VectorizedArray<value_type> > normal = phi.get_normal_vector(q);
      VectorizedArray<value_type> normal_v_plus = val_plus[0]*normal[0];
      VectorizedArray<value_type> normal_v_minus = -val_minus[0]*normal[0];
      for (unsigned int d=1; d<dim; ++d)
      {
        normal_v_plus += val_plus[d] * normal[d];
        normal_v_minus -= val_minus[d] * normal[d];
      }

      VectorizedArray<value_type> lambda = tau_inv*(rho_plus*normal_v_plus+rho_minus*normal_v_minus + tau_plus*val_plus[dim] + tau_minus*val_minus[dim]);
      VectorizedArray<value_type> pres_diff_plus = (val_plus[dim]-lambda)*rho_inv_plus;
      VectorizedArray<value_type> pres_diff_minus = (val_minus[dim]-lambda)*rho_inv_minus;
      for (unsigned int d=0; d<dim; ++d)
      {
        val_plus[d] = pres_diff_plus*normal[d];
        val_minus[d] = -pres_diff_minus*normal[d];
      }
      val_plus[dim] = c_sq_plus * (-rho_plus*normal_v_plus + tau_plus * (lambda - val_plus[dim]));
      val_minus[dim] = c_sq_minus * (-rho_minus*normal_v_minus + tau_minus * (lambda - val_minus[dim]));

      phi.submit_value(val_plus, q);
      phi_neighbor.submit_value(val_minus, q);
    }

    phi.integrate(true,false);
    phi.distribute_local_to_global(dst, 0);

    phi_neighbor.integrate(true,false);
    phi_neighbor.distribute_local_to_global(dst, 0);
  }
}



template <int dim, int fe_degree>
void WaveEquationOperation<dim,fe_degree>::
local_apply_boundary_face (const MatrixFree<dim,value_type> &,
                           std::vector<parallel::distributed::Vector<value_type> >       &dst,
                           const std::vector<parallel::distributed::Vector<value_type> > &src,
                           const std::pair<unsigned int,unsigned int>               &face_range) const
{
  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> phi(this->data, true, 0, 0, true);

  for (unsigned int face=face_range.first; face<face_range.second; face++)
  {
    phi.reinit(face);
    phi.read_dof_values(src, 0);
    phi.evaluate(true,false);
    const VectorizedArray<value_type> rho = phi.read_cell_data(densities);
    const VectorizedArray<value_type> rho_inv = 1./rho;
    const VectorizedArray<value_type> c_sq = phi.read_cell_data(speeds)*phi.read_cell_data(speeds);
    const VectorizedArray<value_type> c = phi.read_cell_data(speeds);
    const VectorizedArray<value_type> tau = 1./phi.read_cell_data(speeds);

    const types::boundary_id boundary_index = this->data.get_boundary_indicator(face);
    const int int_boundary_id = int(boundary_index);

    for (unsigned int q=0; q<phi.n_q_points; ++q)
    {
      Tensor<1,dim,VectorizedArray<value_type> > normal = phi.get_normal_vector(q);
      Tensor<1,dim+1,VectorizedArray<value_type> > val_plus = phi.get_value(q);
      VectorizedArray<value_type> p_plus = val_plus[dim];
      VectorizedArray<value_type> normal_v_plus = val_plus[0] * normal[0];
      for (unsigned int d=1; d<dim; ++d)
        normal_v_plus += val_plus[d] * normal[d];
      Point<dim,VectorizedArray<value_type> > q_point = phi.quadrature_point(q);
      VectorizedArray<value_type> lambda;

      if(int_boundary_id==0) // absorbing boundary
        lambda = tau/(tau+1./c)*p_plus + rho/(tau+1./c)*normal_v_plus;
      else if(int_boundary_id==1) // free boundary
        lambda = VectorizedArray<value_type>();
      else if(int_boundary_id>=2) // dbcs
      {
        for (unsigned int v=0; v<VectorizedArray<value_type>::n_array_elements; ++v)
        {
          Point<dim> point;
          for (unsigned int d=0; d<dim; ++d)
            point[d] = q_point[d][v];
          lambda[v] = dirichlet_boundary_conditions->value(point,int_boundary_id-2);
        }
      }

      for (unsigned int d=0; d<dim; ++d)
        val_plus[d] = (p_plus-lambda)*rho_inv*normal[d];
      val_plus[dim] = c_sq*(-rho*normal_v_plus+tau*(lambda - p_plus));
      phi.submit_value(val_plus,q);
    }
    phi.integrate(true,false);
    phi.distribute_local_to_global(dst, 0);
  }
}



template<int dim, int fe_degree>
void WaveEquationOperation<dim, fe_degree>::
local_apply_mass_matrix(const MatrixFree<dim,value_type>                  &data,
                        std::vector<parallel::distributed::Vector<value_type> >        &dst,
                        const std::vector<parallel::distributed::Vector<value_type> >  &src,
                        const std::pair<unsigned int,unsigned int>    &cell_range) const
{
  internal::InverseMassMatrixData<dim,fe_degree,value_type>& mass_data = mass_matrix_data.get();
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
    {
      mass_data.phi[0].reinit(cell);
      mass_data.phi[0].read_dof_values(src, 0);

      mass_data.inverse.fill_inverse_JxW_values(mass_data.coefficients);
      mass_data.inverse.apply(mass_data.coefficients, dim+1,
                              mass_data.phi[0].begin_dof_values(),
                              mass_data.phi[0].begin_dof_values());

      mass_data.phi[0].set_dof_values(dst,0);
    }
}



template<int dim, int fe_degree>
void WaveEquationOperation<dim, fe_degree>::
apply(const std::vector<parallel::distributed::Vector<value_type> >  &src,
      std::vector<parallel::distributed::Vector<value_type> >        &dst,
      const double                                              &cur_time) const
{
  Timer timer;
  time = cur_time;
  dirichlet_boundary_conditions->set_time(time);
  source_term->set_time(time);

  data.loop (&WaveEquationOperation<dim, fe_degree>::local_apply_domain,
             &WaveEquationOperation<dim, fe_degree>::local_apply_face,
             &WaveEquationOperation<dim, fe_degree>::local_apply_boundary_face,
             this, dst, src);

  computing_times[0] += timer.wall_time();
  timer.restart();

  data.cell_loop(&WaveEquationOperation<dim, fe_degree>::local_apply_mass_matrix,
                 this, dst, dst);

  computing_times[1] += timer.wall_time();
  computing_times[2] += 1.;
}



// explicit instantiations
template class WaveEquationOperation<2,1>;
template class WaveEquationOperation<2,2>;
template class WaveEquationOperation<2,3>;
template class WaveEquationOperation<2,4>;
template class WaveEquationOperation<2,5>;
template class WaveEquationOperation<2,6>;
template class WaveEquationOperation<3,1>;
template class WaveEquationOperation<3,2>;
template class WaveEquationOperation<3,3>;
template class WaveEquationOperation<3,4>;
template class WaveEquationOperation<3,5>;
template class WaveEquationOperation<3,6>;
}


#endif // HAVE_DEAL_II
