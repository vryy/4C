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
#include "../drt_mat/acoustic_sol.H"
#include "acou_ele.H"
#include "acou_sol_ele.H"


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
                      Teuchos::RCP<Function<dim> > source_term,
                      bool sol_in,
                      Teuchos::RCP<Epetra_MultiVector> source_adjoint)
  :
  data(internal::create_matrix_free<dim,value_type>(dof_handler, fe_degree,
                                                    discret->Comm())),
  time(0.),
  computing_times(3),
  source_adjoint_meas(source_adjoint),
  dirichlet_boundary_conditions(boundary_conditions),
  source_term(source_term),
  mass_matrix_data(data),
  mass_matrix_data_solid(data)
{
  this->solid = sol_in;

  densities.resize(data.n_macro_cells()+data.n_macro_ghost_cells());
  speeds.resize(data.n_macro_cells()+data.n_macro_ghost_cells());
  if(this->solid)
    viscs.resize(data.n_macro_cells()+data.n_macro_ghost_cells());

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

      if(this->solid)
      {
        MAT::AcousticSolMat* actmat = static_cast<MAT::AcousticSolMat*>(mat.get());
        densities[i][v] = actmat->Density();
        speeds[i][v] = actmat->SpeedofSound();
        viscs[i][v] = actmat->Viscosity();
      }
      else
      {
        MAT::AcousticMat* actmat = static_cast<MAT::AcousticMat*>(mat.get());
        densities[i][v] = actmat->Density(element_index);
        speeds[i][v] = actmat->SpeedofSound(element_index);
      }
    }
  }

  // only in case of the adjoint run in inverse analysis
  if(source_adjoint_meas!= Teuchos::null)
  {
    // create TableIndices with lengths for table rows and columns
    TableIndices<3> table_indices_ids(data.n_macro_boundary_faces(),VectorizedArray<value_type>::n_array_elements,GeometryInfo<dim>::vertices_per_face);
    TableIndices<4> table_indices_coords(data.n_macro_boundary_faces(),VectorizedArray<value_type>::n_array_elements,GeometryInfo<dim>::vertices_per_face,dim);

    // resize (or init) the tables
    table_node_ids.reinit(table_indices_ids);
    table_node_coords.reinit(table_indices_coords);

    for (unsigned int f=0; f<data.n_macro_boundary_faces(); ++f)
      for (unsigned int v=0; v<VectorizedArray<value_type>::n_array_elements &&
           data.faces[data.n_macro_inner_faces()+f].left_cell[v] !=
               numbers::invalid_unsigned_int; ++v)
      {
        const unsigned int cell_index_non_vectorized = data.faces[data.n_macro_inner_faces()+f].left_cell[v];
        const int element_index = data.get_cell_iterator(cell_index_non_vectorized / VectorizedArray<value_type>::n_array_elements,
                                                         cell_index_non_vectorized % VectorizedArray<value_type>::n_array_elements)->index();
        unsigned int count = 0;
        for(int n=0; n<discret->lColElement(element_index)->NumNode(); ++n)
        {
          unsigned int global_node_id = discret->lColElement(element_index)->NodeIds()[n];
          if(source_adjoint_meas->Map().LID(int(global_node_id))>=0)
          {
            table_node_ids(f,v,count) = global_node_id;
            for(unsigned int d=0; d<dim; ++d)
              table_node_coords(f,v,count,d) = discret->lColElement(element_index)->Nodes()[n]->X()[d];
            count++;
          }
        }
      }

    // in case of acouopt, we need the following quantities:
    densities_grad.resize(data.n_macro_cells()+data.n_macro_ghost_cells());
    speeds_grad.resize(data.n_macro_cells()+data.n_macro_ghost_cells());
    for (unsigned int i=0; i<data.n_macro_cells()+data.n_macro_ghost_cells(); ++i)
    {
      densities_grad[i] = make_vectorized_array<value_type>(0.);
      speeds_grad[i] = make_vectorized_array<value_type>(0.);
    }
  }
  //ConditionalOStream pcout(std::cout, discret->Comm().MyPID() == 0);
  //data.print_memory_consumption(pcout);
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
  if(this->solid==false)
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
    unsigned int ndofs2d = ndofs1d * ndofs1d;

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
            dserror("unknown permutation");
          }
          break;
        }
        case DRT::Element::hex8:
        {
          if(deal_vals_loc[0].distance(baci_vals_loc[0])<1e-10 &&
             deal_vals_loc[1].distance(baci_vals_loc[1])<1e-10 &&
             deal_vals_loc[2].distance(baci_vals_loc[3])<1e-10 &&
             deal_vals_loc[3].distance(baci_vals_loc[2])<1e-10 &&
             deal_vals_loc[4].distance(baci_vals_loc[4])<1e-10 &&
             deal_vals_loc[5].distance(baci_vals_loc[5])<1e-10 &&
             deal_vals_loc[6].distance(baci_vals_loc[7])<1e-10 &&
             deal_vals_loc[7].distance(baci_vals_loc[6])<1e-10)
          {
            // everything is alright
            for (unsigned j=0; j<dofs_per_cell; ++j)
            {
              for (unsigned int d=0; d<dim; ++d)
                phi.begin_dof_values()[d*dofs_per_cell+j][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+j);
              phi.begin_dof_values()[dim*dofs_per_cell+j][v] = acouele->eleinteriorPressnp_(j);
            }
          }
          else if(deal_vals_loc[0].distance(baci_vals_loc[4])<1e-10 &&
                  deal_vals_loc[1].distance(baci_vals_loc[5])<1e-10 &&
                  deal_vals_loc[2].distance(baci_vals_loc[0])<1e-10 &&
                  deal_vals_loc[3].distance(baci_vals_loc[1])<1e-10 &&
                  deal_vals_loc[4].distance(baci_vals_loc[7])<1e-10 &&
                  deal_vals_loc[5].distance(baci_vals_loc[6])<1e-10 &&
                  deal_vals_loc[6].distance(baci_vals_loc[3])<1e-10 &&
                  deal_vals_loc[7].distance(baci_vals_loc[2])<1e-10)
          {
            // negative rotation around x
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              const int ax = i%ndofs1d;
              const int ay = int(i/ndofs1d)%ndofs1d;
              const int az = i/ndofs2d;
              int permute = ax + (ndofs1d-1-ay)*ndofs2d + az*ndofs1d;
              phi.begin_dof_values()[dim*dofs_per_cell+i][v] = acouele->eleinteriorPressnp_(permute);
                for (unsigned int d=0; d<dim; ++d)
                  phi.begin_dof_values()[d*dofs_per_cell+i][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+permute);
            }
          }
          else if(deal_vals_loc[0].distance(baci_vals_loc[5])<1e-10 &&
                  deal_vals_loc[1].distance(baci_vals_loc[6])<1e-10 &&
                  deal_vals_loc[2].distance(baci_vals_loc[1])<1e-10 &&
                  deal_vals_loc[3].distance(baci_vals_loc[2])<1e-10 &&
                  deal_vals_loc[4].distance(baci_vals_loc[4])<1e-10 &&
                  deal_vals_loc[5].distance(baci_vals_loc[7])<1e-10 &&
                  deal_vals_loc[6].distance(baci_vals_loc[0])<1e-10 &&
                  deal_vals_loc[7].distance(baci_vals_loc[3])<1e-10)
          {
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              const int ax = i%ndofs1d;
              const int ay = int(i/ndofs1d)%ndofs1d;
              const int az = i/ndofs2d;
              int permute = ax*ndofs1d + (ndofs1d-1-ay)*ndofs2d + (ndofs1d-1-az);
              phi.begin_dof_values()[dim*dofs_per_cell+i][v] = acouele->eleinteriorPressnp_(permute);
              for (unsigned int d=0; d<dim; ++d)
                phi.begin_dof_values()[d*dofs_per_cell+i][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+permute);
            }
          }
          else if(deal_vals_loc[0].distance(baci_vals_loc[7])<1e-10 &&
                  deal_vals_loc[1].distance(baci_vals_loc[4])<1e-10 &&
                  deal_vals_loc[2].distance(baci_vals_loc[3])<1e-10 &&
                  deal_vals_loc[3].distance(baci_vals_loc[0])<1e-10 &&
                  deal_vals_loc[4].distance(baci_vals_loc[6])<1e-10 &&
                  deal_vals_loc[5].distance(baci_vals_loc[5])<1e-10 &&
                  deal_vals_loc[6].distance(baci_vals_loc[2])<1e-10 &&
                  deal_vals_loc[7].distance(baci_vals_loc[1])<1e-10)
          {
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              const int ax = i%ndofs1d;
              const int ay = int(i/ndofs1d)%ndofs1d;
              const int az = i/ndofs2d;
              int permute = (ndofs1d-1-ax)*ndofs1d + (ndofs1d-1-ay)*ndofs2d + az;
              phi.begin_dof_values()[dim*dofs_per_cell+i][v] = acouele->eleinteriorPressnp_(permute);
              for (unsigned int d=0; d<dim; ++d)
                phi.begin_dof_values()[d*dofs_per_cell+i][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+permute);
            }
          }
          else if(deal_vals_loc[0].distance(baci_vals_loc[6])<1e-10 &&
                  deal_vals_loc[1].distance(baci_vals_loc[7])<1e-10 &&
                  deal_vals_loc[2].distance(baci_vals_loc[2])<1e-10 &&
                  deal_vals_loc[3].distance(baci_vals_loc[3])<1e-10 &&
                  deal_vals_loc[4].distance(baci_vals_loc[5])<1e-10 &&
                  deal_vals_loc[5].distance(baci_vals_loc[4])<1e-10 &&
                  deal_vals_loc[6].distance(baci_vals_loc[1])<1e-10 &&
                  deal_vals_loc[7].distance(baci_vals_loc[0])<1e-10)
          {
            for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              const int ax = i%ndofs1d;
              const int ay = int(i/ndofs1d)%ndofs1d;
              const int az = i/ndofs2d;
              int permute = (ndofs1d-1-ax) + (ndofs1d-1-ay)*ndofs2d + (ndofs1d-1-az)*ndofs1d;
              phi.begin_dof_values()[dim*dofs_per_cell+i][v] = acouele->eleinteriorPressnp_(permute);
              for (unsigned int d=0; d<dim; ++d)
                phi.begin_dof_values()[d*dofs_per_cell+i][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+permute);
            }
          }
          else
          {
            std::cout<<"d "<<deal_vals_loc[0](0)<<" "<<deal_vals_loc[0](1)<<" "<<deal_vals_loc[0](2)<<" b "<<baci_vals_loc[0](0)<<" "<<baci_vals_loc[0](1)<<" "<<baci_vals_loc[0](2)<<std::endl;
            std::cout<<"d "<<deal_vals_loc[1](0)<<" "<<deal_vals_loc[1](1)<<" "<<deal_vals_loc[1](2)<<" b "<<baci_vals_loc[1](0)<<" "<<baci_vals_loc[1](1)<<" "<<baci_vals_loc[1](2)<<std::endl;
            std::cout<<"d "<<deal_vals_loc[2](0)<<" "<<deal_vals_loc[2](1)<<" "<<deal_vals_loc[2](2)<<" b "<<baci_vals_loc[2](0)<<" "<<baci_vals_loc[2](1)<<" "<<baci_vals_loc[2](2)<<std::endl;
            std::cout<<"d "<<deal_vals_loc[3](0)<<" "<<deal_vals_loc[3](1)<<" "<<deal_vals_loc[3](2)<<" b "<<baci_vals_loc[3](0)<<" "<<baci_vals_loc[3](1)<<" "<<baci_vals_loc[3](2)<<std::endl;
            std::cout<<"d "<<deal_vals_loc[4](0)<<" "<<deal_vals_loc[4](1)<<" "<<deal_vals_loc[4](2)<<" b "<<baci_vals_loc[4](0)<<" "<<baci_vals_loc[4](1)<<" "<<baci_vals_loc[4](2)<<std::endl;
            std::cout<<"d "<<deal_vals_loc[5](0)<<" "<<deal_vals_loc[5](1)<<" "<<deal_vals_loc[5](2)<<" b "<<baci_vals_loc[5](0)<<" "<<baci_vals_loc[5](1)<<" "<<baci_vals_loc[5](2)<<std::endl;
            std::cout<<"d "<<deal_vals_loc[6](0)<<" "<<deal_vals_loc[6](1)<<" "<<deal_vals_loc[6](2)<<" b "<<baci_vals_loc[6](0)<<" "<<baci_vals_loc[6](1)<<" "<<baci_vals_loc[6](2)<<std::endl;
            std::cout<<"d "<<deal_vals_loc[7](0)<<" "<<deal_vals_loc[7](1)<<" "<<deal_vals_loc[7](2)<<" b "<<baci_vals_loc[7](0)<<" "<<baci_vals_loc[7](1)<<" "<<baci_vals_loc[7](2)<<std::endl;
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

    for (unsigned int i=0; i<dim+1; ++i)
    {
      dst[i].update_ghost_values();
    }
  }
  else
  {
    FEEvaluation<dim,fe_degree,fe_degree+1,dim+1+dim*dim,value_type> phi(data);
    for (unsigned int j=0; j<phi.dofs_per_cell; ++j)
      phi.submit_dof_value(Tensor<1,dim+1+dim*dim,VectorizedArray<value_type> >(), j);

    unsigned int dofs_per_cell = phi.dofs_per_cell; // i assume, that it is the same for all cells
    unsigned int ndofs1d;
    if(dim==2)
      ndofs1d = std::sqrt(dofs_per_cell);
    else if(dim==3)
      ndofs1d = int(std::pow(dofs_per_cell,1.0/3.0));
    //unsigned int ndofs2d = ndofs1d * ndofs1d;

    unsigned int nodes_per_cell = GeometryInfo< dim >::vertices_per_cell;
    std::vector<Point<dim> > baci_vals_loc(nodes_per_cell);
    std::vector<Point<dim> > deal_vals_loc(nodes_per_cell);

    for (unsigned int i=0; i<data.n_macro_cells(); ++i)
    {
      phi.reinit(i);
      for (unsigned int v=0; v<data.n_components_filled(i); ++v)
      {

        const int element_index = data.get_cell_iterator(i,v)->index();
        DRT::ELEMENTS::AcouSol * acouele = dynamic_cast<DRT::ELEMENTS::AcouSol*>(discret->lColElement(element_index));
        if (acouele == NULL)
          dserror("No acoustic solid element given!");

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
              for(unsigned int d=0; d<dim; ++d)
                for(unsigned int e=0; e<dim; ++e)
                  phi.begin_dof_values()[(d*dim+e+dim+1)*dofs_per_cell+j][v] = acouele->eleinteriorGradVelnp_((d*dim+e)*dofs_per_cell+j);
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
              for(unsigned int d=0; d<dim; ++d)
                for(unsigned int e=0; e<dim; ++e)
                  phi.begin_dof_values()[(d*dim+e+dim+1)*dofs_per_cell+i][v] = acouele->eleinteriorGradVelnp_((d*dim+e)*dofs_per_cell+permute);
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
              for(unsigned int d=0; d<dim; ++d)
                for(unsigned int e=0; e<dim; ++e)
                  phi.begin_dof_values()[(d*dim+e+dim+1)*dofs_per_cell+i][v] = acouele->eleinteriorGradVelnp_((d*dim+e)*dofs_per_cell+permute);

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
              for(unsigned int d=0; d<dim; ++d)
                for(unsigned int e=0; e<dim; ++e)
                  phi.begin_dof_values()[(d*dim+e+dim+1)*dofs_per_cell+i][v] = acouele->eleinteriorGradVelnp_((d*dim+e)*dofs_per_cell+permute);

            }
          }
          else
          {
            std::cout<<"d "<<deal_vals_loc[0](0)<<" "<<deal_vals_loc[0](1)<<" b "<<baci_vals_loc[0](0)<<" "<<baci_vals_loc[0](1)<<" "<<std::endl;
            std::cout<<"d "<<deal_vals_loc[1](0)<<" "<<deal_vals_loc[1](1)<<" b "<<baci_vals_loc[1](0)<<" "<<baci_vals_loc[1](1)<<" "<<std::endl;
            std::cout<<"d "<<deal_vals_loc[2](0)<<" "<<deal_vals_loc[2](1)<<" b "<<baci_vals_loc[2](0)<<" "<<baci_vals_loc[2](1)<<" "<<std::endl;
            std::cout<<"d "<<deal_vals_loc[3](0)<<" "<<deal_vals_loc[3](1)<<" b "<<baci_vals_loc[3](0)<<" "<<baci_vals_loc[3](1)<<" "<<std::endl;
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

    for (unsigned int i=0; i<dim+1+dim*dim; ++i)
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
    ghosted_vector[i].reinit(data.get_dof_handler().locally_owned_dofs(),relevant_dofs, src[0].get_mpi_communicator());
    ghosted_vector[i] = src[i];
    ghosted_vector[i].update_ghost_values();
  }

  unsigned int ndofs1d;
  if(dim==2)
    ndofs1d = std::sqrt(dofs_per_cell);
  else if(dim==3)
    ndofs1d = int(std::pow(dofs_per_cell,1.0/3.0));
  unsigned int ndofs2d = ndofs1d * ndofs1d;

  unsigned int nodes_per_cell = GeometryInfo< dim >::vertices_per_cell;
  std::vector<Point<dim> > baci_vals_loc(nodes_per_cell);
  std::vector<Point<dim> > deal_vals_loc(nodes_per_cell);

  Vector<value_type> local_values(dofs_per_cell);
  if(this->solid)
  for (int i=0; i<discret->NumMyColElements(); ++i)
  {
    typename DoFHandler<dim>::active_cell_iterator cell(&data.get_dof_handler().get_tria(), 0, i, &data.get_dof_handler());
    DRT::ELEMENTS::AcouSol * acouele = dynamic_cast<DRT::ELEMENTS::AcouSol*>(discret->lColElement(i));

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
        // stresses
        for(unsigned d=0; d<dim;  ++d)
          for(unsigned e=0; e<dim; ++e)
          {
            cell->get_interpolated_dof_values(ghosted_vector[dim+1+d*dim+e], local_values);
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              acouele->eleinteriorGradVelnp_((d*dim+e)*dofs_per_cell+j) = local_values[j];
          }
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

        for (unsigned int d=0; d<dim; ++d)
          for(unsigned e=0; e<dim; ++e)
        {
          cell->get_interpolated_dof_values(ghosted_vector[dim+1+d*dim+e], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = j/ndofs1d;
            int permute = (ndofs1d-1-ax)*ndofs1d + ay;
            acouele->eleinteriorGradVelnp_((d*dim+e)*dofs_per_cell+permute) = local_values[j];
          }
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

        for (unsigned int d=0; d<dim; ++d)
          for(unsigned e=0; e<dim; ++e)
        {
          cell->get_interpolated_dof_values(ghosted_vector[dim+1+d*dim+e], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = j/ndofs1d;
            int permute = (ndofs1d-1-ax) + (ndofs1d-1-ay) * ndofs1d;
            acouele->eleinteriorGradVelnp_((d*dim+e)*dofs_per_cell+permute) = local_values[j];
          }
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

        for (unsigned int d=0; d<dim; ++d)
          for(unsigned e=0; e<dim; ++e)
        {
          cell->get_interpolated_dof_values(ghosted_vector[dim+1+d*dim+e], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = j/ndofs1d;
            int permute = (ax) * ndofs1d + (ndofs1d-1-ay);
            acouele->eleinteriorGradVelnp_((d*dim+e)*dofs_per_cell+permute) = local_values[j];
          }
        }
      }
      else
        dserror("unknown permutation");
      break;
    }
    case DRT::Element::hex8:
    {
      if(deal_vals_loc[0].distance(baci_vals_loc[0])<1e-10 &&
         deal_vals_loc[1].distance(baci_vals_loc[1])<1e-10 &&
         deal_vals_loc[2].distance(baci_vals_loc[3])<1e-10 &&
         deal_vals_loc[3].distance(baci_vals_loc[2])<1e-10 &&
         deal_vals_loc[4].distance(baci_vals_loc[4])<1e-10 &&
         deal_vals_loc[5].distance(baci_vals_loc[5])<1e-10 &&
         deal_vals_loc[6].distance(baci_vals_loc[7])<1e-10 &&
         deal_vals_loc[7].distance(baci_vals_loc[6])<1e-10)
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
      else if(deal_vals_loc[0].distance(baci_vals_loc[4])<1e-10 &&
              deal_vals_loc[1].distance(baci_vals_loc[5])<1e-10 &&
              deal_vals_loc[2].distance(baci_vals_loc[0])<1e-10 &&
              deal_vals_loc[3].distance(baci_vals_loc[1])<1e-10 &&
              deal_vals_loc[4].distance(baci_vals_loc[7])<1e-10 &&
              deal_vals_loc[5].distance(baci_vals_loc[6])<1e-10 &&
              deal_vals_loc[6].distance(baci_vals_loc[3])<1e-10 &&
              deal_vals_loc[7].distance(baci_vals_loc[2])<1e-10)
      {
        for (unsigned int d=0; d<dim; ++d)
        {
          cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = int(j/ndofs1d)%ndofs1d;
            const int az = j/ndofs2d;
            int permute = ax + (ndofs1d-1-ay)*ndofs2d + az*ndofs1d;
            acouele->eleinteriorVelnp_(d*dofs_per_cell+permute) = local_values[j];
          }
        }
        cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
        for (unsigned int j=0; j<dofs_per_cell; ++j)
        {
          const int ax = j%ndofs1d;
          const int ay = int(j/ndofs1d)%ndofs1d;
          const int az = j/ndofs2d;
          int permute = ax + (ndofs1d-1-ay)*ndofs2d + az*ndofs1d;
          acouele->eleinteriorPressnp_(permute) = local_values[j];
        }
      }
      else if(deal_vals_loc[0].distance(baci_vals_loc[5])<1e-10 &&
              deal_vals_loc[1].distance(baci_vals_loc[6])<1e-10 &&
              deal_vals_loc[2].distance(baci_vals_loc[1])<1e-10 &&
              deal_vals_loc[3].distance(baci_vals_loc[2])<1e-10 &&
              deal_vals_loc[4].distance(baci_vals_loc[4])<1e-10 &&
              deal_vals_loc[5].distance(baci_vals_loc[7])<1e-10 &&
              deal_vals_loc[6].distance(baci_vals_loc[0])<1e-10 &&
              deal_vals_loc[7].distance(baci_vals_loc[3])<1e-10)
      {
        for (unsigned int d=0; d<dim; ++d)
        {
          cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = int(j/ndofs1d)%ndofs1d;
            const int az = j/ndofs2d;
            int permute = ax*ndofs1d + (ndofs1d-1-ay)*ndofs2d + (ndofs1d-1-az);
            acouele->eleinteriorVelnp_(d*dofs_per_cell+permute) = local_values[j];
          }
        }
        cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
        for (unsigned int j=0; j<dofs_per_cell; ++j)
        {
          const int ax = j%ndofs1d;
          const int ay = int(j/ndofs1d)%ndofs1d;
          const int az = j/ndofs2d;
          int permute = ax*ndofs1d + (ndofs1d-1-ay)*ndofs2d + (ndofs1d-1-az);
          acouele->eleinteriorPressnp_(permute) = local_values[j];
        }
      }
      else if(deal_vals_loc[0].distance(baci_vals_loc[7])<1e-10 &&
              deal_vals_loc[1].distance(baci_vals_loc[4])<1e-10 &&
              deal_vals_loc[2].distance(baci_vals_loc[3])<1e-10 &&
              deal_vals_loc[3].distance(baci_vals_loc[0])<1e-10 &&
              deal_vals_loc[4].distance(baci_vals_loc[6])<1e-10 &&
              deal_vals_loc[5].distance(baci_vals_loc[5])<1e-10 &&
              deal_vals_loc[6].distance(baci_vals_loc[2])<1e-10 &&
              deal_vals_loc[7].distance(baci_vals_loc[1])<1e-10)
      {
        for (unsigned int d=0; d<dim; ++d)
        {
          cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = int(j/ndofs1d)%ndofs1d;
            const int az = j/ndofs2d;
            int permute = (ndofs1d-1-ax)*ndofs1d + (ndofs1d-1-ay)*ndofs2d + az;
            acouele->eleinteriorVelnp_(d*dofs_per_cell+permute) = local_values[j];
          }
        }
        cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
        for (unsigned int j=0; j<dofs_per_cell; ++j)
        {
          const int ax = j%ndofs1d;
          const int ay = int(j/ndofs1d)%ndofs1d;
          const int az = j/ndofs2d;
          int permute = (ndofs1d-1-ax)*ndofs1d + (ndofs1d-1-ay)*ndofs2d + az;
          acouele->eleinteriorPressnp_(permute) = local_values[j];
        }
      }
      else if(deal_vals_loc[0].distance(baci_vals_loc[6])<1e-10 &&
              deal_vals_loc[1].distance(baci_vals_loc[7])<1e-10 &&
              deal_vals_loc[2].distance(baci_vals_loc[2])<1e-10 &&
              deal_vals_loc[3].distance(baci_vals_loc[3])<1e-10 &&
              deal_vals_loc[4].distance(baci_vals_loc[5])<1e-10 &&
              deal_vals_loc[5].distance(baci_vals_loc[4])<1e-10 &&
              deal_vals_loc[6].distance(baci_vals_loc[1])<1e-10 &&
              deal_vals_loc[7].distance(baci_vals_loc[0])<1e-10)
      {
        for (unsigned int d=0; d<dim; ++d)
        {
          cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = int(j/ndofs1d)%ndofs1d;
            const int az = j/ndofs2d;
            int permute = (ndofs1d-1-ax) + (ndofs1d-1-ay)*ndofs2d + (ndofs1d-1-az)*ndofs1d;
            acouele->eleinteriorVelnp_(d*dofs_per_cell+permute) = local_values[j];
          }
        }
        cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
        for (unsigned int j=0; j<dofs_per_cell; ++j)
        {
          const int ax = j%ndofs1d;
          const int ay = int(j/ndofs1d)%ndofs1d;
          const int az = j/ndofs2d;
          int permute = (ndofs1d-1-ax) + (ndofs1d-1-ay)*ndofs2d + (ndofs1d-1-az)*ndofs1d;
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
  else
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
        dserror("unknown permutation");
      break;
    }
    case DRT::Element::hex8:
    {
      if(deal_vals_loc[0].distance(baci_vals_loc[0])<1e-10 &&
         deal_vals_loc[1].distance(baci_vals_loc[1])<1e-10 &&
         deal_vals_loc[2].distance(baci_vals_loc[3])<1e-10 &&
         deal_vals_loc[3].distance(baci_vals_loc[2])<1e-10 &&
         deal_vals_loc[4].distance(baci_vals_loc[4])<1e-10 &&
         deal_vals_loc[5].distance(baci_vals_loc[5])<1e-10 &&
         deal_vals_loc[6].distance(baci_vals_loc[7])<1e-10 &&
         deal_vals_loc[7].distance(baci_vals_loc[6])<1e-10)
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
      else if(deal_vals_loc[0].distance(baci_vals_loc[4])<1e-10 &&
              deal_vals_loc[1].distance(baci_vals_loc[5])<1e-10 &&
              deal_vals_loc[2].distance(baci_vals_loc[0])<1e-10 &&
              deal_vals_loc[3].distance(baci_vals_loc[1])<1e-10 &&
              deal_vals_loc[4].distance(baci_vals_loc[7])<1e-10 &&
              deal_vals_loc[5].distance(baci_vals_loc[6])<1e-10 &&
              deal_vals_loc[6].distance(baci_vals_loc[3])<1e-10 &&
              deal_vals_loc[7].distance(baci_vals_loc[2])<1e-10)
      {
        for (unsigned int d=0; d<dim; ++d)
        {
          cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = int(j/ndofs1d)%ndofs1d;
            const int az = j/ndofs2d;
            int permute = ax + (ndofs1d-1-ay)*ndofs2d + az*ndofs1d;
            acouele->eleinteriorVelnp_(d*dofs_per_cell+permute) = local_values[j];
          }
        }
        cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
        for (unsigned int j=0; j<dofs_per_cell; ++j)
        {
          const int ax = j%ndofs1d;
          const int ay = int(j/ndofs1d)%ndofs1d;
          const int az = j/ndofs2d;
          int permute = ax + (ndofs1d-1-ay)*ndofs2d + az*ndofs1d;
          acouele->eleinteriorPressnp_(permute) = local_values[j];
        }
      }
      else if(deal_vals_loc[0].distance(baci_vals_loc[5])<1e-10 &&
              deal_vals_loc[1].distance(baci_vals_loc[6])<1e-10 &&
              deal_vals_loc[2].distance(baci_vals_loc[1])<1e-10 &&
              deal_vals_loc[3].distance(baci_vals_loc[2])<1e-10 &&
              deal_vals_loc[4].distance(baci_vals_loc[4])<1e-10 &&
              deal_vals_loc[5].distance(baci_vals_loc[7])<1e-10 &&
              deal_vals_loc[6].distance(baci_vals_loc[0])<1e-10 &&
              deal_vals_loc[7].distance(baci_vals_loc[3])<1e-10)
      {
        for (unsigned int d=0; d<dim; ++d)
        {
          cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = int(j/ndofs1d)%ndofs1d;
            const int az = j/ndofs2d;
            int permute = ax*ndofs1d + (ndofs1d-1-ay)*ndofs2d + (ndofs1d-1-az);
            acouele->eleinteriorVelnp_(d*dofs_per_cell+permute) = local_values[j];
          }
        }
        cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
        for (unsigned int j=0; j<dofs_per_cell; ++j)
        {
          const int ax = j%ndofs1d;
          const int ay = int(j/ndofs1d)%ndofs1d;
          const int az = j/ndofs2d;
          int permute = ax*ndofs1d + (ndofs1d-1-ay)*ndofs2d + (ndofs1d-1-az);
          acouele->eleinteriorPressnp_(permute) = local_values[j];
        }
      }
      else if(deal_vals_loc[0].distance(baci_vals_loc[7])<1e-10 &&
              deal_vals_loc[1].distance(baci_vals_loc[4])<1e-10 &&
              deal_vals_loc[2].distance(baci_vals_loc[3])<1e-10 &&
              deal_vals_loc[3].distance(baci_vals_loc[0])<1e-10 &&
              deal_vals_loc[4].distance(baci_vals_loc[6])<1e-10 &&
              deal_vals_loc[5].distance(baci_vals_loc[5])<1e-10 &&
              deal_vals_loc[6].distance(baci_vals_loc[2])<1e-10 &&
              deal_vals_loc[7].distance(baci_vals_loc[1])<1e-10)
      {
        for (unsigned int d=0; d<dim; ++d)
        {
          cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = int(j/ndofs1d)%ndofs1d;
            const int az = j/ndofs2d;
            int permute = (ndofs1d-1-ax)*ndofs1d + (ndofs1d-1-ay)*ndofs2d + az;
            acouele->eleinteriorVelnp_(d*dofs_per_cell+permute) = local_values[j];
          }
        }
        cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
        for (unsigned int j=0; j<dofs_per_cell; ++j)
        {
          const int ax = j%ndofs1d;
          const int ay = int(j/ndofs1d)%ndofs1d;
          const int az = j/ndofs2d;
          int permute = (ndofs1d-1-ax)*ndofs1d + (ndofs1d-1-ay)*ndofs2d + az;
          acouele->eleinteriorPressnp_(permute) = local_values[j];
        }
      }
      else if(deal_vals_loc[0].distance(baci_vals_loc[6])<1e-10 &&
              deal_vals_loc[1].distance(baci_vals_loc[7])<1e-10 &&
              deal_vals_loc[2].distance(baci_vals_loc[2])<1e-10 &&
              deal_vals_loc[3].distance(baci_vals_loc[3])<1e-10 &&
              deal_vals_loc[4].distance(baci_vals_loc[5])<1e-10 &&
              deal_vals_loc[5].distance(baci_vals_loc[4])<1e-10 &&
              deal_vals_loc[6].distance(baci_vals_loc[1])<1e-10 &&
              deal_vals_loc[7].distance(baci_vals_loc[0])<1e-10)
      {
        for (unsigned int d=0; d<dim; ++d)
        {
          cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = int(j/ndofs1d)%ndofs1d;
            const int az = j/ndofs2d;
            int permute = (ndofs1d-1-ax) + (ndofs1d-1-ay)*ndofs2d + (ndofs1d-1-az)*ndofs1d;
            acouele->eleinteriorVelnp_(d*dofs_per_cell+permute) = local_values[j];
          }
        }
        cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
        for (unsigned int j=0; j<dofs_per_cell; ++j)
        {
          const int ax = j%ndofs1d;
          const int ay = int(j/ndofs1d)%ndofs1d;
          const int az = j/ndofs2d;
          int permute = (ndofs1d-1-ax) + (ndofs1d-1-ay)*ndofs2d + (ndofs1d-1-az)*ndofs1d;
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

template<int dim, int fe_degree>
void
WaveEquationOperation<dim,fe_degree>::
compute_gradient_contributions(std::vector<parallel::distributed::Vector<value_type> > &fwnp,
                               std::vector<parallel::distributed::Vector<value_type> > &fwn,
                               std::vector<parallel::distributed::Vector<value_type> > &adnp)
{
  // we need the derivative of the mass matrix with respect to density and sound speed and have to buuild the scalar product with
  // corresponding pressure and velocity values
  FEEvaluation<dim,fe_degree,fe_degree+1,dim,value_type> adjoint_velocity(data);
  FEEvaluation<dim,fe_degree,fe_degree+1,1,value_type>   adjoint_pressure(data);

  FEEvaluation<dim,fe_degree,fe_degree+1,dim,value_type> forward_velocity(data);
  FEEvaluation<dim,fe_degree,fe_degree+1,1,value_type>   forward_pressure(data);

  // build difference vector
  std::vector<parallel::distributed::Vector<value_type> > fw_diff;
  fw_diff.resize(fwnp.size());
  for(unsigned int i=0; i<fwnp.size(); ++i)
  {
    fw_diff[i] = fwnp[i];
    fw_diff[i] -= fwn[i];
  }

  // build adjoint work vector
  std::vector<parallel::distributed::Vector<value_type> > ad;
  ad.resize(adnp.size());
  for(unsigned int i=0; i<adnp.size(); ++i)
    ad[i] = adnp[i];

  for (unsigned int cell=0; cell<data.n_macro_cells()+data.n_macro_ghost_cells(); ++cell)
  {
    // read adjoint solution
    adjoint_velocity.reinit(cell);
    adjoint_velocity.read_dof_values(ad, 0);
    adjoint_velocity.evaluate (true, false, false);

    adjoint_pressure.reinit(cell);
    adjoint_pressure.read_dof_values(ad, dim);
    adjoint_pressure.evaluate(true, false, false);

    // sort the correspondent values
    for (unsigned int q=0; q<adjoint_velocity.n_q_points; ++q)
    {
      const Tensor<1,dim,VectorizedArray<value_type> > adjoint_velocity_value = adjoint_velocity.get_value(q);
      const VectorizedArray<value_type> adjoint_pressure_value = adjoint_pressure.get_value(q);
      adjoint_pressure.submit_value(adjoint_pressure_value,q);
      adjoint_velocity.submit_value(adjoint_velocity_value,q);
    }

    // do integration of adjoint solution
    adjoint_velocity.integrate (true, false);
    adjoint_velocity.distribute_local_to_global (ad, 0);
    adjoint_pressure.integrate(true, false);
    adjoint_pressure.distribute_local_to_global (ad, dim);

    // get the dof values of the forward solutions
    forward_velocity.reinit(cell);
    forward_velocity.read_dof_values(fw_diff,0);
    forward_pressure.reinit(cell);
    forward_pressure.read_dof_values(fw_diff,dim);

    // get the material values
    const VectorizedArray<value_type> rho_fac = 1./densities[cell]; //-1./densities[cell]/densities[cell]/speeds[cell]/speeds[cell];
    const VectorizedArray<value_type> c_fac = 2./speeds[cell]; //-2./speeds[cell]/speeds[cell]/speeds[cell]/densities[cell];

    // get the dof values of integrated adjoint solution and multiply with correspondent forward solutions
    VectorizedArray<value_type> pressure_mass_mat_contrib = VectorizedArray<value_type>();
    for(unsigned int dof=0; dof<adjoint_pressure.dofs_per_cell; ++dof)
    {
      for(unsigned int d=0; d<dim; ++d)
        densities_grad[cell] += 1./densities[cell] * adjoint_velocity.get_dof_value(dof)[d] * forward_velocity.get_dof_value(dof)[d]; // factor is 1
      pressure_mass_mat_contrib += adjoint_pressure.get_dof_value(dof) * forward_pressure.get_dof_value(dof);
    }

    speeds_grad[cell] -= c_fac * pressure_mass_mat_contrib;
    densities_grad[cell] -= rho_fac * pressure_mass_mat_contrib;
  }

  return;
}

template<int dim, int fe_degree>
void WaveEquationOperation<dim,fe_degree>::
write_gradient_contributions(Teuchos::RCP<DRT::DiscretizationHDG> &discret, value_type dt) const
{
  for (unsigned int i=0; i<data.n_macro_cells()+data.n_macro_ghost_cells(); ++i)
  {
    for (unsigned int v=0; v<data.n_components_filled(i); ++v)
    {
      if (data.get_cell_iterator(i,v)->level() != 0)
        dserror("Refined meshes currently not implemented!");

      const int element_index = data.get_cell_iterator(i,v)->index();
      DRT::ELEMENTS::Acou * acouele = dynamic_cast<DRT::ELEMENTS::Acou*>(discret->lColElement(element_index));

      acouele->AddToDensityGradient((densities_grad[i][v])/dt);
      acouele->AddToSoSGradient((speeds_grad[i][v])/dt);
    }
  }
  return;
}

template<int dim, int fe_degree>
void WaveEquationOperation<dim, fe_degree>::
local_apply_domain(const MatrixFree<dim,value_type>                                &data,
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

      pressure.submit_value(c_sq*rhs,q);
      if(this->adjoint_eval==false)
      {
        velocity.submit_value(-rho_inv*pressure_gradient,q);
        pressure.submit_gradient(rho*c_sq*velocity_value,q);
      }
      else
      {
        velocity.submit_value(rho*c_sq*pressure_gradient,q);
        pressure.submit_gradient(-rho_inv*velocity_value,q);
      }
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
  const std::pair<unsigned int,unsigned int>                     &face_range) const
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
    const VectorizedArray<value_type> tau_plus = 1./c_plus/rho_plus;

    phi_neighbor.reinit(face);
    phi_neighbor.read_dof_values(src, 0);
    phi_neighbor.evaluate(true,false);
    const VectorizedArray<value_type> rho_minus = phi_neighbor.read_cell_data(densities);
    const VectorizedArray<value_type> rho_inv_minus = 1./rho_minus;
    const VectorizedArray<value_type> c_minus = phi_neighbor.read_cell_data(speeds);
    const VectorizedArray<value_type> c_sq_minus = c_minus * c_minus;
    const VectorizedArray<value_type> tau_minus = 1./c_minus/rho_minus;

    const VectorizedArray<value_type> tau_inv = 1./(tau_plus + tau_minus);

    AssertDimension(phi.n_q_points, data.get_n_q_points_face(0));

    for (unsigned int q=0; q<phi.n_q_points; ++q)
    {
      Tensor<1,dim+1,VectorizedArray<value_type> > val_plus = phi.get_value(q);
      Tensor<1,dim+1,VectorizedArray<value_type> > val_minus = phi_neighbor.get_value(q);
      Tensor<1,dim,VectorizedArray<value_type> > normal = phi.get_normal_vector(q);
      VectorizedArray<value_type> normal_v_plus = val_plus[0] *normal[0];
      VectorizedArray<value_type> normal_v_minus = -val_minus[0]*normal[0];
      for (unsigned int d=1; d<dim; ++d)
      {
        normal_v_plus += val_plus[d] * normal[d];
        normal_v_minus -= val_minus[d] * normal[d];
      }

      VectorizedArray<value_type> lambda;
      VectorizedArray<value_type> pres_diff_plus;
      VectorizedArray<value_type> pres_diff_minus;
      if(this->adjoint_eval==false)
      {
        lambda = tau_inv*(normal_v_plus + normal_v_minus + tau_plus*val_plus[dim] + tau_minus*val_minus[dim]);
        pres_diff_plus  = (val_plus[dim] - lambda)  * rho_inv_plus;
        pres_diff_minus = (val_minus[dim] - lambda) * rho_inv_minus;
      }
      else
      {
        lambda = tau_inv*(rho_inv_plus*normal_v_plus + rho_inv_minus*normal_v_minus - tau_plus*rho_plus*c_sq_plus*val_plus[dim] - tau_minus*rho_minus*c_sq_minus*val_minus[dim]);
        pres_diff_plus  = -rho_plus*c_sq_plus*val_plus[dim] - lambda ;
        pres_diff_minus = -rho_minus*c_sq_minus*val_minus[dim] - lambda;
      }

      for (unsigned int d=0; d<dim; ++d)
      {
        val_plus[d] = pres_diff_plus*normal[d];
        val_minus[d] = -pres_diff_minus*normal[d];
      }
      if(this->adjoint_eval==false)
      {
        val_plus[dim] = -c_sq_plus * rho_plus * (normal_v_plus + tau_plus * (val_plus[dim] - lambda));
        val_minus[dim] = -c_sq_minus * rho_minus * (normal_v_minus + tau_minus * (val_minus[dim] - lambda));
      }
      else
      {
        val_plus[dim] = -(-rho_inv_plus*normal_v_plus + tau_plus * (c_sq_plus*rho_plus*val_plus[dim] + lambda));
        val_minus[dim] = -(-rho_inv_minus*normal_v_minus + tau_minus * (c_sq_minus*rho_minus*val_minus[dim] + lambda));
      }
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

  // quantities we need in the loop
  Point<dim> point;
  std::vector<value_type> node_values;
  std::vector<std::vector<value_type> > node_coords;
  node_coords.resize(GeometryInfo<dim>::vertices_per_face);
  node_values.resize(GeometryInfo<dim>::vertices_per_face);
  for(unsigned int n=0; n<GeometryInfo<dim>::vertices_per_face; ++n)
    node_coords[n].resize(dim);

  for (unsigned int face=face_range.first; face<face_range.second; face++)
  {
    phi.reinit(face);
    phi.read_dof_values(src, 0);
    phi.evaluate(true,false);
    const VectorizedArray<value_type> rho = phi.read_cell_data(densities);
    const VectorizedArray<value_type> rho_inv = 1./rho;
    const VectorizedArray<value_type> c_sq = phi.read_cell_data(speeds)*phi.read_cell_data(speeds);
    const VectorizedArray<value_type> c = phi.read_cell_data(speeds);
    const VectorizedArray<value_type> tau = 1./phi.read_cell_data(speeds)/phi.read_cell_data(densities);

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
        if(this->adjoint_eval==false)
          lambda = tau/(tau+1./c/rho)*p_plus + 1./(tau+1./c/rho)*normal_v_plus;
        else
          lambda = 1./(tau+1./c/rho)*rho_inv*normal_v_plus - tau*rho*c_sq/(tau+1./c/rho)*p_plus;
      else if(int_boundary_id==1) // monitored
      {
        lambda = 1./tau*normal_v_plus+p_plus;//VectorizedArray<value_type>();

        if(source_adjoint_meas!=Teuchos::null && this->adjoint_eval == true) // second query required for intermediate integration
        {
          for (unsigned int v=0; v<VectorizedArray<value_type>::n_array_elements && data.faces[face].left_cell[v] != numbers::invalid_unsigned_int; ++v)
          {
            for (unsigned int d=0; d<dim; ++d)
              point[d] = q_point[d][v];

            for(unsigned int n=0; n<GeometryInfo<dim>::vertices_per_face; ++n)
            {
              for(unsigned int d=0; d<dim; ++d)
                node_coords[n][d] = table_node_coords(face-data.n_macro_inner_faces(),v,n,d);
              int gid = table_node_ids(face-data.n_macro_inner_faces(),v,n);
              int lid = source_adjoint_meas->Map().LID(gid);
              node_values[n] =  source_adjoint_meas->operator ()(this->timestep_source_number)->operator [](lid);
            }
            lambda[v] -= 1.0/tau[v] * this->evaluate_source_adjoint(point,node_coords,node_values);
          }
        }
      }
      else if(int_boundary_id==2) // monitored and absorbing
      {
        if(this->adjoint_eval==false)
          lambda = tau/(tau+1./c/rho)*p_plus + 1./(tau+1./c/rho)*normal_v_plus;
        else
          lambda = 1./(tau+1./c/rho)*rho_inv*normal_v_plus - tau*rho*c_sq/(tau+1./c/rho)*p_plus;
        if(source_adjoint_meas!=Teuchos::null && this->adjoint_eval == true) // second query required for intermediate integration
        {
          for (unsigned int v=0; v<VectorizedArray<value_type>::n_array_elements && data.faces[face].left_cell[v] != numbers::invalid_unsigned_int; ++v)
          {
            for (unsigned int d=0; d<dim; ++d)
              point[d] = q_point[d][v];

            for(unsigned int n=0; n<GeometryInfo<dim>::vertices_per_face; ++n)
            {
              for(unsigned int d=0; d<dim; ++d)
                node_coords[n][d] = table_node_coords(face-data.n_macro_inner_faces(),v,n,d);
              int gid = table_node_ids(face-data.n_macro_inner_faces(),v,n);
              int lid = source_adjoint_meas->Map().LID(gid);
              node_values[n] =  source_adjoint_meas->operator ()(this->timestep_source_number)->operator [](lid);
            }
            lambda[v] -= 1.0/(tau[v]+1./c[v]/rho[v]) * this->evaluate_source_adjoint(point,node_coords,node_values);
          }
        }
      }
      else if(int_boundary_id==3) // free boundary
        lambda = 1./tau*normal_v_plus+p_plus; // VectorizedArray<value_type>();
      else if(int_boundary_id==4) // dbc from time reversal
      {
        if(source_adjoint_meas!=Teuchos::null)
        {
          for (unsigned int v=0; v<VectorizedArray<value_type>::n_array_elements && data.faces[face].left_cell[v] != numbers::invalid_unsigned_int; ++v)
          {
            for (unsigned int d=0; d<dim; ++d)
              point[d] = q_point[d][v];

            for(unsigned int n=0; n<GeometryInfo<dim>::vertices_per_face; ++n)
            {
              for(unsigned int d=0; d<dim; ++d)
                node_coords[n][d] = table_node_coords(face-data.n_macro_inner_faces(),v,n,d);
              int gid = table_node_ids(face-data.n_macro_inner_faces(),v,n);
              int lid = source_adjoint_meas->Map().LID(gid);
              node_values[n] =  source_adjoint_meas->operator ()(this->timestep_source_number)->operator [](lid);
            }
            lambda[v] = this->evaluate_source_timereversal(point,node_coords,node_values);
          }
        }
      }
      else if(int_boundary_id>=5) // dbcs
      {
        if(this->adjoint_eval==false)
          for (unsigned int v=0; v<VectorizedArray<value_type>::n_array_elements; ++v)
          {
            Point<dim> point;
            for (unsigned int d=0; d<dim; ++d)
              point[d] = q_point[d][v];
            lambda[v] = dirichlet_boundary_conditions->value(point,(int_boundary_id-5)*dim);
          }
        else
          lambda = VectorizedArray<value_type>();
      }

      if(this->adjoint_eval==false)
      {
        for (unsigned int d=0; d<dim; ++d)
          val_plus[d] = (p_plus - lambda)*normal[d]*rho_inv;
        val_plus[dim] = -c_sq*rho*(normal_v_plus - tau*(lambda - p_plus));
      }
      else
      {
        for (unsigned int d=0; d<dim; ++d)
          val_plus[d] = -(lambda + rho*c_sq*p_plus)*normal[d];
        val_plus[dim] = -(-rho_inv*normal_v_plus + tau*(lambda + rho*c_sq*p_plus));
      }
      phi.submit_value(val_plus,q);
    }
    phi.integrate(true,false);
    phi.distribute_local_to_global(dst, 0);
  }
}

template<int dim, int fe_degree>
void WaveEquationOperation<dim, fe_degree>::
local_apply_solid_domain(const MatrixFree<dim,value_type>                                &data,
                         std::vector<parallel::distributed::Vector<value_type> >         &dst,
                         const std::vector<parallel::distributed::Vector<value_type> >   &src,
                         const std::pair<unsigned int,unsigned int>                 &cell_range) const
{
  FEEvaluation<dim,fe_degree,fe_degree+1,dim,value_type> v(data);
  FEEvaluation<dim,fe_degree,fe_degree+1,1,value_type> p(data);
  FEEvaluation<dim,fe_degree,fe_degree+1,dim*dim,value_type> H(data);

  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
  {
    // It is faster to evaluate values of the vector-valued velocity and
    // gradients of the scalar pressure than divergence of velocity and
    // values of pressure
    v.reinit(cell);
    v.read_dof_values(src, 0);
    v.evaluate (true, true, false);

    p.reinit(cell);
    p.read_dof_values(src, dim);
    p.evaluate(true, false, false);

    H.reinit(cell);
    H.read_dof_values(src,dim+1);
    H.evaluate(true, false, false);

    const VectorizedArray<value_type> rho = densities[cell];
    const VectorizedArray<value_type> rho_inv = 1./densities[cell];
    const VectorizedArray<value_type> c_sq = speeds[cell]*speeds[cell];
    const VectorizedArray<value_type> visc = viscs[cell];

    for (unsigned int q=0; q<v.n_q_points; ++q)
    {
      const VectorizedArray<value_type>                p_value = p.get_value(q);
      const Tensor<1,dim,VectorizedArray<value_type> > v_value = v.get_value(q);
      const Tensor<2,dim,VectorizedArray<value_type> > v_gradient = v.get_gradient(q);
      const Tensor<1,dim*dim,VectorizedArray<value_type> > H_value = H.get_value(q);

      Point<dim,VectorizedArray<value_type> > q_points = v.quadrature_point(q);
      Tensor<1,dim,VectorizedArray<value_type> > rhs;
      for (unsigned int d=0; d<dim; ++d)
        for (unsigned int n=0; n<rhs[d].n_array_elements; ++n)
        {
          Point<dim> q_point;
          for (unsigned int d=0; d<dim; ++d)
            q_point[d] = q_points[d][n];
          rhs[d][n] = source_term->value(q_point,d);
        }
      v.submit_value(rho_inv*rhs,q);

      Tensor<2,dim,VectorizedArray<value_type> > muHminp;
      for(unsigned int m=0; m<dim; ++m)
        for(unsigned int n=0; n<dim; ++n)
          if(m==n)
            muHminp[n][n] = visc*H_value[n+dim*m]+p_value;
          else
            muHminp[n][m] = visc*H_value[m+dim*n];

      // contribution to H
      Tensor<1,dim*dim,VectorizedArray<value_type> > help;
      for(unsigned int m=0; m<dim; ++m)
        for(unsigned int n=0; n<dim; ++n)
          help[m+dim*n] = v_gradient[n][m];

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


      if(this->adjoint_eval==false)
      {
        //H.submit_gradient(-contribH,q);
        H.submit_value(help,q);
        v.submit_gradient(-rho_inv*muHminp,q);
        p.submit_gradient(-rho*c_sq*v_value,q);
      }
      else
      {
        dserror("think about it");
      }
    }

    v.integrate (true, true);
    v.distribute_local_to_global (dst, 0);

    p.integrate(false, true);
    p.distribute_local_to_global (dst,dim);

    H.integrate(true, false);
    H.distribute_local_to_global (dst,dim+1);
  }
}



template <int dim, int fe_degree>
void
WaveEquationOperation<dim,fe_degree>::
local_apply_solid_face (const MatrixFree<dim,value_type> &,
                        std::vector<parallel::distributed::Vector<value_type> > &dst,
                        const std::vector<parallel::distributed::Vector<value_type> > &src,
                        const std::pair<unsigned int,unsigned int> &face_range) const
{
  // There is some overhead in the methods in FEEvaluation, so it is faster
  // to combine pressure and velocity in the same object and just combine
  // them at the level of quadrature points
  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim*dim+dim+1,value_type> phi(this->data, true, 0, 0, true);
  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim*dim+dim+1,value_type> phi_neighbor(this->data, false, 0, 0, true);

  for (unsigned int face=face_range.first; face<face_range.second; face++)
  {
    phi.reinit(face);
    phi.read_dof_values(src, 0);
    phi.evaluate(true,false);
    const VectorizedArray<value_type> rho_plus = phi.read_cell_data(densities);
    const VectorizedArray<value_type> rho_inv_plus = 1./rho_plus;
    const VectorizedArray<value_type> c_plus = phi.read_cell_data(speeds);
    const VectorizedArray<value_type> c_sq_plus = c_plus * c_plus;
    const VectorizedArray<value_type> tau_plus = make_vectorized_array(1.0);
    const VectorizedArray<value_type> visc_plus = phi.read_cell_data(viscs);

    phi_neighbor.reinit(face);
    phi_neighbor.read_dof_values(src, 0);
    phi_neighbor.evaluate(true,false);
    const VectorizedArray<value_type> rho_minus = phi_neighbor.read_cell_data(densities);
    const VectorizedArray<value_type> rho_inv_minus = 1./rho_minus;
    const VectorizedArray<value_type> c_minus = phi_neighbor.read_cell_data(speeds);
    const VectorizedArray<value_type> c_sq_minus = c_minus * c_minus;
    const VectorizedArray<value_type> tau_minus = make_vectorized_array(1.0);
    const VectorizedArray<value_type> visc_minus = phi_neighbor.read_cell_data(viscs);

    const VectorizedArray<value_type> tau_inv = 1./(tau_plus + tau_minus);

    AssertDimension(phi.n_q_points, data.get_n_q_points_face(0));

    for(unsigned int q=0; q<phi.n_q_points; ++q)
    {
      Tensor<1,dim*dim+dim+1,VectorizedArray<value_type> > val_plus = phi.get_value(q);
      Tensor<1,dim*dim+dim+1,VectorizedArray<value_type> > val_minus = phi_neighbor.get_value(q);
      Tensor<1,dim,VectorizedArray<value_type> > normal = phi.get_normal_vector(q);

      Tensor<1,dim,VectorizedArray<value_type> > vhat;
      for(unsigned int d=0; d<dim; ++d)
      {
        vhat[d] = 1./2.*(val_plus[d] + val_minus[d]) - tau_inv*(val_plus[dim]-val_minus[dim])*normal[d] ;
        for(unsigned int e=0; e<dim; ++e)
          vhat[d] -= tau_inv * (visc_plus*val_plus[dim+1+d*dim+e]-visc_minus*val_minus[dim+1+d*dim+e])*normal[e];
      }

      // vhat*normal+
      VectorizedArray<value_type> normal_vhat = vhat[0]*normal[0];
      for(unsigned int d=1; d<dim; ++d)
        normal_vhat += vhat[d]*normal[d];

      Tensor<1,dim*dim+dim+1,VectorizedArray<value_type> > submit_plus;
      Tensor<1,dim*dim+dim+1,VectorizedArray<value_type> > submit_minus;

      // 1. contribution to pressure
      submit_plus[dim] = c_sq_plus*rho_plus*normal_vhat;
      submit_minus[dim] = -c_sq_minus*rho_minus*normal_vhat;

      // 2. contribution to velocity
      for(unsigned int d=0; d<dim; ++d)
      {
        submit_plus[d] = -rho_inv_plus * tau_plus * (val_plus[d]-vhat[d]);
        submit_minus[d] = -rho_inv_minus * tau_minus * (val_minus[d]-vhat[d]);
        for(unsigned int e=0; e<dim; ++e)
          if(e==d)
          {
            submit_plus[d] += rho_inv_plus* (visc_plus*val_plus[dim+1+e*dim+d]+val_plus[dim])*normal[e];
            submit_minus[d] -= rho_inv_minus* (visc_minus*val_minus[dim+1+e*dim+d]+val_minus[dim])*normal[e];
          }
          else
          {
            submit_plus[d] += rho_inv_plus* (visc_plus*val_plus[dim+1+d*dim+e])*normal[e];
            submit_minus[d] -= rho_inv_minus* (visc_minus*val_minus[dim+1+d*dim+e])*normal[e];
          }
      }

      // 3. contribution to H
      for(unsigned int d=0; d<dim; ++d)
        for(unsigned int e=0; e<dim; ++e)
        {
          submit_plus[dim+1+d*dim+e] -= (val_plus[d]-vhat[d])*normal[e];
          submit_minus[dim+1+d*dim+e] += (val_minus[d]-vhat[d])*normal[e];
          //submit_plus[dim+1+d*dim+e] += vhat[d]*normal[e];
          //submit_minus[dim+1+d*dim+e] -= vhat[d]*normal[e];
        }


      if(this->adjoint_eval)
        dserror("think about adjoint problem");

      phi.submit_value(submit_plus, q);
      phi_neighbor.submit_value(submit_minus, q);
    }

    phi.integrate(true,false);
    phi.distribute_local_to_global(dst, 0);

    phi_neighbor.integrate(true,false);
    phi_neighbor.distribute_local_to_global(dst, 0);
  }
}

template <int dim, int fe_degree>
void WaveEquationOperation<dim,fe_degree>::
local_apply_solid_boundary_face (const MatrixFree<dim,value_type> &,
                           std::vector<parallel::distributed::Vector<value_type> >       &dst,
                           const std::vector<parallel::distributed::Vector<value_type> > &src,
                           const std::pair<unsigned int,unsigned int>               &face_range) const
{
  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim*dim+dim+1,value_type> phi(this->data, true, 0, 0, true);

  // quantities we need in the loop
  Point<dim> point;
  std::vector<value_type> node_values;
  std::vector<std::vector<value_type> > node_coords;
  node_coords.resize(GeometryInfo<dim>::vertices_per_face);
  node_values.resize(GeometryInfo<dim>::vertices_per_face);
  for(unsigned int n=0; n<GeometryInfo<dim>::vertices_per_face; ++n)
    node_coords[n].resize(dim);

  for (unsigned int face=face_range.first; face<face_range.second; face++)
  {
    phi.reinit(face);
    phi.read_dof_values(src, 0);
    phi.evaluate(true,false);
    const VectorizedArray<value_type> rho = phi.read_cell_data(densities);
    const VectorizedArray<value_type> rho_inv = 1./rho;
    const VectorizedArray<value_type> c_sq = phi.read_cell_data(speeds)*phi.read_cell_data(speeds);
    const VectorizedArray<value_type> c = phi.read_cell_data(speeds);
    const VectorizedArray<value_type> tau = make_vectorized_array(1.0);
    const VectorizedArray<value_type> visc = phi.read_cell_data(viscs);

    const types::boundary_id boundary_index = this->data.get_boundary_indicator(face);
    const int int_boundary_id = int(boundary_index);

    for (unsigned int q=0; q<phi.n_q_points; ++q)
    {
      Tensor<1,dim,VectorizedArray<value_type> > normal = phi.get_normal_vector(q);
      Tensor<1,dim*dim+dim+1,VectorizedArray<value_type> > val_plus = phi.get_value(q);
      Point<dim,VectorizedArray<value_type> > q_point = phi.quadrature_point(q);

      // calculation of vhat dependent on boundary type
      Tensor<1,dim,VectorizedArray<value_type> > vhat;
      if(int_boundary_id==0) // absorbing boundary
      {
        for(unsigned int d=0; d<dim; ++d)
        {
          vhat[d] = tau/(rho*c+tau)*val_plus[d] - 1./(tau+rho*c)*val_plus[dim]*normal[d] ;
          for(unsigned int e=0; e<dim; ++e)
            vhat[d] -= 1./(tau+rho*c) * visc * val_plus[dim+1+d*dim+e] * normal[e];
        }
      }
      else if(int_boundary_id==1) // monitored
      {
        if(this->adjoint_eval==false)
          for(unsigned int d=0; d<dim; ++d)
          {
            vhat[d] = val_plus[d] - 1./tau*val_plus[dim]*normal[d] ;
            for(unsigned int e=0; e<dim; ++e)
              vhat[d] -= 1./tau * visc * val_plus[dim+1+d*dim+e] * normal[e];
          }
        else
          dserror("todo");
      }
      else if(int_boundary_id==2) // monitored and absorbing
      {
        if(this->adjoint_eval==false)
          for(unsigned int d=0; d<dim; ++d)
          {
            vhat[d] = tau/(rho*c+tau)*val_plus[d] - 1./(tau+rho*c)*val_plus[dim]*normal[d] ;
            for(unsigned int e=0; e<dim; ++e)
              vhat[d] -= 1./(tau+rho*c) * visc * val_plus[dim+1+d*dim+e] * normal[e];
          }
        else
          dserror("todo");
      }
      else if(int_boundary_id==3) // free boundary
        for(unsigned int d=0; d<dim; ++d)
        {
          vhat[d] = val_plus[d] - 1./tau*val_plus[dim]*normal[d] ;
          for(unsigned int e=0; e<dim; ++e)
            vhat[d] -= 1./tau * visc * val_plus[dim+1+d*dim+e] * normal[e];
        }
      else if(int_boundary_id==4) // dbc from time reversal
      {
        dserror("todo");
      }
      else if(int_boundary_id>=5) // dbcs
      {
        if(this->adjoint_eval==false)
          for (unsigned int v=0; v<VectorizedArray<value_type>::n_array_elements; ++v)
          {
            Point<dim> point;
            for (unsigned int d=0; d<dim; ++d)
              point[d] = q_point[d][v];
            for(unsigned int d=0; d<dim; ++d)
              vhat[d][v] = dirichlet_boundary_conditions->value(point,(int_boundary_id-5)*dim+d);
          }
        else
          for(unsigned int d=0; d<dim; ++d)
            vhat[d] = VectorizedArray<value_type>();
      }

      // now the actual boundary term evaluation
      Tensor<1,dim*dim+dim+1,VectorizedArray<value_type> > submit_plus;
      if(this->adjoint_eval==false)
      {
        // vhat*normal
        VectorizedArray<value_type> normal_vhat = vhat[0]*normal[0];
        for(unsigned int d=1; d<dim; ++d)
          normal_vhat += vhat[d]*normal[d];

        // 1. contribution to pressure
        submit_plus[dim] = c_sq*rho*normal_vhat;

        // 2. contribution to velocity
        for(unsigned int d=0; d<dim; ++d)
        {
          submit_plus[d] = -rho_inv * tau * (val_plus[d]-vhat[d]);
          for(unsigned int e=0; e<dim; ++e)
            if(e==d)
              submit_plus[d] += rho_inv* (visc*val_plus[dim+1+e*dim+d]+val_plus[dim])*normal[e];
            else
              submit_plus[d] += rho_inv* (visc*val_plus[dim+1+d*dim+e])*normal[e];
        }

        // 3. contribution to H
        for(unsigned int d=0; d<dim; ++d)
          for(unsigned int e=0; e<dim; ++e)
          {
            //submit_plus[dim+1+d*dim+e] += vhat[d]*normal[e];
            submit_plus[dim+1+d*dim+e] -= (val_plus[d]-vhat[d])*normal[e];
          }
      }
      else
      {
        dserror("todo");
      }

      phi.submit_value(submit_plus,q);
    }
    phi.integrate(true,false);
    phi.distribute_local_to_global(dst, 0);
  }
}

template <int dim, int fe_degree>
typename WaveEquationOperationBase<dim>::value_type WaveEquationOperation<dim,fe_degree>::evaluate_source_adjoint(const Point<dim> &p, const std::vector<std::vector<value_type> > nodes, std::vector<value_type> values) const
{
  value_type result = 0.0;
  value_type xyz[dim];
  for(unsigned d=0; d<dim; ++d)
    xyz[d] = p(d);

  if(dim==2 && nodes.size()==2) // quad4 with line2 face element
  {
    value_type node_distance = std::sqrt( (nodes[0][0]-nodes[1][0])*(nodes[0][0]-nodes[1][0])
                                        + (nodes[0][1]-nodes[1][1])*(nodes[0][1]-nodes[1][1]) );
    value_type quad_distance = std::sqrt( (xyz[0]-nodes[0][0])*(xyz[0]-nodes[0][0])
                                        + (xyz[1]-nodes[0][1])*(xyz[1]-nodes[0][1]) );

    result = quad_distance/node_distance * values[1] + (node_distance-quad_distance)/node_distance * values[0];
    result *= -1.0;
    result *= 2.0/node_distance;

    if(node_distance<quad_distance)
    {
      std::cout<<"punkt "<<xyz[0]<<" "<<xyz[1]<<" nodecoords "<<nodes[0][0]<<" "<<nodes[0][1]
                                              <<" "           <<nodes[1][0]<<" "<<nodes[1][1]<<" "
                                              <<" nodevals  " <<values[0]  <<" "<<values[1]<<std::endl;
      std::cout<<"wrong face!!!"<<std::endl<<std::endl;
    }
  }
  else
    dserror("not yet implemented");

  return result;
}

template <int dim, int fe_degree>
typename WaveEquationOperationBase<dim>::value_type WaveEquationOperation<dim,fe_degree>::evaluate_source_timereversal(const Point<dim> &p, const std::vector<std::vector<value_type> > nodes, std::vector<value_type> values) const
{
  value_type result = 0.0;
  value_type xyz[dim];
  for(unsigned d=0; d<dim; ++d)
    xyz[d] = p(d);

  if(dim==2 && nodes.size()==2) // quad4 with line2 face element
  {
    value_type node_distance = std::sqrt( (nodes[0][0]-nodes[1][0])*(nodes[0][0]-nodes[1][0])
                                        + (nodes[0][1]-nodes[1][1])*(nodes[0][1]-nodes[1][1]) );
    value_type quad_distance = std::sqrt( (xyz[0]-nodes[0][0])*(xyz[0]-nodes[0][0])
                                        + (xyz[1]-nodes[0][1])*(xyz[1]-nodes[0][1]) );

    result = quad_distance/node_distance * values[1] + (node_distance-quad_distance)/node_distance * values[0];

    if(node_distance<quad_distance)
    {
      std::cout<<"punkt "<<xyz[0]<<" "<<xyz[1]<<" nodecoords "<<nodes[0][0]<<" "<<nodes[0][1]
                                              <<" "           <<nodes[1][0]<<" "<<nodes[1][1]<<" "
                                              <<" nodevals  " <<values[0]  <<" "<<values[1]<<std::endl;
      std::cout<<"wrong face!!!"<<std::endl<<std::endl;
    }
  }
  else
    dserror("not yet implemented");

  return result;
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
local_apply_solid_mass_matrix(const MatrixFree<dim,value_type>                  &data,
                        std::vector<parallel::distributed::Vector<value_type> >        &dst,
                        const std::vector<parallel::distributed::Vector<value_type> >  &src,
                        const std::pair<unsigned int,unsigned int>    &cell_range) const
{
  internal::InverseMassMatrixDataSolid<dim,fe_degree,value_type>& mass_data = mass_matrix_data_solid.get();
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
    {
      mass_data.phi[0].reinit(cell);
      mass_data.phi[0].read_dof_values(src, 0);

      mass_data.inverse.fill_inverse_JxW_values(mass_data.coefficients);
      mass_data.inverse.apply(mass_data.coefficients, dim*dim+dim+1,
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

  if(this->solid)
    data.loop (&WaveEquationOperation<dim, fe_degree>::local_apply_solid_domain,
               &WaveEquationOperation<dim, fe_degree>::local_apply_solid_face,
               &WaveEquationOperation<dim, fe_degree>::local_apply_solid_boundary_face,
               this, dst, src);
  else
    data.loop (&WaveEquationOperation<dim, fe_degree>::local_apply_domain,
               &WaveEquationOperation<dim, fe_degree>::local_apply_face,
               &WaveEquationOperation<dim, fe_degree>::local_apply_boundary_face,
               this, dst, src);

  computing_times[0] += timer.wall_time();
  timer.restart();

  if(this->solid)
    data.cell_loop(&WaveEquationOperation<dim, fe_degree>::local_apply_solid_mass_matrix,
                 this, dst, dst);
  else
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
