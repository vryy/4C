/*!----------------------------------------------------------------------
\file acou_expl_worker.cpp
\brief Control routine for acoustic explicit time integration.
\level 2

<pre>
\level 2

\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/

#include "acou_expl_worker.H"

#ifdef HAVE_DEAL_II

#include <deal.II/lac/parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>
//#include "fe_evaluation.h"
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/task_info.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/base/timer.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <Epetra_MpiComm.h>

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/acoustic.H"
#include "../drt_mat/acoustic_sol.H"
#include "acou_ele.H"
#include "acou_sol_ele.H"
#include "acou_pml.H"
#include "pat_utils.H"


#define ADERLTS

namespace ACOU
{
namespace internal
{
  template <int dim, typename Number>
  MatrixFree<dim,Number>
  create_matrix_free(const DoFHandler<dim> &dof_handler,
                     const unsigned int     fe_degree,
                     const Epetra_Comm     &comm,
                     const bool             extendedghosting)
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
#ifdef ADERLTS
    additional_data.hold_all_faces_to_owned_cells = extendedghosting;
#endif //ADERLTS

    ConstraintMatrix dummy;
    dummy.close();
    MatrixFree<dim,Number> data;
    data.reinit (dof_handler, dummy, quadrature, additional_data);

    return data;
  }
}

// TODO: also need to have a Mapping for representing curved boundaries
template<int dim, int fe_degree, typename Number>
WaveEquationOperation<dim,fe_degree,Number>::
WaveEquationOperation(const std::vector<const DoFHandler<dim> *> &dof_handlers,
                      Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                      Teuchos::RCP<Function<dim> > boundary_conditions,
                      Teuchos::RCP<Function<dim> > source_term,
                      value_type time_step_in,
                      int sourceno,
                      Teuchos::RCP<PATMonitorManager> monitormanagerin)
  :
  time(0.),
  computing_times(7),
  dirichlet_boundary_conditions(boundary_conditions),
  source_term(source_term),
  monitormanager(monitormanagerin)
//  this->data()
{
  this->time_step = time_step_in;

  //{ init data
  ConstraintMatrix dummy;
  dummy.close();
  std::vector<const ConstraintMatrix *> constraints(dof_handlers.size(),&dummy);

  // Add a second quadrature formula that is used for computing the
  // integrals in post-processing, including the cross terms to the standard
  // DoFHandler.
  std::vector<Quadrature<1> > quadratures(2);
  quadratures[0] = QGauss<1>(fe_degree+1);
  quadratures[1] = QGauss<1>(fe_degree+2);


  typename MatrixFree<dim,value_type>::AdditionalData additional_data;
  additional_data.mpi_communicator = MPI_COMM_WORLD;
  additional_data.tasks_parallel_scheme = MatrixFree<dim,value_type>::AdditionalData::partition_partition;
  additional_data.build_face_info = true;
  additional_data.hold_all_faces_to_owned_cells = true;
  additional_data.mapping_update_flags = (update_gradients | update_JxW_values |
                                          update_quadrature_points | update_normal_vectors |
                                          update_values);

  this->data.reinit(dof_handlers,constraints,quadratures,additional_data);
  //}
  mass_matrix_data.reset(new internal::InverseMassMatrixData<dim,fe_degree,dim+1,value_type>(this->data));
  mass_matrix_data_solid.reset(new internal::InverseMassMatrixData<dim,fe_degree,dim*dim+dim+1,value_type>(this->data));
  mass_matrix_data_pml.reset(new internal::InverseMassMatrixData<dim,fe_degree,2*dim+1,value_type>(this->data));

  source_term_no = sourceno;

  densities.resize(this->data.n_macro_cells()+this->data.n_macro_ghost_cells());
  speeds.resize(this->data.n_macro_cells()+this->data.n_macro_ghost_cells());
  dofpermutations.resize(discret->NumMyColElements());

  // only in case of the adjoint run in inverse analysis
  if(monitormanager != Teuchos::null)
  {
    // in case of acouopt, we need the following quantities:
    densities_grad.resize(this->data.n_macro_cells()+this->data.n_macro_ghost_cells());
    speeds_grad.resize(this->data.n_macro_cells()+this->data.n_macro_ghost_cells());
    for (unsigned int i=0; i<this->data.n_macro_cells()+this->data.n_macro_ghost_cells(); ++i)
    {
      densities_grad[i] = make_vectorized_array<value_type>(0.);
      speeds_grad[i] = make_vectorized_array<value_type>(0.);
    }
  }

  // create everything we need to store the permutations
  {
    // size is 5 for the first dimension since we have 4 cases for quads and 5 for hexs
    const unsigned int dpc = this->data.get_dof_handler().get_fe().dofs_per_cell;
    TableIndices<2> table_indices_permute(9,dpc);
    permutevalues.reinit(table_indices_permute);

    unsigned int ndofs1d;
    if(dim==2)
      ndofs1d = std::sqrt(dpc);
    else if(dim==3)
      ndofs1d = int(std::round(std::pow(dpc,1.0/3.0)));
    unsigned int ndofs2d = ndofs1d * ndofs1d;

    if(dim==2)
    {
      // case 1: no permutation
      for(unsigned int j=0; j<dpc; ++j)
      {
        const int ax = j%ndofs1d;
        const int ay = j/ndofs1d;
        permutevalues(0,j) = j;
        permutevalues(1,j) = (ndofs1d-1-ax)*ndofs1d + ay;
        permutevalues(2,j) = (ndofs1d-1-ax) + (ndofs1d-1-ay) * ndofs1d;
        permutevalues(3,j) = ax*ndofs1d + (ndofs1d-1-ay);
      }
    }
    else if(dim==3)
    {
      for(unsigned int j=0; j<dpc; ++j)
      {
        const int ax = j%ndofs1d;
        const int ay = int(j/ndofs1d)%ndofs1d;
        const int az = j/ndofs2d;
        permutevalues(0,j) = j;
        permutevalues(1,j) = ax + (ndofs1d-1-ay)*ndofs2d + az*ndofs1d;
        permutevalues(2,j) = ax*ndofs1d + (ndofs1d-1-ay)*ndofs2d + (ndofs1d-1-az);
        permutevalues(3,j) = (ndofs1d-1-ax)*ndofs1d + (ndofs1d-1-ay)*ndofs2d + az;
        permutevalues(4,j) = (ndofs1d-1-ax) + (ndofs1d-1-ay)*ndofs2d + (ndofs1d-1-az)*ndofs1d;
        permutevalues(5,j) = ax*ndofs1d + ay*ndofs2d + az;
        permutevalues(6,j) = (ndofs1d-1-ax) + ay*ndofs2d + az;
        permutevalues(7,j) = (ndofs1d-1-ax)*ndofs1d + ay*ndofs2d + (ndofs1d-1-az);
        permutevalues(8,j) = ax + ay*ndofs2d + (ndofs1d-1-az)*ndofs1d;

      }
    }

  }

  //ConditionalOStream pcout(std::cout, discret->Comm().MyPID() == 0);
  //this->data.print_memory_consumption(pcout);
}


template <int dim, int fe_degree, typename Number>
WaveEquationOperation<dim,fe_degree,Number>::~WaveEquationOperation()
{
  mass_matrix_data.reset();
  mass_matrix_data_pml.reset();
  mass_matrix_data_solid.reset();

  /* output of computing time for evaluation and application of the inverse mass matrix
  if (computing_times[2] > 0)
    std::cout << "Computing " << (std::size_t)computing_times[2]
              << " times: evaluate "
              << computing_times[0] << "s, inv mass: " << computing_times[1]
              << "s" << std::endl;
  */
}


template<int dim, int fe_degree, typename Number>
void
WaveEquationOperation<dim,fe_degree,Number>::
read_initial_conditions(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                        std::vector<parallel::distributed::Vector<value_type> > &dst)
                        {

  FEEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> phi(this->data);
  for (unsigned int j=0; j<phi.dofs_per_cell; ++j)
    phi.submit_dof_value(Tensor<1,dim+1,VectorizedArray<value_type> >(), j);

  unsigned int dofs_per_cell = phi.dofs_per_cell; // i assume, that it is the same for all cells

  unsigned int nodes_per_cell = GeometryInfo< dim >::vertices_per_cell;
  std::vector<Point<dim> > baci_vals_loc(nodes_per_cell);
  std::vector<Point<dim> > deal_vals_loc(nodes_per_cell);

  for (unsigned int i=0; i<this->data.n_macro_cells()+this->data.n_macro_ghost_cells(); ++i)
  {
    phi.reinit(i);
    for (unsigned int v=0; v<this->data.n_components_filled(i); ++v)
    {
      const int element_index = this->data.get_cell_iterator(i,v)->index();
      DRT::ELEMENTS::Acou * acouele = dynamic_cast<DRT::ELEMENTS::Acou*>(discret->lColElement(element_index));
      if (acouele == NULL)
        dserror("No acoustic element given!");

      // perform permutation: step 1: get the node coordinates
      for (unsigned int n=0; n<nodes_per_cell; ++n)
      {
        for(int d=0; d<dim; ++d)
        {
          deal_vals_loc[n](d) = this->data.get_cell_iterator(i,v)->vertex(n)(d);
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
          dofpermutations[element_index] = 0;
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
          dofpermutations[element_index] = 1;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            int permute = permutevalues(1,i);
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
          dofpermutations[element_index] = 2;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            int permute = permutevalues(2,i);
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
          dofpermutations[element_index] = 3;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            int permute = permutevalues(3,i);
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
          dofpermutations[element_index] = 0;
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
          dofpermutations[element_index] = 1;
          // negative rotation around x
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            int permute = permutevalues(1,i);
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
          dofpermutations[element_index] = 2;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            int permute = permutevalues(2,i);
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
          dofpermutations[element_index] = 3;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            int permute = permutevalues(3,i);
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
          dofpermutations[element_index] = 4;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            int permute = permutevalues(4,i);
            phi.begin_dof_values()[dim*dofs_per_cell+i][v] = acouele->eleinteriorPressnp_(permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[d*dofs_per_cell+i][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+permute);
          }
        }
        else if(deal_vals_loc[0].distance(baci_vals_loc[0])<1e-10 &&
            deal_vals_loc[1].distance(baci_vals_loc[3])<1e-10 &&
            deal_vals_loc[2].distance(baci_vals_loc[4])<1e-10 &&
            deal_vals_loc[3].distance(baci_vals_loc[7])<1e-10 &&
            deal_vals_loc[4].distance(baci_vals_loc[1])<1e-10 &&
            deal_vals_loc[5].distance(baci_vals_loc[2])<1e-10 &&
            deal_vals_loc[6].distance(baci_vals_loc[5])<1e-10 &&
            deal_vals_loc[7].distance(baci_vals_loc[6])<1e-10)
        {
          dofpermutations[element_index] = 5;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            int permute = permutevalues(5,i);
            phi.begin_dof_values()[dim*dofs_per_cell+i][v] = acouele->eleinteriorPressnp_(permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[d*dofs_per_cell+i][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+permute);
          }
        }
       else if(deal_vals_loc[0].distance(baci_vals_loc[1])<1e-10 &&
            deal_vals_loc[1].distance(baci_vals_loc[0])<1e-10 &&
            deal_vals_loc[2].distance(baci_vals_loc[5])<1e-10 &&
            deal_vals_loc[3].distance(baci_vals_loc[4])<1e-10 &&
            deal_vals_loc[4].distance(baci_vals_loc[2])<1e-10 &&
            deal_vals_loc[5].distance(baci_vals_loc[3])<1e-10 &&
            deal_vals_loc[6].distance(baci_vals_loc[6])<1e-10 &&
            deal_vals_loc[7].distance(baci_vals_loc[7])<1e-10)
        {
          dofpermutations[element_index] = 6;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            int permute = permutevalues(6,i);
            phi.begin_dof_values()[dim*dofs_per_cell+i][v] = acouele->eleinteriorPressnp_(permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[d*dofs_per_cell+i][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+permute);
          }
        }
       else if(deal_vals_loc[0].distance(baci_vals_loc[2])<1e-10 &&
            deal_vals_loc[1].distance(baci_vals_loc[1])<1e-10 &&
            deal_vals_loc[2].distance(baci_vals_loc[6])<1e-10 &&
            deal_vals_loc[3].distance(baci_vals_loc[5])<1e-10 &&
            deal_vals_loc[4].distance(baci_vals_loc[3])<1e-10 &&
            deal_vals_loc[5].distance(baci_vals_loc[0])<1e-10 &&
            deal_vals_loc[6].distance(baci_vals_loc[7])<1e-10 &&
            deal_vals_loc[7].distance(baci_vals_loc[4])<1e-10)
        {
          dofpermutations[element_index] = 7;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            int permute = permutevalues(7,i);
            phi.begin_dof_values()[dim*dofs_per_cell+i][v] = acouele->eleinteriorPressnp_(permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[d*dofs_per_cell+i][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+permute);
          }
        }
        else if(deal_vals_loc[0].distance(baci_vals_loc[3])<1e-10 &&
            deal_vals_loc[1].distance(baci_vals_loc[2])<1e-10 &&
            deal_vals_loc[2].distance(baci_vals_loc[7])<1e-10 &&
            deal_vals_loc[3].distance(baci_vals_loc[6])<1e-10 &&
            deal_vals_loc[4].distance(baci_vals_loc[0])<1e-10 &&
            deal_vals_loc[5].distance(baci_vals_loc[1])<1e-10 &&
            deal_vals_loc[6].distance(baci_vals_loc[4])<1e-10 &&
            deal_vals_loc[7].distance(baci_vals_loc[5])<1e-10)
        {
          dofpermutations[element_index] = 8;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            int permute = permutevalues(8,i);
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

template<int dim, int fe_degree, typename Number>
void
WaveEquationOperation<dim,fe_degree,Number>::write_deal_cell_values(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
    const std::vector<parallel::distributed::Vector<value_type> >   &src) const
{

  const unsigned dofs_per_cell = this->data.get_dof_handler().get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> indices, local_dof_indices (dofs_per_cell);
  for (typename DoFHandler<dim>::active_cell_iterator cell = this->data.get_dof_handler().begin_active();
       cell != this->data.get_dof_handler().end(); ++cell)
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
    ghosted_vector[i].reinit(this->data.get_dof_handler().locally_owned_dofs(),relevant_dofs, src[0].get_mpi_communicator());
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
  for (int i=0; i<discret->NumMyColElements(); ++i)
  {
    typename DoFHandler<dim>::active_cell_iterator cell(&this->data.get_dof_handler().get_triangulation(), 0, i, &this->data.get_dof_handler());
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
    else if(deal_vals_loc[0].distance(baci_vals_loc[0])<1e-10 &&
            deal_vals_loc[1].distance(baci_vals_loc[3])<1e-10 &&
            deal_vals_loc[2].distance(baci_vals_loc[4])<1e-10 &&
            deal_vals_loc[3].distance(baci_vals_loc[7])<1e-10 &&
            deal_vals_loc[4].distance(baci_vals_loc[1])<1e-10 &&
            deal_vals_loc[5].distance(baci_vals_loc[2])<1e-10 &&
            deal_vals_loc[6].distance(baci_vals_loc[5])<1e-10 &&
            deal_vals_loc[7].distance(baci_vals_loc[6])<1e-10)
        {
      for (unsigned int d=0; d<dim; ++d)
      {
        cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
        for (unsigned int j=0; j<dofs_per_cell;  ++j)
        {
          const int ax = j%ndofs1d;
          const int ay = int(j/ndofs1d)%ndofs1d;
          const int az = j/ndofs2d;
          int permute = ax*ndofs1d + ay*ndofs2d + az;

          acouele->eleinteriorVelnp_(d*dofs_per_cell+permute) = local_values[j];
        }
      }
      cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
      for (unsigned int j=0; j<dofs_per_cell; ++j)
      {
        const int ax = j%ndofs1d;
        const int ay = int(j/ndofs1d)%ndofs1d;
        const int az = j/ndofs2d;
        int permute =ax*ndofs1d + ay*ndofs2d + az;

        acouele->eleinteriorPressnp_(permute) = local_values[j];
      }
        }
    else if(deal_vals_loc[0].distance(baci_vals_loc[1])<1e-10 &&
            deal_vals_loc[1].distance(baci_vals_loc[0])<1e-10 &&
            deal_vals_loc[2].distance(baci_vals_loc[5])<1e-10 &&
            deal_vals_loc[3].distance(baci_vals_loc[4])<1e-10 &&
            deal_vals_loc[4].distance(baci_vals_loc[2])<1e-10 &&
            deal_vals_loc[5].distance(baci_vals_loc[3])<1e-10 &&
            deal_vals_loc[6].distance(baci_vals_loc[6])<1e-10 &&
            deal_vals_loc[7].distance(baci_vals_loc[7])<1e-10)
        {
      for (unsigned int d=0; d<dim; ++d)
      {
        cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
        for (unsigned int j=0; j<dofs_per_cell; ++j)
        {
        const int ax = j%ndofs1d;
        const int ay = int(j/ndofs1d)%ndofs1d;
        const int az = j/ndofs2d;
        int permute = (ndofs1d-1-ax) + ay*ndofs2d + az;

        acouele->eleinteriorVelnp_(d*dofs_per_cell+permute) = local_values[j];
        }
      }
      cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
      for (unsigned int j=0; j<dofs_per_cell; ++j)
      {
        const int ax = j%ndofs1d;
        const int ay = int(j/ndofs1d)%ndofs1d;
        const int az = j/ndofs2d;
        int permute =(ndofs1d-1-ax) + ay*ndofs2d + az;

        acouele->eleinteriorPressnp_(permute) = local_values[j];
      }

        }
    else if(deal_vals_loc[0].distance(baci_vals_loc[2])<1e-10 &&
            deal_vals_loc[1].distance(baci_vals_loc[1])<1e-10 &&
            deal_vals_loc[2].distance(baci_vals_loc[6])<1e-10 &&
            deal_vals_loc[3].distance(baci_vals_loc[5])<1e-10 &&
            deal_vals_loc[4].distance(baci_vals_loc[3])<1e-10 &&
            deal_vals_loc[5].distance(baci_vals_loc[0])<1e-10 &&
            deal_vals_loc[6].distance(baci_vals_loc[7])<1e-10 &&
            deal_vals_loc[7].distance(baci_vals_loc[4])<1e-10)
        {
      for (unsigned int d=0; d<dim; ++d)
      {
        cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
        for (unsigned int j=0; j<dofs_per_cell; ++j)
        {
        const int ax = j%ndofs1d;
        const int ay = int(j/ndofs1d)%ndofs1d;
        const int az = j/ndofs2d;
        int permute = (ndofs1d-1-ax)*ndofs1d + ay*ndofs2d + (ndofs1d-1-az);

        acouele->eleinteriorVelnp_(d*dofs_per_cell+permute) = local_values[j];
        }
      }
      cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
      for (unsigned int j=0; j<dofs_per_cell; ++j)
      {
        const int ax = j%ndofs1d;
        const int ay = int(j/ndofs1d)%ndofs1d;
        const int az = j/ndofs2d;
        int permute =(ndofs1d-1-ax)*ndofs1d + ay*ndofs2d + (ndofs1d-1-az);

        acouele->eleinteriorPressnp_(permute) = local_values[j];
      }

        }
      else if(deal_vals_loc[0].distance(baci_vals_loc[3])<1e-10 &&
            deal_vals_loc[1].distance(baci_vals_loc[2])<1e-10 &&
            deal_vals_loc[2].distance(baci_vals_loc[7])<1e-10 &&
            deal_vals_loc[3].distance(baci_vals_loc[6])<1e-10 &&
            deal_vals_loc[4].distance(baci_vals_loc[0])<1e-10 &&
            deal_vals_loc[5].distance(baci_vals_loc[1])<1e-10 &&
            deal_vals_loc[6].distance(baci_vals_loc[4])<1e-10 &&
            deal_vals_loc[7].distance(baci_vals_loc[5])<1e-10)
        {
      for (unsigned int d=0; d<dim; ++d)
      {
        cell->get_interpolated_dof_values(ghosted_vector[d], local_values);
        for (unsigned int j=0; j<dofs_per_cell; ++j)
        {
        const int ax = j%ndofs1d;
        const int ay = int(j/ndofs1d)%ndofs1d;
        const int az = j/ndofs2d;
        int permute = ax + ay*ndofs2d + (ndofs1d-1-az)*ndofs1d;

        acouele->eleinteriorVelnp_(d*dofs_per_cell+permute) = local_values[j];
        }
      }
      cell->get_interpolated_dof_values(ghosted_vector[dim], local_values);
      for (unsigned int j=0; j<dofs_per_cell; ++j)
      {
        const int ax = j%ndofs1d;
        const int ay = int(j/ndofs1d)%ndofs1d;
        const int az = j/ndofs2d;
        int permute = ax + ay*ndofs2d + (ndofs1d-1-az)*ndofs1d;

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


template<int dim, int fe_degree, typename Number>
void
WaveEquationOperationAcousticWavePML<dim,fe_degree,Number>::write_deal_cell_values(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
    const std::vector<parallel::distributed::Vector<value_type> >   &src) const
{

  const unsigned dofs_per_cell = this->data.get_dof_handler().get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> indices, local_dof_indices (dofs_per_cell);
  for (typename DoFHandler<dim>::active_cell_iterator cell = this->data.get_dof_handler().begin_active();
       cell != this->data.get_dof_handler().end(); ++cell)
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
    ghosted_vector[i].reinit(this->data.get_dof_handler().locally_owned_dofs(),relevant_dofs, src[0].get_mpi_communicator());
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
  for (int i=0; i<discret->NumMyColElements(); ++i)
  {
    typename DoFHandler<dim>::active_cell_iterator cell(&this->data.get_dof_handler().get_triangulation(), 0, i, &this->data.get_dof_handler());
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
        for (unsigned int d=0; d<dim; ++d)
        {
          cell->get_interpolated_dof_values(ghosted_vector[dim+1+d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            acouele->eleinteriorAuxiliaryPML_(d*dofs_per_cell+j) = local_values[j];
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
        {
          cell->get_interpolated_dof_values(ghosted_vector[dim+1+d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = j/ndofs1d;
            int permute = (ndofs1d-1-ax)*ndofs1d + ay;
            acouele->eleinteriorAuxiliaryPML_(d*dofs_per_cell+permute) = local_values[j];
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
        {
          cell->get_interpolated_dof_values(ghosted_vector[dim+1+d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = j/ndofs1d;
            int permute = (ndofs1d-1-ax) + (ndofs1d-1-ay) * ndofs1d;
            acouele->eleinteriorAuxiliaryPML_(d*dofs_per_cell+permute) = local_values[j];
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
        {
          cell->get_interpolated_dof_values(ghosted_vector[dim+1+d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = j/ndofs1d;
            int permute = (ax) * ndofs1d + (ndofs1d-1-ay);
            acouele->eleinteriorAuxiliaryPML_(d*dofs_per_cell+permute) = local_values[j];
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
        for (unsigned int d=0; d<dim; ++d)
        {
          cell->get_interpolated_dof_values(ghosted_vector[dim+1+d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            acouele->eleinteriorAuxiliaryPML_(d*dofs_per_cell+j) = local_values[j];
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
        for (unsigned int d=0; d<dim; ++d)
        {
          cell->get_interpolated_dof_values(ghosted_vector[dim+1+d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = int(j/ndofs1d)%ndofs1d;
            const int az = j/ndofs2d;
            int permute = ax + (ndofs1d-1-ay)*ndofs2d + az*ndofs1d;
            acouele->eleinteriorAuxiliaryPML_(d*dofs_per_cell+permute) = local_values[j];
          }
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
        for (unsigned int d=0; d<dim; ++d)
        {
          cell->get_interpolated_dof_values(ghosted_vector[dim+1+d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = int(j/ndofs1d)%ndofs1d;
            const int az = j/ndofs2d;
            int permute = ax*ndofs1d + (ndofs1d-1-ay)*ndofs2d + (ndofs1d-1-az);
            acouele->eleinteriorAuxiliaryPML_(d*dofs_per_cell+permute) = local_values[j];
          }
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
        for (unsigned int d=0; d<dim; ++d)
        {
          cell->get_interpolated_dof_values(ghosted_vector[dim+1+d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = int(j/ndofs1d)%ndofs1d;
            const int az = j/ndofs2d;
            int permute = (ndofs1d-1-ax)*ndofs1d + (ndofs1d-1-ay)*ndofs2d + az;
            acouele->eleinteriorAuxiliaryPML_(d*dofs_per_cell+permute) = local_values[j];
          }
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
        for (unsigned int d=0; d<dim; ++d)
        {
          cell->get_interpolated_dof_values(ghosted_vector[dim+1+d], local_values);
          for (unsigned int j=0; j<dofs_per_cell; ++j)
          {
            const int ax = j%ndofs1d;
            const int ay = int(j/ndofs1d)%ndofs1d;
            const int az = j/ndofs2d;
            int permute = (ndofs1d-1-ax) + (ndofs1d-1-ay)*ndofs2d + (ndofs1d-1-az)*ndofs1d;
            acouele->eleinteriorAuxiliaryPML_(d*dofs_per_cell+permute) = local_values[j];
          }
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

template<int dim, int fe_degree, typename Number>
void
WaveEquationOperation<dim,fe_degree,Number>::
compute_gradient_contributions(std::vector<parallel::distributed::Vector<value_type> > &fwnp,
                               std::vector<parallel::distributed::Vector<value_type> > &fwn,
                               std::vector<parallel::distributed::Vector<value_type> > &adnp)
{
  // we need the derivative of the mass matrix with respect to density and sound speed and have to buuild the scalar product with
  // corresponding pressure and velocity values
  FEEvaluation<dim,fe_degree,fe_degree+1,dim,value_type> adjoint_velocity(this->data);
  FEEvaluation<dim,fe_degree,fe_degree+1,1,value_type>   adjoint_pressure(this->data);

  FEEvaluation<dim,fe_degree,fe_degree+1,dim,value_type> forward_velocity(this->data);
  FEEvaluation<dim,fe_degree,fe_degree+1,1,value_type>   forward_pressure(this->data);

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

  for (unsigned int cell=0; cell<this->data.n_macro_cells()+this->data.n_macro_ghost_cells(); ++cell)
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
        densities_grad[cell] += 1./(densities[cell]) * adjoint_velocity.get_dof_value(dof)[d] * forward_velocity.get_dof_value(dof)[d]; // factor is 1
      pressure_mass_mat_contrib += adjoint_pressure.get_dof_value(dof) * forward_pressure.get_dof_value(dof);
    }
    speeds_grad[cell] -= c_fac * pressure_mass_mat_contrib;
    densities_grad[cell] -= rho_fac * pressure_mass_mat_contrib;
  }

  return;
}

template<int dim, int fe_degree, typename Number>
double
WaveEquationOperation<dim,fe_degree,Number>::get_SoS_gradient(int colid) const
{
  for (unsigned int i=0; i<this->data.n_macro_cells()+this->data.n_macro_ghost_cells(); ++i)
    for (unsigned int v=0; v<this->data.n_components_filled(i); ++v)
      {
        if(this->data.get_cell_iterator(i,v)->index()==colid)
          return speeds_grad[i][v]/this->time_step;
      }

  return 0.0;
}

template<int dim, int fe_degree, typename Number>
double
WaveEquationOperation<dim,fe_degree,Number>::get_density_gradient(int colid) const
{
  for (unsigned int i=0; i<this->data.n_macro_cells()+this->data.n_macro_ghost_cells(); ++i)
    for (unsigned int v=0; v<this->data.n_components_filled(i); ++v)
      if(this->data.get_cell_iterator(i,v)->index()==colid)
        return densities_grad[i][v]/this->time_step;
  return 0.0;
}

template <int dim, int fe_degree, typename Number>
typename WaveEquationOperationBase<dim,Number>::value_type WaveEquationOperation<dim,fe_degree,Number>::evaluate_source_adjoint(const Point<dim> &p, const std::vector<std::vector<value_type> > nodes, std::vector<value_type> values) const
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

//    if(node_distance<quad_distance)
//    {
//      std::cout<<"punkt "<<xyz[0]<<" "<<xyz[1]<<" nodecoords "<<nodes[0][0]<<" "<<nodes[0][1]
//                                              <<" "           <<nodes[1][0]<<" "<<nodes[1][1]<<" "
//                                              <<" nodevals  " <<values[0]  <<" "<<values[1]<<std::endl;
//      std::cout<<"wrong face!!!"<<std::endl<<std::endl;
//    }
  }
  else if(dim==3 && nodes.size()==4) // hex8 with quad4 face element
  {
    value_type distance1 = std::sqrt( (xyz[0]-nodes[0][0])*(xyz[0]-nodes[0][0])
                                    + (xyz[1]-nodes[0][1])*(xyz[1]-nodes[0][1])
                                    + (xyz[2]-nodes[0][2])*(xyz[2]-nodes[0][2]) );
    value_type distance2 = std::sqrt( (xyz[0]-nodes[1][0])*(xyz[0]-nodes[1][0])
                                    + (xyz[1]-nodes[1][1])*(xyz[1]-nodes[1][1])
                                    + (xyz[2]-nodes[1][2])*(xyz[2]-nodes[1][2]) );
    value_type distance3 = std::sqrt( (xyz[0]-nodes[2][0])*(xyz[0]-nodes[2][0])
                                    + (xyz[1]-nodes[2][1])*(xyz[1]-nodes[2][1])
                                    + (xyz[2]-nodes[2][2])*(xyz[2]-nodes[2][2]) );
    value_type distance4 = std::sqrt( (xyz[0]-nodes[3][0])*(xyz[0]-nodes[3][0])
                                    + (xyz[1]-nodes[3][1])*(xyz[1]-nodes[3][1])
                                    + (xyz[2]-nodes[3][2])*(xyz[2]-nodes[3][2]) );
    value_type entiredistance = distance1+distance2+distance3+distance4;

    result = (entiredistance-distance1) * values[0]
           + (entiredistance-distance2) * values[1]
           + (entiredistance-distance3) * values[2]
           + (entiredistance-distance4) * values[3];
    result *= -1.0;
    result /= (3.0*entiredistance);
  }
  else
    dserror("not yet implemented");

  return result;
}

template <int dim, int fe_degree, typename Number>
typename WaveEquationOperationBase<dim,Number>::value_type WaveEquationOperation<dim,fe_degree,Number>::evaluate_source_timereversal(const Point<dim> &p, const std::vector<std::vector<value_type> > nodes, std::vector<value_type> values) const
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

template<int dim, int fe_degree, typename Number>
WaveEquationOperationAcousticWave<dim,fe_degree,Number>::
WaveEquationOperationAcousticWave(const std::vector<const DoFHandler<dim> *> &dof_handlers,
                                  Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                                  Teuchos::RCP<Function<dim> > boundary_conditions,
                                  Teuchos::RCP<Function<dim> > source_term,
                                  value_type time_step_in,
                                  int sourceno,
                                  Teuchos::RCP<PATMonitorManager> monitormanagerin)
  :
  WaveEquationOperation<dim,fe_degree,Number>(dof_handlers,discret,boundary_conditions,source_term,time_step_in,sourceno,monitormanagerin)
{
  for (unsigned int i=0; i<this->data.n_macro_cells()+this->data.n_macro_ghost_cells(); ++i)
  {
    this->densities[i] = make_vectorized_array<value_type>(1.);
    this->speeds[i] = make_vectorized_array<value_type>(1.);
    for (unsigned int v=0; v<this->data.n_components_filled(i); ++v)
    {
      if (this->data.get_cell_iterator(i,v)->level() != 0)
        dserror("Refined meshes currently not implemented!");

      const int element_index = this->data.get_cell_iterator(i,v)->index();

      this->densities[i][v] = discret->lColElement(element_index)->Material()->Parameter()->GetParameter(0,discret->lColElement(element_index)->Id());
      this->speeds[i][v] = discret->lColElement(element_index)->Material()->Parameter()->GetParameter(1,discret->lColElement(element_index)->Id());
    }
  }
}

template<int dim, int fe_degree, typename Number>
void WaveEquationOperation<dim,fe_degree,Number>::
local_apply_domain(const MatrixFree<dim,value_type>                                &data,
                   std::vector<parallel::distributed::Vector<value_type> >         &dst,
                   const std::vector<parallel::distributed::Vector<value_type> >   &src,
                   const std::pair<unsigned int,unsigned int>                 &cell_range) const
{
  FEEvaluation<dim,fe_degree,fe_degree+1,dim,value_type> velocity(data);
  FEEvaluation<dim,fe_degree,fe_degree+1,1,value_type> pressure(data);

  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
    {
      evaluate_cell(velocity,pressure,src,cell);
      velocity.distribute_local_to_global (dst, 0);
      pressure.distribute_local_to_global (dst,dim);
    }
}

template<int dim, int fe_degree, typename Number>
void WaveEquationOperation<dim,fe_degree,Number>::
evaluate_cell(FEEvaluation<dim,fe_degree,fe_degree+1,dim,value_type>        &phi_v,
              FEEvaluation<dim,fe_degree,fe_degree+1,1,value_type>          &phi_p,
              const std::vector<parallel::distributed::Vector<value_type> > &src,
              const unsigned int                                             cell) const
{

    // It is faster to evaluate values of the vector-valued velocity and
    // gradients of the scalar pressure than divergence of velocity and
    // values of pressure
    phi_v.reinit(cell);
    phi_v.read_dof_values(src, 0);
    phi_v.evaluate (true, false, false);

    phi_p.reinit(cell);
    phi_p.read_dof_values(src, dim);
    phi_p.evaluate(false, true, false);

    const VectorizedArray<value_type> rho = this->densities[cell];
    const VectorizedArray<value_type> rho_inv = 1./this->densities[cell];
    const VectorizedArray<value_type> c_sq = this->speeds[cell]*this->speeds[cell];

    for (unsigned int q=0; q<phi_v.n_q_points; ++q)
    {
      const Tensor<1,dim,VectorizedArray<value_type> >
      pressure_gradient = phi_p.get_gradient(q);
      const Tensor<1,dim,VectorizedArray<value_type> >
      velocity_value = phi_v.get_value(q);

      Point<dim,VectorizedArray<value_type> > q_points = phi_v.quadrature_point(q);
      VectorizedArray<value_type> rhs =  make_vectorized_array<value_type>(0.0);
      for (unsigned int n=0; n<rhs.n_array_elements; ++n)
      {
        Point<dim> q_point;
        for (unsigned int d=0; d<dim; ++d)
          q_point[d] = q_points[d][n];
        rhs[n] = this->source_term->value(q_point);
      }

      phi_p.submit_value(c_sq*rhs,q);
      if(this->adjoint_eval==false)
      {
        phi_v.submit_value(-rho_inv*pressure_gradient,q);
        phi_p.submit_gradient(rho*c_sq*velocity_value,q);
      }
      else
      {
        phi_v.submit_value(rho*c_sq*pressure_gradient,q);
        phi_p.submit_gradient(-rho_inv*velocity_value,q);
      }
    }

    phi_v.integrate (true, false);
    phi_p.integrate(true, true);
}


template<int dim, int fe_degree, typename Number>
void WaveEquationOperationAcousticWavePML<dim,fe_degree,Number>::
local_apply_domain(const MatrixFree<dim,value_type>                                &data,
                   std::vector<parallel::distributed::Vector<value_type> >         &dst,
                   const std::vector<parallel::distributed::Vector<value_type> >   &src,
                   const std::pair<unsigned int,unsigned int>                 &cell_range) const
{
  FEEvaluation<dim,fe_degree,fe_degree+1,dim,value_type> velocity(data);
  FEEvaluation<dim,fe_degree,fe_degree+1,1,value_type> pressure(data);
  FEEvaluation<dim,fe_degree,fe_degree+1,dim,value_type> auxiliary(data);

  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
  {
    // It is faster to evaluate values of the vector-valued velocity and
    // gradients of the scalar pressure than divergence of velocity and
    // values of pressure
    velocity.reinit(cell);
    velocity.read_dof_values(src, 0);
    velocity.evaluate (true, true, false);

    pressure.reinit(cell);
    pressure.read_dof_values(src, dim);
    pressure.evaluate(false, true, false);

    auxiliary.reinit(cell);
    auxiliary.read_dof_values(src, dim+1);
    auxiliary.evaluate (true, false, false);

    const VectorizedArray<value_type> rho = this->densities[cell];
    const VectorizedArray<value_type> rho_inv = 1./this->densities[cell];
    const VectorizedArray<value_type> c_sq = this->speeds[cell]*this->speeds[cell];

    Tensor<1,dim,VectorizedArray<value_type> > sigma_values;
    Tensor<1,dim,VectorizedArray<value_type> > eigen_values;

    Tensor<2,dim,VectorizedArray<value_type> > Matrix_A;
    Tensor<3,dim,VectorizedArray<value_type> > eigen_tensors;


    for (unsigned int q=0; q<velocity.n_q_points; ++q)
    {
      const Tensor<1,dim,VectorizedArray<value_type> > pressure_gradient = pressure.get_gradient(q);
      const Tensor<1,dim,VectorizedArray<value_type> > velocity_value    = velocity.get_value(q);
      const Tensor<1,dim,VectorizedArray<value_type> > auxiliary_value   = auxiliary.get_value(q);
      const Tensor<2,dim,VectorizedArray<value_type> > velocity_gradient = velocity.get_gradient(q);

      Point<dim,VectorizedArray<value_type> > q_points = velocity.quadrature_point(q);
      VectorizedArray<value_type> rhs =  make_vectorized_array<value_type>(0.0);
      for (unsigned int n=0; n<rhs.n_array_elements; ++n)
      {
        Point<dim> q_point;
        for (unsigned int d=0; d<dim; ++d)
          q_point[d] = q_points[d][n];
        rhs[n] = this->source_term->value(q_point);
        if (layer_reference[cell][n].size())
          sigma_pml->get_matrix (layer_reference[cell][n], n, q_point, sigma_values,eigen_values, Matrix_A, eigen_tensors);
      }


      // calculate the quota form auxiliary field
      VectorizedArray<value_type>                aux_quota_pressure;
      Tensor<1,dim,VectorizedArray<value_type> > aux_quota_velocity;
      Tensor<1,dim,VectorizedArray<value_type> > aux_quota_auxiliary;

      for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = 0; j < dim; ++j)
          aux_quota_velocity[i] += Matrix_A[i][j] * velocity_value[j];


      aux_quota_pressure  = 0;
      for (unsigned int n = 0; n < dim; ++n)
      {
        aux_quota_pressure += sigma_values[n] * auxiliary_value[n];
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int j = 0; j < dim; ++j)
            aux_quota_auxiliary[n] += eigen_tensors[n][i][j] * velocity_gradient[i][j];
        aux_quota_auxiliary[n] += eigen_values[n] * auxiliary_value[n];
      }

      if(this->adjoint_eval==false)
      {
        velocity.submit_value(-rho_inv*pressure_gradient-aux_quota_velocity,q);
        pressure.submit_value(c_sq*rhs-rho*c_sq*aux_quota_pressure,q);
        pressure.submit_gradient(rho*c_sq*velocity_value,q);
        auxiliary.submit_value(-aux_quota_auxiliary, q);
      }
      else
      {
        velocity.submit_value(rho*c_sq*pressure_gradient-aux_quota_velocity,q);
        pressure.submit_value(c_sq*rhs+rho_inv*aux_quota_pressure,q);
        pressure.submit_gradient(-rho_inv*velocity_value,q);
        auxiliary.submit_value(-aux_quota_auxiliary, q);
      }
    }

    velocity.integrate (true, false);
    velocity.distribute_local_to_global (dst, 0);

    pressure.integrate(true, true);
    pressure.distribute_local_to_global (dst,dim);

    auxiliary.integrate(true, false);
    auxiliary.distribute_local_to_global(dst, dim+1);

  }
}



template <int dim, int fe_degree, typename Number>
void
WaveEquationOperation<dim,fe_degree,Number>::
local_apply_face (const MatrixFree<dim,value_type> &,
  std::vector<parallel::distributed::Vector<value_type> >        &dst,
  const std::vector<parallel::distributed::Vector<value_type> >  &src,
  const std::pair<unsigned int,unsigned int>                     &face_range) const
{
  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> phi(this->data, true, 0, 0, 0, true);
  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> phi_neighbor(this->data, false, 0, 0, 0, true);

  for (unsigned int face=face_range.first; face<face_range.second; face++)
    {
      evaluate_inner_face(phi,phi_neighbor,src,face,1.0);

      phi.distribute_local_to_global(dst, 0);
      phi_neighbor.distribute_local_to_global(dst, 0);
    }
}

template <int dim, int fe_degree, typename Number>
void
WaveEquationOperation<dim,fe_degree,Number>::
evaluate_inner_face(FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type>  &phi,
                    FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type>  &phi_neighbor,
                    const std::vector<parallel::distributed::Vector<value_type> > &src,
                    unsigned int                                                   face,
                    value_type                                                     boundary_fac) const
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

  AssertDimension(phi.n_q_points, this->data.get_n_q_points_face(0));

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
      VectorizedArray<value_type> lambda;
      VectorizedArray<value_type> pres_diff_plus;
      VectorizedArray<value_type> pres_diff_minus;
      if(this->adjoint_eval==false)
      {
        lambda = tau_inv*( normal_v_plus + normal_v_minus
                           + tau_plus * val_plus[dim]
                           + tau_minus * val_minus[dim]);
        pres_diff_plus = (val_plus[dim]-lambda)*rho_inv_plus;
        pres_diff_minus = (val_minus[dim]-lambda)*rho_inv_minus;
      }
      else
      {
        lambda = tau_inv*(rho_inv_plus*normal_v_plus + rho_inv_minus*normal_v_minus - tau_plus*rho_plus*c_sq_plus*val_plus[dim] - tau_minus*rho_minus*c_sq_minus*val_minus[dim]);
        pres_diff_plus  = -rho_plus*c_sq_plus*val_plus[dim] - lambda ;
        pres_diff_minus = -rho_minus*c_sq_minus*val_minus[dim] - lambda;
      }
      for (unsigned int d=0; d<dim; ++d)
      {
        val_plus[d] = boundary_fac*pres_diff_plus*normal[d];
        val_minus[d] = -boundary_fac*pres_diff_minus*normal[d];
      }

      if(this->adjoint_eval==false)
      {
        val_plus[dim] = boundary_fac * c_sq_plus * rho_plus * (-normal_v_plus + tau_plus * (lambda - val_plus[dim]));
        val_minus[dim] = boundary_fac * c_sq_minus * rho_minus * (-normal_v_minus + tau_minus * (lambda - val_minus[dim]));
      }
      else
      {
        val_plus[dim] = -boundary_fac *(-rho_inv_plus*normal_v_plus + tau_plus * (c_sq_plus*rho_plus*val_plus[dim] + lambda));
        val_minus[dim] = -boundary_fac *(-rho_inv_minus*normal_v_minus + tau_minus * (c_sq_minus*rho_minus*val_minus[dim] + lambda));
      }

      phi.submit_value(val_plus, q);
      phi_neighbor.submit_value(val_minus, q);
    }
  phi.integrate(true,false);
  phi_neighbor.integrate(true,false);
}



template <int dim, int fe_degree, typename Number>
void
WaveEquationOperationAcousticWavePML<dim,fe_degree,Number>::
local_apply_face (const MatrixFree<dim,value_type> &,
  std::vector<parallel::distributed::Vector<value_type> >        &dst,
  const std::vector<parallel::distributed::Vector<value_type> >  &src,
  const std::pair<unsigned int,unsigned int>                     &face_range) const
{
  // There is some overhead in the methods in FEEvaluation, so it is faster
  // to combine pressure and velocity in the same object and just combine
  // them at the level of quadrature points
  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> phi(this->data, true, 0, 0, 0, true);
  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> phi_neighbor(this->data, false, 0, 0, 0 , true);

  for (unsigned int face=face_range.first; face<face_range.second; face++)
  {
    phi.reinit(face);
    phi.read_dof_values(src, 0);
    phi.evaluate(true,false);
    const VectorizedArray<value_type> rho_plus = phi.read_cell_data(this->densities);
    const VectorizedArray<value_type> rho_inv_plus = 1./rho_plus;
    const VectorizedArray<value_type> c_plus = phi.read_cell_data(this->speeds);
    const VectorizedArray<value_type> c_sq_plus = c_plus * c_plus;
    const VectorizedArray<value_type> tau_plus = 1./c_plus/rho_plus;

    phi_neighbor.reinit(face);
    phi_neighbor.read_dof_values(src, 0);
    phi_neighbor.evaluate(true,false);
    const VectorizedArray<value_type> rho_minus = phi_neighbor.read_cell_data(this->densities);
    const VectorizedArray<value_type> rho_inv_minus = 1./rho_minus;
    const VectorizedArray<value_type> c_minus = phi_neighbor.read_cell_data(this->speeds);
    const VectorizedArray<value_type> c_sq_minus = c_minus * c_minus;
    const VectorizedArray<value_type> tau_minus = 1./c_minus/rho_minus;

    const VectorizedArray<value_type> tau_inv = 1./(tau_plus + tau_minus);

    AssertDimension(phi.n_q_points, this->data.get_n_q_points_face(0));

    for (unsigned int q=0; q<phi.n_q_points; ++q)
    {
      Point<dim,VectorizedArray<value_type> > q_point = phi.quadrature_point(q);
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
        if(this->adjoint_eval) // second query required for intermediate integration
        {
          lambda = tau_inv*(rho_inv_plus*normal_v_plus + rho_inv_minus*normal_v_minus - tau_plus*rho_plus*c_sq_plus*val_plus[dim] - tau_minus*rho_minus*c_sq_minus*val_minus[dim]);
          if(inner_face_monitored[face].any())
          {
            for (unsigned int v=0; v<VectorizedArray<value_type>::n_array_elements && this->data.faces[face].left_cell[v] != numbers::invalid_unsigned_int; ++v)
            {
              if(inner_face_monitored[face][v])
              {
                std::vector<double> vecpoint(dim);
                double value = 0.0;
                for (unsigned int d=0; d<dim; ++d)
                  vecpoint[d] = q_point[d][v];
                this->monitormanager->EvaluateMonitorValues(this->time,vecpoint,value);
                lambda[v] += tau_inv[v] * value;
              }
            }
          }
        }

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

template <int dim, int fe_degree, typename Number>
void WaveEquationOperationAcousticWave<dim,fe_degree,Number>::
local_apply_domain (const MatrixFree<dim,value_type> &data,
                           std::vector<parallel::distributed::Vector<value_type> >       &dst,
                           const std::vector<parallel::distributed::Vector<value_type> > &src,
                           const std::pair<unsigned int,unsigned int>               &cell_range) const
{
  ACOU::WaveEquationOperation<dim,fe_degree,Number>::local_apply_domain(data,dst,src,cell_range);
}
template <int dim, int fe_degree, typename Number>
void WaveEquationOperationAcousticWave<dim,fe_degree,Number>::
local_apply_face (const MatrixFree<dim,value_type> &data,
                           std::vector<parallel::distributed::Vector<value_type> >       &dst,
                           const std::vector<parallel::distributed::Vector<value_type> > &src,
                           const std::pair<unsigned int,unsigned int>               &face_range) const
{
  ACOU::WaveEquationOperation<dim,fe_degree,Number>::local_apply_face(data,dst,src,face_range);
}
template <int dim, int fe_degree, typename Number>
void WaveEquationOperationAcousticWave<dim,fe_degree,Number>::
local_apply_boundary_face (const MatrixFree<dim,value_type> &data,
                           std::vector<parallel::distributed::Vector<value_type> >       &dst,
                           const std::vector<parallel::distributed::Vector<value_type> > &src,
                           const std::pair<unsigned int,unsigned int>               &face_range) const
{
  ACOU::WaveEquationOperation<dim,fe_degree,Number>::local_apply_boundary_face(data,dst,src,face_range);
}
template <int dim, int fe_degree, typename Number>
void WaveEquationOperationAcousticWave<dim,fe_degree,Number>::
local_apply_mass_matrix (const MatrixFree<dim,value_type> &data,
                           std::vector<parallel::distributed::Vector<value_type> >       &dst,
                           const std::vector<parallel::distributed::Vector<value_type> > &src,
                           const std::pair<unsigned int,unsigned int>               &cell_range) const
{
  ACOU::WaveEquationOperation<dim,fe_degree,Number>::local_apply_mass_matrix(data,dst,src,cell_range);
}


template <int dim, int fe_degree, typename Number>
void WaveEquationOperation<dim,fe_degree,Number>::
local_apply_boundary_face (const MatrixFree<dim,value_type> &,
                           std::vector<parallel::distributed::Vector<value_type> >       &dst,
                           const std::vector<parallel::distributed::Vector<value_type> > &src,
                           const std::pair<unsigned int,unsigned int>               &face_range) const
{
  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> phi(this->data, true, 0, 0, 0, true);
  for (unsigned int face=face_range.first; face<face_range.second; face++)
    {
      evaluate_boundary_face(phi,src,face,1.0);
      phi.distribute_local_to_global(dst, 0);
    }
}

template <int dim, int fe_degree, typename Number>
void WaveEquationOperation<dim,fe_degree,Number>::
evaluate_boundary_face(FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type>  &phi,
                       const std::vector<parallel::distributed::Vector<value_type> > &src,
                       unsigned int                                                   face,
                       value_type                                                     boundary_fac) const
{
  phi.reinit(face);
  phi.read_dof_values(src, 0);
  phi.evaluate(true,false);

  const VectorizedArray<value_type> rho = phi.read_cell_data(this->densities);
  const VectorizedArray<value_type> rho_inv = 1./rho;
  const VectorizedArray<value_type> c = phi.read_cell_data(this->speeds);
  const VectorizedArray<value_type> c_sq = c*c;
  const VectorizedArray<value_type> tau = 1./c*rho_inv;

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
      if(this->adjoint_eval==false)
        lambda = 1./tau*normal_v_plus+p_plus;//VectorizedArray<value_type>();
      else
        lambda = 1./(tau)*rho_inv*normal_v_plus - tau*rho*c_sq/(tau)*p_plus;
      if(this->adjoint_eval == true) // second query required for intermediate integration
      {
        for (unsigned int v=0; v<VectorizedArray<value_type>::n_array_elements && this->data.faces[face].left_cell[v] != numbers::invalid_unsigned_int; ++v)
        {
          std::vector<double> vecpoint(dim);
          double value = 0.0;
          for (unsigned int d=0; d<dim; ++d)
            vecpoint[d] = q_point[d][v];
          this->monitormanager->EvaluateMonitorValues(this->time,vecpoint,value);
          lambda[v] += 1.0/tau[v] * value;
        }
      }
    }
    else if(int_boundary_id==2) // monitored and absorbing
    {
      if(this->adjoint_eval==false)
        lambda = tau/(tau+1./c/rho)*p_plus + 1./(tau+1./c/rho)*normal_v_plus;
      else
      {
        lambda = 1./(tau+1./c/rho)*rho_inv*normal_v_plus - tau*rho*c_sq/(tau+1./c/rho)*p_plus;
        for (unsigned int v=0; v<VectorizedArray<value_type>::n_array_elements && this->data.faces[face].left_cell[v] != numbers::invalid_unsigned_int; ++v)
        {
          std::vector<double> vecpoint(dim);
          double value = 0.0;
          for (unsigned int d=0; d<dim; ++d)
            vecpoint[d] = q_point[d][v];
          this->monitormanager->EvaluateMonitorValues(this->time,vecpoint,value);
          lambda[v] += 1.0/(tau[v]+1./c[v]/rho[v])  * value;
        }
      }
    }
    else if(int_boundary_id==3) // free boundary
    {
      if(this->adjoint_eval==false)
        lambda = 1./tau*normal_v_plus+p_plus; // VectorizedArray<value_type>();
      else
        lambda = 1./(tau)*rho_inv*normal_v_plus - tau*rho*c_sq/(tau)*p_plus;
    }
    else if(int_boundary_id==4) // dbc from time reversal
    {
      for (unsigned int v=0; v<VectorizedArray<value_type>::n_array_elements && this->data.faces[face].left_cell[v] != numbers::invalid_unsigned_int; ++v)
      {
        std::vector<double> vecpoint(dim);
        double value = 0.0;
        for (unsigned int d=0; d<dim; ++d)
          vecpoint[d] = q_point[d][v];
        this->monitormanager->EvaluateMonitorValues(this->time,vecpoint,value);
        lambda[v] = value;
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
          lambda[v] = this->dirichlet_boundary_conditions->value(point,(int_boundary_id-5)*dim);
        }
      else
        lambda = VectorizedArray<value_type>();
    }

    if(this->adjoint_eval==false)
    {
      for (unsigned int d=0; d<dim; ++d)
        val_plus[d] = boundary_fac*(p_plus - lambda)*normal[d]*rho_inv;
      val_plus[dim] = -c_sq*rho*boundary_fac*(normal_v_plus - tau*(lambda - p_plus));
    }
    else
    {
      for (unsigned int d=0; d<dim; ++d)
        val_plus[d] = -boundary_fac*(lambda + rho*c_sq*p_plus)*normal[d];
      val_plus[dim] = -boundary_fac*(-rho_inv*normal_v_plus + tau*(lambda + rho*c_sq*p_plus));
    }
    phi.submit_value(val_plus,q);
  }
  phi.integrate(true,false);
}

template<int dim, int fe_degree, typename Number>
void WaveEquationOperation<dim, fe_degree,Number>::
local_apply_mass_matrix(const MatrixFree<dim,value_type>                  &data,
                        std::vector<parallel::distributed::Vector<value_type> >        &dst,
                        const std::vector<parallel::distributed::Vector<value_type> >  &src,
                        const std::pair<unsigned int,unsigned int>    &cell_range) const
{
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
  {
    mass_matrix_data->phi[0].reinit(cell);
    mass_matrix_data->phi[0].read_dof_values(src, 0);

    mass_matrix_data->inverse.fill_inverse_JxW_values(mass_matrix_data->coefficients);
    mass_matrix_data->inverse.apply(mass_matrix_data->coefficients, dim+1,
                                    mass_matrix_data->phi[0].begin_dof_values(),
                                    mass_matrix_data->phi[0].begin_dof_values());

    mass_matrix_data->phi[0].set_dof_values(dst,0);
  }
}

template<int dim, int fe_degree, typename Number>
void WaveEquationOperationAcousticWavePML<dim, fe_degree,Number>::
local_apply_mass_matrix(const MatrixFree<dim,value_type>                  &data,
                        std::vector<parallel::distributed::Vector<value_type> >        &dst,
                        const std::vector<parallel::distributed::Vector<value_type> >  &src,
                        const std::pair<unsigned int,unsigned int>    &cell_range) const
{
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
  {
    this->mass_matrix_data_pml->phi[0].reinit(cell);
    this->mass_matrix_data_pml->phi[0].read_dof_values(src, 0);

    this->mass_matrix_data_pml->inverse.fill_inverse_JxW_values(this->mass_matrix_data_pml->coefficients);
    this->mass_matrix_data_pml->inverse.apply(this->mass_matrix_data_pml->coefficients, dim+dim+1,
        this->mass_matrix_data_pml->phi[0].begin_dof_values(),
        this->mass_matrix_data_pml->phi[0].begin_dof_values());

    this->mass_matrix_data_pml->phi[0].set_dof_values(dst,0);
  }
}

template<int dim, int fe_degree, typename Number>
void WaveEquationOperation<dim, fe_degree, Number>::
apply(const std::vector<parallel::distributed::Vector<value_type> >  &src,
      std::vector<parallel::distributed::Vector<value_type> >        &dst,
      const double                                              &cur_time,
      const double                                                    &dt) const
{
  for (unsigned int d=0; d<=dim; ++d)
    dst[d] = 0.0;

  Timer timer;
  this->time = cur_time;
  this->dirichlet_boundary_conditions->set_time(this->time);
  this->source_term->set_time(this->time);

  this->data.loop (&WaveEquationOperation<dim, fe_degree, Number>::local_apply_domain,
                   &WaveEquationOperation<dim, fe_degree, Number>::local_apply_face,
                   &WaveEquationOperation<dim, fe_degree, Number>::local_apply_boundary_face,
                   this, dst, src);

  this->computing_times[0] += timer.wall_time();
  timer.restart();

  this->data.cell_loop(&WaveEquationOperation<dim, fe_degree, Number>::local_apply_mass_matrix,
                       this, dst, dst);

  this->computing_times[1] += timer.wall_time();
  this->computing_times[2] += 1.;
}

template<int dim, int fe_degree, typename Number>
void WaveEquationOperationAcousticWave<dim, fe_degree, Number>::
apply(const std::vector<parallel::distributed::Vector<value_type> >  &src,
      std::vector<parallel::distributed::Vector<value_type> >        &dst,
      const double                                              &cur_time,
      const double                                                    &dt) const
{
  Timer timer;
  this->time = cur_time;
  this->dirichlet_boundary_conditions->set_time(this->time);
  this->source_term->set_time(this->time);

  this->data.loop (&WaveEquationOperationAcousticWave<dim, fe_degree, Number>::local_apply_domain,
                   &WaveEquationOperationAcousticWave<dim, fe_degree, Number>::local_apply_face,
                   &WaveEquationOperationAcousticWave<dim, fe_degree, Number>::local_apply_boundary_face,
                   this, dst, src);

  this->computing_times[0] += timer.wall_time();
  timer.restart();

  this->data.cell_loop(&WaveEquationOperationAcousticWave<dim, fe_degree, Number>::local_apply_mass_matrix,
                       this, dst, dst);

  this->computing_times[1] += timer.wall_time();
  this->computing_times[2] += 1.;
}

template<int dim, int fe_degree, typename Number>
void WaveEquationOperation<dim, fe_degree, Number>::
compute_post_gradient(const std::vector<parallel::distributed::Vector<value_type> >  &src,
                      std::vector<parallel::distributed::Vector<value_type> >        &dst,
                      const double                                              &cur_time,
                      const double                                                    &dt) const
{
  Timer timer;
  this->time = cur_time;
  this->dirichlet_boundary_conditions->set_time(this->time);
  this->source_term->set_time(this->time);

  this->data.loop (&WaveEquationOperation<dim, fe_degree, Number>::local_apply_domain,
                   &WaveEquationOperation<dim, fe_degree, Number>::local_apply_face,
                   &WaveEquationOperation<dim, fe_degree, Number>::local_apply_boundary_face,
                   this, dst, src);

  this->computing_times[0] += timer.wall_time();
  timer.restart();

  this->data.cell_loop(&WaveEquationOperation<dim, fe_degree, Number>::local_apply_mass_matrix,
                       this, dst, dst);

  this->computing_times[1] += timer.wall_time();
  this->computing_times[2] += 1.;
}


template<int dim, int fe_degree, typename Number>
WaveEquationOperationAcousticWavePML<dim,fe_degree,Number>::
WaveEquationOperationAcousticWavePML(const std::vector<const DoFHandler<dim> *> &dof_handlers,
                                  Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                                  Teuchos::RCP<Function<dim> > boundary_conditions,
                                  Teuchos::RCP<Function<dim> > source_term,
                                  Teuchos::RCP<AttenuationPML<dim,Number> > sigma_fct,
                                  value_type time_step_in,
                                  int sourceno,
                                  Teuchos::RCP<PATMonitorManager> monitormanagerin,
                                  bool adjoint)
  :
  WaveEquationOperationAcousticWave<dim,fe_degree,Number>(dof_handlers,discret,boundary_conditions,source_term,time_step_in,sourceno,monitormanagerin),
  sigma_pml(sigma_fct)
{
  for (unsigned int i=0; i<this->data.n_macro_cells()+this->data.n_macro_ghost_cells(); ++i)
  {
    layer_reference.push_back(std::vector<std::vector<int> >(1,std::vector<int>()));
    for (unsigned int v=0; v<this->data.n_components_filled(i); ++v)
    {
      layer_reference[i].push_back(std::vector<int>());
      for (unsigned int layer = 0; layer < sigma_fct->get_n_layer(); ++layer)
      {
        bool layer_flag = false;
        for (unsigned int k = 0; k < GeometryInfo<dim>::vertices_per_cell; ++k)
          if (sigma_fct->is_layer_active(layer, this->data.get_cell_iterator(i,v)->vertex(k)))
            layer_flag = true;
        if (layer_flag)
          layer_reference[i][v].push_back(layer);
      }
    }
  }

  if(adjoint)
  {
    inner_face_monitored.resize(this->data.n_macro_inner_faces());

    std::vector<DRT::Condition*> pressmonBC;
    discret->GetCondition("PressureMonitor",pressmonBC);

    for (unsigned int f=0; f<this->data.n_macro_inner_faces(); ++f)
    {
      for (unsigned int v=0; v<VectorizedArray<value_type>::n_array_elements && this->data.faces[f].left_cell[v] != numbers::invalid_unsigned_int; ++v)
      {
        const unsigned int left_cell_index_non_vectorized = this->data.faces[f].left_cell[v];
        const int left_element_index = this->data.get_cell_iterator(left_cell_index_non_vectorized / VectorizedArray<value_type>::n_array_elements,
            left_cell_index_non_vectorized % VectorizedArray<value_type>::n_array_elements)->index();
        const unsigned int right_cell_index_non_vectorized = this->data.faces[f].right_cell[v];
        const int right_element_index = this->data.get_cell_iterator(right_cell_index_non_vectorized / VectorizedArray<value_type>::n_array_elements,
            right_cell_index_non_vectorized % VectorizedArray<value_type>::n_array_elements)->index();

        DRT::Element* leftele = discret->lColElement(left_element_index);
        DRT::Element* rightele = discret->lColElement(right_element_index);

        // find the nodes they share
        std::vector<int> sharednodes;
        for(int n=0; n<leftele->NumNode(); ++n)
        {
          int nodeid = leftele->NodeIds()[n];
          for(int m=0; m<rightele->NumNode(); ++m)
          {
            if(nodeid == rightele->NodeIds()[m])
              sharednodes.push_back(nodeid);
          }
        }
        if(sharednodes.size()!=GeometryInfo<dim>::vertices_per_face)
          dserror("two neighboring elements share less nodes than they should");

        // are all shared nodes part of monitor?
        unsigned int nodepartofmon = 0;
        for(unsigned int i=0; i<sharednodes.size(); ++i)
          for(unsigned int j=0; j<pressmonBC.size(); ++j)
            if(pressmonBC[j]->ContainsNode(int(sharednodes[i])))
              nodepartofmon++;

        unsigned int count = 0;
        if(sharednodes.size() == nodepartofmon)
        {
          for(unsigned int i=0; i<sharednodes.size(); ++i)
          {
            count++;
            inner_face_monitored[f][v] = true;
          }
        }

      }
    }
  }


}

template<int dim, int fe_degree, typename Number>
void
WaveEquationOperationAcousticWavePML<dim,fe_degree,Number>::
read_initial_conditions(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                        std::vector<parallel::distributed::Vector<value_type> > &dst)
                        {

  FEEvaluation<dim,fe_degree,fe_degree+1,dim+dim+1,value_type> phi(this->data);
  for (unsigned int j=0; j<phi.dofs_per_cell; ++j)
    phi.submit_dof_value(Tensor<1,dim+dim+1,VectorizedArray<value_type> >(), j);

  unsigned int dofs_per_cell = phi.dofs_per_cell; // i assume, that it is the same for all cells

  unsigned int nodes_per_cell = GeometryInfo< dim >::vertices_per_cell;
  std::vector<Point<dim> > baci_vals_loc(nodes_per_cell);
  std::vector<Point<dim> > deal_vals_loc(nodes_per_cell);

  for (unsigned int i=0; i<this->data.n_macro_cells()+this->data.n_macro_ghost_cells(); ++i)
  {
    phi.reinit(i);
    for (unsigned int v=0; v<this->data.n_components_filled(i); ++v)
    {
      const int element_index = this->data.get_cell_iterator(i,v)->index();
      DRT::ELEMENTS::Acou * acouele = dynamic_cast<DRT::ELEMENTS::Acou*>(discret->lColElement(element_index));
      if (acouele == NULL)
        dserror("No acoustic element given!");

      // perform permutation: step 1: get the node coordinates
      for (unsigned int n=0; n<nodes_per_cell; ++n)
      {
        for(int d=0; d<dim; ++d)
        {
          deal_vals_loc[n](d) = this->data.get_cell_iterator(i,v)->vertex(n)(d);
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
          this->dofpermutations[element_index] = 0;
          // everything is alright
          for (unsigned j=0; j<dofs_per_cell; ++j)
          {
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[d*dofs_per_cell+j][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+j);
            phi.begin_dof_values()[dim*dofs_per_cell+j][v] = acouele->eleinteriorPressnp_(j);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[(dim+1+d)*dofs_per_cell+j][v] = acouele->eleinteriorAuxiliaryPML_(d*dofs_per_cell+j);
          }
        }
        else if(deal_vals_loc[0].distance(baci_vals_loc[3])<1e-10 &&
            deal_vals_loc[1].distance(baci_vals_loc[0])<1e-10 &&
            deal_vals_loc[2].distance(baci_vals_loc[2])<1e-10 &&
            deal_vals_loc[3].distance(baci_vals_loc[1])<1e-10)
        {
          this->dofpermutations[element_index] = 1;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            int permute = this->permutevalues(1,i);
            phi.begin_dof_values()[dim*dofs_per_cell+i][v] = acouele->eleinteriorPressnp_(permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[d*dofs_per_cell+i][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[(dim+1+d)*dofs_per_cell+i][v] = acouele->eleinteriorAuxiliaryPML_(d*dofs_per_cell+permute);
          }
        }
        else if(deal_vals_loc[0].distance(baci_vals_loc[2])<1e-10 &&
            deal_vals_loc[1].distance(baci_vals_loc[3])<1e-10 &&
            deal_vals_loc[2].distance(baci_vals_loc[1])<1e-10 &&
            deal_vals_loc[3].distance(baci_vals_loc[0])<1e-10)
        {
          this->dofpermutations[element_index] = 2;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            int permute = this->permutevalues(2,i);
            phi.begin_dof_values()[dim*dofs_per_cell+i][v] = acouele->eleinteriorPressnp_(permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[d*dofs_per_cell+i][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[(dim+1+d)*dofs_per_cell+i][v] = acouele->eleinteriorAuxiliaryPML_(d*dofs_per_cell+permute);
          }
        }
        else if(deal_vals_loc[0].distance(baci_vals_loc[1])<1e-10 &&
            deal_vals_loc[1].distance(baci_vals_loc[2])<1e-10 &&
            deal_vals_loc[2].distance(baci_vals_loc[0])<1e-10 &&
            deal_vals_loc[3].distance(baci_vals_loc[3])<1e-10)
        {
          this->dofpermutations[element_index] = 3;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            int permute = this->permutevalues(3,i);
            phi.begin_dof_values()[dim*dofs_per_cell+i][v] = acouele->eleinteriorPressnp_(permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[d*dofs_per_cell+i][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[(dim+1+d)*dofs_per_cell+i][v] = acouele->eleinteriorAuxiliaryPML_(d*dofs_per_cell+permute);
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
          this->dofpermutations[element_index] = 0;
          // everything is alright
          for (unsigned j=0; j<dofs_per_cell; ++j)
          {
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[d*dofs_per_cell+j][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+j);
            phi.begin_dof_values()[dim*dofs_per_cell+j][v] = acouele->eleinteriorPressnp_(j);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[(dim+1+d)*dofs_per_cell+j][v] = acouele->eleinteriorAuxiliaryPML_(d*dofs_per_cell+j);
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
          this->dofpermutations[element_index] = 1;
          // negative rotation around x
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            int permute = this->permutevalues(1,i);
            phi.begin_dof_values()[dim*dofs_per_cell+i][v] = acouele->eleinteriorPressnp_(permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[d*dofs_per_cell+i][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[(dim+1+d)*dofs_per_cell+i][v] = acouele->eleinteriorAuxiliaryPML_(d*dofs_per_cell+permute);
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
          this->dofpermutations[element_index] = 2;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            int permute = this->permutevalues(2,i);
            phi.begin_dof_values()[dim*dofs_per_cell+i][v] = acouele->eleinteriorPressnp_(permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[d*dofs_per_cell+i][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[(dim+1+d)*dofs_per_cell+i][v] = acouele->eleinteriorAuxiliaryPML_(d*dofs_per_cell+permute);
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
          this->dofpermutations[element_index] = 3;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            int permute = this->permutevalues(3,i);
            phi.begin_dof_values()[dim*dofs_per_cell+i][v] = acouele->eleinteriorPressnp_(permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[d*dofs_per_cell+i][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[(dim+1+d)*dofs_per_cell+i][v] = acouele->eleinteriorAuxiliaryPML_(d*dofs_per_cell+permute);
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
          this->dofpermutations[element_index] = 4;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            int permute = this->permutevalues(4,i);
            phi.begin_dof_values()[dim*dofs_per_cell+i][v] = acouele->eleinteriorPressnp_(permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[d*dofs_per_cell+i][v] = acouele->eleinteriorVelnp_(d*dofs_per_cell+permute);
            for (unsigned int d=0; d<dim; ++d)
              phi.begin_dof_values()[(dim+1+d)*dofs_per_cell+i][v] = acouele->eleinteriorAuxiliaryPML_(d*dofs_per_cell+permute);
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

  for (unsigned int i=0; i<dim+dim+1; ++i)
  {
    dst[i].update_ghost_values();
  }

}




// explicit instantiations
template class WaveEquationOperation<2,1,double>;
template class WaveEquationOperation<2,2,double>;
template class WaveEquationOperation<2,3,double>;
template class WaveEquationOperation<2,4,double>;
template class WaveEquationOperation<2,5,double>;
template class WaveEquationOperation<2,6,double>;
template class WaveEquationOperation<3,1,double>;
template class WaveEquationOperation<3,2,double>;
template class WaveEquationOperation<3,3,double>;
template class WaveEquationOperation<3,4,double>;
template class WaveEquationOperation<3,5,double>;
template class WaveEquationOperation<3,6,double>;
template class WaveEquationOperation<2,1,float>;
template class WaveEquationOperation<2,2,float>;
template class WaveEquationOperation<2,3,float>;
template class WaveEquationOperation<2,4,float>;
template class WaveEquationOperation<2,5,float>;
template class WaveEquationOperation<2,6,float>;
template class WaveEquationOperation<3,1,float>;
template class WaveEquationOperation<3,2,float>;
template class WaveEquationOperation<3,3,float>;
template class WaveEquationOperation<3,4,float>;
template class WaveEquationOperation<3,5,float>;
template class WaveEquationOperation<3,6,float>;

template class WaveEquationOperationAcousticWave<2,1,double>;
template class WaveEquationOperationAcousticWave<2,2,double>;
template class WaveEquationOperationAcousticWave<2,3,double>;
template class WaveEquationOperationAcousticWave<2,4,double>;
template class WaveEquationOperationAcousticWave<2,5,double>;
template class WaveEquationOperationAcousticWave<2,6,double>;
template class WaveEquationOperationAcousticWave<3,1,double>;
template class WaveEquationOperationAcousticWave<3,2,double>;
template class WaveEquationOperationAcousticWave<3,3,double>;
template class WaveEquationOperationAcousticWave<3,4,double>;
template class WaveEquationOperationAcousticWave<3,5,double>;
template class WaveEquationOperationAcousticWave<3,6,double>;
template class WaveEquationOperationAcousticWave<2,1,float>;
template class WaveEquationOperationAcousticWave<2,2,float>;
template class WaveEquationOperationAcousticWave<2,3,float>;
template class WaveEquationOperationAcousticWave<2,4,float>;
template class WaveEquationOperationAcousticWave<2,5,float>;
template class WaveEquationOperationAcousticWave<2,6,float>;
template class WaveEquationOperationAcousticWave<3,1,float>;
template class WaveEquationOperationAcousticWave<3,2,float>;
template class WaveEquationOperationAcousticWave<3,3,float>;
template class WaveEquationOperationAcousticWave<3,4,float>;
template class WaveEquationOperationAcousticWave<3,5,float>;
template class WaveEquationOperationAcousticWave<3,6,float>;

template class WaveEquationOperationAcousticWavePML<2,1,double>;
template class WaveEquationOperationAcousticWavePML<2,2,double>;
template class WaveEquationOperationAcousticWavePML<2,3,double>;
template class WaveEquationOperationAcousticWavePML<2,4,double>;
template class WaveEquationOperationAcousticWavePML<2,5,double>;
template class WaveEquationOperationAcousticWavePML<2,6,double>;
template class WaveEquationOperationAcousticWavePML<3,1,double>;
template class WaveEquationOperationAcousticWavePML<3,2,double>;
template class WaveEquationOperationAcousticWavePML<3,3,double>;
template class WaveEquationOperationAcousticWavePML<3,4,double>;
template class WaveEquationOperationAcousticWavePML<3,5,double>;
template class WaveEquationOperationAcousticWavePML<3,6,double>;
template class WaveEquationOperationAcousticWavePML<2,1,float>;
template class WaveEquationOperationAcousticWavePML<2,2,float>;
template class WaveEquationOperationAcousticWavePML<2,3,float>;
template class WaveEquationOperationAcousticWavePML<2,4,float>;
template class WaveEquationOperationAcousticWavePML<2,5,float>;
template class WaveEquationOperationAcousticWavePML<2,6,float>;
template class WaveEquationOperationAcousticWavePML<3,1,float>;
template class WaveEquationOperationAcousticWavePML<3,2,float>;
template class WaveEquationOperationAcousticWavePML<3,3,float>;
template class WaveEquationOperationAcousticWavePML<3,4,float>;
template class WaveEquationOperationAcousticWavePML<3,5,float>;
template class WaveEquationOperationAcousticWavePML<3,6,float>;
}


#endif // HAVE_DEAL_II
