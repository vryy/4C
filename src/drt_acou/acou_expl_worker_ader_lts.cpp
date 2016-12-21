/*!----------------------------------------------------------------------
\file acou_expl_worker_ader_lts.cpp
\brief Control routine for acoustic explicit time integration with ADER LTS

<pre>
\level 3

\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*----------------------------------------------------------------------*/

#include "acou_expl_worker.H"

#ifdef HAVE_DEAL_II

#include <deal.II/lac/parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>
//#include "fe_evaluation.h"
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

int evaluatefacecount = 0;
int evaluateboundaryfacecount = 0;
int evaluatecellcount = 0;
int evaluateupdatecount = 0;

template<int dim, int fe_degree, typename Number>
WaveEquationOperationAcousticWaveADERLTS<dim,fe_degree,Number>::
WaveEquationOperationAcousticWaveADERLTS(const DoFHandler<dim> &dof_handler,
                                         Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                                         Teuchos::RCP<Function<dim> > boundary_conditions,
                                         Teuchos::RCP<Function<dim> > source_term,
                                         int sourceno,
                                         Teuchos::RCP<Epetra_MultiVector> source_adjoint)
  :
  WaveEquationOperation<dim,fe_degree,Number>(dof_handler,discret,boundary_conditions,source_term,sourceno,source_adjoint)
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
      Teuchos::RCP<MAT::Material> mat = discret->lColElement(element_index)->Material();

      MAT::AcousticMat* actmat = static_cast<MAT::AcousticMat*>(mat.get());
      this->densities[i][v] = actmat->Density();
      this->speeds[i][v] = actmat->SpeedofSound();
    }
  }

  timesteps.resize(this->data.n_macro_cells()+this->data.n_macro_ghost_cells());
  neighborcells_i.resize(this->data.n_macro_cells()+this->data.n_macro_ghost_cells());
  adjacentfaces_i.resize(this->data.n_macro_cells()+this->data.n_macro_ghost_cells());
  for(unsigned int i=0; i<neighborcells_i.size(); ++i)
  {
    neighborcells_i[i].resize(2*dim,-1);
    adjacentfaces_i[i].resize(2*dim,-1);
  }
  flux_memory.resize(dim+1);
  this->data.initialize_dof_vector(flux_memory[0]);
  for(unsigned int d=1; d<flux_memory.size(); ++d)
    flux_memory[d] = flux_memory[0];

  const unsigned dofs_per_cell = this->data.get_dof_handler().get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  timesindices.resize(this->data.n_macro_cells()+this->data.n_macro_ghost_cells());

  IndexSet elerowset(flux_memory[0].size()/dofs_per_cell);
  IndexSet elecolset(flux_memory[0].size()/dofs_per_cell);
  for(unsigned int cell=0; cell<this->data.n_macro_cells(); ++cell)
  {
    this->data.get_cell_iterator(cell,0)->get_dof_indices(local_dof_indices);
    elerowset.add_index(local_dof_indices[0]/dofs_per_cell);
    timesindices[cell] = local_dof_indices[0]/dofs_per_cell;
  }
  for(unsigned int cell=this->data.n_macro_cells(); cell<this->data.n_macro_cells()+this->data.n_macro_ghost_cells(); ++cell)
  {
    this->data.get_cell_iterator(cell,0)->get_dof_indices(local_dof_indices);
    elecolset.add_index(local_dof_indices[0]/dofs_per_cell);
    timesindices[cell] = local_dof_indices[0]/dofs_per_cell;
  }

  const Epetra_MpiComm* mpi_comm = dynamic_cast<const Epetra_MpiComm*>(&(discret->Comm()));
  times.reinit(elerowset,elecolset,mpi_comm->Comm());
  oldtimes.reinit(elerowset,elecolset,mpi_comm->Comm());
  times = 0.0;
  oldtimes = 0.0;

  for (unsigned int i=0; i<this->data.n_macro_cells()+this->data.n_macro_ghost_cells(); ++i)
  {
    for (unsigned int v=0; v<this->data.n_components_filled(i); ++v)
    {
      const int element_index = this->data.get_cell_iterator(i,v)->index();
      double h = this->data.get_cell_iterator(i,v)->minimum_vertex_distance();
      //double h = this->data.get_cell_iterator(i,v)->extent_in_direction(0);
      double c = this->speeds[i][v];
      double k = discret->lColElement(element_index)->Degree();
      double acttimestep = 0.1/dim/c/k*h; //0.1/c/k*h; //TODO rausfinden welche CFL Zahl Stabilitaetsgrenze fuer ADER ist
      timesteps[i] = acttimestep;
    }
  }

  for(unsigned int f=0; f<this->data.faces.size(); ++f)
  {
    int iright, ileft;
    iright = this->data.faces[f].right_cell[0];
    ileft  = this->data.faces[f].left_cell[0];

    if(this->data.faces[f].right_cell[0]!=numbers::invalid_unsigned_int && this->data.faces[f].left_cell[0]!=numbers::invalid_unsigned_int)
    {
      // find first free entry
      int free_l = -1;
      for(int i=0; i<2*dim; ++i)
        if(neighborcells_i[ileft][i]==-1)
        {
          free_l=i;
          break;
        }
      int free_r = -1;
      for(int i=0; i<2*dim; ++i)
        if(neighborcells_i[iright][i]==-1)
        {
          free_r=i;
          break;
        }

      neighborcells_i[iright][free_r] = ileft;
      neighborcells_i[ileft][free_l] = iright;
    }
  }

  // fill the adjacent faces vectors
  for(unsigned int f=0; f<this->data.faces.size(); ++f)
  {
    if(this->data.faces[f].left_cell[0]!=numbers::invalid_unsigned_int)
    {
      // find first free entry
      int free_l = -1;
      for(int n=0; n<2*dim; ++n)
        if(adjacentfaces_i[this->data.faces[f].left_cell[0]][n]==-1)
        {
          free_l=n;
          break;
        }
      adjacentfaces_i[this->data.faces[f].left_cell[0]][free_l] = f;
    }
    if(this->data.faces[f].right_cell[0]!=numbers::invalid_unsigned_int)
    {
      int free_r = -1;
      for(int n=0; n<2*dim; ++n)
        if(adjacentfaces_i[this->data.faces[f].right_cell[0]][n]==-1)
        {
          free_r=n;
          break;
        }

      adjacentfaces_i[this->data.faces[f].right_cell[0]][free_r] = f;
    }
  }
}

template<int dim, int fe_degree, typename Number>
void WaveEquationOperationAcousticWaveADERLTS<dim,fe_degree,Number>::
local_apply_aderlts(std::vector<parallel::distributed::Vector<value_type> >             &dst,
                                std::vector<parallel::distributed::Vector<value_type> > &tempsrc,
                                std::vector<parallel::distributed::Vector<value_type> > &src,
                                const value_type                                        eps,
                                const value_type                                        timetoreach) const
{
  for(unsigned int cell=0; cell<this->data.n_macro_cells(); ++cell)
  {
    // check update criterion for all cells in this cell vectorization
    bool evaluate = true;

    // either real update or fill the last time step
    double thisupdatetime = times(timesindices[cell])+timesteps[cell];
    if(thisupdatetime>timetoreach)
      thisupdatetime = timetoreach;

    for(unsigned int n=0; n<2*dim; ++n)
      if(neighborcells_i[cell][n]>=0)
      {
        double neighborupdatetime = times(timesindices[neighborcells_i[cell][n]])+timesteps[neighborcells_i[cell][n]];
        if(neighborupdatetime>timetoreach)
          neighborupdatetime = timetoreach;

        if(thisupdatetime-eps>neighborupdatetime)
          evaluate = false;
      }

    if(evaluate)
    {
      // have to init to allow correct evaluation of next cell
      for(unsigned int d=0; d<dst.size(); ++d)
      {
        tempsrc[d] = 0.;
        dst[d] = 0.;
      }

      // get timings
      const value_type act_time = times(timesindices[cell]);
      value_type act_dt = timesteps[cell];

      // check if we would update further than the desired time level and adapt time step only for this evaluation
      if(act_time+act_dt>timetoreach)
        act_dt = timetoreach-act_time;

      //{  evaluate Psi***_f on adjacent faces
      // get indices of adjacent faces
      std::vector<int> faces(2*dim);
      for(unsigned int n=0; n<2*dim; ++n)
        faces[n] = adjacentfaces_i[cell][n];

      for(unsigned int f=0; f<faces.size(); ++f)
      {
        if(faces[f]<int(this->data.n_macro_inner_faces()) || faces[f]>=int(this->data.n_macro_inner_faces()+this->data.n_macro_boundary_faces())) // inner faces or ghosted faces
        {
          // determinte t1 and t2 for this face
          value_type t1 = times(timesindices[this->data.faces[faces[f]].left_cell[0]]);
          t1 = t1>times(timesindices[this->data.faces[faces[f]].right_cell[0]]) ? t1 : times(timesindices[this->data.faces[faces[f]].right_cell[0]]);

          // when evaluating t2 also consider timetoreach
          value_type t2left = times(timesindices[this->data.faces[faces[f]].left_cell[0]]) + timesteps[this->data.faces[faces[f]].left_cell[0]];
          if(t2left>timetoreach) t2left = timetoreach;
          value_type t2right = times(timesindices[this->data.faces[faces[f]].right_cell[0]])+timesteps[this->data.faces[faces[f]].right_cell[0]];
          if(t2right>timetoreach) t2right = timetoreach;
          value_type t2 = t2left>t2right ? t2right : t2left;

          if(t2-t1<eps)
          {
            // do not evaluate
          }
          else
            evaluateface(faces[f],dst,tempsrc,src,cell,t1,t2,false);
        }
        else // boundary faces
        {
          if(act_dt<eps)
          {
            // do not evaluate
          }
          else
            evaluateboundaryface(faces[f],dst,tempsrc,src,cell,act_time,act_time+act_dt);
        }
      }
      //}

      //{ 4.) evalutae Psi*_e and Psi****_e and do the update
      for(unsigned int d=0; d<dst.size(); ++d)
        tempsrc[d] = 0.;
      evaluatecell(cell,tempsrc,src,act_time,act_time+act_dt,false);
      evaluateupdate(cell,dst,src,tempsrc,act_dt);
      //}

    } // if(evaluate)
    else
    {

    } // else ** if(evaluate)
  }

}

template <int dim, int fe_degree, typename Number>
void
WaveEquationOperationAcousticWaveADERLTS<dim,fe_degree,Number>::
evaluatecell(unsigned int                                                   cell,
             std::vector<parallel::distributed::Vector<value_type> >        &dst,
             const std::vector<parallel::distributed::Vector<value_type> >  &src,
             const value_type                                               t1,
             const value_type                                               t2,
             bool                                                           fluxcomp) const
{
  evaluatecellcount++;

  // this function evaluates Psi**_e und Psi**_en for the cells with index cell in the time interval t1 to t2

  // for calculation of higher spatial derivatives
  //{
  const unsigned int n_q_points = FEEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type>::n_q_points;
  internal::InverseMassMatrixData<dim,fe_degree,dim+1,value_type>& mass_data = this->mass_matrix_data.get();
  FEEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> &phi_eval = mass_data.phi[0];
  //}

  // get all cell quanitites:
  //{
  // velocity and pressure read together
  phi_eval.reinit(cell);
  phi_eval.read_dof_values(src, 0);
  phi_eval.evaluate (true, true, false);

  // and material coefficients
  const VectorizedArray<value_type> rho = this->densities[cell];
  const VectorizedArray<value_type> rho_inv = 1./this->densities[cell];
  const VectorizedArray<value_type> c_sq = this->speeds[cell]*this->speeds[cell];
  value_type act_time = times(timesindices[cell]);
  if(fluxcomp)
    act_time = oldtimes(timesindices[cell]);
  //}

  // create container for quass point pressure and velocity
  // contributions. vcontrib will be set to zero automatically through the
  // Tensor<1,dim> constructor, but pcontrib must be set manually -> this is
  // done in the first loop below
  Tensor<1,dim,VectorizedArray<value_type> > vcontrib[n_q_points];
  VectorizedArray<value_type>                pcontrib[n_q_points];

  // sum over all integration points
  for (unsigned int q=0; q<n_q_points; ++q)
  {
    // get the values at the gauss point
    const Tensor<1,dim+1,VectorizedArray<value_type> > v_and_p = phi_eval.get_value(q);
    const Tensor<1,dim+1,Tensor<1,dim,VectorizedArray<value_type> > > v_and_p_grad = phi_eval.get_gradient(q);

    // contribution from k=0
    for (unsigned int d=0; d<dim; ++d)
      vcontrib[q][d] = (t2-t1)*v_and_p[d];
    pcontrib[q] = (t2-t1)*v_and_p[dim];

    // add contribution from k=1
    vcontrib[q] -= ((t2-act_time)*(t2-act_time)-(t1-act_time)*(t1-act_time))/2.0*rho_inv*v_and_p_grad[dim];
    for(unsigned int d=0; d<dim; ++d)
      pcontrib[q] -= ((t2-act_time)*(t2-act_time)-(t1-act_time)*(t1-act_time))/2.0*c_sq*rho*v_and_p_grad[d][d];

    // add source term contribution from Cauchy-Kovalewski k=1
    //{
    VectorizedArray<value_type> rhs = make_vectorized_array<value_type>(0.0);
    if(this->source_term_no>=0)
      dserror("source term not yet implemented for ADER LTS");
    pcontrib[q] += ((t2-act_time)*(t2-act_time)-(t1-act_time)*(t1-act_time))/2.0*rhs*c_sq;
    //}

    // evaluate phi_1
    Tensor<1,dim+1,VectorizedArray<value_type> > temp;
    for(unsigned int d=0; d<dim; ++d)
    {
      temp[d] = rho_inv*v_and_p_grad[dim][d];
      temp[dim] += c_sq*rho*v_and_p_grad[d][d];
    }
    phi_eval.submit_value(temp,q);
  }

  mass_data.inverse.fill_inverse_JxW_values(mass_data.coefficients);
  double fac = -0.5;
  VectorizedArray<value_type> fac_t;

  // all following contributions can be looped
  for(int k=2; k<=fe_degree; ++k)
  {
    fac /= -(k+1);
    fac_t = std::pow(t2-act_time,value_type(k+1))-std::pow(t1-act_time,value_type(k+1));

    // integrate over element
    phi_eval.integrate(true,false);

    // apply inverse mass matrix
    //{
    mass_data.inverse.apply(mass_data.coefficients, dim+1,
                            phi_eval.begin_dof_values(),
                            phi_eval.begin_dof_values());

    //}

    // evaulate this phi at the gauss points
    phi_eval.evaluate(false,true);

    // sum over all integration points
    for (unsigned int q=0; q<n_q_points; ++q)
    {
      // get the gauss point values
      const Tensor<1,dim+1,Tensor<1,dim,VectorizedArray<value_type> > > phi_gradient = phi_eval.get_gradient(q);

      // calculate contributions
      for(unsigned int d=0; d<dim; ++d)
        vcontrib[q][d] += fac*fac_t*rho_inv*phi_gradient[dim][d];
      for(unsigned int d=0; d<dim; ++d)
        pcontrib[q] += fac*fac_t*c_sq*rho*phi_gradient[d][d];

      // add source term contribution
      //{
      if(this->source_term_no>=0)
        dserror("source term not yet implemented for ADER LTS");
      //}

      // evaluate things phi_k+1 needs
      Tensor<1,dim+1,VectorizedArray<value_type> > temp;
      for(unsigned int d=0; d<dim; ++d)
      {
        temp[d] = rho_inv*phi_gradient[dim][d];
        temp[dim] += c_sq*rho*phi_gradient[d][d];
      }
      phi_eval.submit_value(temp,q);
    }
  }

  // submit what we collected to the field!
  for (unsigned int q=0; q<n_q_points; ++q)
  {
    Tensor<1,dim+1,VectorizedArray<value_type> > v_and_p;
    for (unsigned int d=0; d<dim; ++d)
      v_and_p[d] = vcontrib[q][d];
    v_and_p[dim] = pcontrib[q];
    phi_eval.submit_value(v_and_p, q);
  }
  phi_eval.integrate(true,false);

  // apply [A 0 ; 0 M]^{-1} to the field
  mass_data.inverse.apply(mass_data.coefficients, dim+1,
                          phi_eval.begin_dof_values(),
                          phi_eval.begin_dof_values());

  phi_eval.set_dof_values(dst,0);

  return;
}

template <int dim, int fe_degree, typename Number>
void
WaveEquationOperationAcousticWaveADERLTS<dim,fe_degree,Number>::
evaluateface(int                                                            face,
             std::vector<parallel::distributed::Vector<value_type> >        &dst,
             std::vector<parallel::distributed::Vector<value_type> >        &tempsrc,
             const std::vector<parallel::distributed::Vector<value_type> >  &src,
             unsigned int                                                   callingcell,
             value_type                                                     t1,
             value_type                                                     t2,
             bool                                                           fluxonly) const
{
  evaluatefacecount++;

  evaluatecell(this->data.faces[face].left_cell[0],tempsrc,src,t1,t2,fluxonly);
  evaluatecell(this->data.faces[face].right_cell[0],tempsrc,src,t1,t2,fluxonly);

  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> phi(this->data, true, 0, 0, true);
  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> phi_neighbor(this->data, false, 0, 0, true);

  phi.reinit(face);
  phi.read_dof_values(tempsrc, 0);
  phi.evaluate(true,false);
  const VectorizedArray<value_type> rho_plus = phi.read_cell_data(this->densities);
  const VectorizedArray<value_type> rho_inv_plus = 1./rho_plus;
  const VectorizedArray<value_type> c_plus = phi.read_cell_data(this->speeds);
  const VectorizedArray<value_type> c_sq_plus = c_plus * c_plus;
  const VectorizedArray<value_type> tau_plus = 1./c_plus/rho_plus;
  phi_neighbor.reinit(face);
  phi_neighbor.read_dof_values(tempsrc, 0);
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
      dserror("i don't know yet");

    for (unsigned int d=0; d<dim; ++d)
    {
      val_plus[d] = -pres_diff_plus*normal[d];
      val_minus[d] = pres_diff_minus*normal[d];
    }
    if(this->adjoint_eval==false)
    {
      val_plus[dim] = c_sq_plus * rho_plus * (normal_v_plus + tau_plus * (val_plus[dim] - lambda));
      val_minus[dim] = c_sq_minus * rho_minus * (normal_v_minus + tau_minus * (val_minus[dim] - lambda));
    }
    else
    {
      dserror("i don't know yet");
    }
    phi.submit_value(val_plus, q);
    phi_neighbor.submit_value(val_minus, q);
  }
  phi.integrate(true,false);
  phi_neighbor.integrate(true,false);

  // standard case
  if(!fluxonly)
  {
    if(this->data.faces[face].left_cell[0]==callingcell)
    {
      phi.distribute_local_to_global(dst, 0);
      phi_neighbor.distribute_local_to_global(flux_memory, 0);
    }
    else
    {
      phi.distribute_local_to_global(flux_memory, 0);
      phi_neighbor.distribute_local_to_global(dst, 0);
    }
  }
  else // flux update
  {
    if(this->data.faces[face].left_cell[0]==callingcell)
      phi.distribute_local_to_global(flux_memory, 0);
    else
      phi_neighbor.distribute_local_to_global(flux_memory, 0);
  }

  return;
}


template <int dim, int fe_degree, typename Number>
void
WaveEquationOperationAcousticWaveADERLTS<dim,fe_degree,Number>::
evaluateboundaryface(int                                                            face,
                     std::vector<parallel::distributed::Vector<value_type> >        &dst,
                     std::vector<parallel::distributed::Vector<value_type> >        &tempsrc,
                     const std::vector<parallel::distributed::Vector<value_type> >  &src,
                     int                                                            callingcell,
                     value_type                                                     t1,
                     value_type                                                     t2) const
{
  evaluateboundaryfacecount++;

  evaluatecell(this->data.faces[face].left_cell[0],tempsrc,src,t1,t2,false);

  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> phi(this->data, true, 0, 0, true);

  // quantities we need in the loop
  Point<dim> point;
  std::vector<value_type> node_values(GeometryInfo<dim>::vertices_per_face);
  std::vector<std::vector<value_type> > node_coords(GeometryInfo<dim>::vertices_per_face);
  for(unsigned int n=0; n<GeometryInfo<dim>::vertices_per_face; ++n)
    node_coords[n].resize(dim);

  phi.reinit(face);
  phi.read_dof_values(tempsrc, 0);
  phi.evaluate(true,false);
  const VectorizedArray<value_type> rho = phi.read_cell_data(this->densities);
  const VectorizedArray<value_type> rho_inv = 1./rho;
  const VectorizedArray<value_type> c = phi.read_cell_data(this->speeds);
  const VectorizedArray<value_type> c_sq = phi.read_cell_data(this->speeds)*phi.read_cell_data(this->speeds);
  const VectorizedArray<value_type> tau = 1./phi.read_cell_data(this->speeds)/phi.read_cell_data(this->densities);

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

      if(this->source_adjoint_meas!=Teuchos::null && this->adjoint_eval == true) // second query required for intermediate integration
        dserror("inverse analysis not yet implemented for ader");
    }
    else if(int_boundary_id==2) // monitored and absorbing
    {
      if(this->adjoint_eval==false)
        lambda = tau/(tau+1./c/rho)*p_plus + 1./(tau+1./c/rho)*normal_v_plus;
      else
        lambda = 1./(tau+1./c/rho)*rho_inv*normal_v_plus - tau*rho*c_sq/(tau+1./c/rho)*p_plus;
      if(this->source_adjoint_meas!=Teuchos::null && this->adjoint_eval == true) // second query required for intermediate integration
        dserror("inverse analysis not yet implemented for ader");

    }
    else if(int_boundary_id==3) // free boundary
      lambda = 1./tau*normal_v_plus+p_plus; // VectorizedArray<value_type>();
    else if(int_boundary_id==4) // dbc from time reversal
    {
      if(this->source_adjoint_meas!=Teuchos::null)
      {
        //for (unsigned int v=0; v<n_vectorization && this->data.faces[face].left_cell[v] != numbers::invalid_unsigned_int; ++v)
        int v=0;
        {
          for (unsigned int d=0; d<dim; ++d)
            point[d] = q_point[d][v];

          for(unsigned int n=0; n<GeometryInfo<dim>::vertices_per_face; ++n)
          {
            for(unsigned int d=0; d<dim; ++d)
              node_coords[n][d] = this->table_node_coords(face-this->data.n_macro_inner_faces(),v,n,d);
            int gid = this->table_node_ids(face-this->data.n_macro_inner_faces(),v,n);
            int lid = this->source_adjoint_meas->Map().LID(gid);
            node_values[n] =  this->source_adjoint_meas->operator ()(this->timestep_source_number)->operator [](lid);
          }
          lambda[v] = this->time_step*this->evaluate_source_timereversal(point,node_coords,node_values);
        }
      }
    }
    else if(int_boundary_id>=5) // dbcs
    {
      if(this->adjoint_eval==false)
      {
        //for (unsigned int v=0; v<n_vectorization; ++v)
        int v=0;
        {
          Point<dim> point;
          for (unsigned int d=0; d<dim; ++d)
            point[d] = q_point[d][v];
          lambda[v] = this->time_step*this->dirichlet_boundary_conditions->value(point,(int_boundary_id-5)*dim); // "time integral of dirichlet value"
        }
      }
      else
        lambda = VectorizedArray<value_type>();
    }

    for (unsigned int d=0; d<dim; ++d)
      val_plus[d] = -(p_plus - lambda)*normal[d]*rho_inv;
    val_plus[dim] = c_sq*rho*(normal_v_plus - tau*(lambda - p_plus));

    phi.submit_value(val_plus,q);
  }

  phi.integrate(true,false);
  phi.distribute_local_to_global(dst, 0);

  return;
}

template<int dim, int fe_degree, typename Number>
void WaveEquationOperationAcousticWaveADERLTS<dim,fe_degree,Number>::
evaluateupdate(unsigned int                                                   cell,
               std::vector<parallel::distributed::Vector<value_type> >        &dst,
               std::vector<parallel::distributed::Vector<value_type> >        &src,
               std::vector<parallel::distributed::Vector<value_type> >        &tempsrc,
               value_type                                                     act_dt) const
{
  evaluateupdatecount++;

  const unsigned int n_q_points = FEEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type>::n_q_points;
  internal::InverseMassMatrixData<dim,fe_degree,dim+1,value_type>& mass_data = this->mass_matrix_data.get();
  FEEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> &phi_eval = mass_data.phi[0];
  FEEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> help_eval(this->data); // for memory variable and update of src

  // calculate Psi****_e
  //{
  phi_eval.reinit(cell);
  phi_eval.read_dof_values(tempsrc, 0);
  phi_eval.evaluate (true, true, false);

  // and material coefficients
  const VectorizedArray<value_type> rho = this->densities[cell];
  const VectorizedArray<value_type> rho_inv = 1./this->densities[cell];
  const VectorizedArray<value_type> c_sq = this->speeds[cell]*this->speeds[cell];

  for (unsigned int q=0; q<n_q_points; ++q)
  {
    const Tensor<1,dim+1,Tensor<1,dim,VectorizedArray<value_type> > > v_and_p_grad = phi_eval.get_gradient(q);
    const Tensor<1,dim+1,VectorizedArray<value_type> > v_and_p = phi_eval.get_value(q);

    if(this->source_term_no>=0)
      dserror("not implemented");

    Tensor<1,dim+1,VectorizedArray<value_type> > temp_value;
    for(unsigned int d=0; d<dim; ++d)
      temp_value[d] = rho_inv*v_and_p_grad[dim][d];
    phi_eval.submit_value(temp_value,q);


    Tensor<1,dim+1,Tensor<1,dim,VectorizedArray<value_type> > > temp_gradient;
    for(unsigned int d=0; d<dim; ++d)
      temp_gradient[dim][d] = -rho*c_sq*v_and_p[d];

    phi_eval.submit_gradient(temp_gradient,q);
  }

  phi_eval.integrate (true, true);

  // add memory variable
  help_eval.reinit(cell);
  help_eval.read_dof_values(flux_memory, 0);
  unsigned int dofs_per_cell = phi_eval.dofs_per_cell;
  for (unsigned j=0; j<dofs_per_cell; ++j)
    for (unsigned int d=0; d<dim+1; ++d)
    {
      phi_eval.begin_dof_values()[d*dofs_per_cell+j] += help_eval.begin_dof_values()[d*dofs_per_cell+j];
      help_eval.begin_dof_values()[d*dofs_per_cell+j] = 0.;
    }
  help_eval.set_dof_values(flux_memory, 0); // tell the flux_memory variable, that some of its values are reset

  // add face contribution (stored in dst)
  help_eval.read_dof_values(dst, 0);
  for (unsigned j=0; j<dofs_per_cell; ++j)
    for (unsigned int d=0; d<dim+1; ++d)
    {
      phi_eval.begin_dof_values()[d*dofs_per_cell+j] += help_eval.begin_dof_values()[d*dofs_per_cell+j];
    }

  // apply inverse mass matrix
  mass_data.inverse.apply(mass_data.coefficients, dim+1,
                          phi_eval.begin_dof_values(),
                          phi_eval.begin_dof_values());
  //}

  // update time of the active element
  times(timesindices[cell]) += act_dt;

  // update element degrees of freedom
  help_eval.read_dof_values(src, 0);
  for (unsigned j=0; j<dofs_per_cell; ++j)
    for (unsigned int d=0; d<dim+1; ++d) // velocity and pressure
      help_eval.begin_dof_values()[d*dofs_per_cell+j] =  help_eval.begin_dof_values()[d*dofs_per_cell+j] - phi_eval.begin_dof_values()[d*dofs_per_cell+j];
  help_eval.set_dof_values(src,0);
}

template<int dim, int fe_degree, typename Number>
void WaveEquationOperationAcousticWaveADERLTS<dim, fe_degree, Number>::
apply(const std::vector<parallel::distributed::Vector<value_type> >  &src,
      std::vector<parallel::distributed::Vector<value_type> >  &dst,
      const double                                             &cur_time,
      const double                                             &dt) const
{
  Timer timer;
  this->time = cur_time;
  this->time_step = dt;
  this->dirichlet_boundary_conditions->set_time(this->time);
  this->source_term->set_time(this->time);
  value_type mintime = this->time;
  value_type eps = 1.e-11*this->time_step;
  value_type time_to_reach = this->time + this->time_step; // this is the time level every element should reach to write ouput together at same time level
  parallel::distributed::Vector<value_type> auxtimes;
  auxtimes.reinit(this->times);
  Epetra_MpiComm comm(src[0].get_mpi_communicator());

  std::vector<parallel::distributed::Vector<value_type> > tempsrc(dim+1);
  std::vector<parallel::distributed::Vector<value_type> > constsrc(dim+1);
  for(unsigned int d=0; d<dim+1; ++d)
  {
    tempsrc[d].reinit(src[d],true);
    tempsrc[d] = 0.;
    constsrc[d].reinit(src[d],true);
    constsrc[d] = src[d];
  }

  std::vector<parallel::distributed::Vector<value_type> > srcold(dim+1);
  for(unsigned int d=0; d<dim+1; ++d)
  {
    srcold[d].reinit(src[d],true);
    srcold[d] = 0.;
  }

  int count = 0;
  do{
    //if(!comm.MyPID())
     //std::cout<<"ader time loop: count "<<count<<", time "<<this->time<<", time_step "<<this->time_step<<", mintime "<<std::setprecision(16)<<mintime<<" eps "<<eps<<std::endl;
    count++;

    oldtimes = times;
    oldtimes.update_ghost_values();
    for(unsigned int d=0; d<dim+1; ++d)
    {
      srcold[d] = constsrc[d];
      srcold[d].update_ghost_values();
    }

    local_apply_aderlts(dst,tempsrc,constsrc,eps,time_to_reach);
    times.update_ghost_values();

    // find minimum value in times vector
    auxtimes = 10.*time_to_reach;
    auxtimes.sadd(1.0,-1.0,times);
    value_type inftytime = auxtimes.linfty_norm();
    mintime = 10.*time_to_reach - inftytime;

    // update ghost values in src vector
    for(unsigned int d=0; d<constsrc.size(); ++d)
      constsrc[d].update_ghost_values();

    // calculate the fluxes from the ghosts
    communicate_flux_memory(dst,tempsrc,srcold,eps);

  } while(count< std::numeric_limits<int>::max() && mintime+eps<=this->time+this->time_step); // as long as not all elements are at time+time_step
  if(count==std::numeric_limits<int>::max())
    dserror("something with the update of times went wrong, while loop performed %d repetitions",std::numeric_limits<int>::max());


  for(unsigned int d=0; d<dim+1; ++d)
  {
    dst[d] = constsrc[d];
  }

  // timinig
  this->computing_times[1] += timer.wall_time();
  this->computing_times[2] += 1.;

  std::cout<<"evaluatefacecount " << evaluatefacecount<<std::endl;
  std::cout<<"evaluateboundaryfacecount " << evaluateboundaryfacecount<<std::endl;
  std::cout<<"evaluatecellcount " << evaluatecellcount<<std::endl;
  std::cout<<"evaluateupdatecount " << evaluateupdatecount<<std::endl;


}




template<int dim, int fe_degree, typename Number>
void WaveEquationOperationAcousticWaveADERLTS<dim,fe_degree,Number>::
communicate_flux_memory(std::vector<parallel::distributed::Vector<value_type> >        &dst,
                        std::vector<parallel::distributed::Vector<value_type> >        &tempsrc,
                        const std::vector<parallel::distributed::Vector<value_type> >  &src,
                        const value_type                                                eps) const
{
  for(unsigned int d=0; d<dst.size(); ++d)
  {
    tempsrc[d] = 0.;
    dst[d] = 0.;
  }

  for(unsigned int cell=0; cell<this->data.n_macro_cells(); ++cell)
  {
    // determine ghost neighbors
    std::vector<int> ghostneighbors_i;
    for(unsigned int n=0; n<2*dim; ++n)
      if(neighborcells_i[cell][n]>=int(this->data.n_macro_cells()))
        ghostneighbors_i.push_back(neighborcells_i[cell][n]);

    if(ghostneighbors_i.size()>0) // has ghost neighbors
    {
      // determine the faces inbetween
      std::vector<int> faces(ghostneighbors_i.size());
      for(unsigned int g=0; g<ghostneighbors_i.size(); ++g)
        for(unsigned int n=0; n<2*dim; ++n)
          for(unsigned int m=0; m<2*dim; ++m)
            if(adjacentfaces_i[cell][n]==adjacentfaces_i[ghostneighbors_i[g]][m])
              faces[g] = adjacentfaces_i[cell][n];

      // determine the time interval for which we need flux calculation
      value_type t1 = 0.0;
      value_type t2 = 0.0;
      for(unsigned int g=0; g<ghostneighbors_i.size(); ++g)
      {
        if(times(timesindices[cell])<times(timesindices[ghostneighbors_i[g]])-eps)
        {
          t2 = times(timesindices[ghostneighbors_i[g]]);
          if(times(timesindices[cell])<oldtimes(timesindices[ghostneighbors_i[g]]))
            t1 = oldtimes(timesindices[ghostneighbors_i[g]]);
          else
            t1 = times(timesindices[cell]);

          if(t2-t1>eps)
            evaluateface(faces[g],dst,tempsrc,src,cell,t1,t2,true);
        }
      }
    }
  }

}


// explicit instantiations
template class WaveEquationOperationAcousticWaveADERLTS<2,1,double>;
template class WaveEquationOperationAcousticWaveADERLTS<2,2,double>;
template class WaveEquationOperationAcousticWaveADERLTS<2,3,double>;
template class WaveEquationOperationAcousticWaveADERLTS<2,4,double>;
template class WaveEquationOperationAcousticWaveADERLTS<2,5,double>;
template class WaveEquationOperationAcousticWaveADERLTS<2,6,double>;
template class WaveEquationOperationAcousticWaveADERLTS<3,1,double>;
template class WaveEquationOperationAcousticWaveADERLTS<3,2,double>;
template class WaveEquationOperationAcousticWaveADERLTS<3,3,double>;
template class WaveEquationOperationAcousticWaveADERLTS<3,4,double>;
template class WaveEquationOperationAcousticWaveADERLTS<3,5,double>;
template class WaveEquationOperationAcousticWaveADERLTS<3,6,double>;
template class WaveEquationOperationAcousticWaveADERLTS<2,1,float>;
template class WaveEquationOperationAcousticWaveADERLTS<2,2,float>;
template class WaveEquationOperationAcousticWaveADERLTS<2,3,float>;
template class WaveEquationOperationAcousticWaveADERLTS<2,4,float>;
template class WaveEquationOperationAcousticWaveADERLTS<2,5,float>;
template class WaveEquationOperationAcousticWaveADERLTS<2,6,float>;
template class WaveEquationOperationAcousticWaveADERLTS<3,1,float>;
template class WaveEquationOperationAcousticWaveADERLTS<3,2,float>;
template class WaveEquationOperationAcousticWaveADERLTS<3,3,float>;
template class WaveEquationOperationAcousticWaveADERLTS<3,4,float>;
template class WaveEquationOperationAcousticWaveADERLTS<3,5,float>;
template class WaveEquationOperationAcousticWaveADERLTS<3,6,float>;
}


#endif // HAVE_DEAL_II
