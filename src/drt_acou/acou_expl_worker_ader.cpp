/*!----------------------------------------------------------------------
\file acou_expl_worker_ader.cpp
\brief Control routine for acoustic explicit time integration with ADER
\level 3

<pre>
\level 3

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
//#include <deal.II/matrix_free/fe_evaluation.h>
#include "fe_evaluation.h"
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

template<int dim, int fe_degree, typename Number>
WaveEquationOperationAcousticWaveADER<dim,fe_degree,Number>::
WaveEquationOperationAcousticWaveADER(const DoFHandler<dim> &dof_handler,
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

  /*flux_memory.resize(dim+1);
  this->data.initialize_dof_vector(flux_memory[0]);
  for(unsigned int d=1; d<flux_memory.size(); ++d)
    flux_memory[d] = flux_memory[0];*/
}

template<int dim, int fe_degree, typename Number>
void WaveEquationOperationAcousticWaveADER<dim,fe_degree,Number>::
local_apply_firstader_domain(const MatrixFree<dim,value_type>                                &data,
                             std::vector<parallel::distributed::Vector<value_type> >         &dst,
                             const std::vector<parallel::distributed::Vector<value_type> >   &src,
                             const std::pair<unsigned int,unsigned int>                 &cell_range) const
{
  // for calculation of higher spatial derivatives
  //{
  const unsigned int n_q_points = FEEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type>::n_q_points;
  internal::InverseMassMatrixData<dim,fe_degree,dim+1,value_type>& mass_data = this->mass_matrix_data.get();
  FEEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> &phi_eval = mass_data.phi[0];
  //}

  //FEEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> postpressure(data);
  // cell loop
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
  {
    // get all cell quanitites:
    //{
    // velocity and pressure read together
    phi_eval.reinit(cell);
    phi_eval.read_dof_values(src, 0);
    phi_eval.evaluate (true, true, false);

    //postpressure.reinit(cell);
    //postpressure.read_dof_values(flux_memory,0);
    //postpressure.evaluate(true,false,false);

    // and material coefficients
    const VectorizedArray<value_type> rho = this->densities[cell];
    const VectorizedArray<value_type> rho_inv = 1./this->densities[cell];
    const VectorizedArray<value_type> c_sq = this->speeds[cell]*this->speeds[cell];
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

      //const Tensor<1,dim+1,VectorizedArray<value_type> > postpressure_val = postpressure.get_value(q);

      // contribution from k=0
      for (unsigned int d=0; d<dim; ++d)
        vcontrib[q][d] = this->time_step*v_and_p[d];
      pcontrib[q] = this->time_step*v_and_p[dim];

      // add contribution from k=1
      vcontrib[q] -= this->time_step*this->time_step/2.0*rho_inv*v_and_p_grad[dim];
      //for(unsigned int d=0; d<dim; ++d)
      //  vcontrib[q][d] += this->time_step*this->time_step/2.0*rho_inv*postpressure_val[d];
      for(unsigned int d=0; d<dim; ++d)
        pcontrib[q] -= this->time_step*this->time_step/2.0*c_sq*rho*v_and_p_grad[d][d];
      //pcontrib[q] += this->time_step*this->time_step/2.0*c_sq*rho*postpressure_val[dim];

      // add source term contribution from Cauchy-Kovalewski k=1
      //{
      VectorizedArray<value_type> rhs = make_vectorized_array<value_type>(0.0);
      if(this->source_term_no>=0)
      {
        Point<dim,VectorizedArray<value_type> > q_points = phi_eval.quadrature_point(q);
        for (unsigned int n=0; n<rhs.n_array_elements; ++n)
        {
          double xyz[dim];
          for (unsigned int d=0; d<dim; ++d)
            xyz[d] = q_points[d][n];
          rhs[n] = DRT::Problem::Instance()->Funct(this->source_term_no).Evaluate(1,xyz,this->time); // FIRST component -> actual function of the source term f
        }
      }
      pcontrib[q] += this->time_step*this->time_step/2.0*rhs*c_sq;
      //}

      // evaluate phi_1
      Tensor<1,dim+1,VectorizedArray<value_type> > temp;
      for(unsigned int d=0; d<dim; ++d)
      {
        temp[d] = rho_inv*v_and_p_grad[dim][d];
        temp[dim] += c_sq*rho*v_and_p_grad[d][d];
      }
      phi_eval.submit_value(temp,q);
      //phi_eval.submit_value(-postpressure_val,q);
    }

    mass_data.inverse.fill_inverse_JxW_values(mass_data.coefficients);
    double fac = -this->time_step * this->time_step / 2.;
    // all following contributions can be looped
    for(int k=2; k<=fe_degree+1; ++k)
    {
      fac *= -this->time_step / (k+1);

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
          vcontrib[q][d] += fac*rho_inv*phi_gradient[dim][d];
        for(unsigned int d=0; d<dim; ++d)
          pcontrib[q] += fac*c_sq*rho*phi_gradient[d][d];

        // add source term contribution
        //{
        if(this->source_term_no>=0)
        {
          double facsource = std::abs(fac);

          Point<dim,VectorizedArray<value_type> > q_points = phi_eval.quadrature_point(q);
          int compindex = 2;
          for(int i=2; i<k; ++i)
            compindex += (i-1)*(dim+1)+1;

          // for s=0 the pure time derivative only for pressure
          VectorizedArray<value_type> rhsp = make_vectorized_array<value_type>(0.0);
          for (unsigned int n=0; n<rhsp.n_array_elements; ++n)
          {
            double xyz[dim];
            for (unsigned int d=0; d<dim; ++d)
              xyz[d] = q_points[d][n];
            rhsp[n] = DRT::Problem::Instance()->Funct(this->source_term_no).Evaluate(compindex,xyz,this->time); // (k-1)th time derivative of f
          }
          pcontrib[q] += facsource*rhsp*c_sq;
          compindex++;

          // for the higher s take the following components
          Tensor<1,dim,VectorizedArray<value_type> > rhsv;
          VectorizedArray<value_type> matscal = make_vectorized_array<value_type>(1.0);
          for(int s=1; s<=k-1;++s) // the "s" as from the paper from Dumbser
          {
            facsource *= -1.0;
            if(s%2==0)
              matscal *= c_sq*rho;
            else
              matscal *= rho_inv;

            for (unsigned int n=0; n<rhsp.n_array_elements; ++n)
            {
              double xyz[dim];
              for (unsigned int d=0; d<dim; ++d)
                xyz[d] = q_points[d][n];
              for (unsigned int d=0; d<dim; ++d)
                rhsv[d][n] = DRT::Problem::Instance()->Funct(this->source_term_no).Evaluate(compindex+d,xyz,this->time);
              rhsp[n] = DRT::Problem::Instance()->Funct(this->source_term_no).Evaluate(compindex+dim,xyz,this->time);
            }
            pcontrib[q] += facsource*matscal*rhsp*c_sq;
            vcontrib[q] += facsource*matscal*rhsv*c_sq;
            compindex += dim+1;
          }
        }
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

  } // for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)

}


template <int dim, int fe_degree, typename Number>
void
WaveEquationOperationAcousticWaveADER<dim,fe_degree,Number>::
local_apply_ader_face (const MatrixFree<dim,value_type> &,
  std::vector<parallel::distributed::Vector<value_type> >        &dst,
  const std::vector<parallel::distributed::Vector<value_type> >  &src,
  const std::pair<unsigned int,unsigned int>                     &face_range) const
{
  // basically the same as local_apply_face, but different signs in some places

  // There is some overhead in the methods in FEEvaluation, so it is faster
  // to combine pressure and velocity in the same object and just combine
  // them at the level of quadrature points
  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> phi(this->data, true, 0, 0, 0, true);
  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> phi_neighbor(this->data, false, 0, 0, 0, true);

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
        dserror("i don't know yet");
      }

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
    phi.distribute_local_to_global(dst, 0);

    phi_neighbor.integrate(true,false);
    phi_neighbor.distribute_local_to_global(dst, 0);
  }

}


template <int dim, int fe_degree, typename Number>
void WaveEquationOperationAcousticWaveADER<dim,fe_degree,Number>::
local_apply_ader_boundary_face (const MatrixFree<dim,value_type> &,
                           std::vector<parallel::distributed::Vector<value_type> >       &dst,
                           const std::vector<parallel::distributed::Vector<value_type> > &src,
                           const std::pair<unsigned int,unsigned int>               &face_range) const
{

  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> phi(this->data, true, 0, 0, 0, true);

  // quantities we need in the loop
  Point<dim> point;
  std::vector<value_type> node_values(GeometryInfo<dim>::vertices_per_face);
  std::vector<std::vector<value_type> > node_coords(GeometryInfo<dim>::vertices_per_face);
  for(unsigned int n=0; n<GeometryInfo<dim>::vertices_per_face; ++n)
    node_coords[n].resize(dim);

  for (unsigned int face=face_range.first; face<face_range.second; face++)
  {
    phi.reinit(face);
    phi.read_dof_values(src, 0);
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
          for (unsigned int v=0; v<VectorizedArray<value_type>::n_array_elements && this->data.faces[face].left_cell[v] != numbers::invalid_unsigned_int; ++v)
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
          for (unsigned int v=0; v<VectorizedArray<value_type>::n_array_elements; ++v)
          {
            Point<dim> point;
            for (unsigned int d=0; d<dim; ++d)
              point[d] = q_point[d][v];
            lambda[v] = this->time_step*this->dirichlet_boundary_conditions->value(point,(int_boundary_id-5)*dim); // "time integral of dirichlet value"
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
  }
}

template<int dim, int fe_degree, typename Number>
void WaveEquationOperationAcousticWaveADER<dim,fe_degree,Number>::
local_apply_secondader_domain(const MatrixFree<dim,value_type>                                &data,
                              std::vector<parallel::distributed::Vector<value_type> >         &dst,
                              const std::vector<parallel::distributed::Vector<value_type> >   &src,
                              const std::pair<unsigned int,unsigned int>                 &cell_range) const
{
  FEEvaluation<dim,fe_degree,fe_degree+1,dim,value_type> velocity(data);
  FEEvaluation<dim,fe_degree,fe_degree+1,1,value_type> pressure(data);

  // now: combine face and element stuff
  // cell loop
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
  {
    // get all cell quanitites from VP_CROSS!!
    //{
    // velocity
    velocity.reinit(cell);
    velocity.read_dof_values(src, 0);
    velocity.evaluate (true, false, false);

    // pressure
    pressure.reinit(cell);
    pressure.read_dof_values(src, dim);
    pressure.evaluate(false, true, false);

    // and material coefficients
    const VectorizedArray<value_type> rho = this->densities[cell];
    const VectorizedArray<value_type> rho_inv = 1./this->densities[cell];
    const VectorizedArray<value_type> c_sq = this->speeds[cell]*this->speeds[cell];
    //}

    for (unsigned int q=0; q<velocity.n_q_points; ++q)
    {
      const Tensor<1,dim,VectorizedArray<value_type> > pressure_gradient = pressure.get_gradient(q);
      const Tensor<1,dim,VectorizedArray<value_type> > velocity_value = velocity.get_value(q);

      // add contribution from standard source term (not the terms from Cauchy-Kovalewski)
      //{
      VectorizedArray<value_type> rhs = make_vectorized_array<value_type>(0.0);
      if(this->source_term_no>=0)
      {
        Point<dim,VectorizedArray<value_type> > q_points = velocity.quadrature_point(q);
        for (unsigned int n=0; n<rhs.n_array_elements; ++n)
        {
          double xyz[dim];
          for (unsigned int d=0; d<dim; ++d)
            xyz[d] = q_points[d][n];
          rhs[n] = DRT::Problem::Instance()->Funct(this->source_term_no).Evaluate(0,xyz,this->time+this->time_step)
                   -DRT::Problem::Instance()->Funct(this->source_term_no).Evaluate(0,xyz,this->time); // this is a time integral, so we have to evaluate the function which is the indefinite integral
        }
      }
      pressure.submit_value(-rhs*c_sq,q); // minus because minus in time_integrators.h
      //}

      if(this->adjoint_eval==false)
      {
        velocity.submit_value(rho_inv*pressure_gradient,q);
        pressure.submit_gradient(-rho*c_sq*velocity_value,q);
      }
      else
        dserror("have to think about this");
    }

    velocity.integrate (true, false);
    velocity.distribute_local_to_global (dst, 0);

    pressure.integrate(true, true);
    pressure.distribute_local_to_global (dst,dim);
  }

}

template<int dim, int fe_degree, typename Number>
void WaveEquationOperationAcousticWaveADER<dim, fe_degree,Number>::
local_apply_mass_matrix(const MatrixFree<dim,value_type>                  &data,
                        std::vector<parallel::distributed::Vector<value_type> >        &dst,
                        const std::vector<parallel::distributed::Vector<value_type> >  &src,
                        const std::pair<unsigned int,unsigned int>    &cell_range) const
{
  internal::InverseMassMatrixData<dim,fe_degree,dim+1,value_type>& mass_data = this->mass_matrix_data.get();
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

template<int dim, int fe_degree, typename Number>
void WaveEquationOperationAcousticWaveADER<dim, fe_degree, Number>::
apply(const std::vector<parallel::distributed::Vector<value_type> >  &src,
          std::vector<parallel::distributed::Vector<value_type> >        &dst,
          const double                                                   &cur_time,
          const double                                                   &dt) const
{
  Timer timer;
  this->time = cur_time;
  this->time_step = dt;
  this->dirichlet_boundary_conditions->set_time(this->time);
  this->source_term->set_time(this->time);

  //  calculate postprocessed gradient
  /*for(unsigned d=0; d<flux_memory.size(); ++d)
    flux_memory[d]=0.;
  ACOU::WaveEquationOperation<dim,fe_degree,Number>::compute_post_gradient(src,flux_memory,this->time);*/

  // first ader step
  this->data.cell_loop (&WaveEquationOperationAcousticWaveADER<dim, fe_degree, Number>::local_apply_firstader_domain,
                        this, dst, src);

  // create temp vector holding the dst values as source for boundary actions
  std::vector<parallel::distributed::Vector<value_type> > tempsrc(dim+1);
  for(unsigned int d=0; d<dim+1; ++d)
  {
    tempsrc[d].reinit(dst[d],true);
    tempsrc[d] = dst[d];
    dst[d] = 0.;
  }

  this->data.loop (&WaveEquationOperationAcousticWaveADER<dim, fe_degree, Number>::local_apply_secondader_domain,
                   &WaveEquationOperationAcousticWaveADER<dim, fe_degree, Number>::local_apply_ader_face,
                   &WaveEquationOperationAcousticWaveADER<dim, fe_degree, Number>::local_apply_ader_boundary_face,
                   this, dst, tempsrc);

  // timing
  this->computing_times[0] += timer.wall_time();
  timer.restart();

  // inverse mass matrix
  this->data.cell_loop(&WaveEquationOperationAcousticWaveADER<dim, fe_degree, Number>::local_apply_mass_matrix,
               this, dst, dst);


  // timinig
  this->computing_times[1] += timer.wall_time();
  this->computing_times[2] += 1.;
}


// explicit instantiations

template class WaveEquationOperationAcousticWaveADER<2,1,double>;
template class WaveEquationOperationAcousticWaveADER<2,2,double>;
template class WaveEquationOperationAcousticWaveADER<2,3,double>;
template class WaveEquationOperationAcousticWaveADER<2,4,double>;
template class WaveEquationOperationAcousticWaveADER<2,5,double>;
template class WaveEquationOperationAcousticWaveADER<2,6,double>;
template class WaveEquationOperationAcousticWaveADER<3,1,double>;
template class WaveEquationOperationAcousticWaveADER<3,2,double>;
template class WaveEquationOperationAcousticWaveADER<3,3,double>;
template class WaveEquationOperationAcousticWaveADER<3,4,double>;
template class WaveEquationOperationAcousticWaveADER<3,5,double>;
template class WaveEquationOperationAcousticWaveADER<3,6,double>;
template class WaveEquationOperationAcousticWaveADER<2,1,float>;
template class WaveEquationOperationAcousticWaveADER<2,2,float>;
template class WaveEquationOperationAcousticWaveADER<2,3,float>;
template class WaveEquationOperationAcousticWaveADER<2,4,float>;
template class WaveEquationOperationAcousticWaveADER<2,5,float>;
template class WaveEquationOperationAcousticWaveADER<2,6,float>;
template class WaveEquationOperationAcousticWaveADER<3,1,float>;
template class WaveEquationOperationAcousticWaveADER<3,2,float>;
template class WaveEquationOperationAcousticWaveADER<3,3,float>;
template class WaveEquationOperationAcousticWaveADER<3,4,float>;
template class WaveEquationOperationAcousticWaveADER<3,5,float>;
template class WaveEquationOperationAcousticWaveADER<3,6,float>;
}


#endif // HAVE_DEAL_II
