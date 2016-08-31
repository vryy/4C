/*!----------------------------------------------------------------------
\file acou_expl_worker_sol.cpp
\brief Control routine for acoustic explicit time integration for solids
       or aka elastodynamics
\level 2

<pre>
\maintainer Svenja Schoeder
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
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/acoustic.H"
#include "../drt_mat/acoustic_sol.H"
#include "acou_ele.H"
#include "acou_sol_ele.H"


namespace ACOU
{


template<int dim, int fe_degree, typename Number>
void WaveEquationOperation<dim, fe_degree,Number>::
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



template <int dim, int fe_degree, typename Number>
void
WaveEquationOperation<dim,fe_degree,Number>::
local_apply_solid_face (const MatrixFree<dim,value_type> &,
                        std::vector<parallel::distributed::Vector<value_type> > &dst,
                        const std::vector<parallel::distributed::Vector<value_type> > &src,
                        const std::pair<unsigned int,unsigned int> &face_range) const
{
  // There is some overhead in the methods in FEEvaluation, so it is faster
  // to combine pressure and velocity in the same object and just combine
  // them at the level of quadrature points
  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim*dim+dim+1,value_type> phi(this->data, true, 0, 0, 0, true);
  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim*dim+dim+1,value_type> phi_neighbor(this->data, false, 0, 0, 0, true);

  for (unsigned int face=face_range.first; face<face_range.second; face++)
  {
    phi.reinit(face);
    phi.read_dof_values(src, 0);
    phi.evaluate(true,false);
    const VectorizedArray<value_type> rho_plus = phi.read_cell_data(densities);
    const VectorizedArray<value_type> rho_inv_plus = 1./rho_plus;
    const VectorizedArray<value_type> c_plus = phi.read_cell_data(speeds);
    const VectorizedArray<value_type> c_sq_plus = c_plus * c_plus;
    const VectorizedArray<value_type> tau_plus = make_vectorized_array<value_type>(1.0);
    const VectorizedArray<value_type> visc_plus = phi.read_cell_data(viscs);

    phi_neighbor.reinit(face);
    phi_neighbor.read_dof_values(src, 0);
    phi_neighbor.evaluate(true,false);
    const VectorizedArray<value_type> rho_minus = phi_neighbor.read_cell_data(densities);
    const VectorizedArray<value_type> rho_inv_minus = 1./rho_minus;
    const VectorizedArray<value_type> c_minus = phi_neighbor.read_cell_data(speeds);
    const VectorizedArray<value_type> c_sq_minus = c_minus * c_minus;
    const VectorizedArray<value_type> tau_minus = make_vectorized_array<value_type>(1.0);
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

template <int dim, int fe_degree, typename Number>
void WaveEquationOperation<dim,fe_degree,Number>::
local_apply_solid_boundary_face (const MatrixFree<dim,value_type> &,
                           std::vector<parallel::distributed::Vector<value_type> >       &dst,
                           const std::vector<parallel::distributed::Vector<value_type> > &src,
                           const std::pair<unsigned int,unsigned int>               &face_range) const
{
  FEFaceEvaluation<dim,fe_degree,fe_degree+1,dim*dim+dim+1,value_type> phi(this->data, true, 0, 0, 0, true);

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
    const VectorizedArray<value_type> tau = make_vectorized_array<value_type>(1.0);
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

template<int dim, int fe_degree, typename Number>
void WaveEquationOperation<dim, fe_degree,Number>::
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
}


#endif // HAVE_DEAL_II
