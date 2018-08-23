/*!----------------------------------------------------------------------
\file acou_expl_worker_ader.cpp
\brief Control routine for acoustic explicit time integration with ADER
\level 3

<pre>
\level 3

\maintainer Luca Berardocco
            berardoccoo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
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
WaveEquationOperationAcousticWaveADER(const std::vector<const DoFHandler<dim> *> &dof_handlers,
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
      Teuchos::RCP<MAT::Material> mat = discret->lColElement(element_index)->Material();

      MAT::AcousticMat* actmat = static_cast<MAT::AcousticMat*>(mat.get());
      this->densities[i][v] = actmat->Density();
      this->speeds[i][v] = actmat->SpeedofSound();
    }
  }

  const Teuchos::ParameterList& acouparams = DRT::Problem::Instance()->AcousticParams();
  spectral_evaluation = DRT::INPUT::IntegralValue<bool>(acouparams,"SPECTRAL_EVALUATION");
  use_ader_post = DRT::INPUT::IntegralValue<bool>(acouparams,"USE_ADER_POST");


  // initialize vector for temporary values
  tempsrc.resize(dim+1);
  this->data.initialize_dof_vector(tempsrc[0]);
  for (unsigned int d=1; d<tempsrc.size(); ++d)
    tempsrc[d] = tempsrc[0];
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
  FEEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> phi_eval(this->data);//this->mass_matrix_data->phi[0];
  FEEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type> phi_spectral(this->data,2);
  //}

  // cell loop
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
    {
      integrate_taylor_cauchykovalewski(cell,phi_eval,phi_spectral,src,this->time_step,0.0,0.0,dst);
      phi_eval.set_dof_values(dst,0);
    }
}

template<int dim, int fe_degree, typename Number>
void WaveEquationOperationAcousticWaveADER<dim,fe_degree,Number>::
integrate_taylor_cauchykovalewski(const unsigned int                                              cell,
                                  FEEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type>       &phi_eval,
                                  FEEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type>       &phi_spectral,
                                  const std::vector<parallel::distributed::Vector<value_type> >  &src,
                                  const value_type                                                t2,
                                  const value_type                                                t1,
                                  const value_type                                                te,
                                  const std::vector<parallel::distributed::Vector<value_type> >  &recongraddiv) const
{
  const unsigned int n_q_points = FEEvaluation<dim,fe_degree,fe_degree+1,dim+1,value_type>::n_q_points;

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

  // init cell
  phi_eval.reinit(cell);
  if(spectral_evaluation) phi_spectral.reinit(cell);

  phi_eval.read_dof_values(src, 0);
  phi_eval.evaluate (true, !use_ader_post, false);
  for (unsigned int q=0; q<n_q_points; ++q)
    {
      // contribution from k=0
      const Tensor<1,dim+1,VectorizedArray<value_type> > v_and_p = phi_eval.get_value(q);
      for (unsigned int d=0; d<dim; ++d)
        vcontrib[q][d] = (t2-t1)*v_and_p[d];
      pcontrib[q] = (t2-t1)*v_and_p[dim];

      if (!use_ader_post)
        {
          Tensor<1,dim+1,VectorizedArray<value_type> > v_and_p_next;
          const Tensor<1,dim+1,Tensor<1,dim,VectorizedArray<value_type> > > v_and_p_grad = phi_eval.get_gradient(q);
          v_and_p_next[dim] = VectorizedArray<value_type>();
          for (unsigned int d=0; d<dim; ++d)
            {
              v_and_p_next[d] = -rho_inv*v_and_p_grad[dim][d];
              v_and_p_next[dim] -= c_sq*rho*v_and_p_grad[d][d];
            }
          // add contribution from k=1
          for (unsigned int d=0; d<dim; ++d)
            vcontrib[q][d] += ((t2-te)*(t2-te)-(t1-te)*(t1-te))*0.5*v_and_p_next[d];
          pcontrib[q] +=  ((t2-te)*(t2-te)-(t1-te)*(t1-te))*0.5*v_and_p_next[dim];

          // submit for further evaluation
          if(spectral_evaluation)
          {
            for(unsigned int d=0; d<dim+1; ++d)
              phi_spectral.begin_dof_values()[q+d*n_q_points] = -v_and_p_next[d];
          }
          else
            phi_eval.submit_value(-v_and_p_next,q);
        }
    }

  if (use_ader_post)
    {
      phi_eval.read_dof_values(recongraddiv, 0);
      phi_eval.evaluate(true, false, false);

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          Tensor<1,dim+1,VectorizedArray<value_type> > v_and_p_next = phi_eval.get_value(q);

          // add contribution from k=1
          for (unsigned int d=0; d<dim; ++d)
            vcontrib[q][d] +=  ((t2-te)*(t2-te)-(t1-te)*(t1-te))*0.5*v_and_p_next[d];
          pcontrib[q] +=  ((t2-te)*(t2-te)-(t1-te)*(t1-te))*0.5*v_and_p_next[dim];

          // submit for further evaluation
          if(spectral_evaluation)
          {
            for(unsigned int d=0; d<dim+1; ++d)
              phi_spectral.begin_dof_values()[q+d*n_q_points] = -v_and_p_next[d];
          }
          else
            phi_eval.submit_value(-v_and_p_next,q);
        }
    }

  this->mass_matrix_data->inverse.fill_inverse_JxW_values(this->mass_matrix_data->coefficients);

  // all following contributions can be looped
  double fac = -0.5;
  VectorizedArray<value_type> fac_t;
  for (int k=2; k<=fe_degree; ++k)
    {
      fac /= -(k+1);
      fac_t = std::pow(t2-te,value_type(k+1))-std::pow(t1-te,value_type(k+1));

      if(spectral_evaluation)
        phi_spectral.evaluate(false,true);
      else
      {
        // integrate over element
        phi_eval.integrate(true,false);

        // apply inverse mass matrix
        //{
        this->mass_matrix_data->inverse.apply(this->mass_matrix_data->coefficients, dim+1,
                                              phi_eval.begin_dof_values(),
                                              phi_eval.begin_dof_values());
        //}

        // evaulate this phi at the gauss points
        phi_eval.evaluate(false,true);
      }

      // sum over all integration points
      for (unsigned int q=0; q<n_q_points; ++q)
        {
          // get the gauss point values
          Tensor<1,dim+1,Tensor<1,dim,VectorizedArray<value_type> > > phi_gradient;
          if(spectral_evaluation)
            phi_gradient = phi_spectral.get_gradient(q);
          else
            phi_gradient = phi_eval.get_gradient(q);

          // calculate contributions
          for (unsigned int d=0; d<dim; ++d)
            vcontrib[q][d] += fac*fac_t*rho_inv*phi_gradient[dim][d];
          for (unsigned int d=0; d<dim; ++d)
            pcontrib[q] += fac*fac_t*c_sq*rho*phi_gradient[d][d];

          // evaluate things phi_k+1 needs
          if(spectral_evaluation)
          {
            for(unsigned int d=0; d<dim; ++d)
              phi_spectral.begin_dof_values()[q+d*n_q_points] = -rho_inv*phi_gradient[dim][d];
            phi_spectral.begin_dof_values()[q+dim*n_q_points] = c_sq*rho*phi_gradient[0][0];
            for(unsigned int d=1; d<dim; ++d)
              phi_spectral.begin_dof_values()[q+dim*n_q_points] += c_sq*rho*phi_gradient[d][d];
          }
          else
          {
            Tensor<1,dim+1,VectorizedArray<value_type> > temp;
            for (unsigned int d=0; d<dim; ++d)
              {
                temp[d] = rho_inv*phi_gradient[dim][d];
                temp[dim] += c_sq*rho*phi_gradient[d][d];
              }
            phi_eval.submit_value(temp,q);
          }
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
  this->mass_matrix_data->inverse.apply(this->mass_matrix_data->coefficients, dim+1,
                                        phi_eval.begin_dof_values(),
                                        phi_eval.begin_dof_values());
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
      this->evaluate_inner_face(phi,phi_neighbor,src,face,-1.0);

      phi.distribute_local_to_global(dst, 0);
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

  for (unsigned int face=face_range.first; face<face_range.second; face++)
    {
      this->evaluate_boundary_face(phi,src,face,-1.0);
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
      // get all cell quanitites
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
        pressure.submit_gradient(-rho*c_sq*velocity.get_value(q), q);
        velocity.submit_value(rho_inv*pressure_gradient, q);
      }

      velocity.integrate (true, false);
      velocity.distribute_local_to_global (dst, 0);

      pressure.integrate(false, true);
      pressure.distribute_local_to_global (dst,dim);
    }

}

template<int dim, int fe_degree, typename Number>
void WaveEquationOperationAcousticWaveADER<dim, fe_degree, Number>::
apply_ader(const std::vector<parallel::distributed::Vector<value_type> >  &src,
           std::vector<parallel::distributed::Vector<value_type> >        &dst,
           const double                                                   &cur_time,
           const double                                                   &dt) const
{
  for (unsigned int d=0; d<=dim; ++d)
    {
      dst[d]=0;
      tempsrc[d]=0;
    }
  Timer timer;
  WaveEquationOperation<dim,fe_degree,Number>::apply(src,tempsrc,cur_time,dt);
  this->computing_times[3] += timer.wall_time();

  // first ader step
  timer.restart();
  this->data.cell_loop (&WaveEquationOperationAcousticWaveADER<dim, fe_degree, Number>::local_apply_firstader_domain,
                        this, tempsrc, src);
  this->computing_times[4] += timer.wall_time();
  timer.restart();
  this->data.loop (&WaveEquationOperationAcousticWaveADER<dim, fe_degree, Number>::local_apply_secondader_domain,
                   &WaveEquationOperationAcousticWaveADER<dim, fe_degree, Number>::local_apply_ader_face,
                   &WaveEquationOperationAcousticWaveADER<dim, fe_degree, Number>::local_apply_ader_boundary_face,
                   this, dst, tempsrc);
  this->computing_times[5] += timer.wall_time();

  // inverse mass matrix
  timer.restart();
  this->data.cell_loop(&WaveEquationOperation<dim, fe_degree, Number>::local_apply_mass_matrix,
                       static_cast<const WaveEquationOperation<dim,fe_degree, Number>*>(this), dst, dst);
  this->computing_times[6] += timer.wall_time();

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
