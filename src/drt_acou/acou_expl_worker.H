/*!----------------------------------------------------------------------
\file acou_expl_worker.H
\brief Control routine for acoustic explicit time integration.

<pre>
\level 2

\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/

#ifndef ACOU_EXPL_WORKER_H
#define ACOU_EXPL_WORKER_H

#ifdef HAVE_DEAL_II

#include <deal.II/lac/parallel_vector.h>
#include <deal.II/base/function.h>
#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/table.h>
#include <deal.II/base/thread_local_storage.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>
//#include "fe_evaluation.h"
#include <deal.II/matrix_free/operators.h>
#include <Teuchos_RCP.hpp>
#include "Epetra_MultiVector.h"


namespace DRT
{
  class DiscretizationHDG;
}

namespace ACOU
{
  using namespace dealii;

  // forward declarations
  template <int,typename Number> class AttenuationPML;
  namespace internal
  {
    template <int,int,int,typename> struct InverseMassMatrixData;
  }


  template <int dim, typename Number>
  class WaveEquationOperationBase
  {
  public:

    typedef Number value_type;

    typedef parallel::distributed::Vector<value_type> vector_type;

    virtual void apply (const std::vector<parallel::distributed::Vector<value_type> >  &src,
                        std::vector<parallel::distributed::Vector<value_type> >        &dst,
                        const double                                              &cur_time,
                        const double                                                    &dt = 0.0) const = 0;

    virtual void compute_post_gradient (const std::vector<parallel::distributed::Vector<value_type> >  &src,
                                        std::vector<parallel::distributed::Vector<value_type> >        &dst,
                                        const double                                              &cur_time,
                                        const double                                                    &dt = 0.0) const = 0;

    virtual void initialize_dof_vector(parallel::distributed::Vector<value_type> &src) const = 0;

    virtual void read_initial_conditions(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                                         std::vector<parallel::distributed::Vector<value_type> > &dst) = 0;

    virtual void write_deal_cell_values(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                                        const std::vector<parallel::distributed::Vector<value_type> >   &src) const = 0;

    virtual void compute_gradient_contributions(std::vector<parallel::distributed::Vector<value_type> > &fwnp,
                                                std::vector<parallel::distributed::Vector<value_type> > &fwn,
                                                std::vector<parallel::distributed::Vector<value_type> > &adnp) = 0;

    virtual void write_gradient_contributions(Teuchos::RCP<DRT::DiscretizationHDG> &discret, value_type dt) const = 0;

    virtual int number_of_solutionvectors() = 0;

    void set_timestep_source_number(const unsigned int index)
    {
      timestep_source_number = index;
      return;
    }

    void set_adjoint_eval(const bool adj)
    {
      adjoint_eval = adj;
      return;
    }

  protected:
    unsigned int timestep_source_number;
    bool         adjoint_eval;
  };



  template<int dim, int fe_degree, typename Number>
  class WaveEquationOperation : public WaveEquationOperationBase<dim,Number>
  {
  public:
    typedef typename WaveEquationOperationBase<dim,Number>::value_type value_type;
    static const unsigned int n_vectorization = VectorizedArray<value_type>::n_array_elements;

    WaveEquationOperation(const DoFHandler<dim> &dof_handler,
                          Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                          Teuchos::RCP<Function<dim> > boundary_conditions,
                          Teuchos::RCP<Function<dim> > source_term,
                          int sourceno = 0,
                          Teuchos::RCP<Epetra_MultiVector> source_adjoint = Teuchos::null);

    virtual ~WaveEquationOperation();

    virtual void apply (const std::vector<parallel::distributed::Vector<value_type> >  &src,
                        std::vector<parallel::distributed::Vector<value_type> >        &dst,
                        const double                                              &cur_time,
                        const double                                                    &dt) const = 0;

    virtual void compute_post_gradient (const std::vector<parallel::distributed::Vector<value_type> >  &src,
                                        std::vector<parallel::distributed::Vector<value_type> >        &dst,
                                        const double                                              &cur_time,
                                        const double                                                    &dt = 0.0) const;

    virtual int number_of_solutionvectors() = 0;

    virtual void read_initial_conditions(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                                         std::vector<parallel::distributed::Vector<value_type> > &dst);

    virtual void write_deal_cell_values(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                                        const std::vector<parallel::distributed::Vector<value_type> >   &src) const;

    virtual void compute_gradient_contributions(std::vector<parallel::distributed::Vector<value_type> > &fwnp,
                                                std::vector<parallel::distributed::Vector<value_type> > &fwn,
                                                std::vector<parallel::distributed::Vector<value_type> > &adnp);

    virtual void write_gradient_contributions(Teuchos::RCP<DRT::DiscretizationHDG> &discret, value_type dt) const;

    virtual void initialize_dof_vector(parallel::distributed::Vector<value_type> &src) const
    {
      data.initialize_dof_vector(src);
      return;
    }

    AlignedVector<VectorizedArray<value_type> > densities, speeds, densities_grad, speeds_grad;

    std::vector<unsigned int> dofpermutations;
    Table<2,unsigned int> permutevalues; // each cell has to have the same number of dofs (this does not work for p-adaptivity!)

    MatrixFree<dim,value_type> data;

  protected:
    mutable double                           time;
    mutable double                           time_step;
    mutable std::vector<double>              computing_times;
    mutable int                              source_term_no;

    Table<4,value_type >                     table_node_coords;
    Table<3,unsigned int>                    table_node_ids;
    mutable Teuchos::RCP<Epetra_MultiVector> source_adjoint_meas;

    value_type evaluate_source_adjoint(const Point<dim> &p, const std::vector<std::vector<value_type> > nodes, std::vector<value_type> values) const;
    value_type evaluate_source_timereversal(const Point<dim> &p, const std::vector<std::vector<value_type> > nodes, std::vector<value_type> values) const;

    mutable Teuchos::RCP<Function<dim> > dirichlet_boundary_conditions;
    mutable Teuchos::RCP<Function<dim> > source_term;

    mutable Threads::ThreadLocalStorage<internal::InverseMassMatrixData<dim,fe_degree,dim+1,value_type> > mass_matrix_data;
    mutable Threads::ThreadLocalStorage<internal::InverseMassMatrixData<dim,fe_degree,dim*dim+dim+1,value_type> > mass_matrix_data_solid;
    mutable Threads::ThreadLocalStorage<internal::InverseMassMatrixData<dim,fe_degree,2*dim+1,value_type> > mass_matrix_data_pml;

    virtual void local_apply_domain (const MatrixFree<dim,value_type>                              &data,
                             std::vector<parallel::distributed::Vector<value_type> >        &dst,
                             const std::vector<parallel::distributed::Vector<value_type> >  &src,
                             const std::pair<unsigned int,unsigned int>                &cell_range) const;
    virtual void local_apply_face (const MatrixFree<dim,value_type>       &data,
                           std::vector<parallel::distributed::Vector<value_type> >         &dst,
                           const std::vector<parallel::distributed::Vector<value_type> >   &src,
                           const std::pair<unsigned int,unsigned int>                 &cell_range) const;
    virtual void local_apply_boundary_face (const MatrixFree<dim,value_type>                              &data,
                                    std::vector<parallel::distributed::Vector<value_type> >        &dst,
                                    const std::vector<parallel::distributed::Vector<value_type> >  &src,
                                    const std::pair<unsigned int,unsigned int>              &cell_range) const;

    virtual void local_apply_mass_matrix(const MatrixFree<dim,value_type>                             &data,
                                 std::vector<parallel::distributed::Vector<value_type> >       &dst,
                                 const std::vector<parallel::distributed::Vector<value_type> > &src,
                                 const std::pair<unsigned int,unsigned int>               &cell_range) const;
  };

  template<int dim, int fe_degree, typename Number>
  class WaveEquationOperationAcousticWave : public WaveEquationOperation<dim,fe_degree,Number>
  {
  public:
    typedef typename WaveEquationOperationBase<dim,Number>::value_type value_type;

    WaveEquationOperationAcousticWave(const DoFHandler<dim> &dof_handler,
                          Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                          Teuchos::RCP<Function<dim> > boundary_conditions,
                          Teuchos::RCP<Function<dim> > source_term,
                          int sourceno = 0,
                          Teuchos::RCP<Epetra_MultiVector> source_adjoint = Teuchos::null);

    virtual ~WaveEquationOperationAcousticWave() {}

    virtual void apply (const std::vector<parallel::distributed::Vector<value_type> >  &src,
                        std::vector<parallel::distributed::Vector<value_type> >        &dst,
                        const double                                              &cur_time,
                        const double                                                    &dt) const;

    virtual int number_of_solutionvectors() { return dim+1; }

  protected:
    virtual void local_apply_domain (const MatrixFree<dim,value_type>                              &data,
                             std::vector<parallel::distributed::Vector<value_type> >        &dst,
                             const std::vector<parallel::distributed::Vector<value_type> >  &src,
                             const std::pair<unsigned int,unsigned int>                &cell_range) const;
    virtual void local_apply_face (const MatrixFree<dim,value_type>       &data,
                           std::vector<parallel::distributed::Vector<value_type> >         &dst,
                           const std::vector<parallel::distributed::Vector<value_type> >   &src,
                           const std::pair<unsigned int,unsigned int>                 &cell_range) const;
    virtual void local_apply_boundary_face (const MatrixFree<dim,value_type>                              &data,
                                    std::vector<parallel::distributed::Vector<value_type> >        &dst,
                                    const std::vector<parallel::distributed::Vector<value_type> >  &src,
                                    const std::pair<unsigned int,unsigned int>              &cell_range) const;

    virtual void local_apply_mass_matrix(const MatrixFree<dim,value_type>                             &data,
                                 std::vector<parallel::distributed::Vector<value_type> >       &dst,
                                 const std::vector<parallel::distributed::Vector<value_type> > &src,
                                 const std::pair<unsigned int,unsigned int>               &cell_range) const;
  };

  template<int dim, int fe_degree, typename Number>
  class WaveEquationOperationElasticWave : public WaveEquationOperation<dim,fe_degree,Number>
  {
  public:
    typedef typename WaveEquationOperationBase<dim,Number>::value_type value_type;

    WaveEquationOperationElasticWave(const DoFHandler<dim> &dof_handler,
                          Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                          Teuchos::RCP<Function<dim> > boundary_conditions,
                          Teuchos::RCP<Function<dim> > source_term,
                          int sourceno = 0,
                          Teuchos::RCP<Epetra_MultiVector> source_adjoint = Teuchos::null);

    virtual ~WaveEquationOperationElasticWave() {}

    virtual void apply (const std::vector<parallel::distributed::Vector<value_type> >  &src,
                        std::vector<parallel::distributed::Vector<value_type> >        &dst,
                        const double                                              &cur_time,
                        const double                                                    &dt) const;

    virtual int number_of_solutionvectors() { return dim*dim+dim+1; }

    virtual void read_initial_conditions(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                                         std::vector<parallel::distributed::Vector<value_type> > &dst);

    virtual void write_deal_cell_values(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                                        const std::vector<parallel::distributed::Vector<value_type> >   &src) const;

    AlignedVector<VectorizedArray<value_type> > viscs;

  private:

    void local_apply_mass_matrix(const MatrixFree<dim,value_type>                             &data,
                                 std::vector<parallel::distributed::Vector<value_type> >       &dst,
                                 const std::vector<parallel::distributed::Vector<value_type> > &src,
                                  const std::pair<unsigned int,unsigned int>               &cell_range) const;
    void local_apply_domain(const MatrixFree<dim,value_type>                              &data,
                                   std::vector<parallel::distributed::Vector<value_type> >        &dst,
                                   const std::vector<parallel::distributed::Vector<value_type> >  &src,
                                   const std::pair<unsigned int,unsigned int>                &cell_range) const;
    void local_apply_face (const MatrixFree<dim,value_type>       &data,
                                  std::vector<parallel::distributed::Vector<value_type> >         &dst,
                                  const std::vector<parallel::distributed::Vector<value_type> >   &src,
                                  const std::pair<unsigned int,unsigned int>                 &cell_range) const;
    void local_apply_boundary_face (const MatrixFree<dim,value_type>                              &data,
                                           std::vector<parallel::distributed::Vector<value_type> >        &dst,
                                           const std::vector<parallel::distributed::Vector<value_type> >  &src,
                                           const std::pair<unsigned int,unsigned int>              &cell_range) const;
  };

  template<int dim, int fe_degree, typename Number>
  class WaveEquationOperationAcousticWavePML : public WaveEquationOperationAcousticWave<dim,fe_degree,Number>
  {
  public:
    typedef typename WaveEquationOperationBase<dim,Number>::value_type value_type;

    WaveEquationOperationAcousticWavePML(const DoFHandler<dim> &dof_handler,
                          Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                          Teuchos::RCP<Function<dim> > boundary_conditions,
                          Teuchos::RCP<Function<dim> > source_term,
                          Teuchos::RCP<AttenuationPML<dim,Number> > sigma_fct,
                          int sourceno = 0,
                          Teuchos::RCP<Epetra_MultiVector> source_adjoint = Teuchos::null);

    virtual ~WaveEquationOperationAcousticWavePML() {}

    virtual void read_initial_conditions(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                                         std::vector<parallel::distributed::Vector<value_type> > &dst);

    virtual void write_deal_cell_values(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                                        const std::vector<parallel::distributed::Vector<value_type> >   &src) const;

    virtual int number_of_solutionvectors() { return 1+dim+dim; }

  private:

    std::vector<std::vector<std::vector<int> > >  layer_reference;
    mutable Teuchos::RCP<AttenuationPML<dim,Number> > sigma_pml;

    virtual void local_apply_domain (const MatrixFree<dim,value_type>                              &data,
                             std::vector<parallel::distributed::Vector<value_type> >        &dst,
                             const std::vector<parallel::distributed::Vector<value_type> >  &src,
                             const std::pair<unsigned int,unsigned int>                &cell_range) const;
    void local_apply_face (const MatrixFree<dim,value_type>       &data,
                                  std::vector<parallel::distributed::Vector<value_type> >         &dst,
                                  const std::vector<parallel::distributed::Vector<value_type> >   &src,
                                  const std::pair<unsigned int,unsigned int>                 &cell_range) const;
    virtual void local_apply_mass_matrix(const MatrixFree<dim,value_type>                             &data,
                                 std::vector<parallel::distributed::Vector<value_type> >       &dst,
                                 const std::vector<parallel::distributed::Vector<value_type> > &src,
                                 const std::pair<unsigned int,unsigned int>               &cell_range) const;

    std::vector<std::bitset<WaveEquationOperation<dim,fe_degree,Number>::n_vectorization> > inner_face_monitored;
  };


  template<int dim, int fe_degree, typename Number>
  class WaveEquationOperationAcousticWaveADER : public WaveEquationOperation<dim,fe_degree,Number>
  {
  public:
    typedef typename WaveEquationOperationBase<dim,Number>::value_type value_type;

    WaveEquationOperationAcousticWaveADER(const DoFHandler<dim> &dof_handler,
                          Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                          Teuchos::RCP<Function<dim> > boundary_conditions,
                          Teuchos::RCP<Function<dim> > source_term,
                          int sourceno = 0,
                          Teuchos::RCP<Epetra_MultiVector> source_adjoint = Teuchos::null);

    virtual ~WaveEquationOperationAcousticWaveADER() {}

    virtual void apply (const std::vector<parallel::distributed::Vector<value_type> >  &src,
                        std::vector<parallel::distributed::Vector<value_type> >        &dst,
                        const double                                              &cur_time,
                        const double                                                    &dt) const;

    virtual int number_of_solutionvectors() { return dim+1; }

  private:
    mutable std::vector<parallel::distributed::Vector<value_type> > flux_memory;
    void local_apply_firstader_domain (const MatrixFree<dim,value_type>                              &data,
                                       std::vector<parallel::distributed::Vector<value_type> >        &dst,
                                       const std::vector<parallel::distributed::Vector<value_type> >  &src,
                                       const std::pair<unsigned int,unsigned int>                &cell_range) const;
    void local_apply_ader_face (const MatrixFree<dim,value_type>       &data,
                                std::vector<parallel::distributed::Vector<value_type> >         &dst,
                                const std::vector<parallel::distributed::Vector<value_type> >   &src,
                                const std::pair<unsigned int,unsigned int>                 &cell_range) const;
    void local_apply_ader_boundary_face (const MatrixFree<dim,value_type>                              &data,
                                         std::vector<parallel::distributed::Vector<value_type> >        &dst,
                                         const std::vector<parallel::distributed::Vector<value_type> >  &src,
                                         const std::pair<unsigned int,unsigned int>              &cell_range) const;
    void local_apply_secondader_domain (const MatrixFree<dim,value_type>                              &data,
                                       std::vector<parallel::distributed::Vector<value_type> >        &dst,
                                       const std::vector<parallel::distributed::Vector<value_type> >  &src,
                                       const std::pair<unsigned int,unsigned int>                &cell_range) const;

    void local_apply_mass_matrix(const MatrixFree<dim,value_type>                             &data,
                                 std::vector<parallel::distributed::Vector<value_type> >       &dst,
                                 const std::vector<parallel::distributed::Vector<value_type> > &src,
                                 const std::pair<unsigned int,unsigned int>               &cell_range) const;
  };


  template<int dim, int fe_degree, typename Number>
  class WaveEquationOperationAcousticWaveADERLTS : public WaveEquationOperation<dim,fe_degree,Number>
  {
  public:
    typedef typename WaveEquationOperationBase<dim,Number>::value_type value_type;

    WaveEquationOperationAcousticWaveADERLTS(const DoFHandler<dim> &dof_handler,
                          Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                          Teuchos::RCP<Function<dim> > boundary_conditions,
                          Teuchos::RCP<Function<dim> > source_term,
                          int sourceno = 0,
                          Teuchos::RCP<Epetra_MultiVector> source_adjoint = Teuchos::null);

    virtual ~WaveEquationOperationAcousticWaveADERLTS() {}

    virtual void apply (const std::vector<parallel::distributed::Vector<value_type> >  &src,
                        std::vector<parallel::distributed::Vector<value_type> >        &dst,
                        const double                                              &cur_time,
                        const double                                                    &dt) const;

    virtual int number_of_solutionvectors() { return dim+1; }

  private:

    AlignedVector<value_type> timesteps;
    mutable parallel::distributed::Vector<value_type> times;
    mutable parallel::distributed::Vector<value_type> oldtimes;
    std::vector<std::vector<int> > neighborcells_i;
    std::vector<std::vector<int> > adjacentfaces_i;
    std::vector<int> timesindices;
    mutable std::vector<parallel::distributed::Vector<value_type> > flux_memory;

    void local_apply_aderlts (std::vector<parallel::distributed::Vector<value_type> > &dst,
                              std::vector<parallel::distributed::Vector<value_type> > &tempsrc,
                              std::vector<parallel::distributed::Vector<value_type> > &src,
                              const value_type                                        eps,
                              const value_type                                        timetoreach) const;
    void evaluatecell(unsigned int                                                   cell,
                      std::vector<parallel::distributed::Vector<value_type> >        &dst,
                      const std::vector<parallel::distributed::Vector<value_type> >  &src,
                      const value_type                                               t1,
                      const value_type                                               t2,
                      bool                                                           fluxcomp) const;
    void evaluateface(int                                                            face,
                      std::vector<parallel::distributed::Vector<value_type> >        &dst,
                      std::vector<parallel::distributed::Vector<value_type> >        &tempsrc,
                      const std::vector<parallel::distributed::Vector<value_type> >  &src,
                      unsigned int                                                   callingcell,
                      value_type                                                     t1,
                      value_type                                                     t2,
                      bool                                                           fluxonly) const;
    void evaluateboundaryface(int                                                            face,
                              std::vector<parallel::distributed::Vector<value_type> >        &dst,
                              std::vector<parallel::distributed::Vector<value_type> >        &tempsrc,
                              const std::vector<parallel::distributed::Vector<value_type> >  &src,
                              int                                                            callingcell,
                              value_type                                                     t1,
                              value_type                                                     t2) const;
    void evaluateupdate (unsigned int                                                   cell,
                         std::vector<parallel::distributed::Vector<value_type> >        &dst,
                         std::vector<parallel::distributed::Vector<value_type> >        &src,
                         std::vector<parallel::distributed::Vector<value_type> >        &tempsrc,
                         value_type                                                     act_dt) const;
    void communicate_flux_memory(std::vector<parallel::distributed::Vector<value_type> >        &dst,
                                 std::vector<parallel::distributed::Vector<value_type> >        &tempsrc,
                                 const std::vector<parallel::distributed::Vector<value_type> >  &src,
                                 const value_type                                               eps) const;
  };


  namespace internal
  {
    /// Collect all data for the inverse mass matrix operation in a struct in
    /// order to avoid allocating the memory repeatedly.
    template <int dim, int fe_degree, int component, typename Number>
    struct InverseMassMatrixData
    {
      InverseMassMatrixData(const MatrixFree<dim,Number> &data)
        :
        phi(1, FEEvaluation<dim,fe_degree,fe_degree+1,component,Number>(data)),
        coefficients(FEEvaluation<dim,fe_degree,fe_degree+1,component,Number>::n_q_points),
        inverse(phi[0])
      {}

      // Manually implement the copy operator because CellwiseInverseMassMatrix
      // must point to the object 'phi'
      InverseMassMatrixData(const InverseMassMatrixData &other)
        :
        phi(other.phi),
        coefficients(other.coefficients),
        inverse(phi[0])
      {}

      // For memory alignment reasons, need to place the FEEvaluation object
      // into an aligned vector
      AlignedVector<FEEvaluation<dim,fe_degree,fe_degree+1,component,Number> > phi;
      AlignedVector<VectorizedArray<Number> > coefficients;
      MatrixFreeOperators::CellwiseInverseMassMatrix<dim,fe_degree,component,Number> inverse;
    };
  }
}

#endif // HAVE_DEAL_II

#endif
