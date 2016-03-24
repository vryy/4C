/*!----------------------------------------------------------------------
\file acou_expl_worker.H
\brief Control routine for acoustic explicit time integration.

<pre>
\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
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

  // forward declaration
  namespace internal
  {
    template <int,int,typename> struct InverseMassMatrixData;
    template <int,int,typename> struct InverseMassMatrixDataSolid;
  }


  template <int dim, typename Number>
  class WaveEquationOperationBase
  {
  public:

    typedef Number value_type;

    typedef parallel::distributed::Vector<value_type> vector_type;

    virtual ~WaveEquationOperationBase() {}

    virtual void apply (const std::vector<parallel::distributed::Vector<value_type> >  &src,
                        std::vector<parallel::distributed::Vector<value_type> >        &dst,
                        const double                                              &cur_time) const = 0;

    virtual void applyader (const std::vector<parallel::distributed::Vector<value_type> >  &src,
                            std::vector<parallel::distributed::Vector<value_type> >        &dst,
                            const double                                              &cur_time,
                            const double                                                    &dt) const = 0;

    virtual void initialize_dof_vector(parallel::distributed::Vector<value_type> &src) const = 0;

    virtual void read_initial_conditions(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                                         std::vector<parallel::distributed::Vector<value_type> > &dst) = 0;

    virtual void write_deal_cell_values(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                                        const std::vector<parallel::distributed::Vector<value_type> >   &src) const = 0;

    virtual void compute_gradient_contributions(std::vector<parallel::distributed::Vector<value_type> > &fwnp,
                                                std::vector<parallel::distributed::Vector<value_type> > &fwn,
                                                std::vector<parallel::distributed::Vector<value_type> > &adnp) = 0;

    virtual void write_gradient_contributions(Teuchos::RCP<DRT::DiscretizationHDG> &discret, value_type dt) const = 0;

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
    bool         solid;
  };



  template<int dim, int fe_degree, typename Number>
  class WaveEquationOperation : public WaveEquationOperationBase<dim,Number>
  {
  public:
    typedef typename WaveEquationOperationBase<dim,Number>::value_type value_type;

    WaveEquationOperation(const DoFHandler<dim> &dof_handler,
                          Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                          Teuchos::RCP<Function<dim> > boundary_conditions,
                          Teuchos::RCP<Function<dim> > source_term,
                          bool sol_in,
                          int sourceno = 0,
                          Teuchos::RCP<Epetra_MultiVector> source_adjoint = Teuchos::null);

    virtual ~WaveEquationOperation();

    virtual void apply (const std::vector<parallel::distributed::Vector<value_type> >  &src,
                        std::vector<parallel::distributed::Vector<value_type> >        &dst,
                        const double                                              &cur_time) const;

    virtual void applyader (const std::vector<parallel::distributed::Vector<value_type> >  &src,
                            std::vector<parallel::distributed::Vector<value_type> >        &dst,
                            const double                                              &cur_time,
                            const double                                                    &dt) const;

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
    }

    AlignedVector<VectorizedArray<value_type> > densities, speeds, densities_grad, speeds_grad, viscs;
    std::vector<unsigned int> dofpermutations;

    Table<2,unsigned int> permutevalues; // each cell has to have the same number of dofs (this does not work for p-adaptivity!)

    MatrixFree<dim,value_type> data;

  private:
    mutable double                           time;
    mutable double                           time_step;
    mutable std::vector<double>              computing_times;
    mutable int                              source_term_no;

    Table<4,value_type >                     table_node_coords;
    Table<3,unsigned int>                    table_node_ids;
    mutable Teuchos::RCP<Epetra_MultiVector> source_adjoint_meas;

    //void apply_mass_matrix(const std::vector<parallel::distributed::Vector<value_type> > &src,
    //                       std::vector<parallel::distributed::Vector<value_type> >       &dst);

    void local_apply_mass_matrix(const MatrixFree<dim,value_type>                             &data,
                                 std::vector<parallel::distributed::Vector<value_type> >       &dst,
                                 const std::vector<parallel::distributed::Vector<value_type> > &src,
                                 const std::pair<unsigned int,unsigned int>               &cell_range) const;

    void local_apply_domain (const MatrixFree<dim,value_type>                              &data,
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

    void local_apply_solid_mass_matrix(const MatrixFree<dim,value_type>                             &data,
                                       std::vector<parallel::distributed::Vector<value_type> >       &dst,
                                       const std::vector<parallel::distributed::Vector<value_type> > &src,
                                       const std::pair<unsigned int,unsigned int>               &cell_range) const;
    void local_apply_solid_domain(const MatrixFree<dim,value_type>                              &data,
                                   std::vector<parallel::distributed::Vector<value_type> >        &dst,
                                   const std::vector<parallel::distributed::Vector<value_type> >  &src,
                                   const std::pair<unsigned int,unsigned int>                &cell_range) const;
    void local_apply_solid_face (const MatrixFree<dim,value_type>       &data,
                                  std::vector<parallel::distributed::Vector<value_type> >         &dst,
                                  const std::vector<parallel::distributed::Vector<value_type> >   &src,
                                  const std::pair<unsigned int,unsigned int>                 &cell_range) const;
    void local_apply_solid_boundary_face (const MatrixFree<dim,value_type>                              &data,
                                           std::vector<parallel::distributed::Vector<value_type> >        &dst,
                                           const std::vector<parallel::distributed::Vector<value_type> >  &src,
                                           const std::pair<unsigned int,unsigned int>              &cell_range) const;

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

    value_type evaluate_source_adjoint(const Point<dim> &p, const std::vector<std::vector<value_type> > nodes, std::vector<value_type> values) const;
    value_type evaluate_source_timereversal(const Point<dim> &p, const std::vector<std::vector<value_type> > nodes, std::vector<value_type> values) const;

    mutable Teuchos::RCP<Function<dim> > dirichlet_boundary_conditions;
    mutable Teuchos::RCP<Function<dim> > source_term;

    mutable Threads::ThreadLocalStorage<internal::InverseMassMatrixData<dim,fe_degree,value_type> > mass_matrix_data;
    mutable Threads::ThreadLocalStorage<internal::InverseMassMatrixDataSolid<dim,fe_degree,value_type> > mass_matrix_data_solid;
  };



  namespace internal
  {
    /// Collect all data for the inverse mass matrix operation in a struct in
    /// order to avoid allocating the memory repeatedly.
    template <int dim, int fe_degree, typename Number>
    struct InverseMassMatrixData
    {
      InverseMassMatrixData(const MatrixFree<dim,Number> &data)
        :
        phi(1, FEEvaluation<dim,fe_degree,fe_degree+1,dim+1,Number>(data)),
        coefficients(FEEvaluation<dim,fe_degree,fe_degree+1,dim+1,Number>::n_q_points),
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
      AlignedVector<FEEvaluation<dim,fe_degree,fe_degree+1,dim+1,Number> > phi;
      AlignedVector<VectorizedArray<Number> > coefficients;
      MatrixFreeOperators::CellwiseInverseMassMatrix<dim,fe_degree,dim+1,Number> inverse;
    };

    /// Collect all data for the inverse mass matrix operation in a struct in
    /// order to avoid allocating the memory repeatedly.
    template <int dim, int fe_degree, typename Number>
    struct InverseMassMatrixDataSolid
    {
      InverseMassMatrixDataSolid(const MatrixFree<dim,Number> &data)
        :
        phi(1, FEEvaluation<dim,fe_degree,fe_degree+1,dim+1+dim*dim,Number>(data)),
        coefficients(FEEvaluation<dim,fe_degree,fe_degree+1,dim+1+dim*dim,Number>::n_q_points),
        inverse(phi[0])
      {}

      // Manually implement the copy operator because CellwiseInverseMassMatrix
      // must point to the object 'phi'
      InverseMassMatrixDataSolid(const InverseMassMatrixDataSolid &other)
        :
        phi(other.phi),
        coefficients(other.coefficients),
        inverse(phi[0])
      {}

      // For memory alignment reasons, need to place the FEEvaluation object
      // into an aligned vector
      AlignedVector<FEEvaluation<dim,fe_degree,fe_degree+1,dim+1+dim*dim,Number> > phi;
      AlignedVector<VectorizedArray<Number> > coefficients;
      MatrixFreeOperators::CellwiseInverseMassMatrix<dim,fe_degree,dim+1+dim*dim,Number> inverse;
    };
  }
}

#endif // HAVE_DEAL_II

#endif
