/*!----------------------------------------------------------------------
\file acou_expl_worker.H
\brief Control routine for acoustic explicit time integration.

<pre>
Maintainer: Svenja Schoeder
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
#include <deal.II/matrix_free/matrix_free.h>
#include <Teuchos_RCP.hpp>


namespace DRT
{
  class DiscretizationHDG;
}

namespace ACOU
{
  using namespace dealii;

  template <int dim>
  class WaveEquationOperationBase
  {
  public:
    typedef double value_type;

    virtual ~WaveEquationOperationBase() {}

    virtual void apply (const std::vector<parallel::distributed::Vector<value_type> >  &src,
                        std::vector<parallel::distributed::Vector<value_type> >        &dst,
                        const double                                              &cur_time) const = 0;

    virtual void initialize_dof_vector(parallel::distributed::Vector<value_type> &src) const = 0;

    virtual void read_initial_conditions(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                                         std::vector<parallel::distributed::Vector<value_type> > &dst) const = 0;

    virtual void write_deal_cell_values(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                                        const std::vector<parallel::distributed::Vector<value_type> >   &src) const = 0;
  };



  template<int dim, int fe_degree>
  class WaveEquationOperation : public WaveEquationOperationBase<dim>
  {
  public:
    typedef typename WaveEquationOperationBase<dim>::value_type value_type;

    WaveEquationOperation(const DoFHandler<dim> &dof_handler,
                          Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                          Teuchos::RCP<Function<dim> > boundary_conditions,
                          Teuchos::RCP<Function<dim> > source_term);

    virtual ~WaveEquationOperation();

    virtual void apply (const std::vector<parallel::distributed::Vector<value_type> >  &src,
                        std::vector<parallel::distributed::Vector<value_type> >        &dst,
                        const double                                              &cur_time) const;

    virtual void read_initial_conditions(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                                         std::vector<parallel::distributed::Vector<value_type> > &dst) const;

    virtual void write_deal_cell_values(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
                                        const std::vector<parallel::distributed::Vector<value_type> >   &src) const;

    virtual void initialize_dof_vector(parallel::distributed::Vector<value_type> &src) const
    {
      data.initialize_dof_vector(src);
    }

    AlignedVector<VectorizedArray<value_type> > densities, speeds;

    MatrixFree<dim,value_type> data;
  private:
    mutable double             time;


    mutable std::vector<double>                 computing_times;

    void apply_mass_matrix(const std::vector<parallel::distributed::Vector<value_type> > &src,
                           std::vector<parallel::distributed::Vector<value_type> >       &dst);

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

    mutable Teuchos::RCP<Function<dim> > dirichlet_boundary_conditions;
    mutable Teuchos::RCP<Function<dim> > source_term;
  };
}

#endif // HAVE_DEAL_II

#endif
