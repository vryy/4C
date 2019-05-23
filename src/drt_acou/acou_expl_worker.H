/*!----------------------------------------------------------------------
\file acou_expl_worker.H
\brief Control routine for acoustic explicit time integration.

\level 2

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
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
//#include <deal.II/matrix_free/fe_evaluation.h>
#include "fe_evaluation.h"
#include <deal.II/matrix_free/operators.h>
#include <Teuchos_RCP.hpp>
#include "Epetra_MultiVector.h"


namespace DRT
{
  class DiscretizationHDG;
}

namespace ACOU
{
  class PATMonitorManager;


  using namespace dealii;

  // forward declarations
  template <int, typename Number>
  class AttenuationPML;
  namespace internal
  {
    template <int, int, int, typename>
    struct InverseMassMatrixData;
  }
  template <typename Number>
  class ClusterManager;


  template <int dim, typename Number>
  class WaveEquationOperationBase
  {
   public:
    typedef Number value_type;

    typedef parallel::distributed::Vector<value_type> vector_type;

    virtual void apply(const std::vector<parallel::distributed::Vector<value_type>> &src,
        std::vector<parallel::distributed::Vector<value_type>> &dst, const double &cur_time,
        const double &dt = 0.0) const = 0;

    virtual void apply_ader(const std::vector<parallel::distributed::Vector<value_type>> &,
        std::vector<parallel::distributed::Vector<value_type>> &, const double &,
        const double &) const = 0;

    virtual void compute_post_gradient(
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        std::vector<parallel::distributed::Vector<value_type>> &dst, const double &cur_time,
        const double &dt = 0.0) const = 0;

    virtual void initialize_dof_vector(parallel::distributed::Vector<value_type> &src) const = 0;

    virtual void read_initial_conditions(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
        std::vector<parallel::distributed::Vector<value_type>> &dst) = 0;

    virtual void write_deal_cell_values(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
        const std::vector<parallel::distributed::Vector<value_type>> &src) const = 0;

    virtual void compute_gradient_contributions(
        std::vector<parallel::distributed::Vector<value_type>> &fwnp,
        std::vector<parallel::distributed::Vector<value_type>> &fwn,
        std::vector<parallel::distributed::Vector<value_type>> &adnp) = 0;

    virtual double get_SoS_gradient(int rowid) const = 0;
    virtual double get_density_gradient(int rowid) const = 0;

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

    void set_reduction(const bool red)
    {
      reduction = red;
      return;
    }

    void set_time_step_size(double time_step_in) const
    {
      time_step = time_step_in;
      return;
    }

    value_type get_global_time_step_size() const { return time_step; }

    virtual value_type get_time_step_size(unsigned int) const { return time_step; }

    virtual unsigned int cluster_id(unsigned int) const { return 0; }

    const MatrixFree<dim, value_type> &get_matrix_free() const { return data; }
    MatrixFree<dim, value_type> data;

   protected:
    unsigned int timestep_source_number;
    bool adjoint_eval;
    bool reduction;
    mutable value_type time_step;
  };



  template <int dim, int fe_degree, typename Number>
  class WaveEquationOperation : public WaveEquationOperationBase<dim, Number>
  {
   public:
    typedef typename WaveEquationOperationBase<dim, Number>::value_type value_type;
    static const unsigned int n_vectorization = VectorizedArray<value_type>::n_array_elements;

    WaveEquationOperation(const std::vector<const DoFHandler<dim> *> &dof_handlers,
        Teuchos::RCP<DRT::DiscretizationHDG> &discret,
        Teuchos::RCP<Function<dim>> boundary_conditions, Teuchos::RCP<Function<dim>> source_term,
        value_type time_step_in, int sourceno = 0,
        Teuchos::RCP<PATMonitorManager> monitormanagerin = Teuchos::null);

    virtual ~WaveEquationOperation();

    virtual void apply(const std::vector<parallel::distributed::Vector<value_type>> &src,
        std::vector<parallel::distributed::Vector<value_type>> &dst, const double &cur_time,
        const double &dt) const;

    virtual void apply_ader(const std::vector<parallel::distributed::Vector<value_type>> &,
        std::vector<parallel::distributed::Vector<value_type>> &, const double &,
        const double &) const
    {
      AssertThrow(false, ExcNotImplemented());
    }

    virtual void compute_post_gradient(
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        std::vector<parallel::distributed::Vector<value_type>> &dst, const double &cur_time,
        const double &dt = 0.0) const;

    virtual int number_of_solutionvectors() = 0;

    virtual void read_initial_conditions(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
        std::vector<parallel::distributed::Vector<value_type>> &dst);

    virtual void write_deal_cell_values(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
        const std::vector<parallel::distributed::Vector<value_type>> &src) const;

    virtual void compute_gradient_contributions(
        std::vector<parallel::distributed::Vector<value_type>> &fwnp,
        std::vector<parallel::distributed::Vector<value_type>> &fwn,
        std::vector<parallel::distributed::Vector<value_type>> &adnp);

    virtual double get_SoS_gradient(int colid) const;
    virtual double get_density_gradient(int colid) const;

    virtual void initialize_dof_vector(parallel::distributed::Vector<value_type> &src) const
    {
      this->data.initialize_dof_vector(src);
      return;
    }

    AlignedVector<VectorizedArray<value_type>> densities, speeds, densities_grad, speeds_grad;

    std::vector<unsigned int> dofpermutations;
    Table<2, unsigned int> permutevalues;  // each cell has to have the same number of dofs (this
                                           // does not work for p-adaptivity!)

   protected:
    mutable value_type time;
    mutable std::vector<double> computing_times;
    mutable int source_term_no;

    value_type evaluate_source_adjoint(const Point<dim> &p,
        const std::vector<std::vector<value_type>> nodes, std::vector<value_type> values) const;
    value_type evaluate_source_timereversal(const Point<dim> &p,
        const std::vector<std::vector<value_type>> nodes, std::vector<value_type> values) const;

    mutable Teuchos::RCP<Function<dim>> dirichlet_boundary_conditions;
    mutable Teuchos::RCP<Function<dim>> source_term;

    mutable Teuchos::RCP<PATMonitorManager> monitormanager;

    // Mass matrix data
    mutable std_cxx11::shared_ptr<
        internal::InverseMassMatrixData<dim, fe_degree, dim + 1, value_type>>
        mass_matrix_data;
    mutable std_cxx11::shared_ptr<
        internal::InverseMassMatrixData<dim, fe_degree, dim * dim + dim + 1, value_type>>
        mass_matrix_data_solid;
    mutable std_cxx11::shared_ptr<
        internal::InverseMassMatrixData<dim, fe_degree, 2 * dim + 1, value_type>>
        mass_matrix_data_pml;


    virtual void local_apply_domain(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;

    void evaluate_cell(FEEvaluation<dim, fe_degree, fe_degree + 1, dim, value_type> &phi_v,
        FEEvaluation<dim, fe_degree, fe_degree + 1, 1, value_type> &phi_p,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const unsigned int cell) const;

    virtual void local_apply_face(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;

    void evaluate_inner_face(
        FEFaceEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, value_type> &phi,
        FEFaceEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, value_type> &phi_neighbor,
        const std::vector<parallel::distributed::Vector<value_type>> &src, unsigned int face,
        value_type boundary_fac) const;

    virtual void local_apply_boundary_face(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;

    void evaluate_boundary_face(
        FEFaceEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, value_type> &phi,
        const std::vector<parallel::distributed::Vector<value_type>> &src, unsigned int face,
        value_type boundary_fac) const;

    virtual void local_apply_mass_matrix(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;

    template <int, int, typename>
    friend class WaveEquationOperationAcousticWaveADER;
  };

  template <int dim, int fe_degree, typename Number>
  class WaveEquationOperationAcousticWave : public WaveEquationOperation<dim, fe_degree, Number>
  {
   public:
    typedef typename WaveEquationOperationBase<dim, Number>::value_type value_type;

    WaveEquationOperationAcousticWave(const std::vector<const DoFHandler<dim> *> &dof_handlers,
        Teuchos::RCP<DRT::DiscretizationHDG> &discret,
        Teuchos::RCP<Function<dim>> boundary_conditions, Teuchos::RCP<Function<dim>> source_term,
        value_type time_step_in, int sourceno = 0,
        Teuchos::RCP<PATMonitorManager> monitormanagerin = Teuchos::null);

    virtual ~WaveEquationOperationAcousticWave() {}

    virtual void apply(const std::vector<parallel::distributed::Vector<value_type>> &src,
        std::vector<parallel::distributed::Vector<value_type>> &dst, const double &cur_time,
        const double &dt) const;

    virtual int number_of_solutionvectors() { return dim + 1; }

   protected:
    virtual void local_apply_domain(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;
    virtual void local_apply_face(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;
    virtual void local_apply_boundary_face(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;

    virtual void local_apply_mass_matrix(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;
  };

  template <int dim, int fe_degree, typename Number>
  class WaveEquationOperationElasticWave : public WaveEquationOperation<dim, fe_degree, Number>
  {
   public:
    typedef typename WaveEquationOperationBase<dim, Number>::value_type value_type;

    WaveEquationOperationElasticWave(const std::vector<const DoFHandler<dim> *> &dof_handlers,
        Teuchos::RCP<DRT::DiscretizationHDG> &discret,
        Teuchos::RCP<Function<dim>> boundary_conditions, Teuchos::RCP<Function<dim>> source_term,
        value_type time_step_in, int sourceno = 0,
        Teuchos::RCP<PATMonitorManager> monitormanagerin = Teuchos::null);

    virtual ~WaveEquationOperationElasticWave() {}

    virtual void apply(const std::vector<parallel::distributed::Vector<value_type>> &src,
        std::vector<parallel::distributed::Vector<value_type>> &dst, const double &cur_time,
        const double &dt) const;

    virtual int number_of_solutionvectors() { return dim * dim + dim + 1; }

    virtual void read_initial_conditions(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
        std::vector<parallel::distributed::Vector<value_type>> &dst);

    virtual void write_deal_cell_values(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
        const std::vector<parallel::distributed::Vector<value_type>> &src) const;

    AlignedVector<VectorizedArray<value_type>> viscs;

   private:
    void local_apply_mass_matrix(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;
    void local_apply_domain(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;
    void local_apply_face(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;
    void local_apply_boundary_face(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;
  };

  template <int dim, int fe_degree, typename Number>
  class WaveEquationOperationAcousticWavePML
      : public WaveEquationOperationAcousticWave<dim, fe_degree, Number>
  {
   public:
    typedef typename WaveEquationOperationBase<dim, Number>::value_type value_type;

    WaveEquationOperationAcousticWavePML(const std::vector<const DoFHandler<dim> *> &dof_handlers,
        Teuchos::RCP<DRT::DiscretizationHDG> &discret,
        Teuchos::RCP<Function<dim>> boundary_conditions, Teuchos::RCP<Function<dim>> source_term,
        Teuchos::RCP<AttenuationPML<dim, Number>> sigma_fct, value_type time_step_in,
        int sourceno = 0, Teuchos::RCP<PATMonitorManager> monitormanagerin = Teuchos::null,
        bool adjoint = false, bool redin = false);

    virtual ~WaveEquationOperationAcousticWavePML() {}

    virtual void read_initial_conditions(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
        std::vector<parallel::distributed::Vector<value_type>> &dst);

    virtual void write_deal_cell_values(Teuchos::RCP<DRT::DiscretizationHDG> &discret,
        const std::vector<parallel::distributed::Vector<value_type>> &src) const;

    virtual int number_of_solutionvectors() { return 1 + dim + dim; }

   private:
    std::vector<std::vector<std::vector<int>>> layer_reference;
    mutable Teuchos::RCP<AttenuationPML<dim, Number>> sigma_pml;

    virtual void local_apply_domain(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;
    void local_apply_face(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;
    virtual void local_apply_mass_matrix(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;

    std::vector<std::bitset<WaveEquationOperation<dim, fe_degree, Number>::n_vectorization>>
        inner_face_monitored;
  };


  template <int dim, int fe_degree, typename Number>
  class WaveEquationOperationAcousticWaveADER : public WaveEquationOperation<dim, fe_degree, Number>
  {
   public:
    typedef typename WaveEquationOperationBase<dim, Number>::value_type value_type;

    WaveEquationOperationAcousticWaveADER(const std::vector<const DoFHandler<dim> *> &dof_handlers,
        Teuchos::RCP<DRT::DiscretizationHDG> &discret,
        Teuchos::RCP<Function<dim>> boundary_conditions, Teuchos::RCP<Function<dim>> source_term,
        value_type time_step_in, int sourceno = 0,
        Teuchos::RCP<PATMonitorManager> monitormanagerin = Teuchos::null);

    virtual void apply_ader(const std::vector<parallel::distributed::Vector<value_type>> &src,
        std::vector<parallel::distributed::Vector<value_type>> &dst, const double &cur_time,
        const double &dt) const;

    virtual int number_of_solutionvectors() { return dim + 1; }

   protected:
    mutable std::vector<parallel::distributed::Vector<value_type>> tempsrc;

    void local_apply_firstader_domain(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;

    void local_apply_secondader_domain(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;

    void integrate_taylor_cauchykovalewski(const unsigned int cell,
        FEEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, value_type> &phi_eval,
        FEEvaluation<dim, fe_degree, fe_degree + 1, dim + 1, value_type> &phi_spectral,
        const std::vector<parallel::distributed::Vector<value_type>> &src, const value_type t2,
        const value_type t1, const value_type te,
        const std::vector<parallel::distributed::Vector<value_type>> &recongraddiv) const;

    void local_apply_ader_face(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;
    void local_apply_ader_boundary_face(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;

    bool use_ader_post;
    bool spectral_evaluation;
  };

  template <int dim, int fe_degree, typename Number>
  class WaveEquationOperationAcousticWaveADERLTS
      : public WaveEquationOperationAcousticWaveADER<dim, fe_degree, Number>
  {
    template <typename>
    friend class ClusterManager;

   public:
    typedef typename WaveEquationOperationBase<dim, Number>::value_type value_type;
    static const int dimension = dim;

    WaveEquationOperationAcousticWaveADERLTS(
        const std::vector<const DoFHandler<dim> *> &dof_handlers,
        Teuchos::RCP<DRT::DiscretizationHDG> &discret,
        Teuchos::RCP<Function<dim>> boundary_conditions, Teuchos::RCP<Function<dim>> source_term,
        value_type time_step_in, value_type cfl_number, int sourceno = 0,
        Teuchos::RCP<PATMonitorManager> monitormanagerin = Teuchos::null);

    virtual void apply_ader(const std::vector<parallel::distributed::Vector<value_type>> &src,
        std::vector<parallel::distributed::Vector<value_type>> &dst, const double &cur_time,
        const double &dt) const;

    virtual int number_of_solutionvectors() { return dim + 1; }

    unsigned int cluster_id(unsigned int cell) const
    {
      return cluster_manager.cell_cluster_ids[cell];
    }

    virtual value_type get_time_step_size(unsigned int cell) const
    {
      return cluster_manager.get_cell_time_step(cell);
    }

    void communicate_flux_memory() const
    {
      for (unsigned int d = 0; d < flux_memory.size(); ++d)
        flux_memory[d].compress(VectorOperation::add);
      return;
    }

   private:
    // cluster manager
    ClusterManager<value_type> cluster_manager;

    mutable std::vector<parallel::distributed::Vector<value_type>> flux_memory;

    // we need this frequently for index calculations
    static const unsigned int n_vect = VectorizedArray<value_type>::n_array_elements;

    value_type speed_of_sound(int cell_index, int vect_index) const
    {
      return this->speeds[cell_index][vect_index];
    }


    void evaluate_cells_and_faces_first_ader(
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        std::vector<parallel::distributed::Vector<value_type>> &dst) const;

    void evaluate_cells_second_ader(
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        std::vector<parallel::distributed::Vector<value_type>> &dst) const;

    virtual void local_apply_firstader_domain(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;

    virtual void local_apply_secondader_domain(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;

    void local_apply_dummy_domain(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;

    // overwrite face routines
    virtual void local_apply_ader_face(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;

    virtual void local_apply_ader_boundary_face(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;

    virtual void local_apply_postprocessing_domain(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;

    virtual void local_apply_postprocessing_face(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;

    virtual void local_apply_postprocessing_boundary_face(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;

    void local_apply_postprocessing_mass_matrix(const MatrixFree<dim, value_type> &data,
        std::vector<parallel::distributed::Vector<value_type>> &dst,
        const std::vector<parallel::distributed::Vector<value_type>> &src,
        const std::pair<unsigned int, unsigned int> &cell_range) const;


    void reconstruct_div_grad(const std::vector<parallel::distributed::Vector<value_type>> &src,
        std::vector<parallel::distributed::Vector<value_type>> &dst) const;
  };

  template <typename Number>
  class ClusterManager
  {
   public:
    static const unsigned int n_vect = VectorizedArray<Number>::n_array_elements;

    ClusterManager(const double rel_tol = 1.e-4) : relative_tolerance(rel_tol) {}

    template <typename Operator>
    void perform_time_step(const Operator &op,
        const std::vector<parallel::distributed::Vector<Number>> &src,
        std::vector<parallel::distributed::Vector<Number>> &dst) const;

    template <typename Operator>
    void setup(const Operator &op, const unsigned int max_clusters, const unsigned int max_diff,
        Number cfl_number_in);

    template <typename Operator>
    void setup_mf_index_to_cell_index(const Operator &op);

    //{ routines for WaveEquationOperationADERLTS to ask for:
    Number get_cell_time_step(unsigned int cell) const
    {
      return cluster_timestepmultiples[cell_cluster_ids[cell]] * fastest_time_step;
    }

    bool is_evaluate_cell(unsigned int cell) const { return evaluate_cell[cell]; }
    bool is_update_cell(unsigned int cell) const { return update_cell[cell]; }
    bool is_evaluate_face(unsigned int face) const { return evaluate_face[face]; }
    Number get_t1() const { return t1; }
    Number get_t2() const { return t2; }
    Number get_te(unsigned int cell) const { return cell_timelevels[cell]; }
    Number get_dt() const { return dt; }
    std::bitset<n_vect> get_phi_to_dst(unsigned int face) const { return phi_to_dst[face]; }
    std::bitset<n_vect> get_phi_to_fluxmemory(unsigned int face) const
    {
      return phi_to_fluxmemory[face];
    }
    std::bitset<n_vect> get_phi_neighbor_to_dst(unsigned int face) const
    {
      return phi_neighbor_to_dst[face];
    }
    std::bitset<n_vect> get_phi_neighbor_to_fluxmemory(unsigned int face) const
    {
      return phi_neighbor_to_fluxmemory[face];
    }
    //}

    // vector of length cells with correspondent cluster ids
    std::vector<unsigned int> cell_cluster_ids;

    // vectors of length cells with indicator for neighbors with smaller/bigger time step size
    std::vector<bool> cell_have_faster_neighbor;
    std::vector<bool> cell_have_slower_neighbor;

    // number of clusters
    unsigned int n_clusters;

    // fastest time step
    Number fastest_time_step;

    // courant number
    Number cfl_number;

    // is flux considered
    mutable bool is_fluxmemory_considered;

    // improved gradient and divergence
    mutable std::vector<parallel::distributed::Vector<Number>> improvedgraddiv;

    // state vector
    mutable std::vector<parallel::distributed::Vector<Number>> state;

   private:
    // temporary state for reconstruction
    mutable std::vector<parallel::distributed::Vector<Number>> temporary_recon_state;

    // vector of length clusters with correspondent time levels
    mutable std::vector<Number> cluster_timelevels;

    // vector of length cells with correspondent time levels
    mutable std::vector<Number> cell_timelevels;

    // vector of length clusters with correspondent time step sizes
    std::vector<int> cluster_timestepmultiples;

    // update order (of length updates)
    std::vector<unsigned int> cluster_update_order;

    // update times (of length updates)
    mutable std::vector<std::vector<Number>> cluster_update_times;

    // number of updates
    unsigned int n_updates;

    // number of macro cells (row+column)
    unsigned int n_cells_with_ghosts;

    // number of cells
    unsigned int n_cells;

    // vector with flags for cell evaluation
    mutable std::vector<bool> evaluate_cell;

    // vector with flags for cell evaluation
    mutable std::vector<bool> update_cell;

    // vector with flags for face evaluation
    mutable std::vector<bool> evaluate_face;

    // evaluation times
    mutable Number t1, t2, te, dt;

    // helper times
    mutable Number t1fa, t2fa, t1sl, t2sl, t1sa, t2sa;

    // relative tolerance to check for evaluation
    const Number relative_tolerance;

    // mapping between indices of matrix free thing and the rest of the world
    std::vector<int> mf_index;

    // mapping between indices of matrix free and contiguous cell numbers
    std::vector<std::vector<int>> distr_cell_index;

    // vector for communication of cluster ids
    parallel::distributed::Vector<Number> distr_cluster_ids;

    // masks for each face
    mutable std::vector<std::bitset<n_vect>> phi_to_dst;
    mutable std::vector<std::bitset<n_vect>> phi_neighbor_to_dst;
    mutable std::vector<std::bitset<n_vect>> phi_to_fluxmemory;
    mutable std::vector<std::bitset<n_vect>> phi_neighbor_to_fluxmemory;

    template <typename Operator>
    void update_elements(const Operator &op,
        std::vector<parallel::distributed::Vector<Number>> &dst,
        std::vector<parallel::distributed::Vector<Number>> &local_state,
        const unsigned int actual_cluster, const bool write_to_fluxmemory = true) const;

    void communicate_cluster_ids()
    {
      // write to the communication vector
      distr_cluster_ids = 0.;
      for (unsigned int i = 0; i < n_cells_with_ghosts; ++i)
        for (unsigned int v = 0; v < distr_cell_index[i].size(); ++v)
          distr_cluster_ids[distr_cell_index[i][v]] = cell_cluster_ids[i];

      // communicate
      // distr_cluster_ids.compress(VectorOperation::min);
      distr_cluster_ids.compress(
          VectorOperation::insert);  // TODO: update deal in lnm folder to allow for min evaluation
      distr_cluster_ids.update_ghost_values();

      // bring the communicated values back to the cell_cluster_ids
      for (unsigned int i = 0; i < n_cells_with_ghosts; ++i)
      {
        Number minval = std::numeric_limits<Number>::max();
        for (unsigned int v = 0; v < distr_cell_index[i].size(); ++v)
        {
          if (distr_cluster_ids[distr_cell_index[i][v]] < minval)
            minval = distr_cluster_ids[distr_cell_index[i][v]];
        }
        cell_cluster_ids[i] = minval;
      }
    }
  };


  namespace internal
  {
    /// Collect all data for the inverse mass matrix operation in a struct in
    /// order to avoid allocating the memory repeatedly.
    template <int dim, int fe_degree, int component, typename Number>
    struct InverseMassMatrixData
    {
      InverseMassMatrixData(const MatrixFree<dim, Number> &data)
          : phi(1, FEEvaluation<dim, fe_degree, fe_degree + 1, component, Number>(data)),
            coefficients(
                FEEvaluation<dim, fe_degree, fe_degree + 1, component, Number>::n_q_points),
            inverse(phi[0])
      {
      }

      // Manually implement the copy operator because CellwiseInverseMassMatrix
      // must point to the object 'phi'
      InverseMassMatrixData(const InverseMassMatrixData &other)
          : phi(other.phi), coefficients(other.coefficients), inverse(phi[0])
      {
      }

      // For memory alignment reasons, need to place the FEEvaluation object
      // into an aligned vector
      AlignedVector<FEEvaluation<dim, fe_degree, fe_degree + 1, component, Number>> phi;
      AlignedVector<VectorizedArray<Number>> coefficients;
      MatrixFreeOperators::CellwiseInverseMassMatrix<dim, fe_degree, component, Number> inverse;
    };
  }  // namespace internal
}  // namespace ACOU

#endif  // HAVE_DEAL_II

#endif
