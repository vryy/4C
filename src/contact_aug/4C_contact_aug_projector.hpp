/*----------------------------------------------------------------------------*/
/*! \file
\brief GP projector template

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_AUG_PROJECTOR_HPP
#define FOUR_C_CONTACT_AUG_PROJECTOR_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Mortar
{
  class Element;
}

namespace CONTACT
{
  namespace Aug
  {
    /*--------------------------------------------------------------------------*/
    class ProjectorBase
    {
     protected:
      ProjectorBase(){};

     private:
      static ProjectorBase* get2_d(
          Core::FE::CellType ref_type, Core::FE::CellType tar_type, const bool debug = false);

      template <Core::FE::CellType ref_type>
      static ProjectorBase* get2_d(Core::FE::CellType tar_type, const bool debug = false);

      static ProjectorBase* get3_d(
          Core::FE::CellType ref_type, Core::FE::CellType tar_type, const bool debug = false);

      template <Core::FE::CellType ref_type>
      static ProjectorBase* get3_d(Core::FE::CellType tar_type, const bool debug = false);

     public:
      /// access the singleton pointer of the projector object
      static ProjectorBase* Get(const unsigned probdim, Core::FE::CellType ref_type,
          Core::FE::CellType tar_type, const bool debug = false);

      virtual ~ProjectorBase() = default;

      /** \brief project a point defined on the reference element onto a target
       *  element
       *
       *  \note The auxiliary distance factor is in general NOT the real distance
       *  between the slave and master element, since a non-unit normal vector
       *  is allowed for the projection algorithm. The real distance value is e.g.
       *  given by \f$ d = \| n^{[ref]}(ref_xi) \| alpha \f$ in the case of a
       *  normal defined on the reference element.
       *
       *  \param[in]  ref_ele    reference element
       *  \param[in]  ref_xi     reference parameter coordinates
       *  \param[in]  target_ele target element
       *  \param[out] target_xi  parameter coordinates of the projected point
       *  \param[out] alpha      auxiliary distance factor (see note)
       *
       *  \return TRUE if the local Newton scheme did converge.
       *
       *  \author hiermeier \date 08/17 */
      virtual bool operator()(Mortar::Element& ref_ele, const double* ref_xi,
          Mortar::Element& target_ele, double* target_xi, double& alpha) = 0;

      /// return the relative solution tolerance, i.e. the maximal deviation of
      /// the calculated solution point to the analytical solution
      virtual double get_relative_solution_tolerance() const = 0;

    };  // class ProjectorBase

    /*--------------------------------------------------------------------------*/
    template <class DebugPolicy, unsigned probdim, Core::FE::CellType ref_type,
        Core::FE::CellType tar_type>
    class Projector : public ProjectorBase, public DebugPolicy
    {
      static constexpr unsigned REF_DIM = Core::FE::dim<ref_type>;
      static constexpr unsigned REF_NUMNODES = Core::FE::num_nodes<ref_type>;

      static constexpr unsigned TAR_DIM = Core::FE::dim<tar_type>;
      static constexpr unsigned TAR_NUMNODES = Core::FE::num_nodes<tar_type>;

     public:
      static ProjectorBase* Instance();

     protected:
      /// derived
      bool operator()(Mortar::Element& ref_ele, const double* ref_xi, Mortar::Element& target_ele,
          double* target_xi, double& alpha) override;

      /// derived
      inline double get_relative_solution_tolerance() const override { return rel_sol_tolerance_; }

     private:
      /// constructor
      Projector() : iter_(0), rel_sol_tolerance_(0.0){/* empty */};

      void setup();

      /** \brief Get the jacobian for the GP projection
       *
       *      lmat = [ tarX_{,xi^{1}}, tarX_{,xi^{2}}, -normal(x_ref) ]
       *
       *  \param[out] lmat       : jacobian matrix
       *  \param[out] tar_deriv1 : first derivatives of the nodal shape functions
       *                           at the current target parametric coordinates
       *  \param[in] tar_ele     : reference to the target element
       *  \param[in] tar_coords  : nodal global coordinates of the target element
       *  \param[in] tar_xi      : current parametric coordinates of the target point
       *  \param[in] n_ref       : normal evaluated at the reference point
       *
       *  \author hiermeier \date 07/17 */
      void l_mat_gp(Core::LinAlg::Matrix<probdim, probdim>& lmat,
          Core::LinAlg::Matrix<TAR_DIM, TAR_NUMNODES>& tar_deriv1, Mortar::Element& tar_ele,
          const Core::LinAlg::Matrix<probdim, TAR_NUMNODES>& tar_coords, const double* tar_xi,
          const Core::LinAlg::Matrix<probdim, 1>& n_ref) const;

      /** \brief Get the right-hand-side for the GP projection
       *
       *        rhs = x_tar - x_ref - alpha * n_ref
       *
       *  \param rhs       (out) : calculated right hand side value
       *  \param x_ref      (in) : global position of the reference point
       *  \param n_ref      (in) : normal evaluated at the reference point
       *  \param tar_coords (in) : nodal global coordinates of the target element
       *  \param tar_xi     (in) : current parametric coordinates of the target point
       *  \param alpha      (in) : current distance factor
       *
       *  \return FALSE, if get_global_position failed. Otherwise TRUE.
       *
       *  \author  hiermeier \date 07/17 */
      bool rhs_gp(Core::LinAlg::Matrix<probdim, 1>& rhs,
          const Core::LinAlg::Matrix<probdim, 1>& x_ref,
          const Core::LinAlg::Matrix<probdim, 1>& n_ref, Mortar::Element& target_ele,
          const Core::LinAlg::Matrix<probdim, TAR_NUMNODES>& tar_coords, const double* tar_xi,
          const double& alpha) const;

      /** \brief Get the global position at the parametric coordinates xi
       *
       *  \param[in]  ele    : reference to the considered element
       *  \param[in]  coords : global nodal coordinates
       *  \param[in]  xi     : local parametric coordinates
       *  \param[out] pos    : calculated global position
       *
       *  \return FALSE, if the shape function evaluation failed. Otherwise TRUE.
       *
       *  \author  hiermeier \date 07/17 */
      template <Core::FE::CellType type, unsigned numnodes = Core::FE::num_nodes<type>>
      bool get_global_position(Mortar::Element& ele,
          const Core::LinAlg::Matrix<probdim, numnodes>& coords, const double* xi,
          Core::LinAlg::Matrix<probdim, 1>& pos) const;

     private:
      Core::LinAlg::Matrix<REF_NUMNODES, 1> ref_val_;
      Core::LinAlg::Matrix<probdim, 1> x_ref_;
      Core::LinAlg::Matrix<probdim, 1> n_ref_;

      Core::LinAlg::Matrix<probdim, 1> rhs_;
      Core::LinAlg::Matrix<probdim, probdim> lmat_;
      Core::LinAlg::Matrix<probdim, 1> dx_;

      Core::LinAlg::Matrix<probdim, TAR_NUMNODES> tar_coords_;
      Core::LinAlg::Matrix<TAR_DIM, TAR_NUMNODES> tar_deriv1_;

      unsigned iter_;

      // relative solution tolerance
      double rel_sol_tolerance_;
    };  // class Projector

    /*--------------------------------------------------------------------------*/
    /// empty debugger base class of the projector
    class EmptyProjDebugger
    {
     public:
      inline void write_vector(std::ostream& os, ...) const {};
      inline void writeMatrix(std::ostream& os, ...) const {};
    };

    /*--------------------------------------------------------------------------*/
    /// concrete debugger base class of the projector
    class ProjDebugger
    {
     public:
      inline void write_vector(
          std::ostream& os, unsigned dim, double* vals, const std::string& msg) const
      {
        os << msg << " (vector):\n";
        for (unsigned i = 0; i < dim; ++i)
        {
          os << "(#" << i << "): " << vals[i];
          if (i + 1 < dim) os << ", ";
        }
        os << "\n";
      }

      inline void writeMatrix(std::ostream& os, unsigned rows, unsigned cols, double* vals,
          const std::string& msg) const
      {
        os << msg << " (matrix):\n";
        for (unsigned int i = 0; i < rows; ++i)
        {
          os << "(r#" << i << "): ";
          for (unsigned int j = 0; j < cols; ++j)
          {
            os << vals[i + rows * j];
            if (j + 1 < cols) os << ", ";
          }
          if (i + 1 < rows)
            os << ",\n";
          else
            os << "\n";
        }
      }
    };
  }  // namespace Aug
}  // namespace CONTACT



FOUR_C_NAMESPACE_CLOSE

#endif
