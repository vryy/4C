/*----------------------------------------------------------------------------*/
/*! \file
\brief GP projector template

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef BACI_CONTACT_AUG_PROJECTOR_H
#define BACI_CONTACT_AUG_PROJECTOR_H

#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"
#include "baci_lib_element.H"

BACI_NAMESPACE_OPEN


namespace MORTAR
{
  class MortarElement;
}

namespace CONTACT
{
  namespace AUG
  {
    /*--------------------------------------------------------------------------*/
    class ProjectorBase
    {
     protected:
      ProjectorBase(){};

     private:
      static ProjectorBase* Get2D(
          CORE::FE::CellType ref_type, CORE::FE::CellType tar_type, const bool debug = false);

      template <CORE::FE::CellType ref_type>
      static ProjectorBase* Get2D(CORE::FE::CellType tar_type, const bool debug = false);

      static ProjectorBase* Get3D(
          CORE::FE::CellType ref_type, CORE::FE::CellType tar_type, const bool debug = false);

      template <CORE::FE::CellType ref_type>
      static ProjectorBase* Get3D(CORE::FE::CellType tar_type, const bool debug = false);

     public:
      /// access the singleton pointer of the projector object
      static ProjectorBase* Get(const unsigned probdim, CORE::FE::CellType ref_type,
          CORE::FE::CellType tar_type, const bool debug = false);

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
      virtual bool operator()(MORTAR::MortarElement& ref_ele, const double* ref_xi,
          MORTAR::MortarElement& target_ele, double* target_xi, double& alpha) = 0;

      /// return the relative solution tolerance, i.e. the maximal deviation of
      /// the calculated solution point to the analytical solution
      virtual double getRelativeSolutionTolerance() const = 0;

    };  // class ProjectorBase

    /*--------------------------------------------------------------------------*/
    template <class DebugPolicy, unsigned probdim, CORE::FE::CellType ref_type,
        CORE::FE::CellType tar_type>
    class Projector : public ProjectorBase, public DebugPolicy
    {
      static const unsigned REF_DIM = CORE::FE::dim<ref_type>;
      static const unsigned REF_NUMNODES = CORE::FE::num_nodes<ref_type>;

      static const unsigned TAR_DIM = CORE::FE::dim<tar_type>;
      static const unsigned TAR_NUMNODES = CORE::FE::num_nodes<tar_type>;

     public:
      static ProjectorBase* Instance();

     protected:
      /// derived
      bool operator()(MORTAR::MortarElement& ref_ele, const double* ref_xi,
          MORTAR::MortarElement& target_ele, double* target_xi, double& alpha) override;

      /// derived
      inline double getRelativeSolutionTolerance() const override { return rel_sol_tolerance_; }

     private:
      /// constructor
      Projector() : iter_(0), rel_sol_tolerance_(0.0){/* empty */};

      void Setup();

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
      void LMatGP(CORE::LINALG::Matrix<probdim, probdim>& lmat,
          CORE::LINALG::Matrix<TAR_DIM, TAR_NUMNODES>& tar_deriv1, MORTAR::MortarElement& tar_ele,
          const CORE::LINALG::Matrix<probdim, TAR_NUMNODES>& tar_coords, const double* tar_xi,
          const CORE::LINALG::Matrix<probdim, 1>& n_ref) const;

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
       *  \return FALSE, if GetGlobalPosition failed. Otherwise TRUE.
       *
       *  \author  hiermeier \date 07/17 */
      bool RhsGP(CORE::LINALG::Matrix<probdim, 1>& rhs,
          const CORE::LINALG::Matrix<probdim, 1>& x_ref,
          const CORE::LINALG::Matrix<probdim, 1>& n_ref, MORTAR::MortarElement& target_ele,
          const CORE::LINALG::Matrix<probdim, TAR_NUMNODES>& tar_coords, const double* tar_xi,
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
      template <CORE::FE::CellType type, unsigned numnodes = CORE::FE::num_nodes<type>>
      bool GetGlobalPosition(MORTAR::MortarElement& ele,
          const CORE::LINALG::Matrix<probdim, numnodes>& coords, const double* xi,
          CORE::LINALG::Matrix<probdim, 1>& pos) const;

     private:
      CORE::LINALG::Matrix<REF_NUMNODES, 1> ref_val_;
      CORE::LINALG::Matrix<probdim, 1> x_ref_;
      CORE::LINALG::Matrix<probdim, 1> n_ref_;

      CORE::LINALG::Matrix<probdim, 1> rhs_;
      CORE::LINALG::Matrix<probdim, probdim> lmat_;
      CORE::LINALG::Matrix<probdim, 1> dx_;

      CORE::LINALG::Matrix<probdim, TAR_NUMNODES> tar_coords_;
      CORE::LINALG::Matrix<TAR_DIM, TAR_NUMNODES> tar_deriv1_;

      unsigned iter_;

      // relative solution tolerance
      double rel_sol_tolerance_;
    };  // class Projector

    /*--------------------------------------------------------------------------*/
    /// empty debugger base class of the projector
    class EmptyProjDebugger
    {
     public:
      inline void writeVector(std::ostream& os, ...) const {};
      inline void writeMatrix(std::ostream& os, ...) const {};
    };

    /*--------------------------------------------------------------------------*/
    /// concrete debugger base class of the projector
    class ProjDebugger
    {
     public:
      inline void writeVector(
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
  }  // namespace AUG
}  // namespace CONTACT



BACI_NAMESPACE_CLOSE

#endif  // CONTACT_AUG_PROJECTOR_H
