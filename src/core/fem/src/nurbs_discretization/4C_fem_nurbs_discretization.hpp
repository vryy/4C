/*----------------------------------------------------------------------*/
/*! \file

\brief discretisation with additional knot vectors for nurbs problems
       (isogeometric analysis)

\level 1


*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_FEM_NURBS_DISCRETIZATION_HPP
#define FOUR_C_FEM_NURBS_DISCRETIZATION_HPP

#include "4C_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_discretization_utils.hpp"
#include "4C_fem_nurbs_discretization_control_point.hpp"
#include "4C_fem_nurbs_discretization_knotvector.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::LinAlg
{
  class Solver;
  class MapExtractor;
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace Discret
{
  namespace Nurbs
  {
    /*!
    \brief A class to manage a nurbs discretization in parallel

    Up to now, it's only a standard discretisation extendend by a knotvector.

    The nodes are replaced by control points on construction/input; the control points are
    derived from the node class and are hence managed by the original discretisation

    */
    class NurbsDiscretization : public Core::FE::Discretization
    {
     public:
      /*!
      \brief Standard Constructor

      \param name: name of this nurbs discretization
      \param comm: An epetra comm object associated with this discretization
      \param n_dim: number of space dimensions of this discretization
      */
      NurbsDiscretization(
          const std::string name, Teuchos::RCP<Epetra_Comm> comm, unsigned int n_dim);

      /*!
      \brief Set a knot vector

      Store a the knot vector in the discretization.
      It can then be accessed with the GetKnotVector method.

      \note Knot vectors attached to the discretization have to be
            completely redundant meaning that they are the same on
            each processor. I think its affordable.

      \param knots : The Knotvector class

      \author gammi

      */
      virtual void SetKnotVector(Teuchos::RCP<Discret::Nurbs::Knotvector> knots);

      /*!
      \brief get a pointer to the knotvector from the discretization

      \return knots : The Knotvector class

      \author gammi

      */
      Teuchos::RCP<Discret::Nurbs::Knotvector> GetKnotVector();
      Teuchos::RCP<const Discret::Nurbs::Knotvector> GetKnotVector() const;

      /*!
      \brief return number of knots in each direction

      \param npatch (i)
             the number of the patch we are interested in

      \return  : The number of knots in each direction

      \author gammi

      */
      virtual std::vector<int> Return_n_x_m_x_l(const int npatch)
      {
        return (knots_->Return_n_x_m_x_l(npatch));
      }

      /*!
      \brief return degree in each direction

      \param npatch (i)
             the number of the patch we are interested in

      \return  : The degree in each direction

      \author gammi

      */
      virtual std::vector<int> Return_degree(const int npatch)
      {
        return (knots_->ReturnDegree(npatch));
      }

      /*!
      \brief return the offsets

      \return  : The element offsets of all patches

      \author gammi

      */
      virtual std::vector<int> Return_Offsets() { return (knots_->ReturnOffsets()); }

      /*!
      \brief return number of elements in each direction

      \param npatch (i)
             the number of the patch we are interested in

      \return  : The number of elements in each direction

      \author gammi

      */
      virtual std::vector<int> return_nele_x_mele_x_lele(const int npatch)
      {
        return (knots_->return_nele_x_mele_x_lele(npatch));
      }

     private:
      //! don't want = operator
      NurbsDiscretization operator=(const NurbsDiscretization& old) = delete;

      //! don't want copy constructor
      NurbsDiscretization(const Discret::Nurbs::NurbsDiscretization& old) = delete;

      /*!
      \brief The knot vector

      dimension  u     : nurbs curve   (n)
      dimensions u,v   : nurbs surface (n x m)
      dimensions u,v,w : nurbs volume  (n x m x l)
      */
      Teuchos::RCP<Discret::Nurbs::Knotvector> knots_;

    };  // class NurbsDiscretization
  }     // namespace Nurbs

  namespace UTILS
  {
    class DbcNurbs : public Core::FE::UTILS::Dbc
    {
      using Dbc::do_dirichlet_condition;

     public:
      /// constructor
      DbcNurbs() = default;


     protected:
      void evaluate(const Teuchos::ParameterList& params, const Core::FE::Discretization& discret,
          double time, const Teuchos::RCP<Epetra_Vector>* systemvectors, DbcInfo& info,
          Teuchos::RCP<std::set<int>>* dbcgids) const override;

      void do_dirichlet_condition(const Teuchos::ParameterList& params,
          const Core::FE::Discretization& discret, const Core::Conditions::Condition& cond,
          double time, const Teuchos::RCP<Epetra_Vector>* systemvectors,
          const Epetra_IntVector& toggle,
          const Teuchos::RCP<std::set<int>>* dbcgids) const override;

     private:
      /*!
      \brief Fill mass matrix and rhs vector for evaluation of least squares dirichlet on a boundary

      \param ele          The element that is to be evaluated
      \param knots        element knot vector
      \param lm           reduced location vector of element (DBC DOFs only)
      \param funct        function information (read from the condition)
      \param val          value information (read from the condition)
      \param deg          degree of time derivative needed
      \param time         current time
      \param elemass      element matrix to be filled
      \param elerhs       element right hand side to be filled

      */
      template <Core::FE::CellType distype>

      void fill_matrix_and_rhs_for_ls_dirichlet_boundary(Teuchos::RCP<Core::Elements::Element> ele,
          const std::vector<Core::LinAlg::SerialDenseVector>* knots, const std::vector<int>& lm,
          const std::vector<int>* funct, const std::vector<double>* val, const unsigned deg,
          const double time, Core::LinAlg::SerialDenseMatrix& elemass,
          std::vector<Core::LinAlg::SerialDenseVector>& elerhs,
          const Core::UTILS::FunctionManager& function_manager) const;

      /*!
      \brief Fill mass matrix and rhs vector for evaluation of least squares dirichlet on a domain

      \param ele          The element that is to be evaluated
        \param knots        element knot vector
      \param lm           reduced location vector of element (DBC DOFs only)
      \param funct        function information (read from the condition)
      \param val          value information (read from the condition)
      \param deg          degree of time derivative needed
      \param time         current time
      \param elemass      element matrix to be filled
      \param elerhs       element right hand side to be filled

      */
      template <Core::FE::CellType distype>
      void fill_matrix_and_rhs_for_ls_dirichlet_domain(Teuchos::RCP<Core::Elements::Element> ele,
          const std::vector<Core::LinAlg::SerialDenseVector>* knots, const std::vector<int>& lm,
          const std::vector<int>* funct, const std::vector<double>* val, const unsigned deg,
          const double time, Core::LinAlg::SerialDenseMatrix& elemass,
          std::vector<Core::LinAlg::SerialDenseVector>& elerhs,
          const Core::UTILS::FunctionManager& function_manager) const;

    };  // class DbcNurbs
  }     // namespace UTILS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
