/*----------------------------------------------------------------------*/
/*! \file

\brief discretisation with additional knot vectors for nurbs problems
       (isogeometric analysis)

\level 1


*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_NURBS_DISCRET_HPP
#define FOUR_C_NURBS_DISCRET_HPP

#include "baci_config.hpp"

#include "baci_lib_discret.hpp"
#include "baci_lib_utils_discret.hpp"
#include "baci_nurbs_discret_control_point.hpp"
#include "baci_nurbs_discret_knotvector.hpp"

BACI_NAMESPACE_OPEN

// forward declarations
namespace CORE::LINALG
{
  class Solver;
  class MapExtractor;
  class SparseMatrix;
}  // namespace CORE::LINALG

namespace DRT
{
  namespace NURBS
  {
    /*!
    \brief A class to manage a nurbs discretization in parallel

    Up to now, it's only a standard discretisation extendend by a knotvector.

    The nodes are replaced by control points on construction/input; the control points are
    derived from the node class and are hence managed by the original discretisation

    */
    class NurbsDiscretization : public DRT::Discretization
    {
     public:
      /*!
      \brief Standard Constructor

      \param name (in): name of this nurbs discretization
      \param comm (in): An epetra comm object associated with this discretization
      */
      NurbsDiscretization(const std::string name, Teuchos::RCP<Epetra_Comm> comm);

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
      virtual void SetKnotVector(Teuchos::RCP<DRT::NURBS::Knotvector> knots);

      /*!
      \brief get a pointer to the knotvector from the discretization

      \return knots : The Knotvector class

      \author gammi

      */
      Teuchos::RCP<DRT::NURBS::Knotvector> GetKnotVector();
      Teuchos::RCP<const DRT::NURBS::Knotvector> GetKnotVector() const;

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
      virtual std::vector<int> Return_nele_x_mele_x_lele(const int npatch)
      {
        return (knots_->Return_nele_x_mele_x_lele(npatch));
      }

     private:
      //! don't want = operator
      NurbsDiscretization operator=(const NurbsDiscretization& old) = delete;

      //! don't want copy constructor
      NurbsDiscretization(const DRT::NURBS::NurbsDiscretization& old) = delete;

      /*!
      \brief The knot vector

      dimension  u     : nurbs curve   (n)
      dimensions u,v   : nurbs surface (n x m)
      dimensions u,v,w : nurbs volume  (n x m x l)
      */
      Teuchos::RCP<DRT::NURBS::Knotvector> knots_;

    };  // class NurbsDiscretization
  }     // namespace NURBS

  namespace UTILS
  {
    class DbcNurbs : public DRT::UTILS::Dbc
    {
      using Dbc::DoDirichletCondition;

     public:
      /// constructor
      DbcNurbs() = default;


     protected:
      /** \brief Evaluate the NURBS DBCs
       *
       *  In the case of NURBs the evaluation is split into 4 steps:
       *
       *  (1) Call the base class function and apply the standard DBCs first
       *  (2) Fill the DBC GIDs row set completely ( NURBS + standard DBCs )
       *  (3) Fill new DBC NURBS GID sets ( row + column ) w/o standard DBC
       *  (4) Now, after reading the row and column information of the NURBS DBCs,
       *      start to build and solve the least squares problem
       *
       *  \author hiermeier, vuong \date 01/17 */
      void Evaluate(const DRT::Discretization& discret, double time,
          const Teuchos::RCP<Epetra_Vector>* systemvectors, DbcInfo& info,
          Teuchos::RCP<std::set<int>>* dbcgids) const override;

      /** \brief Determine Dirichlet condition at given time and apply its
       *         values to a system vector
       *  \param cond            The condition object
       *  \param time            Evaluation time
       *  \param systemvector    Vector to apply DBCs to (eg displ. in structure, vel. in fluids)
       *  \param systemvectord   First time derivative of DBCs
       *  \param systemvectordd  Second time derivative of DBCs
       *  \param toggle          Its i-th compononent is set 1 if it has a DBC, otherwise this
       * component remains untouched \param dbcgids         Map containing DOFs subjected to
       * Dirichlet boundary conditions
       *
       * \author vuong */
      void DoDirichletCondition(const DRT::Discretization& discret, const DRT::Condition& cond,
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
      template <CORE::FE::CellType distype>
      void FillMatrixAndRHSForLSDirichletBoundary(Teuchos::RCP<DRT::Element> ele,
          const std::vector<CORE::LINALG::SerialDenseVector>* knots, const std::vector<int>& lm,
          const std::vector<int>* funct, const std::vector<double>* val, const unsigned deg,
          const double time, CORE::LINALG::SerialDenseMatrix& elemass,
          std::vector<CORE::LINALG::SerialDenseVector>& elerhs) const;

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
      template <CORE::FE::CellType distype>
      void FillMatrixAndRHSForLSDirichletDomain(Teuchos::RCP<DRT::Element> ele,
          const std::vector<CORE::LINALG::SerialDenseVector>* knots, const std::vector<int>& lm,
          const std::vector<int>* funct, const std::vector<double>* val, const unsigned deg,
          const double time, CORE::LINALG::SerialDenseMatrix& elemass,
          std::vector<CORE::LINALG::SerialDenseVector>& elerhs) const;

    };  // class DbcNurbs
  }     // namespace UTILS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif  // NURBS_DISCRET_H
