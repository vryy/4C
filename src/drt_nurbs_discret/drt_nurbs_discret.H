/*----------------------------------------------------------------------*/
/*!
\file drt_nurbs_discret.H

\brief discretisation with additional knot vectors for nurbs problems
       (isogeometric analysis)

\level 1

\maintainer Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271

*/
/*----------------------------------------------------------------------*/
#ifndef NURBS_DISCRET_H
#define NURBS_DISCRET_H

#include "../drt_lib/drt_utils_discret.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_nurbs_discret/drt_knotvector.H"
#include "../drt_nurbs_discret/drt_control_point.H"

//! forward declarations
namespace LINALG
{
  class Solver;
  class MapExtractor;
  class SparseMatrix;
}

/*!
\brief DRT: namespace of the BACI discretization module

*/
namespace DRT
{

/*!
\brief NURBS: namespace of all isogeometric analysis tools

*/
namespace NURBS
{

  /*!
  \brief A class to manage a nurbs discretization in parallel

         Up to now, it's only a standard discretisation
         extendend by a knotvector.

         The nodes are replaced by control points on
         construction/input; the control points are
         derived from the node class and are hence
         managed by the original discretisation

  \author gammi (gamnitzer@lnm.mw.tum.de)
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
    \brief Destructor

    */
    virtual ~NurbsDiscretization();

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
        return(knots_->Return_n_x_m_x_l(npatch));
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
        return(knots_->ReturnDegree(npatch));
      }

    /*!
    \brief return the offsets

    \return  : The element offsets of all patches

    \author gammi

    */
    virtual std::vector<int> Return_Offsets()
      {
        return(knots_->ReturnOffsets());
      }

    /*!
    \brief return number of elements in each direction

    \param npatch (i)
           the number of the patch we are interested in

    \return  : The number of elements in each direction

    \author gammi

    */
    virtual std::vector<int> Return_nele_x_mele_x_lele(const int npatch)
      {
        return(knots_->Return_nele_x_mele_x_lele(npatch));
      }

   private:
    //! don't want = operator
    NurbsDiscretization operator = (const NurbsDiscretization& old);
    //! don't want copy constructor
    NurbsDiscretization(const DRT::NURBS::NurbsDiscretization& old);

    //! the number of patches in this nurbs discretisation
    int                                    npatches_;

    //! the knotvector
    //! dimension  u     : nurbs curve   (n)
    //! dimensions u,v   : nurbs surface (n x m)
    //! dimensions u,v,w : nurbs volume  (n x m x l)
    Teuchos::RCP<DRT::NURBS::Knotvector>    knots_;

    //! linear solver to strongly enforce non-constant DirichletBC
    //Teuchos::RCP<LINALG::Solver> dbcsolver_;


  }; // class NurbsDiscretization
} // namespace NURBS
namespace UTILS {
  class DbcNurbs : public DRT::UTILS::Dbc
  {
    using Dbc::DoDirichletCondition;
  public:
    /// constructor
    DbcNurbs() {};

    /// destructor
    virtual ~DbcNurbs() {};

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
    virtual void Evaluate(
        const DRT::DiscretizationInterface &   discret,
        const double &                       time,
        const Teuchos::RCP<Epetra_Vector> *  systemvectors,
        Epetra_Vector &                      toggle,
        Teuchos::RCP<std::set<int> > *       dbcgids) const;

    /** \brief Determine Dirichlet condition at given time and apply its
     *         values to a system vector
     *  \param cond            The condition object
     *  \param time            Evaluation time
     *  \param systemvector    Vector to apply DBCs to (eg displ. in structure, vel. in fluids)
     *  \param systemvectord   First time derivative of DBCs
     *  \param systemvectordd  Second time derivative of DBCs
     *  \param toggle          Its i-th compononent is set 1 if it has a DBC, otherwise this component remains untouched
     *  \param dbcgids         Map containing DOFs subjected to Dirichlet boundary
     *                         conditions
     *
     * \author vuong */
    virtual void DoDirichletCondition(
        const DRT::DiscretizationInterface &   discret,
        const DRT::Condition &               cond,
        const double &                       time,
        const Teuchos::RCP<Epetra_Vector> *  systemvectors,
        const Epetra_Vector &                toggle,
        const Teuchos::RCP<std::set<int> > * dbcgids) const;

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
    template<DRT::Element::DiscretizationType distype>
    void FillMatrixAndRHSForLSDirichletBoundary(
        Teuchos::RCP<DRT::Element>              ele,
        const std::vector<Epetra_SerialDenseVector>* knots,
        const std::vector<int>&                 lm,
        const std::vector<int>*                 funct,
        const std::vector<double>*              val,
        const unsigned                          deg,
        const double                            time,
        Epetra_SerialDenseMatrix&               elemass,
        std::vector<Epetra_SerialDenseVector>&  elerhs) const;

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
    template<DRT::Element::DiscretizationType distype>
    void FillMatrixAndRHSForLSDirichletDomain(
        Teuchos::RCP<DRT::Element>              ele,
        const std::vector<Epetra_SerialDenseVector>* knots,
        const std::vector<int>&                 lm,
        const std::vector<int>*                 funct,
        const std::vector<double>*              val,
        const unsigned                          deg,
        const double                            time,
        Epetra_SerialDenseMatrix&               elemass,
        std::vector<Epetra_SerialDenseVector>&  elerhs) const;

  };  // class DbcNurbs
} // namespace UTILS
} // namespace DRT

#endif  // #ifndef NURBS_DISCRET_H
