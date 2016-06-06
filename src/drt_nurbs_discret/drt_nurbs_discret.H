/*----------------------------------------------------------------------*/
/*!
\file drt_nurbs_discret.H

\brief discretisation with additional knot vectors for nurbs problems
       (isogeometric analysis)

\level 1

\maintainer Anh-Tu Vuong

*/
/*----------------------------------------------------------------------*/
#ifndef NURBS_DISCRET_H
#define NURBS_DISCRET_H

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

  /*!
  \brief Evaluate Dirichlet boundary conditions

  Parameters recognized by this method:
  \code
  params.set("total time",acttime); // current total time
  \endcode

  \param params (in): List of parameters
  \param systemvector (out): Vector holding prescribed Dirichlet values
  \param systemvectord (out): Vector holding 1st time derivative of prescribed Dirichlet values
  \param systemvectordd (out): Vector holding 2nd time derivative prescribed Dirichlet values
  \param toggle (out): Vector containing 1.0 for each Dirichlet dof and 0 for everything else
  \param dbcmapextractor (out): Map extractor containing maps for the DOFs subjected to
                                Dirichlet boundary conditions and the remaining/free DOFs
  */
  virtual void EvaluateDirichlet(Teuchos::ParameterList& params,
                                 Teuchos::RCP<Epetra_Vector> systemvector,
                                 Teuchos::RCP<Epetra_Vector> systemvectord,
                                 Teuchos::RCP<Epetra_Vector> systemvectordd,
                                 Teuchos::RCP<Epetra_Vector> toggle,
                                 Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor=Teuchos::null);

protected:
  /*!
  \brief Determine Dirichlet condition at given time and apply its
         values to a system vector

  \param cond            The condition object
  \param usetime
  \param time            Evaluation time
  \param systemvector    Vector to apply DBCs to (eg displ. in structure, vel. in fluids)
  \param systemvectord   First time derivative of DBCs
  \param systemvectordd  Second time derivative of DBCs
  \param toggle          Its i-th compononent is set 1 if it has a DBC, otherwise this component remains untouched
  \param dbcgids         Map containing DOFs subjected to Dirichlet boundary
  conditions

  */
   void DoNurbsLSDirichletCondition(
    DRT::Condition&             cond,
    const bool                  usetime,
    const double                time,
    Teuchos::RCP<Epetra_Vector> systemvector,
    Teuchos::RCP<Epetra_Vector> systemvectord,
    Teuchos::RCP<Epetra_Vector> systemvectordd);

private:
  //! don't want = operator
  NurbsDiscretization operator = (const NurbsDiscretization& old);
  //! don't want copy constructor
  NurbsDiscretization(const DRT::NURBS::NurbsDiscretization& old);

  /*!
  \brief Find the DOFs subjected to Dirichlet boundary and fill toggle vector and dbc set

  \param cond            The condition object
  \param systemvector    Vector to apply DBCs to (only needed for dofmap here)
  \param systemvectord   First time derivative of DBCs (only needed for dofmap here)
  \param systemvectordd  Second time derivative of DBCs (only needed for dofmap here)
  \param toggle          Its i-th compononent is set 1 if it has a DBC, otherwise this component remains untouched
  \param dbcgids         Map containing DOFs subjected to Dirichlet boundary conditions
  \param dbccolgids      Map containing DOFs (also from column nodes) subjected to Dirichlet boundary conditions

  */
  void FindDBCgidAndToggle(
      DRT::Condition&                     cond,
      const Teuchos::RCP<Epetra_Vector>   systemvector,
      const Teuchos::RCP<Epetra_Vector>   systemvectord,
      const Teuchos::RCP<Epetra_Vector>   systemvectordd,
      Teuchos::RCP<Epetra_Vector>         toggle,
      Teuchos::RCP<std::set<int> >        dbcgids,
      Teuchos::RCP<std::set<int> >        dbccolgids=Teuchos::null);

  /*!
  \brief Fill mass matrix and rhs vector for evaluation of least squares dirichlet on a boundary

  \param ele          The element that is to be evaluated
  \param knots        element knot vector
  \param lm           reduced location vector of element (DBC DOFs only)
  \param funct        function information (read from the condition)
  \param curve        curve information (read from the condition)
  \param val          value information (read from the condition)
  \param usetime      bool indicating if curve is to be evaluated at current time
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
      const std::vector<int>*                 curve,
      const std::vector<double>*              val,
      const bool                              usetime,
      const unsigned                          deg,
      const double                            time,
      Epetra_SerialDenseMatrix&               elemass,
      std::vector<Epetra_SerialDenseVector>&  elerhs);

  /*!
  \brief Fill mass matrix and rhs vector for evaluation of least squares dirichlet on a domain

  \param ele          The element that is to be evaluated
    \param knots        element knot vector
  \param lm           reduced location vector of element (DBC DOFs only)
  \param funct        function information (read from the condition)
  \param curve        curve information (read from the condition)
  \param val          value information (read from the condition)
  \param usetime      bool indicating if curve is to be evaluated at current time
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
      const std::vector<int>*                 curve,
      const std::vector<double>*              val,
      const bool                              usetime,
      const unsigned                          deg,
      const double                            time,
      Epetra_SerialDenseMatrix&               elemass,
      std::vector<Epetra_SerialDenseVector>&  elerhs);

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

} // namespace DRT

#endif  // #ifndef NURBS_DISCRET_H
