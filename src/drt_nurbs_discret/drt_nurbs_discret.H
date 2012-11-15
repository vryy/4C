/*!----------------------------------------------------------------------
\file drt_nurbs_discret.H

\brief discretisation with additional knot vectors for nurbs problems
       (isogeometric analysis)

<pre>
Maintainer: Peter Gamnitzer
            gammi@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef NURBS_DISCRET_H
#define NURBS_DISCRET_H

#include "../drt_lib/drt_discret.H"
#include "../drt_nurbs_discret/drt_knotvector.H"
#include "../drt_nurbs_discret/drt_control_point.H"

/*!
\brief DRT: namespace of the ccarat discretization module

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
  NurbsDiscretization(const string name, RCP<Epetra_Comm> comm);

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
  virtual void SetKnotVector(RCP<DRT::NURBS::Knotvector> knots);

  /*!
  \brief get a pointer to the knotvector from the discretization

  \return knots : The Knotvector class

  \author gammi

  */
  RCP<DRT::NURBS::Knotvector> GetKnotVector();

  /*!
  \brief return number of knots in each direction

  \param npatch (i)
         the number of the patch we are interested in

  \return  : The number of knots in each direction

  \author gammi

  */
  virtual vector<int> Return_n_x_m_x_l(const int npatch)
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
  virtual vector<int> Return_degree(const int npatch)
    {
      return(knots_->ReturnDegree(npatch));
    }

  /*!
  \brief return the offsets

  \return  : The element offsets of all patches

  \author gammi

  */
  virtual vector<int> Return_Offsets()
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
  virtual vector<int> Return_nele_x_mele_x_lele(const int npatch)
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
  RCP<DRT::NURBS::Knotvector>    knots_;

}; // class NurbsDiscretization

} // namespace NURBS

} // namespace DRT

#endif  // #ifndef NURBS_DISCRET_H
