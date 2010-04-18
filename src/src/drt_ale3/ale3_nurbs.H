/*!----------------------------------------------------------------------
\file ale3_nurbs.H

A nurbs implementation of the ale3 element

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ALE
#ifdef CCADISCRET
#ifndef ALE3_NURBS_H
#define ALE3_NURBS_H

#include "ale3.H"
#include "../drt_nurbs_discret/drt_control_point.H"
#include "../drt_fem_general/drt_utils_bspline.H"


/*!
  \brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{
/*!
  \brief ELEMENTS: namespace of elements

*/

namespace ELEMENTS
{

/*!
  \brief NURBS: namespace of all isogeometric analysis tools

*/
namespace NURBS
{



class Ale3Nurbs : public DRT::ELEMENTS::Ale3
{

public:

  /*!
  \brief Standard Constructor

  \param id    : A unique global id
  \param owner : proc id that will own this element
  */
  Ale3Nurbs(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Ale3Nurbs(const Ale3Nurbs& old);

  /*!
  \brief Destructor

  */
  virtual ~Ale3Nurbs();


  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.

  \return my parobject id
  */
  virtual int UniqueParObjectId() const { return ParObject_Ale3_Nurbs;}


  /// Print this element
  void Print(ostream& os) const;


  /*!
  \brief Get shape type of element

  \return nurbs4 or nurbs9

  */
  virtual DiscretizationType Shape() const;


  /*!
  \brief Return number of lines of this element.
  */
  virtual int NumLine() const
  {
    if (NumNode()==27 || NumNode()==8)
    {
      return 12;
    }
    else {
      dserror("Could not determine number of lines");
      return -1;
    }
  }


  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const
  {
    if (NumNode()==27 || NumNode()==8)
    {
      return 6;
    }
    else
    {
      dserror("Could not determine number of surfaces");
      return -1;
    }
  }

private:

};

} // namespace NURBS
} // namespace ELEMENTS
} // namespace DRT


#endif  // #ifndef ALE3_NURBS_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_ALE
