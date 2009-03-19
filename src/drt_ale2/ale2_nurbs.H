/*!----------------------------------------------------------------------
\file ale2_nurbs.H

A nurbs implementation of the ale2 element

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ALE
#ifdef CCADISCRET
#ifndef ALE2_NURBS_H
#define ALE2_NURBS_H

#include "ale2.H"
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

    

class Ale2Nurbs : public DRT::ELEMENTS::Ale2
{

public:

  /*!
  \brief Standard Constructor

  \param id    : A unique global id
  \param owner : proc id that will own this element
  */
  Ale2Nurbs(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Ale2Nurbs(const Ale2Nurbs& old);

  /*!
  \brief Destructor

  */
  virtual ~Ale2Nurbs();


  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  
  \return my parobject id
  */
  virtual int UniqueParObjectId() const { return ParObject_Ale2_Nurbs;}


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
    if (NumNode()==9 || NumNode()==4)
    {
      return 4;
    }
    else {
      dserror("Could not determine number of lines");
      return -1;
    }
  }

  
private:

};

} // namespace NURBS
} // namespace ELEMENTS
} // namespace DRT


#endif  // #ifndef ALE2_NURBS_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_ALE2
