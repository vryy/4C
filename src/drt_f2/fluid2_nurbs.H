/*!----------------------------------------------------------------------
\file fluid2_nurbs.H

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET
#ifndef FLUID2_NURBS_H
#define FLUID2_NURBS_H

#include "fluid2.H"
#include "../drt_lib/drt_control_point.H"
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

    

class Fluid2Nurbs : public DRT::ELEMENTS::Fluid2
{

public:

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  Fluid2Nurbs(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Fluid2Nurbs(const Fluid2Nurbs& old);

  /*!
  \brief Destructor

  */
  virtual ~Fluid2Nurbs();


  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  
  \return my parobject id
  */
  virtual int UniqueParObjectId() const { return ParObject_Fluid2Nurbs;}


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


#endif  // #ifndef FLUID2_NURBS_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID2
