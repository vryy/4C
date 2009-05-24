/*!----------------------------------------------------------------------
\file wall1_nurbs.H

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#ifdef CCADISCRET
#ifndef WALL1_NURBS_H
#define WALL1_NURBS_H

#include "wall1.H"
#include "../drt_nurbs_discret/drt_control_point.H"
#include "../drt_fem_general/drt_utils_bspline.H"
#include "../drt_lib/drt_utils.H"

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

    

class Wall1Nurbs : public DRT::ELEMENTS::Wall1
{

public:

  /*!
  \brief Standard Constructor

  \param id    : A unique global id
  \param owner : proc id that will own this element
  */
  Wall1Nurbs(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Wall1Nurbs(const Wall1Nurbs& old);

  /*!
  \brief Destructor

  */
  virtual ~Wall1Nurbs();

  /*!
  \brief Deep copy this instance of Wall1 and return pointer to it
  */

  virtual DRT::Element* Clone() const;


  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  
  \return my parobject id
  */
  virtual int UniqueParObjectId() const { return ParObject_Wall1Nurbs;}


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


  /*!
  \brief Get vector of RefCountPtrs to the lines of this element
  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();


  /*!
  \brief Get vector of RefCountPtrs to the surfaces of this element
  */
  virtual vector<Teuchos::RCP<DRT::Element> > Surfaces();

  
private:

};

} // namespace NURBS
} // namespace ELEMENTS
} // namespace DRT


#endif  // #ifndef WALL1_NURBS_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_WALL1
