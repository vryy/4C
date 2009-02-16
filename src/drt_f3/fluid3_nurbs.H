/*!----------------------------------------------------------------------
\file fluid3_nurbs.H

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET
#ifndef FLUID3_NURBS_H
#define FLUID3_NURBS_H

#include "fluid3.H"
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

// forward declarations
class Fluid3NurbsSurface;
   
class Fluid3Nurbs : public DRT::ELEMENTS::Fluid3
{

public:

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  */
  Fluid3Nurbs(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Fluid3Nurbs(const Fluid3Nurbs& old);

  /*!
  \brief Destructor

  */
  virtual ~Fluid3Nurbs();

  /*!
  \brief Deep copy this instance of Fluid3Nurbs and return 
  pointer to the copy

  The Clone() method is used from the virtual base class 
  Element in cases where the type of the derived class is 
  unknown and a copy-ctor is needed
  */
  DRT::Element* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  
  \return my parobject id
  */
  virtual int UniqueParObjectId() const { return ParObject_Fluid3Nurbs;}


  /*!
  \brief Get shape type of element

  \return nurbs8 or nurbs27

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
    else 
    {
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

  /*!
  \brief Get vector of RefCountPtrs to the surfaces of this element

  According to Georg, we do NOT store line or surface elements inside 
  the parent element after their creation.
  Reason: if a Redistribute() is performed on the discretization, 
          stored node ids and node pointers owned by these boundary 
	  elements might have become illegal

  -> new surface elements are allocated each time

  */
  virtual vector<RCP<DRT::Element> > Surfaces();
  
private:

};
/*----------------------------------------------------------------------

  \brief Fluid3Surface class. his is a warning. Nurbs surfaces are only 
         valid if the corresponding knotvektor is interpolating.

  \author gammi (gamnitzer@lnm.mw.tum.de), 02/09

  ----------------------------------------------------------------------*/
class Fluid3NurbsSurface : public DRT::ELEMENTS::Fluid3Surface
{

public:

  /*!
  \brief Standard Constructor

  \param id :      A unique global id
  \param owner:    Processor owning this surface
  \param nnode:    Number of nodes attached to this element
  \param nodeids:  global ids of nodes attached to this element
  \param nodes:    the discretizations map of nodes to build ptrs to nodes from
  \param parent:   The parent fluid element of this surface
  \param lsurface: the local surface number of this surface w.r.t. the parent element

  */
  Fluid3NurbsSurface(
    int                    id, 
    int                    owner, 
    int                    nnode, 
    const int*             nodeids,
    DRT::Node**            nodes, 
    DRT::ELEMENTS::Fluid3* parent, 
    const int              lsurface
    );

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Fluid3NurbsSurface(const Fluid3NurbsSurface& old);

  /*!
  \brief Destructor

  */
  virtual ~Fluid3NurbsSurface();

  /*!
  \brief Deep copy this instance of Fluid3NurbsSurface 
  and return pointer to the copy

  The Clone() method is used from the virtual base class 
  Element in cases where the type of the derived class is 
  unknown and a copy-ctor is needed
  */
  DRT::Element* Clone() const;

  /*!
  \brief Return unique ParObject id

  \return my parobject id
  */
  virtual int UniqueParObjectId() const { return ParObject_Fluid3NurbsSurface;}

  /*!
  \brief Get shape type of element

  \return nurbs4 or nurbs9

  */
  virtual DiscretizationType Shape() const;

private:

};

} // namespace NURBS
} // namespace ELEMENTS
} // namespace DRT


#endif  // #ifndef FLUID3_NURBS_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
