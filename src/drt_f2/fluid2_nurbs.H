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
#include "../drt_nurbs_discret/drt_control_point.H"
#include "../drt_fem_general/drt_utils_bspline.H"
#include "../drt_lib/drt_parobjectfactory.H"


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

class Fluid2NurbsType : public DRT::ElementType
{
public:

  std::string Name() const { return "Fluid2NurbsType"; }

  static Fluid2NurbsType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const string eletype,
                                             const string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
  {
    // do nothing. Definition inserted by normal wall element.
  }

private:

  static Fluid2NurbsType instance_;
};


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
  virtual int UniqueParObjectId() const { return Fluid2NurbsType::Instance().UniqueParObjectId(); }

  virtual Fluid2NurbsType & ElementType() const
  { return Fluid2NurbsType::Instance(); }

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
