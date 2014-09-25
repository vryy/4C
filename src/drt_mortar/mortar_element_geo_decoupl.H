/*!----------------------------------------------------------------------
\file mortar_element_geo_decoupl.H
\brief A mortar coupling element with decoupled geometric (=point) and dof (=node) information

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/
#ifndef MORTAR_ELEMENT_GEO_DECOUPL_H
#define MORTAR_ELEMENT_GEO_DECOUPL_H

#include "mortar_element.H"
#include "../drt_meshfree_discret/drt_meshfree_cell.H"
#include "../drt_lib/drt_elementtype.H"


/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms

*/
namespace MORTAR
{

class MortarElementGeoDecouplType : public DRT::ElementType
{
public:

  std::string Name() const { return "MortarElementGeoDecouplType"; }

  static MortarElementGeoDecouplType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

private:

  static MortarElementGeoDecouplType instance_;
};


/*!
\brief A mortar coupling element in which the geometry information (points) is decoupled from the dof information (node)

\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class MortarElementGeoDecoupl : public MORTAR::MortarElement, public DRT::MESHFREE::Cell
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  \param shape (in): shape of this element
  \param numnode (in): Number of nodes adjacent to this element
  \param nodeids (in): ids of nodes adjacent to this element
  \param numpoint (in): Number of points adjacent to this element
  \param pointids (in): ids of points adjacent to this element
  \param isslave (in): flag indicating whether element is slave or master side
  \param isnurbs (in): flag indicating whether element is nurbs element or not
  */
  MortarElementGeoDecoupl(int id, int owner,
                const DRT::Element::DiscretizationType& shape,
                const int numnode,
                const int* nodeids,
                const int numpoint,
                const int* pointids,
                const bool isslave,
                bool isnurbs = false);

  /*!
  \brief Copy Constructor

  Makes a deep copy of this class

  */
  MortarElementGeoDecoupl(const MORTAR::MortarElementGeoDecoupl& old);

  /*!
  \brief Destructor

  */
  virtual ~MortarElementGeoDecoupl() { return; }

  /*!
  \brief Deep copy the derived class and return pointer to it

  */
  virtual DRT::Element* Clone() const;

  /*!
  \brief Return unique ParObject id

  Every class implementing ParObject needs a unique id defined at the
  top of drt_parobject.H

  */
  virtual int UniqueParObjectId() const { return MortarElementGeoDecouplType::Instance().UniqueParObjectId(); }

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Unpack(const std::vector<char>& data);

  virtual DRT::ElementType & ElementType() const
  { return MortarElementGeoDecouplType::Instance(); }

  //! @name Evaluation methods

  /*!
  \brief Evaluate a Neumann boundary condition dummy

  An element derived from this class uses the EvaluateNeumann method to receive commands
  and parameters from some control routine in params and evaluates a Neumann boundary condition
  given in condition

  \note This class implements a dummy of this method that prints a warning and
        returns false.

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : Force vector to be filled by element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(Teuchos::ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL)
  { return 0; }

  //@}

  //! @name Query methods

  /*!
  \brief Get shape type of element

  */
  virtual DiscretizationType Shape() const { return MORTAR::MortarElement::Shape();}

  /*!
  \brief Print this element

  */
  virtual void Print(std::ostream& os) const;

  //@}


protected:

}; // class MortarElementGeoDecoupl

} // namespace MORTAR

// << operator
std::ostream& operator << (std::ostream& os, const MORTAR::MortarElement& ele);

#endif  // #ifndef MORTAR_ELEMENT_GEO_DECOUPL_H
