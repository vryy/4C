/*!----------------------------------------------------------------------
 \file contact_element.H

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
#ifndef CONTACT_ELEMENT_H
#define CONTACT_ELEMENT_H

#include "../drt_mortar/mortar_element.H"

// forward declarations

/*!
 \brief CONTACT: namespace of the baci contact and meshtying algorithms

 */
namespace CONTACT
{

// forward declarations

class CoElementType: public DRT::ElementType
{
public:

  std::string Name() const
  {
    return "CoElementType";
  }

  static CoElementType & Instance()
  {
    return instance_;
  }
  ;

  virtual DRT::ParObject* Create(const std::vector<char> & data);

  virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

  virtual void NodalBlockInformation(DRT::Element * dwele, int & numdf,
      int & dimns, int & nv, int & np);

  virtual void ComputeNullSpace(DRT::Discretization & dis,
      std::vector<double> & ns, const double * x0, int numdf, int dimns);

private:

  static CoElementType instance_;
};

/*!
 \brief A contact element


 \author popp (popp@lnm.mw.tum.de)
 */
class CoElement: public MORTAR::MortarElement
{
public:

  //! @name Constructors and destructors and related methods

  /*!
   \brief Standard Constructor

   \param id    (in): A globally unique element id
   \param owner (in): owner processor of the element
   \param shape (in): shape of this element
   \param numnode (in): Number of nodes to this element
   \param nodeids (in): ids of nodes adjacent to this element
   \param isslave (in): flag indicating whether element is slave or master side
   \param isnurbs (in): flag indicating whether element is nurbs element or not
   */
  CoElement(int id, int owner, const DRT::Element::DiscretizationType& shape,
      const int numnode, const int* nodeids, const bool isslave, bool isnurbs =
          false);

  /*!
   \brief Copy Constructor

   Makes a deep copy of this class

   */
  CoElement(const CONTACT::CoElement& old);

  /*!
   \brief Destructor

   */
  virtual ~CoElement()
  {
    return;
  }

  /*!
   \brief Deep copy the derived class and return pointer to it

   */
  virtual DRT::Element* Clone() const;

  /*!
   \brief Return unique ParObject id

   Every class implementing ParObject needs a unique id defined at the
   top of drt_parobject.H

   */
  virtual int UniqueParObjectId() const
  {
    return CoElementType::Instance().UniqueParObjectId();
  }

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

  virtual CoElementType & ElementType() const
  {
    return CoElementType::Instance();
  }

  //@}

  //! @name Query methods

  /*!
   \brief Get number of degrees of freedom of a certain node

   This CoElement is picky: It cooperates only with CNodes, not with
   standard Node objects!

   */
  int NumDofPerNode(const DRT::Node& node) const;

  /*!
   \brief Print this element

   */
  virtual void Print(std::ostream& os) const;

  //! @name Evaluation methods

  /*!
   \brief Evaluate an element

   An element derived from this class uses the Evaluate method to receive commands
   and parameters from some control routine in params and evaluates element matrices and
   vectors accoring to the command in params.

   \note This class implements a dummy of this method that prints a dserror and
   returns false.

   \param params (in/out)    : ParameterList for communication between control routine
   and elements
   \param discretization (in): A reference to the underlying discretization
   \param lm (in)            : location vector of this element
   \param elemat1 (out)      : matrix to be filled by element depending on commands
   given in params
   \param elemat2 (out)      : matrix to be filled by element depending on commands
   given in params
   \param elevec1 (out)      : vector to be filled by element depending on commands
   given in params
   \param elevec2 (out)      : vector to be filled by element depending on commands
   given in params
   \param elevec3 (out)      : vector to be filled by element depending on commands
   given in params
   \return 0 if successful, negative otherwise
   */
  virtual int Evaluate(Teuchos::ParameterList& params,
      DRT::Discretization& discretization, std::vector<int>& lm,
      Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
      Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
      Epetra_SerialDenseVector& elevec3);

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
      DRT::Discretization& discretization, DRT::Condition& condition,
      std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
      Epetra_SerialDenseMatrix* elemat1 = NULL)
  {
    return 0;
  }

  /*!
   \brief Build element normal derivative at node passed in
   */
  virtual void DerivNormalAtNode(int nid, int& i,
      Epetra_SerialDenseMatrix& elens,
      std::vector<GEN::pairedvector<int, double> >& derivn);

  /*!
   \brief Compute element normal derivative at local coordinate xi
   Caution: This function cannot be called stand-alone! It is
   integrated into the whole nodal normal calculation process.
   */
  virtual void DerivNormalAtXi(double* xi, int& i,
      Epetra_SerialDenseMatrix& elens,
      std::vector<GEN::pairedvector<int, double> >& derivn);

  /*!
   \brief Evaluate derivative J,xi of Jacobian determinant
   */
  virtual void DJacDXi(double* djacdxi, double* xi,
      const LINALG::SerialDenseMatrix& secderiv);

  //@}

private:

  // no new variables as compared to base class

};
// class CoElement
}// namespace CONTACT

// << operator
std::ostream& operator <<(std::ostream& os, const CONTACT::CoElement& ele);

#endif  // #ifndef CONTACT_ELEMENT_H
