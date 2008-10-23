/*!----------------------------------------------------------------------
\file drt_celement.H

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
<\pre>

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef DRT_CELEMENT_H
#define DRT_CELEMENT_H


#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "drt_cnode.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/linalg_serialdensematrix.H"

using namespace std;
using namespace Teuchos;

// forward declarations

/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{

// forward declarations

/*!
\brief A contact element


\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class CElement : public DRT::Element
{
public:

  //! @name Enums and Friends
  enum ShapeType // shape function types recognized by CElement
    {
      lin1D,                  // linear 1D
      lindual1D,              // dual linear 1D
      lindual1D_edge0,
      lindual1D_edge1,
      dual1D_base_for_edge0,
      dual1D_base_for_edge1,
      quad1D,                 // quadratic 1D
      quaddual1D,             // dual quadratic 1D
      quaddual1D_edge0,
      quaddual1D_edge1,
      lin2D,                  // linear 2D
      lindual2D,              // dual linear 2D
      bilin2D,                // bilinear 2D
      bilindual2D,            // dual bilinear 2D
      quad2D,                 // quadratic 2D
      quaddual2D,             // dual quadratic 2D
      serendipity2D,          // serendipity 2D
      serendipitydual2D,      // dual serendipity 2D
      biquad2D,               // biquadratic 2D
      biquaddual2D            // dual biquadratic 2D
    };
      
  //@}

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  \param shape (in): shape of this element
  \param numnode (in): Number of nodes to this element
  \param nodeids (in): ids of nodes adjacent to this element
  \param isslave (in): flag indicating whether element is slave or master side
  */
  explicit CElement(int id, ElementType etype, int owner, 
                    const DRT::Element::DiscretizationType& shape, 
                    const int numnode,
                    const int* nodeids, 
                    const bool isslave);

  /*!
  \brief Copy Constructor

  Makes a deep copy of this class

  */
  CElement(const CONTACT::CElement& old);

  /*!
  \brief Destructor

  */
  virtual ~CElement() { return; }

  /*!
  \brief Deep copy the derived class and return pointer to it

  */
  virtual CONTACT::CElement* Clone() const;


  /*!
  \brief Return unique ParObject id

  Every class implementing ParObject needs a unique id defined at the
  top of drt_parobject.H
  
  */
  virtual int UniqueParObjectId() const { return ParObject_CElement; }

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(vector<char>& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Unpack(const vector<char>& data);

  /*!
  \brief Allocate and return a ElementRegister class

  Assuming the contact elements do not need initialization, we return a 
  dummy base class here.

  */
  virtual RCP<DRT::ElementRegister> ElementRegister() const
  { return rcp(new DRT::ElementRegister(Type())); }

  //@}

  //! @name Query methods

  /*!
  \brief Get shape type of element
  
  */
  virtual DiscretizationType Shape() const { return shape_;}

  /*!
  \brief Return number of lines to this element
  
  */
  virtual int NumLine() const {return 0;}

  /*!
  \brief Return number of surfaces to this element
  
  */
  virtual int NumSurface() const {return 0;}

  /*!
  \brief Get vector of RefCountPtrs to the lines of this element
  
  */
  virtual vector<RCP<DRT::Element> > Lines()
  {
   vector<RCP<DRT::Element> > lines(0);
   return lines;
  }

  /*!
  \brief Get vector of RefCountPtrs to the surfaces of this element
  
  */
  virtual vector<RCP<DRT::Element> > Surfaces()
  {
   vector<RCP<DRT::Element> > surfaces(0);
   return surfaces;
  }

  /*!
  \brief Get number of degrees of freedom of a certain node
  
  This CElement is picky: It cooperates only with CNodes, not with
  standard Node objects!
  
  */
  virtual int NumDofPerNode(const DRT::Node& node) const
  {
    const CONTACT::CNode* cnode = dynamic_cast<const CONTACT::CNode*>(&node);
    if (!cnode) dserror("Node is not a CNode");
    return cnode->NumDof();
  }

  /*!
  \brief Get number of degrees of freedom per element

  For now contact elements do not have degrees of freedom independent of
  the nodes
  
  */
  virtual int NumDofPerElement() const { return 0; }

  /*!
  \brief Print this element

  */
  virtual void Print(ostream& os) const;
  
  /*!
  \brief Return current area

  */
  virtual inline double& Area() { return area_; }
  
  /*!
  \brief Return slave (true) or master status

  */
  virtual inline bool IsSlave() { return isslave_; }
  
  /*!
  \brief Return number of potentially contacting elements
  
  */
  virtual inline int NumSearchElements() const { return (int)searchelements_.size(); }
  
  /*!
  \brief Return global ids of potentially contacting elements
  
  */
  virtual inline vector<int>& SearchElements() {return searchelements_; }
  
  //@}

  //! @name Construction

  //@}

  //! @name Conditions

  //@}

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
  virtual int Evaluate(ParameterList&            params,
                       DRT::Discretization&      discretization,
                       vector<int>&              lm,
                       Epetra_SerialDenseMatrix& elemat1,
                       Epetra_SerialDenseMatrix& elemat2,
                       Epetra_SerialDenseVector& elevec1,
                       Epetra_SerialDenseVector& elevec2,
                       Epetra_SerialDenseVector& elevec3);

  /*!
  \brief Get local coordinates for local node id
  */
  virtual bool LocalCoordinatesOfNode(int lid, double* xi);
    
  /*!
  \brief Get local numbering for global node id
  */
  virtual int GetLocalNodeId(int nid);
  
  /*!
  \brief Build element normal at node passed in
  */
  virtual void BuildNormalAtNode(int nid, int i, Epetra_SerialDenseMatrix& elens);
  
  /*!
  \brief Build element normal derivative at node passed in
  */
  virtual void DerivNormalAtNode(int nid, int i, Epetra_SerialDenseMatrix& elens,
                                 vector<map<int,double> >& derivn);
    
  /*!
  \brief Compute element normal at local coordinate xi
  */
  virtual void ComputeNormalAtXi(double* xi, int i, Epetra_SerialDenseMatrix& elens);
  
  /*!
  \brief Compute element normal derivative at local coordinate xi
  */
  virtual void DerivNormalAtXi(double* xi, int i, Epetra_SerialDenseMatrix& elens,
                               vector<map<int,double> >& derivn);
    
  /*!
  \brief Compute element normal at local coordinate xi
  */
  virtual LINALG::SerialDenseMatrix GetNodalCoords();
  
  /*!
  \brief Evaluate element metrics (local basis vectors)
  */
  virtual void Metrics(double* xi, vector<double>& gxi, vector<double>& geta);
  
  /*!
  \brief Evaluate Jacobian determinant for parameter space integration
  */
  virtual double Jacobian(double* xi);
  
  /*!
  \brief Evaluate derivative J,xi of Jacobian determinant
  */
  virtual void DJacDXi(double* djacdxi,
                       const LINALG::SerialDenseVector& val,
                       const LINALG::SerialDenseMatrix& deriv,
                       const LINALG::SerialDenseMatrix& secderiv,
                       const LINALG::SerialDenseMatrix& coord);
    
    
  /*!
  \brief Compute Jacobian determinant derivative
  */
  virtual void DerivJacobian(double* xi, map<int,double>& derivjac);
  
  /*!
  \brief Compute length/area of the element
  */
  virtual double ComputeArea();
  
  /*!
  \brief Compute length/area derivative of the element
  */
  virtual void DerivArea(map<int,double>& derivarea);
    
  /*!
  \brief A repository for all kinds of 1D/2D shape functions
  */
  virtual void ShapeFunctions(CElement::ShapeType shape,
                              const double* xi,
                              LINALG::SerialDenseVector& val,
                              LINALG::SerialDenseMatrix& deriv);
  
  /*!
  \brief A repository for 1D/2D shape function linearizations
  
  Currently this only contains the dual 1D quadratic shape functions.
  Linearization with respect to the displacements slave displacements
  is necessary here due to the adaption for distorted elements. The
  shape functions are thus dependent of the element geometry !!!
  
  \param derivdual (in): derivative maps to be filled
                         (= derivatives of the dual coefficient matrix Ae) 
  */
  void ShapeFunctionLinearizations(CElement::ShapeType shape,
                                   vector<vector<map<int,double> > >& derivdual);
  /*!
  \brief Evaluate shape functions and derivatives
  */
  virtual bool EvaluateShape(const double* xi, LINALG::SerialDenseVector& val,
                             LINALG::SerialDenseMatrix& deriv, const int valdim);
   
  /*!
  \brief Evaluate dual shape functions and derivatives
  */
  virtual bool EvaluateShapeDual(const double* xi, LINALG::SerialDenseVector& val,
                                 LINALG::SerialDenseMatrix& deriv, const int valdim);
  
  /*!
  \brief Evaluate 2nd derivative of shape functions
  */
  virtual bool Evaluate2ndDerivShape(const double* xi, LINALG::SerialDenseMatrix& val,
                                     const int valdim);
    
  /*!
  \brief Compute directional derivative of dual shape functions
  
  \param derivdual (in): derivative maps to be filled
                         (= derivatives of the dual coefficient matrix Ae) 
  */
  virtual bool DerivShapeDual(vector<vector<map<int,double> > >& derivdual);
  
  /*!
  \brief Interpolate global coordinates for given local element coordinates
  
  This method interpolates global coordinates for a given local element
  coordinate variable using the element node coordinates. For interpolation
  one can choose between shape functions or shape function derivatives!
  
  \param xi (in)        : local element coordinates
  \param inttype (in)   : set to true for shape function usage,
                          set to false for derivative usage
  \param globccord (out): interpolated global coordinates
  */
  virtual bool LocalToGlobal(const double* xi, double* globcoord, bool inttype);

  /*!
  \brief Add CElements to this CElement's potential contact partners
  */
  virtual bool AddSearchElements(const vector<int>& gids);
  
  //@}


private:
  
  DRT::Element::DiscretizationType   shape_;   // shape of this element
  
  bool          isslave_;           // indicating slave or master side
  double        area_;              // element length/area in current configuration
  
  vector<int>   searchelements_;    // global ids of potentially contacting elements 

}; // class CElement
} // namespace CONTACT


// << operator
ostream& operator << (ostream& os, const CONTACT::CElement& ele);




#endif  // #ifndef DRT_CELEMENT_H
#endif  // #ifdef CCADISCRET
