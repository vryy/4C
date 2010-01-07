/*!----------------------------------------------------------------------
\file mortar_element.H
\brief A mortar coupling element

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
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef MORTAR_ELEMENT_H
#define MORTAR_ELEMENT_H

#include "mortar_node.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/linalg_serialdensematrix.H"

using namespace std;
using namespace Teuchos;

// forward declarations

/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms

*/
namespace MORTAR
{

// forward declarations

/*!
\brief A mortar coupling element

\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class MortarElement : public DRT::Element
{
public:

  //! @name Enums and Friends
  enum ShapeType // shape function types recognized by MortarElement
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
  explicit MortarElement(int id, ElementType etype, int owner,
                         const DRT::Element::DiscretizationType& shape,
                         const int numnode,
                         const int* nodeids,
                         const bool isslave);

  /*!
  \brief Copy Constructor

  Makes a deep copy of this class

  */
  MortarElement(const MORTAR::MortarElement& old);

  /*!
  \brief Destructor

  */
  virtual ~MortarElement() { return; }

  /*!
  \brief Deep copy the derived class and return pointer to it

  */
  virtual MORTAR::MortarElement* Clone() const;

  /*!
  \brief Return unique ParObject id

  Every class implementing ParObject needs a unique id defined at the
  top of drt_parobject.H

  */
  virtual int UniqueParObjectId() const { return ParObject_MortarElement; }

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

  This MortarElement is picky: It cooperates only with MortarNodes,
  not with standard Node objects!

  */
  virtual int NumDofPerNode(const DRT::Node& node) const
  {
    const MORTAR::MortarNode* mnode = dynamic_cast<const MORTAR::MortarNode*>(&node);
    if (!mnode) dserror("Node is not a MortarNode");
    return mnode->NumDof();
  }

  /*!
  \brief Get number of degrees of freedom per element

  For now mortar coupling elements do not have degrees of freedom
  independent of the nodes.

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
  \brief Return ansatz type (true = quad) of element

  */
  virtual inline bool IsQuad()
  {
    if (Shape()==quad8 || Shape()==quad9 || Shape()==tri6) return true;
    else return false;
  }

  /*!
  \brief Return number of potentially contacting elements

  */
  virtual inline int NumSearchElements() const { return (int)searchelements_.size(); }

  /*!
  \brief Return global ids of potentially contacting elements

  */
  virtual inline vector<int>& SearchElements() {return searchelements_; }

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
  virtual int EvaluateNeumann(ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL)
  { return 0; }
  
  /*!
  \brief Get local coordinates for local node id
  */
  virtual bool LocalCoordinatesOfNode(int& lid, double* xi);

  /*!
  \brief Get local numbering for global node id
  */
  virtual int GetLocalNodeId(int& nid);
  
  /*!
  \brief Build element normal at node passed in
  */
  virtual void BuildNormalAtNode(int nid, int& i, Epetra_SerialDenseMatrix& elens);

  /*!
  \brief Compute element normal at local coordinate xi
         Caution: This function cannot be called stand-alone! It is
         integrated into the whole nodal normal calculation process.
  */
  virtual void ComputeNormalAtXi(double* xi, int& i, Epetra_SerialDenseMatrix& elens);
  
  /*!
  \brief Compute unit element normal at local coordinate xi
         This function is a real stand-alone function to be called
         for a CElement in order to compute a unit normal at any point.
         Returns the length of the non-unit interpolated normal at xi.
  */
  virtual double ComputeUnitNormalAtXi(double* xi, double* n);

  /*!
  \brief Get nodal reference / spatial coords of current element

  \param isinit (in): true if called for reference coords
  */
  virtual void GetNodalCoords(LINALG::SerialDenseMatrix& coord,
                              bool isinit = false);

  /*!
  \brief Evaluate element metrics (local basis vectors)
  */
  virtual void Metrics(double* xi, vector<double>& gxi, vector<double>& geta);

  /*!
  \brief Evaluate Jacobian determinant for parameter space integration
  */
  virtual double Jacobian(double* xi);

  /*!
  \brief Compute length/area of the element
  */
  virtual double ComputeArea();

  /*!
  \brief A repository for all kinds of 1D/2D shape functions
  */
  virtual void ShapeFunctions(MortarElement::ShapeType shape,
                              const double* xi,
                              LINALG::SerialDenseVector& val,
                              LINALG::SerialDenseMatrix& deriv);

  /*!
  \brief Evaluate shape functions and derivatives
  */
  virtual bool EvaluateShape(const double* xi, LINALG::SerialDenseVector& val,
                             LINALG::SerialDenseMatrix& deriv, const int& valdim);

  /*!
  \brief Evaluate dual shape functions and derivatives
  */
  virtual bool EvaluateShapeDual(const double* xi, LINALG::SerialDenseVector& val,
                                 LINALG::SerialDenseMatrix& deriv, const int& valdim);

  /*!
  \brief Interpolate global coordinates for given local element coordinates

  This method interpolates global coordinates for a given local element
  coordinate variable using the element node coordinates. For interpolation
  one can choose between shape functions or shape function derivatives!

  \param xi (in)        : local element coordinates
  \param inttype (in)   : set to 0 for shape function usage,
                          set to 1 for derivative xi usage
                          set to 2 for derivative eta usage (3D only)
  \param globccord (out): interpolated global coordinates
  */
  virtual bool LocalToGlobal(const double* xi, double* globcoord, int inttype);

  /*!
  \brief Evaluate minimal edge size of CElement

  \param isinit (in): true if called for reference coords
  */
  virtual double MinEdgeSize(bool isinit = false);

  /*!
  \brief Evaluate maximal edge size of CElement

  \param isinit (in): true if called for reference coords
  */
  virtual double MaxEdgeSize(bool isinit = false);

  /*!
  \brief Add MortarElements to this MortarElement's potential contact partners

  This is for the classical node-based brute-force search. We have to
  additionally check, if the given MortarElements have already been added
  to this MortarElement's potential contact partners before, in order to
  avoid double integrations!!

  */
  virtual bool AddSearchElements(const vector<int>& gids);

  /*!
  \brief Add one MortarElement to this MortarElement's potential contact partners

  This is for the element-based brute-force search and for the new
  binary search tree. We do NOT have to additionally check, if the
  given MortarElement has already been added to this MortarCElement's
  potential contact partners before. This cannot happen by construction!

  */
  virtual bool AddSearchElements(const int& gid);
  
  //@}

private:

  DRT::Element::DiscretizationType   shape_;          // shape of this element
  bool                               isslave_;        // indicating slave or master side
  double                             area_;           // element length/area in current configuration
  vector<int>                        searchelements_; // global ids of potentially contacting elements

}; // class MortarElement
} // namespace MORTAR

// << operator
ostream& operator << (ostream& os, const MORTAR::MortarElement& ele);

#endif  // #ifndef MORTAR_ELEMENT_H
#endif  // #ifdef CCADISCRET
