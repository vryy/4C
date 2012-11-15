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
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/
#ifndef MORTAR_ELEMENT_H
#define MORTAR_ELEMENT_H

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_inpar/inpar_mortar.H"


// forward declarations
namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}

/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms

*/
namespace MORTAR
{

class MortarElementType : public DRT::ElementType
{
public:

  std::string Name() const { return "MortarElementType"; }

  static MortarElementType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

private:

  static MortarElementType instance_;
};

/*!
\brief A class containing additional data for mortar elements

This class contains additional information for mortar elements which is
only needed for contact evaluation. Thus, in order to save memory, it is
sufficient to have this data available only on the slave element column map.

\author popp (popp@lnm.mw.tum.de)
*/
class MortarEleDataContainer
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  MortarEleDataContainer();

  /*!
  \brief Destructor

  */
  virtual ~MortarEleDataContainer() { return; }

  /*!
  \brief Pack this class so that it can be communicated

  This function packs the datacontainer. This is only called
  when the class has been initialized and the pointer to this
  class exists.

  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a vector into this class

  This function unpacks the datacontainer. This is only called
  when the class has been initialized and the pointer to this
  class exists.

  */
  virtual void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data);

  //@}

  //! @name Access methods

  /*!
  \brief Return current area

  */
  virtual double& Area() { return area_; }

  /*!
  \brief Return number of potentially contacting elements

  */
  virtual int NumSearchElements() const { return (int)searchelements_.size(); }

  /*!
  \brief Return global ids of potentially contacting elements

  */
  virtual std::vector<int>& SearchElements() {return searchelements_; }

  //@}

protected:

  // don't want = operator and cctor
  MortarEleDataContainer operator = (const MortarEleDataContainer& old);
  MortarEleDataContainer(const MortarEleDataContainer& old);

  double            area_;           // element length/area in current configuration
  std::vector<int>  searchelements_; // global ids of potentially contacting elements

}; // class MortarEleDataContainer


/*!
\brief A mortar coupling element

\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class MortarElement : public DRT::Element
{
public:

  //! @name Enums and Friends

  /*!
  \brief Enum for shape function types recognized by MortarElements

  */
  enum ShapeType
  {
    lin1D,                      // displacements / std LM linear 1D
    quad1D,                     // displacements / std LM quadratic 1D
    lin2D,                      // displacements / std LM linear 2D
    bilin2D,                    // displacements / std LM bilinear 2D
    quad2D,                     // displacements / std LM quadratic 2D
    serendipity2D,              // displacements / std LM serendipity 2D
    biquad2D,                   // displacements / std LM biquadratic 2D

    lindual1D,                  // dual LM linear 1D
    quaddual1D,                 // dual LM quadratic 1D
    lindual2D,                  // dual LM linear 2D
    bilindual2D,                // dual LM bilinear 2D
    quaddual2D,                 // dual LM quadratic 2D
    serendipitydual2D,          // dual LM serendipity 2D
    biquaddual2D,               // dual LM biquadratic 2D

    lin1D_edge0,                // crosspoint LM modification 1D
    lin1D_edge1,                // crosspoint LM modification 1D
    lindual1D_edge0,            // crosspoint LM modification 1D
    lindual1D_edge1,            // crosspoint LM modification 1D
    dual1D_base_for_edge0,      // crosspoint LM modification 1D
    dual1D_base_for_edge1,      // crosspoint LM modification 1D
    quad1D_edge0,               // crosspoint LM modification 1D
    quad1D_edge1,               // crosspoint LM modification 1D
    quaddual1D_edge0,           // crosspoint LM modification 1D
    quaddual1D_edge1,           // crosspoint LM modification 1D

    quad1D_only_lin,            // quad->lin standard LM modification 1D
    quad2D_only_lin,            // quad->lin standard LM modification 2D
    serendipity2D_only_lin,     // quad->lin standard LM modification 2D
    biquad2D_only_lin,          // quad->lin standard LM modification 2D

    quaddual1D_only_lin,        // quad->lin dual LM modification 1D (not yet impl.)
    quaddual2D_only_lin,        // quad->lin dual LM modification 2D
    serendipitydual2D_only_lin, // quad->lin dual LM modification 2D
    biquaddual2D_only_lin,      // quad->lin dual LM modification 2D

    quad1D_modified,            // displacement modification for dual LM quadratic 1D (not yet impl.)
    quad2D_modified,            // displacement modification for dual LM quadratic 2D
    serendipity2D_modified,     // displacement modification for dual LM serendipity 2D
    biquad2D_modified,          // displacement modification for dual LM biquadratic 2D

    quad1D_hierarchical,        // displacement modification for quad->lin dual LM quadratic 1D (not yet impl.)
    quad2D_hierarchical,        // displacement modification for quad->lin dual LM quadratic 2D
    serendipity2D_hierarchical, // displacement modification for quad->lin dual LM serendipity 2D
    biquad2D_hierarchical       // displacement modification for quad->lin dual LM biquadratic 2D
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
  MortarElement(int id, int owner,
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
  virtual int UniqueParObjectId() const { return MortarElementType::Instance().UniqueParObjectId(); }

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
  { return MortarElementType::Instance(); }

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
  \brief Get vector of RCPs to the lines of this element

  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines()
  {
   std::vector<Teuchos::RCP<DRT::Element> > lines(0);
   return lines;
  }

  /*!
  \brief Get vector of RCPs to the surfaces of this element

  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Surfaces()
  {
   std::vector<Teuchos::RCP<DRT::Element> > surfaces(0);
   return surfaces;
  }

  /*!
  \brief Get number of degrees of freedom of a certain node

  This MortarElement is picky: It cooperates only with MortarNodes,
  not with standard Node objects!

  */
  virtual int NumDofPerNode(const DRT::Node& node) const;

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
  \brief Return slave (true) or master status

  */
  virtual bool IsSlave() { return isslave_; }

  /*!
  \brief Change slave (true) or master status

  This changing of contact topology becomes necessary for self contact
  simulations, where slave and master status are assigned dynamically

  */
  virtual bool& SetSlave() { return isslave_; }

  /*!
  \brief Return ansatz type (true = quadratic) of element

  */
  virtual bool IsQuad()
  {
    bool isquad = false;
    switch(Shape())
    {
    case DRT::Element::line2:
    case DRT::Element::tri3:
    case DRT::Element::quad4:
    {
      // do nothing
      break;
    }
    case DRT::Element::line3:
    case DRT::Element::quad8:
    case DRT::Element::quad9:
    case DRT::Element::tri6:
    {
      isquad = true;
      break;
    }
    default:
      dserror("ERROR: Unknown mortar element type identifier");
    }
    return isquad;
  }

  /*!
  \brief Return data container of this element

  This method returns the data container of this mortar element where additional
  mortar specific quantities/information are stored.

  */
  inline MORTAR::MortarEleDataContainer& MoData() { return *modata_; }

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
                       std::vector<int>&         lm,
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
  virtual int EvaluateNeumann(Teuchos::ParameterList&   params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
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
  \brief Compute element unit normal derivative at local coordinate xi
         This function is a real stand-alone function to be called
         for a CoElement in order to compute a unit normal derivative at any point.
  */
  virtual void DerivUnitNormalAtXi(double* xi, std::vector<std::map<int,double> >& derivn);

  /*!
  \brief Get nodal reference / spatial coords of current element

  \param isinit (in): true if called for reference coords
  */
  virtual void GetNodalCoords(LINALG::SerialDenseMatrix& coord,
                              bool isinit = false);

  /*!
  \brief Get nodal spatial coords from previous time step of current element

  \param isinit (in): true if called for reference coords
  */
  virtual void GetNodalCoordsOld(LINALG::SerialDenseMatrix& coord,
                                 bool isinit = false);

  /*!
  \brief Get nodal spatial coords from previous time step of current element

  \param isinit (in): true if called for reference coords
  */
  virtual void GetNodalLagMult(LINALG::SerialDenseMatrix& lagmult,
                               bool isinit = false);

  /*!
  \brief Evaluate element metrics (local basis vectors)
  */
  virtual void Metrics(double* xi, std::vector<double>& gxi, std::vector<double>& geta);

  /*!
  \brief Evaluate Jacobian determinant for parameter space integration
  */
  virtual double Jacobian(double* xi);

  /*!
  \brief Compute Jacobian determinant derivative
  */
  virtual void DerivJacobian(double* xi, std::map<int,double>& derivjac);

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
  \brief A repository for 1D/2D shape function linearizations

  \param derivdual (in): derivative maps to be filled
                         (= derivatives of the dual coefficient matrix Ae)
  */
  void ShapeFunctionLinearizations(MortarElement::ShapeType shape,
                                   std::vector<std::vector<std::map<int,double> > >& derivdual);

  /*!
  \brief Evaluate displacement shape functions and derivatives
  */
  virtual bool EvaluateShape(const double* xi, LINALG::SerialDenseVector& val,
                             LINALG::SerialDenseMatrix& deriv, const int& valdim,
                             bool dualquad3d = false);

  /*!
  \brief Evaluate Lagrange multiplier shape functions and derivatives
  */
  virtual bool EvaluateShapeLagMult(const INPAR::MORTAR::ShapeFcn& lmtype,
                                    const double* xi, LINALG::SerialDenseVector& val,
                                    LINALG::SerialDenseMatrix& deriv, const int& valdim);

  /*!
  \brief Evaluate Lagrange multiplier shape functions and derivatives
  (special version for 3D quadratic mortar with linear Lagrange multipliers)
  */
  virtual bool EvaluateShapeLagMultLin(const INPAR::MORTAR::ShapeFcn& lmtype,
                                       const double* xi, LINALG::SerialDenseVector& val,
                                       LINALG::SerialDenseMatrix& deriv, const int& valdim);

  /*!
  \brief Evaluate 2nd derivative of shape functions
  */
  virtual bool Evaluate2ndDerivShape(const double* xi, LINALG::SerialDenseMatrix& val,
                                     const int& valdim);

  /*!
  \brief Compute directional derivative of dual shape functions

  \param derivdual (in): derivative maps to be filled
                         (= derivatives of the dual coefficient matrix Ae)
  */
  virtual bool DerivShapeDual(std::vector<std::vector<std::map<int,double> > >& derivdual);

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
  \brief Add one MortarElement to this MortarElement's potential contact partners

  This is for the element-based brute-force search and for the new
  binary search tree. We do NOT have to additionally check, if the
  given MortarElement has already been added to this MortarCElement's
  potential contact partners before. This cannot happen by construction!

  */
  virtual bool AddSearchElements(const int& gid);

  /*!
  \brief Initializes the data container of the element

  With this function, the container with mortar specific quantities/information
  is initialized.

  */
  virtual void InitializeDataContainer();

  /*!
  \brief Resets the data container of the element

  With this function, the container with mortar specific quantities/information
  is deleted / reset to Teuchos::null pointer

  */
  virtual void ResetDataContainer();

  //@}

protected:

  DRT::Element::DiscretizationType              shape_;      // shape of this element
  bool                                          isslave_;    // indicating slave or master side

  Teuchos::RCP<MORTAR::MortarEleDataContainer>  modata_;     // additional information

}; // class MortarElement

/*!
\brief A class to perform Gaussian integration on a mortar element

\author popp (popp@lnm.mw.tum.de)
*/

class ElementIntegrator
{
public:

  /*!
  \brief Standard constructor

  */
  ElementIntegrator(DRT::Element::DiscretizationType eletype);

  /*!
  \brief Destructor

  */
  virtual ~ElementIntegrator() {}

  //! @name Access methods

  /*!
  \brief Return number of Gauss points

  */
  int nGP() {return ngp_;}

  /*!
  \brief Return coordinates of a specific GP

  */
  double Coordinate(int& gp, int dir) {return coords_(gp,dir);}

  /*!
  \brief Return weight of a specific GP in 1D/2D CElement

  */
  double Weight(int& gp) {return weights_[gp];}

  //@}

protected:

    // don't want = operator and cctor
    ElementIntegrator operator = (const ElementIntegrator& old);
    ElementIntegrator(const ElementIntegrator& old);

    int                           ngp_;         // number of Gauss points
    Epetra_SerialDenseMatrix      coords_;      // Gauss point coordinates
    std::vector<double>           weights_;     // Gauss point weights

}; // class ElementIntegrator
} // namespace MORTAR

// << operator
ostream& operator << (ostream& os, const MORTAR::MortarElement& ele);

#endif  // #ifndef MORTAR_ELEMENT_H
