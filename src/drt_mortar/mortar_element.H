/*----------------------------------------------------------------------*/
/*!
\file mortar_element.H

\brief A mortar coupling element

\level 2

\maintainer Alexander Seitz

*/
/*----------------------------------------------------------------------*/
#ifndef MORTAR_ELEMENT_H
#define MORTAR_ELEMENT_H

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_inpar/inpar_mortar.H"
#include "../headers/pairedvector.H"
#include "../linalg/linalg_serialdensevector.H"
#include "mortar_node.H"

// forward declarations
namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}
namespace MORTAR{class MortarElementNitscheContainer;};

/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms

*/
namespace MORTAR
{

/*!
\brief A subclass of DRT::ElementType that adds mortar element type specific methods

*/
class MortarElementType : public DRT::ElementType
{
public:

  std::string Name() const { return "MortarElementType"; }

  static MortarElementType & Instance();

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

  /*!
  \brief Return matrix of dual shape function coefficients

  */
  virtual Teuchos::RCP<LINALG::SerialDenseMatrix>& DualShape() {return dualshapecoeff_; }

  /*!
  \brief Return trafo matrix for boundary modification

  */
  virtual Teuchos::RCP<LINALG::SerialDenseMatrix>& Trafo() {return trafocoeff_; }

  /*!
  \brief Return directional derivative of matrix of dual shape function coefficients

  */
  virtual Teuchos::RCP<GEN::pairedvector<int,Epetra_SerialDenseMatrix> >& DerivDualShape() {return derivdualshapecoeff_; }

  /*!
  \brief Reset matrix of dual shape function coefficients and free memory

  */
  virtual void ResetDualShape() { dualshapecoeff_ = Teuchos::null; return; }

  /*!
  \brief Reset directional derivative of matrix of dual shape function coefficients and free memory

  */
  virtual void ResetDerivDualShape() { derivdualshapecoeff_ = Teuchos::null; return; }

  /*!
  \brief Return Parent displacement vector

  */
  virtual std::vector<double>& ParentDisp() {return parentdisp_; }

  /*!
  \brief Return Parent element degrees of freedom

  */
  virtual std::vector<int>& ParentDof() {return parentdofs_; }

  /*!
  \brief Return Parent temperature vector

  */
  virtual std::vector<double>& ParentTemp() {return parenttemp_; }

  //@}

protected:

  // don't want = operator and cctor
  MortarEleDataContainer operator = (const MortarEleDataContainer& old);
  MortarEleDataContainer(const MortarEleDataContainer& old);

  double            area_;           // element length/area in current configuration
  std::vector<int>  searchelements_; // global ids of potentially contacting elements

  // coefficient matrix for dual shape functions
  Teuchos::RCP<LINALG::SerialDenseMatrix> dualshapecoeff_;

  // derivative of coefficient matrix for dual shape functions
  Teuchos::RCP<GEN::pairedvector<int,Epetra_SerialDenseMatrix> > derivdualshapecoeff_;

  // coefficient matrix for boundary trafo
  Teuchos::RCP<LINALG::SerialDenseMatrix> trafocoeff_;

  // Displacement of Parent Element
  std::vector<double> parentdisp_;

  // Displacement Parent element degrees of freedom
  std::vector<int> parentdofs_;

  // Temperature of Parent Element
  std::vector<double> parenttemp_;

}; // class MortarEleDataContainer


/*!
\brief A mortar coupling element

\author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
*/
class MortarElement : public DRT::FaceElement
{
public:

  //! @name Enums and Friends

  /*!
  \brief Enum for shape function types recognized by MortarElements

  */
  enum ShapeType
  {
    p0,                         // displacements / LM constant per element
    lin1D,                      // displacements / std LM linear 1D
    quad1D,                     // displacements / std LM quadratic 1D
    hermite1D,                  // displacements / std LM hermit 1D
    lin2D,                      // displacements / std LM linear 2D
    bilin2D,                    // displacements / std LM bilinear 2D
    quad2D,                     // displacements / std LM quadratic 2D
    serendipity2D,              // displacements / std LM serendipity 2D
    biquad2D,                   // displacements / std LM biquadratic 2D

    lindual1D,                  // dual LM linear 1D
    quaddual1D,                 // dual LM quadratic 1D
    hermite1Ddual,              // dual LM hermit 1D
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
    hermite1D_edge0,            // boundary LM modification 1D
    hermite1D_edge1,            // boundary LM modification 1D
    hermite1Ddual_edge0,        // boundary LM modification 1D
    hermite1Ddual_edge1,        // boundary LM modification 1D

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

  // physical type of mortar element
  // (Scoped enumeration: allows static_cast since C++11)
//  enum struct PhysicalType : int
  enum PhysicalType
  {
    poro = 0,       //poroelastic: (porofluid exists and must be considered in contact/meshtying)
    structure = 1,  //structure
    other = 2       //this should not happen
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
  \param isnurbs (in): flag indicating whether element is nurbs element or not
  */
  MortarElement(int id, int owner,
                const DRT::Element::DiscretizationType& shape,
                const int numnode,
                const int* nodeids,
                const bool isslave,
                bool isnurbs = false);

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
  virtual DRT::Element* Clone() const;

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
  \brief Get vector of Teuchos::RCPs to the lines of this element

  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines()
  {
   std::vector<Teuchos::RCP<DRT::Element> > lines(0);
   return lines;
  }

  /*!
  \brief Get vector of Teuchos::RCPs to the surfaces of this element

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
  virtual void Print(std::ostream& os) const;

  /*!
  \brief Return slave (true) or master status

  */
  virtual bool IsSlave() { return isslave_; }

  /*!
  \brief Return attached status

  */
  virtual bool IsAttached() { return attached_; }
  /*!
  \brief Change slave (true) or master status

  This changing of contact topology becomes necessary for self contact
  simulations, where slave and master status are assigned dynamically

  */
  virtual bool& SetSlave() { return isslave_; }

  /*!
  \brief Return attached status

  */
  virtual bool& SetAttached() { return attached_; }

  /*!
  \brief Return ansatz type (true = quadratic) of element

  */
  virtual bool IsQuad()
  {
    bool isquad = false;
    switch(Shape())
    {
    case DRT::Element::line2:
    case DRT::Element::nurbs2:
    case DRT::Element::tri3:
    case DRT::Element::quad4:
    {
      // do nothing
      break;
    }
    case DRT::Element::line3:
    case DRT::Element::nurbs3:
    case DRT::Element::quad8:
    case DRT::Element::quad9:
    case DRT::Element::nurbs9:
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
  \brief Return spatial dimension

  */
  virtual int Dim()
  {
    switch(Shape())
    {
    case DRT::Element::line2:
    case DRT::Element::nurbs2:
    case DRT::Element::line3:
    case DRT::Element::nurbs3:
    {
      // this is a 2-D problem
      return 2;
      break;
    }
    case DRT::Element::tri3:
    case DRT::Element::quad4:
    case DRT::Element::quad8:
    case DRT::Element::quad9:
    case DRT::Element::nurbs9:
    case DRT::Element::tri6:
    {
      // this is a 3-D problem
      return 3;
      break;
    }
    default:
    {
      dserror("ERROR: Unknown mortar element type identifier");
      return 0;
    }
    }
  }

  /*!
  \brief Return nurbs (true) or not-nurbs (false) status

  */
  virtual bool IsNurbs() { return nurbs_; }

  /*!
  \brief Return data container of this element

  This method returns the data container of this mortar element where additional
  mortar specific quantities/information are stored.

  */
  inline MORTAR::MortarEleDataContainer& MoData() {
      dsassert(modata_!=Teuchos::null, "Mortar data container not set");
      return *modata_; }

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
  virtual int Evaluate(Teuchos::ParameterList&   params,
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
  virtual int EvaluateNeumann(Teuchos::ParameterList& params,
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
  \brief Compute averaged nodal normal at local coordinate xi
  */
  virtual double ComputeAveragedUnitNormalAtXi(double* xi, double* n);

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
  virtual void DerivUnitNormalAtXi(double* xi, std::vector<GEN::pairedvector<int,double> >& derivn);

  /*!
  \brief Get nodal reference / spatial coords of current element

  */
  virtual void GetNodalCoords(LINALG::SerialDenseMatrix& coord);

  /*! \brief Get nodal reference / spatial coords of current element
   *
   *  \author hiermeier \date 03/17 */
  template< unsigned elenumnode >
  inline void GetNodalCoords(LINALG::Matrix<3,elenumnode>& coord)
  {
    LINALG::SerialDenseMatrix sdm_coord( View, coord.A(), 3, 3, elenumnode );
    GetNodalCoords( sdm_coord );
  }

  double inline GetNodalCoords(const int direction, const int node)
  {
  #ifdef DEBUG
    MortarNode* mymrtrnode = dynamic_cast<MortarNode*> (Points()[node]);
    if (!mymrtrnode)
      dserror("ERROR: GetNodalCoords: Null pointer!");
    return mymrtrnode->xspatial()[direction];
  #else
    return static_cast<MortarNode*> (Points()[node])->xspatial()[direction];
  #endif
  }

  /*!
  \brief Get nodal spatial coords from previous time step of current element

  \param isinit (in): true if called for reference coords
  */
  virtual void GetNodalCoordsOld(LINALG::SerialDenseMatrix& coord,
                                 bool isinit = false);

  double inline GetNodalCoordsOld(const int direction, const int node)
  {
  #ifdef DEBUG
    MortarNode* mymrtrnode = dynamic_cast<MortarNode*> (Points()[node]);
    if (!mymrtrnode)
      dserror("ERROR: GetNodalCoords: Null pointer!");
    return mymrtrnode->X()[direction]+mymrtrnode->uold()[direction];
  #else
    return (static_cast<MortarNode*> (Points()[node])->X()[direction]+static_cast<MortarNode*> (Points()[node])->uold()[direction]);
  #endif
  }

  /*!
  \brief Get nodal spatial coords from previous time step of current element

  \param isinit (in): true if called for reference coords
  */
  virtual void GetNodalLagMult(LINALG::SerialDenseMatrix& lagmult,
                               bool isinit = false);

  /*!
  \brief Evaluate element metrics (local basis vectors)
  */
  virtual void Metrics(const double* xi, std::vector<double>& gxi, std::vector<double>& geta);

  /*!
  \brief Evaluate Jacobian determinant for parameter space integration
  */
  virtual double Jacobian( const double* xi );

  /*!
  \brief Compute Jacobian determinant derivative
  */
  virtual void DerivJacobian( const double* xi, GEN::pairedvector<int,double>& derivjac);

  /*!
  \brief Compute length/area of the element
  */
  virtual double ComputeArea();

  /*!
  \brief Compute length/area of the element and its derivative
  */
  virtual double ComputeAreaDeriv(GEN::pairedvector<int,double>& area_deriv);

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
      GEN::pairedvector<int,Epetra_SerialDenseMatrix>& derivdual);

  /*!
  \brief Evaluate displacement shape functions and derivatives
  */
  virtual bool EvaluateShape(const double* xi, LINALG::SerialDenseVector& val,
                             LINALG::SerialDenseMatrix& deriv,
                             const int valdim,
                             bool dualquad3d = false);

  /*! \brief Evaluate displacement shape functions and derivatives
   *
   *  \author hiermeier \date 03/17 */
  template <unsigned elenumnode, unsigned eledim>
  inline bool EvaluateShape(
      const double* xi,
      LINALG::Matrix<elenumnode,1>& val,
      LINALG::Matrix<elenumnode,eledim>& deriv,
      unsigned valdim,
      bool dualquad3d)
  {
    LINALG::SerialDenseVector sdv_val( View, val.A(), elenumnode );
    LINALG::SerialDenseMatrix sdm_deriv( View, deriv.A(),elenumnode, elenumnode, eledim );
    return EvaluateShape( xi, sdv_val, sdm_deriv, valdim, dualquad3d );
  }

  /*!
  \brief Evaluate Lagrange multiplier shape functions and derivatives
  */
  virtual bool EvaluateShapeLagMult(const INPAR::MORTAR::ShapeFcn& lmtype,
                                    const double* xi,
                                    LINALG::SerialDenseVector& val,
                                    LINALG::SerialDenseMatrix& deriv,
                                    const int valdim,
                                    bool boundtrafo = true);

  /*! \brief Evaluate Lagrange multiplier shape functions and derivatives
   *
   *  \author hiermeier \date 03/17 */
  template <unsigned elenumnode, unsigned eledim>
  inline bool EvaluateShapeLagMult(
      INPAR::MORTAR::ShapeFcn lmtype,
      const double* xi,
      LINALG::Matrix<elenumnode,1>& val,
      LINALG::Matrix<elenumnode,eledim>& deriv,
      unsigned valdim,
      bool boundtrafo )
  {
    LINALG::SerialDenseVector sdv_val( View, val.A(), elenumnode );
    LINALG::SerialDenseMatrix sdm_deriv( View, deriv.A(),elenumnode, elenumnode, eledim );
    return EvaluateShapeLagMult( lmtype, xi, sdv_val, sdm_deriv, valdim, boundtrafo );
  }

  /*!
  \brief Evaluate Lagrange multiplier shape functions and derivatives
  (special version for 3D quadratic mortar with linear Lagrange multipliers)
  */
  virtual bool EvaluateShapeLagMultLin(const INPAR::MORTAR::ShapeFcn& lmtype,
                                       const double* xi, LINALG::SerialDenseVector& val,
                                       LINALG::SerialDenseMatrix& deriv, const int valdim);

  /*!
  \brief Evaluate Lagrange multiplier shape functions and derivatives
  (special version for quadratic mortar with element-wise constant Lagrange multipliers)
  */
  virtual bool EvaluateShapeLagMultConst(const INPAR::MORTAR::ShapeFcn& lmtype,
                                       const double* xi, LINALG::SerialDenseVector& val,
                                       LINALG::SerialDenseMatrix& deriv, const int valdim);

  /*! \brief Evaluate Lagrange multiplier shape functions and derivatives
   *  (special version for 3D quadratic mortar with linear Lagrange multipliers)
   *
   *  \author hiermeier \date 03/17 */
  template<unsigned elenumnode, unsigned eledim>
  inline bool EvaluateShapeLagMultLin(
      INPAR::MORTAR::ShapeFcn lmtype,
      const double* xi,
      LINALG::Matrix<elenumnode,1>& val,
      LINALG::Matrix<elenumnode,eledim>& deriv,
      int valdim)
  {
    LINALG::SerialDenseVector sdv_val( View, val.A(), elenumnode );
    LINALG::SerialDenseMatrix sdm_deriv( View, deriv.A(),elenumnode, elenumnode, eledim );
    return EvaluateShapeLagMultLin( lmtype, xi, sdv_val, sdm_deriv, valdim );
  }

  /*!
  \brief Evaluate 2nd derivative of shape functions
  */
  virtual bool Evaluate2ndDerivShape(const double* xi, LINALG::SerialDenseMatrix& val,
                                     const int& valdim);

  template <unsigned elenumnode>
  inline bool Evaluate2ndDerivShape(
      const double* xi,
      LINALG::Matrix<elenumnode,3>& secderiv,
      const int& valdim )
  {
    LINALG::SerialDenseMatrix sdm_secderiv( View, secderiv.A(),elenumnode, elenumnode, 3 );
    return Evaluate2ndDerivShape( xi, sdm_secderiv, valdim );
  }

  /*!
  \brief Compute directional derivative of dual shape functions

  \param derivdual (in): derivative maps to be filled
                         (= derivatives of the dual coefficient matrix Ae)
  */
  virtual bool DerivShapeDual(
      GEN::pairedvector<int,Epetra_SerialDenseMatrix>& derivdual);

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

  */
  virtual double MinEdgeSize();

  /*!
  \brief Evaluate maximal edge size of CElement

  */
  virtual double MaxEdgeSize();

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
  \brief get some information from adjacent elements for hermit smoothing
         status[0] ... type of boundary modification
         status[1] ... numnode for smoothing ele

  */
  virtual bool AdjEleStatus(int* status);

  /*!
  \brief get adjacent nodes for hermit smoothing

  */
  virtual bool HermitEleNodes(DRT::Node** nodes, int status);

  /*!
  \brief get node coordinates from adj elements

  */
  virtual bool AdjNodeCoords(LINALG::SerialDenseMatrix& coord, int inttype);

  /*!
  \brief delete all found meles for this element

  */
  virtual void DeleteSearchElements();

  /*!
  \brief Resets the data container of the element

  With this function, the container with mortar specific quantities/information
  is deleted / reset to Teuchos::null pointer

  */
  virtual void ResetDataContainer();

  // hermit specific
  /*!
  \brief bool for hermit smoothing

  */
  virtual bool& IsHermite() { return hermite_; };

  // nurbs specific:

  /*!
  \brief bool for integration: if true: no integration for this element

  this is only the case if there are more than polynom degree + 1 multiple knot entries

  */
  virtual bool& ZeroSized() { return zero_sized_; };

  /*!
  \brief factor for normal calculation (default 1.0)

  */
  virtual double& NormalFac() { return normalfac_; };

  /*!
  \brief get knot vectors for this mortar element

  */
  virtual std::vector<Epetra_SerialDenseVector>& Knots() { return mortarknots_; };

  /*!
  \brief Get the linearization of the spatial position of the Nodes for this Ele.

         Returns a vector of vector of maps. Outer vector for the nodes,
         inner vector for the spatial dimensions, map for the derivatives.

         Needed to be overloaded by IntElement
  */
  virtual void NodeLinearization(
      std::vector<std::vector<GEN::pairedvector<int, double> > >& nodelin);

  //h.Willmann return physical type of the mortar element
  PhysicalType& PhysType(){return physicaltype_;};

  /*!
  \brief Estimate mesh size and stiffness parameter h/E via Eigenvalues of the trace inequality.
         For Nitsche contact formulations
    */
  void EstimateNitscheTraceMaxEigenvalueCombined();

  /*!
  \brief Estimated mesh size and stiffness parameter h/E via Eigenvalues of the trace inequality.
         For Nitsche contact formulations
    */
  virtual double& TraceHE() {return traceHE_;}

  /*!
  \brief Estimated mesh size and thermal conductivity h/K via Eigenvalues of the trace inequality.
         For Nitsche contact formulations
    */
  virtual double& TraceHCond() {return traceHCond_;}

  /*!
  \brief Get Nitsche data container
    */
  virtual MORTAR::MortarElementNitscheContainer& GetNitscheContainer();
  //@}

protected:

  DRT::Element::DiscretizationType              shape_;        // shape of this element
  bool                                          isslave_;      // indicating slave or master side
  bool                                          attached_;     // bool whether an element contributes to M

  Teuchos::RCP<MORTAR::MortarEleDataContainer>  modata_;       // additional information

  // hermit smoothing
  bool                                          hermite_;        // bool for hermit smoothing

  // nurbs specific:
  bool                                          nurbs_;
  std::vector<Epetra_SerialDenseVector>         mortarknots_;  // mortar element knot vector
  double                                        normalfac_;    // factor for normal orientation
  bool                                          zero_sized_;   // zero-sized element: if true: no integration for this element

  //h.Willmann
  PhysicalType                                  physicaltype_; //physical type

  // approximation of mesh size and stiffness from inverse trace inequality (h/E)
  double traceHE_;
  // approximation of mesh size and stiffness from inverse trace inequality (h/conductivity)
  double traceHCond_;

  // data container for element matrices in Nitsche contact
  Teuchos::RCP<MORTAR::MortarElementNitscheContainer> nitsche_container_;

  /*!
  \brief Protected constructor for use in derived classes that expect standard element constructors

  \param id    (in): A globally unique element id
  \param owner (in): owner processor of the element
  */
  MortarElement(int id, int owner);

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
  int& nGP() {return ngp_;}

  /*!
  \brief Return coordinates of a specific GP

  */
  double& Coordinate(int& gp, int dir) {return coords_(gp,dir);}

  /*!
  \brief Return weight of a specific GP in 1D/2D CElement

  */
  double& Weight(int& gp) {return weights_[gp];}

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
std::ostream& operator << (std::ostream& os, const MORTAR::MortarElement& ele);

#endif  // #ifndef MORTAR_ELEMENT_H
