/*----------------------------------------------------------------------------*/
/*!
\file beam3.H

\brief three dimensional nonlinear Reissner beam element

\level 2

\maintainer Christoph Meier
*/
/*----------------------------------------------------------------------------*/

//header file only included if not yet included!
#ifndef BEAM3_H
#define BEAM3_H


#include "../drt_fem_general/largerotations.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_integration.H"

#include "beam3_base.H"

//#define BEAM3DISCRETELINENEUMANN 0.0 //Flag in order to interpret line Neumann condition as discrete force at given element parameter position DISCRETELINENEUMANN
//Default: off                  //for all elements contained in the corresponding design line!

//! \brief DRT: namespace of the ccarat discretization module
namespace DRT
{

namespace UTILS
{

}

namespace ELEMENTS
{
class Beam3Type : public DRT::ElementType
{
public:

  std::string Name() const { return "Beam3Type"; }

  static Beam3Type & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual int Initialize(DRT::Discretization& dis);

  virtual void NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );



private:

  static Beam3Type instance_;
};

/*!
\brief three dimensional nonlinear corotational Timoshenko beam element according to Crisfield, Vol.2,
section 17.2 ("an interpretation of an element due to simo and vu-quoc"), p. 226 - 231; the current
implementation sets the initial director angles and curvature to zero, which can be done without loss
of generality since any initial rotation cancels out in the equations (17.72) and (17.76). As a consequence
the current nodal directors do not appear in the element, but only the values of the current director
angles minus reference director angles, respectively

The element is capable of using higher order functions from linear to quartic. Please make sure you put the nodes in the right order
in the input file.
LIN2  1---2
LIN3  1---3---2
LIN4  1---4---2---3
LIN5  1---5---2---3---4
LIN6  1---6---2---3---4---5

author cyron (cyron @lnm.mw.tum.de)
*/
class Beam3 : public Beam3Base
{
public:

  //! @name Friends
  friend class Beam3Type;

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  */
  Beam3(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  */
  Beam3(const Beam3& old);

  /*!
  \brief Destructor
  */
  virtual ~Beam3();

  /*!
  \brief Deep copy this instance of Beam3 and return pointer to the copy

  The Clone() method is used by the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed
.
  */
  virtual DRT::Element* Clone() const;

   /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const;

  virtual double ReturnRatioNormForceMoment() const;
  virtual double ReturnNormForce() const;
  virtual double ReturnNormMoment() const;
  virtual double Lcurr() const;
  virtual void TcurrBeam3r(LINALG::Matrix<3,1>& Tcurr1, LINALG::Matrix<3,1>& Tcurr2);
  virtual void TrefBeam3r(LINALG::Matrix<3,1>& Tref1, LINALG::Matrix<3,1>& Tref2);

  /*!
  \brief Return unique ParObject id

  Every class implementing ParObject needs a unique id defined at the
  top of drt_parobject.H
  */
  virtual int UniqueParObjectId() const { return Beam3Type::Instance().UniqueParObjectId(); }

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
  { return Beam3Type::Instance(); }

  //@}

  /*!
  \brief Return number of lines to this element
  */
  virtual int NumLine() const {return 1;}

  /** \brief get number of nodes used for centerline interpolation
   *
   *  \author grill
   *  \date 05/16 */
  virtual inline int NumCenterlineNodes() const { return this->NumNode(); }

  virtual double GetDeltaTheta() const
  {
   return deltatheta_;
  }


  /*!
  \brief Get vector of Teuchos::RCPs to the lines of this element
  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();


  /*!
  \brief Get number of degrees of freedom of a single node
  */
  virtual int NumDofPerNode(const DRT::Node& node) const
  {
    /*note: this is not necessarily the number of DOF assigned to this node by the discretization finally, but only
    *the number of DOF requested for this node by this element; the discretization will finally assign the maximal
    *number of DOF to this node requested by any element connected to this node*/
    return 6;
  }

  /*!
  \brief Get number of degrees of freedom per element not including nodal degrees of freedom
  */
  virtual int NumDofPerElement() const { return 0;}

  /*!
  \brief Print this element
  */
  virtual void Print(std::ostream& os) const;

  /*!
  \brief Get moment of inertia of area Iyy_ of this element
  */
  const double& Iyy() const { return Iyy_;}

  /*!
  \brief Get moment of inertia of area Izz_ of this element
  */
  const double& Izz() const { return Izz_;}

  /*!
  \brief Get polar moment of inertia of area Irr_ of this element
  */
  const double& Irr() const { return Irr_;}

  /*!
  \brief Get jacobi_ factor of first gauss point
  */
  const double& GetJacobi() const { return jacobi_[0];}

  /** \brief get Jacobi factor ds/dxi(xi) at xi \in [-1;1]
   *
   *  \author grill
   *  \date 06/16 */
  inline virtual double GetJacobiFacAtXi(const double& xi) const
  {
    dserror("current implementation only valid for initially straight beams,"
            " where the Jacobi factor is constant for all xi");
    return GetJacobi();
  }

  /** \brief get access to the reference length
   *
   *  \author grill
   *  \date 05/16 */
  virtual inline double RefLength() const{ return 2.0*jacobi_[0]; }

  //! quaternions of triads at GPs during the current iteration step due to Crisfield, Vol. 2, equation (17.65)
  const std::vector<LINALG::Matrix<4,1> >& Qnew() const {return Qnew_;}

  /** \brief return flag whether Hermite polynomials are applied for centerline interpolation
   *
   *  \author grill
   *  \date 05/16 */
  virtual inline bool HermiteCenterlineInterpolation() const { return false; }

  /*!
  \brief Set moment of inertia of area Iyy_ of this element
  */
  void SetIyy(const double& Iyy);

  /*!
  \brief Set moment of inertia of area Iyy_ of this element
  */
  void SetIzz(const double& Izz);

  /*!
  \brief Set moment of inertia of area Iyy_ of this element
  */
  void SetIrr(const double& Irr);

  /*!
  \brief Set cross section area of this element
  */
  void SetCrossSec(const double& crosssec);

  /*!
  \brief Set cross section with shear correction factor of this element
  */
  void SetCrossSecShear(const double& crosssecshear);

  //!computes the number of different random numbers required in each time step for generation of stochastic forces
  int HowManyRandomNumbersINeed();

  //! \brief internal force vector
  const Epetra_SerialDenseVector& InternalForceVector() const {return f_;}
  //! \brief return material forces at Gauss point
  const LINALG::Matrix<3,1>& MatForceGp() const {return Ngp_;}

  //! brief epsilonnormal
  const double& EpsilonSgn() const { return eps_;}

  //! @name Construction

  /*!
  \brief Read input for this element
  */
  bool ReadElement(const std::string&           eletype,
                   const std::string&           distype,
                   DRT::INPUT::LineDefinition*  linedef);

  //@}


  //! @name Evaluation methods


  /*!
  \brief Evaluate an element

  An element derived from this class uses the Evaluate method to receive commands
  and parameters from some control routine in params and evaluates element matrices and
  vectors accoring to the command in params.

  \note This class implements a dummy of this method that prints a warning and
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
  \brief Evaluate a Neumann boundary condition

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
                              Epetra_SerialDenseMatrix* elemat1 = NULL);


  /*!
  \brief Evaluate PTC addition to stiffness for free Brownian motion

  An element derived from this class uses the Evaluate method to receive commands
  and parameters from some control routine in params and evaluates a  statistical Neumann boundary condition
  as used in the problem type STATISTICAL MECHANICS

  \param params (in/out)       : ParameterList for communication between control routine and elements
  \param vector<double> mydisp : current nodal displacement
  \param elemat1 (out)         : artificial damping matrix to be filled by element

  \return 0 if successful, negative otherwise
  */
  template<int nnode> void EvaluatePTC(Teuchos::ParameterList&   params,
                                       Epetra_SerialDenseMatrix& elemat1);

  //@}

  //! @name methods for initilization of the element

  //!sets up from current nodal position all geometric parameters (considering current position as reference configuration)
  template<int nnode>
  void SetUpReferenceGeometry(const std::vector<double>&  xrefe,
                              const std::vector<double>&  rotrefe,
                              const bool                  secondinit = false)
  {
    /*this method initializes geometric variables of the element; the initilization can usually be applied to elements only once;
     *therefore after the first initilization the flag isinit is set to true and from then on this method does not take any action
     *when called again unless it is called on purpose with the additional parameter secondinit. If this parameter is passed into
     *the method and is true the element is initialized another time with respective xrefe and rotrefe;
     *note: the isinit_ flag is important for avoiding reinitialization upon restart. However, it should be possible to conduct a
     *second initilization in principle (e.g. for periodic boundary conditions*/

    //just during initialization of a new crosslinker
    if(!isinit_)
    {
      //store node x and rotational coordinates of this element
      for (int k=0; k<6; k++)
      {
        xactrefe_(k,0) = xrefe[k];
        rotinitrefe_(k,0) = rotrefe[k];
      }
    }

    if(!isinit_ || secondinit)
    {
      isinit_ = true;


      //resize all class STL vectors so that they can each store 1 value at each GP
      jacobi_.resize(nnode-1);
      jacobimass_.resize(nnode);
      jacobinode_.resize(nnode);
      Qref_.resize((nnode));
      Qconv_.resize((nnode-1));
      Qold_.resize((nnode-1));
      Qnew_.resize((nnode-1));
      QoldNode_.resize((nnode));
      QnewNode_.resize((nnode));
      curvconv_.resize((nnode-1));
      curvold_.resize((nnode-1));
      curvnew_.resize((nnode-1));
      thetaconv_.resize((nnode-1));
      thetaold_.resize((nnode-1));
      thetanew_.resize((nnode-1));
      ThetaOldNode_.resize((nnode));
      ThetaNewNode_.resize((nnode));
      thetaprimeconv_.resize((nnode-1));
      thetaprimeold_.resize((nnode-1));
      thetaprimenew_.resize((nnode-1));
      f_.Resize(6*this->NumNode());
      f_.Scale(0.0);

      //create Matrix for the derivates of the shapefunctions at the GP
      LINALG::Matrix<1,nnode> shapefuncderiv;

      //create Matrix for the shapefunctions at the GP
      LINALG::Matrix<1,nnode> funct;

      //Get DiscretizationType
      DRT::Element::DiscretizationType distype = Shape();

      //Get the applied integrationpoints for underintegration
      DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode,gaussunderintegration));

        //Loop through all GPs and calculate jacobi the triads at the GPs
      for(int numgp=0; numgp < gausspoints.nquad; numgp++)
      {
        //Get position xi of GP
        const double xi = gausspoints.qxg[numgp][0];

        //Get derivatives of shapefunctions at GP
        DRT::UTILS::shape_function_1D_deriv1(shapefuncderiv,xi,distype);

        //Get shapefunctions at GP
        DRT::UTILS::shape_function_1D(funct,xi,distype);

        //triad in reference configuration at GP
        LINALG::Matrix<3,3> Tref;

        /*initial triad Tref = [t1,t2,t3] is set in a way for which we don`t have strains in reference configuration*/
        LINALG::Matrix<3,1> dxdxi;

        dxdxi.Clear();
        thetaconv_[numgp].Clear();
        thetaprimeconv_[numgp].Clear();

        //calculate vector dxdxi
        for(int node=0; node<nnode; node++)
        {
          for(int dof=0; dof<3 ; dof++)
          {
            dxdxi(dof) += shapefuncderiv(node) * xrefe[3*node+dof];
            thetaconv_[numgp](dof) += funct(node) * rotrefe[3*node+dof];
            thetaprimeconv_[numgp](dof) += shapefuncderiv(node) * rotrefe[3*node+dof];
          }//for(int dof=0; dof<3 ; dof++)
        }//for(int node=0; node<nnode; node++)



        //Store length factor for every GP
        //note: the length factor jacobi replaces the determinant and refers to the reference configuration by definition
        jacobi_[numgp]= dxdxi.Norm2();

        for (int k=0; k<3; k++)
        {
          //t1 axis points in positive direction along xi and is a unit vector
          Tref(k,0)=dxdxi(k)/jacobi_[numgp];
        }

        //t2 is a unit vector in the x2x3-plane orthogonal to t1
        Tref(0,1) = 0;
        //if t1 is parallel to the x1-axis t2 is set parallel to the x2-axis
        if (Tref(1,0) == 0 && Tref(2,0) == 0)
        {
             Tref(1,1) = 1;
             Tref(2,1) = 0;
        }

        //otherwise t2 is calculated from the scalar product with t1
        else
        {
          //setting t2(0)=0 and calculating other elements by setting scalar product t1 o t2 to zero
          double lin1norm = std::pow(pow(Tref(1,0),2)+pow(Tref(2,0),2),0.5);
          Tref(1,1) = -Tref(2,0)/lin1norm;
          Tref(2,1) =  Tref(1,0)/lin1norm;
        }

        //calculating t3 by crossproduct t1 x t2
        Tref(0,2) = Tref(1,0)*Tref(2,1)-Tref(1,1)*Tref(2,0);
        Tref(1,2) = Tref(2,0)*Tref(0,1)-Tref(2,1)*Tref(0,0);
        Tref(2,2) = Tref(0,0)*Tref(1,1)-Tref(0,1)*Tref(1,0);

        /*the center triad in reference configuration is stored as a quaternion whose equivalent would be the rotation
        * from the identity matrix into the reference configuration*/
        LARGEROTATIONS::triadtoquaternion(Tref,Qconv_[numgp]);

        //Calculating Reference Quternion
        for(int node=0; node<nnode; node++)
          LARGEROTATIONS::triadtoquaternion(Tref,Qref_[node]);



        //the here employed beam element does not need data about the current position of the nodal directors so that
        //initilization of those can be skipped (the nodal displacements handeled in beam3_evaluate.cpp are not the current angles,
        //but only the differences between current angles and angles in reference configuration, respectively. Thus the
        //director orientation in reference configuration cancels out and can be assumed to be zero without loss of
        //generality

        curvconv_[numgp].Clear();
      }//for(int numgp=0; numgp < gausspoints.nquad; numgp++)

      //Now all triads have been calculated and Qold_ and Qnew_ can be updated
      Qold_ = Qconv_;
      Qnew_ = Qconv_;
      QoldNode_[0] = Qconv_[0];
      QoldNode_[1] = Qconv_[0];
      QnewNode_[0] = Qconv_[0];
      QnewNode_[1] = Qconv_[0];

      curvold_ = curvconv_;
      curvnew_ = curvconv_;

      thetaold_ = thetaconv_;
      thetanew_ = thetaconv_;
      ThetaOldNode_[0] = thetaconv_[0];
      ThetaOldNode_[1] = thetaconv_[0];
      ThetaNewNode_[0] = thetaconv_[0];
      ThetaNewNode_[1] = thetaconv_[0];

      thetaprimeold_ = thetaprimeconv_;
      thetaprimenew_ = thetaprimeconv_;

      //Get the applied integrationpoints for exact integration of mass matrix
      DRT::UTILS::IntegrationPoints1D gausspointsmass(MyGaussRule(nnode,gaussexactintegration));

      //Loop through all GPs and calculate jacobi and theta0
      for(int numgp=0; numgp < gausspointsmass.nquad; numgp++)
      {

        //Get position xi of GP
        const double xi = gausspointsmass.qxg[numgp][0];

        //Get derivatives of shapefunctions at GP
        DRT::UTILS::shape_function_1D_deriv1(shapefuncderiv,xi,distype);

        LINALG::Matrix<3,1> dxdximass;

        dxdximass.Clear();
        //calculate dx/dxi and dz/dxi
        for(int node=0; node<nnode; node++)
          for(int dof=0; dof<3; dof++)
            dxdximass(dof)+=shapefuncderiv(node)*xrefe[3*node+dof];


        //Store length factor for every GP
        //note: the length factor jacobi replaces the determinant and refers by definition always to the reference configuration
        jacobimass_[numgp]= dxdximass.Norm2();

      }//for(int numgp=0; numgp < gausspointsmass.nquad; numgp++)


      //compute Jacobi determinant at gauss points for Lobatto quadrature (i.e. at nodes)
      for(int numgp=0; numgp< nnode; numgp++)
      {

        //Get position xi of nodes
        const double xi = -1.0 + 2*numgp / (nnode - 1);

        //Get derivatives of shapefunctions at GP
        DRT::UTILS::shape_function_1D_deriv1(shapefuncderiv,xi,distype);

        LINALG::Matrix<3,1> dxdxi;

        dxdxi.Clear();
        //calculate dx/dxi and dz/dxi
        for(int node=0; node<nnode; node++)
          for(int dof=0; dof<3; dof++)
            dxdxi(dof)+=shapefuncderiv(node)*xrefe[3*node+dof];

        //Store Jacobi determinant for each node (Jacobi determinant refers by definition always to the reference configuration)
        jacobinode_[numgp]= dxdxi.Norm2();

      }//for(int numgp=0; numgp< nnode; numgp++)

      return;

    }//if(!isinit_)

  }//DRT::ELEMENTS::Beam3::SetUpReferenceGeometry()

  //! changes the actual linker length
  void SetReferenceLength(const double& scalefac);

  /** \brief add indices of those DOFs of a given node that are positions
   *
   *  \author grill
   *  \date 07/16 */
  virtual inline void PositionDofIndices(std::vector<int>& posdofs, const DRT::Node& node) const
  {
    posdofs.push_back(0);  posdofs.push_back(1);  posdofs.push_back(2);
    return;
  }

  /** \brief add indices of those DOFs of a given node that are tangents (in the case of Hermite interpolation)
   *
   *  \author grill
   *  \date 07/16 */
  virtual inline void TangentDofIndices(std::vector<int>& tangdofs, const DRT::Node& node) const
  {
    return;
  }

  /** \brief add indices of those DOFs of a given node that are rotation DOFs (non-additive rotation vectors)
   *
   *  \author grill
   *  \date 07/16 */
  virtual inline void RotationVecDofIndices(std::vector<int>& rotvecdofs, const DRT::Node& node) const
  {
    /* note: due to the incremental formulation of beam3, we can not apply the multiplicative update
     *       in the nonlinear solver (see PrePostOperator RotVecUpdater for details) and hence pretend to
     *       not have multiplicative rotation vector DoFs here. Instead, we declare them to be (additive)
     *       1D rotation DoFs (see below) */
    return;
  }

  /** \brief add indices of those DOFs of a given node that are 1D rotation DOFs
   *         (planar rotations are additive, e.g. in case of relative twist DOF of beam3k with rotvec=false)
   *
   *  \author grill
   *  \date 07/16 */
  virtual inline void Rotation1DDofIndices(std::vector<int>& twistdofs, const DRT::Node& node) const
  {
    // see comment in method RotationVecDofIndices
    twistdofs.push_back(3);  twistdofs.push_back(4);  twistdofs.push_back(5);
    return;
  }

  /** \brief add indices of those DOFs of a given node that represent norm of tangent vector
   *         (additive, e.g. in case of beam3k with rotvec=true)
   *
   *  \author grill
   *  \date 07/16 */
  virtual inline void TangentLengthDofIndices(std::vector<int>& tangnormdofs, const DRT::Node& node) const
  {
    return;
  }

  //@}

private:

  //! kind of integration to be performed
  enum IntegrationType
  {
    gaussexactintegration,
    gaussunderintegration,
    lobattointegration
  };

  //! variable saving whether element has already been initialized (then isinit_ == true)
  bool isinit_;

  //! \brief vector storing the internal force vector
  Epetra_SerialDenseVector f_;
  //! \brief material forces at Gauss point
  LINALG::Matrix<3,1> Ngp_;
  //! storing epsilonn(0)
  double eps_;

  //! quaternions of triads at GPs at the end of the preceeding time step due to Crisfield, Vol. 2, equation (17.65)
  std::vector<LINALG::Matrix<4,1> > Qref_;
  //! quaternions of triads at GPs at the end of the preceeding time step due to Crisfield, Vol. 2, equation (17.65)
  std::vector<LINALG::Matrix<4,1> > Qconv_;
  //! quaternions of triads at GPs at the end of the preceeding iteration step due to Crisfield, Vol. 2, equation (17.65)
  std::vector<LINALG::Matrix<4,1> > Qold_;
  //! quaternions of triads at GPs during the current iteration step due to Crisfield, Vol. 2, equation (17.65)
  std::vector<LINALG::Matrix<4,1> > Qnew_;
  //! quaternions of triads at Nodes at the end of the preceeding iteration step due to Crisfield, Vol. 2, equation (17.65)
  std::vector<LINALG::Matrix<4,1> > QoldNode_;
  //! quaternions of triads at Nodes during the current iteration step due to Crisfield, Vol. 2, equation (17.65)
  std::vector<LINALG::Matrix<4,1> > QnewNode_;
  //!curvaturea at the end of the preceeding time step due to Crisfield, Vol. 2, equation (17.72)
  std::vector<LINALG::Matrix<3,1> > curvconv_;
  //!curvaturea at the end of the preceeding iteration step due to Crisfield, Vol. 2, equation (17.72)
  std::vector<LINALG::Matrix<3,1> > curvold_;
  //!curvaturea during the current iteration step due to Crisfield, Vol. 2, equation (17.72)
  std::vector<LINALG::Matrix<3,1> > curvnew_;

  //!angle displacement taken from discretization and interpolated for a specific Gauss point by Crisfield, Vol. 2, equation (17.93); relates to the end of the preceeding time step
  std::vector<LINALG::Matrix<3,1> > thetaconv_;
  //!angle displacement taken from discretization and interpolated for a specific Gauss point by Crisfield, Vol. 2, equation (17.93); relates to the end of the preceeding iteration step
  std::vector<LINALG::Matrix<3,1> > thetaold_;
  //!angle displacement taken from discretization and interpolated for a specific Gauss point by Crisfield, Vol. 2, equation (17.93); relates to current iteration step
  std::vector<LINALG::Matrix<3,1> > thetanew_;
  //!angle displacement at nodes relates to the end of the preceeding iteration step
  std::vector<LINALG::Matrix<3,1> > ThetaOldNode_;
  //!angle displacement at nodes relating to current iteration step
  std::vector<LINALG::Matrix<3,1> > ThetaNewNode_;
  //!derivative of angle displacement taken from discretization with respect to curve parameter at a specific Gauss point; derivative due to Crisfield, Vol. 2, equation (17.93); relates to the end of the preceeding time step
  std::vector<LINALG::Matrix<3,1> > thetaprimeconv_;
  //!derivative of angle displacement taken from discretization with respect to curve parameter at a specific Gauss point; derivative due to Crisfield, Vol. 2, equation (17.93); relates to the end of the preceeding iteration step
  std::vector<LINALG::Matrix<3,1> > thetaprimeold_;
  //!derivative of angle displacement taken from discretization with respect to curve parameter at a specific Gauss point; derivative due to Crisfield, Vol. 2, equation (17.93); relates to the current iteration step
  std::vector<LINALG::Matrix<3,1> > thetaprimenew_;

  //! cross section area
  double crosssec_;

   //! cross section area including shear corretion factor
  double crosssecshear_;

  //!moment of inertia of area around y-axis(orthognoal to beam axis in yz-plane)
  double Iyy_;

  //!moment of inertia of area around z-axis(forming cartesian system with beam- and y-axis)
  double Izz_;

  //!polar moment of inertia of area
  double Irr_;

  //! Change in inclusive angles between filaments
  double deltatheta_;

  //! Current length of crosslinkers
  double lcurr_;

  //! Vector holding value of Jacobi determinant jacobi for each Gauss point for underintegration
  std::vector<double> jacobi_;

  //! Vector holding value of Jacobi determinant jacobi for complete integration of massmatrix
  std::vector<double> jacobimass_;

  //! Vector holding value of Jacobi determinant jacobi at nodes
  std::vector<double> jacobinode_;

  //! \brief vector, which stores the actual x node coordinates of the element; subject to change during runtime
  LINALG::Matrix<6,1> xactrefe_;

  //! \brief vector, which stores the initial rotational node coordinates of the element
  LINALG::Matrix<6,1> rotinitrefe_;


  //force and moments for debugging purpose
  double NormMoment;
  double NormForce;
  double RatioNormForceMoment;

  //! @name Internal calculation methods

  /*!this function performs an update of the rotation (in quaterion form) at the numgp-th Gauss point by the incremental rotation deltatheta,
   *by means of a quaternion product and then computes the respective new triad Tnew at the Gauss point*/
  void updatetriad(const LINALG::Matrix<3,1>& deltatheta, LINALG::Matrix<3,3>& Tnew, const int numgp);

  /*!this function performs an update of the rotation (in quaterion form) at the node by the incremental rotation deltatheta,
    *by means of a quaternion product and then computes the respective new triad Tnew at node*/
  void UpdateNodalTriad(std::vector<double>& disp, std::vector<LINALG::Matrix<3,1> >& Tnew);

  /*!Calculate change in angle from reference configuration */
  void CalcDeltaTheta(std::vector<double>& disp);

  //!updating local curvature according to Crisfield, Vol. 2, pages 209 - 210, basically equation (16.148)
  void updatecurvature(const LINALG::Matrix<3,3>& Tnew, LINALG::Matrix<3,1>& deltatheta,LINALG::Matrix<3,1>& deltathetaprime, const int numgp);

  //!updating local curvature approximately according to Crisfield, Vol. 2, eqs. (17.72) & (17.73) (only possible for linear elements)
  void approxupdatecurvature(const LINALG::Matrix<3,3>& Tnew, LINALG::Matrix<3,1> deltatheta,LINALG::Matrix<3,1> deltathetaprime, const int numgp);

  //!computing stiffness matrix Ksigma1 according to Crisfield, Vol. 2, equation (17.106)
  template<int nnode>
  void computeKsig1(LINALG::Matrix<6*nnode,6*nnode>& Ksig1, const LINALG::Matrix<3,1>& stressn, const LINALG::Matrix<3,1>& stressm, const LINALG::Matrix<1,nnode>& funct, const LINALG::Matrix<1,nnode>& deriv);

  //!computing stiffness matrix Ksigma2 according to Crisfield, Vol. 2, equation (17.107a) and (17.107b)
  template<int nnode>
  void computeKsig2(LINALG::Matrix<6*nnode,6*nnode>& Ksig2, const LINALG::Matrix<3,1>& stressn, const LINALG::Matrix<3,3>& S, const LINALG::Matrix<1,nnode>& funct, const LINALG::Matrix<1,nnode>& deriv);

  //!computing basis of stiffness matrix according to Crisfield, Vol. 2, equation (17.105)
  template<int nnode>
  void computestiffbasis(const LINALG::Matrix<3,3>& Tnew, const LINALG::Matrix<3,1>& Cm, const LINALG::Matrix<3,1>& Cb, const LINALG::Matrix<3,3>& S, LINALG::Matrix<6*nnode,6*nnode>& stiffmatrix, const LINALG::Matrix<1,nnode>& funct, const LINALG::Matrix<1,nnode>& deriv);

  //!calculation of nonlinear stiffness and mass matrix
  template<int nnode>
  void b3_nlnstiffmass( Teuchos::ParameterList&   params,
                        std::vector<double>&      vel,
                        std::vector<double>&      disp,
                        Epetra_SerialDenseMatrix* stiffmatrix,
                        Epetra_SerialDenseMatrix* massmatrix,
                        Epetra_SerialDenseVector* force);

  //!calculation of elastic energy
  template<int nnode>
  void b3_energy( Teuchos::ParameterList&   params,
                  std::vector<double>&      disp,
                  Epetra_SerialDenseVector* intenergy);

  //! lump mass matrix
  template<int nnode>
  void lumpmass(Epetra_SerialDenseMatrix* emass);

  //!determine Gauss rule from required type of integration
  DRT::UTILS::GaussRule1D MyGaussRule(int nnode, IntegrationType integrationtype);

  //@}

  //! @name Methods for Brownian dynamics simulations

   //!compute damping matrix and statistical forces according to fluctuation dissipation theorem
   template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
   void CalcBrownian(Teuchos::ParameterList&    params,  //!<parameter list
                     const std::vector<double>& vel,  //!< element velocity vector
                     const std::vector<double>& disp, //!<element disp vector
                     Epetra_SerialDenseMatrix*  stiffmatrix,  //!< element stiffness matrix
                     Epetra_SerialDenseVector*  force);//!< element internal force vector

   //!computes damping coefficients for employed friction model
   void MyDampingConstants(Teuchos::ParameterList& params,
                           LINALG::Matrix<3,1>&    gamma);

   //!computes rotational damping forces and stiffness
   template<int nnode> //number of nodes
   void MyRotationalDamping(Teuchos::ParameterList&   params,  //!<parameter list
                           const std::vector<double>& vel,  //!< element velocity vector
                           const std::vector<double>& disp, //!<element disp vector
                           Epetra_SerialDenseMatrix*  stiffmatrix,  //!< element stiffness matrix
                           Epetra_SerialDenseVector*  force);//!< element internal force vector

   //!computes translational damping forces and stiffness
   template<int nnode, int ndim, int dof> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node
   void MyTranslationalDamping(Teuchos::ParameterList&  params,  //!<parameter list
                             const std::vector<double>& vel,  //!< element velocity vector
                             const std::vector<double>& disp, //!<element disp vector
                             Epetra_SerialDenseMatrix*  stiffmatrix,  //!< element stiffness matrix
                             Epetra_SerialDenseVector*  force);//!< element internal force vector

   //!computes stochastic translational forces and resulting stiffness
   template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
   void MyStochasticForces(Teuchos::ParameterList&  params,  //!<parameter list
                         const std::vector<double>& vel,  //!< element velocity vector
                         const std::vector<double>& disp, //!<element disp vector
                         Epetra_SerialDenseMatrix*  stiffmatrix,  //!< element stiffness matrix
                         Epetra_SerialDenseVector*  force);//!< element internal force vector

   //!computes stochastic rotational moments and resulting stiffness
   template<int nnode, int randompergauss> //number of nodes, number of random numbers required per Gauss point
   void MyStochasticMoments(Teuchos::ParameterList& params,  //!<parameter list
                         const std::vector<double>& vel,  //!< element velocity vector
                         const std::vector<double>& disp, //!<element disp vector
                         Epetra_SerialDenseMatrix*  stiffmatrix,  //!< element stiffness matrix
                         Epetra_SerialDenseVector*  force);//!< element internal force vector

   //!computes velocity of background fluid and gradient of that velocity at a certain evaluation point in the physical space and adds respective terms to internal forces and damping matrix
   template<int ndim> //number of dimensions of embedding space
   void MyBackgroundVelocity(Teuchos::ParameterList&       params,  //!<parameter list
                             const LINALG::Matrix<ndim,1>& evaluationpoint,  //!<point at which background velocity and its gradient has to be computed
                             LINALG::Matrix<ndim,1>&       velbackground,  //!< velocity of background fluid
                             LINALG::Matrix<ndim,ndim>&    velbackgroundgrad); //!<gradient of velocity of background fluid

   //!shifts nodes so that proper evaluation is possible even in case of periodic boundary conditions
   template<int nnode, int ndim> //number of nodes, number of dimensions
   void NodeShift(Teuchos::ParameterList& params,  //!<parameter list
                  std::vector<double>&    disp); //!<element disp vector

   //@}


  // don't want = operator
  Beam3& operator = (const Beam3& old);


}; // class Beam3



// << operator
std::ostream& operator << (std::ostream& os, const DRT::Element& ele);


} // namespace Elements
} // namespace DRT


#endif  // #ifndef BEAM_H
