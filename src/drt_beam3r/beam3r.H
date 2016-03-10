/*!----------------------------------------------------------------------
\file beam3ii.H

\brief 3D nonlinear Reissner beam element oy type II

3D nonlinear Reissner beam element oy type II (according to "The
interpolation of rotations and its application to finite element models
of geometrically exact rods", Romero 2004)

Attention: For this implementation, prescribed 3D rotation values have no direct physical interpretation so far (Within the
element, the rotation increments from time step to time step are treated as multiplicative updates!).
For 2D rotations, multiplicative and additive increments are identical and rotations can be prescribed without problems.

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>
*----------------------------------------------------------------------*/

//header file only included if not yet included!
#ifndef BEAM3II_H
#define BEAM3II_H

#include "Teuchos_RCP.hpp"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"

#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "Sacado.hpp"
#include "../headers/FAD_utils.H"

//#define BEAM3IIAUTOMATICDIFF  //default: off  //Decide whether automatic differentiation via Sacado is used or not for linearization of f_int

#ifdef BEAM3IIAUTOMATICDIFF
  typedef Sacado::Fad::DFad<double> FADordouble;

  #define MULTIPLICATIVEUPDATES  //default: on  // decide whether multiplicative updates of rotational DOFs are used (in this case, FAD-computed stiffmatrix needs to be transformed)
#else
  typedef double FADordouble;
#endif

//#define BEAM3IICONSTSTOCHFORCE  //Flag in order to hold stochastic forces constant over the element length
                                //and to only provide random numbers for the 3 translational DoFs
                                //(needed in order to compare with beam3eb)

namespace DRT
{

namespace ELEMENTS
{

class Beam3iiType : public DRT::ElementType
{
public:

  std::string Name() const { return "Beam3iiType"; }

  static Beam3iiType & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int         id,
                                             const int         owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual int Initialize(DRT::Discretization& dis);

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static Beam3iiType instance_;
};

/*!
\brief 3D nonlinear Reissner beam element implemented according to the following sources:
Jelenic, Crisfield, 1999, "Geometrically exact 3D beam theory: implementations of a strain-invariant finite element for statics and dynamics",
Crisfield, Jelenic, 1999, "Objectivity of strain measures in the geometrically exact three dimensional beam theory and its finite element implementation",
Romero, 2004, "The interpolation of rotations and its application to finite element models of geometrically exact rods",
Crisfield, 2003, "Non-linear Finite Element Analysis of Solids and Structures", Volume 2

author cyron (cyron @lnm.mw.tum.de)
*/
class Beam3ii : public DRT::Element
{
public:

  //! @name Friends
  friend class Beam3iiType;

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  */
  Beam3ii(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  */
  Beam3ii(const Beam3ii& old);

  /*!
  \brief Destructor
  */
  virtual ~Beam3ii();

  /*!
  \brief Deep copy this instance of Beam3ii and return pointer to the copy

  The Clone() method is used by the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed
.
  */
  virtual DRT::Element* Clone() const;

   /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const;

  /*!
  \brief Return unique ParObject id

  Every class implementing ParObject needs a unique id defined at the
  top of drt_parobject.H
  */
  virtual int UniqueParObjectId() const { return Beam3iiType::Instance().UniqueParObjectId(); }

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

  virtual DRT::ElementType & ElementType() const { return Beam3iiType::Instance(); }

  //@}

  /*!
  \brief Return number of lines to this element
  */
  virtual int NumLine() const {return 1;}

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
   \brief get current tangent at specified node i.e. Tcurr_
   */
   LINALG::Matrix<3,1>  Tcurr(const int NodeID);

  /*!
  \brief get reference tangent at specified node i.e. Tref_
   */
   LINALG::Matrix<3,1> Tref();

  /*!
  \brief Get jacobi_ factor of first gauss point
  */
  const double& GetJacobi() const { return jacobi_[0];}

  /*!
   \brief Get maximal bending curvature occurring in this element
   */
  const double& GetKappaMax() const
  {
    return kappa_max_;
  }

  /*!
  \brief Get current nodal triads as quaternions
  */
  const std::vector<LINALG::Matrix<4,1> >& Qnew() const {return Qnew_;}

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

  //!sets up from current nodal position all geometric parameters (considering current position as reference configuration)
  template<unsigned int nnode>
  void SetUpReferenceGeometry(const std::vector<double>& xrefe,
                              const std::vector<double>& rotrefe,
                              const bool secondinit = false)
  {
    /*this method initializes geometric variables of the element; the initialization can usually be applied to elements only once;
     *therefore after the first initialization the flag isinit is set to true and from then on this method does not take any action
     *when called again unless it is called on purpose with the additional parameter secondinit. If this parameter is passed into
     *the method and is true the element is initialized another time with xrefe;
     *note: the isinit_ flag is important for avoiding reinitialization upon restart. However, it should be possible to conduct a
     *second initialization in principle (e.g. for periodic boundary conditions*/

    for(unsigned int node=0; node<nnode; node++)
      for(int dof=0; dof<3 ; dof++)
      {
        if(fabs(rotrefe[3*node+dof])>1.0e-10)
          dserror("So far, rotrefe is only a dummy and cannot be used in Beam3ii::SetUpReferenceGeometry!");
      }

    if(!isinit_ || secondinit)
    {
      isinit_ = true;

      // ********************************* Initialize/resize variables **************************************

      // first the nodes for the reference triad \Lambda_r of the element are chosen according to eq. (6.2), Crisfield 1999;
      int nodeI = (int)floor(0.5*(NumNode()+1));
      int nodeJ = (int)floor(0.5*(NumNode()+2));

      // The node numbering applied in Crisfield 1999 differs from the order in which nodal quantities are stored in BACI.
      // Therefore we have to apply the following transformation:
      nodeI_=(unsigned int) LARGEROTATIONS::NumberingTrafo(nodeI,NumNode());
      nodeJ_=(unsigned int) LARGEROTATIONS::NumberingTrafo(nodeJ,NumNode());

      //resize and initialize STL vectors for rotational displacements so that they can store one value at each node
      dispthetaconv_.resize(nnode);
      dispthetaold_.resize(nnode);
      dispthetanew_.resize(nnode);
      for(unsigned int i=0; i<nnode; i++)
      {
        dispthetaconv_[i].Clear();
        dispthetaold_[i].Clear();
        dispthetanew_[i].Clear();
      }

      // resize vectors for Jacobi determinants so that they can store one value at each Gauss point
      // (nnode-1) values for under-integration and nnode values for exact integration
      jacobi_.resize(nnode-1);
      kapparef_.resize(nnode-1);
      gammaref_.resize(nnode-1);
      jacobimass_.resize(nnode);
      jacobinode_.resize(nnode);

      f_.Resize(6*nnode);
      f_.Scale(0.0);

      //Get DiscretizationType
      DRT::Element::DiscretizationType distype = Shape();

      // Note: index i refers to the nnode nodes
      //vector whose numgp-th element is a 1xnnode-matrix with all shape functions evaluated at the numgp-th GP
      std::vector<LINALG::Matrix<1,nnode> > I_i(nnode-1);

      //vector whose numgp-th element is a 1xnnode-matrix with the derivatives of all shape functions evaluated at the numgp-th GP
      std::vector<LINALG::Matrix<1,nnode> > I_i_xi(nnode-1);

      //derivative of initial curve in physical space with respect to curve parameter xi \in [-1;1] on element level
      LINALG::Matrix<3,1> drdxi;

      // dummy 3D vector
      LINALG::Matrix<3,1> dummy(true);


      // ****************************** Compute quantities valid for entire element *************************

      // reference quaternion Q_r corresponding to reference triad Lambda_r
      LINALG::Matrix<4,1> Q_r(true);

      // compute reference triad Lambda_r according to (3.9), Jelenic 1999
      // Qnew_ has already been filled with initial values from input file in ReadElement()
      CalcRefQuaternion<double>(Qnew_[nodeI_],Qnew_[nodeJ_],Q_r,dummy);


      // ********************************* Compute nodal quantities ******************************************

      // rotation angles between nodal triads and reference triad according to (3.8), Jelenic 1999
      // use type TMatrix here to enable use of functions Calc_Psi_l and Calc_Psi_l (defined for std::vector<LINALG::TMatrix>)
      std::vector<LINALG::TMatrix<double,3,1> > Psi_li(nnode);

      LINALG::Matrix<1,nnode> shapefuncderiv(true);

      for(unsigned int node=0; node<nnode; node++)
      {
        // compute nodal local rotations according to (3.8), Jelenic 1999
        CalcPsi_li<double>(Qnew_[node],Q_r,Psi_li[node]);

        // compute Jacobi determinant at gauss points for Lobatto quadrature (i.e. at nodes)
        // note: jacobinode_ only used in Statmech methods which are not yet functional for nnode>2 anyway
        // unfortunately, we can't check for Statmech application here (no parameter list available)
        // therefore, just do this for nnode==2 for the time being TODO use flag 'needstatmech_' here as well
        // TODO adapt for case nnode>2: position of Gauss-Lobatto quadrature points no longer at nodes for nnode>2
        if(nnode==2)
        {
          //Get position xi of nodes
          const double xi = -1.0 + 2* (double)node / (double)(nnode - 1);  // TODO this assumes equidistant nodes, right?

          //Get derivatives of shapefunctions at GP
          DRT::UTILS::shape_function_1D_deriv1(shapefuncderiv,xi,distype);

          drdxi.Clear();
          //calculate drdxi
          for(unsigned int node2=0; node2<nnode; node2++)
            for(int dof=0; dof<3; dof++)
              drdxi(dof) += shapefuncderiv(node2) * xrefe[3*node2+dof];

          //Store Jacobi determinant for each node (Jacobi determinant refers by definition always to the reference configuration)
          jacobinode_[node]= drdxi.Norm2();
        }
      }


      // *************************** Compute quantities at GPs for Gauss under-integration ******************
      // *********************************** => numgp = nnode-1 *********************************************

      // these quantities will later be used mainly for calculation of elasticity terms

      //interpolated local relative rotation \Psi^l at a certain Gauss point according to (3.11), Jelenic 1999
      LINALG::Matrix<3,1> Psi_l;
      /*derivative of interpolated local relative rotation \Psi^l with respect to arc-length parameter
       * at a certain Gauss point according to (3.11), Jelenic 1999*/
      LINALG::Matrix<3,1> Psi_l_s;
      // triad at GP
      LINALG::Matrix<3,3> Lambda;

      //Get the applied integrationpoints for underintegration
      DRT::UTILS::IntegrationPoints1D gausspoints_underint(MyGaussRule(nnode,gaussunderintegration));

      // evaluate all shape functions and derivatives with respect to element parameter xi at all specified Gauss points
      EvaluateShapeFunctions<nnode>(gausspoints_underint,I_i,I_i_xi,distype);

      //Loop through all GPs for under-integration and calculate jacobi determinants at the GPs
      for(int numgp=0; numgp < gausspoints_underint.nquad; numgp++)
      {
        drdxi.Clear();

        // calculate vector drdxi
        for(unsigned int node=0; node<nnode; node++)
          for(int dof=0; dof<3 ; dof++)
            drdxi(dof) += I_i_xi[numgp](node) * xrefe[3*node+dof];

        // Store Jacobi determinant with respect to reference configuration
        jacobi_[numgp]= drdxi.Norm2();

        // compute all quantities required for computestrain
        // recycle drdxi here because we need drds (derivative of beam centerline with respect to curve parameter in reference configuration)
        drdxi.Scale(1.0/jacobi_[numgp]);

        Calc_Psi_l<nnode,double>(Psi_li, I_i[numgp], Psi_l);
        Calc_Psi_l_s<nnode,double>(Psi_li, I_i_xi[numgp], jacobi_[numgp], Psi_l_s);
        Calc_Lambda<double>(Psi_l,Q_r,Lambda);

        // compute convected strains gamma and kappa according to Jelenic 1999, eq. (2.12) for reference configuration,
        // i.e. call this function with gammaref=zerovector and kapparef=zerovector
        dummy.Clear();
        computestrain<double>(Psi_l,Psi_l_s,drdxi,Lambda,dummy,dummy,gammaref_[numgp],kapparef_[numgp]);
      }


      // ***************************** Compute quantities at GPs for exact Gauss integration ****************************
      // ****************************************** => numgp = nnode ****************************************************

      // these quantities will later be used mainly for calculation of inertia terms -> named 'mass'
      Qconvmass_.resize(nnode);
      Qnewmass_.resize(nnode);
      wconvmass_.resize(nnode);
      wnewmass_.resize(nnode);
      aconvmass_.resize(nnode);
      anewmass_.resize(nnode);
      rttconvmass_.resize(nnode);
      rttnewmass_.resize(nnode);
      rttmodconvmass_.resize(nnode);
      rttmodnewmass_.resize(nnode);
      rtconvmass_.resize(nnode);
      rtnewmass_.resize(nnode);
      dispconvmass_.resize(nnode);
      dispnewmass_.resize(nnode);
      amodconvmass_.resize(nnode);
      amodnewmass_.resize(nnode);

      // reuse variables for individual shape functions and resize to new numgp=nnode
      I_i.resize(nnode);
      I_i_xi.resize(nnode);

      //Get the applied integrationpoints for exact integration of mass matrix
      DRT::UTILS::IntegrationPoints1D gausspoints_exactint(MyGaussRule(nnode,gaussexactintegration));

      // evaluate all shape functions and derivatives with respect to element parameter xi at all specified Gauss points
      EvaluateShapeFunctions<nnode>(gausspoints_exactint,I_i,I_i_xi,distype);

      //Loop through all GPs for exact integration and compute initial jacobi determinant
      for(int numgp=0; numgp<gausspoints_exactint.nquad; numgp++)
      {
        drdxi.Clear();

        // calculate vector drdxi
        for(unsigned int node=0; node<nnode; node++)
          for(int dof=0; dof<3 ; dof++)
            drdxi(dof) += I_i_xi[numgp](node) * xrefe[3*node+dof];

        //Store Jacobi determinant with respect to reference configuration
        jacobimass_[numgp]= drdxi.Norm2();

        // compute quaternion at this GP and store in Qnewmass_
        Calc_Psi_l<nnode,double>(Psi_li, I_i[numgp], Psi_l);
        Calc_Qgauss<double>(Psi_l,Q_r,Qnewmass_[numgp]);

        // copy Qnewmass_ to Qconvmass_
        Qconvmass_[numgp] = Qnewmass_[numgp];

        wconvmass_[numgp].Clear();
        wnewmass_[numgp].Clear();
        aconvmass_[numgp].Clear();
        anewmass_[numgp].Clear();
        amodconvmass_[numgp].Clear();
        amodnewmass_[numgp].Clear();
        rttconvmass_[numgp].Clear();
        rttnewmass_[numgp].Clear();
        rttmodconvmass_[numgp].Clear();
        rttmodnewmass_[numgp].Clear();
        rtconvmass_[numgp].Clear();
        rtnewmass_[numgp].Clear();
        dispconvmass_[numgp].Clear();
        dispnewmass_[numgp].Clear();
      }


      // **************************************** statmech-specific calculation ********************************
      // TODO not yet functional for nnode>2 => dserror in public methods Beam3ii::Tcurr() and Beam3ii::Tref()

      double norm2=0;
      // For 2 nodes the axial tangents are same for both direction
      Tref_.Clear();
      for(int dof = 0; dof< 3 ; dof++ )
        Tref_(dof) =  xrefe[3+dof] - xrefe[dof];

      norm2 = Tref_.Norm2();
      Tref_.Scale(1/norm2);

    }//if(!isinit_)

    return;
  }//DRT::ELEMENTS::Beam3ii::SetUpReferenceGeometry()

  //@}

  //! @name Construction

  /*!
  \brief Read input for this element
  */
  bool ReadElement(const std::string&          eletype,
                   const std::string&          distype,
                   DRT::INPUT::LineDefinition* linedef);

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
  virtual int EvaluateNeumann(Teuchos::ParameterList&   params,
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
  template<unsigned int nnode>
  void EvaluatePTC(Teuchos::ParameterList&   params,
                   Epetra_SerialDenseMatrix& elemat1);

  //@}


private:

  //! action parameters recognized by Beam3ii
  enum ActionType
  {
    calc_none,
    calc_struct_linstiff,
    calc_struct_nlnstiff,
    calc_struct_internalforce,
    calc_struct_linstiffmass,
    calc_struct_nlnstiffmass,
    calc_struct_nlnstifflmass,  //!< internal force, its stiffness and lumped mass matrix
    calc_struct_stress,
    calc_struct_eleload,
    calc_struct_fsiload,
    calc_struct_update_istep,
    calc_struct_reset_istep, //!< reset elementwise internal variables during iteration to last converged state
    calc_struct_ptcstiff, //calculate artificial stiffness due to PTC solution strategy
    calc_struct_energy //calculate current elastic energy of element
  };

  //! kind of integration to be performed
  enum IntegrationType
  {
    gaussexactintegration,
    gaussunderintegration,
    lobattointegration
  };

  //! variable saving whether element has already been initialized (then isinit_ == true)
  bool isinit_;
  //! variable storing whether Statmech functionalities are needed or not
  bool needstatmech_;
  //! storing epsilon(0)
  double eps_;
  //! vector storing the internal force vector
  Epetra_SerialDenseVector f_;
  //! vector storing internal forces at gauss point in material configuration
  LINALG::Matrix<3,1> Ngp_;
  //! norm of maximal bending curvature occurring in this element
  double kappa_max_;

  //! nodal triads in quaternion form at the end of the preceeding time step
  std::vector<LINALG::Matrix<4,1> > Qconv_;
  //! nodal triads in quaternion form at the end of the preceeding iteration step
  std::vector<LINALG::Matrix<4,1> > Qold_;
  //! nodal triads in quaternion during the current iteration step
  std::vector<LINALG::Matrix<4,1> > Qnew_;
  //! triads at Gauss points for exact integration in quaternion at the end of the preceeding time step (required for computation of angular velocity)
  std::vector<LINALG::Matrix<4,1> > Qconvmass_;
  //! current triads at Gauss points for exact integration in quaternion (required for computation of angular velocity)
  std::vector<LINALG::Matrix<4,1> > Qnewmass_;
  //! spatial angular velocity vector at Gauss points for exact integration at the end of the preceeding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > wconvmass_;
  //! current spatial angular velocity vector at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > wnewmass_;
  //! spatial angular acceleration vector at Gauss points for exact integration at the end of the preceeding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > aconvmass_;
  //! current spatial angular acceleration vector at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > anewmass_;
  //! translational acceleration vector at Gauss points for exact integration at the end of the preceeding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rttconvmass_;
  //! current translational acceleration vector at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rttnewmass_;
  //! modified translational acceleration vector at Gauss points for exact integration at the end of the preceeding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rttmodconvmass_;
  //! modified current translational acceleration vector at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rttmodnewmass_;
  //! translational velocity vector at Gauss points for exact integration at the end of the preceeding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rtconvmass_;
  //! current translational velocity vector at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rtnewmass_;
  //! translational displacement vector at Gauss points for exact integration at the end of the preceeding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > dispconvmass_;
  //! current translational displacement vector at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > dispnewmass_;
  //! modified spatial angular acceleration vector (according to gen-alpha time integration) at Gauss points for exact integration at the end of the preceeding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > amodconvmass_;
  //! current modified spatial angular acceleration vector (according to gen-alpha time integration) at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > amodnewmass_;
  //!rotation angle displacement at nodes; relates to the end of the preceeding time step
  std::vector<LINALG::Matrix<3,1> > dispthetaconv_;
  //!rotation angle displacement at nodes; relates to the end of the preceeding iteration step
  std::vector<LINALG::Matrix<3,1> > dispthetaold_;
  //!rotation angle displacement at nodes; relates to current configuration
  std::vector<LINALG::Matrix<3,1> > dispthetanew_;
  //!rotational pseudovectors at nodes in reference configuration
  std::vector<LINALG::Matrix<3,1> > theta0_;
  //! Matrix holding current tangent at the node connected to Beam3ii element for statmech problem
  LINALG::Matrix<3,1> Tcurr_;
  //! Matrix holding reference tangent at the node connected to Beam3ii element for statmech problem
  LINALG::Matrix<3,1> Tref_;

  //!initial curvature at Gauss points for elasticity (corresponding to \Lambda_0^t \Labmda'_0 in eq. (3.5), Crisfield 1999
  std::vector<LINALG::Matrix<3,1> > kapparef_;

  //!initial axial tension (allways zero) and shear deformation at Gauss points for elasticity (corresponding to \Lambda_0^t rprime_0 - (1,0,0) )
  std::vector<LINALG::Matrix<3,1> > gammaref_;

  //!node I for determination of reference triad, eq. (3.9), (3.10), Jelenic 1999
  unsigned int nodeI_;

  //!node J for determination of reference triad, eq. (3.9), (3.10), Jelenic 1999
  unsigned int nodeJ_;

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

  //! Vector holding value of Jacobi determinant jacobi for each Gauss point for underintegration
  std::vector<double> jacobi_;

  //! Vector holding value of Jacobi determinant jacobi for complete integration of massmatrix
  std::vector<double> jacobimass_;

  //! Vector holding value of Jacobi determinant jacobi at nodes
  std::vector<double> jacobinode_;

  //! kinetic energy of element
  double Ekin_;

  //! internal (elastic) energy of element
  double Eint_;

  //! angular momentum of the element
  LINALG::Matrix<3,1> L_;

  //! linear momentum of the element
  LINALG::Matrix<3,1> P_;

  //!kinetic energy from rotational dofs part1
  double kintorsionenergy_;

  //!kinetic energy from rotational dofs part2
  double kinbendingenergy_;

  //!kinetic energy from translational dofs
  double kintransenergy_;

  //! scale factor for translational inertia terms
  double inertscaletrans_;
  //! scale factor for tangential component of rotational inertia terms
  double inertscalerot1_;
  //! scale factor for remaining components of rotational inertia terms
  double inertscalerot2_;


  //! @name Internal calculation methods

  //!evaluate shape functions at gausspoints
  template<unsigned int nnode>
  void EvaluateShapeFunctions(const DRT::UTILS::IntegrationPoints1D&     gausspoints,
                              std::vector<LINALG::Matrix<1,nnode> >&     I_i,
                              const DRT::Element::DiscretizationType&    distype)       ///< DiscretizationType
  {
    if (I_i.size()!= (unsigned int) gausspoints.nquad)
      dserror("vector for %d individual shape functions has wrong size",gausspoints.nquad);

    for(int numgp=0; numgp < gausspoints.nquad; numgp++)
    {
      I_i[numgp].Clear();

      //Get location of GP in element parameter space xi \in [-1;1]
      const double xi = gausspoints.qxg[numgp][0];

      //evaluate shape functions at xi
      DRT::UTILS::shape_function_1D(I_i[numgp],xi,distype);
    }

    return;
  }

  //!evaluate shape functions and its derivatives at gausspoints
  template<unsigned int nnode>
  void EvaluateShapeFunctions(const DRT::UTILS::IntegrationPoints1D&     gausspoints,
                              std::vector<LINALG::Matrix<1,nnode> >&     I_i,
                              std::vector<LINALG::Matrix<1,nnode> >&     I_i_xi,
                              const DRT::Element::DiscretizationType&    distype)       ///< DiscretizationType
  {
    if (I_i.size()!= (unsigned int) gausspoints.nquad || I_i_xi.size()!= (unsigned int) gausspoints.nquad)
      dserror("vector for %d individual shape functions has wrong size",gausspoints.nquad);

    for(int numgp=0; numgp < gausspoints.nquad; numgp++)
    {
      I_i[numgp].Clear();
      I_i_xi[numgp].Clear();

      //Get location of GP in element parameter space xi \in [-1;1]
      const double xi = gausspoints.qxg[numgp][0];

      //evaluate shape functions and its derivatives at xi
      DRT::UTILS::shape_function_1D(I_i[numgp],xi,distype);
      DRT::UTILS::shape_function_1D_deriv1(I_i_xi[numgp],xi,distype);
    }

    return;
  }

  // compute quaternion corresponding to reference triad Lambda_r according to (3.9), Jelenic 1999
  template<typename T>
  void CalcRefQuaternion(const LINALG::TMatrix<T,4,1>& Q_nodeI,
                    const LINALG::TMatrix<T,4,1>& Q_nodeJ,
                    LINALG::TMatrix<T,4,1>& Q_r,
                    LINALG::TMatrix<T,3,1>& Phi_IJ)
  {
    Q_r.Clear();
    Phi_IJ.Clear();

    //quaternion of relative rotation between node I and J according to (3.10), Jelenic 1999
    LINALG::TMatrix<T,4,1> QIJ(true);
    //quaternion of half relative rotation between node I and J according to (3.9), Jelenic 1999
    LINALG::TMatrix<T,4,1> QIJhalf(true);

    //compute quaternion Qr corresponding to reference triad \Lambda_r according to eq. (3.10) and (3.9), Jelenic 1999
    LARGEROTATIONS::quaternionproduct<T>(Q_nodeJ,LARGEROTATIONS::inversequaternion<T>(Q_nodeI),QIJ);
    LARGEROTATIONS::quaterniontoangle<T>(QIJ,Phi_IJ);

    LINALG::TMatrix<T,3,1> Phi_IJhalf(Phi_IJ);
    Phi_IJhalf.Scale(0.5);
    LARGEROTATIONS::angletoquaternion<T>(Phi_IJhalf,QIJhalf);

    LARGEROTATIONS::quaternionproduct<T>(QIJhalf,Q_nodeI,Q_r);

    return;
  }

  // compute nodal local rotations according to (3.8), Jelenic 1999
  template<typename T>
  void CalcPsi_li(const LINALG::TMatrix<T,4,1>& Q_i,
                  const LINALG::TMatrix<T,4,1>& Q_r,
                  LINALG::TMatrix<T,3,1>& Psi_li)
  {
    // quaternion equivalent to Psi_li
    LINALG::TMatrix<T,4,1>  Q_li(true);

    LARGEROTATIONS::quaternionproduct<T>(Q_i,LARGEROTATIONS::inversequaternion<T>(Q_r),Q_li);
    LARGEROTATIONS::quaterniontoangle<T>(Q_li,Psi_li);

    return;
  }

  //!compute derivative of beam centerline with respect to arc-length parameter s
  template<unsigned int nnode, typename T>
  void Calc_r_s(const LINALG::TMatrix<T,6*nnode,1>& disp_totlag,
                const LINALG::Matrix<1,nnode>& deriv,
                const double& jacobi,
                LINALG::TMatrix<T,3,1>& r_s)
  {
    r_s.Clear();

    for (int dof=0; dof<3; ++dof)
      for (unsigned int node=0; node<nnode; ++node)
        r_s(dof) += deriv(node)*disp_totlag(6*node+dof);

    /*at this point we have computed derivative with respect to the element parameter \xi \in [-1;1];
     *as r' in (2.12) from (3.11) are derivatives with respect to the reference length,
     *we have to divide it by the Jacobi determinant at the respective point*/
    r_s.Scale(1.0/jacobi);

     return;
  }

  //!compute interpolated local relative rotation \Psi^l according to (3.11), Jelenic 1999
  template<unsigned int nnode, typename T>
  void Calc_Psi_l(const std::vector<LINALG::TMatrix<T,3,1> >& Psi_li,
                  const LINALG::Matrix<1,nnode>& func,
                  LINALG::TMatrix<T,3,1>& Psi_l)
  {
    Psi_l.Clear();

    for (int dof=0; dof<3; ++dof)
      for (unsigned int node=0; node<nnode; ++node)
        Psi_l(dof)+= func(node)*(Psi_li[node])(dof);

    return;
  }

  //!compute derivative of interpolated local relative rotation \Psi^l with respect to arc-length parameter s according to (3.11), Jelenic 1999
  template<unsigned int nnode, typename T>
  void Calc_Psi_l_s(const std::vector<LINALG::TMatrix<T,3,1> >& Psi_li,
                  const LINALG::Matrix<1,nnode>& deriv,
                  const double& jacobi,
                  LINALG::TMatrix<T,3,1>& Psi_l_s)
  {
    Psi_l_s.Clear();

    for (int dof=0; dof<3; ++dof)
      for (unsigned int node=0; node<nnode; ++node)
        Psi_l_s(dof)+= deriv(node)*(Psi_li[node])(dof);

    /*at this point we have computed derivative with respect to the element parameter \xi \in [-1;1];
     *as r' in (2.12) from (3.11) are derivatives with respect to the reference length,
     *we have to divide it by the Jacobi determinant at the respective point*/
    Psi_l_s.Scale(1.0/jacobi);

    return;
  }

  //!compute local triad \Lambda from Crisfield 1999, eq. (4.7)
  template <typename T>
  void Calc_Lambda(const LINALG::TMatrix<T,3,1>& Psi_l,
                  const LINALG::TMatrix<T,4,1>& Q_r,
                  LINALG::TMatrix<T,3,3>& Lambda)
  {
    Lambda.Clear();

    LINALG::TMatrix<T,4,1> Ql;
    LINALG::TMatrix<T,4,1> Qgauss;

    //compute relative rotation between triad at Gauss point and reference triad Qr
    LARGEROTATIONS::angletoquaternion(Psi_l,Ql);

    //compute rotation at Gauss point, i.e. the quaternion equivalent to \Lambda(s) in Crisfield 1999, eq. (4.7)
    LARGEROTATIONS::quaternionproduct(Ql,Q_r,Qgauss);

    //compute rotation matrix at Gauss point, i.e. \Lambda(s) in Crisfield 1999, eq. (4.7)
    LARGEROTATIONS::quaterniontotriad(Qgauss,Lambda);

    return;
  }

  //!compute quaternion equivalent to local triad \Lambda from Crisfield 1999, eq. (4.7)
  template <typename T>
  void Calc_Qgauss(const LINALG::TMatrix<T,3,1>& Psi_l,
                   const LINALG::TMatrix<T,4,1>& Q_r,
                   LINALG::TMatrix<T,4,1>& Qgauss)
  {
    Qgauss.Clear();

    LINALG::TMatrix<T,4,1> Ql;

    //compute relative rotation between triad at Gauss point and reference triad Qr
    LARGEROTATIONS::angletoquaternion(Psi_l,Ql);

    //compute rotation at Gauss point, i.e. the quaternion equivalent to \Lambda(s) in Crisfield 1999, eq. (4.7)
    LARGEROTATIONS::quaternionproduct(Ql,Q_r,Qgauss);

    return;
  }

  //!compute convected strain at certain Gauss point with triad rotmat according to Crisfield 1999, eq. (3.4) and eq. (4.9)
  template <typename T>
  void computestrain(const LINALG::TMatrix<T,3,1>& Psi_l,
                     const LINALG::TMatrix<T,3,1>& Psi_l_s,
                     const LINALG::TMatrix<T,3,1>& r_s,
                     const LINALG::TMatrix<T,3,3>& Lambda,
                     const LINALG::Matrix<3,1>& gammaref,
                     const LINALG::Matrix<3,1>& kapparef,
                     LINALG::TMatrix<T,3,1>& gamma,
                     LINALG::TMatrix<T,3,1>& kappa)
  {
    gamma.Clear();
    kappa.Clear();

    //convected strain gamma according to Crisfield 1999, eq. (3.4)
    gamma.MultiplyTN(Lambda,r_s);

    //In contrary to Crisfield 1999, eq. (3.4), the current implementation allows for initial values of the vector gammaref which
    //has also a second and a third component, i.e. it allows for initial shear deformation. This is the case, when the initial triad
    //at the evaluation point is not parallel to the centerline tangent vector at this point. The geometrically exact beam theory
    //does in general allow for such initial triads if they are considered consistently in the reference strains. While it is standard
    //to assume vanishing initial shear strains in the space-continuous setting, the possibility of initial shear strains might be advantageous
    //for the spatially discretized problem: For curved initial geometries, the nodal triad had to be determined such that the resulting
    //interpolated triad at the Gauss point would be tangential to the centerline tangent at this point resulting from the centerline interpolation.
    //In order to avoid this additional effort and to allow for an independent prescription of the nodal triads (e.g. prescribed by an analytical
    //geometry definition), the approach of considering arbitrary initial shear angles at the Gauss points is applied here.
    for (int i=0; i<3; ++i)
      gamma(i) -= gammaref(i);

    //Calculation of material curvature vector according to Crisfield 1999, eq. (4.2) (this equation has been derived
    //for a different beam element formulation but is also valid for the element type considered here),
    //or Jelenic 1999, paragraph on page 153 between NOTE 5 and NOTE 6
    LINALG::TMatrix<T,3,3> Tinv(true);
    Tinv=LARGEROTATIONS::Tinvmatrix<T>(Psi_l);
    //It is important to use the transposed matrix Tinv^T instead of Tinv (these two only differ in one of three terms)
    kappa.MultiplyTN(Tinv,Psi_l_s);

    //mechanically relevant curvature is current curvature minus curvature in reference position
    for (int i=0; i<3; ++i)
      kappa(i) -= kapparef(i);

     return;
  } // DRT::ELEMENTS::Beam3ii::computestrain

  //!computing convected stresses from convected strains vector and returning also convected constitutive matrix between both
  template <typename T>
  void strainstress(const LINALG::TMatrix<T,3,1>& gamma,
                    const LINALG::TMatrix<T,3,1>& kappa,
                    LINALG::TMatrix<T,3,1>& stressN,
                    LINALG::TMatrix<T,3,3>& CN,
                    LINALG::TMatrix<T,3,1>& stressM,
                    LINALG::TMatrix<T,3,3>& CM);

  //!push forward strains, stresses and constitutive matrices to their spatial counterparts by rotation matrix Lambda according to Romero 2004, eq. (3.10)
  template <typename T>
  void pushforward(const LINALG::TMatrix<T,3,3>& Lambda,
                    const LINALG::TMatrix<T,3,1>& stressN,
                    const LINALG::TMatrix<T,3,3>& CN,
                    const LINALG::TMatrix<T,3,1>& stressM ,
                    const LINALG::TMatrix<T,3,3>& CM,
                    LINALG::TMatrix<T,3,1>& stressn,
                    LINALG::TMatrix<T,3,3>& cn,
                    LINALG::TMatrix<T,3,1>& stressm,
                    LINALG::TMatrix<T,3,3>& cm);

  //!compute \tilde{I}^i in (3.18), page 152, Jelenic 1999, for all nodes i at a certain Gauss point
  /*----------------------------------------------------------------------------------------------------------------------*
   |compute \tilde{I}^i in (3.18), page 152, Jelenic 1999, for all nodes i at a certain Gauss point            cyron 04/10|
   *----------------------------------------------------------------------------------------------------------------------*/
  template<unsigned int nnode>
  void computeItilde(const LINALG::TMatrix<double,3,1>&               Psil,
                     std::vector<LINALG::TMatrix<double,3,3> >&       Itilde,
                     const LINALG::TMatrix<double,3,1>&               phiIJ,
                     const LINALG::TMatrix<double,3,3>&               Lambdar,
                     const std::vector<LINALG::TMatrix<double,3,1> >& Psili,
                     const LINALG::Matrix<1,nnode>&                   funct)
  {
    // note: Itilde should always be of type double because it is only needed for computation of analytical stiffness/mass matrix

    //auxiliary matrices for storing intermediate results
    LINALG::TMatrix<double,3,3> auxmatrix;
    auxmatrix.PutScalar(0);
    LINALG::TMatrix<double,3,3> auxmatrix2;
    auxmatrix2.PutScalar(0);
    LINALG::TMatrix<double,3,3> Tinv_Psil = LARGEROTATIONS::Tinvmatrix<double>(Psil);

    //make sure that Itilde has proper dimensions
    Itilde.resize(nnode);

    //compute squared brackets term in (3.18), Jelenic 1999
    LINALG::TMatrix<double,3,3> squaredbrackets;
    squaredbrackets.PutScalar(0);
    for(unsigned int node=0; node<nnode; ++node)
    {
      auxmatrix.Clear();

      auxmatrix = LARGEROTATIONS::Tmatrix<double>(Psili[node]);
      auxmatrix.Scale(funct(node));
      auxmatrix2.Update(-1.0,auxmatrix,1.0);
    }

    squaredbrackets.Multiply(Tinv_Psil,auxmatrix2);

    for(int i=0; i<3; i++)
      squaredbrackets(i,i) += 1;

    //loop through all nodes i
    for (unsigned int node=0; node<nnode; ++node)
    {
      //compute rightmost term in curley brackets in (3.18), Jelenic 1999
      Itilde[node].Clear();
      Itilde[node].Multiply(Tinv_Psil,LARGEROTATIONS::Tmatrix<double>(Psili[node]));
      Itilde[node].Scale(funct(node));

      //if node i is node I then add squared brackets term times v_I
      if(node == nodeI_)
      {
        auxmatrix.Multiply(squaredbrackets,vI(phiIJ));
        Itilde[node] += auxmatrix;
      }

      //if node i is node J then add squared brackets term times v_J
      if(node == nodeJ_)
      {
        auxmatrix.Multiply(squaredbrackets,vJ(phiIJ));
        Itilde[node] += auxmatrix;
      }

      //now the term in the curly brackets has been computed and has to be rotated by \Lambda_r and \Lambda_r^t
      auxmatrix.MultiplyNT(Itilde[node],Lambdar);
      Itilde[node].MultiplyNN(Lambdar,auxmatrix);
    }

     return;
  } // DRT::ELEMENTS::Beam3ii::computeItilde


  /*----------------------------------------------------------------------------------------------------------------------*
   | compute \tilde{I}^i if pre-calculated values are of type FAD                                              grill 02/16|
   *----------------------------------------------------------------------------------------------------------------------*/
  template<unsigned int nnode>
  void computeItilde(const LINALG::TMatrix<FAD,3,1>&               Psil,
                     std::vector<LINALG::TMatrix<double,3,3> >&    Itilde,
                     const LINALG::TMatrix<FAD,3,1>&               phiIJ,
                     const LINALG::TMatrix<FAD,3,3>&               Lambdar,
                     const std::vector<LINALG::TMatrix<FAD,3,1> >& Psili,
                     const LINALG::Matrix<1,nnode>&                funct)
  {
    // note: Itilde should always be of type double because it is only needed for computation of analytical stiffness/mass matrix
    // therefore no templated function used here; cast pre-calculated FAD values to double and call 'regular' function afterwards

    LINALG::TMatrix<double,3,1> Psil_double = FADUTILS::CastToDouble<FAD,3,1>(Psil);
    LINALG::TMatrix<double,3,1> phiIJ_double = FADUTILS::CastToDouble<FAD,3,1>(phiIJ);
    LINALG::TMatrix<double,3,3> Lambdar_double = FADUTILS::CastToDouble<FAD,3,3>(Lambdar);

    std::vector<LINALG::TMatrix<double,3,1> > Psili_double;
    for (unsigned int i=0; i<Psili.size(); ++i)
      Psili_double.push_back(FADUTILS::CastToDouble<FAD,3,1>(Psili[i]));

    computeItilde(Psil_double,Itilde,phiIJ_double,Lambdar_double,Psili_double,funct);

     return;
  } // DRT::ELEMENTS::Beam3ii::computeItilde

  /*----------------------------------------------------------------------------------------------------------------------*
   |compute \tilde{I}^{i'} in (3.19), page 152, Jelenic 1999 for all nodes i at a certain Gauss point          cyron 04/10|
   *----------------------------------------------------------------------------------------------------------------------*/
  template<unsigned int nnode, typename T>
  void computeItildeprime(const LINALG::TMatrix<T,3,1>&               Psil,
                          const LINALG::TMatrix<T,3,1>&               Psilprime,
                          std::vector<LINALG::TMatrix<T,3,3> >&       Itildeprime,
                          const LINALG::TMatrix<T,3,1>&               phiIJ,
                          const LINALG::TMatrix<T,3,3>&               Lambdar,
                          const std::vector<LINALG::TMatrix<T,3,1> >& Psili,
                          const LINALG::Matrix<1,nnode>&              funct,
                          const LINALG::Matrix<1,nnode>&              deriv,
                          const double&                               jacobi)
  {
    //auxiliary matrices for storing intermediate results
    LINALG::TMatrix<T,3,3> auxmatrix;

    //make sure that Itildeprime has proper dimensions
    Itildeprime.resize(nnode);

    //matrix d(T^{-1})/dx
    LINALG::TMatrix<T,3,3> dTinvdx;
    LARGEROTATIONS::computedTinvdx<T>(Psil,Psilprime,dTinvdx);

    //compute T^{~} according to remark subsequent to (3.19), Jelenic 1999
    LINALG::TMatrix<T,3,3> Ttilde;
    Ttilde.PutScalar(0);
    for(unsigned int node=0; node<nnode; ++node)
    {
      auxmatrix = LARGEROTATIONS::Tmatrix<T>(Psili[node]);
      auxmatrix.Scale(funct(node));
      Ttilde += auxmatrix;
    }

    //compute T^{~'} according to remark subsequent to (3.19), Jelenic 1999
    LINALG::TMatrix<T,3,3> Ttildeprime;
    Ttildeprime.PutScalar(0);
    for(unsigned int node=0; node<nnode; ++node)
    {
      auxmatrix = LARGEROTATIONS::Tmatrix<T>(Psili[node]);
      auxmatrix.Scale(deriv(node)/jacobi);
      Ttildeprime += auxmatrix;
    }

    //compute first squared brackets term in (3.18), Jelenic 1999
    LINALG::TMatrix<T,3,3> squaredbrackets;
    squaredbrackets.PutScalar(0);
    squaredbrackets.Multiply(dTinvdx,Ttilde);
    auxmatrix.Multiply(LARGEROTATIONS::Tinvmatrix(Psil),Ttildeprime);
    squaredbrackets += auxmatrix;

    //loop through all nodes i
    for (unsigned int node=0; node<nnode; ++node)
    {

      //compute first term in second squared brackets
      Itildeprime[node] = dTinvdx;
      Itildeprime[node].Scale(funct(node));

      //compute second term in second squared brackets
      auxmatrix.PutScalar(0);
      auxmatrix += LARGEROTATIONS::Tinvmatrix(Psil);
      auxmatrix.Scale(deriv(node)/jacobi);

      //compute second squared brackets
      auxmatrix += Itildeprime[node];

      //compute second squared brackets time T(\Psi^l_j)
      Itildeprime[node].Multiply(auxmatrix,LARGEROTATIONS::Tmatrix<T>(Psili[node]));

      //if node i is node I then add first squared brackets term times v_I
      if(node == nodeI_)
      {
        auxmatrix.Multiply(squaredbrackets,vI(phiIJ));
        Itildeprime[node] -= auxmatrix;
      }

      //if node i is node J then add first squared brackets term times v_J
      if(node == nodeJ_)
      {
        auxmatrix.Multiply(squaredbrackets,vJ(phiIJ));
        Itildeprime[node] -= auxmatrix;
      }

      //now the term in the curly brackets has been computed and has to be rotated by \Lambda_r and \Lambda_r^t
      auxmatrix.MultiplyNT(Itildeprime[node],Lambdar);
      Itildeprime[node].MultiplyNN(Lambdar,auxmatrix);
    }

     return;
  } // DRT::ELEMENTS::Beam3ii::computeItildeprime

  //!compute matrix v_I as outlined in the equations above (3.15) on page 152 of Jelenic 1999
  /*----------------------------------------------------------------------------------------------------------------------*
   |compute matrix v_I as outlined in the equations above (3.15) on page 152 of Jelenic 1999                   cyron 04/10|
   *----------------------------------------------------------------------------------------------------------------------*/
  LINALG::TMatrix<double,3,3> vI(const LINALG::TMatrix<double,3,1>& phiIJ)
  {
    //matrix v_I
    LINALG::TMatrix<double,3,3> result;

    LARGEROTATIONS::computespin(result,phiIJ);
    if(phiIJ.Norm2() == 0)
      result.Scale(0.25);
    else
      result.Scale(tan(phiIJ.Norm2()/4.0)/phiIJ.Norm2());

    for(int i=0; i<3; i++)
      result(i,i) +=1;

    result.Scale(0.5);

    return result;
  } // DRT::ELEMENTS::Beam3ii::vI

   /*----------------------------------------------------------------------------------------------------------------------*
   |compute matrix v_J as outlined in the equations above (3.15) on page 152 of Jelenic 1999                   cyron 04/10|
   *----------------------------------------------------------------------------------------------------------------------*/
  LINALG::TMatrix<double,3,3> vJ(const LINALG::TMatrix<double,3,1>& phiIJ)
  {
    //matrix v_J
    LINALG::TMatrix<double,3,3> result;

    LARGEROTATIONS::computespin(result,phiIJ);
    if(phiIJ.Norm2() == 0)
      result.Scale(-0.25);
    else
      result.Scale(-tan(phiIJ.Norm2()/4.0)/phiIJ.Norm2());

    for(int i=0; i<3; i++)
      result(i,i) +=1;

    result.Scale(0.5);

    return result;
  } // DRT::ELEMENTS::Beam3ii::vJ

  //!calculation of nonlinear stiffness and mass matrix
  template<unsigned int nnode>
  void b3_nlnstiffmass( Teuchos::ParameterList& params,
                        std::vector<double>&      acc,
                        std::vector<double>&      vel,
                        std::vector<double>&      disp,
                        Epetra_SerialDenseMatrix* stiffmatrix,
                        Epetra_SerialDenseMatrix* massmatrix,
                        Epetra_SerialDenseVector* force,
                        Epetra_SerialDenseVector* inertia_force);

  //!calculation of elastic energy
  template<unsigned int nnode>
  void b3_energy(Teuchos::ParameterList& params,
                  std::vector<double>&      disp,
                  Epetra_SerialDenseVector* intenergy);

  //! lump mass matrix
  template<unsigned int nnode>
  void lumpmass(Epetra_SerialDenseMatrix* emass);

  //!determine Gauss rule from required type of integration
  DRT::UTILS::GaussRule1D MyGaussRule(unsigned int nnode, IntegrationType integrationtype);

  //@}

  //! @name Methods for Brownian dynamics simulations

   //!compute damping matrix and statistical forces according to fluctuation dissipation theorem
   template<unsigned int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
   void CalcBrownian(Teuchos::ParameterList&                          params,  //!<parameter list
                     const std::vector<double>&                       vel,  //!< element velocity vector
                     const std::vector<double>&                       disp, //!<element disp vector
                     Epetra_SerialDenseMatrix*                        stiffmatrix,  //!< element stiffness matrix
                     Epetra_SerialDenseVector*                        force,
                     std::vector<LINALG::Matrix<1,nnode> >&           Imass,
                     std::vector<std::vector<LINALG::Matrix<3,3> > >& Itildemass);//!< element internal force vector

   //!computes damping coefficients for employed friction model
   void MyDampingConstants(Teuchos::ParameterList& params,
                           LINALG::Matrix<3,1>&    gamma);

   //!computes rotational damping forces and stiffness
   template<unsigned int nnode> //number of nodes
   void MyRotationalDamping(Teuchos::ParameterList& params,  //!<parameter list
                           const std::vector<double>&                             vel,  //!< element velocity vector
                           const std::vector<double>&                             disp, //!<element disp vector
                           Epetra_SerialDenseMatrix*                              stiffmatrix,//!< element stiffness matrix
                           Epetra_SerialDenseVector*                              force,      //!< element internal force vector
                           const DRT::UTILS::IntegrationPoints1D&                 gausspointsdamping,
                           const std::vector<LINALG::Matrix<1,nnode> >&           Idamping,
                           const std::vector<std::vector<LINALG::Matrix<3,3> > >& Itildedamping,
                           const std::vector<LINALG::Matrix<4,1> >&               Qconvdamping,
                           const std::vector<LINALG::Matrix<4,1> >&               Qnewdamping);

   //!computes translational damping forces and stiffness
   template<unsigned int nnode, int ndim, int dof> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node
   void MyTranslationalDamping(Teuchos::ParameterList&  params,  //!<parameter list
                             const std::vector<double>& vel,  //!< element velocity vector
                             const std::vector<double>& disp, //!<element disp vector
                             Epetra_SerialDenseMatrix*  stiffmatrix,  //!< element stiffness matrix
                             Epetra_SerialDenseVector*  force);//!< element internal force vector

   //!computes stochastic translational forces and resulting stiffness
   template<unsigned int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
   void MyStochasticForces(Teuchos::ParameterList&  params,  //!<parameter list
                         const std::vector<double>& vel,  //!< element velocity vector
                         const std::vector<double>& disp, //!<element disp vector
                         Epetra_SerialDenseMatrix*  stiffmatrix,  //!< element stiffness matrix
                         Epetra_SerialDenseVector*  force);//!< element internal force vector

   //!computes stochastic rotational moments and resulting stiffness
   template<unsigned int nnode, int randompergauss> //number of nodes, number of random numbers required per Gauss point
   void MyStochasticMoments(Teuchos::ParameterList&                             params,  //!<parameter list
                         const std::vector<double>&                             vel,  //!< element velocity vector
                         const std::vector<double>&                             disp, //!<element disp vector
                         Epetra_SerialDenseMatrix*                              stiffmatrix,//!< element stiffness matrix
                         Epetra_SerialDenseVector*                              force,      //!< element internal force vector
                         const DRT::UTILS::IntegrationPoints1D&                 gausspointsdamping,
                         const std::vector<LINALG::Matrix<1,nnode> >&           Idamping,
                         const std::vector<std::vector<LINALG::Matrix<3,3> > >& Itildedamping,
                         const std::vector<LINALG::Matrix<4,1> >&               Qconvdamping,
                         const std::vector<LINALG::Matrix<4,1> >&               Qnewdamping);

   //!computes velocity of background fluid and gradient of that velocity at a certain evaluation point in the physical space and adds respective terms to internal forces and damping matrix
   template<int ndim> //number of dimensions of embedding space
   void MyBackgroundVelocity(Teuchos::ParameterList&       params,  //!<parameter list
                             const LINALG::Matrix<ndim,1>& evaluationpoint,  //!<point at which background velocity and its gradient has to be computed
                             LINALG::Matrix<ndim,1>&       velbackground,  //!< velocity of background fluid
                             LINALG::Matrix<ndim,ndim>&    velbackgroundgrad); //!<gradient of velocity of background fluid

   //!shifts nodes so that proper evaluation is possible even in case of periodic boundary conditions
   template<unsigned int nnode, int ndim> //number of nodes, number of dimensions
   void NodeShift(Teuchos::ParameterList&      params,  //!<parameter list
                 std::vector<double>& disp); //!<element disp vector

   //@}


  // don't want = operator
  Beam3ii& operator = (const Beam3ii& old);


}; // class Beam3ii


// << operator
std::ostream& operator << (std::ostream& os, const DRT::Element& ele);


} // namespace Elements
} // namespace DRT


#endif  // #ifndef BEAM_H
