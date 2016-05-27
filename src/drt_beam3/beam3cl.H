/*----------------------------------------------------------------------------*/
/*!
\file beam3cl.H

\brief 3D nonlinear Reissner beam element of type II with interpolated node positions

\level 3

\maintainer Christoph Meier
*/
/*----------------------------------------------------------------------------*/

//header file only included if not yet included!
#ifndef BEAMCL_H
#define BEAMCL_H


#include "../drt_lib/drt_node.H"
#include "../drt_inpar/inpar_statmech.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_integration.H"

#include "beam3_base.H"

//! \brief DRT: namespace of the ccarat discretization module
namespace DRT
{

  namespace ELEMENTS
  {
    class BeamCLType: public DRT::ElementType
    {
      public:

        std::string Name() const
        {
          return "BeamCLType";
        }

        static BeamCLType & Instance();

        virtual DRT::ParObject* Create(const std::vector<char> & data);

        virtual Teuchos::RCP<DRT::Element> Create(const std::string eletype,
                                                  const std::string eledistype,
                                                  const int id,
                                                  const int owner);

        virtual Teuchos::RCP<DRT::Element> Create(const int id,
                                                  const int owner);

        virtual int Initialize(DRT::Discretization& dis);

        virtual void NodalBlockInformation(DRT::Element * dwele, int & numdf,
                                           int & dimns, int & nv, int & np);

        virtual void ComputeNullSpace(DRT::Discretization & dis,
                                      std::vector<double> & ns,
                                      const double * x0, int numdf, int dimns);

        virtual void SetupElementDefinition(
            std::map<std::string,
                std::map<std::string, DRT::INPUT::LineDefinition> > & definitions);

      private:

        static BeamCLType instance_;
    };

    /*!
     \brief 3D nonlinear Reissner beam element implemented according to the following sources:
     Jelenic, Crisfield, 1999, "Geometrically exact 3D beam theory: implementations of a strain-invariant finite element for statics and dynamics",
     Crisfield, Jelenic, 1998, "Objectivity of strain measures in the geometrically exact three dimensional beam theory and its finite element implementation",
     Romero, 2004, "The interpolation of rotations and its application to finite element models of geometrically exact rods",
     Crisfield, 2003, "Non-linear Finite Element Analysis of Solids and Structures", Volume 2

     The element has been implemented with linear basis functions only so far.

     author cyron (cyron @lnm.mw.tum.de)
     */
    class BeamCL: public Beam3Base
    {
      public:

        //! @name Friends
        friend class BeamCLType;

        //! @name Constructors and destructors and related methods

        /*!
         \brief Standard Constructor

         \param id    (in): A globally unique element id
         \param etype (in): Type of element
         \param owner (in): owner processor of the element
         */
        BeamCL(int id, int owner);

        /*!
         \brief Copy Constructor

         Makes a deep copy of a Element
         */
        BeamCL(const BeamCL& old);

        /*!
         \brief Destructor
         */
        virtual ~BeamCL();

        /*!
         \brief Deep copy this instance of Beam3r and return pointer to the copy

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
        virtual int UniqueParObjectId() const
        {
          return BeamCLType::Instance().UniqueParObjectId();
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

        virtual DRT::ElementType & ElementType() const
        {
          return BeamCLType::Instance();
        }

        //@}

        /*!
         \brief Return number of lines to this element
         */
        virtual int NumLine() const
        {
          return 1;
        }

        /** \brief get number of nodes used for centerline interpolation
         *
         *  \author grill
         *  \date 05/16 */
        virtual inline int NumCenterlineNodes() const { return 2; }

        /*!
         \brief Get vector of RefCountPtrs to the lines of this element
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
        virtual int NumDofPerElement() const
        {
          return 0;
        }

        /*!
         \brief Print this element
         */
        virtual void Print(std::ostream& os) const;
        /*!
         \brief Get moment of inertia of area Iyy_ of this element
         */
        const double& Iyy() const
        {
          return Iyy_;
        }

        /*!
         \brief Get moment of inertia of area Izz_ of this element
         */
        const double& Izz() const
        {
          return Izz_;
        }

        /*!
         \brief Get polar moment of inertia of area Irr_ of this element
         */
        const double& Irr() const
        {
          return Irr_;
        }

        /** \brief get access to the reference length
         *
         *  \author grill
         *  \date 05/16 */
        virtual inline double GetRefLength() const { return 2.0*jacobi_[0]; }

        /*!
         \brief Get crosssection crosssec_ of this element
         */
        const double& Crosssec() const
        {
          return crosssec_;
        }

        /*!
         \brief Get crosssecshear_ of this element
         */
        const double& Crosssecshear() const
        {
          return crosssecshear_;
        }

        const std::vector<double>& MyBindingPosition() const
        {
          return mybindingposition_;
        }

        /** \brief return flag whether Hermite polynomials are applied for centerline interpolation
         *
         *  \author grill
         *  \date 05/16 */
        virtual inline bool HermiteCenterlineInterpolation() const { return false; }

        //@}

        //! @name Construction

        /*!
         \brief Read input for this element

         This class implements a dummy of this method that prints a warning and
         returns false. A derived class would read one line from the input file and
         store all necessary information.

         */
        //virtual bool ReadElement();
        /*!
         \brief Read input for this element
         */
        bool ReadElement(const std::string& eletype, const std::string& distype,
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
        virtual int Evaluate(Teuchos::ParameterList& params,
                             DRT::Discretization& discretization,
                             std::vector<int>& lm,
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
                                    DRT::Discretization& discretization,
                                    DRT::Condition& condition,
                                    std::vector<int>& lm,
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
        template<int nnode>
        void EvaluatePTC(Teuchos::ParameterList& params,
                         Epetra_SerialDenseMatrix& elemat1);

        //!computes the number of different random numbers required in each time step for generation of stochastic forces
        int HowManyRandomNumbersINeed();

        //! \brief internal force vector
        const Epetra_SerialDenseVector& InternalForceVector() const{return f_;}
        //! \brief return material forces at Gauss point
        const LINALG::Matrix<3,1>& MatForceGp() const {return Ngp_;}

        //! brief epsilonnormal
        const double& EpsilonSgn() const {return eps_;}

        //@}

        //! \brief return reference Coordinates of interpolated nodes
        const std::vector<double>& XRef() const {return xrefe_;}

        //! \brief return nodal reference rotations
        const std::vector<double>& RotRef() const {return rotrefe_;}

        //! \brief Get current fictitious nodal triads as quaternions
        const std::vector<LINALG::Matrix<4, 1> >& Qnew() const {return Qnew_;}

        //! \brief Get rotation quaternion previously set by SetRotation()
        const LINALG::Matrix<4,1>& Qrot() const {return Qrot_;}

        //! \brief get length scaling factor
        const double& LScaleFac() const {return lscalefac_;}

        //!\brief Set moment of inertia of area Iyy_ of this element
        void SetIyy(const double& Iyy);

        //! \brief Set moment of inertia of area Iyy_ of this element
        void SetIzz(const double& Izz);

        // !\brief Set moment of inertia of area Iyy_ of this element
        void SetIrr(const double& Irr);

        // !\brief Set cross section area of this element
        void SetCrossSec(const double& crosssec);

        // !\brief Set cross section with shear correction factor of this element
        void SetCrossSecShear(const double& crosssecshear);

        // !\brief Set (Internodal) Binding Position of this element
        void SetBindingPosition(const double& pos1, const double& pos2);

        //!sets up from current nodal position all geometric parameters (considering current position as reference configuration)
        template<int fnnode>
        void SetUpReferenceGeometry(const std::vector<double>& xrefe,
                                    const std::vector<double>& rotrefe,
                                    const bool secondinit = false)
        {
          /*this method initializes geometric variables of the element; the initilization can usually be applied to elements only once;
           *therefore after the first initilization the flag isinit is set to true and from then on this method does not take any action
           *when called again unless it is called on purpose with the additional parameter secondinit. If this parameter is passed into
           *the method and is true the element is initialized another time with respective xrefe and rotrefe;
           *note: the isinit_ flag is important for avoiding reinitialization upon restart. However, it should be possible to conduct a
           *second initilization in principle (e.g. for periodic boundary conditions*/

          if (!isinit_ || secondinit)
          {
            isinit_ = true;
            int nnode = 2 * fnnode;

            f_.Resize(6 * fnnode);
            f_.Scale(0.0);

            /*first the nodes for the reference triad \Lambda_r of the element are chosen according to eq. (6.2), Crisfield 1999;
             *note that the first node of the element in BACI is node 0 so that we need -1 in the end to convert from the notation
             *in Crisfield 1999 to the BACI convention*/
            nodeI_ = (int) floor(0.5 * (0.5 * NumNode() + 1)) - 1;
            nodeJ_ = (int) floor(0.5 * (0.5 * NumNode() + 2)) - 1;

            //resize and initialized STL vectors for interpolated nodes for rotational displacements so that they can store one value at each node
            dispthetaconv_.resize(fnnode);
            dispthetaold_.resize(fnnode);
            dispthetanew_.resize(fnnode);

            //set nodal triads according to input file
            rQnew_.resize(4);
            rQold_.resize(4);
            rQconv_.resize(4);
            Qnew_.resize(2);
            Qold_.resize(2);
            Qconv_.resize(2);

//            //resize and initialized STL vectors for interpolated nodes for rotational displacements so that they can store one value at each node
//            dispthetaconv_.clear();
//            dispthetaold_.clear();
//            dispthetanew_.clear();
//            for (int i = 0; i < fnnode; i++)
//            {
//              dispthetaconv_.push_back(LINALG::Matrix<3, 1>(true));
//              dispthetaold_.push_back(LINALG::Matrix<3, 1>(true));
//              dispthetanew_.push_back(LINALG::Matrix<3, 1>(true));
//            }
//
//            //set nodal tridas according to input file
//            rQnew_.clear();
//            rQold_.clear();
//            rQconv_.clear();
//            for (int i = 0; i < 4; i++)
//            {
//              rQnew_.push_back(LINALG::Matrix<4, 1>(true));
//              rQold_.push_back(LINALG::Matrix<4, 1>(true));
//              rQconv_.push_back(LINALG::Matrix<4, 1>(true));
//            }
//
//            Qnew_.clear();
//            Qold_.clear();
//            Qconv_.clear();
//            for (int i = 0; i < 2; i++)
//            {
//              Qnew_.push_back(LINALG::Matrix<4, 1>(true));
//              Qold_.push_back(LINALG::Matrix<4, 1>(true));
//              Qconv_.push_back(LINALG::Matrix<4, 1>(true));
//            }

            for (int i = 0; i < fnnode; i++)
              for (int j = 0; j < 3; j++)
              {
                dispthetaconv_[i](j) = 0;
                dispthetaold_[i](j) = 0;
                dispthetanew_[i](j) = 0;
              }
            if (!secondinit)
            {
              rdispthetaconv_.resize(nnode);
              rdispthetaold_.resize(nnode);
              rdispthetanew_.resize(nnode);
              for (int i = 0; i < nnode; i++)
                for (int j = 0; j < 3; j++)
                {
                  rdispthetaconv_[i](j) = 0;
                  rdispthetaold_[i](j) = 0;
                  rdispthetanew_[i](j) = 0;
                }
            }

            // When BeamCL is initialized as a crosslinker, it additionally needs all 12 Rotational DOFs
            // Do not overwrite rotational displacements in case of a second initialization (e.g. due to breakage by Periodic Boundary)
            xrefe_.resize(3*fnnode);
            xrefe_ = xrefe;

            if (rotrefe.size() > 6 && !secondinit)
            {
              rotrefe_.resize(3*2*fnnode);
              rotrefe_ = rotrefe;
              for (int i = 0; i < nnode; i++)
                for (int j = 0; j < 3; j++)
                {
                  rdispthetaconv_[i](j) = rotrefe[3 * i + j];
                  rdispthetaold_[i](j) = rotrefe[3 * i + j];
                  rdispthetanew_[i](j) = rotrefe[3 * i + j];
                }
            }

            //resize STL vectors for Jacobi determinants so that they can store one value at each Gauss point
            jacobi_.resize(fnnode - 1);
            kapparef_.resize(fnnode - 1);
            jacobimass_.resize(fnnode);
            jacobinode_.resize(fnnode);

            //create Matrix for the derivates of the shapefunctions at the GP
            LINALG::Matrix<1, fnnode> shapefuncderiv;

            //create Matrix for the shapefunctions at the GP
            LINALG::Matrix<1, fnnode> funct;

            //derivative of curve in physical space with respect to curve parameter xi \in [-1;1] on element level
            LINALG::Matrix<3, 1> drdxi;

            //triad at nodes in reference configuration
            LINALG::Matrix<3, 3> Tref;

            //Get DiscretizationType
            DRT::Element::DiscretizationType distype = Shape();

            //Get the applied integrationpoints for underintegration
            DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(fnnode, gaussunderintegration));

            //Loop through all GPs for underintegration and calculate jacobi determinants and initial curvature at the GPs
            for (int numgp = 0; numgp < gausspoints.nquad; numgp++)
            {
              //Get position xi of GP
              const double xi = gausspoints.qxg[numgp][0];

              //Get derivatives of shapefunctions at GP
              DRT::UTILS::shape_function_1D_deriv1(shapefuncderiv, xi, distype);

              //Get shapefunctions at GP
              DRT::UTILS::shape_function_1D(funct, xi, distype);

              drdxi.Clear();

              //calculate vector dxdxi
              for (int node = 0; node < fnnode; node++)
                for (int dof = 0; dof < 3; dof++)
                  drdxi(dof) += shapefuncderiv(node) * xrefe[3 * node + dof];

              //Store Jacobi determinant with respect to reference configuration
              jacobi_[numgp] = drdxi.Norm2();

              // Calculate appropriate Triads for Reference Configuration
              for (int k = 0; k < 3; k++)
              {
                //t1 axis points in positive direction along xi and is a unit vector
                Tref(k, 0) = drdxi(k) / jacobi_[numgp];
              }

              //t2 is a unit vector in the x2x3-plane orthogonal to t1
              Tref(0, 1) = 0;
              //if t1 is parallel to the x1-axis t2 is set parallel to the x2-axis
              if (Tref(1, 0) == 0 && Tref(2, 0) == 0)
              {
                Tref(1, 1) = 1;
                Tref(2, 1) = 0;
              }

              //otherwise t2 is calculated from the scalar product with t1
              else
              {
                //setting t2(0)=0 and calculating other elements by setting scalar product t1 o t2 to zero
                double lin1norm = pow(pow(Tref(1, 0), 2) + pow(Tref(2, 0), 2),
                    0.5);
                Tref(1, 1) = -Tref(2, 0) / lin1norm;
                Tref(2, 1) = Tref(1, 0) / lin1norm;
              }

              //calculating t3 by crossproduct t1 x t2
              Tref(0, 2) = Tref(1, 0) * Tref(2, 1) - Tref(1, 1) * Tref(2, 0);
              Tref(1, 2) = Tref(2, 0) * Tref(0, 1) - Tref(2, 1) * Tref(0, 0);
              Tref(2, 2) = Tref(0, 0) * Tref(1, 1) - Tref(0, 1) * Tref(1, 0);

              /*the center triad in reference configuration is stored as a quaternion whose equivalent would be the rotation
               * from the identity matrix into the reference configuration*/
              LARGEROTATIONS::triadtoquaternion(Tref, Qconv_[0]);
              Qconv_[1] = Qconv_[0];
              Qnew_ = Qconv_;

              for (int i = 0; i < 4; i++)
              {
                rQnew_[i] = Qconv_[0];
                rQold_[i] = Qconv_[0];
              }

              //compute local rotational vectors phi according to Crisfield 1999,(4.6) in quaternion form
              LINALG::Matrix<4, 1> phi12;
              LARGEROTATIONS::quaternionproduct(Qnew_[1],LARGEROTATIONS::inversequaternion(Qnew_[0]), phi12);

              //according o Crisfield 1999, eq. (4.9), kappa equals the vector corresponding to phi12 divided by the element reference length
              LARGEROTATIONS::quaterniontoangle(phi12, kapparef_[0]);
              kapparef_[0].Scale(0.5 / jacobi_[0]);

            } //for(int numgp=0; numgp < gausspoints.nquad; numgp++)

            //Get the applied integrationpoints for exact integration of mass matrix

            DRT::UTILS::IntegrationPoints1D gausspointsmass(MyGaussRule(fnnode, lobattointegration));

            //vector whose numgp-th element is a 1xfnnode-matrix with all Lagrange polynomial basis functions evaluated at the numgp-th Gauss point
            std::vector<LINALG::Matrix<1, fnnode> > I(fnnode - 1);

            //vector whose numgp-th element is a 1xfnnode-matrix with the derivatives of all Lagrange polynomial basis functions evaluated at fnnode-1 Gauss points for elasticity
            std::vector<LINALG::Matrix<1, fnnode> > Iprime(fnnode - 1);

            //vector whose numgp-th element is a vector with fnnode elements, who represent the 3x3-matrix-shaped interpolation function \tilde{I}^fnnode at fnnode-1 Gauss points for elasticity according to according to (3.18), Jelenic 1999
            std::vector<std::vector<LINALG::Matrix<3, 3> > > Itilde(fnnode - 1);

            //vector whose numgp-th element is a vector with fnnode elements, who represent the 3x3-matrix-shaped interpolation function \tilde{I'}^fnnode at fnnode-1 Gauss points for elasticity according to according to (3.19), Jelenic 1999
            std::vector<std::vector<LINALG::Matrix<3, 3> > > Itildeprime(fnnode - 1);

            //vector with rotation matrices at fnnode-1 Gauss points for elasticity
            std::vector<LINALG::Matrix<3, 3> > Lambda(fnnode - 1);

            //vector whose numgp-th element is a 1xfnnode-matrix with all Lagrange polynomial basis functions evaluated at the fnnode Gauss points for mass matrix
            std::vector<LINALG::Matrix<1, fnnode> > Imass(fnnode);

            //vector whose numgp-th element is a vector with fnnode elements, who represent the 3x3-matrix-shaped interpolation function \tilde{I}^fnnode at the fnnode Gauss points for mass matrix according to according to (3.18), Jelenic 1999
            std::vector<std::vector<LINALG::Matrix<3, 3> > > Itildemass(fnnode);

            //exact Gauss integration requires as many integration points as nodes
            Qconvmass_.resize(fnnode);
            Qnewmass_.resize(fnnode);

            //evaluate basis functions
            evaluatebasisfunctionsandtriads<fnnode>(gausspoints, I, Iprime, Itilde, Itildeprime, Lambda, gausspointsmass, Imass, Itildemass);

            //Loop through all GPs for exact integration and compute initial jacobi determinant
            for (int numgp = 0; numgp < gausspointsmass.nquad; numgp++)
            {

              //Get position xi of GP
              const double xi = gausspointsmass.qxg[numgp][0];

              //Get derivatives of shapefunctions at GP
              DRT::UTILS::shape_function_1D_deriv1(shapefuncderiv, xi, distype);

              drdxi.Clear();
              //calculate dx/dxi and dz/dxi
              for (int node = 0; node < fnnode; node++)
                for (int dof = 0; dof < 3; dof++)
                  drdxi(dof) += shapefuncderiv(node) * xrefe[3 * node + dof];

              //Store Jacobi determinant with respect to reference configuration
              jacobimass_[numgp] = drdxi.Norm2();

              //in evaluatebasisfunctionsandtriads Qnewmass_ was initialized; copy that to Qconvmass_
              Qconvmass_[numgp].PutScalar(0);
              Qconvmass_[numgp] += Qnewmass_[numgp];

            } //for(int numgp=0; numgp < gausspointsmass.nquad; numgp++)

            //compute Jacobi determinant at gauss points for Lobatto quadrature (i.e. at nodes)
            for (int numgp = 0; numgp < fnnode; numgp++)
            {
              //Get position xi of nodes
              const double xi = -1.0 + 2 * numgp / (fnnode - 1);

              //Get derivatives of shapefunctions at GP
              DRT::UTILS::shape_function_1D_deriv1(shapefuncderiv, xi, distype);

              drdxi.Clear();
              //calculate dx/dxi and dz/dxi
              for (int node = 0; node < fnnode; node++)
                for (int dof = 0; dof < 3; dof++)
                  drdxi(dof) += shapefuncderiv(node) * xrefe[3 * node + dof];

              //Store Jacobi determinant for each node (Jacobi determinant refers by definition always to the reference configuration)
              jacobinode_[numgp] = drdxi.Norm2();
            } //for(int numgp=0; numgp< fnnode; numgp++)
          } //if(!isinit_)
          return;

        } //DRT::ELEMENTS::Beam3r::SetUpReferenceGeometry()

        //!\brief Set initial quaternions rQconv_ depending on nodal quaternions from already existing beam elements
        void SetInitialQuaternions(std::vector<LINALG::Matrix<4, 1> >& initquaternions);

        //! changes the actual linker length
        void SetReferenceLength(const double& scalefac);

        //! Set rotation to precalculated value and update additive rotation increment (this is something like a 'prescribed' rotation)
        void SetRotation(LINALG::Matrix<3,1>& theta);

      private:

        //! kind of integration to be performed
        enum IntegrationType
        {
          gaussexactintegration, gaussunderintegration, lobattointegration
        };
        //! variable saving whether element has already been initialized (then isinit_ == true)
        bool isinit_;

        //! reference coordinates of the interpolated nodes
        std::vector<double> xrefe_;
        std::vector<double> rotrefe_;

        //! vector storing the internal force vector
        Epetra_SerialDenseVector f_;
        //! \brief material forces at Gauss point
        LINALG::Matrix<3,1> Ngp_;
        //! storing epsilonn(0)
        double eps_;

        //------------- nodal triads and angle displacements for real nodes

        //! triads of 'real' nodes in quaternion form at the end of the preceeding time step
        std::vector<LINALG::Matrix<4, 1> > rQconv_;
        //! triads of 'real' nodes in quaternion form at the end of the preceeding iteration step
        std::vector<LINALG::Matrix<4, 1> > rQold_;
        //! triads of 'real' nodes in quaternion form during the current iteration step
        std::vector<LINALG::Matrix<4, 1> > rQnew_;
        //!rotation angle displacement at interpolated nodes; relates to the end of the preceeding time step
        std::vector<LINALG::Matrix<3, 1> > rdispthetaconv_;
        //!rotation angle displacement at interpolated nodes; relates to the end of the preceeding iteration step
        std::vector<LINALG::Matrix<3, 1> > rdispthetaold_;
        //!rotation angle displacement at interpolated nodes; relates to current configuration
        std::vector<LINALG::Matrix<3, 1> > rdispthetanew_;

        //------------- nodal triads and angle displacements for interpolated nodes
        //! nodal triads in quaternion form at the end of the preceeding time step
        std::vector<LINALG::Matrix<4, 1> > Qconv_;
        //! nodal triads in quaternion form at the end of the preceeding iteration step
        std::vector<LINALG::Matrix<4, 1> > Qold_;
        //! nodal triads in quaternion during the current iteration step
        std::vector<LINALG::Matrix<4, 1> > Qnew_;
        //! triads at Gauss points for exact integration in quaternion at the end of the preceeding time step (required for computation of angular velocity)
        std::vector<LINALG::Matrix<4, 1> > Qconvmass_;
        //! current triads at Gauss points for exact integration in quaternion (required for computation of angular velocity)
        std::vector<LINALG::Matrix<4, 1> > Qnewmass_;
        //! rotation increment set by SetRotation()
        LINALG::Matrix<4,1> Qrot_;
        //!rotation angle displacement at nodes; relates to the end of the preceeding time step
        std::vector<LINALG::Matrix<3, 1> > dispthetaconv_;
        //!rotation angle displacement at nodes; relates to the end of the preceeding iteration step
        std::vector<LINALG::Matrix<3, 1> > dispthetaold_;
        //!rotation angle displacement at nodes; relates to current configuration
        std::vector<LINALG::Matrix<3, 1> > dispthetanew_;

        //!initial curvature at Gauss points for elasticity (corresponding to \Lambda_0^t \Labmda'_0 in eq. (3.5), Crisfield 1999
        std::vector<LINALG::Matrix<3, 1> > kapparef_;

        //!node I for determination of reference triad, eq. (3.9), (3.10), Jelenic 1999
        int nodeI_;

        //!node J for determination of reference triad, eq. (3.9), (3.10), Jelenic 1999
        int nodeJ_;

        //! cross section area
        double crosssec_;

        //! cross section area including shear corretion factor
        double crosssecshear_;

        //! length scaling factor due to SetReferenceLength() (in order to be able to switch back and forth
        double lscalefac_;

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

        //! Position of the two interpolated nodes relativ to real nodes
        std::vector<double> mybindingposition_;

        //! @name methods for initilization of the element

        //!evaluates basis functions for a given type of Gauss rule  and triads at all Gauss points
        /*----------------------------------------------------------------------------------------------------------------------*
         |evaluates basis functions and their derivatives at Gauss points for elasticity and basis functions only at Gauss points|
         |for mass matrix; furthermore triads are evaluated at all Gauss points; note that derivatives of basis functions are    |
         |required for elasticity only, but not for mass (or damping) matrix                                          cyron 05/10|
         *----------------------------------------------------------------------------------------------------------------------*/
        template<int fnnode>
        void evaluatebasisfunctionsandtriads(const DRT::UTILS::IntegrationPoints1D&            gausspoints,
                                             std::vector<LINALG::Matrix<1, fnnode> >&          I,
                                             std::vector<LINALG::Matrix<1, fnnode> >&          Iprime,
                                             std::vector<std::vector<LINALG::Matrix<3, 3> > >& Itilde,
                                             std::vector<std::vector<LINALG::Matrix<3, 3> > >& Itildeprime,
                                             std::vector<LINALG::Matrix<3, 3> >&               Lambda,
                                             const DRT::UTILS::IntegrationPoints1D&            gausspointsmass,
                                             std::vector<LINALG::Matrix<1, fnnode> >&          Imass,
                                             std::vector<std::vector<LINALG::Matrix<3, 3> > >& Itildemass)
        {
          //declaration of variables

          //quaternion of relative rotation between node I and J according to (3.10), Jelenic 1999
          LINALG::Matrix<4, 1> QIJ;
          //angle of relative rotation between node I and J according to (3.10), Jelenic 1999
          LINALG::Matrix<3, 1> phiIJ;
          //quaternion of half relative rotation between node I and J according to (3.9), Jelenic 1999
          LINALG::Matrix<4, 1> QIJhalf;
          //quaternion of reference triad \Lambda_r according to (3.9), Jelenic 1999
          LINALG::Matrix<4, 1> Qr;
          //rotation quaternion between i-th nodal triadsand refenrece triad according to (3.8), Jelenic 1999
          LINALG::Matrix<4, 1> Qli;
          //rotation angles between nodal triads and refenrece triad according to (3.8), Jelenic 1999
          std::vector<LINALG::Matrix<3, 1> > Psili(fnnode);
          //interpolated local relative rotation \Psi^l at a certain Gauss point according to (3.11), Jelenic 1999
          LINALG::Matrix<3, 1> Psil;
          /*derivative of interpolated local relative rotation \Psi^l at a certain Gauss point according to
           *(3.11), Jelenic 1999, but not with respect to reference length parameter, but with respect to space of integration
           *parameter. Hence, the difference between the variable used in this code and the one used in the
           *equations of Jelenic 1999 is the Jacobi determinant*/
          LINALG::Matrix<3, 1> Psilprime;
          //rotation quaternion between triad at Gauss point and reference triad Qr
          LINALG::Matrix<4, 1> Ql;
          //rotation quaternion at Gauss point
          LINALG::Matrix<4, 1> Qgauss;
          //matrix of reference triad \Lambda_r according to (3.9), Jelenic 1999
          LINALG::Matrix<3, 3> Lambdar;

          //computation of element reference triad and rotations relative to this triad at nodes and Gauss points

          //compute reference rotation quaternion Qr corresponding to triad \Lambda_r according to eq. (3.10) and (3.9), Jelenic 1999
          LARGEROTATIONS::quaternionproduct(Qnew_[nodeJ_], LARGEROTATIONS::inversequaternion(Qnew_[nodeI_]), QIJ);
          LARGEROTATIONS::quaterniontoangle(QIJ, phiIJ);
          phiIJ.Scale(0.5);
          LARGEROTATIONS::angletoquaternion(phiIJ, QIJhalf);
          phiIJ.Scale(2.0);
          LARGEROTATIONS::quaternionproduct(QIJhalf, Qnew_[nodeI_], Qr);

          //compute triad Lambdar corresponding to reference rotation quaternion Qr
          LARGEROTATIONS::quaterniontotriad(Qr, Lambdar);

          //Compute relative rotations \Psi^l_i at all nodes i according to (3.8), Jelenic 1999
          for (int node = 0; node < fnnode; ++node)
          {
            LARGEROTATIONS::quaternionproduct(Qnew_[node], LARGEROTATIONS::inversequaternion(Qr), Qli);
            LARGEROTATIONS::quaterniontoangle(Qli, Psili[node]);
          }

          //evaluation of basis functions and triads at all Gauss points for elasticity
          for (int numgp = 0; numgp < gausspoints.nquad; numgp++)
          {

            //Get location and weight of GP in parameter space
            const double xi = gausspoints.qxg[numgp][0];

            //evaluate Lagrange polynomial basis functions and its derivatives at xi
            DRT::UTILS::shape_function_1D(I[numgp], xi, this->Shape());
            DRT::UTILS::shape_function_1D_deriv1(Iprime[numgp], xi,
                this->Shape());

            //compute local relative rotation \Psi^l and its derivative at current Gauss point according to (3.11), Jelenic 1999
            Psil.PutScalar(0);
            Psilprime.PutScalar(0);
            for (int node = 0; node < fnnode; ++node)
              for (int i = 0; i < 3; i++)
              {
                Psil(i) += I[numgp](node) * Psili[node](i);
                Psilprime(i) += Iprime[numgp](node) * Psili[node](i);
              }

            //compute relative rotation between triad at Gauss point and reference triad Qr
            LARGEROTATIONS::angletoquaternion(Psil, Ql);

            //compute rotation at Gauss point, i.e. the quaternion equivalent to \Lambda(s) in Crisfield 1999, eq. (4.7)
            LARGEROTATIONS::quaternionproduct(Ql, Qr, Qgauss);
            //compute rotation matrix at Gauss point, i.e. \Lambda(s) in Crisfield 1999, eq. (4.7)
            LARGEROTATIONS::quaterniontotriad(Qgauss, Lambda[numgp]);

            //compute at this Gauss point basis functions \tilde{I}^i and \tilde{I}^{i'} in (3.19), page 152, Jelenic 1999, for all nodes
            computeItilde<fnnode>(Psil, Itilde[numgp], phiIJ, Lambdar, Psili,I[numgp]);
            computeItildeprime<fnnode>(Psil, Psilprime, Itildeprime[numgp],phiIJ, Lambdar, Psili, I[numgp], Iprime[numgp]);

          }

          //evaluation of basis functions and triads at all Gauss points for mass matrix
          for (int numgp = 0; numgp < gausspointsmass.nquad; numgp++)
          {
            //Get location and weight of GP in parameter space
            const double xi = gausspointsmass.qxg[numgp][0];

            //evaluate Lagrange polynomial basis functions and its derivatives at xi
            DRT::UTILS::shape_function_1D(Imass[numgp], xi, this->Shape());

            //compute local relative rotation \Psi^l at current Gauss point according to (3.11), Jelenic 1999
            Psil.PutScalar(0);
            for (int node = 0; node < fnnode; ++node)
              for (int i = 0; i < 3; i++)
                Psil(i) += Imass[numgp](node) * Psili[node](i);

            //compute relative rotation between triad at Gauss point and reference triad Qr
            LARGEROTATIONS::angletoquaternion(Psil, Ql);
            Qnewmass_.resize(gausspointsmass.nquad);
            //compute rotation at Gauss point, i.e. the quaternion equivalent to \Lambda(s) in Crisfield 1999, eq. (4.7)
            LARGEROTATIONS::quaternionproduct(Ql, Qr, Qnewmass_[numgp]);

            //compute at this Gauss point basis functions \tilde{I}^i in (3.19), page 152, Jelenic 1999, for all nodes
            computeItilde<fnnode>(Psil, Itildemass[numgp], phiIJ, Lambdar, Psili, Imass[numgp]);

          }

          return;
        } // DRT::ELEMENTS::Beam3r::evaluatebasisfunctionsandtriads

        //Update of the nodal triads of the four real nodes
        template<int nnode>
        void UpdateRealNodalTriads(std::vector<double>& disp)
        {
          //rotational displacement at a certain node between this and last iteration step
          LINALG::Matrix<3, 1> deltatheta;
          //rotational displacement at a certain node between this and last iteration step in quaternion form
          LINALG::Matrix<4, 1> deltaQ;
          for (int node = 0; node < nnode; node++)
          {
            /*rotation increment relative to configuration in last iteration step is difference between current rotation
             *entry in displacement vector minus rotation entry in displacement vector in last iteration step*/
            for (int i = 0; i < 3; i++)
              rdispthetanew_[node](i) = disp[6 * node + 3 + i];

            deltatheta = rdispthetanew_[node];
            deltatheta -= rdispthetaold_[node];

            //compute quaternion from rotation angle relative to last configuration
            LARGEROTATIONS::angletoquaternion(deltatheta, deltaQ);

            //multiply relative rotation with rotation in last configuration to get rotation in new configuration
            LARGEROTATIONS::quaternionproduct(rQold_[node], deltaQ, rQnew_[node]);

            //renormalize quaternion to keep its absolute value one even in case of long simulations and intricate calculations
            rQnew_[node].Scale(1 / rQnew_[node].Norm2());

          }
        } //DRT::ELEMENTS::BeamCL::UpdateRealNodalTriads

        // Interpolate Triads of interpolated Nodes
        void InterpolateFictiveNodalTriads(std::vector<LINALG::Matrix<1, 2> >&               Ibp,
                                           std::vector<LINALG::Matrix<3, 3> >&               rLambdar,
                                           std::vector<LINALG::Matrix<3, 1> >&               rphiIJ,
                                           std::vector<std::vector<LINALG::Matrix<3, 1> > >& rPsili,
                                           std::vector<LINALG::Matrix<3, 1> >&               rPsil)
        {
          const int fnnode = 2;
          std::vector<LINALG::Matrix<4, 1> > rQr(2);
          std::vector<LINALG::Matrix<3, 3> > rLambdaXi(2);
          //quaternion of relative rotation between node I and J
          LINALG::Matrix<4, 1> rQIJ;
          //quaternion of half relative rotation between node I and J
          LINALG::Matrix<4, 1> rQIJhalf;
          //rotation quaternion between triad at Gauss point and reference triad Qr
          LINALG::Matrix<4, 1> rQl;
          LINALG::Matrix<4, 1> rQli;
          //Compute Reference Triad (Lambda_r) for Filament 1 and 2
          for (int filament = 0; filament < 2; filament++)
          {
            LARGEROTATIONS::quaternionproduct(rQnew_[2 * filament + 1], LARGEROTATIONS::inversequaternion(rQnew_[2 * filament]), rQIJ);
            LARGEROTATIONS::quaterniontoangle(rQIJ, rphiIJ[filament]);
            rphiIJ[filament].Scale(0.5);
            LARGEROTATIONS::angletoquaternion(rphiIJ[filament], rQIJhalf);
            rphiIJ[filament].Scale(2.0);
            LARGEROTATIONS::quaternionproduct(rQIJhalf, rQnew_[2 * filament], rQr[filament]);
            LARGEROTATIONS::quaterniontotriad(rQr[filament], rLambdar[filament]);
          }
          //Compute relative rotations \Psi^l_i at all nodes i according to (3.8), Jelenic 1999
          for (int filament = 0; filament < 2; filament++)
          {
            for (int node = 0; node < fnnode; ++node)
            {
              LARGEROTATIONS::quaternionproduct(rQnew_[2 * filament + node], LARGEROTATIONS::inversequaternion(rQr[filament]), rQli);
              LARGEROTATIONS::quaterniontoangle(rQli, rPsili[filament][node]);
            }
          }
          // Compute Triad Lambda at interpolated nodes (f1 and f2)
          rPsil[1].PutScalar(0);
          rPsil[0].PutScalar(0);
          for (int filament = 0; filament < 2; filament++)
            for (int fnode = 0; fnode < 2; fnode++)
              for (int i = 0; i < 3; i++)
                rPsil[filament](i) += Ibp[filament](fnode) * rPsili[filament][fnode](i);

          for (int j = 0; j < 2; j++)
          {
            LARGEROTATIONS::angletoquaternion(rPsil[j], rQl);
            LARGEROTATIONS::quaternionproduct(rQl, rQr[j], Qnew_[j]); // Qxi is Quaternion at f1 and f2
            LARGEROTATIONS::quaterniontotriad(Qnew_[j], rLambdaXi[j]); // Lambda is Triad at f1 and f2
          }
        }

        //@}

        //! @name Internal calculation methods

        //!compute convected strain at a certain Gauss point with triad rotmat according to Crisfield 1999, eq. (3.4) and eq. (4.9)
        void computestrain(const LINALG::Matrix<3, 1>& rprime,
                           const LINALG::Matrix<3, 3>& Lambda,
                           LINALG::Matrix<3, 1>&       gamma,
                           LINALG::Matrix<3, 1>&       kappa);

        //!computing convected stresses from convected strains vector and returing also convected constitutive matrix between both
        void strainstress(const LINALG::Matrix<3, 1>& gamma,
                          const LINALG::Matrix<3, 1>& kappa,
                          LINALG::Matrix<3, 1>& stressN,
                          LINALG::Matrix<3, 3>& CN,
                          LINALG::Matrix<3, 1>& stressM,
                          LINALG::Matrix<3, 3>& CM);

        //!push forward strains, stresses and constitutive matrices to their spatial counterparts by rotation matrix Lambda according to Romero 2004, eq. (3.10)
        void pushforward(const LINALG::Matrix<3, 3>& Lambda,
                         const LINALG::Matrix<3, 1>& stressN,
                         const LINALG::Matrix<3, 3>& CN,
                         const LINALG::Matrix<3, 1>& stressM,
                         const LINALG::Matrix<3, 3>& CM,
                         LINALG::Matrix<3, 1>& stressn,
                         LINALG::Matrix<3, 3>& cn,
                         LINALG::Matrix<3, 1>& stressm,
                         LINALG::Matrix<3, 3>& cm);

        //!compute convected strain at certain Gauss point with triad rotmat according to Crisfield 1999, eq. (3.4) and eq. (4.9)
        void computestrain(const LINALG::Matrix<3, 1> Psil,
                           const LINALG::Matrix<3, 1> Psilprime,
                           const double& xi,
                           const LINALG::Matrix<3, 1>& rprime,
                           const LINALG::Matrix<3, 3>& rotmat,
                           LINALG::Matrix<3, 1>& gamma,
                           LINALG::Matrix<3, 1>& kappa);

        //!compute \tilde{I}^i in (3.18), page 152, Jelenic 1999, for all nodes i at a certaom Gauss point
        /*----------------------------------------------------------------------------------------------------------------------*
         |compute \tilde{I}^i in (3.18), page 152, Jelenic 1999, for all nodes i at a certaom Gauss point            cyron 04/10|
         *----------------------------------------------------------------------------------------------------------------------*/
        template<int fnnode>
        void computeItilde(const LINALG::Matrix<3, 1>& Psil,
                           std::vector<LINALG::Matrix<3, 3> >& Itilde,
                           const LINALG::Matrix<3, 1>& phiIJ,
                           const LINALG::Matrix<3, 3>& Lambdar,
                           const std::vector<LINALG::Matrix<3, 1> >& Psili,
                           const LINALG::Matrix<1, fnnode>& funct)
        {
          //auxiliary matrices for storing intermediate results
          LINALG::Matrix<3, 3> auxmatrix;
          auxmatrix.PutScalar(0);

          //make sure that Itildeprime has proper dimensions
          Itilde.resize(fnnode);

          //compute squared brackets term in (3.18), Jelenic 1999v
          LINALG::Matrix<3, 3> squaredbrackets;
          squaredbrackets.PutScalar(0);
          for (int node = 0; node < fnnode; ++node)
          {
            auxmatrix = LARGEROTATIONS::Tmatrix(Psili[node]);
            auxmatrix.Scale(funct(node));
            squaredbrackets -= auxmatrix;
          }
          for (int i = 0; i < 3; i++)
            squaredbrackets(i, i) += 1;

          //loop through all nodes i
          for (int node = 0; node < fnnode; ++node)
          {

            //compute rightmost term in curley brackets in (3.18), Jelenic 1999
            Itilde[node].PutScalar(0);
            Itilde[node].Multiply(LARGEROTATIONS::Tinvmatrix(Psil), LARGEROTATIONS::Tmatrix(Psili[node]));
            Itilde[node].Scale(funct(node));

            //if node i is node I then add squared bracktets term times v_I
            if (node == nodeI_)
            {
              auxmatrix.Multiply(squaredbrackets, vI(phiIJ));
              Itilde[node] += auxmatrix;
            }

            //if node i is node J then add squared bracktets term times v_J
            if (node == nodeJ_)
            {
              auxmatrix.Multiply(squaredbrackets, vJ(phiIJ));
              Itilde[node] += auxmatrix;
            }

            //now the term in the curley bracktets has been computed and has to be rotated by \Lambda_r and \Lambda_r^t
            auxmatrix.MultiplyNT(Itilde[node], Lambdar);
            Itilde[node].MultiplyNN(Lambdar, auxmatrix);

          }

          return;
        } // DRT::ELEMENTS::BeamCL::computeItilde

        //!compute \tilde{I}^{i'} in (3.19), page 152, Jelenic 1999, for all nodes i at a certaom Gauss point
        /*----------------------------------------------------------------------------------------------------------------------*
         |compute \tilde{I}^{i'} in (3.19), page 152, Jelenic 1999 for all nodes i at a certaom Gauss point          cyron 04/10|
         *----------------------------------------------------------------------------------------------------------------------*/
        template<int fnnode>
        void computeItildeprime(const LINALG::Matrix<3, 1>& Psil,
                                const LINALG::Matrix<3, 1>& Psilprime,
                                std::vector<LINALG::Matrix<3, 3> >& Itildeprime,
                                const LINALG::Matrix<3, 1>& phiIJ,
                                const LINALG::Matrix<3, 3>& Lambdar,
                                const std::vector<LINALG::Matrix<3, 1> >& Psili,
                                const LINALG::Matrix<1, fnnode>& funct,
                                const LINALG::Matrix<1, fnnode>& deriv)
        {
          //auxiliary matrices for storing intermediate results
          LINALG::Matrix<3, 3> auxmatrix;

          //make sure that Itildeprime has proper dimensions
          Itildeprime.resize(fnnode);

          //matrix d(T^{-1})/dx
          LINALG::Matrix<3, 3> dTinvdx;
          LARGEROTATIONS::computedTinvdx(Psil, Psilprime, dTinvdx);

          //compute T^{~} according to remark subsequent to (3.19), Jelenic 1999
          LINALG::Matrix<3, 3> Ttilde;
          Ttilde.PutScalar(0);
          for (int node = 0; node < fnnode; ++node)
          {
            auxmatrix = LARGEROTATIONS::Tmatrix(Psili[node]);
            auxmatrix.Scale(funct(node));
            Ttilde += auxmatrix;
          }

          //compute T^{~'} according to remark subsequent to (3.19), Jelenic 1999
          LINALG::Matrix<3, 3> Ttildeprime;
          Ttildeprime.PutScalar(0);
          for (int node = 0; node < fnnode; ++node)
          {
            auxmatrix = LARGEROTATIONS::Tmatrix(Psili[node]);
            auxmatrix.Scale(deriv(node));
            Ttildeprime += auxmatrix;
          }

          //compute first squared brackets term in (3.18), Jelenic 1999v
          LINALG::Matrix<3, 3> squaredbrackets;
          squaredbrackets.PutScalar(0);
          squaredbrackets.Multiply(dTinvdx, Ttilde);
          auxmatrix.Multiply(LARGEROTATIONS::Tinvmatrix(Psil), Ttildeprime);
          squaredbrackets += auxmatrix;

          //loop through all nodes i
          for (int node = 0; node < fnnode; ++node)
          {

            //compute first term in second squared brackets
            Itildeprime[node] = dTinvdx;
            Itildeprime[node].Scale(funct(node));

            //compute second term in second squared brackets
            auxmatrix.PutScalar(0);
            auxmatrix += LARGEROTATIONS::Tinvmatrix(Psil);
            auxmatrix.Scale(deriv(node));

            //compute second squared brackets
            auxmatrix += Itildeprime[node];

            //compute second squared bracketds time T(\Psi^l_j)
            Itildeprime[node].Multiply(auxmatrix, LARGEROTATIONS::Tmatrix(Psili[node]));

            //if node i is node I then add first squared bracktets term times v_I
            if (node == nodeI_)
            {
              auxmatrix.Multiply(squaredbrackets, vI(phiIJ));
              Itildeprime[node] -= auxmatrix;
            }

            //if node i is node J then add first squared bracktets term times v_J
            if (node == nodeJ_)
            {
              auxmatrix.Multiply(squaredbrackets, vJ(phiIJ));
              Itildeprime[node] -= auxmatrix;
            }

            //now the term in the curley bracktets has been computed and has to be rotated by \Lambda_r and \Lambda_r^t
            auxmatrix.MultiplyNT(Itildeprime[node], Lambdar);
            Itildeprime[node].MultiplyNN(Lambdar, auxmatrix);
          }

          return;
        } // DRT::ELEMENTS::Beam3r::computeItildeprime

        //!compute matrix v_I as outlined in the equations above (3.15) on page 152 of Jelenic 1999
        /*----------------------------------------------------------------------------------------------------------------------*
         |compute matrix v_I as outlined in the equations above (3.15) on page 152 of Jelenic 1999                   cyron 04/10|
         *----------------------------------------------------------------------------------------------------------------------*/
        LINALG::Matrix<3, 3> vI(const LINALG::Matrix<3, 1>& phiIJ)
        {
          //matrix v_I
          LINALG::Matrix<3, 3> result;

          LARGEROTATIONS::computespin(result, phiIJ);
          if (phiIJ.Norm2() == 0)
            result.Scale(0.25);
          else
            result.Scale(tan(phiIJ.Norm2() / 4.0) / phiIJ.Norm2());

          for (int i = 0; i < 3; i++)
            result(i, i) += 1;

          result.Scale(0.5);

          return result;
        } // DRT::ELEMENTS::Beam3r::vI

        //!compute matrix v_J as outlined in the equations above (3.15) on page 152 of Jelenic 1999
        /*----------------------------------------------------------------------------------------------------------------------*
         |compute matrix v_J as outlined in the equations above (3.15) on page 152 of Jelenic 1999                   cyron 04/10|
         *----------------------------------------------------------------------------------------------------------------------*/
        LINALG::Matrix<3, 3> vJ(const LINALG::Matrix<3, 1>& phiIJ)
        {
          //matrix v_J
          LINALG::Matrix<3, 3> result;

          LARGEROTATIONS::computespin(result, phiIJ);
          if (phiIJ.Norm2() == 0)
            result.Scale(-0.25);
          else
            result.Scale(-tan(phiIJ.Norm2() / 4.0) / phiIJ.Norm2());

          for (int i = 0; i < 3; i++)
            result(i, i) += 1;

          result.Scale(0.5);

          return result;
        } // DRT::ELEMENTS::Beam3r::vJ

        //!compute derivative r' (rprime) of line of centroids at Gauss point x with respect to Jelenic 1999 , eq. (2.12)
        /*----------------------------------------------------------------------------------------------------------------------*
         |compute derivative r' (rprime) of line of centroids at Gauss point x with respect to Jelenic 1999, eq. (2.12)        |
         |                                                                                                           cyron 04/10|
         *----------------------------------------------------------------------------------------------------------------------*/
        template<int fnnode, int dof>
        void curvederivative(const std::vector<double>& disp,
                             const LINALG::Matrix<1, fnnode> deriv,
                             LINALG::Matrix<3, 1>& rprime,
                             const double& jacobi)
        {
          //initialize rprime
          rprime.PutScalar(0);

          std::vector<double> rxrefe;
          std::vector<double> xrefe;
          rxrefe.resize(12);
          xrefe.resize(6);
          for (int i = 0; i < 6; i++)
            xrefe[i] = 0;

          for (int node = 0; node < 4; node++)
            for (int d = 0; d < 3; d++)
              rxrefe[node * 3 + d] = Nodes()[node]->X()[d];

          std::vector<LINALG::Matrix<1, 2> > Ibp(2);
          const DiscretizationType dtype = this->Shape();
          for (int filament = 0; filament < 2; filament++)
            DRT::UTILS::shape_function_1D(Ibp[filament],
                mybindingposition_[filament], dtype);

          for (int filament = 0; filament < 2; filament++)
            for (int k = 0; k < 2; k++)
              for (int i = 0; i < 3; i++)
                xrefe[i + 3 * filament] += Ibp[filament](k)
                    * rxrefe[i + 3 * k + 6 * filament];

          for (int i = 0; i < dof; ++i)
            for (int node = 0; node < fnnode; ++node)
              rprime(i) += (xrefe[3 * node + i] + disp[6 * node + i])
                  * deriv(node);

          /*so far we have computed the derivative of the curve with respect to the element parameter \xi \in [-1;1];
           *as r' in (2.12) is the derivative with respect to the reference length, we have to divided it by the Jacobi
           *determinant at the respective point*/
          rprime.Scale(1.0 / jacobi);

          return;
        } // DRT::ELEMENTS::Beam3r::curvederivative

        //!calculation of nonlinear stiffness and mass matrix
        template<int nnode>
        void b3_nlnstiffmass(Teuchos::ParameterList& params,
                             std::vector<double>& vel,
                             std::vector<double>& disp,
                             Epetra_SerialDenseMatrix* stiffmatrix,
                             Epetra_SerialDenseMatrix* massmatrix,
                             Epetra_SerialDenseVector* force);

        //!calculation of elastic energy
        template<int fnnode>
        void b3_energy(Teuchos::ParameterList& params,
                       std::vector<double>& disp,
                       Epetra_SerialDenseVector* intenergy);

        //! lump mass matrix
        template<int fnnode>
        void lumpmass(Epetra_SerialDenseMatrix* emass);

        //!determine Gauss rule from required type of integration
        DRT::UTILS::GaussRule1D MyGaussRule(int fnnode,
                                            IntegrationType integrationtype);

        //@}

        //! @name Methods for Brownian dynamics simulations

        //!compute damping matrix and statistical forces according to fluctuation dissipation theorem
        template<int fnnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
        void CalcBrownian(Teuchos::ParameterList& params, //!<parameter list
                          const std::vector<double>& vel, //!< element velocity vector
                          const std::vector<double>& disp, //!<element disp vector
                          Epetra_SerialDenseMatrix stiffmatrix, //!< element stiffness matrix
                          Epetra_SerialDenseVector force,
                          std::vector<LINALG::Matrix<1, fnnode> >& Imass,
                          std::vector<std::vector<LINALG::Matrix<3, 3> > >& Itildemass); //!< element internal force vector

        //!computes damping coefficients for employed friction model
        void MyDampingConstants(Teuchos::ParameterList& params, LINALG::Matrix<3, 1>& gamma,
                                const INPAR::STATMECH::FrictionModel& frictionmodel);

        //!computes rotational damping forces and stiffness
        template<int fnnode> //number of nodes
        void MyRotationalDamping(Teuchos::ParameterList& params, //!<parameter list
                                 const std::vector<double>& vel, //!< element velocity vector
                                 const std::vector<double>& disp, //!<element disp vector
                                 Epetra_SerialDenseMatrix stiffmatrix, //!< element stiffness matrix
                                 Epetra_SerialDenseVector force, //!< element internal force vector
                                 const DRT::UTILS::IntegrationPoints1D& gausspointsdamping,
                                 const std::vector<LINALG::Matrix<1, fnnode> >& Idamping,
                                 const std::vector<std::vector<LINALG::Matrix<3, 3> > >& Itildedamping,
                                 const std::vector<LINALG::Matrix<4, 1> >& Qconvdamping,
                                 const std::vector<LINALG::Matrix<4, 1> >& Qnewdamping);

        //!computes translational damping forces and stiffness
        template<int fnnode, int ndim, int dof> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node
        void MyTranslationalDamping(Teuchos::ParameterList& params, //!<parameter list
                                    const std::vector<double>& vel, //!< element velocity vector
                                    const std::vector<double>& disp, //!<element disp vector
                                    Epetra_SerialDenseMatrix stiffmatrix, //!< element stiffness matrix
                                    Epetra_SerialDenseVector force); //!< element internal force vector

        //!computes stochastic translational forces and resulting stiffness
        template<int fnnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
        void MyStochasticForces(Teuchos::ParameterList& params, //!<parameter list
            const std::vector<double>& vel, //!< element velocity vector
            const std::vector<double>& disp, //!<element disp vector
            Epetra_SerialDenseMatrix stiffmatrix, //!< element stiffness matrix
            Epetra_SerialDenseVector force); //!< element internal force vector

        //!computes stochastic rotational moments and resulting stiffness
        template<int fnnode, int randompergauss> //number of nodes, number of random numbers required per Gauss point
        void MyStochasticMoments(Teuchos::ParameterList& params, //!<parameter list
                                 const std::vector<double>& vel, //!< element velocity vector
                                 const std::vector<double>& disp, //!<element disp vector
                                 Epetra_SerialDenseMatrix stiffmatrix, //!< element stiffness matrix
                                 Epetra_SerialDenseVector force, //!< element internal force vector
                                 const DRT::UTILS::IntegrationPoints1D& gausspointsdamping,
                                 const std::vector<LINALG::Matrix<1, fnnode> >& Idamping,
                                 const std::vector<std::vector<LINALG::Matrix<3, 3> > >& Itildedamping,
                                 const std::vector<LINALG::Matrix<4, 1> >& Qconvdamping,
                                 const std::vector<LINALG::Matrix<4, 1> >& Qnewdamping);

        //!computes velocity of background fluid and gradient of that velocity at a certain evaluation point in the physical space and adds respective terms to internal forces and damping matrix
        template<int ndim> //number of dimensions of embedding space
        void MyBackgroundVelocity(Teuchos::ParameterList& params, //!<parameter list
                                  const LINALG::Matrix<ndim, 1>& evaluationpoint, //!<point at which background velocity and its gradient has to be computed
                                  LINALG::Matrix<ndim, 1>& velbackground, //!< velocity of background fluid
                                  LINALG::Matrix<ndim, ndim>& velbackgroundgrad); //!<gradient of velocity of background fluid

        //!shifts nodes so that proper evaluation is possible even in case of periodic boundary conditions
        template<int fnnode, int ndim> //number of nodes, number of dimensions
        void NodeShift(Teuchos::ParameterList& params, //!<parameter list
                       std::vector<double>& disp); //!<element disp vector

        //@}

        // don't want = operator
        BeamCL& operator =(const BeamCL& old);

    };
    // class Beam3r

// << operator
    std::ostream& operator <<(std::ostream& os, const DRT::Element& ele);

  } // namespace Elements
} // namespace DRT

#endif  // #ifndef BEAM_H
