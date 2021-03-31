/*----------------------------------------------------------------------------*/
/*! \file
\brief three dimensional total Lagrange truss element

\level 3


*/
/*---------------------------------------------------------------------------*/

// header file only included if not yet included!
#ifndef TRUSS3_H
#define TRUSS3_H

#include "../drt_lib/drt_elementtype.H"
#include "../drt_fem_general/drt_utils_integration.H"

#include "Sacado.hpp"
using FAD = Sacado::Fad::DFad<double>;

namespace STR
{
  namespace ELEMENTS
  {
    class ParamsInterface;
  }
}  // namespace STR

namespace DRT
{
  namespace ELEMENTS
  {
    class Truss3Type : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "Truss3Type"; }

      static Truss3Type& Instance();

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns, const double* x0,
          int numdf, int dimns) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static Truss3Type instance_;
    };

    /*!
     \brief three dimensional total Lagrange truss element

     */
    class Truss3 : public DRT::Element
    {
     public:
      //! @name Friends
      friend class Truss3Type;

      //! @name Constructors and destructors and related methods

      /*!
       \brief Standard Constructor

       \param id    (in): A globally unique element id
       \param etype (in): Type of element
       \param owner (in): owner processor of the element
       */
      Truss3(int id, int owner);

      /*!
       \brief Copy Constructor

       Makes a deep copy of a Element
       */
      Truss3(const Truss3& old);

      ~Truss3() override = default;


      DRT::Element* Clone() const override;

      // TODO: remove once truss3 element is fixed and no longer expects more dofs (6) than it can
      // inherently handle (3)...
      void LocationVector(
          const Discretization& dis, LocationArray& la, bool doDirichlet) const override;


      DiscretizationType Shape() const override;

      int UniqueParObjectId() const override { return Truss3Type::Instance().UniqueParObjectId(); }

      void Pack(DRT::PackBuffer& data) const override;


      void Unpack(const std::vector<char>& data) override;

      DRT::ElementType& ElementType() const override { return Truss3Type::Instance(); }

      //! \brief sets reference coordinates X_ and refenrence length lrefe_ for elements added to
      //! the discretization
      //!
      //! \param xrefe     nodal coordinates in reference frame
      void SetUpReferenceGeometry(const std::vector<double>& xrefe);

      //! scale truss reference length
      void ScaleReferenceLength(double scalefac);

      int NumLine() const override { return 1; }

      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

      int NumDofPerNode(const DRT::Node& node) const override
      {
        /*note: this is not necessarily the number of DOF assigned to this node by the
         *discretization finally, but only the number of DOF requested for this node by this
         *element; the discretization will finally assign the maximal number of DOF to this node
         *requested by any element connected to this node*/
        return 3;
      }

      int NumDofPerElement() const override { return 0; }

      /*!
       \brief Set cross section area of this element
       */
      void SetCrossSec(const double& crosssec);

      /*!
       \brief computes the number of different random numbers required in each time step for
       generation of stochastic forces
       */
      int HowManyRandomNumbersINeed();

      void Print(std::ostream& os) const override;

      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef) override;

      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3) override;

      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseMatrix* elemat1 = nullptr) override;

      /*!
       \brief Evaluate PTC addition to stiffness for free Brownian motion

       An element derived from this class uses the Evaluate method to receive commands
       and parameters from some control routine in params and evaluates a  statistical Neumann
       boundary condition as used in the problem type STATISTICAL MECHANICS

       \param params (in/out)       : ParameterList for communication between control routine and
       elements \param vector<double> mydisp : current nodal displacement \param elemat1 (out) :
       artificial damping matrix to be filled by element

       \return 0 if successful, negative otherwise
       */
      template <int nnode, int ndim, int dof>  // number of nodes, number of dimensions of embedding
                                               // space, number of degrees of freedom per node
      void EvaluatePTC(Teuchos::ParameterList& params, Epetra_SerialDenseMatrix& elemat1);

      //! \brief return reference length
      double Lcurr() const { return lcurr_; }

      //! get internal (elastic) energy of element
      double GetInternalEnergy() const { return eint_; };

      void SetParamsInterfacePtr(const Teuchos::ParameterList& p) override;

      inline bool IsParamsInterface() const override { return (not interface_ptr_.is_null()); }

      Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> ParamsInterfacePtr() override;

      /*! \brief calculate internal force vector and stiffness matrix based
       *         on absolute nodal positions (using kinematic type tr3_totlag)
       *
       *  \author grill
       *  \date 08/17 */
      void CalcInternalForceStiffTotLag(const LINALG::Matrix<6, 1>& nodal_positions_totlag,
          Epetra_SerialDenseVector& forcevec, Epetra_SerialDenseMatrix& stiffmat);

     protected:
      /** \brief get access to the interface
       *
       *  \author hiermeier
       *  \date 04/16 */
      inline STR::ELEMENTS::ParamsInterface& ParamsInterface()
      {
        if (not IsParamsInterface()) dserror("The interface ptr is not set!");
        return *interface_ptr_;
      }

     private:
      //! kind of integration to be performed
      enum IntegrationType
      {
        gaussexactintegration,
        gaussunderintegration,
        lobattointegration
      };

      //! possible kinematic types
      enum KinematicType
      {
        tr3_totlag,
        tr3_engstrain
      };

      /*! \brief interface ptr
       *
       *  data exchange between the element and the time integrator. */
      Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_;

      //! container containing all kind of Truss3 thisandthat
      DRT::Container data_;

      //! variable saving whether element has already been initialized (then isinit_ == true)
      bool isinit_;

      //! reference node position
      LINALG::Matrix<6, 1> X_;

      //! reference tangent position
      std::vector<LINALG::Matrix<3, 1>> trefNode_;

      //! reference tangent position
      LINALG::Matrix<1, 3> diff_disp_ref_;

      //! material type
      int material_;

      //! length in reference configuration
      double lrefe_;

      //! length in reference configuration
      double lcurr_;

      //! Vector holding value of Jacobi determinant jacobi for complete integration of massmatrix
      std::vector<double> jacobimass_;

      //! vector holding value of Jacobi determinant jacobi at nodes
      std::vector<double> jacobinode_;

      //! cross section area
      double crosssec_;

      //! internal energy of element
      double eint_;

      //! Kinematic type
      KinematicType kintype_;

      //! gaussrule_ will be initialized automatically to a 2 point integration rule
      DRT::UTILS::GaussRule1D gaussrule_;

      //! @name Internal calculation methods

      //! determine Gauss rule from required type of integration
      DRT::UTILS::GaussRule1D MyGaussRule(int nnode, IntegrationType integrationtype);

      //! calculation of nonlinear stiffness and mass matrix switching between total lagrange and
      //! enginerring strains
      void t3_nlnstiffmass(Teuchos::ParameterList& params, std::vector<double>& vel,
          std::vector<double>& disp, Epetra_SerialDenseMatrix* stiffmatrix,
          Epetra_SerialDenseMatrix* massmatrix, Epetra_SerialDenseVector* force);

      //! calculation of elastic energy
      void t3_energy(Teuchos::ParameterList& params, std::vector<double>& disp,
          Epetra_SerialDenseVector* intenergy);

      //! calculation of nonlinear stiffness and mass matrix
      void t3_nlnstiffmass_totlag(LINALG::Matrix<1, 6>& DummyDisp,
          Epetra_SerialDenseMatrix& DummyStiffMatrix, Epetra_SerialDenseMatrix* massmatrix,
          Epetra_SerialDenseVector& DummyForce);

      /*!
       * \brief calculate force, nonlinear stiffness and mass matrix using a engineering strain
       * measure.
       *
       * Unlike the fully nonlinear implementation of #t3_nlnstiffmass_totlag, this function uses
       * \f$\varepsilon=\Delta d / d\f$ as strain measure.
       *
       */
      void t3_nlnstiffmass_engstr(const LINALG::Matrix<1, 6>& DummyDisp,
          Epetra_SerialDenseMatrix& DummyStiffMatrix, Epetra_SerialDenseMatrix* massmatrix,
          Epetra_SerialDenseVector& DummyForce);

      //! lump mass matrix
      void t3_lumpmass(Epetra_SerialDenseMatrix* emass);

      //@}

      //! @name Methods for Brownian dynamics simulations

      //! compute damping matrix and statistical forces according to fluctuation dissipation theorem
      template <int nnode, int ndim, int dof,
          int randompergauss>  // number of nodes, number of dimensions of embedding space, number
                               // of degrees of freedom per node, number of random numbers required
                               // per Gauss point
      void
      CalcBrownian(Teuchos::ParameterList& params,     //!< parameter list
          const LINALG::Matrix<1, 6>& DummyVel,        //!< element velocity vector
          const LINALG::Matrix<1, 6>& DummyDisp,       //!< element disp vector
          Epetra_SerialDenseMatrix& DummyStiffMatrix,  //!< element stiffness matrix
          Epetra_SerialDenseVector& DummyForce);       //!< element internal force vector

      //! computes damping coefficients for employed friction model
      void MyDampingConstants(Teuchos::ParameterList& params, LINALG::Matrix<3, 1>& gamma);

      //! computes translational damping forces and stiffness
      template <int nnode, int ndim,
          int dof>  // number of nodes, number of dimensions of embedding space, number of degrees
                    // of freedom per node
      void
      MyTranslationalDamping(Teuchos::ParameterList& params,  //!< parameter list
          const LINALG::Matrix<1, 6>& DummyVel,               //!< element velocity vector
          const LINALG::Matrix<1, 6>& DummyDisp,              //!< element disp vector
          Epetra_SerialDenseMatrix& DummyStiffMatrix,         //!< element stiffness matrix
          Epetra_SerialDenseVector& DummyForce);              //!< element internal force vector

      //! computes stochastic translational forces and resulting stiffness
      template <int nnode, int ndim, int dof,
          int randompergauss>  // number of nodes, number of dimensions of embedding space, number
                               // of degrees of freedom per node, number of random numbers required
                               // per Gauss point
      void
      MyStochasticForces(Teuchos::ParameterList& params,  //!< parameter list
          const LINALG::Matrix<1, 6>& DummyVel,           //!< element velocity vector
          const LINALG::Matrix<1, 6>& DummyDisp,          //!< element disp vector
          Epetra_SerialDenseMatrix& DummyStiffMatrix,     //!< element stiffness matrix
          Epetra_SerialDenseVector& DummyForce);          //!< element internal force vector

      //! computes velocity of background fluid and gradient of that velocity at a certain
      //! evaluation point in the physical space and adds respective terms to internal forces and
      //! damping matrix
      template <int ndim>  // number of dimensions of embedding space
      void MyBackgroundVelocity(Teuchos::ParameterList& params,  //!< parameter list
          const LINALG::Matrix<ndim, 1>& evaluationpoint,  //!< point at which background velocity
                                                           //!< and its gradient has to be computed
          LINALG::Matrix<ndim, 1>& velbackground,          //!< velocity of background fluid
          LINALG::Matrix<ndim, ndim>&
              velbackgroundgrad);  //!< gradient of velocity of background fluid

      //! shifts nodes so that proper evaluation is possible even in case of periodic boundary
      //! conditions
      template <int nnode, int ndim>                  // number of nodes, number of dimensions
      void NodeShift(Teuchos::ParameterList& params,  //!< parameter list
          std::vector<double>& disp);                 //!< element disp vector


      // don't want = operator
      Truss3& operator=(const Truss3& old);
    };

    // << operator
    std::ostream& operator<<(std::ostream& os, const DRT::Element& ele);

  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // #ifndef TRUSS3_H
