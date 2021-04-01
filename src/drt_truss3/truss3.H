/*----------------------------------------------------------------------------*/
/*! \file
\brief three dimensional total Lagrange truss element

\level 3


*/
/*---------------------------------------------------------------------------*/
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
      void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns, const double* x0,
          int numdf, int dimns) override;

      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      static Truss3Type& Instance();

      std::string Name() const override { return "Truss3Type"; }

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

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

      /*!
       \brief Standard Constructor

       \param id    (in): A globally unique element id
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

      //! \brief calculate internal force vector and stiffness matrix based on absolute nodal
      //! positions (using kinematic type tr3_totlag)
      void CalcInternalForceStiffTotLag(const LINALG::Matrix<6, 1>& nodal_positions_totlag,
          Epetra_SerialDenseVector& forcevec, Epetra_SerialDenseMatrix& stiffmat);

      DRT::ElementType& ElementType() const override { return Truss3Type::Instance(); }

      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3) override;

      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseMatrix* elemat1 = nullptr) override;

      //! \tparam nnode    number of nodes
      //! \tparam ndim     number of dimensions of embedding space
      //! \tparam dof      number of degrees of freedom per node
      //! \param params    parameter list
      //! \param elemat1   elemental matrix
      template <int nnode, int ndim, int dof>
      void EvaluatePTC(Teuchos::ParameterList& params, Epetra_SerialDenseMatrix& elemat1);

      //! get internal (elastic) energy of element
      double GetInternalEnergy() const { return eint_; };

      //! computes the number of different random numbers required in each time step for  generation
      //! of stochastic forces
      int HowManyRandomNumbersINeed();

      inline bool IsParamsInterface() const override { return (not interface_ptr_.is_null()); }

      //! return reference length
      double Lcurr() const { return lcurr_; }

      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

      // TODO: remove once truss3 element is fixed and no longer expects more dofs (6) than it can
      // inherently handle (3)...
      void LocationVector(
          const Discretization& dis, LocationArray& la, bool doDirichlet) const override;

      int NumDofPerElement() const override { return 0; }

      int NumDofPerNode(const DRT::Node& node) const override
      {
        /*note: this is not necessarily the number of DOF assigned to this node by the
         *discretization finally, but only the number of DOF requested for this node by this
         *element; the discretization will finally assign the maximal number of DOF to this node
         *requested by any element connected to this node*/
        return 3;
      }

      int NumLine() const override { return 1; }

      void Pack(DRT::PackBuffer& data) const override;

      Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> ParamsInterfacePtr() override;

      void Print(std::ostream& os) const override;

      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef) override;

      //! scale truss reference length
      void ScaleReferenceLength(double scalefac);

      //! set cross section area of this element
      void SetCrossSec(const double& crosssec);

      void SetParamsInterfacePtr(const Teuchos::ParameterList& p) override;

      //! \brief sets reference coordinates X_ and refenrence length lrefe_ for elements added to
      //! the discretization
      //!
      //! \param xrefe     nodal coordinates in reference frame
      void SetUpReferenceGeometry(const std::vector<double>& xrefe);

      DiscretizationType Shape() const override;

      int UniqueParObjectId() const override { return Truss3Type::Instance().UniqueParObjectId(); }

      void Unpack(const std::vector<char>& data) override;

     protected:
      //! get access to the parameter interface
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

      //! compute damping matrix and statistical forces according to fluctuation dissipation theorem
      //!
      //! \tparam nnode             number of nodes
      //! \tparam ndim              number of dimensions of embedding space
      //! \tparam dof               number of degrees of freedom per node
      //! \tparam randompergauss    number of random numbers required per Gauss point
      //! \param params             parameter list
      //! \param DummyVel           element velocity vector
      //! \param DummyDisp          element disp vector
      //! \param DummyStiffMatrix   element stiffness matrix
      //! \param DummyForce         element internal force vector
      template <int nnode, int ndim, int dof, int randompergauss>
      void CalcBrownian(Teuchos::ParameterList& params, const LINALG::Matrix<1, 6>& DummyVel,
          const LINALG::Matrix<1, 6>& DummyDisp, Epetra_SerialDenseMatrix& DummyStiffMatrix,
          Epetra_SerialDenseVector& DummyForce);

      //! determine Gauss rule from required type of integration
      DRT::UTILS::GaussRule1D MyGaussRule(int nnode, IntegrationType integrationtype);

      //!  \brief shifts nodes so that proper evaluation is possible even in case of periodic
      //!  boundary conditions
      //!
      //! \tparam nnode   number of nodes
      //! \tparam ndim    number of dimensions
      //! \param params   parameter list
      //! \param disp     element disp vector
      template <int nnode, int ndim>
      void NodeShift(Teuchos::ParameterList& params, std::vector<double>& disp);

      //! calculation of elastic energy
      void t3_energy(Teuchos::ParameterList& params, std::vector<double>& disp,
          Epetra_SerialDenseVector* intenergy);

      //! lump mass matrix
      void t3_lumpmass(Epetra_SerialDenseMatrix* emass);

      //! calculation of nonlinear stiffness and mass matrix switching between total lagrange
      //! and enginerring strains
      void t3_nlnstiffmass(Teuchos::ParameterList& params, std::vector<double>& vel,
          std::vector<double>& disp, Epetra_SerialDenseMatrix* stiffmatrix,
          Epetra_SerialDenseMatrix* massmatrix, Epetra_SerialDenseVector* force);

      //! \brief calculate force, nonlinear stiffness and mass matrix using a engineering strain
      //! measure.
      //!
      //! Unlike the fully nonlinear implementation of #t3_nlnstiffmass_totlag, this function uses
      //! \f$\varepsilon=\Delta d / d\f$ as strain measure.
      void t3_nlnstiffmass_engstr(const LINALG::Matrix<1, 6>& DummyDisp,
          Epetra_SerialDenseMatrix& DummyStiffMatrix, Epetra_SerialDenseMatrix* massmatrix,
          Epetra_SerialDenseVector& DummyForce);

      //! calculation of nonlinear stiffness and mass matrix
      void t3_nlnstiffmass_totlag(LINALG::Matrix<1, 6>& DummyDisp,
          Epetra_SerialDenseMatrix& DummyStiffMatrix, Epetra_SerialDenseMatrix* massmatrix,
          Epetra_SerialDenseVector& DummyForce);

      //! cross section area
      double crosssec_;

      //! container containing all kind of Truss3 thisandthat
      DRT::Container data_;

      //! reference tangent position
      LINALG::Matrix<1, 3> diff_disp_ref_;

      //! internal energy of element
      double eint_;

      //! gaussrule_ will be initialized automatically to a 2 point integration rule
      DRT::UTILS::GaussRule1D gaussrule_;

      //!  data exchange between the element and the time integrator.
      Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_;

      //! variable saving whether element has already been initialized (then isinit_ == true)
      bool isinit_;

      //! Vector holding value of Jacobi determinant jacobi for complete integration of massmatrix
      std::vector<double> jacobimass_;

      //! vector holding value of Jacobi determinant jacobi at nodes
      std::vector<double> jacobinode_;

      //! Kinematic type
      KinematicType kintype_;

      //! length in reference configuration
      double lcurr_;

      //! length in reference configuration
      double lrefe_;

      //! material type
      int material_;

      //! reference tangent position
      std::vector<LINALG::Matrix<3, 1>> trefNode_;

      //! reference node position
      LINALG::Matrix<6, 1> X_;

      // don't want = operator
      Truss3& operator=(const Truss3& old);
    };

    // << operator
    std::ostream& operator<<(std::ostream& os, const DRT::Element& ele);

  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // #ifndef TRUSS3_H
