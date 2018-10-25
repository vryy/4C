/*!----------------------------------------------------------------------
\file truss3cl.H

\brief three dimensional interpolated total Lagrange hybrid beam-truss element
 (can be connected to beam3 elements)

\level 3

\maintainer Dhrubajyoti Mukherjee
            mukherjee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15270
*----------------------------------------------------------------------*/
// administrative line in order to make a difference between BACI and Ccarat files

// header file only included if not yet included!
#ifndef TRUSS3CL_H
#define TRUSS3CL_H


#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"

#include "../drt_lib/drt_elementtype.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_inpar/inpar_browniandyn.H"


//#define mygaussruleeb intrule_line_6point //define gauss rule; intrule_line_1point -
// intrule_line_10point is implemented

namespace DRT
{
  namespace ELEMENTS
  {
    class Truss3CLType : public DRT::ElementType
    {
     public:
      std::string Name() const { return "Truss3CLType"; }

      static Truss3CLType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

      virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
          const double* x0, int numdf, int dimns);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static Truss3CLType instance_;
    };

    /*!
    \brief three dimensional total Lagrange truss element

    \author cyron (cyron @lnm.mw.tum.de)
    */
    class Truss3CL : public DRT::Element
    {
     public:
      //! @name Friends
      friend class Truss3CLType;

      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id    (in): A globally unique element id
      \param etype (in): Type of element
      \param owner (in): owner processor of the element
      */
      Truss3CL(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element
      */
      Truss3CL(const Truss3CL& old);

      /*!
      \brief Destructor
      */
      virtual ~Truss3CL();

      /*!
      \brief Deep copy this instance of Truss3CL and return pointer to the copy

      The Clone() method is used by the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed
    .
      */
      virtual DRT::Element* Clone() const;

      /*!
      \brief Return the location vector of this element
      */
      // TODO: remove once truss3cl element is fixed and no longer expects more dofs (6) than it can
      // inherently handle (3)...
      virtual void LocationVector(
          const Discretization& dis, LocationArray& la, bool doDirichlet) const;

      /*!
     \brief Get shape type of element
     */
      virtual DiscretizationType Shape() const;


      /*!
      \brief Return unique ParObject id

      Every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H
      */
      virtual int UniqueParObjectId() const { return Truss3CLType::Instance().UniqueParObjectId(); }

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

      virtual DRT::ElementType& ElementType() const { return Truss3CLType::Instance(); }

      //@}

      //! @name methods for initilization of the element

      //! sets reference coordinates X_ and refenrence length lrefe_ for elements added to the
      //! discretization after simulation start
      void SetUpReferenceGeometry(const std::vector<double>& xrefe, const bool secondinit = false);


      //@}

      /*!
      \brief Return number of lines to this element
      */
      virtual int NumLine() const { return 1; }


      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element
      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Lines();


      /*!
      \brief Get number of degrees of freedom of a single node
      */
      virtual int NumDofPerNode(const DRT::Node& node) const
      {
        /*note: this is not necessarily the number of DOF assigned to this node by the
         *discretization finally, but only the number of DOF requested for this node by this
         *element; the discretization will finally assign the maximal number of DOF to this node
         *requested by any element connected to this node*/
        return 3;
      }


      /*!
      \brief Get number of degrees of freedom per element not including nodal degrees of freedom
      */
      virtual int NumDofPerElement() const { return 0; }

      /*!
      \brief Set cross section area of this element
      */
      void SetCrossSec(const double& crosssec);

      /*!
      \brief computes the number of different random numbers required in each time step for
      generation of stochastic forces
      */
      int HowManyRandomNumbersINeed();

      /*!
      \brief Print this element
      */
      virtual void Print(std::ostream& os) const;

      //! \brief return reference Coordinates of interpolated nodes
      const std::vector<double>& XRef() const { return xrefe_; }


      //@}

      //! @name Construction

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
      virtual int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);

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
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1 = NULL);


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
      int EvaluatePTC(Teuchos::ParameterList& params, Epetra_SerialDenseMatrix& elemat1);

      //!\brief return cross section of the truss
      const double& CSec() const { return crosssec_; }

      //! \brief return reference length
      const double& L0() const { return lrefe_; }

      //! \brief internal force vector
      const Teuchos::RCP<Epetra_SerialDenseVector> InternalForceVector() const { return f_; }

      //@}
      const std::vector<double>& MyBindingPosition() const { return mybindingposition_; }

      /*!
       \brief Set (Internodal) Binding Position of this element
       */
      void SetBindingPosition(const double& pos1, const double& pos2);


     private:
      //! action parameters recognized by Truss3CL
      enum ActionType
      {
        calc_none,
        calc_struct_linstiff,
        calc_struct_nlnstiff,
        calc_struct_nlnstiff2,
        calc_struct_internalforce,
        calc_struct_linstiffmass,
        calc_struct_nlnstiffmass,
        calc_struct_nlnstifflmass,  //!< internal force, its stiffness and lumped mass matrix
        calc_struct_stress,
        calc_struct_eleload,
        calc_struct_fsiload,
        calc_struct_update_istep,
        calc_struct_reset_istep,  //!< reset elementwise internal variables
                                  //!< during iteration to last converged state
        postprocess_stress,
        calc_struct_ptcstiff,  // calculate artificial stiffness due to PTC solution strategy
        calc_struct_energy     // calculate current elastic energy of element
      };

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

      //! container containing all kind of Truss3CL thisandthat
      DRT::Container data_;

      //! variable saving whether element has already been initialized (then isinit_ == true)
      bool isinit_;

      //! reference node position
      LINALG::Matrix<6, 1> X_;

      //! material type
      int material_;

      //! length in reference configuration
      double lrefe_;

      //! reference coordinates of the interpolated nodes
      std::vector<double> xrefe_;

      //! Vector holding value of Jacobi determinant jacobi for complete integration of massmatrix
      std::vector<double> jacobimass_;

      //! vector holding value of Jacobi determinant jacobi at nodes
      std::vector<double> jacobinode_;

      //! cross section area
      double crosssec_;

      //! Vector holding value of Jacobi determinant jacobi for each Gauss point for
      //! underintegration
      std::vector<double> jacobi_;

      //! internal force vector
      Teuchos::RCP<Epetra_SerialDenseVector> f_;

      //! Kinematic type
      KinematicType kintype_;

      //! gaussrule_ will be initialized automatically to a 2 point integration rule
      DRT::UTILS::GaussRule1D gaussrule_;

      //! Position of the two interpolated nodes relativ to real nodes
      std::vector<double> mybindingposition_;


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
      void t3_nlnstiffmass_totlag(LINALG::Matrix<6, 1>& fdisp,
          Epetra_SerialDenseMatrix& fstiffmatrix, Epetra_SerialDenseMatrix* massmatrix,
          Epetra_SerialDenseVector& fforce);

      /*!
       * \brief calculate force, nonlinear stiffness and mass matrix using a engineering strain
       * measure.
       *
       * Unlike the fully nonlinear implementation of #t3_nlnstiffmass_totlag, this function uses
       * \f$\varepsilon=\Delta d / d\f$ as strain measure.
       *
       */
      void t3_nlnstiffmass_engstr(const LINALG::Matrix<6, 1>& fdisp,
          Epetra_SerialDenseMatrix& fstiffmatrix, Epetra_SerialDenseMatrix* massmatrix,
          Epetra_SerialDenseVector& fforce);

      //! lump mass matrix
      void t3_lumpmass(Epetra_SerialDenseMatrix* emass);

      //@}

      //! @name Methods for Brownian dynamics simulations

      //! compute damping matrix and statistical forces according to fluctuation dissipation theorem
      template <int fnnode, int ndim, int dof,
          int randompergauss>  // number of nodes, number of dimensions of embedding space, number
                               // of degrees of freedom per node, number of random numbers required
                               // per Gauss point
      void
      CalcBrownian(Teuchos::ParameterList& params,  //!< parameter list
          const LINALG::Matrix<6, 1>& fvel,         //!< element velocity vector
          const std::vector<double>& disp,          //!< element real disp vector
          const LINALG::Matrix<6, 1>& fdisp,        //!< element disp vector
          Epetra_SerialDenseMatrix& fstiffmatrix,   //!< element stiffness matrix
          Epetra_SerialDenseVector& fforce);        //!< element internal force vector

      //! computes damping coefficients for employed friction model
      void MyDampingConstants(Teuchos::ParameterList& params, LINALG::Matrix<3, 1>& gamma);

      //! computes translational damping forces and stiffness
      template <int fnnode, int ndim,
          int dof>  // number of nodes, number of dimensions of embedding space, number of degrees
                    // of freedom per node
      void
      MyTranslationalDamping(Teuchos::ParameterList& params,  //!< parameter list
          const LINALG::Matrix<6, 1>& fvel,                   //!< element velocity vector
          const std::vector<double>& disp,                    //!< element real disp vector
          const LINALG::Matrix<6, 1>& fdisp,                  //!< element fictitious disp vector
          Epetra_SerialDenseMatrix& fstiffmatrix,             //!< element stiffness matrix
          Epetra_SerialDenseVector& fforce);                  //!< element internal force vector

      //! computes stochastic translational forces and resulting stiffness
      template <int fnnode, int ndim, int dof,
          int randompergauss>  // number of nodes, number of dimensions of embedding space, number
                               // of degrees of freedom per node, number of random numbers required
                               // per Gauss point
      void
      MyStochasticForces(Teuchos::ParameterList& params,  //!< parameter list
          const LINALG::Matrix<6, 1>& fdisp,              //!< element disp vector
          Epetra_SerialDenseMatrix& fstiffmatrix,         //!< element stiffness matrix
          Epetra_SerialDenseVector& fforce);              //!< element internal force vector

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

      // computes reference DOFs for TRUSS3CL at t=0, i.e. before first time step (public)
      void ComputeRefDOFs(LINALG::Matrix<6, 1>& Xstart);

      //@}

      // don't want = operator
      Truss3CL& operator=(const Truss3CL& old);


    };  // class Truss3CL



    // << operator
    std::ostream& operator<<(std::ostream& os, const DRT::Element& ele);


  }  // namespace ELEMENTS
}  // namespace DRT



#endif  // #ifndef Truss3CL_H
