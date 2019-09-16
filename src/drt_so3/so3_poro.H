/*----------------------------------------------------------------------*/
/*! \file

\brief implementation of the 3D solid-poro element

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249

\level 2

*----------------------------------------------------------------------*/
#ifndef SO3_PORO_H
#define SO3_PORO_H

#include "../drt_inpar/inpar_structure.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"

// forward declaration
namespace MAT
{
  class StructPoro;
  class FluidPoro;
  class FluidPoroMultiPhase;
}  // namespace MAT

namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    /*!
    \brief A C++ version of a 3 dimensional solid element with modifications for porous media

    A structural 3 dimensional solid displacement element for large deformations
    and (near)-incompressibility.

    */
    template <class so3_ele, DRT::Element::DiscretizationType distype>
    class So3_Poro : public so3_ele
    {
      //! @name Friends
      friend class So_tet4PoroType;
      friend class So_tet10PoroType;
      friend class So_hex8PoroType;
      friend class So_hex27PoroType;

     public:
      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      So3_Poro(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      So3_Poro(const So3_Poro& old);

      /*!
      \brief Destructor

      */
      virtual ~So3_Poro(){};

      //@}

      //! number of element nodes (
      static const int numnod_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      //! number of space dimensions
      static const int numdim_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      //! number of dofs per node
      static const int noddof_ = numdim_;

      //! total dofs per element
      static const int numdof_ = noddof_ * numnod_;

      //! number of strains per node
      static const int numstr_ = (numdim_ * (numdim_ + 1)) / 2;

      //! number of components necessary to store second derivatives
      /*!
       1 component  for nsd=1:  (N,xx)

       3 components for nsd=2:  (N,xx ; N,yy ; N,xy)

       6 components for nsd=3:  (N,xx ; N,yy ; N,zz ; N,xy ; N,xz ; N,yz)
      */
      static const int numderiv2_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;

      //! total gauss points per element
      int numgpt_;

      //! @name Acess methods

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      virtual int UniqueParObjectId() const;

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

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element

      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Lines();

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element

      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Surfaces();

      /*!
      \brief Get vector of Teuchos::RCPs to the volumes of this element

      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Volumes();

      //! @name Access methods

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      virtual DRT::ElementType& ElementType() const;

      //@}

      /*!
      \brief Return value how expensive it is to evaluate this element

      \param double (out): cost to evaluate this element
      */
      virtual double EvaluationCost() { return 10000.0; }

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate element stiffness, mass, internal forces, etc.

      If NULL on input, the controlling method does not expect the element
      to fill these matrices or vectors.

      \return 0 if successful, negative otherwise
      */
      virtual int Evaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  ///< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la,      ///< location array for de-assembly
          Epetra_SerialDenseMatrix& elemat1,    ///< (stiffness-)matrix to be filled by element.
          Epetra_SerialDenseMatrix& elemat2,    ///< (mass-)matrix to be filled by element.
          Epetra_SerialDenseVector& elevec1,    ///< (internal force-)vector to be filled by element
          Epetra_SerialDenseVector& elevec2,    ///< vector to be filled by element
          Epetra_SerialDenseVector& elevec3     ///< vector to be filled by element
      );

      virtual void PreEvaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  ///< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la       ///< location array for de-assembly
      );


      //! initialize the inverse of the jacobian and its determinant in the material configuration
      virtual void InitElement();

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      /*!
      \brief Query names of element data to be visualized using BINIO

      The element fills the provided map with key names of
      visualization data the element wants to visualize AT THE CENTER
      of the element geometry. The values is supposed to be dimension of the
      data to be visualized. It can either be 1 (scalar), 3 (vector), 6 (sym. tensor)
      or 9 (nonsym. tensor)

      Example:
      \code
        // Name of data is 'Owner', dimension is 1 (scalar value)
        names.insert(std::pair<std::string,int>("Owner",1));
        // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
        names.insert(std::pair<std::string,int>("StressesXYZ",6));
      \endcode

      \param names (out): On return, the derived class has filled names with
                          key names of data it wants to visualize and with int dimensions
                          of that data.
      */
      virtual void VisNames(std::map<std::string, int>& names);

      /*!
      \brief Query data to be visualized using BINIO of a given name

      The method is supposed to call this base method to visualize the owner of
      the element.
      If the derived method recognizes a supported data name, it shall fill it
      with corresponding data.
      If it does NOT recognizes the name, it shall do nothing.

      \warning The method must not change size of data

      \param name (in):   Name of data that is currently processed for visualization
      \param data (out):  data to be filled by element if element recognizes the name
      */
      virtual bool VisData(const std::string& name, std::vector<double>& data);

      /*!
      \brief Read input for this element
      */
      virtual bool ReadElement(const std::string& eletype, const std::string& eledistype,
          DRT::INPUT::LineDefinition* linedef);

      //@}

      //! compute porosity at gausspoint
      virtual void ComputePorosity(Teuchos::ParameterList& params, double press, double J, int gp,
          double& porosity, double* dphi_dp, double* dphi_dJ, double* dphi_dJdp, double* dphi_dJJ,
          double* dphi_dpp, bool save);

      //! compute porosity at gauss point at one face of the element
      virtual void ComputeSurfPorosity(Teuchos::ParameterList& params, double press, double J,
          int surfnum, int gp, double& porosity, double* dphi_dp, double* dphi_dJ,
          double* dphi_dJdp, double* dphi_dJJ, double* dphi_dpp, bool save);

      //! return time derivative of reference porosity
      virtual double RefPorosityTimeDeriv();

      //! evaluate Cauchy stress at given point in parameter space
      virtual void GetCauchyAtXi(const LINALG::Matrix<3, 1>& xi, const std::vector<double>& disp,
          const std::vector<double>& pres, const LINALG::Matrix<3, 1>& n,
          const LINALG::Matrix<3, 1>& t, double& sigma_nt, Epetra_SerialDenseMatrix* dsntdd = NULL,
          Epetra_SerialDenseMatrix* dsntdp = NULL, LINALG::Matrix<3, 1>* dsntdn = NULL,
          LINALG::Matrix<3, 1>* dsntdt = NULL, LINALG::Matrix<3, 1>* dsntdpxi = NULL);

     protected:
      //! don't want = operator
      So3_Poro& operator=(const So3_Poro& old);

      /*!
      \brief Evaluate an element

      Evaluate So3_poro element stiffness, mass, internal forces, etc.
      Templated evaluate routine of element matrixes

      If NULL on input, the controlling method does not expect the element
      to fill these matrices or vectors.

      \return 0 if successful, negative otherwise
      */
      virtual int MyEvaluate(
          Teuchos::ParameterList&
              params,  ///< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  ///< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la,      ///< location array for de-assembly
          Epetra_SerialDenseMatrix& elemat1,    ///< (stiffness-)matrix to be filled by element.
          Epetra_SerialDenseMatrix& elemat2,    ///< (mass-)matrix to be filled by element.
          Epetra_SerialDenseVector& elevec1,    ///< (internal force-)vector to be filled by element
          Epetra_SerialDenseVector& elevec2,    ///< vector to be filled by element
          Epetra_SerialDenseVector& elevec3     ///< vector to be filled by element
      );

      // compute porosity at gausspoint and linearization of porosity w.r.t. structural
      // displacements
      virtual void ComputePorosityAndLinearization(Teuchos::ParameterList& params,
          const double& press, const double& J, const int& gp,
          const LINALG::Matrix<numnod_, 1>& shapfct, const LINALG::Matrix<numnod_, 1>* myporosity,
          const LINALG::Matrix<1, numdof_>& dJ_dus, double& porosity,
          LINALG::Matrix<1, numdof_>& dphi_dus);

      // compute porosity at gausspoint and linearization of porosity w.r.t. fluid pressure
      virtual void ComputePorosityAndLinearizationOD(Teuchos::ParameterList& params,
          const double& press, const double& J, const int& gp,
          const LINALG::Matrix<numnod_, 1>& shapfct, const LINALG::Matrix<numnod_, 1>* myporosity,
          double& porosity, double& dphi_dp);

      //! action parameters recognized by so_hex8
      enum ActionType
      {
        none,
        calc_struct_nlnstiff,
        calc_struct_internalforce,
        calc_struct_nlnstiffmass,
        calc_struct_stress,
        calc_struct_multidofsetcoupling,  //!< structure-fluid coupling: internal force, stiffness
                                          //!< for poroelasticity (structural part)
        calc_struct_poroscatracoupling,
        postprocess_stress,
        interpolate_porosity_to_given_point
      };

      //! container containing all kind of solid3 thisandthat
      DRT::Container data_;
      //! vector of inverses of the jacobian in material frame
      std::vector<LINALG::Matrix<numdim_, numdim_>> invJ_;
      //! determinant of Jacobian in material frame
      std::vector<double> detJ_;
      //! vector of coordinates of current integration point in reference coordinates
      std::vector<LINALG::Matrix<numdim_, 1>> xsi_;

      //! Calculate nonlinear stiffness and internal force for poroelasticity problems
      void nlnstiff_poroelast(std::vector<int>& lm,       ///< location matrix
          LINALG::Matrix<numdim_, numnod_>& disp,         //< current displacements
          LINALG::Matrix<numdim_, numnod_>& vel,          //< current velocities
          LINALG::Matrix<numdim_, numnod_>& evelnp,       //< fluid velocity of element
          LINALG::Matrix<numnod_, 1>& epreaf,             //< fluid pressure of element
          LINALG::Matrix<numdof_, numdof_>* stiffmatrix,  //< element stiffness matrix
          LINALG::Matrix<numdof_, numdof_>* reamatrix,    //< element reactive matrix
          LINALG::Matrix<numdof_, 1>* force,              //< element internal force vector
          Teuchos::ParameterList& params                  //< algorithmic parameters e.g. time
      );

      //! Calculate nonlinear stiffness and internal force for poroelasticity problems (pressure
      //! based formulation)
      virtual void nlnstiff_poroelast_presbased(std::vector<int>& lm,  ///< location matrix
          LINALG::Matrix<numdim_, numnod_>& disp,                      //< current displacements
          const std::vector<double>& ephi,  //< current primary variable for poro-multiphase flow
          LINALG::Matrix<numdof_, numdof_>* stiffmatrix,  //< element stiffness matrix
          LINALG::Matrix<numdof_, 1>* force,              //< element internal force vector
          Teuchos::ParameterList& params                  //< algorithmic parameters e.g. time
      );

      //! Calculate coupling terms in nonlinear stiffness and internal force for poroelasticity
      //! problems
      void coupling_poroelast(std::vector<int>& lm,  ///< location matrix
          LINALG::Matrix<numdim_, numnod_>& disp,    //< current displacements
          LINALG::Matrix<numdim_, numnod_>& vel,     //< current velocities
          LINALG::Matrix<numdim_, numnod_>& evelnp,  //< fluid velocity of element
          LINALG::Matrix<numnod_, 1>& epreaf,        //< fluid pressure of element
          LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>*
              stiffmatrix,                                              //< element stiffness matrix
          LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>* reamatrix,  //< element reactive matrix
          LINALG::Matrix<numdof_, 1>* force,  //< element internal force vector
          Teuchos::ParameterList& params      //< algorithmic parameters e.g. time
      );

      //! Calculate coupling terms in nonlinear stiffness and internal force for poroelasticity
      //! problems (pressure based formulation)
      virtual void coupling_poroelast_presbased(std::vector<int>& lm,  ///< location matrix
          LINALG::Matrix<numdim_, numnod_>& disp,                      //< current displacements
          const std::vector<double>& ephi,     // current primary variable for poro-multiphase flow
          Epetra_SerialDenseMatrix& couplmat,  ///< element stiffness matrix
          Teuchos::ParameterList& params       ///< algorithmic parameters e.g. time
      );

      //! Calculate coupling stress for poroelasticity problems
      virtual void couplstress_poroelast(
          LINALG::Matrix<numdim_, numnod_>& disp,    //< current displacements
          LINALG::Matrix<numdim_, numnod_>& evelnp,  //< fluid velocity of element
          LINALG::Matrix<numnod_, 1>& epreaf,        //< fluid pressure of element
          Epetra_SerialDenseMatrix* elestress,       //< stresses at GP
          Epetra_SerialDenseMatrix* elestrain,       //< strains at GP
          Teuchos::ParameterList& params,            //< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress      //< stress output option
      );

      //! Gauss Point Loop evaluating stiffness and force
      void GaussPointLoop(Teuchos::ParameterList& params,
          const LINALG::Matrix<numdim_, numnod_>& xrefe,
          const LINALG::Matrix<numdim_, numnod_>& xcurr,
          const LINALG::Matrix<numdim_, numnod_>& nodaldisp,
          const LINALG::Matrix<numdim_, numnod_>& nodalvel,
          const LINALG::Matrix<numdim_, numnod_>& evelnp, const LINALG::Matrix<numnod_, 1>& epreaf,
          const LINALG::Matrix<numnod_, 1>* porosity_dof, LINALG::Matrix<numdof_, numdof_>& erea_v,
          LINALG::Matrix<numdof_, numdof_>* stiffmatrix, LINALG::Matrix<numdof_, 1>* force);

      //! Gauss Point Loop evaluating stiffness and force
      void GaussPointLoop_presbased(Teuchos::ParameterList& params,
          const LINALG::Matrix<numdim_, numnod_>& xrefe,
          const LINALG::Matrix<numdim_, numnod_>& xcurr,
          const LINALG::Matrix<numdim_, numnod_>& nodaldisp, const std::vector<double>& ephi,
          LINALG::Matrix<numdof_, numdof_>* stiffmatrix, LINALG::Matrix<numdof_, 1>* force);

      //! Gauss Point Loop evaluating stiffness (off diagonal)
      void GaussPointLoopOD_presbased(Teuchos::ParameterList& params,
          const LINALG::Matrix<numdim_, numnod_>& xrefe,
          const LINALG::Matrix<numdim_, numnod_>& xcurr,
          const LINALG::Matrix<numdim_, numnod_>& nodaldisp, const std::vector<double>& ephi,
          Epetra_SerialDenseMatrix& couplmat);

      //! fill stiffness matrix and rhs vector for darcy flow
      void FillMatrixAndVectors(const int& gp, const LINALG::Matrix<numnod_, 1>& shapefct,
          const LINALG::Matrix<numdim_, numnod_>& N_XYZ, const double& J, const double& press,
          const double& porosity, const LINALG::Matrix<numdim_, 1>& velint,
          const LINALG::Matrix<numdim_, 1>& fvelint,
          const LINALG::Matrix<numdim_, numdim_>& fvelder,
          const LINALG::Matrix<numdim_, numdim_>& defgrd_inv,
          const LINALG::Matrix<numstr_, numdof_>& bop,
          const LINALG::Matrix<numdim_, numdim_>& C_inv,
          const LINALG::Matrix<numdim_, 1>& Finvgradp, const LINALG::Matrix<1, numdof_>& dphi_dus,
          const LINALG::Matrix<1, numdof_>& dJ_dus,
          const LINALG::Matrix<numstr_, numdof_>& dCinv_dus,
          const LINALG::Matrix<numdim_, numdof_>& dFinvdus_gradp,
          const LINALG::Matrix<numdim_ * numdim_, numdof_>& dFinvTdus,
          LINALG::Matrix<numdof_, numdof_>& erea_v, LINALG::Matrix<numdof_, numdof_>* stiffmatrix,
          LINALG::Matrix<numdof_, 1>* force, LINALG::Matrix<numstr_, 1>& fstress);

      //! fill stiffness matrix and rhs vector for brinkman flow
      void FillMatrixAndVectorsBrinkman(const int& gp, const double& J, const double& porosity,
          const LINALG::Matrix<numdim_, numdim_>& fvelder,
          const LINALG::Matrix<numdim_, numdim_>& defgrd_inv,
          const LINALG::Matrix<numstr_, numdof_>& bop,
          const LINALG::Matrix<numdim_, numdim_>& C_inv, const LINALG::Matrix<1, numdof_>& dphi_dus,
          const LINALG::Matrix<1, numdof_>& dJ_dus,
          const LINALG::Matrix<numstr_, numdof_>& dCinv_dus,
          const LINALG::Matrix<numdim_ * numdim_, numdof_>& dFinvTdus,
          LINALG::Matrix<numdof_, numdof_>* stiffmatrix, LINALG::Matrix<numdof_, 1>* force,
          LINALG::Matrix<numstr_, 1>& fstress);

      //! fill stiffness matrix and rhs vector for pressure-based formulation
      void FillMatrixAndVectors_presbased(const int& gp, const LINALG::Matrix<numnod_, 1>& shapefct,
          const LINALG::Matrix<numdim_, numnod_>& N_XYZ, const double& J, const double& press,
          const LINALG::Matrix<numstr_, numdof_>& bop,
          const LINALG::Matrix<numdim_, numdim_>& C_inv, const LINALG::Matrix<1, numdof_>& dJ_dus,
          const LINALG::Matrix<numstr_, numdof_>& dCinv_dus,
          const LINALG::Matrix<1, numdof_>& dps_dus, LINALG::Matrix<numdof_, numdof_>* stiffmatrix,
          LINALG::Matrix<numdof_, 1>* force);

      //! fill stiffness matrix and rhs vector for darcy flow (off diagonal terms)
      void FillMatrixAndVectorsOD(const int& gp, const LINALG::Matrix<numnod_, 1>& shapefct,
          const LINALG::Matrix<numdim_, numnod_>& N_XYZ, const double& J, const double& porosity,
          const double& dphi_dp, const LINALG::Matrix<numdim_, 1>& velint,
          const LINALG::Matrix<numdim_, 1>& fvelint,
          const LINALG::Matrix<numdim_, numdim_>& defgrd_inv,
          const LINALG::Matrix<numdim_, 1>& Gradp, const LINALG::Matrix<numstr_, numdof_>& bop,
          const LINALG::Matrix<numdim_, numdim_>& C_inv,
          LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>* stiffmatrix);

      //! fill stiffness matrix (off diagonal terms) -- pressure-based formulation
      void FillMatrixAndVectorsOD_presbased(const int& gp,
          const LINALG::Matrix<numnod_, 1>& shapefct, const LINALG::Matrix<numdim_, numnod_>& N_XYZ,
          const double& J, const LINALG::Matrix<numstr_, numdof_>& bop,
          const LINALG::Matrix<numdim_, numdim_>& C_inv, const std::vector<double>& solpressderiv,
          Epetra_SerialDenseMatrix& couplmat);

      //! fill stiffness matrix and rhs vector for darcy brinkman flow (off diagonal terms)
      void FillMatrixAndVectorsBrinkmanOD(const int& gp, const LINALG::Matrix<numnod_, 1>& shapefct,
          const LINALG::Matrix<numdim_, numnod_>& N_XYZ, const double& J, const double& porosity,
          const double& dphi_dp, const LINALG::Matrix<numdim_, numdim_>& fvelder,
          const LINALG::Matrix<numdim_, numdim_>& defgrd_inv,
          const LINALG::Matrix<numstr_, numdof_>& bop,
          const LINALG::Matrix<numdim_, numdim_>& C_inv,
          LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>* stiffmatrix);

      //! Gauss Point Loop evaluating stiffness (off diagonal)
      void GaussPointLoopOD(Teuchos::ParameterList& params,
          const LINALG::Matrix<numdim_, numnod_>& xrefe,
          const LINALG::Matrix<numdim_, numnod_>& xcurr,
          const LINALG::Matrix<numdim_, numnod_>& nodaldisp,
          const LINALG::Matrix<numdim_, numnod_>& nodalvel,
          const LINALG::Matrix<numdim_, numnod_>& evelnp, const LINALG::Matrix<numnod_, 1>& epreaf,
          LINALG::Matrix<numdof_, (numdim_ + 1) * numnod_>* stiffmatrix);

      //! helper functions to get element vectors from global vector
      void ExtractValuesFromGlobalVector(
          const DRT::Discretization& discretization,       ///< discretization
          const int& dofset,                               ///< number of dofset
          const std::vector<int>& lm,                      ///< location vector
          LINALG::Matrix<numdim_, numnod_>* matrixtofill,  ///< vector field
          LINALG::Matrix<numnod_, 1>* vectortofill,        ///< scalar field
          const std::string state                          ///< state of the global vector
      );

      // private:

      //! Extrapolate Gauss-point values (e.g. stresses) to nodes
      void stress_expol(Epetra_SerialDenseMatrix& stresses, Epetra_MultiVector& expolstresses);

      //! push forward of material stresses to the current, spatial configuration
      void PK2toCauchy(LINALG::Matrix<numstr_, 1>& stress, LINALG::Matrix<numdim_, numdim_>& defgrd,
          LINALG::Matrix<numdim_, numdim_>& cauchystress);

      //! get materials (solid and fluid)
      void GetMaterials();

      //! get materials (solid and fluid) for pressure based formulation
      void GetMaterials_presbased();

      //! Compute Solid-pressure derivative w.r.t. primary variable at GP
      void ComputeSolPressureDeriv(const std::vector<double>& phiAtGP,  ///<< primary variable
          const int numfluidphases,                                     ///<< number of fluid phases
          std::vector<double>& solidpressderiv  ///<< solid pressure derivative at GP
      );

      //! Compute solid pressure at GP
      double ComputeSolPressureAtGP(
          const int totalnumdofpernode,       ///<< total number of multiphase dofs
          const int numfluidphases,           ///<< number of fluid phases
          const std::vector<double>& phiAtGP  ///<< primary variable
      );

      //! recalculate solid pressure at GP in case of volfracs
      double RecalculateSolPressureAtGP(double press, const double porosity,
          const int totalnumdofpernode, const int numfluidphases, const int numvolfrac,
          const std::vector<double>& phiAtGP);

      //! recalculate solid pressure derivative in case of volfracs
      void RecalculateSolPressureDeriv(const std::vector<double>& phiAtGP,
          const int totalnumdofpernode, const int numfluidphases, const int numvolfrac,
          const double press, const double porosity, std::vector<double>& solidpressderiv);

      //! Compute primary variable for multiphase flow at GP
      void ComputePrimaryVariableAtGP(
          const std::vector<double>& ephi,             ///<< primary variable at node
          const int totalnumdofpernode,                ///<< total number of multiphase dofs
          const LINALG::Matrix<numnod_, 1>& shapefct,  ///<< shapefct
          std::vector<double>& phiAtGP                 ///<< primary variable at GP
      );

      //! Compute linearizaton of solid press w.r.t. displacements
      //! only needed if additional volume fractions are present and porosity depends on
      //! Jacobian of deformation gradient
      void ComputeLinearizationOfSolPressWrtDisp(const double fluidpress, const double porosity,
          const int totalnumdofpernode, const int numfluidphases, const int numvolfrac,
          const std::vector<double>& phiAtGP, const LINALG::Matrix<1, numdof_>& dphi_dus,
          LINALG::Matrix<1, numdof_>& dps_dus);

      //! evaluate shape functions and their derivatives at gauss point
      void ComputeShapeFunctionsAndDerivatives(const int& gp, LINALG::Matrix<numnod_, 1>& shapefct,
          LINALG::Matrix<numdim_, numnod_>& deriv, LINALG::Matrix<numdim_, numnod_>& N_XYZ);

      //! Compute Jacobian Determinant and the volume change
      void ComputeJacobianDeterminantVolumeChange(double& J, double& volchange,
          const LINALG::Matrix<numdim_, numdim_>& defgrd,
          const LINALG::Matrix<numdim_, numnod_>& N_XYZ,
          const LINALG::Matrix<numdim_, numnod_>& nodaldisp);

      //! Compute Jacobian Determinant and the volume change and its linearizations
      virtual void ComputeJacobianDeterminantVolumeChangeAndLinearizations(
          double& J,          ///< (o) Jacobian Determinant
          double& volchange,  ///< (o) Change of Volume
          LINALG::Matrix<1, numdof_>&
              dJ_dus,  ///< (o) Linearization of Jacobian Determinant  w.r.t displacements
          LINALG::Matrix<1, numdof_>&
              dvolchange_dus,  ///<  (o) Linearization of Change of Volume  w.r.t displacements
          const LINALG::Matrix<numdim_, numdim_>& defgrd,  ///< (i) deformation gradient
          const LINALG::Matrix<numdim_, numdim_>&
              defgrd_inv,  ///< (i) inverse of deformation gradient
          const LINALG::Matrix<numdim_, numnod_>&
              N_XYZ,  ///< (i) spatial gradient of material shape functions
          const LINALG::Matrix<numdim_, numnod_>& nodaldisp  ///<  (i) nodal displacements
      );

      //! Compute Linearization Of Jacobian
      void ComputeLinearizationOfJacobian(LINALG::Matrix<1, numdof_>& dJ_dus, const double& J,
          const LINALG::Matrix<numdim_, numnod_>& N_XYZ,
          const LINALG::Matrix<numdim_, numdim_>& defgrd_inv);

      //! Compute Linearization some Auxiliary Values in gauss point loop
      void ComputeAuxiliaryValues(const LINALG::Matrix<numdim_, numnod_>& N_XYZ,
          const LINALG::Matrix<numdim_, numdim_>& defgrd_inv,
          const LINALG::Matrix<numdim_, numdim_>& C_inv, const LINALG::Matrix<numdim_, 1>& Gradp,
          LINALG::Matrix<numdim_ * numdim_, numdof_>& dFinvTdus,
          LINALG::Matrix<numdim_, 1>& Finvgradp, LINALG::Matrix<numdim_, numdof_>& dFinvdus_gradp,
          LINALG::Matrix<numstr_, numdof_>& dCinv_dus);

      //! Compute  nonlinear b-operator
      void ComputeBOperator(LINALG::Matrix<numstr_, numdof_>& bop,
          const LINALG::Matrix<numdim_, numdim_>& defgrd,
          const LINALG::Matrix<numdim_, numnod_>& N_XYZ);

      //! Compute deformation gradient
      void ComputeDefGradient(
          LINALG::Matrix<numdim_, numdim_>& defgrd,  ///<<    (i) deformation gradient at gausspoint
          const LINALG::Matrix<numdim_, numnod_>&
              N_XYZ,  ///<<    (i) derivatives of shape functions w.r.t. reference coordinates
          const LINALG::Matrix<numdim_, numnod_>& xcurr  ///<<    (i) current position of gausspoint
      );

      //! Gauss integration rule
      DRT::UTILS::GaussIntegration intpoints_;

      //! flag indicating higher order element
      // bool ishigherorder_;

      //! flag indicating if element has been initialized
      bool init_;

      //! flag for scatra coupling
      bool scatracoupling_;

      //! flag for nurbs
      bool isNurbs_;

      //! weights for nurbs elements
      LINALG::Matrix<numnod_, 1> weights_;
      //! knot vector for nurbs elements
      std::vector<Epetra_SerialDenseVector> myknots_;

      //! corresponding fluid material
      Teuchos::RCP<MAT::FluidPoro> fluidmat_;

      //! corresponding multiphase fluid material
      Teuchos::RCP<MAT::FluidPoroMultiPhase> fluidmultimat_;

      //! own poro structure material
      Teuchos::RCP<MAT::StructPoro> structmat_;

      // int numscal_;

      //! get nodes of element
      DRT::Node** Nodes();

      //! get material of element
      Teuchos::RCP<MAT::Material> Material() const;

      //! get global id of element
      int Id() const;
    };  // class So3_Poro


    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================


  }  // namespace ELEMENTS
}  // namespace DRT


#endif /* SO3_PORO_H */
