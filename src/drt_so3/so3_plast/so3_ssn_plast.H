/*----------------------------------------------------------------------*/
/*! \file
\level 2
*/


/*----------------------------------------------------------------------*
 | definitions                                              seitz 07/13 |
 *----------------------------------------------------------------------*/
#ifndef SO3_SSN_PLAST_H
#define SO3_SSN_PLAST_H

/*----------------------------------------------------------------------*
 | headers                                                  seitz 07/13 |
 *----------------------------------------------------------------------*/
#include "so3_ssn_plast_eletypes.H"
#include "../../drt_thermo/thermo_ele_impl_utils.H"

#include "../../linalg/linalg_serialdensematrix.H"
#include "../../linalg/linalg_serialdensevector.H"

// include header because number of Gauss points are determined
// dependent on distype
#include "../../drt_fem_general/drt_utils_gausspoints.H"

#include "../../drt_fem_general/drt_utils_nurbs_shapefunctions.H"

/*----------------------------------------------------------------------*
 |                                                          seitz 07/13 |
 *----------------------------------------------------------------------*/
namespace DRT
{
  class Discretization;

  namespace UTILS
  {
    struct PlastSsnData;
  }
  namespace ELEMENTS
  {
    /*!
     * \brief EAS technology enhancement types of so_hex8
     *
     * Solid Hex8 has EAS enhancement of GL-strains to avoid locking.
     * New enum for so3_plast. Currently only for hex8
     */
    enum So3Plast_EASType  // with meaningful value for matrix size info
    {
      soh8p_easnone,   //!< no EAS i.e. displacement based with tremendous locking
      soh8p_eassosh8,  //!< related to Solid-Shell, 7 parameters to alleviate
                       //!< in-plane (membrane) locking and main modes for Poisson-locking
      soh8p_easmild,   //!< 9 parameters consisting of modes to alleviate
                       //!< shear locking (bending) and main incompressibility modes
                       //!< (for Solid Hex8)
                       //!< The sosh18 also uses 9 eas parameters. Therefore, we can
                       //!< re-use this EASType here.
      soh8p_easfull,   //!< 21 parameters to prevent almost all locking modes.
                       //!< Equivalent to all 30 parameters to fully complete element
                       //!< with quadratic modes (see Andelfinger 1993 for details) and therefore
                       //!< also suitable for distorted elements.
                       //!< (for Solid Hex8)
      soh18p_eassosh18
    };

    //! A C++ version of a 3 dimensional solid element with modifications von Mises
    //! plasticity using a semi-smooth Newton method
    template <DRT::Element::DiscretizationType distype>
    class So3_Plast : public virtual So_base
    {
      //! @name Friends
      friend class So_hex8PlastType;
      friend class So_hex27PlastType;

     public:
      //@}
      //! @name Constructors and destructors and related methods


      //! Standard Constructor
      So3_Plast(int id,  //!< (i) this element's global id
          int owner      //!< elements owner
      );

      //! Copy Constructor
      //! Makes a deep copy of a Element
      So3_Plast(const So3_Plast& old);

      //! Deep copy this instance of Solid3 and return pointer to the copy
      //!
      //! The Clone() method is used from the virtual base class Element in cases
      //! where the type of the derived class is unknown and a copy-ctor is needed
      DRT::Element* Clone() const;

      //! Destructor
      virtual ~So3_Plast(){};

      //@}

      //! number of element nodes
      static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
      //! number of space dimensions
      static const int nsd_ = 3;
      //! number of dofs per node
      static const int numdofpernode_ = 3;
      //! total dofs per element
      static const int numdofperelement_ = numdofpernode_ * nen_;
      //! number of strains/stresses
      static const int numstr_ = 6;
      //! number of Gauss points per element (value is added in so3_thermo.cpp)
      int numgpt_;
      //! static const is required for fixedsizematrices
      //! TODO maybe more beauty is possible
      static const int numgpt_post = THR::DisTypeToSTRNumGaussPoints<distype>::nquad;


      //! @name Acess methods

      //! Return unique ParObject id
      //!
      //! every class implementing ParObject needs a unique id defined at the top of
      //! this file.
      virtual int UniqueParObjectId() const;

      bool HaveEAS() const { return (eastype_ != soh8p_easnone); };

      //! Pack this class so it can be communicated
      //! Pack and \ref Unpack are used to communicate this element
      virtual void Pack(DRT::PackBuffer& data) const;

      //! Unpack data from a char vector into this class
      //! Pack and \ref Unpack are used to communicate this element
      virtual void Unpack(const std::vector<char>& data);

      //! Get number of degrees of freedom of a certain node in case of multiple
      //! dofsets (implements pure virtual DRT::Element)
      //!
      //! The element decides how many degrees of freedom its nodes must have.
      virtual int NumDofPerNode(const DRT::Node& node) const { return nsd_; };

      //! Get number of degrees of freedom of this element
      virtual int NumDofPerElement() const { return 0; };

      //@}

      //! @name Access methods

      //! Print this element
      void Print(std::ostream& os) const;

      //! return elementtype
      virtual DRT::ElementType& ElementType() const;

      //! return element shape
      virtual DiscretizationType Shape() const { return distype; };

      /*!
      \brief Return number of volumes of this element
      */
      virtual int NumVolume() const;

      /*!
      \brief Return number of surfaces of this element
      */
      virtual int NumSurface() const;

      /*!
      \brief Return number of lines of this element
      */
      virtual int NumLine() const;

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
      //@}

      //! @name Input and Creation

      //! Query names of element data to be visualized using BINIO
      //!
      //! The element fills the provided map with key names of
      //! visualization data the element wants to visualize AT THE CENTER
      //! of the element geometry. The values is supposed to be dimension of the
      //! data to be visualized. It can either be 1 (scalar), 3 (vector), 6 (sym. tensor)
      //! or 9 (nonsym. tensor)
      //!
      //! Example:
      //! \code
      //!   // Name of data is 'Owner', dimension is 1 (scalar value)
      //!   names.insert(std::pair<string,int>("Owner",1));
      //!   // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
      //!   names.insert(std::pair<string,int>("StressesXYZ",6));
      //! \endcode
      //!
      //!  names (out): On return, the derived class has filled names with key
      //!               names of data it wants to visualize and with int dimensions
      //!               of that data.
      virtual void VisNames(std::map<std::string, int>& names);

      //!  Query data to be visualized using BINIO of a given name
      //!
      //! The method is supposed to call this base method to visualize the owner of
      //! the element.
      //! If the derived method recognizes a supported data name, it shall fill it
      //! with corresponding data.
      //! If it does NOT recognizes the name, it shall do nothing.
      //!
      //! warning: the method must not change size of data
      //!
      //!  name (in):   Name of data that is currently processed for visualization
      //! \param data (out):  data to be filled by element if element recognizes the name
      virtual bool VisData(const std::string& name, std::vector<double>& data);

      //! read input for this element
      bool ReadElement(const std::string& eletype,  //!< so3plast(fbar)
          const std::string& eledistype,            //!< hex8,tet4,...
          DRT::INPUT::LineDefinition* linedef       //!< what parameters have to be read
      );

      //@}

      //! @name Evaluation

      //! evaluate an element
      //! evaluate element stiffness, mass, internal forces, etc.
      //!
      //! if NULL on input, the controlling method does not expect the element
      //!  to fill these matrices or vectors.
      //!
      //!  \return 0 if successful, negative otherwise
      int Evaluate(
          Teuchos::ParameterList&
              params,  //!< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  //!< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la,      //!< location array for de-assembly
          Epetra_SerialDenseMatrix&
              elemat1_epetra,  //!< (stiffness-)matrix to be filled by element.
          Epetra_SerialDenseMatrix& elemat2_epetra,  //!< (mass-)matrix to be filled by element.
          Epetra_SerialDenseVector&
              elevec1_epetra,  //!< (internal force-)vector to be filled by element
          Epetra_SerialDenseVector& elevec2_epetra,  //!< vector to be filled by element
          Epetra_SerialDenseVector& elevec3_epetra   //!< vector to be filled by element
      );

      //@}

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a surface Neumann condition on the solid3 element

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param condition (in)     : The condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If NULL on input,

      \return 0 if successful, negative otherwise
      */
      virtual int EvaluateNeumann(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1 = NULL);

      //! init the inverse of the jacobian and its determinant in the material
      //! configuration
      virtual void InitJacobianMapping();

      // get parameter list from ssn_plast_manager
      virtual void ReadParameterList(Teuchos::RCP<Teuchos::ParameterList> plparams);

      /**
       * \brief Evaluate Cauchy stress at given point in parameter space and calculate
       * linearizations
       *
       * \param xi[in]          position in parameter space xi
       * \param disp[in]        vector of displacements
       * \param n[in]           vector n (\f[\bold{n}\f])
       * \param t[in]           vector t (\f[\bold{t}\f])
       * \param sigma_nt[out]   cauchy stress tensor contracted using the vectors n and t
                                (\f[ \bold{\sigma} \cdot \bold{n} \cdot \bold{t} \f])
       * \param DsntDd[out]     derivative of sigma_nt w.r.t. displacements (\f[ \frac{ \mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{d}}
                                \f])
       * \param D2sntDd2[out]   second derivative of sigma_nt w.r.t. displacements (\f[ \frac{
                                \mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}
                                {\mathrm{d} \bold{d}^2} \f])
       * \param D2sntDdDn[out]  second derivative of sigma_nt w.r.t. displacements and vector n (\f[
                                \frac{\mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}
                                {\mathrm{d} \bold{d} \mathrm{d} \bold{n} } \f])
       * \param D2sntDdDt[out]  second derivative of sigma_nt w.r.t. displacements and vector t (\f[
                                \frac{\mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}
                                {\mathrm{d} \bold{d} \mathrm{d} \bold{t} } \f])
       * \param D2sntDdDxi[out] second derivative of sigma_nt w.r.t. displacements and local
                                parameter coordinate xi (\f[ \frac{\mathrm{d}^2 \bold{\sigma} \cdot
                                \bold{n} \cdot \bold{t}} {\mathrm{d} \bold{d} \mathrm{d} \bold{\xi}
                                } \f])
       * \param DsntDn[out]     derivative of sigma_nt w.r.t. vector n (\f[ \frac{ \mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{n}}
                                \f])
       * \param DsntDt[out]     derivative of sigma_nt w.r.t. vector n (\f[ \frac{ \mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{t}}
                                \f])
       * \param DsntDxi[out]    derivative of sigma_nt w.r.t. local parameter coordinate xi (\f[
                                \frac{ \mathrm{d} \bold{\sigma} \cdot \bold{n} \cdot
                                \bold{t}}{\mathrm{d} \bold{\xi}} \f])
       * \param temp[in]        temperature
       * \param DsntDT[out]     derivative of sigma_nt w.r.t. temperature (\f[ \frac{ \mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} T} \f])
       * \param D2sntDdDT[out]  second derivative of sigma_nt w.r.t. displacements and temperature
                                (\f[ \frac{\mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}
                                {\mathrm{d} \bold{d} \mathrm{d} T } \f])
       *
       * \note At the moment this method is only used for the nitsche contact formulation
       */
      void GetCauchyAtXi(const LINALG::Matrix<3, 1>& xi, const std::vector<double>& disp,
          const LINALG::Matrix<3, 1>& n, const LINALG::Matrix<3, 1>& t, double& sigma_nt,
          Epetra_SerialDenseMatrix* DsntDd, Epetra_SerialDenseMatrix* D2sntDd2,
          Epetra_SerialDenseMatrix* D2sntDdDn, Epetra_SerialDenseMatrix* D2sntDdDt,
          Epetra_SerialDenseMatrix* D2sntDdDxi, LINALG::Matrix<3, 1>* DsntDn,
          LINALG::Matrix<3, 1>* DsntDt, LINALG::Matrix<3, 1>* DsntDxi,
          const std::vector<double>* temp = NULL, Epetra_SerialDenseMatrix* DsntDT = NULL,
          Epetra_SerialDenseMatrix* D2sntDdDT = NULL) override;

      void HeatFlux(const std::vector<double>& temp, const std::vector<double>& disp,
          const LINALG::Matrix<nsd_, 1>& xi, const LINALG::Matrix<nsd_, 1>& n, double& q,
          Epetra_SerialDenseMatrix* dq_dT, Epetra_SerialDenseMatrix* dq_dd,
          LINALG::Matrix<nsd_, 1>* dq_dn, LINALG::Matrix<nsd_, 1>* dq_dpxi,
          Epetra_SerialDenseMatrix* d2q_dT_dd, Epetra_SerialDenseMatrix* d2q_dT_dn,
          Epetra_SerialDenseMatrix* d2q_dT_dpxi);

      void HeatFlux(const std::vector<double>& temp, const std::vector<double>& disp,
          const LINALG::Matrix<2, 1>& xi, const LINALG::Matrix<2, 1>& n, double& q,
          Epetra_SerialDenseMatrix* dq_dT, Epetra_SerialDenseMatrix* dq_dd,
          LINALG::Matrix<2, 1>* dq_dn, LINALG::Matrix<2, 1>* dq_dpxi,
          Epetra_SerialDenseMatrix* d2q_dT_dd, Epetra_SerialDenseMatrix* d2q_dT_dn,
          Epetra_SerialDenseMatrix* d2q_dT_dpxi)
      {
        dserror("wrong spatial dimension");
      }


      virtual void SetIsNitscheContactEle(bool val)
      {
        is_nitsche_contact_ = val;
        if (is_nitsche_contact_)
        {
          cauchy_.resize(numgpt_);
          cauchy_deriv_.resize(numgpt_);
          if (tsi_) cauchy_deriv_T_.resize(numgpt_);
        }
      }
      //@}

     protected:
      //! don't want = operator
      So3_Plast& operator=(const So3_Plast& old);

      //! number of plastic variables at each gauss point
      enum PlSpinType
      {
        zerospin = 5,  //!< 5 parameters for zero plastic spin (symmetric traceless tensor)
        plspin = 8,    //!< 8 parameters for zero plastic spin (non-symmetric traceless tensor)
      };

      //! vector of coordinates of current integration point in reference coordinates
      std::vector<LINALG::Matrix<nsd_, 1>> xsi_;
      //! Gauss point weights
      std::vector<double> wgt_;

      //! @name plasticity related stuff


      //! Calculate nonlinear stiffness and mass matrix with condensed plastic matrices
      virtual void nln_stiffmass(std::vector<double>& disp,  // current displacements
          std::vector<double>& vel,                          // current velocities
          std::vector<double>& temp,                         // current temperatures
          LINALG::Matrix<numdofperelement_, numdofperelement_>*
              stiffmatrix,  // element stiffness matrix
          LINALG::Matrix<numdofperelement_, numdofperelement_>* massmatrix,  // element mass matrix
          LINALG::Matrix<numdofperelement_, 1>* force,      // element internal force vector
          LINALG::Matrix<numgpt_post, numstr_>* elestress,  // stresses at GP
          LINALG::Matrix<numgpt_post, numstr_>* elestrain,  // strains at GP
          Teuchos::ParameterList& params,                   // algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,            // stress output option
          const INPAR::STR::StrainType iostrain             // strain output option
      );

      //! Calculate the coupling matrix K_dT for monolithic TSI
      virtual void nln_kdT_tsi(LINALG::Matrix<numdofperelement_, nen_>*
                                   k_dT,  //!< (o): mechanical thermal stiffness term at current gp
          Teuchos::ParameterList& params);

      // Add plastic increment of converged state to plastic history for nonlinear kinematics
      virtual void UpdatePlasticDeformation_nln(PlSpinType spintype);

      //@}

      //! calculate nonlinear B-operator
      virtual void CalculateBop(LINALG::Matrix<numstr_, numdofperelement_>* bop,
          const LINALG::Matrix<nsd_, nsd_>* defgrd, const LINALG::Matrix<nsd_, nen_>* N_XYZ,
          const int gp);

      /// Extrapolate Gauss-point values (e.g. stresses) to nodes and store results in elevectors
      /// todo: Unfortunately, there is no universal extrapolation function in all base elements,
      /// i.e. the hex8 elements would call soh8_expol, hex27 would call soh27_expol...
      /// For now, we copy the extrapolation function here and allow the use only for hex8 elements
      void soh8_expol(LINALG::Matrix<numgpt_post, numstr_>& stresses,  ///< gp stresses
          Epetra_MultiVector& expolstresses);                          ///< nodal stresses

      // EAS element techonolgy *******************************************
      //! Initialize data for EAS (once)
      virtual void EasInit();

      //! setup EAS for each evaluation
      virtual void EasSetup();
      ;

      // evaluate EAS shape functions
      virtual void EasShape(const int gp);

      //! add eas strains to GL strains
      virtual void EasEnhanceStrains();


      //! @name plasticity related stuff

      /*! \brief Calculate the deformation gradient that is consistent
       *         with modified (e.g. EAS) GL strain tensor.
       *         Expensive (two polar decomposition), but required, if
       *         the material evaluation is based on the deformation
       *         gradient rather than the GL strain tensor (e.g. plasticity).
       *
       * \param defgrd_disp  (in)  : displacement-based deformation gradient
       * \param glstrain_mod (in)  : modified GL strain tensor (strain-like Voigt notation)
       * \param defgrd_mod   (out) : consistent modified deformation gradient
       */
      virtual void CalcConsistentDefgrd();


      /*! \brief Evaluate the NCP function and the linearization
       *         and condense the additional degrees of freedom
       *         into the stiffness matrix block.
       *         Required data from the element evaluate is handed in,
       *         with potential EAS matrices or F-bar linearizations.
       *         The template argument "spinype" decides, if
       *         the additional evolution equation for the plastic
       *         spin is linearized and solved for.
       */
      template <int spintype>
      void CondensePlasticity(const LINALG::Matrix<nsd_, nsd_>& defgrd,
          const LINALG::Matrix<nsd_, nsd_>& deltaLp,
          const LINALG::Matrix<numstr_, numdofperelement_>& bop,
          const LINALG::Matrix<nsd_, nen_>* N_XYZ, const LINALG::Matrix<numstr_, 1>* RCG,
          const double detJ_w, const int gp, const double temp, Teuchos::ParameterList& params,
          LINALG::Matrix<numdofperelement_, 1>* force,
          LINALG::Matrix<numdofperelement_, numdofperelement_>* stiffmatrix,
          const Epetra_SerialDenseMatrix* M = NULL, Epetra_SerialDenseMatrix* Kda = NULL,
          std::vector<Epetra_SerialDenseVector>* dHda = NULL, const double* f_bar_factor = NULL,
          const LINALG::Matrix<numdofperelement_, 1>* htensor = NULL);

      void RecoverPlasticityAndEAS(const LINALG::Matrix<numdofperelement_, 1>* res_d,
          const LINALG::Matrix<nen_, 1>* res_t = NULL);

      void RecoverEAS(const LINALG::Matrix<numdofperelement_, 1>* res_d,
          const LINALG::Matrix<nen_, 1>* temp = NULL);

      template <int spintype>
      void RecoverPlasticity(const LINALG::Matrix<numdofperelement_, 1>* res_d, const int gp,
          const double* temp = NULL);

      void ReduceEasStep(const double new_step_length, const double old_step_length);
      void ReducePlasticityStep(
          const double new_step_length, const double old_step_length, const int gp);

      virtual void BuildDeltaLp(const int gp);

      /*! \brief Return if plastic spin is solved for.
       */
      virtual bool HavePlasticSpin();


      //! calculate internal elastic energy
      double CalcIntEnergy(
          std::vector<double>& disp, std::vector<double>& temp, Teuchos::ParameterList& params);

      //@}

      /**
       * \brief Evaluate Cauchy stress at given point in parameter space and calculate
       * linearizations
       *
       * \param xi[in]          position in parameter space xi
       * \param disp[in]        vector of displacements
       * \param n[in]           vector n (\f[\bold{n}\f])
       * \param t[in]           vector t (\f[\bold{t}\f])
       * \param sigma_nt[out]   cauchy stress tensor contracted using the vectors n and t
                                (\f[ \bold{\sigma} \cdot \bold{n} \cdot \bold{t} \f])
       * \param DsntDd[out]     derivative of sigma_nt w.r.t. displacements (\f[ \frac{ \mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{d}}
                                \f])
       * \param D2sntDd2[out]   second derivative of sigma_nt w.r.t. displacements (\f[ \frac{
                                \mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}
                                {\mathrm{d} \bold{d}^2} \f])
       * \param D2sntDdDn[out]  second derivative of sigma_nt w.r.t. displacements and vector n (\f[
                                \frac{\mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}
                                {\mathrm{d} \bold{d} \mathrm{d} \bold{n} } \f])
       * \param D2sntDdDt[out]  second derivative of sigma_nt w.r.t. displacements and vector t (\f[
                                \frac{\mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}
                                {\mathrm{d} \bold{d} \mathrm{d} \bold{t} } \f])
       * \param D2sntDdDxi[out] second derivative of sigma_nt w.r.t. displacements and local
                                parameter coordinate xi (\f[ \frac{\mathrm{d}^2 \bold{\sigma} \cdot
                                \bold{n} \cdot \bold{t}} {\mathrm{d} \bold{d} \mathrm{d} \bold{\xi}
                                } \f])
       * \param DsntDn[out]     derivative of sigma_nt w.r.t. vector n (\f[ \frac{ \mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{n}}
                                \f])
       * \param DsntDt[out]     derivative of sigma_nt w.r.t. vector n (\f[ \frac{ \mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{t}}
                                \f])
       * \param DsntDxi[out]    derivative of sigma_nt w.r.t. local parameter coordinate xi (\f[
                                \frac{ \mathrm{d} \bold{\sigma} \cdot \bold{n} \cdot
                                \bold{t}}{\mathrm{d} \bold{\xi}} \f])
       * \param temp[in]        temperature
       * \param DsntDT[out]     derivative of sigma_nt w.r.t. temperature (\f[ \frac{ \mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} T} \f])
       * \param D2sntDdDT[out]  second derivative of sigma_nt w.r.t. displacements and temperature
                                (\f[ \frac{\mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}
                                {\mathrm{d} \bold{d} \mathrm{d} T } \f])
       *
       * \note At the moment this method is only used for the nitsche contact formulation
       */
      virtual void GetCauchyAtXiPlast(const LINALG::Matrix<3, 1>& xi,
          const std::vector<double>& disp, const LINALG::Matrix<3, 1>& n,
          const LINALG::Matrix<3, 1>& t, double& sigma_nt, Epetra_SerialDenseMatrix* DsntDd,
          Epetra_SerialDenseMatrix* D2sntDd2, Epetra_SerialDenseMatrix* D2sntDdDn,
          Epetra_SerialDenseMatrix* D2sntDdDt, Epetra_SerialDenseMatrix* D2sntDdDxi,
          LINALG::Matrix<3, 1>* DsntDn, LINALG::Matrix<3, 1>* DsntDt, LINALG::Matrix<3, 1>* DsntDxi,
          const std::vector<double>* temp = NULL, Epetra_SerialDenseMatrix* DsntDT = NULL,
          Epetra_SerialDenseMatrix* D2sntDdDT = NULL);

      /**
       * \brief Evaluate Cauchy stress at given point in parameter space and calculate
       * linearizations
       *
       * \param xi[in]          position in parameter space xi
       * \param disp[in]        vector of displacements
       * \param n[in]           vector n (\f[\bold{n}\f])
       * \param t[in]           vector t (\f[\bold{t}\f])
       * \param sigma_nt[out]   cauchy stress tensor contracted using the vectors n and t
                                (\f[ \bold{\sigma} \cdot \bold{n} \cdot \bold{t} \f])
       * \param DsntDd[out]     derivative of sigma_nt w.r.t. displacements (\f[ \frac{ \mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{d}}
                                \f])
       * \param D2sntDd2[out]   second derivative of sigma_nt w.r.t. displacements (\f[ \frac{
                                \mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}
                                {\mathrm{d} \bold{d}^2} \f])
       * \param D2sntDdDn[out]  second derivative of sigma_nt w.r.t. displacements and vector n (\f[
                                \frac{\mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}
                                {\mathrm{d} \bold{d} \mathrm{d} \bold{n} } \f])
       * \param D2sntDdDt[out]  second derivative of sigma_nt w.r.t. displacements and vector t (\f[
                                \frac{\mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}
                                {\mathrm{d} \bold{d} \mathrm{d} \bold{t} } \f])
       * \param D2sntDdDxi[out] second derivative of sigma_nt w.r.t. displacements and local
                                parameter coordinate xi (\f[ \frac{\mathrm{d}^2 \bold{\sigma} \cdot
                                \bold{n} \cdot \bold{t}} {\mathrm{d} \bold{d} \mathrm{d} \bold{\xi}
                                } \f])
       * \param DsntDn[out]     derivative of sigma_nt w.r.t. vector n (\f[ \frac{ \mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{n}}
                                \f])
       * \param DsntDt[out]     derivative of sigma_nt w.r.t. vector n (\f[ \frac{ \mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{t}}
                                \f])
       * \param DsntDxi[out]    derivative of sigma_nt w.r.t. local parameter coordinate xi (\f[
                                \frac{ \mathrm{d} \bold{\sigma} \cdot \bold{n} \cdot
                                \bold{t}}{\mathrm{d} \bold{\xi}} \f])
       * \param temp[in]        temperature
       * \param DsntDT[out]     derivative of sigma_nt w.r.t. temperature (\f[ \frac{ \mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} T} \f])
       * \param D2sntDdDT[out]  second derivative of sigma_nt w.r.t. displacements and temperature
                                (\f[ \frac{\mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}
                                {\mathrm{d} \bold{d} \mathrm{d} T } \f])
       *
       * \note At the moment this method is only used for the nitsche contact formulation
       */
      virtual void GetCauchyAtXiElast(const LINALG::Matrix<3, 1>& xi,
          const std::vector<double>& disp, const LINALG::Matrix<3, 1>& n,
          const LINALG::Matrix<3, 1>& t, double& sigma_nt, Epetra_SerialDenseMatrix* DsntDd,
          Epetra_SerialDenseMatrix* D2sntDd2, Epetra_SerialDenseMatrix* D2sntDdDn,
          Epetra_SerialDenseMatrix* D2sntDdDt, Epetra_SerialDenseMatrix* D2sntDdDxi,
          LINALG::Matrix<3, 1>* DsntDn, LINALG::Matrix<3, 1>* DsntDt, LINALG::Matrix<3, 1>* DsntDxi,
          const std::vector<double>* temp = NULL, Epetra_SerialDenseMatrix* DsntDT = NULL,
          Epetra_SerialDenseMatrix* D2sntDdDT = NULL);

      virtual void OutputStrains(const int gp,
          const INPAR::STR::StrainType iostrain,           // strain output option
          LINALG::Matrix<numgpt_post, numstr_>* elestrain  // strains at GP
      );

      virtual void OutputStress(const int gp,
          const INPAR::STR::StressType iostress,           // strain output option
          LINALG::Matrix<numgpt_post, numstr_>* elestress  // strains at GP
      );

      virtual void Kinematics(const int gp = -1);

      virtual void IntegrateMassMatrix(
          const int gp, LINALG::Matrix<numdofperelement_, numdofperelement_>& mass);

      virtual void IntegrateStiffMatrix(const int gp,
          LINALG::Matrix<numdofperelement_, numdofperelement_>& stiff,
          Epetra_SerialDenseMatrix& Kda);

      virtual void IntegrateForce(const int gp, LINALG::Matrix<numdofperelement_, 1>& force);

      virtual void IntegrateThermoGp(const int gp, Epetra_SerialDenseVector& dHda);


      // algoirthmic parameters
      bool fbar_;

      // plasticity ********************************************
      // Kbb^-1 at each Gauss point for recovery of inner variables
      std::vector<LINALG::SerialDenseMatrix> KbbInv_;

      // Kbd at each Gauss point for recovery of inner variables
      std::vector<LINALG::SerialDenseMatrix> Kbd_;

      // f_b at each Gauss point for recovery of inner variables
      std::vector<LINALG::SerialDenseVector> fbeta_;

      // plastic flow at each Gauss point at last Newton iteration
      std::vector<LINALG::SerialDenseVector> dDp_last_iter_;

      // increment of plastic flow over last Newton step (needed for Newton line search)
      std::vector<LINALG::SerialDenseVector> dDp_inc_;

      PlSpinType plspintype_;

      // line search parameter (old step length)
      double old_step_length_;

      // EAS element technology ************************************
      Teuchos::RCP<LINALG::SerialDenseMatrix> KaaInv_;
      Teuchos::RCP<LINALG::SerialDenseMatrix> Kad_;
      Teuchos::RCP<LINALG::SerialDenseMatrix> KaT_;
      Teuchos::RCP<LINALG::Matrix<numdofperelement_, nen_>> KdT_eas_;
      Teuchos::RCP<LINALG::SerialDenseVector> feas_;
      Teuchos::RCP<std::vector<LINALG::SerialDenseMatrix>> Kba_;
      Teuchos::RCP<LINALG::SerialDenseVector> alpha_eas_;
      Teuchos::RCP<LINALG::SerialDenseVector> alpha_eas_last_timestep_;
      Teuchos::RCP<LINALG::SerialDenseVector> alpha_eas_delta_over_last_timestep_;
      Teuchos::RCP<LINALG::SerialDenseVector> alpha_eas_inc_;
      DRT::ELEMENTS::So3Plast_EASType eastype_;
      int neas_;

      // TSI ******************************************************
      bool tsi_;

      /// derivative of the internal force vector w.r.t. temperature at this GP
      /// derivative w.r.t. gp temperature is sufficient as the gp temperature depends linearily
      /// on the nodal values.
      Teuchos::RCP<std::vector<LINALG::Matrix<numdofperelement_, 1>>> dFintdT_;

      /// derivative of NCP w.r.t. temperatures at gp
      Teuchos::RCP<std::vector<LINALG::SerialDenseVector>> KbT_;

      /// Temperature at each gp in last Newton iteration
      /// this is needed for the recovery of the plastic flow using KbT_
      /// As the gp temperature depends linearly on the nodal temperature dofs
      /// a scalar value is sufficient
      Teuchos::RCP<std::vector<double>> temp_last_;

      // Cauchy stress for Nitsche contact **************************
      bool is_nitsche_contact_;
      std::vector<LINALG::Matrix<numstr_, 1>> cauchy_;
      std::vector<LINALG::Matrix<numstr_, numdofperelement_>> cauchy_deriv_;
      std::vector<LINALG::Matrix<numstr_, nen_>> cauchy_deriv_T_;

      // Static members for faster evaluation **************************
      static std::pair<bool, LINALG::Matrix<nen_, 1>> shapefunct_;
      static std::pair<bool, LINALG::Matrix<nsd_, nen_>> deriv_;
      static std::pair<bool, LINALG::Matrix<nsd_, nsd_>> invJ_;
      static std::pair<bool, double> detJ_;
      static std::pair<bool, LINALG::Matrix<nsd_, nen_>> N_XYZ_;  // N_XYZ = J^-1 * N_rst
      static std::pair<bool, LINALG::Matrix<nsd_, nsd_>>
          defgrd_;  //     deformation gradient consistent with displacements
      static std::pair<bool, LINALG::Matrix<nsd_, nsd_>>
          defgrd_mod_;  // deformation gradient consistent with displacements + element technology
      static std::pair<bool, LINALG::Matrix<nsd_, nsd_>> rcg_;
      static std::pair<bool, LINALG::Matrix<nsd_, nsd_>>
          delta_Lp_;  // plastic velocity increment over this time step
      static std::pair<bool, LINALG::Matrix<numstr_, numdofperelement_>> bop_;
      static std::pair<bool, LINALG::Matrix<numstr_, 1>> pk2_;
      static std::pair<bool, LINALG::Matrix<numstr_, numstr_>> cmat_;

      static std::pair<bool, LINALG::Matrix<nen_, nsd_>> xrefe_;  // X, material coord. of element
      static std::pair<bool, LINALG::Matrix<nen_, nsd_>> xcurr_;  // x, current  coord. of element
      static std::pair<bool, LINALG::Matrix<nen_, nsd_>>
          xcurr_rate_;  // xdot, rate of current  coord. of element
      static std::pair<bool, LINALG::Matrix<nen_, 1>>
          etemp_;  // vector of the current element temperatures

      // Static members for Fbar related stuff **************************
      static std::pair<bool, double> detF_;
      static std::pair<bool, double> detF_0_;
      static std::pair<bool, LINALG::Matrix<nsd_, nsd_>> inv_defgrd_;
      static std::pair<bool, LINALG::Matrix<nsd_, nsd_>> inv_defgrd_0_;
      static std::pair<bool, LINALG::Matrix<nsd_, nen_>> N_XYZ_0_;
      static std::pair<bool, LINALG::Matrix<numstr_, 1>> rcg_vec_;  // strain-like
      static std::pair<bool, double> f_bar_fac_;
      static std::pair<bool, LINALG::Matrix<numdofperelement_, 1>> htensor_;

      // Static members for EAS related stuff **************************
      // transformation matrix T0, maps M-matrix evaluated at origin
      // between local element coords and global coords
      // here we already get the inverse transposed T0
      static std::pair<bool, LINALG::Matrix<numstr_, numstr_>> T0invT_;
      static std::pair<bool, LINALG::Matrix<nsd_, nsd_>> jac_0_;
      static std::pair<bool, double> det_jac_0_;
      static std::pair<bool, LINALG::SerialDenseMatrix> M_eas_;  // EAS matrix M at current GP

      // NURBS-specific
      static std::pair<bool, LINALG::Matrix<nen_, 1>> weights_;
      static std::pair<bool, std::vector<Epetra_SerialDenseVector>> knots_;

      virtual void InvalidGpData()
      {
        shapefunct_.first = false;
        deriv_.first = false;
        invJ_.first = false;
        detJ_.first = false;
        N_XYZ_.first = false;
        defgrd_.first = false;
        defgrd_mod_.first = false;
        rcg_.first = false;
        delta_Lp_.first = false;
        bop_.first = false;
        detF_.first = false;
        f_bar_fac_.first = false;
        htensor_.first = false;
        inv_defgrd_.first = false;
        rcg_vec_.first = false;
        M_eas_.first = false;
        pk2_.first = false;
        cmat_.first = false;
      }

      void InvalidEleData()
      {
        xrefe_.first = false;
        xcurr_.first = false;
        xcurr_rate_.first = false;
        etemp_.first = false;
        detF_0_.first = false;
        inv_defgrd_0_.first = false;
        N_XYZ_0_.first = false;
        T0invT_.first = false;
        jac_0_.first = false;
        det_jac_0_.first = false;
        weights_.first = false;
        knots_.first = false;
      }

      inline const LINALG::Matrix<nen_, 1>& Weights() const
      {
        dsassert(weights_.first, "weights_ not valid");
        return weights_.second;
      }
      inline LINALG::Matrix<nen_, 1>& SetWeights()
      {
        weights_.first = true;
        return weights_.second;
      }

      inline const std::vector<Epetra_SerialDenseVector>& Knots() const
      {
        dsassert(knots_.first, "weights_ not valid");
        return knots_.second;
      }
      inline std::vector<Epetra_SerialDenseVector>& SetKnots()
      {
        knots_.first = true;
        return knots_.second;
      }

      void FillPositionArrays(const std::vector<double>& disp,  // current displacements
          const std::vector<double>& vel,                       // current velocities
          const std::vector<double>& temp                       // current temperatures
      )
      {
        for (int i = 0; i < nen_; ++i)
        {
          for (int d = 0; d < nsd_; ++d)
          {
            xrefe_.second(i, d) = Nodes()[i]->X()[d];
            xcurr_.second(i, d) = Nodes()[i]->X()[d] + disp[i * numdofpernode_ + d];
            if (vel.size()) xcurr_rate_.second(i, d) = vel[i * numdofpernode_ + d];
          }
          if (temp.size()) etemp_.second(i) = temp[i];
        }
        xrefe_.first = true;
        xcurr_.first = true;
        if (vel.size()) xcurr_rate_.first = true;
        if (temp.size()) etemp_.first = true;
      }

      inline const LINALG::Matrix<nen_, nsd_>& Xrefe()
      {
        dsassert(xrefe_.first == true, "xrefe not valid");
        return xrefe_.second;
      }

      inline const LINALG::Matrix<nen_, nsd_>& Xcurr()
      {
        dsassert(xcurr_.first, "xcurr_ not valid");
        return xcurr_.second;
      }

      inline const LINALG::Matrix<nen_, nsd_>& XcurrRate()
      {
        dsassert(xcurr_rate_.first, "xcurr_rate_ not valid");
        return xcurr_rate_.second;
      }

      inline const LINALG::Matrix<nen_, 1>& Temp()
      {
        dsassert(etemp_.first, "etemp not valid");
        return etemp_.second;
      }

      inline const LINALG::Matrix<nen_, 1>& ShapeFunction() const
      {
        dsassert(shapefunct_.first, "shape function not valid");
        return shapefunct_.second;
      }
      inline LINALG::Matrix<nen_, 1>& SetShapeFunction()
      {
        shapefunct_.first = true;
        return shapefunct_.second;
      }

      inline const LINALG::Matrix<nsd_, nen_>& DerivShapeFunction() const
      {
        dsassert(deriv_.first, "deriv shape function not valid");
        return deriv_.second;
      }
      inline LINALG::Matrix<nsd_, nen_>& SetDerivShapeFunction()
      {
        deriv_.first = true;
        return deriv_.second;
      }

      inline const LINALG::Matrix<nsd_, nen_>& DerivShapeFunctionXYZ() const
      {
        dsassert(N_XYZ_.first, "deriv shape function not valid");
        return N_XYZ_.second;
      }
      inline LINALG::Matrix<nsd_, nen_>& SetDerivShapeFunctionXYZ()
      {
        N_XYZ_.first = true;
        return N_XYZ_.second;
      }

      inline const LINALG::Matrix<nsd_, nsd_>& InvJ() const
      {
        dsassert(invJ_.first, "invJ_ not valid");
        return invJ_.second;
      }
      inline LINALG::Matrix<nsd_, nsd_>& SetInvJ()
      {
        invJ_.first = true;
        return invJ_.second;
      }

      inline const double& DetJ() const
      {
        dsassert(detJ_.first, "detJ_ not valid");
        return detJ_.second;
      }
      inline double& SetDetJ()
      {
        detJ_.first = true;
        return detJ_.second;
      }

      inline const LINALG::Matrix<nsd_, nsd_>& Defgrd() const
      {
        dsassert(defgrd_.first, "defgrd_ not valid");
        return defgrd_.second;
      }
      inline LINALG::Matrix<nsd_, nsd_>& SetDefgrd()
      {
        defgrd_.first = true;
        return defgrd_.second;
      }

      inline const LINALG::Matrix<nsd_, nsd_>& DefgrdMod() const
      {
        dsassert(defgrd_mod_.first, "defgrd_mod_ not valid");
        return defgrd_mod_.second;
      }
      inline LINALG::Matrix<nsd_, nsd_>& SetDefgrdMod()
      {
        defgrd_mod_.first = true;
        return defgrd_mod_.second;
      }

      inline const LINALG::Matrix<nsd_, nsd_>& RCG() const
      {
        dsassert(rcg_.first, "rcg_ not valid");
        return rcg_.second;
      }
      inline LINALG::Matrix<nsd_, nsd_>& SetRCG()
      {
        rcg_.first = true;
        return rcg_.second;
      }

      inline const LINALG::Matrix<nsd_, nsd_>& DeltaLp() const
      {
        dsassert(delta_Lp_.first, "delta_Lp_ not valid");
        return delta_Lp_.second;
      }
      inline LINALG::Matrix<nsd_, nsd_>& SetDeltaLp()
      {
        delta_Lp_.first = true;
        return delta_Lp_.second;
      }

      inline const LINALG::Matrix<numstr_, numdofperelement_>& Bop() const
      {
        dsassert(bop_.first, "bop_ not valid");
        return bop_.second;
      }
      inline LINALG::Matrix<numstr_, numdofperelement_>& SetBop()
      {
        bop_.first = true;
        return bop_.second;
      }

      inline const LINALG::Matrix<numstr_, 1>& PK2() const
      {
        dsassert(pk2_.first, "pk2_ not valid");
        return pk2_.second;
      }
      inline LINALG::Matrix<numstr_, 1>& SetPK2()
      {
        pk2_.first = true;
        return pk2_.second;
      }

      inline const LINALG::Matrix<numstr_, numstr_>& Cmat() const
      {
        dsassert(cmat_.first, "cmat_ not valid");
        return cmat_.second;
      }
      inline LINALG::Matrix<numstr_, numstr_>& SetCmat()
      {
        cmat_.first = true;
        return cmat_.second;
      }

      inline const LINALG::Matrix<nsd_, nen_>& DerivShapeFunctionXYZ_0() const
      {
        dsassert(N_XYZ_0_.first, "deriv shape function not valid");
        return N_XYZ_0_.second;
      }
      inline LINALG::Matrix<nsd_, nen_>& SetDerivShapeFunctionXYZ_0()
      {
        N_XYZ_0_.first = true;
        return N_XYZ_0_.second;
      }

      inline const double& DetF() const
      {
        dsassert(detF_.first, "detF_ not valid");
        return detF_.second;
      }
      inline double& SetDetF()
      {
        detF_.first = true;
        return detF_.second;
      }

      inline const double& DetF_0() const
      {
        dsassert(detF_0_.first, "detF_0_ not valid");
        return detF_0_.second;
      }
      inline double& SetDetF_0()
      {
        detF_0_.first = true;
        return detF_0_.second;
      }

      inline const LINALG::Matrix<nsd_, nsd_>& InvDefgrd() const
      {
        dsassert(inv_defgrd_.first, "inv_defgrd_ not valid");
        return inv_defgrd_.second;
      }
      inline LINALG::Matrix<nsd_, nsd_>& SetInvDefgrd()
      {
        inv_defgrd_.first = true;
        return inv_defgrd_.second;
      }

      inline const LINALG::Matrix<nsd_, nsd_>& InvDefgrd_0() const
      {
        dsassert(inv_defgrd_0_.first, "inv_defgrd_0_ not valid");
        return inv_defgrd_0_.second;
      }
      inline LINALG::Matrix<nsd_, nsd_>& SetInvDefgrd_0()
      {
        inv_defgrd_0_.first = true;
        return inv_defgrd_0_.second;
      }

      inline const LINALG::Matrix<nsd_, nsd_>& Jac_0() const
      {
        dsassert(jac_0_.first, "jac_0_ not valid");
        return jac_0_.second;
      }
      inline LINALG::Matrix<nsd_, nsd_>& SetJac_0()
      {
        jac_0_.first = true;
        return jac_0_.second;
      }

      inline const double& DetJac_0() const
      {
        dsassert(det_jac_0_.first, "det_jac_0_ not valid");
        return det_jac_0_.second;
      }
      inline double& SetDetJac_0()
      {
        det_jac_0_.first = true;
        return det_jac_0_.second;
      }

      inline const LINALG::Matrix<numstr_, 1>& RCGvec() const
      {
        dsassert(rcg_vec_.first, "rcg_vec_ not valid");
        return rcg_vec_.second;
      }
      inline LINALG::Matrix<numstr_, 1>& SetRCGvec()
      {
        rcg_vec_.first = true;
        return rcg_vec_.second;
      }

      inline const double& FbarFac() const
      {
        dsassert(f_bar_fac_.first, "f_bar_fac_ not valid");
        return f_bar_fac_.second;
      }
      inline double& SetFbarFac()
      {
        f_bar_fac_.first = true;
        return f_bar_fac_.second;
      }

      inline const LINALG::Matrix<numdofperelement_, 1>& Htensor() const
      {
        dsassert(htensor_.first, "htensor_ not valid");
        return htensor_.second;
      }
      inline LINALG::Matrix<numdofperelement_, 1>& SetHtensor()
      {
        htensor_.first = true;
        return htensor_.second;
      }

      void EvaluateCenter()
      {
        // element coordinate derivatives at centroid
        static LINALG::Matrix<nsd_, nen_> N_rst_0(false);
        DRT::UTILS::shape_function_3D_deriv1(N_rst_0, 0.0, 0.0, 0.0, DRT::Element::hex8);

        // inverse jacobian matrix at centroid
        SetJac_0().Multiply(N_rst_0, Xrefe());
        static LINALG::Matrix<nsd_, nsd_> invJ_0;
        SetDetJac_0() = invJ_0.Invert(Jac_0());
        // material derivatives at centroid
        SetDerivShapeFunctionXYZ_0().Multiply(invJ_0, N_rst_0);

        // deformation gradient and its determinant at centroid
        static LINALG::Matrix<3, 3> defgrd_0(false);
        defgrd_0.MultiplyTT(Xcurr(), DerivShapeFunctionXYZ_0());
        SetDetF_0() = SetInvDefgrd_0().Invert(defgrd_0);
      }

      void SetupFbarGp()
      {
        if (DetF() < 0. || DetF_0() < 0.) dserror("element distortion too large");
        SetFbarFac() = pow(DetF_0() / DetF(), 1. / 3.);
        SetDefgrdMod().Update(SetFbarFac(), Defgrd());
        SetHtensor().Clear();

        for (int n = 0; n < numdofperelement_; n++)
          for (int i = 0; i < 3; i++)
            SetHtensor()(n) += InvDefgrd_0()(i, n % 3) * DerivShapeFunctionXYZ_0()(i, n / 3) -
                               InvDefgrd()(i, n % 3) * DerivShapeFunctionXYZ()(i, n / 3);
      }

      inline const LINALG::Matrix<numstr_, numstr_>& T0invT() const
      {
        dsassert(T0invT_.first, "T0invT_ not valid");
        return T0invT_.second;
      }
      inline LINALG::Matrix<numstr_, numstr_>& SetT0invT()
      {
        T0invT_.first = true;
        return T0invT_.second;
      }

      inline const LINALG::SerialDenseMatrix& M_eas() const
      {
        dsassert(M_eas_.first, "M_eas_ not valid");
        return M_eas_.second;
      }
      inline LINALG::SerialDenseMatrix& SetM_eas()
      {
        M_eas_.first = true;
        return M_eas_.second;
      }

      inline void EvaluateShape(const LINALG::Matrix<3, 1>& xi)
      {
        if (distype == DRT::Element::nurbs27)
          DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv(
              SetShapeFunction(), SetDerivShapeFunction(), xi, Knots(), Weights(), distype);
        else
          DRT::UTILS::shape_function<distype>(xi, SetShapeFunction());
      }

      inline void EvaluateShapeDeriv(const LINALG::Matrix<3, 1>& xi)
      {
        if (distype == DRT::Element::nurbs27)
          DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv(
              SetShapeFunction(), SetDerivShapeFunction(), xi, Knots(), Weights(), distype);
        else
          DRT::UTILS::shape_function_deriv1<distype>(xi, SetDerivShapeFunction());
      }

      void GetNurbsEleInfo(DRT::Discretization* dis = NULL);


    };  // class So3_Plast

    template <DRT::ELEMENTS::So3Plast_EASType eastype>
    struct PlastEasTypeToNumEas
    {
    };
    template <>
    struct PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easmild>
    {
      static const int neas = 9;
    };
    template <>
    struct PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easfull>
    {
      static const int neas = 21;
    };
    template <>
    struct PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_eassosh8>
    {
      static const int neas = 7;
    };
    template <>
    struct PlastEasTypeToNumEas<DRT::ELEMENTS::soh18p_eassosh18>
    {
      static const int neas = 9;
    };
    template <>
    struct PlastEasTypeToNumEas<DRT::ELEMENTS::soh8p_easnone>
    {
      static const int neas = 0;
    };

    int PlastEasTypeToNumEasV(DRT::ELEMENTS::So3Plast_EASType et);

  }  // namespace ELEMENTS
}  // namespace DRT

/*----------------------------------------------------------------------*/
#endif  // SO3_SSN_PLAST_H
