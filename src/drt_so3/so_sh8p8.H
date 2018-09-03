/*----------------------------------------------------------------------*/
/*!
\file so_sh8p8.H
\level 2
\brief
\maintainer Alexander Seitz
*/

/*----------------------------------------------------------------------*/

#ifndef SO_SH8P8_H
#define SO_SH8P8_H

/*----------------------------------------------------------------------*/
#include "so_sh8.H"


// forward declarations
struct _SOH8_DATA;

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    // forward declarations
    class So_sh8;

    class So_sh8p8Type : public So_sh8Type
    {
     public:
      std::string Name() const { return "So_sh8p8Type"; }

      static So_sh8p8Type& Instance();

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
      static So_sh8p8Type instance_;
    };

    /// An incompressible 8-node solid shell element inherited from #DRT::ELEMENTS::so_sh8
    /// utilising a Bochev-stabilised equal-order approach for
    /// tri-linearly Lagrangean interpolated displacement and pressure fields
    ///
    /// <h3>References</h3>
    /// The equal-order Bochev-stabilised approach for incompressible solid
    /// is based on
    /// - [1] C.R. Dohrmann, P.B. Bochev, "A stabilized finite element for the Stokes problem based
    /// on plynomial projections", Int. J. Numer. Fluids, 2000.
    /// - [2] P.B. Bochev, C.R. Dohrmann, "Stabilization of low-order mixed finite elements for the
    /// Stokes equations".
    /// - [3] C. Foerster, "Zur Bochev-Stabilisierung", internal note, 16.7.2008.
    /// - [4] B. Bornemann, "Zum Bochev-stabilisierten Strukturelement fuer
    ///     geometrisch nicht-lineare Verformungen", internal note, 5.11.2008.
    ///
    /// The solid shell/ANS  element technology is based on
    /// - [5] M.A. Frenzel, "An advanced ...", PhD thesis, LNM, TU Muenchen, 2009.
    /// - [6] T. Vu-Quoc, "Optimal solid shells for non-linear analyses
    ///     of multilayer composites", CMAME 2003
    /// - [7] S. Klinkel, Gruttmann, W. Wagner, "A robust non-linear solid shell element
    ///     based on a mixed variational fromulation", Comp. Meth. in Appl. Mech. and Engrg.,
    ///     195:1-3, p. 179-201, 2006.
    ///
    /// \author bborn
    /// \date 03/09
    class So_sh8p8 : public So_sh8
    {
     public:
      /// @name Friends
      //@{
      friend class So_sh8p8Type;
      friend class Soh8Surface;
      friend class Soh8Line;
      //@}

     public:
      /// @name Properties
      //@{

      /// Kind of stabilisation for mixed, equal-order displacement-pressure appraoch
      enum StabilisationType
      {
        stab_affine,     ///< Bochev stabilisation for affine elements in material element domain
                         ///< (default)
        stab_nonaffine,  ///< Bochev stabilisation for non-affine elements in material element
                         ///< domain
        stab_spatialaffine,  ///< Bochev stabilisation for affine elements in spatial element domain
        stab_spatial,  ///< Bochev stabilisation for non-affine elements in spatial element domain
        stab_puredisp  ///< DEBUG ONLY: recover pure displacement-based approach
      };

      /// Kind of ANS anti-locking technique
      enum AnsType
      {
        ans_none = 0,  ///< DEBUG ONLY: ANS switched off
        ans_lateral,   ///< ANS active in t-/out-of-plane/thickness direction, ordinary
                       ///< solid-shell-like
        ans_onspot  ///< ANS active in t-/out-of-plane/thickness direction, specially adjusted at 4
                    ///< spots
      };

      /// way to obtain isochoric material stress response
      enum IsochoricType
      {
        iso_material,  ///< done completely by material model, i.e. material model _is_ isochoric
        iso_enforced   ///< volumetric contribution is split of material stress reponse
      };

      /// Type of Voigt vector
      enum VoigtType
      {
        voigt6_strain,  ///< strain-like vector [ E11, E22, E33, 2*E12, 2*E23, 2*E31 ]
        voigt6_stress,  ///< stress-like vector [ S11, S22, S33, S12, S23, S31 ]
        voigt9          ///< non-symmetric vector [ F11, F22, F33, F12, F23, F31, F13, F32, F21 ]
      };

      // LinearizationType
      enum LinearizationType
      {
        lin_sixth = 0,  ///< minimum linearisation to achieve anti-locking
        lin_half,       ///< half-way linearisation
        lin_one         ///< full linearisation
      };

      //@}

      /// @name Element constants
      //@{

      static const int NUMNOD_ = 8;          ///< number of nodes
      static const int NODDOF_ = 4;          ///< number of DOFs per node
      static const int NODDISP_ = 3;         ///< number of displacements per node
      static const int NODPRES_ = 1;         ///< number of pressures per node
      static const int NUMDOF_ = 32;         ///< total DOFs per element
      static const int NUMDISP_ = 24;        ///< total discrete displacements per element
      static const int NUMDISPSQSYM_ = 300;  ///< (NUMDISP_*NUMDISP_ + NUMDISP_)/2, number
                                             ///< of relevant entries
                                             ///< if symmetric in element displacements
      static const int NUMPRES_ = 8;         ///< total discrete pressures per element
      static const int NUMPRESBRO_ = 1;      ///< total number of discrete, discontinuous/broken
                                             ///< pressures per element
      static const int NUMDFGR_ = 9;         ///< number of deformation gradient components
                                             ///< (deformation gradient is non-symmetric)
      static const int NUMGPT_ = 8;          ///< total gauss points per element
      static const int NUMDIM_ = 3;          ///< number of dimensions, it's 3D

      /// number of ANS sampling points, here 8
      static const int NUMSP_ = 8;
      static const int NUMSP2ND_ = 8;  ///< number of 2nd set of sampling points
      static const int NUMSP3RD_ = 8;  ///< number of 3rd set of sampling points
      /// number of modified ANS strains (E_rt,E_st,E_tt), here 3
      static const int NUMANS_ = 3;

      // number of EAS parameters
      static const int NUMEAS_SOSH8_ = 7;  ///< number of EAS parameters of sosh8-like EAS
      static const int NUMEAS_A_ = 1;      ///< number of EAS parameters of A-type EAS

      // Voigt index arrays
      static const int VOIGT6ROW_[];    ///< 6-Voigt row index of corresponding 2-tensor
      static const int VOIGT6COL_[];    ///< 6-Voigt column index of corresponding 2-tensor
      static const int VOIGT9ROW_[];    ///< 9-Voigt row index of corresponding 2-tensor
      static const int VOIGT9COL_[];    ///< 9-Voigt column index of corresponding 2-tensor
      static const int VOIGT3X3_[];     ///< go from 2-tensor index pair to 9-Voigt index
                                        ///< by [NUMDIM_*i+j] for any i,j=0,1,2
      static const int VOIGT3X3SYM_[];  ///< go from symmetric 2-tensor index pair to 6-Voigt index
                                        ///< by [NUMDIM_*i+j] for any i,j=0,1,2

      // Ordering of element displacement and pressure DOFs
      static const int DISPTODISPPRES_[];  ///< 24 displacement into 32 total DOFs
      static const int PRESTODISPPRES_[];  ///< 8 pressures into 32 total DOFs

      //@}

     public:
      /// @name Constructors and destructors and related methods
      //@{

      /// Standard Constructor
      So_sh8p8(int id,  ///<  A unique global ID
          int owner     ///< elements owning processor
      );

      /// Copy Constructor
      ///
      /// Makes a deep copy of a Element
      So_sh8p8(const So_sh8p8& old);

      /// Deep copy this instance of Solid3 and return pointer to the copy
      ///
      /// The Clone() method is used from the virtual base class Element in cases
      /// where the type of the derived class is unknown and a copy-ctor is needed
      DRT::Element* Clone() const;

      /// Return unique ParObject id
      ///
      /// every class implementing ParObject needs a unique id defined at the
      /// top of this file.
      virtual int UniqueParObjectId() const { return So_sh8p8Type::Instance().UniqueParObjectId(); }

      /// Pack this class so it can be communicated
      ///
      /// \ref Pack and \ref Unpack are used to communicate this element
      virtual void Pack(DRT::PackBuffer& data) const;

      /// Unpack data from a char vector into this class
      ///
      /// \ref Pack and \ref Unpack are used to communicate this element
      virtual void Unpack(const std::vector<char>& data);

      /// Print this element
      void Print(std::ostream& os) const;

      virtual So_sh8p8Type& ElementType() const { return So_sh8p8Type::Instance(); }

      //@}

      /// @name Input and Creation
      //@{

      /// Read input for this element
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef);


      /// Destructor
      virtual ~So_sh8p8();

      //@}

      /// @name Access methods
      //@{

      /// Get number of degrees of freedom of a certain node
      /// (implements pure virtual DRT::Element)
      ///
      /// The element decides how many degrees of freedom its nodes must have.
      /// As this may vary along a simulation, the element can redecide the
      /// number of degrees of freedom per node along the way for each of it's nodes
      /// separately.
      virtual int NumDofPerNode(const DRT::Node& node) const { return 4; }

      /// Get number of degrees of freedom per element
      /// (implements pure virtual DRT::Element)
      ///
      /// The element decides how many element degrees of freedom it has.
      /// It can redecide along the way of a simulation.
      ///
      /// \note Element degrees of freedom mentioned here are dofs that are visible
      /// at the level of the total system of equations. Purely internal
      /// element dofs that are condensed internally should NOT be considered.
      virtual int NumDofPerElement() const { return 0; }

      /// Set ANS type
      void SetANS(const AnsType& newans  ///< ANS to set
      )
      {
        ans_ = newans;
      }

      //@}

      /// @name Evaluation
      //@{

      /// Evaluate an element
      ///
      /// Evaluate so_sh8p8 element stiffness, mass, internal forces, etc.
      ///
      /// \return 0 if successful, negative otherwise
      int Evaluate(Teuchos::ParameterList& params,  ///< (in/out) ParameterList for communication
                                                    ///< between control routine and elements
          DRT::Discretization& discretization,  ///< pointer to discretization for de-assembly
          std::vector<int>& lm,                 ///< (in) location matrix for de-assembly
          Epetra_SerialDenseMatrix& elemat1,  ///< (out) (stiffness-)matrix to be filled by element.
                                              ///< If NULL on input, the controling method does not
                                              ///< expect the element to fill this matrix.
          Epetra_SerialDenseMatrix& elemat2,  ///< (out) (mass-)matrix to be filled by element. If
                                              ///< NULL on input, the controling method does not
                                              ///< expect the element to fill this matrix.
          Epetra_SerialDenseVector& elevec1,  ///< (out) (internal force-)vector to be filled by
                                              ///< element. If NULL on input, the controlling method
                                              ///< does not expect the element to fill this vector
          Epetra_SerialDenseVector& elevec2,  ///< (out)  vector to be filled by element. If NULL on
                                              ///< input, the controlling method does not expect the
                                              ///< element to fill this vector
          Epetra_SerialDenseVector& elevec3  ///< (out) vector to be filled by element. If NULL on
                                             ///< input, the controlling method does not expect the
                                             ///< element to fill this vector
      );


      /// \brief Evaluate a Neumann boundary condition
      ///
      /// this method evaluates a surface Neumann condition on the solid3 element
      ///
      /// \param params (in/out)    : ParameterList for communication between control routine
      ///                             and elements
      /// \param discretization (in): A reference to the underlying discretization
      /// \param condition (in)     : The condition to be evaluated
      /// \param lm (in)            : location vector of this element
      /// \param elevec1 (out)      : vector to be filled by element. If NULL on input,
      ///
      /// \return 0 if successful, negative otherwise
      virtual int EvaluateNeumann(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1 = NULL);

     private:
      /// don't want = operator
      So_sh8p8& operator=(const So_sh8p8& old);

      /// kind of stabilisation, cf. #StabilisationType
      StabilisationType stab_;

      /// kind of ANS, cf. #AnsType
      AnsType ans_;

      /// kind of linearization
      LinearizationType lin_;

      /// kind of isotropic material response handling
      IsochoricType iso_;

      /// @name Evaluate force and stiffness
      //@{

      /// Compute stiffness and mass matrix
      void ForceStiffMass(const std::vector<int>& lm,       ///< location matrix
          const LINALG::Matrix<NUMDISP_, 1>& disp,          ///< current displacements
          const LINALG::Matrix<NUMPRES_, 1>& pres,          ///< current pressures
          const LINALG::Matrix<NUMDISP_, 1>& dispi,         ///< last residual displacements
          const LINALG::Matrix<NUMPRES_, 1>& presi,         //< last residual pressures
          LINALG::Matrix<NUMDISP_, NUMDISP_>* massmatrix,   ///< element mass matrix
          LINALG::Matrix<NUMDISP_, NUMDISP_>* stiffmatrix,  ///< element stiffness matrix
          LINALG::Matrix<NUMDISP_, NUMPRES_>* gradmatrix,   ///< element gradient matrix
          LINALG::Matrix<NUMPRES_, NUMDISP_>* dargmatrix,  ///< element 'transposed' gradient matrix
          LINALG::Matrix<NUMPRES_, NUMPRES_>* stabmatrix,  ///< element stabilisation matrix
          LINALG::Matrix<NUMDISP_, 1>* force,              ///< element internal force vector
          LINALG::Matrix<NUMPRES_, 1>* incomp,             ///< incompressibility residual
          LINALG::Matrix<NUMGPT_, MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
          LINALG::Matrix<NUMGPT_, MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
          double* volume,                                          ///< current element volume
          Teuchos::ParameterList& params,         ///< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress,  ///< stress output option
          const INPAR::STR::StrainType iostrain   ///< strain output option
      );

      /// Return stress at Gauss point
      void Stress(
          LINALG::Matrix<NUMGPT_, MAT::NUM_STRESS_3D>* elestress,  ///< store the stress herein
          const INPAR::STR::StressType iostress,                   ///< stress type
          const int gp,                                            ///< Gauss point index
          const double& detdefgrd,  ///< determinant of (assumed) deformation gradient
          const LINALG::Matrix<NUMDIM_, NUMDIM_>& defgrd,  ///< (assumed) deformation gradient
          const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& glstrain,  ///< Green-Lagrange strain vector
          const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>&
              stress,             ///< (deviatoric) 2nd Piola-Kirchhoff stress vector
          const double& pressure  ///< true pressure
      );

      /// Return strain at Gauss point
      void Strain(
          LINALG::Matrix<NUMGPT_, MAT::NUM_STRESS_3D>* elestrain,  ///< store the strain herein
          const INPAR::STR::StrainType iostrain,                   ///< strain type
          const int gp,                                            ///< Gauss point index
          const double& detdefgrd,  ///< determinant of (assumed) deformation gradient
          const LINALG::Matrix<NUMDIM_, NUMDIM_>& defgrd,  ///< (assumed) deformation gradient
          const LINALG::Matrix<NUMDIM_, NUMDIM_>&
              invdefgrd,  ///< (assumed) inverted deformation gradient
          const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& glstrain  ///< Green-Lagrange strain vector
      );

      /// Recover deformation gradient incoperating assumed natural GL strain
      static void AssDefGrad(double& detdefgrad,  ///< determinat of deformation gradient
          LINALG::Matrix<NUMDIM_, NUMDIM_>&
              defgrad,  ///< deformation gradient \f$[\boldsymbol{F}]\f$
          LINALG::Matrix<NUMDIM_, NUMDIM_>&
              invdefgrad,  ///< inverse deformation gradient \f$[\boldsymbol{F}^{-1}]\f$
          LINALG::Matrix<NUMDIM_, NUMDIM_>&
              rgtstr,  ///< right stretch tensor \f$[\boldsymbol{U}]\f$
          LINALG::Matrix<NUMDIM_, NUMDIM_>&
              defgradD,  ///< pure disp-based deformation gradient \f$[\boldsymbol{F}^d]\f$
          LINALG::Matrix<NUMDIM_, NUMDIM_>&
              rgtstrD,  ///< pure disp-based right stretch tensor \f$[\boldsymbol{U}^d]\f$
          LINALG::Matrix<NUMDIM_, NUMDIM_>& invrgtstrD,  ///< inverted pure disp-based right stretch
                                                         ///< tensor \f$[\boldsymbol{U}^{d-1}]\f$
          const LINALG::Matrix<NUMDIM_, NUMDIM_>&
              Jinv,  ///< inverse of transposed material Jacobi matrix \f$[X_{,\xi}]\f$
          const LINALG::Matrix<NUMDIM_, NUMDIM_>&
              Jac,  ///< transposed material Jacobi matrix \f$[X_{,\xi}]^T\f$
          const LINALG::Matrix<NUMDIM_, NUMDIM_>&
              jac,  ///< transposed spatial Jacobi matrix \f$[x_{,\xi}]^T\f$
          const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>&
              glstrain  ///< material Green-Lagrange strain vector
                        ///< in global Cartesian components \f$[\boldsymbol{E}]\f$
      );

      /// Retrieve shear modulus
      ///
      /// Shear modulus is needed for stabilisation
      double ShearMod() const;

      /// Extrapolate Gauss-point values (e.g. stresses) to nodes and store results in elevectors
      void sosh8p8_expol(LINALG::Matrix<NUMGPT_, MAT::NUM_STRESS_3D>& stresses,  ///< gp stresses
          Epetra_MultiVector& expolstresses                                      ///< nodal stresses
      );

      //@}

      /// determine proportions of element at origin
      void AxialMetricsAtOrigin(
          const LINALG::Matrix<NUMNOD_, NUMDIM_>& xrefe,  ///< (material/reference) element coords
          LINALG::Matrix<NUMDIM_, NUMDIM_>& jac0,         ///< Jacobian at origin
          LINALG::Matrix<NUMDIM_, 1>& metr0               ///< axial metrics at origin
      );

      /// determine metric coefficients in parametric/natural/local co-ordinate system
      static void LocalMetrics(
          const LINALG::Matrix<NUMDIM_, NUMDIM_>& jac, LINALG::Matrix<NUMDIM_, NUMDIM_>& metr);

      /// EXPERIMENTAL: setup of constant ANS data, second set of points
      ///
      /// \sa #sosh8_anssetup()
      void AnsSetup2(
          const LINALG::Matrix<NUMNOD_, NUMDIM_>& xrefe,  ///< material/reference element coords
          const LINALG::Matrix<NUMNOD_, NUMDIM_>& xcurr,  ///< spatial/current element coords
          std::vector<LINALG::Matrix<NUMDIM_, NUMNOD_>>**
              deriv_sp,  ///< derivs eval. at all sampling points
          std::vector<LINALG::Matrix<NUMDIM_, NUMDIM_>>& jac_sps,  ///< jac at all sampling points
          std::vector<LINALG::Matrix<NUMDIM_, NUMDIM_>>&
              jac_cur_sps,  ///< current jac at all sampling points
          LINALG::Matrix<NUMANS_ * NUMSP2ND_, NUMDISP_>& B_ans_loc  ///< modified B
      );

      /// EXPERIMENTAL: setup of constant ANS data, third set of points
      ///
      /// \sa #sosh8_anssetup()
      void AnsSetup3(
          const LINALG::Matrix<NUMNOD_, NUMDIM_>& xrefe,  ///< material/reference element coords
          const LINALG::Matrix<NUMNOD_, NUMDIM_>& xcurr,  ///< spatial/current element coords
          std::vector<LINALG::Matrix<NUMDIM_, NUMNOD_>>**
              deriv_sp,  ///< derivs eval. at all sampling points
          std::vector<LINALG::Matrix<NUMDIM_, NUMDIM_>>& jac_sps,  ///< jac at all sampling points
          std::vector<LINALG::Matrix<NUMDIM_, NUMDIM_>>&
              jac_cur_sps,  ///< current jac at all sampling points
          LINALG::Matrix<NUMANS_ * NUMSP3RD_, NUMDISP_>& B_ans_loc  ///< modified B
      );

      /// @name EAS functions
      //@{

      /// set-up of EAS data
      void EasInit();

      /// retrieve EAS parameters and incremental update of them
      template <int NUMEAS_T>
      static void EasUpdateIncrementally(
          Epetra_SerialDenseMatrix*&
              oldfeas,  ///< EAS constraint \f$f_{EAS}^{k}\f$ of last iteration
          Epetra_SerialDenseMatrix*& oldKaainv,  ///< inverted tangent k_aa^{-1} of last iteration
          Epetra_SerialDenseMatrix*& oldKad,     ///< tangent k_ad of last iteration
          Epetra_SerialDenseMatrix*& oldKap,     ///< tangent k_ap of last iteration
          Teuchos::RCP<Epetra_SerialDenseVector>&
              feas,  ///< current EAS constraint \f$f_{EAS}^{k+1}\f$
          Teuchos::RCP<Epetra_SerialDenseMatrix>& Kaa,  ///< current tangent k_aa
          Teuchos::RCP<Epetra_SerialDenseMatrix>& Kad,  ///< current tangent k_ad
          Teuchos::RCP<Epetra_SerialDenseMatrix>& Kap,  ///< current tangent k_ap
          Epetra_SerialDenseMatrix*& alpha,             ///< EAS parameters
          Teuchos::RCP<Epetra_SerialDenseMatrix>& M,    ///< EAS shape functions
          DRT::Container& data,                         ///< (input) data
          const LINALG::Matrix<NUMDISP_, 1>& dispi,     ///< current residual displacements
          const LINALG::Matrix<NUMPRES_, 1>& presi      ///< current residual pressures
      );

      /// push parametric EAS GL strain to material/reference configuration
      template <int NUMEAS_T>
      static void EasMaterialiseShapeFcts(Teuchos::RCP<Epetra_SerialDenseMatrix>
                                              M,  ///< EAS shape functions in material configuration
          const double& detJ0,  ///< material-to-parameter Jacobian determinant at origin
          const double& detJ,   ///< material-to-parameter Jacobian determinant at GP
          const LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>&
              T0invT,  ///< parameter-to-material transformation for 2-tensors
          const Epetra_SerialDenseMatrix& Mloc  ///< parametric EAS shape function eval. at GP
      );

      /// add EAS strain contribution
      template <int NUMEAS_T>
      static void EasAddStrain(
          LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& glstrain,  ///< Green-Lagrange strain vector
          const Teuchos::RCP<Epetra_SerialDenseMatrix>&
              M,                                 ///< EAS shape functions in material configuration
          const Epetra_SerialDenseMatrix* alpha  ///< EAS parameters
      );

      /// build EAS constraint and its tangents
      template <int NUMEAS_T>
      static void EasConstraintAndTangent(Teuchos::RCP<Epetra_SerialDenseVector>&
                                              feas,  ///< current EAS constraint \f$f_{EAS}^{k+1}\f$
          Teuchos::RCP<Epetra_SerialDenseMatrix>& Kaa,       ///< current tangent k_aa
          Teuchos::RCP<Epetra_SerialDenseMatrix>& Kad,       ///< current tangent k_ad
          Teuchos::RCP<Epetra_SerialDenseMatrix>& Kap,       ///< current tangent k_ap
          const LINALG::Matrix<NUMDIM_, NUMDIM_>& defgradD,  ///< compatible deformation gradient
          const LINALG::Matrix<NUMDIM_, NUMDIM_>&
              invrgtstrD,  ///< compatible inverse right stretch tensor
          const LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>&
              rcgbyrgtstr,           ///< right stretch tensor diff'd w.r.t right CG
          const double& detdefgrad,  ///< determinant of deformation gradient
          const LINALG::Matrix<NUMDFGR_, 1>&
              tinvdefgrad,  ///< vector of transposed inverse deformation gradient
          const LINALG::Matrix<NUMDFGR_, NUMDFGR_>& WmT,  ///< the matrix ( fv . fv^T + Wm )
          const LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>&
              cmat,  ///< (isotropic) elasticity matrix
          const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>&
              stress,                 ///< 2nd Piola-Kirchhoff stress vector
          const double& effpressure,  ///< effective pressure at Gauss point
          const double& detJ_w,       ///< integration factor
          const LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDISP_>&
              cb,  ///< elasticity matrix times B-operator
          const LINALG::Matrix<NUMDFGR_, NUMDISP_>&
              defgradbydisp,  ///< deformation gradient diff'd w.r.t. displacements
          const LINALG::Matrix<NUMPRES_, 1>& prshfct,  ///< effective pressure shape functions
          const Teuchos::RCP<Epetra_SerialDenseMatrix>&
              M  ///< EAS shape functions in material configuration
      );

      /// static condensation
      template <int NUMEAS_T>
      static void EasCondensation(
          LINALG::Matrix<NUMDISP_, 1>* force,               ///< element internal force vector
          LINALG::Matrix<NUMDISP_, NUMDISP_>* stiffmatrix,  ///< element stiffness matrix
          LINALG::Matrix<NUMDISP_, NUMPRES_>* gradmatrix,   ///< element gradient matrix
          LINALG::Matrix<NUMPRES_, 1>* incomp,              ///< incompressibility residual
          LINALG::Matrix<NUMPRES_, NUMDISP_>* dargmatrix,  ///< 'transposed' element gradient matrix
          LINALG::Matrix<NUMPRES_, NUMPRES_>* stabmatrix,  ///< element stabilisation matrix
          Epetra_SerialDenseMatrix*&
              oldfeas,  ///< EAS constraint \f$f_{EAS}^{k}\f$ of last iteration
          Epetra_SerialDenseMatrix*& oldKaainv,  ///< inverted tangent k_aa^{-1} of last iteration
          Epetra_SerialDenseMatrix*& oldKad,     ///< tangent k_ad of last iteration
          Epetra_SerialDenseMatrix*& oldKap,     ///< tangent k_ap of last iteration
          const Teuchos::RCP<Epetra_SerialDenseVector>&
              feas,  ///< current EAS constraint \f$f_{EAS}^{k+1}\f$
          const Teuchos::RCP<Epetra_SerialDenseMatrix>& Kaa,  ///< current tangent k_aa
          const Teuchos::RCP<Epetra_SerialDenseMatrix>& Kad,  ///< current tangent k_ad
          const Teuchos::RCP<Epetra_SerialDenseMatrix>& Kap   ///< current tangent k_ap
      );

      //@}

      /// @name Voigt vector/matrix converters
      //@{

      /// Convert 3x3 matrix of non-sym 2-tensor in 9x1 Voigt vector
      ///
      /// From
      ///\f[
      ///   [\boldsymbol{F}] = \left[\begin{array}{ccc}
      ///   F_{11} & F_{12} & F_{13} \\ F_{21} & F_{22} & F_{23} \\ F_{31} & F_{32} & F_{33}
      ///   \end{array}\right]
      ///\f]
      /// to
      ///\f[
      ///   \tilde{\mathbf{F}} = \left[\begin{array}{c}
      ///   F_{11} \\ F_{22} \\ F_{33} \\ F_{12} \\ F_{23} \\ F_{31} \\ F_{13} \\ F_{32} \\ F_{21}
      ///   \end{array}\right]
      ///\f]
      static void Matrix2TensorToVector9Voigt(
          LINALG::Matrix<NUMDFGR_, 1>& fvct,             ///< (out) 9x1 Voigt vector
          const LINALG::Matrix<NUMDIM_, NUMDIM_>& fmat,  ///< (in) 3x3 matrix
          const bool transpose = false                   ///< use transposed input 3x3 matrix
      );

      /// Convert 3x3 matrix of sym 2-tensor in 9-Voigt vector
      ///
      /// From
      ///\f[
      ///   [\boldsymbol{B}] = \left[\begin{array}{ccc}
      ///   B_{11} & B_{12} & B_{31} \\ B_{12} & B_{22} & B_{23} \\ B_{31} & B_{23} & B_{33}
      ///   \end{array}\right]
      ///\f]
      /// to if \i strain-like vector
      ///\f[
      ///   \tilde{\mathbf{B}} = \left[\begin{array}{c}
      ///   B_{11} \\ B_{22} \\ B_{33} \\ 2 B_{12} \\ 2 B_{23} \\ 2 B_{31}
      ///   \end{array}\right]
      ///\f]
      /// or to if \i stress-like vector
      ///\f[
      ///   \tilde{\mathbf{B}} = \left[\begin{array}{c}
      ///   B_{11} \\ B_{22} \\ B_{33} \\ B_{12} \\ B_{23} \\ B_{31}
      ///   \end{array}\right]
      ///\f]
      static void Matrix2TensorToVector6Voigt(
          LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& bvct,   ///< 6-Voigt vector
          const LINALG::Matrix<NUMDIM_, NUMDIM_>& bmat,  ///< 3x3 sym 2-tensor
          const VoigtType outvoigt6 = voigt6_strain      ///< type of 6-Voigt vector
      );

      /// convert 6-Voigt vector into matrix of 2-tensor
      ///
      /// opposite of #Matrix2TensorToVector6Voigt
      static void Vector6VoigtToMatrix2Tensor(
          LINALG::Matrix<NUMDIM_, NUMDIM_>& bmat,             ///< 3x3 sym 2-tensor
          const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& bvct,  ///< 6-Voigt vector
          const VoigtType invoigt6 = voigt6_strain            ///< type of 6-Voigt vector
      );

      /// Derivative of inverse of non-symmetric 2-tensor with respect to itself
      ///
      /// In index notation:
      ///\f[
      ///   \frac{\partial (F^{-1})_{ij}}{\partial F_{kl}}
      ///   = -(F^{-1})_{ik} \cdot (F^{-1})_{lj}
      ///\f]
      static void InvVector9VoigtDiffByItself(LINALG::Matrix<NUMDFGR_, NUMDFGR_>& invfderf,
          const LINALG::Matrix<NUMDIM_, NUMDIM_>& invfmat,
          const bool transpose = false  ///< use transposed input 3x3 matrix
      );

      /// Derivative of inverse of symmetric 2-tensor with respect to itself
      ///
      /// In index notation:
      ///\f[
      ///   \frac{\partial (F^{-1})_{ij}}{\partial F_{kl}}
      ///   = -1/2*\big( (F^{-1})_{ik} \cdot (F^{-1})_{lj}
      ///                + (F^{-1})_{il} \cdot (F^{-1})_{kj} \big)
      ///\f]
      static void InvVector6VoigtDiffByItself(
          LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& invfderf,  ///< 6x6 derivative
          const LINALG::Matrix<NUMDIM_, NUMDIM_>& invfmat  ///< 3x3 inverse symmetric 2-tensor
      );

      /// 2nd derivative of inverse of symmetric 2-tensor with respect to itself
      ///
      /// In index notation:
      ///\f[
      ///   \begin{array}{rl}
      ///      \dfrac{\partial (C^{-1})^{CD}}{\partial C_{EF} \, \partial C_{GH}}
      ///      = \frac{1}{4} \Big(
      ///       & \big( (C^{-1})^{CG}(C^{-1})^{HE} + (C^{-1})^{CH}(C^{-1})^{GE} \big)(C^{-1})^{FD}
      ///    \\ &  + (C^{-1})^{CE}\big( (C^{-1})^{FG}(C^{-1})^{HD} + (C^{-1})^{FH}(C^{-1})^{GD}
      ///    \big)
      ///    \\ &  + \big( (C^{-1})^{CG}(C^{-1})^{HF} + (C^{-1})^{CH}(C^{-1})^{GF}
      ///    \big)(C^{-1})^{ED}
      ///    \\ &  + (C^{-1})^{CF}\big( (C^{-1})^{EG}(C^{-1})^{HD} + (C^{-1})^{EH}(C^{-1})^{GD}
      ///    \big)
      ///       \Big)
      ///   \end{array}
      ///\f]
      ///
      static void InvVector6VoigtTwiceDiffByItself(
          LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D * MAT::NUM_STRESS_3D>&
              invbvdderb,                              ///< 6x36 matrix
          const LINALG::Matrix<NUMDIM_, NUMDIM_>& ibt  ///< 3x3 inverted symmetric 2-tensor
      );

      /// Derivative of square of symmetric 2-tensor by itself
      ///
      /// In index notation:
      ///\f[
      ///   \frac{\partial F_{im} F_{mj}}{\partial F_{kl}}
      ///   = \delta_{ik}\cdot F_{lj} + \delta_{jl}\cdot F_{ik}
      ///\f]
      static void SqVector6VoigtDiffByItself(LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>&
                                                 sqfderf,  ///< diff. squared matrix
          const LINALG::Matrix<NUMDIM_, NUMDIM_>& fmat,    ///< symmetric 2-tensor
          const VoigtType outvoigt6 = voigt6_strain        ///< 6-Voigt row type of resulting matrix
      );

      /// Derivative of "square" of non-symmetric 2-tensor by itself
      ///
      /// cf. #SqVector6VoigtDiffByItself
      ///
      /// In compact tensor notation
      ///\f[
      ///   \big( \boldsymbol{F}^T \cdot \boldsymbol{F} \big)_{,\boldsymbol{F}}
      ///\f]
      static void SqVector9VoigtDiffByItself(
          LINALG::Matrix<NUMDFGR_, NUMDFGR_>& sqfderf,   ///< diff. squared matrix
          const LINALG::Matrix<NUMDIM_, NUMDIM_>& fmat,  ///< non-symmetric 2-tensor
          const bool transpose = false                   ///< transpose non-symmetric 2-tensor
      );

      /// 2nd derivative of square of symmetric 2-tensor with respect to itself
      ///
      /// In index notation
      ///\f[
      /// \begin{array}{rl}
      ///    \dfrac{\partial (U_{EA} U_{AF})}{\partial U_{GH}\, \partial U_{IJ}}
      ///    = \frac{1}{4}\Big(
      /// &
      ///    \delta_E{}^G\delta^{HI}\delta_F{}^J+\delta_F{}^H\delta_E{}^I\delta^{GJ}
      ///    +\delta_E{}^G\delta^{HJ}\delta_F{}^I+\delta_F{}^H\delta_E{}^J\delta^{GI}
      /// \\ &
      ///    +\delta_E{}^H\delta^{GI}\delta_F{}^J+\delta_F{}^G\delta_E{}^I\delta^{HJ}
      ///    +\delta_E{}^H\delta^{GJ}\delta_F{}^I+\delta_F{}^G\delta_E{}^J\delta^{HI}
      ///     \Big)
      /// \end{array}
      ///\f]
      static void SqVector6VoigtTwiceDiffByItself(
          LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D * MAT::NUM_STRESS_3D>&
              sqfdderf,                                 ///< 2nd derivative
          const LINALG::Matrix<NUMDIM_, NUMDIM_>& fmat  ///< symmetric 2-tensor (not needed)
      );

      /// 2nd derivative of square of symmetric 2-tensor with respect to itself
      /// -- sparse storage version
      ///
      /// cf. #SqVector6VoigtTwiceDiffByItself bit with sparse storage of 6-tensor
      static void SqVector6VoigtTwiceDiffByItself(int* isqfdderf,  ///< indices of non-zero entries
          LINALG::Matrix<MAT::NUM_STRESS_3D, 6>&
              sqfdderf  ///< 2nd derivative data (sparsely stored)
      );

      /// Build 6x9 Voigt matrix for 2-tensor-dot-2-tensor product
      ///
      /// This method builds a matrix $\mathbf{B}$ which allows to perform the
      /// following operation utilising a Voigt matrix.
      ///
      /// The operation in
      /// (1) compact tensor notation
      ///\f[
      ///   \boldsymbol{A} = \boldsymbol{A}^T = \boldsymbol{B}^T \cdot \boldsymbol{C}
      ///\f]
      /// or (2) in index tensor notation
      ///\f[
      ///   A_{ij} = A_{ji} = B_{ki} \cdot C_{kj}
      ///\f]
      /// or (3) in Voigt 9-vector by 6x9-matrix product, i.e.
      ///\f[
      ///   \mathbf{A} = \mathbf{B} \; \tilde{\mathbf{C}}
      ///\f]
      static void Matrix2TensorToMatrix6x9Voigt(
          LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDFGR_>& bm,  ///< (out) 6x9 Voigt matrix
          const LINALG::Matrix<NUMDIM_, NUMDIM_>& bt,        ///< (in) 3x3 matrix
          const bool transpose = true                        ///< transpose entries
      );

      /// Build 6x6-Voigt matrix to multiply 6-Voigt with which is equivalent to
      /// multiplication with non-sym 2-tensor \f$\boldsymbol{F}\f$ from left and right
      /// of a symmetric 2-tensor \f$\boldsymbol{e}\f$
      ///
      /// The operation in
      /// (1) compact tensor notation
      ///\f[
      ///    \boldsymbol{E} = \boldsymbol{F}^T \cdot \boldsymbol{e} \cdot \boldsymbol{F}
      ///\f]
      /// or (2) in index tensor notation
      ///\f[
      ///    E_{AB} = F_A{}^a  \cdot e_{ab}  \cdot F^b{}_B
      ///\f]
      /// or  (3) in Voigt 6-vector by 6x6-matrix product, i.e.
      ///\f[
      ///    \mathbf{E} = \mathbf{B} \; \mathbf{e}
      ///\f]
      ///
      /// This operation is also known as 'pull back' or 'base transformation', respectively.
      static void Matrix2TensorToLeftRightProductMatrix6x6Voigt(
          LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& bm,  ///< (out) 6x6 Voigt matrix
          const LINALG::Matrix<NUMDIM_, NUMDIM_>& bt,  ///< (in) 3x3 matrix of 2-tensor
          const bool transpose,                        ///< 3x3 input matrix is transposed
          const VoigtType rowvoigt6,  ///< 6-Voigt vector layout on rows of 6x6 matrix
          const VoigtType colvoigt6   ///< 6-Voigt vector layout on columns of 6x6 matrix
      );

      //@}

      /// determine (inverse) right stretch tensor
      ///
      /// The right stretch tensor is obtained by polar decomposition
      /// of the right Cauchy-Green 2-tensor in 3 dimensions.
      ///
      /// About:
      ///   Polar decomposition is often applied to the material deformation
      /// tensor F, i.e.
      ///      F = R . U = v . R
      /// in which R is the rotation matrix (two-point tensor), U the material
      /// stretch tensor (refers to the undeformed configuration) and v the
      /// spatial stretch tensor (refers to the deformed configuration).
      ///   This polar decomposition is also applied to the isoparametric
      /// Jacobian tensor (mapping quantities in parameter space to material
      /// configuration). However, the local variables are denoted for the
      /// case of a deformation tensor.
      ///
      /// References:
      /// [1] A. Hoger & D.E. Carlson, "Determination of the stretch and
      ///        rotation in the polar decomposition of the deformation
      ///        gradient", Quart. Appl. Math., 42(2):113-117, 1984.
      /// [2] G.A. Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
      ///        esp. Section 2.6
      static void StretchTensor(double* detut,        ///< determinant of material stretch tensor
          LINALG::Matrix<NUMDIM_, NUMDIM_>* ut,       ///< material stretch tensor
          LINALG::Matrix<NUMDIM_, NUMDIM_>* invut,    ///< inverse material stretch tensor
          const LINALG::Matrix<NUMDIM_, NUMDIM_>& ct  ///< right Cauchy-Green tensor
      );

      /// Spectral decomposition of symmetric 2-tensor calculated
      /// iteratively using Jacobi's method
      ///
      /// References:
      /// [1] I.N. Bronstein & K.A. Semendjajew, "Taschenbuch der
      ///     Mathematik", Teubner, 25.ed, 1991.
      ///     esp. p. 741 ff.
      ///
      /// \return error: 0=success, 1=failure
      static int SymSpectralDecompJacIter(
          LINALG::Matrix<NUMDIM_, NUMDIM_>& ew,  ///< diagional matrix of eigenvalues
          LINALG::Matrix<NUMDIM_, NUMDIM_>& ev,  ///< orthornormal eigenvectors (direction cosines)
          const LINALG::Matrix<NUMDIM_, NUMDIM_>& at,  ///< <i>symmetric</i> input matrix
          const double itertol,                        ///< tolerance
          const int itermax                            ///< maximally allowd iteration steps
      );

      /// @name Local assemble of displacement and pressure quantities
      //@{

      /// Extract 24x1 displacement and 8x1 pressure
      /// of 32x1 displacement-pressure vector
      ///
      /// From 32x1 general element state vector
      ///\f[
      ///   \left[\begin{array}{cccccc}
      ///   \ldots & u_X^k & u_Y^k & u_Z^k & p^k & \ldots
      ///   \end{array}\right]
      ///\f]
      /// to 24x1 element displacements
      ///\f[
      ///   \left[\begin{array}{ccccc}
      ///   \ldots & u_X^k & u_Y^k & u_Z^k & \ldots
      ///   \end{array}\right]
      ///\f]
      /// and to 8x1 element pressures
      ///\f[
      ///   \left[\begin{array}{ccc}
      ///   \ldots & p^k & \ldots
      ///   \end{array}\right]
      ///\f]
      static void ExtractDispAndPres(
          std::vector<double>& mystat,          ///< 32x1 element displacement-pressure vector
          LINALG::Matrix<NUMDISP_, 1>& mydisp,  ///< 24x1 element displacement vector
          LINALG::Matrix<NUMPRES_, 1>& mypres   ///< 8x1 element pressure vector
      );

      /// Build 32x32 element matrix
      ///
      /// Following the pattern given at #ExtractDispAndPres()
      static void BuildElementMatrix(LINALG::Matrix<NUMDOF_, NUMDOF_>* mat,  ///< 32x32 matrix
          const LINALG::Matrix<NUMDISP_, NUMDISP_>*
              matdd,  ///< 24x24 sub-matrix, if NULL then 0s are inserted
          const LINALG::Matrix<NUMDISP_, NUMPRES_>*
              matdp,  ///< 24x8 sub-matrix, if NULL then 0s are inserted
          const LINALG::Matrix<NUMPRES_, NUMDISP_>*
              matpd,  ///< 8x24 sub-matrix, if NULL then transpose of #matdp is inserted or 0s
          const LINALG::Matrix<NUMPRES_, NUMPRES_>*
              matpp  ///< 8x8 sub-matrix, if NULL then 0s are inserted
      );

      /// Build 32x1 element vector
      ///
      /// Following the pattern given at #ExtractDispAndPres()
      static void BuildElementVector(LINALG::Matrix<NUMDOF_, 1>* vct,  ///< 32x1 vector
          const LINALG::Matrix<NUMDISP_, 1>*
              vctd,  ///< 24x1 sub-vector, if NULL then 0s are inserted
          const LINALG::Matrix<NUMPRES_, 1>* vctp  ///< 8x1 sub-vector, if NULL then 0s are inserted
      );

      /// Assemble global volume
      static void AssembleVolume(Teuchos::ParameterList& params,  ///< parameter list for in 'n' out
          const double& elevol                                    ///< current element volume
      );

      //@}

      /// @name Methods for debugging
      //@{

      /// Print files for visualising with Gnuplot
      ///
      /// The routine is for <b>debugging only</b> and should be used
      /// with <b>extreme care</b>. If everything works out well, you achieve
      /// two files: xxx.sosh8p8.gplt and xxx.sosh8p8.txt. The latter
      /// holds the data for the Gnuplot GPLT driver file. Execution
      /// of the GPLT file is going to produce nice graphs of
      /// internal force/incompressibility conditions and tangents
      /// with respect to the element displacements and pressures thereof.
      static void GnuplotOut(Teuchos::ParameterList& params,  ///< parameter list for in 'n' out
          std::vector<double>&
              state,  ///< current state vector, i.e. displacements and pressure DOFs
          LINALG::Matrix<NUMDOF_, 1>&
              resid,  ///< current internal force / incompressibility residual
          LINALG::Matrix<NUMDOF_, NUMDOF_>& tangent  ///< current tangent of inter force WRT state
      );

      //@}

      //! Calculate the STC matrix
      virtual void CalcSTCMatrix(LINALG::Matrix<NUMDOF_, NUMDOF_>& elemat1,
          const INPAR::STR::STC_Scale stc_scaling, const int stc_layer, std::vector<int>& lm,
          DRT::Discretization& discretization, bool calcinverse);
    };  // class So_sh8p8


  }  // namespace ELEMENTS
}  // namespace DRT


/*----------------------------------------------------------------------*/
/* definitions of above declared template functions */
#include "so_sh8p8_eas.H"


#endif  // #ifndef SO_SH8P8_H
