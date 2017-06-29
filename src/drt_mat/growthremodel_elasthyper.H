/*----------------------------------------------------------------------*/
/*!
\file growthremodel_elasthyper.H

\brief This file is used to manage the homogenized constraint mixture during growth and remodeling

\level 3

\maintainer Fabian Braeu
*/
/*----------------------------------------------------------------------*/
#ifndef GROWTHREMODEL_ELASTHYPER_H
#define GROWTHREMODEL_ELASTHYPER_H


#include "../headers/FADmatrix_utils.H"
#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"

// forward declaration due to avoid header definition
namespace MAT
{
  namespace ELASTIC
  {
    class Summand;
    class RemodelFiber;
  }

  // forward declaration
  class GrowthRemodel_ElastHyper;

  namespace PAR
  {
    class GrowthRemodel_ElastHyper : public Parameter
    {

      friend class MAT::GrowthRemodel_ElastHyper;

    public:

      /// standard constructor
      ///
      /// This constructor recursively calls the constructors of the
      /// parameter sets of the hyperelastic summands.
      GrowthRemodel_ElastHyper(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~GrowthRemodel_ElastHyper() { ; }

      /// @name material parameters
      //@{

      /// length of remodelfiber material list
      const int nummat_remodelfiber_;

      /// length of 3d elastin matrix material list
      const int nummat_elastiniso_;

      /// length of membrane elastin matrix material list
      const int nummat_elastinmem_;

      /// the list of remodelfiber material IDs
      const std::vector<int>* matids_remodelfiber_;

      /// the list of 3d elastin matrix material IDs
      const std::vector<int>* matids_elastiniso_;

      /// the list of membrane elastin matrix material IDs
      const std::vector<int>* matids_elastinmem_;

      /// material ID of growth penalty material
      const int matid_penalty_;

      /// initial mass fraction of elastin in constraint mixture
      const double init_w_el_;

      /// material mass density
      const double density_;

      /// circumferential prestretch of elastin matrix
      const double lamb_prestretch_cir_;

      /// axial prestretch of elastin matrix
      const double lamb_prestretch_ax_;

      /// reference wall thickness of the idealized cylindrical aneurysm
      const double t_ref_;

      /// mean blood pressure
      const double p_mean_;

      /// inner radius of the idealized cylindrical aneurysm
      const double ri_;

      /// flag for turning elastin damage on or off: 1: damage on; 0: damage off
      const int damage_;

      /// flag to decide what type of collagen growth is used: 1: anisotropic growth; 0: isotropic growth
      const int growthtype_;

      /// flag to decide what type of local time integration scheme is used: 1: Backward Euler Method; 0: Forward Euler Method
      const int loctimeint_;

      /// Indicator whether Hex or Membrane elements are used ( Membrane: 1, Hex: Everything else )
      const int membrane_;

      /// Indicator whether geometry is a cylinder and that we want to calculate the AXI, CIR and RAD-directions
      /// with the help of the location of the center of each element in the reference configuration
      /// (cylinder_ == 1: cylinder aligned in x-direction; cylinder_ == 2: cylinder aligned in y-direction
      /// cylinder_ == 3: cylinder aligned in z-direction)
      const int cylinder_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class GrowthRemodel_ElastHyper

  }  // namespace PAR

class GrowthRemodel_ElastHyperType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "GrowthRemodel_ElastHyperType"; }

  static GrowthRemodel_ElastHyperType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static GrowthRemodel_ElastHyperType instance_;
};


  /*----------------------------------------------------------------------*/
  class Material;

  class GrowthRemodel_ElastHyper : public So3Material
  {
  public:

    /// construct empty material object
    GrowthRemodel_ElastHyper();

    /// construct the material object given material parameters
    explicit GrowthRemodel_ElastHyper(MAT::PAR::GrowthRemodel_ElastHyper* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const { return GrowthRemodel_ElastHyperType::Instance().UniqueParObjectId(); }

    /// \brief Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by UniqueParObjectId() which will then
    /// identify the exact class on the receiving processor.
    ///
    /// \param data (in/out): char vector to store class information
    virtual void Pack(DRT::PackBuffer& data) const;

    /// \brief Unpack data from a char vector into this class
    ///
    /// The vector data contains all information to rebuild the
    /// exact copy of an instance of a class on a different processor.
    /// The first entry in data has to be an integer which is the unique
    /// parobject id defined at the top of this file and delivered by
    /// UniqueParObjectId().
    ///
    /// \param data (in) : vector storing all data to be unpacked into this
    ///                    instance.
    virtual void Unpack(const std::vector<char>& data);

    //@}

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem)
    { if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
       dserror("element and material kinematics are not compatible");}

    /// material type
    virtual INPAR::MAT::/// axial direction of the cylinder (cylidner axis)
    MaterialType MaterialType() const { return INPAR::MAT::m_growthremodel_elasthyper; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new GrowthRemodel_ElastHyper( *this ) ); }

    /// material mass density
    virtual double Density() const { return params_->density_; }

    /// hyperelastic stress response plus elasticity tensor
    virtual void Evaluate(const LINALG::Matrix<3,3>* defgrd,  ///< Deformation gradient
                          const LINALG::Matrix<6,1>* glstrain,///< Green-Lagrange strain
                          Teuchos::ParameterList& params,     ///< Container for additional information
                          LINALG::Matrix<6,1>* stress,        ///< 2nd Piola-Kirchhoff stresses
                          LINALG::Matrix<6,6>* cmat,          ///< Constitutive matrix
                          const int eleGID);                  ///< Element ID


    /// setup
    virtual void Setup(int const numgp,DRT::INPUT::LineDefinition* linedef);

    /// update
    virtual void Update(LINALG::Matrix<3,3> const& defgrd,    ///< Deformation gradient
                        int const gp,                         ///< Current Gauß-Point
                        Teuchos::ParameterList& params,       ///< Container for additional information
                        int const eleGID);                    ///< Element ID

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// hyperelastic stress response plus elasticity tensor for membrane element (membrane formulation)
    void EvaluateMembrane(LINALG::Matrix<3,3> const& defgrd_glob,  ///< Deformation gradient in global coordinates
                          double & rcg33,                          ///< Adapted stretch in thickness direction
                          Teuchos::ParameterList& params,          ///< Container for additional information
                          LINALG::Matrix<3,3> & pk2M_glob,         ///< 2nd Piola-Kirchhoff stress global coordinates
                          LINALG::Matrix<6,6> & cmat_glob,         ///< Elasticity tensor in global coordinates
                          const int eleGID);                       ///< Element ID

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string,int>& names);

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

  private:
    /// Setup circumferential, radial and axial structural tensor
    void SetupAxiCirRadStructuralTensor(DRT::INPUT::LineDefinition* linedef);

    /// Setup prestretch (optional: setup element axi-, circ-, and rad-directions) for 3D elements
    void SetupGR3D(LINALG::Matrix<3,3> const* const defgrd, ///< Deformation gradient
                   Teuchos::ParameterList& params,          ///< Container for additional information
                   const double dt,                         ///< Time step size
                   const int gp,                            ///< Current Gauß-Point
                   const int eleGID);                       ///< Element ID

    /// Setup prestretch (optional: setup element axi-, circ-, and rad-directions) for 2D elements -> membrane
    void SetupGR2D(Teuchos::ParameterList& params,          ///< Container for additional information
                   const double dt,                         ///< Time step size
                   const int gp);                           ///< Current Gauß-Point

    /// Calculates AXI, CIR and RAD structural tensors and sets new fiber directions in the case of a cylinder
    void SetupAxiCirRadCylinder(LINALG::Matrix<1,3> gprefecoord,
                                double const dt);

    /// Setup anisotropic growth tensors. Here we assume that the growth direction corresponds with
    /// the radial/ thickness direction
    void SetupAnisoGrowthTensors();

    /// Read AXI CIR RAD direction
    void ReadDir(DRT::INPUT::LineDefinition* linedef,
                 std::string specifier,
                 LINALG::Matrix<3,1>& dir);

    /// Evaluate Prestretches
    void EvaluatePrestretch(LINALG::Matrix<3,3> const * const defgrd,
                            int const gp,
                            int const eleGID);

    /// Internal Newton to implicitly solve for current mass density and inelastic remodeling stretch of each fiber family
    void SolveForRhoLambr(LINALG::SerialDenseMatrix & K_T,            ///< Tangent stiffness matrix of internal Newton
                          LINALG::Matrix<3,3> & FgM,                  ///< Growth deformation gradient
                          LINALG::Matrix<3,3> & iFgM,                 ///< Inverse growth deformation gradient
                          LINALG::Matrix<3,3> & dFgdrhoM,             ///< Derivative of growth deformation gradient w.r.t. the mass density
                          LINALG::Matrix<3,3> & diFgdrhoM,            ///< Derivative of inverse growth deformation gradient w.r.t. the mass density
                          LINALG::Matrix<3,3> const * const defgrd,   ///< Deformation gradient
                          double const& dt,                           ///< Time step size
                          int const gp,                               ///< Current Gauß-Point
                          int const eleGID);                          ///< Element ID

    /// Internal Newton to implicitly solve for current mass density and inelastic remodeling stretch of each fiber family
    void SolveFordrhodCdlambrdC(std::vector<LINALG::Matrix<1,6> > & drhodC,     ///< Derivative of collagen mass density w.r.t. right Cauchy Green tensor
                                std::vector<LINALG::Matrix<1,6> > & dlambrdC,   ///< Derivative of inelastic remodel stretch w.r.t. right Cauchy Green tensor
                                LINALG::Matrix<1,6> & sum_drhodC,               ///< Sum of derivatives of all collagen mass density w.r.t. right Cauchy Green tensor
                                LINALG::SerialDenseMatrix & K_T,                ///< Tangent stiffness matrix of internal Newton
                                LINALG::Matrix<3,3> const& iFgM,                ///< Inverse growth deformation gradient
                                LINALG::Matrix<3,3> const * const defgrd,       ///< Deformation gradient
                                double const& dt,                               ///< Time step size
                                int const gp,                                   ///< Current Gauß-Point
                                int const eleGID) const;                        ///< Element ID

    /// Evaluates some kinematic quantities which are used in stress and elasticity tensor calculation
    void EvaluateKinQuantElast(LINALG::Matrix<3,3> const* const defgrd,   ///< Deformation gradient
                               LINALG::Matrix<3,3> const& iFinM,          ///< Inverse inelastic deformation gradient
                               LINALG::Matrix<6,1> & iCinv,               ///< Inverse inelastic right Cauchy-Green tensor
                               LINALG::Matrix<6,1> & iCinCiCinv,          ///< C_{in}^{-1} * C * C_{in}^{-1}
                               LINALG::Matrix<6,1> & iCv,                 ///< Inverse right Cauchy-Green tensor
                               LINALG::Matrix<3,3> & iCinCM,              ///< C_{in}^{-1} * C
                               LINALG::Matrix<3,3> & iFinCeM,             ///< F_{in}^{-1} * C_e
                               LINALG::Matrix<9,1> & CiFin9x1,            ///< C * F_{in}^{-1}
                               LINALG::Matrix<9,1> & CiFinCe9x1,          ///< C * F_{in}^{-1} * C_e
                               LINALG::Matrix<9,1> & CiFiniCe9x1,         ///< C * F_{in}^{-1} * C_e^{-1}
                               LINALG::Matrix<3,1> & prinv,               ///< Principal invariants of elastic right Cauchy-Green tensor
                               int const gp) const;                       ///< Current Gauß-Point

    /// calculate principal invariants
    void InvariantsPrincipal(LINALG::Matrix<3,1> & prinv,                 ///< principal invariants
                             LINALG::Matrix<6,1> const& C_strain) const;  ///< symmetric Cartesian 2-tensor in strain-like 6-Voigt notation

    /// calculate modified invariants
    void InvariantsModified(LINALG::Matrix<3,1> & modinv,             ///< modified invariants
                            LINALG::Matrix<3,1> const& prinv) const;  ///< principal invariants

    /// calculates the derivatives of the hyperelastic laws with respect to the (modified) invariants
    void EvaluateInvariantDerivatives(LINALG::Matrix<3,1> const& prinv,   ///< Principal invariants of right Cauchy-Green tensor
                                      LINALG::Matrix<3,1> & dPIw,         ///< First derivative with respect to invariants weighted with the corresponding volume fraction
                                      LINALG::Matrix<6,1> & ddPIIw,       ///< Second derivative with respect to invariants weighted with the corresponding volume fraction
                                      int const gp,                       ///< Current Gauß-Point
                                      int const eleGID) const;            ///< Element ID

    /// converts the derivatives with respect to the modified invariants in derivatives with respect to principal invariants.
    /// Uses the following conversions:
    ///\f[
    ///  \overline{I}_{\boldsymbol{C}} = J^{-2/3} I_{\boldsymbol{C}},
    ///\f]
    ///\f[
    ///  \overline{II}_{\boldsymbol{C}} = J^{-4/3} II_{\boldsymbol{C}},
    ///\f]
    ///\f[
    ///  J = \sqrt{III_{\boldsymbol{C}}}
    ///\f]
    void ConvertModToPrinc(LINALG::Matrix<3,1> const& prinv,    ///< principal invariants of right Cauchy-Green tensor
                           LINALG::Matrix<3,1> const& dPmodI,   ///< first derivative with respect to modified invariants
                           LINALG::Matrix<6,1> const& ddPmodII, ///< second derivative with respect to modified invariants
                           LINALG::Matrix<3,1> & dPI,           ///< first derivative with respect to invariants
                           LINALG::Matrix<6,1> & ddPII) const;  ///< second derivative with respect to invariants

    /// calculates the isotropic stress and elasticity tensor for coupled configuration
    void EvaluateIsotropicPrincElast(LINALG::Matrix<6,1> & stressisoprinc,    ///< 2nd Piola Kirchhoff stress
                                     LINALG::Matrix<6,6> & cmatisoprinc,      ///< Elasticity tensor
                                     LINALG::Matrix<6,1> const& iCinv,        ///< Inverse inelastic right Cauchy-Green tensor
                                     LINALG::Matrix<6,1> const& iCinCiCinv,   ///< C_{in}^{-1} * C * C_{in}^{-1}
                                     LINALG::Matrix<6,1> const& iCv,          ///< Inverse right Cauchy-Green tensor
                                     LINALG::Matrix<3,1> const& gamma,        ///< Factors for stress calculation
                                     LINALG::Matrix<8,1> const& delta) const; ///< Factors for elasticity tensor calculation

    /// Evaluate derivative of 2nd Piola Kirchhoff stress w.r.t. the inelastic deformation gradient
    void EvaluatedSdiFg(LINALG::Matrix<6,9> & dSdiFg,               ///< Output
                        LINALG::Matrix<3,1> const& gamma,           ///< Factors for stress calculation
                        LINALG::Matrix<8,1> const& delta,           ///< Factors for elasticity tensor calculation
                        LINALG::Matrix<3,3> const& iFinM,           ///< Inverse inelastic deformation gradient
                        LINALG::Matrix<3,3> const& iCinCM,          ///< C_{in}^{-1} * C
                        LINALG::Matrix<6,1> const& iCinv,           ///< Inverse inelastic right Cauchy-Green tensor
                        LINALG::Matrix<9,1> const& CiFin9x1,        ///< C * F_{in}^{-1}
                        LINALG::Matrix<9,1> const& CiFinCe9x1,      ///< C * F_{in}^{-1} * C_e
                        LINALG::Matrix<6,1> const& iCinCiCinv,      ///< C_{in}^{-1} * C * C_{in}^{-1}
                        LINALG::Matrix<9,1> const& CiFiniCe9x1,     ///< C * F_{in}^{-1} * C_e^{-1}
                        LINALG::Matrix<6,1> const& iCv,             ///< Inverse right Cauchy-Green tensor
                        LINALG::Matrix<3,3> const& iFinCeM,         ///< F_{in}^{-1} * C_e
                        int const gp) const;                        ///< Current Gauß-Point

    /// Evaluate additional terms for the elasticity tensor
    /// Additional terms are caused by implicit implementation of the evolution equation
    void EvaluateAdditionalCmat(LINALG::Matrix<6,6> & cmatadd,         ///< Additional elasticity tensor
                                LINALG::Matrix<3,3> const& diFgdrhoM,  ///< Derivative of inverse growth deformation gradient w.r.t. the mass density
                                LINALG::Matrix<1,6> const& sum_drhoC,  ///< Sum over all fiber families of the derivatives of the individual mass density w.r.t. right Cauchy Green tensor
                                LINALG::Matrix<6,9> const& dSdiFg,     ///< Derivative of 2nd Piola Kirchhoff stress w.r.t. the inverse inelastic growth deformation gradient
                                int const gp) const;                   ///< Current Gauß-Point

    /// calculate coefficients gamma and delta from partial derivatives w.r.t. invariants.
    /// The coefficients \f$\gamma_i\f$ and \f$\delta_j\f$ are based
    /// on the summand \f$\Psi(I_{\boldsymbol{C}},II_{\boldsymbol{C}},III_{\boldsymbol{C}})\f$
    /// in which the principal invariants of the right Cauchy-Green tensor \f$\boldsymbol{C}\f$
    /// are the arguments, cf. Holzapfel [1],
    ///\f[
    ///  I_{\boldsymbol{C}} = C_{AA},
    ///  \quad
    ///  II_{\boldsymbol{C}} = 1/2 \big( \mathrm{trace}^2(\boldsymbol{C}) - \mathrm{trace}(\boldsymbol{C}^2) \big),
    ///  \quad
    ///  III_{\boldsymbol{C}} = \det(\boldsymbol{C})
    ///\f]
    ///
    /// cf. Holzapfel [1], p. 216 (6.32) and p. 248
    ///\f[
    ///  \bold{S} = \gamma_1 \ \bold{Id} + \gamma_2 \ \bold{C} + \gamma_3 \ \bold{C}^{-1}
    ///\f]
    ///\f[
    /// \gamma_1 = 2\left( \frac{\partial \Psi}{\partial I_{\boldsymbol{C}}}
    ///          + I_{\boldsymbol{C}}\frac{\partial \Psi}{\partial II_{\boldsymbol{C}}} \right);
    ///\f]
    ///\f[
    /// \gamma_2 = -2\frac{\partial \Psi}{\partial II_{\boldsymbol{C}}};
    ///\f]
    ///\f[
    /// \gamma_3 = 2III_{\boldsymbol{C}} \frac{\partial \Psi}{\partial III_{\boldsymbol{C}}};
    ///\f]
    ///
    /// material constitutive tensor coefficients
    /// cf. Holzapfel [1], p. 261
    ///\f[
    ///  \mathbb{C} = \delta_1 \left( \bold{Id} \otimes \bold{Id} \right) + \delta_2 \left( \bold{Id} \otimes \bold{C} + \bold{C} \otimes \bold{Id} \right)
    ///  + \delta_3 \left( \bold{Id} \otimes \bold{C}^{-1} + \bold{C}^{-1} \otimes \bold{Id} \right) + \delta_4 \left( \bold{C} \otimes \bold{C} \right)
    ///  + \delta_5 \left( \bold{C} \otimes \bold{C}^{-1} + \bold{C}^{-1} \otimes \bold{C} \right) + \delta_6 \left( \bold{C}^{-1} \otimes \bold{C}^{-1} \right)
    ///  + \delta_7 \mathbb{P}
    ///\f]
    ///\f[
    /// \delta_1 = 4\left( \frac{\partial^2 \Psi}{\partial I_{\boldsymbol{C}}^2}
    ///          + 2 Ic \frac{\partial^2 \Psi}{\partial I_{\boldsymbol{C}} \partial II_{\boldsymbol{C}}}
    ///          + \frac{\partial \Psi}{\partial II_{\boldsymbol{C}}}
    ///          + I_{\boldsymbol{C}}^2 \frac{\partial^2 \Psi}{\partial II_{\boldsymbol{C}}^2} \right)
    ///\f]
    ///\f[
    /// \delta_2 = -4\left( \frac{\partial^2 \Psi}{\partial I_{\boldsymbol{C}} \partial II_{\boldsymbol{C}}}
    ///          + I_{\boldsymbol{C}}\frac{\partial^2 \Psi}{\partial II_{\boldsymbol{C}}^2} \right)
    ///\f]
    ///\f[
    /// \delta_3 = 4\left( III_{\boldsymbol{C}} \frac{\partial^2 \Psi}{\partial I_{\boldsymbol{C}} \partial III_{\boldsymbol{C}}}
    ///          + I_{\boldsymbol{C}} III_{\boldsymbol{C}} \frac{\partial^2 \Psi}{\partial IIc \partial III_{\boldsymbol{C}}} \right)
    ///\f]
    ///\f[
    /// \delta_4 = 4\frac{\partial^2 \Psi}{\partial II_{\boldsymbol{C}}^2}
    ///\f]
    ///\f[
    /// \delta_5 = -4 III_{\boldsymbol{C}}\frac{\partial^2 \Psi}{\partial II_{\boldsymbol{C}} \partial III_{\boldsymbol{C}}}
    ///\f]
    ///\f[
    /// \delta_6 = 4\left( III_{\boldsymbol{C}} \frac{\partial \Psi}{\partial III_{\boldsymbol{C}}}
    ///          + III_{\boldsymbol{C}}^2 \frac{\partial^2 \Psi}{\partial III_{\boldsymbol{C}}^2} \right)
    ///\f]
    ///\f[
    /// \delta_7 = -4 III_{\boldsymbol{C}} \frac{\partial \Psi}{\partial III_{\boldsymbol{C}}}
    ///\f]
    ///\f[
    /// \delta_8 = -4 \frac{\partial \Psi}{\partial II_{\boldsymbol{C}}}
    ///\f]
    void CalculateGammaDelta(LINALG::Matrix<3,1> & gamma,             ///< factors for stress calculation
                             LINALG::Matrix<8,1> & delta,             ///< factors for elasticity tensor calculation
                             LINALG::Matrix<3,1> const& prinv,        ///< principal invariants of right Cauchy-Green tensor
                             LINALG::Matrix<3,1> const& dPI,          ///< first derivative with respect to invariants
                             LINALG::Matrix<6,1> const& ddPII) const; ///< second derivative with respect to invariants

    /// Evaluate current individual volume of elastin
    void EvaluateElastinDamage(double const& dt_pre);   ///< time step size of first time step in simulation (initialization and prestretching)

    /// build growth deformation gradient
    void EvaluateGrowthDefGrad(LINALG::Matrix<3,3> & FgM,        ///< Growth deformation gradient
                               LINALG::Matrix<3,3> & iFgM,       ///< Inverse growth deformation gradient
                               LINALG::Matrix<3,3> & dFgdrhoM,   ///< Derivative of growth deformation gradient w.r.t. the mass density
                               LINALG::Matrix<3,3> & diFgdrhoM,  ///< Derivative of inverse growth deformation gradient w.r.t. the mass density
                               const int gp);                    ///< Current Gauß-Point

    /// Evaluates stress and cmat during prestressing
    void EvaluateStressCmatIso(LINALG::Matrix<3,3> const* const defgrd,     ///< Deformation gradient
                               LINALG::Matrix<3,3> const& iFinM,            ///< Inverse inelastic deformation gradient
                               LINALG::Matrix<6,1> & stressiso,             ///< Isotropic stress tensor
                               LINALG::Matrix<6,6> & cmatiso,               ///< Isotropic stiffness matrix
                               LINALG::Matrix<6,9> & dSdiFg,                ///< Derivative of 2nd Piola Kirchhoff stress w.r.t. the inverse inelastic growth deformation gradient
                               int const gp,                                ///< Current Gauß-Point
                               int const eleGID) const;                     ///< Element ID

    /// Contribution of membrane material to 2nd Piola-Kirchhoff stress and elasticity tensor
    void EvaluateStressCmatMembrane(LINALG::Matrix<3,3> const& CM,         ///< Right Cauchy Green tensor
                                    LINALG::Matrix<3,3> const& iFgM,       ///< Inelastic growth deformation gradient
                                    LINALG::Matrix<6,1> & stress,          ///< 2nd Piola-Kirchhoff stress
                                    LINALG::Matrix<6,6> & cmat,            ///< Elasticity tensor
                                    LINALG::Matrix<6,9> & dSdiFg,          ///< Derivative of 2nd Piola Kirchhoff stress w.r.t. the inelastic growth deformation gradient
                                    const int gp,                          ///< Current Gauß-Point
                                    const int eleGID) const;               ///< Element ID

    /// Convert matrix to "stress-like" Voigt notation
    template< typename T >
    inline void MatrixtoStressLikeVoigtNotation(LINALG::TMatrix<T,3,3> const& in,
                                                LINALG::TMatrix<T,6,1> & out) const {
      for(int i=0;i<3;++i)
        out(i) = in(i,i);
      out(3) = 0.5*(in(0,1)+in(1,0));
      out(4) = 0.5*(in(1,2)+in(2,1));
      out(5) = 0.5*(in(0,2)+in(2,0));
    };

    /// Convert matrix to "stress-like" Voigt notation
    inline void MatrixtoStrainLikeVoigtNotation(LINALG::Matrix<3,3> const& in,
                                                LINALG::Matrix<6,1> & out) const {
      for(int i=0;i<3;++i)
        out(i) = in(i,i);
      out(3) = in(0,1)+in(1,0);
      out(4) = in(1,2)+in(2,1);
      out(5) = in(0,2)+in(2,0);
    };

    /// Convert "stress-like" Voigt notation to matrix
    inline void StressLikeVoigtNotationtoMatrix(LINALG::Matrix<6,1> const& in,
                                                LINALG::Matrix<3,3> & out) const {
      for(int i=0;i<3;++i)
        out(i,i) = in(i);
      out(0,1) = out(1,0) = in(3);
      out(1,2) = out(2,1) = in(4);
      out(0,2) = out(2,0) = in(5);
    };

    /// convert non-symmetric 2-tensor to 9x1 vector
    inline void Matrix3x3to9x1(LINALG::Matrix<3,3> const& in,
                               LINALG::Matrix<9,1> & out) const {
      for(int i=0;i<3;++i) out(i) = in(i,i);
      out(3) = in(0,1);
      out(4) = in(1,2);
      out(5) = in(0,2);
      out(6) = in(1,0);
      out(7) = in(2,1);
      out(8) = in(2,0);
    };


    /// My material parameters
    MAT::PAR::GrowthRemodel_ElastHyper* params_;

    /// Map to remodelfiber material summands
    std::vector<Teuchos::RCP<MAT::ELASTIC::RemodelFiber> > potsumrf_;

    /// Map to elastin 3d matrix material summands
    std::vector<Teuchos::RCP<MAT::ELASTIC::Summand> > potsumeliso_;

    /// Map to membrane elastin material summands
    std::vector<Teuchos::RCP<MAT::ELASTIC::Summand> > potsumelmem_;

    /// Map to penalty elastin matrix material summands
    Teuchos::RCP<MAT::ELASTIC::Summand> potsumelpenalty_;

    /// Current individual mass density of elastin
    std::vector<double> cur_rho_el_;

    /// Initial individual mass density of elastin
    std::vector<double> init_rho_el_;

    /// Current volume change induced by growth
    std::vector<double> v_;

    /// Circumferential structural tensor
    LINALG::Matrix<3,3> AcirM_;

    /// Axial structural tensor
    LINALG::Matrix<3,3> AaxM_;

    /// Radial structural tensor
    LINALG::Matrix<3,3> AradM_;

    /// Radial structural tensor in "stress-like" Voigt notation
    LINALG::Matrix<6,1> Aradv_;

    /// Prestretch of elastin matrix in axial, circumferential and radial direction (used for prestressing)
    std::vector<LINALG::Matrix<3,3> > GM_;

    /// Total simulation time
    double t_tot_;

    /// Total number of remodel fibers
    unsigned nr_rf_tot_;

    ///                       (rad_x, axi_x, circ_x)
    /// Cylinder coordinates  (rad_y, axi_y, circ_y)
    ///                       (rad_z, axi_z, circ_z)
    LINALG::Matrix<3,3> radaxicirc_;

    /// Axial coordinate of each Gauß-Point
    std::vector<double> gp_ax_;

    /// Radial coordinate of each Gauß-Point
    std::vector<double> gp_rad_;

    /// Structural tensor of growth direction
    LINALG::Matrix<3,3> AgM_;

    /// Structural tensor of the plane in which all fibers are located
    LINALG::Matrix<3,3> AplM_;

    /// Fraction of 2D material parameter of elastin
    std::vector<double> mue_frac_;

    /// First call of Evaluate()
    std::vector<int> setup_;
  };

}

#endif

