/*----------------------------------------------------------------------*/
/*!
\file growthremodel_elasthyper.H

\brief This file is used to manage the homogenized constraint mixture during growth and remodeling

\level 3

<pre>
\maintainer Fabian Braeu
            braeu@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef GROWTHREMODEL_ELASTHYPER_H
#define GROWTHREMODEL_ELASTHYPER_H



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

      /// length of elastin matrix material list
      const int nummat_elastin_;

      /// the list of remodelfiber material IDs
      const std::vector<int>* matids_remodelfiber_;

      /// the list of elastin matrix material IDs
      const std::vector<int>* matids_elastin_;

      /// material ID of growth penalty material
      const int matid_penalty_;

      /// initial mass fraction of elastin in constraint mixture
      const std::vector<double>* init_w_el_;

      /// material mass density
      const double density_;

      /// time of prestressing
      const double t_prestress_;

      /// circumferential prestretch of elastin matrix
      const double lamb_prestretch_cir_;

      /// axial prestretch of elastin matrix
      const double lamb_prestretch_ax_;

      /// flag for turning elastin damage on or off: 1: damage on; 0: damage off
      const int damage_;

      /// flag to decide what type of collagen growth is used: 1: anisotropic growth; 0: isotropic growth
      const int growthtype_;

      /// flag to decide what type of local time integration scheme is used: 1: Backward Euler Method; 0: Forward Euler Method
      const int loctimeint_;

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
                          const int eleGID                    ///< Element ID
                          );

    /// setup
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    /// update
    virtual void Update(const LINALG::Matrix<3,3> defgrd,     ///< Deformation gradient
                        const int gp,                         ///< Current Gauß-Point
                        Teuchos::ParameterList& params,       ///< Container for additional information
                        const int eleGID);                    ///< Element ID

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string,int>& names);

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

  protected:

    /// Setup circumferential, radial and axial structural tensor
    void SetupAxiCirRadStructuralTensor(DRT::INPUT::LineDefinition* linedef);

    /// Read AXI CIR RAD direction
    void ReadDir(
        DRT::INPUT::LineDefinition* linedef,
        std::string specifier,
        LINALG::Matrix<3,1>& dir);

    virtual void EvaluateKinQuantElast(
        const LINALG::Matrix<3,3>* defgrd,
        const LINALG::Matrix<3,3> iFinel,
        const int gp,
        LINALG::Matrix<6,1>& iCinel,
        LINALG::Matrix<6,1>& iCinelCiCinel,
        LINALG::Matrix<6,1>& iC,
        LINALG::Matrix<3,3>& iCinelC,
        LINALG::Matrix<3,3>& iFinelCe,
        LINALG::Matrix<9,1>& CiFinel,
        LINALG::Matrix<9,1>& CiFinelCe,
        LINALG::Matrix<9,1>& CiFineliCe,
        LINALG::Matrix<3,1>& prinv);

    /// calculate principal invariants
    virtual void InvariantsPrincipal(
      LINALG::Matrix<3,1>& prinv,              ///< principal invariants
      const LINALG::Matrix<6,1>& C_strain      ///< symmetric Cartesian 2-tensor in strain-like 6-Voigt notation
      );

    /// calculate modified invariants
    virtual void InvariantsModified(
      LINALG::Matrix<3,1>& modinv,           ///< modified invariants
      const LINALG::Matrix<3,1>& prinv       ///< principal invariants
      );

    /// calculates the derivatives of the hyperelastic laws with respect to the (modified) invariants
    virtual void EvaluateInvariantDerivatives(
        const LINALG::Matrix<3,1>& prinv,               ///< principal invariants of right Cauchy-Green tensor
        const int& gp,                                  ///< current Gauß-point
        LINALG::Matrix<3,1>& dPIw,                      ///< first derivative with respect to invariants weighted with the corresponding volume fraction
        LINALG::Matrix<6,1>& ddPIIw,                    ///< second derivative with respect to invariants weighted with the corresponding volume fraction
        int eleGID                                      ///< Element ID
        );

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
    virtual void ConvertModToPrinc(
        const LINALG::Matrix<3,1>& prinv,         ///< principal invariants of right Cauchy-Green tensor
        const LINALG::Matrix<3,1>& dPmodI,        ///< first derivative with respect to modified invariants
        const LINALG::Matrix<6,1>& ddPmodII,      ///< second derivative with respect to modified invariants
        LINALG::Matrix<3,1>& dPI,                 ///< first derivative with respect to invariants
        LINALG::Matrix<6,1>& ddPII                ///< second derivative with respect to invariants
        );

    /// calculates the isotropic stress and elasticity tensor
    virtual void EvaluateIsotropicStressCmat(
        LINALG::Matrix<6,1>& stress,              ///< 2nd Piola Kirchhoff stress
        LINALG::Matrix<6,6>& cmat,                ///< elasticity tensor
        LINALG::Matrix<6,1> C_strain,             ///< right Cauchy green in strain-like 6-Voigt Voigt notation
        LINALG::Matrix<6,1> id,                   ///< Cartesian identity 2-tensor I_{AB}
        LINALG::Matrix<6,1> iC_stress,            ///< inverse right Cauchy-Green tensor stress-like 6-Voigt Voigt notation
        LINALG::Matrix<6,6> id4sharp,             ///< Cartesian identity 4-tensor in 6x6-matrix notation (stress-like)
        LINALG::Matrix<3,1> prinv,                ///< principal invariants of right Cauchy-Green tensor
        LINALG::Matrix<3,1> dPI,                  ///< first derivative with respect to invariants
        LINALG::Matrix<6,1> ddPII                 ///< second derivative with respect to invariants
        );

    /// calculates the isotropic stress and elasticity tensor for coupled configuration
    virtual void EvaluateIsotropicPrincElast(
        LINALG::Matrix<6,1>& stressisoprinc,
        LINALG::Matrix<6,6>& cmatisoprinc,
        LINALG::Matrix<6,1> iCg,
        LINALG::Matrix<6,1> iCgCiCg,
        LINALG::Matrix<6,1> iC,
        LINALG::Matrix<3,1> gamma,
        LINALG::Matrix<8,1> delta
        );

    /// Evaluate additional terms for for the elasticity tensor
    /// Additional terms are caused by implicit implementation of the evolution equation
    virtual void EvaluateAdditionalTermsCmat(
        LINALG::Matrix<6,6>& cmatadd,                   ///< additional elasticity tensor
        LINALG::Matrix<3,1> gamma,                      ///< factors for stress calculation
        LINALG::Matrix<8,1> delta,                      ///< factors for elasticity tensor calculation
        double v,                                       ///< current volume change
        int gp,                                         ///< current gp
        LINALG::Matrix<3,3> AgM,                        ///< structural tensor of growth direction
        LINALG::Matrix<1,6> sum_drhoC,                  ///< sum over all fiber families of the derivatives of the individual mass density w.r.t. right Cauchy Green tensor weighted with the corresponding initial mass fraction,k_sig,dt and 1/sig_pre
        LINALG::Matrix<3,3> iFg,
        LINALG::Matrix<3,3> iCgC,
        LINALG::Matrix<3,3> iFgCe,
        LINALG::Matrix<6,1> iCg,
        LINALG::Matrix<9,1> CiFg,
        LINALG::Matrix<9,1> CiFgCe,
        LINALG::Matrix<6,1> iCgCiCg,
        LINALG::Matrix<9,1> CiFgiCe,
        LINALG::Matrix<6,1> iC
    );

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
    virtual void CalculateGammaDelta(
        LINALG::Matrix<3,1>& gamma,             ///< factors for stress calculation
        LINALG::Matrix<8,1>& delta,             ///< factors for elasticity tensor calculation
        const LINALG::Matrix<3,1>& prinv,       ///< principal invariants of right Cauchy-Green tensor
        const LINALG::Matrix<3,1>& dPI,         ///< first derivative with respect to invariants
        const LINALG::Matrix<6,1>& ddPII        ///< second derivative with respect to invariants
        );

    /// Evaluate current individual volume of elastin
    virtual void EvaluateElastinDamage();

    /// add a contribution to the derivative of the second Piola Kirchhoff stress with respect to the growth deformation gradient
    /// in index notation this method does:
    /// dPK2dFginv_IJKL += fac*( A_IK * B_JL + A_JK * B_IL )
    /// Note that the result is symmetric within the first two indices
    /// Thus the first two indices are stored in stress-like Voigt notation.
    /// The 2nd index pair has no symmetries and is therefore stored as 9-vector
    /// Be aware that this is NOT symmetric w.r.t. A and B !!!
    void AddtodPK2diFinel(LINALG::Matrix<6,9>& dPK2diFinel,
        LINALG::Matrix<3,3> A,
        LINALG::Matrix<3,3> B,
        double fac);

    /// This function adds the following contribution to the matrix X_ABCD (9x9)
    /// in index notation this method does:
    /// X_ABCD = fac*( A_AC * B_BD )
    /// There are no symmetries and the matrix is therefore stored in a 9x9 matrix
    /// A and B are not required to be symmetric!!!
    virtual void AddtoMatrixCrossProduct(const double& fac,                    ///< scalar factor
                                         const LINALG::Matrix<3,3>& A,         ///< Matrix in the above mentioned equation
                                         const LINALG::Matrix<3,3>& B,         ///< Matrix in the above mentioned equation
                                         LINALG::Matrix<9,9>& X_ABCD);         ///< Output (see EvaluateInelasticDefGrad)

    /// build growth deformation gradient
    void EvaluateGrowthDefGrad(LINALG::Matrix<3,3>& FgM,        ///< Growth deformation gradient
                               LINALG::Matrix<3,3>& iFgM,       ///< Inverse growth deformation gradient
                               const int gp);                   ///< Current Gauß-Point

    /// convert non-symmetric 2-tensor to 9x1 vector
    void Matrix3x3to9x1(LINALG::Matrix<3,3> A, LINALG::Matrix<9,1>& Out);


    /// evaluates stress and cmat during prestressing
    void EvaluateStressCmat(const LINALG::Matrix<3,3>* defgrd,     ///< Deformation gradient
                            const LINALG::Matrix<3,3> Finel,       ///< Inelastic deformation gradient
                            const int gp,                          ///< Current Gauß-Point
                            LINALG::Matrix<6,1>& stressiso,        ///< Isotropic stress tensor
                            LINALG::Matrix<6,6>& cmatiso,          ///< Isotropic stiffness matrix
                            const int eleGID                       ///< Element ID
                            );

    /// my material parameters
    MAT::PAR::GrowthRemodel_ElastHyper* params_;

    /// map to remodelfiber material summands
    std::vector<Teuchos::RCP<MAT::ELASTIC::RemodelFiber> > potsumrf_;

    /// map to elastin matrix material summands
    std::vector<Teuchos::RCP<MAT::ELASTIC::Summand> > potsumel_;

    /// map to ground matrix material summands
    Teuchos::RCP<MAT::ELASTIC::Summand> potsumpenalty_;

    /// current individual mass density of elastin
    std::vector<std::vector<double> > cur_rho_el_;

    /// initial individual mass density of elastin
    std::vector<double> init_rho_el_;

    /// current volume change induced by growth
    std::vector<double> v_;

    /// circumferential structural tensor
    LINALG::Matrix<3,3> AcirM_;

    /// axial structural tensor
    LINALG::Matrix<3,3> AaxM_;

    /// radial structural tensor
    LINALG::Matrix<3,3> AradM_;

    /// radial structural tensor in "stress-like" Voigt notation
    LINALG::Matrix<6,1> Aradv_;

    /// prestretch of elastin matrix in axial, circumferential and radial direction (used for prestressing)
    std::vector<LINALG::Matrix<3,3> > GM_;

    /// total simulation time
    double t_tot_;

    /// cylinder coordinates (axi, circ, rad)
    std::vector<LINALG::Matrix<1,3> > cylcoords_;

    /// axial coordinate of each Gauß-Point
    std::vector<double> gp_ax_;

    /// radial coordinate of each Gauß-Point
    std::vector<double> gp_rad_;

    /// structural tensor of growth direction
    LINALG::Matrix<3,3> AgM_;

    /// structural tensor of the plane in which all fibers are located
    LINALG::Matrix<3,3> AplM_;


    /// wall thickness in reference configuration
    double t_ref_;

    /// mean blood pressure
    double p_mean_;

  };

}

#endif

