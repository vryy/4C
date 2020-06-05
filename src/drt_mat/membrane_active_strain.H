/*----------------------------------------------------------------------*/
/*! \file
 \brief Active strain membrane material for gastric electromechanics

 The input line should read
 MAT 0 MAT_Membrane_ActiveStrain

 \level 3

 \maintainer Sebastian Brandstaeter


 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                     brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
#ifndef MEMBRANE_ACTIVE_STRAIN_H
#define MEMBRANE_ACTIVE_STRAIN_H

/*----------------------------------------------------------------------*
 | headers                                         brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"
/*----------------------------------------------------------------------*
 | active strain membrane material                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
namespace MAT
{
  // forward declaration
  class Membrane_ActiveStrain;

  namespace PAR
  {
    class Membrane_ActiveStrain : public Parameter
    {
      friend class MAT::Membrane_ActiveStrain;

     public:
      /// standard constructor
      Membrane_ActiveStrain(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~Membrane_ActiveStrain() { ; }

      /// @name material parameters
      //@{

      /// Number of the material that describes the elastic behavior
      const int matid_passive_;

      /// Position of the transmembrane voltage of the SMC cell
      const int scalid_voltage_;

      /// density
      const double density_;

      /// beta (parameter dynamics of the VDCC)
      const double beta1_;

      /// beta (parameter dynamics of the Ca-ions)
      const double beta2_;

      /// voltage level for the activation of the Ca-influx
      const double voltage_threshold_;

      /// alpha (parameter for the intensity of the contraction in fiber 1 direction)
      const double alpha1_;

      /// alpha (parameter for the intensity of the contraction in fiber 2 direction)
      const double alpha2_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();
    };
    // class Membrane_ActiveStrain

  }  // namespace PAR

  class Membrane_ActiveStrainType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "Membrane_ActiveStrainType"; }

    static Membrane_ActiveStrainType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static Membrane_ActiveStrainType instance_;
  };
  // class Membrane_ActiveStrainType

  /*----------------------------------------------------------------------*/
  /// \author brandstaeter
  /// \date 05/18
  // forward declaration
  class Material;

  class Membrane_ActiveStrain : public So3Material
  {
   public:
    /// construct empty material object
    Membrane_ActiveStrain();

    /// construct the material object given material parameters
    explicit Membrane_ActiveStrain(MAT::PAR::Membrane_ActiveStrain* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const
    {
      return Membrane_ActiveStrainType::Instance().UniqueParObjectId();
    }

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
    {
      if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const
    {
      return INPAR::MAT::m_membrane_activestrain;
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new Membrane_ActiveStrain(*this));
    }

    /// material mass density
    virtual double Density() const { return params_->density_; }

    /// setup
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    /// Standard SO3 evaluate (not meant to be used)
    virtual void Evaluate(const LINALG::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const LINALG::Matrix<6, 1>* glstrain,                  ///< Green-Lagrange strain
        Teuchos::ParameterList& params,  ///< Container for additional information
        LINALG::Matrix<6, 1>* stress,    ///< 2nd Piola-Kirchhoff stresses
        LINALG::Matrix<6, 6>* cmat,      ///< Constitutive matrix
        int gp,                          ///< Gauss point
        int eleGID)                      ///< Element ID
    {
      dserror("This a membrane material. Calling So3 evaluate does not make sense.");
    };

    /// hyperelastic stress response plus elasticity tensor for membranes assuming incompressibility
    /// and plane stress
    virtual void Evaluate(LINALG::Matrix<3, 3>& cauchygreen,  ///< right Cauchy-Green tensor
        Teuchos::ParameterList& params,  ///< Container for additional information
        LINALG::Matrix<3, 3>&
            Q_trafo,  ///< Trafo from local membrane orthonormal coordinates to global coordinates
        LINALG::Matrix<3, 1>*
            stress,                  ///< 2nd Piola-Kirchhoff stresses in stress-like voigt notation
        LINALG::Matrix<3, 3>* cmat,  ///< Constitutive matrix
        int gp,                      ///< Gauss point
        int eleGID                   ///< Element GID
    );

    /// Update internal variables
    virtual void Update();

    /// Reset internal variables
    virtual void ResetStep();

    /// Reset internal variables to state in the beginning of the computation (needed for inverse
    /// analysis)
    virtual void ResetAll(int numgp);

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string, int>& names);

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

   private:
    /// My material parameters
    MAT::PAR::Membrane_ActiveStrain* params_;

    /// passive material
    Teuchos::RCP<MAT::So3Material> matpassive_;

    /// (tansmembrane) voltage at every gp
    Teuchos::RCP<std::vector<double>> voltage_;

    /// activation parameter at every gp
    Teuchos::RCP<std::vector<double>> activation_;

    /// indicates if material is initialized
    bool isinit_;

    // setup fiber vectors
    void SetupFiberVectors(int numgp, DRT::INPUT::LineDefinition* linedef);

    // read RAD-AXI-CIR
    void ReadDir(
        DRT::INPUT::LineDefinition* linedef, std::string specifier, LINALG::Matrix<3, 1>& dir);

    // calculate normal direction from FIBER1 and FIBER2
    void SetupNormalDirection();

    // pullback of the tangent from intermediate to reference configuration
    void Pullback4thTensorVoigt(const LINALG::Matrix<2, 2>& defgrd_active_inv_red,
        const LINALG::Matrix<3, 3>& cmat_passive_intermediate,
        LINALG::Matrix<3, 3>* cmat_reference);

    // transform voigt to tensor notation
    void Tensor2x2Indices(int p, int* i, int* j);

    // transform tensor to voigt notation
    void Voigt3Index(int i, int j, int* p);

   protected:
    /// vector of fiber vectors
    std::vector<LINALG::Matrix<3, 1>> fibervecs_;
  };  // class Membrane_ActiveStrain

}  // namespace MAT
// namespace MAT
#endif /* MEMBRANE_ACTIVE_STRAIN_H */
