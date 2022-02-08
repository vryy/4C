/*----------------------------------------------------------------------*/
/*! \file
\brief hyperelastic toolbox for membranes assuming incompressibility and plane stress

The input line should read
MAT 0 MAT_Membrane_ElastHyper NUMMAT 2 MATIDS 1 2 DENS 0

\level 3


*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                           sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
#ifndef MEMBRANE_ELASTHYPER_H
#define MEMBRANE_ELASTHYPER_H

/*----------------------------------------------------------------------*
 | headers                                               sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
#include "elasthyper.H"
#include "membrane_material_interfaces.H"

/*----------------------------------------------------------------------*
 | hyperelastic material for membranes                   sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
namespace MAT
{
  // forward declaration
  namespace ELASTIC
  {
    class Summand;
  }  // namespace ELASTIC

  class Membrane_ElastHyper;

  namespace PAR
  {
    class Membrane_ElastHyper : public MAT::PAR::ElastHyper
    {
      friend class MAT::Membrane_ElastHyper;

     public:
      /// standard constructor
      Membrane_ElastHyper(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~Membrane_ElastHyper() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class Membrane_ElastHyper

  }  // namespace PAR

  class Membrane_ElastHyperType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "Membrane_ElastHyperType"; }

    static Membrane_ElastHyperType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static Membrane_ElastHyperType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Collection of hyperelastic materials
  ///
  /// This collection offers to additively compose a stress response
  /// based on summands defined separately.  This is possible, because
  /// we deal with hyperelastic materials, which are composed
  /// of (Helmholtz free energy density) potentials.  Effectively, we want
  ///\f[
  ///  \Psi(\boldsymbol{C}) = \sum_i \Psi_i(\boldsymbol{C})
  ///\f]
  /// in which the individual \f$\Psi_i\f$ is implemented as #MAT::ELASTIC::Summand.
  ///
  /// Quite often the right Cauchy-Green 2-tensor \f$\boldsymbol{C}\f$
  /// is replaced by its various invariant forms as argument.
  ///
  /// The task of ElastHyper is the evaluation of the
  /// potential energies and their derivatives to obtain the actual
  /// stress response and the elasticity tensor. The storage is located
  /// at the associated member #params_.
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
  /// </ul>
  ///
  /// \author rausch,tk,bborn
  /// \date 05/09

  // forward declaration
  class Material;

  class Membrane_ElastHyper : public ElastHyper, public MAT::MembraneMaterialLocalCoordinates
  {
   public:
    /// construct empty material object
    Membrane_ElastHyper();

    /// construct the material object given material parameters
    explicit Membrane_ElastHyper(MAT::PAR::Membrane_ElastHyper* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const
    {
      return Membrane_ElastHyperType::Instance().UniqueParObjectId();
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

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const
    {
      return INPAR::MAT::m_membrane_elasthyper;
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new Membrane_ElastHyper(*this));
    }

    /// setup
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    void UpdateMembrane(const LINALG::Matrix<3, 3>& defgrd, Teuchos::ParameterList& params,
        const LINALG::Matrix<3, 3>& Q_trafo, int gp, int eleGID) override
    {
      // nothing to do
    }

    void EvaluateMembrane(const LINALG::Matrix<3, 3>& defgrd,
        const LINALG::Matrix<3, 3>& cauchygreen, Teuchos::ParameterList& params,
        const LINALG::Matrix<3, 3>& Q_trafo, LINALG::Matrix<3, 1>& stress,
        LINALG::Matrix<3, 3>& cmat, int gp, int eleGID) override;

    /// evaluate strain energy function
    virtual void StrainEnergy(LINALG::Matrix<3, 3>& cauchygreen,  ///< right Cauchy-Green tensor
        double& psi,                                              ///< Strain energy function
        int gp,                                                   ///< Gauss point
        int eleGID                                                ///< Element GID
    );

   protected:
    /// calculate anisotropic stress and elasticity tensor
    virtual void EvaluateAnisotropicStressCmat(LINALG::Matrix<3, 1>& stress_aniso,
        LINALG::Matrix<3, 3>& cmat_aniso, const LINALG::Matrix<3, 3>& Q_trafo,
        const LINALG::Matrix<3, 1>& rcg, const double& rcg33, Teuchos::ParameterList& params,
        int gp, int eleGID);

    /// vector of fiber vectors
    std::vector<LINALG::Matrix<3, 1>> fibervecs_;
  };

}  // namespace MAT

/*----------------------------------------------------------------------*/
#endif  // MEMBRANE_ELASTHYPER_H
