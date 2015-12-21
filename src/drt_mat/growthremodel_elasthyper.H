/*----------------------------------------------------------------------*/
/*!
\file growthremodel_elasthyper.H

<pre>
Maintainer: Fabian Braeu
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
    class VolGrowthPenalty;
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

      /// length of ground matrix material list
      const int nummat_ground_;

      /// the list of remodelfiber material IDs
      const std::vector<int>* matids_remodelfiber_;

      /// the list of elastin matrix material IDs
      const std::vector<int>* matids_elastin_;

      /// the list of ground matrix material IDs
      const std::vector<int>* matids_ground_;

      /// material ID of growth penalty material
      const int matid_penalty_;

      /// current mass fraction of elastin in constraint mixture
      const std::vector<double>* cur_w_elastin_;

      /// current mass fraction of ground matrix in constraint mixture
      const std::vector<double>* cur_w_ground_;

      /// material mass density
      const double density_;

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      //@}

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growthremodel_elasthyper; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new GrowthRemodel_ElastHyper( *this ) ); }

    /// material mass density
    virtual double Density() const { return params_->density_; }

    /// hyperelastic stress response plus elasticity tensor
    virtual void Evaluate(const LINALG::Matrix<3,3>* defgrd,  ///< Deformation gradient
                          const LINALG::Matrix<6,1>* glstrain,///< Green-Lagrange strain
                          Teuchos::ParameterList& params,     ///< Container for additional information
                          LINALG::Matrix<6,1>* stress,        ///< 2nd Piola-Kirchhoff stresses
                          LINALG::Matrix<6,6>* cmat,         ///< Constitutive matrix
                          const int eleGID                    ///< Element GID
                          );

    /// setup
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    /// update
    virtual void Update();

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string,int>& names);

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

  protected:

    /// @name Voigt index arrays
    //@{
    static const int VOIGT3X3_[3][3];
    static const int VOIGT3X3NONSYM_[3][3];

    //@}

    /// my material parameters
    MAT::PAR::GrowthRemodel_ElastHyper* params_;

    /// map to remodelfiber material summands
    std::vector<Teuchos::RCP<MAT::ELASTIC::RemodelFiber> > potsumrf_;

    /// map to elastin matrix material summands
    std::vector<Teuchos::RCP<MAT::ELASTIC::Summand> > potsumel_;

    /// map to ground matrix material summands
    std::vector<Teuchos::RCP<MAT::ELASTIC::Summand> > potsumgr_;

    /// map to ground matrix material summands
    Teuchos::RCP<MAT::ELASTIC::VolGrowthPenalty> potsumpenalty_;

    /// current mass fractions of elastin in constraint mixture
    std::vector<double> current_w_elastin_;

    /// current mass fractions of ground matrix in constraint mixture
    std::vector<double> current_w_ground_;

  };

}

#endif

