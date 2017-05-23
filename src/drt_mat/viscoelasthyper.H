/*----------------------------------------------------------------------*/
/*!
\file viscoelasthyper.H

\brief
This file contains the viscohyperelastic material toolbox.
The availabele viscous models can be applied to any hyperelastic law
of the Elasthyper toolbox.
The viscous part is summed up with the hyperelastic laws
to build a viscohyperelastic material model.

The input line should read
MAT 0   MAT_ViscoElastHyper   NUMMAT 2 MATIDS 1 2 DENS 0

\level 1

<pre>
\maintainer Anna Birzle
            birzle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef VISCOELASTHYPER_H
#define VISCOELASTHYPER_H



#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"

#include "elasthyper.H"

// forward declaration due to avoid header definition
namespace MAT
{
  namespace ELASTIC
  {
    class Summand;
  }

  // forward declaration
  class ViscoElastHyper;

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// Collection of viscohyperelastic materials
    ///
    /// Storage map of hyperelastic summands.


  class ViscoElastHyper : public MAT::PAR::ElastHyper
//    class ViscoElastHyper : public Parameter
    {

      friend class MAT::ViscoElastHyper;

    public:

      /// standard constructor
      ///
      /// This constructor recursively calls the constructors of the
      /// parameter sets of the hyperelastic summands.
      ViscoElastHyper(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~ViscoElastHyper() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      // elasthyper needs this for stat inverse analysis
      virtual void OptParams(std::map<std::string, int>* pnames){dserror("This material is not implemented in OptParams.");};

      //@}

    };  // class ViscoElastHyper

  }  // namespace PAR

class ViscoElastHyperType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "ViscoElastHyperType"; }

  static ViscoElastHyperType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static ViscoElastHyperType instance_;
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
  class Material;

  class ViscoElastHyper : public MAT::ElastHyper
  {
  public:

    /// construct empty material object
    ViscoElastHyper();

    /// construct the material object given material parameters
    explicit ViscoElastHyper(MAT::PAR::ViscoElastHyper* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const { return ViscoElastHyperType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_viscoelasthyper; }

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem)
    { if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
      dserror("element and material kinematics are not compatible");}

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new ViscoElastHyper( *this ) ); }

    /// Check if history variables are already initialized
    virtual bool Initialized() const
    {
      return isinitvis_&&(histscgcurr_!=Teuchos::null);
      return isinitvis_&&(histstresscurr_!=Teuchos::null);
      return isinitvis_&&(histartstresscurr_!=Teuchos::null);
    }

    /// hyperelastic stress response plus elasticity tensor
    virtual void Evaluate(const LINALG::Matrix<3,3>* defgrd,  ///< Deformation gradient
                          const LINALG::Matrix<6,1>* glstrain,///< Green-Lagrange strain
                          Teuchos::ParameterList& params,     ///< Container for additional information
                          LINALG::Matrix<6,1>* stress,        ///< 2nd Piola-Kirchhoff stresses
                          LINALG::Matrix<6,6>* cmat,
                          const int eleGID);         ///< Constitutive matrix

    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    virtual void ResetAll(const int numgp);
    virtual void SetupTimeVariables(){;};

    virtual void Update();

  protected:
    /// calculates the kinematic quantities and tensors used afterwards for viscous part
    virtual void EvaluateKinQuantVis(
        LINALG::Matrix<6,1> rcg,
        LINALG::Matrix<6,1> scg,
        LINALG::Matrix<6,1>& modrcg,
        LINALG::Matrix<6,1> icg,
        Teuchos::ParameterList& params,
        LINALG::Matrix<3,1> prinv,
        LINALG::Matrix<6,1>& scgrate,
        LINALG::Matrix<6,1>& modrcgrate,
        LINALG::Matrix<7,1>& rateinv,
        LINALG::Matrix<7,1>& modrateinv
      );

    /// calculates the factors associated to the viscous laws
    virtual void EvaluateMuXi(
        LINALG::Matrix<3,1> inv,
        LINALG::Matrix<3,1> modinv,
        LINALG::Matrix<8,1>& mu,
        LINALG::Matrix<8,1>& modmu,
        LINALG::Matrix<33,1>& xi,
        LINALG::Matrix<33,1>& modxi,
        LINALG::Matrix<7,1>& rateinv,
        LINALG::Matrix<7,1>& modrateinv,
        Teuchos::ParameterList& params,
        const int eleGID
        );

    /// calculates the isotropic stress and elasticity tensor for viscous principal configuration
    virtual void EvaluateIsoViscoPrincipal(
        LINALG::Matrix<6,1>& stress,
        LINALG::Matrix<6,6>& cmat,
        LINALG::Matrix<8,1> mu,
        LINALG::Matrix<33,1> xi,
        LINALG::Matrix<6,6> id4sharp,
        LINALG::Matrix<6,1> scgrate
        );

    /// calculates the isotropic stress and elasticity tensor for viscous decoupled configuration
    virtual void EvaluateIsoViscoModified(
        LINALG::Matrix<6,1>& stressisomodisovisco,
        LINALG::Matrix<6,1>& stressisomodvolvisco,
        LINALG::Matrix<6,6>& cmatisomodisovisco,
        LINALG::Matrix<6,6>& cmatisomodvolvisco,
        LINALG::Matrix<3,1> prinv,
        LINALG::Matrix<3,1> modinv,
        LINALG::Matrix<8,1> modmu,
        LINALG::Matrix<33,1> modxi,
        LINALG::Matrix<6,1> rcg,
        LINALG::Matrix<6,1> id2,
        LINALG::Matrix<6,1> icg,
        LINALG::Matrix<6,6> id4,
        LINALG::Matrix<6,1> modrcgrate
        );

    /// calculates the stress and elasticitiy tensor for the viscous Genmax-material
    /// depending on the viscoelastic material isochoric-principal, isochoric-modified
    /// (volumetric and isochoric equal treaded) or anisotropic stresses and elasticity
    /// tensors are added
    virtual void EvaluateViscoGenMax(
        LINALG::Matrix<6,1> stress,
        LINALG::Matrix<6,6> cmat,
        LINALG::Matrix<6,1>& Q,
        LINALG::Matrix<6,6>& cmatq,
        double tau,
        double beta,
        Teuchos::ParameterList& params
        );

    /// @name Flags to specify the viscous formulations
    //@{
    bool isovisco_ ; ///< global indicator for isotropic splitted viscous formulation
    bool viscogenmax_ ; ///< global indicator for viscous contribution according the SLS-Model
    //@}

  private:

    /// viscous history Cauchy-Green-Tensor
    Teuchos::RCP<std::vector<LINALG::Matrix<6,1> > > histscgcurr_; ///< current Cauchy-Green-Tensor
    Teuchos::RCP<std::vector<LINALG::Matrix<6,1> > > histscglast_; ///< Cauchy-Green-Tensor of last converged state
    Teuchos::RCP<std::vector<LINALG::Matrix<6,1> > > histmodrcgcurr_; ///< current decoupled Cauchy-Green-Tensor
    Teuchos::RCP<std::vector<LINALG::Matrix<6,1> > > histmodrcglast_; ///< decoupled Cauchy-Green-Tensor of last converged state
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histstresscurr_; ///< current stress
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histstresslast_; ///< stress of last converged state
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histartstresscurr_; ///< current stress
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histartstresslast_; ///< stress of last converged state

    bool isinitvis_;  ///< indicates if #Initialized routine has been called
  };

}

#endif

