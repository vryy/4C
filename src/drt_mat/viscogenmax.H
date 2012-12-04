/*----------------------------------------------------------------------*/
/*!
\file viscogenmax.H

<pre>
Maintainer: Aline Brunon
            bel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef VISCOGENMAX_H
#define VISCOGENMAX_H



#include "matpar_parameter.H"
#include "material.H"
#include "elasthyper.H"

using namespace Teuchos;

// forward declaration due to avoid header definition
namespace MAT
{
  namespace ELASTIC
  {
    class Summand;
  }

  // forward declaration
  class ViscoGenMax;

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// Collection of hyperelastic materials
    ///
    /// Storage map of hyperelastic summands.
    ///
    /// \author rausch,tk,bborn
    /// \date 05/09

    class ViscoGenMax : public ElastHyper
    {

      friend class MAT::ViscoGenMax;

    public:

      /// standard constructor
      ///
      /// This constructor recursively calls the constructors of the
      /// parameter sets of the hyperelastic summands.
      ViscoGenMax(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~ViscoGenMax() { ; }

      /// @name material parameters
      //@{

      // viscoelastic constant
      // isotropic principal configuration
      double relax_isot_princ_; // relaxation time
      double beta_isot_princ_; // ratio between the two branches of the gen. maxwell model
      // isotropic modified configuration
      double relax_isot_mod_vol_; // relaxation time
      double beta_isot_mod_vol_; // ratio between the two branches of the gen. maxwell model
      double relax_isot_mod_isoc_; // relaxation time
      double beta_isot_mod_isoc_; // ratio between the two branches of the gen. maxwell model
      // anisotropic principal configuration
      double relax_anisot_princ_; // relaxation time
      double beta_anisot_princ_; // ratio between the two branches of the gen. maxwell model
      // anisotropic modified configuration
      double relax_anisot_mod_vol_; // relaxation time
      double beta_anisot_mod_vol_; // ratio between the two branches of the gen. maxwell model
      double relax_anisot_mod_isoc_; // relaxation time
      double beta_anisot_mod_isoc_; // ratio between the two branches of the gen. maxwell model

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// @name Write methods
      //@{
      void SetTauIsoPri(const double relax_isot_princ) { relax_isot_princ_ = relax_isot_princ; }
      void SetTauIsoModVol(const double relax_isot_mod_vol) { relax_isot_mod_vol_ = relax_isot_mod_vol; }
      void SetTauIsoModIsoc(const double relax_isot_mod_isoc) { relax_isot_mod_isoc_ = relax_isot_mod_isoc; }
      void SetTauAnisoPri(const double relax_anisot_princ) { relax_anisot_princ_ = relax_anisot_princ; }
      void SetTauAnisoModVol(const double relax_anisot_mod_vol) { relax_anisot_mod_vol_ = relax_anisot_mod_vol; }
      void SetTauAnisoModIsoc(const double relax_anisot_mod_isoc) { relax_anisot_mod_isoc_ = relax_anisot_mod_isoc; }

      //@}

//      void SetTauIsoPri(const double relax_isot_princ) { tau_isoprinc_ = relax_isot_princ; }
//      void SetTauIsoModVol(const double relax_isot_mod_vol) { tau_isomod_vol_ = relax_isot_mod_vol; }
//      void SetTauIsoModIsoc(const double relax_isot_mod_isoc) { tau_isomod_iso_ = relax_isot_mod_isoc; }
//      void SetTauAnisoPri(const double relax_anisot_princ) { tau_anisoprinc_ = relax_anisot_princ; }
//      void SetTauAnisoModVol(const double relax_anisot_mod_vol) { tau_anisomod_vol_ = relax_anisot_mod_vol; }
//      void SetTauAnisoModIsoc(const double relax_anisot_mod_isoc) { tau_anisomod_iso_ = relax_anisot_mod_isoc; }


   private:

      /// map to materials/potential summands
      std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> > potsum_;

      /// variable including the informations about the hyperelastic part
      std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> > elasthyper_;

      //@}

    };  // class ViscoGenMax

  }  // namespace PAR

class ViscoGenMaxType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "ViscoGenMaxType"; }

  static ViscoGenMaxType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static ViscoGenMaxType instance_;
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

  class ViscoGenMax : public ElastHyper
  {
  public:

    /// construct empty material object
    ViscoGenMax();

    /// construct the material object given material parameters
    explicit ViscoGenMax(MAT::PAR::ViscoGenMax* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const { return ViscoGenMaxType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_viscogenmax; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new ViscoGenMax( *this ) ); }

    /// Check if history variables are already initialized
    virtual bool Initialized() const
    { return isinit_&&(histstressisoprinccurr_!=Teuchos::null);
      return isinit_&&(artstressisoprinccurr_!=Teuchos::null);
      return isinit_&&(histstressisomodisocurr_!=Teuchos::null);
      return isinit_&&(artstressisomodisocurr_!=Teuchos::null);
      return isinit_&&(histstressisomodvolcurr_!=Teuchos::null);
      return isinit_&&(artstressisomodvolcurr_!=Teuchos::null);
      return isinit_&&(histstressanisoprinccurr_!=Teuchos::null);
      return isinit_&&(artstressanisoprinccurr_!=Teuchos::null);
    }

    /// hyperelastic stress response plus elasticity tensor
    virtual void Evaluate(
      const LINALG::Matrix<6,1>& glstrain,  ///< Green-Lagrange strain
      LINALG::Matrix<6,6>& cmat,  ///< elasticity matrix
      LINALG::Matrix<6,1>& stress,  ///< 2nd Piola--Kirchhoff stress
      const int gp,   ///< current Gauss point
      Teuchos::ParameterList& params  ///< parameter list for communication
      );


    virtual void Setup(const int numgp,DRT::INPUT::LineDefinition* linedef);
    /// setup patient-specific AAA stuff
//    void SetupAAA(Teuchos::ParameterList& params);
    //void Setupvisco(const int numgp);

    virtual void Setupvisco(const int numgp);
    virtual void SetupTimeVariables();

    // update internal stress variables
    virtual void Update();

    // reset internal stress variables
    virtual void Reset();

  private:

    /// my material parameters
    MAT::PAR::ViscoGenMax* viscoparams_;

    /// visco history stresses
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histstressisoprinccurr_; ///< current stress
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histstressisoprinclast_; ///< stress of last converged state
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > artstressisoprinccurr_;  ///< current artificial stress
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > artstressisoprinclast_;  ///< artificial stress in last converged state
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histstressisomodisocurr_; ///< current stress
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histstressisomodisolast_; ///< stress of last converged state
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > artstressisomodisocurr_;  ///< current artificial stress
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > artstressisomodisolast_;  ///< artificial stress in last converged state
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histstressisomodvolcurr_; ///< current stress
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histstressisomodvollast_; ///< stress of last converged state
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > artstressisomodvolcurr_;  ///< current artificial stress
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > artstressisomodvollast_;  ///< artificial stress in last converged state
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histstressanisoprinccurr_; ///< current stress
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histstressanisoprinclast_; ///< stress of last converged state
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > artstressanisoprinccurr_;  ///< current artificial stress
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > > artstressanisoprinclast_;  ///< artificial stress in last converged state

    bool isinit_;  ///< indicates if #Initialized routine has been called

    /// flags for the viscous contributions
    bool vis_isoprinc_; ///< indicates that a viscous contribution is added to the isotropic not splitted contribution
    bool vis_isomod_vol_; ///< indicates that a viscous contribution is added to the isotropic splitted contribution on the volumetric part
    bool vis_isomod_iso_; ///< indicates that a viscous contribution is added to the isotropic splitted contribution on the isochoric part
    bool vis_anisoprinc_; ///< indicates that a viscous contribution is added to the anisotropic not splitted contribution
    bool vis_anisomod_vol_; ///< indicates that a viscous contribution is added to the anisotropic splitted contribution on the volumetric part
    bool vis_anisomod_iso_;///< indicates that a viscous contribution is added to the anisotropic splitted contribution on the isochoric part

//    /// viscous constants of the different formulations
//    double tau_isoprinc_;
//    double beta_isoprinc_;
//    double tau_isomod_vol_;
//    double beta_isomod_vol_;
//    double tau_isomod_iso_;
//    double beta_isomod_iso_;
//    double tau_anisoprinc_;
//    double beta_anisoprinc_;
//    double tau_anisomod_vol_;
//    double beta_anisomod_vol_;
//    double tau_anisomod_iso_ ;
//    double beta_anisomod_iso_;

//    void SetTauIsoPri(const double relax_isot_princ) { tau_isoprinc_ = relax_isot_princ; }
//    void SetTauIsoModVol(const double relax_isot_mod_vol) { tau_isomod_vol_ = relax_isot_mod_vol; }
//    void SetTauIsoModIsoc(const double relax_isot_mod_isoc) { tau_isomod_iso_ = relax_isot_mod_isoc; }
//    void SetTauAnisoPri(const double relax_anisot_princ) { tau_anisoprinc_ = relax_anisot_princ; }
//    void SetTauAnisoModVol(const double relax_anisot_mod_vol) { tau_anisomod_vol_ = relax_anisot_mod_vol; }
//    void SetTauAnisoModIsoc(const double relax_anisot_mod_isoc) { tau_anisomod_iso_ = relax_anisot_mod_isoc; }
    };

}

#endif

