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

#ifdef CCADISCRET

#include "matpar_parameter.H"
#include "material.H"
#include "elasthyper.H"

// included in matpar_parameter.H -> matpar_material.H
//#include "../drt_lib/drt_parobjectfactory.H"
// avoided by forward declaration
//#include "../drt_matelast/elast_summand.H"
// included in material.h -> inpar_material.H -> inpar_parameterlist.H
//#include <Teuchos_ParameterList.hpp>

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

      /// provide access to material/summand by its ID
      Teuchos::RCP<const MAT::ELASTIC::Summand> MaterialById(
        const int id  ///< ID to look for in collection of summands
        ) const;

      /// length of material list
//      const int nummat_;

      /// the list of material IDs
//      const std::vector<int>* matids_;

      /// material mass density
//      const double density_;

      /// fiber angle
//      double gamma_;

      // viscoelastic constant
      const double relax_; // relaxation time
      const double beta_; // ratio between the two branches of the gen. maxwell model

      /// initalization modus;
      /// default = -1, Initmodus derived from locsys = 0, Initmodus aligned to locsys = 1
//      const int init_mode_;

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

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

  virtual DRT::ParObject* Create( const vector<char> & data );

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
    virtual void Unpack(const vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_viscogenmax; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new ViscoGenMax( *this ) ); }

    /// number of materials
//    int NumMat() const { return params_->nummat_; }

    /// deliver material ID of index i'th potential summand in collection
//    int MatID(
//      const unsigned index  ///< index i
//      ) const;

    /// provide access to material by its ID
//    Teuchos::RCP<const MAT::ELASTIC::Summand> MaterialById(const int id) const { return params_->MaterialById(id); }

    /// material mass density
//    virtual double Density() const { return params_->density_; }

//    /// a shear modulus equivalent
//    double ShearMod() const;

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

//    /// hyperelastic stress response plus elasticity tensor (Epetra-type interface)
//    void Evaluate(const Epetra_SerialDenseVector* glstrain_e,
//                  Epetra_SerialDenseMatrix* cmat_e,
//                  Epetra_SerialDenseVector* stress_e,
//
//                  //---
//                  const int gp,   ///< current Gauss point
//                  Teuchos::ParameterList& params  ///< parameter list for communication
//                  //---
//                  );

    virtual void Setup(const int numgp,DRT::INPUT::LineDefinition* linedef);
    /// setup patient-specific AAA stuff
//    void SetupAAA(Teuchos::ParameterList& params);
    //void Setupvisco(const int numgp);

    // update internal stress variables
    virtual void Update();

    // reset internal stress variables
    virtual void Reset();

    /// return if anisotropic
//    bool Anisotropic() const { return anisotropic_;}

    /// return a1s
//    LINALG::Matrix<3,1> Geta1() const {return a1_;}

    /// return a2s
//    LINALG::Matrix<3,1> Geta2() const {return a2_;}

    /// Return quick accessible material parameter data
//    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::ViscoGenMax* params_;

    /// visco history stresses
    Teuchos::RCP<vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histstressisoprinccurr_; ///< current stress
    Teuchos::RCP<vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histstressisoprinclast_; ///< stress of last converged state
    Teuchos::RCP<vector<LINALG::Matrix<NUM_STRESS_3D,1> > > artstressisoprinccurr_;  ///< current artificial stress
    Teuchos::RCP<vector<LINALG::Matrix<NUM_STRESS_3D,1> > > artstressisoprinclast_;  ///< artificial stress in last converged state
    Teuchos::RCP<vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histstressisomodisocurr_; ///< current stress
    Teuchos::RCP<vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histstressisomodisolast_; ///< stress of last converged state
    Teuchos::RCP<vector<LINALG::Matrix<NUM_STRESS_3D,1> > > artstressisomodisocurr_;  ///< current artificial stress
    Teuchos::RCP<vector<LINALG::Matrix<NUM_STRESS_3D,1> > > artstressisomodisolast_;  ///< artificial stress in last converged state
    Teuchos::RCP<vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histstressisomodvolcurr_; ///< current stress
    Teuchos::RCP<vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histstressisomodvollast_; ///< stress of last converged state
    Teuchos::RCP<vector<LINALG::Matrix<NUM_STRESS_3D,1> > > artstressisomodvolcurr_;  ///< current artificial stress
    Teuchos::RCP<vector<LINALG::Matrix<NUM_STRESS_3D,1> > > artstressisomodvollast_;  ///< artificial stress in last converged state
    Teuchos::RCP<vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histstressanisoprinccurr_; ///< current stress
    Teuchos::RCP<vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histstressanisoprinclast_; ///< stress of last converged state
    Teuchos::RCP<vector<LINALG::Matrix<NUM_STRESS_3D,1> > > artstressanisoprinccurr_;  ///< current artificial stress
    Teuchos::RCP<vector<LINALG::Matrix<NUM_STRESS_3D,1> > > artstressanisoprinclast_;  ///< artificial stress in last converged state

    bool isinit_;  ///< indicates if #Initialized routine has been called

//    /// Am I anisotropic?
//    bool anisotropic_;
//
//    /// information about HU from medical image data?
//    double HU_;
//    int HUlumen_;
//
//    /// for Gasser AAA-ILT material
//    double normdist_;

//    /// structural tensors in voigt notation for anisotropy
//    LINALG::Matrix<6,1>  A1_;
//    LINALG::Matrix<6,1>  A2_;
//    /// vectors of local coordinate system
//    LINALG::Matrix<3,1>  a1_;
//    LINALG::Matrix<3,1>  a2_;
//    /// mixed structural tensor (asymmetric) \f$a1 \otimes a2\f$
//    LINALG::Matrix<3,3>  A1A2_;
//
//    /// @name Voigt index arrays
//    //@{
//    static const int VOIGT6ROW_[];  ///< 6-Voigt row index of corresponding 2-tensor
//    static const int VOIGT6COL_[];  ///< 6-Voigt column index of corresponding 2-tensor
//    static const int VOIGT3X3SYM_[]; ///< go from symmtric 2-tensor index pair to 6-Voigt index
//                                     ///< by [3*i+j] for any i,j=0,1,2
//    //@}

    };

}

#endif

#endif
