/*!----------------------------------------------------------------------
\file contchainnetw.H
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*----------------------------------------------------------------------*/
#ifndef CONTCHAINNETW_H_
#define CONTCHAINNETW_H_

#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_so3/so_hex8.H"
#include "../drt_so3/so_sh8.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters
    class ContChainNetw : public Parameter
    {
    public:

      /// standard constructor
      ContChainNetw(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~ContChainNetw() { ; }

      /// @name material parameters
      //@{

      const double lambda_;
      const double mue_;
      const double density_;
      const double nchain_;
      const double abstemp_;
      const double contl_l_;
      const double persl_a_;
      const double r0_;
      const double relax_;
      const int initran_;
      const double rembegt_;
      const int updrate_;
      const double difftol_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class ContChainNetw

  }  // namespace PAR

class ContChainNetwType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "ContChainNetwType"; }

  static ContChainNetwType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static ContChainNetwType instance_;
};


  /*----------------------------------------------------------------------*/
  /// Wrapper for Mooney-Rivlin material
  class ContChainNetw : public Material
  {
  public:

    /// construct empty material object
    ContChainNetw();

    /// construct the material object given material parameters
    explicit ContChainNetw(MAT::PAR::ContChainNetw* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ContChainNetwType::Instance().UniqueParObjectId(); }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    virtual void Pack(DRT::PackBuffer& data) const;

    /*!
      \brief Unpack data from a char vector into this class

      The vector data contains all information to rebuild the
      exact copy of an instance of a class on a different processor.
      The first entry in data has to be an integer which is the unique
      parobject id defined at the top of this file and delivered by
      UniqueParObjectId().

      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    virtual void Unpack(const std::vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_contchainnetw; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new ContChainNetw( *this ) ); }

    /// initialize material and remodeling
    void Initialize(const int numgp, const int eleid) ;

    /// evaluate stress and elasticity matrix
    void Evaluate(const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
                                      const int gp,
                                      Teuchos::ParameterList& params,
                                      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> * cmat,
                                      LINALG::Matrix<NUM_STRESS_3D,1> * stress,
                                      int eleId);

    /// return density
    virtual double Density() const {return params_->density_;}

    /// check whether material history is initialized
    bool Initialized() const
      { return isinit_&&(li_!=Teuchos::null);}

    /// return current chain-cell lengths
    RCP<vector<vector<double> > > Getli() const {return li_;}

    /// return initial chain-cell lengths
    RCP<vector<vector<double> > > Getli0() const {return li0_;}

    /// return eigenvalues
    RCP<vector<vector<double> > > Getlambdas() const {return lambda_;}

    /// return chain-cell basis
    RCP<vector<LINALG::Matrix<3,3> > >Getni() const {return ni_;}

    /// evaluate fiber-part of PK2-stress
    LINALG::Matrix<3,3> EvaluateStress(
        const LINALG::Matrix<3,3>& isostress,
        const vector<LINALG::Matrix<3,3> >& Ni,
        const vector<double>& cell_li,
        const vector<double>& cell_Inv,
        const double s_chn_scalar,
        const double stressfree);
    /// evaluate and add fiber-part of PK2-stress
    void UpdateStress(
        LINALG::Matrix<3,3>& stress,
        const vector<LINALG::Matrix<3,3> >& Ni,
        const vector<double>& cell_li,
        const vector<double>& cell_Inv,
        const double s_chn_scalar,
        const double stressfree);

    /// evaluate fiber-part of Cmat
    void EvaluateCmat(LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat,
        const vector<LINALG::Matrix<3,3> >& Ni,
        const vector<double>& cell_li,
        const vector<double>& cell_Inv,
        const double c_chn_scalar,
        const double stressfree);

    /// evaluate structural tensors from chain-cell basis
    vector<LINALG::Matrix<3,3> > EvaluateStructTensors(const int gp);

    /// evaluate fiber invariants
    vector<double> EvaluateInvariants(
        const LINALG::Matrix<3,3> & CG,
        const vector<LINALG::Matrix<3,3> >& Ni);

    /// update fiber basis (remodeling evolution)
    void Update(const Epetra_SerialDenseMatrix& Phi,
        const Epetra_SerialDenseVector& lambda,
        const vector<double> rem_toggle,
        const int gp, const double r0, const double decay);

    /// update fiber basis (remodeling evolution in rate from)
    void UpdateRate(const Epetra_SerialDenseMatrix& Phi,
        const Epetra_SerialDenseVector& lambda,
        const vector<double> rem_toggle,
        const int gp, const double r0, const double decay, const double kappa, const double dt);

    /// print structural tensors to string
    std::string PrintStructTens(const vector<Epetra_SerialDenseMatrix>& Ni);
    /// print chain-cell basis to string
    std::string PrintAnisoVects(const int gp);
    /// print anisotropic part of Cmat to string
    std::string PrintAnisoCmat(const LINALG::Matrix<6,6>& cmat,
        const vector<LINALG::Matrix<3,3> >& Ni,
        const vector<double>& cell_li,
        const vector<double>& cell_Inv,
        const double c_chn_scalar,
        const double stressfree);
    /// print vector to string
    std::string PrintVec(const vector<double> actvec);

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::ContChainNetw* params_;

    /// unit cell dimensions per gp
    Teuchos::RCP<vector<vector<double> > > li_;
    /// unit cell basis in undeformed config. per gp
    Teuchos::RCP<vector<LINALG::Matrix<3,3> > > ni_;
    /// initial unit cell dimensions per gp
    Teuchos::RCP<vector<vector<double> > > li0_;
    /// eigenvalues per gp
    Teuchos::RCP<vector<vector<double> > > lambda_;
    /// current stress
    Teuchos::RCP<vector<LINALG::Matrix<3,3> > > stresses_;

    /// remember time for remodeling evolution
    Teuchos::RCP<vector<double> > mytime_;

    /// flag for history init
    bool isinit_;

  };

  Epetra_SerialDenseMatrix StressVoigt2Mat(Epetra_SerialDenseVector* stress);
  LINALG::Matrix<3,3> StressVoigt2Mat(const LINALG::Matrix<6,1>* stress);
  void StressVoigt2Mat(const LINALG::Matrix<6,1>* stress, LINALG::Matrix<3,3>& mat);
  Epetra_SerialDenseVector StressMat2Voigt(Epetra_SerialDenseMatrix& stressmat);
  LINALG::Matrix<6,1> StressMat2Voigt(LINALG::Matrix<3,3>& stressmat);
  void StressMat2Voigt(LINALG::Matrix<6,1>& stress,const LINALG::Matrix<3,3>& mat);
  Epetra_SerialDenseMatrix StrainVoigt2Mat(const Epetra_SerialDenseVector* strain);
  void ChainOutputToGmsh(const Teuchos::RCP<DRT::Discretization> dis, const double time, const int iter);
  void ChainOutputToTxt(const Teuchos::RCP<DRT::Discretization> dis, const double time, const int iter);
  const vector<double> MatPointCoords(const DRT::Element* actele,const vector<double>& mydisp, int gp=-1);
}

#endif /*CONTCHAINNETW_H_*/
