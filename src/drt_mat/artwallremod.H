/*!----------------------------------------------------------------------
\file artwallremod.H
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*----------------------------------------------------------------------*/
#ifndef ARTWALLREMOD_H_
#define ARTWALLREMOD_H_

#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Teuchos_RCP.hpp>
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
    class ArtWallRemod : public Parameter
    {
    public:

      /// standard constructor
      ArtWallRemod(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~ArtWallRemod() { ; }

      /// @name material parameters
      //@{

      const double kappa_;
      const double mue_;
      const double density_;
      const double k1_;
      const double k2_;
      const double gamma_;
      const int init_;
      double rembegt_;  // this is changed!!!
      const int tensonly_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class ArtWallRemod

  }  // namespace PAR

class ArtWallRemodType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "ArtWallRemodType"; }

  static ArtWallRemodType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static ArtWallRemodType instance_;
};

  /*----------------------------------------------------------------------*/
  /// Wrapper for Mooney-Rivlin material
  class ArtWallRemod : public Material
  {
  public:

    /// construct empty material object
    ArtWallRemod();

    /// construct the material object given material parameters
    explicit ArtWallRemod(MAT::PAR::ArtWallRemod* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ArtWallRemodType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_artwallremod; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new ArtWallRemod( *this ) ); }

    /// evaluate stress and elasticity matrix
    void Evaluate(
            const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
            const int gp,
            Teuchos::ParameterList& params,
            LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> * cmat,
            LINALG::Matrix<NUM_STRESS_3D,1> * stress,
            const LINALG::Matrix<3,3>& defgrd);


    /// return density
    virtual double Density() const {return params_->density_;}

    /// check whether internal variables are initialized
    bool Initialized() const
      { return isinit_&&(remtime_!=Teuchos::null);}

    /// return lambdas
    RCP<vector<std::vector<double> > > Getlambdas() const {return lambda_;}

    /// return phis
    RCP<vector<Epetra_SerialDenseMatrix> > Getphis() const {return phi_;}

    /// return gammas
    RCP<std::vector<double> > Getgammas() const {return gamma_;}

    /// return remtimes
    RCP<std::vector<double> > Getremtimes() const {return remtime_;}

    /// return a1s
    RCP<vector<std::vector<double> > > Geta1() const {return a1_;}

    /// return a2s
    RCP<vector<std::vector<double> > > Geta2() const {return a2_;}

    /*!
      \brief setup material at element input with potentially local fiber cosy

      This anisotropic material consists of fiber angles which can be defined
      at element level via local coordinate systems. This method controls setup
      of corresponding variables and reads e.g. fiber directions into its history
      variables.

      \param numgp (in) : number of gausspoints for this element
      \param eleid (in) : current eleid, used to init random variables
    */

    void Setup(const int numgp, const int eleid, DRT::INPUT::LineDefinition* linedef);

    /*!
      \brief Remodeling driver

      Here the remodeling (reorientation) of fiber directions takes place, such
      as calculating stress eigenvectors and alignment angles

      \param gp (in) : current gp
      \param time (int) : current time
    */
    void Remodel(const int gp, const double time, const LINALG::Matrix<3,3>& defgrd);

    /// evaluate fiber directions from locsys and alignment angle
    void EvaluateFiberVecs(const int gp, const double gamma, const Epetra_SerialDenseMatrix& locsys);

    /// Print vector to string
    std::string PrintVec(const vector<double> actvec);

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// my material parameters
    MAT::PAR::ArtWallRemod* params_;

    /// angle of alignment at gp
    Teuchos::RCP<std::vector<double> > gamma_;
    /// first fiber vector per gp
    Teuchos::RCP<vector<std::vector<double> > > a1_;
    /// second fiber vector per gp
    Teuchos::RCP<vector<std::vector<double> > > a2_;
    /// eigenvectors per gp
    Teuchos::RCP<vector<Epetra_SerialDenseMatrix> > phi_;
    /// eigenvalues per gp
    Teuchos::RCP<vector<std::vector<double> > > lambda_;
    /// current stress
    Teuchos::RCP<vector<Epetra_SerialDenseMatrix> > stresses_;

    /// remember time
    Teuchos::RCP<std::vector<double> > remtime_;

    /// init flag of internal variables
    bool isinit_;

  };

void ArtWallRemodOutputToGmsh(const Teuchos::RCP<DRT::Discretization> dis,
                                        const double time,
                                        const int iter);

void ArtWallRemodOutputToTxt(const Teuchos::RCP<DRT::Discretization> dis,
    const double time,
    const int iter);


}


#endif /*ARTWALLREMOD_H_*/
