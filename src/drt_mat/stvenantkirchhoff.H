/*----------------------------------------------------------------------*/
/*! \file
\brief
St. Venant-Kirchhoff material

\level 2

\maintainer Amadeus Gebauer
*/
/*----------------------------------------------------------------------*/
#ifndef STVENANTKIRCHHOFF_H_
#define STVENANTKIRCHHOFF_H_


#include "matpar_parameter.H"
//#include "material.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Teuchos_ParameterList.hpp>


namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for de St. Venant--Kirchhoff
    class StVenantKirchhoff : public Parameter
    {
     public:
      /// standard constructor
      StVenantKirchhoff(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~StVenantKirchhoff() { ; }

      /// @name material parameters
      //@{

      /// Young's modulus
      const double youngs_;
      /// Possion's ratio
      const double poissonratio_;
      /// mass density
      const double density_;
      /// linear coefficient of thermal expansion
      const double thermexpans_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class StVenantKirchhoff

  }  // namespace PAR

  class StVenantKirchhoffType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "StVenantKirchhoffType"; }

    static StVenantKirchhoffType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static StVenantKirchhoffType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for St.-Venant-Kirchhoff material
  class StVenantKirchhoff : public So3Material
  {
   public:
    /// construct empty material object
    StVenantKirchhoff();

    /// construct the material object given material parameters
    explicit StVenantKirchhoff(MAT::PAR::StVenantKirchhoff* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
    {
      return StVenantKirchhoffType::Instance().UniqueParObjectId();
    }

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

    //! @name Access methods

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_stvenant; }

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem)
    {
      if (!(kinem == INPAR::STR::kinem_linear || kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new StVenantKirchhoff(*this));
    }

    /// Young's modulus
    double Youngs() const { return params_->youngs_; }

    /// Poisson's ratio
    double PoissonRatio() const { return params_->poissonratio_; }

    /// density
    virtual double Density() const { return params_->density_; }

    /// shear modulus
    double ShearMod() const { return 0.5 * params_->youngs_ / (1.0 + params_->poissonratio_); }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    //@}

    //! @name Evaluation methods

    /// evaluates material law
    void Evaluate(const Epetra_SerialDenseVector* glstrain_e, Epetra_SerialDenseMatrix* cmat_e,
        Epetra_SerialDenseVector* stress_e);

    /// evaluate material law
    virtual void Evaluate(const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* glstrain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat,
        const int eleGID);

    /// evaluate strain energy function
    virtual void StrainEnergy(const LINALG::Matrix<6, 1>& glstrain, double& psi, const int eleGID);

    /// evaluate for GEMM
    virtual void EvaluateGEMM(LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* stress,
        LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>* cmat, double* density,
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_m,
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_new,
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_old, LINALG::Matrix<3, 3>* rcg_new,
        LINALG::Matrix<3, 3>* rcg_old, const int eleGID);

    // computes isotropic elasticity tensor in matrix notion for 3d
    void SetupCmat(LINALG::Matrix<6, 6>& cmat);

    // computes isotropic elasticity tensor in matrix notion for 2d
    // plane strain, rotational symmetry
    void SetupCmat2d(Epetra_SerialDenseMatrix* cmat);

    //@}

   private:
    /// my material parameters
    MAT::PAR::StVenantKirchhoff* params_;
  };
}  // namespace MAT

#endif /*STVENANTKIRCHHOFF_H_*/
