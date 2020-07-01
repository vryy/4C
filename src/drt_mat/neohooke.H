/*----------------------------------------------------------------------*/
/*! \file
\brief contains the functions to establish local material law
       stress-strain law for isotropic material for a 3D hex element
       following compressible Neo Hookean material law
       see Holzapfel, Nonlinear Solid Mechanics, pp. 243
       example input line:
       MAT 1 MAT_Struct_NeoHooke  YOUNG 100.0 NUE 0.3 DENS 1.0

       USE ELASTHYPER-Toolbox!!
\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef NEOHOOKE_H_
#define NEOHOOKE_H_


#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for neo-Hooke
    class NeoHooke : public Parameter
    {
     public:
      /// standard constructor
      NeoHooke(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~NeoHooke() { ; }

      /// @name material parameters
      //@{

      /// Young's modulus
      const double youngs_;
      /// Possion's ratio
      const double poissonratio_;
      /// mass density
      const double density_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class NeoHooke

  }  // namespace PAR

  class NeoHookeType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "NeoHookeType"; }

    static NeoHookeType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static NeoHookeType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for Neo-Hookean Material
  class NeoHooke : public Material
  {
   public:
    /// construct empty material object
    NeoHooke();

    /// construct the material object given material parameters
    explicit NeoHooke(MAT::PAR::NeoHooke* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return NeoHookeType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_neohooke; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new NeoHooke(*this)); }

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

    /// evaluate material law
    void Evaluate(const Epetra_SerialDenseVector* glstrain_e, Epetra_SerialDenseMatrix* cmat_e,
        Epetra_SerialDenseVector* stress_e);

    /// evaluate material law
    virtual void Evaluate(const LINALG::Matrix<6, 1>& glstrain, LINALG::Matrix<6, 6>& cmat,
        LINALG::Matrix<6, 1>& stress);

    //@}

   private:
    /// my material parameters
    MAT::PAR::NeoHooke* params_;

    void InverseTensor(
        const Epetra_SerialDenseMatrix& M, Epetra_SerialDenseMatrix& Minv, const double I3);
  };
}  // namespace MAT

#endif /*NEOHOOKE_H_*/
