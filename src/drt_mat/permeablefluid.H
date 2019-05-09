/*----------------------------------------------------------------------*/
/*!
\brief permeable fluid

\level 2

\maintainer Volker Gravemeier
*/
/*----------------------------------------------------------------------*/
#ifndef PERMEABLEFLUID_H
#define PERMEABLEFLUID_H



#include "material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for Permeable fluid
    ///
    /// This object exists only once for each read Newton fluid.
    class PermeableFluid : public Parameter
    {
     public:
      /// standard constructor
      PermeableFluid(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~PermeableFluid() { ; }

      /// @name material parameters
      //@{

      /// problem type: Darcy or Darcy-Stokes
      const std::string* type_;
      /// kinematic or dynamic viscosity
      const double viscosity_;
      /// density
      const double density_;
      /// permeability
      const double permeability_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class PermeableFluid

  }  // namespace PAR

  class PermeableFluidType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "PermeableFluidType"; }

    static PermeableFluidType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static PermeableFluidType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for Permeable fluid material
  ///
  /// This object exists (several times) at every element
  class PermeableFluid : public Material
  {
   public:
    /// construct empty material object
    PermeableFluid();

    /// construct the material object given material parameters
    explicit PermeableFluid(MAT::PAR::PermeableFluid* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
    {
      return PermeableFluidType::Instance().UniqueParObjectId();
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

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_permeable_fluid; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new PermeableFluid(*this)); }

    /// compute reaction coefficient
    double ComputeReactionCoeff() const;

    /// set viscosity (zero for Darcy and greater than zero for Darcy-Stokes)
    double SetViscosity() const;

    /// return type
    std::string Type() const { return *params_->type_; }

    /// return density
    virtual double Density() const { return params_->density_; }

    /// return viscosity
    double Viscosity() const { return params_->viscosity_; }

    /// return permeability
    double Permeability() const { return params_->permeability_; }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// my material parameters
    MAT::PAR::PermeableFluid* params_;
  };

}  // namespace MAT

#endif
