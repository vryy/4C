/*!------------------------------------------------------------------------------------------------*
\brief
Former file of Martin Winklmaier

\level 3

<pre>
\maintainer Martin Kronbichler
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef OPTIMIZATION_DENSITY_H_
#define OPTIMIZATION_DENSITY_H_


#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"



namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// parameters for optimization field
    class TopOptDens : public Parameter
    {
     public:
      /// standard constructor
      TopOptDens(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~TopOptDens() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      double PoroBdDown() const { return poro_bd_down_; };
      double PoroBdUp() const { return poro_bd_up_; };
      double SmearFac() const { return smear_fac_; };

      void UpdateSmearFac(const double new_smearfac)
      {
        smear_fac_ = new_smearfac;
        std::cout << "smoothing factor updated to new value " << smear_fac_ << std::endl;
      };

     private:
      /// @name material parameters
      //@{

      /// lower boundary for pseudo-porosity field
      const double poro_bd_down_;
      /// upper boundary for pseudo-porosity field
      const double poro_bd_up_;
      /// smearing factor for computation of porosity field out of density field
      double smear_fac_;

      //@}

    };  // class TopOptDens

  }  // namespace PAR

  class TopOptDensType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "TopOptDensType"; }

    static TopOptDensType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static TopOptDensType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// wrapper for scalar transport material
  class TopOptDens : public Material
  {
   public:
    /// construct empty material object
    TopOptDens();

    /// construct the material object given material parameters
    explicit TopOptDens(MAT::PAR::TopOptDens* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return TopOptDensType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_opti_dens; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new TopOptDens(*this)); }

    /// return upper boundary of porosity field
    virtual double MaxPoro() const { return params_->PoroBdUp(); }

    /// return upper boundary of porosity field
    virtual double MinPoro() const { return params_->PoroBdDown(); }

    /// return upper boundary of porosity field
    virtual double SmearFactor() const { return params_->SmearFac(); }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// my material parameters
    MAT::PAR::TopOptDens* params_;
  };

}  // namespace MAT

#endif /* OPTIMIZATION_DENSITY_H_ */
