/*----------------------------------------------------------------------*/
/*!
\file maxwell_0d_acinus.H

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef MAXWELL_0D_ACINUS_H
#define MAXWELL_0D_ACINUS_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  
  namespace PAR
  {
    
    /*----------------------------------------------------------------------*/
    /// material parameters for Maxwell 0D acinar material
    ///
    class Maxwell_0d_acinus : public Parameter
    {
    public:
      
      /// standard constructor
      Maxwell_0d_acinus(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~Maxwell_0d_acinus() { }
      
      /// @name material parameters
      //@{
      /// first stiffness of the Maxwell model
      const double stiffness1_;
      /// first stiffness of the Maxwell model
      const double stiffness2_;
      /// first viscosity of the Maxwell model
      const double viscosity1_;
      /// first viscosity of the Maxwell model
      const double viscosity2_;
      //@}
      
      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class Maxwell_0d_acinus
    
  }  // namespace PAR

  class Maxwell_0d_acinusType : public DRT::ParObjectType
  {
  public:

    std::string Name() const { return "maxwell_0d_acinusType"; }
    
    static Maxwell_0d_acinusType & Instance() { return instance_; };
    
    virtual DRT::ParObject* Create( const std::vector<char> & data );
    
  private:
    
    static Maxwell_0d_acinusType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for Maxwell 0D acinar material
  ///
  /// This object exists (several times) at every element
  class Maxwell_0d_acinus : public Material
  {
  public:
    
    /// construct empty material object
    Maxwell_0d_acinus();

    /// construct the material object given material parameters
    explicit Maxwell_0d_acinus(MAT::PAR::Maxwell_0d_acinus* params);
    
    //! @name Packing and Unpacking
    
    /*!
      \brief Return unique ParObject id
      
      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return Maxwell_0d_acinusType::Instance().UniqueParObjectId(); }
    
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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_0d_maxwell_acinus; }
    
    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new Maxwell_0d_acinus( *this ) ); }
    

    /// return density
    virtual double Density() const { return -1; }

    /// return first stiffness of the Maxwell model
    double Stiffness1() const { return params_->stiffness1_;}

    /// return first stiffness of the Maxwell model
    double Stiffness2() const { return params_->stiffness2_;}

    /// return first viscosity of the Maxwell model
    double Viscosity1() const { return params_->viscosity1_;}

    /// return first viscosity of the Maxwell model
    double Viscosity2() const { return params_->viscosity2_;}
    
    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }
    
  private:
    
    /// my material parameters
    MAT::PAR::Maxwell_0d_acinus* params_;
  };
  
}

#endif

