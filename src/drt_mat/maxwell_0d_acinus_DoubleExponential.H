/*----------------------------------------------------------------------*/
/*!
\file maxwell_0d_acinus_DoubleExponential.H

\level 3

\maintainer Lena Yoshihara
*/
/*----------------------------------------------------------------------*/
#ifndef MAXWELL_0D_ACINUS_DOUBLEEXPONENTIAL_H
#define MAXWELL_0D_ACINUS_DOUBLEEXPONENTIAL_H


#include "maxwell_0d_acinus.H"
#include "../drt_red_airways/red_airway.H"


namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for Maxwell 0D acinar material
    ///
    class Maxwell_0d_acinus_DoubleExponential : public Maxwell_0d_acinus
    {
     public:
      /// standard constructor
      Maxwell_0d_acinus_DoubleExponential(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~Maxwell_0d_acinus_DoubleExponential() {}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class Maxwell_0d_acinus_DoubleExponential
  }     // namespace PAR


  class Maxwell_0d_acinusDoubleExponentialType : public Maxwell_0d_acinusType
  {
   public:
    std::string Name() const { return "maxwell_0d_acinusDoubleExponentialType"; }

    static Maxwell_0d_acinusDoubleExponentialType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static Maxwell_0d_acinusDoubleExponentialType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for Maxwell 0D acinar material
  ///
  /// This object exists (several times) at every element
  class Maxwell_0d_acinus_DoubleExponential : public Maxwell_0d_acinus
  {
   public:
    /// construct empty material object
    Maxwell_0d_acinus_DoubleExponential();

    /// construct the material object given material parameters
    Maxwell_0d_acinus_DoubleExponential(MAT::PAR::Maxwell_0d_acinus* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
    {
      return Maxwell_0d_acinusDoubleExponentialType::Instance().UniqueParObjectId();
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
    virtual INPAR::MAT::MaterialType MaterialType() const
    {
      return INPAR::MAT::m_0d_maxwell_acinus_doubleexponential;
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new Maxwell_0d_acinus(*this));
    }

    /*!
      \brief
    */
    virtual void Setup(DRT::INPUT::LineDefinition* linedef);

    /*!
       \brief
     */
    virtual void Evaluate(Epetra_SerialDenseVector& epnp, Epetra_SerialDenseVector& epn,
        Epetra_SerialDenseVector& epnm, Epetra_SerialDenseMatrix& sysmat,
        Epetra_SerialDenseVector& rhs, Teuchos::ParameterList& params, const double NumOfAcini,
        const double Vo, double time, double dt);

   private:
    double e1_01_;
    double e1_lin1_;
    double e1_exp1_;
    double tau1_;

    double e1_02_;
    double e1_lin2_;
    double e1_exp2_;
    double tau2_;
  };

}  // namespace MAT

#endif
