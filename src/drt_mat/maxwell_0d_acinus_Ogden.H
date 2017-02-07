/*----------------------------------------------------------------------*/
/*!
\file maxwell_0d_acinus_Ogden.H

\maintainer Lena Yoshihara

\level 3
*/
/*----------------------------------------------------------------------*/
#ifndef MAXWELL_0D_ACINUS_OGDEN_H
#define MAXWELL_0D_ACINUS_OGDEN_H


#include "maxwell_0d_acinus.H"
#include "../drt_red_airways/red_airway.H"


namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for Maxwell 0D acinar material
    ///
    class Maxwell_0d_acinus_Ogden : public Maxwell_0d_acinus
    {
    public:

      /// standard constructor
      Maxwell_0d_acinus_Ogden(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~Maxwell_0d_acinus_Ogden() { }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// return optimizable parameters
      virtual void OptParams(std::map<std::string, int>* pnames);

      /// enum for mapping between material parameter and entry in the matparams_ vector
      enum matparamnames_ {kappa, beta, first=kappa, last=beta};

  };  // class Maxwell_0d_acinus_Ogden
}  // namespace PAR


class Maxwell_0d_acinusOgdenType : public Maxwell_0d_acinusType
{
  public:

    std::string Name() const { return "maxwell_0d_acinusOgdenType"; }

    static Maxwell_0d_acinusOgdenType & Instance() { return instance_; };

    virtual DRT::ParObject* Create( const std::vector<char> & data );

  private:

    static Maxwell_0d_acinusOgdenType instance_;
  };

/*----------------------------------------------------------------------*/
/// Wrapper for Maxwell 0D acinar material
///
/// This object exists (several times) at every element
class Maxwell_0d_acinus_Ogden : public Maxwell_0d_acinus
{
  public:

    /// construct empty material object
    Maxwell_0d_acinus_Ogden();

    /// construct the material object given material parameters
    Maxwell_0d_acinus_Ogden(MAT::PAR::Maxwell_0d_acinus* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return Maxwell_0d_acinusOgdenType::Instance().UniqueParObjectId(); }


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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_0d_maxwell_acinus_ogden; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new Maxwell_0d_acinus( *this ) ); }

    /*!
      \brief
    */
    virtual void Setup(DRT::INPUT::LineDefinition* linedef);

    /*!
       \brief
     */
    virtual void Evaluate(Epetra_SerialDenseVector& epnp,
                          Epetra_SerialDenseVector& epn,
                          Epetra_SerialDenseVector& epnm,
                          Epetra_SerialDenseMatrix& sysmat,
                          Epetra_SerialDenseVector& rhs,
                          Teuchos::ParameterList & params,
                          const double NumOfAcini,
                          const double Vo,
                          double time,
                          double dt);

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string,int>& names);

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data , int eleID);

    /// Return value of class parameter
    virtual double GetParams(std::string parametername);

    /// Set value of class parameter
    virtual void SetParams(std::string parametername, double new_value);

  private:

     double kappa_;
     double beta_;

  };

}

#endif

