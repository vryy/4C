/*----------------------------------------------------------------------*/
/*! \file
 \brief

This file contains the material for reactive AND chemotactic scalars. It is
in diamond inheritance with MatListReactions and MatListChemotaxis,
which govern the actual doings

\level 3
*----------------------------------------------------------------------*/


#ifndef MATLIST_CHEMOREAC_H_
#define MATLIST_CHEMOREAC_H_



#include "matpar_parameter.H"
#include "material.H"
#include "matlist.H"
#include "matlist_reactions.H"
#include "matlist_chemotaxis.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for list of materials
    class MatListChemoReac : public MatListReactions, public MatListChemotaxis
    {
     public:
      /// standard constructor
      MatListChemoReac(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~MatListChemoReac() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// @name material parameters

    };  // class MatListReactions

  }  // namespace PAR

  class MatListChemoReacType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "MatListChemoReacType"; }

    static MatListChemoReacType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static MatListChemoReacType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for a list of materials
  class MatListChemoReac : public MatListChemotaxis, public MatListReactions
  {
   public:
    /// construct empty material object
    MatListChemoReac();

    /// construct the material object given material parameters
    explicit MatListChemoReac(MAT::PAR::MatListChemoReac* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
    {
      return MatListChemoReacType::Instance().UniqueParObjectId();
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
      return INPAR::MAT::m_matlist_chemoreac;
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new MatListChemoReac(*this));
    }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::MatListChemoReac* Parameter() const { return paramsreachemo_; }

   private:
    /// setup of material map
    virtual void SetupMatMap();

    /// clear everything
    void Clear();

    /// my material parameters
    MAT::PAR::MatListChemoReac* paramsreachemo_;
  };

}  // namespace MAT


#endif /* MATLIST_CHEMOREAC_H_ */
