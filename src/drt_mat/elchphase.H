/*----------------------------------------------------------------------------*/
/*! \file
\brief material that stores a list of ion species in electrolyte solutions

\level 2

\maintainer Christoph Schmidt

*/
/*----------------------------------------------------------------------------*/

#ifndef ELCHPHASE_H
#define ELCHPHASE_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for convection-diffusion
    class ElchPhase : public Parameter
    {
     public:
      /// standard constructor
      ElchPhase(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~ElchPhase() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// provide ids of the individual mat
      const std::vector<int>& MatIds() const { return matids_; }

      /// provide access to phases by its ID
      Teuchos::RCP<MAT::Material> MatById(const int id) const
      {
        if (not local_)
        {
          std::map<int, Teuchos::RCP<MAT::Material>>::const_iterator m = mat_.find(id);

          if (m == mat_.end())
          {
            dserror("Material %d could not be found", id);
            return Teuchos::null;
          }
          else
            return m->second;
        }
        else
          dserror("This is not allowed");

        return Teuchos::null;
      }

      /// @name material parameters
      //@{

      /// porosity
      const double epsilon_;

      /// tortuosity
      const double tortuosity_;

      /// number of materials
      const int nummat_;

      /// the list of material IDs
      const std::vector<int> matids_;

      /// flag for individual materials or only one at global scope
      bool local_;

      //@}

     private:
      /// map to materials (only used for local_==true)
      std::map<int, Teuchos::RCP<MAT::Material>> mat_;

    };  // class ElchPhase

  }  // namespace PAR

  class ElchPhaseType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "ElchPhaseType"; }

    static ElchPhaseType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ElchPhaseType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for the material properties of an ion species in an electrolyte solution
  class ElchPhase : public Material
  {
   public:
    /// construct empty material object
    ElchPhase();

    /// construct the material object given material parameters
    explicit ElchPhase(MAT::PAR::ElchPhase* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ElchPhaseType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_elchphase; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new ElchPhase(*this)); }

    /// return constant porosity
    double Epsilon() const { return params_->epsilon_; }
    /// return constant tortuosity
    double Tortuosity() const { return params_->tortuosity_; }

    int NumMat() const { return params_->nummat_; }

    /// material ID by Index
    int MatID(const unsigned index) const
    {
      if ((int)index < params_->nummat_)
        return params_->matids_.at(index);
      else
      {
        dserror("Index too large");
        return -1;
      }
    }

    /// provide access to material by its ID
    Teuchos::RCP<MAT::Material> MatById(const int id) const
    {
      if (params_->local_)
      {
        std::map<int, Teuchos::RCP<MAT::Material>>::const_iterator m = mat_.find(id);
        if (m == mat_.end())
        {
          dserror("Material %d could not be found", id);
          return Teuchos::null;
        }
        else
          return m->second;
      }
      else  // material is global (stored in material parameters)
        return params_->MatById(id);
    }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// setup of material map
    void SetupMatMap();

    /// clear everything
    void Clear();

    /// my material parameters
    MAT::PAR::ElchPhase* params_;

    /// map to materials
    std::map<int, Teuchos::RCP<MAT::Material>> mat_;
  };

}  // namespace MAT

#endif
