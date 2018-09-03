/*!----------------------------------------------------------------------*/
/*!
\file elchmat.H

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef ELCHMAT_H
#define ELCHMAT_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for list of materials
    class ElchMat : public Parameter
    {
     public:
      /// standard constructor
      ElchMat(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~ElchMat() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// @name material parameters
      //@{

      /// provide ids of the individual phase
      const std::vector<int>& PhaseIds() const { return phaseids_; }

      /// provide access to phases by its ID
      Teuchos::RCP<MAT::Material> PhaseById(const int id) const
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

      /// number of degrees of freedom
      const int numdof_;

      /// number of scalar
      const int numscal_;

      /// length of phase list
      const int numphase_;

      /// the list of material IDs
      const std::vector<int> phaseids_;

      /// flag for individual materials or only one at global scope
      bool local_;

     private:
      /// map to materials (only used for local_==true)
      std::map<int, Teuchos::RCP<MAT::Material>> mat_;

      //@}

    };  // class MatList

  }  // namespace PAR

  class ElchMatType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "ElchMatType"; }

    static ElchMatType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static ElchMatType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for a list of materials
  class ElchMat : public Material
  {
   public:
    /// construct empty material object
    ElchMat();

    /// construct the material object given material parameters
    explicit ElchMat(MAT::PAR::ElchMat* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ElchMatType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_elchmat; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new ElchMat(*this)); }

    // return number of Dof used for this problem type
    int NumDOF() const { return params_->numdof_; }

    // return number of scalars used for this problem type
    int NumScal() const { return params_->numscal_; }

    /// number of materials
    int NumPhase() const { return params_->numphase_; }

    /// material ID by Index
    int PhaseID(const unsigned index) const
    {
      if ((int)index < params_->numphase_)
        return params_->phaseids_.at(index);
      else
      {
        dserror("Index too large");
        return -1;
      }
    }

    /// provide access to material by its ID
    Teuchos::RCP<MAT::Material> PhaseById(const int id) const
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
        return params_->PhaseById(id);
    }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// setup of material map
    void SetupMatMap();

    /// clear everything
    void Clear();

    /// my material parameters
    MAT::PAR::ElchMat* params_;

    /// map to materials
    std::map<int, Teuchos::RCP<MAT::Material>> mat_;
  };

}  // namespace MAT

#endif
