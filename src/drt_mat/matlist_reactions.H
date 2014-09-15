/*!----------------------------------------------------------------------
\file matlist_reactions.cpp

<pre>
Maintainer: Moritz Thon
            thon@mhpc.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-10364
</pre>
*----------------------------------------------------------------------*/

#ifndef MatListReactions_H
#define MatListReactions_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for list of materials
    class MatListReactions : public Parameter
    {
    public:

      /// standard constructor
      MatListReactions(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~MatListReactions() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// @name material parameters
      //@{

      /// provide ids of the individual materials
      const std::vector<int>* MatIds() const {return matids_;}

      /// provide ids of the individual reaction materials
      const std::vector<int>* ReacIds() const {return reacids_;}

      /// provide access to material by its ID
      Teuchos::RCP<MAT::Material> MaterialById(const int id) const
              {
        if (not local_)
        {
          std::map<int,Teuchos::RCP<MAT::Material> >::const_iterator m = mat_.find(id);

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

      /// flag for individual materials or only one at global scope
      bool local_;

      /// length of material list
      const int nummat_;

      /// the list of material IDs
      const std::vector<int>* matids_;

      /// length of reaction list
      const int numreac_;

      /// the list of reaction IDs
      const std::vector<int>* reacids_;

    private:

      /// map to materials (only used for local_==true)
      std::map<int,Teuchos::RCP<MAT::Material> > mat_;

      //@}

    };  // class MatListReactions

  }  // namespace PAR

class MatListReactionsType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "MatListReactionsType"; }

  static MatListReactionsType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static MatListReactionsType instance_;
};

  /*----------------------------------------------------------------------*/
  /// Wrapper for a list of materials
  class MatListReactions : public Material
  {
  public:

    /// construct empty material object
    MatListReactions();

    /// construct the material object given material parameters
    explicit MatListReactions(MAT::PAR::MatListReactions* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return MatListReactionsType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_matlist_reactions; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new MatListReactions( *this ) ); }

    /// number of materials
    int NumMat() const { return params_->nummat_; }

    /// number of reactions
    int NumReac() const { return params_->numreac_; }

    /// material ID by Index
    int MatID(
      const unsigned index
    ) const
    {
      if ((int)index < params_->nummat_)
        return params_->matids_->at(index);
      else
      {
        dserror("Index too large");
        return -1;
      }
    }

    /// reaction ID by Index
    int ReacID(
      const unsigned index
    ) const
    {
      if ((int)index < params_->numreac_)
        return params_->reacids_->at(index);
      else
      {
        dserror("Index too large");
        return -1;
      }
    }

    /// provide access to material by its ID
    Teuchos::RCP<MAT::Material> MaterialById(const int id) const
    {
      if (params_->local_)
      {
      std::map<int,Teuchos::RCP<MAT::Material> >::const_iterator m = mat_.find(id);
      if (m == mat_.end())
      {
        dserror("Material %d could not be found", id);
        return Teuchos::null;
      }
      else
        return m->second;
      }
      else // material is global (stored in material parameters)
        return params_->MaterialById(id);
     }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// setup of material map
    void SetupMatMap();

    /// clear everything
    void Clear();

    /// my material parameters
    MAT::PAR::MatListReactions* params_;

    /// map to materials
    std::map<int,Teuchos::RCP<MAT::Material> > mat_;

  };

}

#endif

