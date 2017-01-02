/*!----------------------------------------------------------------------
\file matlist.H

\level 1

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
*----------------------------------------------------------------------*/
#ifndef MATLIST_H
#define MATLIST_H



#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for list of materials
    class MatList : public Parameter
    {
    public:

      /// standard constructor
      MatList(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~MatList() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// @name material parameters
      //@{

      /// provide ids of the individual materials
      const std::vector<int>* MatIds() const {return matids_;}

      /// provide access to material by its ID
      Teuchos::RCP<MAT::Material> MaterialById(const int id) const;

      std::map<int,Teuchos::RCP<MAT::Material> >* MaterialMapWrite() {return &mat_;}

      /// length of material list
      const int nummat_;

      /// the list of material IDs
      const std::vector<int>* matids_;

      /// flag for individual materials or only one at global scope
      bool local_;

    private:

      /// map to materials (only used for local_==true)
      std::map<int,Teuchos::RCP<MAT::Material> > mat_;

      //@}

    };  // class MatList

  }  // namespace PAR

class MatListType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "MatListType"; }

  static MatListType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static MatListType instance_;
};

  /*----------------------------------------------------------------------*/
  /// Wrapper for a list of materials
  class MatList : public Material
  {
  public:

    /// construct empty material object
    MatList();

    /// construct the material object given material parameters
    explicit MatList(MAT::PAR::MatList* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return MatListType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_matlist; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new MatList( *this ) ); }

    /// number of materials
    int NumMat() const { return params_->nummat_; }

    /// material ID by Index
    int MatID( const unsigned index ) const;

    /// provide access to material by its ID
    virtual Teuchos::RCP<MAT::Material> MaterialById(const int id) const;

    /// Return quick accessible material parameter data
    virtual MAT::PAR::MatList* Parameter() const { return params_; }

  protected:

    ///return pointer to the materials map, which has read-only access.
    const std::map<int,Teuchos::RCP<MAT::Material> >* MaterialMapRead() const {return &mat_;}

    ///return pointer to the materials map, which has read and write access.
    std::map<int,Teuchos::RCP<MAT::Material> >* MaterialMapWrite() {return &mat_;}

  private:

    /// setup of material map
    void SetupMatMap();

    /// clear everything
    void Clear();

    /// my material parameters
    MAT::PAR::MatList* params_;

    /// map to materials
    std::map<int,Teuchos::RCP<MAT::Material> > mat_;

  };

}

#endif

