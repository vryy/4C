/*!
\file dofkey.H

\brief provides a class with a unique identifier for a degree of freedom (DOF)

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */
#ifdef CCADISCRET
#ifndef DOFKEY_H
#define DOFKEY_H


#include "field_enriched.H"


namespace XFEM
{

class DofKeyType : public DRT::ParObjectType
{
public:

  std::string Name() const;

  static DofKeyType & Instance(){ return instance_; };

private:

  static DofKeyType instance_;
};


  /*!
   * \brief A unique identifier for a degree of freedom (DOF)
   *
   * A DofKey identifies a DOF uniquely without generating a global unique integer id.
   * It contains all information that is necessary to know how to treat this DOF.
   *
   * In Haskell one would write:
   * type DofKey = (gid, FieldEnr)
   */
  class DofKey : public DRT::ParObject
  {
    public:
      //! standard constructor
      explicit DofKey(
          const int             gid,      ///< global id
          const FieldEnr&       fieldenr  ///< enriched field
      );

      //! constructor, get values from char array using Unpack member function
      explicit DofKey(
          vector<char>& data      ///< char array with data to unpack
      );

      //! copy constructor
      DofKey(
          const DofKey& other           ///< source
      );

      //! return string representation of this class
      std::string toString() const;

      virtual int UniqueParObjectId() const;

      //! pack all content to char array
      virtual void Pack(DRT::PackBuffer& data) const;

      //! unpack all content from char array
      virtual void Unpack(const vector<char>& data);

      //! return global node or element Id
      int getGid() const;

      //! return enriched field
      const FieldEnr& getFieldEnr() const;

      //! compare: dofkeys are ordered first by the global id and then by their enrichment
      bool operator <(const DofKey& rhs) const;

      //! test for equality: only equal, if gid and enriched field are identical
      bool operator ==(const DofKey& rhs) const;

      //! test for inequality: if any member differs, the objects are not equal
      bool operator !=(const DofKey& rhs) const;

      friend ostream& operator<<(ostream& out, const DofKey& key)
      {
        out << key.toString();
        return out;
      }

    private:

      //! don't want default constructor
      explicit DofKey();

      //! no assignment operator implemented
      DofKey& operator = (const DofKey& old);

      //! specific enriched field
      XFEM::FieldEnr  fieldenr_;

      //! unique id among the choosen type of dof (for example: nodal global id)
      int             gid_;

  };
} // namespace XFEM

#endif  // #ifndef DOFKEY_H
#endif  // #ifdef CCADISCRET
