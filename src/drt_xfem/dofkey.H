/*!
\file dofkey.H

\brief provides a class with a unique identifier for a degree of freedom (DOF)

     The following remark is "from the past" where element data existed on global level

     remark: Only if DLM condensation is turned off, there will be element dofs visible in the global
             DofManager. Of course, irrespective of the DLM condensation setting, there will be
             element dofs on the element level (ElementDofManager), if a stress-based approach is
             chosen to handle the boundary terms.
             So, if there are element dofs, they are not known here, since 'elementalDofs_' in
             fillDofRowDistributionMaps() only knows the different field enrichments for a specific
             element, that is

             e.g.
             element 103 has 28 element dofs ((6 x stess + 1 x disc pressure) x element ansatz (e.g. quad4))

             but 'elementalDofs_' contains only 7 "dof-type entries"
             <ele 103, Tauxx + quad4>
             <ele 103, Tauxy + quad4>
             ...
             <ele 103, DiscPres + quad4>

             That is, element dofs cannot be switched correctly here, since not all element dofs
             are known. They only exist on the element level. They are generated in
             ElementDofManager::ComputeDependentInfo()

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
 */


#ifndef DOFKEY_H
#define DOFKEY_H


#include "field_enriched.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_parobjectfactory.H"


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
          std::vector<char>& data      ///< char array with data to unpack
      );

      //! copy constructor
      DofKey(
          const DofKey& other           ///< source
      );

      //! return std::string representation of this class
      std::string toString() const;

      virtual int UniqueParObjectId() const;

      //! pack all content to char array
      virtual void Pack(DRT::PackBuffer& data) const;

      //! unpack all content from char array
      virtual void Unpack(const std::vector<char>& data);

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

      friend std::ostream& operator<<(std::ostream& out, const DofKey& key)
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
