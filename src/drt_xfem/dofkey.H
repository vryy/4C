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


#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_parobject.H"
#include "field_enriched.H"


namespace XFEM
{
  
  //! Allows to distinguish between nodal and elemental dofs
  enum DofType
  {
    onNode,
    onElem
  };
  
  
  //! create a string from the DofType enum
  inline std::string dofTypeToString(DofType t)
  {
    switch (t)
    {
      case onNode: return "Node";
      case onElem: return "Elem";
      default:
        dserror("undefined behaviour");
        exit(1);
    }
  }
  
  
  /*!
   * \brief A unique identifier for a degree of freedom (DOF)
   * 
   * in Haskell one would write:
   * type DofKey = (gid, FieldEnr)
   * 
   * \tparam doftype to distinguish between element and nodal dofs, now the compiler complains
   *                 when we mix them, since DofKey<onNode> is different from DofKey<onElement>.
   */
  template <DofType doftype>
  class DofKey : public DRT::ParObject
  {
    public:
      //! standard constructor
      explicit DofKey(
          const int             gid,      ///< global id
          const FieldEnr        fieldenr  ///< enriched field
      ) :
        gid_(gid),
        fieldenr_(fieldenr)
        {
            return;
        }

      //! constructor, get values from char array using Unpack member function
      explicit DofKey(
          vector<char>& data      ///< global id
      )
      {
        Unpack(data);
        return;
      }
          
      //! copy constructor
      DofKey(
          const DofKey& other           ///< source
      ) :
        gid_(other.gid_),
        fieldenr_(other.fieldenr_)
      {
          assert(&other != this);
          return;
      }
      
      //! return string representation of this class
      std::string toString() const
      {
        std::stringstream s;
        s << "Dofkey: [" << dofTypeToString(doftype) << " " << setw(3) << gid_ << ", " << fieldenr_.toString() << "]";
        return s.str();
      };
      
      virtual int UniqueParObjectId() const { return ParObject_DofKey;}

      //! pack all content to char array
      virtual void Pack(vector<char>& data) const
      {
        data.resize(0);
        AddtoPack(data,UniqueParObjectId());
        AddtoPack(data,gid_);
        AddtoPack(data,fieldenr_.getField());
        const int label = fieldenr_.getEnrichment().XFEMConditionLabel();
        AddtoPack(data,label);

        const XFEM::Enrichment::EnrType enrtype = fieldenr_.getEnrichment().Type();
        AddtoPack(data,enrtype);
      }
      
      //! unpack all content from char array
      virtual void Unpack(const vector<char>& data)
      {
        int position = 0;
        // extract type
        int type = 0;
        ExtractfromPack(position,data,type);
        if (type != UniqueParObjectId()) dserror("wrong instance type data");
        
        ExtractfromPack(position,data,gid_);
        
        // read field
        int field_int;
        ExtractfromPack(position,data,field_int);
        const PHYSICS::Field field = static_cast<PHYSICS::Field>(field_int);
        
        // read xfem condition label
        int label;
        ExtractfromPack(position,data,label);
        
        // read enrichment type
        int enrtype_int;
        ExtractfromPack(position,data,enrtype_int);
        const XFEM::Enrichment::EnrType enrtype = static_cast<XFEM::Enrichment::EnrType>(enrtype_int);
        
        // check correct reading
        if (position != (int)data.size())
          dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
        
        // create my enriched field
        fieldenr_ = FieldEnr(field,Enrichment(label,enrtype));
        
        return;
      }
      
      //! return global node or element Id
      int getGid() const
      {
        return gid_;
      };
      
      //! return enriched field
      const FieldEnr& getFieldEnr() const
      {
        return fieldenr_;
      };
      
      //! compare: dofkeys are ordered first by the global id and then by their enrichment
      bool operator <(const DofKey& rhs) const
      {
        if (gid_ < rhs.gid_)
          return true;
        else if (gid_ > rhs.gid_)
          return false;
        else
        {
          if (fieldenr_ < rhs.fieldenr_)
            return true;
          else
            return false;
        }
      }
            
      //! test for equality: only equal, if gid and enriched field are identical
      bool operator ==(const DofKey& rhs) const
      {
        if (gid_ == rhs.gid_ and fieldenr_ == rhs.fieldenr_)
          return true;
        else
          return false;
      }
      
      //! test for inequality: if any member differs, the objects are not equal
      bool operator !=(const DofKey& rhs) const
      {
        //        return (!(*this == rhs));
        if (gid_ != rhs.gid_ or fieldenr_ != rhs.fieldenr_)
          return true;
        else
          return false;
      }
              
      friend ostream& operator<<(ostream& out, const DofKey& key)
      {
          out << key.toString();
          return out;
      }
              
    private:
      
      //! don't want default constructor
      explicit DofKey();
      //    :
      //        gid_(-1),
      //        fieldenr_(FieldEnr())
      //    {
      //        dserror("DofKey() -> please don't call me!");
      //        return;
      //    }
      
      //! unique id among the choosen type of dof (for example: nodal global id)
      int             gid_;
      //! specific enriched field
      XFEM::FieldEnr  fieldenr_;
      
  };
  
  //! global Position in Epetra vector 
  typedef int DofPos; 
  
  //! connects each dofkey on Node with a specific position in an solution vector
  typedef std::map<XFEM::DofKey<XFEM::onNode>, XFEM::DofPos > NodalDofPosMap;
  
  //! connects each dofkey on Element with a specific position in an solution vector
  typedef std::map<XFEM::DofKey<XFEM::onElem>, XFEM::DofPos > ElementalDofPosMap;
  
} // namespace XFEM

#endif  // #ifndef DOFKEY_H
#endif  // #ifdef CCADISCRET
