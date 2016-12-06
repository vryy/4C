/*!------------------------------------------------------------------------------------------------*
\file dofkey.cpp

\brief

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
 *------------------------------------------------------------------------------------------------*/


#include "dofkey.H"
#include "enrichment.H"


  XFEM::DofKeyType XFEM::DofKeyType::instance_;


  std::string XFEM::DofKeyType::Name() const
  {
    //return "DofKeyType<doftype>";

    // make the name unique but compiler dependent
    // This is fine since DofKeyType<doftype> is never written to binary io.
    const std::type_info & ti = typeid( *this );
    return ti.name();
  }



  //! standard constructor
  XFEM::DofKey::DofKey(
      const int             gid,      ///< global id
      const FieldEnr&       fieldenr  ///< enriched field
  ) :
    DRT::ParObject(),
    fieldenr_(fieldenr),
    gid_(gid)
    {
        return;
    }

  //! constructor, get values from char array using Unpack member function
  XFEM::DofKey::DofKey(
      std::vector<char>& data      ///< char array with data to unpack
  ) :
    DRT::ParObject()
  {
    Unpack(data);
    return;
  }

  //! copy constructor
  XFEM::DofKey::DofKey(
      const DofKey& other           ///< source
  ) :
    DRT::ParObject(other),
    fieldenr_(other.fieldenr_),
    gid_(other.gid_)

  {
      assert(&other != this);
      return;
  }

  //! return std::string representation of this class
  std::string XFEM::DofKey::toString() const
  {
    std::stringstream s;
    s << "Dofkey: [" << std::setw(3) << gid_ << ", " << fieldenr_.toString() << "]";
    return s.str();
  };

  int XFEM::DofKey::UniqueParObjectId() const
  {
    return XFEM::DofKeyType::Instance().UniqueParObjectId();
  }

  //! pack all content to char array
  void XFEM::DofKey::Pack(DRT::PackBuffer& data) const
  {
    DRT::PackBuffer::SizeMarker sm( data );
    sm.Insert();

    AddtoPack(data,UniqueParObjectId());
    AddtoPack(data,gid_);
    AddtoPack(data,fieldenr_.getField());
    const XFEM::Enrichment::EnrType enrtype = fieldenr_.getEnrichment().Type();
    AddtoPack(data,enrtype);
    const int label = fieldenr_.getEnrichment().XFEMConditionLabel();
    AddtoPack(data,label);
  }

  //! unpack all content from char array
  void XFEM::DofKey::Unpack(const std::vector<char>& data)
  {
    std::vector<char>::size_type position = 0;
    // extract type
    int type = 0;
    ExtractfromPack(position,data,type);
    if (type != UniqueParObjectId()) dserror("wrong instance type data");

    ExtractfromPack(position,data,gid_);

    // read field
    int field;
    ExtractfromPack(position,data,field);

    // read enrichment type
    int enrtype;
    ExtractfromPack(position,data,enrtype);

    // read xfem condition label
    int label;
    ExtractfromPack(position,data,label);

    // check correct reading
    if (position != data.size())
    {
      // this dummy output prevents compilation errors with gcc 5.1 and optimization level -O3
      std::cout << "Detected problem for type " << type << " gid " << gid_ << " field " << field
                << " " << enrtype << " " << label << std::endl;
      dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
    }

    // create my enriched field
    fieldenr_ = FieldEnr(static_cast<PHYSICS::Field>( field ),Enrichment(static_cast<XFEM::Enrichment::EnrType>( enrtype ),label));

    return;
  }

  //! return global node or element Id
  int XFEM::DofKey::getGid() const
  {
    return gid_;
  };

  //! return enriched field
  const XFEM::FieldEnr& XFEM::DofKey::getFieldEnr() const
  {
    return fieldenr_;
  };

  //! compare: dofkeys are ordered first by the global id and then by their enrichment
  bool XFEM::DofKey::operator <(const DofKey& rhs) const
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
  bool XFEM::DofKey::operator ==(const DofKey& rhs) const
  {
    if (gid_ == rhs.gid_ and fieldenr_ == rhs.fieldenr_)
      return true;
    else
      return false;
  }

  //! test for inequality: if any member differs, the objects are not equal
  bool XFEM::DofKey::operator !=(const DofKey& rhs) const
  {
    //        return (!(*this == rhs));
    if (gid_ != rhs.gid_ or fieldenr_ != rhs.fieldenr_)
      return true;
    else
      return false;
  }
