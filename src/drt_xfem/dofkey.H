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

#include "xfem.H"
#include "field_enriched.H"


namespace XFEM
{

/*!
 \brief Allows to distuingish between nodal and elemental dofs
 */
enum DofType
{
    onNode,
    onElement
};

//! create a string from the DofType enum
static std::string dofTypeToString(DofType doftype)
{
    switch (doftype)
    {
        case onNode:    return "onNode";
        case onElement: return "onElement";
        default:
            dserror("unknown doftype");
            exit(1);
    }
}

/*!
 \brief A unique identifier for a degree of freedom (DOF)

 Haskell-esque:
 type FieldEnr = (DofType, Gid, FieldEnr)
 
 \tparam DofType to distinguish between element and nodal dofs, now the compiler complains
         when we mix them, since DofKey<onNode> is different from DofKey<onElement>

 */
template <DofType doftype>
class DofKey
{
public:
    //! empty constructor
    explicit DofKey() :
        gid_(-1),
        fieldenr_(FieldEnr())
    {
        dserror("DofKey() -> please don't call me!");
        return;
    }

    //! standard constructor
    explicit DofKey(
            const int             gid,      ///< global id
            const FieldEnr        fieldenr  ///< enrichmed field
            ) :
                gid_(gid),
                fieldenr_(fieldenr)
    {
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

    
    //! destructor
    ~DofKey()
    {
        return;
    }

    std::string toString() const
    {
        stringstream s;
        s << "Dofkey: " << XFEM::dofTypeToString(doftype) << " with Gid = " << gid_ << ", FieldEnr: " << fieldenr_.toString();
        return s.str();
    };

    int getGid() const
    {
        return gid_;
    };
    
    const FieldEnr& getFieldEnr() const
    {
        return fieldenr_;
    };
    
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
    //!
    bool operator ==(const DofKey& rhs) const
    {
        if (gid_ == rhs.gid_ and fieldenr_ == fieldenr_)
            return true;
        else
            return false;
    }
    
    bool operator !=(const DofKey& rhs) const
    {
//        return (!(*this == rhs));
        if (gid_ != rhs.gid_ or fieldenr_ != fieldenr_)
            return true;
        else
            return false;
    }
    

private:
    //! unique id among the choosen type of dof
    const int       gid_;
    //! specific enrichment field
    const FieldEnr  fieldenr_;

}; // class DofKey

} // namespace XFEM

#endif  // #ifndef DOFKEY_H
#endif  // #ifdef CCADISCRET
