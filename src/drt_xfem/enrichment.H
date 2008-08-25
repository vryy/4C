/*!
 * \file enrichment.H
 * 
 * \brief describes the enrichment class
 *
 * <pre>
 * Maintainer: Axel Gerstenberger
 *             gerstenberger@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15236
 * </pre>
 */
#ifdef CCADISCRET
#ifndef ENRICHMENT_H
#define ENRICHMENT_H


#include "../drt_geometry/vector_definitions.H"
#include "interface.H"

//! namespace for all XFEM related stuff
namespace XFEM
{

  /*!
   * \brief Corresponds to a specific enrichment with unique id and type
   * 
   * Implements a tuple of (iface_label, enrtype) and is sorted in that order
   */
  class Enrichment
  {
    public:
      
      /*! Various kinds of enrichments, which differ in the way,
       *  the enrichment function is evaluated
       */
      enum EnrType
      {
        typeUndefined,   ///< a default value, for which all functions should fail/crash
        typeStandard,    ///< standard enrichment  - corresponds to the standard finite element method
        typeJump,        ///< jump enrichment - Heaviside function from +1 (in Omega+) to -1 (in Omega-)
        typeVoid         ///< void enrichment - Heaviside function from +1 (in Omega+) to  0 (in Omega-)
      };
      
      /*!
       * we want to be able to say, wether we approach froma particular side towards the discontinuity
       * this makes sense for jump or kink enrichments,
       * boundary layers or similar things might not need this
       */
      enum ApproachFrom
      {
        approachUnknown,   ///<  approaching direction is unknown
        approachFromPlus,  ///<  we are coming from Omega+ \note there is a check needed, whether we are in the fluid ore not!!! thin structures will cause problems here!!!!!!
        approachFromMinus  ///<  we are coming from Omega-
      };
      
      //! default constructor
      explicit Enrichment(
      ) :
        xfemconditionlabel_(-1), // means undefined
        type_(typeUndefined)
      {
      //dserror("Enrichment() -> please don't call me!");
      return;
      }
      
      //! constructor
      explicit Enrichment(
          const int      xfemconditionlabel,   ///< unique label identifying the xfem interface
          const EnrType  type                  ///< enrichment type
      ) :
        xfemconditionlabel_(xfemconditionlabel),
        type_(type)
      {
        return;
      }
          
      //! copy constructor
      Enrichment(
          const Enrichment& other
      ) :
        xfemconditionlabel_(other.xfemconditionlabel_),
        type_(other.type_)
      {
        assert(&other != this);
        return;
      }
      
      //! return string representation of the Enrichment
      std::string toString() const;
      
      //! return string representation of an EnrType
      std::string enrTypeToString(const EnrType type) const;
      
      //! return label of the condition that this enrichment belongs to
      int XFEMConditionLabel() const
      {
        return xfemconditionlabel_;
      };
              
      //! return type of enrichment
      EnrType Type() const
      {
        return type_;
      };
      
      //! comparison operator
      bool operator ==(const Enrichment& rhs) const
      {
        if ((xfemconditionlabel_ == rhs.xfemconditionlabel_) and (type_ == rhs.type_))
          return true;
        else
          return false;
      };
      
      //! (redundant?) not equal operator
      bool operator !=(const Enrichment& rhs) const
      {
        //return (!(*this == rhs));
        if ((xfemconditionlabel_ != rhs.xfemconditionlabel_) or (type_ != rhs.type_))
          return true;
        else
          return false;
      }
              
      /*! less than operator
       *  Enrichments are ordered first by their condition label, then by type
       */
      bool operator <(const Enrichment& rhs) const
      {
        if (xfemconditionlabel_ < rhs.xfemconditionlabel_)
          return true;
        else if (xfemconditionlabel_ > rhs.xfemconditionlabel_)
          return false;
        else
        {
          if (type_ < rhs.type_)
            return true;
          else
            return false;
        }
      }
      
      //! calculate modified enrichment value based on input \f$\psi(x) - \psi(x_{\mathrm{n}})\f$
      double ModifiedEnrValue(
          const BlitzVec3&                      actpos,            ///< \f$x\f$
          const BlitzVec3&                      nodalpos,          ///< \f$x_{\mathrm{Node}}\f$
          const XFEM::InterfaceHandle&          ih,                ///< the ubiquitious interface handle with all geometry information
          const XFEM::Enrichment::ApproachFrom  approachdirection  ///< approach direction
      ) const;
      
      //! calculate enrichment value based on input \f$\psi(x)\f$
      double EnrValue(
          const BlitzVec3&                      actpos,            ///< \f$x\f$
          const XFEM::InterfaceHandle&          ih,                ///< the ubiquitious interface handle with all geometry information
          const XFEM::Enrichment::ApproachFrom  approachdirection  ///< approach direction
      ) const;
              
    private:
      //! a global id used to distuinguish several surfaces within one element
      const int xfemconditionlabel_;
      //! the specific type of enrichment
      const EnrType type_;
  };
  
}


#endif  // #ifndef ENRICHMENT_H
#endif  // #ifdef CCADISCRET
