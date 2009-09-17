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


#include "interfacexfsi.H"
#include "../drt_combust/combust_interface.H"

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
        typeUndefined, ///< a default value, for which all functions should fail/crash
        typeStandard,  ///< standard enrichment  - corresponds to the standard finite element method
        typeJump,      ///< jump enrichment - Heaviside function from +1 (in Omega+) to -1 (in Omega-)
        typeVoid,      ///< void enrichment - Heaviside function from +1 (in Omega+) to  0 (in Omega-)
        typeKink       ///< kink enrichment - TODO: discription missing
      };

      /*!
       * we want to be able to say, whether we approach from a particular side towards the discontinuity
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

      //!assignment operator
      Enrichment& operator = (const Enrichment& old) ;

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

      //! calculate enrichment value based on input \f$\psi(x)\f$
      double EnrValue(
          const LINALG::Matrix<3,1>&            actpos,            ///< \f$x\f$
          const XFEM::InterfaceHandle&          ih,                ///< the ubiquitous interface handle with all geometry information
          const XFEM::Enrichment::ApproachFrom  approachdirection  ///< approach direction
      ) const;

      //! compute enrichment value for integration cell (constant!)
      double EnrValueIntCell(
          const GEO::DomainIntCell& cell
      ) const;

      //! compute enrichment value at an embedded interface
      double EnrValueAtInterface(
          const XFEM::Enrichment::ApproachFrom  approachdirection
      ) const;

      /*!
       * \brief get enrichment value at node using the level set function
       *
       * \author henke 07/09
       */
      template<class V1>
      double EnrValueNodeLevelSet(
          const int  inode,
          const V1&  phi
      ) const
      {
          //std::cout << "EnrValueNodeLevelSet()" << std::endl;
          double enrval = 1.0;

          switch (Type())
          {
          case XFEM::Enrichment::typeStandard:
          {
              enrval = 1.0;
              break;
          }
          case XFEM::Enrichment::typeVoid:
          {
              if (phi(inode) < 0.0) {
                  enrval = 0.0;
              } else {
                  enrval = 1.0; // phi == 0 is enriched!
              }
              break;
          }
          case XFEM::Enrichment::typeJump:
          {
              if (phi(inode) < 0.0) {
                  enrval = -1.0;
              } else {
                  enrval = 1.0; // phi == 0 is enriched!
              }
              break;
          }
          case XFEM::Enrichment::typeKink:
          {
              enrval = 0.0;
              dserror("(modified) kink enrichments are zero at all nodes!");
          }
          default:
              dserror("unsupported type of enrichment!");
          }
          return enrval;
      }


      /*!
       * \brief get enrichment value at (Gauss) point using the level set function
       *
       * \author henke 07/09
       */
      template<class V1>
      double EnrValueGaussPointLevelSet(
          const DRT::Element&  ele,
          const V1&            phi,
          const V1&            funct
      ) const
      {
          const size_t numnode = ele.NumNode();

          double enrval = 1.0;

          switch (Type())
          {
          case XFEM::Enrichment::typeStandard:
          {
              enrval = 1.0;
              break;
          }
          case XFEM::Enrichment::typeVoid:
          {
              double phival = 0.0;
              // evaluate phi using shape function values at given point
              for (size_t inode=0; inode<numnode; inode++)
              {
                phival += phi(inode) * funct(inode);
              }
              if (phival < 0.0) {
                  enrval = 0.0;
              } else {
                  enrval = 1.0; // phi == 0 is enriched!
              }
              break;
          }
          case XFEM::Enrichment::typeJump:
          {
              double phival = 0.0;
              // evaluate phi using shape function values at given point
              for (size_t inode=0; inode<numnode; inode++)
              {
                phival += phi(inode) * funct(inode);
              }
              if (phival < 0.0) {
                  enrval = -1.0;
              } else {
                  enrval = 1.0; // phi == 0 is enriched!
              }
              break;
          }
          case XFEM::Enrichment::typeKink:
          {
          // literature:
          // Moës, N., Cloirec, M.,Cartraud, P. and Remacle, J. F.
          // A computational approach to handle complex microstructure geometries:
          // "Computer Methods in Applied Mechanics and Engineering", 192:3163--3177, 2003.
          //
          // psi = sum(abs(phi(i))*N(i))-abs(sum(phi(i)*N(i)), i=0...numnode-1

              double sum1 = 0.0;
              double sum2 = 0.0;
              for (size_t inode=0; inode<numnode; inode++)
              {
                sum1 += fabs(phi(inode)) * funct(inode); //absolute value in C++: fabs() and include math.h
                sum2 += phi(inode) * funct(inode);
              }
              enrval = sum1 - fabs(sum2);
              break;
          }
          default:
              dserror("unsupported enrichment (modified)!");
          }
          return enrval;
      }


      /*!
       * \brief get derivative of enrichment value at (Gauss) point
       *
       * \autor rasthofer 06/09
       */
      template<class V1, class M1>
      std::vector<double> EnrValueGaussPointDerivative(
              const DRT::Element&   ele,
              const V1&             phi,
              const V1&             funct,
              const M1&             derxy
      ) const
      {
          const size_t numnode = ele.NumNode();

          std::vector<double> enrvalderxy(3);
          enrvalderxy[0] = 0.0;
          enrvalderxy[1] = 0.0;
          enrvalderxy[2] = 0.0;

          switch (Type()){
          case XFEM::Enrichment::typeStandard:
          {
              break;
          }
          case XFEM::Enrichment::typeVoid:
          {
              break;
          }
          case XFEM::Enrichment::typeJump:
          {
              break;
          }
          case XFEM::Enrichment::typeKink:
          {
              for (size_t isd=0; isd<3; isd++)
              {
                  double sum1der = 0.0;
                  double sum2der = 0.0;
                  double sum2 = 0.0;

                  for (size_t inode=0; inode<numnode; inode++)
                  {
                      sum1der += fabs(phi(inode)) * derxy(isd,inode);
                      sum2der += phi(inode) * derxy(isd,inode);
                      sum2 += phi(inode) * funct(inode);
                  }

                  if (sum2 < 0)
                  {
                      enrvalderxy[isd] = sum1der + sum2der;
                  }
                  else
                  {
                      enrvalderxy[isd] = sum1der - sum2der;
                  }
              }
              break;
          }
          default:
              dserror("unsupported type of enrichment!");
          }

          return enrvalderxy;
      }

      /*!
       * \brief get second derivative of enrichment value at (Gauss) point
       *
       * \author rasthofer 06/09
       */
      template<class V1, class M2>
      std::vector<double> EnrValueGaussPointDerivative2(
              const DRT::Element&                    ele,
              const V1&                              phi,
              const V1&                              funct,
              const M2&                              derxy2
      ) const
      {
          const size_t numnode = ele.NumNode();

          std::vector<double> enrvalderxy2(6);
          for (int isd=0; isd<6; isd++)
          {
              enrvalderxy2[isd] = 0.0;
          }

          switch (Type()){
          case XFEM::Enrichment::typeStandard:
          {
              break;
          }
          case XFEM::Enrichment::typeVoid:
          {
              break;
          }
          case XFEM::Enrichment::typeJump:
          {
              break;
          }
          case XFEM::Enrichment::typeKink:
          {
              for (size_t isd=0; isd<6; isd++)
              {
                  double sum1der2 = 0.0;
                  double sum2der2 = 0.0;
                  double sum2 = 0.0;

                  for (size_t inode=0; inode<numnode; inode++)
                  {
                      sum1der2 += fabs(phi(inode)) * derxy2(isd,inode);
                      sum2der2 += phi(inode) * derxy2(isd,inode);
                      sum2 += phi(inode) * funct(inode);
                  }

                  if (sum2 < 0)
                  {
                      enrvalderxy2[isd] = sum1der2 + sum2der2;
                  }
                  else
                  {
                      enrvalderxy2[isd] = sum1der2 - sum2der2;
                  }
              }
              break;
          }
          default:
              dserror("unsupported enrichment (modified)!");
          }

        return enrvalderxy2;
      }


      //! calculate modified enrichment value based on input \f$\psi(x) - \psi(x_{\mathrm{n}})\f$
      double ModifiedEnrValue(
          const LINALG::Matrix<3,1>&            actpos,            ///< \f$x\f$
          const LINALG::Matrix<3,1>&            nodalpos,          ///< \f$x_{\mathrm{Node}}\f$
          const XFEM::InterfaceHandle&          ih,                ///< the ubiquitous interface handle with all geometry information
          const XFEM::Enrichment::ApproachFrom  approachdirection  ///< approach direction
      ) const;

    private:

      //! a global id used to distuinguish several surfaces within one element
      int xfemconditionlabel_;
      //! the specific type of enrichment
      EnrType type_;
  };

}


#endif  // #ifndef ENRICHMENT_H
#endif  // #ifdef CCADISCRET
