/*!
\file enrichment_utils.H

\brief describes the enrichment types and classes

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
*/


#ifndef ENRICHMENT_UTILS_H
#define ENRICHMENT_UTILS_H


#include "dof_management_element.H"
#include "field_enriched.H"
#include "xfem_element_utils.H"
#include "../drt_geometry/integrationcell.H"
#include "../drt_combust/combust_defines.H"


namespace COMBUST
{
  class InterfaceHandleCombust;
}

namespace LINALG
{
  class SerialDenseVector;
}

namespace XFEM
{
  /*!
   * \brief pre-calculate enrichment values for an element and store them in this class
   */
  class ElementEnrichmentValues
  {

  public:

    /*!
    * \brief for each given Enrichment in enrset, compute the enrichment value at actpos
    *
    * if one already knows, on which side of the discontinuity one computes something,
    * the information can be given via approachdirection
    */
    explicit ElementEnrichmentValues(
        const DRT::Element&                   ele,
        COMBUST::InterfaceHandleCombust*      ih,               ///< interface information
        const XFEM::ElementDofManager&        dofman,
        const LINALG::Matrix<3,1>&            actpos,           ///< a 3d position
        const bool                            boundary_integral,
        const int                             boundary_label
    );

    /*!
     * \brief compute enrichment value for domain integration cell based on level set function
     *
     * - enrichment value is constant within a domain integration cell (for jump enrichments)
     * - constructor called for combustion problems (jump enrichment, level set)
     *
     * \author henke 06/09
     */
    template<class V1>
    explicit ElementEnrichmentValues(
            const DRT::Element&             ele,    ///< pointer to element
            const XFEM::ElementDofManager&  dofman, ///< pointer to element dof manager
            const GEO::IntCell&             cell,   ///< pointer to integration cell
            const V1&                       phi     ///< vector of level set values
     ) :
       ele_(ele),
       dofman_(dofman)
    {
      enrvals_.clear();
      enrvalnodes_.clear();
      enrvalderxy_.clear();  // not used for jump enrichment (Heaviside function)
      enrvalderxy2_.clear(); // not used for jump enrichment (Heaviside function)

      const std::set<XFEM::Enrichment>& enrset(dofman.getUniqueEnrichments());
      const int numnode = ele_.NumNode();

      //--------------------------------------------------
      // prepare temporary variables for enrichment values
      //--------------------------------------------------
      // enrichment values at nodes
      std::map<int,std::map<XFEM::Enrichment, double> > enrvalnodes;
      // first derivatives of enrichment function
      std::vector<double> enrvalderxy(3);
      enrvalderxy[0] = 0.0;
      enrvalderxy[1] = 0.0; // zero for Heaviside function
      enrvalderxy[2] = 0.0;
      // second derivatives of enrichment function
      std::vector<double> enrvalderxy2(6);
      for (int isd=0; isd<6; isd++)
      {
        enrvalderxy2[isd] = 0.0; // zero for Heaviside function
      }

      //------------------------------------------------------
      // evaluate enrichment values for this integration cell
      //------------------------------------------------------
      for (std::set<XFEM::Enrichment>::const_iterator enriter = enrset.begin(); enriter != enrset.end(); ++enriter)
      {
        const double enrvalcell = enriter->EnrValueIntCell(cell);
        enrvals_[*enriter] = enrvalcell;
        enrvalderxy_[*enriter] = enrvalderxy;
        enrvalderxy2_[*enriter] = enrvalderxy2;
      }

      //------------------------------------------------------------
      // loop over nodes to evaluate enrichment values at the nodes
      //------------------------------------------------------------
      for(int inode=0; inode < numnode; inode++)
      {
        std::map<XFEM::Enrichment, double> enrvalstmp;
        for (std::set<XFEM::Enrichment>::const_iterator enriter = enrset.begin(); enriter != enrset.end(); ++enriter)
        {
          double enrvalnode = 0.0;
          if (enriter->Type()==XFEM::Enrichment::typeJump || enriter->Type()==XFEM::Enrichment::typeVoid)
          {
            enrvalnode = enriter->EnrValueNodeLevelSet(inode,phi);
          }
          enrvalstmp[*enriter] = enrvalnode;
        }
        enrvalnodes[inode] = enrvalstmp;
        enrvalstmp.clear();
      }

      //------------------------------------------------------------------------
      // evaluate final enrichment values for each node H(G(cell)) - H(G(node))
      //------------------------------------------------------------------------
      // literature (p. 1006, penultimate line):
      // Belytschko, T., Moës, N., Usui, S. and Parimi, C.
      // Arbitrary discontinuities in finite elements:
      // "International Journal for Numerical Methods in Engineering", 50:993--1013,2001.

      for(int inode=0; inode < numnode; inode++)
      {
        std::map<XFEM::Enrichment, double> enrvalues;
        for (std::set<XFEM::Enrichment>::const_iterator enriter = enrset.begin(); enriter != enrset.end(); ++enriter)
        {
          enrvalues[*enriter] = enrvals_.find(*enriter)->second - (enrvalnodes.find(inode)->second).find(*enriter)->second;
        }
        enrvalnodes_[inode] = enrvalues;
        enrvalues.clear();
      }

      return;
    }

    /*!
     * \brief compute enrichment value at one side of an embedded interface based on level set function
     *
     * - enrichment value is constant on both sides of the interface (for jump enrichments)
     * - constructor called for evaluation of boundary integrals in combustion problems (jump enrichment, level set)
     *
     * \author henke 09/09
     */
    template<class V1>
    explicit ElementEnrichmentValues(
            const DRT::Element&                   ele,               ///< pointer to element
            const XFEM::ElementDofManager&        dofman,            ///< pointer to element dof manager
            const XFEM::Enrichment::ApproachFrom  approachdirection, ///< side of embedded interface
            const V1&                             phi                ///< vector of level set values
     ) :
       ele_(ele),
       dofman_(dofman)
    {
      enrvals_.clear();
      enrvalnodes_.clear();
      enrvalderxy_.clear();  // 0 for jump enrichment
      enrvalderxy2_.clear(); // 0 for jump enrichment

      const std::set<XFEM::Enrichment>& enrset(dofman.getUniqueEnrichments());
      const int numnode = ele_.NumNode();

      //--------------------------------------------------
      // prepare temporary variables for enrichment values
      //--------------------------------------------------
      // enrichment values at nodes
      std::map<int,std::map<XFEM::Enrichment, double> > enrvalnodes;
      // first derivatives of enrichment function
      std::vector<double> enrvalderxy(3);
      enrvalderxy[0] = 0.0;
      enrvalderxy[1] = 0.0; // zero for Heaviside function
      enrvalderxy[2] = 0.0;
      // second derivatives of enrichment function
      std::vector<double> enrvalderxy2(6);
      for (int isd=0; isd<6; isd++)
      {
        enrvalderxy2[isd] = 0.0; // zero for Heaviside function
      }

      //------------------------------------------------------
      // evaluate enrichment values for this integration cell
      //------------------------------------------------------
      for (std::set<XFEM::Enrichment>::const_iterator enriter = enrset.begin(); enriter != enrset.end(); ++enriter)
      {
        const double enrvalside = enriter->EnrValueAtInterface(approachdirection);
        enrvals_[*enriter] = enrvalside;
        enrvalderxy_[*enriter] = enrvalderxy;
        enrvalderxy2_[*enriter] = enrvalderxy2;
      }

      //------------------------------------------------------------
      // loop over nodes to evaluate enrichment values at the nodes
      //------------------------------------------------------------
      for(int inode=0; inode < numnode; inode++)
      {
        std::map<XFEM::Enrichment, double> enrvalstmp;
        for (std::set<XFEM::Enrichment>::const_iterator enriter = enrset.begin(); enriter != enrset.end(); ++enriter)
        {
          double enrvalnode = 0.0;
          if (enriter->Type()==XFEM::Enrichment::typeJump || enriter->Type()==XFEM::Enrichment::typeVoid)
          {
            enrvalnode = enriter->EnrValueNodeLevelSet(inode,phi);
          }
          enrvalstmp[*enriter] = enrvalnode;
        }
        enrvalnodes[inode] = enrvalstmp;
        enrvalstmp.clear();
      }

      //------------------------------------------------------------------------
      // evaluate final enrichment values for each node H(G(cell)) - H(G(node))
      //------------------------------------------------------------------------
      // literature (p. 1006, penultimate line):
      // Belytschko, T., Moës, N., Usui, S. and Parimi, C.
      // Arbitrary discontinuities in finite elements:
      // "International Journal for Numerical Methods in Engineering", 50:993--1013,2001.

      for(int inode=0; inode < numnode; inode++)
      {
        std::map<XFEM::Enrichment, double> enrvalstmp;
        for (std::set<XFEM::Enrichment>::const_iterator enriter = enrset.begin(); enriter != enrset.end(); ++enriter)
        {
          enrvalstmp[*enriter] = enrvals_.find(*enriter)->second - (enrvalnodes.find(inode)->second).find(*enriter)->second;
          //std::cout << "shifted enrichment value " << enrvalstmp[*enriter] << std::endl;
        }
        enrvalnodes_[inode] = enrvalstmp;
        enrvalstmp.clear();
      }

      return;
    }

    /*!
     * \brief compute enrichment value for domain integration cell
     *
     * constructor called for two phase flow with jump enrichments in pressure and kink enr. in velocity
     * jumps enrichments are called via DomainIntCell, not with level-set values!!!
     * kinks are called using level-set values, called in combust3_sysmat_twophaseflow.H
     *                                                                                                     schott 09/10/10
     */
    template<class V1, class M1, class M2>
    explicit ElementEnrichmentValues(
        const DRT::Element&             ele,
        const XFEM::ElementDofManager&  dofman,
        const V1&                       phi,
        const GEO::DomainIntCell&       cell,
        const V1&                       funct,   ///< vector of shape functions evaluated at Gauss point
        const M1&                       derxy,
        const M2&                       derxy2
    ) :
    ele_(ele),
    dofman_(dofman)
    {
      enrvals_.clear(); // not used for kink enrichment
      enrvalnodes_.clear();
      enrvalderxy_.clear();
      enrvalderxy2_.clear();

      const std::set<XFEM::Enrichment>& enrset(dofman.getUniqueEnrichments());
      const int numnode = ele_.NumNode();

      for(int inode=0; inode<numnode; inode++)
      {
        std::map<XFEM::Enrichment, double> enrvalues;
        for (std::set<XFEM::Enrichment>::const_iterator enriter = enrset.begin(); enriter != enrset.end(); ++enriter)
        {
          double enrval = 0.0;
          if(enriter->Type()==XFEM::Enrichment::typeJump || enriter->Type()==XFEM::Enrichment::typeVoid
              || enriter->Type() == XFEM::Enrichment::typeStandard)
            enrval = enriter->EnrValueIntCell(cell);
          else if(enriter->Type()==XFEM::Enrichment::typeKink)
            enrval = enriter->EnrValueGaussPointLevelSet(ele_,phi,funct);
          else dserror("no known enrichment type");

          double enrvalnode = 0.0;
          // modified enrichment in the case of jump enrichment to maintain the interpolant properties
          // Psi(xGauss)-Psi(xNode)
          if (enriter->Type()==XFEM::Enrichment::typeJump || enriter->Type()==XFEM::Enrichment::typeVoid)
          {
#ifndef COMBUST_SXFEM
            enrvalnode = enriter->EnrValueNodeLevelSet(inode,phi);
#else
            for(int jnode=0; jnode<numnode; jnode++)
            {
              enrvalnode += funct(jnode) * enriter->EnrValueNodeLevelSet(jnode,phi);
            }
#endif
          }
          enrvalues[*enriter] = enrval - enrvalnode;
        }
        enrvalnodes_[inode] = enrvalues;
        enrvalues.clear();
      }

      for (std::set<XFEM::Enrichment>::const_iterator enriter = enrset.begin(); enriter != enrset.end(); ++enriter)
      {
        const std::vector<double> enrvalderxy = enriter->EnrValueGaussPointDerivative(ele_,phi,funct,derxy);
        enrvalderxy_[*enriter] = enrvalderxy;

        const std::vector<double> enrvalderxy2 = enriter->EnrValueGaussPointDerivative2(ele_,phi,funct,derxy2);
        enrvalderxy2_[*enriter] = enrvalderxy2;
      }

      return;
    }


    /*!
     * \brief compute enrichment value at one side of an embedded interface for boundary integration cell
     *
     * constructor called for evaluation of boundary integrals in two-phase-flow problems (jump and kink !!!! enrichment)
     * used for two phase flow with surface tension with kink enrichments in velocity and jump enrichments in pressure
     * call jump enrichments with approachdirection, not with level-set
     * but kink enrichments need the level-set field
     *
     * schott May 12, 2010
     */
    template<class V1, class M1, class M2>
    explicit ElementEnrichmentValues(
        const DRT::Element&                   ele,               ///< pointer to element
        const XFEM::ElementDofManager&        dofman,            ///< pointer to element dof manager
        const XFEM::Enrichment::ApproachFrom  approachdirection, ///< side of embedded interface
        const V1&                             phi,                ///< vector of level set values
        const V1&                             funct,   ///< vector of shape functions evaluated at Gauss point
        const M1&                             derxy,   ///< Matrix of derivatives of shape functions evaluated at Gauss point
        const M2&                             derxy2   ///< second derivatives
    ) :
    ele_(ele),
    dofman_(dofman)
    {
      enrvals_.clear();      // tmp variable for enrvalnodes_
      enrvalnodes_.clear();
      enrvalderxy_.clear();  // 0 for jump enrichment, not 0 for kink enrichment
      enrvalderxy2_.clear(); // 0 for jump enrichment, not 0 for kink enrichment

      const std::set<XFEM::Enrichment>& enrset(dofman.getUniqueEnrichments());
      const int numnode = ele_.NumNode();
      std::map<int,std::map<XFEM::Enrichment, double> > enrvalnodes;

      //------------------------------------------------------
      // evaluate enrichment values for this integration cell
      //------------------------------------------------------
      for (std::set<XFEM::Enrichment>::const_iterator enriter = enrset.begin(); enriter != enrset.end(); ++enriter)
      {
        // get enrichment values all types of enrichment
        double enrvalside = 0.0;

        std::vector<double> enrvalderxy(3);
        enrvalderxy[0] = 0.0;
        enrvalderxy[1] = 0.0;
        enrvalderxy[2] = 0.0;
        std::vector<double> enrvalderxy2(6);
        for (int isd=0; isd<6; isd++)
        {
          enrvalderxy2[isd] = 0.0;
        }

        // kink enrichments
        if(enriter->Type() == XFEM::Enrichment::typeKink)
        {
          enrvalside   = enriter->EnrValueGaussPointLevelSet(ele_,phi,funct);
          enrvalderxy  = enriter->EnrValueGaussPointDerivative(ele_,phi,funct,derxy);
          enrvalderxy2 = enriter->EnrValueGaussPointDerivative2(ele_,phi,funct,derxy2);
        }
        // everything but kink enrichments (jump, void)
        else
        {
          enrvalside = enriter->EnrValueAtInterface(approachdirection);
          // no derivatives needed for jump or void enrichments
#ifdef COMBUST_SXFEM
          enrvalderxy  = enriter->EnrValueGaussPointDerivative(ele_,phi,funct,derxy);
          enrvalderxy2 = enriter->EnrValueGaussPointDerivative2(ele_,phi,funct,derxy2);
#endif
        }

        enrvals_[*enriter] = enrvalside;
        enrvalderxy_[*enriter] = enrvalderxy;
        enrvalderxy2_[*enriter] = enrvalderxy2;
      }

      //------------------------------------------------------------
      // loop over nodes to evaluate enrichment values at the nodes
      //------------------------------------------------------------
      for(int inode=0; inode < numnode; inode++)
      {
        std::map<XFEM::Enrichment, double> enrvalstmp;
        for (std::set<XFEM::Enrichment>::const_iterator enriter = enrset.begin(); enriter != enrset.end(); ++enriter)
        {
          double enrvalnode = 0.0;
          if (enriter->Type()==XFEM::Enrichment::typeJump || enriter->Type()==XFEM::Enrichment::typeVoid)
          {
#ifndef COMBUST_SXFEM
            enrvalnode = enriter->EnrValueNodeLevelSet(inode,phi);
#else
            for(int jnode=0; jnode<numnode; jnode++)
            {
              enrvalnode += funct(jnode) * enriter->EnrValueNodeLevelSet(jnode,phi);
            }
#endif
          }
          // enrvalnode = 0.0 for kink enrichment
          enrvalstmp[*enriter] = enrvalnode;
        }
        enrvalnodes[inode] = enrvalstmp;
        enrvalstmp.clear();
      }

      //------------------------------------------------------------------------
      // evaluate final enrichment values for each node H(G(cell)) - H(G(node))
      //------------------------------------------------------------------------
      // literature (p. 1006, penultimate line):
      // Belytschko, T., Moës, N., Usui, S. and Parimi, C.
      // Arbitrary discontinuities in finite elements:
      // "International Journal for Numerical Methods in Engineering", 50:993--1013,2001.

      for(int inode=0; inode < numnode; inode++)
      {
        std::map<XFEM::Enrichment, double> enrvalstmp;
        for (std::set<XFEM::Enrichment>::const_iterator enriter = enrset.begin(); enriter != enrset.end(); ++enriter)
        {
          enrvalstmp[*enriter] = enrvals_.find(*enriter)->second - (enrvalnodes.find(inode)->second).find(*enriter)->second;
        }
        enrvalnodes_[inode] = enrvalstmp;
        enrvalstmp.clear();
      }

      return;
    }


    /*!
    * \brief for each given Enrichment in enrset, compute the enrichment value at actpos
    *
    * constructor called for two-phase flow problems (kink enrichment)
    *
    * \author rasthofer 06/09
    */
    template<class V1, class M1, class M2>
    explicit ElementEnrichmentValues(
            const DRT::Element&             ele,
            const XFEM::ElementDofManager&  dofman,
            const V1&                       phi,
            const V1&                       funct,   ///< vector of shape functions evaluated at Gauss point
            const M1&                       derxy,
            const M2&                       derxy2
            ) :
              ele_(ele),
              dofman_(dofman)
    {
      enrvals_.clear(); // not used for kink enrichment
      enrvalnodes_.clear();
      enrvalderxy_.clear();
      enrvalderxy2_.clear();

      const std::set<XFEM::Enrichment>& enrset(dofman.getUniqueEnrichments());
      const int numnode = ele_.NumNode();

      for(int inode=0; inode<numnode; inode++)
      {
        std::map<XFEM::Enrichment, double> enrvalues;
        for (std::set<XFEM::Enrichment>::const_iterator enriter = enrset.begin(); enriter != enrset.end(); ++enriter)
        {
          const double enrval = enriter->EnrValueGaussPointLevelSet(ele_,phi,funct);
          double enrvalnode = 0.0;
          // modified enrichment in the case of jump enrichment to maintain the interpolant properties
          // Psi(xGauss)-Psi(xNode)
          if (enriter->Type()==XFEM::Enrichment::typeJump || enriter->Type()==XFEM::Enrichment::typeVoid)
          {
            enrvalnode = enriter->EnrValueNodeLevelSet(inode,phi);
          }
          enrvalues[*enriter] = enrval - enrvalnode;
        }
        enrvalnodes_[inode] = enrvalues;
        enrvalues.clear();
      }

      for (std::set<XFEM::Enrichment>::const_iterator enriter = enrset.begin(); enriter != enrset.end(); ++enriter)
      {
        const std::vector<double> enrvalderxy = enriter->EnrValueGaussPointDerivative(ele_,phi,funct,derxy);
        enrvalderxy_[*enriter] = enrvalderxy;

        const std::vector<double> enrvalderxy2 = enriter->EnrValueGaussPointDerivative2(ele_,phi,funct,derxy2);
        enrvalderxy2_[*enriter] = enrvalderxy2;
      }

      return;
    }


    /*!
    * \brief compute the enriched shape functions
    *
    * remarks:
    * - simplest case: jump or void enrichment
    *   -> no chain rule, since derivative of enrichment function is zero
    * - be carefull: if enrichment function has derivatives not equal to zero, we need the chain rule here!
    *   -> for kink enrichment use ComputeKinkEnrichedNodalShapefunction()
    */
    template<class VI, class VO>
    inline void ComputeEnrichedNodalShapefunction(
        const XFEM::PHYSICS::Field                 field,
        const VI&                                  funct,
        VO&                                        enr_funct
    ) const
    {
      const int* nodeids = ele_.NodeIds();
      const std::size_t numnode = ele_.NumNode();

      std::size_t dofcounter = 0;
      for (std::size_t inode=0; inode < numnode; ++inode)
      {
        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman_.FieldEnrSetPerNode(gid);
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
          enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
          if (enrfield->getField() == field)
          {
            const double enrval = enrvals_.find(enrfield->getEnrichment())->second;
//std::cout << "enrichment value: " << enrval << std::endl;
            enr_funct(dofcounter) = funct(inode) * enrval;
            dofcounter += 1;
          }
        }
      }
      dsassert(dofcounter == dofman_.NumDofPerField(field), "mismatch in information from eledofmanager!");
    }

    /*!
    * \brief compute the modified enriched shape functions
    *
    * derivatives are not computed
    */
    template<class VI, class VO>
    inline void ComputeModifiedEnrichedNodalShapefunction(
        const XFEM::PHYSICS::Field                 field,
        const VI&                                  funct,
        VO&                                        enr_funct
    ) const
    {
      const int* nodeids = ele_.NodeIds();
      const std::size_t numnode = ele_.NumNode();

      std::size_t dofcounter = 0;
      for (std::size_t inode=0; inode < numnode; ++inode)
      {
        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman_.FieldEnrSetPerNode(gid);
        const std::map<XFEM::Enrichment, double> enrvalues = enrvalnodes_.find(inode)->second;

        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
          enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
          if (enrfield->getField() == field)
          {
            const double enrval = enrvalues.find(enrfield->getEnrichment())->second;
            enr_funct(dofcounter) = funct(inode) * enrval;
            dofcounter += 1;
          }
        }
      }
      dsassert(dofcounter == dofman_.NumDofPerField(field), "mismatch in information from eledofmanager!");
    }


    /*!
    * \brief compute the enriched shape functions and the first derivative
    *
    * remarks:
    * - these are not the modified, or shifted, enriched shape functions
    * - the chain rule is not applied here for the first derivative
    * - this is correct for jump or void enrichments, since the derivative of the enrichment
    *   function is zero
    *
    * /!\ be carefull: if enrichment function has derivatives not equal to zero, we need the chain rule here!
    *     -> for kink enrichment use ComputeModifiedEnrichedNodalShapefunction(...,enr_derxy,enr_derxy2)
    */
    template<class VI, class MI, class VO, class MO>
    inline void ComputeEnrichedNodalShapefunction(
        const XFEM::PHYSICS::Field field,
        const VI& funct,
        const MI& derxy,
        VO& enr_funct,
        MO& enr_derxy
    ) const
    {
      const int* nodeids = ele_.NodeIds();
      const std::size_t numnode = ele_.NumNode();

      std::size_t dofcounter = 0;
      for (std::size_t inode=0; inode < numnode; ++inode)
      {
        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman_.FieldEnrSetPerNode(gid);
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
          enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
          if (enrfield->getField() == field)
          {
            const double enrval = enrvals_.find(enrfield->getEnrichment())->second;
            enr_funct(dofcounter) = funct(inode) * enrval;
            for (std::size_t isd = 0; isd < 3; ++isd)
            {
              enr_derxy(isd,dofcounter) = derxy(isd,inode) * enrval;
            }
            dofcounter += 1;
          }
        }
      }
      dsassert(dofcounter == dofman_.NumDofPerField(field), "mismatch in information from eledofmanager!");
    }

    /*!
    * \brief compute the enriched shape functions and the derivatives
    *
    * remarks:
    * - these are not the modified, or shifted, enriched shape functions
    * - the chain rule is not applied here for the first derivative
    * - this is correct for jump or void enrichments, since the derivative of the enrichment
    *   function is zero
    *
    * /!\ be carefull: if enrichment function has derivatives not equal to zero, we need the chain rule here!
    *     -> for kink enrichment use ComputeModifiedEnrichedNodalShapefunction(...,enr_derxy,enr_derxy2)
    */
    template<class VI, class MI1, class MI2, class VO, class MO1, class MO2 >
    inline void ComputeEnrichedNodalShapefunction(
        const XFEM::PHYSICS::Field field,
        const VI& funct,
        const MI1& derxy,
        const MI2& derxy2,
        VO& enr_funct,
        MO1& enr_derxy,
        MO2& enr_derxy2
    ) const
    {
      const int* nodeids = ele_.NodeIds();
      const std::size_t numnode = ele_.NumNode();

      std::size_t dofcounter = 0;
      for (std::size_t inode=0; inode < numnode; ++inode)
      {
        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman_.FieldEnrSetPerNode(gid);
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
          enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
          if (enrfield->getField() == field)
          {
            const double enrval = enrvals_.find(enrfield->getEnrichment())->second;
            enr_funct(dofcounter) = funct(inode) * enrval;
            for (std::size_t isd = 0; isd < 3; ++isd)
            {
              enr_derxy(isd,dofcounter) = derxy(isd,inode) * enrval;
            }
            for (std::size_t isd = 0; isd < 6; ++isd)
            {
              enr_derxy2(isd,dofcounter) = derxy2(isd,inode) * enrval;
            }
            dofcounter += 1;
          }
        }
      }
      dsassert(dofcounter == dofman_.NumDofPerField(field), "mismatch in information from eledofmanager!");
    }

   /*!
    * \brief compute the modified enriched shape functions and derivatives
    *
    * remarks:
    * - these are the modified, or shifted, enriched shape functions
    * - the chain rule is applied here for the derivatives
    * - for jump or void enrichments, zero terms are added, since the derivative of the enrichment
    *   function is zero
    */
    template<class VI, class MI1, class MI2, class VO, class MO1, class MO2 >
    inline void ComputeModifiedEnrichedNodalShapefunction(
        const XFEM::PHYSICS::Field field,
        const VI& funct,
        const MI1& derxy,
        const MI2& derxy2,
        VO& enr_funct,
        MO1& enr_derxy,
        MO2& enr_derxy2
    ) const
    {
      const int* nodeids = ele_.NodeIds();
      const std::size_t numnode = ele_.NumNode();

      std::size_t dofcounter = 0;
      for (std::size_t inode=0; inode < numnode; ++inode)
      {
        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman_.FieldEnrSetPerNode(gid);
        const std::map<XFEM::Enrichment, double> enrvalues = enrvalnodes_.find(inode)->second;

        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
          enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
          if (enrfield->getField() == field)
          {
            const double enrval = enrvalues.find(enrfield->getEnrichment())->second;
//std::cout << "enrval " << enrval << std::endl;
            const std::vector<double> enrvalderxy = enrvalderxy_.find(enrfield->getEnrichment())->second;
//std::cout << "enrvalderxy groesse: " << enrvalderxy.size() << std::endl;
//for (size_t count=0; count < enrvalderxy.size(); ++count)
//  std::cout << "enrvalderxy " << count << "wert " << enrvalderxy[count] << std::endl;

            const std::vector<double> enrvalderxy2 = enrvalderxy2_.find(enrfield->getEnrichment())->second;

            enr_funct(dofcounter) = funct(inode) * enrval;

            for (std::size_t isd = 0; isd < 3; ++isd)
            {
              //product rule
              enr_derxy(isd,dofcounter) = derxy(isd,inode) * enrval + funct(inode) * enrvalderxy[isd];
            }
            enr_derxy2(0,dofcounter) = derxy2(0,inode) * enrval + funct(inode) * enrvalderxy2[0] + derxy(0,inode) * enrvalderxy[0] + derxy(0,inode) * enrvalderxy[0];
            enr_derxy2(1,dofcounter) = derxy2(1,inode) * enrval + funct(inode) * enrvalderxy2[1] + derxy(1,inode) * enrvalderxy[1] + derxy(1,inode) * enrvalderxy[1];
            enr_derxy2(2,dofcounter) = derxy2(2,inode) * enrval + funct(inode) * enrvalderxy2[2] + derxy(2,inode) * enrvalderxy[2] + derxy(2,inode) * enrvalderxy[2];
            enr_derxy2(3,dofcounter) = derxy2(3,inode) * enrval + funct(inode) * enrvalderxy2[3] + derxy(0,inode) * enrvalderxy[1] + derxy(1,inode) * enrvalderxy[0];
            enr_derxy2(4,dofcounter) = derxy2(4,inode) * enrval + funct(inode) * enrvalderxy2[4] + derxy(0,inode) * enrvalderxy[2] + derxy(2,inode) * enrvalderxy[0];
            enr_derxy2(5,dofcounter) = derxy2(5,inode) * enrval + funct(inode) * enrvalderxy2[5] + derxy(1,inode) * enrvalderxy[2] + derxy(2,inode) * enrvalderxy[1];
            dofcounter += 1;
          }
        }
      }
      dsassert(dofcounter == dofman_.NumDofPerField(field), "mismatch in information from eledofmanager!");
    }


    /*!
     * \brief get nodal enrichments values
     *
     * remarks:
     * function for face-based stabilization to extend solution field to opposite nodes
     */
     template<class V>
     inline void GetNodalEnrichmentsFunctionValues(
         const XFEM::PHYSICS::Field field,
         std::vector<int>& num_field_enr_per_node,
         V& nodal_enr
     ) const
     {
       const int* nodeids = ele_.NodeIds();
       const std::size_t numnode = ele_.NumNode();
       num_field_enr_per_node.resize(numnode);

       std::size_t dofcounter = 0;
       for (std::size_t inode=0; inode < numnode; ++inode)
       {
         // set value to zero
         num_field_enr_per_node[inode] = 0;

         const int gid = nodeids[inode];
         const std::set<XFEM::FieldEnr>& enrfieldset = dofman_.FieldEnrSetPerNode(gid);
         const std::map<XFEM::Enrichment, double> enrvalues = enrvalnodes_.find(inode)->second;

         for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
           enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
         {
           if (enrfield->getField() == field)
           {
             const double enrval = enrvalues.find(enrfield->getEnrichment())->second;

             nodal_enr(dofcounter) = enrval;
             num_field_enr_per_node[inode] += 1;

             dofcounter += 1;
           }
         }
       }
       dsassert(dofcounter == dofman_.NumDofPerField(field), "mismatch in information from eledofmanager!");
     }


    template<class V, class M1, size_t numdof>
    void ComputeNormalShapeFunction(
        const V& funct,
        const M1& egradphi,
//#ifdef COLLAPSE_FLAME
//        const LINALG::Matrix<3,1>& normal,
//#endif
        XFEM::ApproxFuncNormalVector<0,numdof>& shp
    ) const
    {
      const int* nodeids = ele_.NodeIds();
      const std::size_t numnode = ele_.NumNode();

      //------------------------------------------------------------------------
      // interpolate gradient field of phi and its derivatives at Gaussian point
      //------------------------------------------------------------------------
      // gradient of phi
      LINALG::Matrix<3,1> gp(true);

      // loop over nodes of this element
      for (std::size_t inode=0;inode<numnode;++inode)
      {
        // compute the normal field n (not normed)
        // grad\phi = N_i*grad\phi_i
        gp(0) += funct(inode)*egradphi(0,inode); // N_i*grad\phi^x_i
        gp(1) += funct(inode)*egradphi(1,inode); // N_i*grad\phi^y_i
        gp(2) += funct(inode)*egradphi(2,inode); // N_i*grad\phi^z_i
      }

      // compute absolute value of gradient of phi
      // |grad\phi| = sqrt( grad\phi_x^2 + grad\phi_y^2 + grad\phi_z^2 )
#ifdef COMBUST_2D
      gp(2) = 0.0;
#endif
      const double normgp = gp.Norm2(); // sqrt(gp(0)*gp(0) + gp(1)*gp(1) + gp(2)*gp(2))
      if (normgp == 0.0) dserror("norm of gradient of phi is zero!");

      //---------------------------------------------------------------
      // compute the normal field and its derivatives at Gaussian point
      //---------------------------------------------------------------
      // remark: definition of normal field based on the level set function \phi:
      //         n = - grad\phi/|grad\phi|
      // compute the normal vector
      LINALG::Matrix<3,1> n = gp;
#ifdef COMBUST_TESTCOUETTEFLOWDECOUPLED
      n(1)=0.0;
#endif
#ifdef COMBUST_TESTCOUETTEFLOW
      n(1)=0.0;
#endif
      n.Scale(-1.0/normgp);
//#ifdef COLLAPSE_FLAME
//      n=normal;
//#endif
      //----------------------------------------
      // compute enriched normal shape functions
      //----------------------------------------
      // loop over nodes of this element
      std::size_t dofcounter = 0;
      for (std::size_t inode=0; inode < numnode; ++inode)
      {
        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman_.FieldEnrSetPerNode(gid);
        const std::map<XFEM::Enrichment, double> enrvalues = enrvalnodes_.find(inode)->second;

        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
            enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
          if (enrfield->getField() == XFEM::PHYSICS::Veln)
          {
            const double enrval = enrvalues.find(enrfield->getEnrichment())->second;

            // normal enriched shape functions
            // N_i*n*\phi_i
            shp.velx.d0.n(dofcounter) = funct(inode)*n(0)*enrval;
            shp.vely.d0.n(dofcounter) = funct(inode)*n(1)*enrval;
            shp.velz.d0.n(dofcounter) = funct(inode)*n(2)*enrval;

            dofcounter += 1;
          }
        }
      }
      dsassert(dofcounter == dofman_.NumDofPerField(XFEM::PHYSICS::Veln), "mismatch in information from eledofmanager!");
    }

   /*!
    * \brief compute the modified enriched shape functions and derivatives for normal enrichment
    *
    * remarks:
    * - these are the modified, or shifted, enriched shape functions
    * - the definition of the normal field is based on the level set field
    *   n = - grad\phi/|grad\phi|
    * - the product rule has to be applied applied here for the derivatives
    *
    * 1. the gradient of \phi 'gp' and its derivatives 'dergp' and 'der2gp' are interpolated
    *    to the Gaussian point from the nodal values 'egradphi' given to the function
    * 2. the normal field 'n' and its derivatives 'dn/dx' and 'd^2n/dxdx' are computed
    * 3. the enriched normal shape functions are computed
    */
    template<class V, class M1, class M2, size_t numdof>
    void ComputeNormalShapeFunction(
        const V& funct,
        const M1& derxy,
        const M2& derxy2,
        const M1& egradphi,
//#ifdef COLLAPSE_FLAME
//        const LINALG::Matrix<3,1>& normal,
//        const LINALG::Matrix<3,1>& point,
//#endif
        XFEM::ApproxFuncNormalVector<2,numdof>& shp
    ) const
    {
      const int* nodeids = ele_.NodeIds();
      const std::size_t numnode = ele_.NumNode();

      //------------------------------------------------------------------------
      // interpolate gradient field of phi and its derivatives at Gaussian point
      //------------------------------------------------------------------------
      // gradient of phi
      LINALG::Matrix<3,1> gp(true);
      // first derivatives of gradient of phi
      LINALG::Matrix<3,3> dergp(true);
      // second derivatives of gradient of phi
      LINALG::Matrix<3,6> der2gp(true);

      // loop over nodes of this element
      for (std::size_t inode=0;inode<numnode;++inode)
      {
        // compute the normal field n (not normed)
        // grad\phi = N_i*grad\phi_i
        gp(0) += funct(inode)*egradphi(0,inode); // N_i*grad\phi^x_i
        gp(1) += funct(inode)*egradphi(1,inode); // N_i*grad\phi^y_i
        gp(2) += funct(inode)*egradphi(2,inode); // N_i*grad\phi^z_i

        // compute first derivatives
        // d grad\phi/dx = dN_i/dx * grad\phi_i
        for (size_t isd=0;isd<3;++isd)
        {
          dergp(isd,0) += derxy(0,inode)*egradphi(isd,inode); // dN_i/dx * grad\phi^isd_i
          dergp(isd,1) += derxy(1,inode)*egradphi(isd,inode); // dN_i/dy * grad\phi^isd_i
          dergp(isd,2) += derxy(2,inode)*egradphi(isd,inode); // dN_i/dz * grad\phi^isd_i
        }

        // compute second derivatives
        // d^2grad\phi / dxdx = d^2N_i/dxdx * grad\phi_i
        for (size_t isd=0;isd<3;++isd)
        {
          der2gp(isd,0) += derxy2(0,inode)*egradphi(isd,inode); // d^2N_i/dxdx * grad\phi^isd_i
          der2gp(isd,1) += derxy2(1,inode)*egradphi(isd,inode); // d^2N_i/dydy * grad\phi^isd_i
          der2gp(isd,2) += derxy2(2,inode)*egradphi(isd,inode); // d^2N_i/dzdz * grad\phi^isd_i
          der2gp(isd,3) += derxy2(3,inode)*egradphi(isd,inode); // d^2N_i/dxdy * grad\phi^isd_i
          der2gp(isd,4) += derxy2(4,inode)*egradphi(isd,inode); // d^2N_i/dxdz * grad\phi^isd_i
          der2gp(isd,5) += derxy2(5,inode)*egradphi(isd,inode); // d^2N_i/dydz * grad\phi^isd_i
        }
      }
      //std::cout << "gp " << gp << std::endl;
      //std::cout << "dergp " << dergp << std::endl;
      //std::cout << "der2gp " << dergp << std::endl;

      // compute absolute value of gradient of phi
      // |grad\phi| = sqrt( grad\phi_x^2 + grad\phi_y^2 + grad\phi_z^2 )
#ifdef COMBUST_2D
      gp(2) = 0.0;
#endif
      const double normgp = gp.Norm2(); // sqrt(gp(0)*gp(0) + gp(1)*gp(1) + gp(2)*gp(2))
      if (normgp == 0.0) dserror("norm of gradient of phi is zero!");

      //---------------------------------------------------------------
      // compute the normal field and its derivatives at Gaussian point
      //---------------------------------------------------------------
      // remark: definition of normal field based on the level set function \phi:
      //         n = - grad\phi/|grad\phi|
      // compute the normal vector
      LINALG::Matrix<3,1> n = gp;
#ifdef COMBUST_TESTCOUETTEFLOWDECOUPLED
      n(1)=0.0;
#endif
#ifdef COMBUST_TESTCOUETTEFLOW
      n(1)=0.0;
#endif
      n.Scale(-1.0/normgp);
//#ifdef COLLAPSE_FLAME
//      n=normal;
//#endif
      // first derivative of normal vector
      LINALG::Matrix<3,3> dern(true);
      // dn_x/dx  dn_x/dy  dn_x/dz
      // dn_y/dx  dn_y/dy  dn_y/dz
      // dn_z/dx  dn_z/dy  dn_z/dz
      for (size_t idim=0;idim<3;++idim) // dimension n_x,n_y,n_z
      {
        for (size_t ider=0;ider<3;++ider) // derivative dx,dy,dz
        {
          dern(idim,ider) = 1/normgp*( -dergp(idim,ider) +
              (gp(0)*dergp(0,ider)+gp(1)*dergp(1,ider)+gp(2)*dergp(2,ider))
              / (normgp*normgp) * gp(idim));
        }
      }
//#ifdef COLLAPSE_FLAME
////std::cout << "normal " << normal << std::endl;
//      const double x = point(0);
////std::cout << "x " << x << std::endl;
//      const double y = point(1);
////std::cout << "y " << y << std::endl;
//      const double r2 = x*x + y*y;
////std::cout << "r2 " << r2 << std::endl;
//      const double r = std::sqrt(r2);
////std::cout << "r " << r << std::endl;
//      const double r2_32 = std::pow(r2,-1.5);
////std::cout << "r2_32 " << r2_32 << std::endl;
//      const double r2_52 = std::pow(r2,-2.5);
////std::cout << "r2_52 " << r2_52 << std::endl;
//      dern(0,0) = -1.0/r + x*x*r2_32;
//      dern(0,1) = x*y*r2_32;
//      dern(0,2) = 0.0;
//      dern(1,0) = y*x*r2_32;
//      dern(1,1) = -1.0/r + y*y*r2_32;
//      dern(1,2) = 0.0;
//      dern(2,0) = 0.0;
//      dern(2,1) = 0.0;
//      dern(2,2) = 0.0;
//#endif
      //std::cout << "dern " << dern << std::endl;

      // second derivative of normal vector
      LINALG::Matrix<3,6> der2n(true);
      // d^2n_x/dxdx  d^2n_x/dydy  d^2n_x/dzdz  d^2n_x/dxdy  d^2n_x/dxdz  d^2n_x/dydz
      // d^2n_y/dxdx  d^2n_y/dydy  d^2n_y/dzdz  d^2n_y/dxdy  d^2n_y/dxdz  d^2n_y/dydz
      // d^2n_z/dxdx  d^2n_z/dydy  d^2n_z/dzdz  d^2n_z/dxdy  d^2n_z/dxdz  d^2n_z/dydz
      // pre-compute some terms
      LINALG::Matrix<3,1> term(true);
      term(0) = gp(0)*dergp(0,0)+gp(1)*dergp(1,0)+gp(2)*dergp(2,0); // grad\phi * d grad\phi/dx
      term(1) = gp(0)*dergp(0,1)+gp(1)*dergp(1,1)+gp(2)*dergp(2,1); // grad\phi * d grad\phi/dy
      term(2) = gp(0)*dergp(0,2)+gp(1)*dergp(1,2)+gp(2)*dergp(2,2); // grad\phi * d grad\phi/dz
      // use these variables to accelerate loop
      //               const double normgp2 = normgp*normgp;
      //               const double normgp3 = normgp*normgp*normgp;
      //               replace terms which do not have 'idim' in them by pre-computed variables
      for (size_t idim=0;idim<3;++idim) // dimension n_x,n_y,n_z
      {
        // d^2n/dxdx
        der2n(idim,0) = (-normgp*normgp*der2gp(idim,0) + dergp(idim,0)*term(0) - 3/(normgp*normgp)*term(0)*term(0)*gp(idim)
            + (dergp(0,0)*dergp(0,0) + gp(0)*der2gp(0,0)
            +  dergp(1,0)*dergp(1,0) + gp(1)*der2gp(1,0)
            +  dergp(2,0)*dergp(2,0) + gp(2)*der2gp(2,0))*gp(idim)
            + term(0)*dergp(idim,0) )/(normgp*normgp*normgp);
        // d^2n/dydy
        der2n(idim,1) = (-normgp*normgp*der2gp(idim,1) + dergp(idim,1)*term(1) - 3/(normgp*normgp)*term(1)*term(1)*gp(idim)
            + (dergp(0,1)*dergp(0,1) + gp(0)*der2gp(0,1)
            +  dergp(1,1)*dergp(1,1) + gp(1)*der2gp(1,1)
            +  dergp(2,1)*dergp(2,1) + gp(2)*der2gp(2,1))*gp(idim)
            + term(1)*dergp(idim,1) )/(normgp*normgp*normgp);
        // d^2n/dzdz
        der2n(idim,2) = (-normgp*normgp*der2gp(idim,2) + dergp(idim,2)*term(2) - 3/(normgp*normgp)*term(2)*term(2)*gp(idim)
            + (dergp(0,2)*dergp(0,2) + gp(0)*der2gp(0,2)
            +  dergp(1,2)*dergp(1,2) + gp(1)*der2gp(1,2)
            +  dergp(2,2)*dergp(2,2) + gp(2)*der2gp(2,2))*gp(idim)
            + term(2)*dergp(idim,2) )/(normgp*normgp*normgp);
        // d^2n/dxdy
        der2n(idim,3) = (-normgp*normgp*der2gp(idim,3) + dergp(idim,0)*term(1) - 3/(normgp*normgp)*term(0)*term(1)*gp(idim)
            + (dergp(0,0)*dergp(0,1) + gp(0)*der2gp(0,3)
            +  dergp(1,0)*dergp(1,1) + gp(1)*der2gp(1,3)
            +  dergp(2,0)*dergp(2,1) + gp(2)*der2gp(2,3))*gp(idim)
            + term(0)*dergp(idim,1) )/(normgp*normgp*normgp);
        // d^2n/dxdz
        der2n(idim,4) = (-normgp*normgp*der2gp(idim,4) + dergp(idim,0)*term(2) - 3/(normgp*normgp)*term(0)*term(2)*gp(idim)
            + (dergp(0,0)*dergp(0,2) + gp(0)*der2gp(0,4)
            +  dergp(1,0)*dergp(1,2) + gp(1)*der2gp(1,4)
            +  dergp(2,0)*dergp(2,2) + gp(2)*der2gp(2,4))*gp(idim)
            + term(0)*dergp(idim,2) )/(normgp*normgp*normgp);
        // d^2n/dydz
        der2n(idim,5) = (-normgp*normgp*der2gp(idim,5) + dergp(idim,1)*term(2) - 3/(normgp*normgp)*term(1)*term(2)*gp(idim)
            + (dergp(0,1)*dergp(0,2) + gp(0)*der2gp(0,5)
            +  dergp(1,1)*dergp(1,2) + gp(1)*der2gp(1,5)
            +  dergp(2,1)*dergp(2,2) + gp(2)*der2gp(2,5))*gp(idim)
            + term(1)*dergp(idim,2) )/(normgp*normgp*normgp);
      }
//#ifdef COLLAPSE_FLAME
//      // second derivative of normal vector
//      der2n.Clear();
//      // d2n_1/dxdx
//      der2n(0,0) = 3.0*x*r2_32 - 3.0*x*x*x*r2_52;
//      // d2n_1/dydy
//      der2n(0,1) = x*r2_32 - 3.0*x*y*y*r2_52;
//      // d2n_1/dxdy = d2n_1/dydx
//      der2n(0,3) = y*r2_32 - 3.0*x*x*y*r2_52;
//
//      // d2n_2/dxdx
//      der2n(1,0) = y*r2_32 - 3.0*x*x*y*r2_52;
//      // d2n_2/dydy
//      der2n(1,1) = 3.0*y*r2_32 - 3.0*y*y*y*r2_52;
//      // d2n_2/dxdy = d2n_2/dydx
//      der2n(1,3) = x*r2_32 - 3.0*x*y*y*r2_52;
//#endif
      //----------------------------------------
      // compute enriched normal shape functions
      //----------------------------------------
      // loop over nodes of this element
      std::size_t dofcounter = 0;
      for (std::size_t inode=0; inode < numnode; ++inode)
      {
        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman_.FieldEnrSetPerNode(gid);
        const std::map<XFEM::Enrichment, double> enrvalues = enrvalnodes_.find(inode)->second;

        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
            enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
          if (enrfield->getField() == XFEM::PHYSICS::Veln)
          {
            const double enrval = enrvalues.find(enrfield->getEnrichment())->second;

            // normal enriched shape functions
            // N_i*n*\phi_i
            shp.velx.d0.n(dofcounter) = funct(inode)*n(0)*enrval;
            shp.vely.d0.n(dofcounter) = funct(inode)*n(1)*enrval;
            shp.velz.d0.n(dofcounter) = funct(inode)*n(2)*enrval;

            // first derivatives of normal enriched shape functions
            // (dN_i/dx*n + N_i*dn/dx)*\phi_i
            shp.velx.dx.n(dofcounter) = (derxy(0,inode)*n(0) + funct(inode)*dern(0,0)) * enrval;
            shp.vely.dx.n(dofcounter) = (derxy(0,inode)*n(1) + funct(inode)*dern(1,0)) * enrval;
            shp.velz.dx.n(dofcounter) = (derxy(0,inode)*n(2) + funct(inode)*dern(2,0)) * enrval;

            shp.velx.dy.n(dofcounter) = (derxy(1,inode)*n(0) + funct(inode)*dern(0,1)) * enrval;
            shp.vely.dy.n(dofcounter) = (derxy(1,inode)*n(1) + funct(inode)*dern(1,1)) * enrval;
            shp.velz.dy.n(dofcounter) = (derxy(1,inode)*n(2) + funct(inode)*dern(2,1)) * enrval;

            shp.velx.dz.n(dofcounter) = (derxy(2,inode)*n(0) + funct(inode)*dern(0,2)) * enrval;
            shp.vely.dz.n(dofcounter) = (derxy(2,inode)*n(1) + funct(inode)*dern(1,2)) * enrval;
            shp.velz.dz.n(dofcounter) = (derxy(2,inode)*n(2) + funct(inode)*dern(2,2)) * enrval;

            // second derivatives of normal enriched shape functions
            // (d^2N_i/dxdx*n + dN_i/dx*dn/dx + dN_i/dx*dn/dx + N_i*d^2n/dxdx)*\phi_i
            shp.velx.dxdx.n(dofcounter) = (derxy2(0,inode)*n(0) + funct(inode)*der2n(0,0) + 2.0*derxy(0,inode)*dern(0,0)) * enrval;
            shp.vely.dxdx.n(dofcounter) = (derxy2(0,inode)*n(1) + funct(inode)*der2n(1,0) + 2.0*derxy(0,inode)*dern(1,0)) * enrval;
            shp.velz.dxdx.n(dofcounter) = (derxy2(0,inode)*n(2) + funct(inode)*der2n(2,0) + 2.0*derxy(0,inode)*dern(2,0)) * enrval;

            shp.velx.dydy.n(dofcounter) = (derxy2(1,inode)*n(0) + funct(inode)*der2n(0,1) + 2.0*derxy(1,inode)*dern(0,1)) * enrval;
            shp.vely.dydy.n(dofcounter) = (derxy2(1,inode)*n(1) + funct(inode)*der2n(1,1) + 2.0*derxy(1,inode)*dern(1,1)) * enrval;
            shp.velz.dydy.n(dofcounter) = (derxy2(1,inode)*n(2) + funct(inode)*der2n(2,1) + 2.0*derxy(1,inode)*dern(2,1)) * enrval;

            shp.velx.dzdz.n(dofcounter) = (derxy2(2,inode)*n(0) + funct(inode)*der2n(0,2) + 2.0*derxy(2,inode)*dern(0,2)) * enrval;
            shp.vely.dzdz.n(dofcounter) = (derxy2(2,inode)*n(1) + funct(inode)*der2n(1,2) + 2.0*derxy(2,inode)*dern(1,2)) * enrval;
            shp.velz.dzdz.n(dofcounter) = (derxy2(2,inode)*n(2) + funct(inode)*der2n(2,2) + 2.0*derxy(2,inode)*dern(2,2)) * enrval;

            shp.velx.dxdy.n(dofcounter) = (derxy2(3,inode)*n(0) + funct(inode)*der2n(0,3) + derxy(0,inode)*dern(0,1) + derxy(1,inode)*dern(0,0)) * enrval;
            shp.vely.dxdy.n(dofcounter) = (derxy2(3,inode)*n(1) + funct(inode)*der2n(1,3) + derxy(0,inode)*dern(1,1) + derxy(1,inode)*dern(1,0)) * enrval;
            shp.velz.dxdy.n(dofcounter) = (derxy2(3,inode)*n(2) + funct(inode)*der2n(2,3) + derxy(0,inode)*dern(2,1) + derxy(1,inode)*dern(2,0)) * enrval;

            shp.velx.dxdz.n(dofcounter) = (derxy2(4,inode)*n(0) + funct(inode)*der2n(0,4) + derxy(0,inode)*dern(0,2) + derxy(2,inode)*dern(0,0)) * enrval;
            shp.vely.dxdz.n(dofcounter) = (derxy2(4,inode)*n(1) + funct(inode)*der2n(1,4) + derxy(0,inode)*dern(1,2) + derxy(2,inode)*dern(1,0)) * enrval;
            shp.velz.dxdz.n(dofcounter) = (derxy2(4,inode)*n(2) + funct(inode)*der2n(2,4) + derxy(0,inode)*dern(2,2) + derxy(2,inode)*dern(2,0)) * enrval;

            shp.velx.dydz.n(dofcounter) = (derxy2(5,inode)*n(0) + funct(inode)*der2n(0,5) + derxy(1,inode)*dern(0,2) + derxy(2,inode)*dern(0,1)) * enrval;
            shp.vely.dydz.n(dofcounter) = (derxy2(5,inode)*n(1) + funct(inode)*der2n(1,5) + derxy(1,inode)*dern(1,2) + derxy(2,inode)*dern(1,1)) * enrval;
            shp.velz.dydz.n(dofcounter) = (derxy2(5,inode)*n(2) + funct(inode)*der2n(2,5) + derxy(1,inode)*dern(2,2) + derxy(2,inode)*dern(2,1)) * enrval;

            dofcounter += 1;
          }
        }
      }
      shp.velx.dydx.n = shp.velx.dxdy.n;
      shp.vely.dydx.n = shp.vely.dxdy.n;
      shp.velz.dydx.n = shp.velz.dxdy.n;

      shp.velx.dzdx.n = shp.velx.dxdz.n;
      shp.vely.dzdx.n = shp.vely.dxdz.n;
      shp.velz.dzdx.n = shp.velz.dxdz.n;

      shp.velx.dzdy.n = shp.velx.dydz.n;
      shp.vely.dzdy.n = shp.vely.dydz.n;
      shp.velz.dzdy.n = shp.velz.dydz.n;

      dsassert(dofcounter == dofman_.NumDofPerField(XFEM::PHYSICS::Veln), "mismatch in information from eledofmanager!");
    }


    /*!
    * \brief For a given situation compute the enriched shape functions
    *
    * simplest case: jump or void enrichment
    * -> no chain rule, since enrichment function derivative is zero
    * remark: if enrichment function has derivatives not equal to zero, we need the chain rule here
    */
    template<class VI, class VO>
    void ComputeEnrichedElementShapefunction(
        const XFEM::PHYSICS::Field field,
        const VI&                  funct,
        VO&                        enr_funct
    ) const
    {
      std::size_t dofcounter = 0;

      const std::set<XFEM::FieldEnr>& enrfieldset(dofman_.getEnrichedFieldsPerEleField(field));
      const DRT::Element::DiscretizationType distype = dofman_.getDisTypePerField(field);
      const std::size_t numvirtualnode = DRT::UTILS::getNumberOfElementNodes(distype);
      dsassert(enrfieldset.size() > 0, "empty enrfieldset not allowed at this point!");
      for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
        enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
      {
        for (std::size_t inode = 0; inode < numvirtualnode; ++inode)
        {
          const double enrval = enrvals_.find(enrfield->getEnrichment())->second;
          enr_funct(dofcounter) = funct(inode) * enrval;
          dofcounter += 1;
        }
      }
      dsassert(dofcounter == dofman_.NumDofPerField(field), "mismatch in information from eledofmanager!");
    }

    /*!
     * \brief compute the enriched shape functions for element dofs
     *
     * remarks:
     * - simplest case: jump or void enrichment
     *   -> no chain rule, since derivative of enrichment function is zero
     * - be carefull: if enrichment function has derivatives not equal to zero, we need the chain rule here!
     *   -> for kink enrichment use ComputeKinkEnrichedNodalShapefunction()
     */
    template<class VI, class VO>
    inline void ComputeShiftedEnrichedElementShapefunction(
        const XFEM::PHYSICS::Field                 field,
        const VI&                                  funct,
        VO&                                        enr_funct
    ) const
    {
      //const int* nodeids = ele_.NodeIds();
      //const std::size_t numnode = ele_.NumNode();

      std::size_t dofcounter = 0;

      const std::set<XFEM::FieldEnr>& enrfieldset(dofman_.getEnrichedFieldsPerEleField(field));
      dsassert(enrfieldset.size() > 0, "empty enrfieldset not allowed at this point!");

      const DRT::Element::DiscretizationType distype = dofman_.getDisTypePerField(field);
      const int numvirtualnode = DRT::UTILS::getNumberOfElementNodes(distype);
      if (numvirtualnode != ele_.NumNode())
        dserror("you can not use the shifted enrichment without correcting the implementation!");

      for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
        enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
      {
        for (int inode = 0; inode < numvirtualnode; ++inode)
        {
          const std::map<XFEM::Enrichment, double> enrvalues = enrvalnodes_.find(inode)->second;
          const double enrval = enrvalues.find(enrfield->getEnrichment())->second;
          //const double enrval = enrvals_.find(enrfield->getEnrichment())->second;
          enr_funct(dofcounter) = funct(inode) * enrval;
          dofcounter += 1;
        }
      }
      dsassert(dofcounter == dofman_.NumDofPerField(field), "mismatch in information from eledofmanager!");
    }


    /*!
    * \brief For a given situation compute the enriched shape functions
    *
    * simplest case: jump or void enrichment
    * -> no chain rule, since enrichment function derivative is zero
    * remark: if enrichment function has derivatives not equal to zero, we need the chain rule here
    */
    template<class VI, class MI1, class VO, class MO1>
    void ComputeEnrichedElementShapefunction(
        const XFEM::PHYSICS::Field field,
        const VI&                  funct,
        const MI1&                 derxy,
        VO&                        enr_funct,
        MO1&                       enr_derxy
    ) const
    {
      std::size_t dofcounter = 0;

      const std::set<XFEM::FieldEnr>& enrfieldset(dofman_.getEnrichedFieldsPerEleField(field));
      const DRT::Element::DiscretizationType distype = dofman_.getDisTypePerField(field);
      const std::size_t numvirtualnode = DRT::UTILS::getNumberOfElementNodes(distype);
      dsassert(enrfieldset.size() > 0, "empty enrfieldset not allowed at this point!");
      for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
        enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
      {
        for (std::size_t inode = 0; inode < numvirtualnode; ++inode)
        {
          const double enrval = enrvals_.find(enrfield->getEnrichment())->second;
          enr_funct(dofcounter) = funct(inode) * enrval;
          for (std::size_t isd = 0; isd < 3; ++isd)
          {
            enr_derxy(isd,dofcounter) = derxy(isd,inode) * enrval;
          }
          dofcounter += 1;
        }
      }
      dsassert(dofcounter == dofman_.NumDofPerField(field), "mismatch in information from eledofmanager!");
    }

  private:

    //! don't want copy constructor
    explicit ElementEnrichmentValues(const ElementEnrichmentValues&);

    //! don't want assignement operator
    ElementEnrichmentValues& operator = (const ElementEnrichmentValues& old);

    //! don't want default constructor
    explicit ElementEnrichmentValues();

    const DRT::Element&                        ele_;
    const XFEM::ElementDofManager&             dofman_;

    //! map holding values for enrichments
    std::map<XFEM::Enrichment, double> enrvals_;
    //! map holding values for enrichments for every node
    std::map<int,std::map<XFEM::Enrichment, double> > enrvalnodes_;

    // maps required for enrichment functions with non-zero derivatives (e.g. kink)
    std::map<XFEM::Enrichment, std::vector<double> > enrvalderxy_;
    std::map<XFEM::Enrichment, std::vector<double> > enrvalderxy2_;

  };


  //! interpolate from nodal vector array to integration point vector using the shape function
  template <class M, class VS>
  LINALG::Matrix<3,1> interpolateVectorFieldToIntPointNormal(
      const M&                       eleVectorField, ///< array with nodal vector values
      const VS&                      shp,            ///< array with nodal shape function
      const size_t                   numnode,        ///< number of parameters
      const DRT::Element*            ele,            ///< pointer to element
      const XFEM::ElementDofManager& dofman          ///< dofmanager of this element
      )
  {
    LINALG::Matrix<3,1> v(true);

    const int* nodeids = ele->NodeIds();

    std::size_t velncounter = 0;
    for (std::size_t inode=0; inode<numnode; ++inode)
    {
      // standard shape functions are identical for all vector components
      // shp.velx.d0.s == shp.vely.d0.s == shp.velz.d0.s
      v(0) += eleVectorField(0,inode)*shp.velx.d0.s(inode);
      v(1) += eleVectorField(1,inode)*shp.vely.d0.s(inode);
      v(2) += eleVectorField(2,inode)*shp.velz.d0.s(inode);

      const int gid = nodeids[inode];
      const std::set<XFEM::FieldEnr>& enrfieldset = dofman.FieldEnrSetPerNode(gid);

      for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
          enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
      {
        if (enrfield->getField() == XFEM::PHYSICS::Veln)
        {
          v(0) += eleVectorField(3,velncounter)*shp.velx.d0.n(velncounter);
          v(1) += eleVectorField(3,velncounter)*shp.vely.d0.n(velncounter);
          v(2) += eleVectorField(3,velncounter)*shp.velz.d0.n(velncounter);

          velncounter += 1;
        }
      }
    }
    if (velncounter != dofman.NumDofPerField(XFEM::PHYSICS::Veln)) dserror("Alles falsch, du Depp!");
    dsassert(velncounter == dofman.NumDofPerField(XFEM::PHYSICS::Veln), "mismatch in information from eledofmanager!");

    return v;
  }


/*! \brief for output, interpolate field to integration cells
 *
 * interpolate field from element node values to cell node values based
 * on a level-set field
 * remark: function used for modified jump enrichment strategy
 */
void InterpolateCellValuesFromElementValuesLevelSet(
  const DRT::Element&                   ele,
  const XFEM::ElementDofManager&        dofman,
  const GEO::DomainIntCell&             cell,
  const std::vector<double>&            phivalues,
  const XFEM::PHYSICS::Field            field,
  const LINALG::SerialDenseMatrix&      elementvalues,
  LINALG::SerialDenseMatrix&            cellvalues
  );

void InterpolateCellValuesFromElementValuesLevelSetNormal(
  const DRT::Element&                   ele,
  const XFEM::ElementDofManager&        dofman,
  const GEO::DomainIntCell&             cell,
  const std::vector<double>&            phivalues,
  const LINALG::Matrix<3,8>&            gradphi,
  const XFEM::PHYSICS::Field            field,
  const LINALG::SerialDenseMatrix&      elementvalues,
  LINALG::SerialDenseMatrix&            cellvalues
  );

/*! \brief for output, interpolate field to integration cells
 *
 * interpolate field from element node values to cell node values based
 * on a level-set field
 * remark: function used for modified kink enrichment strategy
 */
void InterpolateCellValuesFromElementValuesLevelSetKink(
  const DRT::Element&                   ele,
  const XFEM::ElementDofManager&        dofman,
  const GEO::DomainIntCell&             cell,
  const std::vector<double>&            phivalues,
  const XFEM::PHYSICS::Field            field,
  const LINALG::SerialDenseMatrix&      elementvalues,
  LINALG::SerialDenseMatrix&            cellvalues
  );

/*! \brief for output, interpolate field to integration cells
 * schott
 *
 * interpolate field from element node values to cell node values based
 * on a level-set field
 * remark: function used for modified kink enrichment strategy with jumps in pressure and kinks in velocity
 */
void InterpolateCellValuesFromElementValuesLevelSetKinkJump(
  const DRT::Element&                   ele,
  const XFEM::ElementDofManager&        dofman,
  const GEO::DomainIntCell&             cell,
  const std::vector<double>&            phivalues,
  const XFEM::PHYSICS::Field            field,
  const LINALG::SerialDenseMatrix&      elementvalues,
  LINALG::SerialDenseMatrix&            cellvalues
  );


//! for output, compute values on the integration cells
void computeScalarCellNodeValuesFromElementUnknowns(
  const DRT::Element&                   ele,
  COMBUST::InterfaceHandleCombust*      ih,
  const XFEM::ElementDofManager&        dofman,
  const GEO::DomainIntCell&             cell,
  const XFEM::PHYSICS::Field            field,
  const LINALG::SerialDenseVector&      elementvalues,
  LINALG::SerialDenseVector&            cellvalues
  );


//! for output, compute values on the integration cells
void computeTensorCellNodeValuesFromElementUnknowns(
  const DRT::Element&               ele,
  COMBUST::InterfaceHandleCombust*  ih,
  const XFEM::ElementDofManager&    dofman,
  const GEO::DomainIntCell&         cell,
  const XFEM::PHYSICS::Field        field,
  const LINALG::SerialDenseMatrix&  elementvalues,
  LINALG::SerialDenseMatrix&        cellvalues
  );


//! return assembly type depending on the DOF information
XFEM::AssemblyType ComputeAssemblyType(
    const ElementDofManager&  eleDofManager,
    const std::size_t         numnode,
    const int*                nodeids
);

/*! get number of unknowns for this field
 * use the dofmanager in the general case,
 * otherwise use a given number numnode
 */
template<XFEM::AssemblyType>
inline std::size_t getNumParam(
        const XFEM::ElementDofManager&  dofman,
        const XFEM::PHYSICS::Field      field,
        const std::size_t)
{
    return dofman.NumDofPerField(field);
}
template<>
inline std::size_t getNumParam<XFEM::standard_assembly>(
        const XFEM::ElementDofManager&,
        const XFEM::PHYSICS::Field,
        const std::size_t numnode)
{
    return numnode;
}

//! determine number of parameter
template<std::size_t NUMNODE, XFEM::AssemblyType ASSTYPE>
struct NumParam
{
  static std::size_t get(
     const XFEM::ElementDofManager& dofman,
     const XFEM::PHYSICS::Field     field
     )
  {
    return dofman.NumDofPerField(field);
  }
};

/// determine number of parameter for standard_assembly
template<std::size_t NUMNODE> struct NumParam<NUMNODE,XFEM::standard_assembly>
{
  static std::size_t get(
      const XFEM::ElementDofManager& ,
      const XFEM::PHYSICS::Field
      )
  {
    return NUMNODE;
  }
};

}

#endif  // #ifndef ENRICHMENT_UTILS_H
