/*!
\file enrichment_utils.H

\brief describes the enrichment types and classes

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef ENRICHMENT_UTILS_H
#define ENRICHMENT_UTILS_H

#include "interface.H"
#include "../drt_combust/combust_interface.H"
#include "../drt_inpar/inpar_combust.H"
#include "dof_management_element.H"
#include "xfem_enums.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_lib/drt_discret.H"

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
        XFEM::InterfaceHandle*                ih,               ///< interface information
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
            const GEO::DomainIntCell&       cell,   ///< pointer to domain integration cell
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
          //cout << "shifted enrichment value " << enrvalstmp[*enriter] << endl;
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

        if(enriter->Type() != XFEM::Enrichment::typeKink){
          enrvalside = enriter->EnrValueAtInterface(approachdirection);       //version Flo
          //enrvalside = enriter->EnrValueGaussPointLevelSet(ele_,phi,funct); //version Ursula
          // no derivatives needed for not kink enrichments
        }
        else{ // kink enrichment
          enrvalside = enriter->EnrValueGaussPointLevelSet(ele_,phi,funct);
          enrvalderxy = enriter->EnrValueGaussPointDerivative(ele_,phi,funct,derxy);
          enrvalderxy2 = enriter->EnrValueGaussPointDerivative2(ele_,phi,funct,derxy2);
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
            enrvalnode = enriter->EnrValueNodeLevelSet(inode,phi);
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
//cout << "enrichment value: " << enrval << endl;
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
//cout << "enrval " << enrval << endl;
            const std::vector<double> enrvalderxy = enrvalderxy_.find(enrfield->getEnrichment())->second;
//cout << "enrvalderxy groesse: " << enrvalderxy.size() << endl;
//for (size_t count=0; count < enrvalderxy.size(); ++count)
//  cout << "enrvalderxy " << count << "wert " << enrvalderxy[count] << endl;

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
    * \brief For a given situation compute the enriched shape functions
    *
    * simplest case: jump or void enrichment
    * -> no chain rule, since enrichment function derivative is zero
    * TODO: if enrichment function has derivatives not equal to zero, we need the chain rule here
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
    * TODO: if enrichment function has derivatives not equal to zero, we need the chain rule here
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
    //  const RCP<XFEM::InterfaceHandle>           ih_;
    const XFEM::ElementDofManager&             dofman_;

    //! map holding values for enrichments
    std::map<XFEM::Enrichment, double>        enrvals_;
    //! map holding values for enrichments for every node
    std::map<int,std::map<XFEM::Enrichment, double> > enrvalnodes_;

    // maps required for enrichment functions with non-zero derivatives (e.g. kink)
    std::map<XFEM::Enrichment, std::vector<double> > enrvalderxy_;
    std::map<XFEM::Enrichment, std::vector<double> > enrvalderxy2_;

  };


/*! \brief for output, compute values on integration cells
 *
 * used for boundary discretization based XFEM problems
 */
void computeScalarCellNodeValuesFromNodalUnknowns(
  const DRT::Element&                 ele,
  XFEM::InterfaceHandle*  ih,
  const XFEM::ElementDofManager&      dofman,
  const GEO::DomainIntCell&           cell,
  const XFEM::PHYSICS::Field          field,
  const LINALG::SerialDenseVector&    elementvalues,
  LINALG::SerialDenseVector&          cellvalues
  );


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
  XFEM::InterfaceHandle*     ih,
  const XFEM::ElementDofManager&        dofman,
  const GEO::DomainIntCell&             cell,
  const XFEM::PHYSICS::Field            field,
  const LINALG::SerialDenseVector&      elementvalues,
  LINALG::SerialDenseVector&            cellvalues
  );


//! for output, compute values on the integration cells
void computeTensorCellNodeValuesFromElementUnknowns(
  const DRT::Element&  ele,
  XFEM::InterfaceHandle*  ih,
  const XFEM::ElementDofManager& dofman,
  const GEO::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const LINALG::SerialDenseMatrix& elementvalues,
  LINALG::SerialDenseMatrix&     cellvalues
  );

//! for output, compute values on the integration cells
void computeVectorCellNodeValues(
  const DRT::Element&  ele,
  XFEM::InterfaceHandle*  ih,
  const XFEM::ElementDofManager& dofman,
  const GEO::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const LINALG::SerialDenseMatrix& elementvalues,
  LINALG::SerialDenseMatrix&     cellvalues
  );

//! for output, compute values on the integration cells
void computeVectorCellNodeValues(
  const DRT::Element&  ele,
  XFEM::InterfaceHandle*  ih,
  const XFEM::ElementDofManager& dofman,
  const GEO::BoundaryIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const int label,
  const LINALG::SerialDenseMatrix& elementvalues,
  LINALG::SerialDenseMatrix&     cellvalues
  );

/*!
 * \brief calculates the fraction of an element that is covered by fictitious domains.
 * \warning should be done separately for each enrichment, but currently it's not done
 *
 * \return ration between fictitious volume and element volume (can be between 0.0 and 1.0)
 */
double DomainCoverageRatio(
    const DRT::Element&           ele,           ///< the element whose area ratio we want to compute
    const XFEM::InterfaceHandle&  ih             ///< interface information
    );

/*!
 * \brief calculates the fraction of an element that is covered by fictitious domains.
 * \warning should be done separately for each enrichment, but currently it's not done
 *
 * \return ratio between fictitious volume and element volume (can be between 0.0 and 1.0)
 */
std::vector<double> DomainCoverageRatioPerNode(
    const DRT::Element&           ele,           ///< the element whose area ratio we want to compute
    const XFEM::InterfaceHandle&  ih             ///< interface information
    );

/*!
 * \brief calculates the area fraction of an element that is covered by fictitious domains.
 * \warning should be done separately for each enrichment, but currently it's not done
 */
double BoundaryCoverageRatio(
    const DRT::Element&               xele,          ///< the element whose area ratio we want to compute
    const GEO::BoundaryIntCells&      boundaryIntCells,
    const XFEM::InterfaceHandle&      ih             ///< interface information
    );

/*!
 * \brief calculates the volume fraction for each cell of an element .
 *
 * \return ratio between cell volume and element volume (can be between 0.0 and 1.0)
 */
std::vector<double> DomainIntCellCoverageRatio(
        const DRT::Element&           ele,
        const XFEM::InterfaceHandle&  ih
        );

//! return assembly type depending on the DOF information
XFEM::AssemblyType ComputeAssemblyType(
    const ElementDofManager&   eleDofManager,
    const std::size_t          numnode,
    const int*                 nodeids
);

/*! get number of unknowns for this field
 * use the dofmanager in the general case,
 * otherwise use a given number numnode
 */
template<XFEM::AssemblyType>
inline std::size_t getNumParam(
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field     field,
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
#endif  // #ifdef CCADISCRET
