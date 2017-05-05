/*!
\file combust3_utils.H

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/

#ifndef COMBUST3_UTILS_H
#define COMBUST3_UTILS_H


#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_node.H"


namespace COMBUST
{

  template <size_t numparam>
  struct EnrViscs2
  {
    LINALG::Matrix<numparam,1> xx;
    LINALG::Matrix<numparam,1> xy;
    LINALG::Matrix<numparam,1> xz;
    LINALG::Matrix<numparam,1> yx;
    LINALG::Matrix<numparam,1> yy;
    LINALG::Matrix<numparam,1> yz;
    LINALG::Matrix<numparam,1> zx;
    LINALG::Matrix<numparam,1> zy;
    LINALG::Matrix<numparam,1> zz;
  };

  /*!
   * \brief calculate the body force from the nodes condition
   *
   *  get the body force in the nodes of the element
   *  the Neumann condition associated with the nodes is stored in the
   *  array edeadng only if all nodes have a VolumeNeumann condition
   *
   *  \author rasthofer 08/09
   */
  template <DRT::Element::DiscretizationType DISTYPE,
            class M>
  void BodyForce(
      const DRT::Element* ele,
      const double time,
      const std::string neumtype,
      M& edead)
  {
    const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    std::vector<DRT::Condition*> myneumcond;
    int myindex = -1;

    // check whether all nodes have a unique VolumeNeumann condition
    DRT::UTILS::FindElementConditions(ele, "VolumeNeumann", myneumcond);

    for (size_t i = 0; i < myneumcond.size(); ++i)
    {
      const std::string* condtype = myneumcond[i]->Get<std::string>("type");
      if (*condtype == neumtype or (neumtype == "DeadOrLive" and (*condtype == "neum_dead" or *condtype == "neum_live")))
      {
        myindex = (int)i;
        break;
      }
    }

    if (myindex >= 0)
    {

      // get values and switches from the condition
      const std::vector<int>*    onoff = myneumcond[myindex]->Get<std::vector<int> >   ("onoff");
      const std::vector<double>* val   = myneumcond[myindex]->Get<std::vector<double> >("val"  );
      const std::vector<int>*    functions = myneumcond[myindex]->Get<std::vector<int> >("funct");

      // factor given by spatial function
      double functionfac = 1.0;
      int functnum = -1;

      // set this condition to the edeadng array
      for(int isd=0;isd<3;isd++)
      {
        // get factor given by spatial function
        if (functions) functnum = (*functions)[isd];
        else functnum = -1;

        double num = (*onoff)[isd]*(*val)[isd];

        for (std::size_t jnode=0; jnode<numnode; jnode++)
        {
          if (functnum>0)
          {
            // evaluate function at the position of the current node
            functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(isd,(ele->Nodes()[jnode])->X(),time);
          }
          else functionfac = 1.0;

          edead(isd,jnode) = num*functionfac;
        }
      }
    }
    else if(myneumcond.size()==2)
    {

    }
  }


  //! compute largest element diameter for Nitsche stabilization parameter
  template<DRT::Element::DiscretizationType DISTYPE, class M1>
  double getEleDiameter(const M1& xyze)
  {
    double elediam = 0.0;

    // number of nodes of this element
    const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // check all possible connections between nodes of an element
    // node 1 to 2
    //    :
    // node 1 to 8 = numnode
    // node 2 to 3
    //    :
    // node 2 to 8
    //    :
    //    :
    //    :
    // node 7 to 8
    for(size_t i_start=0; i_start< numnode-2; ++i_start)
    {
      for(size_t i_end= i_start+1; i_end < numnode-1; ++i_end)
      {
        LINALG::Matrix<3,1> direction;
        direction.Clear();
        direction(0) = xyze(0, i_start) - xyze(0, i_end);
        direction(1) = xyze(1, i_start) - xyze(1, i_end);
        direction(2) = xyze(2, i_start) - xyze(2, i_end);

        // update elediam
        if (direction.Norm2() > elediam) elediam=direction.Norm2();
      }
    }

    return elediam;
  }

  template <DRT::Element::DiscretizationType DISTYPE,
  class M1, class M2, class M3, class M4>
  void CalcCurvature(
      double&   curvature,          ///< curvature to be computed
      M1&       posXiDomain,        ///< position of gaussian point for evaluating of surface tension terms
      M2&       xyze,               ///< xyz coordinates of element ele
      const M3& grad_ephi_smoothed, ///< smoothed nodal G-function gradient values for this element
      const M4& grad_ephi2_smoothed,///< smoothed nodal G-function second derivatives for this element
      bool set_ephi2_smoothed = false ///< use grad_ephi2_smoothed
  )
  {
    // number space dimensions for 3d combustion element
    const size_t nsd = 3;

    // number of nodes of element
    const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    //==========================reconstruct normal_phi = grad(phi) / |grad(phi)| ===========

    static LINALG::Matrix<nsd,1> grad_phi;
    grad_phi.Clear();

    // get shape functions at gaussian point
    static LINALG::Matrix<numnode,1> funct_gp;
    funct_gp.Clear();
    DRT::UTILS::shape_function_3D(funct_gp,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

    // TODO: already defined above
    static LINALG::Matrix<nsd,numnode> deriv_gp;
    deriv_gp.Clear();
    DRT::UTILS::shape_function_3D_deriv1(deriv_gp,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);

    static LINALG::Matrix<nsd,nsd> xjm_gp;
    xjm_gp.Clear();
    xjm_gp.MultiplyNT(deriv_gp,xyze);

    static LINALG::Matrix<nsd,nsd> xji_gp;
    xji_gp.Clear();
    xji_gp.Invert(xjm_gp);

    // derivatives in xyz-direction
    static LINALG::Matrix<nsd,numnode> deriv_gp_xyz;
    deriv_gp_xyz.Clear();
    deriv_gp_xyz.Multiply(xji_gp,deriv_gp);


    // get second derivatives of phi
    static LINALG::Matrix<9,1> grad_phi2;
    grad_phi2.Clear();

    // loop over nodes
    for(size_t i = 0; i< numnode; i++)
    {
      LINALG::Matrix<nsd,1> nodal_grad_xyz_tmp;
      //      LINALG::Matrix<nsd,1> nodal_grad_xi_tmp;
      nodal_grad_xyz_tmp.Clear();
      //      nodal_grad_xi_tmp.Clear();

      // transform xyz-gradient
      // 1. transformt nodal xyz-gradient (nodal_grad_xyz_tmp) to nodal xi-gradient => nodal_grad_xi_tmp
      // 2. interpolate xi-gradient at Gaussian point          => nodal_grad_xi
      // 3. transform xi-gradient at Gaussian point to xyz-gradient
      // 4. normalize final gradient

      // 1. transform nodal xyz-gradient to nodal xi-gradient
      // get smoothed xyz-gradient at node i
      nodal_grad_xyz_tmp(0) =  grad_ephi_smoothed(0,i);
      nodal_grad_xyz_tmp(1) =  grad_ephi_smoothed(1,i);
      nodal_grad_xyz_tmp(2) =  grad_ephi_smoothed(2,i);

      grad_phi(0) += nodal_grad_xyz_tmp(0) *funct_gp(i);
      grad_phi(1) += nodal_grad_xyz_tmp(1) *funct_gp(i);
      grad_phi(2) += nodal_grad_xyz_tmp(2) *funct_gp(i);

      if (not set_ephi2_smoothed)
      {
        grad_phi2(0) += deriv_gp_xyz(0,i)*nodal_grad_xyz_tmp(0); // ,xx
        grad_phi2(1) += deriv_gp_xyz(1,i)*nodal_grad_xyz_tmp(1); // ,yy
        grad_phi2(2) += deriv_gp_xyz(2,i)*nodal_grad_xyz_tmp(2); // ,zz
        grad_phi2(3) += deriv_gp_xyz(1,i)*nodal_grad_xyz_tmp(0); // ,xy
        grad_phi2(4) += deriv_gp_xyz(2,i)*nodal_grad_xyz_tmp(0); // ,xz
        grad_phi2(5) += deriv_gp_xyz(2,i)*nodal_grad_xyz_tmp(1); // ,yz
        grad_phi2(6) += deriv_gp_xyz(0,i)*nodal_grad_xyz_tmp(1); // ,yx
        grad_phi2(7) += deriv_gp_xyz(0,i)*nodal_grad_xyz_tmp(2); // ,zx
        grad_phi2(8) += deriv_gp_xyz(1,i)*nodal_grad_xyz_tmp(2); // ,zy
      }
      else
      {
        grad_phi2(0) += grad_ephi2_smoothed(0,i)*funct_gp(i); // ,xx
        grad_phi2(1) += grad_ephi2_smoothed(1,i)*funct_gp(i); // ,yy
        grad_phi2(2) += grad_ephi2_smoothed(2,i)*funct_gp(i); // ,zz
        grad_phi2(3) += grad_ephi2_smoothed(3,i)*funct_gp(i); // ,xy
        grad_phi2(4) += grad_ephi2_smoothed(4,i)*funct_gp(i); // ,xz
        grad_phi2(5) += grad_ephi2_smoothed(5,i)*funct_gp(i); // ,yz
        grad_phi2(6) += grad_ephi2_smoothed(6,i)*funct_gp(i); // ,yx
        grad_phi2(7) += grad_ephi2_smoothed(7,i)*funct_gp(i); // ,zx
        grad_phi2(8) += grad_ephi2_smoothed(8,i)*funct_gp(i); // ,zy
      }
    }

    // 3. transform xi-gradient at Gaussian point to xyz-gradient
    //      normal_Phi.Multiply(xji_gp,grad_xi);

    // 4. normalize final gradient
    //normal_Phi.Scale(1.0/normal_Phi.Norm2());

    //=========================================== calculate curvature =========================

    double grad_phi_norm = grad_phi.Norm2();

    // check norm of normal gradient
    if (fabs(grad_phi_norm) < 1.0E-5)// 'ngradnorm' == 0.0
    {
      //std::cout << "grad phi is small -> set to 1.0E12" << grad_phi_norm << std::endl;
      // phi gradient too small -> there must be a local max or min in the level-set field
      // set curvature to a large value (it will be cut off based on the element size)
      curvature = 1.0E12;
    }
    else
    {
      double val = grad_phi_norm*grad_phi_norm*grad_phi_norm;
      double invval = 1.0 / val;
      curvature = -invval*(  grad_phi(0)*grad_phi(0)*grad_phi2(0)
          + grad_phi(1)*grad_phi(1)*grad_phi2(1)
          + grad_phi(2)*grad_phi(2)*grad_phi2(2)  )
          -invval*(  grad_phi(0)*grad_phi(1)*( grad_phi2(3) + grad_phi2(6) )
              + grad_phi(0)*grad_phi(2)*( grad_phi2(4) + grad_phi2(7) )
              + grad_phi(1)*grad_phi(2)*( grad_phi2(5) + grad_phi2(8)) )
              +1.0/grad_phi_norm * ( grad_phi2(0) + grad_phi2(1) + grad_phi2(2) );
    }

    return;

  }


  /*
  \brief calculate rate of strain of (fine-scale) velocity
  \param velderxy   (o) velocity derivatives
 \return computed rate of strain
 */
template <DRT::Element::DiscretizationType DISTYPE>
    double GetStrainRate(const LINALG::Matrix<3,3>& velderxy
                     )
{
  double rateofstrain=0.0;
  const size_t nsd=3;

  // compute (resolved) rate of strain
  //
  //          +-                                 -+ 1
  //          |          /   \           /   \    | -
  //          | 2 * eps | vel |   * eps | vel |   | 2
  //          |          \   / ij        \   / ij |
  //          +-                                 -+
  //

  static LINALG::Matrix<nsd,nsd> two_epsilon;
  for(unsigned rr=0;rr<nsd;++rr)
  {
    for(unsigned mm=0;mm<nsd;++mm)
    {
      two_epsilon(rr,mm) = velderxy(rr,mm) + velderxy(mm,rr);
    }
  }

  for(unsigned rr=0;rr<nsd;++rr)
  {
    for(unsigned mm=0;mm<nsd;++mm)
    {
      rateofstrain += two_epsilon(rr,mm)*two_epsilon(mm,rr);
    }
  }

  return(sqrt(rateofstrain/2.0));
}

} // namespace COMBUST

#endif
