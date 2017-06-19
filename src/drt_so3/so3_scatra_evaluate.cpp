/*!----------------------------------------------------------------------
\file so3_scatra_evaluate.cpp

\brief Solid-scatra elements evaluate

\level 2

<pre>
   \maintainer Thon Moritz
               thon@mhpc.mw.tum.de
               089 - 289-10264
</pre>

*----------------------------------------------------------------------*/

#include "so3_scatra.H"
#include "so_base.H"

#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"

#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lib/drt_element_integration_select.H"

/*----------------------------------------------------------------------*
 |  preevaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Scatra<so3_ele,distype>::PreEvaluate(Teuchos::ParameterList& params,
                                        DRT::Discretization&      discretization,
                                        DRT::Element::LocationArray& la)
{

  if(la.Size()>1)
  {
    //ask for the number of dofs of second dofset (scatra)
    const int numscal = discretization.NumDof(1,Nodes()[0]);

    if (la[1].Size() != numnod_*numscal)
      dserror("calc_struct_nlnstiff: Location vector length for concentrations does not match!");

    if (discretization.HasState(1,"temperature")) //if concentrations were set
    {
      //get reference coordinated of current element
      DRT::Node** nodes = Nodes();
      LINALG::Matrix<numnod_,numdim_> xrefe;  // reference coord. of element

      for (int k=0; k<numnod_; ++k)
      {
        const double* x = nodes[k]->X();
        xrefe(k,0) = x[0];
        xrefe(k,1) = x[1];
        xrefe(k,2) = x[2];
      }

      std::vector<LINALG::Matrix<numnod_,1> > shapefunct(numgpt_);
      std::vector<LINALG::Matrix<numdim_,numnod_> > derivs(numgpt_);
      std::vector<double> detJref(numgpt_);

      if (not (distype == DRT::Element::hex8 or distype == DRT::Element::hex27))
        dserror("The Solidscatra elements are only tested for the Hex8 case. The following should work, but keep you eyes open (especially with the order of the gauß points");

      if (distype == DRT::Element::hex8)
      {
        // (r,s,t) gp-locations of fully integrated linear 8-node Hex
        const double gploc    = 1.0/sqrt(3.0);    // gp sampling point value for linear fct
        const double r[NUMGPT_SOH8] = {-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc,-gploc};
        const double s[NUMGPT_SOH8] = {-gploc,-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc};
        const double t[NUMGPT_SOH8] = {-gploc,-gploc,-gploc,-gploc, gploc, gploc, gploc, gploc};
        // fill up nodal f at each gp
        for (int igp=0; igp<numgpt_; ++igp)
        {
          (shapefunct[igp])(0) = (1.0-r[igp])*(1.0-s[igp])*(1.0-t[igp])*0.125;
          (shapefunct[igp])(1) = (1.0+r[igp])*(1.0-s[igp])*(1.0-t[igp])*0.125;
          (shapefunct[igp])(2) = (1.0+r[igp])*(1.0+s[igp])*(1.0-t[igp])*0.125;
          (shapefunct[igp])(3) = (1.0-r[igp])*(1.0+s[igp])*(1.0-t[igp])*0.125;
          (shapefunct[igp])(4) = (1.0-r[igp])*(1.0-s[igp])*(1.0+t[igp])*0.125;
          (shapefunct[igp])(5) = (1.0+r[igp])*(1.0-s[igp])*(1.0+t[igp])*0.125;
          (shapefunct[igp])(6) = (1.0+r[igp])*(1.0+s[igp])*(1.0+t[igp])*0.125;
          (shapefunct[igp])(7) = (1.0-r[igp])*(1.0+s[igp])*(1.0+t[igp])*0.125;

          // df wrt to r for each node(0..7) at each gp [igp]
          (derivs[igp])(0,0) = -(1.0-s[igp])*(1.0-t[igp])*0.125;
          (derivs[igp])(0,1) =  (1.0-s[igp])*(1.0-t[igp])*0.125;
          (derivs[igp])(0,2) =  (1.0+s[igp])*(1.0-t[igp])*0.125;
          (derivs[igp])(0,3) = -(1.0+s[igp])*(1.0-t[igp])*0.125;
          (derivs[igp])(0,4) = -(1.0-s[igp])*(1.0+t[igp])*0.125;
          (derivs[igp])(0,5) =  (1.0-s[igp])*(1.0+t[igp])*0.125;
          (derivs[igp])(0,6) =  (1.0+s[igp])*(1.0+t[igp])*0.125;
          (derivs[igp])(0,7) = -(1.0+s[igp])*(1.0+t[igp])*0.125;

          // df wrt to s for each node(0..7) at each gp [igp]
          (derivs[igp])(1,0) = -(1.0-r[igp])*(1.0-t[igp])*0.125;
          (derivs[igp])(1,1) = -(1.0+r[igp])*(1.0-t[igp])*0.125;
          (derivs[igp])(1,2) =  (1.0+r[igp])*(1.0-t[igp])*0.125;
          (derivs[igp])(1,3) =  (1.0-r[igp])*(1.0-t[igp])*0.125;
          (derivs[igp])(1,4) = -(1.0-r[igp])*(1.0+t[igp])*0.125;
          (derivs[igp])(1,5) = -(1.0+r[igp])*(1.0+t[igp])*0.125;
          (derivs[igp])(1,6) =  (1.0+r[igp])*(1.0+t[igp])*0.125;
          (derivs[igp])(1,7) =  (1.0-r[igp])*(1.0+t[igp])*0.125;

          // df wrt to t for each node(0..7) at each gp [igp]
          (derivs[igp])(2,0) = -(1.0-r[igp])*(1.0-s[igp])*0.125;
          (derivs[igp])(2,1) = -(1.0+r[igp])*(1.0-s[igp])*0.125;
          (derivs[igp])(2,2) = -(1.0+r[igp])*(1.0+s[igp])*0.125;
          (derivs[igp])(2,3) = -(1.0-r[igp])*(1.0+s[igp])*0.125;
          (derivs[igp])(2,4) =  (1.0-r[igp])*(1.0-s[igp])*0.125;
          (derivs[igp])(2,5) =  (1.0+r[igp])*(1.0-s[igp])*0.125;
          (derivs[igp])(2,6) =  (1.0+r[igp])*(1.0+s[igp])*0.125;
          (derivs[igp])(2,7) =  (1.0-r[igp])*(1.0+s[igp])*0.125;

          LINALG::Matrix<numdim_,numdim_> invJ;
          invJ.Multiply(derivs[igp],xrefe);
          detJref[igp] = invJ.Determinant()*1.0;
        }
      }
      else //all other elements use the standard shape function from DRT::UTILS::shape_function_3D
      {
        const DRT::UTILS::GaussRule3D gaussrule = DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule;
        const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
        for (int igp = 0; igp < numgpt_; ++igp)
        {
          const double r = intpoints.qxg[igp][0];
          const double s = intpoints.qxg[igp][1];
          const double t = intpoints.qxg[igp][2];

          //get shape functions evaluated at current gauß point
          DRT::UTILS::shape_function_3D(shapefunct[igp], r, s, t, distype);

          //get first derivative of shape functions evaluated at current gauß point
          DRT::UTILS::shape_function_3D_deriv1(derivs[igp], r, s, t, distype);

          LINALG::Matrix<numdim_,numdim_> invJ;
          invJ.Multiply(derivs[igp],xrefe);
          detJref[igp] = invJ.Determinant()*intpoints.qwgt[igp];
        }
      }

      /* =========================================================================*/
      // start concentration business
      /* =========================================================================*/
      Teuchos::RCP<std::vector<std::vector<double> > > gpconc =
          Teuchos::rcp(new std::vector<std::vector<double> >(numgpt_,std::vector<double>(numscal,0.0)));

      // check if you can get the scalar state
      Teuchos::RCP<const Epetra_Vector> concnp = discretization.GetState(1,"temperature");

      if (concnp==Teuchos::null)
        dserror("calc_struct_nlnstiff: Cannot get state vector 'temperature' ");

      // extract local values of the global vectors
      Teuchos::RCP<std::vector<double> > myconc = Teuchos::rcp(new std::vector<double>(la[1].lm_.size(),0.0) );

      DRT::UTILS::ExtractMyValues(*concnp,*myconc,la[1].lm_);

      //element vector for k-th scalar
      std::vector<LINALG::Matrix<numnod_,1> > econc(numscal);
      for (int k=0; k<numscal; ++k)
      {
        for (int i=0; i<numnod_; ++i)
        {
          (econc.at(k))(i,0) = myconc->at(numscal*i+k);
        }
      }

      /* =========================================================================*/
      /* ================================================= Loop over Gauss Points */
      /* =========================================================================*/
      // volume of current element in reference configuration
      double volume_ref = 0.0;
      // mass in current element in reference configuration
      std::vector<double> mass_ref(numscal,0.0);

      for (int igp=0; igp<numgpt_; ++igp)
      {
        volume_ref+=detJref[igp];

        //concentrations at current gauß point
        std::vector<double> conc_gp_k(numscal,0.0);

        // shape functions evaluated at current gauß point
        const LINALG::Matrix<numnod_,1> shapefunct_gp = shapefunct.at(igp);

        for (int k=0; k<numscal; ++k)
        {
          // identical shapefunctions for displacements and temperatures
          conc_gp_k.at(k) = shapefunct_gp.Dot(econc.at(k));

          mass_ref.at(k) += conc_gp_k.at(k)*detJref[igp];
        }

        gpconc->at(igp)=conc_gp_k;
      }

      params.set< Teuchos::RCP<std::vector<std::vector<double> > > >("gp_conc",gpconc);

      //compute average concentrations
      for (int k=0; k<numscal; ++k)
      {
        //now mass_ref is the element averaged concentration
        mass_ref.at(k) /= volume_ref;
      }
      Teuchos::RCP<std::vector<std::vector<double> > > avgconc =
                Teuchos::rcp(new std::vector<std::vector<double> >(numgpt_,mass_ref));

      params.set< Teuchos::RCP<std::vector<std::vector<double> > > >("avg_conc",avgconc);

    }

    //If you need a pointer to the scatra material, use these lines:
    //we assume that the second material of the structure is the scatra element material
    //Teuchos::RCP<MAT::Material> scatramat = so3_ele::Material(1);
    //params.set< Teuchos::RCP<MAT::Material> >("scatramat",scatramat);
  }

  //TODO: (thon) actually we do not want this here, since it has nothing to do with scatra specific stuff. But for now we let it be...
  std::vector<double> center = DRT::UTILS::ElementCenterRefeCoords(this);
  Teuchos::RCP<std::vector<double> >xrefe = Teuchos::rcp(new std::vector<double>(center));
  params.set<Teuchos::RCP<std::vector<double> > >("position",xrefe);

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Scatra< so3_ele, distype>::Evaluate(Teuchos::ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    DRT::Element::LocationArray& la,
                                    Epetra_SerialDenseMatrix& elemat1_epetra,
                                    Epetra_SerialDenseMatrix& elemat2_epetra,
                                    Epetra_SerialDenseVector& elevec1_epetra,
                                    Epetra_SerialDenseVector& elevec2_epetra,
                                    Epetra_SerialDenseVector& elevec3_epetra)
{

  PreEvaluate(params,
              discretization,
              la);

  return so3_ele::Evaluate(params,
                    discretization,
                    la[0].lm_,
                    elemat1_epetra,
                    elemat2_epetra,
                    elevec1_epetra,
                    elevec2_epetra,
                    elevec3_epetra);

}


template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8,DRT::Element::hex8>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex27,DRT::Element::hex27>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8fbar,DRT::Element::hex8>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4,DRT::Element::tet4>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet10,DRT::Element::tet10>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6,DRT::Element::wedge6>;

