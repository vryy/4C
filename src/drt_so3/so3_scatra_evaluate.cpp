/*!----------------------------------------------------------------------
\file so3_scatra_evaluate.cpp

<pre>
   Maintainer: Thon Moritz
               thon@mhpc.mw.tum.de
               089 - 289-10264
</pre>

*----------------------------------------------------------------------*/

#include "so3_scatra.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/micromaterial.H"
#include <iterator>

#include "../drt_inpar/inpar_structure.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_globalproblem.H"

//TODO: (thon) delete this header if the lower to-do is done!
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
//#include "Sacado.hpp"

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
      Teuchos::RCP<std::vector<std::vector<double> > > gpconc =
          Teuchos::rcp(new std::vector<std::vector<double> >(numgpt_,std::vector<double>(numscal,0.0)));

      // check if you can get the scalar state
      Teuchos::RCP<const Epetra_Vector> concnp = discretization.GetState(1,"temperature");

      if (concnp==Teuchos::null)
        dserror("calc_struct_nlnstiff: Cannot get state vector 'temperature' ");

      // extract local values of the global vectors
      Teuchos::RCP<std::vector<double> > myconc = Teuchos::rcp(new std::vector<double>(la[1].lm_.size(),0.0) );

      DRT::UTILS::ExtractMyValues(*concnp,*myconc,la[1].lm_);

      std::vector<LINALG::Matrix<numnod_,1> > shapefunct(numgpt_);

      //TODO: (thon) this must be redone, if all base elements finally use the same integration rule!!
      // like simply: DRT::UTILS::shape_function<distype>(xsi_[gp],shapefunct_gp);
      switch(distype)
      {
      case DRT::Element::hex8:
      {
        // (r,s,t) gp-locations of fully integrated linear 8-node Hex
        const double gploc    = 1.0/sqrt(3.0);    // gp sampling point value for linear fct
        const double r[NUMGPT_SOH8] = {-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc,-gploc};
        const double s[NUMGPT_SOH8] = {-gploc,-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc};
        const double t[NUMGPT_SOH8] = {-gploc,-gploc,-gploc,-gploc, gploc, gploc, gploc, gploc};
        // fill up nodal f at each gp

        for (int i=0; i<NUMGPT_SOH8; ++i)
        {
          (shapefunct[i])(0) = (1.0-r[i])*(1.0-s[i])*(1.0-t[i])*0.125;
          (shapefunct[i])(1) = (1.0+r[i])*(1.0-s[i])*(1.0-t[i])*0.125;
          (shapefunct[i])(2) = (1.0+r[i])*(1.0+s[i])*(1.0-t[i])*0.125;
          (shapefunct[i])(3) = (1.0-r[i])*(1.0+s[i])*(1.0-t[i])*0.125;
          (shapefunct[i])(4) = (1.0-r[i])*(1.0-s[i])*(1.0+t[i])*0.125;
          (shapefunct[i])(5) = (1.0+r[i])*(1.0-s[i])*(1.0+t[i])*0.125;
          (shapefunct[i])(6) = (1.0+r[i])*(1.0+s[i])*(1.0+t[i])*0.125;
          (shapefunct[i])(7) = (1.0-r[i])*(1.0+s[i])*(1.0+t[i])*0.125;
        }

        break;
      }
      case DRT::Element::hex27:
      {
        // (r,s,t) gp-locations of fully integrated quadratic Hex 27
        // fill up nodal f at each gp
        const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_hex_27point;
        const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
        for (int igp = 0; igp < intpoints.nquad; ++igp)
        {
          const double r = intpoints.qxg[igp][0];
          const double s = intpoints.qxg[igp][1];
          const double t = intpoints.qxg[igp][2];

          DRT::UTILS::shape_function_3D(shapefunct[igp], r, s, t, DRT::Element::hex27);
        }
        break;
      }
      case DRT::Element::tet4:
      {
        // There is only one gausspoint, so the loop (and the vector) is not really needed.
        for (int gp=0; gp<NUMGPT_SOTET4; gp++)
        {
          (shapefunct[gp])(0) = 0.25;
          (shapefunct[gp])(1) = 0.25;
          (shapefunct[gp])(2) = 0.25;
          (shapefunct[gp])(3) = 0.25;
        }
      break;
      }
      default:
        dserror("Element type not supported!");
        break;
      }

      /* =========================================================================*/
      /* ================================================= Loop over Gauss Points */
      /* =========================================================================*/
      for (int gp=0; gp<numgpt_; ++gp)
      {
        std::vector<double> conck(numscal,0.0);

        // shape functions evaluated at current gauß point
        LINALG::Matrix<numnod_,1> shapefunct_gp = shapefunct.at(gp);

        for (int k=0; k<numscal; ++k)
        {
          //element vector for k-th scalar
          LINALG::Matrix<numnod_,1> econck(true);

          for (int i=0; i<numnod_; ++i)
          {
            econck(i,0) = myconc->at(numscal*i+k);
          }
          // identical shapefunctions for displacements and temperatures
          conck.at(k) = shapefunct_gp.Dot(econck);
        }

        gpconc->at(gp)=conck;
      }

      params.set< Teuchos::RCP<std::vector<std::vector<double> > > >("gp_conc",gpconc);
    }

    DRT::Problem* problem = DRT::Problem::Instance();
    const std::vector<std::string> disnames = problem->GetDisNames();
    int numfield = 0;
    numfield = problem->NumFields();

    for (int i=0; i<numfield; ++i)
    {
      std::string name = problem->GetDis(disnames[i])->Name();
      if (name =="scatra")
      {
        // Get pointer for scatra material in the same element
        Teuchos::RCP<DRT::Discretization> scatradis = Teuchos::null;
        scatradis = problem->GetDis("scatra");
        DRT::Element* scatraele = scatradis->gElement(Id());
        Teuchos::RCP<MAT::Material> scatramat = Teuchos::rcp_dynamic_cast<MAT::Material>(scatraele->Material());
        params.set< Teuchos::RCP<MAT::Material> >("scatramat",scatramat);
      }
    }
  }

  Teuchos::RCP<std::vector<double> >xrefe = Teuchos::rcp(new std::vector<double>(3));
  DRT::Node** nodes = Nodes();
  for (int i=0; i<numnod_; ++i){
      const double* x = nodes[i]->X();
      (*xrefe)[0] +=  x[0]/numnod_;
      (*xrefe)[1] +=  x[1]/numnod_;
      (*xrefe)[2] +=  x[2]/numnod_;

   }
   params.set<Teuchos::RCP<std::vector<double> > >("position",xrefe);
   return;
}
/*----------------------------------------------------------------------*
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Scatra<so3_ele,distype>::MyEvaluate(Teuchos::ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    DRT::Element::LocationArray& la,
                                    Epetra_SerialDenseMatrix& elemat1_epetra,
                                    Epetra_SerialDenseMatrix& elemat2_epetra,
                                    Epetra_SerialDenseVector& elevec1_epetra,
                                    Epetra_SerialDenseVector& elevec2_epetra,
                                    Epetra_SerialDenseVector& elevec3_epetra)
{

  return 0;
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
  // start with "none"
  typename So3_Scatra::ActionType act = So3_Scatra::none;

  // get the required action
  std::string action = params.get<std::string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_multidofsetcoupling")   act = So3_Scatra::calc_struct_multidofsetcoupling;
  else if (action=="postprocess_stress")   act = So3_Scatra::postprocess_stress;
  else if (action=="calc_struct_update_istep") act = So3_Scatra::calc_struct_update_istep;

  // what should the element do
  switch(act)
  {
  //==================================================================================
  // coupling terms in force-vector and stiffness matrix
  case So3_Scatra::calc_struct_multidofsetcoupling:
  {

    MyEvaluate(params,
                      discretization,
                      la,
                      elemat1_epetra,
                      elemat2_epetra,
                      elevec1_epetra,
                      elevec2_epetra,
                      elevec3_epetra);
  }
  break;
  case So3_Scatra::postprocess_stress:
  {
    so3_ele::Evaluate(params,
                          discretization,
                          la[0].lm_,
                          elemat1_epetra,
                          elemat2_epetra,
                          elevec1_epetra,
                          elevec2_epetra,
                          elevec3_epetra);
  }
  break;
  /*case So3_Scatra::calc_struct_update_istep:
  {
    so3_ele::Evaluate(params,
                      discretization,
                      la[0].lm_,
                      elemat1_epetra,
                      elemat2_epetra,
                      elevec1_epetra,
                      elevec2_epetra,
                      elevec3_epetra);
  }
  break;*/
  //==================================================================================
  default:
  {
    //in some cases we need to write/change some data before evaluating

    PreEvaluate(params,
                      discretization,
                      la);

    so3_ele::Evaluate(params,
                      discretization,
                      la[0].lm_,
                      elemat1_epetra,
                      elemat2_epetra,
                      elevec1_epetra,
                      elevec2_epetra,
                      elevec3_epetra);

    MyEvaluate(params,
                      discretization,
                      la,
                      elemat1_epetra,
                      elemat2_epetra,
                      elevec1_epetra,
                      elevec2_epetra,
                      elevec3_epetra);

    break;
  }
  } // action

  return 0;
}


template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8,DRT::Element::hex8>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex27,DRT::Element::hex27>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8fbar,DRT::Element::hex8>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4,DRT::Element::tet4>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet10,DRT::Element::tet10>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6,DRT::Element::wedge6>;

