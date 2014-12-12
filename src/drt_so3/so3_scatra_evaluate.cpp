/*!----------------------------------------------------------------------
\file So3_scatra_evaluate.cpp

<pre>
   Maintainer: Cristobal Bertoglio
               bertoglio@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
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
    //ask for the number ofs dofs of second dofset (fluid)
    const int numdofpernode = discretization.NumDof(1,Nodes()[0]);

    if (la[1].Size() != numnod_*numdofpernode)
      dserror("calc_struct_nlnstiff: Location vector length for velocities does not match!");

    if (discretization.HasState(1,"temperature"))
    {
      // check if you can get the scalar state
      Teuchos::RCP<const Epetra_Vector> tempnp
        = discretization.GetState(1,"temperature");

      if (tempnp==Teuchos::null)
        dserror("calc_struct_nlnstiff: Cannot get state vector 'temperature' ");

      // extract local values of the global vectors
      Teuchos::RCP<std::vector<double> >mytemp = Teuchos::rcp(new std::vector<double>(la[1].lm_.size()) );
      DRT::UTILS::ExtractMyValues(*tempnp,*mytemp,la[1].lm_);

      //number of scalars
      int numscal=mytemp->size()/numnod_;

      Teuchos::RCP<std::vector<double> > meantemp = Teuchos::rcp(new std::vector<double>(numscal));

      //calculate for each scalar the mean scalar value of this element
      for (int i=0; i<numscal;i++)
      {
        double meantempi = 0.0;
        for (int j=0; j<numnod_; ++j){
            meantempi +=  (*mytemp)[i+numscal*j];
        }
        meantemp->at(i) =meantempi/numnod_;
      }

     params.set< Teuchos::RCP<std::vector<double> > >("mean_concentrations",meantemp);
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
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8fbar,DRT::Element::hex8>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4,DRT::Element::tet4>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet10,DRT::Element::tet10>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6,DRT::Element::wedge6>;

