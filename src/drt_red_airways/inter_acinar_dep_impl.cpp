/*----------------------------------------------------------------------*/
/*!
\file inter_acinar_dep_impl.cpp

\brief Internal implementation of RedAcinus element

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>
*/
/*----------------------------------------------------------------------*/



#include "acinus_impl.H"
#include "inter_acinar_dep_impl.H"


#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/maxwell_0d_acinus.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include <fstream>
#include <iomanip>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedInterAcinarDepImplInterface* DRT::ELEMENTS::RedInterAcinarDepImplInterface::Impl(DRT::ELEMENTS::RedInterAcinarDep* red_acinus)
{
  switch (red_acinus->Shape())
  {
  case DRT::Element::line2:
  {
    static InterAcinarDepImpl<DRT::Element::line2>* acinus;
    if (acinus==NULL)
    {
      acinus = new InterAcinarDepImpl<DRT::Element::line2>;
    }
    return acinus;
  }
  default:
    dserror("shape %d (%d nodes) not supported", red_acinus->Shape(), red_acinus->NumNode());
  }
  return NULL;
}



/*----------------------------------------------------------------------*
  | constructor (public)                                    ismail 01/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::InterAcinarDepImpl<distype>::InterAcinarDepImpl()
{

}

/*----------------------------------------------------------------------*
 | evaluate (public)                                       ismail 01/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::InterAcinarDepImpl<distype>::Evaluate(
  RedInterAcinarDep*         ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       discretization,
  std::vector<int>&          lm,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseMatrix&  elemat2_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseVector&  elevec2_epetra,
  Epetra_SerialDenseVector&  elevec3_epetra,
  RCP<MAT::Material> mat)
{

#if 0
  const int   myrank  = discretization.Comm().MyPID();
  // get then nodes connected to this element
  DRT::Node** nodes = ele->Nodes();

  for (int i=0;i<2;i++)
  {
    // get the acinar element connected to the ith (i=0;1) nodes
    DRT::Element ** elems = (nodes[i])->Elements();

    int numOfElems = (nodes[i])->NumElement();
    //RedAirwayType::Instance().UniqueParObjectId()
    for (int j=0;j<numOfElems;j++)
    {
      int generation;

      DRT::ELEMENTS::RedAcinus * actele = dynamic_cast<DRT::ELEMENTS::RedAcinus*>(elems[j]);
      if (actele)
      {
        if (actele->Owner()!=myrank)
        {
          cout<<"+-------------- WOOOOOOW --------------+"<<endl;
        }
        cout<<"Element ("<<actele->Id()<<") rank("<<myrank<<") and HasGlobalElement("<<discretization.HaveGlobalElement(actele->Id())<<")"<<endl;

        vector<int> la;
        // get element location vector, dirichlet flags and ownerships
        discretization.Dof(elems[j],la);
        //        elems[j]->LocationVector(discretization,la,false);

        Epetra_SerialDenseMatrix  elemat_epetra(elemat1_epetra);

        (elems[j])->Evaluate(params,
                         discretization,
                         la,
                         elemat_epetra,
                         elemat2_epetra,
                         elevec1_epetra,
                         elevec2_epetra,
                         elevec3_epetra);

        int index = 0;
        if (i==0)
          index = 1;

        for (int k=0;k<elevec1_epetra.Length();k++)
        {
          elemat1_epetra(index,k) = elemat_epetra(index,k);
        }
        elevec1_epetra.Scale(0.0);
      }
    }
  }


  cout<<"sysmat: "<<elemat1_epetra<<endl;
  fflush(stdout);
#else
  RCP<const Epetra_Vector> sysmat_iad  = discretization.GetState("sysmat_iad");
  RCP<const Epetra_Vector> pn  = discretization.GetState("pn");

  // extract local values from the global vectors
  vector<double> my_sysmat_iad(lm.size());
  DRT::UTILS::ExtractMyValues(*sysmat_iad,my_sysmat_iad,lm);
  vector<double> my_pn(lm.size());
  DRT::UTILS::ExtractMyValues(*pn,my_pn,lm);


#if 0
  double N0 = double(ele->Nodes()[0]->NumElement());
  double N1 = double(ele->Nodes()[1]->NumElement());
#endif

  elemat1_epetra(0,0) =  my_sysmat_iad[0];
  elemat1_epetra(0,1) = -my_sysmat_iad[0];
  elemat1_epetra(1,0) = -my_sysmat_iad[1];
  elemat1_epetra(1,1) =  my_sysmat_iad[1];

  elevec1_epetra.Scale(0.0);
  elevec1_epetra(0)= -my_sysmat_iad[0]*(my_pn[0]-my_pn[1]);
  elevec1_epetra(1)= -my_sysmat_iad[1]*(my_pn[1]-my_pn[0]);
#endif
  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::InterAcinarDepImpl<distype>::Initial(
  RedInterAcinarDep*                     ele,
  Teuchos::ParameterList&                         params,
  DRT::Discretization&                   discretization,
  std::vector<int>&                      lm,
  Teuchos::RCP<const MAT::Material>      material)
{

  RCP<Epetra_Vector> generations   = params.get<RCP<Epetra_Vector> >("generations");

  //--------------------------------------------------------------------
  // get the generation numbers
  //--------------------------------------------------------------------
  //  if(myrank == ele->Owner())
  {
    int    gid = ele->Id();
    double val = -2.0;
    generations->ReplaceGlobalValues(1,&val,&gid);
  }

}//InterAcinarDepImpl::Initial

/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::InterAcinarDepImpl<distype>::Sysmat(
  RedInterAcinarDep*                       ele,
  Epetra_SerialDenseVector&                epnp,
  Epetra_SerialDenseVector&                epn,
  Epetra_SerialDenseVector&                epnm,
  Epetra_SerialDenseMatrix&                sysmat,
  Epetra_SerialDenseVector&                rhs,
  Teuchos::RCP<const MAT::Material>        material,
  ParameterList &                          params,
  double                                   time,
  double                                   dt)
{
}

/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 01/10|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::InterAcinarDepImpl<distype>::EvaluateTerminalBC(
  RedInterAcinarDep*           ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  Epetra_SerialDenseVector&    rhs,
  RCP<MAT::Material>   material)
{

}


/*----------------------------------------------------------------------*
 |  Evaluate the values of the degrees of freedom           ismail 01/10|
 |  at terminal nodes.                                                  |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::InterAcinarDepImpl<distype>::CalcFlowRates(
  RedInterAcinarDep*           ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  Epetra_SerialDenseVector&    elevec1, //a_volumenp,
  Epetra_SerialDenseVector&    elevec2, //a_volume_strain_np,
  std::vector<int>&            lm,
  RCP<MAT::Material>   material)

{

}

/*----------------------------------------------------------------------*
 |  Get the coupled the values on the coupling interface    ismail 07/10|
 |  of the 3D/reduced-D problem                                         |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::InterAcinarDepImpl<distype>::GetCoupledValues(
  RedInterAcinarDep*                   ele,
  Teuchos::ParameterList&      params,
  DRT::Discretization&         discretization,
  std::vector<int>&            lm,
  RCP<MAT::Material>   material)
{

}

