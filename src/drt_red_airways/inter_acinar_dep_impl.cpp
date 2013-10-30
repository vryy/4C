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
  double N0 = double(ele->Nodes()[0]->NumElement());
  double N1 = double(ele->Nodes()[1]->NumElement());
  elemat1_epetra(0,0) =  1.0/(N0-1.0);
  elemat1_epetra(0,1) = -1.0/(N0-1.0);
  elemat1_epetra(1,0) = -1.0/(N1-1.0);
  elemat1_epetra(1,1) =  1.0/(N1-1.0);
  elevec1_epetra.Scale(0.0);

  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)  ismail 01/10|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::InterAcinarDepImpl<distype>::Initial(
  RedInterAcinarDep*                     ele,
  Teuchos::ParameterList&                params,
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
 |  Evaluate the values of the degrees of freedom           ismail 04/13|
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
  const int   myrank  = discretization.Comm().MyPID();

  // get total time
  const double time = params.get<double>("total time");

  // get time-step size
  //  const double dt = params.get<double>("time step size");

  // the number of nodes
  const int numnode = lm.size();
  std::vector<int>::iterator it_vcr;

  RCP<const Epetra_Vector> pnp  = discretization.GetState("pnp");

  if (pnp==Teuchos::null)
    dserror("Cannot get state vectors 'pnp'");

  // extract local values from the global vectors
  std::vector<double> mypnp(lm.size());
  DRT::UTILS::ExtractMyValues(*pnp,mypnp,lm);

  // create objects for element arrays
  Epetra_SerialDenseVector epnp(numnode);

  //get time step size
  //  const double dt = params.get<double>("time step size");

  //get all values at the last computed time step
  for (int i=0;i<numnode;++i)
  {
    // split area and volumetric flow rate, insert into element arrays
    epnp(i)    = mypnp[i];
  }

  // ---------------------------------------------------------------------------------
  // Resolve the BCs
  // ---------------------------------------------------------------------------------

  for(int i = 0; i<ele->NumNode(); i++)
  {
    if (ele->Nodes()[i]->Owner()== myrank)
    {
      if(ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond") )
      {
        std::string Bc;
        double BCin = 0.0;
        if (ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond"))
        {
          DRT::Condition * condition = ele->Nodes()[i]->GetCondition("RedAirwayPrescribedCond");
          // Get the type of prescribed bc
          Bc = *(condition->Get<std::string>("boundarycond"));


          const  std::vector<int>*    curve  = condition->Get<std::vector<int>    >("curve");
          double curvefac = 1.0;
          const  std::vector<double>* vals   = condition->Get<std::vector<double> >("val");
          const std::vector<int>*     functions = condition->Get<std::vector<int> >("funct");

          // -----------------------------------------------------------------
          // Read in the value of the applied BC
          // -----------------------------------------------------------------
          if((*curve)[0]>=0)
          {
            curvefac = DRT::Problem::Instance()->Curve((*curve)[0]).f(time);
            BCin = (*vals)[0]*curvefac;
          }
          else
          {
            dserror("no boundary condition defined!");
            exit(1);
          }

          int functnum = -1;
          if (functions) functnum = (*functions)[0];
          else functnum = -1;

          double functionfac = 0.0;
          if(functnum>0)
          {
            functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(0,(ele->Nodes()[i])->X(),time,NULL);
          }
          BCin += functionfac;

          // -----------------------------------------------------------------------------
          // get the local id of the node to whome the bc is prescribed
          // -----------------------------------------------------------------------------
          int local_id =  discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
          if (local_id< 0 )
          {
            dserror("node (%d) doesn't exist on proc(%d)",ele->Nodes()[i]->Id(),discretization.Comm().MyPID());
            exit(1);
          }
        }
        else
        {

        }

        if (Bc == "pressure" || Bc == "ExponentialPleuralPressure")
        {
          RCP<Epetra_Vector> bcval  = params.get<RCP<Epetra_Vector> >("bcval");
          RCP<Epetra_Vector> dbctog = params.get<RCP<Epetra_Vector> >("dbctog");

          if (bcval==null||dbctog==Teuchos::null)
          {
            dserror("Cannot get state vectors 'bcval' and 'dbctog'");
            exit(1);
          }

          if (Bc == "ExponentialPleuralPressure")
          {
            const double ap =  -977.203;
            const double bp = -3338.290;
            const double cp =    -7.686;
            const double dp =  2034.470;
            const double VFR   = 1240000.0;
            const double TLC   = 4760000.0 - VFR;

            const double lungVolumenp = params.get<double>("lungVolume_np") - VFR;

            const double TLCnp= lungVolumenp/TLC;

            double Pp_np = ap + bp*exp(cp*TLCnp) + dp*TLCnp;

            BCin += Pp_np;
          }

          // set pressure at node i
          int    gid;
          double val;

          gid = lm[i];
          val = BCin;
          bcval->ReplaceGlobalValues(1,&val,&gid);

          gid = lm[i];
          val = 1;
          dbctog->ReplaceGlobalValues(1,&val,&gid);

        }
        else
        {
          dserror("precribed [%s] is not defined for reduced-inter-acinar-linkers",Bc.c_str());
          exit(1);
        }

      }
      else
      {
        // ---------------------------------------------------------------
        // If the node is a terminal node, but no b.c is prescribed to it
        // then a zero output pressure is assumed
        // ---------------------------------------------------------------
        if (ele->Nodes()[i]->NumElement() == 1)
        {
          // -------------------------------------------------------------
          // get the local id of the node to whome the bc is prescribed
          // -------------------------------------------------------------

          int local_id =  discretization.NodeRowMap()->LID(ele->Nodes()[i]->Id());
          if (local_id< 0 )
          {
            dserror("node (%d) doesn't exist on proc(%d)",ele->Nodes()[i],discretization.Comm().MyPID());
            exit(1);
          }

          RCP<Epetra_Vector> bcval  = params.get<RCP<Epetra_Vector> >("bcval");
          RCP<Epetra_Vector> dbctog = params.get<RCP<Epetra_Vector> >("dbctog");

          if (bcval==null||dbctog==Teuchos::null)
          {
            dserror("Cannot get state vectors 'bcval' and 'dbctog'");
            exit(1);
          }


          // set pressure at node i
          int    gid;
          double val;

          gid = lm[i];
          val = 0.0;
          bcval->ReplaceGlobalValues(1,&val,&gid);

          gid = lm[i];
          val = 1;
          dbctog->ReplaceGlobalValues(1,&val,&gid);
        }
      } // END of if there is no BC but the node still is at the terminal
    } // END of if node is available on this processor
  } // End of node i has a condition
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
  RedInterAcinarDep*      ele,
  Teuchos::ParameterList& params,
  DRT::Discretization&    discretization,
  std::vector<int>&       lm,
  RCP<MAT::Material>      material)
{

}

