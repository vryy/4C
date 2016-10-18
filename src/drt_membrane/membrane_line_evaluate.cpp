/*!----------------------------------------------------------------------
\file membrane_line_evaluate.cpp
\brief

\level 3

<pre>
\maintainer Fabian Bräu
            braeu@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

\brief Nonlinear Membrane Finite Element line evaluation

*----------------------------------------------------------------------*/

#include "membrane.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_timecurve.H"

#include "../drt_structure_new/str_elements_paramsinterface.H"


/*----------------------------------------------------------------------*
 |  Integrate a Line Neumann boundary condition (public)   fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::MembraneLine<distype>::EvaluateNeumann(Teuchos::ParameterList&   params,
                                                          DRT::Discretization&      discretization,
                                                          DRT::Condition&           condition,
                                                          std::vector<int>&         lm,
                                                          Epetra_SerialDenseVector& elevec1,
                                                          Epetra_SerialDenseMatrix* elemat1)
{
  // set params interface pointer in the parent element
  ParentElement()->SetParamsInterfacePtr(params);

  // get type of condition
  enum LoadType
  {
    neum_none,
    neum_live
  };
  LoadType ltype;

  const std::string* type = condition.Get<std::string>("type");
  if (*type == "neum_live")
  {
    ltype  = neum_live;
  }
  else dserror("Unknown type of LineNeumann condition");

  // get values and switches from the condition
  const std::vector<int>*    onoff = condition.Get<std::vector<int> >   ("onoff");
  const std::vector<double>* val   = condition.Get<std::vector<double> >("val"  );
  const std::vector<int>*    spa_func = condition.Get<std::vector<int> >("funct");

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  bool usetime = true;

  double time = -1.0;

  if (ParentElement()->IsParamsInterface())
    time = ParentElement()->ParamsInterfacePtr()->GetTotalTime();
  else
    time = params.get("total time",-1.0);

  if (time<0.0) usetime = false;

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < noddof_)
    dserror("Fewer functions or curves defined than the element has dofs.");

  for (int checkdof = noddof_; checkdof < int(onoff->size()); ++checkdof)
  {
    if ((*onoff)[checkdof] != 0)
      dserror("Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
  }

  // find out whether we will use time curves and get the factors
  const std::vector<int>* curve  = condition.Get<std::vector<int> >("curve");
  std::vector<double> curvefacs(noddof_, 1.0);
  for (int i=0; i < noddof_; ++i)
  {
    const int curvenum = (curve) ? (*curve)[i] : -1;
    if (curvenum>=0 && usetime)
      curvefacs[i] = DRT::Problem::Instance()->Curve(curvenum).f(time);
  }

  // element geometry update - currently only material configuration
  LINALG::Matrix<numnod_line_,noddof_> x(true);
  for (int i=0; i<numnod_line_; ++i)
  {
    x(i,0) = Nodes()[i]->X()[0];
    x(i,1) = Nodes()[i]->X()[1];
    x(i,2) = Nodes()[i]->X()[2];
  }

  // allocate vector for shape functions and matrix for derivatives at gp
  LINALG::Matrix<numnod_line_,1> shapefcts(true);
  LINALG::Matrix<1,numnod_line_> derivs(true);

  // integration
  for (int gp=0; gp<intpointsline_.nquad; ++gp)
  {
    // get gausspoints from integration rule
    double xi_gp = intpointsline_.qxg[gp][0];

    // get gauss weight at current gp
    double gpweight = intpointsline_.qwgt[gp];

    // get shape functions and derivatives in the plane of the element
    DRT::UTILS::shape_function_1D(shapefcts,xi_gp,Shape());
    DRT::UTILS::shape_function_1D_deriv1(derivs,xi_gp,Shape());

    switch(ltype)
    {
    case neum_live:
    {
      // uniform load on reference configuration

      // compute dXYZ / dr
      LINALG::Matrix<noddof_,1> dxyzdr(true);
      dxyzdr.MultiplyTT(1.0,x,derivs,0.0);
      // compute line increment dL
      double dL;
      dL=0.0;
      for (int i=0; i<3; ++i)
      {
        dL += dxyzdr(i)*dxyzdr(i);
      }
      dL = sqrt(dL);

      // loop the dofs of a node
      for (int i = 0; i < noddof_; ++i)
      {
        if ((*onoff)[i]) // is this dof activated?
        {
          // factor given by spatial function
          const int functnum = (spa_func) ? (*spa_func)[i] : -1;
          double functfac = 1.0;

          if (functnum > 0)
          {
            // calculate reference position of GP
            LINALG::Matrix<noddof_,1> gp_coord(true);
            gp_coord.MultiplyTN(1.0,x,shapefcts,0.0);

            // write coordinates in another datatype
            double gp_coord2[noddof_];
            for (int k = 0; k < noddof_; k++)
              gp_coord2[k] = gp_coord(k,0);
            const double* coordgpref = &gp_coord2[0]; // needed for function evaluation

            // evaluate function at current gauss point
            functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(i,coordgpref,time,NULL);
          }

          const double fac = (*val)[i] * gpweight * dL * functfac * curvefacs[i];
          for (int node=0; node < numnod_line_; ++node)
          {
            elevec1[noddof_*node+i] += shapefcts(node) * fac;
          }
        }
      }
      break;
    }

    default:
      dserror("Unknown type of LineNeumann load");
    break;
    }

  }

  return 0;
}

template class DRT::ELEMENTS::MembraneLine<DRT::Element::tri3>;
template class DRT::ELEMENTS::MembraneLine<DRT::Element::tri6>;
template class DRT::ELEMENTS::MembraneLine<DRT::Element::quad4>;
template class DRT::ELEMENTS::MembraneLine<DRT::Element::quad9>;
