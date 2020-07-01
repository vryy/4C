/*----------------------------------------------------------------------------*/
/*! \file
\brief Evaluation on the boundary lines of a wall1 element.

\level 1


*/
/*---------------------------------------------------------------------------*/

#include "wall1.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_lib/drt_elements_paramsinterface.H"

#include "../drt_mat/structporo.H"

#include "../drt_nurbs_discret/drt_nurbs_discret.H"

/*----------------------------------------------------------------------*
 |  Integrate a Line Neumann boundary condition (public)      popp 06/13|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Wall1Line::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  // set the interface pointer in the parent element
  ParentElement()->SetParamsInterfacePtr(params);
  // IMPORTANT: The 'neum_orthopressure' case represents a truly nonlinear follower-load
  // acting on the spatial configuration. Therefore, it needs to be linearized. On the
  // contrary, the simplified 'neum_pseudo_orthopressure' option allows for an approximative
  // modeling of an orthopressure load without the need to do any linearization. However,
  // this can only be achieved by referring the 'neum_pseudo_orthopressure' load to the last
  // converged configuration, which introduces an error as compared with 'neum_orthopressure'.
  bool loadlin = (elemat1 != NULL);

  // get type of condition
  enum LoadType
  {
    neum_none,
    neum_live,                  // standard Neumann load
    neum_pseudo_orthopressure,  // pseudo-orthopressure load
    neum_orthopressure
    // orthopressure load
  };

  LoadType ltype = neum_none;
  const std::string* type = condition.Get<std::string>("type");
  if (*type == "neum_live")
    ltype = neum_live;
  else if (*type == "neum_pseudo_orthopressure")
    ltype = neum_pseudo_orthopressure;
  else if (*type == "neum_orthopressure")
    ltype = neum_orthopressure;
  else
    dserror("Unknown type of SurfaceNeumann condition");

  // get values and switches from the condition
  const std::vector<int>* onoff = condition.Get<std::vector<int>>("onoff");
  const std::vector<double>* val = condition.Get<std::vector<double>>("val");
  const std::vector<int>* funct = condition.Get<std::vector<int>>("funct");

  // check total time
  double time = -1.0;
  if (ParentElement()->IsParamsInterface())
    time = ParentElement()->ParamsInterfacePtr()->GetTotalTime();
  else
    time = params.get("total time", -1.0);

  // set number of dofs per node
  const int noddof = NumDofPerNode(*Nodes()[0]);

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < noddof)
    dserror("Fewer functions or curves defined than the element has dofs.");

  // set number of nodes
  const int numnod = NumNode();
  const DiscretizationType distype = Shape();

  // gaussian points
  const DRT::UTILS::GaussRule1D gaussrule = getOptimalGaussrule(distype);
  const DRT::UTILS::IntegrationPoints1D intpoints(gaussrule);

  // allocate vector for shape functions and for derivatives
  LINALG::SerialDenseVector shapefcts(numnod);
  LINALG::SerialDenseMatrix deriv(1, numnod);

  // prepare element geometry 1
  // --> we always need the material configuration
  LINALG::SerialDenseMatrix xye(Wall1::numdim_, numnod);
  for (int i = 0; i < numnod; ++i)
  {
    xye(0, i) = Nodes()[i]->X()[0];
    xye(1, i) = Nodes()[i]->X()[1];
  }

  // prepare element geometry 2
  // --> depending on the type of Neumann condition, we might not need a spatial
  // configuration at all (standard Neumann), we might need the last converged
  // spatial position (pseudo-orthopressure) or the true current geometry (orthopressure).
  LINALG::SerialDenseMatrix xyecurr(Wall1::numdim_, numnod);

  // (1) standard Neumann --> we need only material configuration
  if (ltype == neum_live)
  {
    loadlin = false;  // no linearization needed for load in material configuration
  }

  // (2) pseudo orthopressure --> we need last converged configuration
  else if (ltype == neum_pseudo_orthopressure)
  {
    loadlin = false;  // no linearization needed for load in last converged configuration

    Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
    if (disp == Teuchos::null) dserror("Cannot get state vector 'displacement'");
    std::vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

    for (int i = 0; i < numnod; ++i)
    {
      xyecurr(0, i) = xye(0, i) + mydisp[i * noddof + 0];
      xyecurr(1, i) = xye(1, i) + mydisp[i * noddof + 1];
    }
  }

  // (3) true orthopressure --> we need spatial configuration
  else if (ltype == neum_orthopressure)
  {
    if (!loadlin)
      dserror("No linearization provided for orthopressure load (add 'LOADLIN yes' to input file)");

    Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement new");
    if (disp == Teuchos::null) dserror("Cannot get state vector 'displacement new'");
    std::vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

    for (int i = 0; i < numnod; ++i)
    {
      xyecurr(0, i) = xye(0, i) + mydisp[i * noddof + 0];
      xyecurr(1, i) = xye(1, i) + mydisp[i * noddof + 1];
    }
  }

  // loop over integration points //new
  for (int gpid = 0; gpid < intpoints.nquad; gpid++)
  {
    const double e1 = intpoints.qxg[gpid][0];

    // get shape functions and derivatives in the line
    if (distype == line2 || distype == line3)
    {
      DRT::UTILS::shape_function_1D(shapefcts, e1, distype);
      DRT::UTILS::shape_function_1D_deriv1(deriv, e1, distype);
    }
    else if (distype == nurbs2 || distype == nurbs3)
    {
      DRT::NURBS::NurbsDiscretization* nurbsdis =
          dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

      Teuchos::RCP<DRT::NURBS::Knotvector> knots = (*nurbsdis).GetKnotVector();
      std::vector<Epetra_SerialDenseVector> parentknots(2);
      std::vector<Epetra_SerialDenseVector> boundknots(1);

      double normalfac = 0.0;
      bool zero_size = knots->GetBoundaryEleAndParentKnots(
          parentknots, boundknots, normalfac, ParentMasterElement()->Id(), FaceMasterNumber());

      if (zero_size) return (0);

      Epetra_SerialDenseVector weights(NumNode());
      for (int inode = 0; inode < NumNode(); ++inode)
      {
        DRT::NURBS::ControlPoint* cp = dynamic_cast<DRT::NURBS::ControlPoint*>(Nodes()[inode]);
        weights(inode) = cp->W();
      }

      DRT::NURBS::UTILS::nurbs_get_1D_funct_deriv(
          shapefcts, deriv, e1, boundknots[0], weights, distype);
    }
    else
      dserror("Wrong distype!");

    switch (ltype)
    {
      case neum_live:
      {  // uniform load on reference configuration

        // compute infinitesimal line element dr for integration along the line
        const double dr = w1_substitution(xye, deriv, NULL, numnod);

        double functfac = 1.0;

        // loop the dofs of a node
        for (int i = 0; i < noddof; ++i)
        {
          if ((*onoff)[i])  // is this dof activated?
          {
            // factor given by spatial function
            const int functnum = (funct) ? (*funct)[i] : -1;

            if (functnum > 0)
            {
              // calculate reference position of GP
              LINALG::SerialDenseMatrix gp_coord(1, Wall1::numdim_);
              gp_coord.Multiply('T', 'T', 1.0, shapefcts, xye, 0.0);

              // write coordinates in another datatype
              double gp_coord2[3];  // the position vector has to be given in 3D!!!
              const int numdim = 2;
              for (int k = 0; k < numdim; k++) gp_coord2[k] = gp_coord(0, k);
              for (int k = numdim; k < 3;
                   k++)  // set a zero value for the remaining spatial directions
                gp_coord2[k] = 0.0;
              const double* coordgpref = &gp_coord2[0];  // needed for function evaluation

              // evaluate function at current gauss point
              functfac =
                  DRT::Problem::Instance()->Funct(functnum - 1).Evaluate(i, coordgpref, time);
            }
            else
              functfac = 1.0;

            const double fac = intpoints.qwgt[gpid] * dr * (*val)[i] * functfac;
            for (int node = 0; node < numnod; ++node)
            {
              elevec1[node * noddof + i] += shapefcts[node] * fac;
            }
          }
        }
        break;
      }


      case neum_pseudo_orthopressure:  // pseudo-orthogonal pressure on last converged config.
      case neum_orthopressure:
      {  // orthogonal pressure (nonlinear load) on current config.

        // check for correct input
        if ((*onoff)[0] != 1) dserror("orthopressure on 1st dof only!");
        for (int checkdof = 1; checkdof < noddof; ++checkdof)
        {
          if ((*onoff)[checkdof] != 0) dserror("orthopressure on 1st dof only!");
        }
        double ortho_value = (*val)[0];
        if (!ortho_value) dserror("no orthopressure value given!");

        // outward normal vector (unit vector)
        std::vector<double> unrm(Wall1::numdim_);

        // compute infinitesimal line element dr for integration along the line
        const double dr = w1_substitution(xyecurr, deriv, &unrm, numnod);

        // factor given by spatial function
        const int functnum = (funct) ? (*funct)[0] : -1;
        double functfac = 1.0;
        if (functnum > 0)
        {
          // calculate reference position of GP
          LINALG::SerialDenseMatrix gp_coord(1, Wall1::numdim_);
          gp_coord.Multiply('T', 'T', 1.0, shapefcts, xye, 0.0);

          // write coordinates in another datatype
          double gp_coord2[3];  // the position vector has to be given in 3D!!!
          const int numdim = 2;
          for (int k = 0; k < numdim; k++) gp_coord2[k] = gp_coord(0, k);
          for (int k = numdim; k < 3; k++)  // set a zero value for the remaining spatial directions
            gp_coord2[k] = 0.0;
          const double* coordgpref = &gp_coord2[0];  // needed for function evaluation

          // evaluate function at current gauss point
          functfac = DRT::Problem::Instance()->Funct(functnum - 1).Evaluate(0, coordgpref, time);
        }

        // constant factor for integration
        const double fac = intpoints.qwgt[gpid] * ortho_value * functfac;

        // add load components
        for (int node = 0; node < numnod; ++node)
          for (int j = 0; j < noddof; ++j)
            elevec1[node * noddof + j] += shapefcts[node] * unrm[j] * dr * fac;

        // linearization if needed
        if (loadlin)
        {
          // total number of element DOFs
          int numdof = noddof * numnod;

          // directional derivative of surface
          Epetra_SerialDenseMatrix a_Dnormal(Wall1::numdim_, numdof);

          //******************************************************************
          // compute directional derivative
          //******************************************************************

          // linearization of basis vector
          Epetra_SerialDenseMatrix dg(Wall1::numdim_, numdof);
          for (int node = 0; node < numnod; ++node)
            for (int k = 0; k < Wall1::numdim_; ++k) dg(k, node * noddof + k) = deriv(0, node);

          // linearization of local surface normal vector
          for (int dof = 0; dof < numdof; ++dof)
          {
            a_Dnormal(0, dof) = dg(1, dof);
            a_Dnormal(1, dof) = -dg(0, dof);
          }

          // build surface element load linearization matrix
          // (CAREFUL: Minus sign due to the fact that external forces enter the global
          // residual vector with a minus sign, too! However, the load linaerization is
          // simply added to the global tangent stiffness matrix, thus we explicitly
          // need to set the minus sign here.)
          for (int node = 0; node < numnod; ++node)
            for (int dim = 0; dim < 2; dim++)
              for (int dof = 0; dof < elevec1.M(); dof++)
                (*elemat1)(node * noddof + dim, dof) -= shapefcts[node] * a_Dnormal(dim, dof) * fac;
        }

        break;
      }

      default:
      {
        dserror("Unknown type of SurfaceNeumann load");
        break;
      }
    }
  }

  /*// FD CHECK FOR ORTHOPRESSURE
  switch (ltype)
  {
  case neum_orthopressure:
  {
    // prepare FD check
    Epetra_SerialDenseMatrix fd_deriv(numnod*noddof,numnod*noddof);
    Epetra_SerialDenseMatrix an_deriv = *elemat1;
    Epetra_SerialDenseVector eleforce_ref = elevec1;
    Epetra_SerialDenseVector eleforce_curr(numnod*noddof);
    double eps = 1.0e-8;

    // do FD step for all DOFs
    for (int node=0;node<numnod;++node)
      for (int dof=0;dof<noddof;++dof)
      {
        // move position
        xyecurr(dof,node) += eps;

        // actual FD evaluation
        // loop over integration points //new
        for (int gpid = 0; gpid < intpoints.nquad; gpid++)
        {
          const double e1 = intpoints.qxg[gpid][0];

          // get shape functions and derivatives in the line
          DRT::UTILS::shape_function_1D(shapefcts, e1, distype);
          DRT::UTILS::shape_function_1D_deriv1(deriv, e1, distype);
          double ortho_value = (*val)[0];

          // outward normal vector (unit vector)
          std::vector<double> unrm(Wall1::numdim_);

          // compute infinitesimal line element dr for integration along the line
          const double dr = w1_substitution(xyecurr, deriv, &unrm, numnod);

          // constant factor for integration
          const double fac = intpoints.qwgt[gpid] * dr * ortho_value * curvefac;

          // add load components
          for (int k = 0; k < numnod; ++k)
            for (int j = 0; j < noddof; ++j)
              eleforce_curr[k * noddof + j] += funct[k] * unrm[j] * fac;
        }

        // compute FD linearization
        for (int idx=0;idx<numnod*noddof;++idx)
          fd_deriv(idx,node*noddof+dof) = (eleforce_curr[idx]-eleforce_ref[idx])/eps;

        // unmove position
        xyecurr(dof,node) -= eps;

        // reset eleforce_curr
        for (int k = 0; k < numnod; ++k)
          for (int j = 0; j < noddof; ++j)
            eleforce_curr[k * noddof + j] = 0.0;
      }

    // analyze results
    std::cout << "FD LINEARIZATION \n" << fd_deriv << std::endl;
    std::cout << "ANALYTICAL LINEARIZATION \n" << an_deriv << std::endl;
    //exit(0);

    // overwrite results with FD derivative
    for (int j = 0; j < numnod*noddof; ++j)
      for (int k = 0; k < numnod*noddof; ++k)
        (*elemat1)(j,k) = -fd_deriv(j,k);
  break;
  }
  default:
  {
    // do nothing
    break;
  }
  }*/

  return 0;
}

DRT::UTILS::GaussRule1D DRT::ELEMENTS::Wall1Line::getOptimalGaussrule(
    const DiscretizationType& distype)
{
  DRT::UTILS::GaussRule1D rule = DRT::UTILS::intrule1D_undefined;
  switch (distype)
  {
    case line2:
      rule = DRT::UTILS::intrule_line_2point;
      break;
    case line3:
      rule = DRT::UTILS::intrule_line_3point;
      break;
    case nurbs2:
      rule = DRT::UTILS::intrule_line_2point;
      break;
    case nurbs3:
      rule = DRT::UTILS::intrule_line_3point;
      break;
    default:
      dserror("unknown number of nodes for gaussrule initialization");
      break;
  }
  return rule;
}

// determinant of jacobian matrix

double DRT::ELEMENTS::Wall1Line::w1_substitution(const Epetra_SerialDenseMatrix& xye,
    const Epetra_SerialDenseMatrix& deriv,
    std::vector<double>* unrm,  // unit normal
    const int iel)
{
  /*
   |                                            0 1
   |                                           +-+-+
   |       0 1              0...iel-1          | | | 0
   |      +-+-+             +-+-+-+-+          +-+-+
   |      | | | 1     =     | | | | | 0        | | | .
   |      +-+-+             +-+-+-+-+       *  +-+-+ .
   |                                           | | | .
   |                                           +-+-+
   |                                           | | | iel-1
   |                                           +-+-+
   |
   |       dxyzdrs             deriv^T          xye^T
   |
   |
   |                       +-        -+
   |                        | dx   dy  |
   |      yields   dxydr =  | --   --  |
   |                        | dr   dr  |
   |                       +-        -+
   |
   */
  // compute derivative of parametrization
  double dr = 0.0;
  Epetra_SerialDenseMatrix der_par(1, 2);
  int err = der_par.Multiply('N', 'T', 1.0, deriv, xye, 0.0);
  if (err != 0) dserror("Multiply failed");
  dr = sqrt(der_par(0, 0) * der_par(0, 0) + der_par(0, 1) * der_par(0, 1));
  if (unrm != NULL)
  {
    (*unrm)[0] = 1 / dr * der_par(0, 1);
    (*unrm)[1] = -1 / dr * der_par(0, 0);
  }
  return dr;
}

/*======================================================================*/
int DRT::ELEMENTS::Wall1Line::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elematrix1,
    Epetra_SerialDenseMatrix& elematrix2, Epetra_SerialDenseVector& elevector1,
    Epetra_SerialDenseVector& elevector2, Epetra_SerialDenseVector& elevector3)
{
  const DiscretizationType distype = Shape();

  // set number of dofs per node
  const int noddof = NumDofPerNode(*Nodes()[0]);

  // start with "none"
  DRT::ELEMENTS::Wall1Line::ActionType act = Wall1Line::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "calc_struct_constrarea")
    act = Wall1Line::calc_struct_constrarea;
  else if (action == "calc_struct_centerdisp")
    act = Wall1Line::calc_struct_centerdisp;
  else if (action == "calc_struct_areaconstrstiff")
    act = Wall1Line::calc_struct_areaconstrstiff;
  else
    dserror("Unknown type of action for Wall1_Line");
  // create communicator
  const Epetra_Comm& Comm = discretization.Comm();
  // what the element has to do
  switch (act)
  {
    // just compute the enclosed volume (e.g. for initialization)
    case calc_struct_constrarea:
    {
      if (distype != line2)
      {
        dserror("Area Constraint only works for line2 curves!");
      }
      // We are not interested in volume of ghosted elements
      if (Comm.MyPID() == Owner())
      {
        // element geometry update
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp == Teuchos::null) dserror("Cannot get state vector 'displacement'");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
        const int numnod = NumNode();
        Epetra_SerialDenseMatrix xsrefe(numnod, Wall1::numdim_);  // material coord. of element
        Epetra_SerialDenseMatrix xscurr(numnod, Wall1::numdim_);  // material coord. of element
        for (int i = 0; i < numnod; ++i)
        {
          xsrefe(i, 0) = Nodes()[i]->X()[0];
          xsrefe(i, 1) = Nodes()[i]->X()[1];

          xscurr(i, 0) = xsrefe(i, 0) + mydisp[i * noddof];
          xscurr(i, 1) = xsrefe(i, 1) + mydisp[i * noddof + 1];
        }
        // compute area between line and x-Axis
        double areaele = 0.5 * (xscurr(0, 1) + xscurr(1, 1)) * (xscurr(1, 0) - xscurr(0, 0));
        elevector3[0] = areaele;
      }
    }
    break;
    case calc_struct_centerdisp:
    {
      if (Comm.MyPID() == Owner())
      {
        // element geometry update
        Teuchos::RCP<const Epetra_Vector> disptotal = discretization.GetState("displacementtotal");
        if (disptotal == Teuchos::null) dserror("Cannot get state vector 'displacementtotal'");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disptotal, mydisp, lm);
        const int numnod = NumNode();
        Epetra_SerialDenseMatrix xsrefe(Wall1::numdim_, numnod);  // material coord. of element
        Epetra_SerialDenseMatrix xscurr(Wall1::numdim_, numnod);  // current coord. of element
        for (int i = 0; i < numnod; ++i)
        {
          xsrefe(0, i) = Nodes()[i]->X()[0];
          xsrefe(1, i) = Nodes()[i]->X()[1];

          xscurr(0, i) = xsrefe(0, i) + mydisp[i * noddof];
          xscurr(1, i) = xsrefe(1, i) + mydisp[i * noddof + 1];
        }

        // integration of the displacements over the surface
        const int dim = Wall1::numdim_;
        const DiscretizationType distype = Shape();

        // gaussian points
        const DRT::UTILS::GaussRule1D gaussrule = getOptimalGaussrule(distype);
        const DRT::UTILS::IntegrationPoints1D intpoints(gaussrule);  //

        // allocate vector for shape functions and for derivatives
        LINALG::SerialDenseVector funct(numnod);
        LINALG::SerialDenseMatrix deriv(1, numnod);

        Teuchos::RCP<const Epetra_Vector> dispincr = discretization.GetState("displacementincr");
        std::vector<double> edispincr(lm.size());
        DRT::UTILS::ExtractMyValues(*dispincr, edispincr, lm);

        elevector2[0] = 0;

        for (int gpid = 0; gpid < intpoints.nquad; gpid++)
        {
          const double e1 = intpoints.qxg[gpid][0];  // coordinate of GP

          // get values of shape functions and derivatives in the line at specific GP
          DRT::UTILS::shape_function_1D(funct, e1, distype);
          DRT::UTILS::shape_function_1D_deriv1(deriv, e1, distype);

          double dr = w1_substitution(xscurr, deriv, NULL, numnod);

          elevector2[0] += intpoints.qwgt[gpid] * dr;

          for (int d = 0; d < dim; d++)
          {
            if (gpid == 0) elevector3[d] = 0;

            for (int j = 0; j < numnod; ++j)
            {
              elevector3[d] += funct[j] * intpoints.qwgt[gpid] * edispincr[j * dim + d] * dr;
            }
          }
        }
      }
    }
    break;

    case calc_struct_areaconstrstiff:
    {
      if (distype != line2)
      {
        dserror("Area Constraint only works for line2 curves!");
      }  // element geometry update
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null)
      {
        dserror("Cannot get state vector 'displacement'");
      }
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      const int numnod = NumNode();
      Epetra_SerialDenseMatrix xsrefe(numnod, Wall1::numdim_);  // material coord. of element
      Epetra_SerialDenseMatrix xscurr(numnod, Wall1::numdim_);  // material coord. of element
      for (int i = 0; i < numnod; ++i)
      {
        xsrefe(i, 0) = Nodes()[i]->X()[0];
        xsrefe(i, 1) = Nodes()[i]->X()[1];

        xscurr(i, 0) = xsrefe(i, 0) + mydisp[i * noddof];
        xscurr(i, 1) = xsrefe(i, 1) + mydisp[i * noddof + 1];
      }
      // call submethods
      ComputeAreaConstrStiff(xscurr, elematrix1);
      ComputeAreaConstrDeriv(xscurr, elevector1);
      elevector2 = elevector1;
      // compute area between line and x-Axis
      double areaele = 0.5 * (xscurr(0, 1) + xscurr(1, 1)) * (xscurr(1, 0) - xscurr(0, 0));
      elevector3[0] = areaele;
    }
    break;
    default:
      dserror("Unimplemented type of action for Soh8Surface");
      break;
  }
  return 0;
}

/*----------------------------------------------------------------------*
 * Evaluate method on mutliple dofsets                       vuong 11/12*
 * ---------------------------------------------------------------------*/
int DRT::ELEMENTS::Wall1Line::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, LocationArray& la, Epetra_SerialDenseMatrix& elematrix1,
    Epetra_SerialDenseMatrix& elematrix2, Epetra_SerialDenseVector& elevector1,
    Epetra_SerialDenseVector& elevector2, Epetra_SerialDenseVector& elevector3)
{
  if (la.Size() == 1)
  {
    return Evaluate(params, discretization,
        la[0].lm_,  // location vector is build by the first column of la
        elematrix1, elematrix2, elevector1, elevector2, elevector3);
  }

  const DiscretizationType distype = Shape();

  // start with "none"
  DRT::ELEMENTS::Wall1Line::ActionType act = Wall1Line::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "calc_struct_area_poro")
    act = Wall1Line::calc_struct_area_poro;
  else
    dserror("Unknown type of action for StructuralSurface");

  // what the element has to do
  switch (act)
  {
    case calc_struct_area_poro:
    {
      // get the parent element
      const DRT::Element* parentele = ParentElement();
      const int nenparent = parentele->NumNode();
      // get element location vector and ownerships
      std::vector<int> lmpar;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      parentele->LocationVector(discretization, lmpar, lmowner, lmstride);

      // gaussian points
      const DRT::UTILS::GaussRule1D gaussrule = getOptimalGaussrule(distype);
      // get integration rule
      const DRT::UTILS::IntPointsAndWeights<1> intpoints(gaussrule);

      const int ngp = intpoints.IP().nquad;
      Teuchos::RCP<Epetra_SerialDenseVector> poro = Teuchos::rcp(new Epetra_SerialDenseVector(ngp));
      const int numdim = 2;
      const int numnode = NumNode();
      const int noddof = NumDofPerNode(*(Nodes()[0]));

      // element geometry update
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vector 'displacement'");
      std::vector<double> mydisp(lmpar.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lmpar);

      // update element geometry
      Epetra_SerialDenseMatrix xrefe(numdim, nenparent);  // material coord. of element
      Epetra_SerialDenseMatrix xcurr(numdim, nenparent);  // current  coord. of element

      const DRT::Node* const* nodes = parentele->Nodes();
      for (int i = 0; i < nenparent; ++i)
      {
        const double* x = nodes[i]->X();
        xrefe(0, i) = x[0];
        xrefe(1, i) = x[1];

        xcurr(0, i) = xrefe(0, i) + mydisp[i * noddof + 0];
        xcurr(1, i) = xrefe(1, i) + mydisp[i * noddof + 1];
      }

      // number of degrees of freedom per node of fluid
      const int numdofpernode = 3;

      Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState(1, "fluidvel");
      if (velnp == Teuchos::null) dserror("Cannot get state vector 'fluidvel'");
      // extract local values of the global vectors
      std::vector<double> myvelpres(la[1].lm_.size());
      DRT::UTILS::ExtractMyValues(*velnp, myvelpres, la[1].lm_);

      Epetra_SerialDenseVector mypres(numnode);
      for (int inode = 0; inode < numnode; ++inode)  // number of nodes
      {
        (mypres)(inode, 0) = myvelpres[numdim + (inode * numdofpernode)];
      }

      LINALG::SerialDenseMatrix pqxg;
      Epetra_SerialDenseMatrix derivtrafo;

      DRT::UTILS::BoundaryGPToParentGP<2>(
          pqxg, derivtrafo, intpoints, parentele->Shape(), distype, FaceParentNumber());

      for (int gp = 0; gp < ngp; ++gp)
      {
        // get shape functions and derivatives in the plane of the element
        LINALG::SerialDenseVector funct(nenparent);
        LINALG::SerialDenseMatrix deriv(2, nenparent);
        DRT::UTILS::shape_function_2D(funct, pqxg(gp, 0), pqxg(gp, 1), parentele->Shape());
        DRT::UTILS::shape_function_2D_deriv1(deriv, pqxg(gp, 0), pqxg(gp, 1), parentele->Shape());

        LINALG::SerialDenseVector funct1D(numnode);
        DRT::UTILS::shape_function_1D(funct1D, intpoints.IP().qxg[gp][0], Shape());

        // pressure at integration point
        double press = funct1D.Dot(mypres);

        // get Jacobian matrix and determinant w.r.t. spatial configuration
        //! transposed jacobian "dx/ds"
        LINALG::SerialDenseMatrix xjm(numdim, numdim);
        xjm.Multiply('N', 'T', 1.0, deriv, xcurr, 0.0);
        LINALG::SerialDenseMatrix Jmat(numdim, numdim);
        Jmat.Multiply('N', 'T', 1.0, deriv, xrefe, 0.0);

        double det = 0.0;
        double detJ = 0.0;

        if (numdim == 2)
        {
          det = xjm(0, 0) * xjm(1, 1) - xjm(0, 1) * xjm(1, 0);
          detJ = Jmat(0, 0) * Jmat(1, 1) - Jmat(0, 1) * Jmat(1, 0);
          ;
        }
        else
          dserror("not implemented");

        const double J = det / detJ;

        // get structure material
        Teuchos::RCP<MAT::StructPoro> structmat =
            Teuchos::rcp_static_cast<MAT::StructPoro>(parentele->Material());
        if (structmat == Teuchos::null) dserror("invalid structure material for poroelasticity");
        double porosity = 0.0;
        structmat->ComputeSurfPorosity(params, press, J, FaceParentNumber(), gp, porosity);
      }
    }
    break;
    default:
      dserror("Unimplemented type of action for Soh8Surface");
      break;
  }
  return 0;
}

/*----------------------------------------------------------------------*
 * Compute first derivatives of area                            tk 10/07*
 * with respect to the displacements                                    *
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1Line::ComputeAreaConstrDeriv(
    Epetra_SerialDenseMatrix xscurr, Epetra_SerialDenseVector& elevector)
{
  if (elevector.Length() != 4)
  {
    std::cout << "Length of element Vector: " << elevector.Length() << std::endl;
    dserror("That is not the right size!");
  }
  // implementation of simple analytic solution
  elevector[0] = -xscurr(0, 1) - xscurr(1, 1);
  elevector[1] = xscurr(1, 0) - xscurr(0, 0);
  elevector[2] = xscurr(0, 1) + xscurr(1, 1);
  elevector[3] = xscurr(1, 0) - xscurr(0, 0);
  elevector.Scale(-0.5);
  return;
}

/*----------------------------------------------------------------------*
 * Compute influence of area constraint on stiffness matrix.    tk 10/07*
 * Second derivatives of areas with respect to the displacements        *
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1Line::ComputeAreaConstrStiff(
    Epetra_SerialDenseMatrix xscurr, Epetra_SerialDenseMatrix& elematrix)
{
  elematrix(0, 0) = 0.0;
  elematrix(0, 1) = -0.5;
  elematrix(0, 2) = 0.0;
  elematrix(0, 3) = -0.5;

  elematrix(1, 0) = -0.5;
  elematrix(1, 1) = 0.0;
  elematrix(1, 2) = 0.5;
  elematrix(1, 3) = 0.0;

  elematrix(2, 0) = 0.0;
  elematrix(2, 1) = 0.5;
  elematrix(2, 2) = 0.0;
  elematrix(2, 3) = 0.5;

  elematrix(3, 0) = -0.5;
  elematrix(3, 1) = 0.0;
  elematrix(3, 2) = 0.5;
  elematrix(3, 3) = 0.0;

  elematrix.Scale(-1.0);
  return;
}
