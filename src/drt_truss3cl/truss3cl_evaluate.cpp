/*----------------------------------------------------------------------------*/
/*! \file
\brief three dimensional interpolated total Lagrange hybrid beam-truss element
 (can be connected to beam3 elements)

\level 3

\maintainer Maximilian Grill

*/
/*---------------------------------------------------------------------------*/

#include "truss3cl.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_lib/standardtypes_cpp.H"

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public) mukherjee 01/14|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Truss3CL::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  const int nnode = 4;
  DRT::ELEMENTS::Truss3CL::ActionType act = Truss3CL::calc_none;
  // get the action required
  std::string action = params.get<std::string>("action", "calc_none");
  if (action == "calc_none")
    dserror("No action supplied");
  else if (action == "calc_struct_linstiff")
    act = Truss3CL::calc_struct_linstiff;
  else if (action == "calc_struct_nlnstiff")
    act = Truss3CL::calc_struct_nlnstiff;
  else if (action == "calc_struct_internalforce")
    act = Truss3CL::calc_struct_internalforce;
  else if (action == "calc_struct_linstiffmass")
    act = Truss3CL::calc_struct_linstiffmass;
  else if (action == "calc_struct_nlnstiffmass")
    act = Truss3CL::calc_struct_nlnstiffmass;
  else if (action == "calc_struct_nlnstifflmass")
    act = Truss3CL::calc_struct_nlnstifflmass;
  else if (action == "calc_struct_stress")
    act = Truss3CL::calc_struct_stress;
  else if (action == "calc_struct_eleload")
    act = Truss3CL::calc_struct_eleload;
  else if (action == "calc_struct_fsiload")
    act = Truss3CL::calc_struct_fsiload;
  else if (action == "calc_struct_update_istep")
    act = Truss3CL::calc_struct_update_istep;
  else if (action == "calc_struct_reset_istep")
    act = Truss3CL::calc_struct_reset_istep;
  else if (action == "postprocess_stress")
    act = Truss3CL::postprocess_stress;
  else if (action == "calc_struct_ptcstiff")
    act = Truss3CL::calc_struct_ptcstiff;
  else if (action == "calc_struct_energy")
    act = Truss3CL::calc_struct_energy;
  else
  {
    std::cout << action << std::endl;
    dserror("Unknown type of action for Truss3CL");
  }

  switch (act)
  {
    case Truss3CL::calc_struct_ptcstiff:
    {
      dserror("PTC scheme not fully functional for Truss3CL elements for biopolymer applications!");
      EvaluatePTC<2, 3, 3>(params, elemat1);
    }
    break;
    /*in case that only linear stiffness matrix is required b3_nlstiffmass is called with zero
     displacement and residual values*/
    case Truss3CL::calc_struct_linstiff:
    {
      // only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");
    }
    break;
    // calculate internal energy
    case Truss3CL::calc_struct_energy:
    {
      // need current global displacement and get them from discretization
      // making use of the local-to-global map lm one can extract current displacemnet and residual
      // values for each degree of freedom
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      const int nnode = NumNode();
      switch (nnode)
      {
        case 4:
        {
          t3_energy(params, mydisp, &elevec1);
          break;
        }
        default:
          dserror("Only Line2 Elements implemented.");
      }
    }
    break;

    // nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is
    // required
    case Truss3CL::calc_struct_nlnstiffmass:
    case Truss3CL::calc_struct_nlnstifflmass:
    case Truss3CL::calc_struct_nlnstiff:
    case Truss3CL::calc_struct_internalforce:
    {
      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacement and residual
      // values for each degree of freedom
      //
      // get element displacements
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      // get residual displacements
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (res == Teuchos::null) dserror("Cannot get state vectors 'residual displacement'");
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);

      /*first displacement vector is modified for proper element evaluation in case of periodic
       *boundary conditions; in case that no periodic boundary conditions are to be applied the
       *following code line may be ignored or deleted*/
      NodeShift<nnode, 3>(params, mydisp);

      // only if random numbers for Brownian dynamics are passed to element, get element velocities
      std::vector<double> myvel(lm.size());
      if (params.get<Teuchos::RCP<Epetra_MultiVector>>("RandomNumbers", Teuchos::null) !=
          Teuchos::null)
      {
        Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState("velocity");
        DRT::UTILS::ExtractMyValues(*vel, myvel, lm);
      }

      // for engineering strains instead of total lagrange use t3_nlnstiffmass2
      if (act == Truss3CL::calc_struct_nlnstiffmass)
        t3_nlnstiffmass(params, myvel, mydisp, &elemat1, &elemat2, &elevec1);
      else if (act == Truss3CL::calc_struct_nlnstifflmass)
        t3_nlnstiffmass(params, myvel, mydisp, &elemat1, &elemat2, &elevec1);
      else if (act == Truss3CL::calc_struct_nlnstiff)
      {
        for (int i = 0; i < 24; i++)
          t3_nlnstiffmass(params, myvel, mydisp, &elemat1, NULL, &elevec1);
      }

      else if (act == Truss3CL::calc_struct_internalforce)
        t3_nlnstiffmass(params, myvel, mydisp, NULL, NULL, &elevec1);
    }
    break;
    case calc_struct_update_istep:
    {
      // nothing to do
    }
    break;
    case calc_struct_reset_istep:
    {
      // nothing to do
    }
    break;
    case calc_struct_stress:
    {
      // no stress calculation implemented! Do not crash simulation and just keep quiet!
    }
    break;
    case postprocess_stress:
    {
      // no stress calculation for postprocess. Does not really make sense!
      dserror("No stress output for Truss3CL!");
    }
    break;
    default:
      dserror("Unknown type of action for Truss3CL %d", act);
      break;
  }
  return 0;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public) mukherjee 01/14|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Truss3CL::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  dserror("This method is yet to be configured for bio-polymer networks!");
  // get element displacements
  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp == Teuchos::null) dserror("Cannot get state vector 'displacement'");
  std::vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);


  // find out whether we will use a time curve
  const double time = params.get("total time", -1.0);

  // find out whether we will use a time curve and get the factor
  const std::vector<int>* tmp_funct = condition.Get<std::vector<int>>("funct");
  int functnum = -1;
  // number of the load curve related with a specific line Neumann condition called
  if (tmp_funct) functnum = (*tmp_funct)[0];
  // amplitude of load curve at current time called
  double functfac = 1.0;
  if (functnum >= 0)  // notation for this function similar to Crisfield, Volume 1;
    functfac = DRT::Problem::Instance()->Funct(functnum).EvaluateTime(time);

  // jacobian determinant
  double det = lrefe_ / 2;

  const DiscretizationType distype = this->Shape();

  // gaussian points
  const DRT::UTILS::IntegrationPoints1D intpoints(gaussrule_);

  // declaration of variable in order to store shape function
  Epetra_SerialDenseVector funct(NumNode());

  // get values and switches from the condition

  // onoff is related to the first 6 flags of a line Neumann condition in the input file;
  // value 1 for flag i says that condition is active for i-th degree of freedom
  const std::vector<int>* onoff = condition.Get<std::vector<int>>("onoff");
  // val is related to the 6 "val" fields after the onoff flags of the Neumann condition
  // in the input file; val gives the values of the force as a multiple of the prescribed load curve
  const std::vector<double>* val = condition.Get<std::vector<double>>("val");

  // integration loops
  for (int ip = 0; ip < intpoints.nquad; ++ip)
  {
    // integration points in parameter space and weights
    const double xi = intpoints.qxg[ip][0];
    const double wgt = intpoints.qwgt[ip];

    // evaluation of shape funcitons at Gauss points
    DRT::UTILS::shape_function_1D(funct, xi, distype);

    double fac = 0;
    fac = wgt * det;

    /*load vector ar; regardless of the actual number of degrees of freedom active with respect to
     *this element or certain nodes of it the vector val has always the lengths 6 and in order to
     *deal with possibly different numbers of actually used DOF we always loop through all the 6*/
    double ar[6];
    // loop the dofs of a node

    for (int i = 0; i < 6; ++i) ar[i] = fac * (*onoff)[i] * (*val)[i] * functfac;

    for (int dof = 0; dof < 3; ++dof)
    {
      // computing entries for first node
      elevec1[dof] += funct[0] * ar[dof];
      // computing entries for first node
      elevec1[3 + dof] += funct[1] * ar[dof];
    }

  }  // for (int ip=0; ip<intpoints.nquad; ++ip)

  return 0;
}


/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate PTC damping (public) mukherjee 01/14|
 *----------------------------------------------------------------------------------------------------------*/
template <int fnnode, int ndim, int dof>  // number of nodes, number of dimensions of embedding
                                          // space, number of degrees of freedom per node
int DRT::ELEMENTS::Truss3CL::EvaluatePTC(
    Teuchos::ParameterList& params, Epetra_SerialDenseMatrix& elemat1)
{
  for (int i = 0; i < elemat1.RowDim(); i++)
    for (int j = 0; j < elemat1.ColDim(); j++)
      if (i == j) elemat1(i, j) += params.get<double>("crotptc", 0.0) * 0.5 * lrefe_ / 2;


  return 0;
}  // DRT::ELEMENTS::Truss3CL::EvaluatePTC

/*--------------------------------------------------------------------------------------*
 | calculation of elastic energy                                         mukherjee 01/14|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3CL::t3_energy(
    Teuchos::ParameterList& params, std::vector<double>& disp, Epetra_SerialDenseVector* intenergy)
{
  dserror("This method is yet to be configured for bio-polymer networks!");

  // get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double ym = 0;

  // assignment of material parameters; only St.Venant material is accepted for this truss
  switch (currmat->MaterialType())
  {
    case INPAR::MAT::m_stvenant:  // only linear elastic material supported
    {
      const MAT::StVenantKirchhoff* actmat =
          static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      ym = actmat->Youngs();
    }
    break;
    default:
      dserror("unknown or improper type of material law");
      break;
  }

  // current node position (first entries 0 .. 2 for first node, 3 ..5 for second node)
  LINALG::Matrix<6, 1> xcurr;

  // auxiliary vector for both internal force and stiffness matrix: N^T_(,xi)*N_(,xi)*xcurr
  LINALG::Matrix<6, 1> aux_lcurr;

  // calculate interpolated displacements and velocities of the two fictitious nodes

  // vector containing fictitious nodal displacements
  std::vector<double> fdisp(6);
  for (int i = 0; i < 6; i++) fdisp[i] = 0;
  // Interpolation of Translational displacements at fictitious nodes. Use hermite polynmials
  // Here we need the position of the internodal Binding Spot at time t=0, when the filament beam
  // elements where initialized
  LINALG::Matrix<12, 1> XRefOnlyDisp;
  XRefOnlyDisp.Clear();

  // Obtain cordinates of real nodes
  for (int node = 0; node < 4; node++)
    for (int d = 0; d < 3; d++) XRefOnlyDisp(node * 3 + d) = Nodes()[node]->X()[d];


  // length of first element and second element required for hermite shape function interpolation
  std::vector<double> LengthofElementatRef(2);
  // auxilary vector
  std::vector<double> aux(6);
  aux[0] = XRefOnlyDisp(0) - XRefOnlyDisp(3);
  aux[1] = XRefOnlyDisp(1) - XRefOnlyDisp(4);
  aux[2] = XRefOnlyDisp(2) - XRefOnlyDisp(5);
  aux[3] = XRefOnlyDisp(6) - XRefOnlyDisp(9);
  aux[4] = XRefOnlyDisp(7) - XRefOnlyDisp(10);
  aux[5] = XRefOnlyDisp(8) - XRefOnlyDisp(11);


  for (int filament = 0; filament < 2; filament++)
  {
    LengthofElementatRef[filament] = sqrt(
        pow(aux[3 * filament], 2) + pow(aux[3 * filament + 1], 2) + pow(aux[3 * filament + 2], 2));
  }

  // Evaluate Shape Functions at Binding positions
  std::vector<LINALG::Matrix<1, 4>> Ibp(2);
  const DiscretizationType distype = this->Shape();
  for (int filament = 0; filament < 2; filament++)
  {
    DRT::UTILS::shape_function_hermite_1D(
        Ibp[filament], mybindingposition_[filament], LengthofElementatRef[filament], distype);
  }

  // calculate interpolated fictitious nodal displacaments.
  for (int j = 0; j < 2; j++)  // Number of binding spots
    for (int l = 0; l < 2;
         l++)  // l=0 left real node, l=0 right real node contributing to the binding spot
      for (int r = 0; r < 2; r++)    //  right node contributing to the binding spot
        for (int i = 0; i < 3; i++)  // Number of dimensions
          fdisp[i + 3 * j] += Ibp[j](2 * l + r) * disp[i + 3 * r + 6 * l + 12 * j];

  // strain
  double epsilon;

  // current nodal position (first
  for (int j = 0; j < 3; ++j)
  {
    xcurr(j) = Nodes()[0]->X()[j] + fdisp[j];          // first node
    xcurr(j + 3) = Nodes()[1]->X()[j] + fdisp[3 + j];  // second node
  }

  // computing auxiliary vector aux = 4.0*N^T_{,xi} * N_{,xi} * xcurr
  aux_lcurr(0) = (xcurr(0) - xcurr(3));
  aux_lcurr(1) = (xcurr(1) - xcurr(4));
  aux_lcurr(2) = (xcurr(2) - xcurr(5));
  aux_lcurr(3) = (xcurr(3) - xcurr(0));
  aux_lcurr(4) = (xcurr(4) - xcurr(1));
  aux_lcurr(5) = (xcurr(5) - xcurr(2));

  double lcurr = sqrt(pow(aux_lcurr(0), 2) + pow(aux_lcurr(1), 2) + pow(aux_lcurr(2), 2));


  switch (kintype_)
  {
    case tr3_totlag:
    {
      // calculating Green-Lagrange strain epsilon
      epsilon = 0.5 * (pow(lcurr / lrefe_, 2) - 1.0);

      // W_int = 1/2*E*A*lrefe*\epsilon^2
      (*intenergy)(0) = 0.5 * (ym * crosssec_ * lrefe_ * pow(epsilon, 2));
    }
    break;
    case tr3_engstrain:
    {
      // calculating strain epsilon from node position by scalar product:
      epsilon = (lcurr - lrefe_) / lrefe_;

      // W_int = 1/2*E*A*lrefe*\epsilon^2
      (*intenergy)(0) = 0.5 * (ym * crosssec_ * lrefe_ * pow(epsilon, 2));
    }
    break;
    default:
      dserror("Unknown type kintype_ for Truss3CL");
      break;
  }

  return;
}

/*--------------------------------------------------------------------------------------*
 | switch between kintypes                                               mukherjee 01/14|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3CL::t3_nlnstiffmass(Teuchos::ParameterList& params,
    std::vector<double>& vel, std::vector<double>& disp, Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix, Epetra_SerialDenseVector* force)
{
  const int fnnode = 2;

  // calculate interpolated displacements and velocities of the two fictitious nodes
  /* vector containing fictitious nodal displacements.
   * It is to be noted that the fictitious truss element contains only 6 dofs (3 translational
   * degrees of freedom in each node). However, for ease in calculation and future extension to full
   * 12 dofs for beam element the size of the vector is kept as 12. Therefore, displacement degrees
   * of freedoms are stored in 0,1,2 (1st node) and 6,7,8 (2nd node).
   */
  LINALG::Matrix<6, 1> fdisp;
  fdisp.Clear();

  // vector containing fictitious nodal velocities
  /* It is to be noted that the fictitious truss element contains only 6 dofs (3 translational
   * degrees of freedom in each node). Therefore, displacement degrees of freedoms are stored in
   * 0...2 (1st node) and 3...5 (2nd node).
   */
  LINALG::Matrix<6, 1> fvel;
  fvel.Clear();
  // Evaluate Shape Functions at Binding positions
  /*
   * For interpolating binding sites and displacement quantities, Hermite polynomials need to be
   * implemented.
   */

  // Here we need the position of the internodal Binding Spot at time t=0, when the filament beam
  // elements where initialized

  LINALG::Matrix<12, 1> XRefOnlyDisp;
  XRefOnlyDisp.Clear();
  for (int node = 0; node < 4; node++)  // element has four nodes
    for (int k = 0; k < 3; k++) XRefOnlyDisp(node * 3 + k) = Nodes()[node]->X()[k];

  // length of first element and second element required for hermite shape function interpolation
  std::vector<double> LengthofElementatRef(2);
  // auxiliary vector for norm/length calculations
  std::vector<double> aux(6);

  aux[0] = XRefOnlyDisp(0) - XRefOnlyDisp(3);
  aux[1] = XRefOnlyDisp(1) - XRefOnlyDisp(4);
  aux[2] = XRefOnlyDisp(2) - XRefOnlyDisp(5);
  aux[3] = XRefOnlyDisp(6) - XRefOnlyDisp(9);
  aux[4] = XRefOnlyDisp(7) - XRefOnlyDisp(10);
  aux[5] = XRefOnlyDisp(8) - XRefOnlyDisp(11);


  for (int filament = 0; filament < 2; filament++)
    LengthofElementatRef[filament] = sqrt(
        pow(aux[3 * filament], 2) + pow(aux[3 * filament + 1], 2) + pow(aux[3 * filament + 2], 2));


  // Compute Hermite polynomials in current/spatial position, we need the length of filaments at
  // reference configuration
  std::vector<LINALG::Matrix<1, 4>> Ibp(2);
  const DiscretizationType distype = this->Shape();
  for (int filament = 0; filament < 2; filament++)
    DRT::UTILS::shape_function_hermite_1D(
        Ibp[filament], mybindingposition_[filament], LengthofElementatRef[filament], distype);

  // calculate interpolated fictitious nodal displacaments.
  for (int j = 0; j < 2; j++)  // Number of binding spots
    for (int l = 0; l < 2;
         l++)  // l=0 left real node, l=0 right real node contributing to the binding spot
      for (int r = 0; r < 2; r++)    //  right node contributing to the binding spot
        for (int i = 0; i < 3; i++)  // Number of dimensions
          fdisp(i + 3 * j) += Ibp[j](2 * l + r) * disp[i + 3 * r + 6 * l + 12 * j];

  // calculate interpolated fictitious nodal velocities
  for (int j = 0; j < 2; j++)  // Number of binding spots
    for (int l = 0; l < 2;
         l++)  // l=0 left real node, l=0 right real node contributing to the binding spot
      for (int r = 0; r < 2; r++)    //  right node contributing to the binding spot
        for (int i = 0; i < 3; i++)  // Number of dimensions
          fvel(i + 3 * j) += Ibp[j](2 * l + r) * vel[i + 3 * r + 6 * l + 12 * j];

  // 12x12 Stiffness Matrix of the beam described by the two fictitious nodes
  Epetra_SerialDenseMatrix fstiffmatrix;
  fstiffmatrix.Shape(3 * fnnode, 3 * fnnode);
  fstiffmatrix.Scale(0);
  // 12x1 force vector of the beam described by the two fictitious nodes
  Epetra_SerialDenseVector fforce;
  fforce.Size(3 * fnnode);
  fforce.Scale(0);

  switch (kintype_)
  {
    case tr3_totlag:
    {
      t3_nlnstiffmass_totlag(fdisp, fstiffmatrix, massmatrix, fforce);
    }
    break;
    case tr3_engstrain:
    {
      t3_nlnstiffmass_engstr(fdisp, fstiffmatrix, massmatrix, fforce);
    }
    break;
    default:
      dserror("Unknown type kintype_ for Truss3CL");
      break;
  }

  // internal force vector stored by class variable before application of stochastic excitations
  if (params.get<std::string>("internalforces", "no") == "yes")
    f_ = Teuchos::rcp(new Epetra_SerialDenseVector(fforce));

  /*the following function call applies statistical forces and damping matrix according to the
   * fluctuation dissipation theorem; it is dedicated to the application of Truss3CL elements in the
   * frame of statistical mechanics problems; for these problems a special vector has to be passed
   * to the element packed in the params parameter list; in case that the control routine calling
   * the element does not attach this special vector to params the following method is just doing
   * nothing, which means that for any ordinary problem of structural mechanics it may be ignored*/
  CalcBrownian<2, 3, 3, 3>(params, fvel, disp, fdisp, fstiffmatrix, fforce);

  // TODO: Check why all 24x24 entries in the stiffness matrix and all 24x1 entries in the force
  // vector are non-zero Shouldn't only the original 12x12 matrix entries and 12x1 vector entries be
  // there? On that note, also check whether the NumDofPerNode function of the truss3cl element
  // should return 3 or 6...
  if (stiffmatrix != NULL)
  {
    // Transform 12x12 stiffness matrix for fictitious nodal beam back to 24x24 stiffness matrix of
    // 4-noded beam
    //  Only for degrees of freedom corresponding to nodal displacement
    for (int nodei = 0; nodei < fnnode; nodei++)
      for (int nodej = 0; nodej < fnnode; nodej++)
        for (int Ri = 0; Ri < 2; Ri++)
          for (int Rl = 0; Rl < 2; Rl++)
            for (int Rj = 0; Rj < 2; Rj++)
              for (int Rr = 0; Rr < 2; Rr++)
                for (int i = 0; i < 3; i++)
                  for (int j = 0; j < 3; j++)
                    (*stiffmatrix)(
                        12 * nodei + 6 * Ri + 3 * Rl + i, 12 * nodej + 6 * Rj + 3 * Rr + j) +=
                        Ibp[nodei](2 * Ri + Rl) * Ibp[nodej](2 * Rj + Rr) *
                        fstiffmatrix(3 * nodei + i, 3 * nodej + j);
  }
  if (force != NULL)
  {
    // Tranform 12x1 force vector for fictitious nodal beam to 24x1 force vector of the 4-noded beam
    for (int node = 0; node < fnnode; ++node)
      for (int Ri = 0; Ri < 2; Ri++)
        for (int Rl = 0; Rl < 2; Rl++)
          for (int i = 0; i < 3; i++)
            (*force)(12 * node + 6 * Ri + 3 * Rl + i) +=
                Ibp[node](2 * Ri + Rl) * fforce(3 * node + i);
  }

  return;
}


/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private) mukherjee 01/14|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3CL::t3_nlnstiffmass_totlag(LINALG::Matrix<6, 1>& fdisp,
    Epetra_SerialDenseMatrix& fstiffmatrix, Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector& fforce)
{
  // current node position (first entries 0 .. 2 for first node, 3 ..5 for second node)
  LINALG::Matrix<6, 1> xcurr;

  /*current nodal displacement (first entries 0 .. 2 for first node, 3 ..5 for second node) compared
   * to reference configuration; note: in general this is not equal to the values in fdisp since the
   * fdisp referes to a nodal displacement compared to a reference configuration before the first
   * time step whereas ucurr referes to the displacement with respect to a reference
   * configuration which may have been set up at any point of time during the simulation (usually
   * this is only important if an element has been added to the discretization after the start of
   * the simulation)*/
  LINALG::Matrix<6, 1> ucurr;

  // Green-Lagrange strain
  double epsilon;

  // auxiliary vector for both internal force and stiffness matrix: N^T_(,xi)*N_(,xi)*xcurr
  LINALG::Matrix<6, 1> aux;

  // Vector for storing reference degrees of freedom at t=0 at binding spot or the connecting nodes
  // of the truss element
  LINALG::Matrix<6, 1> Xstart;
  Xstart.Clear();
  ComputeRefDOFs(Xstart);

  // current nodal position
  for (int j = 0; j < 3; ++j)
  {
    xcurr(j) = Xstart(j) + fdisp(j);              // first node
    xcurr(j + 3) = Xstart(3 + j) + fdisp(3 + j);  // second node
  }

  // current displacement = current position - reference position
  ucurr = xcurr;
  ucurr -= X_;

  // computing auxiliary vector aux = N^T_{,xi} * N_{,xi} * xcurr
  aux(0) = (xcurr(0) - xcurr(3));
  aux(1) = (xcurr(1) - xcurr(4));
  aux(2) = (xcurr(2) - xcurr(5));
  aux(3) = (xcurr(3) - xcurr(0));
  aux(4) = (xcurr(4) - xcurr(1));
  aux(5) = (xcurr(5) - xcurr(2));

  // calculating strain epsilon from node position by scalar product:
  // epsilon = (xrefe + 0.5*ucurr)^T * N_{,s}^T * N_{,s} * d
  epsilon = 0;
  epsilon += (X_(0) + 0.5 * ucurr(0)) * (ucurr(0) - ucurr(3));
  epsilon += (X_(1) + 0.5 * ucurr(1)) * (ucurr(1) - ucurr(4));
  epsilon += (X_(2) + 0.5 * ucurr(2)) * (ucurr(2) - ucurr(5));
  epsilon += (X_(3) + 0.5 * ucurr(3)) * (ucurr(3) - ucurr(0));
  epsilon += (X_(4) + 0.5 * ucurr(4)) * (ucurr(4) - ucurr(1));
  epsilon += (X_(5) + 0.5 * ucurr(5)) * (ucurr(5) - ucurr(2));
  epsilon /= lrefe_ * lrefe_;

  // get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double ym = 0;


  // assignment of material parameters; only St.Venant material is accepted for this truss
  switch (currmat->MaterialType())
  {
    case INPAR::MAT::m_stvenant:  // only linear elastic material supported
    {
      const MAT::StVenantKirchhoff* actmat =
          static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      ym = actmat->Youngs();
    }
    break;
    default:
      dserror("unknown or improper type of material law");
      break;
  }

  // computing global internal forces
  for (int i = 0; i < 6; ++i)
  {
    fforce(i) = (ym * crosssec_ * epsilon / lrefe_) * aux(i);
  }

  // computing linear stiffness matrix
  for (int i = 0; i < 3; ++i)
  {
    // stiffness entries for first node
    fstiffmatrix(i, i) = (ym * crosssec_ * epsilon / lrefe_);
    fstiffmatrix(i, 3 + i) = -(ym * crosssec_ * epsilon / lrefe_);
    // stiffness entries for second node
    fstiffmatrix(i + 3, i + 3) = (ym * crosssec_ * epsilon / lrefe_);
    fstiffmatrix(i + 3, i) = -(ym * crosssec_ * epsilon / lrefe_);
  }

  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      fstiffmatrix(i, j) += (ym * crosssec_ / pow(lrefe_, 3)) * aux(i) * aux(j);

  // calculating consistent mass matrix
  if (massmatrix != NULL)
  {
    for (int i = 0; i < 6; ++i)
    {
      (*massmatrix)(i, i) = 1;
    }
  }
  return;
}  // DRT::ELEMENTS::Truss3CL::t3_nlnstiffmass_totlag


/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private) mukherjee 01/14| | engineering strain measure,
 large displacements and rotations                                                |
  *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3CL::t3_nlnstiffmass_engstr(const LINALG::Matrix<6, 1>& fdisp,
    Epetra_SerialDenseMatrix& stiffmatrix, Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector& force)
{
  // current node position (first entries 0 .. 2 for first node, 3 ..5 for second node)
  LINALG::Matrix<6, 1> xcurr;

  // Green-Lagrange strain
  double epsilon;

  // auxiliary vector for both internal force and stiffness matrix: N^T_(,xi)*N_(,xi)*xcurr
  LINALG::Matrix<6, 1> aux;

  // Vector for storing reference degrees of freedom at t=0 at binding spot or the connecting nodes
  // of the truss element
  LINALG::Matrix<6, 1> Xstart;
  Xstart.Clear();
  ComputeRefDOFs(Xstart);

  // current nodal position (first
  for (int j = 0; j < 3; ++j)
  {
    xcurr(j) = Xstart(j) + fdisp(j);              // first node
    xcurr(j + 3) = Xstart(3 + j) + fdisp(3 + j);  // second node
  }

  // computing auxiliary vector aux = 4.0*N^T_{,xi} * N_{,xi} * xcurr
  aux(0) = (xcurr(0) - xcurr(3));
  aux(1) = (xcurr(1) - xcurr(4));
  aux(2) = (xcurr(2) - xcurr(5));
  aux(3) = (xcurr(3) - xcurr(0));
  aux(4) = (xcurr(4) - xcurr(1));
  aux(5) = (xcurr(5) - xcurr(2));

  double lcurr = sqrt(pow(aux(0), 2) + pow(aux(1), 2) + pow(aux(2), 2));

  // calculating strain epsilon from node position by scalar product:
  epsilon = (lcurr - lrefe_) / lrefe_;

  // get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double ym = 0;
  //  double density = 0;

  // assignment of material parameters; only St.Venant material is accepted for this truss
  switch (currmat->MaterialType())
  {
    case INPAR::MAT::m_stvenant:  // only linear elastic material supported
    {
      const MAT::StVenantKirchhoff* actmat =
          static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      ym = actmat->Youngs();
    }
    break;
    default:
      dserror("unknown or improper type of material law");
      break;
  }
  // resulting force scaled by current length
  double forcescalar = (ym * crosssec_ * epsilon) / lcurr;

  // computing global internal forces
  for (int i = 0; i < 6; ++i) force(i) = forcescalar * aux(i);


  // computing linear stiffness matrix
  for (int i = 0; i < 3; ++i)
  {
    // stiffness entries for first node
    stiffmatrix(i, i) = forcescalar;
    stiffmatrix(i, 3 + i) = -forcescalar;
    // stiffness entries for second node
    stiffmatrix(i + 3, i + 3) = forcescalar;
    stiffmatrix(i + 3, i) = -forcescalar;
  }

  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      stiffmatrix(i, j) += (ym * crosssec_ / pow(lcurr, 3)) * aux(i) * aux(j);

  // calculating consistent mass matrix
  if (massmatrix != NULL)
  {
    for (int i = 0; i < 6; ++i)
    {
      (*massmatrix)(i, i) = 1;
    }
  }

  return;
}  // DRT::ELEMENTS::Truss3CL::bt_nlnstiffmass3_engstr


// lump mass matrix
void DRT::ELEMENTS::Truss3CL::t3_lumpmass(Epetra_SerialDenseMatrix* emass)
{
  // lump mass matrix
  if (emass != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (int c = 0; c < (*emass).N(); ++c)  // parse columns
    {
      double d = 0.0;
      for (int r = 0; r < (*emass).M(); ++r)  // parse rows
      {
        d += (*emass)(r, c);  // accumulate row entries
        (*emass)(r, c) = 0.0;
      }
      (*emass)(c, c) = d;  // apply sum of row entries on diagonal
    }
  }
}

/*-----------------------------------------------------------------------------------------------------------*
 | computes damping coefficients per lengthand stores them in a matrix in the following order:
 damping of    | | translation parallel to filament axis, damping of translation orthogonal to
 filament axis, damping of     | | rotation around filament axis (public)       mukherjee   01/14|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Truss3CL::MyDampingConstants(
    Teuchos::ParameterList& params, LINALG::Matrix<3, 1>& gamma)
{
  // translational damping coefficients according to Howard, p. 107, table 6.2;
  gamma(0) = 2 * PI * params.get<double>("ETA", 0.0);
  gamma(1) = 4 * PI * params.get<double>("ETA", 0.0);
  // no rotational damping as no rotaional degrees of freedom
  gamma(2) = 0;

}  // DRT::ELEMENTS::Truss3CL::MyDampingConstants

/*-----------------------------------------------------------------------------------------------------------*
 |computes the number of different random numbers required in each time step for generation of
 stochastic    | |forces; (public)         mukherjee 01/14|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Truss3CL::HowManyRandomNumbersINeed()
{
  /*at each Gauss point one needs as many random numbers as randomly excited degrees of freedom,
   *i.e. three random numbers for the translational degrees of freedom*/
  return (3 * 4);

}  // DRT::ELEMENTS::Truss3CL::HowManyRandomNumbersINeed

/*-----------------------------------------------------------------------------------------------------------*
 |computes velocity of background fluid and gradient of that velocity at a certain evaluation point
 in       | |the physical space                                                         (public)
 mukherjee 01/14|
 *----------------------------------------------------------------------------------------------------------*/
template <int ndim>  // number of dimensions of embedding space
void DRT::ELEMENTS::Truss3CL::MyBackgroundVelocity(
    Teuchos::ParameterList& params,                  //!< parameter list
    const LINALG::Matrix<ndim, 1>& evaluationpoint,  //!< point at which background velocity and its
                                                     //!< gradient has to be computed
    LINALG::Matrix<ndim, 1>& velbackground,          //!< velocity of background fluid
    LINALG::Matrix<ndim, ndim>& velbackgroundgrad)   //!< gradient of velocity of background fluid
{
  /*note: this function is not yet a general one, but always assumes a shear flow, where the
   * velocity of the background fluid is always directed in x-direction. In 3D the velocity
   * increases linearly in z and equals zero for z = 0. In 2D the velocity increases linearly in y
   * and equals zero for y = 0. */

  velbackground.PutScalar(0);
  velbackground(0) = evaluationpoint(ndim - 1) * params.get<double>("CURRENTSHEAR", 0.0);

  if (velbackground(0) != 0)
    dserror("Method MyTranslationalDamping needs to be modified before continuing!");

  velbackgroundgrad.PutScalar(0);
  velbackgroundgrad(0, ndim - 1) = params.get<double>("CURRENTSHEAR", 0.0);
}

/*-----------------------------------------------------------------------------------------------------------*
 | computes translational damping forces and stiffness (public) mukherjee 01/14|
 *----------------------------------------------------------------------------------------------------------*/
template <int fnnode, int ndim,
    int dof>  // number of nodes, number of dimensions of embedding space, number of degrees of
              // freedom per node
inline void
DRT::ELEMENTS::Truss3CL::MyTranslationalDamping(Teuchos::ParameterList& params,  //!< parameter list
    const LINALG::Matrix<6, 1>& fvel,        //!< element velocity vector
    const std::vector<double>& disp,         //!< element real disp vector
    const LINALG::Matrix<6, 1>& fdisp,       //!< element fictitious disp vector
    Epetra_SerialDenseMatrix& fstiffmatrix,  //!< element stiffness matrix
    Epetra_SerialDenseVector& fforce)        //!< element internal force vector
{
  // get time step size
  double dt = params.get<double>("delta time", 0.0);

  // velocity and gradient of background velocity field
  LINALG::Matrix<3, 1> velbackground;
  LINALG::Matrix<3, 3> velbackgroundgrad;

  // evaluation point in physical space corresponding to a certain Gauss point in parameter space
  LINALG::Matrix<ndim, 1> evaluationpoint;

  // damping coefficients for translational and rotatinal degrees of freedom
  LINALG::Matrix<3, 1> gamma(true);
  MyDampingConstants(params, gamma);

  // get vector jacobi with Jacobi determinants at each integration point (gets by default those
  // values required for consistent damping matrix)
  std::vector<double> jacobi(jacobimass_);

  // determine type of numerical integration performed (lumped damping matrix via lobatto
  // integration!)
  IntegrationType integrationtype = gaussexactintegration;

  // get Gauss points and weights for evaluation of damping matrix
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(fnnode, integrationtype));

  // matrix to store basis functions and their derivatives evaluated at a certain Gauss point
  LINALG::Matrix<1, fnnode> funct;
  LINALG::Matrix<1, fnnode> deriv;

  // Vector for storing reference degrees of freedom at t=0 at binding spot or the connecting nodes
  // of the truss element
  LINALG::Matrix<6, 1> Xstart;
  Xstart.Clear();
  ComputeRefDOFs(Xstart);

  // for calculations inside truss element, Lagrange polynomials are used
  for (int gp = 0; gp < gausspoints.nquad; gp++)
  {
    // evaluate basis functions and their derivatives at current Gauss point
    DRT::UTILS::shape_function_1D(funct, gausspoints.qxg[gp][0], Shape());
    DRT::UTILS::shape_function_1D_deriv1(deriv, gausspoints.qxg[gp][0], Shape());

    // compute point in phyiscal space corresponding to Gauss point
    evaluationpoint.PutScalar(0);
    // loop over all line nodes
    for (int i = 0; i < fnnode; i++)
      // loop over all dimensions
      for (int j = 0; j < ndim; j++)
        evaluationpoint(j) += funct(i) * (Xstart(3 * i + j) + fdisp(dof * i + j));

    // compute velocity and gradient of background flow field at evaluationpoint
    MyBackgroundVelocity<ndim>(params, evaluationpoint, velbackground, velbackgroundgrad);
    //    if (velbackground.Norm2()!= 0 && velbackgroundgrad.Norm2()!=0)
    //      dserror("Non-zero background velocity found! Please modify the algorithm accordingly.");

    // compute tangent vector t_{\par} at current Gauss point
    LINALG::Matrix<ndim, 1> tpar(true);
    for (int i = 0; i < fnnode; i++)
      for (int k = 0; k < ndim; k++)
        tpar(k) += deriv(i) * (Xstart(3 * i + k) + fdisp(dof * i + k)) / jacobi[gp];

    // compute velocity vector at this Gauss point
    LINALG::Matrix<ndim, 1> velgp(true);
    for (int i = 0; i < fnnode; i++)
      for (int l = 0; l < ndim; l++) velgp(l) += funct(i) * fvel(dof * i + l);

    // compute matrix product (t_{\par} \otimes t_{\par}) \cdot velbackgroundgrad
    LINALG::Matrix<ndim, ndim> tpartparvelbackgroundgrad(true);
    for (int i = 0; i < ndim; i++)
      for (int j = 0; j < ndim; j++)
        for (int k = 0; k < ndim; k++)
          tpartparvelbackgroundgrad(i, j) += tpar(i) * tpar(k) * velbackgroundgrad(k, j);


    // loop over all line nodes
    for (int i = 0; i < fnnode; i++)
      // loop over lines of matrix t_{\par} \otimes t_{\par}
      for (int k = 0; k < ndim; k++)
        // loop over columns of matrix t_{\par} \otimes t_{\par}
        for (int l = 0; l < ndim; l++)
        {
          fforce(i * dof + k) += funct(i) * jacobi[gp] * gausspoints.qwgt[gp] *
                                 ((k == l) * gamma(1) + (gamma(0) - gamma(1)) * tpar(k) * tpar(l)) *
                                 (velgp(l) - velbackground(l));

          // loop over all column nodes
          for (int j = 0; j < fnnode; j++)
          {
            fstiffmatrix(i * dof + k, j * dof + l) +=
                gausspoints.qwgt[gp] * funct(i) * funct(j) * jacobi[gp] *
                ((k == l) * gamma(1) + (gamma(0) - gamma(1)) * tpar(k) * tpar(l)) / dt;
            fstiffmatrix(i * dof + k, j * dof + l) -=
                gausspoints.qwgt[gp] * funct(i) * funct(j) * jacobi[gp] *
                (velbackgroundgrad(k, l) * gamma(1) +
                    (gamma(0) - gamma(1)) * tpartparvelbackgroundgrad(k, l));

            fstiffmatrix(i * dof + k, j * dof + k) += gausspoints.qwgt[gp] * funct(i) * deriv(j) *
                                                      (gamma(0) - gamma(1)) * tpar(l) *
                                                      (velgp(l) - velbackground(l));
            fstiffmatrix(i * dof + k, j * dof + l) += gausspoints.qwgt[gp] * funct(i) * deriv(j) *
                                                      (gamma(0) - gamma(1)) * tpar(k) *
                                                      (velgp(l) - velbackground(l));
          }
        }
  }

  return;
}  // DRT::ELEMENTS::Truss3CL::MyTranslationalDamping(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes stochastic forces and resulting stiffness (public) mukherjee 01/14|
 *----------------------------------------------------------------------------------------------------------*/
template <int fnnode, int ndim, int dof,
    int randompergauss>  // number of nodes, number of dimensions of embedding space, number of
                         // degrees of freedom per node, number of random numbers required per Gauss
                         // point
inline void
DRT::ELEMENTS::Truss3CL::MyStochasticForces(Teuchos::ParameterList& params,  //!< parameter list
    const LINALG::Matrix<6, 1>& fdisp,       //!< element disp vector
    Epetra_SerialDenseMatrix& fstiffmatrix,  //!< element stiffness matrix
    Epetra_SerialDenseVector& fforce)        //!< element internal force vector
{
  // damping coefficients for three translational and one rotatinal degree of freedom
  LINALG::Matrix<3, 1> gamma(true);
  MyDampingConstants(params, gamma);

  // get vector jacobi with Jacobi determinants at each integration point (gets by default those
  // values required for consistent damping matrix)
  std::vector<double> jacobi(jacobimass_);

  // determine type of numerical integration performed (lumped damping matrix via lobatto
  // integration!)
  IntegrationType integrationtype = gaussexactintegration;

  // get Gauss points and weights for evaluation of damping matrix
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(fnnode, integrationtype));

  // matrix to store basis functions and their derivatives evaluated at a certain Gauss point
  LINALG::Matrix<1, fnnode> funct;
  LINALG::Matrix<1, fnnode> deriv;


  /*get pointer at Epetra multivector in parameter list linking to random numbers for stochastic
   * forces with zero mean and standard deviation (2*kT / dt)^0.5; note carefully: a space between
   * the two subsequal ">" signs is mandatory for the C++ parser in order to avoid confusion with
   * ">>" for streams*/
  Teuchos::RCP<Epetra_MultiVector> randomnumbers =
      params.get<Teuchos::RCP<Epetra_MultiVector>>("RandomNumbers", Teuchos::null);

  // Vector for storing reference degrees of freedom at t=0 at binding spot or the connecting nodes
  // of the truss element
  LINALG::Matrix<6, 1> Xstart;
  Xstart.Clear();
  ComputeRefDOFs(Xstart);

  // For calculations inside truss element, Lagrange polynomials are used.
  for (int gp = 0; gp < gausspoints.nquad; gp++)
  {
    // evaluate basis functions and their derivatives at current Gauss point
    DRT::UTILS::shape_function_1D(funct, gausspoints.qxg[gp][0], Shape());
    DRT::UTILS::shape_function_1D_deriv1(deriv, gausspoints.qxg[gp][0], Shape());

    // compute tangent vector t_{\par} at current Gauss point
    LINALG::Matrix<ndim, 1> tpar(true);
    for (int i = 0; i < fnnode; i++)
      for (int k = 0; k < ndim; k++)
        tpar(k) += deriv(i) * (Xstart(3 * i + k) + fdisp(dof * i + k)) / jacobi[gp];


    // loop over all line nodes
    for (int i = 0; i < fnnode; i++)
      // loop dimensions with respect to lines
      for (int k = 0; k < ndim; k++)
        // loop dimensions with respect to columns
        for (int l = 0; l < ndim; l++)
        {
          fforce(i * dof + k) -=
              funct(i) *
              (sqrt(gamma(1)) * (k == l) + (sqrt(gamma(0)) - sqrt(gamma(1))) * tpar(k) * tpar(l)) *
              (*randomnumbers)[gp * randompergauss + l][LID()] *
              sqrt(jacobi[gp] * gausspoints.qwgt[gp]);

          // loop over all column nodes
          for (int j = 0; j < fnnode; j++)
          {
            fstiffmatrix(i * dof + k, j * dof + k) -=
                funct(i) * deriv(j) * tpar(l) * (*randomnumbers)[gp * randompergauss + l][LID()] *
                sqrt(gausspoints.qwgt[gp] / jacobi[gp]) * (sqrt(gamma(0)) - sqrt(gamma(1)));
            fstiffmatrix(i * dof + k, j * dof + l) -=
                funct(i) * deriv(j) * tpar(k) * (*randomnumbers)[gp * randompergauss + l][LID()] *
                sqrt(gausspoints.qwgt[gp] / jacobi[gp]) * (sqrt(gamma(0)) - sqrt(gamma(1)));
          }
        }
  }
  return;
}  // DRT::ELEMENTS::Truss3CL::MyStochasticForces(.)


/*-----------------------------------------------------------------------------------------------------------*
 | Assemble stochastic and viscous forces and respective stiffness according to fluctuation
 dissipation      | | theorem (public) mukherjee 01/14|
 *----------------------------------------------------------------------------------------------------------*/
template <int fnnode, int ndim, int dof,
    int randompergauss>  // number of nodes, number of dimensions of embedding space, number of
                         // degrees of freedom per node, number of random numbers required per Gauss
                         // point
inline void
DRT::ELEMENTS::Truss3CL::CalcBrownian(Teuchos::ParameterList& params,
    const LINALG::Matrix<6, 1>& fvel,        //!< element velocity vector
    const std::vector<double>& disp,         //!< element real displacement vector
    const LINALG::Matrix<6, 1>& fdisp,       //!< element fictitious displacement vector
    Epetra_SerialDenseMatrix& fstiffmatrix,  //!< element stiffness matrix
    Epetra_SerialDenseVector& fforce)        //!< element internal force vector
{
  // if no random numbers for generation of stochastic forces are passed to the element no Brownian
  // dynamics calculations are conducted
  if (params.get<Teuchos::RCP<Epetra_MultiVector>>("RandomNumbers", Teuchos::null) == Teuchos::null)
    return;

  // add stiffness and forces due to translational damping effects
  MyTranslationalDamping<fnnode, ndim, dof>(params, fvel, disp, fdisp, fstiffmatrix, fforce);

  // add stochastic forces and (if required) resulting stiffness
  MyStochasticForces<fnnode, ndim, dof, randompergauss>(params, fdisp, fstiffmatrix, fforce);

  return;

}  // DRT::ELEMENTS::Truss3CL::CalcBrownian(.)

/*-----------------------------------------------------------------------------------------------------------*
 | shifts nodes so that proper evaluation is possible even in case of periodic boundary conditions;
 if two   | | nodes within one element are separated by a periodic boundary, one of them is shifted
 such that the final | | distance in R^3 is the same as the initial distance in the periodic space;
 the shift affects computation  | | on element level within that very iteration step, only (no
 change in global variables performed)          |                                 | | (public)
 mukherjee 01/14|
 *----------------------------------------------------------------------------------------------------------*/
template <int nnode, int ndim>  // number of nodes, number of dimensions
inline void DRT::ELEMENTS::Truss3CL::NodeShift(Teuchos::ParameterList& params,  //!< parameter list
    std::vector<double>& disp)  //!< element disp vector
{
  dserror("This does not work any longer with new statmech");

  //  /*get number of degrees of freedom per node; note: the following function assumes the same
  //  number of degrees
  //   *of freedom for each element node*/
  //  int numdof = NumDofPerNode(*(Nodes()[0]));
  //  if (nnode==4 && disp.size()==24)
  //      numdof = 6;
  //  double time = params.get<double>("total time",0.0);
  //  double starttime = params.get<double>("STARTTIMEACT",0.0);
  //  double dt = params.get<double>("delta time");
  //  double shearamplitude = params.get<double> ("SHEARAMPLITUDE", 0.0);
  //  int curvenumber = params.get<int> ("CURVENUMBER", -1)-1;
  //  int dbcdispdir = params.get<int> ("DBCDISPDIR", -1)-1;
  //  Teuchos::RCP<std::vector<double> > defvalues = Teuchos::rcp(new std::vector<double>(3,0.0));
  //  Teuchos::RCP<std::vector<double> > periodlength = params.get("PERIODLENGTH", defvalues);
  //  INPAR::STATMECH::DBCType dbctype = params.get<INPAR::STATMECH::DBCType>("DBCTYPE",
  //  INPAR::STATMECH::dbctype_std); bool shearflow = false;
  //  if(dbctype==INPAR::STATMECH::dbctype_shearfixed ||
  //  dbctype==INPAR::STATMECH::dbctype_sheartrans || dbctype==INPAR::STATMECH::dbctype_affineshear)
  //    shearflow = true;
  //
  //  /*only if periodic boundary conditions are in use, i.e. params.get<double>("PeriodLength",0.0)
  //  > 0.0, this
  //   * method has to change the displacement variables*/
  //  if(periodlength->at(0) > 0.0)
  //    //loop through all nodes except for the first node which remains fixed as reference node
  //    for(int i=1;i<nnode;i++)
  //    {
  //      for(int dof= ndim - 1; dof > -1; dof--)
  //      {
  //        /*if the distance in some coordinate direction between some node and the first node
  //        becomes smaller by adding or subtracting
  //         * the period length, the respective node has obviously been shifted due to periodic
  //         boundary conditions and should be shifted
  //         * back for evaluation of element matrices and vectors; this way of detecting shifted
  //         nodes works as long as the element length
  //         * is smaller than half the periodic length*/
  //        if( fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) + periodlength->at(dof) -
  //        (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < fabs(
  //        (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) )
  //        )
  //        {
  //          disp[numdof*i+dof] += periodlength->at(dof);
  //
  //          /*the upper domain surface orthogonal to the z-direction may be subject to shear
  //          Dirichlet boundary condition; the lower surface
  //           *may be fixed by DBC. To avoid problmes when nodes exit the domain through the upper
  //           z-surface and reenter through the lower *z-surface, the shear has to be substracted
  //           from nodal coordinates in that case */
  //          if(shearflow && dof == 2 && curvenumber >= 0 && time>starttime &&
  //          fabs(time-starttime)>dt/1e4)
  //            disp[numdof*i+dbcdispdir] +=
  //            shearamplitude*DRT::Problem::Instance()->Funct(curvenumber).EvaluateTimeDerivative(time);
  //        }
  //
  //        if( fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - periodlength->at(dof) -
  //        (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < fabs(
  //        (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) )
  //        )
  //        {
  //          disp[numdof*i+dof] -= periodlength->at(dof);
  //
  //          /*the upper domain surface orthogonal to the z-direction may be subject to shear
  //          Dirichlet boundary condition; the lower surface
  //           *may be fixed by DBC. To avoid problmes when nodes exit the domain through the lower
  //           z-surface and reenter through the upper *z-surface, the shear has to be added to
  //           nodal coordinates in that case */
  //          if(shearflow && dof == 2 && curvenumber >= 0 && time>starttime &&
  //          fabs(time-starttime)>dt/1e4)
  //            disp[numdof*i+dbcdispdir] -=
  //            shearamplitude*DRT::Problem::Instance()->Funct(curvenumber).EvaluateTimeDerivative(time);
  //        }
  //      }
  //    }


  return;

}  // DRT::ELEMENTS::Truss3CL::NodeShift

/*-----------------------------------------------------------------------------------------------------------*
 | computes reference DOFs for TRUSS3CL at t=0, i.e. before first time step (public) mukherjee
 01/14|
 *----------------------------------------------------------------------------------------------------------*/

/* It is to be stressed out that the reference configuration of the filaments are not same as the
 * reference configuration of the cross-linkers. The filaments exist from the beginning of the
 * time-step. On the other hand, the cross-linkers are added in the middle of the solution. However,
 * we want to evaluate the forces
 * acting on the cross-linkers from the beginning of the simulation. This particular method serves
 * this purpose.*/

inline void DRT::ELEMENTS::Truss3CL::ComputeRefDOFs(LINALG::Matrix<6, 1>& Xstart)
{
  //  Here we need the position of the internodal Binding Spot at time t=0,
  //  when the filament beam elements where initialized
  LINALG::Matrix<1, 12> XRefOnlyDisp;
  XRefOnlyDisp.Clear();

  for (int node = 0; node < 4; node++)  // element has four nodes
    for (int k = 0; k < 3; k++) XRefOnlyDisp(node * 3 + k) = Nodes()[node]->X()[k];

  /* Compute Lagrange polynomials in reference position.
     Since beams are straight at initial configuration, it's a valid assumption. */

  std::vector<LINALG::Matrix<1, 2>> Ibp(2);
  const DiscretizationType distype = this->Shape();

  for (int filament = 0; filament < 2; filament++)
  {
    DRT::UTILS::shape_function_1D(Ibp[filament], mybindingposition_[filament], distype);
  }

  for (int filament = 0; filament < 2; filament++)
    for (int node = 0; node < 2; node++)
      for (int dof = 0; dof < 3; dof++)
      {
        Xstart(dof + 3 * filament) +=
            Ibp[filament](node) * XRefOnlyDisp(dof + 3 * node + 6 * filament);
      }

  return;
}  // DRT::ELEMENTS::Truss3CL::ComputeRefDOFs
