/*!----------------------------------------------------------------------
  \file turbulence_statistics_ccy.cpp

\brief Compute (time and space) averaged values for turbulent flows
       around a rotating cylinder and write them to files.

\level 2
<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include "turbulence_statistics_ccy.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_control_point.H"
#include "../drt_mat/newtonianfluid.H"

/*----------------------------------------------------------------------

                  Standard Constructor (public)

  ---------------------------------------------------------------------*/
FLD::TurbulenceStatisticsCcy::TurbulenceStatisticsCcy(Teuchos::RCP<DRT::Discretization> actdis,
    bool alefluid, Teuchos::RCP<Epetra_Vector> dispnp, Teuchos::ParameterList& params,
    const std::string& statistics_outfilename, const bool withscatra)
    : discret_(actdis),
      dispnp_(dispnp),
      params_(params),
      statistics_outfilename_(statistics_outfilename),
      withscatra_(withscatra),
      numscatradofpernode_(0)
{
  //----------------------------------------------------------------------
  // plausibility check
  int numdim = params_.get<int>("number of velocity degrees of freedom");
  if (numdim != 3)
  {
    dserror("Evaluation of turbulence statistics only for 3d flows!");
  }

  //----------------------------------------------------------------------
  // allocate some vectors
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  meanvelnp_ = LINALG::CreateVector(*dofrowmap, true);

  if (withscatra_)
  {
    meanscanp_ = LINALG::CreateVector(*dofrowmap, true);
    // meanfullphinp_ is initalized in ApplyScatraResults()
  }

  //----------------------------------------------------------------------
  // switches, control parameters, material parameters

  // get the plane normal direction from the parameterlist
  {
    std::string plainstring =
        params_.sublist("TURBULENCE MODEL").get<std::string>("HOMDIR", "not_specified");

    if (plainstring == "z")
    {
      dim_ = 2;
    }
    else
    {
      dserror("homogeneous direction for this flow was specified incorrectly. (need z)");
    }
  }

  // ---------------------------------------------------------------------
  // up to now, there are no records written
  countrecord_ = 0;

  // ---------------------------------------------------------------------
  // compute all planes for sampling

  // available shells of element corners (Nurbs) of elements
  nodeshells_ = Teuchos::rcp(new std::vector<double>);

  // available homogeneous (sampling) shells --- there are
  // numsubdivisions layers per element layer
  shellcoordinates_ = Teuchos::rcp(new std::vector<double>);

  const int numsubdivisions = 5;

  // try to cast discretisation to nurbs variant
  // this tells you what kind of computation of
  // samples is required
  DRT::NURBS::NurbsDiscretization* nurbsdis =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*actdis));

  if (nurbsdis == NULL)
  {
    dserror("Need Nurbs mesh for turbulent flows around a circular cylinder\n");
  }
  else
  {
    // real pointwise control point sampling does not make any sense
    // for Nurbs discretisations since shape functions are not interpolating

    // radial shellcoordinates are determined by the element
    // (cartesian) number in the second knotspan direction and
    // the number of sampling shells in between are added

    // for nurbs discretisations, all vector sizes are already determined
    // by the knotvector size

    // get nurbs dis' knotvector sizes
    std::vector<int> n_x_m_x_l(nurbsdis->Return_n_x_m_x_l(0));

    // get nurbs dis' element numbers
    std::vector<int> nele_x_mele_x_lele(nurbsdis->Return_nele_x_mele_x_lele(0));

    // get the knotvector itself
    Teuchos::RCP<DRT::NURBS::Knotvector> knots = nurbsdis->GetKnotVector();

    // resize and initialise to 0
    {
      (*nodeshells_).resize(nele_x_mele_x_lele[1] + 1);
      (*shellcoordinates_).resize(nele_x_mele_x_lele[1] * (numsubdivisions - 1) + 1);

      std::vector<double>::iterator coord;

      for (coord = (*nodeshells_).begin(); coord != (*nodeshells_).end(); ++coord)
      {
        *coord = 0;
      }
      for (coord = shellcoordinates_->begin(); coord != shellcoordinates_->end(); ++coord)
      {
        *coord = 0;
      }
    }

    // count numbers of nodes in homogeneous shells
    int nodeshellsize = nodeshells_->size();

    std::vector<int> nodeshells_numnodes(nodeshellsize, 0);
    std::vector<int> shellcoordinates_numnodes((*shellcoordinates_).size(), 0);

    // get element map
    const Epetra_Map* elementmap = nurbsdis->ElementRowMap();

    // loop all available elements
    for (int iele = 0; iele < elementmap->NumMyElements(); ++iele)
    {
      // get element pointer
      DRT::Element* const actele = nurbsdis->gElement(elementmap->GID(iele));

      // want to loop all control points of the element,
      // so get the number of points
      const int numnp = actele->NumNode();

      // get the elements control points/nodes
      DRT::Node** nodes = actele->Nodes();

      // acquire weights from nodes
      Epetra_SerialDenseVector weights(numnp);

      for (int inode = 0; inode < numnp; ++inode)
      {
        DRT::NURBS::ControlPoint* cp = dynamic_cast<DRT::NURBS::ControlPoint*>(nodes[inode]);

        weights(inode) = cp->W();
      }

      // get gid, location in the patch
      int gid = actele->Id();

      int patchid = 0;

      std::vector<int> ele_cart_id(3);
      knots->ConvertEleGidToKnotIds(gid, patchid, ele_cart_id);

      // access elements knot span
      std::vector<Epetra_SerialDenseVector> knots(3);
      bool zero_size = (*((*nurbsdis).GetKnotVector())).GetEleKnots(knots, actele->Id());

      // zero sized elements have to be skipped
      if (zero_size)
      {
        continue;
      }

      // get shapefunctions, compute all visualisation point positions
      Epetra_SerialDenseVector nurbs_shape_funct(numnp);

      switch (actele->Shape())
      {
        case DRT::Element::nurbs8:
        case DRT::Element::nurbs27:
        {
          // element local point position
          Epetra_SerialDenseVector uv(3);

          {
            // standard

            //               v
            //              /
            //  w  7       /   8
            //  ^   +---------+
            //  |  /         /|
            //  | /         / |
            // 5|/        6/  |
            //  +---------+   |
            //  |         |   |
            //  |         |   +
            //  |         |  / 4
            //  |         | /
            //  |         |/
            //  +---------+ ----->u
            // 1           2
            // use r-coordinate of point 2 and 4
            // temporary x vector
            std::vector<double> x(3);

            // point 1
            uv(0) = 1.0;
            uv(1) = -1.0;
            uv(2) = -1.0;
            DRT::NURBS::UTILS::nurbs_get_3D_funct(
                nurbs_shape_funct, uv, knots, weights, actele->Shape());
            for (int isd = 0; isd < 3; ++isd)
            {
              double val = 0;
              for (int inode = 0; inode < numnp; ++inode)
              {
                val += (((nodes[inode])->X())[isd]) * nurbs_shape_funct(inode);
              }
              x[isd] = val;
            }

            (*nodeshells_)[ele_cart_id[1]] += sqrt(x[0] * x[0] + x[1] * x[1]);
            (*shellcoordinates_)[ele_cart_id[1] * (numsubdivisions - 1)] +=
                sqrt(x[0] * x[0] + x[1] * x[1]);

            nodeshells_numnodes[ele_cart_id[1]] += 1;
            shellcoordinates_numnodes[ele_cart_id[1] * (numsubdivisions - 1)] += 1;


            for (int rr = 1; rr < numsubdivisions - 1; ++rr)
            {
              uv(1) += 2.0 / (numsubdivisions - 1);

              DRT::NURBS::UTILS::nurbs_get_3D_funct(
                  nurbs_shape_funct, uv, knots, weights, actele->Shape());
              for (int isd = 0; isd < 3; ++isd)
              {
                double val = 0;
                for (int inode = 0; inode < numnp; ++inode)
                {
                  val += (((nodes[inode])->X())[isd]) * nurbs_shape_funct(inode);
                }
                x[isd] = val;
              }
              (*shellcoordinates_)[ele_cart_id[1] * (numsubdivisions - 1) + rr] +=
                  sqrt(x[0] * x[0] + x[1] * x[1]);
              ++(shellcoordinates_numnodes[ele_cart_id[1] * (numsubdivisions - 1) + rr]);
            }


            // set upper point of element, too (only for last layer)
            if (ele_cart_id[1] + 1 == nele_x_mele_x_lele[1])
            {
              // point 8
              uv(0) = 1.0;
              uv(1) = 1.0;
              uv(2) = -1.0;
              DRT::NURBS::UTILS::nurbs_get_3D_funct(
                  nurbs_shape_funct, uv, knots, weights, actele->Shape());
              for (int isd = 0; isd < 3; ++isd)
              {
                double val = 0;
                for (int inode = 0; inode < numnp; ++inode)
                {
                  val += (((nodes[inode])->X())[isd]) * nurbs_shape_funct(inode);
                }
                x[isd] = val;
              }

              (*nodeshells_)[ele_cart_id[1] + 1] += sqrt(x[0] * x[0] + x[1] * x[1]);
              (*shellcoordinates_)[(ele_cart_id[1] + 1) * (numsubdivisions - 1)] +=
                  sqrt(x[0] * x[0] + x[1] * x[1]);
              ++(nodeshells_numnodes[ele_cart_id[1] + 1]);
              ++(shellcoordinates_numnodes[(ele_cart_id[1] + 1) * (numsubdivisions - 1)]);
            }
          }
          break;
        }
        default:
          dserror(
              "Unknown element shape for a nurbs element or nurbs type not valid for turbulence "
              "calculation\n");
      }
    }

    //----------------------------------------------------------------------
    // add contributions from all processors, normalize

    std::vector<double> lnodeplanes(*nodeshells_);
    std::vector<double> lplanecoordinates(*shellcoordinates_);

    std::vector<int> lnodeshells_numnodes(nodeshells_numnodes);
    std::vector<int> lshellcoordinates_numnodes(shellcoordinates_numnodes);

    discret_->Comm().SumAll(&(lnodeplanes[0]), &((*nodeshells_)[0]), nodeshells_->size());
    discret_->Comm().SumAll(
        &(lplanecoordinates[0]), &((*shellcoordinates_)[0]), shellcoordinates_->size());

    discret_->Comm().SumAll(
        &(lnodeshells_numnodes[0]), &(nodeshells_numnodes[0]), nodeshells_numnodes.size());
    discret_->Comm().SumAll(&(lshellcoordinates_numnodes[0]), &(shellcoordinates_numnodes[0]),
        shellcoordinates_numnodes.size());

    {
      (*nodeshells_).resize(nele_x_mele_x_lele[1] + 1);
      (*shellcoordinates_).resize(nele_x_mele_x_lele[1] * (numsubdivisions - 1) + 1);

      for (unsigned rr = 0; rr < (*nodeshells_).size(); ++rr)
      {
        if (fabs(nodeshells_numnodes[rr]) < 1e-9)
        {
          dserror("zero nodes in shell layer %d\n", rr);
        }

        (*nodeshells_)[rr] /= (double)(nodeshells_numnodes[rr]);
      }


      for (unsigned rr = 0; rr < (*shellcoordinates_).size(); ++rr)
      {
        if (fabs(shellcoordinates_numnodes[rr]) < 1e-9)
        {
          dserror("zero nodes in sampling shell layer %d\n", rr);
        }

        (*shellcoordinates_)[rr] /= (double)(shellcoordinates_numnodes[rr]);
      }
    }
  }

  //----------------------------------------------------------------------
  // sort shellcoordinates and nodeshells
  {
    std::set<double, PlaneSortCriterion> shellset;

    std::vector<double>::iterator coord;
    std::set<double, PlaneSortCriterion>::iterator shell;

    {
      for (coord = (*nodeshells_).begin(); coord != (*nodeshells_).end(); ++coord)
      {
        shellset.insert(*coord);
      }

      int rr = 0;
      for (shell = shellset.begin(); shell != shellset.end(); ++shell)
      {
        (*nodeshells_)[rr] = *shell;
        ++rr;
      }
    }

    shellset.clear();

    {
      for (coord = (*shellcoordinates_).begin(); coord != (*shellcoordinates_).end(); ++coord)
      {
        shellset.insert(*coord);
      }

      int rr = 0;
      for (shell = shellset.begin(); shell != shellset.end(); ++shell)
      {
        (*shellcoordinates_)[rr] = *shell;
        ++rr;
      }
    }
  }

  //----------------------------------------------------------------------
  // allocate arrays for sums of in plane mean values

  int size = shellcoordinates_->size();

  // arrays for point based averaging
  // --------------------------------

  // first order moments
  pointsumu_ = Teuchos::rcp(new std::vector<double>);
  pointsumu_->resize(size, 0.0);

  pointsumv_ = Teuchos::rcp(new std::vector<double>);
  pointsumv_->resize(size, 0.0);

  pointsumw_ = Teuchos::rcp(new std::vector<double>);
  pointsumw_->resize(size, 0.0);

  pointsump_ = Teuchos::rcp(new std::vector<double>);
  pointsump_->resize(size, 0.0);

  // now the second order moments
  pointsumuu_ = Teuchos::rcp(new std::vector<double>);
  pointsumuu_->resize(size, 0.0);

  pointsumvv_ = Teuchos::rcp(new std::vector<double>);
  pointsumvv_->resize(size, 0.0);

  pointsumww_ = Teuchos::rcp(new std::vector<double>);
  pointsumww_->resize(size, 0.0);

  pointsumpp_ = Teuchos::rcp(new std::vector<double>);
  pointsumpp_->resize(size, 0.0);

  pointsumuv_ = Teuchos::rcp(new std::vector<double>);
  pointsumuv_->resize(size, 0.0);

  pointsumuw_ = Teuchos::rcp(new std::vector<double>);
  pointsumuw_->resize(size, 0.0);

  pointsumvw_ = Teuchos::rcp(new std::vector<double>);
  pointsumvw_->resize(size, 0.0);

  if (withscatra_)
  {
    pointsumc_ = Teuchos::rcp(new std::vector<double>);
    pointsumc_->resize(size, 0.0);

    pointsumcc_ = Teuchos::rcp(new std::vector<double>);
    pointsumcc_->resize(size, 0.0);

    // pointsumphi_/pointsumphiphi_ are allocated in ApplyScatraResults()
  }

  //----------------------------------------------------------------------
  // initialise output

  Teuchos::RCP<std::ofstream> log;

  if (discret_->Comm().MyPID() == 0)
  {
    std::string s(statistics_outfilename_);
    s.append(".flow_statistics");

    log = Teuchos::rcp(new std::ofstream(s.c_str(), std::ios::out));
    (*log) << "# Statistics for turbulent incompressible flow in a rotating cylinder (first- and "
              "second-order moments)\n\n";

    log->flush();
  }

  // clear statistics
  this->ClearStatistics();

  return;
}  // TurbulenceStatisticsCcy::TurbulenceStatisticsCcy

/*----------------------------------------------------------------------*

                           Destructor

 -----------------------------------------------------------------------*/
FLD::TurbulenceStatisticsCcy::~TurbulenceStatisticsCcy()
{
  return;
}  // TurbulenceStatisticsCcy::~TurbulenceStatisticsCcy()

/*----------------------------------------------------------------------*

       Compute the in-plane mean values of first and second order
       moments for velocities, pressure and Cs are added to global
                            'sum' vectors.

 -----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCcy::DoTimeSample(Teuchos::RCP<Epetra_Vector> velnp,
    Teuchos::RCP<Epetra_Vector> scanp, Teuchos::RCP<Epetra_Vector> fullphinp)
{
  // we have an additional sample
  numsamp_++;

  // meanvelnp is a refcount copy of velnp
  meanvelnp_->Update(1.0, *velnp, 0.0);

  if (withscatra_)
  {
    if (scanp != Teuchos::null)
      meanscanp_->Update(1.0, *scanp, 0.0);
    else
      dserror("Vector scanp is Teuchos::null");

    if (fullphinp != Teuchos::null)
    {
      int err = meanfullphinp_->Update(1.0, *fullphinp, 0.0);
      if (err) dserror("Could not update meanfullphinp_");
    }
    else
      dserror("Vector fullphinp is Teuchos::null");
  }

  //----------------------------------------------------------------------
  // loop planes and calculate pointwise means in each plane
  this->EvaluatePointwiseMeanValuesInPlanes();

  return;
}  // TurbulenceStatisticsCcy::DoTimeSample


/*----------------------------------------------------------------------*

          Compute in plane means of u,u^2 etc. (nodal quantities)

  ----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCcy::EvaluatePointwiseMeanValuesInPlanes()
{
  const int numsubdivisions = 5;

  //----------------------------------------------------------------------
  // sort shellcoordinates and nodeshells
  std::map<double, int, PlaneSortCriterion> countpoints;
  std::map<double, double, PlaneSortCriterion> meanu;
  std::map<double, double, PlaneSortCriterion> meanv;
  std::map<double, double, PlaneSortCriterion> meanw;
  std::map<double, double, PlaneSortCriterion> meanp;
  std::map<double, double, PlaneSortCriterion> meanuu;
  std::map<double, double, PlaneSortCriterion> meanvv;
  std::map<double, double, PlaneSortCriterion> meanww;
  std::map<double, double, PlaneSortCriterion> meanpp;
  std::map<double, double, PlaneSortCriterion> meanuv;
  std::map<double, double, PlaneSortCriterion> meanuw;
  std::map<double, double, PlaneSortCriterion> meanvw;

  std::map<double, double, PlaneSortCriterion> meanc;
  std::map<double, double, PlaneSortCriterion> meancc;

  std::vector<std::map<double, double, PlaneSortCriterion>> meanphi(numscatradofpernode_);
  std::vector<std::map<double, double, PlaneSortCriterion>> meanphiphi(numscatradofpernode_);

  for (std::vector<double>::iterator coord = (*shellcoordinates_).begin();
       coord != (*shellcoordinates_).end(); ++coord)
  {
    double r = *coord;

    meanu.insert(std::pair<double, double>(r, 0.0));
    meanv.insert(std::pair<double, double>(r, 0.0));
    meanw.insert(std::pair<double, double>(r, 0.0));
    meanp.insert(std::pair<double, double>(r, 0.0));

    meanuu.insert(std::pair<double, double>(r, 0.0));
    meanvv.insert(std::pair<double, double>(r, 0.0));
    meanww.insert(std::pair<double, double>(r, 0.0));
    meanpp.insert(std::pair<double, double>(r, 0.0));
    meanuv.insert(std::pair<double, double>(r, 0.0));
    meanuw.insert(std::pair<double, double>(r, 0.0));
    meanvw.insert(std::pair<double, double>(r, 0.0));

    countpoints.insert(std::pair<double, int>(r, 0));

    if (withscatra_)
    {
      meanc.insert(std::pair<double, double>(r, 0.0));
      meancc.insert(std::pair<double, double>(r, 0.0));

      for (int k = 0; k < numscatradofpernode_; ++k)
      {
        meanphi[k].insert(std::pair<double, double>(r, 0.0));
        meanphiphi[k].insert(std::pair<double, double>(r, 0.0));
      }
    }
  }

  // try to cast discretisation to nurbs variant
  // this tells you what kind of computation of
  // samples is required
  DRT::NURBS::NurbsDiscretization* nurbsdis =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*discret_));

  if (nurbsdis == NULL) dserror("Oops. Your discretization is not a NurbsDiscretization.");

  nurbsdis->SetState("velnp", meanvelnp_);

  DRT::NURBS::NurbsDiscretization* scatranurbsdis(NULL);
  if (withscatra_)
  {
    nurbsdis->SetState("scanp", meanscanp_);

    scatranurbsdis = dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*scatradis_));
    if (scatranurbsdis == NULL) dserror("Oops. Your discretization is not a NurbsDiscretization.");

    if (meanfullphinp_ == Teuchos::null)
      dserror("Teuchos::RCP is Teuchos::null");
    else
      scatranurbsdis->SetState("phinp_for_statistics", meanfullphinp_);

    if (not(scatranurbsdis->DofRowMap())->SameAs(meanfullphinp_->Map()))
    {
      scatranurbsdis->DofRowMap()->Print(std::cout);
      meanfullphinp_->Map().Print(std::cout);
      dserror("Global dof numbering in maps does not match");
    }
  }

  // get nurbs dis' knotvector sizes
  std::vector<int> n_x_m_x_l(nurbsdis->Return_n_x_m_x_l(0));

  // get nurbs dis' element numbers
  std::vector<int> nele_x_mele_x_lele(nurbsdis->Return_nele_x_mele_x_lele(0));

  // get the knotvector itself
  Teuchos::RCP<DRT::NURBS::Knotvector> knots = nurbsdis->GetKnotVector();

  // get element map
  const Epetra_Map* elementmap = nurbsdis->ElementRowMap();

  // loop all available elements
  for (int iele = 0; iele < elementmap->NumMyElements(); ++iele)
  {
    // get element pointer
    DRT::Element* const actele = nurbsdis->gElement(elementmap->GID(iele));

    // want to loop all control points of the element,
    // so get the number of points
    const int numnp = actele->NumNode();

    // get the elements control points/nodes
    DRT::Node** nodes = actele->Nodes();

    // acquire weights from nodes
    Epetra_SerialDenseVector weights(numnp);

    for (int inode = 0; inode < numnp; ++inode)
    {
      DRT::NURBS::ControlPoint* cp = dynamic_cast<DRT::NURBS::ControlPoint*>(nodes[inode]);

      weights(inode) = cp->W();
    }
    // get gid, location in the patch
    int gid = actele->Id();

    int patchid = 0;

    std::vector<int> ele_cart_id(3);
    knots->ConvertEleGidToKnotIds(gid, patchid, ele_cart_id);

    // access elements knot span
    std::vector<Epetra_SerialDenseVector> knots(3);
    bool zero_size = (*((*nurbsdis).GetKnotVector())).GetEleKnots(knots, actele->Id());

    // zero sized elements have to be skipped
    if (zero_size)
    {
      continue;
    }

    // get shapefunctions, compute all visualisation point positions
    Epetra_SerialDenseVector nurbs_shape_funct(numnp);

    // extract local values from the global vectors
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;

    actele->LocationVector(*nurbsdis, lm, lmowner, lmstride);

    // extract local values from global vector
    std::vector<double> myvelnp(lm.size());
    DRT::UTILS::ExtractMyValues(*(nurbsdis->GetState("velnp")), myvelnp, lm);

    // create Matrix objects
    LINALG::Matrix<3, 27> evelnp;
    LINALG::Matrix<27, 1> eprenp;

    // insert velocity  into element array
    for (int i = 0; i < 27; ++i)
    {
      const int fi = 4 * i;

      evelnp(0, i) = myvelnp[fi];
      evelnp(1, i) = myvelnp[1 + fi];
      evelnp(2, i) = myvelnp[2 + fi];

      eprenp(i) = myvelnp[3 + fi];
    }

    LINALG::Matrix<1, 27> escanp(true);

    //! scalar at t_(n+1) or t_(n+alpha_F)
    const int nen = 27;  // only quadratic nurbs elements are supported!!
    std::vector<LINALG::Matrix<nen, 1>> ephinp_(numscatradofpernode_);

    if (withscatra_)
    {
      // extract local values from global vector
      std::vector<double> myscanp(lm.size());
      DRT::UTILS::ExtractMyValues(*(nurbsdis->GetState("scanp")), myscanp, lm);

      // insert data into element array (scalar field is stored at pressure dofs)
      for (int i = 0; i < 27; ++i)
      {
        const int fi = 4 * i;
        escanp(0, i) = myscanp[3 + fi];
      }

      // get pointer to corresponding scatra element with identical global id
      DRT::Element* const actscatraele = scatranurbsdis->gElement(gid);
      if (actscatraele == NULL) dserror("could not access transport element with gid %d", gid);

      // extract local values from the global vectors
      std::vector<int> scatralm;
      std::vector<int> scatralmowner;
      std::vector<int> scatralmstride;

      actscatraele->LocationVector(*scatranurbsdis, scatralm, scatralmowner, scatralmstride);

      Teuchos::RCP<const Epetra_Vector> phinp = scatranurbsdis->GetState("phinp_for_statistics");
      if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp' for statistics");
      std::vector<double> myphinp(scatralm.size());
      DRT::UTILS::ExtractMyValues(*phinp, myphinp, scatralm);

      // fill all element arrays
      for (int i = 0; i < nen; ++i)
      {
        for (int k = 0; k < numscatradofpernode_; ++k)
        {
          // split for each transported scalar, insert into element arrays
          ephinp_[k](i, 0) = myphinp[k + (i * numscatradofpernode_)];
        }
      }  // for i
    }

    switch (actele->Shape())
    {
      case DRT::Element::nurbs27:
      {
        LINALG::Matrix<3, 1> vel;

        // element local point position
        Epetra_SerialDenseVector uv(3);

        {
          // standard

          //               v
          //              /
          //  w  7       /   8
          //  ^   +---------+
          //  |  /         /|
          //  | /         / |
          // 5|/        6/  |
          //  +---------+   |
          //  |         |   |
          //  |         |   +
          //  |         |  / 4
          //  |         | /
          //  |         |/
          //  +---------+ ----->u
          // 1           2
          // use r-coordinate of point 1 and 8
          // temporary x vector
          std::vector<double> x(3);

          // point 2
          uv(0) = 1.0;
          uv(1) = -1.0;
          uv(2) = -1.0;
          DRT::NURBS::UTILS::nurbs_get_3D_funct(
              nurbs_shape_funct, uv, knots, weights, actele->Shape());
          for (int isd = 0; isd < 3; ++isd)
          {
            double val = 0;
            for (int inode = 0; inode < numnp; ++inode)
            {
              val += (((nodes[inode])->X())[isd]) * nurbs_shape_funct(inode);
            }
            x[isd] = val;
          }

          const double r = sqrt(x[0] * x[0] + x[1] * x[1]);

          {
            double val = 0;
            for (int inode = 0; inode < numnp; ++inode)
            {
              val += nurbs_shape_funct(inode) * evelnp(0, inode);
            }
            vel(0) = val;

            val = 0;
            for (int inode = 0; inode < numnp; ++inode)
            {
              val += nurbs_shape_funct(inode) * evelnp(1, inode);
            }
            vel(1) = val;

            val = 0;
            for (int inode = 0; inode < numnp; ++inode)
            {
              val += nurbs_shape_funct(inode) * evelnp(2, inode);
            }
            vel(2) = val;

            val = 0;
            for (int inode = 0; inode < numnp; ++inode)
            {
              val += nurbs_shape_funct(inode) * eprenp(inode);
            }
            meanp[r] += val;
            meanpp[r] += val * val;

            if (withscatra_)
            {
              val = 0;
              for (int inode = 0; inode < numnp; ++inode)
              {
                val += nurbs_shape_funct(inode) * escanp(inode);
              }
              meanc[r] += val;
              meancc[r] += val * val;

              for (int k = 0; k < numscatradofpernode_; ++k)
              {
                val = 0;
                for (int inode = 0; inode < numnp; ++inode)
                {
                  val += nurbs_shape_funct(inode) * ((ephinp_[k])(inode));
                }
                (meanphi[k])[r] += val;
                (meanphiphi[k])[r] += val * val;
              }
            }

            std::map<double, int, PlaneSortCriterion>::iterator shell = countpoints.find(r);
            if (shell == countpoints.end())
            {
              dserror("radial coordinate %12.5e was not map\n", r);
            }
            else
            {
              shell->second += 1;
            }

            double uphi = 1.0 / r * (x[0] * vel(1) - x[1] * vel(0));
            double ur = 1.0 / r * (x[0] * vel(0) + x[1] * vel(1));

            meanu[r] += uphi;
            meanv[r] += ur;
            meanw[r] += vel(2);

            meanuu[r] += uphi * uphi;
            meanvv[r] += ur * ur;
            meanww[r] += vel(2) * vel(2);

            meanuv[r] += uphi * ur;
            meanuw[r] += uphi * vel(2);
            meanvw[r] += ur * vel(2);
          }

          for (int rr = 1; rr < numsubdivisions - 1; ++rr)
          {
            uv(1) += 2.0 / (numsubdivisions - 1);

            DRT::NURBS::UTILS::nurbs_get_3D_funct(
                nurbs_shape_funct, uv, knots, weights, actele->Shape());
            for (int isd = 0; isd < 3; ++isd)
            {
              double val = 0;
              for (int inode = 0; inode < numnp; ++inode)
              {
                val += (((nodes[inode])->X())[isd]) * nurbs_shape_funct(inode);
              }
              x[isd] = val;
            }

            const double r = sqrt(x[0] * x[0] + x[1] * x[1]);

            {
              double val = 0;
              for (int inode = 0; inode < numnp; ++inode)
              {
                val += nurbs_shape_funct(inode) * evelnp(0, inode);
              }
              vel(0) = val;

              val = 0;
              for (int inode = 0; inode < numnp; ++inode)
              {
                val += nurbs_shape_funct(inode) * evelnp(1, inode);
              }
              vel(1) = val;

              val = 0;
              for (int inode = 0; inode < numnp; ++inode)
              {
                val += nurbs_shape_funct(inode) * evelnp(2, inode);
              }
              vel(2) = val;

              val = 0;
              for (int inode = 0; inode < numnp; ++inode)
              {
                val += nurbs_shape_funct(inode) * eprenp(inode);
              }
              meanp[r] += val;
              meanpp[r] += val * val;

              if (withscatra_)
              {
                val = 0;
                for (int inode = 0; inode < numnp; ++inode)
                {
                  val += nurbs_shape_funct(inode) * escanp(inode);
                }
                meanc[r] += val;
                meancc[r] += val * val;

                for (int k = 0; k < numscatradofpernode_; ++k)
                {
                  val = 0;
                  for (int inode = 0; inode < numnp; ++inode)
                  {
                    val += nurbs_shape_funct(inode) * ((ephinp_[k])(inode));
                  }
                  (meanphi[k])[r] += val;
                  (meanphiphi[k])[r] += val * val;
                }
              }

              std::map<double, int, PlaneSortCriterion>::iterator shell = countpoints.find(r);
              if (shell == countpoints.end())
              {
                dserror("radial coordinate %12.5e was not map\n", r);
              }
              else
              {
                shell->second += 1;
              }

              double uphi = 1.0 / r * (x[0] * vel(1) - x[1] * vel(0));
              double ur = 1.0 / r * (x[0] * vel(0) + x[1] * vel(1));

              meanu[r] += uphi;
              meanv[r] += ur;
              meanw[r] += vel(2);

              meanuu[r] += uphi * uphi;
              meanvv[r] += ur * ur;
              meanww[r] += vel(2) * vel(2);

              meanuv[r] += uphi * ur;
              meanuw[r] += uphi * vel(2);
              meanvw[r] += ur * vel(2);
            }
          }

          // set upper point of element, too (only for last layer)
          if (ele_cart_id[1] + 1 == nele_x_mele_x_lele[1])
          {
            // point 4
            uv(0) = 1.0;
            uv(1) = 1.0;
            uv(2) = -1.0;
            DRT::NURBS::UTILS::nurbs_get_3D_funct(
                nurbs_shape_funct, uv, knots, weights, actele->Shape());
            for (int isd = 0; isd < 3; ++isd)
            {
              double val = 0;
              for (int inode = 0; inode < numnp; ++inode)
              {
                val += (((nodes[inode])->X())[isd]) * nurbs_shape_funct(inode);
              }
              x[isd] = val;
            }


            const double r = sqrt(x[0] * x[0] + x[1] * x[1]);

            {
              double val = 0;
              for (int inode = 0; inode < numnp; ++inode)
              {
                val += nurbs_shape_funct(inode) * evelnp(0, inode);
              }
              vel(0) = val;

              val = 0;
              for (int inode = 0; inode < numnp; ++inode)
              {
                val += nurbs_shape_funct(inode) * evelnp(1, inode);
              }
              vel(1) = val;

              val = 0;
              for (int inode = 0; inode < numnp; ++inode)
              {
                val += nurbs_shape_funct(inode) * evelnp(2, inode);
              }
              vel(2) = val;

              val = 0;
              for (int inode = 0; inode < numnp; ++inode)
              {
                val += nurbs_shape_funct(inode) * eprenp(inode);
              }
              meanp[r] += val;
              meanpp[r] += val * val;

              if (withscatra_)
              {
                val = 0;
                for (int inode = 0; inode < numnp; ++inode)
                {
                  val += nurbs_shape_funct(inode) * escanp(inode);
                }
                meanc[r] += val;
                meancc[r] += val * val;

                for (int k = 0; k < numscatradofpernode_; ++k)
                {
                  val = 0;
                  for (int inode = 0; inode < numnp; ++inode)
                  {
                    val += nurbs_shape_funct(inode) * ((ephinp_[k])(inode));
                  }
                  (meanphi[k])[r] += val;
                  (meanphiphi[k])[r] += val * val;
                }
              }

              std::map<double, int, PlaneSortCriterion>::iterator shell = countpoints.find(r);
              if (shell == countpoints.end())
              {
                dserror("radial coordinate %12.5e was not map\n", r);
              }
              else
              {
                shell->second += 1;
              }

              double uphi = 1.0 / r * (x[0] * vel(1) - x[1] * vel(0));
              double ur = 1.0 / r * (x[0] * vel(0) + x[1] * vel(1));

              meanu[r] += uphi;
              meanv[r] += ur;
              meanw[r] += vel(2);

              meanuu[r] += uphi * uphi;
              meanvv[r] += ur * ur;
              meanww[r] += vel(2) * vel(2);

              meanuv[r] += uphi * ur;
              meanuw[r] += uphi * vel(2);
              meanvw[r] += ur * vel(2);
            }
          }
        }
        break;
      }
      default:
        dserror(
            "Unknown element shape for a nurbs element or nurbs type not valid for turbulence "
            "calculation\n");
    }
  }  // end element loop

  // clean up
  nurbsdis->ClearState();
  if (scatranurbsdis != NULL)
  {
    scatranurbsdis->ClearState();
  }

  // communicate results among processors
  int size = countpoints.size();
  int rr;

  // collect number of samples
  std::vector<int> lpointcount;
  std::vector<int> pointcount(size);

  for (std::map<double, int, PlaneSortCriterion>::iterator shell = countpoints.begin();
       shell != countpoints.end(); ++shell)
  {
    lpointcount.push_back(shell->second);
  }
  discret_->Comm().SumAll(&(lpointcount[0]), &(pointcount[0]), size);

  // collect number of samples
  std::vector<double> lmeanu;
  std::vector<double> lmeanv;
  std::vector<double> lmeanw;
  std::vector<double> lmeanp;

  std::vector<double> lmeanuu;
  std::vector<double> lmeanvv;
  std::vector<double> lmeanww;
  std::vector<double> lmeanpp;

  std::vector<double> lmeanuv;
  std::vector<double> lmeanuw;
  std::vector<double> lmeanvw;

  std::vector<double> lmeanc;
  std::vector<double> lmeancc;

  std::vector<double> lmeanphi;
  std::vector<double> lmeanphiphi;

  std::vector<double> gmeanu(size);
  std::vector<double> gmeanv(size);
  std::vector<double> gmeanw(size);
  std::vector<double> gmeanp(size);

  std::vector<double> gmeanuu(size);
  std::vector<double> gmeanvv(size);
  std::vector<double> gmeanww(size);
  std::vector<double> gmeanpp(size);

  std::vector<double> gmeanuv(size);
  std::vector<double> gmeanuw(size);
  std::vector<double> gmeanvw(size);

  std::vector<double> gmeanc(size);
  std::vector<double> gmeancc(size);

  std::vector<double> gmeanphi(size);
  std::vector<double> gmeanphiphi(size);

  rr = 0;
  for (std::map<double, double, PlaneSortCriterion>::iterator shell = meanu.begin();
       shell != meanu.end(); ++shell)
  {
    if (fabs(pointcount[rr]) < 1e-6)
    {
      dserror("zero pointcount during computation of averages, layer %d\n", rr);
    }

    lmeanu.push_back(shell->second / pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanu[0]), &((gmeanu)[0]), size);

  rr = 0;
  for (std::map<double, double, PlaneSortCriterion>::iterator shell = meanv.begin();
       shell != meanv.end(); ++shell)
  {
    lmeanv.push_back(shell->second / pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanv[0]), &((gmeanv)[0]), size);

  rr = 0;
  for (std::map<double, double, PlaneSortCriterion>::iterator shell = meanw.begin();
       shell != meanw.end(); ++shell)
  {
    lmeanw.push_back(shell->second / pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanw[0]), &((gmeanw)[0]), size);

  rr = 0;
  for (std::map<double, double, PlaneSortCriterion>::iterator shell = meanp.begin();
       shell != meanp.end(); ++shell)
  {
    lmeanp.push_back(shell->second / pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanp[0]), &((gmeanp)[0]), size);

  rr = 0;
  for (std::map<double, double, PlaneSortCriterion>::iterator shell = meanuu.begin();
       shell != meanuu.end(); ++shell)
  {
    lmeanuu.push_back(shell->second / pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanuu[0]), &((gmeanuu)[0]), size);

  rr = 0;
  for (std::map<double, double, PlaneSortCriterion>::iterator shell = meanvv.begin();
       shell != meanvv.end(); ++shell)
  {
    lmeanvv.push_back(shell->second / pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanvv[0]), &((gmeanvv)[0]), size);

  rr = 0;
  for (std::map<double, double, PlaneSortCriterion>::iterator shell = meanww.begin();
       shell != meanww.end(); ++shell)
  {
    lmeanww.push_back(shell->second / pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanww[0]), &((gmeanww)[0]), size);

  rr = 0;
  for (std::map<double, double, PlaneSortCriterion>::iterator shell = meanpp.begin();
       shell != meanpp.end(); ++shell)
  {
    lmeanpp.push_back(shell->second / pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanpp[0]), &((gmeanpp)[0]), size);

  rr = 0;
  for (std::map<double, double, PlaneSortCriterion>::iterator shell = meanuv.begin();
       shell != meanuv.end(); ++shell)
  {
    lmeanuv.push_back(shell->second / pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanuv[0]), &((gmeanuv)[0]), size);

  rr = 0;
  for (std::map<double, double, PlaneSortCriterion>::iterator shell = meanuw.begin();
       shell != meanuw.end(); ++shell)
  {
    lmeanuw.push_back(shell->second / pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanuw[0]), &((gmeanuw)[0]), size);


  rr = 0;
  for (std::map<double, double, PlaneSortCriterion>::iterator shell = meanvw.begin();
       shell != meanvw.end(); ++shell)
  {
    lmeanvw.push_back(shell->second / pointcount[rr]);
    ++rr;
  }
  discret_->Comm().SumAll(&(lmeanvw[0]), &(gmeanvw[0]), size);

  for (int mm = 0; mm < size; ++mm)
  {
    (*pointsumu_)[mm] += gmeanu[mm];
    (*pointsumv_)[mm] += gmeanv[mm];
    (*pointsumw_)[mm] += gmeanw[mm];
    (*pointsump_)[mm] += gmeanp[mm];

    (*pointsumuu_)[mm] += gmeanuu[mm];
    (*pointsumvv_)[mm] += gmeanvv[mm];
    (*pointsumww_)[mm] += gmeanww[mm];
    (*pointsumpp_)[mm] += gmeanpp[mm];

    (*pointsumuv_)[mm] += gmeanuv[mm];
    (*pointsumuw_)[mm] += gmeanuw[mm];
    (*pointsumvw_)[mm] += gmeanvw[mm];
  }

  if (withscatra_)
  {
    rr = 0;
    for (std::map<double, double, PlaneSortCriterion>::iterator shell = meanc.begin();
         shell != meanc.end(); ++shell)
    {
      lmeanc.push_back(shell->second / pointcount[rr]);
      ++rr;
    }
    discret_->Comm().SumAll(&(lmeanc[0]), &(gmeanc[0]), size);

    rr = 0;
    for (std::map<double, double, PlaneSortCriterion>::iterator shell = meancc.begin();
         shell != meancc.end(); ++shell)
    {
      lmeancc.push_back(shell->second / pointcount[rr]);
      ++rr;
    }
    discret_->Comm().SumAll(&(lmeancc[0]), &(gmeancc[0]), size);

    for (int mm = 0; mm < size; ++mm)
    {
      (*pointsumc_)[mm] += gmeanc[mm];
      (*pointsumcc_)[mm] += gmeancc[mm];
    }

    // safety checks
    if ((size != pointsumphi_->M()) or (size != pointsumphiphi_->M()))
      dserror("Size mismatch: size = %d <-> M = %d", size, pointsumphi_->M());
    if ((numscatradofpernode_ != pointsumphi_->N()) or
        (numscatradofpernode_ != pointsumphiphi_->N()))
      dserror("Size mismatch: numdof = %d <-> N = %d", numscatradofpernode_, pointsumphi_->N());

    // loop all available scatra fields
    for (int k = 0; k < numscatradofpernode_; ++k)
    {
      // sequential reuse of this vector + push back actions below require this:
      lmeanphi.resize(0);
      lmeanphiphi.resize(0);

      rr = 0;
      for (std::map<double, double, PlaneSortCriterion>::iterator shell = meanphi[k].begin();
           shell != meanphi[k].end(); ++shell)
      {
        lmeanphi.push_back(shell->second / pointcount[rr]);
        gmeanphi[rr] = 0.0;  // initialize due to sequential reuse of this vector
        ++rr;
      }
      discret_->Comm().SumAll(&(lmeanphi[0]), &(gmeanphi[0]), size);

      rr = 0;
      for (std::map<double, double, PlaneSortCriterion>::iterator shell = meanphiphi[k].begin();
           shell != meanphiphi[k].end(); ++shell)
      {
        lmeanphiphi.push_back(shell->second / pointcount[rr]);
        gmeanphiphi[rr] = 0.0;  // initialize due to sequential reuse of this vector
        ++rr;
      }
      discret_->Comm().SumAll(&(lmeanphiphi[0]), &(gmeanphiphi[0]), size);

      // insert values for scalar field k
      for (int mm = 0; mm < size; ++mm)
      {
        (*pointsumphi_)(mm, k) += gmeanphi[mm];
        (*pointsumphiphi_)(mm, k) += gmeanphiphi[mm];
      }

    }  // loop over scalar fields k


  }  // if(withscatra)

  return;
}  // TurbulenceStatisticsCcy::EvaluatePointwiseMeanValuesInPlanes()


/*----------------------------------------------------------------------*

       Compute a time average of the mean values over all steps
          since the last output. Dump the result to file.

  ----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCcy::TimeAverageMeansAndOutputOfStatistics(int step)
{
  if (numsamp_ == 0)
  {
    dserror("No samples to do time average");
  }

  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RCP<std::ofstream> log;
  if (discret_->Comm().MyPID() == 0)
  {
    std::string s(statistics_outfilename_);
    s.append(".flow_statistics");

    log = Teuchos::rcp(new std::ofstream(s.c_str(), std::ios::app));
    (*log) << "\n\n\n";
    (*log) << "# Statistics record " << countrecord_;
    (*log) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n";

    (*log) << "#|-------------------";
    (*log) << "---------------------------------------------------";
    (*log) << "--point based (interpolated)-------------------------";
    (*log) << "----------------------------------------------------------|";
    (*log) << "\n";

    (*log) << "#       y        ";
    (*log) << "      u_theta     ";
    (*log) << "       u_r       ";
    (*log) << "       u_z       ";
    (*log) << "        p        ";
    (*log) << " u_theta*u_theta ";
    (*log) << "     u_r*u_r     ";
    (*log) << "     u_z*u_z     ";
    (*log) << "       p*p       ";
    (*log) << "   u_theta*u_r   ";
    (*log) << "   u_theta*u_z   ";
    (*log) << "     u_r*u_z     ";
    if (withscatra_)
    {
      (*log) << "          c          ";
      (*log) << "         c*c         ";

      for (int k = 0; k < numscatradofpernode_; k++)
      {
        (*log) << "          c" << k + 1 << "          ";
        (*log) << "       c" << k + 1 << "*c" << k + 1 << "        ";
      }
    }
    (*log) << "\n";
    (*log) << std::scientific;
    for (unsigned i = 0; i < shellcoordinates_->size(); ++i)
    {
      // y and y+
      (*log) << " " << std::setw(14) << std::setprecision(7) << (*shellcoordinates_)[i];

      // pointwise means
      (*log) << "    " << std::setw(13) << std::setprecision(6) << (*pointsumu_)[i] / numsamp_;
      (*log) << "    " << std::setw(13) << std::setprecision(6) << (*pointsumv_)[i] / numsamp_;
      (*log) << "    " << std::setw(13) << std::setprecision(6) << (*pointsumw_)[i] / numsamp_;
      (*log) << "    " << std::setw(13) << std::setprecision(6) << (*pointsump_)[i] / numsamp_;
      (*log) << "    " << std::setw(13) << std::setprecision(6) << (*pointsumuu_)[i] / numsamp_;
      (*log) << "    " << std::setw(13) << std::setprecision(6) << (*pointsumvv_)[i] / numsamp_;
      (*log) << "    " << std::setw(13) << std::setprecision(6) << (*pointsumww_)[i] / numsamp_;
      (*log) << "    " << std::setw(13) << std::setprecision(6) << (*pointsumpp_)[i] / numsamp_;
      (*log) << "    " << std::setw(13) << std::setprecision(6) << (*pointsumuv_)[i] / numsamp_;
      (*log) << "    " << std::setw(13) << std::setprecision(6) << (*pointsumuw_)[i] / numsamp_;
      (*log) << "    " << std::setw(13) << std::setprecision(6) << (*pointsumvw_)[i] / numsamp_;
      if (withscatra_)
      {
        (*log) << "    " << std::setw(17) << std::setprecision(10) << (*pointsumc_)[i] / numsamp_;
        (*log) << "    " << std::setw(17) << std::setprecision(10) << (*pointsumcc_)[i] / numsamp_;

        for (int k = 0; k < numscatradofpernode_; k++)
        {
          (*log) << "    " << std::setw(17) << std::setprecision(10)
                 << ((*pointsumphi_)(i, k)) / numsamp_;
          (*log) << "    " << std::setw(17) << std::setprecision(10)
                 << ((*pointsumphiphi_)(i, k)) / numsamp_;
        }
      }
      (*log) << "\n";
    }
    log->flush();
  }  // end myrank 0


  // log was written, so increase counter for records
  countrecord_++;

  return;

}  // TurbulenceStatisticsCcy::TimeAverageMeansAndOutputOfStatistics


/*----------------------------------------------------------------------*

                  Reset sums and number of samples to 0

  ----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCcy::ClearStatistics()
{
  // reset the number of samples
  numsamp_ = 0;

  // reset integral and pointwise averages
  for (unsigned i = 0; i < shellcoordinates_->size(); ++i)
  {
    (*pointsumu_)[i] = 0;
    (*pointsumv_)[i] = 0;
    (*pointsumw_)[i] = 0;
    (*pointsump_)[i] = 0;

    (*pointsumuu_)[i] = 0;
    (*pointsumvv_)[i] = 0;
    (*pointsumww_)[i] = 0;
    (*pointsumpp_)[i] = 0;

    (*pointsumuv_)[i] = 0;
    (*pointsumuw_)[i] = 0;
    (*pointsumvw_)[i] = 0;
  }

  meanvelnp_->PutScalar(0.0);

  if (withscatra_)
  {
    // reset integral and pointwise averages
    for (size_t i = 0; i < shellcoordinates_->size(); ++i)
    {
      (*pointsumc_)[i] = 0.0;
      (*pointsumcc_)[i] = 0.0;
    }

    if (meanscanp_ == Teuchos::null)
      dserror("meanscanp_ is Teuchos::null");
    else
      meanscanp_->PutScalar(0.0);

    if (meanfullphinp_ != Teuchos::null)
    {
      meanfullphinp_->PutScalar(0.0);

      // ToDo Is is a good way to initialize everything to zero??
      // Use Shape() instead???
      pointsumphi_->Scale(0.0);
      pointsumphiphi_->Scale(0.0);
    }
  }

  return;
}  // TurbulenceStatisticsCcy::ClearStatistics


/*----------------------------------------------------------------------

Add results from scalar transport fields to statistics

----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCcy::AddScaTraResults(
    Teuchos::RCP<DRT::Discretization> scatradis, Teuchos::RCP<Epetra_Vector> phinp)
{
  if (withscatra_)
  {
    if (scatradis == Teuchos::null)
      dserror("Halleluja.");
    else
      scatradis_ = scatradis;  // now we have access

    // we do not have to cast to a NURBSDiscretization here!
    meanfullphinp_ = LINALG::CreateVector(*(scatradis_->DofRowMap()), true);
    numscatradofpernode_ = scatradis_->NumDof(scatradis_->lRowNode(0));

    // now we know about the number of scatra dofs and can allocate:
    int size = shellcoordinates_->size();
    pointsumphi_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(size, numscatradofpernode_));
    pointsumphiphi_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(size, numscatradofpernode_));

    if (discret_->Comm().MyPID() == 0)
    {
      std::cout << std::endl
                << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
      std::cout << "TurbulenceStatisticsCcy:    added access to ScaTra results" << std::endl;
      std::cout << " | nodeshellsize       = " << nodeshells_->size() << std::endl
                << " | numshellcoordinates = " << size << " (4 subdivisions per element)"
                << std::endl
                << " | numscatradofpernode = " << numscatradofpernode_ << std::endl;
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl
                << std::endl;
    }
  }
  else
  {
    if (discret_->Comm().MyPID() == 0)
    {
      std::cout << "------------------------------------------------------------" << std::endl;
      std::cout << "TurbulenceStatisticsCcy: NO access to ScaTra results !" << std::endl;
      std::cout << "------------------------------------------------------------" << std::endl;
    }
  }

  return;
}  // FLD::TurbulenceStatisticsCcy::AddScaTraResults
