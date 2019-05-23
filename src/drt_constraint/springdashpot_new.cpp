/*!----------------------------------------------------------------------

\brief Methods for spring and dashpot constraints / boundary conditions:

\level 2

\maintainer Amadeus Gebauer

*----------------------------------------------------------------------*/


#include "springdashpot_new.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include <iostream>
#include "../drt_io/io_pstream.H"  // has to go before io.H
#include "../drt_io/io.H"
#include <Epetra_Time.h>

#include "../drt_adapter/adapter_coupling_nonlin_mortar.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_contact/contact_interface.H"
#include "../drt_lib/drt_condition_utils.H"

/*----------------------------------------------------------------------*
 |                                                         pfaller Apr15|
 *----------------------------------------------------------------------*/
UTILS::SpringDashpotNew::SpringDashpotNew(
    Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<DRT::Condition> cond)
    : actdisc_(dis),
      spring_(cond),
      stiff_tens_((*spring_->Get<std::vector<double>>("stiff"))[0]),
      stiff_comp_((*spring_->Get<std::vector<double>>("stiff"))[0]),
      offset_((*spring_->Get<std::vector<double>>("disploffset"))[0]),
      viscosity_((*spring_->Get<std::vector<double>>("visco"))[0]),
      coupling_(spring_->GetInt("coupling id")),
      nodes_(spring_->Nodes()),
      area_(),
      gap0_(),
      gap_(),
      gapdt_(),
      dgap_(),
      normals_(),
      dnormals_(),
      offset_prestr_(),
      offset_prestr_new_(Teuchos::null)
{
  offset_prestr_new_ = Teuchos::rcp(new Epetra_Vector(*actdisc_->DofRowMap()));
  offset_prestr_new_->PutScalar(0.0);

  // set type of this spring
  SetSpringType();

  if (springtype_ != cursurfnormal && coupling_ >= 0)
    dserror(
        "Coupling of spring dashpot to reference surface only possible for DIRECTION "
        "cursurfnormal.");

  if (springtype_ == cursurfnormal && coupling_ == -1)
    dserror("Coupling id necessary for DIRECTION cursurfnormal.");

  if (springtype_ == cursurfnormal && viscosity_ != 0.)
    dserror("Viscous damping not yet implemented for DIRECTION cursurfnormal.");

  // ToDo: delete rest until return statement!
  // get geometry
  std::map<int, Teuchos::RCP<DRT::Element>>& geom = spring_->Geometry();

  // get normal vectors if necessary
  if (springtype_ == cursurfnormal)
  {
    // calculate nodal area
    if (!actdisc_->Comm().MyPID()) IO::cout << "Computing area for spring dashpot condition...\n";
    GetArea(geom);
    InitializeCurSurfNormal();
  }

  // ToDo: do we really need this??
  // initialize prestressing offset
  InitializePrestrOffset();

  return;
}

// NEW version, consistently integrated over element surface!!
/*----------------------------------------------------------------------*
 * Integrate a Surface Robin boundary condition (public)       mhv 08/16|
 * ---------------------------------------------------------------------*/
void UTILS::SpringDashpotNew::EvaluateRobin(Teuchos::RCP<LINALG::SparseMatrix> stiff,
    Teuchos::RCP<Epetra_Vector> fint, const Teuchos::RCP<const Epetra_Vector> disp,
    const Teuchos::RCP<const Epetra_Vector> velo, Teuchos::ParameterList p)
{
  // reset last Newton step
  springstress_.clear();

  const bool assvec = fint != Teuchos::null;
  const bool assmat = stiff != Teuchos::null;

  actdisc_->ClearState();
  actdisc_->SetState("displacement", disp);
  actdisc_->SetState("velocity", velo);
  actdisc_->SetState("offset_prestress", offset_prestr_new_);

  // get values and switches from the condition
  const std::vector<int>* onoff = spring_->Get<std::vector<int>>("onoff");
  const std::vector<double>* springstiff = spring_->Get<std::vector<double>>("stiff");
  const std::vector<double>* dashpotvisc = spring_->Get<std::vector<double>>("visco");
  const std::vector<double>* disploffset = spring_->Get<std::vector<double>>("disploffset");
  const std::string* direction = spring_->Get<std::string>("direction");

  // time-integration factor for stiffness contribution of dashpot, d(v_{n+1})/d(d_{n+1})
  const double time_fac = p.get("time_fac", 0.0);

  Teuchos::ParameterList params;
  params.set("action", "calc_struct_robinforcestiff");
  params.set("onoff", onoff);
  params.set("springstiff", springstiff);
  params.set("dashpotvisc", dashpotvisc);
  params.set("disploffset", disploffset);
  params.set("time_fac", time_fac);
  params.set("direction", direction);

  std::map<int, Teuchos::RCP<DRT::Element>>& geom = spring_->Geometry();

  // if (geom.empty()) dserror("evaluation of condition with empty geometry");
  // no check for empty geometry here since in parallel computations
  // can exist processors which do not own a portion of the elements belonging
  // to the condition geometry
  std::map<int, Teuchos::RCP<DRT::Element>>::iterator curr;
  for (curr = geom.begin(); curr != geom.end(); ++curr)
  {
    // get element location vector and ownerships
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    curr->second->LocationVector(*actdisc_, lm, lmowner, lmstride);

    const int eledim = (int)lm.size();

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2;
    Epetra_SerialDenseVector elevector3;

    elevector1.Size(eledim);
    elevector2.Size(eledim);
    elevector3.Size(eledim);
    elematrix1.Shape(eledim, eledim);

    int err = curr->second->Evaluate(
        params, *actdisc_, lm, elematrix1, elematrix2, elevector1, elevector2, elevector3);
    if (err) dserror("error while evaluating elements");

    if (assvec) LINALG::Assemble(*fint, elevector1, lm, lmowner);
    if (assmat) stiff->Assemble(curr->second->Id(), lmstride, elematrix1, lm, lmowner);

    // save spring stress for postprocessing
    const int numdim = 3;
    const int numdf = 3;
    std::vector<double> stress(numdim, 0.0);

    for (int node = 0; node < curr->second->NumNode(); ++node)
    {
      for (int dim = 0; dim < numdim; dim++) stress[dim] = elevector3[node * numdf + dim];
      springstress_.insert(
          std::pair<int, std::vector<double>>(curr->second->NodeIds()[node], stress));
    }
  } /* end of loop over geometry */

  return;
}


// old version, NOT consistently integrated over element surface!!
/*----------------------------------------------------------------------*
 |                                                         pfaller Mar16|
 *----------------------------------------------------------------------*/
void UTILS::SpringDashpotNew::EvaluateForce(Epetra_Vector& fint,
    const Teuchos::RCP<const Epetra_Vector> disp, const Teuchos::RCP<const Epetra_Vector> velo)
{
  if (disp == Teuchos::null) dserror("Cannot find displacement state in discretization");

  if (springtype_ == cursurfnormal) GetCurNormals(disp);

  // loop nodes of current condition
  const std::vector<int>& nds = *nodes_;
  for (int j = 0; j < (int)nds.size(); ++j)
  {
    // nodes owned by processor
    if (actdisc_->NodeRowMap()->MyGID(nds[j]))
    {
      int gid = nds[j];
      DRT::Node* node = actdisc_->gNode(gid);
      if (!node) dserror("Cannot find global node %d", gid);

      // get nodal values
      const double nodalarea = area_[gid];               // nodal area
      const std::vector<double> normal = normals_[gid];  // normalized nodal normal
      const std::vector<double> offsetprestr =
          offset_prestr_[gid];  // get nodal displacement values of last time step for MULF offset

      const int numdof = actdisc_->NumDof(0, node);
      assert(numdof == 3);
      std::vector<int> dofs = actdisc_->Dof(0, node);

      // initialize
      double gap = 0.;          // displacement
                                //      double gapdt = 0.; // velocity // unused ?!?
      double springstiff = 0.;  // spring stiffness

      // calculation of normals and displacements differs for each spring variant
      switch (springtype_)
      {
        case xyz:  // spring dashpot acts in every surface dof direction
          dserror("You should not be here! Use the new consistent EvaluateRobin routine!!!");
          break;

        case refsurfnormal:  // spring dashpot acts in refnormal direction
          dserror("You should not be here! Use the new consistent EvaluateRobin routine!!!");
          break;

        case cursurfnormal:  // spring dashpot acts in curnormal direction
          // spring displacement
          gap = gap_[gid];
          //        gapdt = gapdt_[gid]; // unused ?!?

          // select spring stiffnes
          springstiff = SelectStiffness(gap);

          // assemble into residual vector
          std::vector<double> out_vec(numdof, 0.);
          for (int k = 0; k < numdof; ++k)
          {
            // force
            const double val =
                -nodalarea * springstiff * (gap - offset_) * normal[k];  // + viscosity_*gapdt
            const int err = fint.SumIntoGlobalValues(1, &val, &dofs[k]);
            if (err) dserror("SumIntoGlobalValues failed!");

            // store spring stress for output
            out_vec[k] = springstiff * (gap - offsetprestr[k] - offset_) * normal[k];
          }
          // add to output
          springstress_.insert(std::pair<int, std::vector<double>>(gid, out_vec));
          break;
      }
    }  // node owned by processor
  }    // loop over nodes

  return;
}


// old version, NOT consistently integrated over element surface!!
/*----------------------------------------------------------------------*
 |                                                         pfaller mar16|
 *----------------------------------------------------------------------*/
void UTILS::SpringDashpotNew::EvaluateForceStiff(LINALG::SparseMatrix& stiff, Epetra_Vector& fint,
    const Teuchos::RCP<const Epetra_Vector> disp, const Teuchos::RCP<const Epetra_Vector> velo,
    Teuchos::ParameterList p)
{
  if (disp == Teuchos::null) dserror("Cannot find displacement state in discretization");

  if (springtype_ == cursurfnormal)
  {
    GetCurNormals(disp);
    stiff.UnComplete();  // sparsity pattern might change
  }

  // time-integration factor for stiffness contribution of dashpot, d(v_{n+1})/d(d_{n+1})
  //  const double time_fac = p.get("time_fac",1.0); // unused for cursurfnormal ?!?

  // loop nodes of current condition
  const std::vector<int>& nds = *nodes_;
  for (int j = 0; j < (int)nds.size(); ++j)
  {
    // nodes owned by processor
    if (actdisc_->NodeRowMap()->MyGID(nds[j]))
    {
      int gid = nds[j];
      DRT::Node* node = actdisc_->gNode(gid);
      if (!node) dserror("Cannot find global node %d", gid);

      // get nodal values
      const double nodalarea = area_[gid];               // nodal area
      const std::vector<double> normal = normals_[gid];  // normalized nodal normal
      const std::vector<double> offsetprestr =
          offset_prestr_[gid];  // get nodal displacement values of last time step for MULF offset

      const int numdof = actdisc_->NumDof(0, node);
      assert(numdof == 3);
      std::vector<int> dofs = actdisc_->Dof(0, node);

      // initialize
      double gap = 0.;          // displacement
      double gapdt = 0.;        // velocity
      double springstiff = 0.;  // spring stiffness

      // calculation of normals and displacements differs for each spring variant
      switch (springtype_)
      {
        case xyz:  // spring dashpot acts in every surface dof direction
          dserror("You should not be here! Use the new consistent EvaluateRobin routine!!!");
          break;

        case refsurfnormal:  // spring dashpot acts in refnormal direction
          dserror("You should not be here! Use the new consistent EvaluateRobin routine!!!");
          break;

        case cursurfnormal:  // spring dashpot acts in curnormal direction
          // spring displacement
          gap = gap_[gid];
          gapdt = gapdt_[gid];

          // select spring stiffnes
          springstiff = SelectStiffness(gap);

          // assemble into residual vector and stiffness matrix
          std::vector<double> out_vec(numdof, 0.);
          for (int k = 0; k < numdof; ++k)
          {
            // force
            const double val =
                -nodalarea * (springstiff * (gap - offset_) + viscosity_ * gapdt) * normal[k];
            const int err = fint.SumIntoGlobalValues(1, &val, &dofs[k]);
            if (err) dserror("SumIntoGlobalValues failed!");

            // stiffness
            std::map<int, double> dgap = dgap_[gid];
            std::vector<GEN::pairedvector<int, double>> dnormal = dnormals_[gid];

            // check if projection exists
            if (!dnormal.empty() && !dgap.empty())
            {
              // linearize gap
              for (std::map<int, double>::iterator i = dgap.begin(); i != dgap.end(); ++i)
              {
                const double dval = -nodalarea * springstiff * (i->second) * normal[k];
                stiff.Assemble(dval, dofs[k], i->first);
              }

              // linearize normal
              for (GEN::pairedvector<int, double>::iterator i = dnormal[k].begin();
                   i != dnormal[k].end(); ++i)
              {
                const double dval = -nodalarea * springstiff * (gap - offset_) * (i->second);
                stiff.Assemble(dval, dofs[k], i->first);
              }
            }
            // store negative value of internal force for output (=reaction force)
            out_vec[k] = -val;
          }
          // add to output
          springstress_.insert(std::pair<int, std::vector<double>>(gid, out_vec));
          break;
      }
    }  // node owned by processor
  }    // loop over nodes

  if (springtype_ == cursurfnormal) stiff.Complete();  // sparsity pattern might have changed

  return;
}

/*----------------------------------------------------------------------*
 |                                                         pfaller Mar16|
 *----------------------------------------------------------------------*/
void UTILS::SpringDashpotNew::ResetNewton()
{
  // all springs
  gap_.clear();
  gapdt_.clear();
  springstress_.clear();

  // only curnormal
  if (springtype_ == cursurfnormal)
  {
    dgap_.clear();
    normals_.clear();
    dnormals_.clear();
  }
}

/*----------------------------------------------------------------------*
 |                                                             mhv 12/15|
 *----------------------------------------------------------------------*/
void UTILS::SpringDashpotNew::ResetPrestress(Teuchos::RCP<Epetra_Vector> dis)
{
  // loop nodes of current condition
  const std::vector<int>& nds = *nodes_;
  for (int j = 0; j < (int)nds.size(); ++j)
  {
    // nodes owned by processor
    if (actdisc_->NodeRowMap()->MyGID(nds[j]))
    {
      int gid = nds[j];
      DRT::Node* node = actdisc_->gNode(gid);
      if (!node) dserror("Cannot find global node %d", gid);

      const int numdof = actdisc_->NumDof(0, node);
      assert(numdof == 3);
      std::vector<int> dofs = actdisc_->Dof(0, node);

      // initialize. calculation of displacements differs for each spring variant
      std::vector<double> uoff(numdof, 0.0);  // displacement vector of condition nodes
      if (springtype_ == refsurfnormal || springtype_ == xyz)
      {
        for (int k = 0; k < numdof; ++k)
        {
          uoff[k] = -(*dis)[dis->Map().LID(dofs[k])];  // extract nodal displacement to be offset
        }

        std::vector<double> offpr = offset_prestr_[gid];
        offset_prestr_.erase(gid);

        for (int k = 0; k < numdof; ++k)
        {
          uoff[k] += offpr[k];  // accumulate displacements to be offset
        }
      }
      offset_prestr_.insert(std::pair<int, std::vector<double>>(gid, uoff));

    }  // node owned by processor
  }    // loop over nodes
}

/*----------------------------------------------------------------------*
 |                                                             mhv 12/15|
 *----------------------------------------------------------------------*/
void UTILS::SpringDashpotNew::SetRestart(Teuchos::RCP<Epetra_Vector> vec)
{
  offset_prestr_new_->Update(1.0, *vec, 0.0);

  return;
}

/*----------------------------------------------------------------------*
 |                                                             mhv 12/15|
 *----------------------------------------------------------------------*/
void UTILS::SpringDashpotNew::SetRestartOld(Teuchos::RCP<Epetra_MultiVector> vec)
{
  // loop nodes of current condition
  const std::vector<int>& nds = *nodes_;
  for (int j = 0; j < (int)nds.size(); ++j)
  {
    // nodes owned by processor
    if (actdisc_->NodeRowMap()->MyGID(nds[j]))
    {
      int gid = nds[j];
      DRT::Node* node = actdisc_->gNode(gid);
      if (!node) dserror("Cannot find global node %d", gid);

      const int numdof = actdisc_->NumDof(0, node);
      assert(numdof == 3);
      std::vector<int> dofs = actdisc_->Dof(0, node);


      if (springtype_ == refsurfnormal || springtype_ == xyz)
      {
        // import spring offset length
        for (std::map<int, std::vector<double>>::iterator i = offset_prestr_.begin();
             i != offset_prestr_.end(); ++i)
        {
          // global id -> local id
          const int lid = vec->Map().LID(i->first);
          // local id on processor
          if (lid >= 0)
          {
            // copy all components of spring offset length vector
            (i->second)[0] = (*(*vec)(0))[lid];
            (i->second)[1] = (*(*vec)(1))[lid];
            (i->second)[2] = (*(*vec)(2))[lid];
          }
        }
      }

    }  // node owned by processor
  }    // loop over nodes
}

/*----------------------------------------------------------------------*
 |                                                         pfaller Jan14|
 *----------------------------------------------------------------------*/
void UTILS::SpringDashpotNew::OutputGapNormal(Teuchos::RCP<Epetra_Vector>& gap,
    Teuchos::RCP<Epetra_MultiVector>& normals, Teuchos::RCP<Epetra_MultiVector>& stress) const
{
  // export gap function
  for (std::map<int, double>::const_iterator i = gap_.begin(); i != gap_.end(); ++i)
  {
    // global id -> local id
    const int lid = gap->Map().LID(i->first);
    // local id on processor
    if (lid >= 0) (*gap)[lid] += i->second;
  }

  // export normal
  for (std::map<int, std::vector<double>>::const_iterator i = normals_.begin(); i != normals_.end();
       ++i)
  {
    // global id -> local id
    const int lid = normals->Map().LID(i->first);
    // local id on processor
    if (lid >= 0)
    {
      // copy all components of normal vector
      (*(*normals)(0))[lid] += (i->second).at(0);
      (*(*normals)(1))[lid] += (i->second).at(1);
      (*(*normals)(2))[lid] += (i->second).at(2);
    }
  }

  // export spring stress
  for (std::map<int, std::vector<double>>::const_iterator i = springstress_.begin();
       i != springstress_.end(); ++i)
  {
    // global id -> local id
    const int lid = stress->Map().LID(i->first);
    // local id on processor
    if (lid >= 0)
    {
      // copy all components of normal vector
      (*(*stress)(0))[lid] += (i->second).at(0);
      (*(*stress)(1))[lid] += (i->second).at(1);
      (*(*stress)(2))[lid] += (i->second).at(2);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |                                                             mhv Dec15|
 *----------------------------------------------------------------------*/
void UTILS::SpringDashpotNew::OutputPrestrOffset(
    Teuchos::RCP<Epetra_Vector>& springprestroffset) const
{
  springprestroffset->Update(1.0, *offset_prestr_new_, 0.0);

  return;
}

/*----------------------------------------------------------------------*
 |                                                             mhv Dec15|
 *----------------------------------------------------------------------*/
void UTILS::SpringDashpotNew::OutputPrestrOffsetOld(
    Teuchos::RCP<Epetra_MultiVector>& springprestroffset) const
{
  // export spring offset length
  for (std::map<int, std::vector<double>>::const_iterator i = offset_prestr_.begin();
       i != offset_prestr_.end(); ++i)
  {
    // global id -> local id
    const int lid = springprestroffset->Map().LID(i->first);
    // local id on processor
    if (lid >= 0)
    {
      // copy all components of spring offset length vector
      (*(*springprestroffset)(0))[lid] = (i->second)[0];
      (*(*springprestroffset)(1))[lid] = (i->second)[1];
      (*(*springprestroffset)(2))[lid] = (i->second)[2];
    }
  }

  return;
}

/*-----------------------------------------------------------------------*
|(private)                                                  pfaller Apr15|
 *-----------------------------------------------------------------------*/
void UTILS::SpringDashpotNew::InitializeCurSurfNormal()
{
  // create MORTAR interface
  mortar_ = Teuchos::rcp(new ADAPTER::CouplingNonLinMortar());

  // create CONTACT elements at interface for normal and gap calculation
  mortar_->SetupSpringDashpot(actdisc_, actdisc_, spring_, coupling_, actdisc_->Comm());

  // create temp vectors for gap initialization
  std::map<int, std::map<int, double>> tmpdgap_;
  std::map<int, std::vector<double>> tmpnormals_;
  std::map<int, std::vector<GEN::pairedvector<int, double>>> tmpdnormals_;

  // empty displacement vector
  Teuchos::RCP<Epetra_Vector> disp;
  disp = LINALG::CreateVector(*(actdisc_->DofRowMap()), true);

  // initialize gap in reference configuration
  mortar_->Interface()->EvaluateDistances(disp, tmpnormals_, tmpdnormals_, gap0_, tmpdgap_);

  return;
}

// ToDo: this function should vanish completely
// obsolete when using new EvaluateRobin function!
/*-----------------------------------------------------------------------*
|(private) adapted from mhv 01/14                           pfaller Apr15|
 *-----------------------------------------------------------------------*/
void UTILS::SpringDashpotNew::GetArea(const std::map<int, Teuchos::RCP<DRT::Element>>& geom)
{
  std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator ele;
  for (ele = geom.begin(); ele != geom.end(); ++ele)
  {
    DRT::Element* element = ele->second.get();

    Teuchos::ParameterList eparams;

    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    element->LocationVector(*(actdisc_), lm, lmowner, lmstride);
    Epetra_SerialDenseMatrix dummat(0, 0);
    Epetra_SerialDenseVector dumvec(0);
    Epetra_SerialDenseVector elevector;
    const int eledim = (int)lm.size();
    elevector.Size(eledim);

    eparams.set("action", "calc_struct_area");
    eparams.set("area", 0.0);
    element->Evaluate(eparams, *(actdisc_), lm, dummat, dummat, dumvec, dumvec, dumvec);

    DRT::Element::DiscretizationType shape = element->Shape();

    double a = eparams.get("area", -1.0);

    // loop over all nodes of the element that share the area
    // do only contribute to my own row nodes
    double apernode = 0.;
    for (int i = 0; i < element->NumNode(); ++i)
    {
      /* here we have to take care to assemble the right stiffness to the nodes!!! (mhv 05/2014):
          we do some sort of "manual" gauss integration here since we have to pay attention to
         assemble the correct stiffness in case of quadratic surface elements*/

      switch (shape)
      {
        case DRT::Element::tri3:
          apernode = a / element->NumNode();
          break;
        case DRT::Element::tri6:
        {
          // integration of shape functions over parameter element surface
          double int_N_cornernode = 0.;
          double int_N_edgemidnode = 1. / 6.;

          int numcornernode = 3;
          int numedgemidnode = 3;

          double weight = numcornernode * int_N_cornernode + numedgemidnode * int_N_edgemidnode;
          double a_inv_weight = a / weight;

          if (i < 3)  // corner nodes
            apernode = int_N_cornernode * a_inv_weight;
          else  // edge mid nodes
            apernode = int_N_edgemidnode * a_inv_weight;
        }
        break;
        case DRT::Element::quad4:
          apernode = a / element->NumNode();
          break;
        case DRT::Element::quad8:
        {
          // integration of shape functions over parameter element surface
          double int_N_cornernode = -1. / 3.;
          double int_N_edgemidnode = 4. / 3.;

          int numcornernode = 4;
          int numedgemidnode = 4;

          double weight = numcornernode * int_N_cornernode + numedgemidnode * int_N_edgemidnode;
          double a_inv_weight = a / weight;

          if (i < 4)  // corner nodes
            apernode = int_N_cornernode * a_inv_weight;
          else  // edge mid nodes
            apernode = int_N_edgemidnode * a_inv_weight;
        }
        break;
        case DRT::Element::quad9:
        {
          // integration of shape functions over parameter element surface
          double int_N_cornernode = 1. / 9.;
          double int_N_edgemidnode = 4. / 9.;
          double int_N_centermidnode = 16. / 9.;

          int numcornernode = 4;
          int numedgemidnode = 4;
          int numcentermidnode = 1;

          double weight = numcornernode * int_N_cornernode + numedgemidnode * int_N_edgemidnode +
                          numcentermidnode * int_N_centermidnode;
          double a_inv_weight = a / weight;

          if (i < 4)  // corner nodes
            apernode = int_N_cornernode * a_inv_weight;
          else if (i == 8)  // center mid node
            apernode = int_N_centermidnode * a / weight;
          else  // edge mid nodes
            apernode = int_N_edgemidnode * a_inv_weight;
        }
        break;
        case DRT::Element::nurbs9:
          dserror(
              "Not yet implemented for Nurbs! To do: Apply the correct weighting of the area per "
              "node!");
          break;
        default:
          dserror("shape type unknown!\n");
          break;
      }

      const int gid = element->Nodes()[i]->Id();
      if (!actdisc_->NodeRowMap()->MyGID(gid)) continue;

      // store area in map (gid, area). erase old value before adding new one
      const double newarea = area_[gid] + apernode;
      area_.erase(gid);
      area_.insert(std::pair<int, double>(gid, newarea));
    }
  }  // for (ele=geom.begin(); ele != geom.end(); ++ele)

  return;
}


/*-----------------------------------------------------------------------*
|(private)                                                    mhv 12/2015|
 *-----------------------------------------------------------------------*/
void UTILS::SpringDashpotNew::InitializePrestrOffset()
{
  offset_prestr_.clear();

  const std::vector<int>& nds = *nodes_;
  for (int j = 0; j < (int)nds.size(); ++j)
  {
    if (actdisc_->NodeRowMap()->MyGID(nds[j]))
    {
      int gid = nds[j];

      DRT::Node* node = actdisc_->gNode(gid);
      if (!node) dserror("Cannot find global node %d", gid);

      int numdof = actdisc_->NumDof(0, node);
      std::vector<int> dofs = actdisc_->Dof(0, node);

      assert(numdof == 3);

      std::vector<double> temp(numdof, 0.0);

      // insert to map
      offset_prestr_.insert(std::pair<int, std::vector<double>>(gid, temp));
    }
  }

  return;
}


/*-----------------------------------------------------------------------*
|(private)                                                  pfaller Apr15|
 *-----------------------------------------------------------------------*/
void UTILS::SpringDashpotNew::GetCurNormals(const Teuchos::RCP<const Epetra_Vector>& disp)
{
  // temp nodal gap
  std::map<int, double> tmpgap;

  // calculate normals and gap using CONTACT elements
  mortar_->Interface()->EvaluateDistances(disp, normals_, dnormals_, tmpgap, dgap_);

  // subtract reference gap from current gap (gap in reference configuration is zero everywhere)
  for (std::map<int, double>::iterator i = tmpgap.begin(); i != tmpgap.end(); ++i)
  {
    std::map<int, double>::iterator j = gap0_.find(i->first);
    if (j == gap0_.end()) gap_[i->first] = i->second;
    //      dserror("The maps of reference gap and current gap are inconsistent.");
    else
      gap_[i->first] = i->second - j->second;
  }

  return;
}

/*-----------------------------------------------------------------------*
|(private)                                                  pfaller Apr15|
 *-----------------------------------------------------------------------*/
void UTILS::SpringDashpotNew::SetSpringType()
{
  // get spring direction from condition
  const std::string* dir = spring_->Get<std::string>("direction");

  if (*dir == "xyz")
    springtype_ = xyz;
  else if (*dir == "refsurfnormal")
    springtype_ = refsurfnormal;
  else if (*dir == "cursurfnormal")
    springtype_ = cursurfnormal;
  else
    dserror(
        "Invalid direction option! Choose DIRECTION xyz, DIRECTION refsurfnormal or DIRECTION "
        "cursurfnormal!");
}
