/*!----------------------------------------------------------------------
\file springdashpot.cpp

\brief Methods for spring and dashpot constraints / boundary conditions:

<pre>
\level 3

\maintainer Martin Pfaller
</pre>

*----------------------------------------------------------------------*/


#include "springdashpot.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include <iostream>
#include "../drt_io/io_pstream.H" // has to go before io.H
#include "../drt_io/io.H"
#include <Epetra_Time.h>

#include "../drt_adapter/adapter_coupling_nonlin_mortar.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_contact/contact_interface.H"
#include "../drt_lib/drt_condition_utils.H"

/*----------------------------------------------------------------------*
 |                                                         pfaller Apr15|
 *----------------------------------------------------------------------*/
UTILS::SpringDashpot::SpringDashpot(
    Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<DRT::Condition> cond):
    actdisc_(dis),
    spring_(cond),
    stiff_tens_(spring_->GetDouble("SPRING_STIFF_TENS")),
    stiff_comp_(spring_->GetDouble("SPRING_STIFF_COMP")),
    offset_(spring_->GetDouble("SPRING_OFFSET")),
    viscosity_(spring_->GetDouble("DASHPOT_VISCOSITY")),
    coupling_(spring_->GetInt("coupling id")),
    nodes_(spring_->Nodes()),
    area_(),
    gap0_(),
    gapnp_(),
    gapn_(),
    gapdt_(),
    dgap_(),
    normals_(),
    dnormals_(),
    offset_prestr_()
{
  // set type of this spring
  SetSpringType();

  if (springtype_ != cursurfnormal && coupling_>=0)
    dserror("Coupling of spring dashpot to reference surface only possible for DIRECTION cursurfnormal.");

  if (springtype_ == cursurfnormal && coupling_==-1)
    dserror("Coupling id necessary for DIRECTION cursurfnormal.");

  // get geometry
  std::map<int,Teuchos::RCP<DRT::Element> >& geom = spring_->Geometry();

  // calculate nodal area
  if (!actdisc_->Comm().MyPID()) IO::cout << "Computing area for spring dashpot condition...\n";
  GetArea(geom);

  // get normal vectors if necessary
  if (springtype_ == cursurfnormal)
    InitializeCurSurfNormal();
  else if (springtype_ == refsurfnormal)
    GetRefNormals(geom);
  // initialize prestressing offset
  InitializePrestrOffset();

  return;
}

/*----------------------------------------------------------------------*
 |                                                             mhv 01/14|
 *----------------------------------------------------------------------*/
void UTILS::SpringDashpot::Evaluate(
    Teuchos::RCP<LINALG::SparseOperator> stiff,
    Teuchos::RCP<Epetra_Vector> fint,
    Teuchos::RCP<Epetra_Vector> disp,
    Teuchos::RCP<Epetra_Vector> velo,
    Teuchos::ParameterList parlist)
{
  if (disp==Teuchos::null) dserror("Cannot find displacement state in discretization");

  // reset last Newton step (common for all spring types)
  gapnp_.clear();
  gapdt_.clear();
  springstress_.clear();

  // get time integrator properties
  const double gamma = parlist.get("scale_gamma",0.0);
  const double beta = parlist.get("scale_beta",1.0);
  const double dt = parlist.get("time_step_size",1.0);
  const double time_scale = gamma/(beta*dt);

  if (springtype_ == cursurfnormal)
  {
    GetCurNormals(disp, parlist);
    stiff->UnComplete(); // sparsity pattern might change
  }

  // loop nodes of current condition
  const std::vector<int>& nds = *nodes_;
  for (int j=0; j<(int)nds.size(); ++j)
  {
    // nodes owned by processor
    if (actdisc_->NodeRowMap()->MyGID(nds[j]))
    {
      int gid = nds[j];
      DRT::Node* node = actdisc_->gNode(gid);
      if (!node) dserror("Cannot find global node %d",gid);

      // get nodal values
      const double nodalarea = area_[gid]; // nodal area
      const std::vector<double> normal = normals_[gid]; // normalized nodal normal
      const std::vector<double> offsetprestr = offset_prestr_[gid]; // get nodal displacement values of last time step for MULF offset

      const int numdof = actdisc_->NumDof(0,node);
      assert (numdof==3);
      std::vector<int> dofs = actdisc_->Dof(0,node);

      // initialize
      double gap = 0.;   // displacement
      double gapdt = 0.; // velocity
      double springstiff = 0.; // spring stiffness

      // assemble into residual vector and stiffness matrix
      std::vector<double> out_vec(numdof, 0.);

      // calculation of normals and displacements differs for each spring variant
      switch (springtype_)
      {
      case all: // spring dashpot acts in every surface dof direction
        // assemble into residual and stiffness matrix
        for (int k=0; k<numdof; ++k)
        {
          if (stiff_tens_ != stiff_comp_)
            dserror("SPRING_STIFF_TENS != SPRING_STIFF_COMP: Different spring moduli for tension and compression not supported "
                "when specifying 'all' as DIRECTION (no ref surface normal information is calculated for that case)! "
                "Only possible for DIRECTION 'refsurfnormal' or 'cursurfnormal'.");

          // extract nodal displacement and velocity
          const double u = (*disp)[disp->Map().LID(dofs[k])]; // displacement
          const double v = (*velo)[velo->Map().LID(dofs[k])]; // velocity

          const double val  = nodalarea*(stiff_tens_*(u-offset_-offsetprestr[k]) + viscosity_*v);
          const double dval = nodalarea*(stiff_tens_ + viscosity_*time_scale);

          const int err = fint->SumIntoGlobalValues(1,&val,&dofs[k]);
          if (err) dserror("SumIntoGlobalValues failed!");
          stiff->Assemble(dval,dofs[k],dofs[k]);
        }
        break;

      case refsurfnormal: // spring dashpot acts in refnormal direction
        for (int k=0; k<numdof; ++k)
        {
          // extract nodal displacement and velocity
          const double u = (*disp)[disp->Map().LID(dofs[k])]; // displacement
          const double v = (*velo)[velo->Map().LID(dofs[k])]; // velocity

          // projection of displacement/velocity onto nodal reference normal
          gap -= (u-offsetprestr[k])*normal[k]; // gap = (u \cdot N)
          gapdt -= v*normal[k]; // gapdt = (v \cdot N)

          // save for output
          gapnp_[gid] = gap;
        }

        // select spring stiffness
        springstiff = SelectStiffness(gap);

        // assemble into residual vector and stiffness matrix
        for (int k=0; k<numdof; ++k)
        {
          // force
          const double val = - nodalarea*(springstiff*(gap-offset_) + viscosity_*gapdt)*normal[k];
          const int err = fint->SumIntoGlobalValues(1,&val,&dofs[k]);
          if (err) dserror("SumIntoGlobalValues failed!");

          // stiffness
          for (int m=0; m<numdof; ++m)
          {
            const double dval = nodalarea*(springstiff + viscosity_*time_scale)*normal[k]*normal[m];
            stiff->Assemble(dval,dofs[k],dofs[m]);
          }

          // store negative value of internal force for output (=reaction force)
          out_vec[k] = - val;
        }
        // add to output
        springstress_.insert(std::pair<int, std::vector<double> >(gid, out_vec));
        break;

      case cursurfnormal: // spring dashpot acts in curnormal direction
        // spring displacement
        gap = gapnp_[gid];
        gapdt = gapdt_[gid];

        // select spring stiffnes
        springstiff = SelectStiffness(gap);

        // assemble into residual vector and stiffness matrix
        for (int k=0; k<numdof; ++k)
        {
          // force
          const double val = - nodalarea*(springstiff*(gap-offsetprestr[k]-offset_) + viscosity_*gapdt)*normal[k];
          const int err = fint->SumIntoGlobalValues(1,&val,&dofs[k]);
          if (err) dserror("SumIntoGlobalValues failed!");

          // stiffness
          std::map<int,double> dgap = dgap_[gid];
          std::vector<GEN::pairedvector<int,double> > dnormal = dnormals_[gid];

          // check if projection exists
          if(!dnormal.empty() && !dgap.empty())
          {
            // linearize gap
            for(std::map<int, double>::iterator i = dgap.begin(); i != dgap.end(); ++i)
            {
              const double dval = -nodalarea*(springstiff*(i->second) + viscosity_*(i->second)/dt)*normal[k];
              stiff->Assemble(dval,dofs[k],i->first);
            }

            // linearize normal
            for(GEN::pairedvector<int,double>::iterator i = dnormal[k].begin(); i != dnormal[k].end(); ++i)
            {
              const double dval = -nodalarea*(springstiff*(gap-offset_) + viscosity_*gapdt)*(i->second);
              stiff->Assemble(dval,dofs[k],i->first);
            }
          }
          else
          {
//            dserror("Projection does not exist for node %d.", gid+1);
          }
          // store negative value of internal force for output (=reaction force)
          out_vec[k] = - val;
        }
        // add to output
        springstress_.insert(std::pair<int, std::vector<double> >(gid, out_vec));
        break;
      }
    } //node owned by processor
  } //loop over nodes

  if (springtype_ == cursurfnormal)
    stiff->Complete(); // sparsity pattern might have changed

  return;
}

/*----------------------------------------------------------------------*
 |                                                             mhv 12/15|
 *----------------------------------------------------------------------*/
void UTILS::SpringDashpot::Reset(
    Teuchos::RCP<Epetra_Vector> dis)
{

  // loop nodes of current condition
  const std::vector<int>& nds = *nodes_;
  for (int j=0; j<(int)nds.size(); ++j)
  {
    // nodes owned by processor
    if (actdisc_->NodeRowMap()->MyGID(nds[j]))
    {
      int gid = nds[j];
      DRT::Node* node = actdisc_->gNode(gid);
      if (!node) dserror("Cannot find global node %d",gid);

      const int numdof = actdisc_->NumDof(0,node);
      assert (numdof==3);
      std::vector<int> dofs = actdisc_->Dof(0,node);

      // initialize. calculation of displacements differs for each spring variant
      std::vector<double> uoff(numdof, 0.0); // displacement vector of condition nodes
      if (springtype_ == refsurfnormal ||
          springtype_ == all)
      {
        for (int k=0; k<numdof; ++k)
        {
          uoff[k] = -(*dis)[dis->Map().LID(dofs[k])]; // extract nodal displacement to be offset
        }

        std::vector<double> offpr = offset_prestr_[gid];
        offset_prestr_.erase(gid);

        for (int k=0; k<numdof; ++k)
        {
          uoff[k] += offpr[k]; // accumulate displacements to be offset
        }

      }
      offset_prestr_.insert(std::pair<int, std::vector<double> >(gid, uoff));

    } //node owned by processor
  } //loop over nodes

}

/*----------------------------------------------------------------------*
 |                                                             mhv 12/15|
 *----------------------------------------------------------------------*/
void UTILS::SpringDashpot::SetRestart(
    Teuchos::RCP<Epetra_MultiVector> vec)
{

  // loop nodes of current condition
  const std::vector<int>& nds = *nodes_;
  for (int j=0; j<(int)nds.size(); ++j)
  {
    // nodes owned by processor
    if (actdisc_->NodeRowMap()->MyGID(nds[j]))
    {
      int gid = nds[j];
      DRT::Node* node = actdisc_->gNode(gid);
      if (!node) dserror("Cannot find global node %d",gid);

      const int numdof = actdisc_->NumDof(0,node);
      assert (numdof==3);
      std::vector<int> dofs = actdisc_->Dof(0,node);


      if (springtype_ == refsurfnormal ||
          springtype_ == all)
      {

        // import spring offset length
        for(std::map<int, std::vector<double> >::iterator i = offset_prestr_.begin(); i != offset_prestr_.end(); ++i)
        {
          // global id -> local id
          const int lid = vec->Map().LID(i->first);
          // local id on processor
          if (lid>=0)
          {
            // copy all components of spring offset length vector
            (i->second)[0] = (*(*vec)(0))[lid];
            (i->second)[1] = (*(*vec)(1))[lid];
            (i->second)[2] = (*(*vec)(2))[lid];

          }
        }

      }

    } //node owned by processor
  } //loop over nodes

}

/*----------------------------------------------------------------------*
 |                                                         pfaller Jan14|
 *----------------------------------------------------------------------*/
void UTILS::SpringDashpot::OutputGapNormal(
    Teuchos::RCP<Epetra_Vector> &gap,
    Teuchos::RCP<Epetra_MultiVector> &normals,
    Teuchos::RCP<Epetra_MultiVector> &stress)
{
  // export gap function
  for(std::map<int, double>::iterator i = gapnp_.begin(); i != gapnp_.end(); ++i)
  {
    // global id -> local id
    const int lid = gap->Map().LID(i->first);
    // local id on processor
    if (lid>=0)
      (*gap)[lid] += i->second;
  }

  // export normal
  for(std::map<int, std::vector<double> >::iterator i = normals_.begin(); i != normals_.end(); ++i)
  {
    // global id -> local id
    const int lid = normals->Map().LID(i->first);
    // local id on processor
    if (lid>=0)
    {
      // copy all components of normal vector
      (*(*normals)(0))[lid] += (i->second)[0];
      (*(*normals)(1))[lid] += (i->second)[1];
      (*(*normals)(2))[lid] += (i->second)[2];
    }
  }

  // export spring stress
  for(std::map<int, std::vector<double> >::iterator i = springstress_.begin(); i != springstress_.end(); ++i)
  {
    // global id -> local id
    const int lid = stress->Map().LID(i->first);
    // local id on processor
    if (lid>=0)
    {
      // copy all components of normal vector
      (*(*stress)(0))[lid] += (i->second)[0];
      (*(*stress)(1))[lid] += (i->second)[1];
      (*(*stress)(2))[lid] += (i->second)[2];
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |                                                             mhv Dec15|
 *----------------------------------------------------------------------*/
void UTILS::SpringDashpot::OutputPrestrOffset(
    Teuchos::RCP<Epetra_MultiVector> &springprestroffset)
{

  // export spring offset length
  for(std::map<int, std::vector<double> >::iterator i = offset_prestr_.begin(); i != offset_prestr_.end(); ++i)
  {
    // global id -> local id
    const int lid = springprestroffset->Map().LID(i->first);
    // local id on processor
    if (lid>=0)
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
|                                                           pfaller May16|
 *-----------------------------------------------------------------------*/
void UTILS::SpringDashpot::Update()
{
  // store current time step
  gapn_ = gapnp_;
}

/*-----------------------------------------------------------------------*
|(private)                                                  pfaller Apr15|
 *-----------------------------------------------------------------------*/
void UTILS::SpringDashpot::InitializeCurSurfNormal()
{
  // create MORTAR interface
  mortar_ = Teuchos::rcp(new ADAPTER::CouplingNonLinMortar());

  // create CONTACT elements at interface for normal and gap calculation
  mortar_->SetupSpringDashpot(actdisc_, actdisc_, spring_, coupling_, actdisc_->Comm());

  // create temp vectors for gap initialization
  std::map<int,std::map<int,double> > tmpdgap_;
  std::map<int, std::vector<double> > tmpnormals_;
  std::map<int,std::vector<GEN::pairedvector<int,double> > > tmpdnormals_;

  // empty displacement vector
  Teuchos::RCP<Epetra_Vector> disp;
  disp = LINALG::CreateVector(*(actdisc_->DofRowMap()), true);

  // initialize gap in reference configuration
  mortar_->Interface()->EvaluateDistances(disp,tmpnormals_,tmpdnormals_,gap0_,tmpdgap_);

  return;
}

/*-----------------------------------------------------------------------*
|(private) adapted from mhv 01/14                           pfaller Apr15|
 *-----------------------------------------------------------------------*/
void UTILS::SpringDashpot::GetArea(const std::map<int,Teuchos::RCP<DRT::Element> >& geom)
{
  std::map<int,Teuchos::RCP<DRT::Element> >::const_iterator ele;
  for (ele=geom.begin(); ele != geom.end(); ++ele)
  {
    DRT::Element* element = ele->second.get();

    Teuchos::ParameterList eparams;

    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    element->LocationVector(*(actdisc_),lm,lmowner,lmstride);
    Epetra_SerialDenseMatrix dummat(0,0);
    Epetra_SerialDenseVector dumvec(0);
    Epetra_SerialDenseVector elevector;
    const int eledim = (int)lm.size();
    elevector.Size(eledim);

    eparams.set("action","calc_struct_area");
    eparams.set("area",0.0);
    element->Evaluate(eparams,*(actdisc_),lm,dummat,dummat,dumvec,dumvec,dumvec);

    DRT::Element::DiscretizationType shape = element->Shape();

    double a = eparams.get("area",-1.0);

    // loop over all nodes of the element that share the area
    // do only contribute to my own row nodes
    double apernode = 0.;
    for (int i=0; i<element->NumNode(); ++i)
    {
      /* here we have to take care to assemble the right stiffness to the nodes!!! (mhv 05/2014):
          we do some sort of "manual" gauss integration here since we have to pay attention to assemble
          the correct stiffness in case of quadratic surface elements*/

      switch(shape)
      {
      case DRT::Element::tri3:
        apernode = a / element->NumNode();
        break;
      case DRT::Element::tri6:
      {
        //integration of shape functions over parameter element surface
        double int_N_cornernode = 0.;
        double int_N_edgemidnode = 1./6.;

        int numcornernode = 3;
        int numedgemidnode = 3;

        double weight = numcornernode * int_N_cornernode + numedgemidnode * int_N_edgemidnode;

        //corner nodes
        if (i==0) apernode = int_N_cornernode * a / weight;
        if (i==1) apernode = int_N_cornernode * a / weight;
        if (i==2) apernode = int_N_cornernode * a / weight;
        //edge mid nodes
        if (i==3) apernode = int_N_edgemidnode * a / weight;
        if (i==4) apernode = int_N_edgemidnode * a / weight;
        if (i==5) apernode = int_N_edgemidnode * a / weight;
      }
      break;
      case DRT::Element::quad4:
        apernode = a / element->NumNode();
        break;
      case DRT::Element::quad8:
      {
        //integration of shape functions over parameter element surface
        double int_N_cornernode = -1./3.;
        double int_N_edgemidnode = 4./3.;

        int numcornernode = 4;
        int numedgemidnode = 4;

        double weight = numcornernode * int_N_cornernode + numedgemidnode * int_N_edgemidnode;

        //corner nodes
        if (i==0) apernode = int_N_cornernode * a / weight;
        if (i==1) apernode = int_N_cornernode * a / weight;
        if (i==2) apernode = int_N_cornernode * a / weight;
        if (i==3) apernode = int_N_cornernode * a / weight;
        //edge mid nodes
        if (i==4) apernode = int_N_edgemidnode * a / weight;
        if (i==5) apernode = int_N_edgemidnode * a / weight;
        if (i==6) apernode = int_N_edgemidnode * a / weight;
        if (i==7) apernode = int_N_edgemidnode * a / weight;
      }
      break;
      case DRT::Element::quad9:
      {
        //integration of shape functions over parameter element surface
        double int_N_cornernode = 1./9.;
        double int_N_edgemidnode = 4./9.;
        double int_N_centermidnode = 16./9.;

        int numcornernode = 4;
        int numedgemidnode = 4;
        int numcentermidnode = 1;

        double weight = numcornernode * int_N_cornernode + numedgemidnode * int_N_edgemidnode + numcentermidnode * int_N_centermidnode;

        //corner nodes
        if (i==0) apernode = int_N_cornernode * a / weight;
        if (i==1) apernode = int_N_cornernode * a / weight;
        if (i==2) apernode = int_N_cornernode * a / weight;
        if (i==3) apernode = int_N_cornernode * a / weight;
        //edge mid nodes
        if (i==4) apernode = int_N_edgemidnode * a / weight;
        if (i==5) apernode = int_N_edgemidnode * a / weight;
        if (i==6) apernode = int_N_edgemidnode * a / weight;
        if (i==7) apernode = int_N_edgemidnode * a / weight;
        //center mid node
        if (i==8) apernode = int_N_centermidnode * a / weight;
      }
      break;
      case DRT::Element::nurbs9:
        dserror("Not yet implemented for Nurbs! To do: Apply the correct weighting of the area per node!");
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
  } // for (ele=geom.begin(); ele != geom.end(); ++ele)

  return;
}

/*-----------------------------------------------------------------------*
|(private) adapted from mhv 01/14                           pfaller Apr15|
 *-----------------------------------------------------------------------*/
void UTILS::SpringDashpot::GetRefNormals(const std::map<int,Teuchos::RCP<DRT::Element> >& geom)
{
  // a vector for all row dofs to hold normals interpolated to the nodes
  Epetra_Vector refnodalnormals(*actdisc_->DofRowMap(),true);

  std::map<int,Teuchos::RCP<DRT::Element> >::const_iterator ele;
  for (ele=geom.begin(); ele != geom.end(); ++ele)
  {
    DRT::Element* element = ele->second.get();

    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    element->LocationVector(*(actdisc_),lm,lmowner,lmstride);
    Epetra_SerialDenseMatrix dummat(0,0);
    Epetra_SerialDenseVector dumvec(0);
    Epetra_SerialDenseVector elevector;
    const int eledim = (int)lm.size();
    elevector.Size(eledim);

    Teuchos::ParameterList eparams2;
    eparams2.set("action","calc_ref_nodal_normals");
    element->Evaluate(eparams2,*(actdisc_),lm,dummat,dummat,elevector,dumvec,dumvec);
    LINALG::Assemble(refnodalnormals,elevector,lm,lmowner);
  }

  const std::vector<int>& nds = *nodes_;
  for (int j=0; j<(int)nds.size(); ++j)
  {
    if (actdisc_->NodeRowMap()->MyGID(nds[j]))
    {
      int gid = nds[j];

      DRT::Node* node = actdisc_->gNode(gid);
      if (!node) dserror("Cannot find global node %d",gid);

      int numdof = actdisc_->NumDof(0,node);
      std::vector<int> dofs = actdisc_->Dof(0,node);

      assert (numdof==3);

      std::vector<double> temp(numdof, 0.0);
      double temp_norm_sq = 0.;
      for (int k=0; k<numdof; ++k)
      {
        temp[k] = refnodalnormals[refnodalnormals.Map().LID(dofs[k])];
        temp_norm_sq += temp[k]*temp[k];
      }

      // calculate vector length
      const double temp_norm = sqrt(temp_norm_sq);

      if (temp_norm == 0.)
        dserror("Nodal normal has length 0.");

      // normalize vector
      for (int k=0; k<numdof; ++k)
        temp[k] /= temp_norm;

      // insert to map
      normals_.insert(std::pair<int, std::vector<double> >(gid, temp));
    }
  }

  return;
}

/*-----------------------------------------------------------------------*
|(private)                                                    mhv 12/2015|
 *-----------------------------------------------------------------------*/
void UTILS::SpringDashpot::InitializePrestrOffset()
{
  offset_prestr_.clear();

  const std::vector<int>& nds = *nodes_;
  for (int j=0; j<(int)nds.size(); ++j)
  {
    if (actdisc_->NodeRowMap()->MyGID(nds[j]))
    {
      int gid = nds[j];

      DRT::Node* node = actdisc_->gNode(gid);
      if (!node) dserror("Cannot find global node %d",gid);

      int numdof = actdisc_->NumDof(0,node);
      std::vector<int> dofs = actdisc_->Dof(0,node);

      assert (numdof==3);

      std::vector<double> temp(numdof, 0.0);

      // insert to map
      offset_prestr_.insert(std::pair<int, std::vector<double> >(gid, temp));

    }
  }

  return;
}


/*-----------------------------------------------------------------------*
|(private)                                                  pfaller Apr15|
 *-----------------------------------------------------------------------*/
void UTILS::SpringDashpot::GetCurNormals(Teuchos::RCP<Epetra_Vector> disp,
                                         Teuchos::ParameterList parlist)
{
  // get current time step size
  const double dt = parlist.get("time_step_size",1.0);

  // reset last newton step (only for curnormal)
  dgap_.clear();
  normals_.clear();
  dnormals_.clear();

  // temp nodal gap
  std::map<int, double> tmpgap;

  // calculate normals and gap using CONTACT elements
  mortar_->Interface()->EvaluateDistances(disp,normals_,dnormals_,tmpgap,dgap_);

  for(std::map<int, double>::iterator i = tmpgap.begin(); i != tmpgap.end(); ++i)
  {
    // subtract reference gap from current gap (gap in reference configuration is zero everywhere)
    std::map<int, double>::iterator j = gap0_.find(i->first);
    if(j == gap0_.end())
      gapnp_[i->first] = i->second;
//      dserror("The maps of reference gap and current gap are inconsistent.");
    else
      gapnp_[i->first] = i->second - j->second;

    // calculate gap velocity via local finite difference (not the best way but also not the worst)
    gapdt_[i->first] = (gapnp_[i->first] - gapn_[i->first])/dt;
  }

  return;
}

/*-----------------------------------------------------------------------*
|(private)                                                  pfaller Apr15|
 *-----------------------------------------------------------------------*/
void UTILS::SpringDashpot::SetSpringType()
{
  // get spring direction from condition
  const std::string* dir = spring_->Get<std::string>("direction");

  if (*dir=="all")
    springtype_ = all;
  else if (*dir=="refsurfnormal")
    springtype_ = refsurfnormal;
  else if (*dir=="cursurfnormal")
    springtype_ = cursurfnormal;
  else
    dserror("Invalid direction option! Choose DIRECTION all, DIRECTION refsurfnormal or DIRECTION cursurfnormal!");
}
