/*!----------------------------------------------------------------------
\file drt_locsys.cpp

\brief Class controlling local coordinate systems on points, lines and
surfaces and supplying all necessary transformation methods for parallel
vectors and matrices.

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*----------------------------------------------------------------------*/

#include "drt_locsys.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_function.H"
#include "../drt_fem_general/largerotations.H"



/*-------------------------------------------------------------------*
 |  ctor (public)                                         meier 06/13|
 *-------------------------------------------------------------------*/
DRT::UTILS::LocsysManager::LocsysManager(DRT::Discretization& discret):
discret_(discret),
dim_(-1),
numlocsys_(-1),
locsyscurvefunct_(false)
{
  // get problem dimension (2D or 3D) and store into dim_
  dim_ = DRT::Problem::Instance()->NDim();

  if (Dim()!= 2 && Dim()!=3) dserror("ERROR: Locsys problem must be 2D or 3D");

  // get node row layout of discretization
  const Epetra_Map* noderowmap = discret_.NodeRowMap();

  // create locsys vector and initialize to -1
  locsystoggle_ = LINALG::CreateVector(*noderowmap,false);
  locsystoggle_->PutScalar(-1.0);

  // check for locsys boundary conditions
  Discret().GetCondition("Locsys",locsysconds_);
  numlocsys_ = (int)locsysconds_.size();
  id_.resize(numlocsys_);
  typelocsys_.resize(numlocsys_);

  for (int i=0; i<NumLocsys(); ++i)
  {
    id_[i] = locsysconds_[i]->Id();
  }

  //First Setup is made in the constructor. If we have no time dependent locsys conditions
  //in our problem, this is the only time where the whole setup routine is conducted.
  Setup(-1.0);
}

/*-------------------------------------------------------------------*
 |  set-up                                                meier 06/13|
 *-------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::Setup(const double time)
{
  // IMPORTANT NOTE:
  // The definition of local coordinate systems only makes sense in
  // combination with Dirichlet boundary conditions. This means that
  // in order to define a boundary condition, both locsys
  // AND Dirichlet condition have to formulated for the same entity
  // (i.e. point, line, surface, volume).

  // LIMITATIONS:
  // - So far locsys only works for 2D and 3D solids and for beam elements of Kirchhoff type
  // - Due to this limitation it's necessary to distinguish between this different element types
  //   by means of there nodal DoFs. If further element types are integrated into locsys
  //   more elaborate criteria might be useful.

  //If we have no time curves or functions in the locsys conditions the whole Setup method is only
  //conducted once in the constructor (where time is set to -1.0).
  if (time>=0.0 and locsyscurvefunct_==false)
    return;

  // get dof row map of discretization
  const Epetra_Map* dofrowmap = discret_.DofRowMap();

  // get node row layout of discretization
  const Epetra_Map* noderowmap = discret_.NodeRowMap();

  // Since also time dependent conditions are possible we clear all local systems in the beginning
  nodalrotvectors_.clear();

  // As for Dirichlet conditions, we keep to a very strict hierarchy
  // for evaluation of the Locsys conditions: Volume locsys conditions
  // are evaluated first, followed by Surface and Line locsys conditions
  // and finally Point locsys conditions. This means that nodes carrying
  // different types of locsys conditions are dominated by the rule "Point
  // above Line above Surface above Volume". When two locsys conditions of
  // the same type are defined for one node, ordering in the input file matters!

  //**********************************************************************
  // read volume locsys conditions
  //**************************+*******************************************
  for (int i=0; i<NumLocsys(); ++i)
  {
    DRT::Condition* currlocsys = locsysconds_[i];

    if (currlocsys->Type() == DRT::Condition::VolumeLocsys)
    {
      typelocsys_[i] = DRT::Condition::VolumeLocsys;

      const std::vector<double>* rotangle = currlocsys->Get<std::vector<double> >("rotangle");
      const std::vector<int>*    curve  = currlocsys->Get<std::vector<int> >("curve");
      const std::vector<int>*    funct  = currlocsys->Get<std::vector<int> >("funct");
      const std::vector<int>*    useUpdatedNodePos  = currlocsys->Get<std::vector<int> >("useupdatednodepos");
      const std::vector<int>*    nodes = currlocsys->Nodes();

      //Check, if we have time dependent locsys conditions (through curves or functions)
      if (((*curve)[0]>0 or (*curve)[1]>0 or (*curve)[2]>0) or
          ((*funct)[0]>0 or (*funct)[1]>0 or (*funct)[2]>0))
        locsyscurvefunct_=true;

      //Here we have the convention that 2D problems "live" in the global xy-plane.
      if (Dim()==2 and ((*rotangle)[0]!=0 or (*rotangle)[1]!=0))
        dserror("For 2D problems (xy-plane) the vector ROTANGLE has to be parallel to the global z-axis!");

      // Check, if the updated node positions shall be used for evaluation of the functions 'funct'
      RCP<const Epetra_Vector> dispnp;
      if (((*useUpdatedNodePos)[0] == 1) && (time >= 0.0)){
        dispnp = Discret().GetState("dispnp");
        if (dispnp == Teuchos::null)
          dserror("Locsys: Cannot find state 'dispnp'! You need to set the state 'dispnp' before calling the locsys setup.");
      }

      //Each component j of the pseudo rotation vector that rotates the global xyz system onto the local system
      //assigned to each node consists of a constant, a time dependent and spatially variable part:
      //currotangle_j(x,t) = rotangle_j * curve_j(t) * funct_j(x)
      LINALG::Matrix<3,1> currotangle;
      currotangle.Clear();

      for (int k=0;k<(int)nodes->size();++k)
      {
        bool havenode = Discret().HaveGlobalNode((*nodes)[k]);
        if (!havenode) continue;

        //Weights of rotations vector due to time curve and spatial function
        for (int j=0;j<3;j++)
        {
          // factor given by time curve
          std::vector<double> curvefac(1,1.0);
          int curvenum = (*curve)[j];
          if (curvenum>=0 && time>=0.0)
          {
            curvefac = DRT::Problem::Instance()->Curve(curvenum).FctDer(time,0);
          }
          // factor given by spatial function
          double functfac = 1.0;
          if ((*funct)[j]>0)
          {
            DRT::Node* node = Discret().gNode((*nodes)[k]);

            // Determine node position, which shall be used for evaluating the function, and evaluate it
            if (((*useUpdatedNodePos)[0] == 1) && (time >= 0.0)){
              // Obtain current displacement for node
              std::vector<int> lm;
              Discret().Dof(node,lm);

              std::vector<double> currDisp;
              currDisp.resize(lm.size());

              DRT::UTILS::ExtractMyValues(*dispnp,currDisp,lm);

              // Calculate current position for node
              double currPos[Dim()];
              const double* xp = node->X();

              for (int i=0; i<Dim(); ++i) {
                currPos[i] = xp[i] + currDisp[i];
              }

              // Evaluate function with current node position
              functfac=(DRT::Problem::Instance()->Funct((*funct)[j]-1)).Evaluate(j, &currPos[0], 0.0, &discret_);
            } else {
              // Evaluate function with reference node position
              functfac=(DRT::Problem::Instance()->Funct((*funct)[j]-1)).Evaluate(j, node->X(), 0.0, &discret_);
            }
          }
          currotangle(j)=(*rotangle)[j]*curvefac[0]*functfac;
        }

        nodalrotvectors_[(*nodes)[k]]=currotangle;

        int indices = (*nodes)[k];
        double values  = i;
        locsystoggle_->ReplaceGlobalValues(1,&values,&indices);
      }
    }
    else if (currlocsys->Type() == DRT::Condition::SurfaceLocsys ||
             currlocsys->Type() == DRT::Condition::LineLocsys ||
             currlocsys->Type() == DRT::Condition::PointLocsys)
    {
      // already done
    }
    else
      dserror("ERROR: Unknown type of locsys condition!");
  }  // end volume locsys
  //**********************************************************************
  // read surface locsys conditions
  //**************************+*******************************************
  for (int i=0; i<NumLocsys(); ++i)
  {
    DRT::Condition* currlocsys = locsysconds_[i];

    if (currlocsys->Type() == DRT::Condition::SurfaceLocsys)
    {
      typelocsys_[i] = DRT::Condition::SurfaceLocsys;

      const std::vector<double>* rotangle = currlocsys->Get<std::vector<double> >("rotangle");
      const std::vector<int>*    curve  = currlocsys->Get<std::vector<int> >("curve");
      const std::vector<int>*    funct  = currlocsys->Get<std::vector<int> >("funct");
      const std::vector<int>*    useUpdatedNodePos  = currlocsys->Get<std::vector<int> >("useupdatednodepos");
      const std::vector<int>*    nodes = currlocsys->Nodes();

      //Check, if we have time dependent locsys conditions (through curves or functions)
      if (((*curve)[0]>0 or (*curve)[1]>0 or (*curve)[2]>0) or
          ((*funct)[0]>0 or (*funct)[1]>0 or (*funct)[2]>0))
        locsyscurvefunct_=true;

      //Here we have the convention that 2D problems "live" in the global xy-plane.
      if (Dim()==2 and ((*rotangle)[0]!=0 or (*rotangle)[1]!=0))
        dserror("For 2D problems (xy-plane) the vector ROTANGLE has to be parallel to the global z-axis!");

      // Check, if the updated node positions shall be used for evaluation of the functions 'funct'
      RCP<const Epetra_Vector> dispnp;
      if (((*useUpdatedNodePos)[0] == 1) && (time >= 0.0)){
        dispnp = Discret().GetState("dispnp");
        if (dispnp == Teuchos::null)
          dserror("Locsys: Cannot find state 'dispnp'! You need to set the state 'dispnp' before calling the locsys setup.");
      }

      //Each component j of the pseudo rotation vector that rotates the global xyz system onto the local system
      //assigned to each node consists of a constant, a time dependent and spatially variable part:
      //currotangle_j(x,t) = rotangle_j * curve_j(t) * funct_j(x)
      LINALG::Matrix<3,1> currotangle;
      currotangle.Clear();

      for (int k=0;k<(int)nodes->size();++k)
      {
        bool havenode = Discret().HaveGlobalNode((*nodes)[k]);
        if (!havenode) continue;

        //Weights of rotations vector due to time curve and spatial function
        for (int j=0;j<3;j++)
        {
          // factor given by time curve
          std::vector<double> curvefac(1,1.0);
          int curvenum = (*curve)[j];
          if (curvenum>=0 && time>=0.0)
          {
            curvefac = DRT::Problem::Instance()->Curve(curvenum).FctDer(time,0);
          }
          // factor given by spatial function
          double functfac = 1.0;
          if ((*funct)[j]>0)
          {
            DRT::Node* node = Discret().gNode((*nodes)[k]);

            // Determine node position, which shall be used for evaluating the function, and evaluate it
            if (((*useUpdatedNodePos)[0] == 1) && (time >= 0.0)){
              // Obtain current displacement for node
              std::vector<int> lm;
              Discret().Dof(node,lm);

              std::vector<double> currDisp;
              currDisp.resize(lm.size());

              DRT::UTILS::ExtractMyValues(*dispnp,currDisp,lm);

              // Calculate current position for node
              double currPos[Dim()];
              const double* xp = node->X();

              for (int i=0; i<Dim(); ++i) {
                currPos[i] = xp[i] + currDisp[i];
              }

              // Evaluate function with current node position
              functfac=(DRT::Problem::Instance()->Funct((*funct)[j]-1)).Evaluate(j, &currPos[0], 0.0, &discret_);
            } else {
              // Evaluate function with reference node position
              functfac=(DRT::Problem::Instance()->Funct((*funct)[j]-1)).Evaluate(j, node->X(), 0.0, &discret_);
            }
          }
          currotangle(j)=(*rotangle)[j]*curvefac[0]*functfac;
        }

        nodalrotvectors_[(*nodes)[k]]=currotangle;

        int indices = (*nodes)[k];
        double values  = i;
        locsystoggle_->ReplaceGlobalValues(1,&values,&indices);
      }
    }
    else if (currlocsys->Type() == DRT::Condition::VolumeLocsys ||
             currlocsys->Type() == DRT::Condition::LineLocsys ||
             currlocsys->Type() == DRT::Condition::PointLocsys)
    {
      // already done
    }
    else
      dserror("ERROR: Unknown type of locsys condition!");
  }  // end surface locsys

  //**********************************************************************
  // read line locsys conditions
  //**********************************************************************
  for (int i=0; i<NumLocsys(); ++i)
  {
    DRT::Condition* currlocsys = locsysconds_[i];

    if (currlocsys->Type() == DRT::Condition::LineLocsys)
    {
      typelocsys_[i] = DRT::Condition::LineLocsys;

      const std::vector<double>* rotangle = currlocsys->Get<std::vector<double> >("rotangle");
      const std::vector<int>*    curve  = currlocsys->Get<std::vector<int> >("curve");
      const std::vector<int>*    funct  = currlocsys->Get<std::vector<int> >("funct");
      const std::vector<int>*    useUpdatedNodePos  = currlocsys->Get<std::vector<int> >("useupdatednodepos");
      const std::vector<int>*    nodes = currlocsys->Nodes();

      //Check, if we have time dependent locsys conditions (through curves or functions)
      if (((*curve)[0]>0 or (*curve)[1]>0 or (*curve)[2]>0) or
          ((*funct)[0]>0 or (*funct)[1]>0 or (*funct)[2]>0))
        locsyscurvefunct_=true;

      //Here we have the convention that 2D problems "live" in the global xy-plane.
      if (Dim()==2 and ((*rotangle)[0]!=0 or (*rotangle)[1]!=0))
        dserror("For 2D problems (xy-plane) the vector ROTANGLE has to be parallel to the global z-axis!");

      // Check, if the updated node positions shall be used for evaluation of the functions 'funct'
      RCP<const Epetra_Vector> dispnp;
      if (((*useUpdatedNodePos)[0] == 1) && (time >= 0.0)){
        dispnp = Discret().GetState("dispnp");
        if (dispnp == Teuchos::null)
          dserror("Locsys: Cannot find state 'dispnp'! You need to set the state 'dispnp' before calling the locsys setup.");
      }

      //Each component j of the pseudo rotation vector that rotates the global xyz system onto the local system
      //assigned to each node consists of a constant, a time dependent and spatially variable part:
      //currotangle_j(x,t) = rotangle_j * curve_j(t) * funct_j(x)
      LINALG::Matrix<3,1> currotangle;
      currotangle.Clear();

      for (int k=0;k<(int)nodes->size();++k)
      {
        bool havenode = Discret().HaveGlobalNode((*nodes)[k]);
        if (!havenode) continue;

        //Weights of rotations vector due to time curve and spatial function
        for (int j=0;j<3;j++)
        {
          // factor given by time curve
          std::vector<double> curvefac(1,1.0);
          int curvenum = (*curve)[j];
          if (curvenum>=0 && time>=0.0)
          {
            curvefac = DRT::Problem::Instance()->Curve(curvenum).FctDer(time,0);
          }
          // factor given by spatial function
          double functfac = 1.0;
          if ((*funct)[j]>0)
          {
            DRT::Node* node = Discret().gNode((*nodes)[k]);

            // Determine node position, which shall be used for evaluating the function, and evaluate it
            if (((*useUpdatedNodePos)[0] == 1) && (time >= 0.0)){
              // Obtain current displacement for node
              std::vector<int> lm;
              Discret().Dof(node,lm);

              std::vector<double> currDisp;
              currDisp.resize(lm.size());

              DRT::UTILS::ExtractMyValues(*dispnp,currDisp,lm);

              // Calculate current position for node
              double currPos[Dim()];
              const double* xp = node->X();

              for (int i=0; i<Dim(); ++i) {
                currPos[i] = xp[i] + currDisp[i];
              }

              // Evaluate function with current node position
              functfac=(DRT::Problem::Instance()->Funct((*funct)[j]-1)).Evaluate(j, &currPos[0], 0.0, &discret_);
            } else {
              // Evaluate function with reference node position
              functfac=(DRT::Problem::Instance()->Funct((*funct)[j]-1)).Evaluate(j, node->X(), 0.0, &discret_);
            }
          }
          currotangle(j)=(*rotangle)[j]*curvefac[0]*functfac;
        }

        nodalrotvectors_[(*nodes)[k]]=currotangle;

        int indices = (*nodes)[k];
        double values  = i;
        locsystoggle_->ReplaceGlobalValues(1,&values,&indices);
      }
    }
    else if (currlocsys->Type() == DRT::Condition::VolumeLocsys ||
             currlocsys->Type() == DRT::Condition::SurfaceLocsys ||
             currlocsys->Type() == DRT::Condition::PointLocsys)
    {
      // already done
    }
    else
      dserror("ERROR: Unknown type of locsys condition!");
  }  // end line locsys

  //**********************************************************************
  // read point locsys conditions
  //**********************************************************************
  for (int i=0; i<NumLocsys(); ++i)
  {
    DRT::Condition* currlocsys = locsysconds_[i];

    if (currlocsys->Type() == DRT::Condition::PointLocsys)
    {
      typelocsys_[i] = DRT::Condition::PointLocsys;

      const std::vector<double>* rotangle = currlocsys->Get<std::vector<double> >("rotangle");
      const std::vector<int>*    curve  = currlocsys->Get<std::vector<int> >("curve");
      const std::vector<int>*    funct  = currlocsys->Get<std::vector<int> >("funct");
      const std::vector<int>*    useUpdatedNodePos  = currlocsys->Get<std::vector<int> >("useupdatednodepos");
      const std::vector<int>*    nodes = currlocsys->Nodes();

      //Check, if we have time dependent locsys conditions (through curves or functions)
      if (((*curve)[0]>0 or (*curve)[1]>0 or (*curve)[2]>0) or
          ((*funct)[0]>0 or (*funct)[1]>0 or (*funct)[2]>0))
        locsyscurvefunct_=true;

      //Here we have the convention that 2D problems "live" in the global xy-plane.
      if (Dim()==2 and ((*rotangle)[0]!=0 or (*rotangle)[1]!=0))
        dserror("For 2D problems (xy-plane) the vector ROTANGLE has to be parallel to the global z-axis!");

      // Check, if the updated node positions shall be used for evaluation of the functions 'funct'
      RCP<const Epetra_Vector> dispnp;
      if (((*useUpdatedNodePos)[0] == 1) && (time >= 0.0)){
        dispnp = Discret().GetState("dispnp");
        if (dispnp == Teuchos::null)
          dserror("Locsys: Cannot find state 'dispnp'! You need to set the state 'dispnp' before calling the locsys setup.");
      }

      //Each component j of the pseudo rotation vector that rotates the global xyz system onto the local system
      //assigned to each node consists of a constant, a time dependent and spatially variable part:
      //currotangle_j(x,t) = rotangle_j * curve_j(t) * funct_j(x)
      LINALG::Matrix<3,1> currotangle;
      currotangle.Clear();

      for (int k=0;k<(int)nodes->size();++k)
      {
        bool havenode = Discret().HaveGlobalNode((*nodes)[k]);
        if (!havenode) continue;

        //Weights of rotations vector due to time curve and spatial function
        for (int j=0;j<3;j++)
        {
          // factor given by time curve
          std::vector<double> curvefac(1,1.0);
          int curvenum = (*curve)[j];
          if (curvenum>=0 && time>=0.0)
          {
            curvefac = DRT::Problem::Instance()->Curve(curvenum).FctDer(time,0);
          }
          // factor given by spatial function
          double functfac = 1.0;
          if ((*funct)[j]>0)
          {
            DRT::Node* node = Discret().gNode((*nodes)[k]);

            // Determine node position, which shall be used for evaluating the function, and evaluate it
            if (((*useUpdatedNodePos)[0] == 1) && (time >= 0.0)){
              // Obtain current displacement for node
              std::vector<int> lm;
              Discret().Dof(node,lm);

              std::vector<double> currDisp;
              currDisp.resize(lm.size());

              DRT::UTILS::ExtractMyValues(*dispnp,currDisp,lm);

              // Calculate current position for node
              double currPos[Dim()];
              const double* xp = node->X();

              for (int i=0; i<Dim(); ++i) {
                currPos[i] = xp[i] + currDisp[i];
              }

              // Evaluate function with current node position
              functfac=(DRT::Problem::Instance()->Funct((*funct)[j]-1)).Evaluate(j, &currPos[0], 0.0, &discret_);
            } else {
              // Evaluate function with reference node position
              functfac=(DRT::Problem::Instance()->Funct((*funct)[j]-1)).Evaluate(j, node->X(), 0.0, &discret_);
            }
          }
          currotangle(j)=(*rotangle)[j]*curvefac[0]*functfac;
        }

        nodalrotvectors_[(*nodes)[k]]=currotangle;

        int indices = (*nodes)[k];
        double values  = i;
        locsystoggle_->ReplaceGlobalValues(1,&values,&indices);
      }
    }
    else if (currlocsys->Type() == DRT::Condition::VolumeLocsys ||
             currlocsys->Type() == DRT::Condition::SurfaceLocsys ||
             currlocsys->Type() == DRT::Condition::LineLocsys)
    {
      // already done
    }
    else
      dserror("ERROR: Unknown type of locsys condition!");
  }  // end point locsys

  if (time < 0.0)
    Print(std::cout);

  // When building the transformation matrix we apply a node-by-node
  // strategy. The global matrix trafo_ will consist of nodal blocks
  // of dimension (numdof)x(numdof). The following code block is designed
  // for 2D and 3D solid elements as well as for beam elements of Kirchhoff
  //type (applying nodal tangents). If special fields are constructed with
  // more than dim geometric dofs, i.e. that have to be transformed, then
  //the following code might have to be modified!

  //**********************************************************************
  // Build transformation matrix trafo_
  //**********************************************************************

  // we need to make sure that two nodes sharing the same dofs are not
  // transformed twice. This is a NURBS/periodic boundary feature.
  Teuchos::RCP<Epetra_Vector> already_processed = LINALG::CreateVector(*dofrowmap,true);
  already_processed->PutScalar(0.0);

  // Perform a check for zero diagonal elements. They will crash the SGS-like preconditioners
  bool sanity_check=false;

  // GIDs of all DoFs subjected to local co-ordinate systems
  std::set<int> locsysdofset;

  trafo_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,3));

  for (int i=0;i<noderowmap->NumMyElements();++i)
  {
    int nodeGID = noderowmap->GID(i);
    DRT::Node* node = Discret().gNode(nodeGID);
    if (!node) dserror("ERROR: Cannot find node with gid %",nodeGID);
    std::vector<int> dofs = Discret().Dof(0,node);
    int numdof = (int)dofs.size();
    int locsysindex = (int)(*locsystoggle_)[i];

    // skip nodes whose dofs have already been processed
    bool skip=false;

    for(int rr=0;rr<numdof;++rr)
    {
      if((*already_processed)[dofrowmap->LID(dofs[rr])]>1e-9)
      {
        skip=true;
      }
    }

    if(skip)
    {
      continue;
    }

    // unity matrix for non-locsys node
    if (locsysindex<0)
    {
      for (int r=0;r<numdof;++r)
        trafo_->Assemble(1.0,dofs[r],dofs[r]);
    }
    // trafo matrix for locsys node
    else
    {
      Epetra_SerialDenseMatrix nodetrafo(numdof,numdof);
      for (int k=0;k<numdof;++k) nodetrafo(k,k)=1.0;

      LINALG::Matrix<3,1> currrotvector = nodalrotvectors_[nodeGID];
      LINALG::Matrix<3,3> currrotationmatrix;

      //Compute rotation matrix out of rotation angle
      LARGEROTATIONS::angletotriad(currrotvector, currrotationmatrix);

      //base vectors of local system
      LINALG::Matrix<3,1> vec1;
      LINALG::Matrix<3,1> vec2;
      LINALG::Matrix<3,1> vec3;

      //The columns of the rotation matrix are the base vectors
      for (int j=0;j<3;j++)
      {
        vec1(j)=currrotationmatrix(j,0);
        vec2(j)=currrotationmatrix(j,1);
        vec3(j)=currrotationmatrix(j,2);
      }

      //Check for zero-diagonal elements
      if(fabs(vec1(0))<1e-9 || fabs(vec2(1))<1e-9 || fabs(vec3(2))<1e-9)
      {
        sanity_check=true;
      }

      if (numdof<6) // for solid elements
      {
        // trafo for 2D case
        if (Dim()==2)
        {
          for (int i=0;i<2;i++)
          {
            nodetrafo(0,i)=vec1(i);
            nodetrafo(1,i)=vec2(i);
          }
        }
        // trafo for 3D case
        if (Dim()==3)
        {
          for (int i=0;i<3;i++)
          {
            nodetrafo(0,i)=vec1(i);
            nodetrafo(1,i)=vec2(i);
            nodetrafo(2,i)=vec3(i);
          }
        }
      }
      else //for Kirchhoff beam elements
      {
        for (int i=0;i<3;i++)
        {
          nodetrafo(3,3+i)=vec1(i);   // In this implementation only the nodal tangents (local dofs 3,4,5),
          nodetrafo(4,3+i)=vec2(i);   // which are necessary for clamped ends with arbitrary spatial orientation,
          nodetrafo(5,3+i)=vec3(i);   // but NOT the nodal positions (local dofs 0,1,2) are transformed.
        }
      }

      // Assemble the rotation of this dofs ('nodetrafo') into the global matrix
      for (int r=0;r<numdof;++r)
      {
        for (int c=0;c<numdof;++c)
        {
          trafo_->Assemble(nodetrafo(r,c),dofs[r],dofs[c]);
        }
      }

      // store the DOF with locsys
      for (int r=0;r<numdof;++r)
        locsysdofset.insert(dofs[r]);

      // node dofs are marked now as already processed
      for (int rr=0;rr<numdof;++rr)
      {
        (*already_processed)[dofrowmap->LID(dofs[rr])]=1.0;
      }
    }
  }

  // complete transformation matrix
  trafo_->Complete();

  // throw warning if transformation matrix has zero diagonal elements since
  // they end up on the diagonal of the system matrix
  if(sanity_check)
  {
    printf("Locsys warning:\n");
    printf("A zero diagonal element on the transformation matrix occured on proc %d.\n",Comm().MyPID());
    printf("This will probably cause a crash in the AZTEC preconditioner.\n");
    printf("Try not to rotate your local coordinate system by 90 degrees \n");
    printf("or more or use the slow version.\n");
  }

  //**********************************************************************
  // Build map holding DOFs linked to nodes with local co-ordinate system
  //**********************************************************************

  // create unique/row map of DOFs subjected to local co-ordinate change
  // transformation matrix for relevent DOFs with local system
  int nummyentries = 0;
  int* myglobalentries = NULL;
  std::vector<int> locsysdofs;
  if (locsysdofset.size() > 0)
  {
    locsysdofs.reserve(locsysdofset.size());
    locsysdofs.assign(locsysdofset.begin(), locsysdofset.end());
    nummyentries = locsysdofs.size();
    myglobalentries = &(locsysdofs[0]);
  }
  locsysdofmap_ = Teuchos::rcp(new Epetra_Map(-1, nummyentries, myglobalentries,
                                     discret_.DofRowMap()->IndexBase(),
                                     discret_.Comm()));
  if (locsysdofmap_ == Teuchos::null) dserror("Creation failed.");

  //The matrix subtrafo_ is used in order to apply the Dirichlet Conditions in a more efficient manner
  subtrafo_ = trafo_->ExtractDirichletRows(*locsysdofmap_);

  /*
  REMARK:
   The most general approach to apply Dirichlet conditions in a rotated, local system would be:
   1) Transform the system into local coordinates by means of
      K \cdot D = F --> \tilde{K} \cdot \tilde{D} = \tilde{F}
      with \tilde{K} = trafo_ \cdot K \cdot trafo_^T, \tilde{F} = trafo_ \cdot F, \tilde{D} = trafo_ \cdot D
   2) Apply Dirichlet conditions in the rotated system
   3) Transform the system back into global coordinates, i.e.
      \tilde{K} \cdot \tilde{D} = \tilde{F} --> K \cdot D = F
      with K = trafo_^T \cdot \tilde{K} \cdot trafo_, F = trafo_^T \cdot \tilde{F}, D = trafo_^T \cdot \tilde{D}

   Nevertheless, we apply a more efficient algorithm which can be shown,to deliver an equivalent system of equations:
   1) Therefore we only apply one left transformation to our system of equations according
      K \cdot D = F --> trafo_ \cdot K \cdot D = trafo_ \cdot F
   2) Afterwards we apply the rotated Dirichlet conditions in an appropriate manner, i.e. we zero the corresponding
      Dirichlet line and than insert the corresponding local base vector vec_i of the assigned local system into the
      corresponding 3*3-block, e.g. if the DoFs of the locsys node are represented by fourth, fifth and sixth column:
      (*,*,*,*,*,*,*,*,*,*,*,*) --> (0,0,0, vec_i^T, 0,0,0,0,0,0)
      We don't invert the left transformation of our system afterwards. This means, that we don't solve the original but
      a algebraic manipulated system of equations, nevertheless we still solve for the original, non-rotated DoFs D.
      However, this is actually no drawback since e.g. zero-diagonal elements resulting from rotated Dirichlet
      conditions would still exist even if we applied the back transformation afterwards.
  */


  // done here
  return;
}

/*----------------------------------------------------------------------*
 |  << operator                                               popp 09/08|
 *----------------------------------------------------------------------*/
std::ostream& operator << (std::ostream& os, const DRT::UTILS::LocsysManager& manager)
{
  manager.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print manager (public)                                   meier 06/13|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::Print(std::ostream& os) const
{
  if (Comm().MyPID()==0)
  {
    os << "\n-------------------------------------DRT::UTILS::LocysManager\n";
    for (int i=0;i<NumLocsys();++i)
    {
      printf("*  *  *  *  *  *  *  *  *  *  *  *  *Locsys entity ID: %1d ",locsysconds_[i]->Id());
      if (TypeLocsys(i)==DRT::Condition::PointLocsys)        printf("Point   ");
      else if (TypeLocsys(i)==DRT::Condition::LineLocsys)    printf("Line    ");
      else if (TypeLocsys(i)==DRT::Condition::SurfaceLocsys) printf("Surface ");
      else if (TypeLocsys(i)==DRT::Condition::VolumeLocsys)  printf("Volume  ");
      else dserror("ERROR: Unknown type of locsys condition!");
      printf("\n");
    }
    os << "-------------------------------------------------------------\n\n";
  }
  return;
}


/*----------------------------------------------------------------------*
 |  Transform system global -> local (public)                 popp 09/08|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateGlobalToLocal(
  Teuchos::RCP<LINALG::SparseMatrix> sysmat,
  Teuchos::RCP<Epetra_Vector> rhs)
{
  // transform rhs vector
  RotateGlobalToLocal(rhs);

  // selective multiplication from left
  Teuchos::RCP<LINALG::SparseMatrix> temp = LINALG::Multiply(*subtrafo_,false,*sysmat,false,true);
  // put transformed rows back into global matrix
  sysmat->Put(*temp,1.0,locsysdofmap_);

  return;
}


/*----------------------------------------------------------------------*
 |  Transform system matrix global -> local (public)       mueller 05/10|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateGlobalToLocal(Teuchos::RCP<LINALG::SparseMatrix> sysmat)
{
  // selective multiplication from left
  Teuchos::RCP<LINALG::SparseMatrix> temp = LINALG::Multiply(*subtrafo_,false,*sysmat,false,true);
  // put transformed rows back into global matrix
  sysmat->Put(*temp,1.0,locsysdofmap_);

  return;
}

/*----------------------------------------------------------------------*
 |  Transform vector global -> local (public)                 popp 09/08|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateGlobalToLocal(Teuchos::RCP<Epetra_Vector> vec, bool offset)
{
  //Add an offset value to the displacement vector. This offsett value is needed for Kirchhoff type
  //beam elements, where tangent vectors and not position vectors are rotated!!!
  if (offset)
  {
    AddOffset(vec, false);
  }

  // y = trafo_ . x  with x = vec
  Epetra_Vector tmp(*vec);
  trafo_->Multiply(false,tmp,*vec);
  return;
}

/*----------------------------------------------------------------------*
 |  Transform result + system local -> global (public)        popp 09/08|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateLocalToGlobal(
  Teuchos::RCP<Epetra_Vector> result,
  Teuchos::RCP<LINALG::SparseMatrix> sysmat,
  Teuchos::RCP<Epetra_Vector> rhs)
{
  // transform result
  RotateLocalToGlobal(result);

  // transform rhs vector
  RotateLocalToGlobal(rhs);

  // transform system matrix
  Teuchos::RCP<LINALG::SparseMatrix> temp;
  Teuchos::RCP<LINALG::SparseMatrix> temp2;

  // We want to keep the SaveGraph() value of sysmat also after transformation.
  // It is not possible to keep ExplicitDirichlet()==true after transformation,
  // so we explicitly set this to false.
  temp = LINALG::Multiply(*sysmat,false,*trafo_,false,false,sysmat->SaveGraph(),true);
  temp2 = LINALG::Multiply(*trafo_,true,*temp,false,false,sysmat->SaveGraph(),true);

  // this is a deep copy (expensive!)
  *sysmat = *temp2;

  return;
}

/*----------------------------------------------------------------------*
 |  Transform vector local -> global (public)                 popp 09/08|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateLocalToGlobal(Teuchos::RCP<Epetra_Vector> vec, bool offset)
{
  Epetra_Vector tmp(*vec);
  trafo_->Multiply(true,tmp,*vec);

  //Remove offset value from the displacement vector. This offset value is needed for Kirchhoff type
  //beam elements, where tangent vectors and not position vectors are rotated!!!
  if (offset)
  {
    AddOffset(vec, true);
  }

  return;

}
/*----------------------------------------------------------------------*
 |  Transform matrix local -> global (public)              mueller 05/10|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateLocalToGlobal(Teuchos::RCP<LINALG::SparseMatrix> sysmat)
{

  Teuchos::RCP<LINALG::SparseMatrix> temp;
  Teuchos::RCP<LINALG::SparseMatrix> temp2;
  temp = LINALG::Multiply(*sysmat,false,*trafo_,false,false,sysmat->SaveGraph(),true);
  temp2 = LINALG::Multiply(*trafo_,true,*temp,false,false,sysmat->SaveGraph(),true);
  *sysmat = *temp2;

  return;
}
/*----------------------------------------------------------------------*
 |  Add displacement offset value                            meier 06/13|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::AddOffset(Teuchos::RCP<Epetra_Vector> vec, bool inverse)
{
  // get dof row map of discretization
  const Epetra_Map* dofrowmap = discret_.DofRowMap();
  for (int i=0; i<discret_.NumMyRowNodes();i++)
  {
    int GID = 0;
    int numdofpn=(discret_.Dof(0,discret_.lRowNode(i))).size();

    //nothing to do for 2D/3D solid elements
    if (numdofpn<6) continue;

    DRT::Node* currnode = discret_.lRowNode(i);
    if (currnode->IsCosserat()==false)
      dserror("No rotational DoFs for node %i initialized! Did you set ROTANGLE  for the nodes in your input file?", (i+1));

    LINALG::Matrix<3,1> rot_angle;
    LINALG::Matrix<3,3> mat_sys;

    for (int j=0;j<3;j++)
    {
      rot_angle(j) = currnode->X()[3+j];
    }

    LARGEROTATIONS::angletotriad(rot_angle, mat_sys);

    for (int k=0;k<3;k++)
    {
      GID = (discret_.Dof(0,discret_.lRowNode(i)))[3 + k];
      if (inverse == false)
      {
        (*vec)[dofrowmap->LID(GID)]+=mat_sys(k,0);
      }
      else
      {
        (*vec)[dofrowmap->LID(GID)]-=mat_sys(k,0);
      }
    }
  }

  return;
}
