/*!----------------------------------------------------------------------
\file drt_locsys.cpp

\brief Class controlling local coordinate systems on points, lines and
surfaces and supplying all necessary transformation methods for parallel
vectors and matrices.

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/

#include "drt_locsys.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_function.H"



/*-------------------------------------------------------------------*
 |  ctor (public)                                         bborn 12/08|
 *-------------------------------------------------------------------*/
DRT::UTILS::LocsysManager::LocsysManager(DRT::Discretization& discret, const bool transformleftonly):
discret_(discret),
transformleftonly_(transformleftonly),
type_(DRT::UTILS::LocsysManager::def),
dim_(-1),
numlocsys_(-1)
{
  Setup();
}



/*-------------------------------------------------------------------*
 |  set-up                                                 popp 09/08|
 *-------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::Setup()
{
  // IMPORTANT NOTE:
  // The definition of local coordinate systems only makes sense in
  // combination with Dirichlet boundary conditions. This means that
  // in order to define a symmetry boundary condition, both locsys
  // AND Dirichlet condition have to formulated for the same entity
  // (i.e. point, line, surface, volume). Due to the coordinate trafo
  // the Dirichlet DOFs in the input file transform: x->n, y->t1, z->t2

  // STILL MISSING:
  // - Testing on fluid problems in 2D and 3D
  // - Check, if D.B.C condition is defined wherever Locsys is defined
  // - Efficiency (at the moment, we use standard methods like LINALG::
  //   Multiply for the transformation of global vectors and matrices)

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
  type_.resize(numlocsys_);
  id_.resize(numlocsys_);
  typelocsys_.resize(numlocsys_);
  normals_.Reshape(numlocsys_,3);
  tangents_.Reshape(numlocsys_,3);
  origins_.Reshape(numlocsys_,3);
  thirddir_.Reshape(numlocsys_,3);

  for (int i=0; i<NumLocsys(); ++i)
  {
    id_[i] = locsysconds_[i]->Id();
    const std::string* type = locsysconds_[i]->Get<std::string>("Type");
    if (*type=="default")
      type_[i] = DRT::UTILS::LocsysManager::def;
    else if (*type=="FunctionEvaluation")
      type_[i] = DRT::UTILS::LocsysManager::functionevaluation;
    else if (*type=="OriginRadialSliding")
      type_[i] = DRT::UTILS::LocsysManager::originradialsliding;
    else dserror("Unknown type of locsys");
  }

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

      const vector<double>* n = currlocsys->Get<vector<double> >("normal");
      const vector<double>* t = currlocsys->Get<vector<double> >("tangent");
      double ln = sqrt((*n)[0]*(*n)[0] + (*n)[1]*(*n)[1] + (*n)[2]*(*n)[2]);
      double lt = sqrt((*t)[0]*(*t)[0] + (*t)[1]*(*t)[1] + (*t)[2]*(*t)[2]);
      const vector<int>* nodes = currlocsys->Nodes();

      // check for sanity of input data
      if (Dim()==2)
      {
        // volume locsys makes no sense for 2D
        dserror("ERROR: Volume locsys definition for 2D problem");
      }
      else
      {
        // normal has to be provided
        if (ln==0.0 || n->size() != 3)
          dserror("ERROR: No normal provided for 3D locsys definition");

        // tangent has to be provided
        if (lt==0.0 || t->size() != 3)
          dserror("ERROR: No tangent provided for 3D locsys definition");

        // tangent has to be orthogonal
        double ndott = (*n)[0]*(*t)[0] + (*n)[1]*(*t)[1] + (*n)[2]*(*t)[2];
        if (abs(ndott)>1.0e-8)
          dserror("ERROR: Locsys normal and tangent not orthogonal");
      }

      // build unit normal and tangent
      for (int k=0;k<3;++k)
      {
        normals_(i,k)  = (*n)[k]/ln;
        tangents_(i,k) = (*t)[k]/lt;
      }

      // build third direction (corkscrew rule)
      thirddir_(i,0) = normals_(i,1)*tangents_(i,2)-normals_(i,2)*tangents_(i,1);
      thirddir_(i,1) = normals_(i,2)*tangents_(i,0)-normals_(i,0)*tangents_(i,2);
      thirddir_(i,2) = normals_(i,0)*tangents_(i,1)-normals_(i,1)*tangents_(i,0);

      // build locsystoggle vector with locsys IDs
      if (type_[i] == DRT::UTILS::LocsysManager::def)
        for (int k=0;k<(int)nodes->size();++k)
        {
          bool havenode = Discret().HaveGlobalNode((*nodes)[k]);
          if (!havenode) continue;

          int indices = (*nodes)[k];
          double values  = i;
          locsystoggle_->ReplaceGlobalValues(1,&values,&indices);
        }
      if (type_[i] == DRT::UTILS::LocsysManager::functionevaluation)
      {
      	// read origin from input file
      	const vector<double>* o = currlocsys->Get<vector<double> >("origin");
      	if(o->size() != 3)
      		dserror("ERROR: No origin provided for locsys definition of type FunctionEvaluation");
      	for (int k=0;k<origins_.N();k++)
      		origins_(i,k) = (*o)[k];
      }
      // special locsys case: originradialsliding
      if (type_[i] == DRT::UTILS::LocsysManager::originradialsliding)
			{
				// read origin from input file
				const vector<double>* o = currlocsys->Get<vector<double> >("origin");
				if(o->size() != 3)
					dserror("ERROR: No origin provided for locsys definition of type OriginRadialSliding");
				// store locsys position (to identify this locsys type later on)
				radslideids_.push_back(i);
				for (int k=0;k<origins_.M();k++)
					origins_(i,k) = (*o)[k];

				for (int k=0;k<(int)nodes->size();++k)
				{
					bool havenode = Discret().HaveGlobalNode((*nodes)[k]);
					if (!havenode) continue;
					int indices = (*nodes)[k];
					double values  = i;
					locsystoggle_->ReplaceGlobalValues(1,&values,&indices);
				}
			}
    }
    else if (currlocsys->Type() == DRT::Condition::SurfaceLocsys ||
             currlocsys->Type() == DRT::Condition::LineLocsys ||
             currlocsys->Type() == DRT::Condition::PointLocsys)
    {
      // do nothing yet
    }
    else
      dserror("ERROR: Unknown type of locsys condition!");
  }

  //**********************************************************************
  // read surface locsys conditions
  //**************************+*******************************************
  for (int i=0; i<NumLocsys(); ++i)
  {
    DRT::Condition* currlocsys = locsysconds_[i];

    if (currlocsys->Type() == DRT::Condition::SurfaceLocsys)
    {
      typelocsys_[i] = DRT::Condition::SurfaceLocsys;

      const vector<double>* n = currlocsys->Get<vector<double> >("normal");
      const vector<double>* t = currlocsys->Get<vector<double> >("tangent");
      double ln = sqrt((*n)[0]*(*n)[0] + (*n)[1]*(*n)[1] + (*n)[2]*(*n)[2]);
      double lt = sqrt((*t)[0]*(*t)[0] + (*t)[1]*(*t)[1] + (*t)[2]*(*t)[2]);
      const vector<int>* nodes = currlocsys->Nodes();

      // check for sanity of input data
      if (Dim()==2)
      {
        // normal has to be provided
        if (ln==0.0 || n->size()!=3)
          dserror("ERROR: No normal provided for 2D locsys definition");

        // normal must lie in x1x2-plane
        if ((*n)[2]!=0.0)
          dserror("ERROR: Invalid normal provided for 2D locsys definition");

        // tangent must not be provided
        if (lt!=0.0)
          dserror("ERROR: Tangent provided for 2D locsys definition");

        // build unit normal and tangent
        normals_(i,0)  = (*n)[0]/ln;
        normals_(i,1)  = (*n)[1]/ln;

        // build unit tangent from 2D orthogonality
        tangents_(i,0)  = -(*n)[1]/ln;
        tangents_(i,1)  =  (*n)[0]/ln;
      }
      else
      {
        // normal has to be provided
        if (ln==0.0 || n->size()!=3)
          dserror("ERROR: No normal provided for 3D locsys definition");

        // tangent has to be provided
        if (lt==0.0 || t->size()!=3)
          dserror("ERROR: No tangent provided for 3D locsys definition");

        // tangent has to be orthogonal
        double ndott = (*n)[0]*(*t)[0] + (*n)[1]*(*t)[1] + (*n)[2]*(*t)[2];
        if (abs(ndott)>1.0e-8)
          dserror("ERROR: Locsys normal and tangent not orthogonal");

        // build unit normal and tangent
        for (int k=0;k<3;++k)
        {
          normals_(i,k)  = (*n)[k]/ln;
          tangents_(i,k) = (*t)[k]/lt;
        }
      }

      // build third direction (corkscrew rule)
      thirddir_(i,0) = normals_(i,1)*tangents_(i,2)-normals_(i,2)*tangents_(i,1);
      thirddir_(i,1) = normals_(i,2)*tangents_(i,0)-normals_(i,0)*tangents_(i,2);
      thirddir_(i,2) = normals_(i,0)*tangents_(i,1)-normals_(i,1)*tangents_(i,0);

      // build locsystoggle vector with locsys IDs
      if (type_[i] == DRT::UTILS::LocsysManager::def)
        for (int k=0;k<(int)nodes->size();++k)
        {
          bool havenode = Discret().HaveGlobalNode((*nodes)[k]);
          if (!havenode) continue;

          int indices = (*nodes)[k];
          double values  = i;
          locsystoggle_->ReplaceGlobalValues(1,&values,&indices);
        }
      if (type_[i] == DRT::UTILS::LocsysManager::functionevaluation)
      {
      	// read origin from input file
      	const vector<double>* o = currlocsys->Get<vector<double> >("origin");
      	if(o->size() != 3)
      		dserror("ERROR: No origin provided for locsys definition of type FunctionEvaluation");

      	for (int k=0;k<origins_.N();k++)
      		origins_(i,k) = (*o)[k];
      }
      if (type_[i] == DRT::UTILS::LocsysManager::originradialsliding)
			{
				// read origin from input file
				const vector<double>* o = currlocsys->Get<vector<double> >("origin");
				if(o->size() != 3)
					dserror("ERROR: No origin provided for locsys definition of type OriginRadialSliding");
				// store locsys position (to identify this locsys type later on)
				radslideids_.push_back(i);
				for (int k=0;k<origins_.N();k++)
					origins_(i,k) = (*o)[k];

				for (int k=0;k<(int)nodes->size();++k)
				{
					bool havenode = Discret().HaveGlobalNode((*nodes)[k]);
					if (!havenode) continue;
					int indices = (*nodes)[k];
					double values  = i;
					locsystoggle_->ReplaceGlobalValues(1,&values,&indices);
				}
			}
    }
    else if (currlocsys->Type() == DRT::Condition::VolumeLocsys ||
             currlocsys->Type() == DRT::Condition::LineLocsys ||
             currlocsys->Type() == DRT::Condition::PointLocsys)
    {
      // do nothing yet
    }
    else
      dserror("ERROR: Unknown type of locsys condition!");
  }

  //**********************************************************************
  // read line locsys conditions
  //**********************************************************************
  for (int i=0; i<NumLocsys(); ++i)
  {
    DRT::Condition* currlocsys = locsysconds_[i];

    if (currlocsys->Type() == DRT::Condition::LineLocsys)
    {
      typelocsys_[i] = DRT::Condition::LineLocsys;

      const vector<double>* n = currlocsys->Get<vector<double> >("normal");
      const vector<double>* t = currlocsys->Get<vector<double> >("tangent");
      double ln = sqrt((*n)[0]*(*n)[0] + (*n)[1]*(*n)[1] + (*n)[2]*(*n)[2]);
      double lt = sqrt((*t)[0]*(*t)[0] + (*t)[1]*(*t)[1] + (*t)[2]*(*t)[2]);
      const vector<int>* nodes = currlocsys->Nodes();

      // check for sanity of input data
      if (Dim()==2)
      {
        // normal has to be provided
        if (ln==0.0)
          dserror("ERROR: No normal provided for 2D locsys definition");

        // normal must lie in x1x2-plane
        if ((*n)[2]!=0.0)
          dserror("ERROR: Invalid normal provided for 2D locsys definition");

        // tangent must not be provided
        if (lt!=0.0)
          dserror("ERROR: Tangent provided for 2D locsys definition");

        // build unit normal and tangent
        normals_(i,0)  = (*n)[0]/ln;
        normals_(i,1)  = (*n)[1]/ln;

        // build unit tangent from 2D orthogonality
        tangents_(i,0)  = -(*n)[1]/ln;
        tangents_(i,1)  =  (*n)[0]/ln;
      }
      else
      {
        // normal has to be provided
        if (ln==0.0)
          dserror("ERROR: No normal provided for 3D locsys definition");

        // tangent has to be provided
        if (lt==0.0)
          dserror("ERROR: No tangent provided for 3D locsys definition");

        // tangent has to be orthogonal
        double ndott = (*n)[0]*(*t)[0] + (*n)[1]*(*t)[1] + (*n)[2]*(*t)[2];
        if (abs(ndott)>1.0e-8)
          dserror("ERROR: Locsys normal and tangent not orthogonal");

        // build unit normal and tangent
        for (int k=0;k<3;++k)
        {
          normals_(i,k)  = (*n)[k]/ln;
          tangents_(i,k) = (*t)[k]/lt;
        }
      }

      // build third direction (corkscrew rule)
      thirddir_(i,0) = normals_(i,1)*tangents_(i,2)-normals_(i,2)*tangents_(i,1);
      thirddir_(i,1) = normals_(i,2)*tangents_(i,0)-normals_(i,0)*tangents_(i,2);
      thirddir_(i,2) = normals_(i,0)*tangents_(i,1)-normals_(i,1)*tangents_(i,0);

      // build locsystoggle vector with locsys IDs
      if (type_[i] == DRT::UTILS::LocsysManager::def)
        for (int k=0;k<(int)nodes->size();++k)
        {
          bool havenode = Discret().HaveGlobalNode((*nodes)[k]);
          if (!havenode) continue;

          int indices = (*nodes)[k];
          double values  = i;
          locsystoggle_->ReplaceGlobalValues(1,&values,&indices);
        }
      if (type_[i] == DRT::UTILS::LocsysManager::functionevaluation)
      {
      	// read origin from input file
      	const vector<double>* o = currlocsys->Get<vector<double> >("origin");
      	if(o->size() != 3)
      		dserror("ERROR: No origin provided for locsys definition of type FunctionEvaluation");
      	for (int k=0;k<origins_.N();k++)
      		origins_(i,k) = (*o)[k];
      }
      if (type_[i] == DRT::UTILS::LocsysManager::originradialsliding)
			{
				// read origin from input file
				const vector<double>* o = currlocsys->Get<vector<double> >("origin");
				if(o->size() != 3)
					dserror("ERROR: No origin provided for locsys definition of type OriginRadialSliding");
				// store locsys position (to identify this locsys type later on)
				radslideids_.push_back(i);
				for (int k=0;k<origins_.N();k++)
					origins_(i,k) = (*o)[k];

				for (int k=0;k<(int)nodes->size();++k)
				{
					bool havenode = Discret().HaveGlobalNode((*nodes)[k]);
					if (!havenode) continue;
					int indices = (*nodes)[k];
					double values  = i;
					locsystoggle_->ReplaceGlobalValues(1,&values,&indices);
				}
			}
    }
    else if (currlocsys->Type() == DRT::Condition::VolumeLocsys ||
             currlocsys->Type() == DRT::Condition::SurfaceLocsys ||
             currlocsys->Type() == DRT::Condition::PointLocsys)
    {
      // already done or do nothing yet
    }
    else
      dserror("ERROR: Unknown type of locsys condition!");
  }

  //**********************************************************************
  // read point locsys conditions
  //**********************************************************************
  for (int i=0; i<NumLocsys(); ++i)
  {
    DRT::Condition* currlocsys = locsysconds_[i];

    if (currlocsys->Type() == DRT::Condition::PointLocsys)
    {
      typelocsys_[i] = DRT::Condition::PointLocsys;

      const vector<double>* n = currlocsys->Get<vector<double> >("normal");
      const vector<double>* t = currlocsys->Get<vector<double> >("tangent");
      double ln = sqrt((*n)[0]*(*n)[0] + (*n)[1]*(*n)[1] + (*n)[2]*(*n)[2]);
      double lt = sqrt((*t)[0]*(*t)[0] + (*t)[1]*(*t)[1] + (*t)[2]*(*t)[2]);
      const vector<int>* nodes = currlocsys->Nodes();

      // check for sanity of input data
      if (Dim()==2)
      {
        // normal has to be provided
        if (ln==0.0)
          dserror("ERROR: No normal provided for 2D locsys definition");

        // normal must lie in x1x2-plane
        if ((*n)[2]!=0.0)
          dserror("ERROR: Invalid normal provided for 2D locsys definition");

        // tangent must not be provided
        if (lt!=0.0)
          dserror("ERROR: Tangent provided for 2D locsys definition");

        // build unit normal and tangent
        normals_(i,0)  = (*n)[0]/ln;
        normals_(i,1)  = (*n)[1]/ln;

        // build unit tangent from 2D orthogonality
        tangents_(i,0)  = -(*n)[1]/ln;
        tangents_(i,1)  =  (*n)[0]/ln;
      }
      else
      {
        // normal has to be provided
        if (ln==0.0)
          dserror("ERROR: No normal provided for 3D locsys definition");

        // tangent has to be provided
        if (lt==0.0)
          dserror("ERROR: No tangent provided for 3D locsys definition");

        // tangent has to be orthogonal
        double ndott = (*n)[0]*(*t)[0] + (*n)[1]*(*t)[1] + (*n)[2]*(*t)[2];
        if (abs(ndott)>1.0e-8)
          dserror("ERROR: Locsys normal and tangent not orthogonal");

        // build unit normal and tangent
        for (int k=0;k<3;++k)
        {
          normals_(i,k)  = (*n)[k]/ln;
          tangents_(i,k) = (*t)[k]/lt;
        }
      }

      // build third direction (corkscrew rule)
      thirddir_(i,0) = normals_(i,1)*tangents_(i,2)-normals_(i,2)*tangents_(i,1);
      thirddir_(i,1) = normals_(i,2)*tangents_(i,0)-normals_(i,0)*tangents_(i,2);
      thirddir_(i,2) = normals_(i,0)*tangents_(i,1)-normals_(i,1)*tangents_(i,0);

      // build locsystoggle vector with locsys IDs
      if (type_[i] == DRT::UTILS::LocsysManager::def)
        for (int k=0;k<(int)nodes->size();++k)
        {
          bool havenode = Discret().HaveGlobalNode((*nodes)[k]);
          if (!havenode) continue;

          int indices = (*nodes)[k];
          double values  = i;
          locsystoggle_->ReplaceGlobalValues(1,&values,&indices);
        }
      if (type_[i] == DRT::UTILS::LocsysManager::functionevaluation)
      {
      	// read origin from input file
      	const vector<double>* o = currlocsys->Get<vector<double> >("origin");
      	if(o->size() != 3)
      		dserror("ERROR: No origin provided for locsys definition of type FunctionEvaluation");
      	for (int k=0;k<origins_.N();k++)
      		origins_(i,k) = (*o)[k];
      }
      if (type_[i] == DRT::UTILS::LocsysManager::originradialsliding)
			{
				// read origin from input file
				const vector<double>* o = currlocsys->Get<vector<double> >("origin");
				if(o->size() != 3)
					dserror("ERROR: No origin provided for locsys definition of type OriginRadialSliding");
				// store locsys position (to identify this locsys type later on)
				radslideids_.push_back(i);
				for (int k=0;k<origins_.N();k++)
					origins_(i,k) = (*o)[k];

				for (int k=0;k<(int)nodes->size();++k)
				{
					bool havenode = Discret().HaveGlobalNode((*nodes)[k]);
					if (!havenode) continue;
					int indices = (*nodes)[k];
					double values  = i;
					locsystoggle_->ReplaceGlobalValues(1,&values,&indices);
				}
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
  }

  Print(cout);

  // When building the transformation matrix we apply a node-by-node
  // strategy. The global matrix trafo_ will consist of nodal blocks
  // of dimension (numdof)x(numdof). To be consistent with as many
  // field types (structure, fluid) as possible, the following rule
  // holds: "For a dim-dimensional problem, the first dim dofs are
  // transformed whereas all other dofs remain untouched. This e.g.
  // accounts for the fact that a node in a 3D fluid discretization
  // has 4 dofs: 3 velocity dofs (to be transformed) and 1 pressure
  // dofs (isotropic). If special fields are constructed with more than
  // dim geometric dofs, i.e. that have to be transformed, then the
  // following code might have to be modified!

  //**********************************************************************
  // Build transformation matrix trafo_
  //**********************************************************************

  // get dof row layout of discretization
  const Epetra_Map* dofrowmap = discret_.DofRowMap();

  // we need to make sure that two nodes sharing the same dofs are not
  // transformed twice. This is a NURBS/periodic boundary feature.
  RCP<Epetra_Vector> already_processed = LINALG::CreateVector(*dofrowmap,true);
  already_processed->PutScalar(0.0);

  // for transformleftonly_ option, perform a check for zero diagonal
  // elements. They will crash the SGS-like preconditioners
  bool sanity_check=false;

  // number of nodes subjected to local co-ordinate systems
  set<int> locsysdofset;

  trafo_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,3));

  for (int i=0;i<noderowmap->NumMyElements();++i)
  {
    int gid = noderowmap->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    vector<int> dofs = Discret().Dof(node);
    int numdof = (int)dofs.size();
    int locsysindex = (int)(*locsystoggle_)[i];

    // skip nodes whos dofs have already been processed
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

    // unity matrix for non-locsys node (programmed late at night)
    if (locsysindex<0)
    {
      for (int r=0;r<numdof;++r)
        trafo_->Assemble(1.0,dofs[r],dofs[r]);
    }
    // trafo matrix for locsys node
    else
    {

      // first create an identity matrix of size numdof --- the default node
      // trafo, applying in particular to pressure dofs in fluid problems
      Epetra_SerialDenseMatrix nodetrafo(numdof,numdof);

      for (int k=0;k<numdof;++k) nodetrafo(k,k)=1.0;

      // ---------------------------------------------------
      //
      // tnb-vectors (tangent, normal, binormal)
      //
      // o in the default case without a spatial function, they are just
      //   copies of the vectors stored in the condition
      // o in the case of spatial functions applied to the local coordinate
      //   system, they are rotated according to the prescribed function
      //
      LINALG::Matrix<3,1> n;
      LINALG::Matrix<3,1> t;
      LINALG::Matrix<3,1> b;

      n.Clear();
      t.Clear();
      b.Clear();

      for (int k=0;k<Dim();++k)
      {
        n(k)=normals_ (locsysindex,k);
        t(k)=tangents_(locsysindex,k);
      }
      if (Dim()==3)
      {
        for (int k=0;k<Dim();++k)
        {
          b(k)=thirddir_(locsysindex,k);
        }
      }

      DRT::Condition* currlocsys = locsysconds_[locsysindex];
      const vector<int>*   funct = currlocsys->Get<vector<int> >("(axis,angle)-funct");
      if(funct)
      {
        // sanity checks
        if(funct->size()!=2)
        {
          dserror("two functions are required. One three-component vector valued (axis), one scalar (angle)\n");
        }

        if((*funct)[0]>0 && (*funct)[1]>0)
        {
          // check number of components in both spatial functions
          if((DRT::Problem::Instance()->Funct((*funct)[0]-1)).NumberComponents()!=3)
          {
            dserror("expecting vector-valued function for rotation axis of local coordinate system, but got %d components\n",(DRT::Problem::Instance()->Funct((*funct)[0]-1)).NumberComponents());
          }
          if((DRT::Problem::Instance()->Funct((*funct)[1]-1)).NumberComponents()!=1)
          {
            dserror("expecting scalar function for rotation angle of local coordinate system, got %d components\n",(DRT::Problem::Instance()->Funct((*funct)[1]-1)).NumberComponents());
          }

          // if necessary, the local coordinate system is rotated according to
          // a given spatial functions (for curved surfaces etc.)
          const double time =0.0;


          const double angle=(DRT::Problem::Instance()->Funct((*funct)[1]-1)).Evaluate(
            0        ,
            node->X(),
            time     ,
            &discret_);
          const double x=(DRT::Problem::Instance()->Funct((*funct)[0]-1)).Evaluate(
            0        ,
            node->X(),
            time     ,
            &discret_);
          const double y=(DRT::Problem::Instance()->Funct((*funct)[0]-1)).Evaluate(
            1        ,
            node->X(),
            time     ,
            &discret_);
          const double z=(DRT::Problem::Instance()->Funct((*funct)[0]-1)).Evaluate(
            2        ,
            node->X(),
            time     ,
            &discret_);

          // set rotation matrix R_
          SetAxisRotation(x,y,z,angle);

          // compute rotated local coordinate system according to precomputed
          // rotation matrix R_
          LocalRotation(n);
          LocalRotation(t);
          LocalRotation(b);
        } // end if(funct)

        if(transformleftonly_)
        {
          // sanity check
          if(fabs(n(0))<1e-9 || fabs(t(1))<1e-9 || fabs(b(2))<1e-9)
          {
            sanity_check=true;
          }
        }
      }

      // case: originradialsliding
      // shift between def and originradialsliding
      bool radslidetoggle = false;
      // look for existing originradialsliding locsys
      for(int k=0;k<(int)radslideids_.size();k++)
        if(locsysindex==radslideids_.at(k))
        {
          currlocsys = locsysconds_[locsysindex];
          radslidetoggle = true;
          continue;
        }

      // ---------------------------------------------------
      // set rotation of this nodes dofs ('nodetrafo')
      /* case: originradialsliding:
       * transformation matrix set up
       * 	- n: standard (outward) normal given by input file
       *  - radialfirstdir: direction given by origin and coordinates of current node
       *  - radialsecdir: binormal to n and radialfirstdir
       */
      if(radslidetoggle)
      {
        vector<double> radialdir;
        vector<double> thirddir;
        radialdir.resize(3, 0.0);
        thirddir.resize(3, 0.0);

        // get nodal coordinates
        const double *nodecoords;
        //lenght of radial vector
        double lrad;

        nodecoords = Discret().gNode(gid)->X();
        // length of radial vector
        lrad = sqrt((nodecoords[0]-origins_(locsysindex,0))*(nodecoords[0]-origins_(locsysindex,0))+
                    (nodecoords[1]-origins_(locsysindex,1))*(nodecoords[1]-origins_(locsysindex,1))+
                    (nodecoords[2]-origins_(locsysindex,2))*(nodecoords[2]-origins_(locsysindex,2)));
        // construct radial unit vector
        for(int l=0;l<3;l++)
          radialdir.at(l) = (nodecoords[l]-origins_(locsysindex,l))/lrad;
        // construct second direction (radialdir x normal)
        thirddir.at(0) = radialdir.at(1)*n(2)-radialdir.at(2)*n(1);
        thirddir.at(1) = radialdir.at(2)*n(0)-radialdir.at(0)*n(2);
        thirddir.at(2) = radialdir.at(0)*n(1)-radialdir.at(1)*n(0);

        // trafo for 2D and 3D case
        nodetrafo(0,0) = n(0);
        nodetrafo(0,1) = n(1);
        nodetrafo(1,0) = radialdir.at(0);
        nodetrafo(1,1) = radialdir.at(1);

        // additional trafo for 3D case
        if (Dim()==3)
        {
          nodetrafo(0,2) = n(2);
          nodetrafo(1,2) = radialdir.at(2);
          nodetrafo(2,0) = thirddir.at(0);
          nodetrafo(2,1) = thirddir.at(1);
          nodetrafo(2,2) = thirddir.at(2);
        }

        // fluid case, 4x4
        if(numdof==4)
        {
          for(int i=0;i<(numdof-1);i++)
          {
            nodetrafo(3,i) = 0.0;
            nodetrafo(i,3) = 0.0;
          }
          nodetrafo(3,3) = 1.0;
        }

        for (int r=0;r<numdof;++r)
          for (int c=0;c<numdof;++c)
            trafo_->Assemble(nodetrafo(r,c),dofs[r],dofs[c]);

        // store the DOF with locsys
        if (transformleftonly_)
          for (int r=0;r<numdof;++r)
            locsysdofset.insert(dofs[r]);
      }
      // build trafo_ conventionally
      else
      {
				// trafo for 2D and 3D case
				nodetrafo(0,0)=n(0);
				nodetrafo(0,1)=n(1);
				nodetrafo(1,0)=t(0);
				nodetrafo(1,1)=t(1);

				// additional trafo for 3D case
				if (Dim()==3)
				{
					nodetrafo(0,2)=n(2);
					nodetrafo(1,2)=t(2);
					nodetrafo(2,0)=b(0);
					nodetrafo(2,1)=b(1);
					nodetrafo(2,2)=b(2);
				}

        // fluid case, 4x4
        if(numdof==4)
        {
          for(int i=0;i<(numdof-1);i++)
          {
            nodetrafo(3,i) = 0.0;
            nodetrafo(i,3) = 0.0;
          }
          nodetrafo(3,3) = 1.0;
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
				if (transformleftonly_)
					for (int r=0;r<numdof;++r)
						locsysdofset.insert(dofs[r]);
      }
    }

    // node dofs are marked now as already processed
    for (int rr=0;rr<numdof;++rr)
    {
      (*already_processed)[dofrowmap->LID(dofs[rr])]=1.0;
    }
  }

  // complete transformation matrix
  trafo_->Complete();

  // throw warning if transformation matrix has zero diagonal elements since
  // they end up on the diagonal of the system matrix in the fast
  // transformleftonly_ case
  if(sanity_check && transformleftonly_)
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
  if (transformleftonly_)
  {
    int nummyentries = 0;
    int* myglobalentries = NULL;
    vector<int> locsysdofs;
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
    if (locsysdofmap_ == null) dserror("Creation failed.");

    subtrafo_ = trafo_->ExtractDirichletRows(*locsysdofmap_);
    //cout << "Subtrafo: nummyrows=" << subtrafo_->EpetraMatrix()->NumMyRows()
    //     << " nummycols=" << subtrafo_->EpetraMatrix()->NumMyCols()
    //     << endl;
  }

  // done here
  return;
}

/*----------------------------------------------------------------------*
 |  << operator                                               popp 09/08|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const DRT::UTILS::LocsysManager& manager)
{
  manager.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print manager (public)                                    popp 09/08|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::Print(ostream& os) const
{
  if (Comm().MyPID()==0)
  {
    os << "\n-------------------------------------DRT::UTILS::LocysManager\n";
    for (int i=0;i<NumLocsys();++i)
    {
      printf("Locsys entity ID: %3d ",locsysconds_[i]->Id());
      if (TypeLocsys(i)==DRT::Condition::PointLocsys)        printf("Point   ");
      else if (TypeLocsys(i)==DRT::Condition::LineLocsys)    printf("Line    ");
      else if (TypeLocsys(i)==DRT::Condition::SurfaceLocsys) printf("Surface ");
      else if (TypeLocsys(i)==DRT::Condition::VolumeLocsys)  printf("Volume  ");
      else dserror("ERROR: Unknown type of locsys condition!");
      printf("Normal: %8.4e %8.4e %8.4e Tangent1: %8.4e %8.4e %8.4e Tangent2: %8.4e %8.4e %8.4e Type: ",
      normals_(i,0),normals_(i,1),normals_(i,2),
      tangents_(i,0),tangents_(i,1),tangents_(i,2),
      thirddir_(i,0),thirddir_(i,1),thirddir_(i,2));
      switch(type_[i])
      {
        case DRT::UTILS::LocsysManager::def:                 printf("default");
        break;
        case DRT::UTILS::LocsysManager::functionevaluation:  printf("functionevaluation");
        break;
        case DRT::UTILS::LocsysManager::originradialsliding: printf("originradialsliding");
        break;
        default:
          dserror("Unknown type of locsys");
        break;
      }
      printf("\n");
    }
    os << "-------------------------------------------------------------\n\n";
  }
  return;
}


/*----------------------------------------------------------------------*
 |  Transform system global -> local (public)                 popp 09/08|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateGlobalToLocal(RCP<LINALG::SparseMatrix> sysmat,
                                                    RCP<Epetra_Vector> rhs)
{
  // transform rhs vector
  RotateGlobalToLocal(rhs);

  // transform system matrix
  if (transformleftonly_)
  {
    // selective multiplication from left
    RCP<LINALG::SparseMatrix> temp = LINALG::Multiply(*subtrafo_,false,*sysmat,false,true);
    // put transformed rows back into global matrix
    sysmat->Put(*temp,1.0,locsysdofmap_);
  }
  else
  {
    RCP<LINALG::SparseMatrix> temp;
    RCP<LINALG::SparseMatrix> temp2;

    // We want to keep the SaveGraph() value of sysmat also after transformation.
    // It is not possible to keep ExplicitDirichlet()==true after transformation,
    // so we explicitly set this to false.
    temp = LINALG::Multiply(*trafo_,false,*sysmat,false,false,sysmat->SaveGraph(),true);  // multiply from left
    temp2 = LINALG::Multiply(*temp,false,*trafo_,true,false,sysmat->SaveGraph(),true);  // multiply from right

    // this is a deep copy (expensive!)
    *sysmat = *temp2;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Transform system matrix global -> local (public)       mueller 05/10|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateGlobalToLocal(RCP<LINALG::SparseMatrix> sysmat)
{
  // transform system matrix
  if (transformleftonly_)
  {
    // selective multiplication from left
    RCP<LINALG::SparseMatrix> temp = LINALG::Multiply(*subtrafo_,false,*sysmat,false,true);
    // put transformed rows back into global matrix
    sysmat->Put(*temp,1.0,locsysdofmap_);
  }
  else
  {
    RCP<LINALG::SparseMatrix> temp;
    RCP<LINALG::SparseMatrix> temp2;

    // We want to keep the SaveGraph() value of sysmat also after transformation.
    // It is not possible to keep ExplicitDirichlet()==true after transformation,
    // so we explicitly set this to false.
    temp = LINALG::Multiply(*trafo_,false,*sysmat,false,false,sysmat->SaveGraph(),true);  // multiply from left
    temp2 = LINALG::Multiply(*temp,false,*trafo_,true,false,sysmat->SaveGraph(),true);  // multiply from right

    // this is a deep copy (expensive!)
    *sysmat = *temp2;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Transform vector global -> local (public)                 popp 09/08|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateGlobalToLocal(RCP<Epetra_Vector> vec)
{
  Epetra_Vector tmp(*vec);
  trafo_->Multiply(false,tmp,*vec);
  return;
}


/*----------------------------------------------------------------------*
 |  Transform result + system local -> global (public)        popp 09/08|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateLocalToGlobal(RCP<Epetra_Vector> result,
                                               RCP<LINALG::SparseMatrix> sysmat,
                                               RCP<Epetra_Vector> rhs)
{
  // transform result
  RotateLocalToGlobal(result);

  // transform rhs vector
  RotateLocalToGlobal(rhs);

  // transform system matrix
  RCP<LINALG::SparseMatrix> temp;
  RCP<LINALG::SparseMatrix> temp2;

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
void DRT::UTILS::LocsysManager::RotateLocalToGlobal(RCP<Epetra_Vector> vec)
{
  Epetra_Vector tmp(*vec);
  trafo_->Multiply(true,tmp,*vec);
  return;
}
/*----------------------------------------------------------------------*
 |  Transform matrix local -> global (public)              mueller 05/10|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateLocalToGlobal(RCP<LINALG::SparseMatrix> sysmat)
{
  RCP<LINALG::SparseMatrix> temp;
  RCP<LINALG::SparseMatrix> temp2;
  temp = LINALG::Multiply(*sysmat,false,*trafo_,false,false,sysmat->SaveGraph(),true);
  temp2 = LINALG::Multiply(*trafo_,true,*temp,false,false,sysmat->SaveGraph(),true);
  *sysmat = *temp2;
  return;
}

