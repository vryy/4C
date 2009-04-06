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
#ifdef CCADISCRET

#include "drt_locsys.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"

using namespace std;
using namespace Teuchos;


/*-------------------------------------------------------------------*
 |  ctor (public)                                         bborn 12/08|
 *-------------------------------------------------------------------*/
DRT::UTILS::LocsysManager::LocsysManager(DRT::Discretization& discret, const bool transformleftonly):
discret_(discret),
transformleftonly_(transformleftonly)
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
  const ParameterList& psize = DRT::Problem::Instance()->ProblemSizeParams();
  dim_=psize.get<int>("DIM");
  if (Dim()!= 2 && Dim()!=3) dserror("ERROR: Locsys problem must be 2D or 3D");
    
  // get node row layout of discretization
  const Epetra_Map* noderowmap = discret_.NodeRowMap();
  
  // create locsys vector and initialize to -1
  locsystoggle_ = LINALG::CreateVector(*noderowmap,false);
  for (int i=0;i<noderowmap->NumMyElements();++i)
    (*locsystoggle_)[i]= -1;
  
  // check for locsys boundary conditions
  Discret().GetCondition("Locsys",locsysconds_);
  numlocsys_ = (int)locsysconds_.size();
  typelocsys_.resize(numlocsys_);
  normals_.Reshape(numlocsys_,3);
  tangents_.Reshape(numlocsys_,3);
  thirddir_.Reshape(numlocsys_,3);
  
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
      for (int k=0;k<(int)nodes->size();++k)
      {
        bool havenode = Discret().HaveGlobalNode((*nodes)[k]);
        if (!havenode) continue;
        
        int indices = (*nodes)[k];
        double values  = i;
        locsystoggle_->ReplaceGlobalValues(1,&values,&indices);
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
      for (int k=0;k<(int)nodes->size();++k)
      {
        bool havenode = Discret().HaveGlobalNode((*nodes)[k]);
        if (!havenode) continue;
        
        int indices = (*nodes)[k];
        double values  = i;
        locsystoggle_->ReplaceGlobalValues(1,&values,&indices);
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
      for (int k=0;k<(int)nodes->size();++k)
      {
        bool havenode = Discret().HaveGlobalNode((*nodes)[k]);
        if (!havenode) continue;
        
        int indices = (*nodes)[k];
        double values  = i;
        locsystoggle_->ReplaceGlobalValues(1,&values,&indices);
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
      for (int k=0;k<(int)nodes->size();++k)
      {
        bool havenode = Discret().HaveGlobalNode((*nodes)[k]);
        if (!havenode) continue;
        
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

  // number of nodes subjected to local co-ordinate systems
  set<int> locsysdofset;
    
  trafo_ = rcp(new LINALG::SparseMatrix(*dofrowmap,3));
  
  for (int i=0;i<noderowmap->NumMyElements();++i)
  {
    int gid = noderowmap->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    vector<int> dofs = Discret().Dof(node);
    int numdof = (int)dofs.size();
    int locsysindex = (int)(*locsystoggle_)[i];
    
    // unity matrix for non-locsys node
    if (locsysindex<0)
    {
      for (int r=0;r<numdof;++r)
        for (int c=0;c<numdof;++c)
          if (r==c) trafo_->Assemble(1.0,dofs[r],dofs[c]);
    }
    
    // trafo matrix for locsys node
    else
    {
      Epetra_SerialDenseMatrix nodetrafo(numdof,numdof);
      
      // first create an identity matrix of size numdof
      for (int k=0;k<numdof;++k) nodetrafo(k,k)=1.0;
      
      // trafo for 2D and 3D case
      nodetrafo(0,0)=normals_(locsysindex,0);
      nodetrafo(0,1)=normals_(locsysindex,1);
      nodetrafo(1,0)=tangents_(locsysindex,0);
      nodetrafo(1,1)=tangents_(locsysindex,1);
      
      // additional trafo for 3D case
      if (Dim()==3)
      {
        nodetrafo(0,2)=normals_(locsysindex,2);
        nodetrafo(1,2)=tangents_(locsysindex,2);
        nodetrafo(2,0)=thirddir_(locsysindex,0);
        nodetrafo(2,1)=thirddir_(locsysindex,1);
        nodetrafo(2,2)=thirddir_(locsysindex,2);
      }
      
      for (int r=0;r<numdof;++r)
        for (int c=0;c<numdof;++c)
          trafo_->Assemble(nodetrafo(r,c),dofs[r],dofs[c]);

      // store the DOF with locsys
      if (transformleftonly_)
        for (int r=0;r<numdof;++r)
          locsysdofset.insert(dofs[r]);

    }
  }

  // complete transformation matrix
  trafo_->Complete();

  //**********************************************************************
  // Build map holding DOFs linked to nodes with local co-ordinate system
  //**********************************************************************

  // create unique/row map of DOFs subjected to local co-ordinate change
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
    locsysdofmap_ = rcp(new Epetra_Map(-1, nummyentries, myglobalentries,
                                       discret_.DofRowMap()->IndexBase(), 
                                       discret_.Comm()));
    if (locsysdofmap_ == null) dserror("Creation failed.");
  }

  // transformation matrix for relevent DOFs with local system
  if (transformleftonly_)
  {
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
    os << "\n***************************************************";
    os << "\n*                       WARNING                   *";
    os << "\n* LocsysManager is still in an experimental stage *";
    os << "\n* and not yet fully tested (e.g. no fluid test!)  *";
    os << "\n***************************************************\n";
    os << "\n-------------------------------------DRT::UTILS::LocysManager\n";
    for (int i=0;i<NumLocsys();++i)
    {
      os << "Locsys ID: " << i << "\t";
      if (TypeLocsys(i)==DRT::Condition::PointLocsys) os << "Point\t";
      else if (TypeLocsys(i)==DRT::Condition::LineLocsys) os << "Line\t";
      else if (TypeLocsys(i)==DRT::Condition::SurfaceLocsys) os << "Surface\t";
      else if (TypeLocsys(i)==DRT::Condition::VolumeLocsys) os << "Volume\t";
      else dserror("ERROR: Unknown type of locsys condition!");
      os << "Normal:  " << normals_(i,0) << " " << normals_(i,1) << " " << normals_(i,2) << "\t";
      os << "Tangent: " << tangents_(i,0) << " " << tangents_(i,1) << " " << tangents_(i,2) << "\n";
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
  // get dof row layout of discretization
  const Epetra_Map* dofrowmap = discret_.DofRowMap();
  
  // transform rhs vector
  if (transformleftonly_)
  {
    RotateGlobalToLocal(rhs);
  }
  else
  {
    RCP<Epetra_Vector> modrhs = LINALG::CreateVector(*dofrowmap);
    trafo_->Multiply(false,*rhs,*modrhs);
    rhs->Update(1.0,*modrhs,0.0);
  }
  
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
    temp = LINALG::Multiply(*trafo_,false,*sysmat,false,true);  // multiply from left
    RCP<LINALG::SparseMatrix> temp2;
    temp2 = LINALG::Multiply(*temp,false,*trafo_,true,true);  // multiply from right
    *sysmat = *temp2;
  }
  
  return;
}


/*----------------------------------------------------------------------*
 |  Transform vector global -> local (public)                 popp 09/08|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateGlobalToLocal(RCP<Epetra_Vector> vec)
{
  // get dof row layout of discretization
  const Epetra_Map* dofrowmap = discret_.DofRowMap();
      
  // transform vec
  RCP<Epetra_Vector> modvec = LINALG::CreateVector(*dofrowmap);
  trafo_->Multiply(false,*vec,*modvec);
  vec->Update(1.0,*modvec,0.0);
    
  return;
}


/*----------------------------------------------------------------------*
 |  Transform result + system local -> global (public)        popp 09/08|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateLocalToGlobal(RCP<Epetra_Vector> result,
                                               RCP<LINALG::SparseMatrix> sysmat,
                                               RCP<Epetra_Vector> rhs)
{
  // get dof row layout of discretization
  const Epetra_Map* dofrowmap = discret_.DofRowMap();
    
  // transform result
  RCP<Epetra_Vector> modres = LINALG::CreateVector(*dofrowmap);
  trafo_->Multiply(true,*result,*modres);
  result->Update(1.0,*modres,0.0);
  
  // transform rhs vector
  RCP<Epetra_Vector> modrhs = LINALG::CreateVector(*dofrowmap);
  trafo_->Multiply(true,*rhs,*modrhs);
  rhs->Update(1.0,*modrhs,0.0);
    
  // transform system matrix
  RCP<LINALG::SparseMatrix> temp;
  RCP<LINALG::SparseMatrix> temp2;
  temp = LINALG::Multiply(*sysmat,false,*trafo_,false,true);
  temp2 = LINALG::Multiply(*trafo_,true,*temp,false,true);
  *sysmat = *temp2;
    
  return;
}

/*----------------------------------------------------------------------*
 |  Transform vector local -> global (public)                 popp 09/08|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateLocalToGlobal(RCP<Epetra_Vector> vec)
{
  // get dof row layout of discretization
  const Epetra_Map* dofrowmap = discret_.DofRowMap();
      
  // transform vec
  RCP<Epetra_Vector> modvec = LINALG::CreateVector(*dofrowmap);
  trafo_->Multiply(true,*vec,*modvec);
  vec->Update(1.0,*modvec,0.0);
    
  return;
}

#endif  // #ifdef CCADISCRET
