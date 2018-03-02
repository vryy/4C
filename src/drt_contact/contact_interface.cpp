/*---------------------------------------------------------------------*/
/*!
\file contact_interface.cpp

\brief One contact interface

\level 2

\maintainer Alexander Seitz

*/
/*---------------------------------------------------------------------*/
#include <Epetra_CrsMatrix.h>
#include <Epetra_FEVector.h>
#include <Epetra_Time.h>
#include "contact_interface.H"
#include "contact_node.H"
#include "contact_element.H"
#include "contact_integrator.H"
#include "contact_interpolator.H"
#include "contact_coupling2d.H"
#include "contact_coupling3d.H"
#include "contact_defines.H"
#include "contact_line_coupling.H"
#include "friction_node.H"
#include "selfcontact_binarytree.H"
#include "../drt_mortar/mortar_binarytree.H"
#include "../drt_mortar/mortar_defines.H"
#include "../drt_mortar/mortar_projector.H"
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_lib/drt_utils_parmetis.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_adapter/adapter_coupling.H"

#include "../drt_mortar/mortar_coupling3d_classes.H"
#include "contact_nitsche_utils.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::IDataContainer::IDataContainer()
    : selfcontact_( false ),
      friction_( false ),
      nonSmoothContact_( false ),
      constr_direction_( INPAR::CONTACT::constr_vague ),
      activenodes_( Teuchos::null ),
      activedofs_( Teuchos::null ),
      inactivenodes_( Teuchos::null ),
      inactivedofs_( Teuchos::null ),
      activen_( Teuchos::null ),
      activet_( Teuchos::null ),
      slipnodes_( Teuchos::null ),
      slipdofs_( Teuchos::null ),
      slipt_( Teuchos::null ),
      nonsmoothnodes_( Teuchos::null ),
      smoothnodes_( Teuchos::null ),
      sdofVertexRowmap_( Teuchos::null ),
      sdofVertexColmap_( Teuchos::null ),
      sdofEdgeRowmap_( Teuchos::null ),
      sdofEdgeColmap_( Teuchos::null ),
      sdofSurfRowmap_( Teuchos::null ),
      sdofSurfColmap_( Teuchos::null ),
      nextendedghosting_( Teuchos::null ),
      eextendedghosting_( Teuchos::null ),
      binarytreeself_( Teuchos::null ),
      cnValues_( Teuchos::null ),
      ctValues_( Teuchos::null ),
      smpairs_( 0 ),
      smintpairs_( 0 ),
      intcells_( 0 )
{

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::CoInterface> CONTACT::CoInterface::Create(
      const int id, const Epetra_Comm& comm, const int dim,
      const Teuchos::ParameterList& icontact,
      const bool selfcontact, INPAR::MORTAR::RedundantStorage redundant)
{
  Teuchos::RCP<MORTAR::IDataContainer> idata_ptr =
      Teuchos::rcp( new CONTACT::IDataContainer() );
  return Teuchos::rcp( new CoInterface( idata_ptr, id, comm, dim,
      icontact, selfcontact, redundant ) );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::CoInterface::CoInterface(
    const Teuchos::RCP<CONTACT::IDataContainer>& idata_ptr )
    : MORTAR::MortarInterface( idata_ptr ),
      idata_ptr_( idata_ptr ),
      idata_( *idata_ptr_ ),
      selfcontact_( idata_.IsSelfContact() ),
      friction_( idata_.IsFriction() ),
      nonSmoothContact_( idata_.IsNonSmoothContact() ),
      constr_direction_( idata_.ConstraintDirection() ),
      activenodes_( idata_.ActiveNodes() ),
      activedofs_( idata_.ActiveDofs() ),
      inactivenodes_( idata_.InActiveNodes() ),
      inactivedofs_( idata_.InActiveDofs() ),
      activen_( idata_.ActiveN() ),
      activet_( idata_.ActiveT() ),
      slipnodes_( idata_.SlipNodes() ),
      slipdofs_( idata_.SlipDofs() ),
      slipt_( idata_.SlipT() ),
      nonsmoothnodes_( idata_.NonSmoothNodes() ),
      smoothnodes_( idata_.SmoothNodes() ),
      sdofVertexRowmap_( idata_.SdofVertexRowmap() ),
      sdofVertexColmap_( idata_.SdofVertexColmap() ),
      sdofEdgeRowmap_( idata_.SdofEdgeRowmap() ),
      sdofEdgeColmap_( idata_.SdofEdgeColmap() ),
      sdofSurfRowmap_( idata_.SdofSurfRowmap() ),
      sdofSurfColmap_( idata_.SdofSurfColmap() ),
      nextendedghosting_( idata_.NExtendedGhosting() ),
      eextendedghosting_( idata_.EExtendedGhosting() ),
      binarytreeself_( idata_.BinaryTreeSelf() ),
      cnValues_( idata_.CnValues() ),
      ctValues_( idata_.CtValues() ),
      smpairs_( idata_.SMIntPairs() ),
      smintpairs_( idata_.SMIntPairs() ),
      intcells_( idata_.IntCells() )
{
  /* do nothing */
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::CoInterface::CoInterface(
    const Teuchos::RCP<MORTAR::IDataContainer>& idata_ptr,
    const int id, const Epetra_Comm& comm,
    const int dim,
    const Teuchos::ParameterList& icontact,
    bool selfcontact,
    INPAR::MORTAR::RedundantStorage redundant)
    : MORTAR::MortarInterface(idata_ptr,id,comm,dim,icontact,redundant),
      idata_ptr_( Teuchos::rcp_dynamic_cast<CONTACT::IDataContainer>( idata_ptr, true ) ),
      idata_( *idata_ptr_ ),
      selfcontact_( idata_.IsSelfContact() ),
      friction_( idata_.IsFriction() ),
      nonSmoothContact_( idata_.IsNonSmoothContact() ),
      constr_direction_( idata_.ConstraintDirection() ),
      activenodes_( idata_.ActiveNodes() ),
      activedofs_( idata_.ActiveDofs() ),
      inactivenodes_( idata_.InActiveNodes() ),
      inactivedofs_( idata_.InActiveDofs() ),
      activen_( idata_.ActiveN() ),
      activet_( idata_.ActiveT() ),
      slipnodes_( idata_.SlipNodes() ),
      slipdofs_( idata_.SlipDofs() ),
      slipt_( idata_.SlipT() ),
      nonsmoothnodes_( idata_.NonSmoothNodes() ),
      smoothnodes_( idata_.SmoothNodes() ),
      sdofVertexRowmap_( idata_.SdofVertexRowmap() ),
      sdofVertexColmap_( idata_.SdofVertexColmap() ),
      sdofEdgeRowmap_( idata_.SdofEdgeRowmap() ),
      sdofEdgeColmap_( idata_.SdofEdgeColmap() ),
      sdofSurfRowmap_( idata_.SdofSurfRowmap() ),
      sdofSurfColmap_( idata_.SdofSurfColmap() ),
      nextendedghosting_( idata_.NExtendedGhosting() ),
      eextendedghosting_( idata_.EExtendedGhosting() ),
      binarytreeself_( idata_.BinaryTreeSelf() ),
      cnValues_( idata_.CnValues() ),
      ctValues_( idata_.CtValues() ),
      smpairs_( idata_.SMIntPairs() ),
      smintpairs_( idata_.SMIntPairs() ),
      intcells_( idata_.IntCells() )
{
  selfcontact_ = selfcontact;
  nonSmoothContact_ = DRT::INPUT::IntegralValue<int>(icontact,"NONSMOOTH_GEOMETRIES");
  constr_direction_ = DRT::INPUT::IntegralValue<INPAR::CONTACT::ConstraintDirection>(icontact,"CONSTRAINT_DIRECTIONS");
  smpairs_ = 0;
  smintpairs_ = 0;
  intcells_ = 0;

  // set frictional contact status
  INPAR::CONTACT::FrictionType ftype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(icontact,"FRICTION");
  if (ftype != INPAR::CONTACT::friction_none)
    friction_ = true;

  // set poro contact
  if (icontact.get<int>("PROBTYPE")==INPAR::CONTACT::poro)
    SetPoroFlag(true);

  // set ehl contact
  if (icontact.get<int>("PROBTYPE")==INPAR::CONTACT::ehl)
    SetEhlFlag(true);

  // check for redundant slave storage
  // (needed for self contact but not wanted for general contact)
//  if ((selfcontact_ or nonSmoothContact_) && redundant != INPAR::MORTAR::redundant_all)
//    dserror("ERROR: We need redundant interface storage (slave+master) for self contact");
  if (!(selfcontact_ or nonSmoothContact_) && redundant == INPAR::MORTAR::redundant_all)
    dserror("ERROR: We do not want redundant interface storage for contact");

  // init extended ghosting for RR loop
  eextendedghosting_ = Teuchos::null;
  nextendedghosting_ = Teuchos::null;

  return;
}

/*----------------------------------------------------------------------*
 |  update master and slave sets (nodes etc.)                farah 10/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::UpdateMasterSlaveSets()
{
  // call mortar function
  MORTAR::MortarInterface::UpdateMasterSlaveSets();

  //********************************************************************
  // DOFS
  //********************************************************************
  // do the same business for dofs
  // (get row and column maps of slave and master dofs seperately)
  if(nonSmoothContact_)
  {
    std::vector<int> sVc; // slave column map
    std::vector<int> sVr; // slave row map
    std::vector<int> sEc; // master column map
    std::vector<int> sEr; // master row map
    std::vector<int> sSc; // master column map
    std::vector<int> sSr; // master row map

    for (int i = 0; i < Discret().NodeColMap()->NumMyElements(); ++i)
    {
      int gid = Discret().NodeColMap()->GID(i);
      DRT::Node* node = Discret().gNode(gid);
      if (!node)
        dserror("ERROR: Cannot find node with gid %", gid);
      CoNode* mrtrnode = dynamic_cast<CoNode*>(node);
      bool isslave = mrtrnode->IsSlave();

      if(isslave)
      {
        // vertex
        if(mrtrnode->IsOnCorner())
        {
          for (int j = 0; j < mrtrnode->NumDof(); ++j)
            sVc.push_back(mrtrnode->Dofs()[j]);

          if (Discret().NodeRowMap()->MyGID(gid))
            for (int j = 0; j < mrtrnode->NumDof(); ++j)
              sVr.push_back(mrtrnode->Dofs()[j]);
        }
        // edge
        else if(mrtrnode->IsOnEdge())
        {
          for (int j = 0; j < mrtrnode->NumDof(); ++j)
            sEc.push_back(mrtrnode->Dofs()[j]);

          if (Discret().NodeRowMap()->MyGID(gid))
            for (int j = 0; j < mrtrnode->NumDof(); ++j)
              sEr.push_back(mrtrnode->Dofs()[j]);
        }
        // surface
        else if(!mrtrnode->IsOnCornerEdge())
        {
          for (int j = 0; j < mrtrnode->NumDof(); ++j)
            sSc.push_back(mrtrnode->Dofs()[j]);

          if (Discret().NodeRowMap()->MyGID(gid))
            for (int j = 0; j < mrtrnode->NumDof(); ++j)
              sSr.push_back(mrtrnode->Dofs()[j]);
        }
        else
        {
          dserror("ERROR: unknown case!");
        }
      }
    }

    sdofVertexRowmap_ = Teuchos::rcp(
        new Epetra_Map(-1, (int) sVr.size(), &sVr[0], 0, Comm()));
    sdofVertexColmap_ = Teuchos::rcp(
        new Epetra_Map(-1, (int) sVc.size(), &sVc[0], 0, Comm()));
    sdofEdgeRowmap_ = Teuchos::rcp(
        new Epetra_Map(-1, (int) sEr.size(), &sEr[0], 0, Comm()));
    sdofEdgeColmap_ = Teuchos::rcp(
        new Epetra_Map(-1, (int) sEc.size(), &sEc[0], 0, Comm()));
    sdofSurfRowmap_ = Teuchos::rcp(
        new Epetra_Map(-1, (int) sSr.size(), &sSr[0], 0, Comm()));
    sdofSurfColmap_ = Teuchos::rcp(
        new Epetra_Map(-1, (int) sSc.size(), &sSc[0], 0, Comm()));
  }

  return;
}

/*----------------------------------------------------------------------*
 |  create and fill cn vector                                farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::SetCnCtValues(const int& iter)
{
  // get cn from the input file
  const double cn = IParams().get<double>("SEMI_SMOOTH_CN");
  const double ct = IParams().get<double>("SEMI_SMOOTH_CT");

  // set all nodal cn-values to the input value
  GetCn() = LINALG::CreateVector(*SlaveRowNodes(),true);
  int err = GetCn()->PutScalar(cn);
  if(err!=0)
    dserror("ERROR: cn definition failed!");

  // set all nodal ct-values to the input value
  if(friction_)
  {
    GetCt() = LINALG::CreateVector(*SlaveRowNodes(),true);
    err = GetCt()->PutScalar(ct);
    if(err!=0)
      dserror("ERROR: cn definition failed!");
  }

  // modification for edge/corner nodes
  for(int i = 0; i<SlaveRowNodes()->NumMyElements();++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %i",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    // calculate characteristic edge length:
    DRT::Element* ele = cnode->Elements()[0];
    CoElement* cele = dynamic_cast<CoElement*>(ele);
    double pos1[3] = {0.0,0.0,0.0};
    double pos2[3] = {0.0,0.0,0.0};
    double vec[3]  = {0.0,0.0,0.0};

    pos1[0] = dynamic_cast<CoNode*>(cele->Nodes()[0])->X()[0];
    pos1[1] = dynamic_cast<CoNode*>(cele->Nodes()[0])->X()[1];
    pos1[2] = dynamic_cast<CoNode*>(cele->Nodes()[0])->X()[2];

    pos2[0] = dynamic_cast<CoNode*>(cele->Nodes()[1])->X()[0];
    pos2[1] = dynamic_cast<CoNode*>(cele->Nodes()[1])->X()[1];
    pos2[2] = dynamic_cast<CoNode*>(cele->Nodes()[1])->X()[2];

    vec[0] = pos1[0]-pos2[0];
    vec[1] = pos1[1]-pos2[1];
    vec[2] = pos1[2]-pos2[2];

    const double length = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if(length<1e-12)
    {
      std::cout << "*** WARNING *** element edge nearly zero" << std::endl;
      continue;
    }

    if(cnode->IsOnEdge())
    {
      GetCnRef()[GetCnRef().Map().LID(cnode->Id())] = cn * (length*length);
      if(friction_)
        GetCtRef()[GetCtRef().Map().LID(cnode->Id())] = ct * (length*length);
    }

    if(cnode->IsOnCorner())
    {
      GetCnRef()[GetCnRef().Map().LID(cnode->Id())] = cn * (length*length*length*length);
      if(friction_)
        GetCtRef()[GetCtRef().Map().LID(cnode->Id())] = ct * (length*length*length*length);
    }
  }


  return;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
std::ostream& operator << (std::ostream& os, const CONTACT::CoInterface& interface)
{
  interface.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print interface (public)                                 mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::Print(std::ostream& os) const
{
  if (Comm().MyPID()==0)
    os << "Contact ";
  MORTAR::MortarInterface::Print(os);

  return;
}

/*----------------------------------------------------------------------*
 |  add contact node (public)                                mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AddCoNode(Teuchos::RCP<CONTACT::CoNode> cnode)
{
  idiscret_->AddNode(cnode);
  return;
}

/*----------------------------------------------------------------------*
 |  add contact element (public)                             mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AddCoElement(Teuchos::RCP<CONTACT::CoElement> cele)
{
  // check for quadratic 2d slave elements to be modified
  if (cele->IsSlave() && (cele->Shape()==DRT::Element::line3))
    quadslave_=true;

  // check for quadratic 3d slave elements to be modified
  if (cele->IsSlave() && (cele->Shape()==DRT::Element::quad8 || cele->Shape()==DRT::Element::tri6))
    quadslave_=true;

  idiscret_->AddElement(cele);
  return;
}


/*----------------------------------------------------------------------*
 |  store the required ghosting within a round              farah 10/13 |
 |  robin iteration for current interface (public)                      |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::RoundRobinExtendGhosting(bool firstevaluation)
{
  std::vector<int> eghosting;
  std::vector<int> nghosting;
  for (int k=0; k<SlaveColElements()->NumMyElements();++k)
  {
    int gid = SlaveColElements()->GID(k);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find ele with gid %i",gid);
    CoElement* sele = dynamic_cast<CoElement*>(ele);

    for (int j=0;j<sele->MoData().NumSearchElements();++j)
    {
      int gid2 = sele->MoData().SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2) dserror("ERROR: Cannot find master element with gid %",gid2);
      CoElement* melement = dynamic_cast<CoElement*>(ele2);

      eghosting.push_back(melement->Id());

      for (int z=0;z<melement->NumNode();++z)
      {
        int gidn = melement->NodeIds()[z];
        nghosting.push_back(gidn);
      }
    }
    //reset found elements
    sele->DeleteSearchElements();
  }

  Teuchos::RCP<Epetra_Map> ecurrghosting = Teuchos::rcp(new Epetra_Map(-1,(int)eghosting.size(),&eghosting[0],0,Comm()));
  Teuchos::RCP<Epetra_Map> ncurrghosting = Teuchos::rcp(new Epetra_Map(-1,(int)nghosting.size(),&nghosting[0],0,Comm()));

  if(firstevaluation)
  {
    eextendedghosting_=Teuchos::rcp(new Epetra_Map(*ecurrghosting));
    nextendedghosting_=Teuchos::rcp(new Epetra_Map(*ncurrghosting));
  }
  else
  {
    eextendedghosting_ = LINALG::MergeMap(eextendedghosting_,ecurrghosting,true);
    nextendedghosting_ = LINALG::MergeMap(nextendedghosting_,ncurrghosting,true);
  }

  return;
}

/*----------------------------------------------------------------------*
 | perform the ownership change within a round robin        farah 10/13 |
 | iteration                                                            |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::RoundRobinChangeOwnership()
{
  // pack/unpack friction nodes is only required for wear problems
  // so we should create a slightly redundant function for the wear-
  // interface and exclude the friction node packing from here

  // get friction type
  INPAR::CONTACT::FrictionType ftype =
    DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(IParams(),"FRICTION");

  // change master-side proc ownership
  // some local variables
  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(Comm().Clone());
  const int myrank  = comm->MyPID();
  const int numproc = comm->NumProc();
  const int torank = (myrank + 1) % numproc;             // to
  const int fromrank = (myrank + numproc - 1) % numproc; // from

  // new eles, new nodes
  std::vector<int> ncol, nrow;
  std::vector<int> ecol, erow;

  //create dummy
  Teuchos::RCP<Epetra_Map> MasterColNodesdummy = Teuchos::rcp(new Epetra_Map(*MasterColNodes()));
  Teuchos::RCP<Epetra_Map> MasterColelesdummy  = Teuchos::rcp(new Epetra_Map(*MasterColElements()));

  //create origin maps
  Teuchos::RCP<Epetra_Map> SCN  = Teuchos::rcp(new Epetra_Map(*SlaveColNodes()));
  Teuchos::RCP<Epetra_Map> SCE  = Teuchos::rcp(new Epetra_Map(*SlaveColElements()));
  Teuchos::RCP<Epetra_Map> SRN  = Teuchos::rcp(new Epetra_Map(*SlaveRowNodes()));
  Teuchos::RCP<Epetra_Map> SRE  = Teuchos::rcp(new Epetra_Map(*SlaveRowElements()));

  // *****************************************
  // Elements
  // *****************************************
  std::vector<char> sdataeles;
  std::vector<char> rdataeles;

  // vector containing all proc ids
  std::vector<int> allproc(numproc);
  for (int i=0; i<numproc; ++i) allproc[i] = i;

  // get exporter
  DRT::Exporter exporter(idiscret_->Comm());

  // create data buffer
  DRT::PackBuffer dataeles;

  // pack data - first just reserving the mem.
  for (int i = 0; i < (int)MasterColelesdummy->NumMyElements(); ++i)
  {
    int gid = MasterColelesdummy->GID(i);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find ele with gid %i",gid);
    MORTAR::MortarElement* mele = dynamic_cast<MORTAR::MortarElement*>(ele);

    mele->Pack(dataeles);

    int ghost=0;
    DRT::ParObject::AddtoPack(dataeles,ghost);
  }

  dataeles.StartPacking();

  // now pack/store
  for (int i = 0; i < (int)MasterColelesdummy->NumMyElements(); ++i)
  {
    int gid = MasterColelesdummy->GID(i);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find ele with gid %i",gid);
    MORTAR::MortarElement* mele = dynamic_cast<MORTAR::MortarElement*>(ele);

    mele->Pack(dataeles);

    // check for ghosting
    int ghost;
    if (mele->Owner()==myrank) ghost=1;
    else                       ghost=0;

    DRT::ParObject::AddtoPack(dataeles,ghost);
  }
  std::swap(sdataeles, dataeles());

  // delete the elements from discretization
  for (int i = 0; i < (int)MasterColelesdummy->NumMyElements(); ++i)
  {
    int gid = MasterColelesdummy->GID(i);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find ele with gid %i",gid);
    MORTAR::MortarElement* mele = dynamic_cast<MORTAR::MortarElement*>(ele);

    // check for ghosting
    if (mele->Owner()==myrank)
    {
      idiscret_->DeleteElement(mele->Id());
    }
  }

  // send the information
  MPI_Request request;
  exporter.ISend(myrank, torank, &(sdataeles[0]), (int)sdataeles.size(), 1234, request);

  // receive the information
  int length = rdataeles.size();
  int tag = -1;
  int from = -1;
  exporter.ReceiveAny(from,tag,rdataeles,length);
  if (tag != 1234 or from != fromrank)
    dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank, myrank, from, myrank);

  // ---- unpack ----
  {
    // Put received nodes into discretization
    std::vector<char>::size_type index = 0;
    while (index < rdataeles.size())
    {
      std::vector<char> data;
      int ghost=-1;
      DRT::ParObject::ExtractfromPack(index,rdataeles,data);
      DRT::ParObject::ExtractfromPack(index,rdataeles,ghost);
      if (ghost==-1) dserror("UNPACK ERROR!!!!!!!!!");

      // this Teuchos::rcp holds the memory of the ele
      Teuchos::RCP<DRT::ParObject> object = Teuchos::rcp(DRT::UTILS::Factory(data),true);
      Teuchos::RCP<MORTAR::MortarElement> ele = Teuchos::rcp_dynamic_cast<MORTAR::MortarElement>(object);
      if (ele == Teuchos::null) dserror("Received object is not an ele");

      // add whether its a row ele
      if(ghost==1)
      {
        ele->SetOwner(myrank);
        idiscret_->AddElement(ele);

        // to new ele
        erow.push_back(ele->Id());
        ecol.push_back(ele->Id());
      }
      else
      {
        ecol.push_back(ele->Id());
      }
    }
  }// end unpack

  // wait for all communication to finish
  exporter.Wait(request);
  comm->Barrier();

  // *****************************************
  // NODES
  // *****************************************
  std::vector<char> sdatanodes;
  std::vector<char> rdatanodes;

  // get exporter
  DRT::Exporter exportern(idiscret_->Comm());

  DRT::PackBuffer datanodes;

  // pack data -- col map --> should prevent further ghosting!
  for (int i = 0; i < (int)MasterColNodesdummy->NumMyElements(); ++i)
  {
    int gid = MasterColNodesdummy->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find ele with gid %i",gid);

    if (ftype==INPAR::CONTACT::friction_none)
    {
      MORTAR::MortarNode* cnode = dynamic_cast<MORTAR::MortarNode*>(node);
      cnode->Pack(datanodes);
    }
    else
    {
      FriNode* cnode = dynamic_cast<FriNode*>(node);
      cnode->Pack(datanodes);
    }

    int ghost=0;
    DRT::ParObject::AddtoPack(datanodes,ghost);
  }

  datanodes.StartPacking();
  for (int i = 0; i < (int)MasterColNodesdummy->NumMyElements(); ++i)
  {
    int gid = MasterColNodesdummy->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find ele with gid %i",gid);

    // check for ghosting
    int ghost;

    if (ftype==INPAR::CONTACT::friction_none)
    {
      MORTAR::MortarNode* cnode = dynamic_cast<MORTAR::MortarNode*>(node);
      cnode->Pack(datanodes);

      if (cnode->Owner()==myrank) ghost=1;
      else                        ghost=0;
    }
    else
    {
      FriNode* cnode = dynamic_cast<FriNode*>(node);
      cnode->Pack(datanodes);

      if (cnode->Owner()==myrank) ghost=1;
      else                        ghost=0;
    }

    DRT::ParObject::AddtoPack(datanodes,ghost);
  }
  std::swap(sdatanodes, datanodes());

  //DELETING
  for (int i = 0; i < (int)MasterColNodesdummy->NumMyElements(); ++i)
  {
    int gid = MasterColNodesdummy->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find ele with gid %i",gid);

    if (ftype==INPAR::CONTACT::friction_none)
    {
      MORTAR::MortarNode* cnode = dynamic_cast<MORTAR::MortarNode*>(node);
      if (cnode->Owner()==myrank)
        idiscret_->DeleteNode(cnode->Id());
    }
    else
    {
      FriNode* cnode = dynamic_cast<FriNode*>(node);
      if (cnode->Owner()==myrank)
        idiscret_->DeleteNode(cnode->Id());
    }
  }

  // ---- send ----
  MPI_Request requestn;
  exportern.ISend(myrank, torank, &(sdatanodes[0]), (int)sdatanodes.size(), 1234, requestn);

  // ---- receive ----
  int lengthn = rdatanodes.size();
  int tagn = -1;
  int fromn = -1;
  exportern.ReceiveAny(fromn,tagn,rdatanodes,lengthn);
  if (tagn != 1234 or fromn != fromrank)
    dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank, myrank, fromn, myrank);

  // ---- unpack ----
  {
    // Put received nodes into discretization
    std::vector<char>::size_type index = 0;
    while (index < rdatanodes.size())
    {
      std::vector<char> data;

      int ghost=-1;
      DRT::ParObject::ExtractfromPack(index,rdatanodes,data);
      DRT::ParObject::ExtractfromPack(index,rdatanodes,ghost);
      if (ghost==-1) dserror("UNPACK ERROR!!!!!!!!!");

      // this Teuchos::rcp holds the memory of the node
      Teuchos::RCP<DRT::ParObject> object = Teuchos::rcp(DRT::UTILS::Factory(data),true);

      if (ftype==INPAR::CONTACT::friction_none)
      {
        Teuchos::RCP<MORTAR::MortarNode> node = Teuchos::rcp_dynamic_cast<MORTAR::MortarNode>(object);
        if (node == Teuchos::null) dserror("Received object is not a node");

        if (ghost==1)
        {
          node->SetOwner(myrank);
          idiscret_->AddNode(node);

          nrow.push_back(node->Id());
          ncol.push_back(node->Id());
        }
        else
        {
          // all others to col
          ncol.push_back(node->Id());
        }
      }
      else // if friction...
      {
        Teuchos::RCP<FriNode> node = Teuchos::rcp_dynamic_cast<FriNode>(object);
        if (node == Teuchos::null) dserror("Received object is not a node");

        if (ghost==1)
        {
          node->SetOwner(myrank);
          idiscret_->AddNode(node);

          nrow.push_back(node->Id());
          ncol.push_back(node->Id());
        }
        else
        {
          // all others to col
          ncol.push_back(node->Id());
        }
      }
    }
  } // end unpack

  // wait for all communication to finish
  exportern.Wait(requestn);
  comm->Barrier();

  //create maps from sending
  Teuchos::RCP<Epetra_Map> noderowmap = Teuchos::rcp(new Epetra_Map(-1,(int)nrow.size(),&nrow[0],0,Comm()));
  Teuchos::RCP<Epetra_Map> nodecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)ncol.size(),&ncol[0],0,Comm()));

  Teuchos::RCP<Epetra_Map> elerowmap = Teuchos::rcp(new Epetra_Map(-1,(int)erow.size(),&erow[0],0,Comm()));
  Teuchos::RCP<Epetra_Map> elecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)ecol.size(),&ecol[0],0,Comm()));

  // Merge s/m column maps for eles and nodes
  Teuchos::RCP<Epetra_Map> colnodesfull = LINALG::MergeMap(nodecolmap,SCN,true);
  Teuchos::RCP<Epetra_Map> colelesfull  = LINALG::MergeMap(elecolmap,SCE,true);

  // Merge s/m row maps for eles and nodes
  Teuchos::RCP<Epetra_Map> rownodesfull = LINALG::MergeMap(noderowmap,SRN,false);
  Teuchos::RCP<Epetra_Map> rowelesfull  = LINALG::MergeMap(elerowmap,SRE,false);

  // to discretization
  // export nodes and elements to the row map
  Discret().ExportRowNodes(*rownodesfull);
  Discret().ExportRowElements(*rowelesfull);

  // export nodes and elements to the col map
  Discret().ExportColumnNodes(*colnodesfull);
  Discret().ExportColumnElements(*colelesfull);

  // ********************************************
  // call the (very) expensive FILLCOMPLETE()!
  // ********************************************
  // make sure discretization is complete
  FillComplete();

  return;
}


/*----------------------------------------------------------------------*
 |  change master ownership clockwise for contact            farah 10/13|
 |  interface without evaluation of the interface                       |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::RoundRobinDetectGhosting()
{
  if (SearchAlg()==INPAR::MORTAR::search_bfele)           EvaluateSearchBruteForce(SearchParam());
  else if (SearchAlg()==INPAR::MORTAR::search_binarytree) EvaluateSearchBinarytree();
  else                                                    dserror("ERROR: Invalid search algorithm");

  // first ghosting for std. distribution
  RoundRobinExtendGhosting(true);

  // Init Maps
  Teuchos::RCP<Epetra_Map> Init_SCN  = Teuchos::rcp(new Epetra_Map(*SlaveColNodes()));
  Teuchos::RCP<Epetra_Map> Init_SCE  = Teuchos::rcp(new Epetra_Map(*SlaveColElements()));
  Teuchos::RCP<Epetra_Map> Init_MCN  = Teuchos::rcp(new Epetra_Map(*MasterColNodes()));
  Teuchos::RCP<Epetra_Map> Init_MCE  = Teuchos::rcp(new Epetra_Map(*MasterColElements()));

  // *************************************
  // start RR loop for current interface
  // *************************************
  // loop over all procs
  if (Comm().NumProc()>1)
    for (int j=0;j<(int)(Comm().NumProc());++j)
    {
      // status output
      if (Comm().MyPID()==0 && j==0) std::cout << "Round-Robin-Iteration #" << j;
      if (Comm().MyPID()==0 && j>0) std::cout << " #" << j;

      // perform the ownership change
      RoundRobinChangeOwnership();

      //build new search tree or do nothing for bruteforce
      if (SearchAlg()==INPAR::MORTAR::search_binarytree)   CreateSearchTree();
      else if (SearchAlg()!=INPAR::MORTAR::search_bfele)   dserror("ERROR: Invalid search algorithm");

      // evaluate interfaces
      if (j<(int)(Comm().NumProc()-1))
      {
        if (SearchAlg()==INPAR::MORTAR::search_bfele)           EvaluateSearchBruteForce(SearchParam());
        else if (SearchAlg()==INPAR::MORTAR::search_binarytree) EvaluateSearchBinarytree();
        else                                                    dserror("ERROR: Invalid search algorithm");

        // other ghostings per iter
        RoundRobinExtendGhosting(false);
      }
      else
      {
        // do nothing -- just switch
      }
    }//end RR

  //NEW VERSION
  eextendedghosting_ = LINALG::MergeMap(eextendedghosting_,Init_SCE,true);
  nextendedghosting_ = LINALG::MergeMap(nextendedghosting_,Init_SCN,true);
  eextendedghosting_ = LINALG::MergeMap(eextendedghosting_,Init_MCE,true);
  nextendedghosting_ = LINALG::MergeMap(nextendedghosting_,Init_MCN,true);

  //finally extend ghosting
  Discret().ExportColumnElements(*eextendedghosting_);
  Discret().ExportColumnNodes(*nextendedghosting_);
  FillComplete();

  // reset extended ghosting maps
  eextendedghosting_ = Teuchos::null;
  nextendedghosting_ = Teuchos::null;

  //build new search tree or do nothing for bruteforce
  if (SearchAlg()==INPAR::MORTAR::search_binarytree)   CreateSearchTree();
  else if (SearchAlg()!=INPAR::MORTAR::search_bfele)   dserror("ERROR: Invalid search algorithm");

  // final output for loop
  if (Comm().MyPID()==0)
    std::cout << " RRL done!" << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 |  change master ownership clockwise for contact            farah 10/13|
 |  interface and evaluate within each iteration                        |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::RoundRobinEvaluate()
{
  // first evaluation with init. parallel redistr.
  Evaluate(0);

  // first ghosting for std. distribution
  RoundRobinExtendGhosting(true);

  // Init Maps
  // TODO: use initial_ele_colmap from abstract_strategy.H !!!
  Teuchos::RCP<Epetra_Map> Init_SCN  = Teuchos::rcp(new Epetra_Map(*SlaveColNodes()));
  Teuchos::RCP<Epetra_Map> Init_SCE  = Teuchos::rcp(new Epetra_Map(*SlaveColElements()));
  Teuchos::RCP<Epetra_Map> Init_MCN  = Teuchos::rcp(new Epetra_Map(*MasterColNodes()));
  Teuchos::RCP<Epetra_Map> Init_MCE  = Teuchos::rcp(new Epetra_Map(*MasterColElements()));

  // *************************************
  // start RR loop for current interface
  // *************************************
  // loop over all procs
  if (Comm().NumProc()>1)
    for (int j=0;j<(int)(Comm().NumProc());++j)
    {
      // status output
      if (Comm().MyPID()==0 && j==0) std::cout << "Round-Robin-Iteration #" << j;
      if (Comm().MyPID()==0 && j>0) std::cout << " #" << j;

      // perform the ownership change
      RoundRobinChangeOwnership();

      //build new search tree or do nothing for bruteforce
      if (SearchAlg()==INPAR::MORTAR::search_binarytree)   CreateSearchTree();
      else if (SearchAlg()!=INPAR::MORTAR::search_bfele)   dserror("ERROR: Invalid search algorithm");

      // evaluate interfaces
      if (j<(int)(Comm().NumProc()-1))
      {
        Evaluate(j+1);

        // other ghostings per iter
        RoundRobinExtendGhosting(false);
      }
      else
      {
        // do nothing -- just switch
      }
    }//end RR

  eextendedghosting_ = LINALG::MergeMap(eextendedghosting_,Init_SCE,true);
  nextendedghosting_ = LINALG::MergeMap(nextendedghosting_,Init_SCN,true);
  eextendedghosting_ = LINALG::MergeMap(eextendedghosting_,Init_MCE,true);
  nextendedghosting_ = LINALG::MergeMap(nextendedghosting_,Init_MCN,true);

  //finally extend ghosting
  Discret().ExportColumnElements(*eextendedghosting_);
  Discret().ExportColumnNodes(*nextendedghosting_);
  FillComplete();

  // reset extended ghosting maps
  eextendedghosting_ = Teuchos::null;
  nextendedghosting_ = Teuchos::null;

  // final output for loop
  if (Comm().MyPID()==0)
    std::cout << " RRL done!" << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 |  redistribute contact interface (public)                   popp 08/10|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::Redistribute(int index)
{
  // make sure we are supposed to be here
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(IParams(),"PARALLEL_REDIST")==INPAR::MORTAR::parredist_none)
    dserror("ERROR: You are not supposed to be here...");

  // some local variables
  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(Comm().Clone());
  const int myrank  = comm->MyPID();
  const int numproc = comm->NumProc();
  Epetra_Time time(*comm);
  std::set<int>::const_iterator iter;

  // vector containing all proc ids
  std::vector<int> allproc(numproc);
  for (int i=0; i<numproc; ++i) allproc[i] = i;

  //**********************************************************************
  // (1) SLAVE splitting in close / non-close parts
  //**********************************************************************
  // perform contact search (still with non-optimal distribution)
  Initialize();
  if (SearchAlg()==INPAR::MORTAR::search_bfele)           EvaluateSearchBruteForce(SearchParam());
  else if (SearchAlg()==INPAR::MORTAR::search_binarytree) EvaluateSearchBinarytree();
  else                                                    dserror("ERROR: Invalid search algorithm");

  // split slave element row map and build redundant vector of
  // all close / non-close slave node ids on all procs
  std::vector<int> closeele, noncloseele;
  std::vector<int> localcns, localfns;

  // loop over all row elements to gather the local information
  for (int i=0; i<SlaveRowElements()->NumMyElements(); ++i)
  {
    // get element
    int gid = SlaveRowElements()->GID(i);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find element with gid %",gid);
    MORTAR::MortarElement* cele = dynamic_cast<MORTAR::MortarElement*>(ele);

    // store element id and adjacent node ids
    int close = cele->MoData().NumSearchElements();
    if (close > 0)
    {
      closeele.push_back(gid);
      for (int k=0;k<cele->NumNode();++k) localcns.push_back(cele->NodeIds()[k]);
    }
    else
    {
      noncloseele.push_back(gid);
      for (int k=0;k<cele->NumNode();++k) localfns.push_back(cele->NodeIds()[k]);
    }
  }

  // loop over all elements to reset candidates / search lists
  // (use standard slave column map)
  for (int i=0;i<SlaveColElements()->NumMyElements();++i)
  {
    int gid = SlaveColElements()->GID(i);
    DRT::Element* ele = Discret().gElement(gid);
    if (!ele) dserror("ERROR: Cannot find ele with gid %i",gid);
    MORTAR::MortarElement* mele = dynamic_cast<MORTAR::MortarElement*>(ele);

    mele->MoData().SearchElements().resize(0);
  }

  // we need an arbitrary preliminary element row map
  Teuchos::RCP<Epetra_Map> scroweles  = Teuchos::rcp(new Epetra_Map(-1,(int)closeele.size(),&closeele[0],0,Comm()));
  Teuchos::RCP<Epetra_Map> sncroweles = Teuchos::rcp(new Epetra_Map(-1,(int)noncloseele.size(),&noncloseele[0],0,Comm()));
  Teuchos::RCP<Epetra_Map> mroweles   = Teuchos::rcp(new Epetra_Map(*MasterRowElements()));

  // check for consistency
  if (scroweles->NumGlobalElements()==0 && sncroweles->NumGlobalElements()==0)
    dserror("ERROR: Redistribute: Both slave sets (close/non-close) are empty");

  //**********************************************************************
  // (2) SPECIAL CASES and output to screen
  //**********************************************************************
  // print element overview
  if (!myrank)
  {
    int cl = scroweles->NumGlobalElements();
    int ncl = sncroweles->NumGlobalElements();
     int ma = mroweles->NumGlobalElements();
    std::cout << "Element overview: " << cl << " / " << ncl << " / " << ma << "  (close-S / non-close-S / M)";
  }

  // print old parallel distribution
  PrintParallelDistribution(index);

  // use simple base class method if there are ONLY close or non-close elements
  // (return value TRUE, because redistribution performed)
  if (scroweles->NumGlobalElements()==0 || sncroweles->NumGlobalElements()==0)
  {
    MORTAR::MortarInterface::Redistribute();
    return true;
  }

  //**********************************************************************
  // (3a) PREPARATIONS decide how many procs are used
  //**********************************************************************
  // first we assume that all procs will be used
  int scproc = numproc;
  int sncproc = numproc;
  int mproc = numproc;

  // minimum number of elements per proc
  int minele = IParams().get<int>("MIN_ELEPROC");

  // calculate real number of procs to be used
  if (minele > 0)
  {
    scproc  = static_cast<int>((scroweles->NumGlobalElements()) / minele);
    sncproc = static_cast<int>((sncroweles->NumGlobalElements()) / minele);
    mproc   = static_cast<int>((mroweles->NumGlobalElements()) / minele);
    if (scroweles->NumGlobalElements() < 2*minele)  scproc = 1;
    if (sncroweles->NumGlobalElements() < 2*minele) sncproc = 1;
    if (mroweles->NumGlobalElements() < 2*minele)   mproc = 1;
    if (scproc > numproc)  scproc = numproc;
    if (sncproc > numproc) sncproc = numproc;
    if (mproc > numproc)   mproc = numproc;
  }

  // print message
  if (!myrank)
  {
    std::cout << "\nRedistributing interface using ZOLTAN.........done!\n";
    std::cout << "Procs used for redistribution: " << scproc << " / " << sncproc << " / " << mproc << " (close-S / non-close-S / M)\n";
  }

  //**********************************************************************
  // (3b) PREPARATIONS build initial node graph
  //**********************************************************************
  // create graph object
  Teuchos::RCP<Epetra_CrsGraph> graph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*SlaveRowNodes(),108,false));

  // loop over all row nodes to fill graph
  for (int k=0;k<SlaveRowNodes()->NumMyElements();++k)
  {
    int gid = SlaveRowNodes()->GID(k);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);

    // find adjacent elements first
    for (int k=0;k<node->NumElement();++k)
    {
      // store adjacent nodes
      DRT::Element* ele = node->Elements()[k];
      int numnode = ele->NumNode();
      std::vector<int> nodeids(numnode);
      for (int n=0;n<numnode;++n) nodeids[n] = ele->NodeIds()[n];

      int err = graph->InsertGlobalIndices(gid,numnode,&nodeids[0]);
      if (err<0) dserror("graph->InsertGlobalIndices returned %d",err);
      if (err==1) dserror("graph->InsertGlobalIndices returned %d",err);
    }
  }

  // fill graph and optimize storage
  graph->FillComplete();
  graph->OptimizeStorage();

  //**********************************************************************
  // (4) CLOSE SLAVE redistribution
  //**********************************************************************
  Teuchos::RCP<Epetra_Map> scrownodes = Teuchos::null;
  Teuchos::RCP<Epetra_Map> sccolnodes = Teuchos::null;

  // build redundant vector of all close slave node ids on all procs
  // (there must not be any double entries in the node lists, thus
  // transform to sets and then back to vectors)
  std::vector<int> globalcns;
  std::set<int> setglobalcns;
  std::vector<int> scnids;
  LINALG::Gather<int>(localcns,globalcns,numproc,&allproc[0],Comm());
  for (int i=0;i<(int)globalcns.size();++i) setglobalcns.insert(globalcns[i]);
  for (iter=setglobalcns.begin();iter!=setglobalcns.end();++iter) scnids.push_back(*iter);

  //**********************************************************************
  // call ZOLTAN for parallel redistribution
  DRT::UTILS::PartUsingParMetis(idiscret_,scroweles,scrownodes,sccolnodes,comm,false,scproc);
  //**********************************************************************

  //**********************************************************************
  // (5) NON-CLOSE SLAVE redistribution
  //**********************************************************************
  Teuchos::RCP<Epetra_Map> sncrownodes = Teuchos::null;
  Teuchos::RCP<Epetra_Map> snccolnodes = Teuchos::null;

  // build redundant vector of all non-close slave node ids on all procs
  // (there must not be any double entries in the node lists, thus
  // transform to sets and then back to vectors)
  std::vector<int> globalfns;
  std::set<int> setglobalfns;
  std::vector<int> sncnids;
  LINALG::Gather<int>(localfns,globalfns,numproc,&allproc[0],Comm());
  for (int i=0;i<(int)globalfns.size();++i) setglobalfns.insert(globalfns[i]);
  for (iter=setglobalfns.begin();iter!=setglobalfns.end();++iter) sncnids.push_back(*iter);

  //**********************************************************************
  // call ZOLTAN for parallel redistribution
  DRT::UTILS::PartUsingParMetis(idiscret_,sncroweles,sncrownodes,snccolnodes,comm,false,sncproc);
  //**********************************************************************

  //**********************************************************************
  // (6) MASTER redistribution
  //**********************************************************************
  Teuchos::RCP<Epetra_Map> mrownodes = Teuchos::null;
  Teuchos::RCP<Epetra_Map> mcolnodes = Teuchos::null;

  RedistributeMasterSide( mrownodes, mcolnodes, mroweles, comm, mproc );

  //**********************************************************************
  // (7) Merge global interface node row and column map
  //**********************************************************************
  // merge slave node row map from close and non-close parts
  Teuchos::RCP<Epetra_Map> srownodes = Teuchos::null;

  //----------------------------------CASE 1: ONE OR BOTH SLAVE SETS EMPTY
  if (scrownodes==Teuchos::null || sncrownodes==Teuchos::null)
  {
     dserror("ERROR: Redistribute: You should not be here");
  }
  //-------------------------------------CASE 2: BOTH SLAVE SETS NON-EMPTY
  else
  {
    // find intersection set of close and non-close nodes
    std::set<int> intersec;
    for (iter=setglobalcns.begin();iter!=setglobalcns.end();++iter)
    {
      std::set<int>::const_iterator found = setglobalfns.find(*iter);
      if (found!=setglobalfns.end()) intersec.insert(*found);
    }

    // build slave node row map
    std::vector<int> mygids(scrownodes->NumMyElements() + sncrownodes->NumMyElements());
    int count = scrownodes->NumMyElements();

    // first get GIDs of input scrownodes
    for (int i=0;i<count;++i) mygids[i] = scrownodes->GID(i);

    // then add GIDs of input sncrownodes (only new ones)
    for (int i=0;i<sncrownodes->NumMyElements();++i)
    {
      // check for intersection gid
      // don't do anything for intersection gids (scrownodes dominates!!!)
      std::set<int>::const_iterator found = intersec.find(sncrownodes->GID(i));
      if (found!=intersec.end()) continue;

      // check for overlap
      if (scrownodes->MyGID(sncrownodes->GID(i)))
        dserror("LINALG::MergeMap: Result map is overlapping");

      // add new GIDs to mygids
      mygids[count]=sncrownodes->GID(i);
      ++count;
    }
    mygids.resize(count);
    sort(mygids.begin(),mygids.end());
    srownodes = Teuchos::rcp(new Epetra_Map(-1,(int)mygids.size(),&mygids[0],0,scrownodes->Comm()));
  }

  // merge interface node row map from slave and master parts
  Teuchos::RCP<Epetra_Map> rownodes = LINALG::MergeMap(srownodes,mrownodes,false);

  // IMPORTANT NOTE:
  // While merging from the two different slave parts of the discretization
  // (close slave, non-close slave) is feasible for the node row map,
  // this is not possible for the node column map. Some necessary
  // information on ghosting at the transition between close and non-close
  // slave region would always be missed! Thus, we reconstruct a
  // suitable slave node column map "by hand" here. This is quite simply
  // done by exporting the initial node graph to the new distribution
  // and by then asking for its column map.

   // create the output graph (with new slave node row map) and export to it
   Teuchos::RCP<Epetra_CrsGraph> outgraph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*srownodes,108,false));
   Epetra_Export exporter(graph->RowMap(),*srownodes);
   int err = outgraph->Export(*graph,exporter,Add);
   if (err<0) dserror("Graph export returned err=%d",err);

  // trash old graph
  graph=Teuchos::null;

  // call fill complete and optimize storage
  outgraph->FillComplete();
  outgraph->OptimizeStorage();

  // get column map from the graph -> build slave node column map
  // (do stupid conversion from Epetra_BlockMap to Epetra_Map)
  const Epetra_BlockMap& bcol = outgraph->ColMap();
  Teuchos::RCP<Epetra_Map> scolnodes = Teuchos::rcp(new Epetra_Map(bcol.NumGlobalElements(),bcol.NumMyElements(),bcol.MyGlobalElements(),0,Comm()));

  // trash new graph
  outgraph=Teuchos::null;

  // merge interface node column map from slave and master parts
  Teuchos::RCP<Epetra_Map> colnodes = LINALG::MergeMap(scolnodes,mcolnodes,false);

  //**********************************************************************
  // (8) Get partitioning information into discretization
  //**********************************************************************
  // build reasonable element maps from the already valid and final node maps
  // (note that nothing is actually redistributed in here)
  Teuchos::RCP<Epetra_Map> roweles  = Teuchos::null;
  Teuchos::RCP<Epetra_Map> coleles  = Teuchos::null;
  Discret().BuildElementRowColumn(*rownodes,*colnodes,roweles,coleles);

  // export nodes and elements to the row map
  Discret().ExportRowNodes(*rownodes);
  Discret().ExportRowElements(*roweles);

  // export nodes and elements to the column map (create ghosting)
  Discret().ExportColumnNodes(*colnodes);
  Discret().ExportColumnElements(*coleles);

  return true;
}

/*----------------------------------------------------------------------*
 | collect distribution data (public)                         popp 10/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::CollectDistributionData(int& loadele, int& crowele)
{
  // loop over proc's column slave elements of the interface
  for (int i=0; i<selecolmap_->NumMyElements();++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1) dserror("ERROR: Cannot find slave element with gid %",gid1);
    CoElement* selement = dynamic_cast<CoElement*>(ele1);

    // bool indicating coupling partners
    bool add = (selement->MoData().NumSearchElements()>0);

    // check if this element has any coupling partners and add
    // element ID to input variable loadele if so
    if (add) loadele += 1;

    // check if - in addition - the active proc owns this element
    // and add element ID to input variable rowele if so
    if (add && selement->Owner()==Comm().MyPID()) crowele += 1;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  create search tree (public)                               popp 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::CreateSearchTree()
{
  // ***WARNING:*** This is commented out here, as idiscret_->SetState()
  // needs all the procs around, not only the interface local ones!
  // if (!lComm()) return;

  // warning
#ifdef MORTARGMSHCTN
  if (Dim()==3 && Comm().MyPID()==0)
  {
    std::cout << "\n******************************************************************\n";
    std::cout << "GMSH output of all contact tree nodes in 3D needs a lot of memory!\n";
    std::cout << "******************************************************************\n";
  }
#endif // #ifdef MORTARGMSHCTN

  // binary tree search
  if (SearchAlg()==INPAR::MORTAR::search_binarytree)
  {
    //*****SELF CONTACT*****
    if (SelfContact())
    {
      // set state in interface to intialize all kinds of quantities
      Teuchos::RCP<Epetra_Vector> zero =Teuchos::rcp(new Epetra_Vector(*idiscret_->DofRowMap()));
      SetState(MORTAR::state_new_displacement,*zero);

      // create fully overlapping map of all contact elements
      Teuchos::RCP<Epetra_Map> elefullmap = LINALG::AllreduceEMap(*idiscret_->ElementRowMap());

      // create binary tree object for self contact search
      // (NOTE THAT SELF CONTACT SEARCH IS NOT YET FULLY PARALLELIZED!)
      binarytreeself_ = Teuchos::rcp(new CONTACT::SelfBinaryTree(Discret(),lComm(),elefullmap,Dim(),SearchParam()));

    }
    //*****TWO BODY CONTACT*****
    else
    {
      // get out of here if not participating in interface
      if (!lComm()) return;

      // create fully overlapping map of all master elements
      // for non-redundant storage (RRloop) we handle the master elements
      // like the slave elements --> melecolmap_
      INPAR::MORTAR::ParallelStrategy strat =
          DRT::INPUT::IntegralValue<INPAR::MORTAR::ParallelStrategy>(IParams(),"PARALLEL_STRATEGY");

      Teuchos::RCP<Epetra_Map> melefullmap = Teuchos::null;
      if (strat==INPAR::MORTAR::roundrobinghost || strat==INPAR::MORTAR::binningstrategy)
      {
        melefullmap = melecolmap_;
      }
      else if (strat==INPAR::MORTAR::roundrobinevaluate)
      {
        melefullmap = melerowmap_;
      }
      else if (strat==INPAR::MORTAR::ghosting_redundant)
      {
        melefullmap = LINALG::AllreduceEMap(*melerowmap_);
      }
      else
        dserror("Chosen parallel strategy not supported!");

      // create binary tree object for contact search and setup tree
      binarytree_ = Teuchos::rcp(new MORTAR::BinaryTree(Discret(),selecolmap_,melefullmap,Dim(),SearchParam(),SearchUseAuxPos()));

      // initialize active contact nodes via binarytree
      // binarytree_->SearchContactInit(binarytree_->Sroot(), binarytree_->Mroot());
    }
  }

  // no binary tree search
  else
  {
    if (SelfContact())
      dserror("ERROR: Binarytree search needed for self contact");
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Initialize Data Container for nodes and elements         farah 02/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::InitializeDataContainer()
{
  // call base class functionality
  MORTAR::MortarInterface::InitializeDataContainer();

  // ==================
  // non-smooth contact:
  // we need this master node data container to create an averaged
  // nodal normal field on the master side for the smoothed cpp
  // normal field!
  if(DRT::INPUT::IntegralValue<int>(IParams(),"CPP_NORMALS") || nonSmoothContact_)
  {
    const Teuchos::RCP<Epetra_Map> masternodes = LINALG::AllreduceEMap(*(MasterRowNodes()));

    for (int i = 0; i < masternodes->NumMyElements(); ++i)
    {
      int gid = masternodes->GID(i);
      DRT::Node* node = Discret().gNode(gid);
      if (!node)
        dserror("ERROR: Cannot find node with gid %i", gid);
      CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(node);
      mnode->InitializeDataContainer();
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  initialize / reset interface for contact                  popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::Initialize()
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over all nodes to reset stuff (fully overlapping column map)
  // (use fully overlapping column map)

  for (int i=0;i<idiscret_->NumMyColNodes();++i)
  {
    CONTACT::CoNode* node = dynamic_cast<CONTACT::CoNode*>(idiscret_->lColNode(i));

    // reset feasible projection and segmentation status
    node->HasProj()    = false;
    node->HasSegment() = false;
  }

  // init normal data in master node data container for cpp calculation
  if(DRT::INPUT::IntegralValue<int>(IParams(),"CPP_NORMALS"))
  {
    for (int i=0;i<MasterColNodes()->NumMyElements();++i)
    {
      int gid = MasterColNodes()->GID(i);
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      CoNode* cnode = dynamic_cast<CoNode*>(node);

      // reset derivative maps of normal vector
      for (int j=0;j<(int)((cnode->CoData().GetDerivN()).size());++j)
        (cnode->CoData().GetDerivN())[j].clear();
      (cnode->CoData().GetDerivN()).resize(0,0);

      // reset derivative maps of tangent vectors
      for (int j=0;j<(int)((cnode->CoData().GetDerivTxi()).size());++j)
        (cnode->CoData().GetDerivTxi())[j].clear();
      (cnode->CoData().GetDerivTxi()).resize(0,0);
      for (int j=0;j<(int)((cnode->CoData().GetDerivTeta()).size());++j)
        (cnode->CoData().GetDerivTeta())[j].clear();
      (cnode->CoData().GetDerivTeta()).resize(0,0);

      for (int j=0;j<(int)((cnode->CoData().GetDerivTangent()).size());++j)
        (cnode->CoData().GetDerivTangent())[j].clear();
      (cnode->CoData().GetDerivTangent()).resize(0,0);
    }
  }

  // loop over all slave nodes to reset stuff (standard column map)
  // (include slave side boundary nodes / crosspoints)
  for (int i=0;i<SlaveColNodesBound()->NumMyElements();++i)
  {
    int gid = SlaveColNodesBound()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    // reset nodal Mortar maps
    // for sts
    cnode->MoData().GetD().clear();
    cnode->MoData().GetM().clear();
    cnode->MoData().GetMmod().clear();
    // for nts
    cnode->MoData().GetDnts().clear();
    cnode->MoData().GetMnts().clear();
    // for lts
    cnode->MoData().GetDlts().clear();
    cnode->MoData().GetMlts().clear();
    // for ltl
    cnode->MoData().GetDltl().clear();
    cnode->MoData().GetMltl().clear();

    // reset derivative maps of normal vector
    for (int j=0;j<(int)((cnode->CoData().GetDerivN()).size());++j)
      (cnode->CoData().GetDerivN())[j].clear();
    (cnode->CoData().GetDerivN()).resize(0,0);

    // reset derivative maps of tangent vectors
    for (int j=0;j<(int)((cnode->CoData().GetDerivTxi()).size());++j)
      (cnode->CoData().GetDerivTxi())[j].clear();
    (cnode->CoData().GetDerivTxi()).resize(0,0);
    for (int j=0;j<(int)((cnode->CoData().GetDerivTeta()).size());++j)
      (cnode->CoData().GetDerivTeta())[j].clear();
    (cnode->CoData().GetDerivTeta()).resize(0,0);

    // reset derivative map of Mortar matrices
    (cnode->CoData().GetDerivD()).clear();
    (cnode->CoData().GetDerivDlts()).clear();
    (cnode->CoData().GetDerivDltl()).clear();
    (cnode->CoData().GetDerivM()).clear();
    (cnode->CoData().GetDerivMnts()).clear();
    (cnode->CoData().GetDerivMlts()).clear();
    (cnode->CoData().GetDerivMltl()).clear();

    // reset nodal weighted gap and derivative
    cnode->CoData().Getg() = 1.0e12;
    cnode->CoData().Getgnts() = 1.0e12;
    cnode->CoData().Getglts() = 1.0e12;
    cnode->CoData().Getgltl()[0] = 1.0e12;
    cnode->CoData().Getgltl()[1] = 1.0e12;
    cnode->CoData().Getgltl()[2] = 1.0e12;
    cnode->MoData().GetDscale() = 0.0;
    (cnode->CoData().GetDerivG()).clear();
    (cnode->CoData().GetDerivGnts()).clear();
    (cnode->CoData().GetDerivGlts()).clear();
    for(int j=0;j<(int)cnode->CoData().GetDerivGltl().size();++j)
      cnode->CoData().GetDerivGltl()[j].clear();
    for(int j=0;j<(int)cnode->CoData().GetDerivJumpltl().size();++j)
      cnode->CoData().GetDerivJumpltl()[j].clear();
//    (cnode->CoData().GetDerivGltl()).resize(0);

    // reset nodal jump
    cnode->CoData().Getjumpltl()[0] = 1.0e12;
    cnode->CoData().Getjumpltl()[1] = 1.0e12;
    cnode->CoData().Getjumpltl()[2] = 1.0e12;

    // hybrid formulation
    cnode->CoData().GetAlphaN() = -1.0;
    cnode->CoData().GetAlpha().clear();

    // reset derivative map of lagrange multipliers
    for (int j=0; j<(int)((cnode->CoData().GetDerivZ()).size()); ++j)
      (cnode->CoData().GetDerivZ())[j].clear();
    (cnode->CoData().GetDerivZ()).resize(0);

    if (friction_)
    {
      FriNode* frinode = dynamic_cast<FriNode*>(cnode);

      // reset SNodes and Mnodes
      frinode->FriData().GetSNodes().clear();
      frinode->FriData().GetMNodes().clear();

      // for gp slip
      if (DRT::INPUT::IntegralValue<int>(IParams(),"GP_SLIP_INCR")==true)
      {
        // reset jump deriv.
        for (int j=0 ;j<(int)((frinode->FriData().GetDerivVarJump()).size()); ++j )
          (frinode->FriData().GetDerivVarJump())[j].clear();

        (frinode->FriData().GetDerivVarJump()).resize(2);

        // reset jumps
        frinode->FriData().jump_var()[0]=0.0;
        frinode->FriData().jump_var()[1]=0.0;
      }
    }

    // just do poro contact relevant stuff!
    if (poro_)
    {
      cnode->CoPoroData().GetnCoup() = 0.0;
      cnode->CoPoroData().GetDerivnCoup().clear();
      cnode->CoPoroData().GetVelDerivnCoup().clear();
      cnode->CoPoroData().GetPresDerivnCoup().clear();
    }

    // just do ehl relevant stuff!
    if (ehl_)
      cnode->CoEhlData().Clear();
  }

  //**********************************************************************
  // In general, it is sufficient to reset search candidates only for
  // all elements in the standard slave column map. However, self contact
  // is an exception here and we need to reset the search candidates of
  // all slave elements in the fully overlapping column map there. This
  // is due to the fact that self contact search is NOT parallelized.
  //**********************************************************************
  if (SelfContact())
  {
    // loop over all elements to reset candidates / search lists
    // (use fully overlapping column map of S+M elements)
    for (int i=0;i<idiscret_->NumMyColElements();++i)
    {
      DRT::Element* ele = idiscret_->lColElement(i);
      MORTAR::MortarElement* mele = dynamic_cast<MORTAR::MortarElement*>(ele);

      mele->MoData().SearchElements().resize(0);

      // dual shape function coefficient matrix
      mele->MoData().ResetDualShape();
      mele->MoData().ResetDerivDualShape();
    }
  }
  else
  {
    // loop over all elements to reset candidates / search lists
    // (use standard slave column map)
    for (int i=0;i<SlaveColElements()->NumMyElements();++i)
    {
      int gid = SlaveColElements()->GID(i);
      DRT::Element* ele = Discret().gElement(gid);
      if (!ele) dserror("ERROR: Cannot find ele with gid %i",gid);
      MORTAR::MortarElement* mele = dynamic_cast<MORTAR::MortarElement*>(ele);

      mele->MoData().SearchElements().resize(0);

      // dual shape function coefficient matrix
      mele->MoData().ResetDualShape();
      mele->MoData().ResetDerivDualShape();
    }
  }

  // clear all Nitsche data
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(imortar_, "ALGORITHM")
      ==INPAR::MORTAR::algorithm_gpts)
    for (int e=0;e<Discret().ElementColMap()->NumMyElements();++e)
      dynamic_cast<MORTAR::MortarElement*>(Discret().gElement(
                Discret().ElementColMap()->GID(e)))->GetNitscheContainer().Clear();

  // reset s/m pairs and intcell counters
  smpairs_    = 0;
  smintpairs_ = 0;
  intcells_   = 0;

  return;
}

/*----------------------------------------------------------------------*
 |  compute element areas (public)                            popp 11/07|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::SetElementAreas()
{
  //**********************************************************************
  // In general, it is sufficient to compute element areas only for
  // all elements in the standard slave column map. However, self contact
  // is an exception here and we need the element areas of all elements
  // (slave and master) in the fully overlapping column map there. At the
  // same time we initialize the element data containers for self contact.
  // This is due to the fact that self contact search is NOT parallelized.
  //**********************************************************************
  if (SelfContact() or DRT::INPUT::IntegralValue<int>(IParams(),"CPP_NORMALS") or
      nonSmoothContact_)
  {
    // loop over all elements to set current element length / area
    // (use fully overlapping column map)
    for (int i=0;i<idiscret_->NumMyColElements();++i)
    {
      MORTAR::MortarElement* element = dynamic_cast<MORTAR::MortarElement*>(idiscret_->lColElement(i));
      element->InitializeDataContainer();
      element->MoData().Area()=element->ComputeArea();
    }
  }
  else
  {
    //refer call back to base class version
    MORTAR::MortarInterface::SetElementAreas();
  }

  return;
}


/*----------------------------------------------------------------------*
 |  pre evaluate to calc normals                            farah 02/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::PreEvaluate(
    const int& step,
    const int& iter)
{
  TEUCHOS_FUNC_TIME_MONITOR( "CONTACT::CoInterface::PreEvaluate" );

  //**********************************************************************
  // search algorithm
  //**********************************************************************
  if (SearchAlg() == INPAR::MORTAR::search_bfele)
    EvaluateSearchBruteForce(SearchParam());
  else if (SearchAlg() == INPAR::MORTAR::search_binarytree)
    EvaluateSearchBinarytree();
  else
    dserror("ERROR: Invalid search algorithm");

  // TODO: maybe we can remove this debug functionality
#ifdef MORTARGMSHCELLS
  // reset integration cell GMSH files
  int proc = Comm().MyPID();
  std::ostringstream filename;
  filename << "o/gmsh_output/cells_" << proc << ".pos";
  FILE* fp = fopen(filename.str().c_str(), "w");
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \"Integration Cells Proc " << proc << "\" {" << std::endl;
  fprintf(fp,gmshfilecontent.str().c_str());
  fclose(fp);
#endif // #ifdef MORTARGMSHCELLS

  // set global vector of cn values
  SetCnCtValues(iter);

  // cpp normals or averaged normal field?
  if(DRT::INPUT::IntegralValue<int>(IParams(),"CPP_NORMALS"))
  {
    // evaluate cpp nodal normals on slave side
    EvaluateCPPNormals();
  }
  else
  {
    // evaluate averaged nodal normals on slave side
    EvaluateNodalNormals();

    // export nodal normals to slave node column map
    // this call is very expensive and the computation
    // time scales directly with the proc number !
    ExportNodalNormals();
  }

  // compute scaling between coupling types
//  if(nonSmoothContact_)
//    ComputeScaling();

  // bye bye
  return;
}


/*----------------------------------------------------------------------*
 |  store nts values into sts data container for assembly   farah 07/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::StoreNTSvalues()
{
  // create iterators for data types
  typedef GEN::pairedvector<int,double>::const_iterator CI;
  typedef std::map<int,double>::const_iterator          CImap;

  // loop over all possibly non smooth nodes
  for(int i=0; i<SlaveRowNodes()->NumMyElements();++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    // check if integration is done
    if (cnode->MoData().GetDnts().size()<1)
      continue;

    // if nonsmooth contact is activated and the node is no corner node continue
    // if non-smooth contact is not activated go on
    if(!cnode->IsOnCorner() and nonSmoothContact_)
      continue;

    //-------------------------------------------------------------------------------------
    // store D matrix entries
    // resize pairedvector to nts size
    if ((int)cnode->MoData().GetD().size()==0)
      cnode->MoData().GetD().resize(cnode->MoData().GetDnts().size());

    for (CI p=cnode->MoData().GetDnts().begin();p!=cnode->MoData().GetDnts().end();++p)
      cnode->MoData().GetD()[p->first] += (p->second);

    //-------------------------------------------------------------------------------------
    // store M matrix entries
    for (CImap p=cnode->MoData().GetMnts().begin();p!=cnode->MoData().GetMnts().end();++p)
      cnode->MoData().GetM()[p->first] += (p->second);

    //-------------------------------------------------------------------------------------
    // store weighted gap
    cnode->CoData().Getg() = cnode->CoData().Getgnts();

    //-------------------------------------------------------------------------------------
    // store weighted gap linearization
    for (CImap p=cnode->CoData().GetDerivGnts().begin();p!=cnode->CoData().GetDerivGnts().end();++p)
      cnode->CoData().GetDerivG()[p->first] += (p->second);

    //-------------------------------------------------------------------------------------
    // store D deriv
    // --> No D linearization!

    //-------------------------------------------------------------------------------------
    // store M deriv
    {
      // Mortar M derivatives
      std::map<int,std::map<int,double> >& mntsderiv = cnode->CoData().GetDerivMnts();

      // get sizes and iterator start
      int mastersize = (int)mntsderiv.size();
      std::map<int,std::map<int,double> >::iterator mntscurr = mntsderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l=0;l<mastersize;++l)
      {
        int mgid = mntscurr->first;
        ++mntscurr;

        // Mortar matrix M derivatives
        std::map<int,double>&thismderivnts    = cnode->CoData().GetDerivMnts()[mgid];
        std::map<int,double>&thismderivmortar = cnode->CoData().GetDerivM()[mgid];

        int mapsize = (int)(thismderivnts.size());

        std::map<int,double>::iterator mcolcurr = thismderivnts.begin();

        // loop over all directional derivative entries
        for (int c=0;c<mapsize;++c)
        {
          thismderivmortar[mcolcurr->first] += (mcolcurr->second);
          ++mcolcurr;
        }

        // check for completeness of DerivM-Derivatives-iteration
        if (mcolcurr!=thismderivnts.end())
          dserror("ERROR: StoreNTS: Not all derivative entries of DerivM considered!");
      }
    }
  }// end node loop


  return;
}


/*----------------------------------------------------------------------*
 |  store lts values into sts data container for assembly   farah 07/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::StoreLTSvalues()
{
  // create iterators for data types
  typedef GEN::pairedvector<int,double>::const_iterator CI;
  typedef std::map<int,double>::const_iterator          CImap;

  // loop over all possibly non smooth nodes
  for(int i=0; i<SlaveRowNodes()->NumMyElements();++i)
  {
    double msum = 0.0;
    double ssum = 0.0;

    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    // check if this is an edge or a corner and nonsmooth contact is activated
    if(!cnode->IsOnEdge() and nonSmoothContact_)
      continue;
    if(cnode->IsOnCorner() and nonSmoothContact_)
      continue;

    // check if integration is done
    if (cnode->MoData().GetDlts().size()<1)
      continue;

    //-------------------------------------------------------------------------------------
    // store D matrix entries
    // resize pairedvector to nts size
    if ((int)cnode->MoData().GetD().size()==0)
      cnode->MoData().GetD().resize(cnode->MoData().GetDlts().size() + cnode->MoData().GetDltl().size());

    for (CI p=cnode->MoData().GetDlts().begin();p!=cnode->MoData().GetDlts().end();++p)
    {
      cnode->MoData().GetD()[p->first] += (p->second);
      ssum += (p->second);
    }

    //-------------------------------------------------------------------------------------
    // store M matrix entries
    for (CImap p=cnode->MoData().GetMlts().begin();p!=cnode->MoData().GetMlts().end();++p)
    {
      cnode->MoData().GetM()[p->first] += (p->second);
      msum += (p->second);
    }

    //-------------------------------------------------------------------------------------
    // store weighted gap
    cnode->CoData().Getg() = cnode->CoData().Getglts();

    //-------------------------------------------------------------------------------------
    // store weighted gap linearization
    for (CImap p=cnode->CoData().GetDerivGlts().begin();p!=cnode->CoData().GetDerivGlts().end();++p)
      cnode->CoData().GetDerivG()[p->first] += (p->second);

    //-------------------------------------------------------------------------------------
    // store D deriv
    {
      // Mortar M derivatives
      std::map<int,std::map<int,double> >& mntsderiv = cnode->CoData().GetDerivDlts();

      // get sizes and iterator start
      int mastersize = (int)mntsderiv.size();
      std::map<int,std::map<int,double> >::iterator mntscurr = mntsderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l=0;l<mastersize;++l)
      {
        int mgid = mntscurr->first;
        ++mntscurr;

        // Mortar matrix M derivatives
        std::map<int,double>&thismderivnts    = cnode->CoData().GetDerivDlts()[mgid];
        std::map<int,double>&thismderivmortar = cnode->CoData().GetDerivD()[mgid];

        int mapsize = (int)(thismderivnts.size());

        std::map<int,double>::iterator mcolcurr = thismderivnts.begin();

        // loop over all directional derivative entries
        for (int c=0;c<mapsize;++c)
        {
          thismderivmortar[mcolcurr->first] += (mcolcurr->second);
          ++mcolcurr;
        }

        // check for completeness of DerivM-Derivatives-iteration
        if (mcolcurr!=thismderivnts.end())
          dserror("ERROR: StoreNTS: Not all derivative entries of DerivM considered!");
      }
    }

    //-------------------------------------------------------------------------------------
    // store M deriv
    {
      // Mortar M derivatives
      std::map<int,std::map<int,double> >& mntsderiv = cnode->CoData().GetDerivMlts();

      // get sizes and iterator start
      int mastersize = (int)mntsderiv.size();
      std::map<int,std::map<int,double> >::iterator mntscurr = mntsderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l=0;l<mastersize;++l)
      {
        int mgid = mntscurr->first;
        ++mntscurr;

        // Mortar matrix M derivatives
        std::map<int,double>&thismderivnts    = cnode->CoData().GetDerivMlts()[mgid];
        std::map<int,double>&thismderivmortar = cnode->CoData().GetDerivM()[mgid];

        int mapsize = (int)(thismderivnts.size());

        std::map<int,double>::iterator mcolcurr = thismderivnts.begin();

        // loop over all directional derivative entries
        for (int c=0;c<mapsize;++c)
        {
          thismderivmortar[mcolcurr->first] += (mcolcurr->second);
          ++mcolcurr;
        }

        // check for completeness of DerivM-Derivatives-iteration
        if (mcolcurr!=thismderivnts.end())
          dserror("ERROR: StoreNTS: Not all derivative entries of DerivM considered!");
      }
    }
//    std::cout << "ssum = " << ssum << "  msum = " << msum << "  balance= " << ssum-msum << std::endl;
    if( abs(ssum-msum)>1e-12 )
      dserror("ERROR: no slave master balance!");

  }// end node loop


  return;
}



/*----------------------------------------------------------------------*
 |  store lts values into sts data container for assembly   farah 07/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::StoreLTLvalues()
{
  dserror("ERROR: StoreLTLvalues() is outdated!");
  return;
//  // create iterators for data types
//  typedef GEN::pairedvector<int,double>::const_iterator CI;
//  typedef std::map<int,double>::const_iterator          CImap;
//
//  // loop over all possibly non smooth nodes
//  for(int i=0; i<SlaveRowNodes()->NumMyElements();++i)
//  {
//    int gid = SlaveRowNodes()->GID(i);
//    DRT::Node* node = idiscret_->gNode(gid);
//    if (!node)
//      dserror("ERROR: Cannot find node with gid %",gid);
//    CoNode* cnode = dynamic_cast<CoNode*>(node);
//
//    if(!cnode->IsOnEdge())
//      continue;
//    if(cnode->IsOnCorner())
//      continue;
//
//    // check if integration is done
//    if (cnode->MoData().GetDltl().size()<1)
//      continue;
//
//    //-------------------------------------------------------------------------------------
//    // store D matrix entries
//    // resize pairedvector to nts size
//    if ((int)cnode->MoData().GetD().size()==0)
//      cnode->MoData().GetD().resize(cnode->MoData().GetDltl().size());
//
//    for (CI p=cnode->MoData().GetDltl().begin();p!=cnode->MoData().GetDltl().end();++p)
//      cnode->MoData().GetD()[p->first] += (p->second);
//
//    //-------------------------------------------------------------------------------------
//    // store M matrix entries
//    for (CImap p=cnode->MoData().GetMltl().begin();p!=cnode->MoData().GetMltl().end();++p)
//      cnode->MoData().GetM()[p->first] += (p->second);
//
//    //-------------------------------------------------------------------------------------
//    // store weighted gap
//    cnode->CoData().Getg() = cnode->CoData().Getgltl();
//
//    //-------------------------------------------------------------------------------------
//    // store weighted gap linearization
//    for (CImap p=cnode->CoData().GetDerivGltl().begin();p!=cnode->CoData().GetDerivGltl().end();++p)
//      cnode->CoData().GetDerivG()[p->first] += (p->second);
//
//    //-------------------------------------------------------------------------------------
//    // store D deriv
//    {
//      // Mortar M derivatives
//      std::map<int,std::map<int,double> >& mntsderiv = cnode->CoData().GetDerivDltl();
//
//      // get sizes and iterator start
//      int mastersize = (int)mntsderiv.size();
//      std::map<int,std::map<int,double> >::iterator mntscurr = mntsderiv.begin();
//
//      /********************************************** LinMMatrix **********/
//      // loop over all master nodes in the DerivM-map of the current LM slave node
//      for (int l=0;l<mastersize;++l)
//      {
//        int mgid = mntscurr->first;
//        ++mntscurr;
//
//        // Mortar matrix M derivatives
//        std::map<int,double>&thismderivnts    = cnode->CoData().GetDerivDltl()[mgid];
//        std::map<int,double>&thismderivmortar = cnode->CoData().GetDerivD()[mgid];
//
//        int mapsize = (int)(thismderivnts.size());
//
//        std::map<int,double>::iterator mcolcurr = thismderivnts.begin();
//
//        // loop over all directional derivative entries
//        for (int c=0;c<mapsize;++c)
//        {
//          thismderivmortar[mcolcurr->first] += (mcolcurr->second);
//          ++mcolcurr;
//        }
//
//        // check for completeness of DerivM-Derivatives-iteration
//        if (mcolcurr!=thismderivnts.end())
//          dserror("ERROR: StoreNTS: Not all derivative entries of DerivM considered!");
//      }
//    }
//
//    //-------------------------------------------------------------------------------------
//    // store M deriv
//    {
//      // Mortar M derivatives
//      std::map<int,std::map<int,double> >& mntsderiv = cnode->CoData().GetDerivMltl();
//
//      // get sizes and iterator start
//      int mastersize = (int)mntsderiv.size();
//      std::map<int,std::map<int,double> >::iterator mntscurr = mntsderiv.begin();
//
//      /********************************************** LinMMatrix **********/
//      // loop over all master nodes in the DerivM-map of the current LM slave node
//      for (int l=0;l<mastersize;++l)
//      {
//        int mgid = mntscurr->first;
//        ++mntscurr;
//
//        // Mortar matrix M derivatives
//        std::map<int,double>&thismderivnts    = cnode->CoData().GetDerivMltl()[mgid];
//        std::map<int,double>&thismderivmortar = cnode->CoData().GetDerivM()[mgid];
//
//        int mapsize = (int)(thismderivnts.size());
//
//        std::map<int,double>::iterator mcolcurr = thismderivnts.begin();
//
//        // loop over all directional derivative entries
//        for (int c=0;c<mapsize;++c)
//        {
//          thismderivmortar[mcolcurr->first] += (mcolcurr->second);
//          ++mcolcurr;
//        }
//
//        // check for completeness of DerivM-Derivatives-iteration
//        if (mcolcurr!=thismderivnts.end())
//          dserror("ERROR: StoreNTS: Not all derivative entries of DerivM considered!");
//      }
//    }
//  }// end node loop
//  return;
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces                         farah 10/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AddLTLforcesFric(Teuchos::RCP<Epetra_FEVector> feff)
{
  double slaveforce = 0.0;
  double masterforce = 0.0;

  const double penalty    = IParams().get<double>("PENALTYPARAM");
  const double penaltytan = IParams().get<double>("PENALTYPARAMTAN");
  const double frcoeff    = IParams().get<double>("FRCOEFF");

  double oldtraction[3] = {0.0 , 0.0 , 0.0};

  typedef GEN::pairedvector<int,double>::const_iterator CI;
  typedef std::map<int,double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j=0;j<snoderowmap_->NumMyElements();++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    double x = cnode->FriData().tractionoldLTL()[0] * cnode->FriData().tractionoldLTL()[0];
    double y = cnode->FriData().tractionoldLTL()[1] * cnode->FriData().tractionoldLTL()[1];
    double z = cnode->FriData().tractionoldLTL()[2] * cnode->FriData().tractionoldLTL()[2];
    double tracvalue = sqrt(x+y+z);

    if(tracvalue>1e-8)
    {
      oldtraction[0] = cnode->FriData().tractionoldLTL()[0];
      oldtraction[1] = cnode->FriData().tractionoldLTL()[1];
      oldtraction[2] = cnode->FriData().tractionoldLTL()[2];
      break;
    }
  }

  // maybe the old traction is here zero (first contact step)
  // loop over all slave nodes
  for (int j=0;j<snoderowmap_->NumMyElements();++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    // check if this is active node
    if(cnode->CoData().Getgltl()[0]<1e8 and
       cnode->CoData().Getgltl()[1]<1e8 and
       cnode->CoData().Getgltl()[2]<1e8)
    {
      // normal force
      double fn[3] = {0.0 , 0.0 , 0.0};
      for(int dim = 0; dim<Dim();++dim)
        fn[dim] = -penalty * cnode->CoData().Getgltl()[dim];

      // f trial tangential
      double ftrial[3] = {0.0 , 0.0 , 0.0};
      for(int dim = 0; dim<Dim();++dim)
        ftrial[dim] = oldtraction[dim] - penaltytan * cnode->CoData().Getjumpltl()[dim];

      // trial norm
      double trialnorm = sqrt(ftrial[0]*ftrial[0] + ftrial[1]*ftrial[1] + ftrial[2]*ftrial[2]);

      // maxtrac
      double maxtrac = sqrt(fn[0]*fn[0] + fn[1]*fn[1] + fn[2]*fn[2]);

      // real traction
      double ftan[3] = {0.0 , 0.0 , 0.0};

      if(trialnorm - frcoeff*maxtrac <= 0.0)
      {
        for(int dim = 0; dim<Dim();++dim)
          ftan[dim] = ftrial[dim];
      }
      else
      {
        double coeff = frcoeff * maxtrac/trialnorm;
        for(int dim = 0; dim<Dim();++dim)
          ftan[dim] = coeff* ftrial[dim];
      }

      // store
      cnode->FriData().traction()[0] = ftan[0];
      cnode->FriData().traction()[1] = ftan[1];
      cnode->FriData().traction()[2] = ftan[2];

      // ASSEMBLE
      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetDltl()).size()>0)
      {
        GEN::pairedvector<int,double>  map = cnode->MoData().GetDltl();

        for (CI p=map.begin();p!=map.end();++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("ERROR: Cannot find node with gid");
          CoNode* csnode = dynamic_cast<CoNode*>(snode);

          for(int dim = 0; dim<Dim();++dim)
          {
            double value[1] = {(p->second) * ftan[dim]};
            const int ltlid[1] = {csnode->Dofs()[dim]};
            int err = feff->SumIntoGlobalValues(1,ltlid,value);
            if (err<0)
              dserror("stop");
            slaveforce += value[0];
          }
        }
      }
      else
      {
        dserror("ERROR: no d matrix entries available for ltlt contact");
      }

      /**************************************************** M-matrix ******/
      if ((cnode->MoData().GetMltl()).size()>0)
      {
        std::map<int,double>  map = cnode->MoData().GetMltl();

        for (CImap p=map.begin();p!=map.end();++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("ERROR: Cannot find node with gid");
          CoNode* csnode = dynamic_cast<CoNode*>(snode);

          for(int dim = 0; dim<Dim();++dim)
          {
            double value[1] = {-(p->second) * ftan[dim]};
            const int ltlid[1] = {csnode->Dofs()[dim]};
            int err = feff->SumIntoGlobalValues(1,ltlid,value);
            if (err<0)
              dserror("stop");

            masterforce += value[0];
          }
        }
      }
      else
      {
        dserror("ERROR: no m matrix entries available for ltlt contact");
      }

      break;
    }
  }

//  std::cout << "--------------------------------FRICTION:"<<std::endl;
//  std::cout << "LTL slaveforce  = " << slaveforce << std::endl;
//  std::cout << "LTL masterforce = " << masterforce << std::endl;
//  std::cout << "LTL sum         = " << slaveforce+masterforce << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces                         farah 10/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AddLTLstiffnessFric(Teuchos::RCP<LINALG::SparseMatrix> kteff)
{
  const double penalty    = IParams().get<double>("PENALTYPARAM");
  const double penaltytan = IParams().get<double>("PENALTYPARAMTAN");
  const double frcoeff    = IParams().get<double>("FRCOEFF");

  typedef GEN::pairedvector<int,double>::const_iterator CI;
  typedef std::map<int,double>::const_iterator CImap;

  double oldtraction[3] = {0.0 , 0.0 , 0.0};

  // loop over all slave nodes
  for (int j=0;j<snoderowmap_->NumMyElements();++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    double x = cnode->FriData().tractionoldLTL()[0] * cnode->FriData().tractionoldLTL()[0];
    double y = cnode->FriData().tractionoldLTL()[1] * cnode->FriData().tractionoldLTL()[1];
    double z = cnode->FriData().tractionoldLTL()[2] * cnode->FriData().tractionoldLTL()[2];
    double tracvalue = sqrt(x+y+z);

    if(tracvalue>1e-8)
    {
      oldtraction[0] = cnode->FriData().tractionoldLTL()[0];
      oldtraction[1] = cnode->FriData().tractionoldLTL()[1];
      oldtraction[2] = cnode->FriData().tractionoldLTL()[2];
      break;
    }
  }

  // maybe the old traction is here zero (first contact step)
  // loop over all slave nodes
  for (int j=0;j<snoderowmap_->NumMyElements();++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    // check if this is active node
    if(cnode->CoData().Getgltl()[0]<1e8 and
       cnode->CoData().Getgltl()[1]<1e8 and
       cnode->CoData().Getgltl()[2]<1e8)
    {
      // state
      bool stick = true;
      GEN::pairedvector<int,double> coefflin(100);

      // normal force
      double fn[3] = {0.0 , 0.0 , 0.0};
      for(int dim = 0; dim<Dim();++dim)
        fn[dim] = -penalty * cnode->CoData().Getgltl()[dim];

      // f trial tangential
      double ftrial[3] = {0.0 , 0.0 , 0.0};
      for(int dim = 0; dim<Dim();++dim)
        ftrial[dim] = oldtraction[dim] - penaltytan * cnode->CoData().Getjumpltl()[dim];

      double coeff = 0.0;

      // trial norm
      double trialnorm = sqrt(ftrial[0]*ftrial[0] + ftrial[1]*ftrial[1] + ftrial[2]*ftrial[2]);

      // maxtrac
      double maxtrac = sqrt(fn[0]*fn[0] + fn[1]*fn[1] + fn[2]*fn[2]);

      // real traction
      double ftan[3] = {0.0 , 0.0 , 0.0};

      if(trialnorm - frcoeff*maxtrac <= 0.0)
      {
        stick = true;
        for(int dim = 0; dim<Dim();++dim)
          ftan[dim] = ftrial[dim];
      }
      else
      {
        stick = false;
        coeff = frcoeff * maxtrac/trialnorm;
        for(int dim = 0; dim<Dim();++dim)
          ftan[dim] = coeff * ftrial[dim];

        GEN::pairedvector<int,double> fn_x(100);
        GEN::pairedvector<int,double> fn_y(100);
        GEN::pairedvector<int,double> fn_z(100);

        GEN::pairedvector<int,double> ft_x(100);
        GEN::pairedvector<int,double> ft_y(100);
        GEN::pairedvector<int,double> ft_z(100);

        for (CImap pp=cnode->CoData().GetDerivGltl()[0].begin();pp!=cnode->CoData().GetDerivGltl()[0].end();++pp)
          fn_x[pp->first] -= penalty * (pp->second);
        for (CImap pp=cnode->CoData().GetDerivGltl()[1].begin();pp!=cnode->CoData().GetDerivGltl()[1].end();++pp)
          fn_y[pp->first] -= penalty * (pp->second);
        for (CImap pp=cnode->CoData().GetDerivGltl()[2].begin();pp!=cnode->CoData().GetDerivGltl()[2].end();++pp)
          fn_z[pp->first] -= penalty * (pp->second);

        for (CImap pp=cnode->CoData().GetDerivJumpltl()[0].begin();pp!=cnode->CoData().GetDerivJumpltl()[0].end();++pp)
          ft_x[pp->first] -= penaltytan * (pp->second);
        for (CImap pp=cnode->CoData().GetDerivJumpltl()[1].begin();pp!=cnode->CoData().GetDerivJumpltl()[1].end();++pp)
          ft_y[pp->first] -= penaltytan * (pp->second);
        for (CImap pp=cnode->CoData().GetDerivJumpltl()[2].begin();pp!=cnode->CoData().GetDerivJumpltl()[2].end();++pp)
          ft_z[pp->first] -= penaltytan * (pp->second);

        GEN::pairedvector<int,double> maxtraclin(100);
        for (CI pp=fn_x.begin();pp!=fn_x.end();++pp)
          maxtraclin[pp->first] -= 0.5 * (1.0/maxtrac) * (pp->second) * 2.0 * fn[0] * pp->second;
        for (CI pp=fn_y.begin();pp!=fn_y.end();++pp)
          maxtraclin[pp->first] -= 0.5 * (1.0/maxtrac) * (pp->second) * 2.0 * fn[1] * pp->second;
        for (CI pp=fn_z.begin();pp!=fn_z.end();++pp)
          maxtraclin[pp->first] -= 0.5 * (1.0/maxtrac) * (pp->second) * 2.0 * fn[2] * pp->second;

        GEN::pairedvector<int,double> trialnormlin(100);
        for (CI pp=ft_x.begin();pp!=ft_x.end();++pp)
          trialnormlin[pp->first] -= 0.5 * (1.0/trialnorm) * (pp->second) * 2.0 * ftrial[0] * pp->second;
        for (CI pp=ft_y.begin();pp!=ft_y.end();++pp)
          trialnormlin[pp->first] -= 0.5 * (1.0/trialnorm) * (pp->second) * 2.0 * ftrial[1] * pp->second;
        for (CI pp=ft_z.begin();pp!=ft_z.end();++pp)
          trialnormlin[pp->first] -= 0.5 * (1.0/trialnorm) * (pp->second) * 2.0 * ftrial[2] * pp->second;

        for (CI pp=maxtraclin.begin();pp!=maxtraclin.end();++pp)
          coefflin[pp->first] += frcoeff * pp->second * (1.0/trialnorm);

        for (CI pp=trialnormlin.begin();pp!=trialnormlin.end();++pp)
          coefflin[pp->first] -= frcoeff * pp->second * (1.0/(trialnorm*trialnorm)) * maxtrac;
      }


      std::map<int,std::map<int,double> >& dderiv = cnode->CoData().GetDerivDltl();

      // get sizes and iterator start
      int slavesize = (int)dderiv.size();
      std::map<int,std::map<int,double> >::iterator scurr = dderiv.begin();

      /********************************************** LinDMatrix **********/
      // loop over all DISP slave nodes in the DerivD-map of the current LM slave node
      for (int k=0;k<slavesize;++k)
      {
        int sgid = scurr->first;
        ++scurr;

        DRT::Node* snode = idiscret_->gNode(sgid);
        if (!snode) dserror("ERROR: Cannot find node with gid %",sgid);
        CoNode* csnode = dynamic_cast<CoNode*>(snode);

        // Mortar matrix D derivatives
        std::map<int,double>& thisdderiv = cnode->CoData().GetDerivDltl()[sgid];
        int mapsize = (int)(thisdderiv.size());

        // inner product D_{jk,c} * z_j for index j
        for (int prodj=0;prodj<Dim();++prodj)
        {
          int row = csnode->Dofs()[prodj];
          std::map<int,double>::iterator scolcurr = thisdderiv.begin();

          // loop over all directional derivative entries
          for (int c=0;c<mapsize;++c)
          {
            int col = scolcurr->first;
            double val = ftan[prodj] * (scolcurr->second);
            ++scolcurr;

            kteff->FEAssemble(val,row,col);
          }

          // check for completeness of DerivD-Derivatives-iteration
          if (scolcurr!=thisdderiv.end())
            dserror("ERROR: AssembleLinDM: Not all derivative entries of DerivD considered!");
        }
      }

      // Mortar matrix D and M derivatives
      std::map<int,std::map<int,double> >& mderiv = cnode->CoData().GetDerivMltl();

      // get sizes and iterator start
      int mastersize = (int)mderiv.size();
      std::map<int,std::map<int,double> >::iterator mcurr = mderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l=0;l<mastersize;++l)
      {
        int mgid = mcurr->first;
        ++mcurr;

        DRT::Node* mnode = idiscret_->gNode(mgid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",mgid);
        CoNode* cmnode = dynamic_cast<CoNode*>(mnode);

        // Mortar matrix M derivatives
        std::map<int,double>&thismderiv = cnode->CoData().GetDerivMltl()[mgid];
        int mapsize = (int)(thismderiv.size());

        // inner product M_{jl,c} * z_j for index j
        for (int prodj=0;prodj<Dim();++prodj)
        {
          int row = cmnode->Dofs()[prodj];
          std::map<int,double>::iterator mcolcurr = thismderiv.begin();

          // loop over all directional derivative entries
          for (int c=0;c<mapsize;++c)
          {
            int col = mcolcurr->first;
            double val = ftan[prodj] * (mcolcurr->second);
            ++mcolcurr;

            kteff->FEAssemble(-val,row,col);
          }

          // check for completeness of DerivM-Derivatives-iteration
          if (mcolcurr!=thismderiv.end())
            dserror("ERROR: AssembleLinDM: Not all derivative entries of DerivM considered!");
        }
      }

      // ****************************************************************
      // **************************************************** stick state
      // ****************************************************************
      if(stick)
      {
        /**************************************************** D-matrix ******/
        if ((cnode->MoData().GetDltl()).size()>0)
        {
          GEN::pairedvector<int,double>  map = cnode->MoData().GetDltl();

          for (CI p=map.begin();p!=map.end();++p)
          {
            // node id
            int gid3 = p->first;
            DRT::Node* snode = idiscret_->gNode(gid3);
            if (!snode) dserror("ERROR: Cannot find node with gid");
            CoNode* csnode = dynamic_cast<CoNode*>(snode);

            for(int dim = 0; dim<Dim();++dim)
            {
              for (CImap pp=cnode->CoData().GetDerivJumpltl()[dim].begin();pp!=cnode->CoData().GetDerivJumpltl()[dim].end();++pp)
              {
                double value[1] = {penaltytan*(p->second) * (pp->second)};
                kteff->FEAssemble(value[0],csnode->Dofs()[dim],pp->first);
              }
            }
          }
        }
        else
        {
          dserror("ERROR: no d matrix entries available for ltlt contact");
        }
        /**************************************************** D-matrix ******/
        if ((cnode->MoData().GetMltl()).size()>0)
        {
          std::map<int,double>  map = cnode->MoData().GetMltl();

          for (CImap p=map.begin();p!=map.end();++p)
          {
            // node id
            int gid3 = p->first;
            DRT::Node* snode = idiscret_->gNode(gid3);
            if (!snode) dserror("ERROR: Cannot find node with gid");
            CoNode* csnode = dynamic_cast<CoNode*>(snode);

            for(int dim = 0; dim<Dim();++dim)
            {
              for (CImap pp=cnode->CoData().GetDerivJumpltl()[dim].begin();pp!=cnode->CoData().GetDerivJumpltl()[dim].end();++pp)
              {
                double value[1] = {-penaltytan*(p->second) * (pp->second)};
                kteff->FEAssemble(value[0],csnode->Dofs()[dim],pp->first);
              }
            }
          }
        }
        else
        {
          dserror("ERROR: no m matrix entries available for ltlt contact");
        }
      }
      // ****************************************************************
      // **************************************************** slip state
      // ****************************************************************
      else
      {
        /**************************************************** D-matrix ******/
        if ((cnode->MoData().GetDltl()).size()>0)
        {
          GEN::pairedvector<int,double>  map = cnode->MoData().GetDltl();

          for (CI p=map.begin();p!=map.end();++p)
          {
            // node id
            int gid3 = p->first;
            DRT::Node* snode = idiscret_->gNode(gid3);
            if (!snode) dserror("ERROR: Cannot find node with gid");
            CoNode* csnode = dynamic_cast<CoNode*>(snode);

            for(int dim = 0; dim<Dim();++dim)
            {
              for (CImap pp=cnode->CoData().GetDerivJumpltl()[dim].begin();pp!=cnode->CoData().GetDerivJumpltl()[dim].end();++pp)
              {
                double value[1] = {penaltytan*coeff*(p->second) * (pp->second)};
                kteff->FEAssemble(value[0],csnode->Dofs()[dim],pp->first);
              }
            }
          }
        }
        else
        {
          dserror("ERROR: no d matrix entries available for ltlt contact");
        }

        /**************************************************** D-matrix ******/
        if ((cnode->MoData().GetDltl()).size()>0)
        {
          GEN::pairedvector<int,double>  map = cnode->MoData().GetDltl();

          for (CI p=map.begin();p!=map.end();++p)
          {
            // node id
            int gid3 = p->first;
            DRT::Node* snode = idiscret_->gNode(gid3);
            if (!snode) dserror("ERROR: Cannot find node with gid");
            CoNode* csnode = dynamic_cast<CoNode*>(snode);

            for(int dim = 0; dim<Dim();++dim)
            {
              for (CI pp=coefflin.begin();pp!=coefflin.end();++pp)
              {
                double value[1] = {-penaltytan*ftan[dim]*(p->second) * (pp->second)};
                kteff->FEAssemble(value[0],csnode->Dofs()[dim],pp->first);
              }
            }
          }
        }
        else
        {
          dserror("ERROR: no d matrix entries available for ltlt contact");
        }

        /**************************************************** D-matrix ******/
        if ((cnode->MoData().GetMltl()).size()>0)
        {
          std::map<int,double>  map = cnode->MoData().GetMltl();

          for (CImap p=map.begin();p!=map.end();++p)
          {
            // node id
            int gid3 = p->first;
            DRT::Node* snode = idiscret_->gNode(gid3);
            if (!snode) dserror("ERROR: Cannot find node with gid");
            CoNode* csnode = dynamic_cast<CoNode*>(snode);

            for(int dim = 0; dim<Dim();++dim)
            {
              for (CImap pp=cnode->CoData().GetDerivJumpltl()[dim].begin();pp!=cnode->CoData().GetDerivJumpltl()[dim].end();++pp)
              {
                double value[1] = {-penaltytan*coeff*(p->second) * (pp->second)};
                kteff->FEAssemble(value[0],csnode->Dofs()[dim],pp->first);
              }
            }
          }
        }
        else
        {
          dserror("ERROR: no m matrix entries available for ltlt contact");
        }
        if ((cnode->MoData().GetMltl()).size()>0)
        {
          std::map<int,double>  map = cnode->MoData().GetMltl();

          for (CImap p=map.begin();p!=map.end();++p)
          {
            // node id
            int gid3 = p->first;
            DRT::Node* snode = idiscret_->gNode(gid3);
            if (!snode) dserror("ERROR: Cannot find node with gid");
            CoNode* csnode = dynamic_cast<CoNode*>(snode);

            for(int dim = 0; dim<Dim();++dim)
            {
              for (CI pp=coefflin.begin();pp!=coefflin.end();++pp)
              {
                double value[1] = {penaltytan*ftan[dim]*(p->second) * (pp->second)};
                kteff->FEAssemble(value[0],csnode->Dofs()[dim],pp->first);
              }
            }
          }
        }
        else
        {
          dserror("ERROR: no m matrix entries available for ltlt contact");
        }
      }

    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Add nts penalty forces master                           farah 11/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AddNTSforcesMaster(Teuchos::RCP<Epetra_FEVector> feff)
{
  const double penalty = IParams().get<double>("PENALTYPARAM");

  typedef GEN::pairedvector<int,double>::const_iterator CI;
  typedef std::map<int,double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j=0;j<mnoderowmap_->NumMyElements();++j)
  {
    int gid = mnoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    // only for corners
    if(!cnode->IsOnCorner())
      continue;

    // is gap is in contact
    if(cnode->CoData().Getgnts()<1e-12)
    {
      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetDnts()).size()>0)
      {
        GEN::pairedvector<int,double>  map = cnode->MoData().GetDnts();

        for (CI p=map.begin();p!=map.end();++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("ERROR: Cannot find node with gid");
          CoNode* csnode = dynamic_cast<CoNode*>(snode);

          for(int dim = 0; dim<Dim();++dim)
          {
            double value[1] = {penalty*(p->second) * cnode->CoData().Getgnts()*cnode->MoData().n()[dim]};
            const int ltlid[1] = {csnode->Dofs()[dim]};
            int err = feff->SumIntoGlobalValues(1,ltlid,value);
            if (err<0)
              dserror("stop");
          }
        }
      }
      else
      {
        dserror("ERROR: no d matrix entries available for ltlt contact");
      }

      /**************************************************** M-matrix ******/
      if ((cnode->MoData().GetMnts()).size()>0)
      {
        std::map<int,double>  map = cnode->MoData().GetMnts();

        for (CImap p=map.begin();p!=map.end();++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("ERROR: Cannot find node with gid");
          CoNode* csnode = dynamic_cast<CoNode*>(snode);

          for(int dim = 0; dim<Dim();++dim)
          {
            double value[1] = {-penalty*(p->second) * cnode->CoData().Getgnts()*cnode->MoData().n()[dim]};
            const int ltlid[1] = {csnode->Dofs()[dim]};
            int err = feff->SumIntoGlobalValues(1,ltlid,value);
            if (err<0)
              dserror("stop");
          }
        }
      }
      else
      {
        dserror("ERROR: no m matrix entries available for ltlt contact");
      }
    }

  }
  return;
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces master                  farah 11/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AddLTSforcesMaster(Teuchos::RCP<Epetra_FEVector> feff)
{
  const double penalty = IParams().get<double>("PENALTYPARAM");

  typedef GEN::pairedvector<int,double>::const_iterator CI;
  typedef std::map<int,double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j=0;j<mnoderowmap_->NumMyElements();++j)
  {
    int gid = mnoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    // only for edges (without corner)
    if(!cnode->IsOnCornerEdge())
      continue;
    if(cnode->IsOnCorner())
      continue;

    // scale penalty
    double penaltyLts = penalty * cnode->CoData().Kappa();

    // is gap is in contact
    if(cnode->CoData().Getglts()<1e-12)
    {
      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetDlts()).size()>0)
      {
        GEN::pairedvector<int,double>  map = cnode->MoData().GetDlts();

        for (CI p=map.begin();p!=map.end();++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("ERROR: Cannot find node with gid");
          CoNode* csnode = dynamic_cast<CoNode*>(snode);

          for(int dim = 0; dim<Dim();++dim)
          {
            double value[1] = {penaltyLts*(p->second) * cnode->CoData().Getglts()*cnode->MoData().n()[dim]};
            const int ltlid[1] = {csnode->Dofs()[dim]};
            int err = feff->SumIntoGlobalValues(1,ltlid,value);
            if (err<0)
              dserror("stop");
          }
        }
      }
      else
      {
        dserror("ERROR: no d matrix entries available for ltlt contact");
      }

      /**************************************************** M-matrix ******/
      if ((cnode->MoData().GetMlts()).size()>0)
      {
        std::map<int,double>  map = cnode->MoData().GetMlts();

        for (CImap p=map.begin();p!=map.end();++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("ERROR: Cannot find node with gid");
          CoNode* csnode = dynamic_cast<CoNode*>(snode);

          for(int dim = 0; dim<Dim();++dim)
          {
            double value[1] = {-penaltyLts*(p->second) * cnode->CoData().Getglts()*cnode->MoData().n()[dim]};
            const int ltlid[1] = {csnode->Dofs()[dim]};
            int err = feff->SumIntoGlobalValues(1,ltlid,value);
            if (err<0)
              dserror("stop");
          }
        }
      }
      else
      {
        dserror("ERROR: no m matrix entries available for ltlt contact");
      }
    }

  }
  return;
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces                         farah 10/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AddLTLforces(Teuchos::RCP<Epetra_FEVector> feff)
{
  double slaveforce  = 0.0;
  double masterforce = 0.0;

  // gap = g_n * n
  // D/M = sval/mval
  const double penalty = IParams().get<double>("PENALTYPARAM");

  typedef GEN::pairedvector<int,double>::const_iterator CI;
  typedef std::map<int,double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j=0;j<snoderowmap_->NumMyElements();++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    // check if this is valid node
    if(cnode->CoData().Getgltl()[0]<1e8 and
       cnode->CoData().Getgltl()[1]<1e8 and
       cnode->CoData().Getgltl()[2]<1e8)
    {
      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetDltl()).size()>0)
      {
        GEN::pairedvector<int,double>  map = cnode->MoData().GetDltl();

        for (CI p=map.begin();p!=map.end();++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("ERROR: Cannot find node with gid");
          CoNode* csnode = dynamic_cast<CoNode*>(snode);

          for(int dim = 0; dim<Dim();++dim)
          {
            double value[1] = {penalty*(p->second) * cnode->CoData().Getgltl()[dim]};
            const int ltlid[1] = {csnode->Dofs()[dim]};
            int err = feff->SumIntoGlobalValues(1,ltlid,value);
            if (err<0)
              dserror("stop");
            slaveforce +=value[0];
          }
        }
      }
      else
      {
        dserror("ERROR: no d matrix entries available for ltlt contact");
      }

      /**************************************************** M-matrix ******/
      if ((cnode->MoData().GetMltl()).size()>0)
      {
        std::map<int,double>  map = cnode->MoData().GetMltl();

        for (CImap p=map.begin();p!=map.end();++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("ERROR: Cannot find node with gid");
          CoNode* csnode = dynamic_cast<CoNode*>(snode);

          for(int dim = 0; dim<Dim();++dim)
          {
            double value[1] = {-penalty*(p->second) * cnode->CoData().Getgltl()[dim]};
            const int ltlid[1] = {csnode->Dofs()[dim]};
            int err = feff->SumIntoGlobalValues(1,ltlid,value);
            if (err<0)
              dserror("stop");
            masterforce +=value[0];
          }
        }
      }
      else
      {
        dserror("ERROR: no m matrix entries available for ltlt contact");
      }
    }
  }

//  std::cout << "LTL slaveforce  = " << slaveforce << std::endl;
//  std::cout << "LTL masterforce = " << masterforce << std::endl;
//  std::cout << "LTL sum         = " << slaveforce+masterforce << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces                         farah 11/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AddLTSstiffnessMaster(Teuchos::RCP<LINALG::SparseMatrix> kteff)
{
  const double penalty = IParams().get<double>("PENALTYPARAM");

  typedef GEN::pairedvector<int,double>::const_iterator CI;
  typedef std::map<int,double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j=0;j<mnoderowmap_->NumMyElements();++j)
  {
    int gid = mnoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    // only for edges (without corner)
    if(!cnode->IsOnCornerEdge())
      continue;
    if(cnode->IsOnCorner())
      continue;

    // scale penalty
    double penaltyLts = penalty * cnode->CoData().Kappa();

    // is gap is in contact
    if(cnode->CoData().Getglts()<1e-12)
    {
      double lm[3] = {0.0,0.0,0.0};
      lm[0] = penaltyLts * cnode->CoData().Getglts() * cnode->MoData().n()[0];
      lm[1] = penaltyLts * cnode->CoData().Getglts() * cnode->MoData().n()[1];
      lm[2] = penaltyLts * cnode->CoData().Getglts() * cnode->MoData().n()[2];

      std::map<int,std::map<int,double> >& dderiv = cnode->CoData().GetDerivDlts();

      // get sizes and iterator start
      int slavesize = (int)dderiv.size();
      std::map<int,std::map<int,double> >::iterator scurr = dderiv.begin();

      /********************************************** LinDMatrix **********/
      // loop over all DISP slave nodes in the DerivD-map of the current LM slave node
      for (int k=0;k<slavesize;++k)
      {
        int sgid = scurr->first;
        ++scurr;

        DRT::Node* snode = idiscret_->gNode(sgid);
        if (!snode) dserror("ERROR: Cannot find node with gid %",sgid);
        CoNode* csnode = dynamic_cast<CoNode*>(snode);

        // Mortar matrix D derivatives
        std::map<int,double>& thisdderiv = cnode->CoData().GetDerivDlts()[sgid];
        int mapsize = (int)(thisdderiv.size());

        // inner product D_{jk,c} * z_j for index j
        for (int prodj=0;prodj<Dim();++prodj)
        {
          int row = csnode->Dofs()[prodj];
          std::map<int,double>::iterator scolcurr = thisdderiv.begin();

          // loop over all directional derivative entries
          for (int c=0;c<mapsize;++c)
          {
            int col = scolcurr->first;
            double val = lm[prodj] * (scolcurr->second);
            ++scolcurr;

            kteff->FEAssemble(-val,row,col);
          }

          // check for completeness of DerivD-Derivatives-iteration
          if (scolcurr!=thisdderiv.end())
            dserror("ERROR: AssembleLinDM: Not all derivative entries of DerivD considered!");
        }
      }

      // Mortar matrix D and M derivatives
      std::map<int,std::map<int,double> >& mderiv = cnode->CoData().GetDerivMlts();

      // get sizes and iterator start
      int mastersize = (int)mderiv.size();
      std::map<int,std::map<int,double> >::iterator mcurr = mderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l=0;l<mastersize;++l)
      {
        int mgid = mcurr->first;
        ++mcurr;

        DRT::Node* mnode = idiscret_->gNode(mgid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",mgid);
        CoNode* cmnode = dynamic_cast<CoNode*>(mnode);

        // Mortar matrix M derivatives
        std::map<int,double>&thismderiv = cnode->CoData().GetDerivMlts()[mgid];
        int mapsize = (int)(thismderiv.size());

        // inner product M_{jl,c} * z_j for index j
        for (int prodj=0;prodj<Dim();++prodj)
        {
          int row = cmnode->Dofs()[prodj];
          std::map<int,double>::iterator mcolcurr = thismderiv.begin();

          // loop over all directional derivative entries
          for (int c=0;c<mapsize;++c)
          {
            int col = mcolcurr->first;
            double val = lm[prodj] * (mcolcurr->second);
            ++mcolcurr;

            kteff->FEAssemble(val,row,col);
          }

          // check for completeness of DerivM-Derivatives-iteration
          if (mcolcurr!=thismderiv.end())
            dserror("ERROR: AssembleLinDM: Not all derivative entries of DerivM considered!");
        }
      }

      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetDlts()).size()>0)
      {
        GEN::pairedvector<int,double>  map = cnode->MoData().GetDlts();

        for (CI p=map.begin();p!=map.end();++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("ERROR: Cannot find node with gid");
          CoNode* csnode = dynamic_cast<CoNode*>(snode);

          for(int dim = 0; dim<Dim();++dim)
          {
            // gap linearization
            for (CImap pp=cnode->CoData().GetDerivGlts().begin();pp!=cnode->CoData().GetDerivGlts().end();++pp)
            {
              double value[1] = {-penaltyLts*(p->second) * (pp->second) * cnode->MoData().n()[dim]};
              kteff->FEAssemble(value[0],csnode->Dofs()[dim],pp->first);
            }
            // normal linearization
            for (CI pp=cnode->CoData().GetDerivN()[dim].begin();pp!=cnode->CoData().GetDerivN()[dim].end();++pp)
            {
              double value[1] = {-penaltyLts*(p->second) * (pp->second) * cnode->CoData().Getglts()};
              kteff->FEAssemble(value[0],csnode->Dofs()[dim],pp->first);
            }
          }
        }
      }
      else
      {
        dserror("ERROR: no d matrix entries available for ltlt contact");
      }
      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetMlts()).size()>0)
      {
        std::map<int,double>  map = cnode->MoData().GetMlts();

        for (CImap p=map.begin();p!=map.end();++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("ERROR: Cannot find node with gid");
          CoNode* csnode = dynamic_cast<CoNode*>(snode);

          for(int dim = 0; dim<Dim();++dim)
          {
            for (CImap pp=cnode->CoData().GetDerivGlts().begin();pp!=cnode->CoData().GetDerivGlts().end();++pp)
            {
              double value[1] = {penaltyLts*(p->second) * (pp->second)* cnode->MoData().n()[dim]};
              kteff->FEAssemble(value[0],csnode->Dofs()[dim],pp->first);
            }
            // normal linearization
            for (CI pp=cnode->CoData().GetDerivN()[dim].begin();pp!=cnode->CoData().GetDerivN()[dim].end();++pp)
            {
              double value[1] = {penaltyLts*(p->second) * (pp->second)* cnode->CoData().Getglts()};
              kteff->FEAssemble(value[0],csnode->Dofs()[dim],pp->first);
            }
          }
        }
      }
      else
      {
        dserror("ERROR: no m matrix entries available for ltlt contact");
      }

    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces                         farah 11/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AddNTSstiffnessMaster(Teuchos::RCP<LINALG::SparseMatrix> kteff)
{
  const double penalty = IParams().get<double>("PENALTYPARAM");

  typedef GEN::pairedvector<int,double>::const_iterator CI;
  typedef std::map<int,double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j=0;j<mnoderowmap_->NumMyElements();++j)
  {
    int gid = mnoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    // only for corners
    if(!cnode->IsOnCorner())
      continue;

    // is gap is in contact
    if(cnode->CoData().Getgnts()<1e-12)
    {
      double lm[3] = {0.0,0.0,0.0};
      lm[0] = penalty * cnode->CoData().Getgnts() * cnode->MoData().n()[0];
      lm[1] = penalty * cnode->CoData().Getgnts() * cnode->MoData().n()[1];
      lm[2] = penalty * cnode->CoData().Getgnts() * cnode->MoData().n()[2];

      // Mortar matrix D and M derivatives
      std::map<int,std::map<int,double> >& mderiv = cnode->CoData().GetDerivMnts();

      // get sizes and iterator start
      int mastersize = (int)mderiv.size();
      std::map<int,std::map<int,double> >::iterator mcurr = mderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l=0;l<mastersize;++l)
      {
        int mgid = mcurr->first;
        ++mcurr;

        DRT::Node* mnode = idiscret_->gNode(mgid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",mgid);
        CoNode* cmnode = dynamic_cast<CoNode*>(mnode);

        // Mortar matrix M derivatives
        std::map<int,double>&thismderiv = cnode->CoData().GetDerivMnts()[mgid];
        int mapsize = (int)(thismderiv.size());

        // inner product M_{jl,c} * z_j for index j
        for (int prodj=0;prodj<Dim();++prodj)
        {
          int row = cmnode->Dofs()[prodj];
          std::map<int,double>::iterator mcolcurr = thismderiv.begin();

          // loop over all directional derivative entries
          for (int c=0;c<mapsize;++c)
          {
            int col = mcolcurr->first;
            double val = lm[prodj] * (mcolcurr->second);
            ++mcolcurr;

            kteff->FEAssemble(val,row,col);
          }

          // check for completeness of DerivM-Derivatives-iteration
          if (mcolcurr!=thismderiv.end())
            dserror("ERROR: AssembleLinDM: Not all derivative entries of DerivM considered!");
        }
      }

      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetDnts()).size()>0)
      {
        GEN::pairedvector<int,double>  map = cnode->MoData().GetDnts();

        for (CI p=map.begin();p!=map.end();++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("ERROR: Cannot find node with gid");
          CoNode* csnode = dynamic_cast<CoNode*>(snode);

          for(int dim = 0; dim<Dim();++dim)
          {
            // gap linearization
            for (CImap pp=cnode->CoData().GetDerivGnts().begin();pp!=cnode->CoData().GetDerivGnts().end();++pp)
            {
              double value[1] = {-penalty*(p->second) * (pp->second) * cnode->MoData().n()[dim]};
              kteff->FEAssemble(value[0],csnode->Dofs()[dim],pp->first);
            }
            // normal linearization
            for (CI pp=cnode->CoData().GetDerivN()[dim].begin();pp!=cnode->CoData().GetDerivN()[dim].end();++pp)
            {
              double value[1] = {-penalty*(p->second) * (pp->second) * cnode->CoData().Getgnts()};
              kteff->FEAssemble(value[0],csnode->Dofs()[dim],pp->first);
            }
          }
        }
      }
      else
      {
        dserror("ERROR: no d matrix entries available for ltlt contact");
      }
      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetMnts()).size()>0)
      {
        std::map<int,double>  map = cnode->MoData().GetMnts();

        for (CImap p=map.begin();p!=map.end();++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("ERROR: Cannot find node with gid");
          CoNode* csnode = dynamic_cast<CoNode*>(snode);

          for(int dim = 0; dim<Dim();++dim)
          {
            for (CImap pp=cnode->CoData().GetDerivGnts().begin();pp!=cnode->CoData().GetDerivGnts().end();++pp)
            {
              double value[1] = {penalty*(p->second) * (pp->second)* cnode->MoData().n()[dim]};
              kteff->FEAssemble(value[0],csnode->Dofs()[dim],pp->first);
            }
            // normal linearization
            for (CI pp=cnode->CoData().GetDerivN()[dim].begin();pp!=cnode->CoData().GetDerivN()[dim].end();++pp)
            {
              double value[1] = {penalty*(p->second) * (pp->second)* cnode->CoData().Getgnts()};
              kteff->FEAssemble(value[0],csnode->Dofs()[dim],pp->first);
            }
          }
        }
      }
      else
      {
        dserror("ERROR: no m matrix entries available for ltlt contact");
      }

    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Add line to line penalty forces                         farah 10/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AddLTLstiffness(Teuchos::RCP<LINALG::SparseMatrix> kteff)
{
  const double penalty = IParams().get<double>("PENALTYPARAM");

  typedef GEN::pairedvector<int,double>::const_iterator CI;
  typedef std::map<int,double>::const_iterator CImap;

  // loop over all slave nodes
  for (int j=0;j<snoderowmap_->NumMyElements();++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    // check if this is valid node
    if(cnode->CoData().Getgltl()[0]<1e8 and
       cnode->CoData().Getgltl()[1]<1e8 and
       cnode->CoData().Getgltl()[2]<1e8)
    {
      double lm[3] = {0.0,0.0,0.0};
      lm[0] = penalty * cnode->CoData().Getgltl()[0];
      lm[1] = penalty * cnode->CoData().Getgltl()[1];
      lm[2] = penalty * cnode->CoData().Getgltl()[2];

      std::map<int,std::map<int,double> >& dderiv = cnode->CoData().GetDerivDltl();

      // get sizes and iterator start
      int slavesize = (int)dderiv.size();
      std::map<int,std::map<int,double> >::iterator scurr = dderiv.begin();

      /********************************************** LinDMatrix **********/
      // loop over all DISP slave nodes in the DerivD-map of the current LM slave node
      for (int k=0;k<slavesize;++k)
      {
        int sgid = scurr->first;
        ++scurr;

        DRT::Node* snode = idiscret_->gNode(sgid);
        if (!snode) dserror("ERROR: Cannot find node with gid %",sgid);
        CoNode* csnode = dynamic_cast<CoNode*>(snode);

        // Mortar matrix D derivatives
        std::map<int,double>& thisdderiv = cnode->CoData().GetDerivDltl()[sgid];
        int mapsize = (int)(thisdderiv.size());

        // inner product D_{jk,c} * z_j for index j
        for (int prodj=0;prodj<Dim();++prodj)
        {
          int row = csnode->Dofs()[prodj];
          std::map<int,double>::iterator scolcurr = thisdderiv.begin();

          // loop over all directional derivative entries
          for (int c=0;c<mapsize;++c)
          {
            int col = scolcurr->first;
            double val = lm[prodj] * (scolcurr->second);
            ++scolcurr;

            kteff->FEAssemble(-val,row,col);
          }

          // check for completeness of DerivD-Derivatives-iteration
          if (scolcurr!=thisdderiv.end())
            dserror("ERROR: AssembleLinDM: Not all derivative entries of DerivD considered!");
        }
      }

      // Mortar matrix D and M derivatives
      std::map<int,std::map<int,double> >& mderiv = cnode->CoData().GetDerivMltl();

      // get sizes and iterator start
      int mastersize = (int)mderiv.size();
      std::map<int,std::map<int,double> >::iterator mcurr = mderiv.begin();

      /********************************************** LinMMatrix **********/
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l=0;l<mastersize;++l)
      {
        int mgid = mcurr->first;
        ++mcurr;

        DRT::Node* mnode = idiscret_->gNode(mgid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",mgid);
        CoNode* cmnode = dynamic_cast<CoNode*>(mnode);

        // Mortar matrix M derivatives
        std::map<int,double>&thismderiv = cnode->CoData().GetDerivMltl()[mgid];
        int mapsize = (int)(thismderiv.size());

        // inner product M_{jl,c} * z_j for index j
        for (int prodj=0;prodj<Dim();++prodj)
        {
          int row = cmnode->Dofs()[prodj];
          std::map<int,double>::iterator mcolcurr = thismderiv.begin();

          // loop over all directional derivative entries
          for (int c=0;c<mapsize;++c)
          {
            int col = mcolcurr->first;
            double val = lm[prodj] * (mcolcurr->second);
            ++mcolcurr;

            kteff->FEAssemble(val,row,col);
          }

          // check for completeness of DerivM-Derivatives-iteration
          if (mcolcurr!=thismderiv.end())
            dserror("ERROR: AssembleLinDM: Not all derivative entries of DerivM considered!");
        }
      }

      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetDltl()).size()>0)
      {
        GEN::pairedvector<int,double>  map = cnode->MoData().GetDltl();

        for (CI p=map.begin();p!=map.end();++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("ERROR: Cannot find node with gid");
          CoNode* csnode = dynamic_cast<CoNode*>(snode);

          for(int dim = 0; dim<Dim();++dim)
          {
            for (CImap pp=cnode->CoData().GetDerivGltl()[dim].begin();pp!=cnode->CoData().GetDerivGltl()[dim].end();++pp)
            {
              double value[1] = {-penalty*(p->second) * (pp->second)};
              kteff->FEAssemble(value[0],csnode->Dofs()[dim],pp->first);
            }
          }
        }
      }
      else
      {
        dserror("ERROR: no d matrix entries available for ltlt contact");
      }
      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetMltl()).size()>0)
      {
        std::map<int,double>  map = cnode->MoData().GetMltl();

        for (CImap p=map.begin();p!=map.end();++p)
        {
          // node id
          int gid3 = p->first;
          DRT::Node* snode = idiscret_->gNode(gid3);
          if (!snode) dserror("ERROR: Cannot find node with gid");
          CoNode* csnode = dynamic_cast<CoNode*>(snode);

          for(int dim = 0; dim<Dim();++dim)
          {
            for (CImap pp=cnode->CoData().GetDerivGltl()[dim].begin();pp!=cnode->CoData().GetDerivGltl()[dim].end();++pp)
            {
              double value[1] = {penalty*(p->second) * (pp->second)};
              kteff->FEAssemble(value[0],csnode->Dofs()[dim],pp->first);
            }
          }
        }
      }
      else
      {
        dserror("ERROR: no m matrix entries available for ltlt contact");
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  post evaluate to scale calculated terms                 farah 02/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::PostEvaluate(const int step, const int iter)
{
  // nonsmooth contact
  if(nonSmoothContact_)
  {
    // store lts into mortar data container
    StoreLTSvalues();

    // store nts into mortar data container
    StoreNTSvalues();

    // FINAL TEST:



  }
  // node assemble of terms
  else
  {
    // decide which type of coupling should be evaluated
    INPAR::MORTAR::AlgorithmType algo =
        DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(imortar_, "ALGORITHM");

    // smooth contact
    switch(algo)
    {
    //*********************************
    // Mortar Coupling (STS)    (2D/3D)
    // Gauß-Point-To-Segment (GPTS)
    //*********************************
    case INPAR::MORTAR::algorithm_mortar:
    case INPAR::MORTAR::algorithm_gpts:
    {
      // already stored
      return;
      break;
    }
    //*********************************
    // Line-to-Segment Coupling (3D)
    //*********************************
    case INPAR::MORTAR::algorithm_lts:
    {
      // store lts into mortar data container
      StoreLTSvalues();
      break;
    }
    //*********************************
    // Node-to-Segment Coupling (2D/3D)
    //*********************************
    case INPAR::MORTAR::algorithm_nts:
    {
      // store nts into mortar data container
      StoreNTSvalues();
      break;
    }
    //*********************************
    // line-to-line Coupling (3D)
    //*********************************
    case INPAR::MORTAR::algorithm_ltl:
    {
      return;
      break;
    }
    //*********************************
    // Node-to-Line Coupling (3D)
    //*********************************
    case INPAR::MORTAR::algorithm_ntl:
    {
      dserror("ERROR: not yet implemented!");
      break;
    }
    //*********************************
    // Segment-to-Line Coupling (3D)
    //*********************************
    case INPAR::MORTAR::algorithm_stl:
    {
      // store lts into mortar data container
      StoreLTSvalues();
      break;
    }
    //*********************************
    // Default case
    //*********************************
    default:
    {
      dserror("ERROR: Unknown discr. type for constraints!");
      break;
    }
    }
  }

#ifdef MORTARGMSHCELLS
  // finish integration cell GMSH files
  int proc = Comm().MyPID();
  std::ostringstream filename;
  filename << "o/gmsh_output/cells_" << proc << ".pos";
  FILE* fp = fopen(filename.str().c_str(), "a");
  std::stringstream gmshfilecontent2;
  gmshfilecontent2 << "};" << std::endl;
  fprintf(fp,gmshfilecontent2.str().c_str());
  fclose(fp);

  // construct unique filename for gmsh output
  // first index = time step index
  std::ostringstream newfilename;
  newfilename << "o/gmsh_output/cells_";
  if (step<10) newfilename << 0 << 0 << 0 << 0;
  else if (step<100) newfilename << 0 << 0 << 0;
  else if (step<1000) newfilename << 0 << 0;
  else if (step<10000) newfilename << 0;
  else if (step>99999) dserror("Gmsh output implemented for a maximum of 99.999 time steps");
  newfilename << step;

  // second index = Newton iteration index
  newfilename << "_";
  if (iter<10) newfilename << 0;
  else if (iter>99) dserror("Gmsh output implemented for a maximum of 99 iterations");
  newfilename << iter << "_p" << proc << ".pos";

  // rename file
  rename(filename.str().c_str(),newfilename.str().c_str());
#endif // #ifdef MORTARGMSHCELLS

  return;
}

/*----------------------------------------------------------------------*
 |  scale calculated entries for nts, mortar etc...         farah 02/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::ScaleTerms()
{
  // only for 3d contact
  if(Dim()==2)
    return;

  // scale ltl terms for point and line contact
  ScaleTermsLTL();

  // bye
  return;
}

/*----------------------------------------------------------------------*
 |  scale calculated entries for lts contact                farah 09/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::ScaleTermsLTL()
{
  std::cout << "ScaleTerms LTL" << std::endl;
  dserror("ERROR: ScaleTermsLTL() is outdated");

  return;
//  // create iterators for data types
//  typedef GEN::pairedvector<int,double>::const_iterator CI;
//  typedef std::map<int,double>::const_iterator          CImap;
//
//  // loop over all possibly non smooth nodes
//  for(int i=0; i<snoderowmap_->NumMyElements();++i)
//  {
//    int gid = snoderowmap_->GID(i);
//    DRT::Node* node = idiscret_->gNode(gid);
//    if (!node)
//      dserror("ERROR: Cannot find node with gid %",gid);
//    CoNode* cnode = dynamic_cast<CoNode*>(node);
//
//    // only for edge nodes
//    if(!cnode->IsOnEdge())
//      continue;
//
//    // get scale factor alpha
//    const double alpha = cnode->CoData().GetAlphaN();
//    if(alpha<-1e-12)
//      continue;
//
//    std::cout << "PERFORM SCALING" << std::endl;
//
//    // check if integration is done
////    if (cnode->MoData().GetD().size()<1)
////      continue;
//
//    // ############################################################
//    //                   create aux terms
//    // ############################################################
//    std::map<int, double> derivg_lts =
//        cnode->CoData().GetDerivGlts();
//
//    std::map<int, std::map<int, double> >derivm_lts =
//        cnode->CoData().GetDerivMlts();
//
//    std::map<int, std::map<int, double> > derivd_lts =
//        cnode->CoData().GetDerivDlts();
//
//    double g_lts = cnode->CoData().Getglts();
//
//    std::map<int, double> m_lts =
//        cnode->MoData().GetMlts();
//
//    GEN::pairedvector<int, double> d_lts =
//        cnode->MoData().GetDlts();
//
//    // ############################################################
//    //                   scale lts terms:
//    // ############################################################
//    {
//      //-------------------------------------------------------------------------------------
//      // scale M matrix entries for mortar data
//      for (CImap p=m_lts.begin();p!=m_lts.end();++p)
//        cnode->MoData().GetMlts()[p->first] = alpha*(p->second);
//
//      //-------------------------------------------------------------------------------------
//      // scale weighted gap
//      cnode->CoData().Getglts() = alpha * g_lts;
//
//      //-------------------------------------------------------------------------------------
//      // scale M deriv
//      {
//        // Mortar M derivatives
//        std::map<int,std::map<int,double> >& mderiv = cnode->CoData().GetDerivMlts();
//
//        // get sizes and iterator start
//        int mastersize = (int)mderiv.size();
//        std::map<int,std::map<int,double> >::iterator mcurr = mderiv.begin();
//
//        /********************************************** LinMMatrix **********/
//        // loop over all master nodes in the DerivM-map of the current LM slave node
//        for (int l=0;l<mastersize;++l)
//        {
//          int mgid = mcurr->first;
//          ++mcurr;
//
//          // Mortar matrix M derivatives
//          std::map<int,double>&thismderiv = cnode->CoData().GetDerivMlts()[mgid];
//          std::map<int,double>&auxderiv   = derivm_lts[mgid];
//
//          int mapsize = (int)(thismderiv.size());
//
//          // inner product M_{jl,c} * z_j for index j
//          for (int prodj=0;prodj<Dim();++prodj)
//          {
//            std::map<int,double>::iterator mcolcurr    = thismderiv.begin();
//            std::map<int,double>::iterator mcolcurraux = auxderiv.begin();
//
//            // loop over all directional derivative entries
//            for (int c=0;c<mapsize;++c)
//            {
//              mcolcurr->second = alpha * (mcolcurraux->second);
//              ++mcolcurr;
//            }
//
//            // check for completeness of DerivM-Derivatives-iteration
//            if (mcolcurr!=thismderiv.end())
//              dserror("ERROR: ScaleTerms: Not all derivative entries of DerivM considered!");
//          }
//        }
//      }
//
//      //-------------------------------------------------------------------------------------
//      // scale alpha lin
//      {
//        for (CI p=cnode->CoData().GetAlpha().begin();p!=cnode->CoData().GetAlpha().end();++p)
//        {
//          for (CImap pp=m_lts.begin();pp!=m_lts.end();++pp)
//          {
//            (cnode->CoData().GetDerivMlts()[pp->first])[p->first] += (p->second) * (pp->second);
//          }
//        }
//      }
//
//      //-------------------------------------------------------------------------------------
//      // scale gap derivative entries for mortar data
//      for (CImap p=derivg_lts.begin();p!=derivg_lts.end();++p)
//        cnode->CoData().GetDerivGlts()[p->first] = (p->second)*alpha;
//
//      //-------------------------------------------------------------------------------------
//      // scale gap derivative with lin alpha entries for mortar data
//      for (CI p=cnode->CoData().GetAlpha().begin();p!=cnode->CoData().GetAlpha().end();++p)
//      {
//        cnode->CoData().GetDerivGlts()[p->first] += g_lts *(p->second);
//      }
//    }// end lts terms
//
//    // ############################################################
//    // scale ltl terms and add to lts terms:
//    // ############################################################
//    {
//      // get scale factor alpha
//      const double alphaNew = 1.0 - cnode->CoData().GetAlphaN();
//
//      if(cnode->MoData().GetDlts().size()>1)
//        dserror("ERROR: D matrix must be diagonal!");
//
//      //-------------------------------------------------------------------------------------
//      const double dvalue = d_lts.begin()->second;
//
//      //-------------------------------------------------------------------------------------
//      // scale M matrix entries for mortar data --> Ntilde * D * (1-a)
//      for (CImap p=cnode->MoData().GetMltl().begin();p!=cnode->MoData().GetMltl().end();++p)
//        cnode->MoData().GetMlts()[p->first] += alphaNew*(p->second)*dvalue;
//
//      //-------------------------------------------------------------------------------------
//      // scale weighted gap
//      cnode->CoData().Getglts() += alphaNew *  cnode->CoData().Getgltl() * dvalue;
//
//      //-------------------------------------------------------------------------------------
//      // Mortar matrix D and M derivatives
//      std::map<int,std::map<int,double> >& dderiv = cnode->CoData().GetDerivDlts();
//
//      // get sizes and iterator start
//      int slavesize = (int)dderiv.size();
////      if (slavesize>1)
////        dserror("wrong dlin dimension");
//
//      std::map<int,std::map<int,double> >::iterator scurr = dderiv.begin();
//
//      /********************************************** LinDMatrix **********/
//      // loop over all DISP slave nodes in the DerivD-map of the current LM slave node
//      for (int k=0;k<slavesize;++k)
//      {
//        int sgid = scurr->first;
//        ++scurr;
//
//        DRT::Node* snode = idiscret_->gNode(sgid);
//        if (!snode) dserror("ERROR: Cannot find node with gid %",sgid);
//
//        // Mortar matrix D derivatives
//        std::map<int,double>& thisdderivD        = derivd_lts[sgid];
//        std::map<int,double>& thismderivmortarM  = cnode->CoData().GetDerivMlts()[sgid];
//
//        int mapsize = (int)(thisdderivD.size());
//
//        // inner product D_{jk,c} * z_j for index j
////        for (int prodj=0;prodj<Dim();++prodj)
////        {
//          std::map<int,double>::iterator scolcurrD = thisdderivD.begin();
//
//          // loop over all directional derivative entries
//          for (int c=0;c<mapsize;++c)
//          {
//            // derivG
//            cnode->CoData().GetDerivGlts()[scolcurrD->first] += alphaNew*(scolcurrD->second)*cnode->CoData().Getgltl();
//
//            //derivM with delta D
//            for (CImap p=cnode->MoData().GetMltl().begin();p!=cnode->MoData().GetMltl().end();++p)
//              thismderivmortarM[scolcurrD->first] += alphaNew*(scolcurrD->second) * p->second;
//            ++scolcurrD;
//          }
//
//          // check for completeness of DerivD-Derivatives-iteration
//          if (scolcurrD!=thisdderivD.end())
//            dserror("ERROR: ScaleTerms: Not all derivative entries of DerivD considered!");
//        //}
//      }
//
//
//      //-------------------------------------------------------------------------------------
//      // scale M LTL deriv
//      {
//        // Mortar M derivatives
//        std::map<int,std::map<int,double> >& mntsderiv = cnode->CoData().GetDerivMltl();
//
//        // get sizes and iterator start
//        int mastersize = (int)mntsderiv.size();
//        std::map<int,std::map<int,double> >::iterator mntscurr = mntsderiv.begin();
//
//        /********************************************** LinMMatrix **********/
//        // loop over all master nodes in the DerivM-map of the current LM slave node
//        for (int l=0;l<mastersize;++l)
//        {
//          int mgid = mntscurr->first;
//          ++mntscurr;
//
//          // Mortar matrix M derivatives
//          std::map<int,double>&thismderivnts    = cnode->CoData().GetDerivMltl()[mgid];
//          std::map<int,double>&thismderivmortar = cnode->CoData().GetDerivMlts()[mgid];
//
//          int mapsize = (int)(thismderivnts.size());
//
//          // inner product M_{jl,c} * z_j for index j
////          for (int prodj=0;prodj<Dim();++prodj)
////          {
//            std::map<int,double>::iterator mcolcurr = thismderivnts.begin();
//
//            // loop over all directional derivative entries
//            for (int c=0;c<mapsize;++c)
//            {
//              thismderivmortar[mcolcurr->first] += alphaNew * (mcolcurr->second) * dvalue;
//              ++mcolcurr;
//            }
//
//            // check for completeness of DerivM-Derivatives-iteration
//            if (mcolcurr!=thismderivnts.end())
//              dserror("ERROR: ScaleTerms: Not all derivative entries of DerivM considered!");
////          }
//        }
//      }
//
//
//      //-------------------------------------------------------------------------------------
//      // scale alpha lin TODO check this
//      {
//        for (CI p=cnode->CoData().GetAlpha().begin();p!=cnode->CoData().GetAlpha().end();++p)
//        {
//          for (CImap pp=cnode->MoData().GetMltl().begin();pp!=cnode->MoData().GetMltl().end();++pp)
//          {
//            (cnode->CoData().GetDerivMlts()[pp->first])[p->first] -= (p->second) * (pp->second) * dvalue;
//          }
//        }
//      }
//
//      //-------------------------------------------------------------------------------------
//      // scale gap derivative entries for mortar data
//      for (CImap p=cnode->CoData().GetDerivGltl().begin();p!=cnode->CoData().GetDerivGltl().end();++p)
//        cnode->CoData().GetDerivGlts()[p->first] += alphaNew*(p->second)*dvalue;
//
//      //-------------------------------------------------------------------------------------
//      // scale gap derivative with lin alpha entries for mortar data
//      for (CI p=cnode->CoData().GetAlpha().begin();p!=cnode->CoData().GetAlpha().end();++p)
//      {
//        cnode->CoData().GetDerivGlts()[p->first] -= cnode->CoData().Getgltl()* dvalue *(p->second);
//      }
//    }// end scale nts
//  }

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate coupling type segment-to-segment coupl          farah 02/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::EvaluateSTS(
    const Epetra_Map& selecolmap,
    const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  MORTAR::MortarInterface::EvaluateSTS(selecolmap,mparams_ptr);
  return;
//  // loop over all slave col elements
//  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
//  {
//    int gid1 = selecolmap_->GID(i);
//    DRT::Element* ele1 = idiscret_->gElement(gid1);
//    if (!ele1)
//      dserror("ERROR: Cannot find slave element with gid %", gid1);
//    MORTAR::MortarElement* selement = dynamic_cast<MORTAR::MortarElement*>(ele1);
//
//    // loop over all slave nodes of this element to check for active nodes
//    bool eval = true;
//    for(int k = 0; k<selement->NumNode(); ++k)
//    {
//      CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(selement->Nodes()[k]);
//      if(cnode->Active())
//        eval=false;
//    }
//    if(!eval)
//      continue;
//
//    // skip zero-sized nurbs elements (slave)
//    if (selement->ZeroSized())
//      continue;
//
//    // empty vector of master element pointers
//    std::vector<MORTAR::MortarElement*> melements;
//
//    // loop over the candidate master elements of sele_
//    // use slave element's candidate list SearchElements !!!
//    for (int j = 0; j < selement->MoData().NumSearchElements(); ++j)
//    {
//      int gid2 = selement->MoData().SearchElements()[j];
//      DRT::Element* ele2 = idiscret_->gElement(gid2);
//      if (!ele2)
//        dserror("ERROR: Cannot find master element with gid %", gid2);
//      MORTAR::MortarElement* melement = dynamic_cast<MORTAR::MortarElement*>(ele2);
//
//      // skip zero-sized nurbs elements (master)
//      if (melement->ZeroSized())
//        continue;
//
//      melements.push_back(melement);
//    }
//
//    // concrete coupling evaluation routine
//    MortarCoupling(selement,melements,mparams_ptr);
//  }

  return;
}

/*----------------------------------------------------------------------*
 |  protected evaluate routine                               farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::EvaluateCoupling(
    const Epetra_Map& selecolmap,
    const Epetra_Map* snoderowmap,
    const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  // ask if nonsmooth contact is activated!
  if(nonSmoothContact_)
  {
    // 2D: only STS and nts has to be performed
    if(Dim()==2)
    {
      //********************************************************************
      // 1) perform coupling (projection + overlap detection for sl/m pairs)
      // 2) integrate Mortar matrix M and weighted gap g
      // 3) compute directional derivative of M and g and store into nodes
      //    (only for contact setting)
      //********************************************************************
      EvaluateSTS(selecolmap,mparams_ptr);

      //********************************************************************
      // 1) perform coupling (find closest point between to lines)
      // 2) evaluate gap and shape functions at this point
      // 3) compute directional derivative of entries and store into nodes
      //    (only for contact setting)
      //********************************************************************
      EvaluateNTS();

      //********************************************************************
      // NTN is a special case of NTS and an additional implementation is
      // not required!
      //********************************************************************
    }
    else if(Dim() == 3)
    {
      //********************************************************************
      // TODO: remove this hack!
      // HACK: LTL is not yet included in nonsmooth contact framework!
      // However, we want to test the LTL code separately. Thus, the "if"-
      // statement is included:
      // decide which type of coupling should be evaluated
      //********************************************************************
      INPAR::MORTAR::AlgorithmType algo =
          DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(imortar_, "ALGORITHM");
      if(algo == INPAR::MORTAR::algorithm_ltl)
      {
        EvaluateLTL();
        return;
      }

      //********************************************************************
      // 1) try to project slave nodes onto master elements
      // 2) evaluate shape functions at projected positions
      // 3) compute directional derivative of M and g and store into nodes
      //********************************************************************
      EvaluateNTS();

      //********************************************************************
      // 1) perform coupling (find closest point between to lines)
      // 2) evaluate gap and shape functions at this point
      // 3) compute directional derivative of entries and store into nodes
      //********************************************************************
      EvaluateLTL();

      //********************************************************************
      // 1) perform coupling (projection + line clipping edge surface pairs)
      // 2) integrate Mortar matrices D + M and weighted gap g
      // 3) compute directional derivative of D + M and g and store into nodes
      //********************************************************************
      EvaluateLTS();

      //********************************************************************
      // 1) perform coupling (projection + overlap detection for sl/m pairs)
      // 2) integrate Mortar matrix M and weighted gap g
      // 3) compute directional derivative of M and g and store into nodes
      //********************************************************************
      EvaluateSTS(selecolmap,mparams_ptr);

      //********************************************************************
      // perform LTS steps for master edges
      //********************************************************************
//      EvaluateLTSMaster();

      //********************************************************************
      // perform NTS steps for master edges
      //********************************************************************
//      EvaluateNTSMaster();

      //********************************************************************
      // NTN is a special case of NTS and an additional implementation is
      // not required!
      //********************************************************************
    }
    else
    {
      dserror("ERROR: Wrong dimension!");
    }
  }
  else
  {
    //********************************************************************
    // call base routine for standard mortar/nts evaluation
    //********************************************************************
    MORTAR::MortarInterface::EvaluateCoupling(selecolmap,snoderowmap,mparams_ptr);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Check and initialize corner/edge contact                 farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::InitializeCornerEdge()
{
  //return if nonsmooth contact is activated
  if(nonSmoothContact_)
    return;

  // call base function
  MORTAR::MortarInterface::InitializeCornerEdge();

  return;
}


/*----------------------------------------------------------------------*
 |  stuff for non-smooth contact geometries                 farah 02/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::DetectNonSmoothGeometries()
{
  dserror("ERROR: outdated!");

  std::vector<int> nonsmoothnodegids(0);
  std::vector<int> smoothnodegids(0);

  // loop over slave nodes to find nodes and eles
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find node with gid %", gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: Node ownership inconsistency!");

    if(cnode->NumElement()<2)
    {
      nonsmoothnodegids.push_back(cnode->Id());
      continue;
    }

    int loc = 0;
    double normalsele0[3] = {0.0, 0.0, 0.0};
    double normalsele1[3] = {0.0, 0.0, 0.0};

    // build normal at node for 1. ele
    CoElement* sele0 = dynamic_cast<CoElement*>(cnode->Elements()[0]);
    Epetra_SerialDenseMatrix elen0(6,1);
    sele0->BuildNormalAtNode(cnode->Id(),loc,elen0);

    // create 1. unit normal
    normalsele0[0] = elen0(0,0) / elen0(4,0);
    normalsele0[1] = elen0(1,0) / elen0(4,0);
    normalsele0[2] = elen0(2,0) / elen0(4,0);

    // build normal at node for 2. ele
    loc = 0;
    CoElement* sele1 = dynamic_cast<CoElement*>(cnode->Elements()[1]);
    Epetra_SerialDenseMatrix elen1(6,1);
    sele1->BuildNormalAtNode(cnode->Id(),loc,elen1);

    // create 2. unit normal
    normalsele1[0] = elen1(0,0) / elen1(4,0);
    normalsele1[1] = elen1(1,0) / elen1(4,0);
    normalsele1[2] = elen1(2,0) / elen1(4,0);

    double dot = normalsele0[0] * normalsele1[0] + normalsele0[1] * normalsele1[1] + normalsele0[2] * normalsele1[2];
    double R = abs(dot); // this is the abs ( cos(phi) ) --> 1.0 if parallel; 0.0 if orthogonal

    // critical node
    if(R<0.8) // circa 90 - 11.5 = 78.5
      nonsmoothnodegids.push_back(cnode->Id());
    else
      smoothnodegids.push_back(cnode->Id());
  }

  // create maps
  nonsmoothnodes_ = LINALG::CreateMap(nonsmoothnodegids,Comm());
  smoothnodes_    = LINALG::CreateMap(smoothnodegids,   Comm());

  // debug output
  std::cout << "# nonsmooth nodes: " << nonsmoothnodes_->NumGlobalElements() << std::endl;
  std::cout << "#    smooth nodes: " << smoothnodes_->NumGlobalElements()    << std::endl;

  // bye bye
  return;
}

/*----------------------------------------------------------------------*
 |  cpp to edge + Lin                                       farah 11/16 |
 *----------------------------------------------------------------------*/
double CONTACT::CoInterface::ComputeNormalNodeToEdge(
    MORTAR::MortarNode& snode,
    MORTAR::MortarElement& mele,
    double* normal,
    std::vector<GEN::pairedvector<int,double> >& normaltonodelin)
{
  // define tolerance
  const double tol = 1e-8;
  double dist = 1e12;
  int nrow = mele.NumNode();

  CoNode* node1 = dynamic_cast<CoNode*>(mele.Nodes()[0]);
  CoNode* node2 = dynamic_cast<CoNode*>(mele.Nodes()[1]);

  double length1=sqrt(node1->MoData().EdgeTangent()[0] * node1->MoData().EdgeTangent()[0] + node1->MoData().EdgeTangent()[1] * node1->MoData().EdgeTangent()[1] + node1->MoData().EdgeTangent()[2] * node1->MoData().EdgeTangent()[2]);
  double length2=sqrt(node2->MoData().EdgeTangent()[0] * node2->MoData().EdgeTangent()[0] + node2->MoData().EdgeTangent()[1] * node2->MoData().EdgeTangent()[1] + node2->MoData().EdgeTangent()[2] * node2->MoData().EdgeTangent()[2]);

  if(length1<1e-12 or length2<1e-12)
    return dist;

  // calc angle between tangents
  double t1[3] = {0.0,0.0,0.0};
  double t2[3] = {0.0,0.0,0.0};

  t1[0] = node1->MoData().EdgeTangent()[0];
  t1[1] = node1->MoData().EdgeTangent()[1];
  t1[2] = node1->MoData().EdgeTangent()[2];

  t2[0] = node2->MoData().EdgeTangent()[0];
  t2[1] = node2->MoData().EdgeTangent()[1];
  t2[2] = node2->MoData().EdgeTangent()[2];

  double test = t1[0]*t2[0] + t1[1]*t2[1] + t1[2]*t2[2];
  if(test<tol)
    dserror("ERROR: tangents have wrong direction!");


  double f  = 0.0;
  double df = 0.0;
  double xi[1] = {0.0};

  // newton loop
  for (int k = 0; k < MORTARMAXITER; ++k)
  {
//    std::cout << "k= " << k << std::endl;
    //**********************************************
    //  F CALCULATION                             //
    //**********************************************
    LINALG::SerialDenseVector sval(nrow);
    LINALG::SerialDenseMatrix sderiv(nrow,1);
    mele.EvaluateShape(xi,sval,sderiv,nrow);

    // tangent part
    double tangent[3] = {0.0,0.0,0.0};
    tangent[0] += sval[0] * dynamic_cast<CoNode*>(mele.Nodes()[0])->MoData().EdgeTangent()[0];
    tangent[1] += sval[0] * dynamic_cast<CoNode*>(mele.Nodes()[0])->MoData().EdgeTangent()[1];
    tangent[2] += sval[0] * dynamic_cast<CoNode*>(mele.Nodes()[0])->MoData().EdgeTangent()[2];

    tangent[0] += sval[1] * dynamic_cast<CoNode*>(mele.Nodes()[1])->MoData().EdgeTangent()[0];
    tangent[1] += sval[1] * dynamic_cast<CoNode*>(mele.Nodes()[1])->MoData().EdgeTangent()[1];
    tangent[2] += sval[1] * dynamic_cast<CoNode*>(mele.Nodes()[1])->MoData().EdgeTangent()[2];

    double tangentSlave = 0.0;
    tangentSlave =
        tangent[0] * snode.xspatial()[0] +
        tangent[1] * snode.xspatial()[1] +
        tangent[2] * snode.xspatial()[2];

    // master part
    double master[3] = {0.0,0.0,0.0};
    master[0] += sval[0] * node1->xspatial()[0];
    master[1] += sval[0] * node1->xspatial()[1];
    master[2] += sval[0] * node1->xspatial()[2];

    master[0] += sval[1] * node2->xspatial()[0];
    master[1] += sval[1] * node2->xspatial()[1];
    master[2] += sval[1] * node2->xspatial()[2];

    double tangentMaster = 0.0;
    tangentMaster =
        tangent[0] * master[0] +
        tangent[1] * master[1] +
        tangent[2] * master[2];

    f=tangentSlave - tangentMaster;
    if (abs(f) <= MORTARCONVTOL)
      break;
    //**********************************************
    //   F GRADIENT CALCULATION                   //
    //**********************************************
    // lin tangent part
    double lintangent[3] = {0.0,0.0,0.0};
    lintangent[0] += sderiv(0,0) * node1->MoData().EdgeTangent()[0];
    lintangent[1] += sderiv(0,0) * node1->MoData().EdgeTangent()[1];
    lintangent[2] += sderiv(0,0) * node1->MoData().EdgeTangent()[2];

    lintangent[0] += sderiv(1,0) * node2->MoData().EdgeTangent()[0];
    lintangent[1] += sderiv(1,0) * node2->MoData().EdgeTangent()[1];
    lintangent[2] += sderiv(1,0) * node2->MoData().EdgeTangent()[2];

    double lintangentSlave = 0.0;
    lintangentSlave =
        lintangent[0] * snode.xspatial()[0] +
        lintangent[1] * snode.xspatial()[1] +
        lintangent[2] * snode.xspatial()[2];

    // lin master part
    double linmaster[3] = {0.0,0.0,0.0};
    linmaster[0] += sderiv(0,0) * node1->xspatial()[0];
    linmaster[1] += sderiv(0,0) * node1->xspatial()[1];
    linmaster[2] += sderiv(0,0) * node1->xspatial()[2];

    linmaster[0] += sderiv(1,0) * node2->xspatial()[0];
    linmaster[1] += sderiv(1,0) * node2->xspatial()[1];
    linmaster[2] += sderiv(1,0) * node2->xspatial()[2];

    double lintangentMaster = 0.0;
    lintangentMaster =
        lintangent[0] * master[0] +
        lintangent[1] * master[1] +
        lintangent[2] * master[2];

    double tangentlinMaster = 0.0;
    tangentlinMaster =
        tangent[0] * linmaster[0] +
        tangent[1] * linmaster[1] +
        tangent[2] * linmaster[2];

    df = lintangentSlave - lintangentMaster - tangentlinMaster;
    if(abs(df)<1e-12)
      dserror("ERROR: df zero");
    xi[0] += -f/df;
  }

  //**********************************************
  //   CHECK XI                                 //
  //**********************************************
  if(-1.0-tol>xi[0] or xi[0]>1.0+tol)
    return dist;

  //**********************************************
  //   LINEARIZATION   df                       //
  //**********************************************

  LINALG::SerialDenseVector sval(nrow);
  LINALG::SerialDenseMatrix sderiv(nrow,1);
  mele.EvaluateShape(xi,sval,sderiv,nrow);

  // tangent part
  double tangent[3] = {0.0,0.0,0.0};
  tangent[0] += sval[0] * dynamic_cast<CoNode*>(mele.Nodes()[0])->MoData().EdgeTangent()[0];
  tangent[1] += sval[0] * dynamic_cast<CoNode*>(mele.Nodes()[0])->MoData().EdgeTangent()[1];
  tangent[2] += sval[0] * dynamic_cast<CoNode*>(mele.Nodes()[0])->MoData().EdgeTangent()[2];

  tangent[0] += sval[1] * dynamic_cast<CoNode*>(mele.Nodes()[1])->MoData().EdgeTangent()[0];
  tangent[1] += sval[1] * dynamic_cast<CoNode*>(mele.Nodes()[1])->MoData().EdgeTangent()[1];
  tangent[2] += sval[1] * dynamic_cast<CoNode*>(mele.Nodes()[1])->MoData().EdgeTangent()[2];

  // master part
  double master[3] = {0.0,0.0,0.0};
  master[0] += sval[0] * node1->xspatial()[0];
  master[1] += sval[0] * node1->xspatial()[1];
  master[2] += sval[0] * node1->xspatial()[2];

  master[0] += sval[1] * node2->xspatial()[0];
  master[1] += sval[1] * node2->xspatial()[1];
  master[2] += sval[1] * node2->xspatial()[2];

  // lin tangent part
  double lintangent[3] = {0.0,0.0,0.0};
  lintangent[0] += sderiv(0,0) * dynamic_cast<CoNode*>(mele.Nodes()[0])->MoData().EdgeTangent()[0];
  lintangent[1] += sderiv(0,0) * dynamic_cast<CoNode*>(mele.Nodes()[0])->MoData().EdgeTangent()[1];
  lintangent[2] += sderiv(0,0) * dynamic_cast<CoNode*>(mele.Nodes()[0])->MoData().EdgeTangent()[2];

  lintangent[0] += sderiv(1,0) * dynamic_cast<CoNode*>(mele.Nodes()[1])->MoData().EdgeTangent()[0];
  lintangent[1] += sderiv(1,0) * dynamic_cast<CoNode*>(mele.Nodes()[1])->MoData().EdgeTangent()[1];
  lintangent[2] += sderiv(1,0) * dynamic_cast<CoNode*>(mele.Nodes()[1])->MoData().EdgeTangent()[2];

  // lin master part
  double linmaster[3] = {0.0,0.0,0.0};
  linmaster[0] += sderiv(0,0) * node1->xspatial()[0];
  linmaster[1] += sderiv(0,0) * node1->xspatial()[1];
  linmaster[2] += sderiv(0,0) * node1->xspatial()[2];

  linmaster[0] += sderiv(1,0) * node2->xspatial()[0];
  linmaster[1] += sderiv(1,0) * node2->xspatial()[1];
  linmaster[2] += sderiv(1,0) * node2->xspatial()[2];


  //**********************************************
  //   LINEARIZATION    f                       //
  //**********************************************
  typedef GEN::pairedvector<int,double>::const_iterator _CI;

  std::vector<GEN::pairedvector<int,double> > linT(3,100);// added all sizes

  for (_CI p=node1->CoData().GetDerivTangent()[0].begin();p!=node1->CoData().GetDerivTangent()[0].end();++p)
    linT[0][p->first] += sval[0] * p->second;
  for (_CI p=node1->CoData().GetDerivTangent()[1].begin();p!=node1->CoData().GetDerivTangent()[1].end();++p)
    linT[1][p->first] += sval[0] * p->second;
  for (_CI p=node1->CoData().GetDerivTangent()[2].begin();p!=node1->CoData().GetDerivTangent()[2].end();++p)
    linT[2][p->first] += sval[0] * p->second;

  for (_CI p=node2->CoData().GetDerivTangent()[0].begin();p!=node2->CoData().GetDerivTangent()[0].end();++p)
    linT[0][p->first] += sval[1] * p->second;
  for (_CI p=node2->CoData().GetDerivTangent()[1].begin();p!=node2->CoData().GetDerivTangent()[1].end();++p)
    linT[1][p->first] += sval[1] * p->second;
  for (_CI p=node2->CoData().GetDerivTangent()[2].begin();p!=node2->CoData().GetDerivTangent()[2].end();++p)
    linT[2][p->first] += sval[1] * p->second;

  std::vector<GEN::pairedvector<int,double> > linXsl(3,100);// added all sizes
  linXsl[0][snode.Dofs()[0]] += 1.0;
  linXsl[1][snode.Dofs()[1]] += 1.0;
  linXsl[2][snode.Dofs()[2]] += 1.0;

  std::vector<GEN::pairedvector<int,double> > linXm(3,100);// added all sizes
  linXm[0][node1->Dofs()[0]] += sval[0];
  linXm[1][node1->Dofs()[1]] += sval[0];
  linXm[2][node1->Dofs()[2]] += sval[0];

  linXm[0][node2->Dofs()[0]] += sval[1];
  linXm[1][node2->Dofs()[1]] += sval[1];
  linXm[2][node2->Dofs()[2]] += sval[1];

  GEN::pairedvector<int,double> linf(100);// added all sizes
  for (_CI p=linT[0].begin();p!=linT[0].end();++p)
    linf[p->first] += snode.xspatial()[0] * p->second;
  for (_CI p=linT[1].begin();p!=linT[1].end();++p)
    linf[p->first] += snode.xspatial()[1] * p->second;
  for (_CI p=linT[2].begin();p!=linT[2].end();++p)
    linf[p->first] += snode.xspatial()[2] * p->second;

  for (_CI p=linXsl[0].begin();p!=linXsl[0].end();++p)
    linf[p->first] += tangent[0] * p->second;
  for (_CI p=linXsl[1].begin();p!=linXsl[1].end();++p)
    linf[p->first] += tangent[1] * p->second;
  for (_CI p=linXsl[2].begin();p!=linXsl[2].end();++p)
    linf[p->first] += tangent[2] * p->second;

  for (_CI p=linT[0].begin();p!=linT[0].end();++p)
    linf[p->first] -= master[0] * p->second;
  for (_CI p=linT[1].begin();p!=linT[1].end();++p)
    linf[p->first] -= master[1] * p->second;
  for (_CI p=linT[2].begin();p!=linT[2].end();++p)
    linf[p->first] -= master[2] * p->second;

  for (_CI p=linXm[0].begin();p!=linXm[0].end();++p)
    linf[p->first] -= tangent[0] * p->second;
  for (_CI p=linXm[1].begin();p!=linXm[1].end();++p)
    linf[p->first] -= tangent[1] * p->second;
  for (_CI p=linXm[2].begin();p!=linXm[2].end();++p)
    linf[p->first] -= tangent[2] * p->second;

  GEN::pairedvector<int,double> linXi(100);// added all sizes
  for (_CI p=linf.begin();p!=linf.end();++p)
    linXi[p->first] -= p->second / df;

  //**********************************************
  //   CALC NORMAL                              //
  //**********************************************
  double auxnormal[3] = {0.0, 0.0, 0.0};
  auxnormal[0] = snode.xspatial()[0] - master[0];
  auxnormal[1] = snode.xspatial()[1] - master[1];
  auxnormal[2] = snode.xspatial()[2] - master[2];

  // remove numerical artifacts
//  if(abs(auxnormal[0])<1e-12)
//    auxnormal[0] = 0.0;
//  if(abs(auxnormal[1])<1e-12)
//    auxnormal[1] = 0.0;
//  if(abs(auxnormal[2])<1e-12)
//    auxnormal[2] = 0.0;

  // calc distance
  dist =
      sqrt(auxnormal[0]*auxnormal[0] + auxnormal[1]*auxnormal[1] + auxnormal[2]* auxnormal[2]);

  if(abs(dist)<1e-12)
    return 1e12;

  //*******************************************
  // Lin:
  std::vector<GEN::pairedvector<int,double> > auxlin(3,100);// added all sizes

  // xslave
   for(int k = 0; k<3; ++k)
     (auxlin[k])[snode.Dofs()[k]] += 1.0;

   // xmaster n1
   for(int k = 0; k<3; ++k)
     (auxlin[k])[node1->Dofs()[k]] -= sval[0];
   // xmaster n2
   for(int k = 0; k<3; ++k)
     (auxlin[k])[node2->Dofs()[k]] -= sval[1];

   for (_CI p=linXi.begin();p!=linXi.end();++p)
   {
     for(int k = 0; k<3; ++k)
     {
       (auxlin[k])[p->first] -= sderiv(0,0) * node1->xspatial()[k] * p->second;
       (auxlin[k])[p->first] -= sderiv(1,0) * node2->xspatial()[k] * p->second;
     }
   }

   normal[0] = auxnormal[0];
   normal[1] = auxnormal[1];
   normal[2] = auxnormal[2];

   //******************************
   // Orientation check:
   double slavebasednormal[3] = {0.0, 0.0, 0.0};
   int nseg = snode.NumElement();
   DRT::Element** adjeles = snode.Elements();

   // we need to store some stuff here
   //**********************************************************************
   // elens(0,i): x-coord of element normal
   // elens(1,i): y-coord of element normal
   // elens(2,i): z-coord of element normal
   // elens(3,i): id of adjacent element i
   // elens(4,i): length of element normal
   // elens(5,i): length/area of element itself
   //**********************************************************************
   Epetra_SerialDenseMatrix elens(6,nseg);
   MORTAR::MortarElement* adjmrtrele = dynamic_cast<MORTAR::MortarElement*> (adjeles[0]);

   // build element normal at current node
   // (we have to pass in the index i to be able to store the
   // normal and other information at the right place in elens)
   int i = 0;
   adjmrtrele->BuildNormalAtNode(snode.Id(),i,elens);

   // add (weighted) element normal to nodal normal n
   for (int j=0;j<3;++j)
     slavebasednormal[j]+=elens(j,0)/elens(4,0);

   // create unit normal vector
   const double length = sqrt(slavebasednormal[0]*slavebasednormal[0]+slavebasednormal[1]*slavebasednormal[1]+slavebasednormal[2]*slavebasednormal[2]);
   if (abs(length)<1e-12)
   {
     dserror("ERROR: Nodal normal length 0, node ID %i",snode.Id());
   }
   else
   {
     for (int j=0;j<3;++j)
       slavebasednormal[j]/=length;
   }

   const double dotprod =
        -(normal[0]/dist*slavebasednormal[0] + normal[1]/dist*slavebasednormal[1] + normal[2]/dist*slavebasednormal[2]);

   if(dotprod < - 1e-12)
   {
     // get the cpp normal
     normal[0] = -normal[0];
     normal[1] = -normal[1];
     normal[2] = -normal[2];

     for(int j = 0; j<3 ;++j)
       for (_CI p=auxlin[j].begin();p!=auxlin[j].end();++p)
         (normaltonodelin[j])[p->first] -= (p->second);
   }
   else
   {
     // linearization
     for(int j = 0; j<3 ;++j)
       for (_CI p=auxlin[j].begin();p!=auxlin[j].end();++p)
         (normaltonodelin[j])[p->first] += (p->second);
   }

  return dist;
}

/*----------------------------------------------------------------------*
 |  cpp to node + Lin                                       farah 01/16 |
 *----------------------------------------------------------------------*/
double CONTACT::CoInterface::ComputeNormalNodeToNode(
    MORTAR::MortarNode& snode,
    MORTAR::MortarNode& mnode,
    double* normal,
    std::vector<GEN::pairedvector<int,double> >& normaltonodelin)
{
  const int dim = Dim();

  // distance between node and surface
  double gdist = 1e12;
  double gnormal[3] = {0.0, 0.0, 0.0};
  std::vector<GEN::pairedvector<int,double> > glin(3,1000);
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  double dist = 1e12;
  double auxnormal[3] = {0.0, 0.0, 0.0};

  // loop over found master nodes
  std::vector<GEN::pairedvector<int,double> > auxlin(3,1000);

  // calc vector
  auxnormal[0] = snode.xspatial()[0] - mnode.xspatial()[0];
  auxnormal[1] = snode.xspatial()[1] - mnode.xspatial()[1];
  auxnormal[2] = snode.xspatial()[2] - mnode.xspatial()[2];

  // remove numerical artifacts
  if(abs(auxnormal[0])<1e-12)
    auxnormal[0] = 0.0;
  if(abs(auxnormal[1])<1e-12)
    auxnormal[1] = 0.0;
  if(abs(auxnormal[2])<1e-12)
    auxnormal[2] = 0.0;

  // calc distance
  dist =
      sqrt(auxnormal[0]*auxnormal[0] + auxnormal[1]*auxnormal[1] + auxnormal[2]* auxnormal[2]);

  // if nodes lying on each other: continue to next master node
  if(abs(dist)<1e-12)
    return dist;

  //*******************************************
  // Lin:
  // xslave
   for(int k = 0; k<dim; ++k)
     (auxlin[k])[snode.Dofs()[k]] += 1.0;

   //xmaster
   for(int k = 0; k<dim; ++k)
     (auxlin[k])[mnode.Dofs()[k]] -= 1.0;

  // get normal
  gdist = dist;

  // normalize vector
  gnormal[0] = auxnormal[0];///dist;
  gnormal[1] = auxnormal[1];///dist;
  gnormal[2] = auxnormal[2];///dist;

  // linearization
  glin = auxlin;

  // get the cpp normal
  normal[0] = gnormal[0];
  normal[1] = gnormal[1];
  normal[2] = gnormal[2];

  //******************************
  // Orientation check:
  double slavebasednormal[3] = {0.0, 0.0, 0.0};
  int nseg = snode.NumElement();
  DRT::Element** adjeles = snode.Elements();

  // we need to store some stuff here
  //**********************************************************************
  // elens(0,i): x-coord of element normal
  // elens(1,i): y-coord of element normal
  // elens(2,i): z-coord of element normal
  // elens(3,i): id of adjacent element i
  // elens(4,i): length of element normal
  // elens(5,i): length/area of element itself
  //**********************************************************************
  Epetra_SerialDenseMatrix elens(6,nseg);
  MORTAR::MortarElement* adjmrtrele = dynamic_cast<MORTAR::MortarElement*> (adjeles[0]);

  // build element normal at current node
  // (we have to pass in the index i to be able to store the
  // normal and other information at the right place in elens)
  int i = 0;
  adjmrtrele->BuildNormalAtNode(snode.Id(),i,elens);

  // add (weighted) element normal to nodal normal n
  for (int j=0;j<3;++j)
    slavebasednormal[j]+=elens(j,0)/elens(4,0);

  // create unit normal vector
  const double length = sqrt(slavebasednormal[0]*slavebasednormal[0]+slavebasednormal[1]*slavebasednormal[1]+slavebasednormal[2]*slavebasednormal[2]);
  if (abs(length)<1e-12)
  {
    dserror("ERROR: Nodal normal length 0, node ID %i",snode.Id());
  }
  else
  {
    for (int j=0;j<3;++j)
      slavebasednormal[j]/=length;
  }

  const double dotprod =
      - (normal[0]*slavebasednormal[0] + normal[1]*slavebasednormal[1] + normal[2]*slavebasednormal[2]);

  if(dotprod < - 1e-12)
  {
    // get the cpp normal
    normal[0] = -normal[0];
    normal[1] = -normal[1];
    normal[2] = -normal[2];

    for(int j = 0; j<dim ;++j)
      for (CI p=glin[j].begin();p!=glin[j].end();++p)
        (normaltonodelin[j])[p->first] -= (p->second);
  }
  else
  {
    // linearization
    for(int j = 0; j<dim ;++j)
      for (CI p=glin[j].begin();p!=glin[j].end();++p)
        (normaltonodelin[j])[p->first] += (p->second);
  }


  return gdist;
}

/*----------------------------------------------------------------------*
 |  evaluate closest point normals                          farah 08/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::EvaluateCPPNormals()
{
  // Build averaged normal field on physically smooth surface
  // loop over proc's master nodes of the interface
  // use row map and export to column map later
  for (int i = 0; i < MasterRowNodes()->NumMyElements(); ++i)
  {
    int gid = MasterRowNodes()->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find node with gid %", gid);
    CoNode* mrtrnode = dynamic_cast<CoNode*>(node);

    // build averaged normal at each master node
    mrtrnode->BuildAveragedNormal();

    // build tangent
    if(mrtrnode->IsOnEdge())
      mrtrnode->BuildAveragedEdgeTangent();
  }

  // export nodal normals
  ExportMasterNodalNormals();

  // loop over slave nodes
  for (int i = 0; i < SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find node with gid %", gid);
    MORTAR::MortarNode* mrtrnode = dynamic_cast<MORTAR::MortarNode*>(node);

    if (mrtrnode->Owner() != Comm().MyPID())
      dserror("ERROR: Node ownership inconsistency!");

    // vector with possible contacting master eles/nodes
    std::vector<MORTAR::MortarElement*> meles;
    std::vector<MORTAR::MortarNode*>    mnodes;

    // fill vector with possibly contacting meles
    FindMEles(*mrtrnode,meles);

    // fallback solution if no mele is available
    if(meles.size()<1)// or !mrtrnode->IsOnCornerEdge())
    {
      CoNode* cnode = dynamic_cast<CoNode*>(mrtrnode);
      cnode->BuildAveragedNormal();
      continue;
    }


    // Here we have all found master elements for one slave node.
    // distance for cpp
    double normaltoline[3] = {0.0, 0.0, 0.0};
    std::vector<GEN::pairedvector<int,double> > normaltolineLin(3,1); // 1 dummy

    // Now, calculate distance between node and master line
    double dist = ComputeCPPNormal(
        *mrtrnode,
        meles,
        normaltoline,
        normaltolineLin);

    // if no projection was posible
    if(dist>1e11)
    {
      CoNode* cnode = dynamic_cast<CoNode*>(mrtrnode);
      cnode->BuildAveragedNormal();
      continue;
    }

    // set the normal and its lineratization
    SetCPPNormal(
        *mrtrnode,
        normaltoline,
        normaltolineLin);
  }

  // export slave normals
  ExportNodalNormals();

  // bye bye
  return;
}


/*----------------------------------------------------------------------*
 |  export master nodal normals (protected)                  farah 08/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::ExportMasterNodalNormals()
{
  std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > triad;

  std::map<int,std::vector<int> > n_x_key;
  std::map<int,std::vector<int> > n_y_key;
  std::map<int,std::vector<int> > n_z_key;
  std::map<int,std::vector<int> > txi_x_key;
  std::map<int,std::vector<int> > txi_y_key;
  std::map<int,std::vector<int> > txi_z_key;
  std::map<int,std::vector<int> > teta_x_key;
  std::map<int,std::vector<int> > teta_y_key;
  std::map<int,std::vector<int> > teta_z_key;

  std::map<int,std::vector<double> > n_x_val;
  std::map<int,std::vector<double> > n_y_val;
  std::map<int,std::vector<double> > n_z_val;
  std::map<int,std::vector<double> > txi_x_val;
  std::map<int,std::vector<double> > txi_y_val;
  std::map<int,std::vector<double> > txi_z_val;
  std::map<int,std::vector<double> > teta_x_val;
  std::map<int,std::vector<double> > teta_y_val;
  std::map<int,std::vector<double> > teta_z_val;

  GEN::pairedvector<int,double>::iterator iter;

  const Teuchos::RCP<Epetra_Map> masternodes = LINALG::AllreduceEMap(*(mnoderowmap_));

  // build info on row map
  for(int i=0; i<mnoderowmap_->NumMyElements();++i)
  {
    int gid = mnoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);

    // fill nodal matrix
    Teuchos::RCP<Epetra_SerialDenseMatrix> loc = Teuchos::rcp(new Epetra_SerialDenseMatrix(3,3));
    (*loc)(0,0) = cnode->MoData().n()[0];
    (*loc)(1,0) = cnode->MoData().n()[1];
    (*loc)(2,0) = cnode->MoData().n()[2];
    (*loc)(0,1) = cnode->CoData().txi()[0];
    (*loc)(1,1) = cnode->CoData().txi()[1];
    (*loc)(2,1) = cnode->CoData().txi()[2];
    (*loc)(0,2) = cnode->CoData().teta()[0];
    (*loc)(1,2) = cnode->CoData().teta()[1];
    (*loc)(2,2) = cnode->CoData().teta()[2];

    triad[gid] = loc;

    // fill nodal derivative vectors
    std::vector<GEN::pairedvector<int,double> >& derivn    = cnode->CoData().GetDerivN();
    std::vector<GEN::pairedvector<int,double> >& derivtxi  = cnode->CoData().GetDerivTxi();
    std::vector<GEN::pairedvector<int,double> >& derivteta = cnode->CoData().GetDerivTeta();

    for(iter=derivn[0].begin();iter!=derivn[0].end();++iter)
    {
      n_x_key[gid].push_back(iter->first);
      n_x_val[gid].push_back(iter->second);
    }
    for(iter=derivn[1].begin();iter!=derivn[1].end();++iter)
    {
      n_y_key[gid].push_back(iter->first);
      n_y_val[gid].push_back(iter->second);
    }
    for(iter=derivn[2].begin();iter!=derivn[2].end();++iter)
    {
      n_z_key[gid].push_back(iter->first);
      n_z_val[gid].push_back(iter->second);
    }

    for(iter=derivtxi[0].begin();iter!=derivtxi[0].end();++iter)
    {
      txi_x_key[gid].push_back(iter->first);
      txi_x_val[gid].push_back(iter->second);
    }
    for(iter=derivtxi[1].begin();iter!=derivtxi[1].end();++iter)
    {
      txi_y_key[gid].push_back(iter->first);
      txi_y_val[gid].push_back(iter->second);
    }
    for(iter=derivtxi[2].begin();iter!=derivtxi[2].end();++iter)
    {
      txi_z_key[gid].push_back(iter->first);
      txi_z_val[gid].push_back(iter->second);
    }

    for(iter=derivteta[0].begin();iter!=derivteta[0].end();++iter)
    {
      teta_x_key[gid].push_back(iter->first);
      teta_x_val[gid].push_back(iter->second);
    }
    for(iter=derivteta[1].begin();iter!=derivteta[1].end();++iter)
    {
      teta_y_key[gid].push_back(iter->first);
      teta_y_val[gid].push_back(iter->second);
    }
    for(iter=derivteta[2].begin();iter!=derivteta[2].end();++iter)
    {
      teta_z_key[gid].push_back(iter->first);
      teta_z_val[gid].push_back(iter->second);
    }
  }

  // communicate from master node row to column map
  DRT::Exporter ex(*mnoderowmap_,*masternodes,Comm());
  ex.Export(triad);

  ex.Export(n_x_key);
  ex.Export(n_x_val);
  ex.Export(n_y_key);
  ex.Export(n_y_val);
  ex.Export(n_z_key);
  ex.Export(n_z_val);

  ex.Export(txi_x_key);
  ex.Export(txi_x_val);
  ex.Export(txi_y_key);
  ex.Export(txi_y_val);
  ex.Export(txi_z_key);
  ex.Export(txi_z_val);

  ex.Export(teta_x_key);
  ex.Export(teta_x_val);
  ex.Export(teta_y_key);
  ex.Export(teta_y_val);
  ex.Export(teta_z_key);
  ex.Export(teta_z_val);

  // extract info on column map
  for(int i=0; i<masternodes->NumMyElements();++i)
  {
    // only do something for ghosted nodes
    int gid = masternodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
    int linsize = cnode->GetLinsize()+(int)(n_x_key[gid].size());

    if (cnode->Owner()==Comm().MyPID())
      continue;

    // extract info
    Teuchos::RCP<Epetra_SerialDenseMatrix> loc = triad[gid];
    cnode->MoData().n()[0]    = (*loc)(0,0);
    cnode->MoData().n()[1]    = (*loc)(1,0);
    cnode->MoData().n()[2]    = (*loc)(2,0);
    cnode->CoData().txi()[0]  = (*loc)(0,1);
    cnode->CoData().txi()[1]  = (*loc)(1,1);
    cnode->CoData().txi()[2]  = (*loc)(2,1);
    cnode->CoData().teta()[0] = (*loc)(0,2);
    cnode->CoData().teta()[1] = (*loc)(1,2);
    cnode->CoData().teta()[2] = (*loc)(2,2);

    // extract derivative info
    std::vector<GEN::pairedvector<int,double> >& derivn    = cnode->CoData().GetDerivN();
    std::vector<GEN::pairedvector<int,double> >& derivtxi  = cnode->CoData().GetDerivTxi();
    std::vector<GEN::pairedvector<int,double> >& derivteta = cnode->CoData().GetDerivTeta();

    for (int k=0;k<(int)(derivn.size());++k)
      derivn[k].clear();
    derivn.resize(3,linsize);
    for (int k=0;k<(int)(derivtxi.size());++k)
      derivtxi[k].clear();
    derivtxi.resize(3,linsize);
    for (int k=0;k<(int)(derivteta.size());++k)
      derivteta[k].clear();
    derivteta.resize(3,linsize);

    cnode->CoData().GetDerivN()[0].resize(linsize);
    cnode->CoData().GetDerivN()[1].resize(linsize);
    cnode->CoData().GetDerivN()[2].resize(linsize);

    cnode->CoData().GetDerivTxi()[0].resize(linsize);
    cnode->CoData().GetDerivTxi()[1].resize(linsize);
    cnode->CoData().GetDerivTxi()[2].resize(linsize);

    cnode->CoData().GetDerivTeta()[0].resize(linsize);
    cnode->CoData().GetDerivTeta()[1].resize(linsize);
    cnode->CoData().GetDerivTeta()[2].resize(linsize);

    for (int k=0;k<(int)(n_x_key[gid].size());++k)
      (cnode->CoData().GetDerivN()[0])[n_x_key[gid][k]] = n_x_val[gid][k];
    for (int k=0;k<(int)(n_y_key[gid].size());++k)
      (cnode->CoData().GetDerivN()[1])[n_y_key[gid][k]] = n_y_val[gid][k];
    for (int k=0;k<(int)(n_z_key[gid].size());++k)
      (cnode->CoData().GetDerivN()[2])[n_z_key[gid][k]] = n_z_val[gid][k];

    for (int k=0;k<(int)(txi_x_key[gid].size());++k)
      (cnode->CoData().GetDerivTxi()[0])[txi_x_key[gid][k]] = txi_x_val[gid][k];
    for (int k=0;k<(int)(txi_y_key[gid].size());++k)
      (cnode->CoData().GetDerivTxi()[1])[txi_y_key[gid][k]] = txi_y_val[gid][k];
    for (int k=0;k<(int)(txi_z_key[gid].size());++k)
      (cnode->CoData().GetDerivTxi()[2])[txi_z_key[gid][k]] = txi_z_val[gid][k];

    for (int k=0;k<(int)(teta_x_key[gid].size());++k)
      (cnode->CoData().GetDerivTeta()[0])[teta_x_key[gid][k]] = teta_x_val[gid][k];
    for (int k=0;k<(int)(teta_y_key[gid].size());++k)
      (cnode->CoData().GetDerivTeta()[1])[teta_y_key[gid][k]] = teta_y_val[gid][k];
    for (int k=0;k<(int)(teta_z_key[gid].size());++k)
      (cnode->CoData().GetDerivTeta()[2])[teta_z_key[gid][k]] = teta_z_val[gid][k];
  }

  // free memory
  triad.clear();

  n_x_key.clear();
  n_y_key.clear();
  n_z_key.clear();
  txi_x_key.clear();
  txi_y_key.clear();
  txi_z_key.clear();
  teta_x_key.clear();
  teta_y_key.clear();
  teta_z_key.clear();

  n_x_val.clear();
  n_y_val.clear();
  n_z_val.clear();
  txi_x_val.clear();
  txi_y_val.clear();
  txi_z_val.clear();
  teta_x_val.clear();
  teta_y_val.clear();
  teta_z_val.clear();

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate nodal normals (public)                          farah 02/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::EvaluateAveragedNodalNormals()
{
  dserror("ERROR: outdated function!");
  //safety
  if(smoothnodes_ == Teuchos::null)
    dserror("ERROR: map of non smooth nodes is wrong!");

  // loop over proc's slave nodes of the interface
  // use row map and export to column map later
  // (use boundary map to include slave side boundary nodes)
  for (int i = 0; i < smoothnodes_->NumMyElements(); ++i)
  {
    int gid = smoothnodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find node with gid %", gid);
    CoNode* mrtrnode = dynamic_cast<CoNode*>(node);

    // build averaged normal at each slave node
    mrtrnode->BuildAveragedNormal();
  }

  return;
}


/*----------------------------------------------------------------------*
 |  calc scaling for hybrid formulation                     farah 09/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::ComputeScaling()
{
  // ltl scaling only for 3D setting
  if(Dim() == 2)
    return;

  // compute ltl scaling for point and line contact
  ComputeScalingLTL();

  // bye
  return;
}


/*----------------------------------------------------------------------*
 |  calc scaling for hybrid formulation                     farah 01/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::ComputeScalingLTL()
{
  std::cout << "ComputeScalingLTL" << std::endl;

  // define iterator for linerization
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  // angle in degree
  const double minAngle = IParams().get<double>("HYBRID_ANGLE_MIN");
  const double maxAngle = IParams().get<double>("HYBRID_ANGLE_MAX");

  // check
  if(minAngle<0.0 or maxAngle < 0.0)
    dserror("ERROR: invalid hybrid angle!");
  if(minAngle>=maxAngle)
    dserror("ERROR: invalid hybrid angle!");

  // angle in rad
  const double alphaMin = minAngle*(M_PI/180.0); // min angle in rad
  const double alphaMax = maxAngle*(M_PI/180.0); // max angle in rad

  //======================================================
  // search slave edges and master edges and store them

  // guarantee uniquness of slave edges
  std::set<std::pair<int,int> > donebeforeS;

  // loop over slave elements
  for(int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1)
      dserror("ERROR: Cannot find slave element with gid %", gid1);
    CoElement* selement = dynamic_cast<CoElement*>(ele1);

    // empty vector of slave element pointers
    std::vector<Teuchos::RCP<MORTAR::MortarElement> > lineElementsS;

    if(selement->Shape() == DRT::Element::quad4)
    {
      for(int j = 0; j< 4 ; ++j)
      {
        int nodeIds[2]  = {0,0};
        int nodeLIds[2] = {0,0};

        if(j == 0)
        {
          nodeIds[0] = selement->NodeIds()[0];
          nodeIds[1] = selement->NodeIds()[1];

          nodeLIds[0] = 0;
          nodeLIds[1] = 1;
        }
        else if(j == 1)
        {
          nodeIds[0] = selement->NodeIds()[1];
          nodeIds[1] = selement->NodeIds()[2];

          nodeLIds[0] = 1;
          nodeLIds[1] = 2;
        }
        else if(j == 2)
        {
          nodeIds[0] = selement->NodeIds()[2];
          nodeIds[1] = selement->NodeIds()[3];

          nodeLIds[0] = 2;
          nodeLIds[1] = 3;
        }
        else if(j == 3)
        {
          nodeIds[0] = selement->NodeIds()[3];
          nodeIds[1] = selement->NodeIds()[0];

          nodeLIds[0] = 3;
          nodeLIds[1] = 0;
        }

        // check if both nodes on edge geometry
        bool node0Edge = dynamic_cast<MORTAR::MortarNode*>(selement->Nodes()[nodeLIds[0]])->IsOnEdge();
        bool node1Edge = dynamic_cast<MORTAR::MortarNode*>(selement->Nodes()[nodeLIds[1]])->IsOnEdge();

        if(!node0Edge or !node1Edge)
          continue;

        //create pair
        std::pair<int,int> actIDs   = std::pair<int,int>(nodeIds[0],nodeIds[1]);
        std::pair<int,int> actIDstw = std::pair<int,int>(nodeIds[1],nodeIds[0]);

        // check if processed before
        std::set<std::pair<int,int> >::iterator iter   = donebeforeS.find(actIDs);
        std::set<std::pair<int,int> >::iterator itertw = donebeforeS.find(actIDstw);

         // if not then create ele
         if (iter == donebeforeS.end() and itertw == donebeforeS.end())
         {
           // add to set of processed nodes
           donebeforeS.insert(actIDs);
           donebeforeS.insert(actIDstw);

           // create line ele:
           Teuchos::RCP<MORTAR::MortarElement> lineEle = Teuchos::rcp(
                       new MORTAR::MortarElement(
                           j,
                           selement->Owner(),
                           DRT::Element::line2,
                           2,
                           nodeIds,
                           false));

           // get nodes
           DRT::Node* nodes[2] = {selement->Nodes()[nodeLIds[0]],selement->Nodes()[nodeLIds[1]]};
           lineEle->BuildNodalPointers(nodes);

           // init data container for dual shapes
           lineEle->InitializeDataContainer();

           // push back into vector
           lineElementsS.push_back(lineEle);
         }
      } // end edge loop
    }
    else
      dserror("ERROR: LTL only for quad4!");

    // guarantee uniquness of master edges
    std::set<std::pair<int,int> > donebeforeM;

    // empty vector of slave element pointers
    std::vector<Teuchos::RCP<MORTAR::MortarElement> > lineElementsM;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int k = 0; k < selement->MoData().NumSearchElements(); ++k)
    {
      int gid2 = selement->MoData().SearchElements()[k];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2)
        dserror("ERROR: Cannot find master element with gid %", gid2);
      CoElement* melement = dynamic_cast<CoElement*>(ele2);

      if(melement->Shape() == DRT::Element::quad4)
      {
        for(int j = 0; j< 4 ; ++j)
        {
          int nodeIds[2]  = {0,0};
          int nodeLIds[2] = {0,0};

          if(j == 0)
          {
            nodeIds[0] = melement->NodeIds()[0];
            nodeIds[1] = melement->NodeIds()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if(j == 1)
          {
            nodeIds[0] = melement->NodeIds()[1];
            nodeIds[1] = melement->NodeIds()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if(j == 2)
          {
            nodeIds[0] = melement->NodeIds()[2];
            nodeIds[1] = melement->NodeIds()[3];

            nodeLIds[0] = 2;
            nodeLIds[1] = 3;
          }
          else if(j == 3)
          {
            nodeIds[0] = melement->NodeIds()[3];
            nodeIds[1] = melement->NodeIds()[0];

            nodeLIds[0] = 3;
            nodeLIds[1] = 0;
          }

          // check if both nodes on edge geometry
          bool node0Edge = dynamic_cast<MORTAR::MortarNode*>(melement->Nodes()[nodeLIds[0]])->IsOnEdge();
          bool node1Edge = dynamic_cast<MORTAR::MortarNode*>(melement->Nodes()[nodeLIds[1]])->IsOnEdge();

          if(!node0Edge or !node1Edge)
            continue;

          //create pair
          std::pair<int,int> actIDs   = std::pair<int,int>(nodeIds[0],nodeIds[1]);
          std::pair<int,int> actIDstw = std::pair<int,int>(nodeIds[1],nodeIds[0]);

          // check if processed before
          std::set<std::pair<int,int> >::iterator iter = donebeforeM.find(actIDs);
          std::set<std::pair<int,int> >::iterator itertw = donebeforeM.find(actIDstw);

           // if not then create ele
           if (iter == donebeforeM.end() and itertw == donebeforeM.end())
           {
             // add to set of processed nodes
             donebeforeM.insert(actIDs);
             donebeforeM.insert(actIDstw);

             // create line ele:
             Teuchos::RCP<MORTAR::MortarElement> lineEle = Teuchos::rcp(
                         new MORTAR::MortarElement(
                             j,
                             melement->Owner(),
                             DRT::Element::line2,
                             2,
                             nodeIds,
                             false));

             // get nodes
             DRT::Node* nodes[2] = {melement->Nodes()[nodeLIds[0]],melement->Nodes()[nodeLIds[1]]};
             lineEle->BuildNodalPointers(nodes);

             // init data container for dual shapes
             lineEle->InitializeDataContainer();

             // push back into vector
             lineElementsM.push_back(lineEle);
           }
        } // end edge loop
      }
      else
        dserror("ERROR: LTL only for quad4!");
    }// end found mele loop

    // loop over slave edges
    for( int s = 0; s<(int)lineElementsS.size(); ++s)
    {
      double gR = 1e12;
      bool parallel = false;

      // loop over master edges
      for( int m = 0; m<(int)lineElementsM.size(); ++m)
      {
        // 1. search slave master pair with cpp point
        LineToLineCouplingPoint3d coup(
            *idiscret_,
            3,
            IParams(),
            lineElementsS[s],
            lineElementsM[m]);

        parallel = coup.CheckParallelity();
        if(parallel)
        {
          continue;
//          dserror("ERROR: edges parallel");
          break;
        }
        // create empty points
        double sxi[1] = {0.0};
        double mxi[1] = {0.0};
        GEN::pairedvector<int,double> dsxi(3*lineElementsM[m]->NumNode() + 3*lineElementsS[s]->NumNode());
        GEN::pairedvector<int,double> dmxi(3*lineElementsM[m]->NumNode() + 3*lineElementsS[s]->NumNode());

        coup.LineIntersection(sxi, mxi, dsxi, dmxi);
        bool valid = coup.CheckIntersection(sxi,mxi);

        // no valid intersection
        if(!valid)
        {
          continue;
        }
        // found valid intersection
        else
        {
          std::cout << "VALID INTERSECTION" << std::endl;
          // 2. calc current angle for this element pair
          GEN::pairedvector<int,double> linAngle(100 + 3*lineElementsM[m]->NumNode() + 3*lineElementsS[s]->NumNode());
          gR = coup.CalcCurrentAngle(linAngle);

          // get vector of nodes
          std::vector<CoNode*> cnodes;
          for(int nn = 0; nn<lineElementsS[s]->NumNode(); ++nn)
            cnodes.push_back(dynamic_cast<CoNode*>(lineElementsS[s]->Nodes()[nn]));

          //======================================================
          // calc alpha:
          double alpha = 0.0;

          std::cout << "gR= " << gR << std::endl;

          // ltl line contact
          if(gR<=alphaMin)
          {
            std::cout << "LTS contact" << std::endl;
            alpha = 1.0;
          }
          // ltl hybrid
          else if(gR>alphaMin and gR<alphaMax)
          {
            std::cout << "current transition angle in degree = " << gR*(180.0/M_PI) << std::endl;
            alpha=0.5*(1.0 - cos(M_PI*(gR-alphaMax)/(alphaMin-alphaMax)));
            std::cout << "alpha = " << alpha << std::endl;
          }
          // ltl point contact
          else
          {
            std::cout << "LTL contact" << std::endl;
            alpha = 0.0;
          }

          // store to data container
          for(int nn = 0; nn<lineElementsS[s]->NumNode(); ++nn)
            cnodes[nn]->CoData().GetAlphaN() = alpha;

          //======================================================
          // lin alpha:
          if(gR<=alphaMin)
          {
            // nothing: pure ltl line contact
          }
          // hybrid linearization:
          else if (gR>alphaMin and gR<alphaMax)
          {
            // clear old data
            for(int nn = 0; nn<lineElementsS[s]->NumNode(); ++nn)
              if(cnodes[nn]->CoData().GetAlpha().size() == 0)
                cnodes[nn]->CoData().GetAlpha().resize(1000);

            const double fac = (M_PI/(2.0*(alphaMin-alphaMax))) * sin(M_PI*(gR-alphaMax)/(alphaMin-alphaMax));
            for (CI p=linAngle.begin();p!=linAngle.end();++p)
            {
              for(int nn = 0; nn<lineElementsS[s]->NumNode(); ++nn)
                (cnodes[nn]->CoData().GetAlpha())[p->first] = fac * (p->second);
            }
          }
          else
          {
            // nothing: pure ltl point contact
          }
          break;
        }// valid projection

      }
    }
  }// end slave loop








//  //======================================================
//  // loop over smooth slave nodes
//  for (int i = 0; i < smoothnodes_->NumMyElements(); ++i)
//  {
//    int gid = smoothnodes_->GID(i);
//    DRT::Node* node = idiscret_->gNode(gid);
//    if (!node)
//      dserror("ERROR: Cannot find node with gid %", gid);
//    CoNode* cnode = dynamic_cast<CoNode*>(node);
//
//    if (cnode->Owner() != Comm().MyPID())
//      dserror("ERROR: Node ownership inconsistency!");
//
//    // MORTAR
//    cnode->CoData().GetAlphaN() = -1.0;
//  }
//
//  //======================================================
//  // loop over non smooth slave nodes
//  for (int i = 0; i < nonsmoothnodes_->NumMyElements(); ++i)
//  {
//    int gid = nonsmoothnodes_->GID(i);
//    DRT::Node* node = idiscret_->gNode(gid);
//    if (!node)
//      dserror("ERROR: Cannot find node with gid %", gid);
//    CoNode* cnode = dynamic_cast<CoNode*>(node);
//
//    if (cnode->Owner() != Comm().MyPID())
//      dserror("ERROR: Node ownership inconsistency!");
//
//    // clear old data
//    if(cnode->CoData().GetAlpha().size() == 0)
//      cnode->CoData().GetAlpha().resize(1000);
//
//    // get cpp normal
//    double normalcpp[3] = {0.0, 0.0, 0.0};
//    normalcpp[0] = cnode->MoData().n()[0];
//    normalcpp[1] = cnode->MoData().n()[1];
//    normalcpp[2] = cnode->MoData().n()[2];
//
//    // NTS evtl...
//    // calc element normal
//    const int eles = cnode->NumElement();
//    if(eles>2)
//      dserror("ERROR: element number to high!");
//
//    // calc smallest angle between cpp and nele
//    double gR = 1e12;
//    GEN::pairedvector<int,double> ddotcppfinal(1000);
//
//    for (int n = 0; n<eles; ++n)
//    {
//      // create element normal:
//      int loc = 0;
//      CoElement* seleaux = dynamic_cast<CoElement*>(cnode->Elements()[n]);
//      Epetra_SerialDenseMatrix elens(6,1);
//      seleaux->BuildNormalAtNode(cnode->Id(),loc,elens);
//
//      double normalsele[3] = {0.0, 0.0, 0.0};
//      normalsele[0] = elens(0,0) / elens(4,0);
//      normalsele[1] = elens(1,0) / elens(4,0);
//      normalsele[2] = elens(2,0) / elens(4,0);
//
//      // calculate angle between cpp and elenormal
//      double dotcpp = normalsele[0] * normalcpp[0] + normalsele[1] * normalcpp[1] + normalsele[2] * normalcpp[2];
//      double Rl = acos(dotcpp); // current angle in rad
//
//      if(Rl<0.0)
//        dserror("ERROR: angle less than 0.0");
//
//      //======================================================
//      // lin angle
//      // lin Rl
//      double fac = (-1.0/(sqrt(1.0-dotcpp*dotcpp)));
//
//      // init lin
//      GEN::pairedvector<int,double> ddotcpp(1000);
//      std::vector<GEN::pairedvector<int,double> > dnele(3,1000);
//
//      // deriv unit normal!
//      seleaux->DerivNormalAtNode(cnode->Id(),loc,elens,dnele);
//
//      // linearization of cpp normal
//      for (CI p=cnode->CoData().GetDerivN()[0].begin();p!=cnode->CoData().GetDerivN()[0].end();++p)
//        (ddotcpp)[p->first] += (p->second) * normalsele[0];
//      for (CI p=cnode->CoData().GetDerivN()[1].begin();p!=cnode->CoData().GetDerivN()[1].end();++p)
//        (ddotcpp)[p->first] += (p->second) * normalsele[1];
//      for (CI p=cnode->CoData().GetDerivN()[2].begin();p!=cnode->CoData().GetDerivN()[2].end();++p)
//        (ddotcpp)[p->first] += (p->second) * normalsele[2];
//
//      // linearization of ele normal
//      for (CI p=dnele[0].begin();p!=dnele[0].end();++p)
//        (ddotcpp)[p->first] += (p->second) * normalcpp[0];
//      for (CI p=dnele[1].begin();p!=dnele[1].end();++p)
//        (ddotcpp)[p->first] += (p->second) * normalcpp[1];
//      for (CI p=dnele[2].begin();p!=dnele[2].end();++p)
//        (ddotcpp)[p->first] += (p->second) * normalcpp[2];
//
//      // multiply with fac: -asin(dotcpp)
//      for (CI p=ddotcpp.begin();p!=ddotcpp.end();++p)
//        (ddotcpp)[p->first] = fac * (p->second);
//
//
//      // safe smaller angle and linearization
//      if(Rl<gR)
//      {
//        gR=Rl;
//        ddotcppfinal = ddotcpp;
//      }
//    }// end ele loop
//
//
//    //======================================================
//    // calc alpha:
//    double alpha = 0.0;
//    if(gR<=alphaMin)
//    {
//      alpha = 1.0;
//    }
//    else if(gR>alphaMin and gR<alphaMax)
//    {
//      std::cout << "current transition angle in degree = " << gR*(180.0/3.141592) << std::endl;
//      alpha=0.5*(1.0 - cos(3.141592*(gR-alphaMax)/(alphaMin-alphaMax)));
//      std::cout << "alpha = " << alpha << std::endl;
//    }
//    else
//    {
//      alpha = 0.0;
//    }
//
//    // store to data container
//    cnode->CoData().GetAlphaN() = alpha;
//
//    //======================================================
//    // lin alpha:
//    if(gR<=alphaMin)
//    {
//      // nothing: pure mortar
//    }
//    // hybrid linearization:
//    else if (gR>alphaMin and gR<alphaMax)
//    {
//      const double fac = (3.141592/(2.0*(alphaMin-alphaMax))) * sin(3.141592*(gR-alphaMax)/(alphaMin-alphaMax));
//      for (CI p=ddotcppfinal.begin();p!=ddotcppfinal.end();++p)
//        (cnode->CoData().GetAlpha())[p->first] = fac * (p->second);
//    }
//    else
//    {
//      // nothing: pure nts
//    }
//  }// end node loop

  // bye bye
  return;
}

/*----------------------------------------------------------------------*
 |  scale normals for hybrid formulation                    farah 05/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::ScaleNormals()
{
  dserror("ERROR: outdated!");

  if(Dim()==2)
    ScaleNormals2D();
  else if(Dim() == 3)
    ScaleNormals3D();
  else
    dserror("ERROR: Wrong dimension!");
}

/*----------------------------------------------------------------------*
 |  scale normals for hybrid formulation                    farah 05/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::ScaleNormals2D()
{
  // define iterator for paired vector
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  //======================================================
  // loop over non smooth slave nodes
  for (int i = 0; i < nonsmoothnodes_->NumMyElements(); ++i)
  {
    int gid = nonsmoothnodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find node with gid %", gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: Node ownership inconsistency!");

    // get cpp normal
    double normalcpp[3]  = {0.0, 0.0, 0.0};

    // get cpp normal lin.
    std::vector<GEN::pairedvector<int,double> > dcppaux(3,1000);

    dcppaux = cnode->CoData().GetDerivN();
    normalcpp[0] = cnode->MoData().n()[0];
    normalcpp[1] = cnode->MoData().n()[1];
    normalcpp[2] = cnode->MoData().n()[2];

    if(cnode->Active())
    {
      std::cout << "cnode->CoData().GetAlphaN()] = " << cnode->CoData().GetAlphaN()  << std::endl;
      std::cout << "normalcpp[0] = " << normalcpp[0]  << std::endl;
      std::cout << "normalcpp[1] = " << normalcpp[1]  << std::endl;
      std::cout << "normalcpp[2] = " << normalcpp[2]  << std::endl;
    }


    // calc element normal
    const int eles = cnode->NumElement();
    if(eles>2)
      dserror("ERROR: element number to high!");

    double gR = 1e12;
    int localid=-1;

    for (int n = 0; n<eles; ++n)
    {
      // create element normal:
      int loc = 0;
      CoElement* seleaux = dynamic_cast<CoElement*>(cnode->Elements()[n]);
      Epetra_SerialDenseMatrix elens(6,1);
      seleaux->BuildNormalAtNode(cnode->Id(),loc,elens);

      double normalsele[3] = {0.0, 0.0, 0.0};
      normalsele[0] = elens(0,0) / elens(4,0);
      normalsele[1] = elens(1,0) / elens(4,0);
      normalsele[2] = elens(2,0) / elens(4,0);

      // calculate angle between cpp and elenormal
      double dotcpp = normalsele[0] * normalcpp[0] + normalsele[1] * normalcpp[1] + normalsele[2] * normalcpp[2];
      double Rl = acos(dotcpp); // current angle in rad

      if(cnode->Active())
      {
        std::cout << "normalsele[0] = " << normalsele[0]  << std::endl;
        std::cout << "normalsele[1] = " << normalsele[1]  << std::endl;
        std::cout << "normalsele[2] = " << normalsele[2]  << std::endl;

        std::cout << "dotcpp= " << dotcpp << std::endl;
        std::cout << "Rl= " << Rl << std::endl;
      }


      if(Rl<0.0)
        dserror("ERROR: angle less than 0.0");



      // safe smaller angle and linearization
      if(Rl<gR)
      {
        gR=Rl;
        localid = n;
      }
    }// end ele loop

    // create new normals:
    const double alpha = cnode->CoData().GetAlphaN();

    // create element normal:
    int loc = 0;
    CoElement* seleaux = dynamic_cast<CoElement*>(cnode->Elements()[localid]);
    Epetra_SerialDenseMatrix elens(6,1);
    seleaux->BuildNormalAtNode(cnode->Id(),loc,elens);

    double normalsele[3] = {0.0, 0.0, 0.0};
    normalsele[0] = elens(0,0) / elens(4,0);
    normalsele[1] = elens(1,0) / elens(4,0);
    normalsele[2] = elens(2,0) / elens(4,0);

    cnode->MoData().n()[0] = alpha * normalsele[0] + (1.0 - alpha) * normalcpp[0];
    cnode->MoData().n()[1] = alpha * normalsele[1] + (1.0 - alpha) * normalcpp[1];
    cnode->MoData().n()[2] = alpha * normalsele[2] + (1.0 - alpha) * normalcpp[2];


    // deriv unit normal!
    std::vector<GEN::pairedvector<int,double> > dnele(3,1000);
    seleaux->DerivNormalAtNode(cnode->Id(),loc,elens,dnele);

    cnode->CoData().GetDerivN()[0].clear();
    cnode->CoData().GetDerivN()[1].clear();
    cnode->CoData().GetDerivN()[2].clear();

    // lin ele normal * alpha
    for (CI p=dnele[0].begin();p!=dnele[0].end();++p)
      (cnode->CoData().GetDerivN()[0])[p->first] += alpha * (p->second);
    for (CI p=dnele[1].begin();p!=dnele[1].end();++p)
      (cnode->CoData().GetDerivN()[1])[p->first] += alpha * (p->second);
    for (CI p=dnele[2].begin();p!=dnele[2].end();++p)
      (cnode->CoData().GetDerivN()[2])[p->first] += alpha * (p->second);

    // lin cpp normal * alpha
    for (CI p=dcppaux[0].begin();p!=dcppaux[0].end();++p)
      (cnode->CoData().GetDerivN()[0])[p->first] += (1.0 -alpha) * (p->second);
    for (CI p=dcppaux[1].begin();p!=dcppaux[1].end();++p)
      (cnode->CoData().GetDerivN()[1])[p->first] += (1.0 -alpha) * (p->second);
    for (CI p=dcppaux[2].begin();p!=dcppaux[2].end();++p)
      (cnode->CoData().GetDerivN()[2])[p->first] += (1.0 -alpha) * (p->second);

    // lin alpha
    for (CI p=cnode->CoData().GetAlpha().begin();p!=cnode->CoData().GetAlpha().end();++p)
    {
      (cnode->CoData().GetDerivN()[0])[p->first] += (p->second) * (normalsele[0] - normalcpp[0]);
      (cnode->CoData().GetDerivN()[1])[p->first] += (p->second) * (normalsele[1] - normalcpp[1]);
      (cnode->CoData().GetDerivN()[2])[p->first] += (p->second) * (normalsele[2] - normalcpp[2]);
    }

    // 2D Tangent!
    if (cnode->NumDof()==2)
    {
      // simple definition for txi
      cnode->CoData().txi()[0] = -cnode->MoData().n()[1];
      cnode->CoData().txi()[1] =  cnode->MoData().n()[0];
      cnode->CoData().txi()[2] =  0.0;

      // teta is z-axis
      cnode->CoData().teta()[0] = 0.0;
      cnode->CoData().teta()[1] = 0.0;
      cnode->CoData().teta()[2] = 1.0;

      cnode->CoData().GetDerivTxi()[0].clear();
      cnode->CoData().GetDerivTxi()[1].clear();
      cnode->CoData().GetDerivTxi()[2].clear();


      for (CI p=cnode->CoData().GetDerivN()[1].begin();p!=cnode->CoData().GetDerivN()[1].end();++p)
        (cnode->CoData().GetDerivTxi()[0])[p->first] -= (p->second);
      for (CI p=cnode->CoData().GetDerivN()[0].begin();p!=cnode->CoData().GetDerivN()[0].end();++p)
        (cnode->CoData().GetDerivTxi()[1])[p->first] += (p->second);
    }
    else
      dserror("ERROR: only 2D");



    double length = sqrt(cnode->MoData().n()[0] *cnode->MoData().n()[0] +cnode->MoData().n()[1]*cnode->MoData().n()[1] + cnode->MoData().n()[2]*cnode->MoData().n()[2]);

    if(cnode->Active())
    {
      std::cout << "cnode->MoData().n()[0]= " << cnode->MoData().n()[0] << std::endl;
      std::cout << "cnode->MoData().n()[1]= " << cnode->MoData().n()[1] << std::endl;
      std::cout << "cnode->MoData().n()[2]= " << cnode->MoData().n()[2] << std::endl;
      std::cout << "length= " << length << std::endl;
    }


  }

  return;
}

/*----------------------------------------------------------------------*
 |  scale normals for hybrid formulation                    farah 05/16 |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::ScaleNormals3D()
{
  dserror("ERROR: not yet implemented!");
}


/*----------------------------------------------------------------------*
 |  cpp to line based on averaged nodal normal field        farah 08/16 |
 *----------------------------------------------------------------------*/
double CONTACT::CoInterface::ComputeCPPNormal2D(
    MORTAR::MortarNode& mrtrnode,
    std::vector<MORTAR::MortarElement*> meles,
    double* normal,
    std::vector<GEN::pairedvector<int,double> >& normaltolineLin)
{
  // define tolerance
  const double tol        = 1e-8;
  const double validAngle = 20;

  // distance between node and surface
  double gdist  = 1e12;  // distance
  double gnormal[3] = {0.0, 0.0, 0.0};
  std::vector<GEN::pairedvector<int,double> > glin(3,1); // 1 dummy
  std::set<int> donebeforeMasterCorner;

  bool nodeOnNode = false;    // flag for node on node (corner on corner) setting
  bool pathdependent = false; // flag if we have to check path from last converged check
  double vect[3] = {0.0 ,0.0, 0.0}; // patch

  // calc trajectory and node-to-node distance for corner nodes
  if(mrtrnode.IsOnCorner())
  {
    pathdependent = true;
    CONTACT::CoNode& coNode = dynamic_cast<CONTACT::CoNode&>(mrtrnode);
    if(coNode.Active())
      pathdependent = false;

    // calculate path
    // check trajectory from considered node
    double Posn[3]  = {0.0 ,0.0, 0.0};
    double Posnp[3] = {0.0 ,0.0, 0.0};
    Posn[0] = mrtrnode.X()[0] + mrtrnode.uold()[0];
    Posn[1] = mrtrnode.X()[1] + mrtrnode.uold()[1];
    Posn[2] = mrtrnode.X()[2] + mrtrnode.uold()[2];
    Posnp[0] = mrtrnode.xspatial()[0];
    Posnp[1] = mrtrnode.xspatial()[1];
    Posnp[2] = mrtrnode.xspatial()[2];

    vect[0] = Posnp[0] - Posn[0];
    vect[1] = Posnp[1] - Posn[1];
    vect[2] = Posnp[2] - Posn[2];

    double lvec = sqrt(vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2]);
    if(lvec<1e-12)
    {
      pathdependent = false;
    }
    else
    {
      vect[0] /=lvec;
      vect[1] /=lvec;
      vect[2] /=lvec;
    }

    // Compute normal to node:
    // loop over found eles
    for (size_t ele = 0; ele< meles.size();++ele)
    {
      // get linsize
      int linsize = 0;
      for (int i = 0;i<meles[ele]->NumNode();++i)
      {
        DRT::Node* node = meles[ele]->Nodes()[i];
        if (!node)
          dserror("ERROR: Cannot find master node");
        CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(node);
        linsize += mnode->GetLinsize();

        // if master node is also corner node
        if(mnode->IsOnCorner())
        {
          std::set<int>::iterator iter = donebeforeMasterCorner.find(mnode->Id());

          // if not then create ele
          if (iter != donebeforeMasterCorner.end())
            continue;

          // set master corner node id
          donebeforeMasterCorner.insert(mnode->Id());

          // aux variables
          double dist = 1e12;
          double auxnormal[3] = {0.0, 0.0, 0.0};
          std::vector<GEN::pairedvector<int,double> > auxlin(3,linsize + 1 + meles[ele]->NumNode());

          // compute distance between corners
          dist= ComputeNormalNodeToNode(mrtrnode,
              *mnode,
              auxnormal,
              auxlin);

          // if nodes lying on each other
          if(abs(dist)<1e-12)
          {
            nodeOnNode = true;
            continue;
          }

          // angle between trajectory and normal
          if(pathdependent)
          {
            double auxl  = sqrt(auxnormal[0]*auxnormal[0] + auxnormal[1]*auxnormal[1]);
            if(auxl<1e-12)
              continue;
            double angle = acos(-vect[0] * auxnormal[0]/auxl - vect[1] * auxnormal[1]/auxl );

            angle = 180 * (angle/3.14159265359);
            if(abs(angle) > validAngle)
              continue;
          }

          // get closest valid distance
          if (abs(dist) < abs(gdist))
          {
            gdist = dist;
            gnormal[0] = auxnormal[0];
            gnormal[1] = auxnormal[1];
            gnormal[2] = auxnormal[2];
            glin = auxlin;
          }
        }
      } // end node loop
    } // end element loop
  }

  // loop over found eles
  for (size_t ele = 0; ele< meles.size();++ele)
  {
    bool cornerele = false;
    // get linsize
    int linsize = 0;
    for (int i = 0;i<meles[ele]->NumNode();++i)
    {
      DRT::Node* node = meles[ele]->Nodes()[i];
      if (!node) dserror("ERROR: Cannot find master node");
      CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(node);
      linsize += mnode->GetLinsize();

      if(mnode->IsOnCorner())
        cornerele = true;
    }

    double xi[2] = {0.0 , 0.0};
    double dist = 1e12;
    double auxnormal[3] = {0.0, 0.0, 0.0};
    std::vector<GEN::pairedvector<int,double> > auxlin(3,linsize + 1 + meles[ele]->NumNode());

    // check for nonsmooth mele
    if(cornerele and !nodeOnNode)
    {
      // perform CPP to find normals based on element normals
      MORTAR::MortarProjector::Impl(*meles[ele])->ProjectSNodeByMNormalLin(
          mrtrnode,
          *meles[ele],
          xi,
          auxnormal,
          dist,
          auxlin);
    }
    // compute normal with averaged nodal normal field from master surface
    else
    {
      // perform CPP to find normals based on averaged nodal normal field
      MORTAR::MortarProjector::Impl(*meles[ele])->ProjectSNodeByMNodalNormalLin(
          mrtrnode,
          *meles[ele],
          xi,
          auxnormal,
          dist,
          auxlin);
    }

    // check if found parameter space coordinate is within element domain
    if(meles[ele]->Shape() == DRT::Element::line2 or
       meles[ele]->Shape() == DRT::Element::line3)
    {
      if(-1.0-tol>xi[0] or xi[0]>1.0+tol)
        continue;
    }
    else
    {
      dserror("ERROR: Unknown ele type!");
    }

    // angle between trajectory and normal
    if(pathdependent)
    {
      double auxl  = sqrt(auxnormal[0]*auxnormal[0] + auxnormal[1]*auxnormal[1]);
      if(auxl<1e-12)
        continue;
      double angle = acos(-vect[0] * auxnormal[0]/auxl - vect[1] * auxnormal[1]/auxl );

      angle = 180 * (angle/3.14159265359);
      if(abs(angle) > validAngle)
        continue;
    }

    // get closest valid distance
    if (abs(dist) < abs(gdist))
    {
      gdist = dist;
      gnormal[0] = auxnormal[0];
      gnormal[1] = auxnormal[1];
      gnormal[2] = auxnormal[2];
      glin = auxlin;
    }
  }// end mele loop

  // get the cpp normal
  normal[0] = gnormal[0];
  normal[1] = gnormal[1];
  normal[2] = gnormal[2];
  normaltolineLin = glin;

  return gdist;
}

/*----------------------------------------------------------------------*
 |  cpp to line based on averaged nodal normal field        farah 08/16 |
 *----------------------------------------------------------------------*/
double CONTACT::CoInterface::ComputeCPPNormal3D(
    MORTAR::MortarNode& mrtrnode,
    std::vector<MORTAR::MortarElement*> meles,
    double* normal,
    std::vector<GEN::pairedvector<int,double> >& normaltolineLin)
{
  // define tolerance
  const double tol = 1e-8;

  // distance between node and surface
  bool pathdependent = true;
  const double validAngle = 5.0;
  double gdist  = 1e12;  // distance
  double gnormal[3] = {0.0, 0.0, 0.0};
  std::vector<GEN::pairedvector<int,double> > glin(3,1); // 1 dummy

  //******************************************************
  //             CALC TRAJECTORY
  //******************************************************
  double Posn[3]  = {0.0 ,0.0, 0.0};
  double Posnp[3] = {0.0 ,0.0, 0.0};
  double vect[3]  = {0.0 ,0.0, 0.0};

  Posn[0] = mrtrnode.X()[0] + mrtrnode.uold()[0];
  Posn[1] = mrtrnode.X()[1] + mrtrnode.uold()[1];
  Posn[2] = mrtrnode.X()[2] + mrtrnode.uold()[2];
  Posnp[0] = mrtrnode.xspatial()[0];
  Posnp[1] = mrtrnode.xspatial()[1];
  Posnp[2] = mrtrnode.xspatial()[2];

  vect[0] = Posnp[0] - Posn[0];
  vect[1] = Posnp[1] - Posn[1];
  vect[2] = Posnp[2] - Posn[2];

  double lvec = sqrt(vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2]);
  if(lvec<1e-12)
  {
    pathdependent = false;
  }
  else
  {
    vect[0] /=lvec;
    vect[1] /=lvec;
    vect[2] /=lvec;
  }

  if(dynamic_cast<CONTACT::CoNode&>(mrtrnode).Active())
    pathdependent = false;

  //******************************************************
  //             COMPUTE NORMAL TO SURFACE
  //******************************************************
  // loop over found eles for all geometrical nodes
  for (size_t ele = 0; ele< meles.size();++ele)
  {
    double xi[2] = {0.0 , 0.0};
    double dist = 1e12;
    double auxnormal[3] = {0.0, 0.0, 0.0};
    std::vector<GEN::pairedvector<int,double> > auxlin(3,1000);

    // perform CPP to find normals
    bool success = MORTAR::MortarProjector::Impl(*meles[ele])->ProjectSNodeByMNodalNormalLin(
        mrtrnode,
        *meles[ele],
        xi,
        auxnormal,
        dist,
        auxlin);

    // newton not converged
    if(!success)
      continue;

    // check if found parameter space coordinate is within element domain
    if(meles[ele]->Shape()  == DRT::Element::quad4 or
        meles[ele]->Shape() == DRT::Element::quad8 or
        meles[ele]->Shape() == DRT::Element::quad9)
    {
      if(-1.0-tol>xi[0] or xi[0]>1.0+tol or
         -1.0-tol>xi[1] or xi[1]>1.0+tol)
        continue;
    }
    else if(meles[ele]->Shape() == DRT::Element::tri3 or
            meles[ele]->Shape() == DRT::Element::tri6)
    {
      if(xi[0]<0.0-tol or xi[1]<0.0-tol or
         xi[0]>1.0+tol or xi[1]>1.0+tol or
         xi[0]+xi[1]>1.0+2*tol)
        continue;
    }
    else
    {
      dserror("ERROR: Unknown ele type!");
    }

    // angle between trajectory and normal
    if(pathdependent)
    {
      double auxl  = sqrt(auxnormal[0]*auxnormal[0] + auxnormal[1]*auxnormal[1] + auxnormal[2]*auxnormal[2]);
      if(auxl<1e-12)
        continue;
      double angle = acos(-vect[0] * auxnormal[0]/auxl - vect[1] * auxnormal[1]/auxl -vect[2] * auxnormal[2]/auxl);

      angle = 180 * (angle/3.14159265359);
      if(abs(angle) > validAngle)
        continue;
    }

    if (dist < gdist)
    {
      gdist = dist;
      gnormal[0] = auxnormal[0];
      gnormal[1] = auxnormal[1];
      gnormal[2] = auxnormal[2];
      glin = auxlin;
    }
  }// end mele loop
  //******************************************************
  //             COMPUTE NORMAL TO LINE
  //******************************************************
  if(mrtrnode.IsOnCornerEdge()) // only for edge or corner nodes possible
  {
    // guarantee uniquness
    std::set<std::pair<int,int> > donebefore;

    // calc
    for (size_t ele = 0; ele< meles.size();++ele)
    {
      // loop over master edges -> match node number for quad4
      for(int j = 0; j< meles[ele]->NumNode() ; ++j)
      {
        int nodeIds[2]  = {0,0};
        int nodeLIds[2] = {0,0};

        if(meles[ele]->Shape() == DRT::Element::quad4)
        {
          if(j == 0)
          {
            nodeIds[0] = meles[ele]->NodeIds()[0];
            nodeIds[1] = meles[ele]->NodeIds()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if(j == 1)
          {
            nodeIds[0] = meles[ele]->NodeIds()[1];
            nodeIds[1] = meles[ele]->NodeIds()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if(j == 2)
          {
            nodeIds[0] = meles[ele]->NodeIds()[2];
            nodeIds[1] = meles[ele]->NodeIds()[3];

            nodeLIds[0] = 2;
            nodeLIds[1] = 3;
          }
          else if(j == 3)
          {
            nodeIds[0] = meles[ele]->NodeIds()[3];
            nodeIds[1] = meles[ele]->NodeIds()[0];

            nodeLIds[0] = 3;
            nodeLIds[1] = 0;
          }
          else
            dserror("ERROR: loop counter and edge number do not match!");
        }

        // check if both nodes on edge geometry
        bool node0Edge =
            dynamic_cast<MORTAR::MortarNode*>(meles[ele]->Nodes()[nodeLIds[0]])->IsOnCornerEdge();
        bool node1Edge =
            dynamic_cast<MORTAR::MortarNode*>(meles[ele]->Nodes()[nodeLIds[1]])->IsOnCornerEdge();

        if(!node0Edge or !node1Edge)
          continue;

        //create pair
        std::pair<int,int> actIDs   = std::pair<int,int>(nodeIds[0],nodeIds[1]);
        std::pair<int,int> actIDstw = std::pair<int,int>(nodeIds[1],nodeIds[0]);

        // check if processed before
        std::set<std::pair<int,int> >::iterator iter   = donebefore.find(actIDs);
        std::set<std::pair<int,int> >::iterator itertw = donebefore.find(actIDstw);

         // if not then create ele
         if (iter == donebefore.end() and itertw == donebefore.end() )
         {
           // add to set of processed nodes
           donebefore.insert(actIDs);
           donebefore.insert(actIDstw);

           // create line ele:
           Teuchos::RCP<MORTAR::MortarElement> lineEle = Teuchos::rcp(
                       new MORTAR::MortarElement(
                           j,
                           meles[ele]->Owner(),
                           DRT::Element::line2,
                           2,
                           nodeIds,
                           false));

           // get nodes
           DRT::Node* nodes[2] =
           {meles[ele]->Nodes()[nodeLIds[0]],meles[ele]->Nodes()[nodeLIds[1]]};
           lineEle->BuildNodalPointers(nodes);

           // init data container for dual shapes
           lineEle->InitializeDataContainer();

           // call cpp function for edge to edge

           double dist = 1e12;
           double auxnormal[3] = {0.0, 0.0, 0.0};
           std::vector<GEN::pairedvector<int,double> > auxlin(3,100 + 1 + meles[ele]->NumNode());

           // compute distance between node and edge
           dist= ComputeNormalNodeToEdge(mrtrnode,
               *lineEle,
               auxnormal,
               auxlin);

           // angle between trajectory and normal
           if(pathdependent)
           {
             double auxl  = sqrt(auxnormal[0]*auxnormal[0] + auxnormal[1]*auxnormal[1] + auxnormal[2]*auxnormal[2]);
             if(auxl<1e-12)
               continue;
             double angle = acos(-vect[0] * auxnormal[0]/auxl - vect[1] * auxnormal[1]/auxl -vect[2] * auxnormal[2]/auxl);

             angle = 180 * (angle/3.14159265359);
             if(abs(angle) > validAngle)
               continue;
           }

           if (dist <= gdist+tol)
           {
             std::cout << "CLOSE TO EDGE!!!" << std::endl;
             gdist = dist;
             gnormal[0] = auxnormal[0];
             gnormal[1] = auxnormal[1];
             gnormal[2] = auxnormal[2];
             glin = auxlin;
           }

         }
      } // end mele node loop

    }
  }

  //******************************************************
  //             COMPUTE NORMAL TO NODE
  //******************************************************
  if(mrtrnode.IsOnCorner()) // only for corner nodes possible
  {
    std::set<int> donebeforeMasterCorner;

    for (size_t ele = 0; ele< meles.size();++ele)
    {
      // get linsize
      int linsize = 0;
      for (int i = 0;i<meles[ele]->NumNode();++i)
      {
        DRT::Node* node = meles[ele]->Nodes()[i];
        if (!node)
          dserror("ERROR: Cannot find master node");
        CONTACT::CoNode* mnode = dynamic_cast<CONTACT::CoNode*>(node);
        linsize += mnode->GetLinsize();

        // if master node is also corner node
        if(mnode->IsOnCorner())
        {
          std::set<int>::iterator iter = donebeforeMasterCorner.find(mnode->Id());

          // if not then create ele
          if (iter != donebeforeMasterCorner.end())
            continue;

          // set master corner node id
          donebeforeMasterCorner.insert(mnode->Id());

          double dist = 1e12;
          double auxnormal[3] = {0.0, 0.0, 0.0};
          std::vector<GEN::pairedvector<int,double> > auxlin(3,linsize + 1 + meles[ele]->NumNode());

          // compute distance between corners
          dist= ComputeNormalNodeToNode(mrtrnode,
              *mnode,
              auxnormal,
              auxlin);

          // angle between trajectory and normal
          if(pathdependent)
          {
            double auxl  = sqrt(auxnormal[0]*auxnormal[0] + auxnormal[1]*auxnormal[1] + auxnormal[2]*auxnormal[2]);
            if(auxl<1e-12)
              continue;
            double angle = acos(-vect[0] * auxnormal[0]/auxl - vect[1] * auxnormal[1]/auxl -vect[2] * auxnormal[2]/auxl);

            angle = 180 * (angle/3.14159265359);
            if(abs(angle) > validAngle)
              continue;
          }

          if (dist < gdist)
          {
            gdist = dist;
            gnormal[0] = auxnormal[0];
            gnormal[1] = auxnormal[1];
            gnormal[2] = auxnormal[2];
            glin = auxlin;
          }
        }
      }
    }
  }

  //******************************************************
  //             FINAL STORAGE
  //******************************************************

  // get the cpp normal
  normal[0] = gnormal[0];
  normal[1] = gnormal[1];
  normal[2] = gnormal[2];
  normaltolineLin = glin;

  // bye bye
  return gdist;
}

/*----------------------------------------------------------------------*
 |  cpp to line based on averaged nodal normal field        farah 05/16 |
 *----------------------------------------------------------------------*/
double CONTACT::CoInterface::ComputeCPPNormal(
    MORTAR::MortarNode& mrtrnode,
    std::vector<MORTAR::MortarElement*> meles,
    double* normal,
    std::vector<GEN::pairedvector<int,double> >& normaltolineLin)
{
  // define distance
  double gdist = 1e12;

  //===================================================================
  //===================================================================
  //                           2D case
  //===================================================================
  //===================================================================
  if (Dim()==2)
  {
    gdist= ComputeCPPNormal2D(mrtrnode,
                              meles,
                              normal,
                              normaltolineLin);
  }
  //===================================================================
  //===================================================================
  //                           3D case
  //===================================================================
  //===================================================================
  else if(Dim()==3)
  {
    gdist= ComputeCPPNormal3D(mrtrnode,
                              meles,
                              normal,
                              normaltolineLin);
  }
  //===================================================================
  //===================================================================
  //                           Invalid
  //===================================================================
  //===================================================================
  else
  {
    dserror("ERROR: invalid dimension!");
  }

  // return distance
  return gdist;
}

/*----------------------------------------------------------------------*
 |  set cpp normal                                           farah 01/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::SetCPPNormal(
    MORTAR::MortarNode& snode,
    double* normal,
    std::vector<GEN::pairedvector<int,double> >& normallin)
{
  CoNode& cnode = dynamic_cast<CoNode&>(snode);

  const double length = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  if(length<1e-12)
    dserror("ERROR: normal length is zero!!!");

  // negative sign because it is a master normal!
  cnode.MoData().n()[0] = -normal[0]/length;
  cnode.MoData().n()[1] = -normal[1]/length;
  cnode.MoData().n()[2] = -normal[2]/length;

//  if (cnode.IsOnEdge())
//    std::cout << "normal =  " << cnode.MoData().n()[0] << "  "<< cnode.MoData().n()[1] << "  "<< cnode.MoData().n()[2] << std::endl;

  // prepare nodal storage maps for derivative
  if ((int)cnode.CoData().GetDerivN().size()==0)
    cnode.CoData().GetDerivN().resize(3,normallin[0].size()*3);
  if ((int)cnode.CoData().GetDerivTxi().size()==0)
    cnode.CoData().GetDerivTxi().resize(3,normallin[0].size()*3);
  if ((int)cnode.CoData().GetDerivTeta().size()==0)
    cnode.CoData().GetDerivTeta().resize(3,normallin[0].size()*3);

  // init tangent length
  double ltxi = -1.0;

  //------------------------------------------------------------------
  // 2D Tangent!
  if (cnode.NumDof()==2)
  {
    // simple definition for txi
    cnode.CoData().txi()[0] = -cnode.MoData().n()[1];
    cnode.CoData().txi()[1] =  cnode.MoData().n()[0];
    cnode.CoData().txi()[2] =  0.0;

    // teta is z-axis
    cnode.CoData().teta()[0] = 0.0;
    cnode.CoData().teta()[1] = 0.0;
    cnode.CoData().teta()[2] = 1.0;
  }
  // 3D Tangent!
  else
  {
    if (abs(cnode.MoData().n()[0])>1.0e-6 || abs(cnode.MoData().n()[1])>1.0e-6 )
    {
      cnode.CoData().txi()[0]=-cnode.MoData().n()[1];
      cnode.CoData().txi()[1]=cnode.MoData().n()[0];
      cnode.CoData().txi()[2]=0.0;
    }
    else
    {
      cnode.CoData().txi()[0]=0.0;
      cnode.CoData().txi()[1]=-cnode.MoData().n()[2];
      cnode.CoData().txi()[2]=cnode.MoData().n()[1];
    }

    ltxi = sqrt(cnode.CoData().txi()[0]*cnode.CoData().txi()[0]+cnode.CoData().txi()[1]*cnode.CoData().txi()[1]+cnode.CoData().txi()[2]*cnode.CoData().txi()[2]);
    if(ltxi<1e-12)
      dserror("ERROR: tangent txi length is zero!!!");
    for (int j=0;j<3;++j)
      cnode.CoData().txi()[j]/=ltxi;

    // teta follows from corkscrew rule (teta = n x txi)
    cnode.CoData().teta()[0] = cnode.MoData().n()[1]*cnode.CoData().txi()[2]-cnode.MoData().n()[2]*cnode.CoData().txi()[1];
    cnode.CoData().teta()[1] = cnode.MoData().n()[2]*cnode.CoData().txi()[0]-cnode.MoData().n()[0]*cnode.CoData().txi()[2];
    cnode.CoData().teta()[2] = cnode.MoData().n()[0]*cnode.CoData().txi()[1]-cnode.MoData().n()[1]*cnode.CoData().txi()[0];
  }


  //------------------------------------------------------------------
  typedef GEN::pairedvector<int,double>::const_iterator CI;

  for (CI p=normallin[0].begin();p!=normallin[0].end();++p)
    (cnode.CoData().GetDerivN()[0])[p->first] -= (p->second);
  for (CI p=normallin[1].begin();p!=normallin[1].end();++p)
    (cnode.CoData().GetDerivN()[1])[p->first] -= (p->second);
  for (CI p=normallin[2].begin();p!=normallin[2].end();++p)
    (cnode.CoData().GetDerivN()[2])[p->first] -= (p->second);

  // normalize directional derivative
  // (length differs for weighted/unweighted case bot not the procedure!)
  // (be careful with reference / copy of derivative maps!)
  typedef GEN::pairedvector<int,double>::const_iterator CI;
  GEN::pairedvector<int,double>& derivnx = cnode.CoData().GetDerivN()[0];
  GEN::pairedvector<int,double>& derivny = cnode.CoData().GetDerivN()[1];
  GEN::pairedvector<int,double>& derivnz = cnode.CoData().GetDerivN()[2];
  GEN::pairedvector<int,double> cderivnx = cnode.CoData().GetDerivN()[0];
  GEN::pairedvector<int,double> cderivny = cnode.CoData().GetDerivN()[1];
  GEN::pairedvector<int,double> cderivnz = cnode.CoData().GetDerivN()[2];
  const double nxnx = cnode.MoData().n()[0] * cnode.MoData().n()[0];
  const double nxny = cnode.MoData().n()[0] * cnode.MoData().n()[1];
  const double nxnz = cnode.MoData().n()[0] * cnode.MoData().n()[2];
  const double nyny = cnode.MoData().n()[1] * cnode.MoData().n()[1];
  const double nynz = cnode.MoData().n()[1] * cnode.MoData().n()[2];
  const double nznz = cnode.MoData().n()[2] * cnode.MoData().n()[2];

  // build a vector with all keys from x,y,z maps
  // (we need this in order not to miss any entry!)
  std::vector<int> allkeysn;
  for (CI p=derivnx.begin();p!=derivnx.end();++p)
  {
    bool found = false;
    for (int j=0;j<(int)allkeysn.size();++j)
      if ((p->first)==allkeysn[j]) found = true;
    if (!found) allkeysn.push_back(p->first);

  }
  for (CI p=derivny.begin();p!=derivny.end();++p)
  {
    bool found = false;
    for (int j=0;j<(int)allkeysn.size();++j)
      if ((p->first)==allkeysn[j]) found = true;
    if (!found) allkeysn.push_back(p->first);

  }
  for (CI p=derivnz.begin();p!=derivnz.end();++p)
  {
    bool found = false;
    for (int j=0;j<(int)allkeysn.size();++j)
      if ((p->first)==allkeysn[j]) found = true;
    if (!found) allkeysn.push_back(p->first);
  }

  // normalize x-components
  for (int j=0;j<(int)allkeysn.size();++j)
  {
    double val = cderivnx[allkeysn[j]];
    derivnx[allkeysn[j]] = (val-nxnx*val-nxny*cderivny[allkeysn[j]]-nxnz*cderivnz[allkeysn[j]])/length;
  }

  // normalize y-components
  for (int j=0;j<(int)allkeysn.size();++j)
  {
    double val = cderivny[allkeysn[j]];
    derivny[allkeysn[j]] = (val-nxny*cderivnx[allkeysn[j]]-nyny*val-nynz*cderivnz[allkeysn[j]])/length;
  }

  // normalize z-components
  for (int j=0;j<(int)allkeysn.size();++j)
  {
    double val = cderivnz[allkeysn[j]];
    derivnz[allkeysn[j]] = (val-nxnz*cderivnx[allkeysn[j]]-nynz*cderivny[allkeysn[j]]-nznz*val)/length;
  }

  //------------------------------------------------------------------
  // 2D Tangent!
  if (cnode.NumDof()==2)
  {
    for (CI p=cnode.CoData().GetDerivN()[1].begin();p!=cnode.CoData().GetDerivN()[1].end();++p)
      (cnode.CoData().GetDerivTxi()[0])[p->first] -= (p->second);
    for (CI p=cnode.CoData().GetDerivN()[0].begin();p!=cnode.CoData().GetDerivN()[0].end();++p)
      (cnode.CoData().GetDerivTxi()[1])[p->first] += (p->second);
  }
  // 3D Tangent!
  else
  {
    // unnormalized tangent derivative txi
    // use definitions for txi from BuildAveragedNormal()
    if (abs(cnode.MoData().n()[0])>1.0e-6 || abs(cnode.MoData().n()[1])>1.0e-6)
    {
      GEN::pairedvector<int,double>& derivtxix = cnode.CoData().GetDerivTxi()[0];
      GEN::pairedvector<int,double>& derivtxiy = cnode.CoData().GetDerivTxi()[1];

      for (CI p=derivny.begin();p!=derivny.end();++p)
        derivtxix[p->first] -= (p->second);

      for (CI p=derivnx.begin();p!=derivnx.end();++p)
        derivtxiy[p->first] += (p->second);
    }
    else
    {
      GEN::pairedvector<int,double>& derivtxiy = cnode.CoData().GetDerivTxi()[1];
      GEN::pairedvector<int,double>& derivtxiz = cnode.CoData().GetDerivTxi()[2];

      for (CI p=derivnz.begin();p!=derivnz.end();++p)
        derivtxiy[p->first] -= (p->second);

      for (CI p=derivny.begin();p!=derivny.end();++p)
        derivtxiz[p->first] += (p->second);
    }

    if(ltxi<1e-12)
      dserror("ERROR: tangent txi length is zero!!!");

    // normalize txi directional derivative
    // (identical to normalization of normal derivative)
    typedef GEN::pairedvector<int,double>::const_iterator CI;
    GEN::pairedvector<int,double>& derivtxix = cnode.CoData().GetDerivTxi()[0];
    GEN::pairedvector<int,double>& derivtxiy = cnode.CoData().GetDerivTxi()[1];
    GEN::pairedvector<int,double>& derivtxiz = cnode.CoData().GetDerivTxi()[2];
    GEN::pairedvector<int,double> cderivtxix = cnode.CoData().GetDerivTxi()[0];
    GEN::pairedvector<int,double> cderivtxiy = cnode.CoData().GetDerivTxi()[1];
    GEN::pairedvector<int,double> cderivtxiz = cnode.CoData().GetDerivTxi()[2];
    const double txtx = cnode.CoData().txi()[0] * cnode.CoData().txi()[0];
    const double txty = cnode.CoData().txi()[0] * cnode.CoData().txi()[1];
    const double txtz = cnode.CoData().txi()[0] * cnode.CoData().txi()[2];
    const double tyty = cnode.CoData().txi()[1] * cnode.CoData().txi()[1];
    const double tytz = cnode.CoData().txi()[1] * cnode.CoData().txi()[2];
    const double tztz = cnode.CoData().txi()[2] * cnode.CoData().txi()[2];

    // build a vector with all keys from x,y,z maps
    // (we need this in order not to miss any entry!)
    std::vector<int> allkeyst;
    for (CI p=derivtxix.begin();p!=derivtxix.end();++p)
    {
      bool found = false;
      for (int j=0;j<(int)allkeyst.size();++j)
        if ((p->first)==allkeyst[j]) found = true;
      if (!found) allkeyst.push_back(p->first);

    }
    for (CI p=derivtxiy.begin();p!=derivtxiy.end();++p)
    {
      bool found = false;
      for (int j=0;j<(int)allkeyst.size();++j)
        if ((p->first)==allkeyst[j]) found = true;
      if (!found) allkeyst.push_back(p->first);

    }
    for (CI p=derivtxiz.begin();p!=derivtxiz.end();++p)
    {
      bool found = false;
      for (int j=0;j<(int)allkeyst.size();++j)
        if ((p->first)==allkeyst[j]) found = true;
      if (!found) allkeyst.push_back(p->first);
    }

    // normalize x-components
    for (int j=0;j<(int)allkeyst.size();++j)
    {
      double val = cderivtxix[allkeyst[j]];
      derivtxix[allkeyst[j]] = (val-txtx*val-txty*cderivtxiy[allkeyst[j]]-txtz*cderivtxiz[allkeyst[j]])/ltxi;
    }

    // normalize y-components
    for (int j=0;j<(int)allkeyst.size();++j)
    {
      double val =cderivtxiy[allkeyst[j]];
      derivtxiy[allkeyst[j]] = (val-txty*cderivtxix[allkeyst[j]]-tyty*val-tytz*cderivtxiz[allkeyst[j]])/ltxi;
    }

    // normalize z-components
    for (int j=0;j<(int)allkeyst.size();++j)
    {
      double val = cderivtxiz[allkeyst[j]];
      derivtxiz[allkeyst[j]] = (val-txtz*cderivtxix[allkeyst[j]]-tytz*cderivtxiy[allkeyst[j]]-tztz*val)/ltxi;
    }

    // get normalized tangent derivative teta
    // use corkscrew rule from BuildAveragedNormal()
    GEN::pairedvector<int,double>& derivtetax = cnode.CoData().GetDerivTeta()[0];
    GEN::pairedvector<int,double>& derivtetay = cnode.CoData().GetDerivTeta()[1];
    GEN::pairedvector<int,double>& derivtetaz = cnode.CoData().GetDerivTeta()[2];

    for (CI p=derivnx.begin();p!=derivnx.end();++p)
    {
      derivtetay[p->first] -= cnode.CoData().txi()[2]*(p->second);
      derivtetaz[p->first] += cnode.CoData().txi()[1]*(p->second);
    }
    for (CI p=derivny.begin();p!=derivny.end();++p)
    {
      derivtetax[p->first] += cnode.CoData().txi()[2]*(p->second);
      derivtetaz[p->first] -= cnode.CoData().txi()[0]*(p->second);
    }
    for (CI p=derivnz.begin();p!=derivnz.end();++p)
    {
      derivtetax[p->first] -= cnode.CoData().txi()[1]*(p->second);
      derivtetay[p->first] += cnode.CoData().txi()[0]*(p->second);
    }
    for (CI p=derivtxix.begin();p!=derivtxix.end();++p)
    {
      derivtetay[p->first] += cnode.MoData().n()[2]*(p->second);
      derivtetaz[p->first] -= cnode.MoData().n()[1]*(p->second);
    }
    for (CI p=derivtxiy.begin();p!=derivtxiy.end();++p)
    {
      derivtetax[p->first] -= cnode.MoData().n()[2]*(p->second);
      derivtetaz[p->first] += cnode.MoData().n()[0]*(p->second);
    }
    for (CI p=derivtxiz.begin();p!=derivtxiz.end();++p)
    {
      derivtetax[p->first] += cnode.MoData().n()[1]*(p->second);
      derivtetay[p->first] -= cnode.MoData().n()[0]*(p->second);
    }

    // OLD VERSION:

//    if (abs(cnode.MoData().n()[0])>1.0e-6 || abs(cnode.MoData().n()[1])>1.0e-6 )
//    {
//      for (CI p=cnode.CoData().GetDerivN()[1].begin();p!=cnode.CoData().GetDerivN()[1].end();++p)
//        (cnode.CoData().GetDerivTxi()[0])[p->first] -= (p->second);
//      for (CI p=cnode.CoData().GetDerivN()[0].begin();p!=cnode.CoData().GetDerivN()[0].end();++p)
//        (cnode.CoData().GetDerivTxi()[1])[p->first] += (p->second);
//    }
//    else
//    {
//      for (CI p=cnode.CoData().GetDerivN()[2].begin();p!=cnode.CoData().GetDerivN()[2].end();++p)
//        (cnode.CoData().GetDerivTxi()[1])[p->first] -= (p->second);
//      for (CI p=cnode.CoData().GetDerivN()[1].begin();p!=cnode.CoData().GetDerivN()[1].end();++p)
//        (cnode.CoData().GetDerivTxi()[2])[p->first] += (p->second);
//    }
//
//    double ltxi = sqrt(cnode.CoData().txi()[0]*cnode.CoData().txi()[0]+cnode.CoData().txi()[1]*cnode.CoData().txi()[1]+cnode.CoData().txi()[2]*cnode.CoData().txi()[2]);
//    for (int j=0;j<3;++j) cnode.CoData().txi()[j]/=ltxi;

//     //teta follows from corkscrew rule (teta = n x txi)
//    cnode.CoData().teta()[0] = cnode.MoData().n()[1]*cnode.CoData().txi()[2]-cnode.MoData().n()[2]*cnode.CoData().txi()[1];
//    cnode.CoData().teta()[1] = cnode.MoData().n()[2]*cnode.CoData().txi()[0]-cnode.MoData().n()[0]*cnode.CoData().txi()[2];
//    cnode.CoData().teta()[2] = cnode.MoData().n()[0]*cnode.CoData().txi()[1]-cnode.MoData().n()[1]*cnode.CoData().txi()[0];
//
//    for (CI p=cnode.CoData().GetDerivN()[1].begin();p!=cnode.CoData().GetDerivN()[1].end();++p)
//      (cnode.CoData().GetDerivTeta()[0])[p->first] += (p->second) * cnode.CoData().txi()[2];
//    for (CI p=cnode.CoData().GetDerivTxi()[2].begin();p!=cnode.CoData().GetDerivTxi()[2].end();++p)
//      (cnode.CoData().GetDerivTeta()[0])[p->first] += (p->second) * cnode.MoData().n()[1];
//    for (CI p=cnode.CoData().GetDerivN()[2].begin();p!=cnode.CoData().GetDerivN()[2].end();++p)
//      (cnode.CoData().GetDerivTeta()[0])[p->first] -= (p->second) * cnode.CoData().txi()[1];
//    for (CI p=cnode.CoData().GetDerivTxi()[1].begin();p!=cnode.CoData().GetDerivTxi()[1].end();++p)
//      (cnode.CoData().GetDerivTeta()[0])[p->first] -= (p->second) * cnode.MoData().n()[2];
//
//    for (CI p=cnode.CoData().GetDerivN()[2].begin();p!=cnode.CoData().GetDerivN()[2].end();++p)
//      (cnode.CoData().GetDerivTeta()[1])[p->first] += (p->second) * cnode.CoData().txi()[0];
//    for (CI p=cnode.CoData().GetDerivTxi()[0].begin();p!=cnode.CoData().GetDerivTxi()[0].end();++p)
//      (cnode.CoData().GetDerivTeta()[1])[p->first] += (p->second) * cnode.MoData().n()[2];
//    for (CI p=cnode.CoData().GetDerivN()[0].begin();p!=cnode.CoData().GetDerivN()[0].end();++p)
//      (cnode.CoData().GetDerivTeta()[1])[p->first] -= (p->second) * cnode.CoData().txi()[2];
//    for (CI p=cnode.CoData().GetDerivTxi()[2].begin();p!=cnode.CoData().GetDerivTxi()[2].end();++p)
//      (cnode.CoData().GetDerivTeta()[1])[p->first] -= (p->second) * cnode.MoData().n()[0];
//
//    for (CI p=cnode.CoData().GetDerivN()[0].begin();p!=cnode.CoData().GetDerivN()[0].end();++p)
//      (cnode.CoData().GetDerivTeta()[2])[p->first] += (p->second) * cnode.CoData().txi()[1];
//    for (CI p=cnode.CoData().GetDerivTxi()[1].begin();p!=cnode.CoData().GetDerivTxi()[1].end();++p)
//      (cnode.CoData().GetDerivTeta()[2])[p->first] += (p->second) * cnode.MoData().n()[0];
//    for (CI p=cnode.CoData().GetDerivN()[1].begin();p!=cnode.CoData().GetDerivN()[1].end();++p)
//      (cnode.CoData().GetDerivTeta()[2])[p->first] -= (p->second) * cnode.CoData().txi()[0];
//    for (CI p=cnode.CoData().GetDerivTxi()[0].begin();p!=cnode.CoData().GetDerivTxi()[0].end();++p)
//      (cnode.CoData().GetDerivTeta()[2])[p->first] -= (p->second) * cnode.MoData().n()[1];
  }

  return;
}


/*----------------------------------------------------------------------*
 |  export nodal normals (public)                             popp 11/10|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::ExportNodalNormals() const
{
  // create empty data objects
  std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > triad;

  std::map<int,std::vector<int> > n_x_key;
  std::map<int,std::vector<int> > n_y_key;
  std::map<int,std::vector<int> > n_z_key;
  std::map<int,std::vector<int> > txi_x_key;
  std::map<int,std::vector<int> > txi_y_key;
  std::map<int,std::vector<int> > txi_z_key;
  std::map<int,std::vector<int> > teta_x_key;
  std::map<int,std::vector<int> > teta_y_key;
  std::map<int,std::vector<int> > teta_z_key;

  std::map<int,std::vector<double> > n_x_val;
  std::map<int,std::vector<double> > n_y_val;
  std::map<int,std::vector<double> > n_z_val;
  std::map<int,std::vector<double> > txi_x_val;
  std::map<int,std::vector<double> > txi_y_val;
  std::map<int,std::vector<double> > txi_z_val;
  std::map<int,std::vector<double> > teta_x_val;
  std::map<int,std::vector<double> > teta_y_val;
  std::map<int,std::vector<double> > teta_z_val;

  std::map<int,double>::iterator iter;
  GEN::pairedvector<int,double>::iterator _iter;

  // build info on row map
  for(int i=0; i<snoderowmapbound_->NumMyElements();++i)
  {
    int gid = snoderowmapbound_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    // fill nodal matrix
    Teuchos::RCP<Epetra_SerialDenseMatrix> loc = Teuchos::rcp(new Epetra_SerialDenseMatrix(3,3));
    (*loc)(0,0) = cnode->MoData().n()[0];
    (*loc)(1,0) = cnode->MoData().n()[1];
    (*loc)(2,0) = cnode->MoData().n()[2];
    (*loc)(0,1) = cnode->CoData().txi()[0];
    (*loc)(1,1) = cnode->CoData().txi()[1];
    (*loc)(2,1) = cnode->CoData().txi()[2];
    (*loc)(0,2) = cnode->CoData().teta()[0];
    (*loc)(1,2) = cnode->CoData().teta()[1];
    (*loc)(2,2) = cnode->CoData().teta()[2];

    triad[gid] = loc;

    // fill nodal derivative vectors
    std::vector<GEN::pairedvector<int,double> >& derivn    = cnode->CoData().GetDerivN();
    std::vector<GEN::pairedvector<int,double> >& derivtxi  = cnode->CoData().GetDerivTxi();
    std::vector<GEN::pairedvector<int,double> >& derivteta = cnode->CoData().GetDerivTeta();

    for(_iter=derivn[0].begin();_iter!=derivn[0].end();++_iter)
    {
      n_x_key[gid].push_back(_iter->first);
      n_x_val[gid].push_back(_iter->second);
    }
    for(_iter=derivn[1].begin();_iter!=derivn[1].end();++_iter)
    {
      n_y_key[gid].push_back(_iter->first);
      n_y_val[gid].push_back(_iter->second);
    }
    for(_iter=derivn[2].begin();_iter!=derivn[2].end();++_iter)
    {
      n_z_key[gid].push_back(_iter->first);
      n_z_val[gid].push_back(_iter->second);
    }

    for(_iter=derivtxi[0].begin();_iter!=derivtxi[0].end();++_iter)
    {
      txi_x_key[gid].push_back(_iter->first);
      txi_x_val[gid].push_back(_iter->second);
    }
    for(_iter=derivtxi[1].begin();_iter!=derivtxi[1].end();++_iter)
    {
      txi_y_key[gid].push_back(_iter->first);
      txi_y_val[gid].push_back(_iter->second);
    }
    for(_iter=derivtxi[2].begin();_iter!=derivtxi[2].end();++_iter)
    {
      txi_z_key[gid].push_back(_iter->first);
      txi_z_val[gid].push_back(_iter->second);
    }

    for(_iter=derivteta[0].begin();_iter!=derivteta[0].end();++_iter)
    {
      teta_x_key[gid].push_back(_iter->first);
      teta_x_val[gid].push_back(_iter->second);
    }
    for(_iter=derivteta[1].begin();_iter!=derivteta[1].end();++_iter)
    {
      teta_y_key[gid].push_back(_iter->first);
      teta_y_val[gid].push_back(_iter->second);
    }
    for(_iter=derivteta[2].begin();_iter!=derivteta[2].end();++_iter)
    {
      teta_z_key[gid].push_back(_iter->first);
      teta_z_val[gid].push_back(_iter->second);
    }
  }


  // communicate from slave node row to column map
  DRT::Exporter& ex = idata_.Exporter();

  ex.Export(triad);

  ex.Export(n_x_key);
  ex.Export(n_x_val);
  ex.Export(n_y_key);
  ex.Export(n_y_val);
  ex.Export(n_z_key);
  ex.Export(n_z_val);

  ex.Export(txi_x_key);
  ex.Export(txi_x_val);
  ex.Export(txi_y_key);
  ex.Export(txi_y_val);
  ex.Export(txi_z_key);
  ex.Export(txi_z_val);

  ex.Export(teta_x_key);
  ex.Export(teta_x_val);
  ex.Export(teta_y_key);
  ex.Export(teta_y_val);
  ex.Export(teta_z_key);
  ex.Export(teta_z_val);

  // extract info on column map
  for(int i=0; i<snodecolmapbound_->NumMyElements();++i)
  {
    // only do something for ghosted nodes
    int gid = snodecolmapbound_->GID(i);
    if (snoderowmapbound_->MyGID(gid)) continue;

    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);
    int linsize = cnode->GetLinsize()+(int)(n_x_key[gid].size());

    // extract info
    Teuchos::RCP<Epetra_SerialDenseMatrix> loc = triad[gid];
    cnode->MoData().n()[0]    = (*loc)(0,0);
    cnode->MoData().n()[1]    = (*loc)(1,0);
    cnode->MoData().n()[2]    = (*loc)(2,0);
    cnode->CoData().txi()[0]  = (*loc)(0,1);
    cnode->CoData().txi()[1]  = (*loc)(1,1);
    cnode->CoData().txi()[2]  = (*loc)(2,1);
    cnode->CoData().teta()[0] = (*loc)(0,2);
    cnode->CoData().teta()[1] = (*loc)(1,2);
    cnode->CoData().teta()[2] = (*loc)(2,2);

    // extract derivative info
    std::vector<GEN::pairedvector<int,double> >& derivn    = cnode->CoData().GetDerivN();
    std::vector<GEN::pairedvector<int,double> >& derivtxi  = cnode->CoData().GetDerivTxi();
    std::vector<GEN::pairedvector<int,double> >& derivteta = cnode->CoData().GetDerivTeta();

    for (int k=0;k<(int)(derivn.size());++k)
      derivn[k].clear();
    derivn.resize(3,linsize);
    for (int k=0;k<(int)(derivtxi.size());++k)
      derivtxi[k].clear();
    derivtxi.resize(3,linsize);
    for (int k=0;k<(int)(derivteta.size());++k)
      derivteta[k].clear();
    derivteta.resize(3,linsize);

    for (int k=0;k<(int)(n_x_key[gid].size());++k)
      (cnode->CoData().GetDerivN()[0])[n_x_key[gid][k]] = n_x_val[gid][k];
    for (int k=0;k<(int)(n_y_key[gid].size());++k)
      (cnode->CoData().GetDerivN()[1])[n_y_key[gid][k]] = n_y_val[gid][k];
    for (int k=0;k<(int)(n_z_key[gid].size());++k)
      (cnode->CoData().GetDerivN()[2])[n_z_key[gid][k]] = n_z_val[gid][k];

    for (int k=0;k<(int)(txi_x_key[gid].size());++k)
      (cnode->CoData().GetDerivTxi()[0])[txi_x_key[gid][k]] = txi_x_val[gid][k];
    for (int k=0;k<(int)(txi_y_key[gid].size());++k)
      (cnode->CoData().GetDerivTxi()[1])[txi_y_key[gid][k]] = txi_y_val[gid][k];
    for (int k=0;k<(int)(txi_z_key[gid].size());++k)
      (cnode->CoData().GetDerivTxi()[2])[txi_z_key[gid][k]] = txi_z_val[gid][k];

    for (int k=0;k<(int)(teta_x_key[gid].size());++k)
      (cnode->CoData().GetDerivTeta()[0])[teta_x_key[gid][k]] = teta_x_val[gid][k];
    for (int k=0;k<(int)(teta_y_key[gid].size());++k)
      (cnode->CoData().GetDerivTeta()[1])[teta_y_key[gid][k]] = teta_y_val[gid][k];
    for (int k=0;k<(int)(teta_z_key[gid].size());++k)
      (cnode->CoData().GetDerivTeta()[2])[teta_z_key[gid][k]] = teta_z_val[gid][k];
  }

  // free memory
  triad.clear();

  n_x_key.clear();
  n_y_key.clear();
  n_z_key.clear();
  txi_x_key.clear();
  txi_y_key.clear();
  txi_z_key.clear();
  teta_x_key.clear();
  teta_y_key.clear();
  teta_z_key.clear();

  n_x_val.clear();
  n_y_val.clear();
  n_z_val.clear();
  txi_x_val.clear();
  txi_y_val.clear();
  txi_z_val.clear();
  teta_x_val.clear();
  teta_y_val.clear();
  teta_z_val.clear();

  /*std::cout << "---------- nodal normals ----------------------------------------------------------" << std::endl;

  // print nodal normals
  for (int p=0;p<Comm().NumProc();++p)
  {
    // one proc after the other
    if (p==Comm().MyPID())
    {
      std::cout << "\n*****\nPROC " << p << "\n*****" << std::endl;
      for(int i=0; i<snodecolmapbound_->NumMyElements();++i)
      {
        int gid = snodecolmapbound_->GID(i);
        DRT::Node* node = idiscret_->gNode(gid);
        if (!node) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* cnode = dynamic_cast<CoNode*>(node);

        // print normal and tangents at each slave node
        std::cout << "Proc: " << p << " Node: " << gid << " Owner: " << cnode->Owner()
             << " Normal: " << cnode->MoData().n()[0]
             << " " << cnode->MoData().n()[1] << " " << cnode->MoData().n()[2] << std::endl;
        std::cout << "Proc: " << p << " Node: " << gid << " Owner: " << cnode->Owner()
             << " TXi: " << cnode->CoData().txi()[0]
             << " " << cnode->CoData().txi()[1] << " " << cnode->CoData().txi()[2] << std::endl;
        std::cout << "Proc: " << p << " Node: " << gid << " Owner: " << cnode->Owner()
             << " TEta: " << cnode->CoData().teta()[0]
             << " " << cnode->CoData().teta()[1] << " " << cnode->CoData().teta()[2] << std::endl;

        // print linearizations at each slave node
        std::cout << "Proc: " << p << " Node: " << gid  << " Owner: " << cnode->Owner() << " LinN: ";
        for (_iter=cnode->CoData().GetDerivN()[0].begin();_iter!=cnode->CoData().GetDerivN()[0].end();++_iter)
          std::cout << "\n" << _iter->first << "\t" << _iter->second;
        std::cout << std::endl;
        std::cout << "Proc: " << p << " Node: " << gid  << " Owner: " << cnode->Owner() << " LinTxi: ";
        for (_iter=cnode->CoData().GetDerivTxi()[0].begin();_iter!=cnode->CoData().GetDerivTxi()[0].end();++_iter)
          std::cout << "\n" << _iter->first << "\t" << _iter->second;
        std::cout << std::endl;
        std::cout << "Proc: " << p << " Node: " << gid  << " Owner: " << cnode->Owner() << " LinTeta: ";
        for (_iter=cnode->CoData().GetDerivTeta()[0].begin();_iter!=cnode->CoData().GetDerivTeta()[0].end();++_iter)
          std::cout << "\n" << _iter->first << "\t" << _iter->second;
        std::cout << std::endl;
      }
      std::cout << std::endl << std::endl;
    }

    // barrier
    Comm().Barrier();
  }*/

  return;
}

/*----------------------------------------------------------------------*
 |  Search for potentially contacting sl/ma pairs (public)    popp 10/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::EvaluateSearchBinarytree()
{
  // ***WARNING:*** This is commented out here, as UpdateMasterSlaveSets()
  // needs all the procs around, not only the interface local ones!
  // if (!lComm()) return true;

  // *********************************************************************
  // Possible versions for self contact:
  // *********************************************************************
  //
  // 1) Combined Update and Search
  // -> In this case we have to call SearchContactCombined(), which
  //    does both top-down update (where necessary) and search. Then
  //    the dynamics master/slave assignment routine UpdateMasterSlaveSets()
  //    is called and the new slave nodes' data containers are initialized.
  //
  // 2) Separate Update and Search
  // -> In this case we have to call SearchContactSeparate(), which
  //    does both bottom-up update (on whole interface) and search. Then
  //    the dynamics master/slave assignment routine UpdateMasterSlaveSets()
  //    is called and the new slave nodes' data containers are initialized.
  //
  // *********************************************************************
  if (SelfContact())
  {
    // evaluate search itself
    binarytreeself_->EvaluateSearch();

    // update master/slave sets of interface
    UpdateMasterSlaveSets();

    // initialize node data container
    // (include slave side boundary nodes / crosspoints)
    for (int i=0; i<SlaveColNodesBound()->NumMyElements(); ++i)
    {
      int gid = SlaveColNodesBound()->GID(i);
      DRT::Node* node = Discret().gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %i",gid);
      MORTAR::MortarNode* mnode = dynamic_cast<MORTAR::MortarNode*>(node);

      // initialize container if not yet initialized before
      mnode->InitializeDataContainer();
    }

    // no initialization of element data container as this would
    // possibly destroy the information on search elements again
    // (this was already done in SetElementAreas())
  }

  // *********************************************************************
  // Possible versions for 2-body contact:
  // *********************************************************************
  //
  // 1) Combined Update and Contact Search
  // -> In this case we only have to call SearchContactCombined(), which
  //    does both top-down update (where necessary) and search.
  //
  // 2) Separate Update and Contact Search
  // -> In this case we have to explicitly call and updating routine, i.e.
  //    UpdateTreeTopDown() or UpdateTreeBottomUp() before calling the
  //    search routine SearchContactSeparate(). Of course, the bottom-up
  //    update makes more sense here. For very large contact problems,
  //    this version is preferable and thus chosen as default.
  //
  // *********************************************************************
  else
  {
    // call mortar routine
    MORTAR::MortarInterface::EvaluateSearchBinarytree();
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  evaluate coupling type segment-to-line coupl             farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::EvaluateSTL()
{
  // check
  if(Dim() == 2)
    dserror("ERROR: LTS algorithm only for 3D simulations!");

  // counter
  int count = 0;

  // loop over slave elements
  for(int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1)
      dserror("ERROR: Cannot find slave element with gid %", gid1);
    CoElement* selement = dynamic_cast<CoElement*>(ele1);

    // guarantee uniquness
    std::set<std::pair<int,int> > donebefore;

    // loop over found meles
    for (int j = 0; j < selement->MoData().NumSearchElements(); ++j)
    {
      int gid2 = selement->MoData().SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2)
        dserror("ERROR: Cannot find master element with gid %", gid2);
      CoElement* melement = dynamic_cast<CoElement*>(ele2);

      if(melement->Shape() == DRT::Element::quad4)
      {
        for(int j = 0; j< 4 ; ++j)
        {
          int nodeIds[2]  = {0,0};
          int nodeLIds[2] = {0,0};

          if(j == 0)
          {
            nodeIds[0] = melement->NodeIds()[0];
            nodeIds[1] = melement->NodeIds()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if(j == 1)
          {
            nodeIds[0] = melement->NodeIds()[1];
            nodeIds[1] = melement->NodeIds()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if(j == 2)
          {
            nodeIds[0] = melement->NodeIds()[2];
            nodeIds[1] = melement->NodeIds()[3];

            nodeLIds[0] = 2;
            nodeLIds[1] = 3;
          }
          else if(j == 3)
          {
            nodeIds[0] = melement->NodeIds()[3];
            nodeIds[1] = melement->NodeIds()[0];

            nodeLIds[0] = 3;
            nodeLIds[1] = 0;
          }

          //create pair
          std::pair<int,int> actIDs   = std::pair<int,int>(nodeIds[0],nodeIds[1]);
          std::pair<int,int> actIDstw = std::pair<int,int>(nodeIds[1],nodeIds[0]);

          // check if processed before
          std::set<std::pair<int,int> >::iterator iter  = donebefore.find(actIDs);
          std::set<std::pair<int,int> >::iterator itertw= donebefore.find(actIDstw);

           // if not then create ele
           if (iter == donebefore.end() and itertw == donebefore.end())
           {
             // add to set of processed nodes
             donebefore.insert(actIDs);
             donebefore.insert(actIDstw);

             // create line ele:
             Teuchos::RCP<MORTAR::MortarElement> lineEle = Teuchos::rcp(
                         new MORTAR::MortarElement(
                             j,
                             melement->Owner(),
                             DRT::Element::line2,
                             2,
                             nodeIds,
                             false));

             // get nodes
             DRT::Node* nodes[2] = {melement->Nodes()[nodeLIds[0]],melement->Nodes()[nodeLIds[1]]};
             lineEle->BuildNodalPointers(nodes);

             // init data container for dual shapes
             lineEle->InitializeDataContainer();

             std::vector<CoElement* > seleElements;
             seleElements.push_back(selement);

             // create coupling object
             LineToSurfaceCoupling3d coup(
                 *idiscret_,
                 3,
                 IParams(),
                 *melement,
                 lineEle,
                 seleElements,
                 LineToSurfaceCoupling3d::stl);

             // perform evaluate!
             coup.EvaluateCoupling();

             // count for safety
             count++;
           }
        } // end edge loop
      }
      else
        dserror("ERROR: LTS only for quad4!");

    }// end found mele loop
  }// end slave ele loop

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate coupling type ndoe-to-segment coupl             farah 11/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::EvaluateNTSMaster()
{
  // create one interpolator instance which is valid for all nodes!
  Teuchos::RCP<NTS::CoInterpolator> interpolator =
      Teuchos::rcp(new NTS::CoInterpolator(IParams(),Dim()));

  // guarantee uniquness
  std::set<int> donebefore;

  // loop over all slave col elements
  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1)
      dserror("ERROR: Cannot find slave element with gid %", gid1);
    CoElement* selement = dynamic_cast<CoElement*>(ele1);

    // skip zero-sized nurbs elements (slave)
    if (selement->ZeroSized())
      continue;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j = 0; j < selement->MoData().NumSearchElements(); ++j)
    {
      int gid2 = selement->MoData().SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2)
        dserror("ERROR: Cannot find master element with gid %", gid2);
      CoElement* melement = dynamic_cast<CoElement*>(ele2);

      // skip zero-sized nurbs elements (master)
      if (melement->ZeroSized())
        continue;

      for(int n=0;n<(int)melement->NumNode();++n)
      {
        MORTAR::MortarNode* cnode = dynamic_cast<MORTAR::MortarNode*>(melement->Nodes()[n]);

        std::set<int>::iterator iter   = donebefore.find(cnode->Id());

         // if not evaluated before
         if (iter == donebefore.end())
         {
           std::vector<MORTAR::MortarElement*> dummy;
           MORTAR::MortarElement* sele = dynamic_cast<MORTAR::MortarElement*>(selement);
           dummy.push_back(sele);

           // call interpolation functions
           bool success = interpolator->Interpolate(*cnode,dummy);
           if(success)
             donebefore.insert(cnode->Id());
         }

      } // node loop
    }// found contact eles
  } // sele loop

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate coupling type line-to-segment coupl             farah 11/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::EvaluateLTSMaster()
{
  // check
  if(Dim() == 2)
    dserror("ERROR: LTS algorithm only for 3D simulations!");

  // counter
  int count = 0;

  // guarantee uniquness
  std::set<std::pair<int,int> > donebefore;

  // loop over slave elements
  for(int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1)
      dserror("ERROR: Cannot find slave element with gid %", gid1);
    CoElement* selement = dynamic_cast<CoElement*>(ele1);

    // ele check
    if(selement->Shape() != DRT::Element::quad4 and
       selement->Shape() != DRT::Element::tri3)
      dserror("ERROR: LTS algorithm only for tri3/quad4!");

    // empty vector of master element pointers
    std::vector<Teuchos::RCP<MORTAR::MortarElement> >  lineElements;
    std::vector<CoElement* > meleElements;

    // compute slave normal
    double slaveN[3]    = {0.0, 0.0, 0.0};
    double loccenter[2] = {0.0, 0.0 };

    DRT::Element::DiscretizationType dt = selement->Shape();
    if (dt == DRT::Element::tri3 || dt == DRT::Element::tri6)
    {
      loccenter[0] = 1.0 / 3.0;
      loccenter[1] = 1.0 / 3.0;
    }
    else if (dt == DRT::Element::quad4 || dt == DRT::Element::quad8
          || dt == DRT::Element::quad9)
    {
      loccenter[0] = 0.0;
      loccenter[1] = 0.0;
    }
    else
      dserror("ERROR: AuxiliaryPlane called for unknown element type");

    // we then compute the unit normal vector at the element center
    selement->ComputeUnitNormalAtXi(loccenter, slaveN);

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j = 0; j < selement->MoData().NumSearchElements(); ++j)
    {
      int gid2 = selement->MoData().SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2)
        dserror("ERROR: Cannot find master element with gid %", gid2);
      CoElement* melement = dynamic_cast<CoElement*>(ele2);

      // check orientation
      // compute slave normal
      double masterN[3] = {0.0, 0.0, 0.0};
      double loccenterM[2] = { 0.0, 0.0 };

      DRT::Element::DiscretizationType dt = melement->Shape();
      if (dt == DRT::Element::tri3 || dt == DRT::Element::tri6)
      {
        loccenterM[0] = 1.0 / 3.0;
        loccenterM[1] = 1.0 / 3.0;
      }
      else if (dt == DRT::Element::quad4 || dt == DRT::Element::quad8
            || dt == DRT::Element::quad9)
      {
        loccenterM[0] = 0.0;
        loccenterM[1] = 0.0;
      }
      else
        dserror("ERROR: AuxiliaryPlane called for unknown element type");

      // we then compute the unit normal vector at the element center
      melement->ComputeUnitNormalAtXi(loccenterM, masterN);

      double scaprod = slaveN[0] * masterN[0] + slaveN[1] * masterN[1] + slaveN[2] * masterN[2];

      // tolerance for line clipping
      const double sminedge = selement->MinEdgeSize();
      const double mminedge = melement->MinEdgeSize();
      const double tol = 0.001 * std::min(sminedge, mminedge);
      if(abs(scaprod)<tol)
        continue;

      // if orientation is okay
      meleElements.push_back(melement);
    }

    // no valid maste elements?
    if(meleElements.size()<1)
      continue;

    for(int m = 0; m<(int)meleElements.size();++m)
    {
      // loop over slave edges -> match node number for tri3/quad4
      for(int j = 0; j< meleElements[m]->NumNode() ; ++j)
      {
        int nodeIds[2]  = {0,0};
        int nodeLIds[2] = {0,0};

        if(meleElements[m]->Shape() == DRT::Element::quad4)
        {
          if(j == 0)
          {
            nodeIds[0] = meleElements[m]->NodeIds()[0];
            nodeIds[1] = meleElements[m]->NodeIds()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if(j == 1)
          {
            nodeIds[0] = meleElements[m]->NodeIds()[1];
            nodeIds[1] = meleElements[m]->NodeIds()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if(j == 2)
          {
            nodeIds[0] = meleElements[m]->NodeIds()[2];
            nodeIds[1] = meleElements[m]->NodeIds()[3];

            nodeLIds[0] = 2;
            nodeLIds[1] = 3;
          }
          else if(j == 3)
          {
            nodeIds[0] = meleElements[m]->NodeIds()[3];
            nodeIds[1] = meleElements[m]->NodeIds()[0];

            nodeLIds[0] = 3;
            nodeLIds[1] = 0;
          }
          else
            dserror("ERROR: loop counter and edge number do not match!");
        }
        else if(meleElements[m]->Shape() == DRT::Element::tri3)
        {
          if(j == 0)
          {
            nodeIds[0] = meleElements[m]->NodeIds()[0];
            nodeIds[1] = meleElements[m]->NodeIds()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if(j == 1)
          {
            nodeIds[0] = meleElements[m]->NodeIds()[1];
            nodeIds[1] = meleElements[m]->NodeIds()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if(j == 2)
          {
            nodeIds[0] = meleElements[m]->NodeIds()[2];
            nodeIds[1] = meleElements[m]->NodeIds()[0];

            nodeLIds[0] = 2;
            nodeLIds[1] = 0;
          }
          else
            dserror("ERROR: loop counter and edge number do not match!");
        }

        // check if both nodes on edge geometry
        bool node0Edge =
            dynamic_cast<MORTAR::MortarNode*>(meleElements[m]->Nodes()[nodeLIds[0]])->IsOnEdge();
        bool node1Edge =
            dynamic_cast<MORTAR::MortarNode*>(meleElements[m]->Nodes()[nodeLIds[1]])->IsOnEdge();

        if(nonSmoothContact_ and (!node0Edge or !node1Edge))
          continue;

        //create pair
        std::pair<int,int> actIDs   = std::pair<int,int>(nodeIds[0],nodeIds[1]);
        std::pair<int,int> actIDstw = std::pair<int,int>(nodeIds[1],nodeIds[0]);

        // check if processed before
        std::set<std::pair<int,int> >::iterator iter   = donebefore.find(actIDs);
        std::set<std::pair<int,int> >::iterator itertw = donebefore.find(actIDstw);

         // if not then create ele
         if (iter == donebefore.end() and itertw == donebefore.end() )
         {
           // add to set of processed nodes
           donebefore.insert(actIDs);
           donebefore.insert(actIDstw);

           // create line ele:
           Teuchos::RCP<MORTAR::MortarElement> lineEle = Teuchos::rcp(
                       new MORTAR::MortarElement(
                           j,
                           meleElements[m]->Owner(),
                           DRT::Element::line2,
                           2,
                           nodeIds,
                           false));

           // get nodes
           DRT::Node* nodes[2] = {meleElements[m]->Nodes()[nodeLIds[0]],meleElements[m]->Nodes()[nodeLIds[1]]};
           lineEle->BuildNodalPointers(nodes);

           // init data container for dual shapes
           lineEle->InitializeDataContainer();

           // push back into vector
           lineElements.push_back(lineEle);

           std::vector<CoElement*> dummy;
           dummy.push_back(selement);

           // create coupling object
           LineToSurfaceCoupling3d coup(
               *idiscret_,
               3,
               IParams(),
               *meleElements[m],
               lineEle,
               dummy,
               LineToSurfaceCoupling3d::lts);


           // perform evaluate!
           coup.EvaluateCoupling();

           // count for safety
           count++;
         }
      } // end edge loop
    }


    // loop over all created line elements
//    for(int l = 0; l<(int)lineElements.size(); ++l)
//    {
//      // create coupling object
//      LineToSurfaceCoupling3d coup(
//          *idiscret_,
//          3,
//          IParams(),
//          *selement,
//          lineElements[l],
//          meleElements,
//          LineToSurfaceCoupling3d::lts);
//
//      // perform evaluate!
//      coup.EvaluateCoupling();
//    }
  }//slave ele loop

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate coupling type line-to-segment coupl             farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::EvaluateLTS()
{
  // check
  if(Dim() == 2)
    dserror("ERROR: LTS algorithm only for 3D simulations!");

  // counter
  int count = 0;

  // guarantee uniquness
  std::set<std::pair<int,int> > donebefore;

  // loop over slave elements
  for(int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1)
      dserror("ERROR: Cannot find slave element with gid %", gid1);
    CoElement* selement = dynamic_cast<CoElement*>(ele1);

    // ele check
    if(selement->Shape() != DRT::Element::quad4 and
       selement->Shape() != DRT::Element::tri3)
      dserror("ERROR: LTS algorithm only for tri3/quad4!");

    // empty vector of master element pointers
    std::vector<Teuchos::RCP<MORTAR::MortarElement> >  lineElements;
    std::vector<CoElement* > meleElements;

    // compute slave normal
    double slaveN[3]    = {0.0, 0.0, 0.0};
    double loccenter[2] = {0.0, 0.0 };

    DRT::Element::DiscretizationType dt = selement->Shape();
    if (dt == DRT::Element::tri3 || dt == DRT::Element::tri6)
    {
      loccenter[0] = 1.0 / 3.0;
      loccenter[1] = 1.0 / 3.0;
    }
    else if (dt == DRT::Element::quad4 || dt == DRT::Element::quad8
          || dt == DRT::Element::quad9)
    {
      loccenter[0] = 0.0;
      loccenter[1] = 0.0;
    }
    else
      dserror("ERROR: AuxiliaryPlane called for unknown element type");

    // we then compute the unit normal vector at the element center
    selement->ComputeUnitNormalAtXi(loccenter, slaveN);

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j = 0; j < selement->MoData().NumSearchElements(); ++j)
    {
      int gid2 = selement->MoData().SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2)
        dserror("ERROR: Cannot find master element with gid %", gid2);
      CoElement* melement = dynamic_cast<CoElement*>(ele2);

      // check orientation
      // compute slave normal
      double masterN[3] = {0.0, 0.0, 0.0};
      double loccenterM[2] = { 0.0, 0.0 };

      DRT::Element::DiscretizationType dt = melement->Shape();
      if (dt == DRT::Element::tri3 || dt == DRT::Element::tri6)
      {
        loccenterM[0] = 1.0 / 3.0;
        loccenterM[1] = 1.0 / 3.0;
      }
      else if (dt == DRT::Element::quad4 || dt == DRT::Element::quad8
            || dt == DRT::Element::quad9)
      {
        loccenterM[0] = 0.0;
        loccenterM[1] = 0.0;
      }
      else
        dserror("ERROR: AuxiliaryPlane called for unknown element type");

      // we then compute the unit normal vector at the element center
      melement->ComputeUnitNormalAtXi(loccenterM, masterN);

      double scaprod = slaveN[0] * masterN[0] + slaveN[1] * masterN[1] + slaveN[2] * masterN[2];

      // tolerance for line clipping
      const double sminedge = selement->MinEdgeSize();
      const double mminedge = melement->MinEdgeSize();
      const double tol = 0.001 * std::min(sminedge, mminedge);
      if(abs(scaprod)<tol)
        continue;

      // if orientation is okay
      meleElements.push_back(melement);
    }

    // no valid maste elements?
    if(meleElements.size()<1)
      continue;

    // loop over slave edges -> match node number for tri3/quad4
    for(int j = 0; j< selement->NumNode() ; ++j)
    {
      int nodeIds[2]  = {0,0};
      int nodeLIds[2] = {0,0};

      if(selement->Shape() == DRT::Element::quad4)
      {
        if(j == 0)
        {
          nodeIds[0] = selement->NodeIds()[0];
          nodeIds[1] = selement->NodeIds()[1];

          nodeLIds[0] = 0;
          nodeLIds[1] = 1;
        }
        else if(j == 1)
        {
          nodeIds[0] = selement->NodeIds()[1];
          nodeIds[1] = selement->NodeIds()[2];

          nodeLIds[0] = 1;
          nodeLIds[1] = 2;
        }
        else if(j == 2)
        {
          nodeIds[0] = selement->NodeIds()[2];
          nodeIds[1] = selement->NodeIds()[3];

          nodeLIds[0] = 2;
          nodeLIds[1] = 3;
        }
        else if(j == 3)
        {
          nodeIds[0] = selement->NodeIds()[3];
          nodeIds[1] = selement->NodeIds()[0];

          nodeLIds[0] = 3;
          nodeLIds[1] = 0;
        }
        else
          dserror("ERROR: loop counter and edge number do not match!");
      }
      else if(selement->Shape() == DRT::Element::tri3)
      {
        if(j == 0)
        {
          nodeIds[0] = selement->NodeIds()[0];
          nodeIds[1] = selement->NodeIds()[1];

          nodeLIds[0] = 0;
          nodeLIds[1] = 1;
        }
        else if(j == 1)
        {
          nodeIds[0] = selement->NodeIds()[1];
          nodeIds[1] = selement->NodeIds()[2];

          nodeLIds[0] = 1;
          nodeLIds[1] = 2;
        }
        else if(j == 2)
        {
          nodeIds[0] = selement->NodeIds()[2];
          nodeIds[1] = selement->NodeIds()[0];

          nodeLIds[0] = 2;
          nodeLIds[1] = 0;
        }
        else
          dserror("ERROR: loop counter and edge number do not match!");
      }

      // check if both nodes on edge geometry
      bool node0Edge =
          dynamic_cast<MORTAR::MortarNode*>(selement->Nodes()[nodeLIds[0]])->IsOnEdge();
      bool node1Edge =
          dynamic_cast<MORTAR::MortarNode*>(selement->Nodes()[nodeLIds[1]])->IsOnEdge();

      if(nonSmoothContact_ and (!node0Edge or !node1Edge))
        continue;

      //create pair
      std::pair<int,int> actIDs   = std::pair<int,int>(nodeIds[0],nodeIds[1]);
      std::pair<int,int> actIDstw = std::pair<int,int>(nodeIds[1],nodeIds[0]);

      // check if processed before
      std::set<std::pair<int,int> >::iterator iter   = donebefore.find(actIDs);
      std::set<std::pair<int,int> >::iterator itertw = donebefore.find(actIDstw);

       // if not then create ele
       if (iter == donebefore.end() and itertw == donebefore.end() )
       {
         // add to set of processed nodes
         donebefore.insert(actIDs);
         donebefore.insert(actIDstw);

         // create line ele:
         Teuchos::RCP<MORTAR::MortarElement> lineEle = Teuchos::rcp(
                     new MORTAR::MortarElement(
                         j,
                         selement->Owner(),
                         DRT::Element::line2,
                         2,
                         nodeIds,
                         false));

         // get nodes
         DRT::Node* nodes[2] = {selement->Nodes()[nodeLIds[0]],selement->Nodes()[nodeLIds[1]]};
         lineEle->BuildNodalPointers(nodes);

         // init data container for dual shapes
         lineEle->InitializeDataContainer();

         // push back into vector
         lineElements.push_back(lineEle);

         // count for safety
         count++;
       }
    } // end edge loop

    // loop over all created line elements
    for(int l = 0; l<(int)lineElements.size(); ++l)
    {
      // create coupling object
      LineToSurfaceCoupling3d coup(
          *idiscret_,
          3,
          IParams(),
          *selement,
          lineElements[l],
          meleElements,
          LineToSurfaceCoupling3d::lts);

      // perform evaluate!
      coup.EvaluateCoupling();
    }
  }//slave ele loop

  // bye bye
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate coupling type line-to-line coupl                farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::EvaluateLTL()
{
  // check
  if(Dim() == 2)
    dserror("ERROR: LTL algorithm only for 3D simulations!");

  // guarantee uniquness of slave edges
  std::set<std::pair<int,int> > donebeforeS;

  // loop over slave elements
  for(int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1)
      dserror("ERROR: Cannot find slave element with gid %", gid1);
    CoElement* selement = dynamic_cast<CoElement*>(ele1);

    // empty vector of slave element pointers
    std::vector<Teuchos::RCP<MORTAR::MortarElement> > lineElementsS;

    if(selement->Shape() == DRT::Element::quad4)
    {
      for(int j = 0; j< 4 ; ++j)
      {
        int nodeIds[2]  = {0,0};
        int nodeLIds[2] = {0,0};

        if(j == 0)
        {
          nodeIds[0] = selement->NodeIds()[0];
          nodeIds[1] = selement->NodeIds()[1];

          nodeLIds[0] = 0;
          nodeLIds[1] = 1;
        }
        else if(j == 1)
        {
          nodeIds[0] = selement->NodeIds()[1];
          nodeIds[1] = selement->NodeIds()[2];

          nodeLIds[0] = 1;
          nodeLIds[1] = 2;
        }
        else if(j == 2)
        {
          nodeIds[0] = selement->NodeIds()[2];
          nodeIds[1] = selement->NodeIds()[3];

          nodeLIds[0] = 2;
          nodeLIds[1] = 3;
        }
        else if(j == 3)
        {
          nodeIds[0] = selement->NodeIds()[3];
          nodeIds[1] = selement->NodeIds()[0];

          nodeLIds[0] = 3;
          nodeLIds[1] = 0;
        }

        // check if both nodes on edge geometry
        bool node0Edge = dynamic_cast<MORTAR::MortarNode*>(selement->Nodes()[nodeLIds[0]])->IsOnEdge();
        bool node1Edge = dynamic_cast<MORTAR::MortarNode*>(selement->Nodes()[nodeLIds[1]])->IsOnEdge();

        if(!node0Edge or !node1Edge)
          continue;

        //create pair
        std::pair<int,int> actIDs   = std::pair<int,int>(nodeIds[0],nodeIds[1]);
        std::pair<int,int> actIDstw = std::pair<int,int>(nodeIds[1],nodeIds[0]);

        // check if processed before
        std::set<std::pair<int,int> >::iterator iter   = donebeforeS.find(actIDs);
        std::set<std::pair<int,int> >::iterator itertw = donebeforeS.find(actIDstw);

         // if not then create ele
         if (iter == donebeforeS.end() and itertw == donebeforeS.end())
         {
           // add to set of processed nodes
           donebeforeS.insert(actIDs);
           donebeforeS.insert(actIDstw);

           // create line ele:
           Teuchos::RCP<MORTAR::MortarElement> lineEle = Teuchos::rcp(
                       new MORTAR::MortarElement(
                           j,
                           selement->Owner(),
                           DRT::Element::line2,
                           2,
                           nodeIds,
                           false));

           // get nodes
           DRT::Node* nodes[2] = {selement->Nodes()[nodeLIds[0]],selement->Nodes()[nodeLIds[1]]};
           lineEle->BuildNodalPointers(nodes);

           // init data container for dual shapes
           lineEle->InitializeDataContainer();

           // push back into vector
           lineElementsS.push_back(lineEle);
         }
      } // end edge loop
    }
    else if(selement->Shape() == DRT::Element::tri3)
    {
      for(int j = 0; j< 3 ; ++j)
      {
        int nodeIds[2]  = {0,0};
        int nodeLIds[2] = {0,0};

        if(j == 0)
        {
          nodeIds[0] = selement->NodeIds()[0];
          nodeIds[1] = selement->NodeIds()[1];

          nodeLIds[0] = 0;
          nodeLIds[1] = 1;
        }
        else if(j == 1)
        {
          nodeIds[0] = selement->NodeIds()[1];
          nodeIds[1] = selement->NodeIds()[2];

          nodeLIds[0] = 1;
          nodeLIds[1] = 2;
        }
        else if(j == 2)
        {
          nodeIds[0] = selement->NodeIds()[2];
          nodeIds[1] = selement->NodeIds()[0];

          nodeLIds[0] = 2;
          nodeLIds[1] = 0;
        }

        // check if both nodes on edge geometry
        bool node0Edge = dynamic_cast<MORTAR::MortarNode*>(selement->Nodes()[nodeLIds[0]])->IsOnEdge();
        bool node1Edge = dynamic_cast<MORTAR::MortarNode*>(selement->Nodes()[nodeLIds[1]])->IsOnEdge();

        if(!node0Edge or !node1Edge)
          continue;

        //create pair
        std::pair<int,int> actIDs   = std::pair<int,int>(nodeIds[0],nodeIds[1]);
        std::pair<int,int> actIDstw = std::pair<int,int>(nodeIds[1],nodeIds[0]);

        // check if processed before
        std::set<std::pair<int,int> >::iterator iter   = donebeforeS.find(actIDs);
        std::set<std::pair<int,int> >::iterator itertw = donebeforeS.find(actIDstw);

         // if not then create ele
         if (iter == donebeforeS.end() and itertw == donebeforeS.end())
         {
           // add to set of processed nodes
           donebeforeS.insert(actIDs);
           donebeforeS.insert(actIDstw);

           // create line ele:
           Teuchos::RCP<MORTAR::MortarElement> lineEle = Teuchos::rcp(
                       new MORTAR::MortarElement(
                           j,
                           selement->Owner(),
                           DRT::Element::line2,
                           2,
                           nodeIds,
                           false));

           // get nodes
           DRT::Node* nodes[2] = {selement->Nodes()[nodeLIds[0]],selement->Nodes()[nodeLIds[1]]};
           lineEle->BuildNodalPointers(nodes);

           // init data container for dual shapes
           lineEle->InitializeDataContainer();

           // push back into vector
           lineElementsS.push_back(lineEle);
         }
      } // end edge loop
    }
    else
      dserror("ERROR: LTL only for quad4 and tri3!");

    // guarantee uniquness of master edges
    std::set<std::pair<int,int> > donebeforeM;

    // empty vector of slave element pointers
    std::vector<Teuchos::RCP<MORTAR::MortarElement> > lineElementsM;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int k = 0; k < selement->MoData().NumSearchElements(); ++k)
    {
      int gid2 = selement->MoData().SearchElements()[k];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2)
        dserror("ERROR: Cannot find master element with gid %", gid2);
      CoElement* melement = dynamic_cast<CoElement*>(ele2);

      if(melement->Shape() == DRT::Element::quad4)
      {
        for(int j = 0; j< 4 ; ++j)
        {
          int nodeIds[2]  = {0,0};
          int nodeLIds[2] = {0,0};

          if(j == 0)
          {
            nodeIds[0] = melement->NodeIds()[0];
            nodeIds[1] = melement->NodeIds()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if(j == 1)
          {
            nodeIds[0] = melement->NodeIds()[1];
            nodeIds[1] = melement->NodeIds()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if(j == 2)
          {
            nodeIds[0] = melement->NodeIds()[2];
            nodeIds[1] = melement->NodeIds()[3];

            nodeLIds[0] = 2;
            nodeLIds[1] = 3;
          }
          else if(j == 3)
          {
            nodeIds[0] = melement->NodeIds()[3];
            nodeIds[1] = melement->NodeIds()[0];

            nodeLIds[0] = 3;
            nodeLIds[1] = 0;
          }

          // check if both nodes on edge geometry
          bool node0Edge = dynamic_cast<MORTAR::MortarNode*>(melement->Nodes()[nodeLIds[0]])->IsOnEdge();
          bool node1Edge = dynamic_cast<MORTAR::MortarNode*>(melement->Nodes()[nodeLIds[1]])->IsOnEdge();

          if(!node0Edge or !node1Edge)
            continue;

          //create pair
          std::pair<int,int> actIDs   = std::pair<int,int>(nodeIds[0],nodeIds[1]);
          std::pair<int,int> actIDstw = std::pair<int,int>(nodeIds[1],nodeIds[0]);

          // check if processed before
          std::set<std::pair<int,int> >::iterator iter   = donebeforeM.find(actIDs);
          std::set<std::pair<int,int> >::iterator itertw = donebeforeM.find(actIDstw);

           // if not then create ele
           if (iter == donebeforeM.end() and itertw == donebeforeM.end())
           {
             // add to set of processed nodes
             donebeforeM.insert(actIDs);
             donebeforeM.insert(actIDstw);

             // create line ele:
             Teuchos::RCP<MORTAR::MortarElement> lineEle = Teuchos::rcp(
                         new MORTAR::MortarElement(
                             j,
                             melement->Owner(),
                             DRT::Element::line2,
                             2,
                             nodeIds,
                             false));

             // get nodes
             DRT::Node* nodes[2] = {melement->Nodes()[nodeLIds[0]],melement->Nodes()[nodeLIds[1]]};
             lineEle->BuildNodalPointers(nodes);

             // init data container for dual shapes
             lineEle->InitializeDataContainer();

             // push back into vector
             lineElementsM.push_back(lineEle);
           }
        } // end edge loop
      }
      else if(melement->Shape() == DRT::Element::tri3)
      {
        for(int j = 0; j< 3 ; ++j)
        {
          int nodeIds[2]  = {0,0};
          int nodeLIds[2] = {0,0};

          if(j == 0)
          {
            nodeIds[0] = melement->NodeIds()[0];
            nodeIds[1] = melement->NodeIds()[1];

            nodeLIds[0] = 0;
            nodeLIds[1] = 1;
          }
          else if(j == 1)
          {
            nodeIds[0] = melement->NodeIds()[1];
            nodeIds[1] = melement->NodeIds()[2];

            nodeLIds[0] = 1;
            nodeLIds[1] = 2;
          }
          else if(j == 2)
          {
            nodeIds[0] = melement->NodeIds()[2];
            nodeIds[1] = melement->NodeIds()[0];

            nodeLIds[0] = 2;
            nodeLIds[1] = 0;
          }

          // check if both nodes on edge geometry
          bool node0Edge = dynamic_cast<MORTAR::MortarNode*>(melement->Nodes()[nodeLIds[0]])->IsOnEdge();
          bool node1Edge = dynamic_cast<MORTAR::MortarNode*>(melement->Nodes()[nodeLIds[1]])->IsOnEdge();

          if(!node0Edge or !node1Edge)
            continue;

          //create pair
          std::pair<int,int> actIDs   = std::pair<int,int>(nodeIds[0],nodeIds[1]);
          std::pair<int,int> actIDstw = std::pair<int,int>(nodeIds[1],nodeIds[0]);

          // check if processed before
          std::set<std::pair<int,int> >::iterator iter   = donebeforeM.find(actIDs);
          std::set<std::pair<int,int> >::iterator itertw = donebeforeM.find(actIDstw);

           // if not then create ele
           if (iter == donebeforeM.end() and itertw == donebeforeM.end())
           {
             // add to set of processed nodes
             donebeforeM.insert(actIDs);
             donebeforeM.insert(actIDstw);

             // create line ele:
             Teuchos::RCP<MORTAR::MortarElement> lineEle = Teuchos::rcp(
                         new MORTAR::MortarElement(
                             j,
                             melement->Owner(),
                             DRT::Element::line2,
                             2,
                             nodeIds,
                             false));

             // get nodes
             DRT::Node* nodes[2] = {melement->Nodes()[nodeLIds[0]],melement->Nodes()[nodeLIds[1]]};
             lineEle->BuildNodalPointers(nodes);

             // init data container for dual shapes
             lineEle->InitializeDataContainer();

             // push back into vector
             lineElementsM.push_back(lineEle);
           }
        } // end edge loop
      }
      else
        dserror("ERROR: LTL only for quad4 and tri3!");
    }// end found mele loop

    // loop over slave edges
    for( int s = 0; s<(int)lineElementsS.size(); ++s)
    {
      // loop over master edges
      for( int m = 0; m<(int)lineElementsM.size(); ++m)
      {
        // create coupling object
        LineToLineCouplingPoint3d coup(
            *idiscret_,
            3,
            IParams(),
            lineElementsS[s],
            lineElementsM[m]);

        // perform evaluate!
        coup.EvaluateCoupling();
      }
    }

  }// end slave loop

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate coupling type node-to-segment coupl             farah 02/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::EvaluateNTS()
{
  // create one interpolator instance which is valid for all nodes!
  Teuchos::RCP<NTS::CoInterpolator> interpolator =
      Teuchos::rcp(new NTS::CoInterpolator(IParams(),Dim()));

  // loop over slave nodes
  for (int i = 0; i < snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find node with gid %", gid);
    MORTAR::MortarNode* mrtrnode = dynamic_cast<MORTAR::MortarNode*>(node);

    if(!mrtrnode->IsOnCorner() and nonSmoothContact_)
      continue;

    if (mrtrnode->Owner() != Comm().MyPID())
      dserror("ERROR: Node ownership inconsistency!");

    // vector with possible contacting master eles
    std::vector<MORTAR::MortarElement*> meles;

    // fill vector with possibly contacting meles
    FindMEles(*mrtrnode,meles);

    // skip calculation if no meles vector is empty
    if(meles.size() < 1)
      continue;

    // call interpolation functions
    interpolator->Interpolate(*mrtrnode,meles);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Integrate matrix M and gap g on slave/master overlaps     popp 11/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::MortarCoupling(
    MORTAR::MortarElement* sele,
    std::vector<MORTAR::MortarElement*> mele,
    const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr)
{
  // do stuff before the actual coupling is going to be evaluated
  PreMortarCoupling(sele,mele,mparams_ptr);

  // increase counter of slave/master pairs
  smpairs_ += (int)mele.size();

  // check if quadratic interpolation is involved
  bool quadratic = false;
  if (sele->IsQuad())
    quadratic = true;
  for (int m=0;m<(int)mele.size();++m)
    if (mele[m]->IsQuad())
      quadratic = true;

  // *********************************************************************
  // do interface coupling within a new class
  // (projection slave and master, overlap detection, integration and
  // linearization of the Mortar matrix M)
  // ************************************************************** 2D ***
  if (Dim()==2)
  {
    // *************************************************** linear 2D ***
    // ************************************************ quadratic 2D ***
    // neither quadratic interpolation nor mixed linear and quadratic
    // interpolation need any special treatment in the 2d case

    // create CoCoupling2dManager
    CONTACT::CoCoupling2dManager coup(Discret(),Dim(),quadratic,IParams(),sele,mele);
    // evaluate
    coup.EvaluateCoupling(mparams_ptr);

    // increase counter of slave/master integration pairs and intcells
    smintpairs_ += (int)mele.size();
    intcells_   += (int)mele.size();
  }
  // ************************************************************** 3D ***
  else if (Dim()==3)
  {
    // *************************************************** linear 3D ***
    if (!quadratic)
    {
      // create CoCoupling3dManager
      CONTACT::CoCoupling3dManager coup(Discret(),Dim(),quadratic,IParams(),sele,mele);
      // evaluate
      coup.EvaluateCoupling(mparams_ptr);

      // increase counter of slave/master integration pairs and intcells
      smintpairs_ += (int)mele.size();
      intcells_   += coup.IntegrationCells();
    }

    // ************************************************** quadratic 3D ***
    else
    {
      //create Coupling3dQuadManager
      CONTACT::CoCoupling3dQuadManager coup(Discret(),Dim(),quadratic,IParams(),sele,mele);
      // evaluate
      coup.EvaluateCoupling(mparams_ptr);
    } // quadratic
  } // 3D
  else
    dserror("ERROR: Dimension for Mortar coupling must be 2D or 3D!");
  // *********************************************************************

  // do stuff after the coupling evaluation
  PostMortarCoupling(sele,mele,mparams_ptr);

  return true;
}

/*----------------------------------------------------------------------*
 |  Integrate penalty scaling factor kapp (public)            popp 11/09|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::IntegrateKappaPenalty(CONTACT::CoElement& sele)
{
  // create correct integration limits
  double sxia[2] = {0.0, 0.0};
  double sxib[2] = {0.0, 0.0};
  if (sele.Shape()==DRT::Element::tri3 || sele.Shape()==DRT::Element::tri6)
  {
    // parameter space is [0,1] for triangles
    sxib[0] = 1.0; sxib[1] = 1.0;
  }
  else
  {
    // parameter space is [-1,1] for quadrilaterals
    sxia[0] = -1.0; sxia[1] = -1.0;
    sxib[0] =  1.0; sxib[1] =  1.0;
  }

  // ************************************************** quadratic 3D ***
  if (Dim()==3 && sele.IsQuad())
  {
    // get LM interpolation and testing type
    INPAR::MORTAR::LagMultQuad lmtype =
      DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(IParams(),"LM_QUAD");

    // build linear integration elements from quadratic CElements
    std::vector<Teuchos::RCP<MORTAR::IntElement> > sauxelements(0);
    SplitIntElements(sele,sauxelements);

    // different options for mortar integration
    if (lmtype == INPAR::MORTAR::lagmult_quad || lmtype == INPAR::MORTAR::lagmult_lin)
    {
      // do the element integration of kappa and store into gap
      int nrow = sele.NumNode();
      Teuchos::RCP<Epetra_SerialDenseVector> gseg = Teuchos::rcp(new Epetra_SerialDenseVector(nrow));

      // create a CONTACT integrator instance with correct NumGP and Dim
      CONTACT::CoIntegrator integrator(imortar_,sele.Shape(),Comm());
      integrator.IntegrateKappaPenalty(sele,sxia,sxib,gseg);

      // do the assembly into the slave nodes
      integrator.AssembleG(Comm(),sele,*gseg);
    }

    else if (lmtype == INPAR::MORTAR::lagmult_pwlin)
    {
      // integrate each int element seperately
      for (int i=0;i<(int)sauxelements.size();++i)
      {
        // do the int element integration of kappa and store into gap
        int nrow = sauxelements[i]->NumNode();
        Teuchos::RCP<Epetra_SerialDenseVector> gseg = Teuchos::rcp(new Epetra_SerialDenseVector(nrow));

        // create a CONTACT integrator instance with correct NumGP and Dim
        CONTACT::CoIntegrator integrator(imortar_,sauxelements[i]->Shape(),Comm());
        integrator.IntegrateKappaPenalty(sele,*(sauxelements[i]),sxia,sxib,gseg);

        // do the assembly into the slave nodes
        integrator.AssembleG(Comm(),*(sauxelements[i]),*gseg);
      }
    }

    else
    {
      dserror("ERROR: IntegrateKappaPenalty: Invalid case for 3D mortar contact LM interpolation");
    }
  }

  // *************************************************** other cases ***
  else
  {
    // do the element integration of kappa and store into gap
    int nrow = sele.NumNode();
    Teuchos::RCP<Epetra_SerialDenseVector> gseg = Teuchos::rcp(new Epetra_SerialDenseVector(nrow));

    // create a CONTACT integrator instance with correct NumGP and Dim
    CONTACT::CoIntegrator integrator(imortar_,sele.Shape(),Comm());
    integrator.IntegrateKappaPenalty(sele,sxia,sxib,gseg);

    // do the assembly into the slave nodes
    integrator.AssembleG(Comm(),sele,*gseg);
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate relative movement (jump) of a slave node     gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::EvaluateRelMov(const Teuchos::RCP<Epetra_Vector> xsmod,
                                          const Teuchos::RCP<LINALG::SparseMatrix> dmatrixmod,
                                          const Teuchos::RCP<LINALG::SparseMatrix> doldmod)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  if (friction_ == false)
    dserror("Error in CoInterface::EvaluateRelMov(): Only evaluated for frictional contact");

  // parameters
  double pp = IParams().get<double>("PENALTYPARAM");

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);
    double cn = GetCnRef()[GetCnRef().Map().LID(cnode->Id())];

    // get some informatiom form the node
    double gap = cnode->CoData().Getg();

    int dim = cnode->NumDof();

    // compute normal part of Lagrange multiplier
    double nz = 0.0;
    for (int k=0;k<3;++k)
      nz += cnode->MoData().n()[k] * cnode->MoData().lm()[k];

    std::vector <double> jump(dim);
    for(int dim=0;dim<Dim();dim++)
      jump[dim] = 0;

    double lmuzawan = 0.0;
    for (int k=0;k<dim;++k)
      lmuzawan += cnode->MoData().lmuzawa()[k]*cnode->MoData().n()[k];

    double kappa = cnode->CoData().Kappa();

    // evaluate jump (relative displacement) of this node
    // only when the node is going to be active, otherwise,
    // this value isn't needed.
    bool activeinfuture = false;

    if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(IParams(),"STRATEGY")== INPAR::CONTACT::solution_penalty)
    {
      if (-gap >= 0) activeinfuture = true;
    }
    else if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(IParams(),"STRATEGY")== INPAR::CONTACT::solution_lagmult and
             DRT::INPUT::IntegralValue<int>(IParams(),"SEMI_SMOOTH_NEWTON")!=1)
    {
      if (-gap >= 0) activeinfuture = true;
    }
    else if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(IParams(),"STRATEGY")== INPAR::CONTACT::solution_lagmult and
             DRT::INPUT::IntegralValue<int>(IParams(),"SEMI_SMOOTH_NEWTON")==1)
    {
      if((nz - cn*gap > 0) or cnode->Active()) activeinfuture = true;
    }
    else if (DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(IParams(),"STRATEGY")== INPAR::CONTACT::solution_uzawa)
    {
      if(lmuzawan - kappa * pp * gap >= 0) activeinfuture = true;
    }
    else
      dserror("Error in Interface::EvaluateRelMov(): Solution strategy not known!");

    if(activeinfuture==true)
    {
      GEN::pairedvector<int,double>& dmap = cnode->MoData().GetD();
      GEN::pairedvector<int,double>& dmapold = cnode->FriData().GetDOld();

      std::set <int> snodes = cnode->FriData().GetSNodes();

      // check if there are entries in the old D map
      if(dmapold.size()< 1)
        dserror("Error in Interface::EvaluateRelMov(): No old D-Map!");

      std::map<int,double>::iterator colcurr;
      std::set<int>::iterator scurr;

      // loop over all slave nodes with an entry adjacent to this node
      for (scurr=snodes.begin(); scurr != snodes.end(); scurr++)
      {
        int gid = *scurr;
        DRT::Node* snode = idiscret_->gNode(gid);
        if (!snode) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* csnode = dynamic_cast<CoNode*>(snode);

        double dik = dmap[csnode->Id()];
        double dikold = dmapold[csnode->Id()];

        std::map<int,double>::iterator mcurr;

        for (int dim=0;dim<csnode->NumDof();++dim)
        {
          int locid = (xsmod->Map()).LID(csnode->Dofs()[dim]);
          jump[dim]-=(dik-dikold)*(*xsmod)[locid];
        }
      } //  loop over adjacent slave nodes

      std::map<int,double>& mmap = cnode->MoData().GetM();
      std::map<int,double>& mmapold = cnode->FriData().GetMOld();

      const std::set <int>& mnodescurrent = cnode->FriData().GetMNodes();
      const std::set <int>& mnodesold = cnode->FriData().GetMNodesOld();

      // check if there are entries in the M map
      if(mmap.size()< 1)
        dserror("Error in Interface::EvaluateRelMov(): No M-Map!");

      // check if there are entries in the old M map
      if(mmapold.size()< 1)
        dserror("Error in Interface::EvaluateRelMov(): No old M-Map!");

      if(mnodesold.size() <1)
        dserror ("Error in Interface::EvaluateRelMov(): No old M-Set!");

      std::set <int> mnodes;
      std::set<int>::iterator mcurr;

      for (mcurr=mnodescurrent.begin(); mcurr != mnodescurrent.end(); mcurr++)
        mnodes.insert(*mcurr);

      for (mcurr=mnodesold.begin(); mcurr != mnodesold.end(); mcurr++)
        mnodes.insert(*mcurr);

      // loop over all master nodes (find adjacent ones to this slip node)
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* cmnode = dynamic_cast<CoNode*>(mnode);

        double mik = mmap[cmnode->Id()];
        double mikold = mmapold[cmnode->Id()];

        std::map<int,double>::iterator mcurr;

        for (int dim=0;dim<cnode->NumDof();++dim)
        {
           jump[dim]+= (mik-mikold)*(cmnode->xspatial()[dim]);
        }
      } //  loop over master nodes

      // write it to nodes
      for(int dim=0;dim<Dim();dim++)
        cnode->FriData().jump()[dim] = jump[dim];

      // linearization of jump vector

      // reset derivative map of jump
      for (int j=0; j<(int)((cnode->FriData().GetDerivJump()).size()); ++j)
        (cnode->FriData().GetDerivJump())[j].clear();
      (cnode->FriData().GetDerivJump()).resize(0);

      /*** 01  **********************************************************/

      if(dmatrixmod==Teuchos::null)
      {
        // loop over according slave nodes
        for (scurr=snodes.begin(); scurr != snodes.end(); scurr++)
        {
          int gid = *scurr;
          DRT::Node* snode = idiscret_->gNode(gid);
          if (!snode) dserror("ERROR: Cannot find node with gid %",gid);
          CoNode* csnode = dynamic_cast<CoNode*>(snode);

          double dik = dmap[csnode->Id()];
          double dikold=dmapold[csnode->Id()];

          for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
          {
            int col = csnode->Dofs()[dimrow];
            double val = -(dik-dikold);
            if (abs(val)>1e-14)
              cnode->AddDerivJumpValue(dimrow,col,val);
          }
        }
      }
      // in the 3D quadratic case, the values are obtained from the
      // global matrices Dmod and Doldmod
      else
      {
        // loop over dimension of the node
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int NumEntries = 0;
          int NumEntriesOld = 0;
          std::vector<double> Values((dmatrixmod->EpetraMatrix())->MaxNumEntries());
          std::vector<int> Indices((dmatrixmod->EpetraMatrix())->MaxNumEntries());
          std::vector<double> ValuesOld((dmatrixmod->EpetraMatrix())->MaxNumEntries());
          std::vector<int> IndicesOld((dmatrixmod->EpetraMatrix())->MaxNumEntries());

          // row
          int row = cnode->Dofs()[dim];

          // extract entries of this row from matrix
          int err = (dmatrixmod->EpetraMatrix())->ExtractGlobalRowCopy(row,(dmatrixmod->EpetraMatrix())->MaxNumEntries(),NumEntries,&Values[0], &Indices[0]);
          if (err) dserror("ExtractMyRowView failed: err=%d", err);

          int errold = (doldmod->EpetraMatrix())->ExtractGlobalRowCopy(row,(doldmod->EpetraMatrix())->MaxNumEntries(),NumEntriesOld,&ValuesOld[0], &IndicesOld[0]);
          if (errold) dserror("ExtractMyRowView failed: err=%d", err);

          // loop over entries of this vector
          for (int j=0;j<NumEntries;++j)
          {
            double ValueOld=0;
            bool found = false;

            // find value with the same index in vector of Dold
            for (int k=0;k<NumEntriesOld;++k)
            {
              if (Indices[k]==Indices[j])
              {
                ValueOld = ValuesOld[k];
                found = true;
                break;
              }
            }

            if(found==false or abs(ValueOld) < 1e-12)
              dserror("Error in EvaluareRelMov(): No old D value exists");

            // write to node
            cnode->AddDerivJumpValue(dim,Indices[j],(Values[j]-ValueOld));
          }
        }
      }

      /*** 02  **********************************************************/
      // loop over according master nodes
      for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
      {
        int gid = *mcurr;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* cmnode = dynamic_cast<CoNode*>(mnode);

        double mik = mmap[cmnode->Id()];
        double mikold=mmapold[cmnode->Id()];

        for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
        {
            int col = cmnode->Dofs()[dimrow];
            double val = (mik-mikold);
            if (abs(val)>1e-14)
              cnode->AddDerivJumpValue(dimrow,col,val);
        }
      }

      /*** 03 ***********************************************************/
      // we need the Lin(D-matrix) entries of this node
      std::map<int,std::map<int,double> >& ddmap = cnode->CoData().GetDerivD();
      std::map<int,std::map<int,double> >::iterator dscurr;

      // loop over all slave nodes in the DerivM-map of the stick slave node
      for (dscurr=ddmap.begin();dscurr!=ddmap.end();++dscurr)
      {
        int gid = dscurr->first;
        DRT::Node* snode = idiscret_->gNode(gid);
        if (!snode) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* csnode = dynamic_cast<CoNode*>(snode);

        // compute entry of the current stick node / slave node pair
        std::map<int,double>& thisdmmap = cnode->CoData().GetDerivD(gid);

        // loop over all entries of the current derivative map
        for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
        {
          int col = colcurr->first;

          // loop over dimensions
          for(int dim=0;dim<cnode->NumDof();++dim)
          {
            int locid = (xsmod->Map()).LID(csnode->Dofs()[dim]);
            double val =-colcurr->second*(*xsmod)[locid];
            if (abs(val)>1e-14)
              cnode->AddDerivJumpValue(dim,col,val);
          }
        }
      }

      /*** 04 ***********************************************************/
      // we need the Lin(M-matrix) entries of this node
      std::map<int,std::map<int,double> >& dmmap = cnode->CoData().GetDerivM();
      std::map<int,std::map<int,double> >::iterator dmcurr;

      // loop over all master nodes in the DerivM-map of the stick slave node
      for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
      {
        int gid = dmcurr->first;
        DRT::Node* mnode = idiscret_->gNode(gid);
        if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
        CoNode* cmnode = dynamic_cast<CoNode*>(mnode);
        double* mxi = cmnode->xspatial();

        // compute entry of the current stick node / master node pair
        std::map<int,double>& thisdmmap = cnode->CoData().GetDerivM(gid);

        // loop over all entries of the current derivative map
        for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
        {
          int col = colcurr->first;

          // loop over dimensions
          for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
          {
            double val =colcurr->second*mxi[dimrow];
            if (abs(val)>1e-14)
              cnode->AddDerivJumpValue(dimrow,col,val);
          }
        }
      }

      if (constr_direction_==INPAR::CONTACT::constr_xyz)
        for (int j=0; j<Dim(); j++)
          if (cnode->DbcDofs()[j]==true)
          {
            cnode->FriData().jump()[j]=0.;
            cnode->FriData().GetDerivJump()[j].clear();
          }

    } // active nodes
  } // loop over slave nodes
  return;
}

/*----------------------------------------------------------------*
 |  Assemble slave coordinates (xs)                 gitterle 10/09|
 *----------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleSlaveCoord(Teuchos::RCP<Epetra_Vector>& xsmod)
{

  // loop over all slave nodes
  for (int j=0; j<snoderowmap_->NumMyElements(); ++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    int dim = cnode->NumDof();

    Epetra_SerialDenseVector xspatial(dim);
    std::vector<int> dof(dim);
    std::vector<int> owner(dim);

    for( int k=0; k<dim; ++k )
    {
      xspatial(k) = cnode->xspatial()[k];
      dof[k] = cnode->Dofs()[k];
      owner[k] = cnode->Owner();
    }

    // do assembly
    LINALG::Assemble(*xsmod, xspatial, dof, owner);
  }
  return;
}


/*----------------------------------------------------------------------*
 |  calculate nodal distances (public)                     pfaller Jan15|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::EvaluateDistances(
    const Teuchos::RCP<const Epetra_Vector>& vec,
    std::map<int, std::vector<double> >& mynormals,
    std::map<int,std::vector<GEN::pairedvector<int,double> > >& dmynormals,
    std::map<int,double>& mygap,
    std::map<int,std::map<int,double> >& dmygap)
{
  SetState(MORTAR::state_new_displacement,*vec);
  Initialize();

  // interface needs to be complete
  if (!Filled() && Comm().MyPID() == 0)
    dserror("ERROR: FillComplete() not called on interface %", id_);

  // get out of here if not participating in interface
  if (!lComm())
    return;

  // create an interpolator instance
  Teuchos::RCP<NTS::CoInterpolator> interpolator =
      Teuchos::rcp(new NTS::CoInterpolator(imortar_,Dim()));

  // create normals
  PreEvaluate(-1,-1); // dummy values

  // loop over proc's slave elements of the interface for integration
  // use standard column map to include processor's ghosted elements
  Comm().Barrier();

  for (int i = 0; i < selecolmap_->NumMyElements(); ++i)
  {
    int gid1 = selecolmap_->GID(i);
    DRT::Element* ele1 = idiscret_->gElement(gid1);
    if (!ele1)
      dserror("ERROR: Cannot find slave element with gid %", gid1);
    CONTACT::CoElement* selement = dynamic_cast<CONTACT::CoElement*>(ele1);

    if(selement->MoData().NumSearchElements()<1)
    {
      std::cout << "WARNING: No elements found!" << std::endl;
      continue;
    }

    // skip zero-sized nurbs elements (slave)
    if (selement->ZeroSized())
      continue;

    // empty vector of master element pointers
    std::vector<CONTACT::CoElement*> melements;

    // loop over the candidate master elements of sele_
    // use slave element's candidate list SearchElements !!!
    for (int j = 0; j < selement->MoData().NumSearchElements(); ++j)
    {
      int gid2 = selement->MoData().SearchElements()[j];
      DRT::Element* ele2 = idiscret_->gElement(gid2);
      if (!ele2)
        dserror("ERROR: Cannot find master element with gid %", gid2);
      CONTACT::CoElement* melement = dynamic_cast<CONTACT::CoElement*>(ele2);

      // skip zero-sized nurbs elements (master)
      if (melement->ZeroSized())
        continue;

      melements.push_back(melement);
    }

    //**************************************************************
    //                loop over all Slave nodes
    //**************************************************************
    for(int snodes = 0; snodes<selement->NumNode() ;++snodes)
    {
      CONTACT::CoNode* mynode = dynamic_cast<CONTACT::CoNode*>(selement->Nodes()[snodes]);

      // skip this node if already considered
      if (mynode->HasProj())
        continue;

      //                store node normals
      //**************************************************************
//      int gid = snoderowmapbound_->GID(snodes);
      int gid = mynode->Id();

      int numdofs = mynode->NumDof();
      std::vector<double> temp(numdofs, 0.0);
      for (int kk = 0; kk < numdofs; kk++)
      {
        temp[kk] = mynode->MoData().n()[kk];
      }
      mynormals.insert(std::pair<int, std::vector<double> >(gid, temp));
      dmynormals.insert(std::pair<int, std::vector<GEN::pairedvector<int,double> > >(gid, mynode->CoData().GetDerivN()));

      //**************************************************************
      double sxi[2] = {0.0, 0.0};

      if(selement->Shape() == DRT::Element::quad4 or
         selement->Shape() == DRT::Element::quad8 or
         selement->Shape() == DRT::Element::quad9)
      {
        // TODO (pfaller): switch case
        if(snodes==0)       {sxi[0] = -1; sxi[1] = -1;}
        else if(snodes==1)  {sxi[0] =  1; sxi[1] = -1;}
        else if(snodes==2)  {sxi[0] =  1; sxi[1] =  1;}
        else if(snodes==3)  {sxi[0] = -1; sxi[1] =  1;}
        else if(snodes==4)  {sxi[0] =  0; sxi[1] = -1;}
        else if(snodes==5)  {sxi[0] =  1; sxi[1] =  0;}
        else if(snodes==6)  {sxi[0] =  0; sxi[1] =  1;}
        else if(snodes==7)  {sxi[0] = -1; sxi[1] =  0;}
        else if(snodes==8)  {sxi[0] =  0; sxi[1] =  0;}
        else dserror("ERORR: wrong node LID");
      }
      else if (selement->Shape() == DRT::Element::tri3 or
               selement->Shape() == DRT::Element::tri6)
      {
        if(snodes==0)       {sxi[0] = 0;    sxi[1] = 0;}
        else if(snodes==1)  {sxi[0] = 1;    sxi[1] = 0;}
        else if(snodes==2)  {sxi[0] = 0;    sxi[1] = 1;}
        else if(snodes==3)  {sxi[0] = 0.5;  sxi[1] = 0;}
        else if(snodes==4)  {sxi[0] = 0.5;  sxi[1] = 0.5;}
        else if(snodes==5)  {sxi[0] = 0;    sxi[1] = 0.5;}
        else dserror("ERORR: wrong node LID");
      }
      else
      {
        dserror("ERROR: Chosen element type not supported for NTS!");
      }

      //**************************************************************
      //                loop over all Master Elements
      //**************************************************************
      // create vectors to store projection information for several master elements in case projection is not unique
      std::vector<double> gap_vec;
      std::vector<std::map<int,double> > dgap_vec;

      for (int nummaster=0;nummaster<(int)melements.size();++nummaster)
      {
        // project Gauss point onto master element
        double mxi[2]    = {0.0, 0.0};
        double projalpha =  0.0;
        bool is_projected = MORTAR::MortarProjector::Impl(*selement,*melements[nummaster])->ProjectGaussPoint3D(
            *selement,sxi,*melements[nummaster],mxi,projalpha);

        bool is_on_mele = true;

        // check GP projection
        DRT::Element::DiscretizationType dt = melements[nummaster]->Shape();
        const double tol = 1e-8;
        if (dt==DRT::Element::quad4 || dt==DRT::Element::quad8 || dt==DRT::Element::quad9)
        {
          if (mxi[0]<-1.0-tol || mxi[1]<-1.0-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol)
          {
            is_on_mele=false;
          }
        }
        else
        {
          if (mxi[0]<-tol || mxi[1]<-tol || mxi[0]>1.0+tol || mxi[1]>1.0+tol || mxi[0]+mxi[1]>1.0+2*tol)
          {
            is_on_mele=false;
          }
        }

        // node on mele?
        if (is_on_mele && is_projected)
        {
          // store information of projection so that this node is not considered again
          mynode->HasProj() = true;

          int ndof = 3;
          int ncol = melements[nummaster]->NumNode();
          LINALG::SerialDenseVector mval(ncol);
          LINALG::SerialDenseMatrix mderiv(ncol,2);
          melements[nummaster]->EvaluateShape(mxi,mval,mderiv,ncol,false);

//          int linsize    = mynode->GetLinsize();
          int linsize = 100;
          double gpn[3]  = {0.0, 0.0, 0.0};
          //**************************************************************

          // evalute the GP slave coordinate derivatives --> no entries
          std::vector<GEN::pairedvector<int,double> > dsxi(2,0);
          std::vector<GEN::pairedvector<int,double> > dmxi(2,4*linsize+ncol*ndof);

          (*interpolator).DerivXiGP3D(*selement, *melements[nummaster],sxi,mxi, dsxi, dmxi, projalpha);
          (*interpolator).nwGap3D(*mynode, *melements[nummaster], mval, mderiv, dmxi, gpn);

          // store linearization for node
          std::map<int,double> dgap = mynode->CoData().GetDerivGnts(); // (dof,value)

          // store gap information
          gap_vec.push_back(mynode->CoData().Getgnts());
          dgap_vec.push_back(dgap);

          // reset nodal weighted gap and derivative
          mynode->CoData().Getgnts() = 1.0e12;
          (mynode->CoData().GetDerivGnts()).clear();
        }//End hit ele
      }//End Loop over all Master Elements

      if (gap_vec.size()>0)
      {
        // find projection with smallest absoluate value of gap
        std::vector<double>::iterator iter_min = std::min_element(gap_vec.begin(), gap_vec.end(), abs_compare);
        const int i_min = std::distance(gap_vec.begin(), iter_min);

        // save to map at GID
        mygap.insert(std::pair<int, double>(gid, gap_vec[i_min]));
        dmygap.insert(std::pair<int, std::map<int,double> >(gid, dgap_vec[i_min]));
      }
    }
  }

  Comm().Barrier();

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate L2 Norm of tangential contact conditions     gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::EvaluateTangentNorm(double& cnormtan)
{
  // friction coefficient
  double frcoeff = IParams().get<double>("FRCOEFF");

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    // get some information from node
    double* n = cnode->MoData().n();
    int dim = cnode->NumDof();

    // tangential plane
    Epetra_SerialDenseMatrix tanplane(dim,dim);
    if (dim==3)
    {
      tanplane(0,0)= 1-(n[0]*n[0]);
      tanplane(0,1)=  -(n[0]*n[1]);
      tanplane(0,2)=  -(n[0]*n[2]);
      tanplane(1,0)=  -(n[1]*n[0]);
      tanplane(1,1)= 1-(n[1]*n[1]);
      tanplane(1,2)=  -(n[1]*n[2]);

      tanplane(2,0)=  -(n[2]*n[0]);
      tanplane(2,1)=  -(n[2]*n[1]);
      tanplane(2,2)= 1-(n[2]*n[2]);
    }
    else if (dim==2)
    {
      tanplane(0,0)= 1-(n[0]*n[0]);
      tanplane(0,1)=  -(n[0]*n[1]);

      tanplane(1,0)=  -(n[1]*n[0]);
      tanplane(1,1)= 1-(n[1]*n[1]);
    }
    else
      dserror("Error in AssembleTangentForces: Unknown dimension.");

    // jump vector
    Epetra_SerialDenseMatrix jumpvec(dim,1);
    for (int i=0;i<dim;i++)
      jumpvec(i,0) = cnode->FriData().jump()[i];

    // evaluate jump in tangential direction
    Epetra_SerialDenseMatrix jumptan(dim,1);
    jumptan.Multiply('N','N',1,tanplane,jumpvec,0.0);

    // force vector
    Epetra_SerialDenseMatrix forcevec(dim,1);
    for (int i=0;i<dim;i++)
      forcevec(i,0) = cnode->MoData().lm()[i];

    // evaluate force in normal direction
    double forcen = 0.0;
    for (int k=0;k<dim;++k)
      forcen += forcevec(k,0) * n[k];

    // norm of constraint violation for stick nodes
    if (cnode->Active()== true and cnode->FriData().Slip()==false)
    {
      for (int j=0;j<dim;j++)
        cnormtan += jumptan(j,0)*jumptan(j,0);
    }
    // norm of constraint violation for slip nodes
    else if (cnode->Active()== true and cnode->FriData().Slip()==true)
    {
      double part1 = 0.0;
      double jumpnorm = 0.0;

      for (int j=0;j<dim;j++)
      {
        jumpnorm += jumptan(j,0)*jumptan(j,0);
        part1    += jumptan(j,0)*forcevec(j,0);
      }

      jumpnorm = sqrt(jumpnorm);
      cnormtan += (part1 - frcoeff*forcen*jumpnorm)*(part1 - frcoeff*forcen*jumpnorm);
    }
  } // loop over slave nodes

  // get cnorm from all procs
  double sumcnormtanallprocs=0.0;
  Comm().SumAll(&cnormtan,&sumcnormtanallprocs,1);
  cnormtan=sumcnormtanallprocs;

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate regularized normal forces (nodes)                popp 05/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleRegNormalForces(bool& localisincontact,
                                                 bool& localactivesetchange)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // penalty parameter
  double pp = IParams().get<double>("PENALTYPARAM");

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    int dim = cnode->NumDof();
    double gap = cnode->CoData().Getg();

    // modified gap for zero initial gap
    // (if gap is below zero, it is explicitly set to zero)
    //double modgap = cnode->CoData().Getg();
    //if (abs(modgap) < 1.0e-10) modgap=0.0;

    double kappa = cnode->CoData().Kappa();

    double lmuzawan = 0.0;
    for (int k=0;k<dim;++k)
      lmuzawan += cnode->MoData().lmuzawa()[k]*cnode->MoData().n()[k];

#ifdef CONTACTFDPENALTYKC1
    // set lagrangian multipliers explicitly to constant
    // and corresponding derivatives to zero

    for( int j=0;j<dim;++j)
      cnode->MoData().lm()[j] = i*j;

    cnode->CoData().GetDerivZ().clear();

    continue;
#endif

    //********************************************************************
    // Decision on active /  inactive nodes (regularization)
    //
    // CASE 1: Penalty approach
    // A node is activated if its weighted gap is negative or deactivated
    // if its gap is equal zero or positive.
    // -> the regularization reads: lambda_n = kappa * pp * < -gap >
    //
    // CASE 2: Uzawa augmented Lagrange approach
    // A node is activated if its Lagrange multiplier, stemming from the
    // last Uzawa Lagrange multiplier AND the current regularization is
    // negative or deactivated if its LM is equal zero or positive.
    // -> the regularization reads: lambda_n = < lmuzawa_n - kappa * pp * gap >
    //
    // As the Uzawa Lagrange multipliers are zero in the penalty approach,
    // the two cases can formally be treted identically, see below.
    // We do not need an explicit separation of cases!
    //
    //********************************************************************

    // Activate/Deactivate node and notice any change
    if( (cnode->Active() == false) && (lmuzawan - kappa * pp * gap >= 0) )
    {
        cnode->Active() = true;
        localactivesetchange = true;

        //std::cout << "node #" << gid << " is now active (";
        //for( int j=0; j<dim; j++)
        //  std::cout << " " << cnode->Dofs()[j] << " ";
        //std::cout << ") gap=" << gap << std::endl;
    }

    else if( (cnode->Active() == true) && (lmuzawan - kappa * pp * gap < 0) )
    {
        cnode->Active() = false;
        localactivesetchange = true;

        //std::cout << "node #" << gid << " is now inactive, gap=" << gap << std::endl;
    }
    //********************************************************************

    // Compute derivZ-entries with the Macauley-Bracket
    // of course, this is only done for active constraints in order
    // for linearization and r.h.s to match!
    if( cnode->Active()==true )
    {

//      std::cout << "GID " << gid << std::endl;
//      std::cout << "LMUZAWAN " << lmuzawan << std::endl;
//      std::cout << "GAP " << gap << std::endl;

      localisincontact = true;

      double* normal = cnode->MoData().n();

      // compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->MoData().lm()[j] = (lmuzawan - kappa * pp * gap) * normal[j];

      // compute derivatives of lagrange multipliers and store into node

      // contribution of derivative of weighted gap
      std::map<int,double>& derivg = cnode->CoData().GetDerivG();
      std::map<int,double>::iterator gcurr;

      // contribution of derivative of normal
      std::vector<GEN::pairedvector<int,double> >& derivn = cnode->CoData().GetDerivN();
      GEN::pairedvector<int,double>::iterator ncurr;

      for( int j=0;j<dim;++j)
      {
        for( gcurr = derivg.begin(); gcurr != derivg.end(); ++gcurr )
          cnode->AddDerivZValue(j, gcurr->first, - kappa * pp * (gcurr->second) * normal[j]);
        for( ncurr = (derivn[j]).begin(); ncurr != (derivn[j]).end(); ++ncurr )
          cnode->AddDerivZValue(j, ncurr->first, - kappa * pp * gap * ncurr->second);
        for( ncurr = (derivn[j]).begin(); ncurr != (derivn[j]).end(); ++ncurr )
          cnode->AddDerivZValue(j, ncurr->first, + lmuzawan * ncurr->second);
      }
    }

    // be sure to remove all LM-related stuff from inactive nodes
    else
    {
      // clear lagrange multipliers
      for( int j=0;j<dim;++j) cnode->MoData().lm()[j] = 0.0;

      // clear derivz
      cnode->CoData().GetDerivZ().clear();

    } // Macauley-Bracket
  } // loop over slave nodes

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate regularized tangential forces                gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleRegTangentForcesPenalty()
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // penalty parameter in tangential direction
  double ppnor   = IParams().get<double>("PENALTYPARAM");
  double pptan   = IParams().get<double>("PENALTYPARAMTAN");
  double frcoeff = IParams().get<double>("FRCOEFF");

  INPAR::CONTACT::FrictionType ftype =
    DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(IParams(),"FRICTION");

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    // get some informatiom form the node
    double gap = cnode->CoData().Getg();
    int dim = cnode->NumDof();
    double kappa = cnode->CoData().Kappa();
    double* n = cnode->MoData().n();

    // Lagrange multiplier from Uzawa algorithm
    Epetra_SerialDenseMatrix lmuzawa(dim,1);
    for (int k=0;k<dim;++k)
      lmuzawa(k,0) = cnode->MoData().lmuzawa()[k];

    // Lagrange multiplier in normal direction
    double lmuzawan = 0.0;
    for (int k=0;k<dim;++k)
      lmuzawan += cnode->MoData().lmuzawa()[k]*cnode->MoData().n()[k];

    // tangential plane
    Epetra_SerialDenseMatrix tanplane(dim,dim);
    if (dim==3)
     {
       tanplane(0,0)= 1-(n[0]*n[0]);
       tanplane(0,1)=  -(n[0]*n[1]);
       tanplane(0,2)=  -(n[0]*n[2]);
       tanplane(1,0)=  -(n[1]*n[0]);
       tanplane(1,1)= 1-(n[1]*n[1]);
       tanplane(1,2)=  -(n[1]*n[2]);

       tanplane(2,0)=  -(n[2]*n[0]);
       tanplane(2,1)=  -(n[2]*n[1]);
       tanplane(2,2)= 1-(n[2]*n[2]);
     }
     else if (dim==2)
     {
       tanplane(0,0)= 1-(n[0]*n[0]);
       tanplane(0,1)=  -(n[0]*n[1]);

       tanplane(1,0)=  -(n[1]*n[0]);
       tanplane(1,1)= 1-(n[1]*n[1]);
     }
     else
       dserror("Error in AssembleTangentForces: Unknown dimension.");

    // evaluate traction
    Epetra_SerialDenseMatrix jumpvec(dim,1);

    for (int i=0;i<dim;i++)
      jumpvec(i,0) = cnode->FriData().jump()[i];

    // evaluate kappa.pptan.jumptan
    Epetra_SerialDenseMatrix temptrac(dim,1);
    temptrac.Multiply('N','N',kappa*pptan,tanplane,jumpvec,0.0);

    // fill vector tractionold
    std::vector<double> tractionold(dim);
    for (int i=0;i<dim;i++)
      tractionold[i] = cnode->FriData().tractionold()[i];

    // Evaluate trailtraction (tractionold+temptrac in penalty case)
    std::vector<double> trailtraction(dim);
    double magnitude = 0;
    for (int i=0;i<dim;i++)
    {
      trailtraction[i]=tractionold[i]+temptrac(i,0);
      magnitude += (trailtraction[i]*trailtraction[i]);
    }

    // evaluate magnitude of trailtraction
    magnitude = sqrt(magnitude);

    // evaluate maximal tangential traction
    double maxtantrac = frcoeff*(lmuzawan - kappa * ppnor * gap);

    if(cnode->Active()==false)
    {
      // do nothing
      cnode->FriData().Slip() = false;
    }
    else if (cnode->Active()==true &&
        ((abs(maxtantrac) - magnitude >= 0) or ftype==INPAR::CONTACT::friction_stick))
    {
      //std::cout << "Node " << gid << " is stick" << std::endl;
      cnode->FriData().Slip() = false;

      // in the stick case, traction is trailtraction
      for (int i=0;i<dim;i++)
        cnode->FriData().traction()[i]=trailtraction[i];

      // compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->MoData().lm()[j] = n[j]*(- kappa * ppnor * gap) + trailtraction[j];
    }
    else
    {
      //std::cout << "Node " << gid << " is slip" << std::endl;
      cnode->FriData().Slip() = true;

      // in the slip case, traction is evaluated with a return map algorithm
      for (int i=0;i<dim;i++)
        cnode->FriData().traction()[i]=maxtantrac/magnitude*trailtraction[i];

      // compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->MoData().lm()[j] = n[j]*(- kappa * ppnor * gap)+maxtantrac/magnitude*trailtraction[j];
    }

    // linearization of contact forces (lagrange multipliers)
    // this consists the linearization of the tangential part,
    // the normal part was already done in AssembleRegNormalTraction

    // stick nodes
    if(cnode->Active() == true && cnode->FriData().Slip() == false)
    {
      /***************************************** tanplane.deriv(jump) ***/
      std::vector<std::map<int,double> >& derivjump = cnode->FriData().GetDerivJump();
      std::map<int,double>::iterator      colcurr;
      GEN::pairedvector<int,double>::iterator _colcurr;

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
          {
              int col = colcurr->first;
              double val =pptan*kappa*(colcurr->second)*tanplane(dimrow,dim);
              cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      /**************************************** deriv(tanplane).jump  ***/
      std::vector<GEN::pairedvector<int,double> >& derivn = cnode->CoData().GetDerivN();

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        // loop over all entries of the current derivative map
        for (_colcurr=derivn[dimrow].begin();_colcurr!=derivn[dimrow].end();++_colcurr)
        {
          for (int dim =0;dim<cnode->NumDof();++dim)
          {
            int col = _colcurr->first;
            double val =-pptan*kappa*(_colcurr->second)*n[dim]*cnode->FriData().jump()[dim];
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      // loop over dimensions
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        // loop over all entries of the current derivative map
        for (_colcurr=derivn[dim].begin();_colcurr!=derivn[dim].end();++_colcurr)
        {
          for (int dimrow =0;dimrow<cnode->NumDof();++dimrow)
          {
            int col = _colcurr->first;
            double val =-pptan*kappa*(_colcurr->second)*n[dimrow]*cnode->FriData().jump()[dim];
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }
    }
    // slip nodes
    else if (cnode->Active() == true && cnode->FriData().Slip()== true)
    {
      /******************** tanplane.deriv(jump).maxtantrac/magnidude ***/

      std::vector<std::map<int,double> >& derivjump = cnode->FriData().GetDerivJump();
      std::map<int,double>::iterator      colcurr;
      GEN::pairedvector<int,double>::iterator _colcurr;

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
          {
              int col = colcurr->first;
              double val =pptan*kappa*(colcurr->second)*tanplane(dimrow,dim)*maxtantrac/magnitude;
              cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      /******************** deriv(tanplane).jump.maxtantrac/magnitude ***/
      std::vector<GEN::pairedvector<int,double> >& derivn = cnode->CoData().GetDerivN();

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        // loop over all entries of the current derivative map
        for (_colcurr=derivn[dimrow].begin();_colcurr!=derivn[dimrow].end();++_colcurr)
        {
          for (int dim =0;dim<cnode->NumDof();++dim)
          {
            int col = _colcurr->first;
            double val =-pptan*kappa*(_colcurr->second)*n[dim]*cnode->FriData().jump()[dim]*maxtantrac/magnitude;
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }
      // loop over dimensions
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        // loop over all entries of the current derivative map
        for (_colcurr=derivn[dim].begin();_colcurr!=derivn[dim].end();++_colcurr)
        {
          for (int dimrow =0;dimrow<cnode->NumDof();++dimrow)
          {
            int col = _colcurr->first;
            double val =-pptan*kappa*(_colcurr->second)*n[dimrow]*cnode->FriData().jump()[dim]*maxtantrac/magnitude;
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      /******************** tanplane.jump.deriv(maxtantrac)/magnitude ***/
      std::map<int,double>& derivg = cnode->CoData().GetDerivG();
      std::map<int,double>::iterator gcurr;

      for( int j=0;j<cnode->NumDof();++j)
      {
        for( gcurr = derivg.begin(); gcurr != derivg.end(); ++gcurr )
        {
          cnode->AddDerivZValue(j, gcurr->first, - frcoeff*kappa * ppnor * (gcurr->second)*trailtraction[j]/magnitude);
        }
      }

      /******************** tanplane.jump.maxtantrac/deriv(magnitude) ***/
      // vector double temp
      std::vector <double> temp(cnode->NumDof());
      for (int dim=0;dim<cnode->NumDof();++dim)
        temp[dim] = -maxtantrac/(magnitude*magnitude)*trailtraction[dim];

      // loop over dimensions
      for (int dimout=0;dimout<cnode->NumDof();++dimout)
      {
        double traction = 0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          traction += tanplane(dimout,dim)*cnode->FriData().jump()[dim]*kappa*pptan;

        traction+= tractionold[dimout];

        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
          {
            int col = colcurr->first;
            double val = tanplane(dimout,dim)*pptan*kappa*(colcurr->second)*traction/magnitude;

            for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
            {
              double val1 = val*temp[dimrow];
              cnode->AddDerivZValue(dimrow,col,val1);
            }
          }
        }
      }

      // loop over dimensions
      for (int dimout=0;dimout<cnode->NumDof();++dimout)
      {
        double traction = 0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          traction += tanplane(dimout,dim)*cnode->FriData().jump()[dim]*kappa*pptan;

        traction+=tractionold[dimout];

        // loop over all entries of the current derivative map
        for (_colcurr=derivn[dimout].begin();_colcurr!=derivn[dimout].end();++_colcurr)
        {
          int col = _colcurr->first;

          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            double val =-_colcurr->second*n[dim]*cnode->FriData().jump()[dim]*traction/magnitude*pptan*kappa;
            for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
            {
              double val1 = val*temp[dimrow];
              cnode->AddDerivZValue(dimrow,col,val1);
            }
          }
        }
      }

      // loop over dimensions
      for (int dimout=0;dimout<cnode->NumDof();++dimout)
      {
        double traction = 0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          traction += tanplane(dimout,dim)*cnode->FriData().jump()[dim]*kappa*pptan;

        traction += tractionold[dimout];

        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (_colcurr=derivn[dim].begin();_colcurr!=derivn[dim].end();++_colcurr)
          {
            int col = _colcurr->first;

            double val =-_colcurr->second*n[dimout]*cnode->FriData().jump()[dim]*traction/magnitude*pptan*kappa;

            for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
            {
               double val1 = val*temp[dimrow];
               cnode->AddDerivZValue(dimrow,col,val1);
            }
          }
        }
      }
    } // if Slip == true
    else
    {
      // clear tractions
      for( int j=0;j<dim;++j) cnode->MoData().lm()[j] = 0;
      // clear derivz
      cnode->CoData().GetDerivZ().clear();
    }
  } // loop over active nodes
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate regularized tangential forces                gitterle 10/09|
 |  (Uzawa Aug. Lagr.)                                                  |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleRegTangentForcesUzawa()
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // penalty parameter in tangential direction
  double ppnor = IParams().get<double>("PENALTYPARAM");
  double pptan = IParams().get<double>("PENALTYPARAMTAN");
  double frcoeff = IParams().get<double>("FRCOEFF");

  INPAR::CONTACT::FrictionType ftype =
    DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(IParams(),"FRICTION");

  // loop over all slave row nodes on the current interface
  for (int i=0; i<SlaveRowNodes()->NumMyElements(); ++i)
  {
    int gid = SlaveRowNodes()->GID(i);
    DRT::Node* node = Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    // get some informatiom form the node
    double gap = cnode->CoData().Getg();
    int dim = cnode->NumDof();
    double kappa = cnode->CoData().Kappa();
    double* n = cnode->MoData().n();

    // Lagrange multiplier from Uzawa algorithm
    Epetra_SerialDenseMatrix lmuzawa(dim,1);
    for (int k=0;k<dim;++k)
      lmuzawa(k,0) = cnode->MoData().lmuzawa()[k];

    // Lagrange multiplier in normal direction
    double lmuzawan = 0.0;
    for (int k=0;k<dim;++k)
      lmuzawan += cnode->MoData().lmuzawa()[k]*cnode->MoData().n()[k];

    // tangential plane
    Epetra_SerialDenseMatrix tanplane(dim,dim);
    if (dim==3)
    {
      tanplane(0,0)= 1-(n[0]*n[0]);
      tanplane(0,1)=  -(n[0]*n[1]);
      tanplane(0,2)=  -(n[0]*n[2]);
      tanplane(1,0)=  -(n[1]*n[0]);
      tanplane(1,1)= 1-(n[1]*n[1]);
      tanplane(1,2)=  -(n[1]*n[2]);

      tanplane(2,0)=  -(n[2]*n[0]);
      tanplane(2,1)=  -(n[2]*n[1]);
      tanplane(2,2)= 1-(n[2]*n[2]);
    }
    else if (dim==2)
    {
      tanplane(0,0)= 1-(n[0]*n[0]);
      tanplane(0,1)=  -(n[0]*n[1]);

      tanplane(1,0)=  -(n[1]*n[0]);
      tanplane(1,1)= 1-(n[1]*n[1]);
    }
    else
      dserror("Error in AssembleTangentForces: Unknown dimension.");

    // Lagrange multiplier in tangential direction
    Epetra_SerialDenseMatrix lmuzawatan(dim,1);
    lmuzawatan.Multiply('N','N',1,tanplane,lmuzawa,0.0);

    // evaluate traction
    Epetra_SerialDenseMatrix jumpvec(dim,1);

    for (int i=0;i<dim;i++)
      jumpvec(i,0) = cnode->FriData().jump()[i];

    // evaluate kappa.pptan.jumptan
    Epetra_SerialDenseMatrix temptrac(dim,1);
    temptrac.Multiply('N','N',kappa*pptan,tanplane,jumpvec,0.0);

    // Evaluate trailtraction
    std::vector<double> trailtraction(dim);
    double magnitude = 0;
    for (int i=0;i<dim;i++)
    {
      trailtraction[i]=lmuzawatan(i,0)+temptrac(i,0);
      magnitude += (trailtraction[i]*trailtraction[i]);
    }

    // evaluate magnitude of trailtraction
    magnitude = sqrt(magnitude);

    // evaluate maximal tangential traction
    double maxtantrac = frcoeff*(lmuzawan - kappa * ppnor * gap);

    if(cnode->Active()==false)
    {
    }
    else if (cnode->Active()==true && ((abs(maxtantrac) - magnitude >= 0)or ftype==INPAR::CONTACT::friction_stick))    {
      //std::cout << "Node " << gid << " is stick" << std::endl;
      cnode->FriData().Slip() = false;

      // compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->MoData().lm()[j] = n[j]*(lmuzawan - kappa * ppnor * gap)+trailtraction[j];
    }
    else
    {
      //std::cout << "Node " << gid << " is slip" << std::endl;
      cnode->FriData().Slip() = true;

      // compute lagrange multipliers and store into node
      for( int j=0;j<dim;++j)
        cnode->MoData().lm()[j] = n[j]*(lmuzawan - kappa * ppnor * gap)+trailtraction[j]*maxtantrac/magnitude;
    }

    // linearization of contact forces (lagrange multipliers)
    // this consists the linearization of the tangential part,
    // the normal part was already done in AssembleRegNormalTraction

    // stick nodes
    if(cnode->Active() == true && cnode->FriData().Slip() == false)
    {
      /***************************************** tanplane.deriv(jump) ***/
      std::vector<std::map<int,double> >& derivjump = cnode->FriData().GetDerivJump();
      std::map<int,double>::iterator      colcurr;
      GEN::pairedvector<int,double>::iterator _colcurr;

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
          {
            int col = colcurr->first;
            double val =pptan*kappa*(colcurr->second)*tanplane(dimrow,dim);
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      /******************************* deriv(tanplane).(lmuzawa+jump) ***/
      std::vector<GEN::pairedvector<int,double> >& derivn = cnode->CoData().GetDerivN();

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        // loop over all entries of the current derivative map
        for (_colcurr=derivn[dimrow].begin();_colcurr!=derivn[dimrow].end();++_colcurr)
        {
          for (int dim =0;dim<cnode->NumDof();++dim)
          {
            int col = _colcurr->first;
            double val =-pptan*kappa*(_colcurr->second)*n[dim]*(cnode->FriData().jump()[dim]);
            val = val - (_colcurr->second)*n[dim]*(cnode->MoData().lmuzawa()[dim]);
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      // loop over dimensions
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        // loop over all entries of the current derivative map
        for (_colcurr=derivn[dim].begin();_colcurr!=derivn[dim].end();++_colcurr)
        {
          for (int dimrow =0;dimrow<cnode->NumDof();++dimrow)
          {
            int col = _colcurr->first;
            double val =-pptan*kappa*(_colcurr->second)*n[dimrow]*(cnode->FriData().jump()[dim]);
            val = val-(_colcurr->second)*n[dimrow]*(cnode->MoData().lmuzawa()[dim]);
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }
    }

    // slip nodes
    else if (cnode->Active() == true && cnode->FriData().Slip()== true)
    {
      /***************************************** tanplane.deriv(jump) ***/
      std::vector<std::map<int,double> >& derivjump = cnode->FriData().GetDerivJump();
      std::map<int,double>::iterator      colcurr;
      GEN::pairedvector<int,double>::iterator _colcurr;

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
          {
            int col = colcurr->first;
            double val =pptan*kappa*(colcurr->second)*tanplane(dimrow,dim)*maxtantrac/magnitude;
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      /******************************* deriv(tanplane).(lmuzawa+jump) ***/
      std::vector<GEN::pairedvector<int,double> >& derivn = cnode->CoData().GetDerivN();

      // loop over dimensions
      for (int dimrow=0;dimrow<cnode->NumDof();++dimrow)
      {
        // loop over all entries of the current derivative map
        for (_colcurr=derivn[dimrow].begin();_colcurr!=derivn[dimrow].end();++_colcurr)
        {
          for (int dim =0;dim<cnode->NumDof();++dim)
          {
            int col = _colcurr->first;
            double val =-pptan*kappa*(_colcurr->second)*n[dim]*cnode->FriData().jump()[dim];
            val = (val - (_colcurr->second)*n[dim]*(cnode->MoData().lmuzawa()[dim]))*maxtantrac/magnitude;
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      // loop over dimensions
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        // loop over all entries of the current derivative map
        for (_colcurr=derivn[dim].begin();_colcurr!=derivn[dim].end();++_colcurr)
        {
          for (int dimrow =0;dimrow<cnode->NumDof();++dimrow)
          {
            int col = _colcurr->first;
            double val =-pptan*kappa*(_colcurr->second)*n[dimrow]*cnode->FriData().jump()[dim];
            val = (val-(_colcurr->second)*n[dimrow]*(cnode->MoData().lmuzawa()[dim]))*maxtantrac/magnitude;
            cnode->AddDerivZValue(dimrow,col,val);
          }
        }
      }

      /******************** tanplane.jump.deriv(maxtantrac)/magnitude ***/
      std::map<int,double>& derivg = cnode->CoData().GetDerivG();
      std::map<int,double>::iterator gcurr;

      for( int j=0;j<cnode->NumDof();++j)
      {
        for( gcurr = derivg.begin(); gcurr != derivg.end(); ++gcurr )
        {
          cnode->AddDerivZValue(j,gcurr->first,- frcoeff*kappa*ppnor*(gcurr->second)*trailtraction[j]/magnitude);
        }
      }

      for( int j=0;j<cnode->NumDof();++j)
      {
        for( _colcurr = (derivn[j]).begin(); _colcurr != (derivn[j]).end(); ++_colcurr )
        {
          for( int k=0;k<cnode->NumDof();++k)
          {
            double val = frcoeff*(_colcurr->second)*lmuzawa(j,0)*trailtraction[k]/magnitude;
            cnode->AddDerivZValue(k,_colcurr->first,val);
          }
        }
      }

      /******************** tanplane.jump.maxtantrac/deriv(magnitude) ***/
      // vector double temp
      std::vector <double> temp(cnode->NumDof());
      for (int dim=0;dim<cnode->NumDof();++dim)
        temp[dim] = -maxtantrac/(magnitude*magnitude)*trailtraction[dim];

      // loop over dimensions
      for (int dimout=0;dimout<cnode->NumDof();++dimout)
      {
        double traction = 0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          traction += tanplane(dimout,dim)*(lmuzawa(dim,0)+cnode->FriData().jump()[dim]*kappa*pptan);

        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
          {
            int col = colcurr->first;
            double val = tanplane(dimout,dim)*pptan*kappa*(colcurr->second)*traction/magnitude;

            for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
            {
              double val1 = val*temp[dimrow];
              cnode->AddDerivZValue(dimrow,col,val1);
            }
          }
        }
      }

      // loop over dimensions
      for (int dimout=0;dimout<cnode->NumDof();++dimout)
      {
        double traction = 0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          traction += tanplane(dimout,dim)*(lmuzawa(dim,0)+cnode->FriData().jump()[dim]*kappa*pptan);

        // loop over all entries of the current derivative map
        for (_colcurr=derivn[dimout].begin();_colcurr!=derivn[dimout].end();++_colcurr)
        {
          int col = _colcurr->first;

          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            double val =-_colcurr->second*n[dim]*(lmuzawa(dim,0)+cnode->FriData().jump()[dim]*pptan*kappa)*traction/magnitude;
            for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
            {
              double val1 = val*temp[dimrow];
              cnode->AddDerivZValue(dimrow,col,val1);
            }
          }
        }
      }

      // loop over dimensions
      for (int dimout=0;dimout<cnode->NumDof();++dimout)
      {
        double traction = 0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          traction += tanplane(dimout,dim)*(lmuzawa(dim,0)+cnode->FriData().jump()[dim]*kappa*pptan);

        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map
          for (_colcurr=derivn[dim].begin();_colcurr!=derivn[dim].end();++_colcurr)
          {
            int col = _colcurr->first;
            double val =-_colcurr->second*n[dimout]*(lmuzawa(dim,0)+cnode->FriData().jump()[dim]*pptan*kappa)*traction/magnitude;

            for(int dimrow=0;dimrow<cnode->NumDof();++dimrow)
            {
              double val1 = val*temp[dimrow];
              cnode->AddDerivZValue(dimrow,col,val1);
            }
          }
        }
      }
    } // if Slip == true
    else
    {
      // clear tractions
      for( int j=0;j<dim;++j) cnode->MoData().lm()[j] = 0;
      // clear derivz
      cnode->CoData().GetDerivZ().clear();
    }
  } // loop over active nodes
  return;
}

/*----------------------------------------------------------------------*
 |  Assemble derivatives of lagrange multipliers              popp 05/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinZ(LINALG::SparseMatrix& linzglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // loop over all slave nodes (row map)
  for (int i=0; i<snoderowmap_->NumMyElements(); ++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node)
      dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleLinZ: Node ownership inconsistency!");

    // derivz is the std::vector<map> we want to assemble
    std::vector<std::map<int,double> >& derivz = cnode->CoData().GetDerivZ();

    if ( (int) derivz.size()>0 )
    {
      int rowsize = cnode->NumDof();
      int colsize = (int) derivz[0].size();

      // consistency check
      for (int j=0; j<rowsize-1; ++j)
        if ((int)derivz[j].size() != (int)derivz[j+1].size())
          dserror("ERROR: AssembleLinZ: Column dim. of nodal derivz-map is inconsistent!");

      std::map<int,double>::iterator colcurr;

      // loop over dofs
      for ( int k=0; k<rowsize; ++k )
      {
        int row = cnode->Dofs()[k]; // row index equals global dof index of this #i node's dof k
        int l = 0;

        // loop over all directional derivative entries using the map iterator
        for( colcurr = derivz[k].begin(); colcurr != derivz[k].end(); ++colcurr )
        {
          int col = colcurr->first; // col index equals global id of directional derivative component ,l
          double val = colcurr->second;
          linzglobal.Assemble(val, row, col);
          l++;
        }

        if( l != colsize )
          dserror("ERROR: AssembleLinZ: l = %i but colsize = %i",k,colsize);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix with nodal tangents or/and normals         popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleTN(Teuchos::RCP<LINALG::SparseMatrix> tglobal,
                                      Teuchos::RCP<LINALG::SparseMatrix> nglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (activenodes_==Teuchos::null)
    return;

  if (Dim() != 2 && Dim() != 3)
    dserror("ERROR: Dim() must be either 2 or 3!");

// loop over all active slave nodes of the interface
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleTN: Node ownership inconsistency!");

    if (tglobal != Teuchos::null)
    {
      if (constr_direction_==INPAR::CONTACT::constr_xyz)
      {
        if (Dim()==2)
        {
          // prepare assembly
          std::vector<int> lmrowT(cnode->NumDof());
          std::vector<int> lmrowownerT(cnode->NumDof());
          std::vector<int> lmcol(cnode->NumDof());

          Epetra_SerialDenseMatrix Tnode(cnode->NumDof(),cnode->NumDof());
          for (int i=0; i<cnode->NumDof(); ++i)
          {
            lmrowT[i] = cnode->Dofs()[i];
            lmrowownerT[i] = cnode->Owner();
            for (int j=0; j<cnode->NumDof(); ++j)
            {
              lmcol[j] = cnode->Dofs()[j];
              Tnode(i,j)= cnode->CoData().txi()[i]*cnode->CoData().txi()[j];
            }
          }
          tglobal->Assemble(-1,Tnode,lmrowT,lmrowownerT,lmcol);
        }
        else if (Dim()==3)
        {
          std::vector<int> lmrowT(cnode->NumDof());
          std::vector<int> lmrowownerT(cnode->NumDof());
          std::vector<int> lmcol(cnode->NumDof());

          Epetra_SerialDenseMatrix Tnode(cnode->NumDof(),cnode->NumDof());

          for (int i=0; i<cnode->NumDof(); ++i)
          {
            lmrowT[i]=cnode->Dofs()[i];
            lmrowownerT[i]=cnode->Owner();
            for (int j=0; j<cnode->NumDof(); ++j)
            {
              lmcol[j] = cnode->Dofs()[j];
              Tnode(i,j) = cnode->CoData().txi()[i]*cnode->CoData().txi()[j]
                                                                          +cnode->CoData().teta()[i]*cnode->CoData().teta()[j];
            }
          }
          tglobal->Assemble(-1,Tnode,lmrowT,lmrowownerT,lmcol);
        }
        else
          dserror("ERROR: Dim() must be either 2D or 3D");

      }
      else
      {
        if (Dim()==2)
        {
          // prepare assembly
          int colsize = cnode->NumDof();
          std::vector<int> lmrowT(1);
          std::vector<int> lmrowownerT(1);
          std::vector<int> lmcol(colsize);

          lmrowT[0] = activet_->GID(i);
          lmrowownerT[0] = cnode->Owner();

          /**************************************************** T-matrix ******/
          Epetra_SerialDenseMatrix Tnode(1,colsize);

          for (int j=0;j<colsize;++j)
          {
            lmcol[j] = cnode->Dofs()[j];
            Tnode(0,j) = cnode->CoData().txi()[j];
          }

          // assemble into matrix of normal vectors T
          tglobal->Assemble(-1,Tnode,lmrowT,lmrowownerT,lmcol);
        }

        else if (Dim()==3)
        {
          // prepare assembly
          int colsize = cnode->NumDof();
          std::vector<int> lmrowT(2);
          std::vector<int> lmrowownerT(2);
          std::vector<int> lmcol(colsize);

          lmrowT[0] = activet_->GID(2*i);
          lmrowT[1] = activet_->GID(2*i+1);
          lmrowownerT[0] = cnode->Owner();
          lmrowownerT[1] = cnode->Owner();

          /**************************************************** T-matrix ******/
          Epetra_SerialDenseMatrix Tnode(2,colsize);

          for (int j=0;j<colsize;++j)
          {
            lmcol[j] = cnode->Dofs()[j];
            Tnode(0,j) = cnode->CoData().txi()[j];
            Tnode(1,j) = cnode->CoData().teta()[j];
          }

          // assemble into matrix of normal vectors T
          tglobal->Assemble(-1,Tnode,lmrowT,lmrowownerT,lmcol);
        }
        else
          dserror("ERROR: Dim() must be either 2D or 3D");
      }
    }

    if (nglobal != Teuchos::null)
    {
      // nodal normal
      double* nodalnormal = cnode->MoData().n();

      int row = activen_->GID(i);

      // add normal to corresponding row in global matrix
      for (int k = 0; k < cnode->NumDof(); ++k)
        nglobal->Assemble(nodalnormal[k], row, cnode->Dofs()[k]); //use the first dof for normal direction!!!
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble matrix S containing gap g~ derivatives           popp 02/09|
 |  PS: "AssembleS" is an outdated name which could make                |
 |  you confused.                                                       |
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleS(LINALG::SparseMatrix& sglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (activenodes_==Teuchos::null)
    return;

  // loop over all active slave nodes of the interface
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleS: Node ownership inconsistency!");

    // prepare assembly
    std::map<int,double>& dgmap = cnode->CoData().GetDerivG();
    std::map<int,double>::iterator colcurr;
    int row = activen_->GID(i);

    for (colcurr=dgmap.begin();colcurr!=dgmap.end();++colcurr)
    {
      int col = colcurr->first;
      double val = colcurr->second;

      // do not assemble zeros into s matrix
      if (constr_direction_==INPAR::CONTACT::constr_xyz)
      {
        for (int j=0; j<cnode->NumDof(); j++)
          if (abs(val*cnode->MoData().n()[j])>1.0e-12)
              sglobal.Assemble(val*cnode->MoData().n()[j],cnode->Dofs()[j],col);
      }
      else
        if (abs(val)>1.0e-12) sglobal.Assemble(val,row,col);
    }
  } //for (int i=0;i<activenodes_->NumMyElements();++i)

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble tangent deriv or/and normal deriv matrix         popp 05/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleTNderiv(Teuchos::RCP<LINALG::SparseMatrix> tderivglobal,
                                           Teuchos::RCP<LINALG::SparseMatrix> nderivglobal,
                                           bool usePoroLM)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (activenodes_==Teuchos::null)
    return;

  if (Dim() != 2 && Dim() != 3)
    dserror("ERROR: Dim() must be either 2 or 3!");

  // loop over all active slave nodes of the interface
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID()) //move this check into debug?
      dserror("ERROR: AssembleTNderiv: Node ownership inconsistency!");

    if (tderivglobal != Teuchos::null) //assemble tangential derivs?
    {
      std::vector<GEN::pairedvector<int,double> >& dtmap = cnode->CoData().GetDerivTxi();
      GEN::pairedvector<int,double>::iterator colcurr;

      for (int dim = 0; dim < Dim() -1; ++dim) //for both tangents
      {
        if (dim == 1) //just for 3d case, 2nd tangent
          dtmap = cnode->CoData().GetDerivTeta();
        int colsize = (int)dtmap[0].size();
        int mapsize = (int)dtmap.size();

        int row = activet_->GID((Dim()-1)*i+dim);

        if (Dim() == 2 && mapsize==3)
          mapsize=2; //??

        for (int j=0;j<mapsize-1;++j) //move this check into debug?
          if ((int)dtmap[j].size() != (int)dtmap[j+1].size())
          {
            std::cout << "size j = " << dtmap[j].size() << "   size j+1 = " << dtmap[j+1].size() << std::endl;
            dserror("ERROR: AssembleTNderiv: Column dim. of nodal DerivT-map is inconsistent!");
          }

        // begin assembly of Tderiv-matrix
        //std::cout << std::endl << "->Assemble P for Node ID: " << cnode->Id() << std::endl;

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val;
            if (!usePoroLM)
              val= cnode->MoData().lm()[j]*(colcurr->second);
            else
              val= cnode->CoPoroData().poroLM()[j]*(colcurr->second);
            //std::cout << "lm[" << j << "]=" << cnode->MoData().lm()[j] << " deriv=" << colcurr->second << std::endl;
            //std::cout << "Assemble P: " << row << " " << col << " " << val << std::endl;
            // do not assemble zeros into P matrix

            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int i=0; i<cnode->NumDof(); ++i)
              {
                if (abs(val)>1.0e-12)
                {
                  double t;
                  if (dim == 0)
                    t = cnode->CoData().txi()[i];
                  else
                    t = cnode->CoData().teta()[i];
                  tderivglobal->Assemble(val*t,cnode->Dofs()[i],col);
                }
              }
            }
            else
              if (abs(val)>1.0e-12) tderivglobal->Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleTNderiv: k = %i but colsize = %i",k,colsize);
        }
      }
    }

    if (nderivglobal != Teuchos::null) //assemble normal derivs?
    {
      std::vector<GEN::pairedvector<int,double> >& dnmap = cnode->CoData().GetDerivN();
      GEN::pairedvector<int,double>::iterator colcurr;

        int colsize = (int)dnmap[0].size();
        int mapsize = (int)dnmap.size();

        int row = activen_->GID(i);

        if (Dim() == 2 && mapsize==3)
          mapsize=2; //??

        for (int j=0;j<mapsize-1;++j) //move this check into debug?
          if ((int)dnmap[j].size() != (int)dnmap[j+1].size())
            dserror("ERROR: AssembleTNderiv: Column dim. of nodal DerivN-map is inconsistent!");

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dnmap[j].begin();colcurr!=dnmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val;
            if (!usePoroLM)
              val= cnode->MoData().lm()[j]*(colcurr->second);
            else
              val= cnode->CoPoroData().poroLM()[j]*(colcurr->second);

            // do not assemble zeros into P matrix

            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int i=0; i<cnode->NumDof(); ++i)
                if (abs(val)>1.0e-12)
                  nderivglobal->Assemble(val*cnode->MoData().n()[i],cnode->Dofs()[i],col);
            }
            else
              if (abs(val)>1.0e-12) nderivglobal->Assemble(val,row,col);
            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleTNderiv: k = %i but colsize = %i",k,colsize);
        }
    }

  } //for (int i=0;i<activenodes_->NumMyElements();++i)

  return;
}


/*----------------------------------------------------------------------*
 |  Assemble matrices LinD containing fc derivatives          popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinD(LINALG::SparseMatrix& lindglobal,
                                       bool usePoroLM)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  /**********************************************************************/
  // NEW VERSION (09/2010): No more communication, thanks to FE_MATRIX!
  /**********************************************************************/
  // we have: D_jk,c with j = Lagrange multiplier slave dof
  //                 with k = Displacement slave dof
  //                 with c = Displacement slave or master dof
  // we compute (LinD)_kc = D_jk,c * z_j
  /**********************************************************************/

  // loop over all LM slave nodes (row map)
  for (int j=0;j<snoderowmap_->NumMyElements();++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);
    int dim = cnode->NumDof();

    // Mortar matrix D and M derivatives
    std::map<int,std::map<int,double> >& dderiv = cnode->CoData().GetDerivD();

    // current Lagrange multipliers
    double* lm;
    if (!usePoroLM)
      lm = cnode->MoData().lm();
    else
      lm = cnode->CoPoroData().poroLM();

    // get sizes and iterator start
    int slavesize = (int)dderiv.size();
    std::map<int,std::map<int,double> >::iterator scurr = dderiv.begin();

    /********************************************** LinDMatrix **********/
    // loop over all DISP slave nodes in the DerivD-map of the current LM slave node
    for (int k=0;k<slavesize;++k)
    {
      int sgid = scurr->first;
      ++scurr;

      DRT::Node* snode = idiscret_->gNode(sgid);
      if (!snode) dserror("ERROR: Cannot find node with gid %",sgid);
      CoNode* csnode = dynamic_cast<CoNode*>(snode);

      // Mortar matrix D derivatives
      std::map<int,double>& thisdderiv = cnode->CoData().GetDerivD()[sgid];
      int mapsize = (int)(thisdderiv.size());

      // inner product D_{jk,c} * z_j for index j
      for (int prodj=0;prodj<dim;++prodj)
      {
        int row = csnode->Dofs()[prodj];
        std::map<int,double>::iterator scolcurr = thisdderiv.begin();

        // loop over all directional derivative entries
        for (int c=0;c<mapsize;++c)
        {
          int col = scolcurr->first;
          double val = lm[prodj] * (scolcurr->second);
          ++scolcurr;

          // owner of LM slave node can do the assembly, although it actually
          // might not own the corresponding rows in lindglobal (DISP slave node)
          // (FE_MATRIX automatically takes care of non-local assembly inside!!!)
          //std::cout << "Assemble LinD: " << row << " " << col << " " << val << std::endl;
          if (abs(val)>1.0e-12) lindglobal.FEAssemble(val,row,col);
        }

        // check for completeness of DerivD-Derivatives-iteration
        if (scolcurr!=thisdderiv.end())
          dserror("ERROR: AssembleLinDM: Not all derivative entries of DerivD considered!");
      }
    }

    // check for completeness of DerivD-Slave-iteration
    if (scurr!=dderiv.end())
      dserror("ERROR: AssembleLinDM: Not all DISP slave entries of DerivD considered!");
    /******************************** Finished with LinDMatrix **********/
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Assemble matrices LinM containing fc derivatives          popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinM(LINALG::SparseMatrix& linmglobal,
                                       bool usePoroLM)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  /**********************************************************************/
  // NEW VERSION (09/2010): No more communication, thanks to FE_MATRIX!
  /**********************************************************************/
  // we have: M_jl,c with j = Lagrange multiplier slave dof
  //                 with l = Displacement master dof
  //                 with c = Displacement slave or master dof
  // we compute (LinM)_lc = M_jl,c * z_j
  /**********************************************************************/

  // loop over all LM slave nodes (row map)
  for (int j=0;j<snoderowmap_->NumMyElements();++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);
    int dim = cnode->NumDof();

    // Mortar matrix D and M derivatives
    std::map<int,std::map<int,double> >& mderiv = cnode->CoData().GetDerivM();

    // current Lagrange multipliers
    double* lm;
    if (!usePoroLM)
      lm = cnode->MoData().lm();
    else
      lm = cnode->CoPoroData().poroLM();

    // get sizes and iterator start
    int mastersize = (int)mderiv.size();
    std::map<int,std::map<int,double> >::iterator mcurr = mderiv.begin();

    /********************************************** LinMMatrix **********/
    // loop over all master nodes in the DerivM-map of the current LM slave node
    for (int l=0;l<mastersize;++l)
    {
      int mgid = mcurr->first;
      ++mcurr;

      DRT::Node* mnode = idiscret_->gNode(mgid);
      if (!mnode) dserror("ERROR: Cannot find node with gid %",mgid);
      CoNode* cmnode = dynamic_cast<CoNode*>(mnode);

      // Mortar matrix M derivatives
      std::map<int,double>&thismderiv = cnode->CoData().GetDerivM()[mgid];
      int mapsize = (int)(thismderiv.size());

      // inner product M_{jl,c} * z_j for index j
      for (int prodj=0;prodj<dim;++prodj)
      {
        int row = cmnode->Dofs()[prodj];
        std::map<int,double>::iterator mcolcurr = thismderiv.begin();

        // loop over all directional derivative entries
        for (int c=0;c<mapsize;++c)
        {
          int col = mcolcurr->first;
          double val = lm[prodj] * (mcolcurr->second);
          ++mcolcurr;

          // owner of LM slave node can do the assembly, although it actually
          // might not own the corresponding rows in lindglobal (DISP slave node)
          // (FE_MATRIX automatically takes care of non-local assembly inside!!!)
          //std::cout << "Assemble LinM: " << row << " " << col << " " << val << std::endl;
          if (abs(val)>1.0e-12) linmglobal.FEAssemble(-val,row,col);
        }

        // check for completeness of DerivM-Derivatives-iteration
        if (mcolcurr!=thismderiv.end())
          dserror("ERROR: AssembleLinDM: Not all derivative entries of DerivM considered!");
      }
    }

    // check for completeness of DerivM-Master-iteration
    if (mcurr!=mderiv.end())
      dserror("ERROR: AssembleLinDM: Not all master entries of DerivM considered!");
    /******************************** Finished with LinMMatrix **********/
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Assemble matrices LinDM containing fc derivatives        farah 02/16|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinDM(LINALG::SparseMatrix& lindglobal,
                                       LINALG::SparseMatrix& linmglobal,
                                       bool usePoroLM)
{
  // call both sub functions
  AssembleLinD(lindglobal,usePoroLM);
  AssembleLinM(linmglobal,usePoroLM);

  return;
}


/*----------------------------------------------------------------------*
 |  Assemble normal weighted gap                              popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleG(Epetra_Vector& gglobal)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleG: Node ownership inconsistency!");

    /**************************************************** g-vector ******/
    if (cnode->CoData().Getg()!=0.0)
    {
      double gap = cnode->CoData().Getg();

      // std::cout << "Node ID: " << cnode->Id() << " HasProj: " << cnode->HasProj()
      //      << " IsActive: " << cnode->Active() << " Gap: " << gap << std::endl;

      // check if this inactive node has a feasible projection
      // else, it cannot be in contact and weighted gap should be positive
      // (otherwise wrong results possible for g~ because of non-positivity
      // of dual shape functions!!!)
      //******************************************************************
      // TODO: This is only necessary for quadratic LM shape functions!
      // By the way, it makes the method slightly inconsistent
      // (e.g. patch tests with slave side being wider than master side).
      // However, we are able to solve many problems with this little trick.
      // But not all problems, e.g. dropping edge problems would still fail!!!
      // To solve this dilemma, we need a clever modification of the LM shape
      // functions such that they have positive integral values on the
      // "projecting" element part. Once we have this, the following trick
      // can (and should) also be removed in order to make the method
      // consistent again! (08/2013)
      //******************************************************************

      if(!nurbs_) // only for Lagrange elements
      {
        bool node_has_quad_element = false;
        for (int i=0; i<cnode->NumElement(); i++)
        {
          if (dynamic_cast<MORTAR::MortarElement*>(cnode->Elements()[i])->IsQuad()==true)
          {
            node_has_quad_element=true;
            break;
          }
        }
        if (!cnode->HasProj() && !cnode->Active() && node_has_quad_element)
        {
          gap = 1.0e12;
          cnode->CoData().Getg()=gap;
        }
      }

      if (constr_direction_==INPAR::CONTACT::constr_xyz)
      {
        Epetra_SerialDenseVector gnode(Dim());
        std::vector<int> lm(Dim());
        std::vector<int> lmowner(Dim());
        for (int i=0; i<Dim(); i++)
        {
          gnode(i) = gap * cnode->MoData().n()[i];
          lm[i] = cnode->Dofs()[i];
          lmowner[i]=cnode->Owner();
        }
        LINALG::Assemble(gglobal,gnode,lm,lmowner);
      }
      else
      {
        static Epetra_SerialDenseVector gnode(1);
        static std::vector<int> lm(1);
        static std::vector<int> lmowner(1);

        gnode(0)   = gap;
        lm[0]      = cnode->Id();
        lmowner[0] = cnode->Owner();
        LINALG::Assemble(gglobal,gnode,lm,lmowner);
      }

    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Assemble inactive right hand side                    hiermeier 08/13|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleInactiverhs(Epetra_Vector& inactiverhs)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // FIXME It's possible to improve the performance, if only recently active nodes of the inactive node set,
  // i.e. nodes, which were active in the last iteration, are considered. Since you know, that the lagrange
  // multipliers of former inactive nodes are still equal zero.

  Teuchos::RCP<Epetra_Map> inactivenodes  = LINALG::SplitMap(*snoderowmap_, *activenodes_);
  Teuchos::RCP<Epetra_Map> inactivedofs   = LINALG::SplitMap(*sdofrowmap_, *activedofs_);

  static std::vector<int> lm_gid(Dim());
  static std::vector<int> lm_owner(Dim());
  static Epetra_SerialDenseVector lm_i(Dim());

  for (int i=0;i<inactivenodes->NumMyElements();++i)
  {
    int gid = inactivenodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleInactiverhs: Node ownership inconsistency!");
    if (Dim() == 2)
    {
      // calculate the tangential rhs
      for (int j=0; j<Dim(); ++j)
      {
        lm_owner[j] =   cnode->Owner();
        lm_i[j]     = - cnode->MoData().lm()[j];    // already negative rhs!!!
      }
      lm_gid[0]   = inactivedofs->GID(2*i);
      lm_gid[1]   = inactivedofs->GID(2*i+1);

      LINALG::Assemble(inactiverhs, lm_i, lm_gid, lm_owner);
    }
    else if (Dim() == 3)
    {
      // calculate the tangential rhs
      for (int j=0; j<Dim(); ++j)
      {
        lm_owner[j] =   cnode->Owner();
        lm_i[j]     = - cnode->MoData().lm()[j];    // already negative rhs!!!
      }
      lm_gid[0] = inactivedofs->GID(3*i);
      lm_gid[1] = inactivedofs->GID(3*i+1);
      lm_gid[2] = inactivedofs->GID(3*i+2);

      LINALG::Assemble(inactiverhs, lm_i, lm_gid, lm_owner);
    }
  }
}

/*----------------------------------------------------------------------*
 |  Assemble tangential right-hand side                  hiermeier 08/13|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleTangrhs(Epetra_Vector& tangrhs)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  static std::vector<int> lm_gid(Dim()-1);
  static std::vector<int> lm_owner(Dim()-1);
  static Epetra_SerialDenseVector lm_t(Dim()-1);

  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleTangrhs: Node ownership inconsistency!");
    if (constr_direction_==INPAR::CONTACT::constr_xyz)
    {
      if (Dim()==2)
      {
        std::vector<int> lm_gid(2);
        std::vector<int> lm_owner(2);
        Epetra_SerialDenseVector lm_t(2);
        for (int i=0; i<Dim(); ++i)
        {
          lm_gid[i] = cnode->Dofs()[i];
          lm_owner[i]=cnode->Owner();
          lm_t[i]=0.;
          for (int j=0; j<Dim(); ++j)
          {
            lm_t[i] -= cnode->CoData().txi()[i] * cnode->CoData().txi()[j] * cnode->MoData().lm()[j];
          }
        }
        LINALG::Assemble(tangrhs, lm_t, lm_gid, lm_owner);
      }
      else if (Dim()==3)
      {
        std::vector<int> lm_gid(3);
        std::vector<int> lm_owner(3);
        Epetra_SerialDenseVector lm_t(3);

        for (int i=0; i<Dim(); ++i)
        {
          lm_gid[i] = cnode->Dofs()[i];
          lm_owner[i]=cnode->Owner();
          lm_t[i]=0.;
          for (int j=0; j<Dim(); ++j)
          {
            lm_t[i] -= cnode->CoData().txi()[i] * cnode->CoData().txi()[j] * cnode->MoData().lm()[j];
            lm_t[i] -= cnode->CoData().teta()[i] * cnode->CoData().teta()[j] * cnode->MoData().lm()[j];
          }
        }
        LINALG::Assemble(tangrhs, lm_t, lm_gid, lm_owner);

      }
    }
    else
    {
      if (Dim()==2)
      {
        lm_gid[0]   = activet_->GID(i);
        lm_owner[0] = cnode->Owner();

        lm_t[0] = 0.0;
        for (int j=0; j<Dim(); ++j)
          lm_t[0] -= cnode->CoData().txi()[j] * cnode->MoData().lm()[j];   // already negative rhs!!!

        LINALG::Assemble(tangrhs, lm_t, lm_gid, lm_owner);
      }
      else if (Dim()==3)
      {
        lm_gid[0]   = activet_->GID(2*i);   // even
        lm_gid[1]   = activet_->GID(2*i+1);   // odd
        lm_owner[0] = cnode->Owner();
        lm_owner[1] = cnode->Owner();

        // calculate the tangential rhs
        lm_t[0] = 0.0;
        lm_t[1] = 0.0;
        for (int j=0; j<Dim(); ++j)
        {
          lm_t[0] -= cnode->CoData().txi()[j]  * cnode->MoData().lm()[j];    // already negative rhs!!!
          lm_t[1] -= cnode->CoData().teta()[j] * cnode->MoData().lm()[j];   // already negative rhs!!!
        }
        LINALG::Assemble(tangrhs, lm_t, lm_gid, lm_owner);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 |  Assemble matrix LinStick with tangential+D+M derivatives  mgit 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinStick(LINALG::SparseMatrix& linstickLMglobal,
                                          LINALG::SparseMatrix& linstickDISglobal,
                                          Epetra_Vector& linstickRHSglobal)
{
  // FIXGIT: Assemble LinStick is containing a matrix for the de-
  // rivatives of the Lagrange multipliers. This is according to Hüeber.
  // Because of worse convergence, this is not implemented, but the
  // code is commented after the algorithm.

  // get out of here if not participating in interface
  if (!lComm())
    return;

  // create map of stick nodes
  Teuchos::RCP<Epetra_Map> sticknodes = LINALG::SplitMap(*activenodes_,*slipnodes_);
  Teuchos::RCP<Epetra_Map> stickt = LINALG::SplitMap(*activet_,*slipt_);

  // nothing to do if no stick nodes
  if (sticknodes->NumMyElements()==0)
    return;

  // information from interface contact parameter list
  INPAR::CONTACT::FrictionType ftype
  = DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(IParams(),"FRICTION");
  double frcoeff_in = IParams().get<double>("FRCOEFF"); // the friction coefficient from the input
  bool gp_slip      = DRT::INPUT::IntegralValue<int>(IParams(),"GP_SLIP_INCR");
  bool frilessfirst = DRT::INPUT::IntegralValue<int>(IParams(),"FRLESS_FIRST");

  double frcoeff=0.; // the friction coefficient actually used
  bool consistent = false;

#if defined(CONSISTENTSTICK) && defined(CONSISTENTSLIP)
  dserror("It's not reasonable to activate both, the consistent stick and slip branch, "
      "because both together will lead again to an inconsistent formulation!");
#endif
#ifdef CONSISTENTSTICK
  consistent = true;
#endif

  // loop over all stick nodes of the interface
  for (int i=0;i<sticknodes->NumMyElements();++i)
  {
    int gid = sticknodes->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleLinStick: Node ownership inconsistency!");

    // get friction coefficient for this node
    // in case of TSI, the nodal temperature influences the local friction coefficient
    frcoeff = cnode->FrCoeff(frcoeff_in);

    double cn_input   = GetCnRef()[GetCnRef().Map().LID(cnode->Id())];
    double ct_input   = GetCtRef()[GetCtRef().Map().LID(cnode->Id())];

    // more information from node
    double* n    = cnode->MoData().n();
    double* txi  = cnode->CoData().txi();
    double* teta = cnode->CoData().teta();
    double* z    = cnode->MoData().lm();

    // evaluation of specific components of entries to assemble
    double znor = 0.0;
    double ztxi = 0.0;
    double zteta = 0.0;

    for (int j=0;j<Dim();j++)
    {
      znor  += n[j]*z[j];
      ztxi  += txi[j]*z[j];
      zteta += teta[j]*z[j];
    }

    // prepare assembly, get information from node
    std::vector<GEN::pairedvector<int,double> > dnmap = cnode->CoData().GetDerivN();
    std::vector<GEN::pairedvector<int,double> > dtximap = cnode->CoData().GetDerivTxi();
    std::vector<GEN::pairedvector<int,double> > dtetamap = cnode->CoData().GetDerivTeta();
    std::map<int,double> dscmap;

    // iterator for maps
    std::map<int,double>::iterator colcurr;
    GEN::pairedvector<int,double>::iterator _colcurr;

    // row number of entries
    std::vector<int> row (Dim()-1);
    if (Dim()==2)
    {
      row[0] = stickt->GID(i);
    }
    else if (Dim()==3)
    {
      row[0] = stickt->GID(2*i);
      row[1] = stickt->GID(2*i)+1;
    }
    else
      dserror("ERROR: AssemblelinStick: Dimension not correct");

    //****************************************************************
    // CONSISTENT TREATMENT OF CASE FRCOEFF=0 (FRICTIONLESS)
    //****************************************************************
    // popp 08/2012
    //
    // There is a problem with our frictional nonlinear complementarity
    // function when applied to the limit case frcoeff=0 (frictionless).
    // In this case, the simple frictionless sliding condition should
    // be consistently recovered, which unfortunately is not the case.
    // This fact is well-known (see PhD thesis S. Hüeber) and now
    // taken care of by a special treatment as can be seen below
    //
    // seitz 11/2013
    // In case of a consistent treatment of the stick condition (i.e. using
    // the full linearization as in Diss. Gitterle eq.4.21) the case of
    // a vanishing frictional bound (frcoeff*znor==0) needs to be treated
    // as frictionless contact as well.
    // The FRLESS_FIRST option might end up here, depending on if nodes newly
    // in contact are initialized as stick or slip
    //
    //****************************************************************
    if (frcoeff*znor==0.0
        || (cnode->CoData().ActiveOld()==false
            && frilessfirst))
    {
      //**************************************************************
      // calculation of matrix entries of linearized slip condition
      //**************************************************************

      // 1) Entries from differentiation with respect to LM
      /******************************************************************/

      // loop over the dimension
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        int col = cnode->Dofs()[dim];
        double valtxi = txi[dim];

        double valteta = 0.0;
        if (Dim()==3) valteta = teta[dim];

        if (constr_direction_==INPAR::CONTACT::constr_xyz)
        {
          for (int j=0; j<Dim(); j++)
          {
            if (abs(valtxi*txi[j])>1.e-12)
              linstickLMglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
            if (Dim()==3)
              if (abs(valteta*teta[j])>1.e-12)
                linstickLMglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
          }
        }
        else
        {
          if (abs(valtxi)>1.0e-12) linstickLMglobal.Assemble(valtxi,row[0],col);
          if (Dim()==3)
            if (abs(valteta)>1.0e-12) linstickLMglobal.Assemble(valteta,row[1],col);
        }
      }

      // 2) Entries on right hand side
      /******************************************************************/
      if (constr_direction_==INPAR::CONTACT::constr_xyz)
      {
        Epetra_SerialDenseVector rhsnode(Dim());
        std::vector<int> lm(Dim());
        std::vector<int> lmowner(Dim());
        for (int j=0; j<Dim(); j++)
        {
          lm[j]=cnode->Dofs()[j];
          lmowner[j]=cnode->Owner();
          rhsnode[j] -= ztxi*txi[j];
          if (Dim()==3)
            rhsnode[j] -= zteta*teta[j];
        }
        LINALG::Assemble(linstickRHSglobal,rhsnode,lm,lmowner);
      }
      else
      {
        Epetra_SerialDenseVector rhsnode(Dim()-1);
        std::vector<int> lm(Dim()-1);
        std::vector<int> lmowner(Dim()-1);

        lm[0] = cnode->Dofs()[1];
        lmowner[0] = cnode->Owner();

        rhsnode[0] = -ztxi;   // already negative rhs!!!

        if(Dim()==3)
        {
          rhsnode[1] = -zteta;    // already negative rhs!!!

          lm[1] = cnode->Dofs()[2];
          lmowner[1] = cnode->Owner();
        }
        LINALG::Assemble(linstickRHSglobal,rhsnode,lm,lmowner);
      }

      // 3) Entries from differentiation with respect to displacements
      /******************************************************************/

      // loop over dimensions
      for (int j=0;j<Dim();++j)
      {
        // loop over all entries of the current derivative map (txi)
        for (_colcurr=dtximap[j].begin();_colcurr!=dtximap[j].end();++_colcurr)
        {
          int col = _colcurr->first;
          double val = (_colcurr->second)*z[j];

          // do not assemble zeros into matrix
          if (constr_direction_==INPAR::CONTACT::constr_xyz)
          {
            for (int j=0; j<Dim(); j++)
              if (abs(val*txi[j])>1.e-12)
                linstickDISglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
          }
          else
            if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[0],col);
        }

        if (Dim()==3)
        {
          // loop over all entries of the current derivative map (teta)
          for (_colcurr=dtetamap[j].begin();_colcurr!=dtetamap[j].end();++_colcurr)
          {
            int col = _colcurr->first;
            double val = (_colcurr->second)*z[j];

            // do not assemble zeros into s matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(val*teta[j])>1.e-12)
                  linstickDISglobal.Assemble(val*teta[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[1],col);
          }
        }
      }
    }
    else if (consistent && ftype == INPAR::CONTACT::friction_coulomb)
    {
      std::map<int,double>& dgmap = cnode->CoData().GetDerivG();

      // check for Dimension of derivative maps
      for (int j=0;j<Dim()-1;++j)
        if ((int)dnmap[j].size() != (int)dnmap[j+1].size())
          dserror("ERROR: AssembleLinStick: Column dim. of nodal DerivN-map is inconsistent!");

      for (int j=0;j<Dim()-1;++j)
        if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
          dserror("ERROR: AssembleLinStick: Column dim. of nodal DerivTxi-map is inconsistent!");

      if (Dim()==3)
      {
        for (int j=0;j<Dim()-1;++j)
          if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
            dserror("ERROR: AssembleLinStick: Column dim. of nodal DerivTeta-map is inconsistent!");
      }

      double cn = cn_input;
      double ct = ct_input;

      double& wgap = cnode->CoData().Getg();

      // evaluation of specific components of entries to assemble
      double jumptxi = 0;
      double jumpteta = 0;

      double* jump = cnode->FriData().jump();

      // choose slip increment
      if (gp_slip)
      {
        jumptxi=cnode->FriData().jump_var()[0];

        if (Dim()==3)
          jumpteta=cnode->FriData().jump_var()[1];
      }
      else
      {
        // more information from node
        for (int i=0;i<Dim();i++)
        {
          jumptxi += txi[i]*jump[i];
          jumpteta += teta[i]*jump[i];
        }
      }

      // check for dimensions
      if(Dim()==2 and (jumpteta != 0.0))
        dserror ("ERROR: AssembleLinStick: jumpteta must be zero in 2D");

      //**************************************************************
      // calculation of matrix entries of linearized stick condition
      //**************************************************************

      // 1) Entries from differentiation with respect to LM
      /******************************************************************/

      // loop over the dimension
      for (int dim=0;dim<cnode->NumDof();++dim)
      {
        double valtxi = 0.0;
        double valteta = 0.0;
        int col = cnode->Dofs()[dim];
        valtxi = -frcoeff * ct * jumptxi * n[dim];

        if (Dim()==3)
        {
          valteta = -frcoeff * ct * jumpteta * n[dim];
        }

        // do not assemble zeros into matrix
        if (constr_direction_==INPAR::CONTACT::constr_xyz)
        {
          for (int j=0; j<Dim(); j++)
          {
            if (abs(valtxi*txi[j])>1.e-12)
              linstickLMglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
            if (Dim()==3)
              if (abs(valteta*teta[j])>1.e-12)
                linstickLMglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
          }
        }
        else
        {
          if (abs(valtxi)>1.0e-12) linstickLMglobal.Assemble(valtxi,row[0],col);
          if (Dim()==3)
            if (abs(valteta)>1.0e-12) linstickLMglobal.Assemble(valteta,row[1],col);
        }
      }

      // Entries on right hand side ****************************
      if (constr_direction_==INPAR::CONTACT::constr_xyz)
      {
        Epetra_SerialDenseVector rhsnode(Dim());
        std::vector<int> lm(Dim());
        std::vector<int> lmowner(Dim());

        for (int j=0; j<Dim(); j++)
        {
          lm[j]=cnode->Dofs()[j];
          lmowner[j]=cnode->Owner();
          rhsnode(j) += frcoeff*(znor-cn*wgap)*ct*jumptxi*txi[j];
          if (Dim()==3)
            rhsnode(j) += frcoeff*(znor-cn*wgap)*ct*jumpteta*teta[j];
        }
        LINALG::Assemble(linstickRHSglobal,rhsnode,lm,lmowner);
      }
      else
      {
        Epetra_SerialDenseVector rhsnode(Dim()-1);
        std::vector<int> lm(Dim()-1);
        std::vector<int> lmowner(Dim()-1);
        rhsnode(0) = frcoeff * (znor - cn * wgap) * ct * jumptxi;

        lm[0] = cnode->Dofs()[1];
        lmowner[0] = cnode->Owner();
        if (Dim()==3)
        {
          rhsnode(1) = frcoeff * (znor - cn * wgap) * ct * jumpteta;
          lm[1] = cnode->Dofs()[2];
          lmowner[1] = cnode->Owner();
        }
        LINALG::Assemble(linstickRHSglobal,rhsnode,lm,lmowner);
      }


      // 3) Entries from differentiation with respect to displacements
      /******************************************************************/

      if (gp_slip)
      {
        std::vector<std::map<int,double> > derivjump_ = cnode->FriData().GetDerivVarJump();

        //txi
        for (colcurr=derivjump_[0].begin();colcurr!=derivjump_[0].end();++colcurr)
        {
          int col = colcurr->first;
          double valtxi = - frcoeff * (znor - cn * wgap) * ct * (colcurr->second);

          if (constr_direction_==INPAR::CONTACT::constr_xyz)
          {
            for (int j=0; j<Dim(); j++)
              if (abs(valtxi*txi[j])>1.e-12)
                linstickDISglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
          }
          else
            if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);
        }
        //teta
        for (colcurr=derivjump_[1].begin();colcurr!=derivjump_[1].end();++colcurr)
        {
          int col = colcurr->first;
          double valteta = - frcoeff * (znor - cn * wgap) * ct * (colcurr->second);

          if (constr_direction_==INPAR::CONTACT::constr_xyz)
          {
            for (int j=0; j<Dim(); j++)
              if (abs(valteta*teta[j])>1.e-12)
                linstickDISglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
          }
          else
            if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
        }

        // ... old slip
        // get linearization of jump vector
        std::vector<std::map<int,double> > derivjump = cnode->FriData().GetDerivJump();

        // loop over dimensions
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // linearization of normal direction *****************************************
          // loop over all entries of the current derivative map
          for (_colcurr=dnmap[dim].begin();_colcurr!=dnmap[dim].end();++_colcurr)
          {
            int col = _colcurr->first;
            double valtxi=0.0;
            valtxi = - frcoeff * z[dim] * _colcurr->second * ct * jumptxi;

            // do not assemble zeros into matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(valtxi*txi[j])>1.e-12)
                  linstickDISglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);

            if (Dim()==3)
            {
              double valteta=0.0;
              valteta = - frcoeff * z[dim] * _colcurr->second * ct * jumpteta;

              // do not assemble zeros into matrix
              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim(); j++)
                  if (abs(valteta*teta[j])>1.e-12)
                    linstickDISglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
              }
              else
                if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
            }
          }
        } // loop over all dimensions
      }
      else // std slip
      {
        // ... old slip
        // get linearization of jump vector
        std::vector<std::map<int,double> > derivjump = cnode->FriData().GetDerivJump();

        // loop over dimensions
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          // loop over all entries of the current derivative map (jump)
          for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
          {
            int col = colcurr->first;

            double valtxi=0.0;
            valtxi = - frcoeff * (znor - cn * wgap) * ct * txi[dim] * (colcurr->second);

            // do not assemble zeros into matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(valtxi*txi[j])>1.e-12)
                  linstickDISglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);

            if (Dim()==3)
            {
              double valteta=0.0;
              valteta = - frcoeff * (znor - cn * wgap) * ct * teta[dim] * (colcurr->second);

              // do not assemble zeros into matrix
              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim(); j++)
                  if (abs(valteta*teta[j])>1.e-12)
                    linstickDISglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
              }
              else
                if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
            }
          }

          // linearization first tangential direction *********************************
          // loop over all entries of the current derivative map (txi)
          for (_colcurr=dtximap[dim].begin();_colcurr!=dtximap[dim].end();++_colcurr)
          {
            int col = _colcurr->first;
            double valtxi=0.0;
            valtxi = - frcoeff*(znor-cn*wgap) * ct * jump[dim] * _colcurr->second;

            // do not assemble zeros into matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(valtxi*txi[j])>1.e-12)
                  linstickDISglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);
          }
          // linearization second tangential direction *********************************
          if (Dim()==3)
          {
            // loop over all entries of the current derivative map (teta)
            for (_colcurr=dtetamap[dim].begin();_colcurr!=dtetamap[dim].end();++_colcurr)
            {
              int col = _colcurr->first;
              double valteta=0.0;
              valteta = - frcoeff * (znor-cn*wgap) * ct * jump[dim] * _colcurr->second;

              // do not assemble zeros into matrix
              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim(); j++)
                  if (abs(valteta*teta[j])>1.e-12)
                    linstickDISglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
              }
              else
                if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
            }
          }

          // linearization of normal direction *****************************************
          // loop over all entries of the current derivative map
          for (_colcurr=dnmap[dim].begin();_colcurr!=dnmap[dim].end();++_colcurr)
          {
            int col = _colcurr->first;
            double valtxi=0.0;
            valtxi = - frcoeff * z[dim] * _colcurr->second * ct * jumptxi;

            // do not assemble zeros into matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(valtxi*txi[j])>1.e-12)
                  linstickDISglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);

            if (Dim()==3)
            {
              double valteta=0.0;
              valteta = - frcoeff * z[dim] * _colcurr->second * ct * jumpteta;

              // do not assemble zeros into matrix
              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim(); j++)
                  if (abs(valteta*teta[j])>1.e-12)
                    linstickDISglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
              }
              else
                if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
            }
          }
        } // loop over all dimensions
      }

        // linearization of weighted gap**********************************************
        // loop over all entries of the current derivative map fixme
        for (colcurr=dgmap.begin();colcurr!=dgmap.end();++colcurr)
        {
          int col = colcurr->first;
          double valtxi=0.0;
          valtxi = frcoeff * colcurr->second * ct * cn * jumptxi;

          // do not assemble zeros into matrix
          if (constr_direction_==INPAR::CONTACT::constr_xyz)
          {
            for (int j=0; j<Dim(); j++)
              if (abs(valtxi*txi[j])>1.e-12)
                linstickDISglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
          }
          else
            if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);

          if (Dim()==3)
          {
            double valteta=0.0;
            valteta = frcoeff * colcurr->second * ct * cn * jumpteta;

            // do not assemble zeros into matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(valteta*teta[j])>1.e-12)
                  linstickDISglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
          }
        }
      }
      else // not consistent stick
      {
        for (int j=0;j<Dim()-1;++j)
          if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
            dserror("ERROR: AssembleLinStick: Column dim. of nodal DerivTxi-map is inconsistent!");

        if (Dim()==3)
        {
          for (int j=0;j<Dim()-1;++j)
            if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
              dserror("ERROR: AssembleLinStick: Column dim. of nodal DerivTeta-map is inconsistent!");
        }

        // evaluation of specific components of entries to assemble
        double jumptxi=0;
        double jumpteta=0;

        // more information from node
        double* jump = cnode->FriData().jump();

        if (gp_slip)
        {
          jumptxi=cnode->FriData().jump_var()[0];

          if (Dim()==3)
            jumpteta=cnode->FriData().jump_var()[1];
        }
        else
        {
          for (int i=0;i<Dim();i++)
          {
            jumptxi += txi[i]*jump[i];
            jumpteta += teta[i]*jump[i];
          }
        }

        // check for dimensions
        if(Dim()==2 and (jumpteta != 0.0))
          dserror ("ERROR: AssembleLinStick: jumpteta must be zero in 2D");

        // Entries on right hand side
        /************************************************ (-utxi, -uteta) ***/
        if (constr_direction_==INPAR::CONTACT::constr_xyz)
        {
          Epetra_SerialDenseVector rhsnode(Dim());
          std::vector<int> lm(Dim());
          std::vector<int> lmowner(Dim());
          for (int j=0; j<Dim(); j++)
          {
            lm[j]=cnode->Dofs()[j];
            lmowner[j]=cnode->Owner();
            rhsnode(j) -= jumptxi*txi[j];
            if (Dim()==3)
              rhsnode(j) -= jumpteta*teta[j];
          }
          LINALG::Assemble(linstickRHSglobal,rhsnode,lm,lmowner);
        }
        else
        {
          Epetra_SerialDenseVector rhsnode(Dim()-1);
          std::vector<int> lm(Dim()-1);
          std::vector<int> lmowner(Dim()-1);

          // modification to stabilize the convergence of the lagrange multiplier incr (hiermeier 08/13)
          if (abs(jumptxi)<1e-15)
            rhsnode(0) = 0.0;
          else
            rhsnode(0) = -jumptxi;

          lm[0] = cnode->Dofs()[1];
          lmowner[0] = cnode->Owner();

          if (Dim()==3)
          {
            // modification to stabilize the convergence of the lagrange multiplier incr (hiermeier 08/13)
            if (abs(jumpteta)<1e-15)
              rhsnode(1) = 0.0;
            else
              rhsnode(1) = -jumpteta;

            lm[1] = cnode->Dofs()[2];
            lmowner[1] = cnode->Owner();
          }
          LINALG::Assemble(linstickRHSglobal,rhsnode,lm,lmowner);
        }


        // The routine "ApplyDirichlet" in SaddlepointSolve can only set ones on the diagonal
        // if there has already been a diagonal entry in the sparse matrix
        for (int j=0; j<Dim(); j++)
        {
          if (cnode->DbcDofs()[j]==true)
          {
            linstickLMglobal.Assemble(1.e-12,cnode->Dofs()[j],cnode->Dofs()[j]);
          }
        }

        // Entries from differentiation with respect to displacements
        /*** 1 ************************************** tangent.deriv(jump) ***/
        if (gp_slip)
        {
          std::map<int,double> derivjump1 = cnode->FriData().GetDerivVarJump()[0];
          std::map<int,double> derivjump2 = cnode->FriData().GetDerivVarJump()[1];

          for (colcurr=derivjump1.begin();colcurr!=derivjump1.end();++colcurr)
          {
            int col = colcurr->first;
            double valtxi = colcurr->second;

            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(valtxi*txi[j])>1.e-12)
                  linstickDISglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);
          }

          if(Dim()==3)
          {
            for (colcurr=derivjump2.begin();colcurr!=derivjump2.end();++colcurr)
            {
              int col = colcurr->first;
              double valteta = colcurr->second;

              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim(); j++)
                  if (abs(valteta*teta[j])>1.e-12)
                    linstickDISglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
              }
              else
                if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
            }
          }
        }
        else // std. slip
        {
          // get linearization of jump vector
          std::vector<std::map<int,double> > derivjump = cnode->FriData().GetDerivJump();

          if (derivjump.size()<1)
            dserror ("AssembleLinStick: Derivative of jump is not exiting!");

          // loop over dimensions
          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            // loop over all entries of the current derivative map (jump)
            for (colcurr=derivjump[dim].begin();colcurr!=derivjump[dim].end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi = txi[dim]*colcurr->second;

              // do not assemble zeros into matrix
              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim(); j++)
                  if (abs(valtxi*txi[j])>1.e-12)
                    linstickDISglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
              }
              else
                if (abs(valtxi)>1.0e-12) linstickDISglobal.Assemble(valtxi,row[0],col);

              if(Dim()==3)
              {
                double valteta = teta[dim]*colcurr->second;

                if (constr_direction_==INPAR::CONTACT::constr_xyz)
                {
                  for (int j=0; j<Dim(); j++)
                    if (abs(valteta*teta[j])>1.e-12)
                      linstickDISglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
                }
                else
                  if (abs(valteta)>1.0e-12) linstickDISglobal.Assemble(valteta,row[1],col);
              }
            }
          }

          /*** 2 ************************************** deriv(tangent).jump ***/
          // loop over dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (_colcurr=dtximap[j].begin();_colcurr!=dtximap[j].end();++_colcurr)
            {
              int col = _colcurr->first;
              double val = jump[j]*_colcurr->second;

              // do not assemble zeros into s matrix
              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim(); j++)
                  if (abs(val*txi[j])>1.e-12)
                    linstickDISglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
              }
              else
                if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[0],col);
            }

            if(Dim()==3)
            {
              // loop over all entries of the current derivative map (teta)
              for (_colcurr=dtetamap[j].begin();_colcurr!=dtetamap[j].end();++_colcurr)
              {
                int col = _colcurr->first;
                double val = jump[j]*_colcurr->second;

                // do not assemble zeros into matrix
                if (constr_direction_==INPAR::CONTACT::constr_xyz)
                {
                  for (int j=0; j<Dim(); j++)
                    if (abs(val*teta[j])>1.e-12)
                      linstickDISglobal.Assemble(val*teta[j],cnode->Dofs()[j],col);
                }
                else
                  if (abs(val)>1.0e-12) linstickDISglobal.Assemble(val,row[1],col);
              }
            }
          }
        }
    }
  }
  return;
}

/*---------------------------------------------------------------------*
 | Assemble matrix LinSlip with tangential+D+M derivatives  mgit 02/09 |
 *---------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleLinSlip(LINALG::SparseMatrix& linslipLMglobal,
    LINALG::SparseMatrix& linslipDISglobal,
    Epetra_Vector& linslipRHSglobal)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no slip nodes
  if (slipnodes_->NumMyElements()==0)
    return;

  // information from interface contact parameter list
  INPAR::CONTACT::FrictionType ftype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(IParams(),"FRICTION");
  double frbound    = IParams().get<double>("FRBOUND");
  double frcoeff_in = IParams().get<double>("FRCOEFF"); // the friction coefficient from the input
  bool gp_slip      = DRT::INPUT::IntegralValue<int>(IParams(),"GP_SLIP_INCR");
  bool frilessfirst = DRT::INPUT::IntegralValue<int>(IParams(),"FRLESS_FIRST");

  // the friction coefficient adapted by every node (eg depending on the local temperature)
  double frcoeff=0.;

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
  // Coulomb Friction
  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
  if (ftype == INPAR::CONTACT::friction_coulomb)
  {
    // loop over all slip nodes of the interface
    for (int i=0;i<slipnodes_->NumMyElements();++i)
    {
      int gid = slipnodes_->GID(i);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      FriNode* cnode = dynamic_cast<FriNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleLinSlip: Node ownership inconsistency!");

      // get friction coefficient for this node
      // in case of TSI, the nodal temperature influences the local friction coefficient
      frcoeff = cnode->FrCoeff(frcoeff_in);

      double cn_input = GetCnRef()[GetCnRef().Map().LID(cnode->Id())];
      double ct_input = GetCtRef()[GetCtRef().Map().LID(cnode->Id())];

      // prepare assembly, get information from node
      std::vector<GEN::pairedvector<int,double> > dnmap = cnode->CoData().GetDerivN();
      std::vector<GEN::pairedvector<int,double> > dtximap = cnode->CoData().GetDerivTxi();
      std::vector<GEN::pairedvector<int,double> > dtetamap = cnode->CoData().GetDerivTeta();

      double cn = cn_input;
      double ct = ct_input;

      // check for Dimension of derivative maps
      for (int j=0;j<Dim()-1;++j)
        if ((int)dnmap[j].size() != (int)dnmap[j+1].size())
          dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivTxi-map is inconsistent!");

      for (int j=0;j<Dim()-1;++j)
        if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
          dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivTxi-map is inconsistent!");

      if (Dim()==3)
      {
        for (int j=0;j<Dim()-1;++j)
          if ((int)dtximap[j].size() != (int)dtximap[j+1].size())
            dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivTeta-map is inconsistent!");
      }

      // more information from node
      double* n    = cnode->MoData().n();
      double* txi  = cnode->CoData().txi();
      double* teta = cnode->CoData().teta();
      double* z    = cnode->MoData().lm();
      double wgap  = cnode->CoData().Getg();

      // iterator for maps
      GEN::pairedvector<int,double>::iterator colcurr;
      std::map<int,double>::iterator _colcurr;

      // row number of entries
      std::vector<int> row (Dim()-1);
      if (Dim()==2)
      {
        row[0] = slipt_->GID(i);
      }
      else if (Dim()==3)
      {
        row[0] = slipt_->GID(2*i);
        row[1] = slipt_->GID(2*i)+1;
      }
      else
        dserror("ERROR: AssemblelinSlip: Dimension not correct");

      // evaluation of specific components of entries to assemble
      double znor = 0.0;
      double ztxi = 0.0;
      double zteta = 0.0;
      double jumptxi = 0.0;
      double jumpteta = 0.0;
      double euclidean = 0.0;
      double* jump = cnode->FriData().jump();

      if (gp_slip)
      {
        jumptxi=cnode->FriData().jump_var()[0];

        if (Dim()==3)
          jumpteta=cnode->FriData().jump_var()[1];

        for (int i=0;i<Dim();i++)
        {
          znor  += n[i]*z[i];
          ztxi  += txi[i]*z[i];
          zteta += teta[i]*z[i];
        }
      }
      else // std. slip
      {
        for (int i=0;i<Dim();i++)
        {
          znor += n[i]*z[i];
          ztxi += txi[i]*z[i];
          zteta += teta[i]*z[i];
          jumptxi += txi[i]*jump[i];
          jumpteta += teta[i]*jump[i];
        }
      }

      // evaluate euclidean norm ||vec(zt)+ct*vec(jumpt)||
      std::vector<double> sum1 (Dim()-1,0);
      sum1[0] = ztxi+ct*jumptxi;
      if (Dim()==3) sum1[1] = zteta+ct*jumpteta;
      if (Dim()==2) euclidean = abs(sum1[0]);
      if (Dim()==3) euclidean = sqrt(sum1[0]*sum1[0]+sum1[1]*sum1[1]);

      // check of dimensions
      if(Dim()==2 and (zteta != 0.0 or jumpteta != 0.0))
        dserror ("ERROR: AssemblelinSlip: zteta and jumpteta must be zero in 2D");

      //****************************************************************
      // CONSISTENT TREATMENT OF CASE FRCOEFF=0 (FRICTIONLESS)
      //****************************************************************
      // popp 08/2012
      //
      // There is a problem with our frictional nonlinear complementarity
      // function when applied to the limit case frcoeff=0 (frictionless).
      // In this case, the simple frictionless sliding condition should
      // be consistently recovered, which unfortunately is not the case.
      // This fact is well-known (see PhD thesis S. Hüeber) and now
      // taken care of by a special treatment as can be seen below
      //
      // seitz11/13
      // Like that FRLESS_FIRST results in a frictionless constraint.
      // Nodes with a vanishing slip and tangential Lagrange multiplier
      // are treated as frictionless as well.
      //
      //****************************************************************
      if (frcoeff==0.0
          || (cnode->CoData().ActiveOld()==false
              && frilessfirst)
          || euclidean==0.0)
      {
        //**************************************************************
        // calculation of matrix entries of linearized slip condition
        //**************************************************************

        // 1) Entries from differentiation with respect to LM
        /******************************************************************/

        // loop over the dimension
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cnode->Dofs()[dim];
          double valtxi = txi[dim];

          double valteta = 0.0;
          if (Dim()==3) valteta = teta[dim];

          // do not assemble zeros into matrix
          if (constr_direction_==INPAR::CONTACT::constr_xyz)
          {
            for (int j=0; j<Dim(); j++)
            {
              if (abs(valtxi*txi[j])>1.e-12)
                linslipLMglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
              if (Dim()==3)
                if (abs(valteta*teta[j])>1.e-12)
                  linslipLMglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
            }
          }
          else
          {
            if (abs(valtxi)>1.0e-12) linslipLMglobal.Assemble(valtxi,row[0],col);
            if (Dim()==3)
              if (abs(valteta)>1.0e-12) linslipLMglobal.Assemble(valteta,row[1],col);
          }
        }

        // 2) Entries on right hand side
        /******************************************************************/
        if (constr_direction_==INPAR::CONTACT::constr_xyz)
        {
          Epetra_SerialDenseVector rhsnode(Dim());
          std::vector<int> lm(Dim());
          std::vector<int> lmowner(Dim());
          for (int j=0; j<Dim(); j++)
          {
            lm[j]=cnode->Dofs()[j];
            lmowner[j]=cnode->Owner();
            rhsnode[j] -= ztxi*txi[j];
            if (Dim()==3)
              rhsnode[j] -= zteta*teta[j];
          }
          LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);
        }
        else
        {
          Epetra_SerialDenseVector rhsnode(Dim()-1);
          std::vector<int> lm(Dim()-1);
          std::vector<int> lmowner(Dim()-1);

          lm[0] = cnode->Dofs()[1];
          lmowner[0] = cnode->Owner();

          rhsnode[0] = -ztxi;   // already negative rhs!!!

          if(Dim()==3)
          {
            rhsnode[1] = -zteta;    // already negative rhs!!!

            lm[1] = cnode->Dofs()[2];
            lmowner[1] = cnode->Owner();
          }
          LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);
        }

        // 3) Entries from differentiation with respect to displacements
        /******************************************************************/

        // loop over dimensions
        for (int j=0;j<Dim();++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (colcurr->second)*z[j];

            // do not assemble zeros into matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(val*txi[j])>1.e-12)
                  linslipDISglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
          }

          if (Dim()==3)
          {
            // loop over all entries of the current derivative map (teta)
            for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double val = (colcurr->second)*z[j];

              // do not assemble zeros into s matrix
              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim(); j++)
                  if (abs(val*teta[j])>1.e-12)
                    linslipDISglobal.Assemble(val*teta[j],cnode->Dofs()[j],col);
              }
              else
                if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[1],col);
            }
          }
        }
      }

      //****************************************************************
      // STANDARD TREATMENT OF CASE FRCOEFF!=0 (FRICTIONAL)
      //****************************************************************
      else
      {
        //**************************************************************
        // calculation of matrix entries of linearized slip condition
        //**************************************************************

        // 1) Entries from differentiation with respect to LM
        /******************************************************************/

        // loop over the dimension
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          double valtxi = 0.0;
          int col = cnode->Dofs()[dim];

          double valtxi0 = euclidean*txi[dim];
          double valtxi1 = ((ztxi+ct*jumptxi)/euclidean*ztxi)*txi[dim];
          double valtxi3 = (zteta+ct*jumpteta)/euclidean*ztxi*teta[dim];

#ifdef CONSISTENTSLIP
          valtxi0 = valtxi0 / (znor - cn * wgap);
          valtxi1 = valtxi1 / (znor - cn * wgap);
          valtxi3 = valtxi3 / (znor - cn * wgap);

          // Additional term
          valtxi0 -= euclidean * ztxi / pow(znor -cn * wgap,2.0) * n[dim];

          double valtxi2 = -frcoeff * txi[dim];
#else
          double valtxi2 = -frcoeff*(znor-cn*wgap)*txi[dim]-frcoeff*(ztxi+ct*jumptxi)*n[dim];
#endif
          valtxi = valtxi0 + valtxi1 + valtxi2 + valtxi3;

          double valteta = 0.0;
          if (Dim()==3)
          {
            double valteta0 = euclidean*teta[dim];
            double valteta1 = ((ztxi+ct*jumptxi)/euclidean*zteta)*txi[dim];
            double valteta3 = (zteta+ct*jumpteta)/euclidean*zteta*teta[dim];
#ifdef CONSISTENTSLIP
            valteta0 = valteta0 / (znor - cn * wgap);
            valteta1 = valteta1 / (znor - cn * wgap);
            valteta3 = valteta3 / (znor - cn * wgap);

            // Additional term
            valteta0 -= euclidean * zteta / pow(znor -cn * wgap,2.0) * n[dim];

            double valteta2 = -frcoeff * teta[dim];
#else
            double valteta2 = -frcoeff*(znor-cn*wgap)*teta[dim]-frcoeff*(zteta+ct*jumpteta)*n[dim];
#endif
            valteta = valteta0 + valteta1 + valteta2 + valteta3;
          }

          // do not assemble zeros into matrix
          if (constr_direction_==INPAR::CONTACT::constr_xyz)
          {
            for (int j=0; j<Dim(); j++)
            {
              if (abs(valtxi*txi[j])>1.e-12)
                linslipLMglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
              if (Dim()==3)
                if (abs(valteta*teta[j])>1.e-12)
                  linslipLMglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
            }
          }
          else
          {
            if (abs(valtxi)>1.0e-12) linslipLMglobal.Assemble(valtxi,row[0],col);
            if (Dim()==3)
              if (abs(valteta)>1.0e-12) linslipLMglobal.Assemble(valteta,row[1],col);
          }
        }

        // 2) Entries on right hand side
        /******************************************************************/
        if (constr_direction_==INPAR::CONTACT::constr_xyz)
        {
          Epetra_SerialDenseVector rhsnode(Dim());
          std::vector<int> lm(Dim());
          std::vector<int> lmowner(Dim());

#ifdef CONSISTENTSLIP
          double valuetxi1 = -(euclidean)*ztxi / (znor - cn * wgap) + frcoeff*(ztxi+ct*jumptxi);
#else
          double valuetxi1 = -(euclidean)*ztxi+(frcoeff*(znor-cn*wgap))*(ztxi+ct*jumptxi);
#endif

          for (int j=0; j<Dim(); j++)
          {
            lm[j]=cnode->Dofs()[j];
            lmowner[j]=cnode->Owner();
            rhsnode(j) +=valuetxi1*txi[j];
          }

          if(Dim()==3)
          {
#ifdef CONSISTENTSLIP
            double valueteta1 = -(euclidean)*zteta / (znor - cn * wgap) + frcoeff * (zteta + ct * jumpteta);
#else
            double valueteta1 = -(euclidean)*zteta+(frcoeff*(znor-cn*wgap))*(zteta+ct*jumpteta);
#endif

            for (int j=0; j<Dim(); j++)
              rhsnode(j) += valueteta1*teta[j];
          }
          LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);
        }
        else
        {
          Epetra_SerialDenseVector rhsnode(Dim()-1);
          std::vector<int> lm(Dim()-1);
          std::vector<int> lmowner(Dim()-1);
#ifdef CONSISTENTSLIP
          double valuetxi1 = -(euclidean)*ztxi / (znor - cn * wgap) + frcoeff*(ztxi+ct*jumptxi);
#else
          double valuetxi1 = -(euclidean)*ztxi+(frcoeff*(znor-cn*wgap))*(ztxi+ct*jumptxi);
#endif
          rhsnode(0) = valuetxi1;
          lm[0] = cnode->Dofs()[1];
          lmowner[0] = cnode->Owner();

          if(Dim()==3)
          {
#ifdef CONSISTENTSLIP
            double valueteta1 = -(euclidean)*zteta / (znor - cn * wgap) + frcoeff * (zteta + ct * jumpteta);
#else
            double valueteta1 = -(euclidean)*zteta+(frcoeff*(znor-cn*wgap))*(zteta+ct*jumpteta);
#endif
            rhsnode(1) = valueteta1;

            lm[1] = cnode->Dofs()[2];
            lmowner[1] = cnode->Owner();
          }
          LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);
        }

        // 3) Entries from differentiation with respect to displacements
        /******************************************************************/
        std::map<int,double> derivjump1, derivjump2;    // for gp slip
        std::vector<std::map<int,double> > derivjump ;  // for dm slip

        /*** 01  ********* -Deriv(euclidean).ct.tangent.deriv(u)*ztan ***/
        if (gp_slip)
        {
          derivjump1 = cnode->FriData().GetDerivVarJump()[0];
          derivjump2 = cnode->FriData().GetDerivVarJump()[1];

          for (_colcurr=derivjump1.begin();_colcurr!=derivjump1.end();++_colcurr)
          {
            int col = _colcurr->first;
            double valtxi1 = (ztxi+ct*jumptxi)/euclidean*ct*_colcurr->second*ztxi;
            double valteta1 = (ztxi+ct*jumptxi)/euclidean*ct*_colcurr->second*zteta;

            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
              {
                if (abs(valtxi1*txi[j])>1.e-12)
                  linslipDISglobal.Assemble(valtxi1*txi[j],cnode->Dofs()[j],col);
                if (abs(valteta1*teta[j])>1.e-12)
                  linslipDISglobal.Assemble(valteta1*teta[j],cnode->Dofs()[j],col);
              }
            }
            else
            {
              if (abs(valtxi1)>1.0e-12) linslipDISglobal.Assemble(valtxi1,row[0],col);
              if (abs(valteta1)>1.0e-12) linslipDISglobal.Assemble(valteta1,row[1],col);
            }
          }

          if (Dim()==3)
          {
            for (_colcurr=derivjump2.begin();_colcurr!=derivjump2.end();++_colcurr)
            {
              int col = _colcurr->first;
              double valtxi2 = (zteta+ct*jumpteta)/euclidean*ct*_colcurr->second*ztxi;
              double valteta2 = (zteta+ct*jumpteta)/euclidean*ct*_colcurr->second*zteta;

              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim(); j++)
                {
                  if (abs(valtxi2*txi[j])>1.e-12)
                    linslipDISglobal.Assemble(valtxi2*txi[j],cnode->Dofs()[j],col);
                  if (abs(valteta2*teta[j])>1.e-12)
                    linslipDISglobal.Assemble(valteta2*teta[j],cnode->Dofs()[j],col);
                }
              }
              else
              {
                if (abs(valtxi2)>1.0e-12) linslipDISglobal.Assemble(valtxi2,row[0],col);
                if (abs(valteta2)>1.0e-12) linslipDISglobal.Assemble(valteta2,row[1],col);
              }
            }
          }
        }
        else // std. slip increment
        {
          // get linearization of jump vector
          derivjump = cnode->FriData().GetDerivJump();

          // loop over dimensions
          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            // loop over all entries of the current derivative map (jump)
            for (_colcurr=derivjump[dim].begin();_colcurr!=derivjump[dim].end();++_colcurr)
            {
              int col = _colcurr->first;

              double valtxi1 = (ztxi+ct*jumptxi)/euclidean*ct*txi[dim]*_colcurr->second*ztxi;
              double valteta1 = (ztxi+ct*jumptxi)/euclidean*ct*txi[dim]*_colcurr->second*zteta;
              double valtxi2 = (zteta+ct*jumpteta)/euclidean*ct*teta[dim]*_colcurr->second*ztxi;
              double valteta2 = (zteta+ct*jumpteta)/euclidean*ct*teta[dim]*_colcurr->second*zteta;

  #ifdef CONSISTENTSLIP
              valtxi1   = valtxi1  / (znor - cn * wgap);
              valteta1  = valteta1 / (znor - cn * wgap);
              valtxi2   = valtxi2  / (znor - cn * wgap);
              valteta2  = valteta2 / (znor - cn * wgap);
  #endif

              // do not assemble zeros into matrix
              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim();j++)
                {
                  if (abs((valtxi1+valtxi2)*txi[j])>1.e-12)
                    linslipDISglobal.Assemble((valtxi1+valtxi2)*txi[j],cnode->Dofs()[j],col);
                  if (Dim()==3)
                    if (abs((valteta1+valteta2)*teta[j])>1.e-12)
                      linslipDISglobal.Assemble((valteta1+valteta2)*teta[j],cnode->Dofs()[j],col);
                }
              }
              else
              {
                if (abs(valtxi1)>1.0e-12) linslipDISglobal.Assemble(valtxi1,row[0],col);
                if (abs(valteta1)>1.0e-12) linslipDISglobal.Assemble(valteta1,row[1],col);
                if (abs(valtxi2)>1.0e-12) linslipDISglobal.Assemble(valtxi2,row[0],col);
                if (abs(valteta2)>1.0e-12) linslipDISglobal.Assemble(valteta2,row[1],col);
              }
            }

  #ifdef CONSISTENTSLIP
            /*** Additional Terms ***/
            // normal derivative
            for (colcurr=dnmap[dim].begin();colcurr!=dnmap[dim].end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi   = - euclidean * ztxi * z[dim] / pow(znor - cn * wgap, 2.0) * (colcurr->second);
              double valteta  = - euclidean * zteta * z[dim] / pow(znor - cn * wgap, 2.0) * (colcurr->second);

              // do not assemble zeros into s matrix
              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim(); j++)
                {
                  if (abs(valtxi*txi[j])>1.e-12)
                    linslipDISglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
                  if (Dim()==3)
                    if (abs(valteta*teta[j])>1.e-12)
                      linslipDISglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
                }
              }
              else
              {
                if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
                if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
              }
            }
#endif
          }
        }

#ifdef CONSISTENTSLIP
        /*** Additional Terms ***/
        // wgap derivative
        std::map<int,double>& dgmap = cnode->CoData().GetDerivG();

        for (colcurr=dgmap.begin(); colcurr!=dgmap.end(); ++colcurr)
        {
          int col = colcurr->first;
          double valtxi  = + euclidean * ztxi  / pow(znor - cn * wgap, 2.0) * cn * (colcurr->second);
          double valteta = + euclidean * zteta / pow(znor - cn * wgap, 2.0) * cn * (colcurr->second);

          //do not assemble zeros into matrix
          if (constr_direction_==INPAR::CONTACT::constr_xyz)
          {
            for (int j=0; j<Dim(); j++)
            {
              if (abs(valtxi*txi[j])>1.e-12)
                linslipDISglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
              if (Dim()==3)
                if (abs(valteta*teta[j])>1.e-12)
                  linslipDISglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
            }
          }
          else
          {
            if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
            if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
          }
        }
#endif


        /*** 02 ***************** frcoeff*znor*ct*tangent.deriv(jump) ***/
        if (gp_slip)
        {
          for (_colcurr=derivjump1.begin();_colcurr!=derivjump1.end();++_colcurr)
          {
            int col = _colcurr->first;
            double valtxi = -frcoeff*(znor-cn*wgap)*ct*_colcurr->second;
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(valtxi*txi[j])>1.e-12)
                  linslipDISglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
          }

          if (Dim()==3)
          {
            for (_colcurr=derivjump2.begin();_colcurr!=derivjump2.end();++_colcurr)
            {
              int col = _colcurr->first;
              double valteta = -frcoeff*(znor-cn*wgap)*ct*_colcurr->second;

              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim(); j++)
                  if (abs(valteta*teta[j])>1.e-12)
                    linslipDISglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
              }
              else
                if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }
          }
        }
        else
        {
          // loop over dimensions
          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            // loop over all entries of the current derivative map (jump)
            for (_colcurr=derivjump[dim].begin();_colcurr!=derivjump[dim].end();++_colcurr)
            {
              int col = _colcurr->first;

              //std::cout << "val " << colcurr->second << std::endl;
#ifdef CONSISTENTSLIP
              double valtxi = - frcoeff * ct * txi[dim] * colcurr->second;
              double valteta = - frcoeff * ct * teta[dim] * colcurr->second;
#else
              double valtxi = (-1)*(frcoeff*(znor-cn*wgap))*ct*txi[dim]*_colcurr->second;
              double valteta = (-1)*(frcoeff*(znor-cn*wgap))*ct*teta[dim]*_colcurr->second;
#endif
              // do not assemble zeros into matrix
              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim(); j++)
                  if (abs(valtxi*txi[j])>1.e-12)
                    linslipDISglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
              }
              else
                if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);

              if (Dim()==3)
              {
                if (constr_direction_==INPAR::CONTACT::constr_xyz)
                {
                  for (int j=0; j<Dim(); j++)
                    if (abs(valteta*teta[j])>1.e-12)
                      linslipDISglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
                }
                else
                  if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
              }
            }
          }
        }
        /*** 1 ********************************* euclidean.deriv(T).z ***/
        // loop over dimensions
        for (int j=0;j<Dim();++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = euclidean*(colcurr->second)*z[j];

#ifdef CONSISTENTSLIP
            val = val / (znor - cn * wgap);
#endif

            // do not assemble zeros into s matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(val*txi[j])>1.e-12)
                  linslipDISglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
          }

          if (Dim()==3)
          {
            // loop over all entries of the current derivative map (teta)
            for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double val = euclidean*(colcurr->second)*z[j];

#ifdef CONSISTENTSLIP
              val = val / (znor - cn * wgap);
#endif

              // do not assemble zeros into s matrix
              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim(); j++)
                  if (abs(val*teta[j])>1.e-12)
                    linslipDISglobal.Assemble(val*teta[j],cnode->Dofs()[j],col);
              }
              else
                if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[1],col);
            }
          }
        }

        /*** 2 ********************* deriv(euclidean).deriv(T).z.ztan ***/
        // loop over dimensions
        for (int j=0;j<Dim();++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double valtxi = (ztxi+ct*jumptxi)/euclidean*(colcurr->second)*z[j]*ztxi;
            double valteta = (ztxi+ct*jumptxi)/euclidean*(colcurr->second)*z[j]*zteta;

#ifdef CONSISTENTSLIP
            valtxi  = valtxi / (znor - cn*wgap);
            valteta   = valteta / (znor - cn*wgap);
#endif

            // do not assemble zeros into matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(valtxi*txi[j])>1.e-12)
                  linslipDISglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);

            if (Dim()==3)
            {
              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim(); j++)
                  if (abs(valteta*teta[j])>1.e-12)
                    linslipDISglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
              }
              else
                if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }
          }

          if(Dim()==3)
          {
            // 3D loop over all entries of the current derivative map (teta)
            for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi = (zteta+ct*jumpteta)/euclidean*(colcurr->second)*z[j]*ztxi;
              double valteta = (zteta+ct*jumpteta)/euclidean*(colcurr->second)*z[j]*zteta;

#ifdef CONSISTENTSLIP
              valtxi  = valtxi / (znor - cn*wgap);
              valteta = valteta / (znor - cn*wgap);
#endif

              // do not assemble zeros into matrix
              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim(); j++)
                  if (abs(valtxi*txi[j])>1.e-12)
                    linslipDISglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
              }
              else
                if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);

              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim(); j++)
                  if (abs(valteta*teta[j])>1.e-12)
                    linslipDISglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
              }
              else
                if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }
          }
        }

        /*** 3 ****************** deriv(euclidean).deriv(T).jump.ztan ***/
        if (gp_slip)
        {
          //!!!!!!!!!!!!!!! DO NOTHING !!!!!!!
        }
        else
        {
          // loop over dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
            {
              int col = colcurr->first;
              double valtxi = (ztxi+ct*jumptxi)/euclidean*ct*(colcurr->second)*jump[j]*ztxi;
              double valteta = (ztxi+ct*jumptxi)/euclidean*ct*(colcurr->second)*jump[j]*zteta;

#ifdef CONSISTENTSLIP
              valtxi  = valtxi / (znor - cn*wgap);
              valteta   = valteta / (znor - cn*wgap);
#endif

              // do not assemble zeros into s matrix
              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim(); j++)
                  if (abs(valtxi*txi[j])>1.e-12)
                    linslipDISglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
                for (int j=0; j<Dim(); j++)
                  if (abs(valteta*teta[j])>1.e-12)
                    linslipDISglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
              }
              else
              {
                if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
                if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
              }
            }

            if(Dim()==3)
            {
              // loop over all entries of the current derivative map (teta)
              for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
              {
                int col = colcurr->first;
                double valtxi = (zteta+ct*jumpteta)/euclidean*ct*(colcurr->second)*jump[j]*ztxi;
                double valteta = (zteta+ct*jumpteta)/euclidean*ct*(colcurr->second)*jump[j]*zteta;

#ifdef CONSISTENTSLIP
                valtxi  = valtxi / (znor - cn*wgap);
                valteta   = valteta / (znor - cn*wgap);
#endif

                // do not assemble zeros into matrix
                if (constr_direction_==INPAR::CONTACT::constr_xyz)
                {
                  for (int j=0; j<Dim(); j++)
                    if (abs(valtxi*txi[j])>1.e-12)
                      linslipDISglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
                  for (int j=0; j<Dim(); j++)
                    if (abs(valteta*teta[j])>1.e-12)
                      linslipDISglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
                }
                else
                {
                  if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
                  if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
                }
              }
            }
          }
        }
      /*** 4 ************************** (frcoeff*znor).deriv(T).z ***/
        // loop over all dimensions
        for (int j=0;j<Dim();++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
          {
            int col = colcurr->first;
#ifdef CONSISTENTSLIP
            double val = - frcoeff * (colcurr->second)*z[j];
#else
            double val = (-1)*(frcoeff*(znor-cn*wgap))*(colcurr->second)*z[j];
#endif
            // do not assemble zeros into matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(val*txi[j])>1.e-12)
                  linslipDISglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
          }

          if(Dim()==3)
          {
            // loop over all entries of the current derivative map (teta)
            for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
            {
              int col = colcurr->first;
#ifdef CONSISTENTSLIP
              double val = - frcoeff * (colcurr->second) * z[j];
#else
              double val = (-1.0)*(frcoeff*(znor-cn*wgap))*(colcurr->second)*z[j];
#endif
              // do not assemble zeros into matrix
              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim(); j++)
                  if (abs(val*teta[j])>1.e-12)
                    linslipDISglobal.Assemble(val*teta[j],cnode->Dofs()[j],col);
              }
              else
                if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[1],col);
            }
          }
        }

        /*** 5 *********************** (frcoeff*znor).deriv(T).jump ***/
        if (gp_slip)
        {
          //!!!!!!!!!!!!!!! DO NOTHING !!!!!!!!!!
        }
        else
        {
          // loop over all dimensions
          for (int j=0;j<Dim();++j)
          {
            // loop over all entries of the current derivative map (txi)
            for (colcurr=dtximap[j].begin();colcurr!=dtximap[j].end();++colcurr)
            {
              int col = colcurr->first;
#ifdef CONSISTENTSLIP
              double val = - frcoeff * ct * (colcurr->second) * jump[j];
#else
              double val = (-1)*(frcoeff*(znor-cn*wgap))*ct*(colcurr->second)*jump[j];
#endif
              // do not assemble zeros into matrix
              if (constr_direction_==INPAR::CONTACT::constr_xyz)
              {
                for (int j=0; j<Dim(); j++)
                  if (abs(val*txi[j])>1.e-12)
                    linslipDISglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
              }
              else
                if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[0],col);
            }

            if(Dim()==3)
            {
              // loop over all entries of the current derivative map (teta)
              for (colcurr=dtetamap[j].begin();colcurr!=dtetamap[j].end();++colcurr)
              {
                int col = colcurr->first;
#ifdef CONSISTENTSLIP
                double val = - frcoeff * ct * (colcurr->second) * jump[j];
#else
                double val = (-1)*(frcoeff*(znor-cn*wgap))*ct*(colcurr->second)*jump[j];
#endif
              // do not assemble zeros into s matrix
                if (constr_direction_==INPAR::CONTACT::constr_xyz)
                {
                  for (int j=0; j<Dim(); j++)
                    if (abs(val*teta[j])>1.e-12)
                      linslipDISglobal.Assemble(val*teta[j],cnode->Dofs()[j],col);
                }
                else
                  if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row[1],col);
              }
            }
          }
        }

#ifndef CONSISTENTSLIP
        /*** 6 ******************* -frcoeff.Deriv(n).z(ztan+ct*utan) ***/
        // loop over all dimensions
        for (int j=0;j<Dim();++j)
        {
          // loop over all entries of the current derivative map
          for (colcurr=dnmap[j].begin();colcurr!=dnmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double valtxi = (-1)*(ztxi+ct*jumptxi)*frcoeff*(colcurr->second)*z[j];
            double valteta = (-1)*(zteta+ct*jumpteta)*frcoeff*(colcurr->second)*z[j];

            // do not assemble zeros into s matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(valtxi*txi[j])>1.e-12)
                  linslipDISglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
              for (int j=0; j<Dim(); j++)
                if (abs(valteta*teta[j])>1.e-12)
                  linslipDISglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
            }
            else
            {
              if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
              if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
            }
          }
        }

        /*** 7 ****************** frcoeff*cn*deriv (g).(ztan+ct*utan) ***/
        // prepare assembly
        std::map<int,double>& dgmap = cnode->CoData().GetDerivG();

        // loop over all entries of the current derivative map
        for (_colcurr=dgmap.begin();_colcurr!=dgmap.end();++_colcurr)
        {
          int col = _colcurr->first;
          double valtxi = frcoeff*cn*(_colcurr->second)*(ztxi+ct*jumptxi);
          double valteta = frcoeff*cn*(_colcurr->second)*(zteta+ct*jumpteta);

          // do not assemble zeros into matrix
          if (constr_direction_==INPAR::CONTACT::constr_xyz)
          {
            for (int j=0; j<Dim(); j++)
              if (abs(valtxi*txi[j])>1.e-12)
                linslipDISglobal.Assemble(valtxi*txi[j],cnode->Dofs()[j],col);
            for (int j=0; j<Dim(); j++)
              if (abs(valteta*teta[j])>1.e-12)
                linslipDISglobal.Assemble(valteta*teta[j],cnode->Dofs()[j],col);
          }
          else
          {
            if (abs(valtxi)>1.0e-12) linslipDISglobal.Assemble(valtxi,row[0],col);
            if (abs(valteta)>1.0e-12) linslipDISglobal.Assemble(valteta,row[1],col);
          }
        }
#endif
      } // if (frcoeff==0.0)
    } // loop over all slip nodes of the interface
  } // Coulomb friction

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
  // Tresca Friction
  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
  if (ftype == INPAR::CONTACT::friction_tresca)
  {
    // loop over all slip nodes of the interface
    for (int i=0; i<slipnodes_->NumMyElements();++i)
    {
      int gid = slipnodes_->GID(i);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node) dserror("ERROR: Cannot find node with gid %",gid);
      FriNode* cnode = dynamic_cast<FriNode*>(node);

      if (cnode->Owner() != Comm().MyPID())
        dserror("ERROR: AssembleLinSlip: Node ownership inconsistency!");

      double ct = GetCtRef()[GetCtRef().Map().LID(cnode->Id())];

      // preparation of assembly
      // get Deriv N and calculate DerivD form DerivN

      // only for 2D so far, in this case calculation is very easy
      // dty =  dnx
      // dtx = -dny
      // FIXGIT: in the future DerivD will be called directly form node

      std::vector<GEN::pairedvector<int,double> > dnmap = cnode->CoData().GetDerivN();

      // iterator
      GEN::pairedvector<int,double>::iterator _colcurr;
      std::map<int,double>::iterator colcurr;

      std::vector <std::map<int,double> > dtmap(Dim());

      for (_colcurr=dnmap[0].begin(); _colcurr!=dnmap[0].end(); _colcurr++)
        dtmap[1].insert(std::pair<int,double>(_colcurr->first,_colcurr->second));

      for (_colcurr=dnmap[1].begin(); _colcurr!=dnmap[1].end(); _colcurr++)
        dtmap[0].insert(std::pair<int,double>(_colcurr->first,(-1)*_colcurr->second));

      // get more information from node
      double* jump = cnode->FriData().jump();
      double* txi = cnode->CoData().txi();
      double* xi = cnode->xspatial();
      double* z = cnode->MoData().lm();
      int row = slipt_->GID(i);

      int colsize = (int)dtmap[0].size();
      int mapsize = (int)dtmap.size();

      for (int j=0;j<mapsize-1;++j)
        if ((int)dtmap[j].size() != (int)dtmap[j+1].size())
          dserror("ERROR: AssembleLinSlip: Column dim. of nodal DerivT-map is inconsistent!");

      // calculation of parts of the complementary function
      double ztan    = txi[0]*z[0] + txi[1]*z[1];
      double jumptan = txi[0]*jump[0] + txi[1]*jump[1];
      //double temp = ztan + ct*jumptan;
      //double epk = frbound/abs(temp);
      //double Fpk = ztan*temp/(frbound*abs(temp));
      //double Mpk = epk*(1-Fpk);
      //double fac = 1/(abs(ztan+ct*jumptan))*1/(1-Mpk)*(-1);

      // calculation of |ztan+ct*utan|
      double sum = 0;
      int prefactor = 1;
      for (int dim = 0;dim < Dim();dim++)
        sum += txi[dim]*z[dim]+ct*txi[dim]*jump[dim];

      // calculate |sum| and prefactor
      if (sum < 0)
      {
        sum = -sum;
        prefactor = (-1);
      }

      //****************************************************************
      // CONSISTENT TREATMENT OF CASE FRBOUND=0 (FRICTIONLESS)
      //****************************************************************
      // popp 08/2012
      //
      // There is a problem with our frictional nonlinear complementarity
      // function when applied to the limit case frbound=0 (frictionless).
      // In this case, the simple frictionless sliding condition should
      // be consistently recovered, which unfortunately is not the case.
      // This fact is well-known (see PhD thesis S. Hüeber) and now
      // taken care of by a special treatment as can be seen below
      //
      //****************************************************************
      if (frbound==0.0
          || (cnode->CoData().ActiveOld()==false
              && frilessfirst))
      {
        //**************************************************************
        // calculation of matrix entries of linearized slip condition
        //**************************************************************

        // 1) Entries from differentiation with respect to LM
        /******************************************************************/

        // loop over the dimension
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cnode->Dofs()[dim];
          double valtxi = txi[dim];

          // do not assemble zeros into matrix
          if (abs(valtxi)>1.0e-12) linslipLMglobal.Assemble(valtxi,row,col);
        }

        // 2) Entries on right hand side
        /******************************************************************/

        Epetra_SerialDenseVector rhsnode(1);
        std::vector<int> lm(1);
        std::vector<int> lmowner(1);

        rhsnode(0) = -ztan;
        lm[0] = cnode->Dofs()[1];
        lmowner[0] = cnode->Owner();

        LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);

        // 3) Entries from differentiation with respect to displacements
        /******************************************************************/

        // loop over dimensions
        for (int j=0;j<Dim();++j)
        {
          // loop over all entries of the current derivative map (txi)
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (colcurr->second)*z[j];

            // do not assemble zeros into matrix
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }
      }

      //****************************************************************
      // STANDARD TREATMENT OF CASE FRBOUND!=0 (FRICTIONAL)
      //****************************************************************
      else
      {
        //****************************************************************
        // calculation of matrix entries of the linearized slip condition
        //****************************************************************

        // 1) Entries from differentiation with respect to LM
        /**************** (Deriv(abs)*ztan+|ztan+ct*jumptan|-frbound).tan ***/

        // loop over the dimension
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cnode->Dofs()[dim];
          double val = (prefactor*ztan+sum-frbound)*txi[dim];

          // do not assemble zeros into matrix
          if (constr_direction_==INPAR::CONTACT::constr_xyz)
          {
            for (int j=0; j<Dim(); j++)
              if (abs(val*txi[j])>1.e-12)
                linslipLMglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
          }
          else
            if (abs(val)>1.0e-12) linslipLMglobal.Assemble(val,row,col);
        }

        // 2) Entries on right hand side
        /************ -C + entries from writing Delta(z) as z(k+1)-z(k) ***/

        // -C and remaining terms
        double value1= -(abs(ztan+ct*jumptan))*ztan+frbound*(ztan+ct*jumptan);

        if (constr_direction_==INPAR::CONTACT::constr_xyz)
        {
          Epetra_SerialDenseVector rhsnode(Dim());
          std::vector<int> lm(Dim());
          std::vector<int> lmowner(Dim());

          for (int j=0; j<Dim(); j++)
          {
            lm[j] = cnode->Dofs()[j];
            lmowner[j] = cnode->Owner();
            rhsnode(j) = value1*txi[j];
          }

          LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);
        }
        else
        {
          Epetra_SerialDenseVector rhsnode(1);
          rhsnode(0) = value1;

          std::vector<int> lm(1);
          std::vector<int> lmowner(1);

          lm[0] = cnode->Dofs()[1];
          lmowner[0] = cnode->Owner();

          LINALG::Assemble(linslipRHSglobal,rhsnode,lm,lmowner);
        }

        // 3) Entries from differentiation with respect to displacements
        /***************************** -Deriv(abs)*ct*tan.(D-Dn-1)*ztan ***/

        // we need the nodal entries of the D-matrix and the old one
        double D= cnode->MoData().GetD()[cnode->Id()];
        double Dold= cnode->FriData().GetDOld()[cnode->Id()];

        if (abs(Dold)<0.0001)
          dserror ("Error:No entry for Dold");

        // loop over all derivative maps (=dimensions)
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cnode->Dofs()[dim];
          double val = prefactor*(-1)*ct*txi[dim]*(D-Dold)*ztan;
          //std::cout << "01 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

          // do not assemble zeros into matrix
          if (constr_direction_==INPAR::CONTACT::constr_xyz)
          {
            for (int j=0; j<Dim(); j++)
              if (abs(val*txi[j])>1.e-12)
                linslipDISglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
          }
          else
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }

        /***************************** -Deriv(abs)*ct*tan.(M-Mn-1)*ztan ***/

        // we need the nodal entries of the M-matrix and the old one
        std::map<int,double>& mmap = cnode->MoData().GetM();
        std::map<int,double>& mmapold = cnode->FriData().GetMOld();

        // create a set of nodes including nodes according to M entries
        // from current and previous time step
        std::set <int> mnodes;

        // iterator
        std::set<int>::iterator mcurr;

        std::set <int> mnodescurrent = cnode->FriData().GetMNodes();
        std::set <int> mnodesold = cnode->FriData().GetMNodesOld();

        for (mcurr=mnodescurrent.begin(); mcurr != mnodescurrent.end(); mcurr++)
          mnodes.insert(*mcurr);

        for (mcurr=mnodesold.begin(); mcurr != mnodesold.end(); mcurr++)
          mnodes.insert(*mcurr);

        // loop over all master nodes (find adjacent ones to this stick node)
        for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
        {
          int gid = *mcurr;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          FriNode* cmnode = dynamic_cast<FriNode*>(mnode);

          double mik = mmap[cmnode->Id()];
          double mikold = mmapold[cmnode->Id()];

          // compute linstick-matrix entry of the current active node / master node pair
          // loop over all derivative maps (=dimensions)
          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            int col = cmnode->Dofs()[dim];
            double val = prefactor*(+1)*ct*txi[dim]*(mik-mikold)*ztan;
            //std::cout << "02 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

            // do not assemble zeros into matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(val*txi[j])>1.e-12)
                  linslipDISglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }

        /************************************** frbound*ct*tan.(D-Dn-1) ***/

        // loop over all derivative maps (=dimensions)
        for (int dim=0;dim<cnode->NumDof();++dim)
        {
          int col = cnode->Dofs()[dim];
          double val = frbound*ct*txi[dim]*(D-Dold);
          //std::cout << "03 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

          // do not assemble zeros into matrix
          if (constr_direction_==INPAR::CONTACT::constr_xyz)
          {
            for (int j=0; j<Dim(); j++)
              if (abs(val*txi[j])>1.e-12)
                linslipDISglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
          }
          else
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }

        /********************************** -frbound*ct*tan.(M-Mn-1).xm ***/

        // loop over all master nodes
        for (mcurr=mnodes.begin(); mcurr != mnodes.end(); mcurr++)
        {
          int gid = *mcurr;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          FriNode* cmnode = dynamic_cast<FriNode*>(mnode);

          double mik = mmap[cmnode->Id()];
          double mikold = mmapold[cmnode->Id()];

          // loop over all derivative maps (=dimensions)
          for (int dim=0;dim<cnode->NumDof();++dim)
          {
            int col = cmnode->Dofs()[dim];
            double val = frbound*(-1)*ct*txi[dim]*(mik-mikold);
            //std::cout << "04 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;
            // do not assemble zeros into matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(val*txi[j])>1.e-12)
                  linslipDISglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }

        /************************************ |ztan+ct*utan|.DerivT.z ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = sum*(colcurr->second)*z[j];
            //std::cout << "1 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

            // do not assemble zeros into s matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(val*txi[j])>1.e-12)
                  linslipDISglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);

            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /*********************************** Deriv(abs)*DerivT.z*ztan ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*(colcurr->second)*z[j]*ztan;
            //std::cout << "2 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

            // do not assemble zeros into matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(val*txi[j])>1.e-12)
                  linslipDISglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);

            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /******************************* Deriv(abs)*DerivT.jump+*ztan ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*ct*(colcurr->second)*jump[j]*ztan;
            //std::cout << "3 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

            // do not assemble zeros into s matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(val*txi[j])>1.e-12)
                  linslipDISglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);

            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /*************************** -Deriv(abs).ct.tan.DerivD.x*ztan ***/

        // we need the dot product t*x of this node
        double tdotx = 0.0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          tdotx += txi[dim]*xi[dim];

        // prepare assembly
        std::map<int,double>& ddmap = cnode->CoData().GetDerivD()[gid];

        // loop over all entries of the current derivative map
        for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
        {
          int col = colcurr->first;
          double val = (-1)*prefactor*ct*tdotx*colcurr->second*ztan;
          //std::cout << "4 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

          // do not assemble zeros into matrix
          if (constr_direction_==INPAR::CONTACT::constr_xyz)
          {
            for (int j=0; j<Dim(); j++)
              if (abs(val*txi[j])>1.e-12)
                linslipDISglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
          }
          else
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }

        /**************************** Deriv(abs).ct.tan.DerivM.x*ztan ***/

        // we need the Lin(M-matrix) entries of this node
        std::map<int,std::map<int,double> >& dmmap = cnode->CoData().GetDerivM();
        std::map<int,std::map<int,double> >::iterator dmcurr;

        // loop over all master nodes in the DerivM-map of the active slave node
        for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
        {
          int gid = dmcurr->first;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          FriNode* cmnode = dynamic_cast<FriNode*>(mnode);
          double* mxi = cmnode->xspatial();

          // we need the dot product ns*xm of this node pair
          double tdotx = 0.0;
          for (int dim=0;dim<cnode->NumDof();++dim)
            tdotx += txi[dim]*mxi[dim];

          // compute entry of the current active node / master node pair
          std::map<int,double>& thisdmmap = cnode->CoData().GetDerivM(gid);

          // loop over all entries of the current derivative map
          for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
          {
            int col = colcurr->first;
            double val = prefactor*ct*tdotx*colcurr->second*ztan;
            //std::cout << "5 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

            // do not assemble zeros into matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(val*txi[j])>1.e-12)
                  linslipDISglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }

        /****************************************** -frbound.DerivT.z ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (-1)*frbound*(colcurr->second)*z[j];
            //std::cout << "6 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

            // do not assemble zeros into s matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(val*txi[j])>1.e-12)
                  linslipDISglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);

            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /************************************ -frbound.ct.DerivT.jump ***/

        // loop over all derivative maps (=dimensions)
        for (int j=0;j<mapsize;++j)
        {
          int k=0;

          // loop over all entries of the current derivative map
          for (colcurr=dtmap[j].begin();colcurr!=dtmap[j].end();++colcurr)
          {
            int col = colcurr->first;
            double val = (-1)*frbound*ct*(colcurr->second)*jump[j];
            //std::cout << "7 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

            // do not assemble zeros into s matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(val*txi[j])>1.e-12)
                  linslipDISglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);

            ++k;
          }

          if (k!=colsize)
            dserror("ERROR: AssembleLinSlip: k = %i but colsize = %i",k,colsize);
        }

        /************************************* +frbound.ct.T.DerivD.x ***/

        // we need the dot product t*x of this node
        tdotx = 0.0;
        for (int dim=0;dim<cnode->NumDof();++dim)
          tdotx += txi[dim]*xi[dim];

        // loop over all entries of the current derivative map
        for (colcurr=ddmap.begin();colcurr!=ddmap.end();++colcurr)
        {
          int col = colcurr->first;
          double val = (-1)*(-1)*frbound*ct*tdotx*colcurr->second;
          //std::cout << "8 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;

          // do not assemble zeros into matrix
          if (constr_direction_==INPAR::CONTACT::constr_xyz)
          {
            for (int j=0; j<Dim(); j++)
              if (abs(val*txi[j])>1.e-12)
                linslipDISglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
          }
          else
            if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
        }

        /********************************  -frbound.ct.T.DerivM.x ******/

        // loop over all master nodes in the DerivM-map of the active slave node
        for (dmcurr=dmmap.begin();dmcurr!=dmmap.end();++dmcurr)
        {
          int gid = dmcurr->first;
          DRT::Node* mnode = idiscret_->gNode(gid);
          if (!mnode) dserror("ERROR: Cannot find node with gid %",gid);
          FriNode* cmnode = dynamic_cast<FriNode*>(mnode);
          double* mxi = cmnode->xspatial();

          // we need the dot product ns*xm of this node pair
          double tdotx = 0.0;
          for (int dim=0;dim<cnode->NumDof();++dim)
            tdotx += txi[dim]*mxi[dim];

          // compute entry of the current active node / master node pair
          std::map<int,double>& thisdmmap = cnode->CoData().GetDerivM(gid);

          // loop over all entries of the current derivative map
          for (colcurr=thisdmmap.begin();colcurr!=thisdmmap.end();++colcurr)
          {
            int col = colcurr->first;
            double val = (-1)*frbound*ct*tdotx*colcurr->second;
            //std::cout << "9 GID " << gid << " row " << row << " col " << col << " val " << val << std::endl;


            // do not assemble zeros into matrix
            if (constr_direction_==INPAR::CONTACT::constr_xyz)
            {
              for (int j=0; j<Dim(); j++)
                if (abs(val*txi[j])>1.e-12)
                  linslipDISglobal.Assemble(val*txi[j],cnode->Dofs()[j],col);
            }
            else
              if (abs(val)>1.0e-12) linslipDISglobal.Assemble(val,row,col);
          }
        }
      }
    }
  }// Tresca friction

  return;
}

/*----------------------------------------------------------------------*
 |  build active set (nodes / dofs)                           popp 02/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::BuildActiveSet(bool init)
{
  // define local variables
  std::vector<int> mynodegids(0);
  std::vector<int> mynodegidsInactive(0);
  std::vector<int> mydofgids(0);
  std::vector<int> mydofgidsInactive(0);
  std::vector<int> myslipnodegids(0);
  std::vector<int> myslipdofgids(0);
  std::vector<int> mymnodegids(0);
  std::vector<int> mymdofgids(0);

  // loop over all slave nodes
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);
    const int numdof = cnode->NumDof();

    // *******************************************************************
    // INITIALIZATION OF THE ACTIVE SET (t=0)
    // *******************************************************************
    // This is given by the CoNode member variable IsInitActive(), which
    // has been introduced via the contact conditions in the input file.
    // Thus, if no design line has been chosen to be active at t=0,
    // the active node set will be empty at t=0. Yet, if one or more
    // design lines have been specified as "Slave" AND "Active" then
    // the corresponding CoNodes are put into an initial active set!
    // This yields a very flexible solution for contact initialization.
    // *******************************************************************
    if (init)
    {
      // flag for initialization of init active nodes with nodal gaps
      bool initcontactbygap = DRT::INPUT::IntegralValue<int>(IParams(),"INITCONTACTBYGAP");
      // value
      double initcontactval = IParams().get<double>("INITCONTACTGAPVALUE");

      // Either init contact by definition or by gap
      if(cnode->IsInitActive() and initcontactbygap)
        dserror("Init contact either by definition in condition or by gap!");

      // check if node is initially active or, if initialization with nodal, gap,
      // the gap is smaller than the prescribed value
      if (
          cnode->IsInitActive() or (initcontactbygap and cnode->CoData().Getg() < initcontactval)
//          sqrt(cnode->X()[0]*cnode->X()[0]+cnode->X()[1]*cnode->X()[1])>65.
//          (sqrt(cnode->X()[0]*cnode->X()[0]+cnode->X()[1]*cnode->X()[1])>12. && abs(cnode->X()[2]-10.5)<1.e-12)
//          || abs(cnode->X()[2]-0.5)<1.e-12
         )
      {
/*
      // **********************************************************************************
      // ***                     CONSISTENT INITIALIZATION STATE                        ***
      // **********************************************************************************
      // hiermeier 08/2013
      //
      // The consistent nonlinear complementarity function C_tau is for the special initialization case
      //
      //                zn == 0   and   gap == 0
      //
      // in conjunction with the Coulomb friction model no longer unique and a special treatment is
      // necessary. For this purpose we identify the critical nodes here and set the slip-state to true.
      // In the next step the tangential part of the Lagrange multiplier vector will be set to zero. Hence,
      // we treat these nodes as frictionless nodes! (see AssembleLinSlip)
      INPAR::CONTACT::FrictionType ftype =
          DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(IParams(),"FRICTION");
      if (ftype == INPAR::CONTACT::friction_coulomb)
      {
      dynamic_cast<FriNode*>(cnode)->FriData().InconInit() = true;
      myslipnodegids.push_back(cnode->Id());

      for (int j=0;j<numdof;++j)
        myslipdofgids.push_back(cnode->Dofs()[j]);
      }

*/
        cnode->Active()=true;
        mynodegids.push_back(cnode->Id());

        for (int j=0;j<numdof;++j)
          mydofgids.push_back(cnode->Dofs()[j]);
      }

      // check if frictional node is initially in slip state
      if (friction_)
      {
        // do nothing: we always assume STICK at t=0
      }
    }

    // *******************************************************************
    // RE-BUILDING OF THE ACTIVE SET
    // *******************************************************************
    else
    {
      // check if node is active
      if (cnode->Active())
      {
        mynodegids.push_back(cnode->Id());

        for (int j=0;j<numdof;++j)
          mydofgids.push_back(cnode->Dofs()[j]);
      }
      else
      {
        mynodegidsInactive.push_back(cnode->Id());

        for (int j=0;j<numdof;++j)
          mydofgidsInactive.push_back(cnode->Dofs()[j]);
      }

      // check if frictional node is in slip state
      if (friction_)
      {
        if (dynamic_cast<FriNode*>(cnode)->FriData().Slip())
        {
          myslipnodegids.push_back(cnode->Id());

          for (int j=0;j<numdof;++j)
            myslipdofgids.push_back(cnode->Dofs()[j]);
        }
      }
    }
  }

  // create active node map and active dof map
  activenodes_ = LINALG::CreateMap(mynodegids,Comm());
  activedofs_  = LINALG::CreateMap(mydofgids,Comm());
  inactivenodes_ = LINALG::CreateMap(mynodegidsInactive,Comm());
  inactivedofs_  = LINALG::CreateMap(mydofgidsInactive,Comm());

  if (friction_)
  {
    // create slip node map and slip dof map
    slipnodes_ = LINALG::CreateMap(myslipnodegids,Comm());
    slipdofs_  = LINALG::CreateMap(myslipdofgids,Comm());
  }

  // split active dofs and slip dofs
  SplitActiveDofs();

  return true;
}

/*----------------------------------------------------------------------*
 |  split active dofs into Ndofs, Tdofs and slipTdofs         popp 02/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CoInterface::SplitActiveDofs()
{
  // get out of here if active set is empty
  if (activenodes_==Teuchos::null)
  {
    activen_ = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
    activet_ = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
    slipt_ = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
    return true;
  }

  else if (activenodes_->NumGlobalElements()==0)
  {
    activen_ = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
    activet_ = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
    slipt_ = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
    return true;
  }

  // define local variables
  int countN=0;
  int countT=0;
  std::vector<int> myNgids(activenodes_->NumMyElements());
  std::vector<int> myTgids((Dim()-1)*activenodes_->NumMyElements());

  // dimension check
  double dimcheck =(activedofs_->NumGlobalElements())/(activenodes_->NumGlobalElements());
  if (dimcheck != Dim()) dserror("ERROR: SplitActiveDofs: Nodes <-> Dofs dimension mismatch!");

  // loop over all active row nodes
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);
    const int numdof = cnode->NumDof();

    // add first dof to Nmap
    myNgids[countN] = cnode->Dofs()[0];
    ++countN;

    // add remaining dofs to Tmap
    for (int j=1;j<numdof;++j)
    {
      myTgids[countT] = cnode->Dofs()[j];
      ++countT;
    }
  }

  // resize the temporary vectors
  myNgids.resize(countN);
  myTgids.resize(countT);

  // communicate countN and countT among procs
  int gcountN, gcountT;
  Comm().SumAll(&countN,&gcountN,1);
  Comm().SumAll(&countT,&gcountT,1);

  // check global dimensions
  if ((gcountN+gcountT)!=activedofs_->NumGlobalElements())
    dserror("ERROR: SplitActiveDofs: Splitting went wrong!");

  // create Nmap and Tmap objects
  activen_ = Teuchos::rcp(new Epetra_Map(gcountN,countN,&myNgids[0],0,Comm()));
  activet_ = Teuchos::rcp(new Epetra_Map(gcountT,countT,&myTgids[0],0,Comm()));

  // *******************************************************************
  // FRICTION - EXTRACTING TANGENTIAL DOFS FROM SLIP DOFS
  // *******************************************************************

  // get out of here if there is no friction
  if(friction_==false)
    return true;

  // get out of here if slip set is empty
  if (slipnodes_==Teuchos::null)
  {
    slipt_ = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
    return true;
  }

  if (slipnodes_->NumGlobalElements()==0)
  {
    slipt_ = Teuchos::rcp(new Epetra_Map(0,0,Comm()));
    return true;
  }

  // define local variables
  int countslipT=0;
  std::vector<int> myslipTgids((Dim()-1)*slipnodes_->NumMyElements());

  // dimension check
  dimcheck =(slipdofs_->NumGlobalElements())/(slipnodes_->NumGlobalElements());
  if (dimcheck != Dim()) dserror("ERROR: SplitActiveDofs: Nodes <-> Dofs dimension mismatch!");

  // loop over all slip row nodes
  for (int i=0;i<slipnodes_->NumMyElements();++i)
  {
    int gid = slipnodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);
    const int numdof = cnode->NumDof();

    // add dofs to slipTmap
    for (int j=1;j<numdof;++j)
    {
      myslipTgids[countslipT] = cnode->Dofs()[j];
      ++countslipT;
    }
  }

  // resize the temporary vectors
  myslipTgids.resize(countslipT);

  // communicate countslipT among procs
  int gcountslipT;
  Comm().SumAll(&countslipT,&gcountslipT,1);

  // create Tslipmap objects
  slipt_   = Teuchos::rcp(new Epetra_Map(gcountslipT,countslipT,&myslipTgids[0],0,Comm()));

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::CoInterface::GetForceOfNode(
    LINALG::Matrix<3,1>& nodal_force,
    const Epetra_Vector& force,
    const DRT::Node& node ) const
{
  std::vector<int> dofs;
  idiscret_->Dof(&node,dofs);

  // reset nodal force vector
  std::fill( nodal_force.A(), nodal_force.A()+3, 0.0 );
  const double* f_vals = force.Values();

  if ( dofs.size()>3 )
    dserror("The interface node seems to have more than 3 DOFs!");

  for ( unsigned i=0; i<dofs.size(); ++i )
  {
    const int dof = dofs[i];
    const int f_lid = force.Map().LID( dof );
    if ( f_lid == -1 )
      dserror("Couldn't find the nodal DOF %d in the global force vector!",
          dof );

    nodal_force(i,0) = f_vals[f_lid];
  }
}

/*----------------------------------------------------------------------*
 |  Calculate angular interface moments                  hiermeier 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::EvalResultantMoment(
    const Epetra_Vector& fs,
    const Epetra_Vector& fm,
    LINALG::SerialDenseMatrix* conservation_data_ptr) const
{
  static int step = 0;

  // get out of here if not participating in interface
  if (!lComm()) return;

  double lresMoSl[3] = {0.0,0.0,0.0};   // local slave moment
  double gresMoSl[3] = {0.0,0.0,0.0};   // global slave moment
  double lresMoMa[3] = {0.0,0.0,0.0};   // local master momemnt
  double gresMoMa[3] = {0.0,0.0,0.0};   // global master moment
  double gbalMo[3]   = {0.0,0.0,0.0};   // global moment balance

  double lresFSl[3] = {0.0,0.0,0.0};  // local slave force
  double gresFSl[3] = {0.0,0.0,0.0};  // global slave force
  double lresFMa[3] = {0.0,0.0,0.0};  // local master force
  double gresFMa[3] = {0.0,0.0,0.0};  // global master force
  double gbalF[3]   = {0.0,0.0,0.0};  // global force balance

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  LINALG::Matrix<3,1> nforce(true);
  for (int i=0;i<snoderowmap_->NumMyElements();++i)
  {
    int gid = snoderowmap_->GID(i);
    CoNode* snode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));

    GetForceOfNode(nforce, fs, *snode);

    if (Dim()==2)
      lresMoSl[2] += snode->xspatial()[0] * nforce(1,0) - snode->xspatial()[1]*nforce(0,0);
    else
    {
      lresMoSl[0] += snode->xspatial()[1] * nforce(2,0) - snode->xspatial()[2] * nforce(1,0);
      lresMoSl[1] += snode->xspatial()[2] * nforce(0,0) - snode->xspatial()[0] * nforce(2,0);
      lresMoSl[2] += snode->xspatial()[0] * nforce(1,0) - snode->xspatial()[1] * nforce(0,0);
    }
    for ( unsigned k=0; k<3; ++k )
      lresFSl[k] += nforce(k,0);
  }
  Comm().SumAll(&lresMoSl[0],&gresMoSl[0],3);
  Comm().SumAll(&lresFSl[0],&gresFSl[0],3);

  // loop over proc's master nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<mnoderowmap_->NumMyElements();++i)
  {
    int gid = mnoderowmap_->GID(i);
    CoNode* mnode = dynamic_cast<CoNode*>(idiscret_->gNode(gid));

    GetForceOfNode(nforce, fm, *mnode);

    if (Dim()==2)
      lresMoMa[2] += mnode->xspatial()[0] * nforce(1,0) - mnode->xspatial()[1]*nforce(0,0);
    else
    {
      lresMoMa[0] += mnode->xspatial()[1] * nforce(2,0) - mnode->xspatial()[2] * nforce(1,0);
      lresMoMa[1] += mnode->xspatial()[2] * nforce(0,0) - mnode->xspatial()[0] * nforce(2,0);
      lresMoMa[2] += mnode->xspatial()[0] * nforce(1,0) - mnode->xspatial()[1] * nforce(0,0);
    }
    for ( unsigned k=0; k<3; ++k )
      lresFMa[k] += nforce(k,0);
  }
  Comm().SumAll(&lresMoMa[0],&gresMoMa[0],3);
  Comm().SumAll(&lresFMa[0],&gresFMa[0],3);

  for (int d=0;d<3;++d)
  {
    gbalMo[d] = gresMoSl[d] + gresMoMa[d];
    gbalF[d] = gresFSl[d] + gresFMa[d];

//    if (abs(balMo[d])>1.0e-11) dserror("Conservation of angular momentum is not fulfilled!");
  }
  if (Comm().MyPID()==0)
  {
    std::cout << "SLAVE:   " <<
        " [" << std::setw(14) << std::setprecision(5) << std::scientific << gresMoSl[0] <<
        ", " << std::setw(14) << std::setprecision(5) << std::scientific << gresMoSl[1] <<
        ", " << std::setw(14) << std::setprecision(5) << std::scientific << gresMoSl[2] <<
        "]"  << std::endl;
    std::cout << "Master:  " <<
        " [" << std::setw(14) << std::setprecision(5) << std::scientific << gresMoMa[0] <<
        ", " << std::setw(14) << std::setprecision(5) << std::scientific << gresMoMa[1] <<
        ", " << std::setw(14) << std::setprecision(5) << std::scientific << gresMoMa[2] <<
        "]"  << std::endl;
    std::cout << "Balance: " <<
        " [" << std::setw(14) << std::setprecision(5) << std::scientific << gbalMo[0] <<
        ", " << std::setw(14) << std::setprecision(5) << std::scientific << gbalMo[1] <<
        ", " << std::setw(14) << std::setprecision(5) << std::scientific << gbalMo[2] <<
        "]"  << std::endl;
  }

  if ( conservation_data_ptr )
  {
    LINALG::SerialDenseMatrix& conservation_data = *conservation_data_ptr;
    conservation_data.Zero();
    if (conservation_data.M() < 18 )
      dserror("conservation_data length is too short!");

    std::copy( gresFSl, gresFSl+3, conservation_data.A() );
    std::copy( gresFMa, gresFMa+3, conservation_data.A()+3 );
    std::copy( gbalF, gbalF+3, conservation_data.A()+6 );

    std::copy( gresMoSl, gresMoSl+3, conservation_data.A()+9 );
    std::copy( gresMoMa, gresMoMa+3, conservation_data.A()+12 );
    std::copy( gbalMo, gbalMo+3, conservation_data.A()+15 );
  }

  ++step;
  return;
}


/*---------------------------------------------------------------------------*
 |  Assemble normal coupling weighted condition for poro contact   ager 07/14|
 *--------------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleNCoup(Epetra_Vector& gglobal)
{
  // get out of here if not participating in interface
  if (!lComm()) return;

  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* mrtnode = dynamic_cast<CoNode*>(node);

    if (mrtnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleDMG: Node ownership inconsistency!");

    /**************************************************** nCoup-vector ******/
    if (mrtnode->CoPoroData().GetnCoup()!=0.0)
    {
      double nCoup = mrtnode->CoPoroData().GetnCoup();

      Epetra_SerialDenseVector gnode(1);
      std::vector<int> lm(1);
      std::vector<int> lmowner(1);

      gnode(0) = nCoup;
      lm[0] =  activen_->GID(i);

      lmowner[0] = mrtnode->Owner();

      LINALG::Assemble(gglobal,gnode,lm,lmowner);
    }
  }

  return;
}

/*---------------------------------------------------------------------*
 |  Assemble linearisation of normal coupling                          |
 |          weighted condition for poro contact              ager 07/14|
 *--------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleNCoupLin(LINALG::SparseMatrix& sglobal, ADAPTER::Coupling& coupfs,
    bool AssembleVelocityLin)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // nothing to do if no active nodes
  if (activenodes_==Teuchos::null)
    return;

  Teuchos::RCP<const Epetra_Map> MasterDofMap_full;
  Teuchos::RCP<const Epetra_Map> PermSlaveDofMap_full;

  if (AssembleVelocityLin)
  {
    //store map on all processors, simple but expensive
    MasterDofMap_full = LINALG::AllreduceEMap(*coupfs.MasterDofMap());
    PermSlaveDofMap_full = LINALG::AllreduceEMap(*coupfs.PermSlaveDofMap());
  }

  for (int i=0;i<activenodes_->NumMyElements();++i)
  {
    int gid = activenodes_->GID(i);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      dserror("ERROR: AssembleS: Node ownership inconsistency!");

    std::map<int,double>::iterator colcurr;
    int row = activen_->GID(i);

    std::map<int,double>& dgmap = cnode->CoPoroData().GetDerivnCoup();
    if (AssembleVelocityLin)
    {//Assign fluid velocity linearization to matrix

      dgmap = cnode->CoPoroData().GetVelDerivnCoup();
        for (colcurr=dgmap.begin();colcurr!=dgmap.end();++colcurr)
        {
          int col = PermSlaveDofMap_full->GID(MasterDofMap_full->LID(colcurr->first));
          double val = colcurr->second;

          // do not assemble zeros into s matrix
          if (abs(val)>1.0e-12) sglobal.Assemble(val,row,col);
        }
      //Assign pressure linearization to matrix
      dgmap = cnode->CoPoroData().GetPresDerivnCoup();
        for (colcurr=dgmap.begin();colcurr!=dgmap.end();++colcurr)
        {
          int col = PermSlaveDofMap_full->GID(MasterDofMap_full->LID(colcurr->first))+Dim();
          double val = colcurr->second;

          // do not assemble zeros into s matrix
          if (abs(val)>1.0e-12) sglobal.Assemble(val,row,col);
        }
    }
    else
    {//Assign skeleton displacement linearization to matrix
      for (colcurr=dgmap.begin();colcurr!=dgmap.end();++colcurr)
      {
        int col = colcurr->first;
        double val = colcurr->second;

        // do not assemble zeros into s matrix
        if (abs(val)>1.0e-12) sglobal.Assemble(val,row,col);
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const CONTACT::CoInterface&
CONTACT::CoInterface::GetMaSharingRefInterface() const
{
  return dynamic_cast<const CoInterface&>(idata_.GetMaSharingRefInterface());
}

/*------------------------------------------------------------------------*
 | Derivative of D-matrix multiplied with a slave dof vector              |
 *------------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleCoupLinD(LINALG::SparseMatrix& CoupLin,
                           const Teuchos::RCP<Epetra_Vector> x)
{
  // get out of here if not participating in interface
  if (!lComm())
    return;

  // we have: D_jk,c with j = Slave dof
  //                 with k = Displacement slave dof
  //                 with c = Displacement slave or master dof
  // we compute (LinD)_kc = D_jk,c * x_j

  for (int j=0;j<snoderowmap_->NumMyElements();++j)
  {

    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    // Mortar matrix D derivatives
    std::map<int,std::map<int,double> >& dderiv = cnode->CoData().GetDerivD();

    // get sizes and iterator start
    int slavesize  = (int)dderiv.size();
    std::map<int,std::map<int,double> >::iterator scurr = dderiv.begin();

    /********************************************** LinDMatrix **********/
    // loop over all DISP slave nodes in the DerivD-map of the current slave node
    for (int k=0;k<slavesize;++k)
    {
      int sgid = scurr->first;
      ++scurr;

      DRT::Node* snode = idiscret_->gNode(sgid);
      if (!snode) dserror("ERROR: Cannot find node with gid %",sgid);
      CoNode* csnode = dynamic_cast<CoNode*>(snode);   //current slave node

      // Mortar matrix D derivatives
      std::map<int,double>& thisdderiv = cnode->CoData().GetDerivD()[sgid];
      int mapsize = (int)(thisdderiv.size());

      if(cnode->NumDof() != csnode->NumDof())
        dserror("ERROR: Mortar Nodes on interface must have same number of dofs!");

      // inner product D_{jk,c} * z_j for index j
      for (int prodj=0;prodj<(cnode->NumDof());++prodj)
      {
        int row = csnode->Dofs()[prodj];
        std::map<int,double>::iterator scolcurr = thisdderiv.begin();

        // loop over all directional derivative entries
        for (int c=0;c<mapsize;++c)
        {
          int col = scolcurr->first;
          int slavedofgid = cnode->Dofs()[prodj];

          int slavedoflid = x->Map().LID(slavedofgid);
          if (slavedoflid<0)
            dserror("invalid slave dof lid");
          double val = (*x)[slavedoflid]*(scolcurr->second);

          ++scolcurr;


          //************   ASSEMBLY INTO THE MATRIX    **********************
          if (abs(val)>1.0e-12)
            CoupLin.FEAssemble(val,row,col);
        }
        // check for completeness of DerivD-Derivatives-iteration
        if (scolcurr!=thisdderiv.end())
          dserror("ERROR: AssembleCoupLin: Not all derivative entries of Lin(D*z_s) considered!");
      }
    }

    if (scurr!=dderiv.end())
      dserror("ERROR: AssembleCoupLin: Not all connected slave nodes to Lin(D*z_s) considered!");

  }

  return;
}

/*-----------------------------------------------------------------------------------*
 | Derivative of transposed M-matrix multiplied with a slave dof vector  seitz 01/18 |
 *-----------------------------------------------------------------------------------*/
void CONTACT::CoInterface::AssembleCoupLinM(LINALG::SparseMatrix& CoupLin,
                             const Teuchos::RCP<Epetra_Vector> x)
{
  if (!lComm())
    return;

  // we have: M_jl,c with j = Slave dof
  //                 with l = Displacement master dof
  //                 with c = Displacement slave or master dof
  // we compute (CoupLin)_lc = M_jl,c * x_j

  // loop over all slave nodes (row map)
  for (int j=0;j<snoderowmap_->NumMyElements();++j)
  {
    int gid = snoderowmap_->GID(j);
    DRT::Node* node = idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %",gid);
    CoNode* cnode = dynamic_cast<CoNode*>(node);

    // Mortar matrix M derivatives
    std::map<int,std::map<int,double> >& mderiv = cnode->CoData().GetDerivM();

    // get sizes and iterator start
    int mastersize = (int)mderiv.size();
    std::map<int,std::map<int,double> >::iterator mcurr = mderiv.begin();

    /********************************************** LinMMatrix **********/
    // loop over all master nodes in the DerivM-map of the current LM slave node
    for (int l=0;l<mastersize;++l)
    {
      int mgid = mcurr->first;
      ++mcurr;

      DRT::Node* mnode = idiscret_->gNode(mgid);
      if (!mnode) dserror("ERROR: Cannot find node with gid %",mgid);
      CoNode* cmnode = dynamic_cast<CoNode*>(mnode);

      // Mortar matrix M derivatives
      std::map<int,double>&thismderiv = cnode->CoData().GetDerivM()[mgid];
      int mapsize = (int)(thismderiv.size());

      if(cnode->NumDof() != cmnode->NumDof())
        dserror("ERROR: Mortar Nodes on interface must have same number of dofs!");

      // inner product M_{jl,c} * z_j for index j
      for (int prodj=0;prodj<(cmnode->NumDof());++prodj)
      {
        int row = cmnode->Dofs()[prodj];
        std::map<int,double>::iterator mcolcurr = thismderiv.begin();

        // loop over all directional derivative entries
        for (int c=0;c<mapsize;++c)
        {
          int col = mcolcurr->first;
          int slavedofgid = cnode->Dofs()[prodj];

          int slavedoflid = x->Map().LID(slavedofgid);
          double val = (*x)[slavedoflid]*(mcolcurr->second);

          ++mcolcurr;

          ///************   ASSEMBLY INTO THE MATRIX    **********************
          if (abs(val)>1.0e-12)
            CoupLin.FEAssemble(val,row,col);
        }

        // check for completeness of DerivM-Derivatives-iteration
        if (mcolcurr!=thismderiv.end())
          dserror("ERROR: AssembleCoupLin: Not all derivative entries of DerivM considered!");
      }
    }//loop over all master nodes, connected to this slave node => Lin(M*z_s) finished

    // check for completeness of DerivM-Master-iteration
    if (mcurr!=mderiv.end())
      dserror("ERROR: AssembleCoupLin: Not all master entries of Lin(M*z_s) considered!");

    /*********************************************************************/
    /*******************       Finished with Lin(M*z_s)         **********/
    /*********************************************************************/
  }
}

/*----------------------------------------------------------------------*
 | Store nodal quant. to old ones (last conv. time step)  gitterle 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::CoInterface::StoreToOld(
    MORTAR::StrategyBase::QuantityType type)
{
    // loop over all slave row nodes on the current interface
    for (int j = 0; j < SlaveColNodes()->NumMyElements(); ++j)
    {
      int gid = SlaveColNodes()->GID(j);
      DRT::Node* node = idiscret_->gNode(gid);
      if (!node)
        dserror("ERROR: Cannot find node with gid %", gid);

      switch (type)
      {
      case MORTAR::StrategyBase::dm:
      {
        // store D and M entries
        dynamic_cast<FriNode*>(node)->StoreDMOld();
        break;
      }
      case MORTAR::StrategyBase::pentrac:
      {
        // store penalty tractions to old ones
        dynamic_cast<FriNode*>(node)->StoreTracOld();
        break;
      }
      case MORTAR::StrategyBase::n_old:
      {
        dynamic_cast<CoNode*>(node)->StoreOldNormal();
        break;
      }
      default:
        dserror("ERROR: StoreDMToNodes: Unknown state std::string variable!");
        break;
      } // switch
  }
  return;
}
