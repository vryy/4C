/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_validate.cpp

\brief validate a given .dat-file

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bauer
            089 - 289-15252
</pre>

Validate a given BACI input file (after all preprocessing steps)

*/
/*----------------------------------------------------------------------*/
#ifdef D_EXODUS
#ifdef CCADISCRET

#include "pre_exodus_validate.H"
#include "pre_exodus_soshextrusion.H" //just temporarly for gmsh-plot


/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles; // this CCARAT struct is needed by ReadConditions()


using namespace std;
using namespace Teuchos;
using namespace EXODUS;


void EXODUS::ValidateInputFile(const string datfile)
{
  // read and check the provided header file
  //cout << "checking BACI input file       --> "<<datfile<< endl;

  // do some dirty tricks in order to keep ReadConditions() running
  // (compare with ntainp_ccadiscret() )
  char* datfilename = (char*) datfile.c_str();
  allfiles.inputfile_name = datfilename;
  sprintf(allfiles.outputfile_name, "%s.err",datfile.c_str());
  if ((allfiles.out_err = fopen(allfiles.outputfile_name,"w"))==NULL)
  {
    printf("Opening of output file .err failed\n");
  }

  // communication
#ifdef PARALLEL
  int myrank = 0;
  int nproc  = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if ((nproc>1) && (myrank==0)) dserror("Using more than one processor is not supported.");
  RefCountPtr<Epetra_Comm> comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  RefCountPtr<Epetra_Comm> comm = rcp(new Epetra_SerialComm());
#endif

  //create a problem instance and a DatFileReader
  Teuchos::RCP<DRT::Problem> problem = DRT::Problem::Instance();
  DRT::INPUT::DatFileReader reader(datfile, comm, 0,false);
  reader.Activate();

  // validate dynamic and solver sections
  cout<<"...Read parameters"<<endl;
  problem->ReadParameter(reader);

  // validate all condition definitions
  cout<<"...";
  problem->ReadConditions(reader);

  // materials cannot be checked via problem->ReadMaterial()
  // since the filters use a dummy definition for this method

  // do not read the different fields (discretizations) here,
  // since RAM might be a problem for huge problems!

  // the input file seems to be valid
  cout<<"...OK"<<endl<<endl;

  // clean up
  problem->Done();

  return;
}

void EXODUS::ValidateMeshElementJacobians(Mesh& mymesh)
{
  if (mymesh.GetNumDim() != 3) dserror("Element Validation only for 3 Dimensions");

  map<int,RCP<ElementBlock> > myebs = mymesh.GetElementBlocks();
  map<int,RCP<ElementBlock> >::iterator i_eb;

  for(i_eb=myebs.begin(); i_eb!=myebs.end(); ++i_eb){
    RCP<ElementBlock> eb = i_eb->second;
    const DRT::Element::DiscretizationType distype = PreShapeToDrt(eb->GetShape());
    ValidateElementJacobian(mymesh,distype,eb);
    // full check at all gausspoints
    int invalid_dets = ValidateElementJacobian_fullgp(mymesh,distype,eb);
    if (invalid_dets > 0) cout << invalid_dets << " negative Jacobian determinants in EB of shape " << ShapeToString(eb->GetShape()) << endl;
  }
  return;
}

void EXODUS::ValidateElementJacobian(Mesh& mymesh, const DRT::Element::DiscretizationType distype, RCP<ElementBlock> eb)
{
  // use one point gauss rule to calculate jacobian at element center
  DRT::UTILS::GaussRule3D integrationrule_1point = DRT::UTILS::intrule3D_undefined;
  switch(distype)
  {
  case DRT::Element::hex8: case DRT::Element::hex20: case DRT::Element::hex27:
      integrationrule_1point = DRT::UTILS::intrule_hex_1point;
      break;
  case DRT::Element::tet4: case DRT::Element::tet10:
      integrationrule_1point = DRT::UTILS::intrule_tet_1point;
      break;
  case DRT::Element::wedge6: case DRT::Element::wedge15:
      integrationrule_1point = DRT::UTILS::intrule_wedge_1point;
      break;
  case DRT::Element::pyramid5:
      integrationrule_1point = DRT::UTILS::intrule_pyramid_1point;
      break;
  default:
      // cout<<"No Validation for this kind of Element implemented! Good luck!"<<endl;
      return;
  }
  const DRT::UTILS::IntegrationPoints3D  intpoints = getIntegrationPoints3D(integrationrule_1point);
  const int iel = eb->GetEleNodes(0).size();
  // shape functions derivatives
  const int NSD = 3;
  Epetra_SerialDenseMatrix    deriv(NSD, iel);
  DRT::UTILS::shape_function_3D_deriv1(deriv,intpoints.qxg[0][0],intpoints.qxg[0][1],intpoints.qxg[0][2],distype);

  // go through all elements
  RCP<map<int,vector<int> > > eleconn = eb->GetEleConn();
  map<int,vector<int> >::iterator i_ele;
  for(i_ele=eleconn->begin();i_ele!=eleconn->end();++i_ele){
    if(!PositiveEle(i_ele->second,mymesh,deriv)){
      i_ele->second = RewindEle(i_ele->second,distype);
      // double check
      if(!PositiveEle(i_ele->second,mymesh,deriv)) dserror("Could not determine a proper rewinding");
    }
  }

  return;
}

int EXODUS::ValidateElementJacobian_fullgp(Mesh& mymesh, const DRT::Element::DiscretizationType distype, RCP<ElementBlock> eb)
{
  DRT::UTILS::GaussRule3D integrationrule = DRT::UTILS::intrule3D_undefined;
  switch(distype)
  {
  case DRT::Element::hex8:
      integrationrule = DRT::UTILS::intrule_hex_8point;
      break;
  case DRT::Element::hex20:
      integrationrule = DRT::UTILS::intrule_hex_27point;
      break;
  case DRT::Element::hex27:
      integrationrule = DRT::UTILS::intrule_hex_27point;
      break;
  case DRT::Element::tet4:
      integrationrule = DRT::UTILS::intrule_tet_4point;
      break;
  case DRT::Element::tet10:
      integrationrule = DRT::UTILS::intrule_tet_10point;
      break;
  case DRT::Element::wedge6: case DRT::Element::wedge15:
      integrationrule = DRT::UTILS::intrule_wedge_6point;
      break;
  case DRT::Element::pyramid5:
      integrationrule = DRT::UTILS::intrule_pyramid_8point;
      break;
  case DRT::Element::line2:
      return 0;
      break;
  default:
    integrationrule = DRT::UTILS::intrule3D_undefined;
    break;
  }
  const DRT::UTILS::IntegrationPoints3D  intpoints = getIntegrationPoints3D(integrationrule);
  const int iel = eb->GetEleNodes(0).size();
  // shape functions derivatives
  const int NSD = 3;
  Epetra_SerialDenseMatrix    deriv(NSD, iel);

  // go through all elements
  int invalids = 0;
  RCP<map<int,vector<int> > > eleconn = eb->GetEleConn();
  map<int,vector<int> >::iterator i_ele;
  for(i_ele=eleconn->begin();i_ele!=eleconn->end();++i_ele){
    for (int igp = 0; igp < intpoints.nquad; ++igp) {
      DRT::UTILS::shape_function_3D_deriv1(deriv,intpoints.qxg[igp][0],intpoints.qxg[igp][1],intpoints.qxg[igp][2],distype);
      if (PositiveEle(i_ele->second,mymesh,deriv) == false){
        invalids++;
      }
    }
  }

  return invalids;
}


bool EXODUS::PositiveEle(const vector<int>& nodes,const Mesh& mymesh,const Epetra_SerialDenseMatrix& deriv)
{
  const int iel = deriv.N();
  const int NSD = deriv.M();
  LINALG::SerialDenseMatrix xyze(deriv.M(),iel);
  for (int inode=0; inode<iel; inode++)
  {
    const vector<double> x = mymesh.GetNode(nodes.at(inode));
    xyze(0,inode) = x[0];
    xyze(1,inode) = x[1];
    xyze(2,inode) = x[2];
  }
  // get Jacobian matrix and determinant
  // actually compute its transpose....
  LINALG::SerialDenseMatrix xjm(NSD,NSD);
  xjm.Multiply('N','T',1.0,deriv,xyze,0.0);
  const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                     xjm(0,1)*xjm(1,2)*xjm(2,0)+
                     xjm(0,2)*xjm(1,0)*xjm(2,1)-
                     xjm(0,2)*xjm(1,1)*xjm(2,0)-
                     xjm(0,0)*xjm(1,2)*xjm(2,1)-
                     xjm(0,1)*xjm(1,0)*xjm(2,2);
  if (abs(det) < 1E-16) dserror("ZERO JACOBIAN DETERMINANT");

  if (det < 0.0) return false;
  else return true;
}

int EXODUS::EleSaneSign(const vector<int>& nodes,const map<int,vector<double> >& nodecoords)
{
  const int iel = nodes.size();
  // to be even stricter we test the Jacobian at every Node, not just at the gausspoints
  LINALG::SerialDenseMatrix local_nodecoords(iel,3);
  DRT::Element::DiscretizationType distype;
  switch(iel)
  {
  case 8: // hex8
    local_nodecoords(0,0) = -1.; local_nodecoords(0,1) = -1.; local_nodecoords(0,2) = -1.;
    local_nodecoords(1,0) =  1.; local_nodecoords(1,1) = -1.; local_nodecoords(1,2) = -1.;
    local_nodecoords(2,0) =  1.; local_nodecoords(2,1) =  1.; local_nodecoords(2,2) = -1.;
    local_nodecoords(3,0) = -1.; local_nodecoords(3,1) =  1.; local_nodecoords(3,2) = -1.;
    local_nodecoords(4,0) = -1.; local_nodecoords(4,1) = -1.; local_nodecoords(4,2) =  1.;
    local_nodecoords(5,0) =  1.; local_nodecoords(5,1) = -1.; local_nodecoords(5,2) =  1.;
    local_nodecoords(6,0) =  1.; local_nodecoords(6,1) =  1.; local_nodecoords(6,2) =  1.;
    local_nodecoords(7,0) = -1.; local_nodecoords(7,1) =  1.; local_nodecoords(7,2) =  1.;
    distype = DRT::Element::hex8;
    break;
  case 6: // wedge6
    local_nodecoords(0,0) = 0.; local_nodecoords(0,1) = 0.; local_nodecoords(0,2) = -1.;
    local_nodecoords(1,0) = 1.; local_nodecoords(1,1) = 0.; local_nodecoords(1,2) = -1.;
    local_nodecoords(2,0) = 0.; local_nodecoords(2,1) = 1.; local_nodecoords(2,2) = -1.;
    local_nodecoords(3,0) = 0.; local_nodecoords(3,1) = 0.; local_nodecoords(3,2) =  1.;
    local_nodecoords(4,0) = 1.; local_nodecoords(4,1) = 0.; local_nodecoords(4,2) =  1.;
    local_nodecoords(5,0) = 0.; local_nodecoords(5,1) = 1.; local_nodecoords(5,2) =  1.;
    distype = DRT::Element::wedge6;
    break;
  default:
    dserror("No Element Sanity Check for this distype");
    break;
  }
  // shape functions derivatives
  const int NSD = 3;
  Epetra_SerialDenseMatrix    deriv(NSD, iel);

  LINALG::SerialDenseMatrix xyze(deriv.M(),iel);
  for (int inode=0; inode<iel; inode++)
  {
    const vector<double> x = nodecoords.find(nodes[inode])->second;
    xyze(0,inode) = x[0];
    xyze(1,inode) = x[1];
    xyze(2,inode) = x[2];
  }
  // get Jacobian matrix and determinant
  // actually compute its transpose....
  LINALG::SerialDenseMatrix xjm(NSD,NSD);
  int n_posdet = 0;
  int n_negdet = 0;

  for (int i = 0; i < iel; ++i) {
    DRT::UTILS::shape_function_3D_deriv1(deriv,local_nodecoords(i,0),local_nodecoords(i,1),local_nodecoords(i,2),distype);
    xjm.Multiply('N','T',1.0,deriv,xyze,0.0);
    const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                       xjm(0,1)*xjm(1,2)*xjm(2,0)+
                       xjm(0,2)*xjm(1,0)*xjm(2,1)-
                       xjm(0,2)*xjm(1,1)*xjm(2,0)-
                       xjm(0,0)*xjm(1,2)*xjm(2,1)-
                       xjm(0,1)*xjm(1,0)*xjm(2,2);
    if (abs(det) < 1E-16) dserror("ZERO JACOBIAN DETERMINANT");
    if (det<0){
      ++n_negdet;
    }
    else ++n_posdet;
  }

  if (n_posdet==iel && n_negdet==0) return 1;
  else if (n_posdet==0 && n_negdet==iel) return -1;
  else return 0;

  return 0;
}

vector<int> EXODUS::RewindEle(vector<int> old_nodeids, const DRT::Element::DiscretizationType distype)
{
  vector<int> new_nodeids(old_nodeids.size());
  // rewinding of nodes to arrive at mathematically positive element
  switch(distype)
  {
  case DRT::Element::tet4:{
    new_nodeids[0] = old_nodeids[0];
    new_nodeids[1] = old_nodeids[2];
    new_nodeids[2] = old_nodeids[1];
    new_nodeids[3] = old_nodeids[3];
    break;
  }
  case DRT::Element::tet10:{
    new_nodeids[0] = old_nodeids[0];
    new_nodeids[1] = old_nodeids[2];
    new_nodeids[2] = old_nodeids[1];
    new_nodeids[3] = old_nodeids[3];
    new_nodeids[4] = old_nodeids[6];
    new_nodeids[5] = old_nodeids[5];
    new_nodeids[6] = old_nodeids[4];
    new_nodeids[7] = old_nodeids[7];
    new_nodeids[8] = old_nodeids[8];
    new_nodeids[9] = old_nodeids[9];
    break;
  }
  case DRT::Element::hex8:{
    new_nodeids[0] = old_nodeids[4];
    new_nodeids[1] = old_nodeids[5];
    new_nodeids[2] = old_nodeids[6];
    new_nodeids[3] = old_nodeids[7];
    new_nodeids[4] = old_nodeids[0];
    new_nodeids[5] = old_nodeids[1];
    new_nodeids[6] = old_nodeids[2];
    new_nodeids[7] = old_nodeids[3];
    break;
  }
  case DRT::Element::wedge6:{
    new_nodeids[0] = old_nodeids[3];
    new_nodeids[1] = old_nodeids[4];
    new_nodeids[2] = old_nodeids[5];
    new_nodeids[3] = old_nodeids[0];
    new_nodeids[4] = old_nodeids[1];
    new_nodeids[5] = old_nodeids[2];
    break;
  }
  case DRT::Element::pyramid5:{
    new_nodeids[1] = old_nodeids[3];
    new_nodeids[3] = old_nodeids[1];
    // the other nodes can stay the same
    new_nodeids[0] = old_nodeids[0];
    new_nodeids[2] = old_nodeids[2];
    new_nodeids[4] = old_nodeids[4];
    break;
  }
  default: dserror("no rewinding scheme for this type of element");
  }
  return new_nodeids;
}

#endif
#endif
