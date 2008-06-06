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


void EXODUS::ValidateInputFile(string datfile)
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
  DRT::INPUT::DatFileReader reader(datfile.c_str(), comm, 0,false);
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
    switch(distype)
    {
    case DRT::Element::hex8: 
      ValidateElementJacobian(mymesh,distype,eb); break;
    case DRT::Element::tet4: case DRT::Element::tet10:
      ValidateElementJacobian(mymesh,distype,eb); break;
    case DRT::Element::wedge6: 
      ValidateElementJacobian(mymesh,distype,eb); break;
    case DRT::Element::pyramid5:
      ValidateElementJacobian(mymesh,distype,eb); break;
    default:
        cout << "Warning: No ElementJacobian Validation for this distype: " << DRT::DistypeToString(distype) << endl;
    }
  }
  return;
}

void EXODUS::ValidateElementJacobian(Mesh& mymesh, const DRT::Element::DiscretizationType distype, RCP<ElementBlock> eb)
{
  // use one point gauss rule to calculate tau at element center
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
      dserror("invalid discretization type");
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
    }
  }

  return;
}

bool EXODUS::PositiveEle(vector<int>& nodes,Mesh& mymesh, Epetra_SerialDenseMatrix& deriv)
{
  const int iel = deriv.N();
  const int NSD = deriv.M();
  LINALG::SerialDenseMatrix xyze(deriv.M(),iel);
  for (int inode=0; inode<iel; inode++)
  {
    const vector<double> x = mymesh.GetNodeExo(nodes.at(inode));
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
  
  bool rewind = false;
  if (det < 0.0) rewind = true;
  else rewind = false;
 
  return rewind;
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
