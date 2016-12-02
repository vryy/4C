/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_validate.cpp

\brief validate a given .dat-file

\maintainer Martin Kronbichler

\level 1

Validate a given BACI input file (after all preprocessing steps)

*/
/*----------------------------------------------------------------------*/

#include "pre_exodus_validate.H"
#include "pre_exodus_soshextrusion.H" //just temporarly for gmsh-plot
#include "../drt_io/io_control.H"     //for writing to the error file
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_inputreader.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::ValidateInputFile(const Teuchos::RCP<Epetra_Comm> comm, const std::string datfile)
{
  // read and check the provided header file
  //std::cout << "checking BACI input file       --> "<<datfile<< std::endl;

  // access our problem instance
  DRT::Problem* problem = DRT::Problem::Instance();

  // create a DatFileReader
  DRT::INPUT::DatFileReader reader(datfile,comm, 0);

  // read and validate dynamic and solver sections
  std::cout<<"...Read parameters"<<std::endl;
  problem->ReadParameter(reader);

  // read and validate all material definitions
  std::cout<<"...Read materials"<<std::endl;
  problem->ReadMaterials(reader);

  // do NOT allocate the different fields (discretizations) here,
  // since RAM might be a problem for huge problems!
  // But, we have to perform at least the problem-specific setup since
  // some reading procedures depend on the number of fields (e.g., ReadKnots())
  std::cout<<"...Read field setup"<<std::endl;
  problem->ReadFields(reader,false);  // option false is important here!

  // read and validate all condition definitions
  std::cout<<"...";
  problem->ReadTimeFunctionResult(reader);
  problem->ReadConditions(reader);

  // input of materials of cloned fields (if needed)
  problem->ReadCloningMaterialMap(reader);

  // read all knot information for isogeometric analysis
  // and add it to the (derived) nurbs discretization
  problem->ReadKnots(reader);

  // inform user about unused/obsolete section names being found
  // and force him/her to correct the input file accordingly
  reader.PrintUnknownSections();

  // the input file seems to be valid
  std::cout<<"...OK"<<std::endl<<std::endl;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::ValidateMeshElementJacobians(Mesh& mymesh)
{
  if (mymesh.GetNumDim() != 3) dserror("Element Validation only for 3 Dimensions");

  std::map<int,Teuchos::RCP<ElementBlock> > myebs = mymesh.GetElementBlocks();
  std::map<int,Teuchos::RCP<ElementBlock> >::iterator i_eb;

  for(i_eb=myebs.begin(); i_eb!=myebs.end(); ++i_eb){
    Teuchos::RCP<ElementBlock> eb = i_eb->second;
    const DRT::Element::DiscretizationType distype = PreShapeToDrt(eb->GetShape());
    // check and rewind if necessary
    ValidateElementJacobian(mymesh,distype,eb);
    // full check at all gausspoints
    int invalid_dets = ValidateElementJacobian_fullgp(mymesh,distype,eb);
    if (invalid_dets > 0) std::cout << invalid_dets << " negative Jacobian determinants in EB of shape " << ShapeToString(eb->GetShape()) << std::endl;
  }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EXODUS::ValidateElementJacobian(Mesh& mymesh, const DRT::Element::DiscretizationType distype, Teuchos::RCP<ElementBlock> eb)
{
  // use one point gauss rule to calculate jacobian at element center
  DRT::UTILS::GaussRule3D integrationrule_1point = DRT::UTILS::intrule3D_undefined;
  switch(distype)
  {
  case DRT::Element::hex8: case DRT::Element::hex20:
      integrationrule_1point = DRT::UTILS::intrule_hex_1point;
      break;
  case DRT::Element::hex27:
    integrationrule_1point = DRT::UTILS::intrule_hex_27point;  // one point is not enough for hex27!!
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
  // do nothing for 2D, 1D and 0D elements
  case DRT::Element::quad4: case DRT::Element::quad8: case DRT::Element::quad9:
  case DRT::Element::tri3  : case DRT::Element::tri6:
  case DRT::Element::line2: case DRT::Element::line3:
  case DRT::Element::point1:
      return;
  default:
      dserror("Unknown element type, validation failed!");
      break;
  }
  const DRT::UTILS::IntegrationPoints3D  intpoints(integrationrule_1point);
  const int iel = eb->GetEleNodes(0).size();
  // shape functions derivatives
  const int NSD = 3;
  Epetra_SerialDenseMatrix    deriv(NSD, iel);

  // go through all elements
  Teuchos::RCP<std::map<int,std::vector<int> > > eleconn = eb->GetEleConn();
  std::map<int,std::vector<int> >::iterator i_ele;
  int numrewindedeles=0;
  for(i_ele=eleconn->begin();i_ele!=eleconn->end();++i_ele){
    int rewcount=0;
    for (int igp = 0; igp < intpoints.nquad; ++igp) {
      DRT::UTILS::shape_function_3D_deriv1(deriv,intpoints.qxg[igp][0],intpoints.qxg[igp][1],intpoints.qxg[igp][2],distype);
      if(!PositiveEle(i_ele->first,i_ele->second,mymesh,deriv))
      {
        // rewind the element nodes
        if (rewcount==0)
        {
          i_ele->second = RewindEle(i_ele->second,distype);
          // write info to error file
          FILE* errfile = DRT::Problem::Instance()->ErrorFile()->Handle();
          fprintf(errfile,"GAUSS POINT %d: REWINDED ELEMENT %d\n",igp, i_ele->first);
          fflush(errfile);
          numrewindedeles++;
        }
        // double check
        if(!PositiveEle(i_ele->first,i_ele->second,mymesh,deriv))
          dserror("No proper rewinding for element id %d at gauss point %d",i_ele->first,igp);
        rewcount++;
      }
    }
  }
  if (numrewindedeles>0)
    std::cout<<"...Successfully rewinded "<<numrewindedeles<<" elements. For details see *.err file"<<std::endl;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int EXODUS::ValidateElementJacobian_fullgp(Mesh& mymesh, const DRT::Element::DiscretizationType distype, Teuchos::RCP<ElementBlock> eb)
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
  // do nothing for 2D, 1D and 0D elements
  case DRT::Element::quad4: case DRT::Element::quad8: case DRT::Element::quad9:
  case DRT::Element::tri3  : case DRT::Element::tri6:
  case DRT::Element::line2: case DRT::Element::line3:
  case DRT::Element::point1:
      return 0;
  default:
      dserror("Unknown element type, validation failed!");
      break;
  }
  const DRT::UTILS::IntegrationPoints3D  intpoints(integrationrule);
  const int iel = eb->GetEleNodes(0).size();
  // shape functions derivatives
  const int NSD = 3;
  Epetra_SerialDenseMatrix    deriv(NSD, iel);

  // go through all elements
  int invalids = 0;
  Teuchos::RCP<std::map<int,std::vector<int> > > eleconn = eb->GetEleConn();
  std::map<int,std::vector<int> >::iterator i_ele;
  for(i_ele=eleconn->begin();i_ele!=eleconn->end();++i_ele){
    for (int igp = 0; igp < intpoints.nquad; ++igp) {
      DRT::UTILS::shape_function_3D_deriv1(deriv,intpoints.qxg[igp][0],intpoints.qxg[igp][1],intpoints.qxg[igp][2],distype);
      if (PositiveEle(i_ele->first,i_ele->second,mymesh,deriv) == false){
        invalids++;
      }
    }
  }

  return invalids;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool EXODUS::PositiveEle(const int& eleid, const std::vector<int>& nodes,const Mesh& mymesh,const Epetra_SerialDenseMatrix& deriv)
{
  const int iel = deriv.N();
  const int NSD = deriv.M();
  LINALG::SerialDenseMatrix xyze(deriv.M(),iel);
  for (int inode=0; inode<iel; inode++)
  {
    const std::vector<double> x = mymesh.GetNode(nodes.at(inode));
    xyze(0,inode) = x[0];
    xyze(1,inode) = x[1];
    xyze(2,inode) = x[2];
  }
  // get Jacobian matrix and determinant
  // actually compute its transpose....
  if (NSD==3)
  {
    LINALG::SerialDenseMatrix xjm(NSD,NSD);
    xjm.Multiply('N','T',1.0,deriv,xyze,0.0);
    LINALG::Matrix<3,3> jac(xjm.A(),true);
    const double det =  jac.Determinant();

    if (abs(det) < 1E-16)
        dserror("ZERO JACOBIAN DETERMINANT FOR ELEMENT %d: DET = %f",eleid,det);

    if (det < 0.0)
    {
      // write info to the error log file
      FILE* errfile = DRT::Problem::Instance()->ErrorFile()->Handle();
      fprintf(errfile,"NEGATIVE JACOBIAN DETERMINANT FOR ELEMENT %d: DET = %f -> REWIND\n",eleid,det);
      fflush(errfile);
      return false;
    }
  }
  return true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int EXODUS::EleSaneSign(const std::vector<int>& nodes,const std::map<int,std::vector<double> >& nodecoords)
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
    const std::vector<double> x = nodecoords.find(nodes[inode])->second;
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

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<int> EXODUS::RewindEle(std::vector<int> old_nodeids, const DRT::Element::DiscretizationType distype)
{
  std::vector<int> new_nodeids(old_nodeids.size());
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
  case DRT::Element::hex27:{
    // nodes 1 - 20 can stay the same (no rewinding for hex20)
    for (int i = 0; i < 20; i++)
    {
      new_nodeids[i] = old_nodeids[i];
    }
    // rewind the nodes on the center of the 6 sides
    // and the center node of the actual hex27 element
    new_nodeids[20] = old_nodeids[21];
    new_nodeids[21] = old_nodeids[25];
    new_nodeids[22] = old_nodeids[24];
    new_nodeids[23] = old_nodeids[26];
    new_nodeids[24] = old_nodeids[23];
    new_nodeids[25] = old_nodeids[22];
    new_nodeids[26] = old_nodeids[20];
    break;
  }
  default: dserror("no rewinding scheme for this type of element");
  break;
  }
  return new_nodeids;
}
