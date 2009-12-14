/*!
\file io_gmsh.cpp

\brief simple element print library for Gmsh

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/

#ifdef CCADISCRET

#include "io_gmsh.H"
#include "io_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/linalg_utils.H"  // LINALG::Export
#include "../drt_geometry/intersection_service.H"


std::string IO::GMSH::distypeToGmshElementHeader(
    const DRT::Element::DiscretizationType distype)
{
  switch (distype)
  {
    case DRT::Element::hex8:    return "H";  break;
    case DRT::Element::hex20:   return "H";  break;
    case DRT::Element::hex27:   return "H";  break;
    case DRT::Element::tet4:    return "S";  break;
    case DRT::Element::tet10:   return "S";  break;
    case DRT::Element::point1:  return "P";  break;
    case DRT::Element::quad4:   return "Q";  break;
    case DRT::Element::quad8:   return "Q";  break;
    case DRT::Element::quad9:   return "Q";  break;
    case DRT::Element::tri3:    return "T";  break;
    case DRT::Element::tri6:    return "T";  break;
    case DRT::Element::line2:   return "L";  break;
    case DRT::Element::line3:   return "L2"; break;
    default:
      dserror("distypeToGmshElementHeader: distype not supported for printout!");
  }
  return "xxx";
}

int IO::GMSH::distypeToGmshNumNode(
    const DRT::Element::DiscretizationType distype   ///< element shape
)
{
  switch (distype)
  {
    case DRT::Element::hex8:   return 8;  break;
    case DRT::Element::hex20:  return 8;  break;
    case DRT::Element::hex27:  return 8;  break;
    case DRT::Element::tet4:   return 4;  break;
    case DRT::Element::tet10:  return 4;  break;
    case DRT::Element::point1: return 1;  break;
    case DRT::Element::quad4:  return 4;  break;
    case DRT::Element::quad8:  return 4;  break;
    case DRT::Element::quad9:  return 4;  break;
    case DRT::Element::tri3:   return 3;  break;
    case DRT::Element::tri6:   return 3;  break;
    case DRT::Element::line2:  return 2;  break;
    case DRT::Element::line3:  return 3;  break;
    default:
      dserror("distypeToGmshNumNode: distype not supported for printout!");
  }
  return -1;
}

/*------------------------------------------------------------------------------------------------*
 | write scalar field to Gmsh postprocessing file                                     henke 12/09 |
 *------------------------------------------------------------------------------------------------*/
void IO::GMSH::ScalarFieldToGmsh(
    const Teuchos::RCP<DRT::Discretization> discret,
    const Teuchos::RCP<Epetra_Vector>       scalarfield_row,
    std::ostream&                           s
)
{
#ifdef PARALLEL
  // tranform solution vector from DofRowMap to DofColMap
  const Teuchos::RCP<const Epetra_Vector> scalarfield = DRT::UTILS::GetColVersionOfRowVector(discret,scalarfield_row);
#else
  const Teuchos::RCP<const Epetra_Vector> scalarfield = scalarfield_row;
#endif

  // loop all row elements on this processor
  for (int iele=0; iele<discret->NumMyRowElements(); ++iele)
  {
    const DRT::Element* ele = discret->lRowElement(iele);
    const DRT::Element::DiscretizationType distype = ele->Shape();
    const int numnode = distypeToGmshNumNode(distype);

    LINALG::SerialDenseMatrix xyze(3,numnode);

    const DRT::Node*const* nodes = ele->Nodes();
    for(int inode = 0; inode < numnode; ++inode)
    {
      xyze(0,inode) = nodes[inode]->X()[0];
      xyze(1,inode) = nodes[inode]->X()[1];
      xyze(2,inode) = nodes[inode]->X()[2];
    }

    s << "S"; // scalar field indicator
    s << distypeToGmshElementHeader(distype);

    // write node coordinates to Gmsh stream
    CoordinatesToStream(xyze, distype, s);

    std::vector<int> lm;
    std::vector<int> lmowner;
    ele->LocationVector(*discret, lm, lmowner);

    // extract local values from the global vector
    Epetra_SerialDenseVector myscalarfield(lm.size());
    DRT::UTILS::ExtractMyValues(*scalarfield, myscalarfield, lm);

    // write scalar field to Gmsh stream
    ScalarFieldToStream(myscalarfield, distype, s);

    s << "\n";
  }
}

/*------------------------------------------------------------------------------------------------*
 | write dof-based vector field to Gmsh postprocessing file                           henke 12/09 |
 *------------------------------------------------------------------------------------------------*/
void IO::GMSH::VectorFieldDofBasedToGmsh(
    const Teuchos::RCP<DRT::Discretization> discret,
    const Teuchos::RCP<Epetra_Vector>       vectorfield_row,
    std::ostream&                           s
)
{
#ifdef PARALLEL
  // tranform solution vector from DofRowMap to DofColMap
  const Teuchos::RCP<const Epetra_Vector> vectorfield = DRT::UTILS::GetColVersionOfRowVector(discret,vectorfield_row);
#else
  const Teuchos::RCP<const Epetra_Vector> vectorfield = vectorfield_row;
#endif

  // loop all row elements on this processor
  for (int iele=0; iele<discret->NumMyRowElements(); ++iele)
  {
    const DRT::Element* ele = discret->lRowElement(iele);
    const DRT::Element::DiscretizationType distype = ele->Shape();
    const int numnode = distypeToGmshNumNode(distype);

    LINALG::SerialDenseMatrix xyze(3,numnode);

    const DRT::Node*const* nodes = ele->Nodes();
    for(int inode = 0; inode < numnode; ++inode)
    {
      xyze(0,inode) = nodes[inode]->X()[0];
      xyze(1,inode) = nodes[inode]->X()[1];
      xyze(2,inode) = nodes[inode]->X()[2];
    }

    s << "V"; // vector field indicator
    s << distypeToGmshElementHeader(distype);

    // write node coordinates to Gmsh stream
    CoordinatesToStream(xyze, distype, s);

    std::vector<int> lm;
    std::vector<int> lmowner;
    ele->LocationVector(*discret, lm, lmowner);

    // extract local values from the global vector
    Epetra_SerialDenseVector myvectorfield(lm.size());
    DRT::UTILS::ExtractMyValues(*vectorfield, myvectorfield, lm);

    // write vector field to Gmsh stream
    // remark: only the first 3 components are written to the Gmsh postprocessing file
    //         -> e.g. pressure in velnp_ fluid state vector is not written ('myvectorfield': velx,vely,velz,pressure)
    VectorFieldToStream(myvectorfield, distype, s);

    s << "\n";
  }
}

/*------------------------------------------------------------------------------------------------*
 | write node-based vector field to Gmsh postprocessing file                          henke 12/09 |
 *------------------------------------------------------------------------------------------------*/
void IO::GMSH::VectorFieldNodeBasedToGmsh(
    const Teuchos::RCP<DRT::Discretization> discret,
    const Teuchos::RCP<Epetra_MultiVector>  vectorfield_row,
    std::ostream&                           s
)
{
#ifdef PARALLEL
  // tranform solution vector from NodeRowMap to NodeColMap
  // remark: DRT::UTILS::GetColVersionOfRowVector() does only work for Epetra_Vectors on DofRowMap
  //         something similar is done in COMBUST::FlameFront::ProcessFlameFront, although not for Epetra_MultiVectors
  const Teuchos::RCP<Epetra_MultiVector> vectorfield = rcp(new Epetra_MultiVector(*discret->NodeColMap(),3,true));
  LINALG::Export(*vectorfield_row,*vectorfield);
#else
  const Teuchos::RCP<Epetra_MultiVector> vectorfield = vectorfield_row;
#endif

  // loop all row elements on this processor
  for (int iele=0; iele<discret->NumMyRowElements(); ++iele)
  {
    const DRT::Element* ele = discret->lRowElement(iele);
    const DRT::Element::DiscretizationType distype = ele->Shape();
    const int numnode = distypeToGmshNumNode(distype);

    LINALG::SerialDenseMatrix xyze(3,numnode);

    const DRT::Node*const* nodes = ele->Nodes();
    for(int inode = 0; inode < numnode; ++inode)
    {
      xyze(0,inode) = nodes[inode]->X()[0];
      xyze(1,inode) = nodes[inode]->X()[1];
      xyze(2,inode) = nodes[inode]->X()[2];
    }

    s << "V"; // vector field indicator
    s << distypeToGmshElementHeader(distype);

    // write node coordinates to Gmsh stream
    CoordinatesToStream(xyze, distype, s);

    // extract local values from the global vector
    Epetra_SerialDenseMatrix myvectorfield(3,numnode);
    DRT::UTILS::ExtractMyNodeBasedValues(ele, myvectorfield, vectorfield, 3);

    // write vector field to Gmsh stream
    VectorFieldToStream(myvectorfield, distype, s);

    s << "\n";
  }
}

void IO::GMSH::ScalarToStream(
    const double                           scalar,
    const DRT::Element::DiscretizationType distype,
    std::ostream&                          s
    )
{
  s.setf(ios::scientific,ios::floatfield);
  s.precision(12);

  const int numnode = distypeToGmshNumNode(distype);

  // values
  s << "{";
  for (int i = 0; i<numnode; ++i)
  {
    s << scalar;
    if (i < numnode-1)
    {
      s << ",";
    }
  };
  s << "};";
}

void IO::GMSH::elementAtInitialPositionToStream(
    const double scalar,
    const DRT::Element* ele,
    std::ostream& s)
{
  const DRT::Node*const* nodes = ele->Nodes();

  const DRT::Element::DiscretizationType distype = ele->Shape();
  const int numnode = distypeToGmshNumNode(distype);

  s.setf(ios::scientific,ios::floatfield);
  s.precision(12);

  s << "S" << distypeToGmshElementHeader(distype) << "(";
  for (int i = 0; i<numnode; ++i)
  {
    const DRT::Node* node = nodes[i];
    const double* x = node->X();
    s << x[0] << ",";
    s << x[1] << ",";
    s << x[2];
    if (i < numnode-1)
    {
      s << ",";
    }
  };
  s << ")";
  // values
  ScalarToStream(scalar, distype, s);
  s << "\n";
}


std::string IO::GMSH::elementAtInitialPositionToString(
    const double scalar,
    const DRT::Element* ele)
{
  std::ostringstream s;
  elementAtInitialPositionToStream(scalar, ele, s);
  return s.str();
}


void IO::GMSH::elementAtCurrentPositionToStream(
    const double                            scalar,
    const DRT::Element*                     ele,
    const map<int,LINALG::Matrix<3,1> >&    currentelepositions,
    std::ostream&                           s
    )
{
  IO::GMSH::cellWithScalarToStream(
      ele->Shape(),
      scalar,
      GEO::getCurrentNodalPositions(ele,currentelepositions),
      s);
}


std::string IO::GMSH::elementAtCurrentPositionToString(
    const double                            scalar,
    const DRT::Element*                     ele,
    const map<int,LINALG::Matrix<3,1> >&    currentelepositions)
{
  std::ostringstream s;
  IO::GMSH::elementAtCurrentPositionToStream(
      scalar,
      ele,
      currentelepositions,
      s);
  return s.str();
}


std::string IO::GMSH::text3dToString(
    const LINALG::Matrix<3,1>&            xyz,      ///< 3d Position of text
    const std::string&                    text,     ///< text to be printed
    const int                             fontsize  ///< font size
    )
{
  std::ostringstream s;

  s << "T3";
  // coordinates
  s << "(";
  s << scientific << xyz(0) <<",";
  s << scientific << xyz(1) <<",";
  s << scientific << xyz(2) <<",";
  s << fontsize << ")";
  s << "{\"" << text <<"\"};";
  s << "\n";
  return s.str();
}

void IO::GMSH::disToStream(
    const std::string&                      text,
    const double                            scalar,
    const Teuchos::RCP<DRT::Discretization> dis,
    std::ostream&                           s
    )
{
  s << "View \" " << text << " Elements \" {\n";
  for (int i=0; i<dis->NumMyRowElements(); ++i)
  {
    const DRT::Element* actele = dis->lRowElement(i);
    IO::GMSH::elementAtInitialPositionToStream(scalar, actele, s);
  };
  s << "};\n";
}

std::string IO::GMSH::disToString(
    const std::string& text,
    const double scalar,
    const Teuchos::RCP<DRT::Discretization> dis)
{
  std::ostringstream s;
  disToStream(text, scalar, dis, s);
  return s.str();
}

void IO::GMSH::disToStream(
    const std::string&                          text,
    const double                                scalar,
    const Teuchos::RCP<DRT::Discretization>     dis,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    std::ostream&                               s)
{
  s << "View \" " << text << " Elements \" {\n";

  for (int i=0; i<dis->NumMyColElements(); ++i)
  {
    const DRT::Element* actele = dis->lColElement(i);
    IO::GMSH::cellWithScalarToStream(actele->Shape(),
        scalar, GEO::getCurrentNodalPositions(actele,currentpositions), s);
  };
  s << "};\n";
}

std::string IO::GMSH::disToString(
    const std::string&                          text,
    const double                                scalar,
    const Teuchos::RCP<DRT::Discretization>     dis,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions)
{
  std::ostringstream s;
  disToStream(text, scalar, dis, currentpositions, s);
  return s.str();
}

std::string IO::GMSH::GetNewFileNameAndDeleteOldFiles(
    const std::string&   filename_base,
    const int&           actstep,           ///< generate filename for this step
    const int&           step_diff,         ///< how many steps are kept
    const bool           screen_out,
    const int            pid
    )
{
  std::ostringstream filename;
  std::ostringstream filenamedel;
  const std::string filebase(DRT::Problem::Instance()->OutputControlFile()->FileName());

  std::ostringstream pid_stream;
  pid_stream << ".p" << std::setw(2) << setfill('0') << pid;

  filename    << filebase << "." << filename_base << "_" << std::setw(5) << setfill('0') << actstep           << pid_stream.str() << ".pos";
  filenamedel << filebase << "." << filename_base << "_" << std::setw(5) << setfill('0') << actstep-step_diff << pid_stream.str() << ".pos";
  std::remove(filenamedel.str().c_str());
  if (screen_out) std::cout << "writing " << left << std::setw(60) <<filename.str()<<"...";
  return filename.str();
}

#endif // #ifdef CCADISCRET
