/*!
 \brief simple element print library for Gmsh (mostly for debuging, only sequential)

 <pre>
 Maintainer: Axel Gerstenberger
             gerstenberger@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236
 </pre>
 */
#ifndef GMSH_H
#define GMSH_H
#ifdef CCADISCRET

#include <blitz/array.h>
#include <string>
#include <vector>

#include "../drt_xfem/integrationcell.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_discret.H"

namespace IO
{
  //! routines to transform results (elements, nodes, integrationcells, whole discretizations...) to Gmsh output
  namespace GMSH
  {

    inline std::string distypeToGmshElementHeader(
        const DRT::Element::DiscretizationType distype)
    {
      switch (distype)
      {
        case DRT::Element::point1:
          return "P";
          break;
        case DRT::Element::quad4:
          return "Q";
          break;
        case DRT::Element::quad8:
          return "Q";
          break;
        case DRT::Element::quad9:
          return "Q";
          break;
        case DRT::Element::tri3:
          return "T";
          break;
        case DRT::Element::tri6:
          return "T";
          break;
        case DRT::Element::hex8:
          return "H";
          break;
        case DRT::Element::hex20:
          return "H";
          break;
        case DRT::Element::hex27:
          return "H";
          break;
        case DRT::Element::tet4:
          return "S";
          break;
        case DRT::Element::tet10:
          return "S";
          break;
        case DRT::Element::line2:
          return "L";
          break;
        case DRT::Element::line3:
          return "L2";
          break;
        default:
          dserror("distypeToGmshElementHeader: distype not supported for printout!");
      }
      return "xxx";
    }

    inline int distypeToGmshNumNode(
        const DRT::Element::DiscretizationType distype)
    {
      switch (distype)
      {
        case DRT::Element::point1:
          return 1;
          break;
        case DRT::Element::quad4:
          return 4;
          break;
        case DRT::Element::quad8:
          return 4;
          break;
        case DRT::Element::quad9:
          return 4;
          break;
        case DRT::Element::tri3:
          return 3;
          break;
        case DRT::Element::tri6:
          return 3;
          break;
        case DRT::Element::hex8:
          return 8;
          break;
        case DRT::Element::hex20:
          return 8;
          break;
        case DRT::Element::hex27:
          return 8;
          break;
        case DRT::Element::tet4:
          return 4;
          break;
        case DRT::Element::tet10:
          return 4;
          break;
        case DRT::Element::line2:
          return 2;
          break;
        case DRT::Element::line3:
          return 3;
          break;
        default:
          dserror("distypeToGmshNumNode: distype not supported for printout!");
      }
      return -1;
    }

    /*!
     \brief  fill array with current nodal positions

     \return array with element nodal positions
     */
    inline blitz::Array<double,2> getCurrentNodalPositions(
        const DRT::Element* ele, ///< element with nodal pointers
        const std::map<int,blitz::TinyVector<double,3> >& currentcutterpositions ///< current positions of all cutter nodes
    )
    {
      const int numnode = ele->NumNode();
      blitz::Array<double,2> xyze(3, numnode);
      const DRT::Node*const* nodes = ele->Nodes();
      for (int inode = 0; inode < numnode; ++inode)
      {
        const blitz::TinyVector<double,3> x = currentcutterpositions.find(nodes[inode]->Id())->second;
        xyze(0, inode) = x(0);
        xyze(1, inode) = x(1);
        xyze(2, inode) = x(2);
      }
      return xyze;
    }

    template<class M> std::string CoordinatesToString(
        const M& coord,
        const DRT::Element::DiscretizationType distype)
    {
      stringstream pos_array_string;

      const int numnode = distypeToGmshNumNode(distype);

      // coordinates
      pos_array_string << "(";
      for (int inen = 0; inen < numnode; ++inen)
      {
        // print position
        pos_array_string << scientific << coord(0, inen) << ",";
        pos_array_string << scientific << coord(1, inen) << ",";
        pos_array_string << scientific << coord(2, inen);
        if (inen < numnode-1)
        {
          pos_array_string << ",";
        }
      };
      pos_array_string << ")";
      return pos_array_string.str();
    }

    std::string ScalarToString(
        const double scalar,
        const DRT::Element::DiscretizationType distype);

    template<class V> std::string ScalarFieldToString(
        const V& scalarfield,
        const DRT::Element::DiscretizationType distype)
    {
      stringstream pos_array_string;

      const int numnode = distypeToGmshNumNode(distype);

      // values
      pos_array_string << "{";
      for (int inode = 0; inode<numnode; ++inode)
      {
        pos_array_string << scientific << scalarfield(inode);
        if (inode < numnode-1)
        {
          pos_array_string << ",";
        }
      };
      pos_array_string << "};";
      return pos_array_string.str();
    }

    template<class M> std::string VectorFieldToString(
        const M& vectorfield,
        const DRT::Element::DiscretizationType distype)
    {
      stringstream pos_array_string;

      const int numnode = distypeToGmshNumNode(distype);
      const int nsd = 3;

      // values
      pos_array_string << "{";
      for (int inode = 0; inode<numnode; ++inode)
      {
        for (int isd = 0; isd < nsd; ++isd)
        {
          pos_array_string << scientific << vectorfield(isd, inode);
          if (isd < nsd-1)
            pos_array_string << ",";
        }
        if (inode < numnode-1)
        {
          pos_array_string << ",";
        }
      };
      pos_array_string << "};";
      return pos_array_string.str();
    }
    
    template<class M> std::string TensorFieldToString(
        const M& tensorfield,
        const DRT::Element::DiscretizationType distype)
    {
      stringstream pos_array_string;

      const int numnode = distypeToGmshNumNode(distype);
      const int numTensorEntries = 9;

      // values
      pos_array_string << "{";
      for (int inode = 0; inode<numnode; ++inode)
      {
        for (int ientry = 0; ientry < numTensorEntries; ++ientry)
        {
          pos_array_string << scientific << tensorfield(ientry, inode);
          if (ientry < numTensorEntries-1)
            pos_array_string << ",";
        }
        if (inode < numnode-1)
        {
          pos_array_string << ",";
        }
      };
      pos_array_string << "};";
      return pos_array_string.str();
    }

    std::string elementToString(
        const double scalar,
        const DRT::Element* ele);

    template<class M> std::string cellWithScalarToString(
        const DRT::Element::DiscretizationType distype,
        const double scalar,
        const M& xyze)
    {
      stringstream gmsh_ele_line;
      // header for this element
      gmsh_ele_line << "S" << distypeToGmshElementHeader(distype);
      // coordinates
      gmsh_ele_line << CoordinatesToString(xyze, distype);
      // values
      gmsh_ele_line << ScalarToString(scalar, distype);

      return gmsh_ele_line.str();
    }

    template<class V, class M> std::string cellWithScalarFieldToString(
        const DRT::Element::DiscretizationType distype,
        const V& scalarfield,
        const M& xyze)
    {
      stringstream gmsh_ele_line;
      // header for this element
      gmsh_ele_line << "S" << distypeToGmshElementHeader(distype);
      // coordinates
      gmsh_ele_line << CoordinatesToString(xyze, distype);
      // values
      gmsh_ele_line << ScalarFieldToString(scalarfield, distype);

      return gmsh_ele_line.str();
    }

    template<class M1, class M2> std::string cellWithVectorFieldToString(
        const DRT::Element::DiscretizationType distype,
        const M1& vectorfield,
        const M2& xyze)
    {
      stringstream gmsh_ele_line;
      // header for this element
      gmsh_ele_line << "V" << distypeToGmshElementHeader(distype);
      // coordinates
      gmsh_ele_line << CoordinatesToString(xyze, distype);
      // values
      gmsh_ele_line << VectorFieldToString(vectorfield, distype);

      return gmsh_ele_line.str();
    }
    
    template<class M1, class M2> std::string cellWithTensorFieldToString(
        const DRT::Element::DiscretizationType distype,
        const M1& tensorfield,
        const M2& xyze)
    {
      stringstream gmsh_ele_line;
      // header for this element
      gmsh_ele_line << "T" << distypeToGmshElementHeader(distype);
      // coordinates
      gmsh_ele_line << CoordinatesToString(xyze, distype);
      // values
      gmsh_ele_line << TensorFieldToString(tensorfield, distype);

      return gmsh_ele_line.str();
    }

    //! print discretization in initial configuration (t = 0)
    std::string disToString(
        const std::string& s,
        const double scalar,
        const Teuchos::RCP<DRT::Discretization> dis);

    //! print discretization in current configuration (t > 0)
    std::string disToString(
        const std::string& s,
        const double scalar,
        const Teuchos::RCP<DRT::Discretization> dis,
        std::map<int,blitz::TinyVector<double,3> > currentpositions);

    std::string disToString(const std::string& s,
        const double scalar,
        const Teuchos::RCP<DRT::Discretization> dis,
        const std::map<int,XFEM::DomainIntCells >& elementDomainIntCellsMap,
        const std::map<int,XFEM::BoundaryIntCells >& elementBoundaryIntCellsMap);

    std::string getConfigString(const int numview);

  }
}

#endif // #ifdef CCADISCRET
#endif /*GMSH_H*/
