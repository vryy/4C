/*!
 * \file io_gmsh.H
 * \brief simple element print library for Gmsh (mostly for debugging, only output per processor (no communication)
 *
 * <pre>
 * Maintainer: Axel Gerstenberger
 *             gerstenberger@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15236
 * </pre>
 *
 * For a complete description, see the section about post processing formations in the Gmsh documentation
 */
#ifndef GMSH_H
#define GMSH_H
#ifdef CCADISCRET

#include "../drt_lib/drt_discret.H"

namespace IO
{
  //! routines to transform results (elements, nodes, integrationcells, whole discretizations...) to Gmsh output
  namespace GMSH
  {
    //! for each DRT::Element::DiscretizationType find the appropriate Gmsh file type
    std::string distypeToGmshElementHeader(
        const DRT::Element::DiscretizationType distype   ///< element shape
        );

    //! for each DRT::Element::DiscretizationType find the apropriate number of element nodes for Gmsh output
    int distypeToGmshNumNode(
        const DRT::Element::DiscretizationType distype   ///< element shape
      );

    //! take an array (3,numnode) and translate it to the coordinate section of en Gmsh element entry
    template<class M> std::string CoordinatesToString(
        const M& coord,
        const DRT::Element::DiscretizationType distype   ///< element shape
    )
    {
      std::stringstream pos_array_string;
      pos_array_string.setf(ios::scientific,ios::floatfield);
      pos_array_string.precision(12);

      const int numnode = distypeToGmshNumNode(distype);

      // coordinates
      pos_array_string << "(";
      for (int inen = 0; inen < numnode; ++inen)
      {
        // print position
        pos_array_string << coord(0, inen) << ",";
        pos_array_string << coord(1, inen) << ",";
        pos_array_string << coord(2, inen);
        if (inen < numnode-1)
        {
          pos_array_string << ",";
        }
      };
      pos_array_string << ")";
      return pos_array_string.str();
    }

    //! take a double value and translate it to the value section of an Gmsh element entry -> gives a constant value over the element
    std::string ScalarToString(
        const double scalar,
        const DRT::Element::DiscretizationType distype      ///< element shape
    );

    //! take an array (numnode) and translate it to the scalar value section of an Gmsh element entry
    template<class V> std::string ScalarFieldToString(
        const V& scalarfield,
        const DRT::Element::DiscretizationType distype      ///< element shape
    )
    {
      std::stringstream s;
      s.setf(ios::scientific,ios::floatfield);
      s.precision(12);

      const int numnode = distypeToGmshNumNode(distype);

      // values
      s << "{";
      for (int inode = 0; inode<numnode; ++inode)
      {
        s << scalarfield(inode);
        if (inode < numnode-1)
        {
          s << ",";
        }
      };
      s << "};";
      return s.str();
    }

    //! take an array (3,numnode) and translate it to the vector value section of an Gmsh element entry
    template<class M> std::string VectorFieldToString(
        const M& vectorfield,
        const DRT::Element::DiscretizationType distype      ///< element shape
    )
    {
      std::stringstream s;
      s.setf(ios::scientific,ios::floatfield);
      s.precision(12);

      const int numnode = distypeToGmshNumNode(distype);
      const int nsd = 3;

      // values
      s << "{";
      for (int inode = 0; inode<numnode; ++inode)
      {
        for (int isd = 0; isd < nsd; ++isd)
        {
          s << scientific << vectorfield(isd, inode);
          if (isd < nsd-1)
            s << ",";
        }
        if (inode < numnode-1)
        {
          s << ",";
        }
      };
      s << "};";
      return s.str();
    }

    //! take an array (9,numnode) and translate it to the tensor value section of an Gmsh element entry
    template<class M> std::string TensorFieldToString(
        const M& tensorfield,
        const DRT::Element::DiscretizationType distype      ///< element shape
    )
    {
      std::stringstream s;
      s.setf(ios::scientific,ios::floatfield);
      s.precision(12);

      const int numnode = distypeToGmshNumNode(distype);
      const int numTensorEntries = 9;

      // values
      s << "{";
      for (int inode = 0; inode<numnode; ++inode)
      {
        for (int ientry = 0; ientry < numTensorEntries; ++ientry)
        {
          s << scientific << tensorfield(ientry, inode);
          if (ientry < numTensorEntries-1)
            s << ",";
        }
        if (inode < numnode-1)
        {
          s << ",";
        }
      };
      s << "};";
      return s.str();
    }

    //! take an entire DRT::Element and print it with constant scalar value at its initial pusition
    std::string elementAtInitialPositionToString(
        const double scalar,
        const DRT::Element* ele);

    //! take an entire DRT::Element and print it with constant scalar value at the given pusition
    std::string elementAtCurrentPositionToString(
        const double                            scalar,                 ///< scalar value for the entire element
        const DRT::Element*                     ele,                    ///< element to print
        const map<int,LINALG::Matrix<3,1> >&    currentelepositions     ///< nodal position array
    );

    //! take an array (numnode) and translate it to the scalar value section of an Gmsh element entry
    template<class M> std::string cellWithScalarToString(
        const DRT::Element::DiscretizationType distype,   ///< element shape
        const double scalar,                              ///< scalar value for the entire element
        const M& xyze                                     ///< position array (3, numnode)
    )
    {
      std::stringstream gmsh_ele_line;
      // header for this element
      gmsh_ele_line << "S" << distypeToGmshElementHeader(distype);
      // coordinates
      gmsh_ele_line << CoordinatesToString(xyze, distype);
      // values
      gmsh_ele_line << ScalarToString(scalar, distype);

      return gmsh_ele_line.str();
    }

    //! take an value array (numnode) and a position array (3,numnode) and translate it to an element section of a Gmsh postprocessing file
    template<class V, class M> std::string cellWithScalarFieldToString(
        const DRT::Element::DiscretizationType  distype,        ///< element shape
        const V&                                scalarfield,    ///< scalar field in the element
        const M&                                xyze            ///< position array (3, numnode)
    )
    {
      std::stringstream gmsh_ele_line;
      // header for this element
      gmsh_ele_line << "S" << distypeToGmshElementHeader(distype);
      // coordinates
      gmsh_ele_line << CoordinatesToString(xyze, distype);
      // values
      gmsh_ele_line << ScalarFieldToString(scalarfield, distype);

      return gmsh_ele_line.str();
    }

    //! take an value array (3,numnode) and a position array (3,numnode) and translate it to an element section of a Gmsh postprocessing file
    template<class M1, class M2> std::string cellWithVectorFieldToString(
        const DRT::Element::DiscretizationType  distype,        ///< element shape
        const M1&                               vectorfield,    ///< vector field in the element
        const M2&                               xyze            ///< position array
        )
    {
      std::stringstream gmsh_ele_line;
      // header for this element
      gmsh_ele_line << "V" << distypeToGmshElementHeader(distype);
      // coordinates
      gmsh_ele_line << CoordinatesToString(xyze, distype);
      // values
      gmsh_ele_line << VectorFieldToString(vectorfield, distype);

      return gmsh_ele_line.str();
    }

    //! take an value array (9,numnode) and a position array (3,numnode) and translate it to an element section of a Gmsh postprocessing file
    template<class M1, class M2> std::string cellWithTensorFieldToString(
        const DRT::Element::DiscretizationType  distype,        ///< element shape
        const M1&                               tensorfield,    ///< tensor field in the element
        const M2&                               xyze            ///< position array
        )
    {
      std::stringstream gmsh_ele_line;
      // header for this element
      gmsh_ele_line << "T" << distypeToGmshElementHeader(distype);
      // coordinates
      gmsh_ele_line << CoordinatesToString(xyze, distype);
      // values
      gmsh_ele_line << TensorFieldToString(tensorfield, distype);

      return gmsh_ele_line.str();
    }

    /// print a piece of text at a given position
    std::string text3dToString(
        const LINALG::Matrix<3,1>&            xyz,      ///< 3d Position of text
        const std::string                     text,     ///< text to be printed
        const int                             fontsize  ///< font size
        );

    //! print discretization in initial configuration (t = 0)
    std::string disToString(
        const std::string& s,
        const double scalar,
        const Teuchos::RCP<DRT::Discretization> dis);

    //! print discretization in current configuration (t > 0)
    std::string disToString(
        const std::string&                         s,
        const double                               scalar,
        const Teuchos::RCP<DRT::Discretization>    dis,
        const std::map<int,LINALG::Matrix<3,1> >&  currentpositions);

  }
}

#endif // #ifdef CCADISCRET
#endif /*GMSH_H*/
