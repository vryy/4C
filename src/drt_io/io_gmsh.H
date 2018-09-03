/*!
 * \file io_gmsh.H
 * \brief simple element print library for Gmsh
 *
 * Useful for debugging, only output per processor (no communication supported)
 *
 * For a complete description, see the section about post processing formats in the Gmsh
documentation.
 *
 * Most routines come in pairs: one gives a std::string as return argument, another pipes
 * the same std::string directly to an output stream. The latter avoid numerous std::string copies
 * and should be much faster when large discretizations are piped into a file.
 *

<pre>

\level 2
\maintainer Martin Kronbichler
             kronbichler@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15235
</pre>

 */
#ifndef GMSH_H
#define GMSH_H

#include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

namespace IO
{
  //! routines to transform results (elements, nodes, integrationcells, whole discretizations...) to
  //! Gmsh output
  namespace GMSH
  {
    //! for each DRT::Element::DiscretizationType find the appropriate Gmsh file type
    inline std::string distypeToGmshElementHeader(
        const DRT::Element::DiscretizationType distype  ///< element shape
    )
    {
      switch (distype)
      {
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
        case DRT::Element::line2:
          return "L";
          break;
        case DRT::Element::line3:
          return "L2";
          break;
        case DRT::Element::wedge6:
          return "I";
          break;
        case DRT::Element::wedge15:
          return "I";
          break;
        default:
          dserror("distypeToGmshElementHeader: distype not supported for printout!");
      }
      return "xxx";
    }

    //! for each DRT::Element::DiscretizationType find the appropriate number of element nodes for
    //! Gmsh output
    inline int distypeToGmshNumNode(
        const DRT::Element::DiscretizationType distype  ///< element shape
    )
    {
      switch (distype)
      {
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
        case DRT::Element::line2:
          return 2;
          break;
        case DRT::Element::line3:
          return 3;
          break;
        case DRT::Element::wedge6:
          return 6;
          break;
        case DRT::Element::wedge15:
          return 6;
          break;
        default:
          dserror("distypeToGmshNumNode: distype not supported for printout!");
      }
      return -1;
    }

    //! write scalar field to Gmsh postprocessing file
    void ScalarFieldToGmsh(const Teuchos::RCP<DRT::Discretization> discret,  ///< discretization
        const Teuchos::RCP<const Epetra_Vector> scalarfield,  ///< scalar field to output
        std::ostream& s                                       ///< output stream
    );

    //! write scalar field to Gmsh postprocessing file
    //  Can we replace the ScalarFieldToGmsh function with this (nds=0)?
    void ScalarFieldDofBasedToGmsh(
        const Teuchos::RCP<DRT::Discretization> discret,      ///< discretization
        const Teuchos::RCP<const Epetra_Vector> scalarfield,  ///< scalar field to output
        const int nds,                                        ///< dofset
        std::ostream& s                                       ///< output stream
    );

    //! write scalar field to Gmsh postprocessing file
    void ScalarElementFieldToGmsh(
        const Teuchos::RCP<DRT::Discretization> discret,      ///< discretization
        const Teuchos::RCP<const Epetra_Vector> scalarfield,  ///< scalar field to output
        std::ostream& s                                       ///< output stream
    );

    //! write dof-based vector field to Gmsh postprocessing file
    //! when writing the given vectorfield is displacement, displacenodes can be set true
    //! then it writes the nodal coordinates at the present time step.
    void VectorFieldDofBasedToGmsh(
        const Teuchos::RCP<DRT::Discretization> discret,      ///< discretization
        const Teuchos::RCP<const Epetra_Vector> vectorfield,  ///< vector field to output
        std::ostream& s,                                      ///< output stream
        const int nds = 0,  ///< number of dofset associated with vector field
        bool displacenodes = false);

    //! write dof-based vector field to Gmsh postprocessing file
    void VectorFieldMultiVectorDofBasedToGmsh(
        const Teuchos::RCP<const DRT::Discretization> discret,     ///< discretization
        const Teuchos::RCP<const Epetra_MultiVector> vectorfield,  ///< vector field to output
        std::ostream& s,                                           ///< output stream
        const int nds = 0  //< which dof-set to use from vector
    );

    //! write node-based vector field to Gmsh postprocessing file
    void VectorFieldNodeBasedToGmsh(
        const Teuchos::RCP<const DRT::Discretization> discret,     ///< discretization
        const Teuchos::RCP<const Epetra_MultiVector> vectorfield,  ///< vector field to output
        std::ostream& s                                            ///< output stream
    );

    //! write dof-based vector field to Gmsh postprocessing file
    void SurfaceVectorFieldDofBasedToGmsh(
        const Teuchos::RCP<DRT::Discretization> discret,      ///< discretization
        const Teuchos::RCP<const Epetra_Vector> vectorfield,  ///< vector field to output
        std::map<int, LINALG::Matrix<3, 1>>& currpos,
        std::ostream& s,  ///< output stream
        const int nsd, const int numdofpernode);

    //! write dof-based velocity / pressure field to Gmsh postprocessing file
    //!
    void VelocityPressureFieldDofBasedToGmsh(
        const Teuchos::RCP<DRT::Discretization> discret,      ///< discretization
        const Teuchos::RCP<const Epetra_Vector> vectorfield,  ///< vector field to output
        const std::string field,                              ///< "velocity" or "pressure"
        std::ostream& s,                                      ///< output stream
        const int nds = 0                                     ///< which dof-set to use from vector
    );

    //! write node-based scalar field to Gmsh postprocessing file
    void ScalarFieldNodeBasedToGmsh(
        const Teuchos::RCP<const DRT::Discretization> discret,  ///< discretization
        const Teuchos::RCP<const Epetra_Vector> scalarfield,    ///< scalar field to output
        std::ostream& s                                         ///< output stream
    );

    //! take an array (3,numnode) and translate it to the coordinate section of a Gmsh element entry
    template <class M>
    inline void CoordinatesToStream(const M& coord,      ///< position array (3, numnode)
        const DRT::Element::DiscretizationType distype,  ///< element shape
        std::ostream& s                                  ///< output stream
    )
    {
      s.setf(std::ios::scientific, std::ios::floatfield);
      s.precision(12);

      const int numnode = distypeToGmshNumNode(distype);

      // coordinates
      s << "(";
      for (int inen = 0; inen < numnode; ++inen)
      {
        // print position
        s << coord(0, inen) << ",";
        s << coord(1, inen) << ",";
        s << coord(2, inen);
        if (inen < numnode - 1)
        {
          s << ",";
        }
      };
      s << ")";
      return;
    }

    //! take an array (3,numnode) and translate it to the coordinate section of a Gmsh element entry
    template <class M>
    inline void CoordinatesToStream2D(const M& coord,    ///< position array (3, numnode)
        const DRT::Element::DiscretizationType distype,  ///< element shape
        std::ostream& s)
    {
      s.setf(std::ios::scientific, std::ios::floatfield);
      s.precision(12);

      const int numnode = distypeToGmshNumNode(distype);

      // coordinates
      s << "(";
      for (int inen = 0; inen < numnode; ++inen)
      {
        // print position 2D
        s << coord(0, inen) << ",";
        s << coord(1, inen) << ",";
        s << "0.000000000000e+00";
        if (inen < numnode - 1)
        {
          s << ",";
        }
      };
      s << ")";
      return;
    }

    //! take a double value and translate it to the value section of an Gmsh element entry
    //! -> gives a constant value over the element; used to plot the mesh without any physical field
    inline void ScalarToStream(
        const double scalar,  ///< constant (arbitrary) value assigned to elements
        const DRT::Element::DiscretizationType distype,  ///< element shape
        std::ostream& s)
    {
      s.setf(std::ios::scientific, std::ios::floatfield);
      s.precision(12);

      const int numnode = distypeToGmshNumNode(distype);

      // values
      s << "{";
      for (int i = 0; i < numnode; ++i)
      {
        s << scalar;
        if (i < numnode - 1)
        {
          s << ",";
        }
      };
      s << "};";
    }

    //! take a scalar value at a point and translate it into Gmsh postprocessing format
    void ScalarToStream(const LINALG::Matrix<3, 1>& pointXYZ,  ///< coordinates of point
        const double scalarvalue,                              ///< scalar value at this point
        std::ostream& s                                        ///< stream
    );

    //! take a scalar value at a point and translate it into Gmsh postprocessing format
    void VectorToStream(const LINALG::Matrix<3, 1>& pointXYZ,  ///< coordinates of point
        const LINALG::Matrix<3, 1>& vectorvalue,               ///< vector at this point
        std::ostream& s                                        ///< stream
    );

    //! take an array (numnode) and translate it to the scalar value section of an Gmsh element
    //! entry
    template <class V>
    inline void ScalarFieldToStream(const V& scalarfield,
        const DRT::Element::DiscretizationType distype,  ///< element shape
        std::ostream& s)
    {
      s.setf(std::ios::scientific, std::ios::floatfield);
      s.precision(12);

      const int numnode = distypeToGmshNumNode(distype);

      // values
      s << "{";
      for (int inode = 0; inode < numnode; ++inode)
      {
        s << scalarfield(inode);
        if (inode < numnode - 1)
        {
          s << ",";
        }
      };
      s << "};";
    }

    //! take an array (3,numnode) and translate it to the vector value section of an Gmsh element
    //! entry
    template <class M>
    inline void VectorFieldToStream(const M& vectorfield,  ///< vector value array (3, numnode)
        const DRT::Element::DiscretizationType distype,    ///< element shape
        std::ostream& s)
    {
      s.setf(std::ios::scientific, std::ios::floatfield);
      s.precision(12);

      const int numnode = distypeToGmshNumNode(distype);
      const int nsd = 3;

      // values
      s << "{";
      for (int inode = 0; inode < numnode; ++inode)
      {
        for (int isd = 0; isd < nsd; ++isd)
        {
          s << std::scientific << vectorfield(isd, inode);
          if (isd < nsd - 1) s << ",";
        }
        if (inode < numnode - 1)
        {
          s << ",";
        }
      };
      s << "};";
    }


    //! take an array (3,numnode) and translate it to the vector value section of an Gmsh element
    //! entry
    template <class M>
    inline void VectorFieldToStream2D(const M& vectorfield,  ///< vector value array (3, numnode)
        const DRT::Element::DiscretizationType distype,      ///< element shape
        std::ostream& s)
    {
      s.setf(std::ios::scientific, std::ios::floatfield);
      s.precision(12);

      const int numnode = distypeToGmshNumNode(distype);

      // values
      s << "{";
      for (int inode = 0; inode < numnode; ++inode)
      {
        // 2D
        s << std::scientific << vectorfield(0, inode) << ",";
        s << std::scientific << vectorfield(1, inode) << ",";
        s << "0.000000000000e+00";
        if (inode < numnode - 1)
        {
          s << ",";
        }
      };
      s << "};";
    }

    //! take an array (9,numnode) and translate it to the tensor value section of an Gmsh element
    //! entry
    template <class M>
    inline void TensorFieldToStream(const M& tensorfield,
        const DRT::Element::DiscretizationType distype,  ///< element shape
        std::ostream& s)
    {
      s.setf(std::ios::scientific, std::ios::floatfield);
      s.precision(12);

      const int numnode = distypeToGmshNumNode(distype);
      const int numTensorEntries = 9;

      // values
      s << "{";
      for (int inode = 0; inode < numnode; ++inode)
      {
        for (int ientry = 0; ientry < numTensorEntries; ++ientry)
        {
          s << std::scientific << tensorfield(ientry, inode);
          if (ientry < numTensorEntries - 1) s << ",";
        }
        if (inode < numnode - 1)
        {
          s << ",";
        }
      };
      s << "};";
    }

    //! take an entire DRT::Element and print it with constant scalar value at its initial position
    void elementAtInitialPositionToStream(
        const double scalar, const DRT::Element* ele, std::ostream& s);

    //! take an entire DRT::Element and print it with constant scalar value at its initial position
    std::string elementAtInitialPositionToString(const double scalar, const DRT::Element* ele);

    //! take an entire DRT::Element and print it with constant scalar value at the given position
    void elementAtCurrentPositionToStream(
        const double scalar,      ///< scalar value for the entire element
        const DRT::Element* ele,  ///< element to print
        const std::map<int, LINALG::Matrix<3, 1>>& currentelepositions,  ///< nodal position array
        std::ostream& s);

    //! take an entire DRT::Element and print it with constant scalar value at the given position
    std::string elementAtCurrentPositionToString(
        const double scalar,      ///< scalar value for the entire element
        const DRT::Element* ele,  ///< element to print
        const std::map<int, LINALG::Matrix<3, 1>>& currentelepositions  ///< nodal position array
    );

    //! take an array (numnode) and translate it to the scalar value section of an Gmsh element
    //! entry
    template <class M>
    void cellWithScalarToStream(const DRT::Element::DiscretizationType distype,  ///< element shape
        const double scalar,  ///< scalar value for the entire element
        const M& xyze,        ///< position array (3, numnode)
        std::ostream& s)
    {
      s << "S";  // scalar field indicator
      s << distypeToGmshElementHeader(distype);
      CoordinatesToStream(xyze, distype, s);
      ScalarToStream(scalar, distype, s);
      s << "\n";
    }

    //! take an array (numnode) and translate it to the scalar value section of an Gmsh element
    //! entry
    template <class M>
    std::string cellWithScalarToString(
        const DRT::Element::DiscretizationType distype,  ///< element shape
        const double scalar,                             ///< scalar value for the entire element
        const M& xyze                                    ///< position array (3, numnode)
    )
    {
      std::ostringstream s;
      cellWithScalarToStream(distype, scalar, xyze, s);
      return s.str();
    }

    //! take an value array (numnode) and a position array (3,numnode) and translate it to an
    //! element section of a Gmsh postprocessing file
    template <class V, class M>
    inline void cellWithScalarFieldToStream(
        const DRT::Element::DiscretizationType distype,  ///< element shape
        const V& scalarfield,                            ///< scalar field in the element
        const M& xyze,                                   ///< position array (3, numnode)
        std::ostream& s)
    {
      s << "S";  // scalar field indicator
      s << distypeToGmshElementHeader(distype);
      CoordinatesToStream(xyze, distype, s);
      ScalarFieldToStream(scalarfield, distype, s);
      s << "\n";
    }

    //! take an value array (numnode) and a position array (3,numnode) and translate it to an
    //! element section of a Gmsh postprocessing file
    template <class V, class M>
    inline std::string cellWithScalarFieldToString(
        const DRT::Element::DiscretizationType distype,  ///< element shape
        const V& scalarfield,                            ///< scalar field in the element
        const M& xyze                                    ///< position array (3, numnode)
    )
    {
      std::ostringstream s;
      cellWithScalarFieldToStream(distype, scalarfield, xyze, s);
      return s.str();
    }

    //! take an value array (3,numnode) and a position array (3,numnode) and translate it to an
    //! element section of a Gmsh postprocessing file
    template <class M1, class M2>
    inline void cellWithVectorFieldToStream(
        const DRT::Element::DiscretizationType distype,  ///< element shape
        const M1& vectorfield,  ///< vector field in the element (3, numnode)
        const M2& xyze,         ///< position array (3, numnode)
        std::ostream& s)
    {
      s << "V";  // vector field indicator
      s << distypeToGmshElementHeader(distype);
      CoordinatesToStream(xyze, distype, s);
      VectorFieldToStream(vectorfield, distype, s);
      s << "\n";
    }

    //! take an value array (3,numnode) and a position array (3,numnode) and translate it to an
    //! element section of a Gmsh postprocessing file
    template <class M1, class M2>
    inline std::string cellWithVectorFieldToString(
        const DRT::Element::DiscretizationType distype,  ///< element shape
        const M1& vectorfield,  ///< vector field in the element (3, numnode)
        const M2& xyze          ///< position array (3, numnode)
    )
    {
      std::ostringstream s;
      cellWithVectorFieldToStream(distype, vectorfield, xyze, s);
      return s.str();
    }

    //! take an value array (9,numnode) and a position array (3,numnode) and translate it to an
    //! element section of a Gmsh postprocessing file
    template <class M1, class M2>
    inline void cellWithTensorFieldToStream(
        const DRT::Element::DiscretizationType distype,  ///< element shape
        const M1& tensorfield,  ///< tensor field in the element (9,numnode)
        const M2& xyze,         ///< position array (3, numnode)
        std::ostream& s)
    {
      s << "T";  // tensor field indicator
      s << distypeToGmshElementHeader(distype);
      CoordinatesToStream(xyze, distype, s);
      TensorFieldToStream(tensorfield, distype, s);
      s << "\n";
    }

    //! take an value array (9,numnode) and a position array (3,numnode) and translate it to an
    //! element section of a Gmsh postprocessing file
    template <class M1, class M2>
    inline std::string cellWithTensorFieldToString(
        const DRT::Element::DiscretizationType distype,  ///< element shape
        const M1& tensorfield,  ///< tensor field in the element (9,numnode)
        const M2& xyze          ///< position array (3, numnode)
    )
    {
      std::ostringstream s;
      cellWithTensorFieldToStream(distype, tensorfield, xyze, s);
      return s.str();
    }

    /// print a piece of text at a given position
    std::string text3dToString(const LINALG::Matrix<3, 1>& xyz,  ///< 3d Position of text
        const std::string& text,                                 ///< text to be printed
        const int fontsize                                       ///< font size
    );

    //! print discretization in initial configuration (t = 0)
    void disToStream(const std::string& text, const double scalar,
        const Teuchos::RCP<DRT::Discretization> dis, std::ostream& s);

    //! print discretization in initial configuration (t = 0)
    std::string disToString(
        const std::string& text, const double scalar, const Teuchos::RCP<DRT::Discretization> dis);

    //! print discretization in current configuration (t > 0)
    void disToStream(const std::string& text, const double scalar,
        const Teuchos::RCP<DRT::Discretization> dis,
        const std::map<int, LINALG::Matrix<3, 1>>& currentpositions, std::ostream& s);

    //! print discretization in current configuration (t > 0)
    std::string disToString(const std::string& text, const double scalar,
        const Teuchos::RCP<DRT::Discretization> dis,
        const std::map<int, LINALG::Matrix<3, 1>>& currentpositions);

    std::string GetNewFileNameAndDeleteOldFiles(const std::string& filename_base,
        const int& actstep,    ///< generate filename for this step
        const int& step_diff,  ///< how many steps are kept
        const bool screen_out,
        const int pid = 0  ///< my processor id
    );

    //! open Gmsh output file to add data
    std::string GetFileName(const std::string& filename_base,
        const int& actstep,  ///< generate filename for this step
        const bool screen_out,
        const int pid = 0  ///< my processor id
    );
  }  // namespace GMSH
}  // namespace IO

#endif /*GMSH_H*/
