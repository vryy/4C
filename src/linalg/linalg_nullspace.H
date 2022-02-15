/*----------------------------------------------------------------------*/
/*! \file

\brief Calcultion of nullspace

\level 0

*----------------------------------------------------------------------*/
#ifndef LINALG_NULLSPACE_H
#define LINALG_NULLSPACE_H

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_elementtype.H"

namespace LINALG
{
  class Nullspace
  {
   public:
    static void ComputeNullSpace(DRT::Discretization& dis,
        Teuchos::RCP<std::vector<double>> nullspace, int numdf = 0, int dimns = 0)
    {
      // get number of dofs to check or allocate nullspace storage
      const int lrows = dis.DofRowMap(0)->NumMyElements();

      // if nullspace is null then do everything from scratch
      if (nullspace->size() == 0)
      {
        int nv, np = 0;  // dummy - not needed here

        // downwinding needs nodal block information, compute it
        if (dis.NumMyRowElements())
        {
          // We assume that all elements are of equal type
          DRT::Element* dwele = dis.lRowElement(0);
          dwele->ElementType().NodalBlockInformation(dwele, numdf, dimns, nv, np);
        }

        // communicate data to procs without row element
        int ldata[2] = {numdf, dimns};
        int gdata[2] = {0, 0};
        dis.Comm().MaxAll(&ldata[0], &gdata[0], 2);
        numdf = gdata[0];
        dimns = gdata[1];

        // allocate storage for vector holding the nullspace data
        nullspace->resize(dimns * lrows, 0.0);
      }
      // check if nullspace has correct size
      else if (((ssize_t)(nullspace->size()) != (dimns * lrows) or numdf == 0))
      {
        dserror("nullspace does not have correct size or numdf or dimns are zero");
      }

      // check whether numdf and dimns are consistent among all procs
      const int numproc = dis.Comm().NumProc();
      int sumnumdf;
      dis.Comm().SumAll(&numdf, &sumnumdf, 1);
      if (sumnumdf != numdf * numproc) dserror("numdf not consistent among procs");
      int sumdimns;
      dis.Comm().SumAll(&dimns, &sumdimns, 1);
      if (sumdimns != dimns * numproc) dserror("dimns not consistent among procs");

      // check if dimns is possible
      if (dimns > 10) dserror("Nullspace size only upto 10 supported");

      // compute nullspace for simple case: vector of ones
      if (dimns == 1 && numdf == 1)
      {
        for (int i = 0; i < lrows; ++i) (*nullspace)[i] = 1.0;
        return;
      }

      // for rigid body rotations
      // compute nodal center of the discretization
      double x0send[3] = {0.0, 0.0, 0.0};
      for (int i = 0; i < dis.NumMyRowNodes(); ++i)
        for (int j = 0; j < 3; ++j) x0send[j] += dis.lRowNode(i)->X()[j];
      double x0[3];
      dis.Comm().SumAll(x0send, x0, 3);
      for (int i = 0; i < 3; ++i) x0[i] /= dis.NumGlobalNodes();

      // let the elementtype compute the nullspace
      DRT::Element* dwele = dis.lRowElement(0);
      dwele->ElementType().ComputeNullSpace(dis, *nullspace, x0, numdf, dimns);
    }
  };
}  // namespace LINALG

#endif
