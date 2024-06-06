/*-----------------------------------------------------------------------*/
/*! \file
\brief Contact interface capable of TSI


\level 3
*/
/*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Header                                                    seitz 08/15|
 *----------------------------------------------------------------------*/
#include "4C_contact_tsi_interface.hpp"

#include "4C_contact_defines.hpp"
#include "4C_contact_element.hpp"
#include "4C_contact_friction_node.hpp"
#include "4C_contact_interface.hpp"
#include "4C_contact_node.hpp"
#include "4C_inpar_mortar.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_mortar_dofset.hpp"
#include "4C_mortar_element.hpp"
#include "4C_mortar_node.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                            seitz 08/15|
 *----------------------------------------------------------------------*/
CONTACT::TSIInterface::TSIInterface(
    const Teuchos::RCP<Mortar::InterfaceDataContainer>& interfaceData_ptr, const int id,
    const Epetra_Comm& comm, const int dim, const Teuchos::ParameterList& icontact,
    bool selfcontact)
    : CONTACT::Interface(interfaceData_ptr, id, comm, dim, icontact, selfcontact)
{
  return;
}
void CONTACT::TSIInterface::AssembleLinStick(Core::LinAlg::SparseMatrix& linstickLMglobal,
    Core::LinAlg::SparseMatrix& linstickDISglobal, Core::LinAlg::SparseMatrix& linstickTEMPglobal,
    Epetra_Vector& linstickRHSglobal)
{
  CONTACT::Interface::AssembleLinStick(linstickLMglobal, linstickDISglobal, linstickRHSglobal);

#ifdef CONSISTENTSTICK

  // do the additional thermal linearizations
  // be aware, that there is another contribution to the linDIS part, since the temperature used
  // in the calculation of the friction coefficient may be the master side temperature, which has
  // a displacement derivative due to the projection


  // create map of stick nodes
  Teuchos::RCP<Epetra_Map> sticknodes = Core::LinAlg::SplitMap(*activenodes_, *slipnodes_);
  Teuchos::RCP<Epetra_Map> stickt = Core::LinAlg::SplitMap(*activet_, *slipt_);

  // nothing to do if no stick nodes
  if (sticknodes->NumMyElements() == 0) return;

  // information from interface contact parameter list
  Inpar::CONTACT::FrictionType ftype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(interface_params(), "FRICTION");
  if (ftype != Inpar::CONTACT::friction_coulomb) FOUR_C_THROW("only coulomb friction for CTSI");

  double frcoeff_in =
      interface_params().get<double>("FRCOEFF");  // the friction coefficient from the input
  double cn = interface_params().get<double>("SEMI_SMOOTH_CN");

  // some things that are not implemented
  bool gp_slip = Core::UTILS::IntegralValue<int>(interface_params(), "GP_SLIP_INCR");
  bool frilessfirst = Core::UTILS::IntegralValue<int>(interface_params(), "FRLESS_FIRST");
  if (gp_slip || frilessfirst)
    FOUR_C_THROW("this fancy option for the contact algorithm is not implemented for TSI");

  // consistent equation is:
  // mu(d,T)*(z_n-c_n*g)ut = 0

  typedef std::map<int, double>::const_iterator _CI;
  // loop over all stick nodes of the interface
  for (int i = 0; i < sticknodes->NumMyElements(); ++i)
  {
    int gid = sticknodes->GID(i);
    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      FOUR_C_THROW("AssembleLinStick: Node ownership inconsistency!");

    const Core::LinAlg::Matrix<3, 1> n(cnode->MoData().n(), true);
    const Core::LinAlg::Matrix<3, 1> lm(cnode->MoData().lm(), true);
    const double lm_n = n.Dot(lm);
    const double wgap = cnode->Data().Getg();
    Core::LinAlg::Matrix<3, 1> txi(cnode->Data().txi(), true);
    Core::LinAlg::Matrix<3, 1> teta(cnode->Data().teta(), true);
    Core::LinAlg::Matrix<3, 1> jump(cnode->FriData().jump(), true);
    double jump_txi = jump.Dot(txi);
    double jump_teta = jump.Dot(teta);

    // row number of entries
    std::vector<int> row(Dim() - 1);
    if (Dim() == 2)
      FOUR_C_THROW(
          "2-Dimensional thermo-contact not implemented 'cause there's no 2-D thermo element");
    else if (Dim() == 3)
    {
      row[0] = stickt->GID(2 * i);
      row[1] = stickt->GID(2 * i) + 1;
    }
    else
      FOUR_C_THROW("AssemblelinStick: Dimension not correct");

    std::map<int, double> dfrdT, dfrdD;
    cnode->derivFrCoeffTemp(frcoeff_in, dfrdT, dfrdD);

    double fac = lm_n - cn * wgap;
    for (_CI p = dfrdT.begin(); p != dfrdT.end(); ++p)
    {
      if (constr_direction_ == Inpar::CONTACT::constr_xyz)
        for (int j = 0; j < Dim(); j++)
        {
          linstickTEMPglobal.Assemble(
              fac * jump_txi * p->second * txi(j), cnode->Dofs()[j], p->first);
          linstickTEMPglobal.Assemble(
              fac * jump_teta * p->second * teta(j), cnode->Dofs()[j], p->first);
        }
      else
      {
        linstickTEMPglobal.Assemble(fac * jump_txi * p->second, row[0], p->first);
        linstickTEMPglobal.Assemble(fac * jump_teta * p->second, row[1], p->first);
      }
    }
    for (_CI p = dfrdD.begin(); p != dfrdD.end(); ++p)
    {
      if (constr_direction_ == Inpar::CONTACT::constr_xyz)
        for (int j = 0; j < Dim(); j++)
        {
          linstickDISglobal.Assemble(
              fac * jump_txi * p->second * txi(j), cnode->Dofs()[j], p->first);
          linstickDISglobal.Assemble(
              fac * jump_teta * p->second * teta(j), cnode->Dofs()[j], p->first);
        }
      else
      {
        linstickDISglobal.Assemble(fac * jump_txi * p->second, row[0], p->first);
        linstickDISglobal.Assemble(fac * jump_teta * p->second, row[1], p->first);
      }
    }
  }
#else
  // if we don't do the inconsistent stick linearization, it reduces to a mesh-tying
  // problem and hence it is independent of the friction coefficient and therefore
  // independent of the temperature. So do nothing here.
#endif

  return;
}
void CONTACT::TSIInterface::AssembleLinSlip(Core::LinAlg::SparseMatrix& linslipLMglobal,
    Core::LinAlg::SparseMatrix& linslipDISglobal, Core::LinAlg::SparseMatrix& linslipTEMPglobal,
    Epetra_Vector& linslipRHSglobal)
{
  CONTACT::Interface::AssembleLinSlip(linslipLMglobal, linslipDISglobal, linslipRHSglobal);

  // do the additional thermal linearizations
  // be aware, that there is another contribution to the linDIS part, since the temperature used
  // in the calculation of the friction coefficient may be the master side temperature, which has
  // a displacement derivative due to the projection


#ifdef CONSISTENTSLIP
  FOUR_C_THROW("CONSISTENTSLIP not implemented for CTSI");
#endif

  // create map of stick nodes
  Teuchos::RCP<Epetra_Map> slipnodes = slipnodes_;
  Teuchos::RCP<Epetra_Map> slipt = slipt_;

  // nothing to do if no stick nodes
  if (slipnodes->NumMyElements() == 0) return;

  // information from interface contact parameter list
  Inpar::CONTACT::FrictionType ftype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(interface_params(), "FRICTION");
  if (ftype != Inpar::CONTACT::friction_coulomb) FOUR_C_THROW("only coulomb friction for CTSI");

  if (Dim() != 3) FOUR_C_THROW("CTSI only for 3D");

  double frcoeff_in =
      interface_params().get<double>("FRCOEFF");  // the friction coefficient from the input
  double ct_input = interface_params().get<double>("SEMI_SMOOTH_CT");
  double cn_input = interface_params().get<double>("SEMI_SMOOTH_CN");

  // some things that are not implemented
  bool gp_slip = Core::UTILS::IntegralValue<int>(interface_params(), "GP_SLIP_INCR");
  bool frilessfirst = Core::UTILS::IntegralValue<int>(interface_params(), "FRLESS_FIRST");
  if (gp_slip || frilessfirst)
    FOUR_C_THROW("this fancy option for the contact algorithm is not implemented for TSI");

  // consistent equation is:
  // || z_t^tr || z_t - \mu (z_n -c_n*wgap) z_t^tr = 0

  typedef std::map<int, double>::const_iterator _CI;
  // loop over all stick nodes of the interface
  for (int i = 0; i < slipnodes->NumMyElements(); ++i)
  {
    int gid = slipnodes->GID(i);
    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    FriNode* cnode = dynamic_cast<FriNode*>(node);

    if (cnode->Owner() != Comm().MyPID())
      FOUR_C_THROW("AssembleLinStick: Node ownership inconsistency!");

    const Core::LinAlg::Matrix<3, 1> n(cnode->MoData().n(), true);
    const Core::LinAlg::Matrix<3, 1> lm(cnode->MoData().lm(), true);
    const double lm_n = n.Dot(lm);
    const double wgap = cnode->Data().Getg();
    Core::LinAlg::Matrix<3, 1> txi(cnode->Data().txi(), true);
    Core::LinAlg::Matrix<3, 1> teta(cnode->Data().teta(), true);
    Core::LinAlg::Matrix<3, 1> jump(cnode->FriData().jump(), true);
    const double jump_txi = jump.Dot(txi);
    const double jump_teta = jump.Dot(teta);
    const double lm_txi = lm.Dot(txi);
    const double lm_teta = lm.Dot(teta);

    // row number of entries
    std::vector<int> row(Dim() - 1);
    if (Dim() == 2)
      FOUR_C_THROW(
          "2-Dimensional thermo-contact not implemented 'cause there's no 2-D thermo element");
    else if (Dim() == 3)
    {
      row[0] = slipt->GID(2 * i);
      if (slipt->GID(2 * i) != cnode->Dofs()[1]) FOUR_C_THROW("why");
      row[1] = slipt->GID(2 * i) + 1;
      if (slipt->GID(2 * i) + 1 != cnode->Dofs()[2]) FOUR_C_THROW("why2");
    }
    else
      FOUR_C_THROW("AssemblelinStick: Dimension not correct");

    std::map<int, double> dfrdT, dfrdD;
    cnode->derivFrCoeffTemp(frcoeff_in, dfrdT, dfrdD);

    double cn = cn_input;
    double ct = ct_input;

    double fac = lm_n - cn * wgap;
    for (_CI p = dfrdT.begin(); p != dfrdT.end(); ++p)
    {
      if (constr_direction_ == Inpar::CONTACT::constr_xyz)
        for (int j = 0; j < Dim(); ++j)
        {
          linslipTEMPglobal.Assemble(
              -fac * (lm_txi + ct * jump_txi) * p->second * txi(j), cnode->Dofs()[j], p->first);
          linslipTEMPglobal.Assemble(
              -fac * (lm_teta + ct * jump_teta) * p->second * teta(j), cnode->Dofs()[j], p->first);
        }
      else
      {
        linslipTEMPglobal.Assemble(-fac * (lm_txi + ct * jump_txi) * p->second, row[0], p->first);
        linslipTEMPglobal.Assemble(-fac * (lm_teta + ct * jump_teta) * p->second, row[1], p->first);
      }
    }
    for (_CI p = dfrdD.begin(); p != dfrdD.end(); ++p)
    {
      if (constr_direction_ == Inpar::CONTACT::constr_xyz)
        for (int j = 0; j < Dim(); ++j)
        {
          linslipDISglobal.Assemble(
              -fac * (lm_txi + ct * jump_txi) * p->second * txi(j), cnode->Dofs()[j], p->first);
          linslipDISglobal.Assemble(
              -fac * (lm_teta + ct * jump_teta) * p->second * teta(j), cnode->Dofs()[j], p->first);
        }
      else
      {
        linslipDISglobal.Assemble(-fac * (lm_txi + ct * jump_txi) * p->second, row[0], p->first);
        linslipDISglobal.Assemble(-fac * (lm_teta + ct * jump_teta) * p->second, row[1], p->first);
      }
    }
  }

  return;
}


void CONTACT::TSIInterface::AssembleLinConduct(Core::LinAlg::SparseMatrix& linConductDISglobal,
    Core::LinAlg::SparseMatrix& linConductTEMPglobal,
    Core::LinAlg::SparseMatrix& linConductThermoLMglobal,
    Core::LinAlg::SparseMatrix& linConductContactLMglobal)
{
  // nothing to do if no active contact nodes
  if (activenodes_->NumMyElements() == 0) return;

  // heat transfer factors
  const double alpha_s = interface_params().get<double>("HEATTRANSSLAVE");
  const double alpha_m = interface_params().get<double>("HEATTRANSMASTER");
  const double beta_bar = alpha_s * alpha_m / (alpha_s + alpha_m);
  const double delta_bar = alpha_s / (alpha_s + alpha_m);

  assemble_dual_mass_lumped(linConductThermoLMglobal, linConductDISglobal);

  AssembleLinDM_X(&linConductDISglobal, nullptr, -delta_bar, LinDM_Diss, activenodes_);
  AssembleDM_linDiss(
      &linConductDISglobal, nullptr, &linConductContactLMglobal, nullptr, -delta_bar);

  AssembleDM_LMn(-beta_bar, &linConductTEMPglobal);
  assemble_lin_l_mn_dm_temp(-beta_bar, &linConductDISglobal, &linConductContactLMglobal);

  return;
}


void CONTACT::TSIInterface::assemble_dual_mass_lumped(
    Core::LinAlg::SparseMatrix& dualMassGlobal, Core::LinAlg::SparseMatrix& linDualMassGlobal)
{
  // loop over proc's slave nodes of the interface for assembly
  // use standard row map to assemble each node only once
  for (int i = 0; i < activenodes_->NumMyElements(); ++i)
  {
    int gid = activenodes_->GID(i);
    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    CONTACT::Node* conode = dynamic_cast<CONTACT::Node*>(node);

    if (conode->Owner() != Comm().MyPID())
      FOUR_C_THROW("AssembleDualMass: Node ownership inconsistency!");

    double thermo_lm = conode->TSIData().ThermoLM();
    std::map<int, std::map<int, double>>& derivDualMass = conode->Data().GetDerivD();

    if (Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(interface_params(), "LM_QUAD") !=
        Inpar::Mortar::lagmult_const)
    {
      /**********************************************dual mass matrix ******/
      if (conode->MoData().GetD().size() > 0)
      {
        const Core::Gen::Pairedvector<int, double>& dualMassmap = conode->MoData().GetD();
        Core::Gen::Pairedvector<int, double>::const_iterator colcurr;

        for (colcurr = dualMassmap.begin(); colcurr != dualMassmap.end(); ++colcurr)
        {
          Core::Nodes::Node* knode = Discret().gNode(colcurr->first);
          if (!knode) FOUR_C_THROW("node not found");
          Node* kcnode = dynamic_cast<Node*>(knode);
          if (!kcnode) FOUR_C_THROW("node not found");

          // create the mass matrix
          dualMassGlobal.FEAssemble(colcurr->second, conode->Dofs()[0], kcnode->Dofs()[0]);
        }
      }

      for (std::map<int, std::map<int, double>>::const_iterator a = derivDualMass.begin();
           a != derivDualMass.end(); ++a)
      {
        int sgid = a->first;
        Core::Nodes::Node* snode = idiscret_->gNode(sgid);
        if (!snode) FOUR_C_THROW("Cannot find node with gid %", sgid);
        Node* csnode = dynamic_cast<Node*>(snode);

        for (std::map<int, double>::const_iterator b = a->second.begin(); b != a->second.end(); ++b)
          // val               row                col
          linDualMassGlobal.FEAssemble(b->second * thermo_lm, csnode->Dofs()[0], b->first);
      }
    }

    else  // Inpar::Mortar::lagmult_const
    {
      if (conode->NumElement() != 1)
        FOUR_C_THROW(
            "some inconsistency: for lagmult_const every slave node may only have one element "
            "attached (it's the center-node!)");

      CONTACT::Element* coele = dynamic_cast<CONTACT::Element*>(conode->Elements()[0]);
      if (!coele) FOUR_C_THROW("this should be a contact element");

      Core::Gen::Pairedvector<int, double> derivArea(coele->num_node() * Dim());
      double area = coele->compute_area_deriv(derivArea);

      dualMassGlobal.FEAssemble(area, conode->Dofs()[0], conode->Dofs()[0]);

      for (Core::Gen::Pairedvector<int, double>::const_iterator p = derivArea.begin();
           p != derivArea.end(); ++p)
        linDualMassGlobal.FEAssemble(p->second * thermo_lm, conode->Dofs()[0], p->first);
    }
  }
  return;
}

void CONTACT::TSIInterface::AssembleLinDM_X(Core::LinAlg::SparseMatrix* linD_X,
    Core::LinAlg::SparseMatrix* linM_X, const double fac, const LinDmXMode mode,
    const Teuchos::RCP<Epetra_Map> node_rowmap)
{
  // get out if there's nothing to do
  if (linD_X == nullptr && linM_X == nullptr) return;

  // no dissipation without friction
  if (!friction_ && mode == LinDM_Diss) return;

  const double dt = interface_params().get<double>("TIMESTEP");

  // loop over all LM slave nodes (row map)
  for (int j = 0; j < node_rowmap->NumMyElements(); ++j)
  {
    int gid = node_rowmap->GID(j);
    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    // Mortar matrix D and M derivatives
    std::map<int, std::map<int, double>>& dderiv = cnode->Data().GetDerivD();
    std::map<int, std::map<int, double>>& mderiv = cnode->Data().GetDerivM();

    // current Lagrange multipliers
    double lm = 0.;
    switch (mode)
    {
      case LinDM_ThermoLM:
        lm = cnode->TSIData().ThermoLM();
        break;
      case LinDM_Diss:
      {
        CONTACT::FriNode* frnode = dynamic_cast<CONTACT::FriNode*>(cnode);
        if (frnode == nullptr)
          lm = 0.;
        else
        {
          double dval = 1.;
          if (cnode->MoData().GetD().size() == 0)
            continue;
          else
            cnode->MoData().GetD()[cnode->Id()];
          const Core::LinAlg::Matrix<3, 1> lmc(cnode->MoData().lm(), true);
          const Core::LinAlg::Matrix<3, 1> n(cnode->MoData().n(), true);
          const Core::LinAlg::Matrix<3, 1> jump(frnode->FriData().jump(), true);
          double diss = (-lmc.Dot(jump) + lmc.Dot(n) * jump.Dot(n)) / (dt * dval);
          lm = diss;
        }
        break;
      }
      case linDM_ContactLMnormal:
      {
        const Core::LinAlg::Matrix<3, 1> contact_LM(cnode->MoData().lm(), true);
        const Core::LinAlg::Matrix<3, 1> n(cnode->MoData().n(), true);
        lm = contact_LM.Dot(n);
        break;
      }
      default:
        FOUR_C_THROW("mode not implemented");
        break;
    }

    // get sizes and iterator start
    int slavesize = (int)dderiv.size();
    int mastersize = (int)mderiv.size();
    std::map<int, std::map<int, double>>::iterator scurr = dderiv.begin();
    std::map<int, std::map<int, double>>::iterator mcurr = mderiv.begin();

    /********************************************** LinDMatrix **********/
    // loop over all DISP slave nodes in the DerivD-map of the current LM slave node
    if (linD_X)
    {
      for (int k = 0; k < slavesize; ++k)
      {
        int sgid = scurr->first;
        ++scurr;

        Core::Nodes::Node* snode = idiscret_->gNode(sgid);
        if (!snode) FOUR_C_THROW("Cannot find node with gid %", sgid);
        Node* csnode = dynamic_cast<Node*>(snode);

        // Mortar matrix D derivatives
        std::map<int, double>& thisdderiv = cnode->Data().GetDerivD()[sgid];
        int mapsize = (int)(thisdderiv.size());

        // inner product D_{jk,c} * z_j for index j
        int row = csnode->Dofs()[0];
        std::map<int, double>::iterator scolcurr = thisdderiv.begin();

        // loop over all directional derivative entries
        for (int c = 0; c < mapsize; ++c)
        {
          int col = scolcurr->first;
          double val = lm * (scolcurr->second);
          ++scolcurr;

          // owner of LM slave node can do the assembly, although it actually
          // might not own the corresponding rows in lindglobal (DISP slave node)
          // (FE_MATRIX automatically takes care of non-local assembly inside!!!)
          // std::cout << "Assemble LinD: " << row << " " << col << " " << val << std::endl;
          if (abs(val) > 1.0e-12) linD_X->FEAssemble(val * fac, row, col);
        }

        // check for completeness of DerivD-Derivatives-iteration
        if (scolcurr != thisdderiv.end())
          FOUR_C_THROW("AssembleLinDM: Not all derivative entries of DerivD considered!");
      }

      // check for completeness of DerivD-Slave-iteration
      if (scurr != dderiv.end())
        FOUR_C_THROW("AssembleLinDM: Not all DISP slave entries of DerivD considered!");
    } /******************************** Finished with LinDMatrix **********/

    /********************************************** LinMMatrix **********/
    if (linM_X)
    {
      // loop over all master nodes in the DerivM-map of the current LM slave node
      for (int l = 0; l < mastersize; ++l)
      {
        int mgid = mcurr->first;
        ++mcurr;

        Core::Nodes::Node* mnode = idiscret_->gNode(mgid);
        if (!mnode) FOUR_C_THROW("Cannot find node with gid %", mgid);
        Node* cmnode = dynamic_cast<Node*>(mnode);

        // Mortar matrix M derivatives
        std::map<int, double>& thismderiv = cnode->Data().GetDerivM()[mgid];
        int mapsize = (int)(thismderiv.size());

        int row = cmnode->Dofs()[0];
        std::map<int, double>::iterator mcolcurr = thismderiv.begin();

        // loop over all directional derivative entries
        for (int c = 0; c < mapsize; ++c)
        {
          int col = mcolcurr->first;
          double val = lm * (mcolcurr->second);
          ++mcolcurr;

          // owner of LM slave node can do the assembly, although it actually
          // might not own the corresponding rows in lindglobal (DISP slave node)
          // (FE_MATRIX automatically takes care of non-local assembly inside!!!)
          // std::cout << "Assemble LinM: " << row << " " << col << " " << val << std::endl;
          if (abs(val) > 1.0e-12) linM_X->FEAssemble(-val * fac, row, col);
        }

        // check for completeness of DerivM-Derivatives-iteration
        if (mcolcurr != thismderiv.end())
          FOUR_C_THROW("AssembleLinDM: Not all derivative entries of DerivM considered!");
      }

      // check for completeness of DerivM-Master-iteration
      if (mcurr != mderiv.end())
        FOUR_C_THROW("AssembleLinDM: Not all master entries of DerivM considered!");
    } /******************************** Finished with LinMMatrix **********/
  }

  return;
}

void CONTACT::TSIInterface::AssembleDM_linDiss(Core::LinAlg::SparseMatrix* d_LinDissDISP,
    Core::LinAlg::SparseMatrix* m_LinDissDISP, Core::LinAlg::SparseMatrix* d_LinDissContactLM,
    Core::LinAlg::SparseMatrix* m_LinDissContactLM, const double fac)
{
  // get out if there's nothing to do
  if (d_LinDissDISP == nullptr && m_LinDissDISP == nullptr && d_LinDissContactLM == nullptr &&
      m_LinDissContactLM == nullptr)
    return;

  // there's no dissipation without friction
  if (!friction_) return;

  typedef std::map<int, double>::const_iterator _cim;
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _cip;

  const double dt = interface_params().get<double>("TIMESTEP");

  // loop over all LM slave nodes (row map)
  for (int j = 0; j < activenodes_->NumMyElements(); ++j)
  {
    int gid = activenodes_->GID(j);
    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);
    FriNode* fnode = dynamic_cast<FriNode*>(cnode);
    // if this is not a friction node, there is no dissipation so go to the next one
    if (fnode == nullptr) continue;

    // get nodal normal
    const Core::LinAlg::Matrix<3, 1> n(cnode->MoData().n(), true);

    // projection into tangential plane = 1 - n \otimes n
    Core::LinAlg::Matrix<3, 3> tang_proj(true);
    for (int i = 0; i < 3; ++i) tang_proj(i, i) = 1.;
    tang_proj.MultiplyNT(-1., n, n, 1.);

    // get D entry of this node
    // remember: D is diagonal
    int id = cnode->Id();
    const std::map<int, double>& derivD = cnode->Data().GetDerivD(id);
    const double dval = cnode->MoData().GetD()[cnode->Id()];

    // get nodal values
    Core::LinAlg::Matrix<3, 1> jump(fnode->FriData().jump());
    const Core::LinAlg::Matrix<3, 1> lm(cnode->MoData().lm());
    Core::LinAlg::Matrix<3, 1> lm_t;
    lm_t.Multiply(tang_proj, lm);
    const double lm_n = lm.Dot(n);
    const double jump_n = jump.Dot(n);
    Core::LinAlg::Matrix<3, 1> jump_tan;
    jump_tan.Multiply(tang_proj, jump);


    // linearization w.r.t. displacements
    if (d_LinDissDISP != nullptr || m_LinDissDISP != nullptr)
    {
      // get nodal jump and deriv
      const std::vector<std::map<int, double>>& derivJump = fnode->FriData().GetDerivJump();
      const std::vector<Core::Gen::Pairedvector<int, double>>& derivN = cnode->Data().GetDerivN();

      // calculate derivative of Dissipation *******************************
      std::map<int, double> derivDiss;
      for (int i = 0; i < 3; ++i)
      {
        for (_cim p = derivJump[i].begin(); p != derivJump[i].end(); ++p)
          derivDiss[p->first] -= ((lm(i) - lm_n * n(i)) * p->second) / (dt * dval);
        for (_cip p = derivN[i].begin(); p != derivN[i].end(); ++p)
          derivDiss[p->first] += ((lm_n * jump(i) + jump_n * lm(i)) * p->second) / (dt * dval);
        for (_cim p = derivD.begin(); p != derivD.end(); ++p)
          derivDiss[p->first] +=
              (-lm.Dot(jump) + lm.Dot(n) * jump.Dot(n)) / (dt * dval * dval) * (-p->second);
      }

      // put everything together*******************************************
      /**************************************************** D-matrix ******/
      if (d_LinDissDISP != nullptr)
        if ((cnode->MoData().GetD()).size() > 0)
          for (_cip k = cnode->MoData().GetD().begin(); k != cnode->MoData().GetD().end(); ++k)
          {
            Core::Nodes::Node* knode = Discret().gNode(k->first);
            if (!knode) FOUR_C_THROW("node not found");
            CONTACT::Node* kcnode = dynamic_cast<CONTACT::Node*>(knode);
            for (_cim currDeriv = derivDiss.begin(); currDeriv != derivDiss.end(); ++currDeriv)
              d_LinDissDISP->FEAssemble(
                  fac * k->second * currDeriv->second, kcnode->Dofs()[0], currDeriv->first);
          }
      /**************************************************** M-matrix ******/
      if (m_LinDissDISP != nullptr)
        if ((cnode->MoData().GetM()).size() > 0)
          for (_cim k = cnode->MoData().GetM().begin(); k != cnode->MoData().GetM().end(); ++k)
          {
            Core::Nodes::Node* knode = Discret().gNode(k->first);
            if (!knode) FOUR_C_THROW("node not found");
            CONTACT::Node* kcnode = dynamic_cast<CONTACT::Node*>(knode);
            for (_cim currDeriv = derivDiss.begin(); currDeriv != derivDiss.end(); ++currDeriv)
              m_LinDissDISP->FEAssemble(
                  fac * k->second * currDeriv->second, kcnode->Dofs()[0], currDeriv->first);
          }
    }  // linearization w.r.t. displacements

    // linearization wrt contact Lagrange multiplier
    if (d_LinDissContactLM != nullptr)
      // put everything together*******************************************
      /**************************************************** D-matrix ******/
      if ((cnode->MoData().GetD()).size() > 0)
        for (_cip k = cnode->MoData().GetD().begin(); k != cnode->MoData().GetD().end(); ++k)
        {
          Core::Nodes::Node* knode = Discret().gNode(k->first);
          if (!knode) FOUR_C_THROW("node not found");
          CONTACT::Node* kcnode = dynamic_cast<CONTACT::Node*>(knode);
          for (int d = 0; d < 3; ++d)
            d_LinDissContactLM->FEAssemble(
                -fac * k->second * jump_tan(d) / (dt * dval), kcnode->Dofs()[0], cnode->Dofs()[d]);
        }

    /**************************************************** M-matrix ******/
    if (m_LinDissContactLM)
      if ((cnode->MoData().GetM()).size() > 0)
        for (_cim k = cnode->MoData().GetM().begin(); k != cnode->MoData().GetM().end(); ++k)
        {
          Core::Nodes::Node* knode = Discret().gNode(k->first);
          if (!knode) FOUR_C_THROW("node not found");
          CONTACT::Node* kcnode = dynamic_cast<CONTACT::Node*>(knode);
          for (int d = 0; d < 3; ++d)
            m_LinDissContactLM->FEAssemble(
                -fac * k->second * jump_tan(d) / (dt * dval), kcnode->Dofs()[0], cnode->Dofs()[d]);
        }
  }  // loop over all active LM slave nodes (row map)
}

void CONTACT::TSIInterface::assemble_lin_l_mn_dm_temp(
    const double fac, Core::LinAlg::SparseMatrix* lin_disp, Core::LinAlg::SparseMatrix* lin_lm)
{
  // get out if there's nothing to do
  if (lin_disp == nullptr) FOUR_C_THROW("called to assemble something but didn't provide a matrix");

  typedef std::map<int, double>::const_iterator _cim;
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _cip;
  typedef std::map<int, std::map<int, double>>::const_iterator _cimm;

  // loop over all LM slave nodes (row map)
  for (int j = 0; j < activenodes_->NumMyElements(); ++j)
  {
    int gid = activenodes_->GID(j);
    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    const Core::LinAlg::Matrix<3, 1> n(cnode->MoData().n(), true);
    const Core::LinAlg::Matrix<3, 1> lm(cnode->MoData().lm(), true);
    const double lm_n = lm.Dot(n);

    for (_cimm k = cnode->Data().GetDerivD().begin(); k != cnode->Data().GetDerivD().end(); ++k)
    {
      Core::Nodes::Node* knode = Discret().gNode(k->first);
      if (!knode) FOUR_C_THROW("Cannot find node with gid %", gid);
      double temp_k = dynamic_cast<Node*>(knode)->TSIData().Temp();
      for (_cim l = k->second.begin(); l != k->second.end(); ++l)
        if (abs(l->second) > 1.e-12)
          lin_disp->FEAssemble(fac * lm_n * temp_k * l->second, cnode->Dofs()[0], l->first);
    }
    for (_cimm k = cnode->Data().GetDerivM().begin(); k != cnode->Data().GetDerivM().end(); ++k)
    {
      Core::Nodes::Node* knode = Discret().gNode(k->first);
      if (!knode) FOUR_C_THROW("Cannot find node with gid %", gid);
      double temp_k = dynamic_cast<Node*>(knode)->TSIData().Temp();
      for (_cim l = k->second.begin(); l != k->second.end(); ++l)
        if (abs(l->second) > 1.e-12)
          lin_disp->FEAssemble(-fac * lm_n * temp_k * l->second, cnode->Dofs()[0], l->first);
    }

    for (_cip k = cnode->MoData().GetD().begin(); k != cnode->MoData().GetD().end(); ++k)
    {
      Core::Nodes::Node* knode = Discret().gNode(k->first);
      if (!knode) FOUR_C_THROW("Cannot find node with gid %", gid);
      Node* cnodek = dynamic_cast<Node*>(knode);
      double temp_k = cnodek->TSIData().Temp();
      for (int d = 0; d < 3; ++d)
        for (_cip l = cnode->Data().GetDerivN()[d].begin(); l != cnode->Data().GetDerivN()[d].end();
             ++l)
          if (abs(l->second) > 1.e-12)
            lin_disp->FEAssemble(
                lm(d) * k->second * l->second * fac * temp_k, cnode->Dofs()[0], l->first);

      for (int d = 0; d < 3; ++d)
        lin_lm->FEAssemble(n(d) * temp_k * fac * k->second, cnode->Dofs()[0], cnode->Dofs()[d]);
    }

    for (_cim k = cnode->MoData().GetM().begin(); k != cnode->MoData().GetM().end(); ++k)
    {
      Core::Nodes::Node* knode = Discret().gNode(k->first);
      if (!knode) FOUR_C_THROW("Cannot find node with gid %", gid);
      Node* cnodek = dynamic_cast<Node*>(knode);
      double temp_k = cnodek->TSIData().Temp();
      for (int d = 0; d < 3; ++d)
        for (_cip l = cnode->Data().GetDerivN()[d].begin(); l != cnode->Data().GetDerivN()[d].end();
             ++l)
          if (abs(l->second) > 1.e-12)
            lin_disp->FEAssemble(
                -lm(d) * k->second * l->second * fac * temp_k, cnode->Dofs()[0], l->first);

      for (int d = 0; d < 3; ++d)
        lin_lm->FEAssemble(-n(d) * temp_k * fac * k->second, cnode->Dofs()[0], cnode->Dofs()[d]);
    }
  }

  return;
}

void CONTACT::TSIInterface::AssembleDM_LMn(const double fac, Core::LinAlg::SparseMatrix* DM_LMn)
{
  // get out if there's nothing to do
  if (DM_LMn == nullptr) FOUR_C_THROW("called to assemble something but didn't provide a matrix");

  typedef std::map<int, double>::const_iterator _cim;
  typedef Core::Gen::Pairedvector<int, double>::const_iterator _cip;

  // loop over all LM slave nodes (row map)
  for (int j = 0; j < activenodes_->NumMyElements(); ++j)
  {
    int gid = activenodes_->GID(j);
    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    const Core::LinAlg::Matrix<3, 1> n(cnode->MoData().n(), true);
    const Core::LinAlg::Matrix<3, 1> lm(cnode->MoData().lm(), true);
    const double lm_n = lm.Dot(n);

    for (_cip k = cnode->MoData().GetD().begin(); k != cnode->MoData().GetD().end(); ++k)
      if (abs(k->second) > 1.e-12)
      {
        Core::Nodes::Node* knode = Discret().gNode(k->first);
        if (!knode) FOUR_C_THROW("node not found");
        CONTACT::Node* kcnode = dynamic_cast<CONTACT::Node*>(knode);
        DM_LMn->FEAssemble(fac * lm_n * k->second, cnode->Dofs()[0], kcnode->Dofs()[0]);
      }

    for (_cim k = cnode->MoData().GetM().begin(); k != cnode->MoData().GetM().end(); ++k)
      if (abs(k->second) > 1.e-12)
      {
        Core::Nodes::Node* knode = Discret().gNode(k->first);
        if (!knode) FOUR_C_THROW("node not found");
        CONTACT::Node* kcnode = dynamic_cast<CONTACT::Node*>(knode);
        DM_LMn->FEAssemble(-fac * lm_n * k->second, cnode->Dofs()[0], kcnode->Dofs()[0]);
      }
  }
  return;
}


void CONTACT::TSIInterface::AssembleInactive(Core::LinAlg::SparseMatrix* linConductThermoLM)
{
  // get out if there's nothing to do
  if (linConductThermoLM == nullptr)
    FOUR_C_THROW("called to assemble something but didn't provide a matrix");

  // inactive nodes
  Teuchos::RCP<Epetra_Map> inactivenodes = Core::LinAlg::SplitMap(*snoderowmap_, *activenodes_);

  // loop over all LM slave nodes (row map)
  for (int j = 0; j < inactivenodes->NumMyElements(); ++j)
  {
    int gid = inactivenodes->GID(j);
    Core::Nodes::Node* node = idiscret_->gNode(gid);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
    Node* cnode = dynamic_cast<Node*>(node);

    if (cnode->Active()) FOUR_C_THROW("something with active and inactive maps is not in order");

    linConductThermoLM->FEAssemble(1., cnode->Dofs()[0], cnode->Dofs()[0]);
  }
}
/*----------------------------------------------------------------------*
 |  initialize / reset interface for tsi                     seitz 08/15|
 *----------------------------------------------------------------------*/
void CONTACT::TSIInterface::Initialize()
{
  CONTACT::Interface::Initialize();

  // loop over all nodes to reset stuff (fully overlapping column map)
  // (use fully overlapping column map)

  for (int i = 0; i < idiscret_->NumMyColNodes(); ++i)
  {
    CONTACT::Node* node = dynamic_cast<CONTACT::Node*>(idiscret_->lColNode(i));
    node->initialize_tsi_data_container(
        imortar_.get<double>("TEMP_REF"), imortar_.get<double>("TEMP_DAMAGE"));
    node->TSIData().Clear();
  }
  return;
}

FOUR_C_NAMESPACE_CLOSE
