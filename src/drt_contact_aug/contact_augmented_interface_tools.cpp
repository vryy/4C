/*---------------------------------------------------------------------*/
/*! \file
\brief Tools for the augmented contact interface evaluation.

\level 2

\maintainer Matthias Mayr
*/
/*---------------------------------------------------------------------*/
#include "contact_augmented_interface.H"
#include "contact_integrator_utils.H"
#include "../drt_contact/contact_integrator.H"
#include "../drt_contact/contact_defines.H"
#include "../drt_contact/contact_node.H"
#include "../drt_contact/contact_paramsinterface.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_mortar/mortar_dofset.H"
#include "../drt_mortar/mortar_integrator.H"
#include "../drt_mortar/mortar_defines.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_contact.H"

#include "../drt_structure_new/str_model_evaluator_data.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename T0, typename T1, typename T2>
CONTACT::AUG::Interface::FiniteDifference<T0, T1, T2>::FiniteDifference(
    Interface& interface, const double delta)
    : inter_(interface), stored_actiontype_(MORTAR::eval_none), ref_x_(true), delta_(delta)
{ /* empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename T0, typename T1, typename T2>
void CONTACT::AUG::Interface::FiniteDifference<T0, T1, T2>::GPCheck1stOrder(
    CONTACT::ParamsInterface& cparams)
{
  if (not inter_.lComm()) return;

  Epetra_Map& snode_row_map = *(inter_.snoderowmap_);
  Epetra_Map& mnode_row_map = *(inter_.mnoderowmap_);

  std::map<int, T0> ref;
  T0 new_ref;
  std::map<int, T1> ref_deriv1st;

  std::set<int> empty_node_gids;
  GetReference(snode_row_map, ref, ref_deriv1st, empty_node_gids);

  // nothing to do, if all nodes have no linearization
  if (empty_node_gids.size() >= static_cast<unsigned>(snode_row_map.NumMyElements())) return;

  RefEvaluate(cparams);

  const unsigned dim = inter_.Dim();

  Teuchos::RCP<Epetra_Map> fullmaps[2] = {Teuchos::null, Teuchos::null};

  fullmaps[0] = GetFullMap(snode_row_map);
  fullmaps[1] = GetFullMap(mnode_row_map);

  for (unsigned m = 0; m < 2; ++m)
  {
    for (unsigned fd = 0; fd < dim * fullmaps[m]->NumMyElements(); ++fd)
    {
      // store warnings for this finite difference
      int w = 0;

      // Initialize
      inter_.Initialize();

      // now get the node we want to apply the FD scheme to
      int gid = fullmaps[m]->GID(fd / dim);
      CoNode* snode = dynamic_cast<CoNode*>(inter_.idiscret_->gNode(gid));
      if (!snode) dserror("ERROR: Cannot find slave node with gid %", gid);

      const int sdof = snode->Dofs()[fd % dim];

      std::cout << "\nDERIVATIVE FOR " << (m == 0 ? "S" : "M") << "-NODE # " << gid
                << " DOF: " << sdof << std::endl;

      // do step forward (modify nodal displacement)
      DoPerturbation(*snode, fd);

      // *******************************************************************
      // contents of Evaluate()
      // *******************************************************************
      PerturbedEvaluate(cparams);

      // compute finite difference derivative
      for (int k = 0; k < snode_row_map.NumMyElements(); ++k)
      {
        const int kgid = snode_row_map.GID(k);

        DRT::Node* knode = inter_.idiscret_->gNode(kgid);
        if (!knode) dserror("ERROR: Cannot find node with gid %", kgid);
        CoNode* kcnode = dynamic_cast<CoNode*>(knode);

        if (GetDeriv1st(kcnode->AugData()).empty())
        {
          continue;
        }

        const T0& ref_vals = GetValues(kcnode->AugData());
        ResizeValues(new_ref);
        std::copy(&ref_vals, &ref_vals + GetNumVectors(ref_vals), &new_ref);

        for (unsigned v = 0; v < GetNumVectors(new_ref); ++v)
        {
          // print results (derivatives) to screen
          const int dof_gid = AtRef(new_ref, v).first;
          const double ref_val = AtRef(ref[kgid], v).second;

          const int var_node_gid = sdof / dim;
          const CoNode& var_node =
              dynamic_cast<const CoNode&>(*inter_.idiscret_->gNode(var_node_gid));
          const std::string sm(var_node.IsSlave() ? "s" : "m");

          if (std::abs(AtRef(new_ref, v).second - ref_val) > 1e-12)
          {
            const double fd_approx = (AtRef(new_ref, v).second - ref_val) / delta_;

            const double ref_d_var = AtDeriv1st(ref_deriv1st[kgid], v, sdof);
            const double dev = fd_approx - ref_d_var;

            /* output legend:
             * v    = vector component
             * kgid = node holding the derivatives
             * var_gid = global id of the varied dof (row ID in tangential matrix)
             * var_gid/dim = global id of the node holding the varied dof
             * s_dof = linearized dof (col ID in the tangential matrix */
            std::cout << "[" << v << "]"
                      << "{" << dof_gid << "," << kgid << "}"
                      << "(" << sdof << "(" << sdof / dim << "[" << sm << "]):";

            std::cout << "   fd=" << std::setw(14) << std::setprecision(5) << std::scientific
                      << fd_approx;
            std::cout << "   lin=" << std::setw(14) << std::setprecision(5) << std::scientific
                      << ref_d_var << "   DEVIATION= " << std::setw(14) << std::setprecision(5)
                      << std::scientific << dev;

            std::cout << "   REL-ERROR [%]= ";
            if (fd_approx != 0.0)
              std::cout << std::setw(14) << std::setprecision(5) << std::scientific
                        << abs(dev / fd_approx) * 100;
            else
              std::cout << " undef.";

            if (abs(dev) > 1e-4)
            {
              std::cout << " ***** WARNING ***** ";
              w++;
            }
            else if (abs(dev) > 1e-5)
            {
              std::cout << " ***** warning ***** ";
              w++;
            }

            std::cout << std::endl;
          }
        }
      }

      // undo finite difference modification
      UnDoPerturbation(*snode);

      std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** "
                << std::endl;
    }
  }

  FinalEvaluate(cparams);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename T0, typename T1, typename T2>
void CONTACT::AUG::Interface::FiniteDifference<T0, T1, T2>::GPCheck2ndOrder(
    CONTACT::ParamsInterface& cparams)
{
  if (not inter_.lComm()) return;

  Epetra_Map& node_row_map = *(inter_.snoderowmap_);

  std::map<int, T1> ref_deriv1st;
  T1 new_deriv1st;
  std::map<int, T2> ref_deriv2nd;

  RefEvaluate(cparams);

  std::set<int> empty_node_gids;
  GetReference(node_row_map, ref_deriv1st, ref_deriv2nd, empty_node_gids);

  // nothing to do, if all nodes have no linearization
  if (empty_node_gids.size() >= static_cast<unsigned>(node_row_map.NumMyElements())) return;

  const unsigned dim = inter_.Dim();

  Teuchos::RCP<Epetra_Map> fullmap = GetFullMap(node_row_map);

  for (unsigned fd = 0; fd < dim * fullmap->NumMyElements(); ++fd)
  {
    // store warnings for this finite difference
    int w = 0;

    // Initialize
    inter_.Initialize();

    // now get the node we want to apply the FD scheme to
    int gid = fullmap->GID(fd / dim);
    CoNode* snode = dynamic_cast<CoNode*>(inter_.idiscret_->gNode(gid));
    if (!snode) dserror("ERROR: Cannot find slave node with gid %", gid);

    const int sdof = snode->Dofs()[fd % dim];

    std::cout << "\nDERIVATIVE FOR S-NODE # " << gid << " DOF: " << sdof << std::endl;

    // do step forward (modify nodal displacement)
    DoPerturbation(*snode, fd);

    // *******************************************************************
    // contents of Evaluate()
    // *******************************************************************
    PerturbedEvaluate(cparams);

    // compute finite difference derivative
    for (int k = 0; k < node_row_map.NumMyElements(); ++k)
    {
      const int kgid = node_row_map.GID(k);

      // clear the calculated new r.h.s. map
      unsigned deriv1st_cap = GetCapacity(ref_deriv1st[kgid]);
      GEN_DATA::reset(deriv1st_cap, new_deriv1st);

      DRT::Node* knode = inter_.idiscret_->gNode(kgid);
      if (!knode) dserror("ERROR: Cannot find node with gid %", kgid);
      CoNode* kcnode = dynamic_cast<CoNode*>(knode);

      if (GetDeriv2nd(kcnode->AugData()).empty())
      {
        continue;
      }

      const T1& deriv1st = GetDeriv1st(kcnode->AugData());
      GEN_DATA::copy(deriv1st, new_deriv1st);

      for (unsigned v = 0; v < GetNumVectors(new_deriv1st); ++v)
      {
        for (auto& new_d_var : AtDeriv1st(new_deriv1st, v))
        {
          const int var_gid = new_d_var.first;
          const double ref_d_var = AtDeriv1st(ref_deriv1st[kgid], v, var_gid);

          const int var_node_gid = var_gid / dim;
          const CoNode& var_node =
              dynamic_cast<const CoNode&>(*inter_.idiscret_->gNode(var_node_gid));
          const std::string sm(var_node.IsSlave() ? "s" : "m");

          const std::string ai(kcnode->Active() ? "a" : "i");

          if (std::abs(new_d_var.second - AtDeriv1st(ref_deriv1st[kgid], v, var_gid)) > 1e-12)
          {
            const double fd_approx = (new_d_var.second - ref_d_var) / delta_;

            const double ref_dd_var = AtDeriv2nd(ref_deriv2nd[kgid], v, var_gid, sdof);
            const double dev = fd_approx - ref_dd_var;

            /* output legend:
             * v    = vector component
             * kgid = node holding the derivatives
             * var_gid = global id of the varied dof (row ID in tangential matrix)
             * var_gid/dim = global id of the node holding the varied dof
             * s_dof = linearized dof (col ID in the tangential matrix */
            std::cout << "[" << v << "]"
                      << "{" << kgid << "[" << ai << "]}"
                      << "(" << var_gid << "(" << var_gid / dim << "[" << sm << "])"
                      << "," << sdof
                      << ") :"
                         "   fd="
                      << std::setw(14) << std::setprecision(5) << std::scientific << fd_approx
                      << "   lin=" << std::setw(14) << std::setprecision(5) << std::scientific
                      << ref_dd_var << "   DEVIATION= " << std::setw(14) << std::setprecision(5)
                      << std::scientific << dev << "   REL-ERROR [%]= " << std::setw(14)
                      << std::setprecision(5) << std::scientific << abs(dev / fd_approx) * 100;

            if (abs(dev) > 1e-4)
            {
              std::cout << " ***** WARNING ***** ";
              w++;
            }
            else if (abs(dev) > 1e-5)
            {
              std::cout << " ***** warning ***** ";
              w++;
            }

            std::cout << std::endl;
          }
        }
      }
    }

    // undo finite difference modification
    UnDoPerturbation(*snode);

    std::cout << " ******************** GENERATED " << w << " WARNINGS ***************** "
              << std::endl;
  }

  FinalEvaluate(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename T0, typename T1, typename T2>
void CONTACT::AUG::Interface::FiniteDifference<T0, T1, T2>::DoPerturbation(
    CoNode& snode, const int fd_gid)
{
  const int dim = inter_.Dim();

  std::copy(snode.xspatial(), snode.xspatial() + 3, ref_x_.A());

  if (fd_gid % dim == 0)
  {
    snode.xspatial()[0] += delta_;
  }
  else if (fd_gid % dim == 1)
  {
    snode.xspatial()[1] += delta_;
  }
  else
  {
    snode.xspatial()[2] += delta_;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename T0, typename T1, typename T2>
void CONTACT::AUG::Interface::FiniteDifference<T0, T1, T2>::UnDoPerturbation(CoNode& snode)
{
  std::copy(ref_x_.A(), ref_x_.A() + 3, snode.xspatial());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename T0, typename T1, typename T2>
void CONTACT::AUG::Interface::FiniteDifference<T0, T1, T2>::RefEvaluate(
    CONTACT::ParamsInterface& cparams)
{
  stored_actiontype_ = cparams.GetActionType();

  SetActionType(MORTAR::eval_force_stiff, cparams);

  Evaluate(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename T0, typename T1, typename T2>
void CONTACT::AUG::Interface::FiniteDifference<T0, T1, T2>::PerturbedEvaluate(
    CONTACT::ParamsInterface& cparams)
{
  // compute element areas
  inter_.SetElementAreas();

  Evaluate(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename T0, typename T1, typename T2>
void CONTACT::AUG::Interface::FiniteDifference<T0, T1, T2>::FinalEvaluate(
    CONTACT::ParamsInterface& cparams)
{
  SetActionType(stored_actiontype_, cparams);

  Evaluate(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename T0, typename T1, typename T2>
void CONTACT::AUG::Interface::FiniteDifference<T0, T1, T2>::GetReference(
    const Epetra_Map& node_row_map, std::map<int, T0>& ref_map, std::map<int, T1>& ref_deriv1st,
    std::set<int>& empty_node_gids)
{
  int* mygids = node_row_map.MyGlobalElements();
  for (int i = 0; i < node_row_map.NumMyElements(); ++i)
  {
    const int gid = mygids[i];
    DRT::Node* node = inter_.idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);

    CoNode& cnode = static_cast<CoNode&>(*node);

    const T1& deriv1st = GetDeriv1st(cnode.AugData());

    if (deriv1st.empty())
    {
      empty_node_gids.insert(gid);
      continue;
    }

    const T0& ref = GetValues(cnode.AugData());
    if (AtRef(ref, 0).first == 0 and AtRef(ref, 0).second == 0.0)
    {
      empty_node_gids.insert(gid);
      continue;
    }


    T0& ref_vals = ref_map[gid];
    ResizeValues(ref_vals);

    std::copy(&ref, &ref + GetNumVectors(ref), &ref_vals);

    T1& ref_vars = ref_deriv1st[gid];
    GEN_DATA::copy(deriv1st, ref_vars);

    std::cout << "1-st order derivatives for node " << gid << "\n";
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename T0, typename T1, typename T2>
void CONTACT::AUG::Interface::FiniteDifference<T0, T1, T2>::GetReference(
    const Epetra_Map& node_row_map, std::map<int, T1>& ref_deriv1st,
    std::map<int, T2>& ref_deriv2nd, std::set<int>& empty_node_gids)
{
  int* mygids = node_row_map.MyGlobalElements();
  for (int i = 0; i < node_row_map.NumMyElements(); ++i)
  {
    const int gid = mygids[i];
    DRT::Node* node = inter_.idiscret_->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);

    CoNode& cnode = static_cast<CoNode&>(*node);

    const T2& deriv2nd = GetDeriv2nd(cnode.AugData());

    if (deriv2nd.empty())
    {
      empty_node_gids.insert(gid);
      continue;
    }

    const T1& deriv1st = GetDeriv1st(cnode.AugData());

    T1& ref_vals = ref_deriv1st[gid];
    GEN_DATA::copy(deriv1st, ref_vals);

    T2& ref_lins = ref_deriv2nd[gid];
    GEN_DATA::copy(deriv2nd, ref_lins);

    std::cout << "2-nd order derivatives for node " << gid << "\n";
    bool symm_check = INTEGRATOR::TestForSymmetry(ref_lins);
    if (not symm_check)
    {
      std::cout << "\nTest symmetry again with acitivated output ...\n";
      INTEGRATOR::TestForSymmetry(ref_lins, true);
      //      std::cout << "\n\n______________________________________\n";
      //      std::cout << "Print the 2-nd derivatives:\n";
      //      Print( ref_lins );
      //      dserror( "Symmetry check has not been passed!" );
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename T0, typename T1, typename T2>
void CONTACT::AUG::Interface::FiniteDifference<T0, T1, T2>::SetActionType(
    enum MORTAR::ActionType actiontype, CONTACT::ParamsInterface& cparams) const
{
  STR::MODELEVALUATOR::ContactData& cparams_mutable =
      dynamic_cast<STR::MODELEVALUATOR::ContactData&>(cparams);

  cparams_mutable.SetActionType(actiontype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename T0, typename T1, typename T2>
void CONTACT::AUG::Interface::FiniteDifference<T0, T1, T2>::Evaluate(
    CONTACT::ParamsInterface& cparams)
{
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr = Teuchos::rcpFromRef(cparams);

  // Initialize
  inter_.Initialize();

  // compute element areas
  inter_.SetElementAreas();

  // *******************************************************************
  // We have to do both evaluate steps here
  // *******************************************************************
  // evaluate averaged weighted gap
  inter_.Evaluate(cparams_ptr);
  // evaluate remaining entities and linearization
  inter_.RedEvaluate(cparams_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename T0, typename T1, typename T2>
Teuchos::RCP<Epetra_Map> CONTACT::AUG::Interface::FiniteDifference<T0, T1, T2>::GetFullMap(
    const Epetra_Map& rowmap) const
{
  return LINALG::AllreduceEMap(rowmap);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::pair<int, double>& CONTACT::AUG::Interface::FD_Debug::GetValues(
    const NodeDataContainer& data) const
{
  return data.Get_Debug();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const CONTACT::AUG::Deriv1stMap& CONTACT::AUG::Interface::FD_Debug::GetDeriv1st(
    const NodeDataContainer& data) const
{
  return data.GetDeriv1st_Debug();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const CONTACT::AUG::Deriv2ndMap& CONTACT::AUG::Interface::FD_Debug::GetDeriv2nd(
    const NodeDataContainer& data) const
{
  return data.GetDeriv2nd_Debug();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::vector<std::pair<int, double>>& CONTACT::AUG::Interface::FD_DebugVec::GetValues(
    const NodeDataContainer& data) const
{
  return data.Get_DebugVec();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const CONTACT::AUG::Deriv1stVecMap& CONTACT::AUG::Interface::FD_DebugVec::GetDeriv1st(
    const NodeDataContainer& data) const
{
  return data.GetDeriv1st_DebugVec();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const CONTACT::AUG::Deriv2ndVecMap& CONTACT::AUG::Interface::FD_DebugVec::GetDeriv2nd(
    const NodeDataContainer& data) const
{
  return data.GetDeriv2nd_DebugVec();
}

template class CONTACT::AUG::Interface::FiniteDifference<std::pair<int, double>,
    CONTACT::AUG::Deriv1stMap, CONTACT::AUG::Deriv2ndMap>;
template class CONTACT::AUG::Interface::FiniteDifference<std::vector<std::pair<int, double>>,
    CONTACT::AUG::Deriv1stVecMap, CONTACT::AUG::Deriv2ndVecMap>;
