#pragma once

#include <numeric>


#include <ranges>
#include "pnhda.h"
#include "utils.h"

namespace convmaxhda
{

/**
 *
 * @tparam CONC conclist type
 * @return true iff \a l1 union \a l2 == \a r1 union \a r2
 * taken as multisets
 */
template <class CONC>
bool conc_union_equiv(const CONC& l1,const CONC& l2,
                      const CONC& r1,const CONC& r2)
{
  if (l1.size() +  l2.size() != r1.size() + r2.size())
    return false;

  auto l1it = l1.begin();
  auto l2it = l2.begin();
  auto r1it = r1.begin();
  auto r2it = r2.begin();
  auto l1e = l1.end();
  auto l2e = l2.end();
  auto r1e = r1.end();
  auto r2e = r2.end();

  while ((l1it != l1e || l2it != l2e) && (r1it != r1e || r2it != r2e))
  {
    auto& lit = (l2it == l2e) || ((l1it != l1e) && (*l1it < *l2it)) ? l1it : l2it;
    auto& rit = (r2it == r2e) || ((r1it != r1e) && (*r1it < *r2it)) ? r1it : r2it;

    // Advance both if equal
    if (*lit == *rit)
    {
      ++lit;
      ++rit;
    }
    else
      break;
  }
  return l1it == l1e && l2it == l2e && r1it == r1e && r2it == r2e;
}

/**
 *
 * @tparam CONC conclist type
 * @return true iff \a l1 union l2 subseteq r1 union r2
 * taken as multisets
 */
template <class CONC>
bool conc_union_inclusion(const CONC& l1,const CONC& l2,
                      const CONC& r1,const CONC& r2)
{
  if (l1.size() +  l2.size() > r1.size() + r2.size())
    return false;

  auto l1it = l1.begin();
  auto l2it = l2.begin();
  auto r1it = r1.begin();
  auto r2it = r2.begin();
  auto l1e = l1.end();
  auto l2e = l2.end();
  auto r1e = r1.end();
  auto r2e = r2.end();

  while ((l1it != l1e || l2it != l2e) && (r1it != r1e || r2it != r2e))
  {
    auto& lit = (l2it == l2e) || ((l1it != l1e) && (*l1it < *l2it)) ? l1it : l2it;
    auto& rit = (r2it == r2e) || ((r1it != r1e) && (*r1it < *r2it)) ? r1it : r2it;

    // Advance both if equal
    if (*lit == *rit)
    {
      ++lit;
      ++rit;
    }
    else if (*lit > *rit) // Skip on the right
      ++rit;
    else
      break; //Element left too small
  }
  return l1it == l1e && l2it == l2e;
}

/**
 * @tparam CONC conclist type
 * @return true iff \a l subseteq \a r
 * taken as multiset
 */
template <class CONC>
bool conc_inclusion(const CONC& l, const CONC& r)
{
  if (l.size() > r.size())
    return false;

  auto lit = l.begin();
  auto rit = r.begin();
  auto le = l.end();
  auto re = r.end();

  while (lit != le && rit != re)
  {
    // Advance both if equal
    if (*lit == *rit)
    {
      ++lit;
      ++rit;
    }
    else if (*lit > *rit) // Skip on the right
      ++rit;
    else
      break; // Element left too small
  }
  return lit == le;
}


/**
 *
 * @tparam CONC
 * @param ltc_src
 * @param lsc
 * @param ltc_dst
 * @param rtc_src
 * @param rsc
 * @param rtc_dst
 * @return true iff left is "covered" by right: at least as large interface, at most as large to corner
 */
template <class CONC>
bool edge_conc_sub(const CONC& ltc0_src, const CONC& ltc1_src,
                   const CONC& lsc,
                   const CONC& ltc0_dst, const CONC& ltc1_dst,
                   const CONC& rtc0_src, const CONC& rtc1_src,
                   const CONC& rsc,
                   const CONC& rtc0_dst, const CONC& rtc1_dst)
{
  return conc_inclusion(lsc, rsc)
         && conc_inclusion(rtc0_src, ltc0_src)
         && conc_inclusion(rtc1_src, ltc1_src)
         && conc_inclusion(rtc0_dst, ltc0_dst)
         && conc_inclusion(rtc1_dst, ltc1_dst);
}

/**
 *
 * @tparam HDAT Base type of the HDA
 * Class representing HDAs which were converted from
 * stadard HDAs. The transitions between max-cells are labelled
 * with the shared conclist, as well as the lower and upper
 * face maps that have to be taken in the source and
 * destination cell to attain the shared cell.
 */
template <class HDAT>
class basic_max_conv_hda
{
public:
  /**
   * @addtogroup PN definitions
   * @{
   */
  using pnames_t = HDAT::pnames_t;
  using pnamesid_t = HDAT::pnamesid_t;
  using tnames_t = HDAT::tnames_t;
  using tnamesid_t = HDAT::tnamesid_t;
  using marking_t = HDAT::marking_t;
  using transition_t = HDAT::transition_t;
  //using all_trans_t = HDAT::all_trans_t;
  using place2id_t = HDAT::place2id_t;
  using trans2id_t = HDAT::trans2id_t;
  /** @} */

  /**
   * @addtogroup PNHDA definitions
   * @{
   */
  using cid_t = HDAT::cid_t;
  using conc_t = HDAT::conc_t;
  using concid_t = HDAT::concid_t;
  using markingid_t = HDAT::markingid_t;
  using config_t = HDAT::config_t; // The meaning is different: initial marking and max-conclist
  using fm_t = HDAT::fm_t;
  //using get_dim = HDAT::get_dim;
  struct edge_t
  { // todo I do not think that this is the optimal way to do it
    int edge_type;
    std::array<concid_t, 2> to_corner_src; // 0 -> unstarted, 1 -> terminated
    concid_t shared_conc; // This one is actually superfluous, can be deduced from the others
    std::array<concid_t, 2> to_corner_dst; // 0 -> unstarted, 1 -> terminated
    cid_t dst;
  };
  /** @} */

private:
  bimap<conc_t, concid_t> square_ess_; /** Bimap for reachable conclists (include those necessary or edges */
  bimap<marking_t, markingid_t> mark_ess_; /** Bimap for reachable markings */

  bimap<config_t, cid_t> X_ext_; /** bimap between configurations and cells (one-to-one) */

  /** Todo inefficient */
  std::map<cid_t, std::vector<edge_t>> edges_; /** Instead of facemaps, edges are used */

  std::set<cid_t> init_set_; /** The \bot set */
  std::set<cid_t> final_set_; /** The \top set */

  place2id_t place2id_; /** Copied labeling from the underlying petri net for states */
  trans2id_t trans2id_; /** Copied labeling from the underlying petri net for transitions */

  unsigned nmarkings_;

  void convert_me_(const auto& hda)
  {
    // Get all the max-cells and rename them
    auto all_cells_ismax = std::vector<bool>(hda.X_ext().get_fwd().size(), true);
    for (const auto& m : {&hda.delta_0(), &hda.delta_1()})
      for (const auto& t : *m)
        for (const auto& tt : t.second)
        {
          if (tt.second >= all_cells_ismax.size())
            all_cells_ismax.resize(2*tt.second, true);
          all_cells_ismax[tt.second] = false;
        }
    // Do not rename
    // Copy all the necessary conclists and markings
    // They do not need to be renamed
    for (const auto& [conf, id]: hda.X_ext().get_fwd())
    {
      if (!all_cells_ismax[id])
        continue;
      X_ext_.try_emplace(conf, id);
      mark_ess_.try_emplace(conf.first, hda.mark_ess().at(conf.first));
      square_ess_.try_emplace(conf.second, hda.square_ess().at(conf.second));
    }

    // Check for shared subcells
    // Type 0 transitions are asym, Type 1 are sym
    // We need to all unordered pairs
    // Either its type 1 -> sym
    // Or type 0 and we can derive who is source and who is dst
    // This only works when the hda is full
    // I.e. it needs !!l0it <-> !!l1it
    struct cell_pair_iterator
    {
      using iter_t = std::set<fm_t>::iterator;
      using iter_opt_t = std::optional<iter_t>;

      cid_t left;
      cid_t right;
      iter_opt_t l0e;
      iter_opt_t l0it;
      iter_opt_t l1it;
      iter_opt_t l1e;
      iter_opt_t r0it;
      iter_opt_t r0e;
      iter_opt_t r1it;
      iter_opt_t r1e;
      iter_opt_t lit;
      iter_opt_t le;
      iter_opt_t rit;
      iter_opt_t re;
      bool left_orig = true;
      bool right_orig = false;
      bool llower = true;
      bool rlower = true;

      cell_pair_iterator(cid_t l, cid_t r, const HDAT& hda)
        : left{l}// l < r ? l : r}
        , right{r}// l < r ? r : l}
      {
        //std::cout << "cstr " << l << "; " << r << std::endl;
        // Get the iterators if available
        if (auto ite = hda.delta_0().find(left); ite != hda.delta_0().end())
        {
          l0it = ite->second.begin();
          l0e = ite->second.end();
        }
        if (auto ite = hda.delta_1().find(left); ite != hda.delta_1().end())
        {
          l1it = ite->second.begin();
          l1e = ite->second.end();
        }
        if (auto ite = hda.delta_0().find(right); ite != hda.delta_0().end())
        {
          r0it = ite->second.begin();
          r0e = ite->second.end();
        }
        if (auto ite = hda.delta_1().find(right); ite != hda.delta_1().end())
        {
          r1it = ite->second.begin();
          r1e = ite->second.end();
        }
        lit = l0it;
        le = l0e;
        rit = r0it;
        re = r0e;

        // Forward to ok position
        // Can not do left_orig if right is empty
        if (rit == re)
        {
          // Forward to right orig
          left_orig = false;
          right_orig = true;
          llower = true;
          rlower = true;
          // Can not do right_orig if left is empty
          if (lit == le)
          {
            left_orig = false;
            right_orig = false;
            llower = false;
            rlower = false; // No successors at all
          }
        }
      }

      struct ret_val
      {
        cid_t left; /** left cid */
        cid_t right; /** right cid */
        std::optional<tnamesid_t> ltr; /** concerned transition left */
        std::optional<tnamesid_t> rtr; /** concerned transition right */
        bool llower; /** Lower or upper face */
        bool rlower; /** Lower or upper face */
      };

      /**
       * @brief Returns informations about the next successor couple
       */
      ret_val get() const
      {
        //cid_t ll, rr;
        //std::optional<tnamesid_t> lltr, rrtr;
        //bool lllower, rrlower;

        if (left_orig)
        {
          assert(llower);
          assert(!right_orig);
          //ll = left;
          //rr = (*rit)->second;
          //lltr = std::nullopt;
          //rrtr = (*rit)->first;
          //lllower = llower;
          //rrlower = rlower;
          return {left, (*rit)->second, std::nullopt, (*rit)->first, llower, rlower};
        }
        else if (right_orig)
        {
          assert(rlower);
          assert(!left_orig);
          return {(*lit)->second, right, (*lit)->first, std::nullopt, llower, rlower};
          //ll = (*lit)->second;
          //rr = right;
          //lltr = (*lit)->first;
          //rrtr = std::nullopt;
          //lllower = llower;
          //rrlower = rlower;
        }
        else
        {
          assert(lit && lit != le);
          assert(rit && rit != re);
          return {(*lit)->second, (*rit)->second, (*lit)->first, (*rit)->first, llower, rlower};
          //ll = (*lit)->second;
          //rr = (*rit)->second;
          //lltr = (*lit)->first;
          //rrtr = (*rit)->first;
          //lllower = llower;
          //rrlower = rlower;
        }
        //if (ll < rr)
        //  return {ll, rr, lltr, rrtr, lllower, rrlower};
        //else
        //  return {rr, ll, rrtr, lltr, rrlower, lllower};
      }

      operator bool() const
      {
        // Attention these are optional
        // Either we iterate over orig -> only one need to be non-empty
        // Or we iterate over coupled -> both need to be non-empty
        return ((!!lit && lit != le) && (!!rit && rit != re))
               || (left_orig && (!!rit && rit != re))
               || (right_orig && (!!lit && lit != le));
      }

      bool right_incr_()
      {
        // Regular advance
        if ((rit != re) && (++*rit != re))
          return true;

        // Switch to upper
        if (rlower)
        {
          rlower = false;
          rit = r1it;
          re = r1e;
          return rit != re;
        }
        // Can't advance
        return false;
      }
      bool left_incr_()
      {
        // Regular advance
        if ((lit != le) && (++*lit != le))
          return true;

        // Switch to upper
        if (llower)
        {
          llower = false;
          lit = l1it;
          le = l1e;
          return lit != le;
        }
        // Can't advance
        return false;
      }

      bool next()
      {
        if (!(*this))
          return false; // Already exhausted

        // We are iterating over left orig
        if (left_orig)
        {
          if (right_incr_())
            return true;
          // Switch to right orig
          left_orig = false;
          right_orig = true;
          llower = true;
          rlower = true;
          lit = l0it;
          le = l0e;
          rit = r0it;
          re = r0e;
          return !!(*this);
        }
        if (right_orig)
        {
          if (left_incr_())
            return true;
          // Switch to combined iteration
          left_orig = false;
          right_orig = false;
          llower = true;
          rlower = true;
          lit = l0it;
          le = l0e;
          rit = r0it;
          re = r0e;
          return !!(*this);
        }

        // Early exit, advance only one
        left_orig = false;
        right_orig = false;
        llower = false;
        rlower = false;
        lit = l1e;
        le = l1e;
        rit = r1e;
        re = r1e;

        return false;

        //// Going down the dimensions
        //if (right_incr_())
        //  return true;
        //// rit was exhausted -> Try advance left
        //if (left_incr_())
        //{
        //  rit = r0it;
        //  re = r0e;
        //  rlower = true;
        //  return !!(*this);
        //}
        //// Left and right are exhausted
        //assert(lit == le); // Optional not present or exhausted
        //assert(rit == re);
        //assert(!llower && !rlower);
        //assert(!left_orig && !right_orig);
        //return false;
      }
    };

    auto share_cell =
      [&, stack = std::vector<cell_pair_iterator>(), matched = std::set<cid_t>(),
       seen = std::set<std::pair<cid_t, cid_t>>()](cid_t left, cid_t right) mutable
    {
      // Reset matching
      matched.clear();
      seen.clear();
      // Cells that have already been matched and do not need further investigation
      stack.emplace_back(left, right, hda);
      //if (left < right)
      //  stack.emplace_back(left, right, hda);
      //else
      //  stack.emplace_back(right, left, hda);

      auto decl_seen = [&seen](cid_t l, cid_t r) -> void
      {
        seen.emplace(l < r ? l : r, l < r ? r : l);
      };
      auto check_seen = [&seen](cid_t l, cid_t r) -> bool
      {
        return (bool) seen.count({l < r ? l : r, l < r ? r : l});
      };

      // 0 -> unstarted, 1 -> terminated
      std::array<conc_t, 2> lconc;
      std::array<conc_t, 2> rconc;

      /** Returns (was_already_marked, revisits premarked)
      * Due to "pseudo auto-concurrency" we can not simply
      * skip cells that are already marked as they might have
      * different predescessors
      */
      auto mark_me = [&](cid_t c, auto&& mark_me_) -> std::pair<bool, bool>
        {
          if (matched.insert(c).second)
          {
            bool remarks = false;
            // A new cell, so we need to recurse
            for (const auto& m : {&hda.delta_0(), &hda.delta_1()}) // Take ptr, otherwise it is a copy!!
              if (auto it = m->find(c); it != m->end())
                for (const auto& t : it->second)
                  remarks |= mark_me_(t.second, mark_me_).first;
            return {false, remarks};
          }
          return {true, true};
        };

      while (!stack.empty())
      {
        // Check if we look at the same cells
        const auto& curr = stack.back();
        // Check if current was already explored and skip
        if (check_seen(curr.left, curr.right))
        {
          stack.pop_back();
          if (!stack.empty())
            stack.back().next(); // Advance
          continue;
        }

        if (curr.left == curr.right)
        {
          // We have found a match
          // Mark all matched cells
          auto [m_covered, m_covering] = mark_me(curr.left, mark_me);
          //std::cout << m_covered << m_covering << std::endl;
          // There is a bug somewhere that causes inclusions to be missed
          m_covered = true;
          m_covering = true;

          // Reconstruct the transitions
          for (auto i : {0, 1})
          {
            lconc[i].clear();
            rconc[i].clear();
          }
          bool llower_ex = true;
          bool rlower_ex = true;

          // last elem
          auto last = stack.back(); // Do not touch "curr" from here on, only last
          // We have finished exploring this pair -> "see it"
          decl_seen(last.left, last.right);
          // Remove it
          stack.pop_back();
          // Take into account all transitions
          for (const auto& t : stack)
          {
            const auto& tn = t.get();
            llower_ex &= tn.llower;
            rlower_ex &= tn.rlower;
            // add transition that NEED to be unstarted or terminated
            if (tn.ltr)
              lconc[!tn.llower].push_back(*tn.ltr);
            if (tn.rtr)
              rconc[!tn.rlower].push_back(*tn.rtr);
          }

          // We need to sort them, they may be out of order
          for (auto i : {0, 1})
          {
            std::ranges::sort(lconc[i]);
            std::ranges::sort(rconc[i]);
          }

          // The type depends on the lower_ex of the other
          auto get_conc = [&](const conc_t& cc) -> concid_t
            {
              auto new_id = square_ess_.get_rev().rbegin()->first + 1;
              return square_ess_.try_emplace(cc, new_id).first->second;
            };
          auto get_concp = [&](const std::array<conc_t, 2>& c) -> std::array<concid_t, 2>
            {
              return {get_conc(c[0]), get_conc(c[1])};
            };
          // Shared conc
          concid_t scid = get_conc(hda.square_ess().at(hda.X_ext().at(last.left).second)); // We need the conc of the shared cell

          // Before inserting check if covered
          // If the marking was already covered we need to check
          // if there exists a "covering" transition
          bool covered = false;
          if (m_covered)
          {
            auto itl = edges_.find(left);
            if (itl != edges_.end())
            {
              for (const auto& e : itl->second)
              {
                if (right == e.dst
                    && edge_conc_sub(lconc[0], lconc[1],
                                     square_ess_.at(scid),
                                     rconc[0], rconc[1],
                                     square_ess_.at(e.to_corner_src[0]), square_ess_.at(e.to_corner_src[1]),
                                     square_ess_.at(e.shared_conc),
                                     square_ess_.at(e.to_corner_dst[0]), square_ess_.at(e.to_corner_dst[1])))
                {
                  covered = true;
                  break;
                }
              }
            }
            // Since edges are symmetric
            // We can inverse left and right
            auto itr = edges_.find(right);
            if (!covered && itr != edges_.end())
            {
              for (const auto& e : itr->second)
              {
                if (left == e.dst
                    && edge_conc_sub(rconc[0], rconc[1],
                                     square_ess_.at(scid),
                                     lconc[0], lconc[1],
                                     square_ess_.at(e.to_corner_src[0]), square_ess_.at(e.to_corner_src[1]),
                                     square_ess_.at(e.shared_conc),
                                     square_ess_.at(e.to_corner_dst[0]), square_ess_.at(e.to_corner_dst[1])))
                {
                  covered = true;
                  break;
                }
              }
            }

          }
          if (m_covering && !covered) // If this trans is covered, all those that it would cover are already treated
          {
            // If an edge is covered by this one -> remove it
            auto itl = edges_.find(left);
            if (itl != edges_.end())
            {
              // Make sure that iterators remain valid
              // -> backward traversal by idx
              // The edges themselves are not ordered
              for (auto ite = itl->second.rbegin(); ite != itl->second.rend(); ++ite)
              {
                if (right == ite->dst
                    && edge_conc_sub(square_ess_.at(ite->to_corner_src[0]), square_ess_.at(ite->to_corner_src[1]),
                                     square_ess_.at(ite->shared_conc),
                                     square_ess_.at(ite->to_corner_dst[0]), square_ess_.at(ite->to_corner_dst[1]),
                                     lconc[0], lconc[1],
                                     square_ess_.at(scid),
                                     rconc[0], rconc[1]))
                {
                  // Remove me here
                  // reverse iter -> ok for removal
                  *ite = itl->second.back();
                  itl->second.pop_back();
                }
              }
            }
            // Since edges are symmetric
            // We can inverse left and right
            auto itr = edges_.find(right);
            if (itr != edges_.end())
            {
              // Make sure that iterators remain valid
              // -> backward traversal by idx
              // The edges themselves are not ordered
              for (auto ite = itr->second.rbegin(); ite != itr->second.rend(); ++ite)
              {
                if (left == ite->dst
                    && edge_conc_sub(square_ess_.at(ite->to_corner_src[0]), square_ess_.at(ite->to_corner_src[1]),
                                     square_ess_.at(ite->shared_conc),
                                     square_ess_.at(ite->to_corner_dst[0]), square_ess_.at(ite->to_corner_dst[1]),
                                     rconc[0], rconc[1],
                                     square_ess_.at(scid),
                                     lconc[0], lconc[1]))
                {
                  // Remove me here
                  // reverse iter -> ok for removal
                  *ite = itr->second.back();
                  itr->second.pop_back();
                }
              }
            }
          }

          if (!covered)
          {
            // left to right
            auto& eleft = edges_[left];
            eleft.emplace_back((int) rlower_ex, get_concp(lconc), scid, get_concp(rconc), right);
            // right to left
            // Right to left can either be constructed afterwards
            // Or is implicit
            //auto& eright = edges_[right];
            //eright.emplace_back((int) llower_ex, get_conc(rconc), scid, get_conc(lconc), left);
          }
          // done
          // -> Advance the parent
          if (!stack.empty())
            stack.back().next();
          continue;
        }
        // Push the successor or pop if empty
        if (curr)
        {
          const auto& nxt = curr.get();
          //if (!matched.count(nxt.left) && !matched.count(nxt.right)) This is too strong
          stack.emplace_back(nxt.left, nxt.right, hda);
        } else {
          // We have exhausted
          decl_seen(stack.back().left, stack.back().right);
          stack.pop_back();
          if (!stack.empty())
            stack.back().next(); // Advance
        }
      }

      stack.clear();
      matched.clear();
      seen.clear();
    };

    for (auto it1 = X_ext_.get_fwd().begin(); it1 != X_ext_.get_fwd().end(); ++it1)
    {
      auto it2 = it1;
      for (++it2; it2 != X_ext_.get_fwd().end(); ++it2)
        share_cell(it1->second, it2->second);
    }
  }

  std::string get_csv_base_header_() const { return "ncells,nconclist,nmarkings,nEdges,maxdim,L(dim;ncell)"; }

    std::string get_csv_base_() const
    {
      auto ss = std::ostringstream();

      ss << X_ext_.size() << ',' << square_ess_.size() << ',' << mark_ess_.size() << ','
         << accu_size_sec(edges_) << ',';

      // dimension of cell and how many
      auto celldet = std::map<uint, uint>();

      for (const auto& [conf, _] : X_ext_.get_fwd())
      {
        concid_t cid = conf.second;
        const auto& concl = square_ess_.at(cid);
        uint dim = concl.size();//get_dim(concl);

        auto it = celldet.try_emplace(dim, 0).first;
        ++(it->second);
      }

      ss << celldet.rbegin()->first << ",[";

      for (const auto& [dim, ndim] : celldet)
        ss << '[' << dim << ';' << ndim << "];";
      ss << "]";

      return ss.str();
    }


    std::string to_dot_(const std::string& fine_opt) const
    {
      // bool legend_left = fine_opt.contains("legend=left;");
      std::string dot_dir = "TB";
      if (fine_opt.contains("LR;"))
        dot_dir = "LR";

      uint nCols = 3;
      if (auto sidx = fine_opt.find("nCols"); sidx != std::string::npos)
      {
        auto fidx = fine_opt.find(";", sidx);
        if (fidx == std::string::npos)
          throw std::invalid_argument("Option does not terminate with ;");
        std::cerr << fine_opt.substr(sidx + 6, fidx - sidx - 6) << std::endl;
        nCols = std::stol(fine_opt.substr(sidx + 6, fidx - sidx - 6));
      }
      //int print_ep = 0;
      //if (auto sidx = fine_opt.find("printEp"); sidx != std::string::npos)
      //{
      //  auto fidx = fine_opt.find(";", sidx);
      //  if (fidx != sidx + 9)
      //    throw std::invalid_argument("Option does not terminate with ; or does not use 0,1 to specify value");
      //  std::cerr << fine_opt.substr(sidx + 8, sidx+9) << std::endl;
      //  print_ep = std::stoi(fine_opt.substr(sidx + 8, sidx + 9));
      //}


      std::string rnodes = "rank=source;\n";
      std::string rcnodes = "rank=source;\n";

      auto ss = std::ostringstream();

      ss << "digraph automaton {\n"
         << "rankdir=" << dot_dir << ";  // Ensure top-to-bottom layout\n"
         << "// Subgraph for node identification\n";

      // Already define all nodes
      // Put all the actual states
      ss << "node [shape = circle];\n"
         << "nodesep=0.15;\n"
         << "ranksep=0.35;\n"
         << rnodes;


      for (const auto& p : X_ext_.get_rev())
      {
        auto ns = "N" + std::to_string(p.first);
        ss << ns << " [label = \"" << ns << "\"];\n";
      }
      // Continue with the subgraph
      ss << "subgraph cluster_comments {\n"
         << rcnodes << "node [shape=none];\n"
         << "edge [style=invis];\n"
         << "nodesep=0.15;\n"
         << "ranksep=0.35;\n";
      // Create all the nodes
      // todo factorize me
      auto fmtnode = [&](cid_t cid)
      {
        auto [markid, concid] = X_ext_.at(cid);
        const auto& mark = mark_ess_.at(markid);
        const auto& conc = square_ess_.at(concid);

        auto concstream = std::ostringstream();
        auto concstr = std::string();
        auto markstream = std::ostringstream();
        auto markstr = std::string();

        markstream.str("");
        for (auto m : mark)
          markstream << m << ',';
        markstr = markstream.str();
        if (!markstr.empty())
          markstr.pop_back();
        markstr = "[" + markstr + "]";

        concstream.str("");
        for (auto tid : conc)
          concstream << trans2id_.at(tid) << ',';
        concstr = concstream.str();
        if (!concstr.empty())
          concstr.pop_back();
        concstr = "[" + concstr + "]";

        return markstr + "x" + concstr;
      };
      auto fullfmtnode = [&](cid_t cid)
      { return "<TD ALIGN=\"LEFT\">N" + std::to_string(cid) + ": " + fmtnode(cid) + "</TD>"; };

      auto cnNodes = std::vector<std::string>(nCols, "<TR>");
      {
        std::size_t i = 0;
        for (auto it = X_ext_.get_rev().begin(); it != X_ext_.get_rev().end(); ++i, ++it)
        {
          cnNodes.at(i % nCols) += fullfmtnode(it->first);
        }
        for (; i % nCols != 0; ++i)
          cnNodes.at(i % nCols) += "<TD ALIGN=\"LEFT\"> </TD>";
        for (auto& s : cnNodes)
          s += "</TR>";
      }
      ss << "CN0 [label=<<TABLE BORDER=\"0\" CELLBORDER=\"0\" CELLSPACING=\"0\">";
      for (const auto& s : cnNodes)
        ss << s;
      ss << "</TABLE>>];\n";

      ss << "}\n";
      // And all the transitions
      auto fmtedge = [&](const auto& e)
        {
          auto es = std::ostringstream();
          //es << "[style = " << (e.edge_type == 0 ? "solid" : "dashed") << ", label = \"[(";
          es << " [style = solid, label = \"[(";
          es << join(square_ess().at(e.to_corner_src[0])) << ")(";
          es << join(square_ess().at(e.to_corner_src[1])) << ");(";
          es << join(square_ess().at(e.shared_conc)) << ");(";
          es << join(square_ess().at(e.to_corner_dst[0])) << ")(";
          es << join(square_ess().at(e.to_corner_dst[1])) << ")]\"]";
          return es.str();
        };
      for (const auto& [src, ev] : edges_)
        for (const auto& e : ev)
          {
            bool invert = !square_ess_.at(e.to_corner_dst[1]).empty() && square_ess_.at(e.to_corner_src[1]).empty();
            if (!invert)
            {
              auto srcs = "N" + std::to_string(src);
              auto dsts = "N" + std::to_string(e.dst);
              ss << srcs << " -> " << dsts << fmtedge(e) << ";\n";
            }
            else
            {
              auto ep = edge_t{(int) square_ess_.at(e.to_corner_src[1]).empty(), e.to_corner_dst, e.shared_conc, e.to_corner_src, src};
              auto srcsp = "N" + std::to_string(e.dst);
              auto dstsp = "N" + std::to_string(src);
              ss << srcsp << " -> " << dstsp << fmtedge(ep) << ";\n";
            }
          }
      ss << "}\n";
      return ss.str();
    }

public:
  /**
   * Constructs the max-cell HDA from a regular, explicit HDA
   * Only for full (i.e. not partial) HDAs
   * @param hda The hda to transform
   */

  basic_max_conv_hda(const auto& hda)
    : place2id_(hda.place2id())
    , trans2id_(hda.trans2id())
    , nmarkings_(0)
  {
    nmarkings_ = std::ranges::count_if(hda.X_ext().get_fwd(), [&hda](const auto& p)
      {return hda.square_ess().at(p.first.second).size() == 0; });
    convert_me_(hda);
  }

  /**
 * \addtogroup All getters
 * @{
 */
  const bimap<conc_t, concid_t>& square_ess() const { return square_ess_; }
  const bimap<marking_t, markingid_t>& mark_ess() const { return mark_ess_; }

  const bimap<config_t, cid_t>& X_ext() const { return X_ext_; } // map between configurations and cells (one-to-one)

  // Keep the labeling
  const place2id_t& place2id() const { return place2id_; }
  const trans2id_t& trans2id() const { return trans2id_; }
  /** @} */

  const auto& edges() const {return edges_;}

  unsigned nmarkings() const { return nmarkings_; }

  /**
 * Function dispatching to the different representation functionalities
 * @param opt ["str", "dot"] are currently supported
 * @return The corresponding string
 */
  std::string format(const std::string& opt) const
  {
    if (opt.starts_with("str"))
    {
      auto ss = std::ostringstream();
      ss << *this;
      return ss.str();
    }
    else if (opt.starts_with("dot"))
    {
      return to_dot_(opt.substr(3));
    }
    else
      throw std::invalid_argument("unknown option: " + opt);
  }

  friend std::ostream& operator<<(std::ostream& os, const basic_max_conv_hda& maxhda)
  {
    auto concstream = std::ostringstream();
    auto concstr = std::string();
    auto markstream = std::ostringstream();
    auto markstr = std::string();
    // Cells by increasing conclist and marking
    auto cc = config_t();
    for (const auto& [conc, concid] : maxhda.square_ess_.get_fwd())
    {
      cc.second = concid;
      concstream.str("");
      for (auto tid : conc)
        concstream << maxhda.trans2id_.at(tid) << ',';
      concstr = concstream.str();
      if (!concstr.empty())
        concstr.pop_back();
      concstr = "[" + concstr + "]";

      for (const auto& [mark, markid] : maxhda.mark_ess_.get_fwd())
      {
        cc.first = markid;
        markstream.str("");
        for (auto m : mark)
          markstream << m << ',';
        markstr = markstream.str();
        if (!markstr.empty())
          markstr.pop_back();
        markstr = "[" + markstr + "]";
        if (maxhda.X_ext_.get_fwd().contains(cc))
        {
          // Actual cell
          auto cidx = maxhda.X_ext_.at(cc);
          os << cidx << " : " << markstr << " x " << concstr << '\n';
          os << "edges:\n";
          if (auto it = maxhda.edges().find(cidx); it!=maxhda.edges().end())
          {
            for (const auto& t : it->second)
            {
              os << t.edge_type << ";(";
              os << join(maxhda.square_ess().at(t.to_corner_src[0])) << ")|(";
              os << join(maxhda.square_ess().at(t.to_corner_src[1])) << ");(";
              os << join(maxhda.square_ess().at(t.shared_conc)) << ")|(";
              os << join(maxhda.square_ess().at(t.to_corner_dst[0])) << ")|(";
              os << join(maxhda.square_ess().at(t.to_corner_dst[1])) << ");";
              os << t.dst << '\n';
            }
          }
        }
      }
    }
    return os;
  }

  /**
 * Get basic statistics about the HDA as csv line
 * @param want_header If true, the first line is the header
 */
  std::string get_csv_stats(bool want_header) const
  {
    if (want_header)
      return get_csv_base_header_() + "\n" + get_csv_base_();
    return get_csv_base_();
  }
  /**
   * Parametrizable version of getting statistics in csv form about the HDA.
   * Will append the string obtained from \a gen to the basic statistics
   * @tparam CSVGEN Class of the gen function
   * @param want_header Return header?
   * @param add_header Additional header for the infos generated by gen
   * @param gen The function taking an HDA instance and returning a csv compatible string
   * @return csv compatible string containing the addtional information
   */
  template <class CSVGEN>
  std::string get_csv_stats(bool want_header, const std::string& add_header, const CSVGEN& gen)
  {
    if (want_header)
      return get_csv_base_header_() + add_header + "\n" + get_csv_base_() + gen(*this);
    return get_csv_base_() + gen(*this);
  }

};

class standard_max_conv_hda : public basic_max_conv_hda<pnhda::standard_pn_hda>
{
};

class weighted_max_conv_hda : public basic_max_conv_hda<pnhda::weighted_pn_hda>
{
};

}