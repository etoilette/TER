#pragma once

#pragma once

#include <algorithm>
#include <functional>
#include <iterator>
#include <memory>
#include <optional>
#include <ranges>
#include <set>
#include <unordered_set>
#include <vector>
#include <z3++.h>
#include "pn.h"
#include "pnhda.h"
#include "utils.h"

namespace naivemaxhda
{
  /**
   * Standard definition of a max-HDA derived from a Petri Net
   * @tparam PNHDAT Structures defining all the types
   */
  template <class PNHDAT>
  class basic_max_hda : public PNHDAT
  {
  public:
    constexpr static inline bool ismax = true;
    /**
     * @addtogroup PN definitions
     * @{
     */
    using pnames_t = PNHDAT::pnames_t;
    using pnamesid_t = PNHDAT::pnamesid_t;
    using tnames_t = PNHDAT::tnames_t;
    using tnamesid_t = PNHDAT::tnamesid_t;
    using marking_t = PNHDAT::marking_t;
    using transition_t = PNHDAT::transition_t;
    using all_trans_t = PNHDAT::all_trans_t;
    using place2id_t = PNHDAT::place2id_t;
    using trans2id_t = PNHDAT::trans2id_t;
    /** @} */

    /**
     * @addtogroup PNHDA definitions
     * @{
     */
    using cid_t = PNHDAT::cid_t;
    using conc_t = PNHDAT::conc_t;
    using concid_t = PNHDAT::concid_t;
    using markingid_t = PNHDAT::markingid_t;
    using config_t = PNHDAT::config_t;
    using fm_t = PNHDAT::fm_t;
    /** @} */

  protected:
    // used to get a unique id per cell
    cid_t count = 0; // cell id 0 is invalid, it must be incremented before

    bimap<conc_t, concid_t> square_ess_; /** Bimap for reachable conclists */
    bimap<marking_t, markingid_t> mark_ess_; /** Bimap for reachable markings */

    bimap<config_t, cid_t> X_ext_; /** bimap between configurations and cells (one-to-one) */

    std::map<cid_t, std::set<cid_t>> next_; /** Set of immediately reachable cells from the current one (eg: next max-cells) */

    std::set<cid_t> init_set_; /** The \bot set */
    std::set<cid_t> final_set_; /** The \top set */

    place2id_t place2id_; /** Copied labeling from the underlying petri net for states */
    trans2id_t trans2id_; /** Copied labeling from the underlying petri net for transitions */

    /**
     * Iterator type to iterate over transitions firing within a conclist
     */
    struct transition_iterator {
    private:
      // state type
      struct transition_state {
        tnamesid_t tid;
        uint32_t max_conc;
        uint32_t curr_conc;
      }; // struct transition_state

      // current state
      std::vector<transition_state> data;

      // add a Transition: called by constructor
      void add(tnamesid_t tid, uint32_t max_conc) {
        data.emplace_back(tid, max_conc, 0U);
      }

    public:
      transition_iterator(const conc_t& c) {
        transition_state s(0, 0, 0);
        for (const tnamesid_t& t: c) {
          if ((s.tid != t) and (s.max_conc > 0)) {
            add(s.tid, s.max_conc);
            s.tid = t;
            s.max_conc = 0;
          }
          s.tid = t;
          ++s.max_conc;
        }
        if (s.max_conc > 0) {
          add(s.tid, s.max_conc);
        }
      }

      bool next() {
        auto it = data.begin();
        while (it != data.end()) {
          if (it->curr_conc != it->max_conc) {
            ++it->curr_conc;
            return true;
          }
          it->curr_conc = 0;
          ++it;
        }
        return false;
      }
      // todo create an incremental less costly version
      /**
       * Fire all current concurrent transitions described by the current states
       * @param m The initial marking to start firing the transitions
       * @param pn The Petri Net used to fire the transitions
       * @return the resulting marking
       * @note We assert the given marking can fire all of those transitions...
       */
      marking_t fire(marking_t m, const pn::standard_pn_t& pn)
      {
        auto mp = m;
        fire_me(mp, pn);
        return mp;
      }

      void fire_me(marking_t& m, const pn::standard_pn_t& pn)
      {
        for (const auto& s: data) {
          for (uint32_t i = 0; i < s.curr_conc; i++) {
            pn.try_fire_me(pn.all_trans().find(s.tid)->second, m);
          }
        }
      }

      /**
       * Create the conclist with all current concurrent transitions described by the current states
       * @return the resulting conclist
       */
      conc_t create_conclist() {
        conc_t c;
        for (const auto& s: data) {
          for (uint32_t i = 0; i < s.curr_conc; i++) {
            c.emplace_back(s.tid);
          }
        }
        // Should be the case by default as the tid are pushed back from a sorted vector...
        // std::sort(c.begin(), c.end());
        return c;
      }
    }; // struct transition_iterator

  public:
    /**
     * Construct a max-HDA from the given petri net
     * @param pn The Petri Net to convert
     * @param version Which algorithm version to use
     */
    basic_max_hda(const std::pair<const pn::standard_pn_t* const, int>& p)
      : place2id_(p.first->place2id()),
        trans2id_(p.first->trans2id()) {
      convert_me_(*(p.first), p.second);
    }

    /**
     * Add links from a parent cell to a range of child cell (inherited from a false-max-cell)
     */
    void add_link(cid_t src, cid_t dst)
    {
      if (src > 0) // 0 is initial
        next_[src].insert(dst);
    }

    /**
     * Add a new cell using the config (marking, conc) as a successor of parentId cell
     * Copy conc and marking
     * ¸\pre The cell does not exist before
     */
    cid_t add_cell(cid_t parentId, const marking_t& marking, const conc_t& conc)
    {
      cid_t nc = ++count;
      config_t conf = config_t();
      conf.first = mark_ess_.try_emplace(marking, mark_ess_.size()).first->second;
      conf.second = square_ess_.try_emplace(conc, square_ess_.size()).first->second;

      if (!X_ext_.try_emplace(conf, nc).second)
        throw std::runtime_error("Adding existing cell");

      add_link(parentId, nc);

      return nc;
    }

    /**
     * Declares the cell s to be a subscell of m.
     * No check is performed, s is erased and replaced in links by m
     */
    void declare_subcell(cid_t s, cid_t m)
    {

      // Replace targets
      for (auto src_it = next_.begin(); src_it != next_.end(); ++src_it)
        if (src_it->second.erase(s))
          src_it->second.insert(m);
      // add next of s to m
      for (const auto& d : next_[s])
        next_[m].insert(d);
      // Remove s
      next_.erase(s);
      // Also erase the cell
      X_ext_.erase_value(s);
      // At some point we should clean up square_ess_ and mark_ess_
      // as they might be polluted with "dead" objects
    }

    /**
     * Remove all unused markings and conclists
     */
    void clean_me_()
    {
      std::unordered_set<markingid_t> unused_mark_;
      std::unordered_set<concid_t> unused_conc_;

      for (const auto& p : mark_ess().get_fwd())
        unused_mark_.insert(p.second);
      for (const auto& p : square_ess().get_fwd())
        unused_conc_.insert(p.second);

      for (const auto& p : X_ext_.get_fwd())
      {
        unused_mark_.erase(p.first.first);
        unused_conc_.erase(p.first.second);
      }

      for (markingid_t m : unused_mark_)
        mark_ess_.erase_value(m);
      for (concid_t c : unused_conc_)
        square_ess_.erase_value(c);

      //for (auto it = mark_ess_.get_fwd().begin(); it != mark_ess_.get_fwd().end(); ++it)
      //  if (!used_mark_.count(it->second))
      //    mark_ess_.erase_value(it->second);
      //for (auto it = square_ess_.get_fwd().begin(); it != square_ess_.get_fwd().end(); ++it)
      //  if (!used_conc_.count(it->second))
      //    square_ess_.erase_value(it->second);
    }

    /**
     * Validate a max conc step (verify it is max)
     * We assume m is already a copy
     */
    static
    bool validate_max_conc_step(const pn::standard_pn_t& pn, marking_t& m)
    {
      for (tnamesid_t t = 0; t < pn.all_trans().size(); t++) {
        if (pn::standard_pn_t::try_start_me(pn.all_trans().find(t)->second, m)) {
          return false;
        }
      }
      return true;
    }

    /**
     * Generate all of the max conc step for the current marking
     *
     * Pass marking and current by copy
     * TODO: The current version use a brute force recursive algorithm. It must be optimised!
     */
    static
    std::vector<conc_t>
    gen_max_conc_step(const pn::standard_pn_t& pn, marking_t marking,
                      conc_t current = {}, tnamesid_t id = 0)
    {
      if (id == pn.all_trans().size()) {
        if (validate_max_conc_step(pn, marking)) {
          return { current };
        }
        return { };
      }

      std::vector<conc_t> out = std::vector<conc_t>();

      for (auto& c : gen_max_conc_step(pn, marking, current, id + 1)) {
        if (not c.empty()) {
          out.push_back(c);
        }
      }

      while (pn::standard_pn_t::try_start_me(pn.all_trans().find(id)->second, marking)) {
        current.emplace_back(id);

        for (auto& c : gen_max_conc_step(pn, marking, current, id + 1)) {
          if (not c.empty()) {
            out.push_back(c);
          }
        }
      }
      return out;
    }

    std::string format(const std::string& opt) const
    {
      if (opt.starts_with("str"))
      {
        auto ss = std::ostringstream();
        ss << *this;
        return ss.str();
      }
      if (opt.starts_with("dot"))
      {
        return to_dot_(opt.substr(3));
      }
      throw std::invalid_argument("unknown option: " + opt);
    }

    std::string get_csv_base_header_() const { return "ncells,nconclist,nmarkings,nlinks,maxdim,L(dim;ncell)"; }
    std::string get_csv_base_() const
    {
      auto ss = std::ostringstream();

      ss << X_ext_.size() << ',' << square_ess_.size() << ',' << mark_ess_.size() << ',' << next_.size() << ',';

      // dimension of cell and how many
      auto celldet = std::map<uint, uint>();

      for (const auto& [conf, _] : X_ext_.get_fwd())
      {
        concid_t cid = conf.second;
        const auto& concl = square_ess_.at(cid);
        uint dim = concl.size();

        auto it = celldet.try_emplace(dim, 0).first;
        ++(it->second);
      }

      ss << celldet.rbegin()->first << ",[";

      for (const auto& [dim, ndim] : celldet)
        ss << '[' << dim << ';' << ndim << "];";
      ss << "]";

      return ss.str();
    }

    std::string get_csv_stats(bool want_header) const
    {
      if (want_header)
        return get_csv_base_header_() + "\n" + get_csv_base_();
      return get_csv_base_();
    }

    template <class CSVGEN>
    std::string get_csv_stats(bool want_header, const std::string& add_header, const CSVGEN& gen)
    {
      if (want_header)
        return get_csv_base_header_() + add_header + "\n" + get_csv_base_() + gen(*this);
      return get_csv_base_() + gen(*this);
    }

    friend std::ostream& operator<<(std::ostream& os, const basic_max_hda& maxhda)
    {
      auto cc = config_t();
      for (const auto& [conc, concid] : maxhda.square_ess_.get_fwd())
      {
        cc.second = concid;
        auto concstr = "[" + join((conc | std::views::transform([&](const auto& tid) -> tnames_t{return maxhda.trans2id().at(tid);}))) + "]";

        for (const auto& [mark, markid] : maxhda.mark_ess_.get_fwd())
        {
          cc.first = markid;
          auto markstr = "[" + join(mark) + "]";
          if (maxhda.X_ext_.get_fwd().contains(cc))
          {
            auto cidx = maxhda.X_ext_.at(cc);
            os << cidx << " : " << markstr << " x " << concstr << '\n';
            os << "links:\n";
            if (auto it_linkset = maxhda.next().find(cidx); it_linkset != maxhda.next().end())
            {
              for (const auto& dstidx : it_linkset->second)
                os << cidx << " -> " << dstidx << '\n';
            }
          }
        }
      }
      return os;
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

      std::string rnodes = "rank=source;\n";
      std::string rcnodes = "rank=source;\n";

      auto ss = std::ostringstream();

      ss << "digraph automaton {\n"
         << "rankdir=" << dot_dir << ";  // Ensure top-to-bottom layout\n"
         << "// Subgraph for node identification\n";

      // Already define all nodes
      // Put all the actual states
      ss << "node [shape = circle];\n"
         << "nodesep=0.15; // Reverting back to standard\n"
         << "ranksep=0.35; // Reverting back to standard\n"
         << rnodes;


      for (const auto& p : X_ext_.get_rev())
      {
        auto ns = "N" + std::to_string(p.first);
        ss << ns << " [label = \"" << ns << "\"];\n";
      }
      // Continue with the subgraph
      ss << "subgraph cluster_comments {\n"
         << rcnodes << "node [shape=none];\n"
         << "edge [style=invis];  // Make edges invisible\n"
         << "nodesep=0.15; // Less horizontal distance\n"
         << "ranksep=0.35; // Less vertical distance\n";
      // Create all the nodes
      // todo factorize me
      auto fmtnode = [&](cid_t cid)
      {
        auto [markid, concid] = X_ext_.at(cid);
        const auto& mark = mark_ess_.at(markid);
        const auto& conc = square_ess_.at(concid);

        return "[" + join(mark) + "]" + "x" + "["
               + join(conc | std::views::transform(
                 [&](const auto& tid){return trans2id().at(tid);})) + "]";
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

      ss << "} // End subgraph\n";
      // And all the transitions
      for (const auto& [src, dstset] : next_)
        for (const auto& dst : dstset)
          ss << "N" << std::to_string(src) << " -> N" << std::to_string(dst) << ";\n";

      ss << "}\n";
      return ss.str();
    }

    /**
     * \addtogroup All getters
     * @{
     */
    const bimap<conc_t, concid_t>& square_ess() const { return square_ess_; }
    const bimap<marking_t, markingid_t>& mark_ess() const { return mark_ess_; }

    const bimap<config_t, cid_t>& X_ext() const { return X_ext_; } // map between configurations and cells (one-to-one)

    const std::map<cid_t, std::set<cid_t>>& next() const { return next_; } // Set of immediately reachable cells from a given cell-id (Partially filled during the max-HDA construction)

    const std::set<cid_t>& init_set() const { return init_set_; } // The \bot set
    const std::set<cid_t>& final_set() const { return final_set_; } // The \top set

    // Keep the labeling from the PN
    const place2id_t& place2id() const { return place2id_; }
    const trans2id_t& trans2id() const { return trans2id_; }
    /** @} */

  private:
    void convert_me_(const pn::standard_pn_t& pn, int version)
    {
      switch (version)
      {
      case 1: return convert_me_1_(pn);
      case 2: return convert_me_2_(pn);
      default: throw std::runtime_error("unknown version");
      }
    }

    // Definition of additional structures
    struct cell_iterator_t
    {
      marking_t m;
      cid_t parent;
      static bool
      cmp(const cell_iterator_t& l, const cell_iterator_t& r)
      {
        return l.m < r.m;
      }
    };

    struct opt_context_t
    {
      struct transvar_t
      {
        z3::expr texpr;
        int max_mult;
        int eff_mult;
      };
      struct mark_expr_t
      {
        z3::expr expr;
        int diff;
      };
      struct conc_context_t
      {
        z3::context ctx; // Context is not a regular variable
        std::map<tnamesid_t, transvar_t> transvars; // trans -> variable and
        // (max multiplicity, remaining multiplicity)
        std::map<unsigned, mark_expr_t> mark_expr; // affected place -> expr; sorted by place
      };
      std::map<concid_t, conc_context_t> conc2ctx; // Mapping each conclist to its context
    };
    // Implementation of subcell implementation
    bool is_subcell_impl_(const marking_t& m, const conc_t& c,
                          const marking_t& mp, const conc_t& cp,
                          opt_context_t& opt_context, const pn::standard_pn_t& pn )
    {
      using transvar_t = opt_context_t::transvar_t;
      using mark_expr_t = opt_context_t::mark_expr_t;
      // using conc_context_t = opt_context_t::conc_context_t;
      // std::cerr << "Testing [" << join(m) << "]x[" << join(c)
      //           <<  "] included in [" << join(mp) << "]x[" << join(cp)
      //           << "]\n";

      // Find the context
      auto [it, ins] = opt_context.conc2ctx.try_emplace(square_ess_.at(cp));
      if (ins)
      {
        // We need to create the context
        auto& c_ctx = it->second;
        // c_ctx.ctx = z3::context(); No assignment for context, default cstr is ok

        // Create all transition variables
        for (auto it_tr = cp.begin(); it_tr != cp.end(); ++it_tr)
        {
          // Create the var
          auto [itvar, insvar] = c_ctx.transvars.try_emplace(
            *it_tr,
            transvar_t(c_ctx.ctx.int_const(trans2id_[*it_tr].c_str()), 1, 0)
            );
          if (!insvar)
            throw std::runtime_error("Can not yet exist");
          const auto& texpr = itvar->second.texpr; // For later
          // Count the multiplicty
          auto it_tr2 = it_tr;
          ++it_tr2;
          while (it_tr2 != cp.end() && *it_tr2 == *it_tr)
          {
            ++itvar->second.max_mult;
            ++it_tr;
            ++it_tr2;
          }
          // Modify the expression for pre and postplaces
          // Loop over pre and post places:
          const auto& t = pn.all_trans().at(*it_tr);
          for (const auto& p : std::get<0>(t)) // places of preplaceset
          {
            // todo Improvement : Ignore if in pre AND post
            auto [itmexpr, ins] = c_ctx.mark_expr.try_emplace(p,
                mark_expr_t{-texpr, 0});
            if (!ins)
              itmexpr->second.expr = itmexpr->second.expr - texpr;
          }
          for (const auto& p : std::get<1>(t)) // places of postplaceset
          {
            // todo Improvement : Ignore if in pre AND post
            auto [itmexpr, ins] = c_ctx.mark_expr.try_emplace(p,
                mark_expr_t{texpr, 0});
            if (!ins)
              itmexpr->second.expr = itmexpr->second.expr + texpr;
          }
        }
      } // Done with setup
      // Return if same
      if (m == mp && c == cp)
        //std::cerr << "Selfcomp\n";
        return true;

      // If they differ in an unaffected place -> false
      // Walk through the map at the same time as the marking
      auto& ctx = it->second;
      const unsigned N = m.size();
      auto itmpla = ctx.mark_expr.begin();
      for (unsigned idx = 0; idx < N; ++idx)
      {
        if (m[idx] != mp[idx])
        {
          // Advance the map iterator
          while (itmpla != ctx.mark_expr.end()
                 && itmpla->first < idx)
          {
            itmpla->second.diff = 0; // All prior plaes had no diff
            ++itmpla;
          }
          if (itmpla == ctx.mark_expr.end()
              || itmpla->first != idx)
            return false;
          // Remember the diff between the markings for later
          itmpla->second.diff = (int) m[idx] - (int) mp[idx];
          // Next place
          ++itmpla;
        }
      }
      // set remaining to zero as well
      for (; itmpla != ctx.mark_expr.end(); ++itmpla)
        itmpla->second.diff = 0;
      // Second precheck: compute the remaining conclist
      // cp \setminus c
      // compute the effective multiplicity
      auto itcp = ctx.transvars.begin();
      for (auto itc = c.begin(); itc != c.end(); ++itc)
      {
        int tmultc = 1;
        auto itc2 = itc;
        ++itc2;
        while (itc2 != c.end() && *itc2 == *itc)
        {
          ++tmultc;
          ++itc;
          ++itc2;
        }
        // Compare to availability in cp
        //auto itcp = ctx.transvars.find(*itc);
        //if (itcp == ctx.transvars.end() || itcp->second.max_mult < tmultc)
        //  return false;
        // Search for correct entry
        while (itcp != ctx.transvars.end()
               && itcp->first < *itc)
        {
          itcp->second.eff_mult = itcp->second.max_mult;
          ++itcp;
        }
        // Check if found
        if (itcp == ctx.transvars.end() // Not found end
            || itcp->first != *itc  // Not found other
            || itcp->second.max_mult < tmultc) // not enough
          return false;
        // All good ->
        // Update available multiplicity of events
        itcp->second.eff_mult = itcp->second.max_mult - tmultc;
        // Used and advance
        ++itcp;
      }
      // Set rest
      for (; itcp != ctx.transvars.end(); ++itcp)
        itcp->second.eff_mult = itcp->second.max_mult;
      // Minimal condition passed
      // We need to actually solve the optimization prob
      auto s = z3::solver(ctx.ctx);
      // Limit the multiplicity
      for (const auto& p : ctx.transvars)
      {
        s.add(0 <= p.second.texpr);
        s.add(p.second.texpr <= p.second.eff_mult);
      }
      // the firing must correpond  to the diff
      for (const auto& p : ctx.mark_expr)
        s.add(p.second.expr == p.second.diff);

      //std::cerr << s.to_smt2() << std::endl;
      switch (s.check())
      {
      case z3::unsat: return false;
      case z3::sat: return true;
      default: throw std::logic_error("Must be SAT or UNSAT"); break;
      }
    }; // todo seperate precheck, heuristics and solver construction


    void convert_me_1_(const pn::standard_pn_t& pn)
    {
      /**
       * Check if (m,c) is a subcell of (mp,cp)
       */
      auto is_subcell_
        = [&, opt_context = opt_context_t()](const marking_t& m, const conc_t& c,
                                          const marking_t& mp, const conc_t& cp) mutable
        {
          return is_subcell_impl_(m, c, mp, cp, opt_context, pn);
        };


      auto stack = std::vector<cell_iterator_t>();

      auto push = [&](auto&& m, cid_t p)
        {
          //std::cerr << "push " << join(m) << "; " << p << std::endl;
          stack.emplace_back(std::forward<decltype(m)>(m), p);
        };
      auto pop = [&]()
        {
          auto r = std::move(stack.back());
          //std::cerr << "poph " << join(r.m) << "; " << r.parent << std::endl;
          stack.pop_back();
          return r;
        };

      push(pn.init_marking(), count);

      // Main alg
      auto subpairs = std::vector<std::pair<cid_t, cid_t>>();
      while (!stack.empty())
      {
        auto [nm, pid] = pop();

        auto all_MCS = gen_max_conc_step(pn, nm);
        for (const auto& amcs : all_MCS)
        {
          // Check if subcell
          if (std::ranges::any_of(X_ext_.get_fwd(),
              [&](const auto& p)
              {
                // (nm, amcs) is not yet registered -> No self-equiv possible
                return is_subcell_(nm, amcs, mark_ess_.at(p.first.first),
                                   square_ess_.at(p.first.second));
              }))
            continue; // Skip this cell

          // Generate new cell an link
          auto nc = add_cell(pid, nm, amcs);

          // Check if any existing cells are covered
          subpairs.clear();
          for (const auto& p : X_ext_.get_fwd())
          {
            if (p.second == nc)
              continue; // itself
            if (is_subcell_(mark_ess_.at(p.first.first), square_ess_.at(p.first.second),
                            nm, amcs))
              subpairs.emplace_back(p.second, nc);
          }
          for (const auto& [s,m] : subpairs)
              declare_subcell(s, m);


          // Push all corner markings
          auto it = transition_iterator(amcs);
          while (it.next()) // First is first no trans
            push(it.fire(nm, pn), nc);
        }
      }
      clean_me_(); // Remove temp markings and conclists
    }


    void convert_me_2_(const pn::standard_pn_t& pn)
    {
      /**
       * Check if (m,c) is a subcell of (mp,cp)
       */
      auto is_subcell_
        = [&, opt_context = opt_context_t()](const marking_t& m, const conc_t& c,
                                             const marking_t& mp, const conc_t& cp) mutable
        {
          return is_subcell_impl_(m, c, mp, cp, opt_context, pn);
        };

      /**
       * Check if the marking m can be generated by (mp, cp)
       */
      auto is_generated_ =
        [&] (const marking_t& m,
             const marking_t& mp, const conc_t& cp)
          {
            return is_subcell_(m, conc_t(), mp, cp);
          };

      auto stack = std::set<cell_iterator_t, decltype(cell_iterator_t::cmp)*>(&cell_iterator_t::cmp);

      auto try_push = [&](auto&& m, cid_t p) -> auto
        {
          //std::cerr << "push " << join(m) << "; " << p << std::endl;
          return stack.emplace(std::forward<decltype(m)>(m), p);
        };
      auto try_push_c = [&](auto&& c) -> auto
        {
          //std::cerr << "push " << join(m) << "; " << p << std::endl;
          return stack.emplace(std::forward<decltype(c)>(c));
        };
      auto pop = [&]()
        {
          auto it = stack.begin();
          auto r = std::move(*it);
          stack.erase(it);
          return r;
        };

      try_push(pn.init_marking(), count);

      // Main alg
      auto subpairs = std::vector<std::pair<cid_t, cid_t>>();
      while (!stack.empty())
      {
        auto [nm, pid] = pop();

        auto all_MCS = gen_max_conc_step(pn, nm);
        for (const auto& amcs : all_MCS)
        {
          // Check if subcell
          if (std::ranges::any_of(X_ext_.get_fwd(),
              [&](const auto& p)
              {
                // (nm, amcs) is not yet registered -> No self-equiv possible
                return is_subcell_(nm, amcs, mark_ess_.at(p.first.first),
                                   square_ess_.at(p.first.second));
              }))
            continue; // Skip this cell

          // Generate new cell an link
          auto nc = add_cell(pid, nm, amcs);

          // Check if any existing cells are covered
          subpairs.clear();
          for (const auto& p : X_ext_.get_fwd())
          {
            if (p.second == nc)
              continue; // itself
            if (is_subcell_(mark_ess_.at(p.first.first), square_ess_.at(p.first.second),
                            nm, amcs))
              subpairs.emplace_back(p.second, nc);
          }
          for (const auto& [s,m] : subpairs)
              declare_subcell(s, m);


          // Push all corner markings
          auto it = transition_iterator(amcs);
          while (it.next()) // First is first no trans
          {
            auto next_cell_it = cell_iterator_t(it.fire(nm, pn), nc);
            if (stack.count(next_cell_it) == 0)
            {
              // Only use if not generated
              //try_push_c(std::move(next_cell_it));

              if (std::ranges::none_of(X_ext_.get_fwd(),
                [&](const auto& p)
                  {
                    return p.second != nc
                           && is_generated_(next_cell_it.m,
                                        mark_ess_.at(p.first.first), square_ess_.at(p.first.second));
                  }))
                  try_push_c(std::move(next_cell_it));
            }
          }
        }
      }
      clean_me_(); // Remove temp markings and conclists
    }

  };

  /**
   * @addtogroup maxHDA instances
   * Define all (partial) HDA types by suitable parametrization
   * @{
   */
  class standard_pn_maxhda : public basic_max_hda<pnhda::basic_pnhda_types<pn::basic_pn_types>>
  {
    using basic_max_hda<pnhda::basic_pnhda_types<pn::basic_pn_types>>::basic_max_hda;
  };
}