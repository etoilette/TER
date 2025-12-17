#pragma once
#pragma once

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <iostream>
#include <set>
#include <sstream>
#include <vector>
#include <ranges>
#include <numeric>

#include "pn.h"


namespace pnhda
{

  /**
   * Struct containing the invariant HDA types
   * @tparam BPNT They are dervied from the underlying petri net
   */
  template <class BPNT>
  struct basic_pnhda_types : public BPNT
  {
    using tnamesid_t = BPNT::tnamesid_t;

    using cid_t = uint;

    using conc_t = std::vector<tnamesid_t>; /** tnames will be used for events;
    The event order is represented by keeping this vector sorted. */
    using concid_t = uint; /** Identifier for conclists*/

    using markingid_t = uint; /** Identifier for marking */
    using config_t = std::pair<markingid_t, concid_t>; /** A configuration is a pair
    (marking, conclist) */

    using fm_t = std::pair<tnamesid_t, cid_t>; /** Helper for face-maps*/

    /**
     *Dimension of a cell / conclist
     * -> Should have a size method
     */
    //uint get_dim(const conc_t& conc) const noexcept { return conc.size(); }
  };

  /**
   * Standard definition of a HDA derived from a Petri Net
   * @tparam PNHDAT Structures defining all the types
   */
  template <class PNHDAT>
  class basic_pn_hda : public PNHDAT
  {
  public:
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
    //using get_dim = PNHDAT::get_dim;
    /** @} */

  protected:
    bimap<conc_t, concid_t> square_ess_; /** Bimap for reachable conclists */
    bimap<marking_t, markingid_t> mark_ess_; /** Bimap for reachable markings */

    bimap<config_t, cid_t> X_ext_; /** bimap between configurations and cells (one-to-one) */

    // These are sets due to symmetry
    std::map<cid_t, std::set<fm_t>> delta_0_; /** Lower face map of a cell (All "unstartings") */
    std::map<cid_t, std::set<fm_t>> delta_1_; /** Upper face map of a cell (All "terminations") */

    std::set<cid_t> init_set_; /** The \bot set */
    std::set<cid_t> final_set_; /** The \top set */

    place2id_t place2id_; /** Copied labeling from the underlying petri net for states */
    trans2id_t trans2id_; /** Copied labeling from the underlying petri net for transitions */


    /**
     * Transforms the given petri net into a (p)HDA
     */
    template <class PN>
    void convert_me_(const PN& pn)
    {
      // Construct initial cells, then recurse

      // Basic version, each instance has its own marking and
      // conclist -> Use some on the fly construction
      struct estate_iterator
      {
        marking_t m_; // Current marking
        conc_t c_; // Current conclist
        cid_t cid_;

        all_trans_t::const_iterator it_start_; // Next "starting" transition to explore
        all_trans_t::const_iterator it_start_end_; // Nahh
        // The vector moves -> can not use iterators
        std::size_t idx_term_; // Next "terminating" transition to explore
        std::size_t idx_term_end_; // Nahh

        const PN* spn_ = nullptr;


        operator bool() const { return it_start_ == it_start_end_ && idx_term_ == idx_term_end_; }

        /// Returns the next reachable configuration
        /// @return Return std::nullopt if the iterator is exhausted,
        /// otherwise it returns wether a transition has been started (0)
        /// or terminated (1), which transition is concerned and
        /// the iterator over this new state
        std::optional<std::tuple<int, tnamesid_t, estate_iterator>> next()
        {
          if (!!*this)
            return std::nullopt;

          auto res = std::optional<std::tuple<int, tnamesid_t, estate_iterator>>();
          res.emplace();
          auto& res_es = std::get<2>(*res);
          res_es.spn_ = spn_;
          res_es.m_ = m_;
          res_es.c_ = c_;

          while (it_start_ != it_start_end_)
          {
            // Try if current transition is fireable
            // else increment and try again
            bool succ = PN::try_start_me(it_start_->second, res_es.m_);
            auto se_idx = it_start_->first;
            ++it_start_;
            if (succ)
            {
              // Update conclist; Add started event
              res_es.c_.push_back(se_idx);
              // Sort
              // todo, do this better
              std::sort(res_es.c_.begin(), res_es.c_.end());
              // Update the iterators

              // starting transitions and terminating transition
              // go from front to back
              res_es.it_start_ = spn_->all_trans().begin();
              res_es.it_start_end_ = spn_->all_trans().end();
              res_es.idx_term_ = 0;
              res_es.idx_term_end_ = res_es.c_.size();
              std::get<0>(*res) = 0;
              std::get<1>(*res) = se_idx;
              return res;
            }
          }
          // We exhausted up transitions
          // So now we will terminate
          // Quick-fix a-priori semantics
          // (Instead of (super-a-priori semantics)
          // A transition may only terminate if it does not
          // impact the inhibitor places of an transition
          // that is currently executing
          if constexpr (std::derived_from<PNHDAT, pn::weighted_pni_types>)
          {
            while (true)
            {
              if (idx_term_ == idx_term_end_)
                return std::nullopt;

              // Here comes the hotfix
              // This is _very_ inefficient todo
              // Loop through all transitions that still exist
              // and check if their inhibitor places are now populated
              const unsigned Nt = c_.size();
              bool okT = true;
              for (uint i = 0; i < Nt; ++i)
              {
                if (i == idx_term_)
                  continue; // Disregard the one to be terminated
                const auto& new_marks = spn_->all_trans().at(c_[idx_term_]).second;
                // newMarks == weighted_place_set_t
                auto it_n = new_marks.begin();
                auto i_inhib = spn_->all_trans().at(c_[i]).first.first;
                // i_inhib == place_set_t
                auto it_i = i_inhib.begin();

                // Loop through both
                while (it_i != i_inhib.end() && it_n != new_marks.end())
                {
                  if (*it_i == it_n->first && it_n->second) // todo zero weight allowed??
                  {
                    okT = false;
                    break;
                  }
                  if (*it_i < it_n->first)
                    ++it_i;
                  else
                    ++it_n;
                }
                if (!okT)
                  break;
              }
              if (!okT)
              {
                // Try next
                ++idx_term_;
                continue;
              }

              // Terminate the current and advance
              // Update marking
              PN::terminate_me(spn_->all_trans().at(c_[idx_term_]), res_es.m_);
              // Remove from conclist
              res_es.c_.erase(res_es.c_.begin() + idx_term_);

              // Update the rest of the res
              res_es.it_start_ = spn_->all_trans().begin();
              res_es.it_start_end_ = spn_->all_trans().end();
              res_es.idx_term_ = 0;
              res_es.idx_term_end_ = res_es.c_.size();
              std::get<0>(*res) = 1;
              std::get<1>(*res) = c_[idx_term_];
              // Increment the current once everything is done
              ++idx_term_;
              return res;
            }
          }
          else
          {
            if (idx_term_ == idx_term_end_)
              return std::nullopt;
            // Terminate the current and advance
            // Update marking
            PN::terminate_me(spn_->all_trans().at(c_[idx_term_]), res_es.m_);
            // Remove from conclist
            res_es.c_.erase(res_es.c_.begin() + idx_term_);
            // Update the rest of the res
            res_es.it_start_ = spn_->all_trans().begin();
            res_es.it_start_end_ = spn_->all_trans().end();
            res_es.idx_term_ = 0;
            res_es.idx_term_end_ = res_es.c_.size();
            std::get<0>(*res) = 1;
            std::get<1>(*res) = c_[idx_term_];
            // Increment the current once everything is done
            ++idx_term_;
            return res;
          }
        }
      }; // estate iterator

      auto stack = std::vector<estate_iterator>();

      auto id_one = [&](estate_iterator& it)
      {
        auto cc = config_t();
        cc.first = mark_ess_.try_emplace(it.m_, mark_ess_.size()).first->second;
        cc.second = square_ess_.try_emplace(it.c_, square_ess_.size()).first->second;

        auto [itX, ins] = X_ext_.try_emplace(cc, X_ext_.size());
        it.cid_ = itX->second;
        return ins;
      };

      {
        // Push the first
        stack.push_back(estate_iterator());
        stack.back().spn_ = &pn;
        stack.back().m_ = pn.init_marking();
        stack.back().c_ = conc_t();
        stack.back().it_start_ = pn.all_trans().begin();
        stack.back().it_start_end_ = pn.all_trans().end();
        stack.back().idx_term_ = 0;
        stack.back().idx_term_end_ = stack.back().c_.size();
        id_one(stack.back());
        // Set it as initial
        init_set_.insert(stack.back().cid_);
      }

      while (!stack.empty())
      {

        auto ne = stack.back().next();
        if (!ne)
        {
          stack.pop_back();
          continue;
        }
        // Compute the corresponding cell
        auto& [dir, tid, eit] = *ne;
        bool is_new = id_one(eit);
        // Set the face maps
        if (dir)
          delta_1_[stack.back().cid_].emplace(tid, eit.cid_);
        else
          delta_0_[eit.cid_].emplace(tid, stack.back().cid_);

        if (is_new)
          stack.push_back(std::move(eit));
      }
    }

    std::string get_csv_base_header_() const { return "ncells,nconclist,nmarkings,nd0,nd1,maxdim,L(dim;ncell)"; }

    std::string get_csv_base_() const
    {
      auto ss = std::ostringstream();

      ss << X_ext_.size() << ',' << square_ess_.size() << ',' << mark_ess_.size() << ','
      << (std::accumulate(delta_0_.begin(), delta_0_.end(), (size_t)0,
        [](auto lhs, const auto& ev)
          {return lhs + ev.second.size(); })) << ','
      << (std::accumulate(delta_1_.begin(), delta_1_.end(), (size_t)0,
        [](auto lhs, const auto& ev)
          {return lhs + ev.second.size(); })) << ',';

      // dimension of cell and how many
      auto celldet = std::map<uint, uint>();

      for (const auto& [conf, _] : X_ext_.get_fwd())
      {
        concid_t cid = conf.second;
        const auto& concl = square_ess_.at(cid);
        uint dim = concl.size(); //PNHDAT::get_dim(concl);

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

      // Dummy for initial node
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

      ss << "} // End subgraph\n";
      // And all the transitions
      for (const auto& p : delta_0_)
      {
        // "Unstarting" transition
        // key cell is tgt
        auto tgt = "N" + std::to_string(p.first);
        for (const auto& fm : p.second)
        {
          const auto& lbl = trans2id_.at(fm.first);
          auto src = "N" + std::to_string(fm.second);
          ss << src << " -> " << tgt << " [label = \"" << lbl << "+" << "\"];\n";
        }
      }
      for (const auto& p : delta_1_)
      {
        // "Terminating" transition
        // key cell is src
        auto src = "N" + std::to_string(p.first);
        for (const auto& fm : p.second)
        {
          const auto& lbl = trans2id_.at(fm.first);
          auto tgt = "N" + std::to_string(fm.second);
          ss << src << " -> " << tgt << " [label = \"" << lbl << "-" << "\"];\n";
        }
      }

      ss << "}\n";
      return ss.str();
    }

  public:
    /**
     * Constructs a partial HDA from the given petri net \a spn
     * @param spn The net to convert
     */
    basic_pn_hda(const auto& spn) : place2id_(spn.place2id()), trans2id_(spn.trans2id()) { convert_me_(spn); }


    /**
     * \addtogroup All getters
     * @{
     */
    const bimap<conc_t, concid_t>& square_ess() const { return square_ess_; }
    const bimap<marking_t, markingid_t>& mark_ess() const { return mark_ess_; }

    const bimap<config_t, cid_t>& X_ext() const { return X_ext_; } // map between configurations and cells (one-to-one)

    const std::map<cid_t, std::set<fm_t>>& delta_0() const { return delta_0_; } // All "unstartings" from cell cid_t
    const std::map<cid_t, std::set<fm_t>>& delta_1() const { return delta_1_; } // All terminations from cell cid_t

    const std::set<cid_t>& init_set() const { return init_set_; } // The \bot set
    const std::set<cid_t>& final_set() const { return final_set_; } // The \top set

    // Keep the labeling from the PN
    const place2id_t& place2id() const { return place2id_; }
    const trans2id_t& trans2id() const { return trans2id_; }
    /** @} */

    /**
     * Get basic statistics about the HDA as csv line
     * @param want_header If true, the first line is the header
     */
    std::string get_csv_stats(bool want_header) const
    {
      if (want_header)
        return get_csv_base_header_() + "\n" + get_csv_base_();
      else
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
      else
        return get_csv_base_() + gen(*this);
    }

    /**
     * Function dispatching to the different representation functionalities
     * @param opt ["str", "dot"] are currently supported
     * @return The corresponding string
     */
    std::string format(const std::string& opt) const
    {
      if (opt.starts_with("str_max"))
      {

        auto all_max_cells = get_max_cells();

        auto clist2cset = [](const conc_t& c)
          {
            std::map<tnamesid_t, size_t> cset; // trans Id / Multiplicity

            for (const auto& t: c)
              if (!cset.contains(t))
                cset.emplace_hint(cset.end(), t, std::count(c.cbegin(), c.cend(), t));
            return cset;
          };

        auto ss = std::ostringstream();

        ss << "NumMaxCells: " << all_max_cells.size() << '\n';
        ss << "Transitions: ";
        for (const auto& [t_name, t_id] : trans2id_.get_fwd())
          ss << '(' << t_id << ",\"" << t_name << "\"),";
        ss << '\n';
        ss << "<Cells>\n";
        for (auto cidx : all_max_cells)
        {
          const auto& cf = X_ext().at(cidx);
          const auto& mark = mark_ess().at(cf.first);
          const auto& conc = square_ess().at(cf.second);
          auto cset = clist2cset(conc);

          ss << cidx << ": [";
          for (auto m : mark)
            ss << m << ',';
          ss << "] x [";
          for (const auto& [t_id, mult] : cset)
            ss << '(' << t_id << ", " << mult << "),";
          ss << "]\n";
        }

        ss << "</Cells>\n";

        return ss.str();
      }
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

    friend std::ostream& operator<<(std::ostream& os, const basic_pn_hda& spnhda)
    {
      auto concstream = std::ostringstream();
      auto concstr = std::string();
      auto markstream = std::ostringstream();
      auto markstr = std::string();
      // Cells by increasing conclist and marking
      auto cc = config_t();
      for (const auto& [conc, concid] : spnhda.square_ess_.get_fwd())
      {
        cc.second = concid;
        concstream.str("");
        for (auto tid : conc)
          concstream << spnhda.trans2id_.at(tid) << ',';
        concstr = concstream.str();
        if (!concstr.empty())
          concstr.pop_back();
        concstr = "[" + concstr + "]";

        for (const auto& [mark, markid] : spnhda.mark_ess_.get_fwd())
        {
          cc.first = markid;
          markstream.str("");
          for (auto m : mark)
            markstream << m << ',';
          markstr = markstream.str();
          if (!markstr.empty())
            markstr.pop_back();
          markstr = "[" + markstr + "]";
          if (spnhda.X_ext_.get_fwd().contains(cc))
          {
            // Actual cell
            auto cidx = spnhda.X_ext_.at(cc);
            os << cidx << " : " << markstr << " x " << concstr << '\n';
            os << "d0:\n";
            if (auto it_d0 = spnhda.delta_0_.find(cidx); it_d0 != spnhda.delta_0_.end())
            {
              for (const auto& [tid, tgtcid] : it_d0->second)
                os << '(' << spnhda.trans2id_.at(tid) << ", " << tgtcid << ')' << ';';
              os << '\n';
            }
            os << "d1:\n";
            if (auto it_d1 = spnhda.delta_1_.find(cidx); it_d1 != spnhda.delta_1_.end())
            {
              for (const auto& [tid, tgtcid] : it_d1->second)
                os << '(' << spnhda.trans2id_.at(tid) << ", " << tgtcid << ')' << ';';
              os << '\n';
            }
          }
        }
      }
      return os;
    }

    auto get_max_cells() const
    {
      // Take all cells, remove all cells that are the lower or upper face of some other
      auto all_cells = std::vector<bool>(X_ext().get_fwd().size(), true);
      for (const auto& t : delta_0())
        for (const auto& tt : t.second)
          all_cells[tt.second] = false;
      for (const auto& t : delta_1())
        for (const auto& tt : t.second)
          all_cells[tt.second] = false;

      auto r = std::vector<cid_t>();
      for (cid_t i = 0; i < (cid_t) all_cells.size(); ++i)
        if (all_cells[i])
          r.push_back(i);
      return r;
    }

  }; // basic_pn_hda

  /**
   * @addtogroup pHDA instances
   * Define all (partial) HDA types by suitable parametrization
   * @{
   */
  class standard_pn_hda : public basic_pn_hda<basic_pnhda_types<pn::basic_pn_types>>
  {
  };
  class weighted_pn_hda : public basic_pn_hda<basic_pnhda_types<pn::weighted_pn_types>>
  {
  };
  class weighted_pni_hda : public basic_pn_hda<basic_pnhda_types<pn::weighted_pni_types>>
  {
  };
  /** @} */

} // namespace pnhda
