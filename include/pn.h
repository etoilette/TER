#pragma once

#include <algorithm>
#include <vector>

#include "bimap.h"

#include <symmetri/types.h>

namespace pn
{
  // Namespace regrouping all petri net functionnality

  /**
   * The arbitrary order on events
   * For convenience we use the shortlex order
   */
  struct short_lex
  {
    bool operator()(const std::string& lhs, const std::string& rhs) const noexcept
    {
      return lhs.size() < rhs.size() || (lhs.size() == rhs.size() && lhs < rhs);
    }
  };

  /**
   * Types that remain unchanged no matter the type of net
   */
  struct invar_pn_types
  {
    using pnames_t = std::string; /** Names for places */
    using pnamesid_t = uint; /** Place identifier */

    using tnames_t = std::string; /** Names for transitions */
    using tnamesid_t = uint; /** Transition identifier */

    using marking_t = std::vector<uint>; /** Marking of a PT-net; Number
    of tokens per place: pnamesid_t -> num tokens */

    using place_set_t = std::vector<pnamesid_t>; /** Use a _sorted_
    vector as a set of places */

    using place2id_t = bimap<pnames_t, pnamesid_t, short_lex, std::less<pnamesid_t>>; /** Bimap for identifaction */
    using trans2id_t = bimap<tnames_t, tnamesid_t, short_lex, std::less<pnamesid_t>>; /** Bimap for identifaction */
  };


  /**
   * Basic uncolored, unweighted petri net without extensions
   */
  struct basic_pn_types : public invar_pn_types
  {
    using transition_t = std::pair<place_set_t, place_set_t>; /** A transition is defined by
                         the input and output places */

    using all_trans_t = std::map<tnamesid_t, transition_t>; /** All the
    transitions in a PN */
  };

  /**
   * Types invariant for weighted petri net with inhibitor arcs
   */
  struct invar_weighted_pn_types : invar_pn_types
  {
    using weighted_place_t = std::pair<pnamesid_t, uint>; /** A weighted place */
    using weighted_place_set_t = std::vector<weighted_place_t>; /** A set of weighted places;
    Vector is sorted with respect to the places */
  };

  /**
   * Types for weighted petri nets
   */
  struct weighted_pn_types : public invar_weighted_pn_types
  {
    using transition_t = std::pair<weighted_place_set_t, weighted_place_set_t>; /** A transition
                         is defined by the weighted input and output places */

    using all_trans_t = std::map<tnamesid_t, transition_t>; /** Type to store all
    the transitions in a WPN */
  };

  /**
   * Types for weighted petri nets with inhibitor arcs
   */
  struct weighted_pni_types : invar_weighted_pn_types
  {

    using transition_t =
      std::pair<std::pair<place_set_t, weighted_place_set_t>, weighted_place_set_t>; /** A transition is defined by the
                                                (inhibitor places, weighted input) and output places */

    using all_trans_t = std::map<tnamesid_t, transition_t>; /** All the transitions
    in a  WPNI */

    static constexpr uint inhib_value = 999u; /** Not many passers / tools
    support inhibitor arcs, so we choose a special value to designate them */
  };


  /**
   * Configurable class to represent all nets
   * @tparam PNT Defines all necessary types
   * @tparam TSM Function defining how to initiate firing
   * @tparam TM Function defining how to terminate firing
   */
  template <class PNT, class TSM, class TM>
  class basic_pn
  {
  public:
    using pnames_t = PNT::pnames_t;
    using pnamesid_t = PNT::pnamesid_t;
    using tnames_t = PNT::tnames_t;
    using tnamesid_t = PNT::tnamesid_t;
    using marking_t = PNT::marking_t;
    using transition_t = PNT::transition_t;
    using all_trans_t = PNT::all_trans_t;

    using place2id_t = PNT::place2id_t;
    using trans2id_t = PNT::trans2id_t;

  protected:
    place2id_t place2id_; /** bidirectional map between places and their id */
    trans2id_t trans2id_; /** bidirectional map between transitions and their id */

    all_trans_t all_trans_;

    marking_t init_m_;

    //
    static inline TSM tsm_; /** (Lambda)-function implementing try_start_me_.
    That is try to fire the first part of a transition given a current marking */
    static inline TM tm_; /** (Lambda)-function implementing terminate_me_.
    That is to fire the second part of a transition given a current marking */


  public:
    basic_pn() = default;
    basic_pn(const TSM& tsm, const TM& tm)
    {
      tsm_ = tm;
      tm_ = tm;
    }

    // Getters
    const place2id_t& place2id() const noexcept { return place2id_; }
    const trans2id_t& trans2id() const noexcept { return trans2id_; }
    const all_trans_t& all_trans() const noexcept { return all_trans_; }
    const marking_t& init_marking() const noexcept { return init_m_; }

    /**
     * Tries to do start the transition \a t WITHOUT terminating it
     * (i.e. it only removes tokens from input places)
     * and changes \a m accordingly
     * If it can not be started, it leaves m unchanged
     * @param t transition id to start
     * @param m marking as in-out
     * @return true iff successfully fired
     * @pre It is the users responsibility to make sure
     * that the transition is compatible with the marking
     */
    static bool try_start_me(const transition_t& tt, marking_t& m) { return tsm_(tt, m); }

    /**
     * Terminates the transition \a t WITHOUT starting it
     * (i.e. it only adds tokens to target places)
     * @param t transition id to start
     * @param m marking as in-out
     * @pre It is the users responsibility to make sure
     * that the transition is compatible with the marking
     */
    static void terminate_me(const transition_t& tt, marking_t& m) { return tm_(tt, m); }

    /**
     * Tries to fire a given transition \a trans and change the marking \a m accordingly.
     * If it can not be fired, leave m unchanged
     * @param t transition id to fire
     * @param m marking as in-out
     * @return true iff successfully fired
     * @pre It is the users responsibility to make sure
     * that the transition is compatible with the marking
     */
    static bool try_fire_me(const transition_t& tt, marking_t& m)
    {
      if (try_start_me(tt, m))
      {
        terminate_me(tt, m);
        return true;
      }
      return false;
    }

    /**
     * @defgroup try_firing
     * @ingroup try_firing
     * Tries to fire transition \a t and change the marking \a m accordingly.
     * If it can not be fired, leave m unchanged
     * @param t transition id to fire
     * @param m marking as in-out
     * @return true iff successfully fired
     * @{
     */
    bool try_fire(tnamesid_t t, marking_t& m)
    {
      const auto& tt = all_trans_[t];
      return try_fire_me(tt, m);
    }
    bool try_fire(const tnames_t& t, marking_t& m) { return try_fire(trans2id_.at(t), m); }
    /** @} */
  };

  // Concrete implementation of starting and terminating functions

  /**
   * Implementing a firing tentative of a normal, unweighted petri net
   * @param tt basic pn transition
   * @param m basic pn marking
   * @return true iff the transition was fired
   */
  auto try_start_me_standard = [](const basic_pn_types::transition_t& tt, basic_pn_types::marking_t& m)
  {
    auto it = tt.first.begin();
    for (; it != tt.first.end(); ++it)
    {
      if (!m[*it])
        break;
      --m[*it];
    }
    if (it != tt.first.end())
    {
      // Undo changes
      // Current pos was not modified
      assert(m[*it] == 0);
      if (it == tt.first.begin())
        return false;
      while (true)
      {
        --it;
        ++m[*it];
        if (it == tt.first.begin())
          break;
      }
      return false;
    }
    return true;
  };
  /**
   * Terminates the transition \a t WITHOUT starting it
   * (i.e. it only adds tokens to target places)
   * @param t transition id to start
   * @param m marking as in-out
   * @pre It is the users responsibility to make sure
   * that the transition is compatible with the marking
   */
  auto terminate_me_standard = [](const basic_pn_types::transition_t& tt, basic_pn_types::marking_t& m)
  {
    // Successful firing
    for (auto it2 = tt.second.begin(); it2 != tt.second.end(); ++it2)
      ++m[*it2];
  };

  /**
   * Class for standard petri nets.
   * Neither weighted nor inhibited, derives from a suitable parametrised basic_pn class
   */
  class standard_pn_t
      : public basic_pn<basic_pn_types, decltype(try_start_me_standard), decltype(terminate_me_standard)>
  {
    using bpn = basic_pn<basic_pn_types, decltype(try_start_me_standard), decltype(terminate_me_standard)>;

  public:
    using place_set_t = basic_pn_types::place_set_t;

    /**
     * Transforms a parsed PN into the standard representation
     */
    standard_pn_t(const symmetri::Net& net, const symmetri::Marking& marking) : bpn()
    {

      // Loop once to get all names of transitions and places
      auto all_p = std::set<pnames_t, short_lex>{};
      auto all_t = std::set<tnames_t, short_lex>{};

      for (const auto& [nt, t] : net)
      {
        all_t.insert(nt);
        for (const auto& vec : {t.first, t.second}) // Impacted before and after
          for (const auto& [pn, _] : vec) // All place / token combos
            all_p.insert(pn);
      }
      // Initial marking
      for (const auto& [pn, _] : marking)
        all_p.insert(pn);

      // Set the indices
      auto idx = -1u;
      for (const auto& pn : all_p)
        place2id_.insert(pn, ++idx);
      idx = -1u;
      for (const auto& pt : all_t)
        trans2id_.insert(pt, ++idx);

      // Set the initial marking
      // Ignored sanity check: All tokens should have the same index
      init_m_ = marking_t(all_p.size(), 0);
      for (const auto& [pn, _] : marking)
        ++init_m_[place2id_[pn]];

      // Translate all the transitions
      for (const auto& [nt, t] : net)
      {
        auto nti = trans2id_[nt];
        auto [it, ins] =
          all_trans_.emplace(std::piecewise_construct, std::forward_as_tuple(nti),
                             std::forward_as_tuple(place_set_t(t.first.size()), place_set_t(t.second.size())));
        assert(ins);
        for (int i : {0, 1})
        {
          auto& int_rep = i == 0 ? it->second.first : it->second.second;
          auto& ext_rep = i == 0 ? t.first : t.second;
          int_rep.clear();
          for (const auto& [pn, _] : ext_rep)
            int_rep.push_back(place2id_[pn]);
          std::sort(int_rep.begin(), int_rep.end());
        }
      }
    }
  };

  // Implementations for weighted PN

  /**
   * Implementing a firing tentative of a weighted petri net
   * @param wt weighted pn transition
   * @param m basic pn marking
   * @return true iff the transition was fired
   * Actual implem followed by a wrapper
   */
  bool try_start_me_weighted_(const weighted_pn_types::weighted_place_set_t& wt, weighted_pn_types::marking_t& m)
  {
    auto it = wt.begin();
    for (; it != wt.end(); ++it)
    {
      const auto& [place, weight] = *it;
      if (m[place] < weight)
        break;
      m[place] -= weight;
    }
    if (it != wt.end())
    {
      // Undo changes
      // Current pos was not modified
      assert(m[it->first] < it->second);
      if (it == wt.begin())
        return false;
      while (true)
      {
        --it;
        const auto& [place, weight] = *it;
        m[place] += weight;
        if (it == wt.begin())
          break;
      }
      return false;
    }
    return true;
  }
  auto try_start_me_weighted = [](const weighted_pn_types::transition_t& tt, weighted_pn_types::marking_t& m)
  { return try_start_me_weighted_(tt.first, m); };

  /**
   * Terminates the transition \a t WITHOUT starting it
   * (i.e. it only adds tokens to target places)
   * @defgroup terminate_me
   * @ingroup terminate_me
   * @param t transition id to start
   * @param m marking as in-out
   * @pre It is the users responsibility to make sure
   * that the transition is compatible with the marking
   * @{
   */
  void terminate_me_weighted_(const weighted_pn_types::weighted_place_set_t& wt, weighted_pn_types::marking_t& m)
  {
    // Successful firing
    for (auto it2 = wt.begin(); it2 != wt.end(); ++it2)
    {
      const auto& [place, weight] = *it2;
      m[place] += weight;
    }
  };

  auto terminate_me_weighted = [](const weighted_pn_types::transition_t& tt, weighted_pn_types::marking_t& m)
  {
    // Successful firing
    return terminate_me_weighted_(tt.second, m);
  };
  /** @} */

  /**
   * Class for weighted petri nets.
   * Weighted but not inhibited, derives from a suitable parametrised basic_pn class
   */
  class weighted_pn_t
      : public basic_pn<weighted_pn_types, decltype(try_start_me_weighted), decltype(terminate_me_weighted)>
  {
    using wpn = basic_pn<weighted_pn_types, decltype(try_start_me_weighted), decltype(terminate_me_weighted)>;

  public:
    using weighted_place_t = weighted_pn_types::weighted_place_t;
    using weighted_place_set_t = weighted_pn_types::weighted_place_set_t;

    /**
     * Transforms a PN into this representation
     */
    weighted_pn_t(const symmetri::Net& net, const symmetri::Marking& marking) : wpn()
    {

      // Loop once to get all names of transitions and places
      auto all_p = std::set<pnames_t, short_lex>{};
      auto all_t = std::set<tnames_t, short_lex>{};

      for (const auto& [nt, t] : net)
      {
        all_t.insert(nt);
        for (const auto& vec : {t.first, t.second}) // Impacted before and after
          for (const auto& [pn, _] : vec) // All place / token combos
            all_p.insert(pn);
      }
      // Initial marking
      for (const auto& [pn, _] : marking)
        all_p.insert(pn);

      // Set the indices
      auto idx = -1u;
      for (const auto& pn : all_p)
        place2id_.insert(pn, ++idx);
      idx = -1u;
      for (const auto& pt : all_t)
        trans2id_.insert(pt, ++idx);

      // Set the initial marking
      // Ignored sanity check: All tokens should have the same index
      init_m_ = marking_t(all_p.size(), 0);
      for (const auto& [pn, _] : marking)
        ++init_m_[place2id_[pn]];

      // Translate all the transitions
      // Weighted transitions are representing by having the token n times
      // -> Sort it first into a set
      auto thistrans = std::pair<std::map<pnamesid_t, uint>, std::map<pnamesid_t, uint>>();
      for (const auto& [nt, t] : net)
      {
        auto nti = trans2id_[nt];
        auto [itfinal, insfinal] =
          all_trans_.emplace(std::piecewise_construct, std::forward_as_tuple(nti),
                             std::forward_as_tuple(weighted_place_set_t(0), weighted_place_set_t(0)));
        assert(insfinal);

        // Compute internal representation
        for (int i : {0, 1})
        {
          auto& int_rep = i == 0 ? thistrans.first : thistrans.second;
          auto& ext_rep = i == 0 ? t.first : t.second;
          auto& int_rep_final = i == 0 ? itfinal->second.first : itfinal->second.second;

          int_rep.clear();
          for (const auto& [pn, _] : ext_rep)
          {
            auto thisit = int_rep.try_emplace(place2id_[pn], 0u).first;
            ++(thisit->second);
          }
          // Insert
          int_rep_final.reserve(int_rep.size());
          for (const auto& [place, weight] : int_rep)
            int_rep_final.emplace_back(place, weight);
        }
      }
    }
  };

  // Implementations for weighted and inhibited PN

  /**
   * Implementing a firing tentative of a weighted and inhibited petri net
   * @param tt weighted and inhibited pn transition
   * @param m basic pn marking
   * @return true iff the transition was fired
   */
  auto try_start_me_weighted_inhib = [](const weighted_pni_types::transition_t& tt, weighted_pni_types::marking_t& m)
  {
    // First test for inhibiting places
    if (std::any_of(tt.first.first.cbegin(), tt.first.first.cend(), [&m](const auto& idx) -> bool { return m[idx]; }))
      return false; // Was inhibited
    // Check the weighted part
    return try_start_me_weighted_(tt.first.second, m);
  };
  /**
   * Terminates the transition \a t WITHOUT starting it
   * (i.e. it only adds tokens to target places)
   * @param t transition id to start
   * @param m marking as in-out
   * @pre It is the users responsibility to make sure
   * that the transition is compatible with the marking
   */
  auto terminate_me_weighted_inhib = [](const weighted_pni_types::transition_t& tt, weighted_pn_types::marking_t& m)
  {
    // Successful firing
    // Inhibition is only important for starting
    return terminate_me_weighted_(tt.second, m);
  };

  /**
  * Class for standard petri nets.
  * Neither weighted nor inhibited, derives from a suitable parametrised basic_pn class
  */
  class weighted_pni_t : public basic_pn<weighted_pni_types, decltype(try_start_me_weighted_inhib),
                                         decltype(terminate_me_weighted_inhib)>
  {
    using wpni =
      basic_pn<weighted_pni_types, decltype(try_start_me_weighted_inhib), decltype(terminate_me_weighted_inhib)>;

  public:
    using place_set_t = weighted_pni_types::place_set_t;
    using weighted_place_t = weighted_pni_types::weighted_place_t;
    using weighted_place_set_t = weighted_pni_types::weighted_place_set_t;

    /**
    * Transforms a PN into this representation
    */
    weighted_pni_t(const symmetri::Net& net, const symmetri::Marking& marking) : wpni()
    {

      // Loop once to get all names of transitions and places
      auto all_p = std::set<pnames_t, short_lex>{};
      auto all_t = std::set<tnames_t, short_lex>{};

      for (const auto& [nt, t] : net)
      {
        all_t.insert(nt);
        for (const auto& vec : {t.first, t.second}) // Impacted before and after
          for (const auto& [pn, _] : vec) // All place / token combos
            all_p.insert(pn);
      }
      // Initial marking
      for (const auto& [pn, _] : marking)
        all_p.insert(pn);

      // Set the indices
      auto idx = -1u;
      for (const auto& pn : all_p)
        place2id_.insert(pn, ++idx);
      idx = -1u;
      for (const auto& pt : all_t)
        trans2id_.insert(pt, ++idx);

      // Set the initial marking
      // Ignored sanity check: All tokens should have the same index
      init_m_ = marking_t(all_p.size(), 0);
      for (const auto& [pn, _] : marking)
        ++init_m_[place2id_[pn]];

      // Translate all the transitions
      // Weighted transitions are representing by having the token n times
      // -> Sort it first into a set
      auto int_rep = std::map<pnamesid_t, uint>();
      for (const auto& [nt, t] : net)
      {
        auto nti = trans2id_[nt];
        auto [itfinal, insfinal] = all_trans_.emplace(
          std::piecewise_construct, std::forward_as_tuple(nti),
          std::forward_as_tuple(std::make_pair(place_set_t(0), weighted_place_set_t(0)), weighted_place_set_t(0)));
        assert(insfinal);

        // Compute internal representation
        // Now "input" and "output" are structurally different
        // Input
        int_rep.clear();
        for (const auto& [pn, _] : t.first)
        {
          auto thisit = int_rep.try_emplace(place2id_[pn], 0u).first;
          ++(thisit->second);
        }
        for (const auto& [place, weight] : int_rep)
        {
          if (weight == weighted_pni_types::inhib_value)
            itfinal->second.first.first.push_back(place);
          else
            itfinal->second.first.second.emplace_back(place, weight);
        }
        // Output
        int_rep.clear();
        for (const auto& [pn, _] : t.second)
        {
          auto thisit = int_rep.try_emplace(place2id_[pn], 0u).first;
          ++(thisit->second);
        }
        itfinal->second.second.reserve(int_rep.size());
        for (const auto& [place, weight] : int_rep)
          itfinal->second.second.emplace_back(place, weight);
      }
    }
  };


} // namespace pn
