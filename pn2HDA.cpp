#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <symmetri/parsers.h>
#include "pn.h"
#include "pnhda.h"

template<class PNHDA>
void
act(PNHDA& pnhda, const std::string& act)
{
  const auto print_act = std::set<std::string>{"str", "str_max", "dot"};

  auto get_avg_occ = [](const auto& pnhda)
    {
      constexpr size_t nBins = 9;

      auto nS = pnhda.mark_ess().get_fwd().begin()->first.size();
      auto res = std::vector<double>(nBins, 0); // Last bin is for size

      auto toBin = [nBins, &nS, nsD = (double) nS](size_t count)
        {
          count -= (count == nS);
          return (std::size_t) std::floor(((double) count / nsD) * nBins);
        };

      res.back() = (double) nS;
      for (const auto& [c, cid] : pnhda.X_ext().get_fwd())
      {
        const auto& mid = c.first;
        const auto& m = pnhda.mark_ess().at(mid);

        res[toBin(nS - std::count(m.begin(), m.end(), 0))] += 1;
      }

      auto rs = std::ostringstream();
      rs << ",[";
      for (const auto& mr : res)
        rs << mr << ';';
      rs << ']';

      return rs.str();
    };

  auto count_max_cell = [](const auto& pnhda)
    {
      auto all_max_cells = pnhda.get_max_cells();

      // Dimension to occurrence
      auto cell_count = std::map<unsigned, unsigned>{};
      for (auto idx : all_max_cells)
      {
        const auto& cf = pnhda.X_ext().at(idx);
        const auto& conc = pnhda.square_ess().at(cf.second);
        cell_count[conc.size()] += 1;
      }

      auto rs = std::ostringstream();
      rs << ",[";
      for (const auto& mr : cell_count)
        rs << '[' << mr.first << ';' << mr.second << "];";
      rs << "]";
      return rs.str();
    };

  auto all_f = [&](const auto& pnhda)
    {
      return get_avg_occ(pnhda) + count_max_cell(pnhda);
    };

  if (std::any_of(print_act.begin(), print_act.end(),
                  [&act](const auto& s)->bool{return act.starts_with(s); }))
    std::cout << pnhda.format(act) << std::endl;
  else if (act == "csv")
    std::cout << pnhda.get_csv_stats(false, ",mOcc,MaxCells", all_f) << std::endl;
  else if (act == "csv_header")
    std::cout << pnhda.get_csv_stats(true, ",mOcc,MaxCells", all_f) << std::endl;
  else
    throw std::invalid_argument("unknown action: " + act);

}

int main(int argc, char *argv[])
{
  if (argc == 2 && (strcmp(argv[1], "-h") || strcmp(argv[1], "--help")))
  {
    std::cout << "pnhda type file act\n"
              << "Translate the given pnml of type\n"
                 "\"standard\": Standard PT net\n"
                 "\"weighted\": Weighted PT net\n"
                 "\"wpni\": Weighted PT net with inhibitor arcs. Inhibitor arcs are marked\n"
                 "with a special weight, 999 by default\n"
                 "file is the path to the pnml file\n"
                 "and finally act describes what to do:\n"
                 "\"base\": Output basic informations of the Petri Net\n"
                 "\"str\": Generate a text-based representation of the HDA\n"
                 "\"str_max\": Generate a text-based representation of the maximal representation of the HDA\n"
                 "\"dot\": Generate the dot file for the corresponding ST-automaton\n";
    return 0;
  }
  if (argc != 4)
    throw std::invalid_argument("Wrong number of arguments");

  // Parse the pnml file
  auto [n, m] = symmetri::readPnml({std::string(argv[2])});

  auto type_s = std::string(argv[1]);
  auto act_s = std::string(argv[3]);

  if (type_s == "standard")
  {
    auto spn = pn::standard_pn_t(n, m);
    if (act_s == "base")
      std::cout << spn << std::endl;
    else
    {
      auto spnhda = pnhda::standard_pn_hda(spn);
      act(spnhda, argv[3]);
    }
    return 0;
  }
  if (type_s == "weighted")
  {
    auto spn = pn::weighted_pn_t(n, m);
    if (act_s == "base")
      std::cout << spn << std::endl;
    else
    {
      auto spnhda = pnhda::weighted_pn_hda(spn);
      act(spnhda, argv[3]);
    }
    return 0;
  }
  if (type_s == "wpni")
  {
    auto spn = pn::weighted_pni_t(n, m);

    if (act_s == "base")
      std::cout << spn << std::endl;
    else
    {
      auto spnhda = pnhda::weighted_pni_hda(spn);
      act(spnhda, argv[3]);
    }
    return 0;
  }

  throw std::invalid_argument("unknown type: " + type_s);
}
