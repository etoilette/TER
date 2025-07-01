#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <symmetri/parsers.h>
#include "convmaxhda.h"
#include "pn.h"
#include "pnhda.h"

template <class MAXHDA>
void
act(MAXHDA& maxhda, const std::string& act)
{
  const auto print_act = std::set<std::string>{"str", "dot"};

  if (std::any_of(print_act.begin(), print_act.end(),
                  [&act](const auto& s)->bool{return act.starts_with(s); }))
    std::cout << maxhda.format(act) << std::endl;
  else if (act == "csv")
    std::cout << maxhda.get_csv_stats(false) << std::endl;
  else if (act == "csv_header")
    std::cout << maxhda.get_csv_stats(true) << std::endl;
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
                 "\"weighted\": Weighted PT net.\n"
                 "file is the path to the pnml file\n"
                 "and finally act describes what to do:\n"
                 "\"str\": Generate a text-based representation of the maxHDA\n"
                 "\"dot\": Generate the dot file for the corresponding "
                 "ST-automaton\n"
                 "Note that we currently first translate the net into an explicit "
                 "HDA, then reduce it.";
    return 0;
  }
  if (argc != 4)
    throw std::invalid_argument("Wrong number of arguments");

  // Parse the pnml file
  auto [n, m] = symmetri::readPnml({std::string(argv[2])});

  auto type_s = std::string(argv[1]);

  if (type_s == "standard")
  {
    auto spn = pn::standard_pn_t(n, m);
    auto spnhda = pnhda::standard_pn_hda(spn);
    auto maxhda = convmaxhda::standard_max_conv_hda(spnhda);
    act(maxhda, argv[3]);
  }
  else if (type_s == "weighted")
  {
    auto spn = pn::weighted_pn_t(n, m);
    auto spnhda = pnhda::weighted_pn_hda(spn);
    auto maxhda = convmaxhda::weighted_max_conv_hda(spnhda);
    act(maxhda, argv[3]);
  }
  else
    throw std::invalid_argument("unknown type: " + type_s);

  return 0;
}
