#include <algorithm>
#include <cmath>
#include <cstring>
#include <chrono>
#include <iostream>
#include <symmetri/parsers.h>
#include "convmaxhda.h"
#include "pn.h"
#include "pnhda.h"
#include "standardmaxhda.h"

struct voider
{
  std::string operator()(const auto& mhda) const
  {
    return "";
  }
};

template <class MAXHDA, class DATA = voider>
void
act(MAXHDA& maxhda, const std::string& act, const std::string& header,
    DATA&& data)
{
  const auto print_act = std::set<std::string>{"str", "dot"};

  if (std::any_of(print_act.begin(), print_act.end(),
                  [&act](const auto& s)->bool{return act.starts_with(s); }))
    std::cout << maxhda.format(act) << std::endl;
  else if (act == "csv")
    std::cout << maxhda.get_csv_stats(false) << std::endl;
  else if (act == "csv_header")
    std::cout << maxhda.get_csv_stats(true) << std::endl;
  else if (act == "csv_time")
  {
    if (std::is_same_v<DATA, voider>)
      throw std::runtime_error("csv with add. info can not use "
                               "default args");
    std::cout << maxhda.get_csv_stats(false, header, data) << std::endl;
  }
  else if (act == "csv_time_header")
  {
    if (std::is_same_v<DATA, voider>)
      throw std::runtime_error("csv with add. info and header "
                               "can not use default args");

    std::cout << maxhda.get_csv_stats(true, header, data) << std::endl;
  }
  else
    throw std::invalid_argument("unknown action: " + act);
}

struct timer_
{
  std::chrono::time_point<std::chrono::high_resolution_clock> t0_;

  void start()
  {
    t0_ = std::chrono::high_resolution_clock::now();
  }
  [[nodiscard]] size_t time() const
  {
    return std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::high_resolution_clock::now() - t0_).count();
  }
};

int main(int argc, char *argv[])
{
  if (argc == 2 && (strcmp(argv[1], "-h") || strcmp(argv[1], "--help")))
  {
    std::cout << "pnhda type file act\n"
              << "This tool provides different strategies to "
                 "obtain maxHDAs from pnml files. Finally different post-processes "
                 "can be choosen.\n"
                 "\"standard\": Standard PT net converted to detailed max HDA\n"
                 "\"standardmaxhda_dfs\": Direct translation to maxHDA, "
                 "first algo, dfs\n"
                 "\"standardmaxhda_bfs\": Direct translation to maxHDA, "
                 "first algo, bfs\n"
                 "\"standardmaxhda2_dfs\": Direct translation to maxHDA, "
                 "second algo, dfs\n"
                 "\"standardmaxhda2_bfs\": Direct translation to maxHDA, "
                 "second algo, bfs\n"
                 "file is the path to the pnml file\n"
                 "and finally act describes what to do:\n"
                 "\"str\": Generate a text-based representation of the maxHDA\n"
                 "\"dot\": Generate the dot file for the corresponding "
                 "ST-automaton\n"
                 "\"csv\": Produces standard csv data\n"
                 "\"csv_header\": Produces standard csv header and data\n"
                 "\"csv_time\": Produces standard csv data plus timing info\n"
                 "\"csv_time_header\": Produces standard csv header and data "
                 "for standard and timing info\n";
    return 0;
  }
  if (argc != 4)
    throw std::invalid_argument("Wrong number of arguments");

  // Parse the pnml file
  auto [n, m] = symmetri::readPnml({std::string(argv[2])});

  auto type_s = std::string(argv[1]);

  auto timer = timer_();

  if (type_s == "standard")
  {
    auto spn = pn::standard_pn_t(n, m);

    timer.start();
    auto spnhda = pnhda::standard_pn_hda(spn);
    auto tCstr = timer.time();

    timer.start();
    auto maxhda = convmaxhda::standard_max_conv_hda(spnhda);
    auto tConv = timer.time();

    auto d = [tCstrs=std::to_string(tCstr), tConvs = std::to_string(tConv)](const auto&)
      {return "," + tCstrs + "," + tConvs; };

    act(maxhda, argv[3], ",tCstr,tConv", d);
  }
  else if (type_s == "standardmaxhda_dfs")
  {
    auto spn = pn::standard_pn_t(n, m);

    timer.start();
    auto spnmaxhda =
      standardmaxhda::standard_pn_maxhda(std::make_pair(&spn, 0));
    auto tCstr = timer.time();

    auto d = [tCstrs=std::to_string(tCstr)](const auto&)
      {return "," + tCstrs; };

    act(spnmaxhda, argv[3], ",tCstr", d);
  }
  else if (type_s == "standardmaxhda_bfs")
  {
    auto spn = pn::standard_pn_t(n, m);

    timer.start();
    auto spnmaxhda =
      standardmaxhda::standard_pn_maxhda(std::make_pair(&spn, 1));
    auto tCstr = timer.time();

    auto d = [tCstrs=std::to_string(tCstr)](const auto&)
      {return "," + tCstrs; };

    act(spnmaxhda, argv[3], ",tCstr", d);
  }
  else if (type_s == "standardmaxhda2_dfs")
  {
    auto spn = pn::standard_pn_t(n, m);

    timer.start();
    auto spnmaxhda =
      standardmaxhda::standard_pn_maxhda(std::make_pair(&spn, 2));
    auto tCstr = timer.time();

    auto d = [tCstrs=std::to_string(tCstr)](const auto&)
      {return "," + tCstrs; };

    act(spnmaxhda, argv[3], ",tCstr", d);
  }
  else if (type_s == "standardmaxhda2_bfs")
  {
    auto spn = pn::standard_pn_t(n, m);

    timer.start();
    auto spnmaxhda =
      standardmaxhda::standard_pn_maxhda(std::make_pair(&spn, 3));
    auto tCstr = timer.time();

    auto d = [tCstrs=std::to_string(tCstr)](const auto&)
    {return ","+tCstrs; };

    act(spnmaxhda, argv[3], ",tCstr", d);
  }
  else
    throw std::invalid_argument("unknown type: " + type_s);

  return 0;
}
