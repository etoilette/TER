#pragma once

#include <sstream>
#include <ranges>
#include <numeric>

template <class CONT>
std::string join(const CONT& c, const std::string& sep = ",")
{
  auto out = std::ostringstream();
  for (const auto& e : c)
    out << e << ",";

  auto outs = out.str();
  if (outs.empty())
    return outs;
  else
    return outs.substr(0, outs.size() - sep.size());
}

std::size_t
accu_size_sec(const auto& c)
{
  return std::accumulate(c.begin(), c.end(), (size_t)0,
             [](auto lhs, const auto& ev)
             {return lhs + ev.second.size(); });
}