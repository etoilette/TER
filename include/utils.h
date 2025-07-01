#pragma once

#include <sstream>

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