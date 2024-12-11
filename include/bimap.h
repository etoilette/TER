#pragma once
#include <map>
#include <stdexcept>

template <typename A, typename B, typename CA = std::less<A>, typename CB = std::less<B>>
class bimap
{
  using fw_t = std::map<A, B, CA>;
  using fb_t = std::map<B, A, CB>;
  fw_t forward_map_;
  fb_t reverse_map_;

  public:

  auto size() const noexcept
  {
    assert(forward_map_.size() == reverse_map_.size());
    return forward_map_.size();
  }

  bool empty() const
  {
    assert(forward_map_.size() == reverse_map_.size());
    return forward_map_.empty();
  }


  void insert(const A& a, const B& b)
  {
    forward_map_[a] = b;
    reverse_map_[b] = a;
  }

  template< class... Args >
  auto try_emplace( const A& a, Args&&... args )
  {
    auto res = forward_map_.try_emplace( a, std::forward<Args>(args)... );
    if (res.second)
      reverse_map_.emplace( res.first->second, res.first->first );
    return res;
  }
  template< class... Args >
  auto try_emplace( const B& b, Args&&... args )
  {
    auto res = reverse_map_.try_emplace( b, std::forward<Args>(args)... );
    if (res.second)
      forward_map_.emplace( res.first->second, res.first->first );
    return res;
  }


  A& operator[](const B& key) { return reverse_map_.at(key); }

  B& operator[](const A& key) { return forward_map_.at(key); }

  const A& at(const B& key) const { return reverse_map_.at(key); }

  const B& at(const A& key) const { return forward_map_.at(key); }

  auto find(const B& key) { return reverse_map_.find(key); }
  auto find(const A& key) { return forward_map_.find(key); }
  auto find(const B& key) const { return reverse_map_.find(key); }
  auto find(const A& key) const { return forward_map_.find(key); }

  // Check if value exists
  bool contains_value(const A& value) const { return forward_map_.contains(value); }
  bool contains_value(const B& value) const { return reverse_map_.contains(value); }

  // Erase by key
  void erase_key(const A& key)
  {
    if (auto it = forward_map_.find(key); it != forward_map_.end())
    {
      assert(reverse_map_.contains(it->second));
      reverse_map_.erase(it->second);
      forward_map_.erase(it);
    }
  }

  // Erase by value
  void erase_value(const B& value)
  {
    if (auto it = reverse_map_.find(value); it != reverse_map_.end())
    {
      assert(forward_map_.contains(it->second));
      forward_map_.erase(it->second);
      reverse_map_.erase(it);
    }
  }

  // Clear the bimap
  void clear()
  {
    forward_map_.clear();
    reverse_map_.clear();
  }

  const auto& get_fwd() const noexcept { return forward_map_; }
  const auto& get_rev() const noexcept { return reverse_map_; }


};
