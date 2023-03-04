#pragma once

// Standard includes
#include <ostream>
#include <unordered_map>

// Empirical includes
#include "emp/bits/BitSet.hpp"

// SignalGP includes
#include "sgp/EventLibrary.hpp"

/// Basic event type!
/// - Just contains a tag.
template<size_t TAG_SIZE>
struct Event : public sgp::BaseEvent {
  using tag_t = emp::BitSet<TAG_SIZE>;
  tag_t tag;

  Event(size_t _id, const tag_t& _tag)
    : BaseEvent(_id), tag(_tag)
  { ; }

  tag_t& GetTag() { return tag; }
  const tag_t& GetTag() const { return tag; }

  void Print(std::ostream& os) const {
    os << "{id:" << GetID() << ",tag:";
    tag.Print(os);
    os << "}";
  }
};

/// Message event type
/// - contains a tag and data
template<size_t TAG_SIZE>
struct MessageEvent : public Event<TAG_SIZE> {
  using tag_t = typename Event<TAG_SIZE>::tag_t;
  using data_t = std::unordered_map<int, double>;
  data_t data;

  MessageEvent(size_t _id, const tag_t& _tag, const data_t& _data=data_t())
    : Event<TAG_SIZE>(_id, _tag), data(_data)
  { ; }

  data_t & GetData() { return data; }
  const data_t & GetData() const { return data; }
};