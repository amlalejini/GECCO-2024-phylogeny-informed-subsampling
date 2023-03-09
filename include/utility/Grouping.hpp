#pragma once

#include <utility>
#include <algorithm>

#include "emp/base/vector.hpp"

namespace utils {

class GroupManager;

class Grouping {
public:
  friend class GroupManager;
protected:
  size_t group_id = 0;
  size_t group_size = 0;
  emp::vector<size_t> member_ids;

public:
  Grouping() = default;

  void Resize(size_t size, size_t id=0) {
    group_size = size;
    member_ids.resize(group_size, id);
  }

  size_t GetGroupID() const {
    return group_id;
  }

  size_t GetSize() const {
    return member_ids.size();
  }

  void SetGroupID(size_t id) {
    group_id = id;
  }

  emp::vector<size_t>& GetMembers() {
    return member_ids;
  }

  const emp::vector<size_t>& GetMembers() const {
    return member_ids;
  }

};


class GroupManager {
public:



protected:

  emp::vector<size_t> possible_ids;       ///< List of all possible entity IDS (things to be grouped)
  emp::vector<size_t> member_group_assignments;  ///< Group ID of each entity.
  emp::vector<Grouping> groupings;

  std::function<void(void)> assign_groupings;

  bool dirty_groups = true;

public:

  void SetPossibleIDs(const emp::vector<size_t>& ids) {
    // Update possible ids
    possible_ids = ids;
    // Resize and clear out exiting group assignments
    member_group_assignments.resize(ids.size(), 0);
    std::fill(
      member_group_assignments.begin(),
      member_group_assignments.end(),
      0
    );

    dirty_groups = true;
  }

  const emp::vector<size_t>& GetPossibleIDs() const {
    return possible_ids;
  }

  size_t GetNumGroups() const {
    return groupings.size();
  }

  const Grouping& GetGroup(size_t group_id) const {
    return groupings[group_id];
  }

  /// @brief Configure SingleGrouping mode.
  //  All members are assigned to a single group.
  void SetSingleGroupMode() {
    // Resize groupings to 1
    groupings.resize(1);
    // Set group ID of single group to 0
    auto& group = groupings.back();
    group.SetGroupID(0);
    // Resize group to # of members
    group.Resize(possible_ids.size(), 0);
    // Copy member ids into group
    auto& group_member_ids = group.GetMembers();
    std::copy(
      possible_ids.begin(),
      possible_ids.end(),
      group_member_ids.begin()
    );
    // Update group assignments for each member
    std::fill(
      member_group_assignments.begin(),
      member_group_assignments.end(),
      0
    );
    // Group is static, so do nothing on group update.
    assign_groupings = [this]() {
      emp_assert(groupings.size() == 1);
      emp_assert(groupings.back().GetMembers().size() == possible_ids.size());
      /* Do nothing! */
    };
  }

  void UpdateGroupings() {
    assign_groupings();
  }

};

}