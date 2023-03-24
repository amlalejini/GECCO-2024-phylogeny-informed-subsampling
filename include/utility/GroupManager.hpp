#pragma once

#include <utility>
#include <algorithm>

#include "emp/base/vector.hpp"
#include "emp/math/Random.hpp"

namespace utils {

class GroupManager;

class Grouping {
public:
  friend class GroupManager;
protected:
  size_t group_id = 0;
  // size_t group_size = 0;
  emp::vector<size_t> member_ids;

public:
  Grouping() = default;

  void Resize(size_t size, size_t id=0) {
    // group_size = size;
    member_ids.resize(size, id);
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
  emp::Random& random;
  emp::vector<size_t> possible_ids;       ///< List of all possible entity IDS (things to be grouped)
  emp::vector<size_t> member_group_assignments;  ///< Group ID of each entity.
  emp::vector<Grouping> groupings;

  std::function<void(void)> assign_groupings;

  bool dirty_groups = true;

public:

  GroupManager(emp::Random& rnd) : random(rnd) { ; }

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

  size_t GetMemberGroupID(size_t member_id) const {
    emp_assert(member_id < member_group_assignments.size());
    return member_group_assignments[member_id];
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

  void SetCohortsMode(size_t num_cohorts) {
    emp_assert(num_cohorts > 0);
    // Compute cohort sizes
    const size_t total_members = possible_ids.size();
    const size_t base_cohort_size = (size_t)(total_members / num_cohorts);
    size_t leftover_members = total_members - (base_cohort_size * num_cohorts);
    groupings.resize(num_cohorts);
    for (size_t i = 0; i < groupings.size(); ++i) {
      size_t group_size = base_cohort_size;
      if (leftover_members > 0) {
        ++group_size;
        --leftover_members;
      }
      auto& group = groupings[i];
      group.SetGroupID(i);
      group.Resize(group_size, 0);
    }

    assign_groupings = [this]() {
      // Shuffle all possible test ids
      emp::Shuffle(random, possible_ids);
      // Assign to cohorts in shuffled order
      size_t cur_pos = 0;
      for (size_t cohort_id = 0; cohort_id < groupings.size(); ++cohort_id) {
        auto& cohort = groupings[cohort_id];
        for (size_t member_i = 0; member_i < cohort.GetSize(); ++member_i) {
          const size_t member_id = possible_ids[cur_pos];
          cohort.member_ids[member_i] = member_id;
          member_group_assignments[member_id] = cohort.GetGroupID();
          ++cur_pos;
        }
      }
      emp_assert(cur_pos == possible_ids.size());
    };
  }

  void SetDownSampleMode(size_t sample_size) {
    // Initialize single group to hold down-sample
    groupings.resize(1);
    auto& group = groupings.back();
    group.SetGroupID(0);
    group.Resize(sample_size, 0);
    // Configure assign_groupings
    assign_groupings = [this, sample_size]() {
      emp_assert(groupings.size() == 1);
      // Get focal group
      auto& group = groupings.back();
      emp_assert(group.member_ids.size() == sample_size);
      // Shuffle all possible ids
      emp::Shuffle(random, possible_ids);
      // Assign to down-sample in shuffled order
      for (size_t i = 0; i < sample_size; ++i) {
        const size_t member_id = possible_ids[i];
        group.member_ids[i] = member_id;
        member_group_assignments[i] = 0;
      }
      for (size_t i = sample_size; i < possible_ids.size(); ++i) {
        const size_t member_id = possible_ids[i];
        member_group_assignments[member_id] = 1;
      }
    };
  }

  void UpdateGroupings() {
    assign_groupings();
  }

};

}