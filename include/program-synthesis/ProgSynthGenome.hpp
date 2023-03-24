#pragma once

namespace psynth {

template<typename PROGRAM_T>
struct ProgSynthGenome {
  using this_t = ProgSynthGenome<PROGRAM_T>;
  using program_t = PROGRAM_T;

  program_t program;

  ProgSynthGenome(const program_t& p) : program(p) {}
  ProgSynthGenome(const this_t&) = default;
  ProgSynthGenome(this_t&&) = default;

  bool operator==(const this_t& other) const {
    return program == other.program;
  }

  bool operator!=(const this_t& other) const {
    return !(*this == other);
  }

  bool operator<(const this_t& other) const {
    return program < other.program;
  }

  const program_t& GetProgram() const { return program; }
  program_t& GetProgram() { return program; }
};

}