#pragma once

#include <iostream>
#include <utility>

#include "emp/base/vector.hpp"
#include "emp/bits/BitSet.hpp"
#include "emp/tools/string_utils.hpp"

#include "sgp/cpu/lfunprg/LinearFunctionsProgram.hpp"
#include "sgp/cpu/linprg/Instruction.hpp"


namespace psynth {

namespace internal {

template<typename TAG_T>
void PrintTagsJSON(
  std::ostream& out,
  const emp::vector<TAG_T>& tags
) {
  out << "[";
  for (size_t tag_i = 0; tag_i < tags.size(); ++tag_i) {
    if (tag_i) out << ",";
    out << "\"" << tags[tag_i] << "\"";
  }
  out << "]";
}

template<typename INST_LIB_T, typename INST_T>
void PrintInstructionJSON(
  std::ostream& out,
  const INST_LIB_T& ilib,
  const INST_T& inst
) {
  size_t inst_id = inst.GetID();
  const std::string& inst_name = ilib.GetName(inst_id);
  const auto& inst_args = inst.GetArgs();

  out << "{";
  out << "\"name\":\"" << inst_name << "\",";
  out << "\"tags\":";
  PrintTagsJSON(out, inst.GetTags());
  out << ",";
  out << "\"args\":[";
  for (size_t arg_i = 0; arg_i < inst_args.size(); ++arg_i) {
    if (arg_i) out << ",";
    out << "\"" << inst_args[arg_i] << "\"";
  }
  out << "]";
  out << "}";
}

}

template<typename INST_LIB_T, typename TAG_T, typename INST_ARG_T>
void PrintProgramJSON(
  std::ostream& out,
  const sgp::cpu::lfunprg::LinearFunctionsProgram<TAG_T, INST_ARG_T>& program,
  const INST_LIB_T& ilib
) {
  out << "{";
  out << "\"functions\":[";
  for (size_t func_i = 0; func_i < program.GetSize(); ++func_i) {
    auto& function = program[func_i];
    if (func_i) out << ",";
    out << "{";
    out << "\"id\":" << func_i << ",";
    out << "\"tags\":";
    internal::PrintTagsJSON(out, function.GetTags());
    out << ",";
    out << "\"instructions\":[";
    for (size_t inst_i = 0; inst_i < function.GetSize(); ++inst_i) {
      if (inst_i) out << ",";
      internal::PrintInstructionJSON(
        out,
        ilib,
        function[inst_i]
      );
    }
    out << "]"; // End instructions
    out << "}"; // End function
  }
  out << "]"; // End functions
  out << "}";
}

template<size_t TAG_WIDTH>
emp::BitSet<TAG_WIDTH> FromString_BitSet(
  const std::string& bit_str
) {
  using tag_t = emp::BitSet<TAG_WIDTH>;
  tag_t tag;
  tag.Clear();
  for (size_t i = 0; i < bit_str.size(); ++i) {
    if (i >= tag.GetSize()) {
      break;
    }
    if (bit_str[i] == '1') {
      tag.Set(tag.GetSize() - i - 1, true);
    }
  }
  return tag;
}



/// Load linear functions program from print format
/// Example:
template<typename INST_LIB_T, size_t TAG_WIDTH>
sgp::cpu::lfunprg::LinearFunctionsProgram<emp::BitSet<TAG_WIDTH>, int> LoadLinearFunctionsProgram_PrintFormat(
  std::istream& input,
  const INST_LIB_T& inst_lib
) {
  using program_t = sgp::cpu::lfunprg::LinearFunctionsProgram<emp::BitSet<TAG_WIDTH>, int>;
  using tag_t = emp::BitSet<TAG_WIDTH>;
  program_t program;
  std::string cur_line;
  emp::vector<std::string> line_components;
  while (!input.eof()) {
    std::getline(input, cur_line);
    // Remove all whitespace
    emp::remove_whitespace(cur_line);
    // Handle comments
    cur_line = emp::string_get(cur_line, "#");
    // If line is empty, skip.
    if (cur_line == emp::empty_string()) {
      continue;
    }
    // Collect tags
    emp::vector<tag_t> tags;
    // Grab tags
    std::string tags_str = emp::string_pop(
      cur_line,
      ")"
    );
    // Pop () off tags list
    emp::remove_chars(tags_str, "()");
    emp::slice(tags_str, line_components, ',');
    for (const auto& tag_str : line_components) {
      tags.emplace_back(FromString_BitSet<TAG_WIDTH>(tag_str));
    }
    // Function definition or instruction definition?
    emp::slice(cur_line, line_components, '-');
    if (emp::to_lower(line_components[0]) == "fn" && line_components.size() > 1) {
      // Function definition
      program.PushFunction(tags);
    } else {
      // Instruction definition
      // Get instruction name
      std::string inst_name = emp::string_get(cur_line, "[");
      // Isolate and parse instruction arguments
      emp::string_pop(cur_line, "[");
      emp::remove_chars(cur_line, "[]");
      emp::vector<int> args;
      emp::slice(cur_line, line_components, ',');
      for (const auto& arg_str : line_components) {
        args.emplace_back(emp::from_string<int>(arg_str));
      }
      program.PushInst(
        inst_lib,
        inst_name,
        args,
        tags
      );
    }
  }
  return program;
}

/// Load linear functions program from print format
/// - Each program is ID'd
template<typename INST_LIB_T, size_t TAG_WIDTH>
emp::vector<
  std::pair<
    sgp::cpu::lfunprg::LinearFunctionsProgram<emp::BitSet<TAG_WIDTH>, int>,
    std::string
  >
> LoadLinearFunctionsPrograms_PrintFormat(
  std::istream& input,
  const INST_LIB_T& inst_lib
) {
  using program_t = sgp::cpu::lfunprg::LinearFunctionsProgram<emp::BitSet<TAG_WIDTH>, int>;
  using tag_t = emp::BitSet<TAG_WIDTH>;

  emp::vector< std::pair<program_t, std::string> > programs;

  std::string cur_line;
  emp::vector<std::string> line_components;
  while (!input.eof()) {
    std::getline(input, cur_line);

    // Remove all whitespace
    emp::remove_whitespace(cur_line);

    // Handle comments
    cur_line = emp::string_get(cur_line, "#");

    // If line is empty, skip.
    if (cur_line == emp::empty_string()) {
      continue;
    }

    // New program, collect id
    if (emp::has_prefix(cur_line, "!")) {
      emp::string_pop_fixed(cur_line, 1);
      programs.emplace_back(
        std::make_pair(
          program_t(),
          cur_line
        )
      );
      continue;
    }

    // Program wasn't denoted by id
    if (programs.size() == 0) {
      programs.emplace_back(
        std::make_pair(
          program_t(),
          "unknown"
        )
      );
    }

    // Collect tags
    emp::vector<tag_t> tags;
    // Grab tags
    std::string tags_str = emp::string_pop(
      cur_line,
      ")"
    );
    // Pop () off tags list
    emp::remove_chars(tags_str, "()");
    emp::slice(tags_str, line_components, ',');
    for (const auto& tag_str : line_components) {
      tags.emplace_back(FromString_BitSet<TAG_WIDTH>(tag_str));
    }
    // Function definition or instruction definition?
    emp::slice(cur_line, line_components, '-');
    if (emp::to_lower(line_components[0]) == "fn" && line_components.size() > 1) {
      // Function definition
      programs.back().first.PushFunction(tags);
    } else {
      // Instruction definition
      // Get instruction name
      std::string inst_name = emp::string_get(cur_line, "[");
      // Isolate and parse instruction arguments
      emp::string_pop(cur_line, "[");
      emp::remove_chars(cur_line, "[]");
      emp::vector<int> args;
      emp::slice(cur_line, line_components, ',');
      for (const auto& arg_str : line_components) {
        args.emplace_back(emp::from_string<int>(arg_str));
      }
      programs.back().first.PushInst(
        inst_lib,
        inst_name,
        args,
        tags
      );
    }
  }
  return programs;
}


/// Load linear functions program from JSON format
void LoadLinearFunctionsProgram_JSONFormat(
  std::istream& input
) {
  // TODO - implement JSON format load
  emp_assert(false);
}

}