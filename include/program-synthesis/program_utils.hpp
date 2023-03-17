#pragma once

#include <iostream>

#include "emp/base/vector.hpp"
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

/// Load linear functions program from print format
void LoadLinearFunctionsProgram_PrintFormat(
  std::istream& input
) {
  // TODO
}

/// Load linear functions program from JSON format
void LoadLinearFunctionsProgram_JSONFormat(
  std::istream& input
) {
  // TODO
}

}