#pragma once

#include <utility>

namespace psynth {

/// A TestResult records a program's problem-specific "score" on a test case.
struct TestResult {
  bool has_response=false;            ///< Did the program produce a response?
  bool is_correct=false;              ///< Is output correct?
  double score = 0.0;                 ///< Typically value between 0 and 1 indicating "closeness" to similarity (used for partial credit if incorrect).

  TestResult(
    bool _responded=false,
    bool _correct=false,
    double _score=0.0
  ) :
    has_response(_responded),
    is_correct(_correct),
    score(_score)
  { ; }

  TestResult& operator=(const TestResult&) = default;
  TestResult& operator=(TestResult&&) = default;
};

}