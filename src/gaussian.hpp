#ifndef _gaussian_hpp_INCLUDED
#define _gaussian_hpp_INCLUDED

#include <vector>
#include <optional>

namespace CaDiCaL {

class Internal; // forward declaration

// Simple Gauss-Jordan elimination over GF(2) for XOR clauses.
// This code is derived in spirit from CryptoMiniSat's Gaussian
// elimination implementation (MIT licensed).
// Simplified Gaussian elimination solver for XOR clauses. Inspired by
// CryptoMiniSat's implementation but stripped down for CaDiCaL.
class Gaussian {
public:
  struct Xor {
    std::vector<int> vars; // variable indices starting at 1
    bool rhs = false;      // right-hand side
  };

  void add_clause(const Xor &clause);
  // Solve the system. Returns empty optional if unsatisfiable.
  std::optional<std::vector<bool>> eliminate() const;

  // Propagate given the current assignments in 'internal'.
  // Returns true if a conflict was detected.
  bool propagate(Internal *internal);

  bool empty() const { return equations.empty(); }
  void clear() { equations.clear(); }

private:
  std::vector<Xor> equations;
};

} // namespace CaDiCaL

#endif // _gaussian_hpp_INCLUDED
