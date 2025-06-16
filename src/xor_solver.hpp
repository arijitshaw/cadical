#ifndef _xor_solver_hpp_INCLUDED
#define _xor_solver_hpp_INCLUDED

#include <vector>
#include <optional>

namespace CaDiCaL {

// Simple Gauss-Jordan elimination over GF(2) for XOR clauses.
// This code is derived in spirit from CryptoMiniSat's Gaussian
// elimination implementation (MIT licensed).
class XORSolver {
public:
  struct Equation {
    std::vector<int> vars; // variable indices starting at 1
    bool rhs = false;      // right-hand side
  };

  void add_equation(const Equation &eq);
  // Solve the system. Returns empty optional if no solution exists.
  std::optional<std::vector<bool>> solve() const;

  bool empty() const { return equations.empty(); }
  void clear() { equations.clear(); }

private:
  std::vector<Equation> equations;
};

} // namespace CaDiCaL

#endif // _xor_solver_hpp_INCLUDED
