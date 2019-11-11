#include "qif"

namespace qif::lp {

bool   Defaults::presolve  = true;
string Defaults::msg_level = MsgLevel::OFF;
string Defaults::method    = Method::AUTO;
string Defaults::solver    = Solver::AUTO;

} // namespace qif::lp
