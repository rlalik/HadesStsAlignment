#include "lib.hpp"

auto main() -> int
{
  auto const lib = library {};

  return lib.name == "HadesStsAlignment" ? 0 : 1;
}
