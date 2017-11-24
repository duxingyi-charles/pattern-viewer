#ifndef STRING_SUPPORT
#define STRING_SUPPORT

#include <string>
#include <functional>

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief The Hasher class computes the hash value of a string.
 */
class Hasher : public std::hash<std::string> {
public:
  typedef std::hash<std::string>::result_type   HashValue;

  /**
   * @brief StringToHashValue converts a string into a single integer value.
   * This is extremely useful for avoiding to occupy too much memory.
   * @param str The input string to convert.
   * @return The corresponding unique value associated to the string str.
   */
  static HashValue StringToHashValue(const std::string &str) {
    Hasher hasher;
    return hasher(str);
  }
};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // STRING_SUPPORT

