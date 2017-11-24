#ifndef TIME_SUPPORT_H
#define TIME_SUPPORT_H

#include <chrono>
#include <string>
#include <sstream>

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief days represents the time in number of days.
 */
typedef std::chrono::duration<int, std::ratio<86400,1>> days;

/**
 * @brief The TimeElements struct represents a time as a set of standard elements, like seconds.
 */
struct TimeElements {
  days                      d;    ///< Days.
  std::chrono::hours        h;    ///< Hours.
  std::chrono::minutes      m;    ///< Minutes.
  std::chrono::seconds      s;    ///< Seconds.
  std::chrono::milliseconds ms;   ///< Milliseconds.
  std::chrono::microseconds us;   ///< Microseconds.
  std::chrono::nanoseconds  ns;   ///< Nanoseconds.
};

class TimeSupport {
public:
  /**
   * @brief SplitTimeElements extracts elements from time.
   * @param duration The time to split.
   * @param timeElements Elements to be extracted.
   */
  template < typename Rep, typename Period >
  static void SplitTimeElements(std::chrono::duration<Rep,Period> duration, TimeElements &timeElements) {
    timeElements.d = std::chrono::duration_cast<days>(duration);
    duration -= timeElements.d;
    timeElements.h = std::chrono::duration_cast<std::chrono::hours>(duration);
    duration -= timeElements.h;
    timeElements.m = std::chrono::duration_cast<std::chrono::minutes>(duration);
    duration -= timeElements.m;
    timeElements.s = std::chrono::duration_cast<std::chrono::seconds>(duration);
    duration -= timeElements.s;
    timeElements.ms = std::chrono::duration_cast<std::chrono::milliseconds>(duration);
    duration -= timeElements.ms;
    timeElements.us = std::chrono::duration_cast<std::chrono::microseconds>(duration);
    duration -= timeElements.us;
    timeElements.ns = duration;
  }

  /**
   * @brief TimeElementsToString converts a time into a string representation.
   * @param timeElements The time elements to convert.
   * @param timeStr The output string.
   */
  static void TimeElementsToString(const TimeElements &timeElements, std::string &timeStr) {
    std::stringstream stream;
    if (timeElements.d.count() > 0)
      stream << timeElements.d.count() << "d.";
    if (timeElements.h.count() > 0 || timeElements.d.count() > 0)
      stream << std::setw(2) << std::setfill('0') << timeElements.h.count() << "h.";
    if (timeElements.m.count() > 0 || timeElements.h.count() > 0 || timeElements.d.count() > 0)
      stream << std::setw(2) << std::setfill('0') << timeElements.m.count() << "m.";
    if (timeElements.s.count() > 0 || timeElements.m.count() > 0 || timeElements.h.count() > 0 || timeElements.d.count() > 0)
      stream << std::setw(2) << std::setfill('0') << timeElements.s.count() << "s.";
    stream << std::setw(3) << std::setfill('0') << timeElements.ms.count() << "ms.";
    stream << std::setw(3) << std::setfill('0') << timeElements.us.count() << "us.";
    if (timeElements.ns.count() > 0)
      stream << std::setw(3) << std::setfill('0') << timeElements.ns.count() << "ns.";
    timeStr = stream.str();
  }

  /**
   * @brief TimeToString converts duration into a string representation.
   * @param duration The time to convert.
   * @param timeStr The output string.
   */
  template < typename Rep, typename Period >
  static void TimeToString(std::chrono::duration<Rep,Period> duration, std::string &timeStr) {
    TimeElements timeElements;
    SplitTimeElements(duration, timeElements);
    TimeElementsToString(timeElements, timeStr);
  }
};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // TIME_SUPPORT_H
