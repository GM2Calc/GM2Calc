#pragma once
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

const char PATH_SEPARATOR =
#ifdef _WIN32
   '\\';
#else
   '/';
#endif

namespace gm2calc {
namespace test {

/**
 * Reads real numbers from a file line by line into a vector of vectors of T.
 *
 * @param filename file name
 * @tparam T data type
 *
 * @return vector of vectors of doubles.
 */
template <typename T>
std::vector<std::vector<T>>
read_from_file(const std::string& filename)
{
   std::vector<std::vector<T>> data;
   std::string line;
   std::ifstream fstr(filename);

   if (!fstr.is_open()) {
      throw std::runtime_error("Cannot open file: " + filename);
   }

   while (std::getline(fstr, line)) {
      std::istringstream iss(line);
      std::vector<T> vec{std::istream_iterator<T>(iss),
                         std::istream_iterator<T>()};
      data.emplace_back(vec);
   }

   return data;
}

} // namespace test
} // namespace gm2calc
