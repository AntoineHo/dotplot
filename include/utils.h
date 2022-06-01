#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <iostream>
#include <iterator>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <sstream>
#include <fstream> // Necessary to read/write from/to files
#include <list>
#include <cmath>

/**
 * \brief Function that reads command line and stores arguments
 * \param &set_to_sort : reference of set of string
 * \param &sorted_lengths : reference of a map string:int
 * \returns a map<string, int>
 **/
std::list<std::string> sort_list_by_list(
  std::set<std::string> &set_to_sort,
  std::list<std::string> &sorted_list,
  bool reverse
);

/**
 * \brief Function that sorts keys from a map by value
 * \param &map : reference of a map string:int
 * \returns a list of sorted keys
 **/
std::list<std::string> sort_keys_by_values(
  std::map<std::string,uint> &map
);

#endif
