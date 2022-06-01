#include "utils.h"

typedef std::pair<std::string, int> pair;

std::list<std::string> sort_list_by_list(
  std::set<std::string> &set_to_sort,
  std::list<std::string> &sorted_list,
  bool reverse
){
  std::string contig;
  std::list<std::string> output;

  if (reverse) {
    std::list<std::string>::reverse_iterator it;
    for (it = sorted_list.rbegin(); it != sorted_list.rend(); it++) {
      contig = *it;
      if (set_to_sort.count(contig)) {
        output.push_back(contig);
      }
    }
  } else {
    std::list<std::string>::iterator it;
    for (it = sorted_list.begin(); it != sorted_list.end(); it++) {
      contig = *it;
      if (set_to_sort.count(contig)) {
        output.push_back(contig);
      }
    }
  }

  return output;
}

std::list<std::string> sort_keys_by_values(std::map<std::string,uint> &map) {

  // create an empty vector of pairs
  std::list<std::string> output;
  std::vector<pair> vec;
  // copy key-value pairs from the map to the vector
  std::copy(map.begin(), map.end(), std::back_inserter<std::vector<pair>>(vec));
  // sort the vector by increasing the order of its pair's second value
  // if the second value is equal, order by the pair's first value
  std::sort(vec.begin(), vec.end(), [](const pair &l, const pair &r) {
    if (l.second != r.second) {
      return l.second < r.second;
    }
  return l.first < r.first;
  });

  // print the vector
  for (auto const &pair: vec) {
    //std::cout << '{' << pair.first << "," << pair.second << '}' << std::endl;
    output.push_back(pair.first);
  }
  return output;
}
