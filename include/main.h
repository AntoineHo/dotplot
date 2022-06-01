#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

// LOCAL MODULES
#include "svg.h"
using namespace svg;
#include "classes.h"

// Includes
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
using namespace boost;
namespace po = boost::program_options;

typedef std::pair<std::string, std::string> stringpair;
typedef std::pair<uint,uint> intpair;
typedef std::pair<double,double> doublepair;
typedef std::pair<std::string,uint> intstringpair;
typedef std::pair<std::string,double> doublestringpair;


//#include <map> // Necessary to use maps
//#include <string> // Necessary to use strings
//#include <iostream> // Preprocessor indication > read/write from/to console
//#include <vector> // Necessary to use dynamically allocated tables
//#include <sstream> // Necessary to use std::to_string
//#include <algorithm>    // std::sort
//#include <cstring> // Necessary to use c_string functions ex: strtok
//#include <cstdlib> // Necessary to use strtol
//#include <unistd.h> // Necessary for page_size
//#include <cmath>

//#include "classes.h"
//

// Defines
#define MAX_STRING_SIZE 1024 // Defined string size for reading .paf file


/**
 * \brief Function that reads command line and stores arguments
 * \param argc : int
 * \param argv : pointer to char array
 * \param vm : address of po::variables_map
 * \return int (1 fail, 0 pass)
 **/
int parse_arguments(int argc, char *argv[], po::variables_map &vm);


/**
 * \brief Function that reads a .paf file and stores alignments
 * \param file_name : std::string
 * \param alignments : address of a vector of alignements
 * \param queries : address of a vector of strings
 * \param targets : address of a vector of strings
 * \param lengths : address of a map <std::string, int>
 * \param vm : address of po::variables_map
 * \return int (1 fail, 0 pass)
 **/
int parse_paf(
  std::string filename,
  std::vector<Alignment> &alignments,
  std::set<std::string> &queries,
  std::set<std::string> &targets,
  std::map<std::string, uint> &lengths,
  po::variables_map &vm
);


/**
 * \brief Function that creates a .svg file using .paf alignments
 * \param alignments : pointer to a vector of alignements
 * \param queries : reference of a vector containing queries
 * \param targets : reference of a vector containing targets
 * \param vm : reference of a variables_map containing parameters
 * \param vm : reference of a variables_map containing parameters
 * \return int (1 fail, 0 pass)
 **/
int make_svg(
  std::vector<Alignment> &alignments,
  std::list<std::string> &queries,
  std::list<std::string> &targets,
  std::map<std::string,uint> &lengths,
  po::variables_map &vm
);

/**
 * \brief Function that checks whether we should keep an alignment
 * \params ql,tl,mq,bl : integers
 * \param identity : float
 * \param vm : address of a variable map containing parameters
 * \return bool (true pass, false fail)
 **/
bool check_alignment(int ql, int tl, int mq, int bl, float identity, po::variables_map &vm);

/**
 * \brief Function that reads command line and stores arguments
 * \param &queries : reference of vector of string
 * \param &targets : reference of vector of string
 * \param &lengths : reference of a map string:int
 * \param &doc : reference to a SVG document
 * \param factor_x : double
 * \param factor_y : double
 * \param x_dim : uint
 * \param y_dim : uint
 * \param x_size : uint
 * \param y_size : uint
 * \param offset : uint
 * \param fontsize : int
 * \returns map of coordinates of axis origin for each query target pair
 **/
std::map<stringpair, doublepair> draw_axis(
  std::list<std::string> &queries,
  std::list<std::string> &targets,
  std::map<std::string, uint> &lengths,
  Document &doc,
  double factor_x,
  double factor_y,
  uint x_dim,
  uint y_dim,
  uint x_size,
  uint y_size,
  uint offset,
  int fontsize
);

/**
 * \brief Function that reads command line and stores arguments
 * \param &alignments : reference of vector of alignments
 * \param &lengths : reference of a map string:int
 * \param &coords : reference of map of stringpair to intpair
 * \param &doc : reference to a SVG document
 * \param factor_x : double
 * \param factor_y : double
 * \param sw : float (stroke width)
 * \returns map of coordinates of axis origin for each query target pair
 **/
void draw_alignments(
  std::vector<Alignment> &alignments,
  std::map<std::string,uint> &lengths,
  std::map<stringpair,doublepair> &coords,
  Document &doc,
  double factor_x,
  double factor_y,
  float sw
);

/**
 * \brief Function that returns a ceiling division (always round up)
 * \param i : int
 * \param r : int
 * \returns int
 **/
int ceildiv(uint i, uint r);

#endif
