#include "main.h"

// FUNCTIONS
int main(int argc, char* argv[]) {

  int ret;
  po::variables_map vm;

  // parse arguments
  ret = parse_arguments(argc, argv, vm); // return 0 for proceed & 1 for help
  if (ret == 1) { return 0; } else if (ret > 1) { return 1; }

  // prepare to read and store data
  std::string filename = vm["paf"].as<std::string>();
  std::cout << "Input .paf file is: " << filename << std::endl;
  std::vector<Alignment> alignments;
  std::set<std::string> queries; // store queries names
  std::set<std::string> targets; // store targets names
  std::map<std::string, uint> lengths; // map of lengths of contigs

  // read .paf, count queries and store alignements info
  std::cout << "Reading .paf file..." << std::endl;
  ret = parse_paf(filename, alignments, queries, targets, lengths, vm);
  if (ret == 1) {
    std::cout << "ERROR: Could not read .paf file." << std::endl;
    return 1;
  }

  /*
  std::cout << "Printing lengths" << std::endl;
  for (auto const &pair: lengths) {
    std::cout << '{' << pair.first << "," << pair.second << '}' << std::endl;
  }
  */

  // Sort contigs by lengths
  /*
  std::map<std::string,int> sorted_lengths = sort_lengths(lengths);
  std::cout << "Printing sorted lengths" << std::endl;
  for (auto const &pair: sorted_lengths) {
    std::cout << '{' << pair.first << "," << pair.second << '}' << std::endl;
  }
  */

  // Sort queries and targets by lengths
  std::list<std::string> sorted_contigs = sort_keys_by_values(lengths);
  std::list<std::string> sorted_queries = sort_list_by_list(
    queries, sorted_contigs, true
  );
  std::list<std::string> sorted_targets = sort_list_by_list(
    targets, sorted_contigs, true
  );

  //make_axis()
  /*
  std::cout << "Printing sorted queries" << std::endl;
  std::list<std::string>::iterator it;
  for (it = sorted_queries.begin(); it != sorted_queries.end(); it++) {
    std::cout << *it << std::endl;
  }
  */

  //std::cout << "---" << std::endl;
  /* // Example of reversed for loop through map<std::string, int>
  std::map<std::string,int>::reverse_iterator rit;
  for (rit = sorted_lengths.rbegin(); rit != sorted_lengths.rend(); rit++) {
    std::cout << rit->first << " - " << rit->second << std::endl;
  }
  */

  /* // Example for loop through map<std::string, int>
  std::map<std::string,int>::iterator it;
  for (it = lengths.begin(); it != lengths.end(); it++) {
    std::cout << it->first << " - " << it->second << std::endl;
  }
  */

   /* // DEBUG
  for (int i(0); i < alignments.size(); i++) {
    alignments[i].print();
  }
  */

  std::cout << "Making .svg" << std::endl;

  ret = make_svg(alignments, sorted_queries, sorted_targets, lengths, vm);
  if (ret == 1) {
    std::cout << "ERROR: Could not read .paf file." << std::endl;
    return 1;
  }

  return 0;
}





int parse_arguments(int argc, char* argv[], po::variables_map &vm) {
  // Declare the supported options.
  po::options_description desc("Allowed options");
  try {
    std::string default_output = "output";
    desc.add_options()
        ("help,h", "Prints help message")
        ("paf,p", po::value<std::string>(), "Input file (.paf format)")
        ("output,o", po::value<std::string>()->default_value(default_output), "prefix of output")
        ("qlen,q", po::value<int>()->default_value(100000), "Minimum query length")
        ("tlen,t", po::value<int>()->default_value(100000), "Minimum target length")
        ("alen,a", po::value<int>()->default_value(10000), "Minimum alignment length")
        ("mapq,m", po::value<int>()->default_value(0), "Minimum mapping quality")
        ("stroke,s", po::value<float>()->default_value(1.0), "Stroke width")
        ("idy,i", po::value<float>()->default_value(0.31), "Minimum alignment identity")
        ("xsize,x", po::value<int>()->default_value(4096), "x-size. Default")
        ("ysize,y", po::value<int>()->default_value(4096), "y-size. Default")
        ("fontsize,f", po::value<int>()->default_value(12), "font size")
        ("offset", po::value<int>()->default_value(80), "offset of origin")
    ;

    // no argument in command line
    if (argc == 1) {
        std::cout << desc << std::endl;
        return 1;
    }

    // store command line
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // if help
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    // error catch
  } catch ( const std::exception& e ) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  return 0;
}







int parse_paf(
  std::string filename,
  std::vector<Alignment> &alignments,
  std::set<std::string> &queries,
  std::set<std::string> &targets,
  std::map<std::string, uint> &lengths,
  po::variables_map &vm
) {

  std::ifstream paf (filename);
  std::vector<std::string> s;
  std::string line;

  std::string qn;
  uint qs;
  uint qe;
  uint ql;
  std::string tn;
  uint ts;
  uint te;
  uint tl;
  uint bl;
  uint nm;
  uint mq;
  bool sd;

  int id = 0;
  if (paf.is_open()) {
    while ( getline(paf, line) ) {

      boost::split(s, line, boost::is_any_of("\t"));

      // Fill Alignments
      qn = s[0]; tn = s[5];
      ql = std::stoi(s[1]); qs = std::stoi(s[2]); qe = std::stoi(s[3]);
      tl = std::stoi(s[6]); ts = std::stoi(s[7]); te = std::stoi(s[8]);
      nm = std::stoi(s[9]); bl = std::stoi(s[10]); mq = std::stoi(s[11]);
      if (s[4] == "+") { sd = true; } else { sd = false; }

      float identity = (float)nm/(float)bl;
      //std::cout << "Identity is " << identity << std::endl;
      if (check_alignment(ql, tl, mq, bl, identity, vm)) {

        // Add query and target to sets
        queries.insert(qn);
        targets.insert(tn);

        // Add length of query and target if not already in map
        // count in a map can only return 1
        if ( lengths.count(qn) == 0 ) { lengths[qn] = ql; }
        if ( lengths.count(tn) == 0 ) { lengths[tn] = tl; }

        // Add alignment
        alignments.push_back(
          Alignment(id, qn, ql, qs, qe, sd, tn, tl, ts, te, nm, bl, mq)
        );
        id++;
      }
    }
  }

  paf.close();

  return 0;

}

bool check_alignment(int ql, int tl, int mq, int bl, float identity, po::variables_map &vm) {
  if (ql < vm["q"].as<int>()) { return false; }
  if (tl < vm["t"].as<int>()) { return false; }
  if (mq < vm["m"].as<int>()) { return false; }
  if (bl < vm["a"].as<int>()) { return false; }
  if (identity < vm["i"].as<float>()) { return false; }
  return true;
}







int make_svg(
  std::vector<Alignment> &alignments,
  std::list<std::string> &queries,
  std::list<std::string> &targets,
  std::map<std::string,uint> &lengths,
  //std::map<stringpair,intpair> &coords,
  po::variables_map &vm
){

  // Create output file
  //std::string filename = vm["p"].as<std::string>();
  //std::string ext = ".dotplot.svg";
  std::string output = vm["o"].as<std::string>() + ".dotplot.svg";
  std::cout << "Output is: " << output << std::endl;

  /* previous way of doing it
  // The reduction factor is the size of a minimum alignment which should be a "pixel sized"
  // If the minimum alignment size is 10Kb then a pixel should be 10Kb
  // We will set the minimum size of an axis to 4 pixels so here 40Kb
  int reduction_factor = vm["r"].as<int>();
  */

  double x_size = (double) vm["x"].as<int>(); // + offset;
  double y_size = (double) vm["y"].as<int>(); // + offset;

  uint x_dim = 0;
  uint y_dim = 0;
  std::list<std::string>::iterator it;
  for (it = queries.begin(); it != queries.end(); it++ ) {
    x_dim += lengths[*it];
  }
  for (it = targets.begin(); it != targets.end(); it++ ) {
    y_dim += lengths[*it];
  }

  /* previous way
  // The condition (x_dim % reduction_factor != 0) returns 1 if different from 0
  uint adjusted_x = x_dim / reduction_factor + (x_dim % reduction_factor != 0);
  uint adjusted_y = y_dim / reduction_factor + (y_dim % reduction_factor != 0);
  */

  double factor_x = x_size / x_dim; // fraction used to adjust x coordinates
  double factor_y = y_size / y_dim; // fraction used to adjust y coordinates

  int offset = vm["of"].as<int>();
  x_size += 2*offset;
  y_size += 2*offset;

  /* //DEBUG
  //std::cout << "factor_x = x_size / x_dim" << std::endl;
  //std::cout << factor_x << " = " << x_size <<" / " << x_dim << std::endl;
  //std::cout << "factor_y = y_size / y_dim" << std::endl;
  //std::cout << factor_y << " = " << y_size <<" / " << y_dim << std::endl;
  */

  // Create a layout
  //std::cout << "Dimensions: (" << x_dim << ", " << y_dim << ")" << std::endl;
  Dimensions dimensions(x_size, y_size);
  Document doc(output, Layout(dimensions, Layout::BottomLeft));

  // Create axis
  int fontsize = vm["f"].as<int>();
  std::map<stringpair,doublepair> coords = draw_axis(
    queries, targets, lengths, doc, factor_x, factor_y,
    x_dim, y_dim, x_size, y_size, offset, fontsize
  );

  // Draw alignments
  float sw = vm["s"].as<float>();
  draw_alignments(alignments, lengths, coords, doc, factor_x, factor_y, sw);


  doc.save();

  return 0;

}








std::map<stringpair,doublepair> draw_axis(
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
){

  std::map<stringpair,doublepair> coords;
  std::vector<double> x_starts;
  std::vector<double> y_starts;
  std::vector<std::string> y_names;
  std::vector<std::string> x_names;
  bool first_read = true;

  //std::stringstream ss;
  double cumulative_x = 0.0 + offset;  // starts with a 20 offset (origin is at (20, 20))


  // Iterate through queries x targets to pass through each pair once
  std::list<std::string>::iterator xit;
  std::list<std::string>::iterator yit;

  std::string query;
  std::string target;
  for (xit = queries.begin(); xit != queries.end(); xit++) {

    // Compute x position of the pair (0,0) anchor at south west
    query = *xit;
    uint x_length = lengths[query];
    double adjusted_x = factor_x * x_length;
    x_starts.push_back(cumulative_x);
    x_names.push_back(query);


    double cumulative_y = 0.0 + offset;  // starts with a 20 offset (origin is at (20, 20))
    for (yit = targets.begin(); yit != targets.end(); yit++) {

      // Compute y position of the pair (0,0) anchor at south west
      target = *yit;
      uint y_length = lengths[target];
      double adjusted_y = factor_y * y_length;

      // Get adjusted coordinates of origin for the pair query-target
      stringpair qry_tgt = {query, target};
      doublepair orig = {cumulative_x, cumulative_y};
      coords[qry_tgt] = orig;

      if (first_read) {
        y_starts.push_back(cumulative_y);
        y_names.push_back(target);
      }

      // Increments y-axis origin of the pair at the end of previous contig
      cumulative_y += adjusted_y;
    }

    // If we read the targets once already do not push this anymore
    if (first_read) {
      first_read = false;
      // last pair // I think next two lines are an error
      //y_starts.push_back(cumulative_y);
      //y_names.push_back(target);
    }

    // Increments x-axis origin of the pair at the end of previous contig
    cumulative_x += adjusted_x;
  }

  // Last pair // I think next two lines are an error
  //x_starts.push_back(cumulative_x);
  //x_names.push_back(query);

  // draw axis lines & contig names
  std::vector<double>::iterator it;
  uint index = 0;
  uint cumulative_length = 0;
  for (it = x_starts.begin(); it != x_starts.end(); it++) {

    // Draw vertical line
    double xs = *it; // xs is already offset by +20.0
    Point origin = Point(xs, offset - 5.0); // starts an (origin is at (offset, offset))
    Point xend = Point(xs, y_size - offset); // - 2*offset because we added 2 offsets
    Line x_axis(origin, xend, Stroke(0.5, Color::Black));
    doc << x_axis;

    int text_offset = 25.0;
    //if (index % 2 != 0) { text_offset += 10.0; }

    // Print name of contig below x axis
    std::string query = x_names[index];
    double adjusted_x_length = (double)lengths[query] * factor_x;
    double center_x_text = xs + adjusted_x_length/2.0;
    double center_y_text = offset - text_offset;
    double rotation_value = 90.0;
    origin = Point(center_x_text, center_y_text);
    Text query_name_text(
      origin, query, Fill(Color::Black), Font((double)fontsize),
      rotation_value, center_x_text, center_y_text, Stroke(0.0)
    ); // 0.0 is rotation, center_x_text is x coord and 5.0 is y coord for center of text
    doc << query_name_text;

    // Add a tick line at the end showing total length span
    // with a tick label
    cumulative_length += lengths[query];
    std::stringstream tick_val_ss;
    tick_val_ss << std::floor( ((double)cumulative_length / 1000000)*100.0 ) / 100.0 << "M";
    std::string tick_val_string = tick_val_ss.str();

    int tick_offset = 25.0;
    //if (index % 2 == 0) { tick_offset += 10.0; }

    center_y_text = offset - tick_offset;
    center_x_text = xs+adjusted_x_length;
    origin = Point(center_x_text, center_y_text);
    rotation_value = 90.0;
    Text tick_value(
      origin, tick_val_string, Fill(Color::Black), Font((double)fontsize),
      rotation_value, center_x_text, center_y_text, Stroke(0.0)
    );
    doc << tick_value;

    index += 1;
  }

  index = 0;
  cumulative_length = 0;
  for (it = y_starts.begin(); it != y_starts.end(); it++) {
    double ys = *it; // ys is already offset
    Point origin = Point(offset - 5.0, ys);
    Point yend = Point(x_size-offset, ys);  // -offset because we added 2 offsets
    Line y_axis(origin, yend, Stroke(0.5, Color::Black));
    doc << y_axis;

    int text_offset = 25.0;
    //if (index % 2 != 0) { text_offset += 10.0; }

    // Print name of contig on left of y axis
    std::string target = y_names[index];
    double adjusted_y_length = (double)lengths[target] * factor_y;
    double center_y_text = ys + adjusted_y_length/2.0;
    double center_x_text = offset - text_offset;
    double rotation_value = 0.0;
    origin = Point(center_x_text, center_y_text);  // should be at x = 10.0
    Text target_name_text(
      origin, target, Fill(Color::Black), Font((double)fontsize),
      rotation_value, center_x_text, center_y_text, Stroke(0.0)
    );  // 90.0 is rotation, 10.0 is x coord and center_y_text is y coord for center of text
    doc << target_name_text;

    // Add a tick line at the end showing total length span
    // with a tick label
    cumulative_length += lengths[target];
    std::stringstream tick_val_ss;
    tick_val_ss << std::floor( ((double)cumulative_length / 1000000)*100.0 ) / 100.0 << "M";
    std::string tick_val_string = tick_val_ss.str();

    int tick_offset = 40.0;
    center_x_text = offset - tick_offset;
    center_y_text = ys+adjusted_y_length;
    origin = Point(center_x_text, center_y_text);  // here text needs to be more offset or goes into the axis
    rotation_value = 0.0;
    Text tick_value(
      origin, tick_val_string, Fill(Color::Black), Font((double)fontsize),
      rotation_value, center_x_text, center_y_text, Stroke(0.0)
    );
    doc << tick_value;

    index += 1;
  }

  // Add a 0.0 tick name
  Text origin_text(
    Point(offset - 15.0, offset - 15.0), "0.00M", Fill(Color::Black), Font((double)fontsize),
    0.0, offset - 15.0, offset - 15.0, Stroke(0.0)
  );  // 90.0 is rotation, 10.0 is x coord and center_y_text is y coord for center of text
  doc << origin_text;

  return coords;
}


void draw_alignments(
  std::vector<Alignment> &alignments,
  std::map<std::string,uint> &lengths,
  std::map<stringpair,doublepair> &coords,
  Document &doc,
  double factor_x,
  double factor_y,
  float sw
){

  // For each alignment
  for (auto &aln : alignments) {
    stringpair qry_tgt = {aln.get_qn(), aln.get_tn()};
    doublepair origin = coords[qry_tgt];
    //std::cout << "Alignement on axis: (" << origin.first << "," << origin.second << ")" << std::endl;

    double xs = 0.0; double xe = 0.0;
    double ys = 0.0; double ye = 0.0;

    bool sd = aln.get_sd();
    uint qs = aln.get_qs(); uint qe = aln.get_qe();
    uint ts = 0; uint te = 0;

    if (sd) {
    	//uint ql = aln.get_ql(); uint tl = aln.get_tl();
    	ts = aln.get_ts(); te = aln.get_te();
    } else {
    	ts = aln.get_te(); te = aln.get_ts();
    }

    double identity = (double) aln.get_nm() / aln.get_bl();

    double adjusted_qs = factor_x * qs;
    double adjusted_qe = factor_x * qe;
    double adjusted_ts = factor_y * ts;
    double adjusted_te = factor_y * te;

    xs = origin.first + adjusted_qs;
    xe = origin.first + adjusted_qe;
    ys = origin.second + adjusted_ts;
    ye = origin.second + adjusted_te;

    Point ori = Point(xs, ys);
    Point end = Point(xe, ye);

    if (identity < 0.5) {
      doc << Line (ori, end, Stroke(sw, Color::Red));
    } else if (identity < 0.7) {
      doc << Line (ori, end, Stroke(sw, Color::Blue));
    } else if (identity < 0.9) {
      doc << Line (ori, end, Stroke(sw, Color::Green));
    } else {
      doc << Line (ori, end, Stroke(sw, Color::Black));
    }
  }

}

// UNUSED
int ceildiv(uint i, uint r) {
  return i / r + (i % r != 0);
}
