#ifndef CLASSES_H_INCLUDED
#define CLASSES_H_INCLUDED

#include "utils.h"

// DECLARATION OF DEFINES
#define MAX_STRING_SIZE 1024 // Defined string size for reading .paf file

//CLASSES
class Alignment {
public: // constructor
  Alignment(uint pid, std::string pqn, uint pql, uint pqs,
    uint pqe, bool psd, std::string ptn, uint ptl,
    uint pts, uint pte, uint pnm, uint pbl, uint pmq);
  Alignment();
  void print();

  // Parameters
  uint get_id(); // index
  std::string get_qn(); // query name
  uint get_qs(); // query start
  uint get_qe(); // query end
  uint get_ql(); // query length
  std::string get_tn();
  uint get_ts(); // target start
  uint get_te(); // target end
  uint get_tl(); // target length
  bool get_sd(); //strand + (true) / - (false)
  uint get_mq(); // mapq
  uint get_nm(); // num matches
  uint get_bl(); // aligned length
  float get_identity(); // identity

  // FOR SETTING
  //void set_parameters(std::string p_name, int p_alignedLength, int p_startPos,
  //  int p_endPos, int p_chrStartPos, int p_chrEndPos, bool p_relativeStrand,
  //  int p_mappingQuality);
  void set_id(uint pid);
  void set_qn(std::string pqn);
  void set_qs(uint pqs);
  void set_qe(uint pqe);
  void set_ql(uint pql);
  void set_tn(std::string ptn);
  void set_ts(uint pts);
  void set_te(uint pte);
  void set_tl(uint ptl);
  void set_bl(uint pbl);
  void set_nm(uint pnm);
  void set_mq(uint pmq);
  void set_sd(bool psd);

private:
  int id;
  std::string qn;
  int qs;
  int qe;
  int ql;
  std::string tn;
  int ts;
  int te;
  int tl;
  int bl;
  int nm;
  int mq;
  bool sd;
  float identity;
};

/*
class Axis {
public: // constructor
  Axis(int pid, int pql, int pqs);

  // Parameters
  int get_id(); // index
  int get_x();
  int get_y();
  int* get_coords(); // index

  // Set parameters
  void set_id(int pid);
  void set_x(int px);
  void set_y(int py);

  void draw();

private:
  int id;
  int x;
  int y;
};
*/


#endif
