#include "classes.h"

Alignment::Alignment(uint pid, std::string pqn, uint pql, uint pqs,
  uint pqe, bool psd, std::string ptn, uint ptl, uint pts,
  uint pte, uint pnm, uint pbl, uint pmq) : id(pid), qn(pqn), ql(pql),
  qs(pqs), qe(pqe), sd(psd), tn(ptn), tl(ptl), ts(pts),
  te(pte), nm(pnm), bl(pbl), mq(pmq)
{
  float identity = (float)pnm / (float)pbl;
}

uint Alignment::get_id() { return id; }
std::string Alignment::get_qn() { return qn; }
uint Alignment::get_ql() { return ql; }
uint Alignment::get_qs() { return qs; }
uint Alignment::get_qe() { return qe; }
std::string Alignment::get_tn() { return tn; }
uint Alignment::get_tl() { return tl; }
uint Alignment::get_ts() { return ts; }
uint Alignment::get_te() { return te; }
bool Alignment::get_sd() { return sd; }
uint Alignment::get_bl() { return bl; }
uint Alignment::get_mq() { return mq; }
uint Alignment::get_nm() { return nm; }
float Alignment::get_identity() { return identity; }

void Alignment::set_id(uint pid) { id = pid; }
void Alignment::set_qn(std::string pqn) { qn = pqn; }
void Alignment::set_ql(uint pql) { ql = pql; }
void Alignment::set_qs(uint pqs) { qs = pqs; }
void Alignment::set_qe(uint pqe) { qe = pqe; }
void Alignment::set_tn(std::string ptn) { tn = ptn; }
void Alignment::set_tl(uint ptl) { tl = ptl; }
void Alignment::set_ts(uint pts) { ts = pts; }
void Alignment::set_te(uint pte) { te = pte; }
void Alignment::set_sd(bool psd) { sd = psd; }
void Alignment::set_bl(uint pbl) { bl = pbl; }
void Alignment::set_mq(uint pmq) { mq = pmq; }
void Alignment::set_nm(uint pnm) { nm = pnm; }

void Alignment::print() {
  std::cout << "Alignment_" << id << "(" << std::endl;
  std::cout << "Query = " << qn << std::endl;
  std::cout << "Query length = " << ql << std::endl;
  std::cout << "Query start = " << qs << std::endl;
  std::cout << "Query end = " << qe << std::endl;
  if (sd == true) {
    std::cout << "Strand = +" << std::endl;
  } else {
    std::cout << "Strand = -" << std::endl;
  }
  std::cout << "Target = " << tn << std::endl;
  std::cout << "Target length = " << tl << std::endl;
  std::cout << "Target start = " << ts << std::endl;
  std::cout << "Target end = " << te << std::endl;
  std::cout << "Number of matches = " << nm << std::endl;
  std::cout << "Block length = " << bl << std::endl;
  std::cout << "Mapping quality = " << mq << std::endl;
  std::cout << "Identity = " << identity << std::endl;
  std::cout << ")" << std::endl;
}

/*
Axis::Axis(int pid, int px, int py) : id(pid), x(px), y(py) {} // Each axis contains
int Axis::get_id() { return id; }
int Axis::get_x() { return x; }
int Axis::get_y() { return y; }

int* Axis::get_coords() {
  int *coords = new int[2];
  coords[0] = x; coords[1] = y;
  return coords;
}
void Axis::set_id(int pid) {id = pid;};
void Axis::set_x(int px) {x = px;};
void Axis::set_y(int py) {y = py;};

void Axis::draw() {
  // 1. Draw axis based on their coordinates

}
*/
