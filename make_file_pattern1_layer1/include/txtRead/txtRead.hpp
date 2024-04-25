#pragma once
#include <cstddef>
#include <fstream>
#include <iterator>
#include <point_in/point.hpp>
#include <string_view>
#include <vec_in/CellVec.hpp>
#include <vec_in/CellVecAdd.hpp>
#include <vector>
/*
 * usage:
 * 1.
 *
 */
class txtRead {
public:
  void Open(std::string_view filename); // 打开文件
  cell get_master() { return master; }
  CellVec get_mental() { return mental; }
  CellVecAdd get_die() { return die; }
  void Read();
  void test();
  void Write(std::string_view filename);
  std::vector<std::size_t> Search(Point<2, double> po); // 以顶点搜索单元
  std::vector<std::size_t> Searchy(double up,
                                   double down); // 以纵坐标范围搜索单元
  std::vector<std::size_t> Searchx(double left,
                                   double right); // 以纵坐标范围搜索单元
  void Replace_layer(std::size_t la);             // 交换层
  void Replace_space(double w1, double w2, double edge, double edge_space);
  std::size_t Same_number(std::vector<std::size_t> a,
                          std::vector<std::size_t> b);
  void Replace_init();

private:
  // 用于Read()
  void ReadCell(cell &ce);
  void ReadMaster() { ReadCell(master); }
  void ReadMental();
  void ReadDie();
  void ReadCorner();
  // 用于Weite()
  void WriteCell(cell ce);
  void WriteCell(cellAdd ce);
  void WriteMaster();
  void WriteMental();
  void WriteDie();
  void WriteCorner();

private:
  std::ifstream ifs;
  cell master{};
  CellVec mental{};
  CellVecAdd die{};
  cell corner{};
  std::ofstream ofs;
};

// 类成员函数的实现
//
// 打开文件
void txtRead::Open(std::string_view filename) {
  ifs.open(filename.data());
  if (!ifs.is_open()) {
    throw std::runtime_error{"file not open"};
  }
}
// Read()系列
void txtRead::ReadCell(cell &ce) {
  std::string buffer;
  std::size_t n_point = 0;
  Point<2, double> point;
  ifs >> buffer;
  ce.add_name(buffer);
  ifs >> n_point;
  ce.add_number(n_point);
  for (auto i : std::ranges::iota_view{0uz, n_point}) {
    ifs >> point;
    ce += point;
  }
}
void txtRead::ReadMental() {
  std::size_t n_cell = 0;
  ifs >> n_cell;
  mental.set_number(n_cell);
  for (auto i : std::ranges::iota_view{0uz, n_cell}) {
    cell ce;
    ReadCell(ce);
    mental += ce;
  }
}
void txtRead::ReadDie() {
  std::size_t n_cell = 0;
  double value = 0;
  ifs >> n_cell;
  die.set_number(n_cell);
  for (auto i : std::ranges::iota_view{0uz, n_cell}) {
    cellAdd ce;
    ReadCell(ce);
    ifs >> value;
    ce.add_value(value);
    die += ce;
  }
}
void txtRead::ReadCorner() {
  double n = 0;
  Point<2, double> point;
  ifs >> n;
  corner.add_number(n);
  for (auto i : std::ranges::iota_view{0uz, n}) {
    ifs >> point;
    corner += point;
  }
}
void txtRead::Read() {
  ReadMaster();
  ReadMental();
  ReadDie();
  ReadCorner();
}

//
void txtRead::test() {
  ReadMaster();
  master.show();
  std::cout << std::endl;
  ReadMental();
  mental.show();
  ReadDie();
  die.show();
  ReadCorner();
  corner.show();
}
// Write（）系列
void txtRead::Write(std::string_view filename) {
  ofs.open(filename.data());
  WriteMaster();
  WriteMental();
  WriteDie();
  WriteCorner();
}
void txtRead::WriteCell(cell ce) {
  ofs << ce.get_name() << "\n";
  ofs << ce.get_number() << "\n";
  for (auto i : ce) {
    for (auto j : i) {
      ofs << j << " ";
    }
    ofs << "\n";
  }
  ofs << "\n";
}
void txtRead::WriteCell(cellAdd ce) {
  ofs << ce.get_name() << "\n";
  ofs << ce.get_number() << "\n";
  for (auto i : ce) {
    for (auto j : i) {
      ofs << j << " ";
    }
    ofs << "\n";
  }
  ofs << ce.get_value() << "\n\n";
}
void txtRead::WriteMaster() { WriteCell(master); }
void txtRead::WriteMental() {
  ofs << mental.get_number() << "\n";
  for (auto i : mental) {
    WriteCell(i);
  }
  ofs << "\n";
}
void txtRead::WriteDie() {
  ofs << die.get_number() << "\n";
  for (auto i : die) {
    WriteCell(i);
  }
}
void txtRead::WriteCorner() { WriteCell(corner); }

// Search()系列
std::vector<std::size_t> txtRead::Search(Point<2, double> po) {
  std::vector<std::size_t> number{};
  for (auto k = 0; auto i : die) {
    for (auto j : i) {
      if (po == j) {
        number.push_back(k);
      }
    }
    k++;
  }
  return number;
}
std::vector<std::size_t> txtRead::Searchy(double up, double down) {
  std::vector<std::size_t> number{};
  for (auto k = 0; auto i : die) {
    int number_of_point = 0;
    for (auto j : i) {
      // if (*(j.begin() + 1) <= up && *(j.begin() + 1) >= down) {
      if (j.data()[1] <= up + 0.001 && j.data()[1] >= down - 0.001) {
        number_of_point++;
      }
      // point j
      if (number_of_point == i.size()) {
        number.push_back(k);
      }
    }
    k++;
  }
  return number;
}
std::vector<std::size_t> txtRead::Searchx(double left, double right) {
  std::vector<std::size_t> number{};
  for (auto k = 0; auto i : die) {
    int number_of_point = 0;
    for (auto j : i) {
      if (j.data()[0] <= right + 0.001 && j.data()[0] >= left - 0.001) {
        number_of_point++;
      }
      // point j
      if (number_of_point == i.size()) {
        number.push_back(k);
      }
    }
    k++;
  }
  return number;
}

// Replace系列
void txtRead::Replace_layer(std::size_t la) {
  const double high_base = 0.4665;
  const double thick_layer = 0.187;
  for (auto &i : master) {
    i.data()[1] = i.data()[1] + (la - 2) * thick_layer;
  }

  for (auto &i : mental) {
    if (i.get_name() != "botleft" && i.get_name() != "botright") {
      for (auto &j : i) {
        j.data()[1] = j.data()[1] + (la - 2) * thick_layer;
      }
    }
  }

  std::vector<std::size_t> a, b;
  a = Searchy(high_base + thick_layer, high_base);
  b = Searchy(high_base + thick_layer + (la - 2) * thick_layer,
              high_base + (la - 2) * thick_layer);

  for (auto i : a) {
    for (auto &j : *(die.begin() + i)) {
      j.data()[1] = j.data()[1] + (la - 2) * thick_layer;
    }
  }

  for (auto i : b) {
    for (auto &j : *(die.begin() + i)) {
      j.data()[1] = j.data()[1] - (la - 2) * thick_layer;
    }
  }
}
void txtRead::Replace_space(double w1, double w2, double edge,
                            double edge_space) {
  /*
   * 目前写到水平方向扩展，刚在cell中加了pop——back函数，后面用于重新编写介质
   * */

  // 主导体w1
  master.data()[0].data()[0] -= (w1 - 0.045) / 2;
  master.data()[5].data()[0] -= (w1 - 0.045) / 2;
  master.data()[6].data()[0] -= (w1 - 0.045) / 2;
  master.data()[7].data()[0] -= (w1 - 0.045) / 2;
  master.data()[1].data()[0] += (w1 - 0.045) / 2;
  master.data()[2].data()[0] += (w1 - 0.045) / 2;
  master.data()[3].data()[0] += (w1 - 0.045) / 2;
  master.data()[4].data()[0] += (w1 - 0.045) / 2;
  // 右边导体
  mental.data()[3].data()[0].data()[0] += (edge - 0.0225 + (w1 - 0.045) / 2);
  mental.data()[3].data()[5].data()[0] += (edge - 0.0225 + (w1 - 0.045) / 2);
  mental.data()[3].data()[6].data()[0] += (edge - 0.0225 + (w1 - 0.045) / 2);
  mental.data()[3].data()[7].data()[0] += (edge - 0.0225 + (w1 - 0.045) / 2);
  mental.data()[3].data()[1].data()[0] +=
      (edge - 0.0225 + (w1 - 0.045) / 2 + (w2 - 0.045));
  mental.data()[3].data()[2].data()[0] +=
      (edge - 0.0225 + (w1 - 0.045) / 2 + (w2 - 0.045));
  mental.data()[3].data()[3].data()[0] +=
      (edge - 0.0225 + (w1 - 0.045) / 2 + (w2 - 0.045));
  mental.data()[3].data()[4].data()[0] +=
      (edge - 0.0225 + (w1 - 0.045) / 2 + (w2 - 0.045));
  // 左边导体
  mental.data()[2].data()[1].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045));
  mental.data()[2].data()[2].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045));
  mental.data()[2].data()[3].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045));
  mental.data()[2].data()[4].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045));
  mental.data()[2].data()[0].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045) + (w1 - 0.045));
  mental.data()[2].data()[5].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045) + (w1 - 0.045));
  mental.data()[2].data()[6].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045) + (w1 - 0.045));
  mental.data()[2].data()[7].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045) + (w1 - 0.045));
  // 最右导体
  mental.data()[4].data()[0].data()[0] +=
      (w1 - 0.045) / 2 + (w2 - 0.045) + (edge - 0.0225) + (edge_space - 0.045);
  mental.data()[4].data()[5].data()[0] +=
      (w1 - 0.045) / 2 + (w2 - 0.045) + (edge - 0.0225) + (edge_space - 0.045);
  mental.data()[4].data()[6].data()[0] +=
      (w1 - 0.045) / 2 + (w2 - 0.045) + (edge - 0.0225) + (edge_space - 0.045);
  mental.data()[4].data()[7].data()[0] +=
      (w1 - 0.045) / 2 + (w2 - 0.045) + (edge - 0.0225) + (edge_space - 0.045);
  mental.data()[4].data()[1].data()[0] += (w1 - 0.045) / 2 + 2 * (w2 - 0.045) +
                                          (edge - 0.0225) +
                                          (edge_space - 0.045);
  mental.data()[4].data()[2].data()[0] += (w1 - 0.045) / 2 + 2 * (w2 - 0.045) +
                                          (edge - 0.0225) +
                                          (edge_space - 0.045);
  mental.data()[4].data()[3].data()[0] += (w1 - 0.045) / 2 + 2 * (w2 - 0.045) +
                                          (edge - 0.0225) +
                                          (edge_space - 0.045);
  mental.data()[4].data()[4].data()[0] += (w1 - 0.045) / 2 + 2 * (w2 - 0.045) +
                                          (edge - 0.0225) +
                                          (edge_space - 0.045);

  // 介质
  // 待调试代码
  // 1
  auto die2 = die;
  std::size_t a;
  // 第一列
  // 1
  a = Same_number(Searchx(-9, -0.10792), Searchy(0.3955, 0.3705));
  die2[a][1] = mental[2][0];
  die2[a][2] = mental[2][7];
  // 2
  a = Same_number(Searchx(-9, -0.110305), Searchy(0.4065, 0.3955));
  die2[a][1] = mental[2][7];
  die2[a][2] = mental[2][6];
  // 3
  a = Same_number(Searchx(-9, -0.111355), Searchy(0.4665, 0.4065));
  die2[a][1] = mental[2][6];
  die2[a][2] = mental[2][5];

  // 第二列
  // 1
  a = Same_number(Searchx(-0.07208, -0.01792), Searchy(0.3955, 0.3705));
  die2[a][0] = mental[2][1];
  die2[a][3] = mental[2][2];
  die2[a][1] = master[0];
  die2[a][2] = master[7];
  // 2
  a = Same_number(Searchx(-0.069695, -0.020305), Searchy(0.4065, 0.3955));
  die2[a][0] = mental[2][2];
  die2[a][3] = mental[2][3];
  die2[a][1] = master[7];
  die2[a][2] = master[6];
  // 3
  a = Same_number(Searchx(-0.068645, -0.021355), Searchy(0.4665, 0.4065));
  die2[a][0] = mental[2][3];
  die2[a][3] = mental[2][4];
  die2[a][1] = master[6];
  die2[a][2] = master[5];

  // 第三列
  // 1
  a = Same_number(Searchx(0.01792, 0.04958), Searchy(0.3955, 0.3705));
  die2[a][0] = master[1];
  die2[a][3] = master[2];
  die2[a][1] = mental[3][0];
  die2[a][2] = mental[3][7];
  // 2
  a = Same_number(Searchx(0.020305, 0.047195), Searchy(0.4065, 0.3955));
  die2[a][0] = master[2];
  die2[a][3] = master[3];
  die2[a][1] = mental[3][7];
  die2[a][2] = mental[3][6];
  // 3
  a = Same_number(Searchx(0.021355, 0.046145), Searchy(0.4665, 0.4065));
  die2[a][0] = master[3];
  die2[a][3] = master[4];
  die2[a][1] = mental[3][6];
  die2[a][2] = mental[3][5];

  // 第四列
  // 1
  a = Same_number(Searchx(0.08542, 0.13958), Searchy(0.3955, 0.3705));
  die2[a][0] = mental[3][1];
  die2[a][3] = mental[3][2];
  die2[a][1] = mental[4][0];
  die2[a][2] = mental[4][7];
  // 2
  a = Same_number(Searchx(0.087805, 0.137195), Searchy(0.4065, 0.3955));
  die2[a][0] = mental[3][2];
  die2[a][3] = mental[3][3];
  die2[a][1] = mental[4][7];
  die2[a][2] = mental[4][6];
  // 3
  a = Same_number(Searchx(0.088855, 0.136145), Searchy(0.4665, 0.4065));
  die2[a][0] = mental[3][3];
  die2[a][3] = mental[3][4];
  die2[a][1] = mental[4][6];
  die2[a][2] = mental[4][5];

  // 第五列
  // 1
  a = Same_number(Searchx(0.17542, 9), Searchy(0.3955, 0.3705));
  die2[a][0] = mental[4][1];
  die2[a][3] = mental[4][2];
  // 2
  a = Same_number(Searchx(0.177805, 9), Searchy(0.4065, 0.3955));
  die2[a][0] = mental[4][2];
  die2[a][3] = mental[4][3];
  // 3
  a = Same_number(Searchx(0.178855, 9), Searchy(0.4665, 0.4065));
  die2[a][0] = mental[4][3];
  die2[a][3] = mental[4][4];

  die = die2;
  /*
    a=Same_number(Searchx(-0.6125,-0.118924),Searchy(0.5545,0.5445));
    die.data()[a].data()[1].data()[0]=mental.data()[2].data()[0].data()[0];
    die.data()[a].data()[1].data()[1]=mental.data()[2].data()[0].data()[1];
    die.data()[a].data()[2].data()[0]=mental.data()[2].data()[3].data()[0];
    die.data()[a].data()[2].data()[1]=mental.data()[2].data()[3].data()[1];
  */
}

// a=Same_number(Search(Point<2,double>()),Search(Point<2,double>()));

std::size_t txtRead::Same_number(std::vector<std::size_t> a,
                                 std::vector<std::size_t> b) {
  std::size_t n = 999;
  for (auto i : a) {
    for (auto j : b) {
      if (i == j) {
        n = i;
        return n;
      }
    }
  }
  return n;
}

void txtRead::Replace_init() {
  std::size_t a;
  a = Same_number(Searchx(-9, 9), Searchy(0.4875, 0.4665));
  cellAdd ce;

  ce += die[a].get_name();
  ce.add_number(12);
  ce += Point<2>{-0.6125, 0.4665};
  ce += mental[2][5];
  ce += mental[2][4];

  ce += master[5];
  ce += master[4];

  ce += mental[3][5];
  ce += mental[3][4];

  ce += mental[4][5];
  ce += mental[4][4];

  ce += Point<2>{0.68, 0.4665};
  ce += Point<2>{0.68, 0.4875};
  ce += Point<2>{-0.6125, 0.4875};
  ce.add_value(die[a].get_value());
  die[a] = ce;
}
