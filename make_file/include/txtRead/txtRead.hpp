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

  /*
if (edge - 4.5 < 1e-6 && edge_space - 4.5 < 1e-6) {

  // 主导体w1
  master.data()[0].data()[0] -= (w1 - 0.045) / 2;
  master.data()[3].data()[0] -= (w1 - 0.045) / 2;
  master.data()[1].data()[0] += (w1 - 0.045) / 2;
  master.data()[2].data()[0] += (w1 - 0.045) / 2;
} else if (edge - 4.5 < 1e-6) {
  // 主导体
  master.data()[0].data()[0] -= (w1 - 0.045) / 2;
  master.data()[3].data()[0] -= (w1 - 0.045) / 2;
  master.data()[1].data()[0] += (w1 - 0.045) / 2;
  master.data()[2].data()[0] += (w1 - 0.045) / 2;
  // 左边导体
  mental.data()[2].data()[1].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045));
  mental.data()[2].data()[2].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045));
  mental.data()[2].data()[0].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045) + (w1 - 0.045));
  mental.data()[2].data()[3].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045) + (w1 - 0.045));
} else if (edge_space - 4.5 < 1e-6) {
  // 主导体w1
  master.data()[0].data()[0] -= (w1 - 0.045) / 2;
  master.data()[3].data()[0] -= (w1 - 0.045) / 2;
  master.data()[1].data()[0] += (w1 - 0.045) / 2;
  master.data()[2].data()[0] += (w1 - 0.045) / 2;
  // 右边导体
  mental.data()[3].data()[0].data()[0] += (edge - 0.0225 + (w1 - 0.045) / 2);
  mental.data()[3].data()[3].data()[0] += (edge - 0.0225 + (w1 - 0.045) / 2);
  mental.data()[3].data()[1].data()[0] +=
      (edge - 0.0225 + (w1 - 0.045) / 2 + (w2 - 0.045));
  mental.data()[3].data()[2].data()[0] +=
      (edge - 0.0225 + (w1 - 0.045) / 2 + (w2 - 0.045));
} else {
  // 主导体w1
  master.data()[0].data()[0] -= (w1 - 0.045) / 2;
  master.data()[3].data()[0] -= (w1 - 0.045) / 2;
  master.data()[1].data()[0] += (w1 - 0.045) / 2;
  master.data()[2].data()[0] += (w1 - 0.045) / 2;
  // 右边导体
  mental.data()[3].data()[0].data()[0] += (edge - 0.0225 + (w1 - 0.045) / 2);
  mental.data()[3].data()[3].data()[0] += (edge - 0.0225 + (w1 - 0.045) / 2);
  mental.data()[3].data()[1].data()[0] +=
      (edge - 0.0225 + (w1 - 0.045) / 2 + (w2 - 0.045));
  mental.data()[3].data()[2].data()[0] +=
      (edge - 0.0225 + (w1 - 0.045) / 2 + (w2 - 0.045));
  // 左边导体
  mental.data()[2].data()[1].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045));
  mental.data()[2].data()[2].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045));
  mental.data()[2].data()[0].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045) + (w1 - 0.045));
  mental.data()[2].data()[3].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045) + (w1 - 0.045));
  // 最右导体
  mental.data()[4].data()[0].data()[0] += (w1 - 0.045) / 2 + (w2 - 0.045) +
                                          (edge - 0.0225) +
                                          (edge_space - 0.045);
  mental.data()[4].data()[3].data()[0] += (w1 - 0.045) / 2 + (w2 - 0.045) +
                                          (edge - 0.0225) +
                                          (edge_space - 0.045);
  mental.data()[4].data()[1].data()[0] += (w1 - 0.045) / 2 +
                                          2 * (w2 - 0.045) + (edge - 0.0225) +
                                          (edge_space - 0.045);
  mental.data()[4].data()[2].data()[0] += (w1 - 0.045) / 2 +
                                          2 * (w2 - 0.045) + (edge - 0.0225) +
                                          (edge_space - 0.045);
}*/

  // 主导体w1
  master.data()[0].data()[0] -= (w1 - 0.045) / 2;
  master.data()[3].data()[0] -= (w1 - 0.045) / 2;
  master.data()[1].data()[0] += (w1 - 0.045) / 2;
  master.data()[2].data()[0] += (w1 - 0.045) / 2;
  // 右边导体
  mental.data()[3].data()[0].data()[0] += (edge - 0.0225 + (w1 - 0.045) / 2);
  mental.data()[3].data()[3].data()[0] += (edge - 0.0225 + (w1 - 0.045) / 2);
  mental.data()[3].data()[1].data()[0] +=
      (edge - 0.0225 + (w1 - 0.045) / 2 + (w2 - 0.045));
  mental.data()[3].data()[2].data()[0] +=
      (edge - 0.0225 + (w1 - 0.045) / 2 + (w2 - 0.045));
  // 左边导体
  mental.data()[2].data()[1].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045));
  mental.data()[2].data()[2].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045));
  mental.data()[2].data()[0].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045) + (w1 - 0.045));
  mental.data()[2].data()[3].data()[0] -=
      ((w1 - 0.045) / 2 + (edge_space - 0.045) + (w1 - 0.045));
  // 最右导体
  mental.data()[4].data()[0].data()[0] +=
      (w1 - 0.045) / 2 + (w2 - 0.045) + (edge - 0.0225) + (edge_space - 0.045);
  mental.data()[4].data()[3].data()[0] +=
      (w1 - 0.045) / 2 + (w2 - 0.045) + (edge - 0.0225) + (edge_space - 0.045);
  mental.data()[4].data()[1].data()[0] += (w1 - 0.045) / 2 + 2 * (w2 - 0.045) +
                                          (edge - 0.0225) +
                                          (edge_space - 0.045);
  mental.data()[4].data()[2].data()[0] += (w1 - 0.045) / 2 + 2 * (w2 - 0.045) +
                                          (edge - 0.0225) +
                                          (edge_space - 0.045);

  // 介质
  // 待调试代码
  // 1
  std::size_t a;
  a = Same_number(Searchx(-9, -0.114519), Searchy(0.6535, 0.5545));
  die.data()[a].data()[2].data()[0] =
      mental.data()[2].data()[0].data()[0] - 0.005009;
  die.data()[a].data()[3].data()[0] =
      mental.data()[2].data()[3].data()[0] - 0.005009;
  // 2
  a = Same_number(Searchx(-0.120499, -0.10951), Searchy(0.6535, 0.5545));
  die.data()[a].data()[1].data()[0] = mental.data()[2].data()[0].data()[0];
  die.data()[a].data()[2].data()[0] = mental.data()[2].data()[3].data()[0];
  die.data()[a].data()[0].data()[0] =
      mental.data()[2].data()[0].data()[0] - 0.005009;
  die.data()[a].data()[3].data()[0] =
      mental.data()[2].data()[3].data()[0] - 0.005009;
  // 3
  a = Same_number(Searchx(-0.07049, -0.059501), Searchy(0.6535, 0.5545));
  die.data()[a].data()[0].data()[0] = mental.data()[2].data()[1].data()[0];
  die.data()[a].data()[3].data()[0] = mental.data()[2].data()[2].data()[0];
  die.data()[a].data()[1].data()[0] =
      mental.data()[2].data()[1].data()[0] + 0.005009;
  die.data()[a].data()[2].data()[0] =
      mental.data()[2].data()[2].data()[0] + 0.005009;
  // 4
  a = Same_number(Searchx(-0.065481, -0.024519), Searchy(0.6535, 0.5545));
  die.data()[a].data()[0].data()[0] =
      mental.data()[2].data()[1].data()[0] + 0.005009;
  die.data()[a].data()[5].data()[0] =
      mental.data()[2].data()[2].data()[0] + 0.005009;
  die.data()[a].data()[3].data()[0] = master.data()[0].data()[0] - 0.005009;
  die.data()[a].data()[4].data()[0] = master.data()[3].data()[0] - 0.005009;
  // 5
  a = Same_number(Searchx(-0.030499, -0.01951), Searchy(0.6535, 0.5545));
  die.data()[a].data()[1].data()[0] = master.data()[0].data()[0];
  die.data()[a].data()[2].data()[0] = master.data()[3].data()[0];
  die.data()[a].data()[0].data()[0] = master.data()[0].data()[0] - 0.005009;
  die.data()[a].data()[3].data()[0] = master.data()[3].data()[0] - 0.005009;
  // 6
  a = Same_number(Searchx(0.01951, 0.030499), Searchy(0.6535, 0.5545));
  die.data()[a].data()[0].data()[0] = master.data()[1].data()[0];
  die.data()[a].data()[3].data()[0] = master.data()[2].data()[0];
  die.data()[a].data()[1].data()[0] = master.data()[1].data()[0] + 0.005009;
  die.data()[a].data()[2].data()[0] = master.data()[2].data()[0] + 0.005009;
  // 7
  a = Same_number(Searchx(0.024519, 0.042981), Searchy(0.6535, 0.5545));
  die.data()[a].data()[0].data()[0] = master.data()[1].data()[0] + 0.005009;
  die.data()[a].data()[5].data()[0] = master.data()[2].data()[0] + 0.005009;
  die.data()[a].data()[3].data()[0] =
      mental.data()[3].data()[0].data()[0] - 0.005009;
  die.data()[a].data()[4].data()[0] =
      mental.data()[3].data()[3].data()[0] - 0.005009;
  // 8
  a = Same_number(Searchx(0.037001, 0.04799), Searchy(0.6535, 0.5545));
  die.data()[a].data()[0].data()[0] =
      mental.data()[3].data()[0].data()[0] - 0.005009;
  die.data()[a].data()[3].data()[0] =
      mental.data()[3].data()[3].data()[0] - 0.005009;
  die.data()[a].data()[1].data()[0] = mental.data()[3].data()[0].data()[0];
  die.data()[a].data()[2].data()[0] = mental.data()[3].data()[3].data()[0];
  // 9
  a = Same_number(Searchx(0.08701, 0.097999), Searchy(0.6535, 0.5545));
  die.data()[a].data()[0].data()[0] = mental.data()[3].data()[1].data()[0];
  die.data()[a].data()[3].data()[0] = mental.data()[3].data()[2].data()[0];
  die.data()[a].data()[1].data()[0] =
      mental.data()[3].data()[1].data()[0] + 0.005009;
  die.data()[a].data()[2].data()[0] =
      mental.data()[3].data()[2].data()[0] + 0.005009;
  // 10
  a = Same_number(Searchx(0.092019, 0.132981), Searchy(0.6535, 0.5545));
  die.data()[a].data()[0].data()[0] =
      mental.data()[3].data()[1].data()[0] + 0.005009;
  die.data()[a].data()[5].data()[0] =
      mental.data()[3].data()[2].data()[0] + 0.005009;
  die.data()[a].data()[3].data()[0] =
      mental.data()[4].data()[0].data()[0] - 0.005009;
  die.data()[a].data()[4].data()[0] =
      mental.data()[4].data()[3].data()[0] - 0.005009;
  // 11
  a = Same_number(Searchx(0.127001, 0.13799), Searchy(0.6535, 0.5545));
  die.data()[a].data()[0].data()[0] =
      mental.data()[4].data()[0].data()[0] - 0.005009;
  die.data()[a].data()[3].data()[0] =
      mental.data()[4].data()[3].data()[0] - 0.005009;
  die.data()[a].data()[1].data()[0] = mental.data()[4].data()[0].data()[0];
  die.data()[a].data()[2].data()[0] = mental.data()[4].data()[3].data()[0];
  // 12
  a = Same_number(Searchx(0.17701, 0.187999), Searchy(0.6535, 0.5545));
  die.data()[a].data()[0].data()[0] = mental.data()[4].data()[1].data()[0];
  die.data()[a].data()[3].data()[0] = mental.data()[4].data()[2].data()[0];
  die.data()[a].data()[1].data()[0] =
      mental.data()[4].data()[1].data()[0] + 0.005009;
  die.data()[a].data()[2].data()[0] =
      mental.data()[4].data()[2].data()[0] + 0.005009;
  // 13-----------------
  a = Same_number(Searchx(0.182019, 9), Searchy(0.6535, 0.5545));
  die.data()[a].data()[0].data()[0] =
      mental.data()[4].data()[1].data()[0] + 0.005009;
  die.data()[a].data()[4].data()[0] =
      mental.data()[4].data()[2].data()[0] + 0.005009;

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
