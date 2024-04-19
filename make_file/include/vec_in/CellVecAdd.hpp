#pragma once
#include <cell_in/cellAdd.hpp>
#include <vec_in/CellVec.hpp>
// #include <vec_in/CellVec.hpp>

class CellVecAdd : public CellVec {
  // like CellVec,but replace cell by cellAdd
public:
  CellVecAdd &operator+=(cellAdd cell) {
    this->m_data.push_back(cell);
    return *this;
  }
  void show() {
    std::cout << m_number << std::endl;
    for (auto i : m_data) {
      i.show();
      std::cout << std::endl;
    }
  }
  std::vector<cellAdd>::iterator begin() { return m_data.begin(); }
  std::vector<cellAdd>::iterator end() { return m_data.end(); }
  CellVecAdd &operator=(CellVecAdd Ce) {
    this->m_number = Ce.m_number;
    for (auto i : Ce) {
      *this += i;
    }
    return *this;
  }
  std::vector<cellAdd> data() { return m_data; }

private:
  std::vector<cellAdd> m_data{};
};
