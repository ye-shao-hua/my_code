#pragma once
#include <cell_in/cell.hpp>
#include <vector>

class CellVec {
  // CellVec class ,a component by Cell
public:
  CellVec() = default;
  void set_number(std::size_t number) { m_number = number; }
  std::size_t get_number() { return m_number; }
  CellVec &operator+=(cell cell) {
    this->m_data.push_back(cell);
    return *this;
  }
  CellVec &operator=(CellVec Ce) {
    this->m_number = Ce.m_number;
    for (auto i : Ce) {
      *this += i;
    }
    return *this;
  }
  std::vector<cell>::iterator begin() { return m_data.begin(); }
  std::vector<cell>::iterator end() { return m_data.end(); }
  void show() {
    std::cout << m_number << std::endl;
    for (auto i : m_data) {
      i.show();
      std::cout << std::endl;
    }
  }

protected:
  std::vector<cell> m_data{};
  std::size_t m_number{};
};
