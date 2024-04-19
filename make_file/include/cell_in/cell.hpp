#pragma once
#include <cstddef>
#include <point_in/point.hpp>
#include <string>
#include <vector>

class cell {
  // cell class
public:
  cell() = default;
  void reserve(std::size_t m_PointNumber) { m_data.reserve(m_PointNumber); }
  void push_back(Point<2, double> point) { m_data.push_back(point); }
  void pop_back() { m_data.pop_back(); }

  cell &operator+=(Point<2, double> point) {
    this->push_back(point);
    return *this;
  }
  cell &operator+=(std::string name) {
    this->m_name += name;
    return *this;
  }
  cell &operator=(cell ce) {
    this->m_name = ce.m_name;
    this->m_PointNumber = ce.m_PointNumber;
    for (auto i : ce) {
      *this += i;
    }
    return *this;
  }
  std::size_t size() { return m_data.size(); }
  void add_name(std::string st) { m_name += st; }
  void add_number(std::size_t si) { m_PointNumber = si; }
  std::string get_name() { return m_name; }
  std::size_t get_number() { return m_PointNumber; }
  std::vector<Point<2, double>>::iterator begin() { return m_data.begin(); }
  std::vector<Point<2, double>>::iterator end() { return m_data.end(); }
  void show() {
    std::cout << m_name << std::endl;
    std::cout << m_PointNumber << std::endl;
    for (auto i : m_data) {
      i.show();
      std::cout << std::endl;
    }
  }
  std::vector<Point<2, double>> &data() { return m_data; }

protected:
  std::vector<Point<2, double>> m_data{};
  std::string m_name{};
  std::size_t m_PointNumber{};
};
