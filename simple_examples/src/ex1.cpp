//Serialisation de simples variables dans des fichiers

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <iostream>
#include <fstream>

using namespace boost::archive;

void save()
{
  std::ofstream file("archive.ser");
  text_oarchive oa(file);
  int i = 1;
  int j = 42;
  oa << i;
  oa << j;
}

void load()
{
  std::ifstream file("archive.ser");
  text_iarchive ia(file);
  int i = 0;
  int j = 0;
  ia >> i;
  ia >> j;
  std::cout << "i = " << i << ", j = " << j << '\n';
}

int main()
{
  save();
  load();
}