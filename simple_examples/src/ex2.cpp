//Serialisation de simples variables dans des stringstream

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <iostream>
#include <sstream>

using namespace boost::archive;

std::stringstream ss;

void save()
{
  text_oarchive oa(ss);
  int i = 1;
  int j = 42;
  oa << i;
  oa << j;
}

void load()
{
  text_iarchive ia(ss);
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