//Serialisation d'un classe attributs simples variables
//Pas de modifications de la classe, ajout de serialize() avec comme argument une instance du type à serialiser.

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <iostream>
#include <sstream>

using namespace boost::archive;

std::stringstream ss;

struct animal
{
  int legs_;

  animal(){}
  animal(int legs) {
    legs_ = legs;
  }
  int legs() const { return legs_; }
};

template <typename Archive>
void serialize(Archive &ar, animal &a, const unsigned int version)
{
  ar & a.legs_;
}

void save()
{
  text_oarchive oa(ss);
  animal a(4);
  oa << a;
}

void load()
{
  text_iarchive ia(ss);
  animal a;
  ia >> a;
  std::cout  << "#legs = " << a.legs() << '\n';
}

int main()
{
  save();
  load();
}