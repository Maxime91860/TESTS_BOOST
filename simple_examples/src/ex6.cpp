//Serialisation d'une instance de classe pointée

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <iostream>
#include <sstream>

std::stringstream ss;

class animal
{
public:
  animal() = default;
  animal(int legs, std::string name) : legs_(legs), name_(std::move(name)) {}
  int legs() const { return legs_; }
  const std::string &name() const { return name_; }

private:
  friend class boost::serialization::access;

  template <typename Archive>
  void serialize(Archive &ar, const unsigned int version) { 
    ar & legs_;
    ar & name_;
  }

  int legs_;
  std::string name_;
};

void save()
{
  boost::archive::text_oarchive oa{ss};
  animal *a = new animal(4,"cat");
  oa << a;
  std::cout << std::hex << a << '\n';
  delete a;
}

void load()
{
  boost::archive::text_iarchive ia{ss};
  animal *a;
  ia >> a;
  std::cout << std::hex << a << '\n';
  std::cout << std::dec << "# legs = " << a->legs() << ", name = " <<  a->name() << '\n';
  delete a;
}

int main()
{
  save();
  load();
}