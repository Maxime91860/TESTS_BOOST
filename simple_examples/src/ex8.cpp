//Serialisation d'une instance de classe pointée

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
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

  void add(int i){ 
    list_.push_back(i); 
  }

  int list(int i){
    return list_[i];
  }

  int list_size () { 
    return list_.size();
  }

  void display () {
    std::cout << "# legs = " << legs_ << ", name = " <<  name_ << '\n';
    for (int i; i < list_.size(); i++){
      std::cout << list_[i] << '\n';  
    }
  }


private:
  friend class boost::serialization::access;

  template <typename Archive>
  void serialize(Archive &ar, const unsigned int version) { 
    ar & legs_;
    ar & name_;
    ar & list_;
  }

  int legs_;
  std::string name_;
  std::vector<int> list_;
};

void save()
{
  boost::archive::text_oarchive oa(ss);
  animal *a = new animal(4,"cat");
  a->add(6);
  a->add(7);
  a->add(10);
  oa << a;
  delete a;
}

void load()
{
  boost::archive::text_iarchive ia(ss);
  animal *a;
  ia >> a;

  a->display();

  delete a;
}

int main()
{
  save();
  load();
}