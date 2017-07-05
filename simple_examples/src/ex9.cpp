//Serialization of vectors

#include <vector>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <iostream>
// #include <sstream>
#include <fstream>

// std::stringstream ss;

class A
{
	public :
		A() {}

		A(int size_tab) { 
			m_size_tab = size_tab;
		}

		void display (){
			std::cout << "size_tab = " << m_size_tab << '\n';
			for (int i = 0; i < m_v.size(); i++){
				if (i != (m_v.size() - 1)){
					std::cout << m_v[i] << ", ";
				} 
				else {
					std::cout << m_v[i] << "\n";	
				}
			}
		}

		void push_double (double x){
			m_v.push_back(x);
		}

	private :

		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive &ar, const unsigned int version) {
			ar & m_size_tab;
			ar & m_v;
		}

		int m_size_tab;
		std::vector<double> m_v;
};


void save()
{
	// boost::archive::text_oarchive oa(ss);
	std::ofstream file("vector.ser");
	boost::archive::text_oarchive oa(file);

	int size = 5;
	A *a = new A(size);

	a->push_double(5.00);
	a->push_double(10.00);
	a->push_double(3.14);

	oa << a;
}

void load()
{
	// boost::archive::text_iarchive ia(ss);
	std::ifstream file("vector.ser");
	boost::archive::text_iarchive ia(file);

	A *a;
	ia >> a;
	a -> display();
}

int main()
{
	save();
	load();
}