

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
// #include <fstream>
#include <sstream>

// class Test;
// const std::ostream& operator<<(const std::ostream &strm, const Test &test);
class Test
{
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize (Archive & ar, const unsigned int version)
	{
		ar & a;
		ar & b;
	}

	int a;
	int b;

	friend  std::ostream& operator<<( std::ostream &strm, const Test &test);

public:
	Test(int a, int b)
	{
		this->a = a;
		this->b = b;
	}

	Test(){}	
};

 std::ostream& operator<<( std::ostream &strm, const Test &test) {
  return strm << "Test(" << test.a << "," << test.b << ")";
}

int main(int argc, char const *argv[])
{
	
	Test var1 = Test(17,9);
	Test var2 = Test(18,10);

	std::cout << var1 << std::endl;
	std::cout << var2 << std::endl;

	//Ouverture du stream de sortie et de l'archive de sortie
	// std::ofstream ofs("files_serialized/test.ser");
	std::stringstream ss;
	boost::archive::text_oarchive oa(ss);


	//Ecriture de l'objet dans l'archive
	oa << var1;
	// oa << var2;

	// {
	// 	//Ouverture du stream d'entrée et de l'archive d'entrée
		// std::ifstream ifs("files_serialized/test.ser");
		boost::archive::text_iarchive ia(ss);

		Test var3;
	// 	Test var4;

		ia >> var3;
	// 	// ia >> var4;

		std::cout << var3 << std::endl;

	// }
	return 0;
}