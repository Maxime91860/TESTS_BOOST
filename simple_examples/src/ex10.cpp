//Serialization of pointers
//05/07/2017 - Maxime Kermarquer

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <iostream>

// #include <fstream>
#include <sstream>
std::stringstream ss;


class A
{
	public :
		A() {}

		A(int size){
			m_size = size;
			tab_1D = new int[size];
			for (int i = 0; i < size; i++){
				tab_1D[i] = std::rand()%20;
			}

			tab_2D = new int*[size];
			for (int i = 0; i < size; i++){
				tab_2D[i] = new int[size];
				for (int j = 0; j < size; j++){
					tab_2D[i][j] = std::rand()%20;
				}
			}
		}

		void display(){
			std::cout << "tab_1D = \n";
			display_tab1D(tab_1D, m_size);

			std::cout << "tab_2D = \n";
			for (int i = 0; i < m_size; i++)
			{
				std::cout << std::hex << tab_2D[i] << ": ";
				std::cout << std::dec;
				display_tab1D(tab_2D[i], m_size);
			}

		}

	private :
		int m_size;
		int *tab_1D;
		int **tab_2D;

		void display_tab1D (int	*tab, int size){
			for (int i = 0; i < size; i++){
				if (i != (size - 1)){
					std::cout << tab[i] << ", ";
				} 
				else {
					std::cout << tab[i] << "\n";	
				}
			}
		}

		friend class boost::serialization::access;
		template <typename Archive>
		void serialize(Archive &ar, const unsigned int version){
			ar & m_size;
			
			//Allocation for the deserialization
			if(Archive::is_loading::value){
				tab_1D = new int[m_size];
			}

			// for (int i = 0; i < m_size; i++){
			// 	ar & tab_1D[i];
			// }
			
			// This statement is equivalent to the commented loop above
			ar & boost::serialization::make_array <int> (tab_1D, m_size);

			//Allocation for the deserialization
			if(Archive::is_loading::value){
				tab_2D = new int*[m_size];
				for (int i = 0; i < m_size; i++){
					tab_2D[i] = new int[m_size];
				}				
			}

			for (int i = 0; i < m_size; i++){
				ar & boost::serialization::make_array <int> (tab_2D[i], m_size);				
			}

		}

};

void save ()
{
	
	// std::ofstream file("pointers.ser");
	// boost::archive::text_oarchive oa(file);

	boost::archive::text_oarchive oa(ss);

	A *a1 = new A(5);
	std::cout << "---- a1 ----\n";
	a1->display();
	std::cout << "-------------\n\n";

	oa << a1;
	delete a1;
}

void load ()
{
	// std::ifstream file("pointers.ser");
	// boost::archive::text_iarchive ia(file);

	boost::archive::text_iarchive ia(ss);

	A *a2;
	ia >> a2;
	std::cout << "---- a2 ----\n";
	a2->display();
	std::cout << "-------------\n";
	delete a2;
}

int main()
{
	save();
	load();
}