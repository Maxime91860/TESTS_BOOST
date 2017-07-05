#include <iostream>
// #include <fstream>
#include <sstream>

// include input and output archivers
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

// include this header to serialize vectors
#include <boost/serialization/vector.hpp>

using namespace std;

class gps_position
{
    friend std::ostream & operator<<(std::ostream &os, const gps_position &gp);
    friend class boost::serialization::access;
    int degrees;
    int minutes;
    float seconds;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */){
        ar & degrees & minutes & seconds;
    }
public:
    // every serializable class needs a constructor
    gps_position(){};
    gps_position(int _d, int _m, float _s) : 
        degrees(_d), minutes(_m), seconds(_s)
    {}
};
std::ostream & operator<<(std::ostream &os, const gps_position &gp)
{
    return os << ' ' << gp.degrees << "°" << gp.minutes << '\'' << gp.seconds << '\"';
}

int main()
{

    std::vector<double> v;
    v.push_back (1);
    v.push_back (2);
    v.push_back (3.4);
    v.push_back (5.6);
    v.push_back (5.014625687652);
    v.push_back (18);
    v.push_back (213054);



    // serialize vector

    // std::ofstream ofs(test);
    std::stringstream ofs;
    boost::archive::text_oarchive oa(ofs);
    // std::cout << "Avant sérialisation dans l'archive" << "\n";
    oa << v;


    std::cout << ofs << ofs.str() << std::endl;

    //Conversion stringstream to char*
    std::string vect_serialized = ofs.str();
    std::cout << "vect_serialized = " << vect_serialized << std::endl;

    int taille_buffer = vect_serialized.length() + 1;
    char* vect_serialized_char = new char[taille_buffer];
    strcpy(vect_serialized_char, vect_serialized.c_str()); //vect_serialized.c_str() renvoie une const char*

    std::cout << "taille_buffer = " << taille_buffer << std::endl;
    std::cout << "vect_serialized_char = " << vect_serialized_char << std::endl;


    std::vector<double> v2;

    //Conversion char* to stringstrem

    // load serialized vector into vector 2

    // std::ifstream ifs(test);
    boost::archive::text_iarchive ia(ofs);
    ia & v2;


    std::cout << "Affichage v2 :\n";
    for (int i = 0; i < v2.size(); ++i)
    {
        std::cout << v2[i] << "\n";
    }


    return 0;
}