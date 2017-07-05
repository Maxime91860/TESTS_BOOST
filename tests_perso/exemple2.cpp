#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>

// archives Boost
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

// pour la sérialisation de std::vector
#include <boost/serialization/vector.hpp>

using namespace boost::archive;

/**
 * classe basique à sérialiser
 */
class User{
public:
    User(){}
    std::string nom, prenom;
    int num, age;

    void display(){
        std::cout << nom << " " << prenom << " " << num << " " << age << std::endl;
    }

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version){
        ar & nom & prenom & num & age;
    }
};

// enregistrer une collection dans un fichier
template <class T>
void saveIntoFile(std::vector<User>& d, const char* file){
    std::ofstream ofile(file);
    T ar(ofile);

    ar << d;
    ofile.close();
}

// charger une collection depuis un fichier
template <class T>
void getFromFile(std::vector<User>& d, const char* file){
    std::ifstream ifile(file);
    T ar(ifile);

    ar >> d;
    ifile.close();
}

int main (){
    // créer un tableau d'objets
    std::vector<User> d;
    for (int i = 0; i<10; i++){
        User u;
        u.age = rand()%30;
        u.num = i;
        
        std::ostringstream ss1, ss2;
        ss1 << "nom" << i;
        u.nom = ss1.str();

        ss2 << "prenom"<< i;
        u.prenom =  ss2.str();
        d.push_back(u);

        u.display() ;
    }

    // sauver le tableau d'objets dans un fichier
    saveIntoFile<text_oarchive>(d, "out.txt");

    std::cout << std::endl << std::endl;
    d.clear();

    // relire les données depuis le fichier
    getFromFile<text_iarchive>(d, "out.txt");
    for (int i = 0; i<d.size(); i++)
        d[i].display();
        
    return 0;
}