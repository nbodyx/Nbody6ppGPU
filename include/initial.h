//==================================================//
// Initial paraments interface for common programs  //
// By l.wang========================================//
//==================================================//
//Notes: Used for initialize parameters in program==//
//       beginning. Can use command line input,=====//
//       user interactive input, input from last====//
//       config file.===============================//
//Public member functions:==========================//
// (this)add(str name,str intro,int/double/str value) : add new item.
//    (T)get<T>(str name): get item value to type T.
//  (str)getintro(str name): get item introduction
// (bool)initcin(str name): initial an item with cin.
// (bool)initall_cin(): initial all items with cin.
// (bool)initall_argv(int argc, char* argv[],int begin(1)):
//       initial all items with argv. read data from argv[begin-1].
// (bool)initial_iconf(): initial all items from config.
// (bool)initial_oconf(): output all items name and value to config.
// (bool)initial(int argc, char* argv[], int begin(1), bool print(1), bool kill(0)):
//       initial function, print=1 will print information, kill=1 will terminate if error exist.
// (int) set(str name): set currentname=name for <<, return 1(int),2(double),3(string).
// (int) getnextargc(): get next input order, can be used to set as beginning for initial for second pars_initial class.
//Operator: =: copy all items, config file name.
//Operator: <<: out value of currentname item.
//==================================================//
//--06/26/12 --:---creat--------------------------------//
//--07/01/12 18:58-Add synthesized assignment operator--//
//--07/26/12 15:02-Rebuild the code---------------------//

#ifndef initial_h
#define initial_h

#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <utility>
#include <vector>

typedef std::map<std::string,std::pair<std::string,std::string> > mapitem;
typedef std::map<std::string,std::pair<std::string,std::string> >::iterator mapitemitr;
typedef std::pair<std::string,std::string>* pairitr;

class pars_initial {
 public:
  // Add items
  template <typename T>
  pars_initial& add(const std::string& name, const std::string& intro, const T& data);
  
  // Get item value===================================//
  template <typename T>
  T get(const std::string& name);

  double getd(const std::string& name); 
  int geti(const std::string& name); 
  float getf(const std::string& name);
  std::string gets(const std::string& name);
  bool getb(const std::string& name);

  //item check, if exist, return true=================//
  bool itemcheck(const std::string& name);
  
  // Get item introduction============================//
  std::string getintro(const std::string& name);
  
  //print all items===================================//
  void print()  { pars_print("current"); }
    
  // initial function ================================//
  void initial(const int& argc, char* argv[], bool info_print=true, bool kill_if_unsuccess=false);

  // Use for multi initial conditions=================//
  //return next unread argv position==================//
  int getnextargc();
    
  // Constructor======================================//
  pars_initial(): argvleft(1),config(".config") {}
  pars_initial(const std::string& name): argvleft(1),config(name) {}
  pars_initial(const std::string& name,int offset): argvleft(offset),config(name) {}
  // Synthesized Assignment operator==================//
  pars_initial& operator=(const pars_initial& pio);

  // Destructor=======================================//
  ~pars_initial() {}

private:

  // --- from cin ---
  bool init_cin();

  // --- from argv ---
  // return true: success; false: fail
  bool init_argv(const int &argc, char* argv[]);

  // --- from config ---
  bool init_conf();

  // Save config data=================================//
  bool out_conf();

  // Initial from iostream============================//
  bool pars_print(const std::string& extraintro);

  // members==========================================//
  int argvleft;         //begining index of unread argv[]
  std::vector<std::string> namelist; //arguments list
  std::string config;   //Name of config file
  mapitem items;        //string items
};

template <typename T>
pars_initial& pars_initial::add(const std::string& name, const std::string& intro, const T& data)
{
  std::stringstream sstr;
  items[name].first=intro;
  sstr.str("");
  sstr.clear();
  sstr<<data;
  sstr>>items[name].second;
  namelist.push_back(name);
  return *this;
}

template <typename T>
T pars_initial::get(const std::string& name) {
  if(itemcheck(name)) {
    std::stringstream sstr;
    T tmp;
    sstr.str("");
    sstr.clear();
    sstr<<items[name].second;
    sstr>>tmp;
    return tmp;
  }
  return 0;
}

#endif
