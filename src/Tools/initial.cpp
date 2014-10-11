//==================================================//
// Source for class pars_initial====================//
//==================================================//

#include "initial.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <map>
#include <utility>
#include <fstream>
#include <vector>

//check whether name is inside map items============//
bool pars_initial::itemcheck(const std::string& name){
  if(!items.count(name)) {
    std::cerr<<"Error: "<<name<<" not found in arguments list!"<<std::endl;
    return false;
  }
  return true;
}

double pars_initial::getd(const std::string& name) {
  return get<double>(name);
}

int pars_initial::geti(const std::string& name) {
  return get<int>(name);
}

float pars_initial::getf(const std::string& name) {
  return get<float>(name);
}

std::string pars_initial::gets(const std::string& name) {
  return get<std::string>(name);
}  

bool pars_initial::getb(const std::string& name) {
  return get<bool>(name);
}  

std::string pars_initial::getintro(const std::string& name) {
  if(itemcheck(name)) return items[name].first;
  return "";
}

bool pars_initial::pars_print(const std::string& extraintro) {
  for(int i=1; i<50;i++) std::cout<<"-";
  std::cout<<std::endl;
  std::cout<<"Parameters\tDefination\t"<<extraintro<<std::endl;
  for(int i=0; i<namelist.size();++i)  {
    if(!itemcheck(namelist[i])) return false;
    pairitr it=&items[namelist[i]];
    std::cout<<namelist[i]<<"\t"<<it->first<<"\t"<<it->second<<std::endl;
  }
  for(int i=1; i<50;i++) std::cout<<"-";
  std::cout<<std::endl;
  return true;
}

//private===========================================//
bool pars_initial::init_cin() {
  for(int i=0; i<namelist.size();++i)  {
    if(!itemcheck(namelist[i])) return false;
    pairitr it=&items[namelist[i]];
    std::string check;
    std::cout<<namelist[i]<<" ("<<it->first<<")\nCurrent value: "<<it->second<<".   Change it? (y or n)"<<std::endl;
    while (std::cin>>check)
    {
      if (std::cin && check[0]=='n' && check.length()==1) break;
      else if (std::cin && check[0]=='y' && check.length()==1) {
        std::cout<<"Input new "<<namelist[i]<<": ";
        std::cin>>it->second;
        break;
      }
      else std::cerr<<"Input uncorrect, please choose y or n: ";
    }
  }
  return true;
}

bool pars_initial::init_argv(const int& argc, char* argv[]) {
  if (argc-argvleft >= namelist.size())  {
    for(int i=0; i<namelist.size();++i)  {
      if(!itemcheck(namelist[i])) return false;
      items[namelist[i]].second=argv[argvleft];
      argvleft++;
    }
    return true;
  }
  else std::cerr<<"Error: Input parameters incomplele, (input)"<<argc-argvleft<<" < (required)"<<namelist.size()<<std::endl;
  return false;
}

bool pars_initial::init_conf() {
  std::ifstream configi(config.c_str());
  if (!configi.is_open())  {
    std::cerr<<"No config file \""<<config<<"\""<<std::endl;
    return false;
  }
  else {
    std::string name;
    while(true) {
      configi>>name;
      if(configi.eof()) break;
      if(items.count(name)) configi>>items[name].second;
      else {
        std::cerr<<"Warning!: miss match parameter name: "<<name<<std::endl;
        configi>>name;
      }
    }
  }
  configi.close();
  return true;
}

bool pars_initial::out_conf() {
  std::ofstream configo(config.c_str());
  for(int i=0; i<namelist.size();++i)  {
    if(!itemcheck(namelist[i])) return false;
    configo<<namelist[i]<<"\t"<<items[namelist[i]].second<<std::endl;
  }
  return true;
}

void pars_initial::initial(const int& argc, char * argv[], bool info_print, bool kill_if_unsuccess) {
  std::string inflag;
  std::string checkopt;
  if (argc-argvleft==0) {
    std::cout<<"Reset parameters (y/n?)";
    while (std::cin>>inflag) {
      if (inflag[0]=='y' && inflag.length()==1) {
        if (!init_conf()&&kill_if_unsuccess) exit(1);
        if (!init_cin()&&kill_if_unsuccess) exit(1);
        out_conf();
        break;
      }
      else if (inflag[0]=='n' && inflag.length()==1) {
        if(!init_conf()&&kill_if_unsuccess) exit(1);
        break;
      }
      else std::cerr<<"Error: Input uncorrect, please choose y or n: ";
    }
  }
  else {
    checkopt=argv[argvleft];
    if (checkopt=="-l" && checkopt.length()==2) {
      argvleft++;
      if(!init_conf()&&kill_if_unsuccess) exit(1);
    }
    else if (checkopt=="-h" && checkopt.length()==2) {
      pars_print("Default");
      exit(1);
    }
    else if (checkopt=="-f" && checkopt.length()==2) {
      config=argv[argvleft+1];
      argvleft +=2;
      if(!init_conf()&&kill_if_unsuccess) exit(1);
    }
    else {
      if(!init_argv(argc,argv)&&kill_if_unsuccess) exit(1);
      out_conf();
    }
  }
  if (info_print) pars_print("Current");
}

int pars_initial::getnextargc() {
  return argvleft;
}

pars_initial& pars_initial::operator=(const pars_initial& pio) {
  argvleft=pio.argvleft;
  namelist=pio.namelist;
  config=pio.config;
  items=pio.items;
  return *this;
}
