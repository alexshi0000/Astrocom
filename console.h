#include <stdlib.h>
#ifndef console
#define console
#include "science.h"
#include "prototype.h"
#ifndef _GLIBCXX_NO_ASSERT
#include <cassert>
#endif
#include <cctype>
#include <cerrno>
#include <cfloat>
#include <ciso646>
#include <climits>
#include <clocale>
#include <cmath>
#include <csetjmp>
#include <csignal>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#if __cplusplus >= 201103L
#include <ccomplex>
#include <cfenv>
#include <cinttypes>
#include <cstdbool>
#include <cstdint>
#include <ctgmath>
#include <cwchar>
#include <cwctype>
#endif

// C++
#include <algorithm>
#include <bitset>
#include <complex>
#include <deque>
#include <exception>
#include <fstream>
#include <functional>
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <iterator>
#include <limits>
#include <list>
#include <locale>
#include <map>
#include <memory>
#include <new>
#include <numeric>
#include <ostream>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <typeinfo>
#include <utility>
#include <valarray>
#include <vector>

#if __cplusplus >= 201103L
#include <array>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <forward_list>
#include <future>
#include <initializer_list>
#include <mutex>
#include <random>
#include <ratio>
#include <regex>
#include <scoped_allocator>
#include <system_error>
#include <thread>
#include <tuple>
#include <typeindex>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#endif
using namespace std;
	//this header file contains the console data
void tokenize_current_process(int i, int j);
string last = "nothing";
vector <string> tokens;
string processed = "";
string display   = "~console:  ";
bool debug_mode = false;
	//leave double space for caret
int caret = 0;	// the caret

void set_time_step(int value){
	time_step = value;
}

void set_delay(int value){
	delay = value;
}

string get_display_text(){
	return display;
}

string substring(string s, int i, int j){
	return s.substr(i,(j-i));
}

void display_help_options(){
	for(int i = 0; i < 28; i++)
		cout << endl;
	cout << "Commands:" << endl;
	cout << "- quit() or exit()" << endl; 				//
	cout << "- simulate(TEST_ID)" << endl;   				//
	cout << "- add(m,x,y,z,i,j,k,)" << endl;   			//
	cout << "- set_width(double)" << endl;      		//
	cout << "- set_ignore_width(double)" << endl;       //
	cout << "- set_color(R,G,B)" << endl;     			//
	cout << "- set_color_all(R,G,B)" << endl; 			//
	cout << "- set_time_step(t)" << endl;               //
	cout << "- set_delay(d)" << endl;                   //
	cout << "- set_radius(r)" << endl;                  //
	cout << "- set_theta(TH)" << endl;                  //
	cout << "- set_time_unit(s,d,y)" << endl;           //
	cout << "- set_visualization(v)" << endl;           //
	cout << "- set_debug_lower_bound(lo)" << endl;
	cout << "- set_debug_upper_bound(hi)" << endl;
	cout << "- clear()"<< endl;
	cout << "- clear_particles()" << endl;
	cout << "- save()" << endl;
	cout << "- load()" << endl;
	cout << "- debug()" << endl;
	cout << "    - debug_tree()" << endl;
	cout << "    - peek(id)" << endl;
	cout << "    - find(id,R,G,B)" << endl;
}

void process(){
	#pragma omp task
	{
		//do something, called after enter is pressed
		cout << ">> " << processed << endl;
		if(processed.compare("quit()") == 0 || processed.compare("exit()") == 0)
			exit(0);
		else if(processed.compare("clear()") == 0){
			for(int i = 0; i < 28; i++)
				cout << endl;
		}
		else if(processed.compare("debug()") == 0)
				//debug mode must be entered in order to turn on debug tree which is a viewtype
			debug_mode = !debug_mode;
		else if(debug_mode && processed.compare("debug_tree()") == 0)
			debug_tree = !debug_tree;
		else if(processed.compare("help()") == 0)
			display_help_options();
		else{
			try{
				if(processed.find("set_time_unit") != string::npos){
					time_unit = substring(processed, processed.find("(")+1, processed.find(")"));
				}
				else if(processed.find("add") != string::npos){
					tokenize_current_process(processed.find("("), processed.find(")"));
					double mass = atof(tokens.at(0).c_str());
					double x = atof(tokens.at(1).c_str()), y = atof(tokens.at(2).c_str()), z = atof(tokens.at(3).c_str());
					double v_i = atof(tokens.at(4).c_str()), v_j = atof(tokens.at(5).c_str()), v_k = atof(tokens.at(6).c_str());
					cout << mass <<" "<< x <<" "<< y <<" "<< z <<" "<< v_i <<" "<< v_j <<" "<< v_k << endl;
					new Particle(mass,x,y,z,v_i,v_j,v_k);
				}
				else if(processed.find("set_width") != string::npos){
					string curr = substring(processed, processed.find("(")+1, processed.find(")"));
					root -> width = atof(curr.c_str());
					root -> height = root -> width;
					root -> length = root -> width;
					width = root -> width;
				}
				else if(processed.find("simulate") != string::npos){
					int id = stoi(substring( processed, processed.find("(")+1, processed.find(")")));
					if(id == 1)
						solar_system_test();
					else if(id == 2)
						center_of_mass_test();
					else if(id == 3)
						galaxy_test();
					else if(id == 4)
						galaxy_test0();
					else if(id == 5)
						cube_test();
					else if(id == 6)
						galaxy_test1();
					else if(id == 7)
						galaxy_test2();
					else
						cout << ">> no simulation with corresponding id found" << endl;
				}
				else if(processed.find("set_color_all") != string::npos){
					tokenize_current_process(processed.find("("), processed.find(")"));
					if(field.size() > 0){
						float R = atof(tokens.at(0).c_str());
						float G = atof(tokens.at(1).c_str());
						float B = atof(tokens.at(2).c_str());
						for(int i = 0; i < field.size(); i++)
							field.at(i) -> set_color(R,G,B);
					}
				}
				else if(processed.find("set_color") != string::npos){
						//same as the above but only changes color for the lastest particle
					tokenize_current_process(processed.find("("), processed.find(")"));
					if(field.size() > 0){
						Particle *p = field.at(field.size()-1);
						float R = atof(tokens.at(0).c_str());
						float G = atof(tokens.at(1).c_str());
						float B = atof(tokens.at(2).c_str());
						p -> set_color(R,G,B);
					}
				}
				else if(processed.find("set_ignore_width") != string::npos){
					ignore_width = atof(substring(processed, processed.find("(")+1, processed.find(")")).c_str());
					for(int i = 0; i < field.size(); i++)
						field.at(i) -> ignore = false;
							//reset the ingore flag to false
				}
				else if(processed.find("set_time_step") != string::npos){
					time_step = std::atof(substring(processed, processed.find("(")+1, processed.find(")")).c_str());
				}
				else if(processed.find("set_radius") != string::npos){
					particle_radius = std::atof(substring(processed, processed.find("(")+1, processed.find(")")).c_str());
				}
				else if(processed.find("set_delay") != string::npos){
					delay = std::stoi(substring(processed, processed.find("(")+1, processed.find(")")).c_str());
						//integer value
				}
				else if(processed.find("set_theta") != string::npos){
					theta = std::atof(substring(processed, processed.find("(")+1, processed.find(")")).c_str());
				}
				else if(processed.find("set_debug_upper_bound") != string::npos){
					debug_upper_bound = std::stoi(substring(processed, processed.find("(")+1, processed.find(")")).c_str());
				}
				else if(processed.find("set_debug_lower_bound") != string::npos){
					debug_lower_bound = std::stoi(substring(processed, processed.find("(")+1, processed.find(")")).c_str());
				}
				else if(processed.find("clear_particles") != string::npos){
					field.clear();
					N = 0;
				}
			} catch (int e){
				cout << ">> invalid command" << endl;
			}
		}
	}
}

void tokenize_current_process(int i, int j){
	//tokenizer (by comma)
	tokens.clear();
		//clear
	string s = substring(processed,i+1,j);
	stringstream ss(s);
	string string_builder = "";
	char curr;
	while(ss >> curr){
		string_builder += curr;
		if(ss.peek() == ','){
			//deliminator
			tokens.push_back(string_builder);
			string_builder = "";
			ss.ignore();
		}
	}
	if(string_builder.compare("") != 0)			//still more
		tokens.push_back(string_builder);
}

void glut_special_func(int key, int x, int y){
	//caret movement for ease of use
	if(key == GLUT_KEY_RIGHT || key == GLUT_KEY_LEFT){
		if(key == GLUT_KEY_RIGHT && caret < processed.length()) {
			caret++;
		}
		else if (key == GLUT_KEY_LEFT && caret > 0){
			caret--;
		}
		if(caret != processed.length() && caret >= 0)
			display = "~console: " + substring(processed,0,caret) + " " + substring(processed,caret,processed.length());
		else
			display = "~console: " + processed + " ";
		display[caret+10] = '|';
	}
	else if(key == GLUT_KEY_UP && last.compare("nothing") != 0){
		processed = last;
		last = "nothing";
		caret = processed.length();
		display = "~console: " + processed + " ";
		display[caret+10] = '|';
	}
}

void do_something(unsigned char key, int x, int y){
	//keyboard functions
	if(key == 127){
		//backspace
		if(caret == processed.length()){
			processed = substring(processed,0,processed.length()-1);
			display = "~console: " + processed + " ";
			caret = processed.length();
			display[caret+10] = '|';
		}
		else if(caret > 0 && caret < processed.length()){
			display = "~console: " + substring(processed, 0, caret-1) + " " + substring(processed, caret, processed.length());
			processed = substring(processed, 0, caret-1) + substring(processed, caret, processed.length());
			caret--;
			display[caret+10] = '|';
		}
	}
	else if(key == 13){
			//enter is pressed
		process();
			//handle the process
		last = processed;
		processed = "";
		display = "~console:  ";
		caret = 0;
		display[caret+10] = '|';
			//reset for next query
	}
	else if(key == 45) {
		zoom /= 1.1;
	}
	else if(key == 43) {
		zoom *= 1.1;
	}
	else{
		if(caret == processed.length()){
			processed += key;
			display = "~console: " + processed + " ";
			caret = processed.length();
		}
		else{
			display = "~console: " + substring(processed,0,caret) + (char)key + " "+ substring(processed,caret,processed.length());
			processed = substring(processed,0,caret) + (char)key + substring(processed,caret,processed.length());
			caret++;
		}
		display[caret+10] = '|';
	}
}

#endif
