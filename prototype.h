#include <stdlib.h>
#ifndef prototype
#define prototype
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
#include "science.h"

void user_def(){
		//do nothing, user definations
	cout << "Welcome to astrocom" << endl;
	cout << "Introduction:" << endl;
	cout << "- n-body simulation software" << endl;
	cout << "- scientific data visualization" << endl;
	cout << "- efficient octree algorithm" << endl;
	cout << "- cpu multithread parallelism support" << endl;
	cout << "- type help() to begin" << endl;
		//user defined simulation only
	theta = 0.5;
	time_step = 250000;
	width = 20*AU;
	delay = 10;
}

void solar_system_test(){
	//test for accuracy by using this function and comparing earth x == 0 at T = 365.25 days
	Particle *sun = new Particle(1.989e30,0,0,0,0,0,0);	//sun
	Particle *earth = new Particle(5.972e24,AU,0,0,0,2.978589e4,0); //earth
	Particle *venus = new Particle(5.972e24,108.2e9,0,0,0,35e3,0); //venus
	Particle *mars = new Particle(6.39e23,227.9e9,0,0,0,24.1e3,0); //mars
	sun -> set_color(1,0.98,0.75);
	earth -> set_color(0.1,1,0.1);
	venus -> set_color(0.8,0.8,0.15);
	mars -> set_color(0.9,0.3,0.1);
	particle_radius = 0.005;
	theta = 0.5;
	time_step = 250000;
	width = 20*AU;
	delay = 10;
}

void center_of_mass_test(){
	//tests how accurate the COM APROX is with four very massive bodies
	Particle *a = new Particle(1e31,AU,AU,0,0,-6e4,0);
	Particle *b = new Particle(1e31,-AU,AU,0,6e4,0,0);
	Particle *c = new Particle(1e31,AU,-AU,0,-6e4,0,0);
	Particle *d = new Particle(1e31,-AU,-AU,0,0,6e4,0);
	a -> set_color(1,0,0);
	b -> set_color(0,1,0);
	c -> set_color(0,0,1);
	d -> set_color(1,1,0);
	theta = 100;
	time_step = 25000;
	width = 20*AU;
	delay = 10;
}

void disk_test(bool bulge, unsigned int PARTICLE_COUNT, double radius, double init_vel, float red, float green, float blue){
	if(bulge){
		Particle *q = new Particle(1e33,0,0,0,0,0,0);
		q->set_color(0,0,0);
		q->black_hole = true;
	}
	for(unsigned int i = 0; i < PARTICLE_COUNT; i++){
		double x = ((double)rand()/(double)RAND_MAX) * 2.0 * radius - radius;
		double y = ((double)rand()/(double)RAND_MAX) * 2.0 * sqrt(radius*radius - x*x) - sqrt(radius*radius - x*x);
		double z = 0;
		if(bulge)
			z = ((double)rand()/(double)RAND_MAX) * 0.8 * sqrt(radius*radius - x*x) - 0.4 * sqrt(radius*radius - x*x);
		else
			z = ((double)rand()/(double)RAND_MAX) * 0.2 * sqrt(radius*radius - x*x) - 0.1 * sqrt(radius*radius - x*x);
		double m = y / x;
		m = -1/m;
		double b = y - x*m;
		if(y > 0){
			Vect *velocity = new Vect(-1, (x-1)*m + b - y, 0);
			velocity -> scalar(init_vel);
			Particle *q = new Particle(2e23, x, y, z, velocity -> i, velocity -> j, velocity -> k);
			q -> set_color(red,green,blue);
		}
		else{
			Vect *velocity = new Vect(1, (x+1)*m + b - y, 0);
			velocity -> scalar(init_vel);
			Particle *q = new Particle(2e23, x, y, z, velocity -> i, velocity -> j, velocity -> k);
			q -> set_color(red,green,blue);
		}
	}
}

void spiral_arms(unsigned int PARTICLE_COUNT, double radius, double init_vel, float red, float green, float blue){
	for(unsigned int i = 0; i < PARTICLE_COUNT; i++){
		double x = ((double)rand()/(double)RAND_MAX) * 2.0 * radius - radius;
		double y = ((double)rand()/(double)RAND_MAX) * 2.0 * sqrt(radius*radius - x*x) - sqrt(radius*radius - x*x);
			if(abs(y) < 70*AU || abs(x) < 70*AU || (y > 0 && x > 0) || (y < 0 && x < 0)){
				i--;
				continue;
			}
		double z = ((double)rand()/(double)RAND_MAX) * 0.02 * radius - radius * 0.01;
		double m = y / x;
		m = -1/m;
		double b = y - x*m;
		if(y > 0){
			Vect *velocity = new Vect(-1, (x-1)*m + b - y, 0);
			velocity -> scalar(init_vel);
			Particle *q = new Particle(2e24, x, y, z, velocity -> i, velocity -> j, velocity -> k);
			q -> set_color(red,green,blue);
		}
		else{
			Vect *velocity = new Vect(1, (x+1)*m + b - y, 0);
			velocity -> scalar(init_vel);
			Particle *q = new Particle(2e24, x, y, z, velocity -> i, velocity -> j, velocity -> k);
			q -> set_color(red,green,blue);
		}
	}
}

void kuzmin_disk(bool bulge, double mass, int PARTICLE_COUNT, double radius, float red, float green, float blue, double x1, double y1, double z1){
	if(bulge){
		Particle *q = new Particle(1e33,0+x1,0+y1,0+z1,0,0,0);
		q->set_color(0,0,0);
		q->black_hole = true;
	}
	for(unsigned int i = 0; i < PARTICLE_COUNT; i++){
		double x = ((double)rand()/(double)RAND_MAX) * 2.0 * radius - radius;
		double y = ((double)rand()/(double)RAND_MAX) * 2.0 * sqrt(radius*radius - x*x) - sqrt(radius*radius - x*x);
		double z = 0;
		if(bulge)
			z = ((double)rand()/(double)RAND_MAX) * 0.1 * sqrt(radius*radius - x*x) - 0.1 * sqrt(radius*radius - x*x);
		else
			z = ((double)rand()/(double)RAND_MAX) * 0.05 * sqrt(radius*radius - x*x) - 0.05 * sqrt(radius*radius - x*x);
		double m = y / x;
		m = -1/m;
		double b = y - x*m;
		if(y > 0){
			Vect *velocity = new Vect(-1, (x-1)*m + b - y, 0);
			velocity -> scalar(
				sqrt((GC * (1e33)) /
					sqrt(x*x + y*y + z*z))
			);
			Particle *q = new Particle(2e19*mass, x+x1, y+y1, z+z1, velocity -> i, velocity -> j, velocity -> k);
			q -> set_color(red,green,blue);
		}
		else{
			Vect *velocity = new Vect(1, (x+1)*m + b - y, 0);
			velocity -> scalar(
				sqrt((GC * (1e33)) /
					sqrt(x*x + y*y + z*z))
			);
			Particle *q = new Particle(2e19*mass, x+x1, y+y1, z+z1, velocity -> i, velocity -> j, velocity -> k);
			q -> set_color(red,green,blue);
		}
	}
}

double gen_rand(double lo, double hi){
	return (hi-lo) * ((double)rand()/(double)RAND_MAX) + lo;
}

void spheriod(Particle *center, int particle_count, double mass, double radius, double c, double x_shift, double y_shift, double z_shift,
					   double h, double f, double l, double angle){	//galactic center, particle cound, radius, flattening constant, position vector, rotational vector, angle of rotation
	for(int i = 0; i < particle_count; i++){
		double x_pos = gen_rand(-1,1);
		double y_min = -sqrt(1-pow(x_pos,2));
		double y_max = sqrt(1-pow(x_pos,2));
		double y_pos = gen_rand(y_min,y_max);
		double z_max = sqrt(c*(1-pow(x_pos,2)-pow(y_pos,2)));
		double z_min = -sqrt(c*(1-pow(x_pos,2)-pow(y_pos,2)));
		double z_pos = gen_rand(z_min,z_max);
		x_pos *= radius;
		y_pos *= radius;
		z_pos *= radius;
		double m = y_pos / x_pos;
		m = -1/m;
		double b = y_pos - x_pos*m;
		Vect *velocity = NULL;
		if(y_pos > 0){
			velocity = new Vect(-1e2*AU,(x_pos-1e2*AU)*m + b - y_pos, 0);
			velocity->scalar(sqrt(
				(double)(GC*center->m)	/
				(double)(sqrt(pow(x_pos,2) + pow(y_pos,2) + pow(z_pos,2)))
			));
		}
		else{
			velocity = new Vect(1e2*AU,(x_pos+1e2*AU)*m + b - y_pos, 0);
			velocity->scalar(sqrt(
				(double)(GC*center->m)	/
				(double)(sqrt(pow(x_pos,2) + pow(y_pos,2) + pow(z_pos,2)))
			));
		}
		new Particle(mass, x_pos+x_shift, y_pos+y_shift, z_pos+z_shift, velocity->i, velocity->j, velocity->k);
	}
}

void galaxy_test2(){

}

void galaxy_test1(){
	Particle *center = new Particle(1e37,0,0,0,0,0,0);
	center->black_hole = true;
	double radius = 1e8*AU;
	double mass = 1e30;
	for(double i = 0.0002; i <= 1; i=i*1.15){
		spheriod(center,200,mass,radius, i, 0,0,0 ,0,0,0, 0);
		radius = radius * 0.97;
		mass *= 0.92;
	}
	particle_radius = 0.0035;
	alpha = 0.2333;
	render_quality = 4;
	theta = 1000;
	time_step = 5e9;
	width = 2e12*AU;
	delay = 10;
	ignore_width = 1e12*AU;
}

void galaxy_test0(){
	//kuzmin_disk(true,1e3,500,1e4*AU,1,1,1,0,0,0);
	kuzmin_disk(true,5e8,100000,2e4*AU,1,1,1,0,0,0);
	particle_radius = 0.002;
	alpha = 0.333;
	render_quality = 4;
	theta = 0.5;
	time_step = 2e9;
	width = 6e7*AU;
	delay = 20;
	ignore_width = 5e7*AU;
	debug_lower_bound = 0;
	debug_upper_bound = 10;
	time_unit = "y";
}

void galaxy_test(){
	/*disk_test(true, 3000, 50*AU, 1e5, 1, 0, 0);
	disk_test(false, 3500, 100*AU, 0.75e5, 0, 1, 0);
	spiral_arms(4000, 200*AU, 0.51e5, 0, 0, 1);*/

	disk_test(true, 1250, 50*AU, 1e5, 1, 0.1, 0.1);
	disk_test(false, 2750, 100*AU, 0.75e5, 0.1, 1, 0.1);
	spiral_arms(1250, 200*AU, 0.51e5, 0.1, 0.1, 1);
	particle_radius = 0.004;
	render_quality = 3;

	theta = 1000000;
	time_step = 5500000;
	width = 1e4*AU;
	delay = 10;
	ignore_width = 1e3*AU;
	debug_lower_bound = 0;
	debug_upper_bound = 20;
}

void cube(double offset_x, double offset_y){
	for(double i = -AU*36; i < AU*36; i+=AU*6){
		for(double j = -AU*36; j < AU*36; j+=AU*6){
			for(double k = -AU*36; k < AU*36; k+=AU*6){
				new Particle(1e31,i+offset_x,j+offset_y,k,0,0,0);
				field.at(field.size()-1) -> set_color((i+AU*28)/(AU*56), (j+AU*28)/(AU*56), (k+AU*28)/(AU*56));
			}
		}
	}
	theta = 10000;
	time_step = 25000;
	delay = 10;
	particle_radius = 0.005;
	width = 1e5*AU;
	ignore_width = 1e4*AU;		//ignore with must be set to value below the width
	debug_lower_bound = 0;
	debug_upper_bound = 20;
}

void cube_test(){
	cube(-100*AU,-10*AU);
	cube(100*AU,10*AU);
}
#endif
