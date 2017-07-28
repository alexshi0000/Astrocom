#ifndef science				//this is our global header file
#define science
#include <bits/stdc++.h>
using namespace std;

//mathematical
#define GC (6.673e-11)		//gravitational constant
#define PI (3.14159265359)	//pi

//conversions
#define SD (1.0/86400.0)	//seconds to days
#define PC (3.086e+16)		//parsec
#define AU (1.496e+11) 		//one astronomical unit in meters(m)
#define LY (9.461e+15)		//one lightyear
#define LC (299792458)		//speed of light in m/s

//default settings
#define DTS (25000) 		//default time per calculation
#define DVR (1)				//default view zoom
#define DFW (20*AU)			//default width
#define DPR (0.003)			//default particle radius
#define DIW (-1)			//default igonre width	
#define DLB (0)             //default lower bound for recursion tree
#define DUB (1000)          //default upper bound for recursion tree
#define DLY	(30)			//default delay
#define DTH (0.5)   		
							//	default opening angle constant, default value is 0.5 because the particles are sufficiently far enough to not be included 
							//	in APROXIMATION from node center of mass. However, TH -> MAX may still yield an accurate theoretical result
#define DTU "s"				//default time unit
#define DRQ (5)
#define DAV (1.0)			//default alpha value

typedef unsigned long long ull;

class Vect;
class Particle;
class Node;
vector <Particle*> field;

Node *root = nullptr; 		//the parent of all octants in 3d space
ull N      = 0; 			//N bodies
ull T      = 0; 			//time of simulation in seconds
ull OCTANT = 0;
ull VECTOR = 0;
float mx;					//mouse x and y positions
float my;
bool locked     = false;	//locked camera
bool debug_tree = false;

//misc
float  zoom              = DVR;
float  render_quality    = DRQ;
float  alpha             = DAV;
double time_step         = DTS;
double theta             = DTH;
double width             = DFW;
double ignore_width      = DIW;
double particle_radius   = DPR;
int    delay             = DLY;
int    debug_lower_bound = DLB;
int    debug_upper_bound = DUB;
string time_unit         = DTU;

class Point {
	//point in 4d space
	public:
		double x, y, z, w;
		Point(double x, double y, double z){
			this -> x = x;
			this -> y = y;
			this -> z = z;
		}
		void rotate(double u, double v, double w, double angle){
			//rotate on axis defined by the vector units i,j,k representing x,y,z axis
			double x = this -> x, y = this -> y, z = this -> z;
			double cos_p = cos(angle); 
			double sin_p = sin(angle);
			double x_prime = u*(u*x + v*y + w*z)*(1-cos_p) + x*cos_p + (-w*y + v*z) * sin_p;
			double y_prime = v*(u*x + v*y + w*z)*(1-cos_p) + y*cos_p + ( w*x - u*z) * sin_p;
			double z_prime = w*(u*x + v*y + w*z)*(1-cos_p) + z*cos_p + (-v*x + u*y) * sin_p;
			this -> x = x_prime;
			this -> y = y_prime;
			this -> z = z_prime;
		}
};

class Vect{
	public:
		double i, j, k, m;
		~Vect(){
			VECTOR--;
		}
		Vect(double i, double j, double k){
				//r3 implementation
			this -> i = i;
			this -> j = j;
			this -> k = k;
			this -> m = sqrt(i*i + j*j + k*k);
			VECTOR++;
		}
		void scalar(double m){
				//scaler multiplication of a non-zero directional vec
			if(this -> m == 0){
				//cout << "warning:: cannot apply scalar multiplication to a vector of magnitude 0" << endl;
				//prevents the multiplication of ZERO vectors and is a added safety to precision issues involving coordinates
				return;
			}
			double coeff = m / this -> m;
			this -> m = m;
			this -> i = this -> i * coeff;
			this -> j = this -> j * coeff;
			this -> k = this -> k * coeff;
		}
		static Vect* add(Vect *v, Vect *u){
				//addition of two vectors
			Vect *result = new Vect(v -> i + u -> i, v -> j + u -> j, v -> k + u -> k);
			return result;
		}
		static Vect* unit(Vect *u){
				//returns a unit vector
			Vect *unitVect = new Vect(u -> i, u -> j, u -> k);
			unitVect -> scalar(1);
			return unitVect;
		}
		static Vect* crossProduct(Vect *u, Vect *v){									
				//returns a vec pointer that is perpendicular to both u and v, the normal
			Vect *resultant = new Vect(0,0,0);
			resultant -> i = u -> j * v -> k - u -> k * v -> j;
			resultant -> j = v -> i * u -> k - u -> i * v -> k;
			resultant -> k = u -> i * v -> j - v -> i * u -> j;
			return resultant;
		}
		static double dotProduct(Vect *u, Vect *v){
				//dot operation
			return u -> i * v -> i + u -> j * v -> j + u -> k * v -> k;
		}
};

class Particle{
	public:
		bool black_hole;
		bool ignore;
		float color[3];
		double m,x,y,z;
		Vect *f,*v,*a;
		Particle(double m, double x, double y, double z, double i, double j, double k){
				//mass, x, y, z, the vector components
			black_hole = false;
			ignore     = false;
			this -> m = m;
			this -> x = x;
			this -> y = y;
			this -> z = z;
			this -> v = new Vect(i,j,k);
			this -> f = new Vect(0,0,0);
			this -> a = new Vect(0,0,0);
				//to be updated set to ZERO vector
			set_color(1,1,1);
			field.push_back(this);
			N++;
		}
		void set_color(float r, float g, float b){
			this -> color[0] = r;
			this -> color[1] = g;
			this -> color[2] = b;
		} 
};

class Node{
	public:
		double width, height, length, x, y, z, m, cmx, cmy, cmz;
			//since the node is representing a 3d quad, width, length and height are the same in magnitude
		Node *a, *b, *c, *d, *e, *f, *g, *h;
		Particle *p;
		Node(double width, double height, double length, double x, double y, double z){ 
			this -> width = width;
			this -> height = height;
			this -> length = length;
			this -> x = x;
			this -> y = y;
			this -> z = z;
			this -> m = 0;
			this -> cmx = 0;
			this -> cmy = 0;
			this -> cmz = 0;
			p = NULL, a = NULL, b = NULL, c = NULL, d = NULL, 
			e = NULL, f = NULL, g = NULL, h = NULL;
				//no octants formed nor any particle added 
			OCTANT++;
				//include ocatant count
		}	
		~Node(){
			OCTANT--;
				//remove octant count
			delete a, delete b, delete c, delete d, delete e, delete f, delete g, delete h;
			a = NULL, b = NULL, c = NULL, d = NULL, 
			e = NULL, f = NULL, g = NULL, h = NULL;
				//set to null to prevent undefined behaviour, pointers will be compared. subtrees are deleted recursively
		}
		void update_mass(){
			this -> m = 0, cmx = 0, cmy = 0, cmz = 0;
				//reset the total and center of masses
			if(a -> p == NULL){
				m += a -> m;
				cmx += a -> cmx * a -> m;
				cmy += a -> cmy * a -> m;
			}
			else{
				m += a -> p -> m;
				cmx += a -> p -> x * a -> p -> m;
				cmy += a -> p -> y * a -> p -> m;
			} // since this function is only called on an internal node, we do not need to check for existance of octants
			if(b -> p == NULL){
				m += b -> m;
				cmx += b -> cmx * b -> m;
				cmy += b -> cmy * b -> m;
			}
			else{
				m += b -> p -> m;
				cmx += b -> p -> x * b -> p -> m;
				cmy += b -> p -> y * b -> p -> m;
			}
			if(c -> p == NULL){
				m += c -> m;
				cmx += c -> cmx * c -> m;
				cmy += c -> cmy * c -> m;
			}
			else{
				m += c -> p -> m;
				cmx += c -> p -> x * c -> p -> m;
				cmy += c -> p -> y * c -> p -> m;
			}
			if(d -> p == NULL){
				m += d -> m;
				cmx += d -> cmx * d -> m;
				cmy += d -> cmy * d -> m;
			}
			else{
				m += d -> p -> m;
				cmx += d -> p -> x * d -> p -> m;
				cmy += d -> p -> y * d -> p -> m;
			}
			if(e -> p == NULL){
				m += e -> m;
				cmx += e -> cmx * e -> m;
				cmy += e -> cmy * e -> m;
			}
			else{
				m += e -> p -> m;
				cmx += e -> p -> x * e -> p -> m;
				cmy += e -> p -> y * e -> p -> m;
			}
			if(f -> p == NULL){
				m += f -> m;
				cmx += f -> cmx * f -> m;
				cmy += f -> cmy * f -> m;
			}
			else{
				m += f -> p -> m;
				cmx += f -> p -> x * f -> p -> m;
				cmy += f -> p -> y * f -> p -> m;
			}
			if(g -> p == NULL){
				m += g -> m;
				cmx += g -> cmx * g -> m;
				cmy += g -> cmy * g -> m;
			}
			else{
				m += g -> p -> m;
				cmx += g -> p -> x * g -> p -> m;
				cmy += g -> p -> y * g -> p -> m;
			}
			if(h -> p == NULL){
				m += h -> m;
				cmx += h -> cmx * h -> m;
				cmy += h -> cmy * h -> m;
			}
			else{
				m += h -> p -> m;
				cmx += h -> p -> x * h -> p -> m;
				cmy += h -> p -> y * h -> p -> m;
			}
			cmx = cmx / this -> m;
			cmy = cmy / this -> m; 
		}
};

#endif