#include <bits/stdc++.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/gl.h>
#include "science.h"	//includes the science file written for utility classes and functions
#include "prototype.h"	//testing file data including n body simulations
#include "console.h"	//console and ui features
	
using namespace std;

bool internal(Node *focus){
		//returns if the node is internal
	if(focus -> a == NULL && focus -> b == NULL && focus -> c == NULL && focus -> d == NULL
    && focus -> e == NULL && focus -> f == NULL && focus -> g == NULL && focus -> h == NULL)
		return false;
	return true;
}

void insert_particle(Particle *body, Node *focus);
void expand_octant(Node *octant){
		//initialize 8 octants for the given node
	double width = octant -> width, height = octant -> height, length = octant -> length;
	double x = octant -> x, y = octant -> y, z = octant -> z;
		//dimension and position management
	octant -> a = new Node( width/2, height/2, length/2, x - width/4, y+height/4, z+length/4);
	octant -> b = new Node( width/2, height/2, length/2, x + width/4, y+height/4, z+length/4);
	octant -> c = new Node( width/2, height/2, length/2, x - width/4, y-height/4, z+length/4);
	octant -> d = new Node( width/2, height/2, length/2, x + width/4, y-height/4, z+length/4);
	octant -> e = new Node( width/2, height/2, length/2, x - width/4, y+height/4, z-length/4);
	octant -> f = new Node( width/2, height/2, length/2, x + width/4, y+height/4, z-length/4);
	octant -> g = new Node( width/2, height/2, length/2, x - width/4, y-height/4, z-length/4);
	octant -> h = new Node( width/2, height/2, length/2, x + width/4, y-height/4, z-length/4);
}
	
void insert_to_octant(Particle *body, Node *focus){
	//choose the branch of the octree to traverse, this is ONLY a re-direction function
	if(abs(body -> x) > root -> width || abs(body -> y) > root -> height || abs(body -> z) > root -> length)
		cout <<"warning:: particle out of bounds root dimensions"<< endl;
			//particle has left bounds set by the simulation

	if(body -> x >= focus -> x){
		if(body -> y >= focus -> y){
			if(body -> z >= focus -> z)
				insert_particle(body,focus -> b);
			else
				insert_particle(body,focus -> f);
		}
		else{
			if(body -> z >= focus -> z)
				insert_particle(body,focus -> d);
			else 
				insert_particle(body,focus -> h);
		}
	}
	else{
		if(body -> y >= focus -> y){
			if(body -> z >= focus -> z)
				insert_particle(body,focus -> a);
			else
				insert_particle(body,focus -> e);
		}
		else{
			if(body -> z >= focus -> z)
				insert_particle(body, focus -> c);
			else
				insert_particle(body, focus -> g);
		}
	}
}

void insert_particle(Particle *body, Node *focus){
	if(focus -> p == NULL && !internal(focus)){
		focus -> p = body;
	}
	else if(internal(focus)){
		insert_to_octant(body,focus);
		focus -> update_mass();
			//after recursively added the body, we need to update te mass from the bottom of the tree
	}
	else{
		expand_octant(focus);
			//creates eight branches that will represent the expansion
		Particle *body_a = body;
		Particle *body_b = focus -> p;			
			//reference two seperate particles and add in recursively
		focus -> p = NULL;
			//focus pointer points to NULL, body_b has reference to original particle
		insert_to_octant(body_a,focus);
		insert_to_octant(body_b,focus);
		focus -> update_mass();
			//update mass;
	}
}

void build_tree(double w, double h, double l){
		//w,h,l must be the sane in magnitude
	delete root;
		//root and its subtrees will be deleted recursively 
	root = new Node(w,h,l,0,0,0);
	for(int i = 0; i < field.size(); i++)
		insert_particle(field.at(i),root);	
			//add particles in the field to the octree
}

bool in_octant(Node *n, Particle *p){
	//checks if the given aprticle is in the given node, if so it is advised to delete mass
	return 
	   p -> x <= n -> x + n -> width/2 && 
	   p -> x >= n -> x - n -> width/2 && 
	   p -> y <= n -> y + n -> width/2 && 
	   p -> y >= n -> y - n -> width/2 && 
	   p -> z <= n -> z + n -> width/2 && 
	   p -> z >= n -> z - n -> width/2;
}

void particle_update_traversal(Node *focus, Particle *body);
void particle_update(){
	#pragma omp parallel for num_threads(2)
	for(int i = 0; i < field.size(); i++){
		Particle *p = field.at(i);
		if(sqrt(pow(p->x-root->x,2) + pow(p->y-root->y,2) + pow(p->z-root->z,2)) > ignore_width && ignore_width > 0){
				//particle is out of focus and will no longer be visualized or calculated
			p -> ignore = true;
			field.erase(field.begin()+i);	//remove this line if particles should not be removed for crossing bounds
			continue;
		}
		if(p -> black_hole)
			continue;
		delete p -> f;
		delete p -> a;
		p -> f = new Vect(0,0,0);
		particle_update_traversal(root,p);
			//traverse the tree to calculate net force of particle p
		p -> a = new Vect(p -> f -> i, p -> f -> j, p -> f -> k);
			//acceleration is in the same direction as the net force
		p -> a -> scalar( (p -> f -> m / p -> m) * time_step);
			//force magnitude / mass multiplied by the timestep, TS
		double v_i = p -> v -> i, 
			   v_j = p -> v -> j, 
			   v_k = p -> v -> k;
		delete p -> v;
		Vect *curr_vel = new Vect(v_i,v_j,v_k);
		p -> v = Vect::add(curr_vel, p -> a);
		delete curr_vel;
			//vector addition of copy velocity and acceleration
		p -> x = p -> v -> i * time_step + p -> x; 
		p -> y = p -> v -> j * time_step + p -> y;
		p -> z = p -> v -> k * time_step + p -> z;
			//update the velocity and displacement vect
	}
}

void particle_update_traversal(Node *focus, Particle *body){
	if(!internal(focus) && focus -> p != NULL && focus -> p != body){
		if(body -> x != focus -> p -> x || body -> y != focus -> p -> y || body -> z != focus -> p -> z){
			Vect *gravity = new Vect(focus -> p -> x - body -> x, focus -> p -> y - body -> y, focus -> p -> z - body -> z);
			gravity -> scalar( 
			    ( GC * body -> m * focus -> p -> m ) /
				( pow(focus -> p -> x - body -> x, 2) + pow(focus -> p -> y - body -> y, 2) + pow(focus -> p -> z - body -> z, 2) )
			);
				//newtonian universal law of gravity
			Vect *temp = new Vect(body -> f -> i, body -> f -> j, body -> f -> k);
			delete body -> f;
			body -> f = Vect::add(temp, gravity);
			delete temp;
			delete gravity;
				//garbage collection
		}
	}
	else if(
			((focus -> width) / 
			sqrt(
				pow(body -> x - focus -> x, 2) +
				pow(body -> y - focus -> y, 2) +
				pow(body -> z - focus -> z, 2)
			)) < theta)
		{
			// if s/d < theta, use center of mass and total mass aproximations 
		double temp_mass = focus -> m;
		double temp_cmx = focus -> cmx;
		double temp_cmy = focus -> cmy;
		double temp_cmz = focus -> cmz;
			//take a temporary total mass of octant for newtonian calculation below
		if(in_octant(focus,body)){
			temp_mass = focus -> m - body -> m;
			temp_cmx = (temp_cmx*focus -> m - body -> x * body -> m)/temp_mass;
			temp_cmy = (temp_cmy*focus -> m - body -> y * body -> m)/temp_mass;
			temp_cmz = (temp_cmz*focus -> m - body -> z * body -> m)/temp_mass;
				//removing influence of body on node focus center of mass (com) in case on body in octant
		}
		//if the octant contains the body, changes must be made to ensure the center of mass and mass sum does not already include the particle body
		Vect *gravity = new Vect(temp_cmx - body -> x, temp_cmy - body -> y, temp_cmz - body -> z);
		gravity -> scalar( 
		    ( GC * body -> m * temp_mass ) /
			( pow(temp_cmx - body -> x, 2) + pow(temp_cmy - body -> y, 2) + pow(temp_cmz - body -> z, 2) )
		);
		Vect *temp = new Vect(body -> f -> i, body -> f -> j, body -> f -> k);
			//temp variable points to a copy of foce vec
		delete body -> f;
		body -> f = Vect::add(temp, gravity);
		delete temp;
		delete gravity;
	}
	else{
		if(focus -> a != NULL)
			particle_update_traversal(focus -> a, body);
		if(focus -> b != NULL)
			particle_update_traversal(focus -> b, body);
		if(focus -> c != NULL)
			particle_update_traversal(focus -> c, body);
		if(focus -> d != NULL)
			particle_update_traversal(focus -> d, body);
		if(focus -> e != NULL)
			particle_update_traversal(focus -> e, body);
		if(focus -> f != NULL)
			particle_update_traversal(focus -> f, body);
		if(focus -> g != NULL)
			particle_update_traversal(focus -> g, body);
		if(focus -> h != NULL)
			particle_update_traversal(focus -> h, body);
		//recursively visit more branches to update force from
	}
}

void display_text(string s, float x, float y){
		//gl display text s coordinates x and y
    glRasterPos2f(x, y);
	void * font = GLUT_BITMAP_HELVETICA_12;
	for (string::iterator i = s.begin(); i != s.end(); ++i){
		char c = *i;
		glColor3d(1.0, 0.0, 0.0);
		glutBitmapCharacter(font, c);
	}
}

void display_debug_tree_traversal_util(Node *sub_root, int layer_lower_bound, int layer_upper_bound, int layer_current){
		//visualizes the octree
	if(sub_root == NULL || layer_current > layer_upper_bound)
		return;
	if(layer_current <= layer_upper_bound && layer_current >= layer_lower_bound){
		glVertex3f((sub_root -> x + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z + sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z + sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z + sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z + sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z + sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z + sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z + sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z + sub_root -> width / 2) / (root -> width * zoom));
			//front face
		glVertex3f((sub_root -> x + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z - sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z - sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z - sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z - sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z - sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z - sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z - sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z - sub_root -> width / 2) / (root -> width * zoom));
			//back face
		glVertex3f((sub_root -> x + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z + sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z - sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z + sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z - sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z + sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y + sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z - sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z + sub_root -> width / 2) / (root -> width * zoom));
		glVertex3f((sub_root -> x - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> y - sub_root -> width / 2) / (root -> width * zoom), 
				   (sub_root -> z - sub_root -> width / 2) / (root -> width * zoom));
			//sides
	}
	display_debug_tree_traversal_util(sub_root -> a, layer_lower_bound, layer_upper_bound, layer_current+1);
	display_debug_tree_traversal_util(sub_root -> b, layer_lower_bound, layer_upper_bound, layer_current+1);
	display_debug_tree_traversal_util(sub_root -> c, layer_lower_bound, layer_upper_bound, layer_current+1);
	display_debug_tree_traversal_util(sub_root -> d, layer_lower_bound, layer_upper_bound, layer_current+1);
	display_debug_tree_traversal_util(sub_root -> e, layer_lower_bound, layer_upper_bound, layer_current+1);
	display_debug_tree_traversal_util(sub_root -> f, layer_lower_bound, layer_upper_bound, layer_current+1);
	display_debug_tree_traversal_util(sub_root -> g, layer_lower_bound, layer_upper_bound, layer_current+1);
	display_debug_tree_traversal_util(sub_root -> h, layer_lower_bound, layer_upper_bound, layer_current+1);
}

void display_debug_tree(){
	display_debug_tree_traversal_util(root, debug_lower_bound, debug_upper_bound, 0);
}

struct sys_time{
		//utility container made for fps calculations
	public:
		static unsigned int frame_count;
		static double prev_time;
		static double curr_time;
		static void update();
};

unsigned int sys_time::frame_count = 1;
double sys_time::prev_time = 0;
double sys_time::curr_time = clock();

void sys_time::update(){
	sys_time::prev_time = sys_time::curr_time;
	sys_time::curr_time = clock();
}

void display_fps(float x, float y){
	//display the framerate. analyze performance of the algorithm
	glColor3f(0.99,0.99,0.99);
	sys_time::update();
	int fps = (int)( CLOCKS_PER_SEC / (sys_time::curr_time - sys_time::prev_time) + 0.001);
	display_text("framerate: "+to_string(fps),x,y);
}

void additive_blend(){
//TODO
}	//adds extra alpha blending based on the visualization algorithm

void display_func(){
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glMatrixMode(GL_PROJECTION); 
		//set matrix to modelview
	glRotatef(my*200,1.0,0.0,0.0);
    glRotatef(mx*200,0.0,1.0,0.0);
    	//rotations
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		//clear
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
		//additive blending if required
	if(debug_tree){
    	glBegin(GL_LINES);
	    	glColor3f(0.2,0.9,0.3);
	    	display_debug_tree();
    		//debug mode is enabled
    	glEnd();
    }
    else{
	    for(int i = 0; i < field.size(); i++){
	    	Particle *render = field.at(i);
	    	if(render -> color[0] < 0.001 && render -> color[1] < 0.001 && render -> color[2] < 0.001)
	    		continue;
	    	float x1,y1,z1;
			x1 = render -> x / (root -> width * zoom);	
			y1 = render -> y / (root -> width * zoom);
			z1 = render -> z / (root -> width * zoom);
				//divided by camera dimensions
			glColor4f(render -> color[0], render -> color[1], render -> color[2], alpha);
			if(render -> ignore)
				glColor4f(0,0,0,0);
			glPushMatrix ();
		        glTranslatef    (x1, y1, z1);
		        glutSolidSphere (particle_radius, render_quality, render_quality);
	   		glPopMatrix ();
	    }
	}
    //glDisable(GL_BLEND);
	glLoadIdentity();
		//set matrix back to default when displaying text
	glColor3f(0.99,0.99,0.99);
	if(time_unit.compare("d") == 0)
    	display_text("time: "+to_string((unsigned long long)(T*SD))+" days", 0.65, 0.95);
    else if(time_unit.compare("y") == 0)
    	display_text("time: "+to_string((unsigned long long)((T*SD)/365.25))+" years", 0.65, 0.95);
    else if(time_unit.compare("s") == 0)
    	display_text("time: "+to_string(T)+" seconds", 0.65, 0.95);
    else
    	display_text("time: unknown unit", 0.65, 0.95);
    glColor3f(0.99,0.99,0.99);
    display_text("view res: "+to_string((int)(zoom*100))+"%", -0.98, 0.96);
    	//view 
    glColor3f(0.99,0.99,0.99);
    display_text("treenodes: "+to_string(OCTANT), 0.65, 0.92);
    glColor3f(0.99,0.99,0.99);
    display_text("particles: "+to_string(N), 0.65, 0.89);
    glColor3f(0.99,0.99,0.99);
    display_text("vec objects: "+to_string(VECTOR), 0.65, 0.86);
    	//printing some essential information
    glColor3f(0.99,0.99,0.99);
    if(debug_mode)
    	display_text("debug mode: true", 0.65, 0.83);
    else
    	display_text("debug mode: false", 0.65, 0.83);
  		//display debug status
    display_fps(0.65,0.8);
    glColor3f(0.99,0.99,0.99);
    display_text(get_display_text(), -0.9,-0.95);
    glutSwapBuffers();
    	//swap the second buffer
}

void physics(int data){
	if(field.size() < 0)
		return;
	glutTimerFunc(delay, physics, -1); 
	build_tree(width, width, width);
    particle_update();
    T += time_step;
    glutPostRedisplay();
}

void scroll_func(int button, int state, int x, int y){
   	if ((button == 3) || (button == 4)){
		if (button == 3 && state == GLUT_UP) 
			zoom *= 1.1;
		if (state == GLUT_DOWN && button == 4)  
			zoom /= 1.1;
   	} 
   	if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
   		locked = !locked;
}

void passive_func(int x, int y){
    //passive_func init tools
	if(!locked){
		mx = x;
	    my = y;
	    my = (my/glutGet(GLUT_WINDOW_WIDTH)) - 0.5f;
	    mx = -((mx/glutGet(GLUT_WINDOW_HEIGHT)) - 0.5f);
		//converting into usable coordinates	
	}
} 

void reshape(int x, int y){
	//reshaping the window
    if (y == 0 || x == 0) return;   
    glMatrixMode(GL_PROJECTION);  
    glLoadIdentity(); 
    gluPerspective(39.0,(GLdouble)x/(GLdouble)y,0.6,21.0);
    glMatrixMode(GL_MODELVIEW);
    glViewport(0,0,x,y);  
}

int main(int argc, char **argv){
	srand(time(NULL)); 
	user_def();
	if(field.size() >= 0)
		build_tree(width, width, width);
		//build octree 
	glutInit(&argc,argv);			
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(700,30);
	glutInitWindowSize(980,980);
	glutCreateWindow("astrocom");
	glutDisplayFunc(display_func);
	glutPassiveMotionFunc(passive_func);
	glutMouseFunc(scroll_func);
	glutTimerFunc(delay, physics, -1); 
	glutKeyboardFunc(do_something);
	glutSpecialFunc(glut_special_func);
	glutReshapeFunc(reshape);
	glutMainLoop();
	return 0;
}