#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include "header.h"
#include "bitmap_image.hpp"
#include "config.h"

#define pi (2*acos(0.0))
#define ANGLE_OF_ROTATION 3.0       // in degrees

#define FLOOR_AMB 0.25
#define FLOOR_DIFF 0.25
#define FLOOR_SPEC 0.25
#define FLOOR_REC 0.25
#define FLOOR_SHININESS 10

using namespace std;

int drawaxes;

string input_filepath = "scene.txt";
string output_filepath = "output";

extern int recursion_level;
int pixels_along_dim = 0;
int object_count = 0;
int point_light_count = 0;
int spot_light_count = 0;

double viewAngle;
double windowHeight = 500;
double windowWidth = 500;

int bmp_image_count = 0;

extern vector<Object*> objects;
extern vector<PointLight> point_lights;
extern vector<SpotLight> spot_lights;

extern Vector3D eye;            /// position of the camera
Vector3D u;                     /// up direction
Vector3D r;                     /// right direction
Vector3D l;                     /// look direction

void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
}

void capture()
{
    bitmap_image img(pixels_along_dim, pixels_along_dim);

    double planeDistance = (windowHeight/2)/tan(viewAngle*pi/360);
    Vector3D topLeft = eye + l*planeDistance - r*(windowWidth/2) + u*(windowHeight/2);

    double du = (double) windowWidth/pixels_along_dim;
    double dv = (double) windowHeight/pixels_along_dim;

    topLeft = topLeft + r*(du/2) - u*(dv/2);

    int nearest;
    double t, tMin;

    for(int i=0; i<pixels_along_dim; i++)
    {
        for(int j=0; j<pixels_along_dim; j++)
        {
            Vector3D currentPixel = topLeft + r*(i*du) - u*(j*dv);
            Ray ray(eye, currentPixel - eye);       
            Color dummy;

            nearest = -1;
            tMin = INF;

            for(int i=0; i<objects.size(); i++)
            {
                t = objects[i]->intersect(ray, dummy, 0);

                if(t > 0 && t < tMin)
                {
                    tMin = t;
                    nearest = i;
                }
            }

            if(nearest != -1)
            {
                Color colorContainer;
                tMin = objects[nearest]->intersect(ray, colorContainer, 1);
                img.set_pixel(i, j, (int) round(colorContainer.red * 255), (int) round(colorContainer.green * 255), (int) round(colorContainer.blue * 255));
            }
            else
            {
                img.set_pixel(i, j, 0, 0, 0);
            }
        }
    }
    img.save_image(output_filepath + to_string(++bmp_image_count) + ".bmp");
    cout << "Image " << bmp_image_count << " captured" << endl;
}


void keyboardListener(unsigned char key, int x,int y)
{
	switch(key){
        case '0':
            capture();
            break;
		case '1':
            l = l*cos(ANGLE_OF_ROTATION * pi/180) + (u^l)*sin(ANGLE_OF_ROTATION * pi/180);
            r = l^u;
			break;
        case '2':
            l = l*cos(-ANGLE_OF_ROTATION * pi/180) + (u^l)*sin(-ANGLE_OF_ROTATION * pi/180);
            r = l^u;
			break;
        case '3':
            l = l*cos(ANGLE_OF_ROTATION * pi/180) + (r^l)*sin(ANGLE_OF_ROTATION * pi/180);
            u = r^l;
            break;
        case '4':
            l = l*cos(-ANGLE_OF_ROTATION * pi/180) + (r^l)*sin(-ANGLE_OF_ROTATION * pi/180);
            u = r^l;
            break;
        case '5':
            u = u*cos(ANGLE_OF_ROTATION * pi/180) + (l^u)*sin(ANGLE_OF_ROTATION * pi/180);
            r = l^u;
            break;
        case '6':
            u = u*cos(-ANGLE_OF_ROTATION * pi/180) + (l^u)*sin(-ANGLE_OF_ROTATION * pi/180);
            r = l^u;
            break;
		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y)
{
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			eye = eye - l;
			break;
		case GLUT_KEY_UP:		// up arrow key
			eye = eye + l;
			break;
		case GLUT_KEY_RIGHT:
			eye = eye + r;
			break;
		case GLUT_KEY_LEFT:
			eye = eye - r;
			break;

		case GLUT_KEY_PAGE_UP:
		    eye = eye + u;
			break;
		case GLUT_KEY_PAGE_DOWN:
		    eye = eye - u;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y)
{
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}


void display()
{

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();
	
    gluLookAt(eye.x, eye.y, eye.z,  eye.x + l.x, eye.y + l.y, eye.z + l.z,   u.x, u.y, u.z);

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects
	drawAxes();

    for(int i=0; i<objects.size(); i++){
        objects[i]->draw();
    }

    for(int i=0; i<point_lights.size(); i++){
        point_lights[i].draw();
    }

    for(int i=0; i<spot_lights.size(); i++){
        spot_lights[i].draw();
    }

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate()
{
	glutPostRedisplay();
}

void init()
{

    viewAngle = 80;

	//codes for initialization
	drawaxes=1;

    // initializing camera position
	eye.x = 100;
	eye.y = 100;
	eye.z = 0;

    // initializing up direction
	u.x = 0;
	u.y = 0;
	u.z = 1;

    // initializing down direction
	r.x = -1/sqrt(2);
	r.y = 1/sqrt(2);
	r.z = 0;

	// initializing look direction
	l = u ^ r;

	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(viewAngle,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

void loadData()
{
    ifstream sceneFile;
    sceneFile.open(input_filepath);

    if(!sceneFile.is_open())
    {
        cout << "Error opening file" << endl;
        exit(1);
    }

    sceneFile >> recursion_level;
    sceneFile >> pixels_along_dim;
    sceneFile >> object_count;

    string objectName;

    Object* obj = NULL;

    // read objects
    for(int i=0; i<object_count; i++)
    {
        sceneFile >> objectName;

        if(objectName == "sphere")
        {
            double center_x, center_y, center_z, radius;
            double color_r, color_g, color_b;
            double ambient, diffuse, specular, reflection, shininess;

            sceneFile >> center_x >> center_y >> center_z;
            sceneFile >> radius;
            sceneFile >> color_r >> color_g >> color_b;
            sceneFile >> ambient >> diffuse >> specular >> reflection;
            sceneFile >> shininess;

            Vector3D center(center_x, center_y, center_z);
            obj = new Sphere(center, radius);

            obj->setColor(color_r, color_g, color_b);
            obj->setCoefficients(ambient, diffuse, specular, reflection);
            obj->setShine(shininess);
        }
        else if(objectName == "triangle")
        {
            double x1, y1, z1, x2, y2, z2, x3, y3, z3;
            double color_r, color_g, color_b;
            double ambient, diffuse, specular, reflection, shininess;

            sceneFile >> x1 >> y1 >> z1;
            sceneFile >> x2 >> y2 >> z2;
            sceneFile >> x3 >> y3 >> z3;
            sceneFile >> color_r >> color_g >> color_b;
            sceneFile >> ambient >> diffuse >> specular >> reflection;
            sceneFile >> shininess;

            Vector3D a(x1, y1, z1);
            Vector3D b(x2, y2, z2);
            Vector3D c(x3, y3, z3);

            obj = new Triangle(a, b, c);
            obj->setColor(color_r, color_g, color_b);
            obj->setCoefficients(ambient, diffuse, specular, reflection);
            obj->setShine(shininess);
        }
        else if(objectName == "general")
        {
            double A, B, C, D, E, F, G, H, I, J;
            double cube_ref_X, cube_ref_Y, cube_ref_Z;
            double length, width, height;
            double color_r, color_g, color_b;
            double ambient, diffuse, specular, reflection, shininess;

            sceneFile >> A >> B >> C >> D >> E >> F >> G >> H >> I >> J;
            sceneFile >> cube_ref_X >> cube_ref_Y >> cube_ref_Z >> length >> width >> height;
            sceneFile >> color_r >> color_g >> color_b;
            sceneFile >> ambient >> diffuse >> specular >> reflection;
            sceneFile >> shininess;

            Vector3D cube_reference_point(cube_ref_X, cube_ref_Y, cube_ref_Z);
            
            obj  = new GeneralQuadricSurface(A, B, C, D, E, F, G, H, I, J, cube_reference_point, length, width, height);
            obj->setColor(color_r, color_g, color_b);
            obj->setCoefficients(ambient, diffuse, specular, reflection);
            obj->setShine(shininess);
        }

        else{
            cout << "Invalid Object" << endl;
            exit(1);
        }

        objects.push_back(obj);
    }

    sceneFile >> point_light_count;
    double light_x, light_y, light_z;
    double light_color_r, light_color_g, light_color_b;

    // read in point lights
    for(int i = 0; i < point_light_count; i++)
    {
        sceneFile >> light_x >> light_y >> light_z;
        sceneFile >> light_color_r >> light_color_g >> light_color_b;

        Vector3D light_pos(light_x, light_y, light_z);

        point_lights.push_back(PointLight(light_pos, light_color_r, light_color_g, light_color_b));
    }

    sceneFile >> spot_light_count;
    double spot_light_x, spot_light_y, spot_light_z;
    double spot_light_dir_x, spot_light_dir_y, spot_light_dir_z;
    double spot_light_color_r, spot_light_color_g, spot_light_color_b;
    double spot_light_angle;

    for(int i = 0; i < spot_light_count; i++)
    {
        sceneFile >> spot_light_x >> spot_light_y >> spot_light_z;
        sceneFile >> spot_light_color_r >> spot_light_color_g >> spot_light_color_b;
        sceneFile >> spot_light_dir_x >> spot_light_dir_y >> spot_light_dir_z;
        sceneFile >> spot_light_angle;

        Vector3D spot_light_pos(spot_light_x, spot_light_y, spot_light_z);
        Vector3D spot_light_dir(spot_light_dir_x, spot_light_dir_y, spot_light_dir_z);

        spot_lights.push_back(SpotLight(spot_light_pos, spot_light_dir, spot_light_color_r, spot_light_color_g, spot_light_color_b, spot_light_angle));
    }


    sceneFile.close();

    // creating floor
    obj = new Floor(1000, 20);
    obj->setCoefficients(FLOOR_AMB, FLOOR_DIFF, FLOOR_SPEC, FLOOR_REC);
    obj->setShine(FLOOR_SHININESS);
    objects.push_back(obj);
}

void clearObjects()
{
    for(int i = 0; i < objects.size(); i++)
    {
        delete objects[i];
    }
    objects.clear();
}

void clearPointLights()
{
    point_lights.clear();
}

void clearSpotLights()
{
    spot_lights.clear();
}

int main(int argc, char **argv)
{
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("Ray Tracing");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

    if((atexit(clearObjects) != 0) || (atexit(clearPointLights) != 0) || (atexit(clearSpotLights) != 0))
    {
        cout << "Error: Could not register exit function" << endl;
        exit(1);
    }

    loadData();

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
