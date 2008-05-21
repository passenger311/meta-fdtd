#define DEBUG

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "Scene.h"
#include "errorclasses.h"

using namespace std;

int main(int argc, char **argv) {

	Scene scene;
	try {
		scene.readPriorityCommandLineArguments(argc, argv);
		if (scene.bQuit) return 0;
	}
	catch (Exception* ex) {
		cerr << "Error while reading command-line arguments: ";
		cerr << ex->getMessage() << endl;
		cerr << "Use option --help for further information" << endl;
		delete ex;
		return 1;
	}
	try {
		scene.readSceneFile();
	}
	catch (Exception* ex) {
		cerr << "Error while reading scene-file: ";
		cerr << ex->getMessage() << endl;
		delete ex;
		return 1;
	}
	try {
		scene.readCommandLineArguments(argc, argv);
	}
	catch (Exception* ex) {
		cerr << "Error while reading command-line arguments: ";
		cerr << ex->getMessage() << endl;
		cerr << "Use option --help for further information" << endl;
		delete ex;
		return 1;
	}
	scene.buildScene();
	try {
		scene.generateOutput();
	}
	catch (Exception* ex) {
		cerr << "Error processing Output: ";
		cerr << ex->getMessage() << endl;
		delete ex;
		return 1;
	}

	return 0;
	
/*	CCylinder cyl(7,17);
	vec3 trslt(-15,-15,-10);
	vec3 rotax(1,1,1);
	double angle = 35;
	double sf = 1;
	CSimpleTransform trf(&cyl,trslt,rotax,angle,sf);

	vec3 pos(5,5,12);
	//frame boxfrm(pos,4,4,4);
	CEllipsoid box(pos,12,3,4);
	
	CCylinder cyl2(4,35);
	CLogicAndNotObject bobj(&cyl,&cyl2);
	
	scene.objects.push_back(&bobj);
	scene.objects.push_back(&trf);
	scene.objects.push_back(&box); */
	
//	vector<vec2> pl;
//	pl.push_back(vec2(4.0,4.0));
//	pl.push_back(vec2(7.0,0.0));
//	pl.push_back(vec2(6.0,-1.0));
//	pl.push_back(vec2(6.0,-2.0));
//	pl.push_back(vec2(5.0,-1.5));
//	pl.push_back(vec2(5.0,0.0));
//	pl.push_back(vec2(2.5,-1.5));
//	pl.push_back(vec2(4.0,-2.0));
//	pl.push_back(vec2(4.0,-3.0));
//	pl.push_back(vec2(1.0,-4.0));
//	pl.push_back(vec2(-4.0,4.0));
//	pl.push_back(vec2(0.0,2.0));
//	pl.push_back(vec2(1.0,2.0));
//	CSimplePrism test(pl,4.0);
//	vector<vec2> pl;
//	pl.push_back(vec2(2.0,0.0));
//	pl.push_back(vec2(2.0,-5.0));
//	pl.push_back(vec2(3.5,-3.0));
//	pl.push_back(vec2(3.5,-6.0));
//	pl.push_back(vec2(3.5,-10.0));
//	pl.push_back(vec2(-3.5,-10.0));
//	pl.push_back(vec2(-3.5,-6.0));
//	pl.push_back(vec2(-3.5,-3.0));
//	pl.push_back(vec2(-2.0,-5.0));
//	pl.push_back(vec2(-2.0,0.0));
//	CBezierPrism test(pl,4.0,10);
//	scene.objects.push_back(&test);
//
//	vec3 origin(0,0,-1);
//	vec3 endpoint(20,20,1);
//	frame gridframe(origin,endpoint);
//	
//	GridBox gridbox(&scene, gridframe, 21, 21, 1);
//	gridbox.iSubGriddingDivX = 10;
//	gridbox.iSubGriddingDivY = gridbox.iSubGriddingDivX;
//	gridbox.iSubGriddingDivZ = gridbox.iSubGriddingDivX;
//	
//	gridbox.bNoSubgridding = false;
//	gridbox.bAlwaysSubgridding = false;
//	
//	gridbox.bYeeGrid = true;
//	
//	FileHandler fhd;
//	fhd.sFile = "/home/tobias/gridsetup/output_test.vtk";
//	fhd.bYeeGrid = true;
//	fhd.gbGridBox = &gridbox;
//	
//	gridbox.generateOutput(&fhd);
//		
//	
//	return 0;
}
