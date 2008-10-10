#ifndef SCENE_H_
#define SCENE_H_

#include <vector>
#include "errorclasses.h"
#include "objects/CObject.h"
#include <string>
#include <map>

using namespace std;

class Scene
{
public:
	vector<CObject*> objects;

public:
	Scene() { }
	~Scene() { }
};

#endif /*SCENE_H_*/
