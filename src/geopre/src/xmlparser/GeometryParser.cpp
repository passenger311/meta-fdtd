#include "GeometryParser.h"

void GeometryParser::startElement(const char* name, const char** atts)
{
	if (sScene->exFileReading != NULL)
		return;
	try {
		static int elemcounter = 0;
		bool bObjLink = false;
		elemcounter++;
		dictmap attrmap;
		  
		if (atts) {
			int i;
			for (i=0; atts[i]; i+=2) {
				attrmap[atts[i]] = atts[i+1];
			}
		}
		CObject* newobj = NULL;
		ArgumentReader* reader = sScene->ptReader;
		vec3 defpos(0,0,0);
		vec3 defscale(1,1,1);
		double defradius = 1;
		double defheight = 1;
		double defscalefactor = 1;
		if ( strcmp(name,"Sphere") == 0) {
			vec3 pos = reader->readVector3(attrmap["position"],&defpos);
			double rad = reader->readDouble(attrmap["radius"],&defradius);
			newobj = new CSphere(pos,rad);
		} else if ( strcmp(name,"Cylinder") == 0) {
			double height = reader->readDouble(attrmap["height"],&defheight);
			double radius = reader->readDouble(attrmap["radius"],&defradius);
			newobj = new CCylinder(radius,height);
		} else if ( strcmp(name,"Box") == 0) {
			vec3 startpos = reader->readVector3(attrmap["position"],&defpos);
			vec3 endpos = reader->readVector3(attrmap["endpoint"],&defpos);
			if (attrmap["size"] != NULL) {
				vec3 size = reader->readVector3(attrmap["size"]);
				endpos = startpos + size;
			}
			frame box(startpos,endpos);
			newobj = new CBox(box);
		} else if ( strcmp(name,"Ellipsoid") == 0) {
			vec3 pos = reader->readVector3(attrmap["position"],&defpos);
			vec3 scale = reader->readVector3(attrmap["scale"],&defscale);
			newobj = new CEllipsoid(pos,scale[VX],scale[VY],scale[VZ]);
		} else if ( strcmp(name,"ConvexPrism") == 0) {
			double height = reader->readDouble(attrmap["height"],&defheight);
			vector<vec2>* points = reader->readPointlist2(attrmap["points"]);
			newobj = new CConvexPrism(*points,height);
			delete points;
		} else if ( strcmp(name,"SimplePrism") == 0) {
			double height = reader->readDouble(attrmap["height"],&defheight);
			vector<vec2>* points = reader->readPointlist2(attrmap["points"]);
			newobj = new CSimplePrism(*points,height);
			delete points;
		} else if ( strcmp(name,"BezierPrism") == 0) {
			double height = reader->readDouble(attrmap["height"],&defheight);
			int defsteps = 10;
			int steps = reader->readInteger(attrmap["steps"],&defsteps);
			vector<vec2>* points = reader->readPointlist2(attrmap["points"]);
			newobj = new CBezierPrism(*points,height,steps);
			delete points;
		} else if ( strcmp(name,"SimpleTransform") == 0) {
			vec3 trans = reader->readVector3(attrmap["translation"],&defpos);
			vec3 defrotaxis(0,0,1);
			vec3 rotaxis = reader->readVector3(attrmap["rotationaxis"],&defrotaxis);
			double scale = reader->readDouble(attrmap["scale"],&defscalefactor);
			double defrotangle = 0;
			double rotang = reader->readDouble(attrmap["rotationangle"],&defrotangle);
			newobj = new CSimpleTransform(NULL,trans,rotaxis,rotang,scale);
		} else if ( strcmp(name,"SimpleRotationZ") == 0) {
			newobj = new CSimpleRotationZ();
		} else if ( strcmp(name,"LogicAndNot") == 0) {
			newobj = new CLogicAndNotObject(NULL, NULL);
		} else if ( strcmp(name,"LogicAnd") == 0) {
			newobj = new CLogicAndObject(NULL, NULL);
		} else if ( strcmp(name,"LogicOr") == 0) {
			newobj = new CLogicOrObject(NULL, NULL);
		} else if ( strcmp(name,"LogicXOr") == 0) {
			newobj = new CLogicXOrObject(NULL, NULL);
		} else if ( strcmp(name,"ObjectLink") == 0) {
			newobj = reader->readObjectReference(attrmap["ref"]);
			bObjLink = true;
		} else {
			// TODO: error handling
			printf ("Invalid object: %s\n", name);
			movCurrentObjects.push_back(NULL);
			return;
		}
		if (mDepth > 1) {
			if (!movCurrentObjects.back()->addSubObject(newobj))
				// TODO: error handling
				printf ("No more sub objects allowed for \"%s\"\n", movCurrentObjects.back()->name.c_str());
		} else if (mDepth == 1) {
			sScene->objects.push_back(newobj);
		}
		if (!bObjLink) {
			if (attrmap["id"] != NULL)	
				newobj->name = reader->readIdentifier(attrmap["id"]);
			else {
				char oname[50];
				sprintf(oname, "object%s%i", name, elemcounter); 
				newobj->name = oname;
			}
		}
		momObjects[newobj->name.c_str()] = newobj;
#ifdef DEBUG
		printf ("########start %s (%i)\n",name,mDepth);
		if (attrmap.size() > 0) {
			for (dictmap::iterator iter = attrmap.begin(); iter != attrmap.end(); iter++) {
				if (iter->second != NULL)
					printf("attribute: %s = \"%s\"\n", iter->first, iter->second);	
			}	
		}
		printf("%s\n", newobj->name.c_str());
		printf ("########end\n");
#endif
		movCurrentObjects.push_back(newobj);
	} catch (Exception* ex) {
		sScene->exFileReading = ex;
	}
}

//CObject* 
//GeometryParser::getCObject(const char* name, dictmap attributes)
//{
//	if ( strcmp(name,"Sphere") == 0) {
//		if (attributes["position"] == NULL || attributes["radius"] == NULL) {
//			printf ("Attribute(s) missing in object \"%s\"\n", name);
//			return NULL;	
//		}
//		vec3 pos = readVector(attributes["position"]);
//		double rad = atof(attributes["radius"]);
//		return new CSphere(pos,rad);
//	} else {
//		return NULL;	
//	}
//}

void 
GeometryParser::endElement(const char*)
{
	movCurrentObjects.pop_back();
}


void 
GeometryParser::charData(const XML_Char *s, int len)
{
#ifdef DEBUG
  const int leadingSpace = skipWhiteSpace(s);
  if (len==0 || len==leadingSpace)
  	return;  // called with whitespace between elements

/* write out the user data bracketed by ()*/
  putchar('(');
  fwrite(s, len, 1, stdout);
  puts(")");
#endif
}
