#include "GridBoxParser.h"
#include "../errorclasses.h"

GridBoxParser::GridBoxParser(expatppNesting *parent)
 : expatppNesting(parent)
{
	m_bDefaultGridSet = false;
}

GridBoxParser::~GridBoxParser()
{
}


//TODO: error numbering
void GridBoxParser::startElement(const char* name, const char** atts)
{
	if (sScene->exFileReading != NULL)
		return;
	try {
		dictmap attrmap;
		if (atts) {
			int i;
			for (i=0; atts[i]; i+=2) {
				attrmap[atts[i]] = atts[i+1];
			}
		}
		ArgumentReader* reader = sScene->ptReader;
		if ( strcmp(name,"Grid") == 0) {
			if (m_bDefaultGridSet && attrmap["id"] == NULL) {
				throw new ValueParseException(name, "Grid-object has no id attribute", 4001);
			} else if (!m_bDefaultGridSet) {
				attrmap["id"] = "default";
				m_bDefaultGridSet = true;
			}
			string gid = attrmap["id"]; 
			if (attrmap["points"] == NULL)
				throw new ValueParseException(gid, "Grid-object has no points attribute", 4001);
			if (attrmap["size"] == NULL)
				throw new ValueParseException(gid, "Grid-object has no size attribute", 4001);
			sScene->selectGrid(gid,false);
			GridBox* grid = sScene->ptCurrentGrid;
			vec3 size = reader->readVector3(attrmap["size"]);
			if (attrmap["origin"] != NULL)
				grid->frBBox.position_start = reader->readVector3(attrmap["origin"]);
			grid->frBBox.position_end = grid->frBBox.position_start + size;
			grid->frBBox.calcDimensions();
			vec3 points = reader->readVector3(attrmap["points"]);
			if (points[VX] < 1 || points[VY] < 1 || points[VZ] < 1)
				throw new ValueParseException(attrmap["points"],"invalid value for points parameter, provide three positive integer values",4005);
			grid->iCellsX = (int) points[VX]; 	
			grid->iCellsY = (int) points[VY]; 	
			grid->iCellsZ = (int) points[VZ]; 	
			if (attrmap["subgridding"] != NULL) {
				if ( strcmp(attrmap["subgridding"],"ENABLED") == 0) {
					grid->bAlwaysSubgridding = false;
					grid->bNoSubgridding = false;
				} else if ( strcmp(attrmap["subgridding"],"DISABLED") == 0) {
					grid->bAlwaysSubgridding = false;
					grid->bNoSubgridding = true;
				} else if ( strcmp(attrmap["subgridding"],"ALWAYS") == 0) {
					grid->bAlwaysSubgridding = true;
					grid->bNoSubgridding = false;
				} else
					throw new ValueParseException(attrmap["subgridding"],"invalid value for subgridding parameter",4005);
				
			}
			if (attrmap["subgriddingdivisions"] != NULL) {
				vec3 sgdivs = reader->readVector3(attrmap["subgriddingdivisions"]);
				if (sgdivs[VX] < 1 || sgdivs[VY] < 1 || sgdivs[VZ] < 1)
					throw new ValueParseException(attrmap["subgriddingdivisions"],"invalid value for subgriddingdivisions parameter, provide three positive integer values",4005);
				grid->iSubGriddingDivX = (int) sgdivs[VX]; 	
				grid->iSubGriddingDivY = (int) sgdivs[VY]; 	
				grid->iSubGriddingDivZ = (int) sgdivs[VZ]; 	
			}
			if (attrmap["yeegrid"] != NULL) {
				bool btrue = true;
				grid->bYeeGrid = reader->readBoolean(attrmap["yeegrid"],&btrue);
			}
			if (attrmap["pointframescale"] != NULL) {
				grid->dPointframeScale = reader->readDouble(attrmap["pointframescale"]);
			}
			if (attrmap["pointframedivisions"] != NULL) {
				vec3 pfdivs = reader->readVector3(attrmap["pointframedivisions"]);
				if (pfdivs[VX] < 1 || pfdivs[VY] < 1 || pfdivs[VZ] < 1)
					throw new ValueParseException(attrmap["pointframedivisions"],"invalid value for pointframedivisions parameter, provide three positive integer values",4005);
				grid->iPointframeDivisionsX = (int) pfdivs[VX]; 	
				grid->iPointframeDivisionsY = (int) pfdivs[VY]; 	
				grid->iPointframeDivisionsZ = (int) pfdivs[VZ]; 	
			}
			if (attrmap["pointframefacepoints"] != NULL) {
				grid->bPointframeFacepoints = reader->readBoolean(attrmap["pointframefacepoints"]);
			}
#ifdef DEBUG
			printf ("########start %s (%i)\n",name,mDepth);
			if (attrmap.size() > 0) {
				for (dictmap::iterator iter = attrmap.begin(); iter != attrmap.end(); iter++) {
					if (iter->second != NULL)
						printf("attribute: %s = \"%s\"\n", iter->first, iter->second);	
				}	
			}
			printf("%s\n", grid->sName.c_str());
			printf ("########end\n");
#endif
		} else {
			throw new FileParsingException(name, "invalid object", 4001);
		}
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
GridBoxParser::endElement(const char*)
{

}


void 
GridBoxParser::charData(const XML_Char *s, int len)
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
