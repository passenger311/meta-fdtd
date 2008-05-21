#include "OutputParser.h"
#include "../errorclasses.h"
#include "../constants.h"

OutputParser::OutputParser(expatppNesting *parent)
 : expatppNesting(parent)
{
	m_bDefaultOLSet = false;
	m_bDefaultOFSet = false;
}

OutputParser::~OutputParser()
{
}


//TODO: error numbering
void OutputParser::startElement(const char* name, const char** atts)
{
	if (sScene->exFileReading != NULL)
		return;
	try {
		static int elemcounter = 0;
		elemcounter++;
		dictmap attrmap;
		if (atts) {
			int i;
			for (i=0; atts[i]; i+=2) {
				attrmap[atts[i]] = atts[i+1];
			}
		}
		ArgumentReader* reader = sScene->ptReader;
		if ( strcmp(name,"OutputList") == 0) {
			if (m_sCurrentOL != "") {
				throw new FileParsingException(name, "OutputList-objects are not allowed to have OutputList-objects as child-nodes", 5031);
			}
			if (attrmap["id"] == NULL) {
				if (m_bDefaultOLSet) 
					throw new ValueParseException(name, "OutputList-object has no id attribute", 4001);
				else {
					attrmap["id"] = "default";
				}
			}
			if (strcmp(attrmap["id"],"default")==0)
				m_bDefaultOLSet = true;
			m_sCurrentOL = attrmap["id"]; 
			sScene->selectOutputList(m_sCurrentOL,false);
			OutputList* ptol = sScene->ptCurrentOutputList;
			string nullstr = "";
			ptol->sParameter = reader->readIdentifier(attrmap["parameter"],&nullstr);
			ptol->sGrid = reader->readIdentifier(attrmap["grid"],&nullstr);
			ptol->sFormat = reader->readIdentifier(attrmap["format"],&nullstr);
			ptol->sObjectCut = reader->readIdentifier(attrmap["objectcut"],&nullstr);
			ptol->bYeeGrid_isset = (attrmap["yeegrid"] != NULL);
			bool btrue = true;
			ptol->bYeeGrid = reader->readBoolean(attrmap["yeegrid"],&btrue);
			ptol->bUseCluster_isset = (attrmap["cluster"] != NULL);
			bool bfalse = false;
			ptol->bUseCluster = reader->readBoolean(attrmap["cluster"],&bfalse);
#ifdef DEBUG
			printf ("########start %s (%i)\n",name,mDepth);
			if (attrmap.size() > 0) {
				for (dictmap::iterator iter = attrmap.begin(); iter != attrmap.end(); iter++) {
					if (iter->second != NULL)
						printf("attribute: %s = \"%s\"\n", iter->first, iter->second);	
				}	
			}
			printf("%s\n", ptol->sName.c_str());
			printf ("########end\n");
#endif
		} else if ( strcmp(name,"OutputFile") == 0) {
			if (mDepth > 2) {
				throw new FileParsingException(name, "OutputFile-objects must be nested in OutputList-objects", 5033);
			}
			char oname[50];
			if (attrmap["id"] == NULL) {
				if (mDepth > 1) {
					sprintf(oname, "object%s%i", name, elemcounter); 
					attrmap["id"] = oname; 
				} else if (m_bDefaultOFSet) { 
					throw new ValueParseException(name, "OutputFile-object has no id attribute", 4001);
				} else {
					attrmap["id"] = "default";
				}
			}
			if (strcmp(attrmap["id"],"default")==0)
				m_bDefaultOFSet = true;
			sScene->selectOutputFile(attrmap["id"],false);
			OutputList* ptol = sScene->ptCurrentOutputList;
			OutputFile* ptof = sScene->ptCurrentOutputFile;
			if (mDepth == 2)
				ptol->vOutputFiles.push_back(ptof);
//			string stt = ptol->sParameter;
			ptof->sFile = reader->readIdentifier(attrmap["file"]);
			ptof->sParameter = reader->readIdentifier(attrmap["parameter"],&(ptol->sParameter));
			ptof->sGrid = reader->readIdentifier(attrmap["grid"],&(ptol->sGrid));
			ptof->sFormat = reader->readIdentifier(attrmap["format"],&(ptol->sFormat));
			ptof->sObjectCut = reader->readIdentifier(attrmap["objectcut"],&(ptol->sObjectCut));
			bool byee = ptol->bYeeGrid_isset ? ptol->bYeeGrid : true;
			ptof->bYeeGrid = reader->readBoolean(attrmap["yeegrid"],&byee);
			bool bclus = ptol->bUseCluster_isset ? ptol->bUseCluster : false;
			ptof->bUseCluster = reader->readBoolean(attrmap["cluster"],&bclus);
#ifdef DEBUG
			printf ("########start %s (%i)\n",name,mDepth);
			if (attrmap.size() > 0) {
				for (dictmap::iterator iter = attrmap.begin(); iter != attrmap.end(); iter++) {
					if (iter->second != NULL)
						printf("attribute: %s = \"%s\"\n", iter->first, iter->second);	
				}	
			}
			printf("%s\n", ptof->sName.c_str());
			printf ("########end\n");
#endif
		} else {
			throw new FileParsingException(name, "invalid object", 5001);
		}
	} catch (Exception *ex) {
		sScene->exFileReading = ex;
	}
}

void 
OutputParser::endElement(const char* elname)
{
	if (strcmp(elname,"OutputList") == 0)
		m_sCurrentOL = "";
}


void 
OutputParser::charData(const XML_Char *s, int len)
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
