#ifndef ALGEBRA3EXTENSION_H_
#define ALGEBRA3EXTENSION_H_

#include "algebra3aux.h"

#define SIGN(A) (A > 0) ? 1 : ((A < 0) ? -1 : 0)

class frame {
	protected:
	public:
		vec3 position_start;
		vec3 position_end;
		double mindimension;	
		double maxdimension;	
	
	public:
		frame() { }
		frame(vec3& pos_start,vec3& pos_end)
		{
			position_start = vec3(MIN(pos_start[VX],pos_end[VX]),MIN(pos_start[VY],pos_end[VY]),MIN(pos_start[VZ],pos_end[VZ]));
			position_end = vec3(MAX(pos_start[VX],pos_end[VX]),MAX(pos_start[VY],pos_end[VY]),MAX(pos_start[VZ],pos_end[VZ]));
			calcDimensions();	
		}
		frame(vec3* pointlist,int count)
		{
			getFromPointlist(pointlist, count);
		}
		frame(vec3& pos_start,double xsize,double ysize, double zsize)
		{
			vec3 pos_end = pos_start+vec3(xsize,ysize,zsize); 
			position_start = vec3(MIN(pos_start[VX],pos_end[VX]),MIN(pos_start[VY],pos_end[VY]),MIN(pos_start[VZ],pos_end[VZ]));
			position_end = vec3(MAX(pos_start[VX],pos_end[VX]),MAX(pos_start[VY],pos_end[VY]),MAX(pos_start[VZ],pos_end[VZ]));
			calcDimensions();	
		}
		
		inline void getFromPointlist(vec3* pointlist,int count)
		{
			position_start = pointlist[0];
			position_end = pointlist[0];
			for (int i=1; i<count; i++) {
				if (pointlist[i][VX] > position_end[VX])
					position_end[VX] = pointlist[i][VX]; 
				if (pointlist[i][VX] < position_start[VX])
					position_start[VX] = pointlist[i][VX]; 
				if (pointlist[i][VY] > position_end[VY])
					position_end[VY] = pointlist[i][VY]; 
				if (pointlist[i][VY] < position_start[VY])
					position_start[VY] = pointlist[i][VY]; 
				if (pointlist[i][VZ] > position_end[VZ])
					position_end[VZ] = pointlist[i][VZ]; 
				if (pointlist[i][VZ] < position_start[VZ])
					position_start[VZ] = pointlist[i][VZ]; 
			}
			calcDimensions();	
		}
		
		inline bool pointInFrame(vec3& point)
		{
			return (position_start[VX]<=point[VX] && point[VX]<=position_end[VX] &&
				position_start[VY]<=point[VY] && point[VY]<=position_end[VY] && 
				position_start[VZ]<=point[VZ] && point[VZ]<=position_end[VZ]); 
		}
		
		inline bool frameIntersection(const frame& ff)
		{
			return ((position_start[VX]<=ff.position_start[VX] && ff.position_start[VX]<=position_end[VX]) || 
				(position_start[VX]<=ff.position_end[VX] && ff.position_end[VX]<=position_end[VX]) ||
				(ff.position_start[VX]<=position_start[VX] && position_end[VX]<=ff.position_end[VX]))&&
((position_start[VY]<=ff.position_start[VY] && ff.position_start[VY]<=position_end[VY]) || 
				(position_start[VY]<=ff.position_end[VY] && ff.position_end[VY]<=position_end[VY]) ||
				(ff.position_start[VY]<=position_start[VY] && position_end[VY]<=ff.position_end[VY]))&&
((position_start[VZ]<=ff.position_start[VZ] && ff.position_start[VZ]<=position_end[VZ]) || 
				(position_start[VZ]<=ff.position_end[VZ] && ff.position_end[VZ]<=position_end[VZ]) ||
				(ff.position_start[VZ]<=position_start[VZ] && position_end[VZ]<=ff.position_end[VZ]));
		}
		
		inline bool frameIntersectionZ(const frame& ff)
		{
			return (position_start[VZ]<=ff.position_start[VZ] && ff.position_start[VZ]<=position_end[VZ]) || 
				(position_start[VZ]<=ff.position_end[VZ] && ff.position_end[VZ]<=position_end[VZ]) ||
				(ff.position_start[VZ]<=position_start[VZ] && position_end[VZ]<=ff.position_end[VZ]);
		}
		
		inline frame combineWith(frame& ff)
		{
			frame combined;
			combined.position_start[VX] = MIN(position_start[VX],ff.position_start[VX]);
			combined.position_start[VY] = MIN(position_start[VY],ff.position_start[VY]);
			combined.position_start[VZ] = MIN(position_start[VZ],ff.position_start[VZ]);
			combined.position_end[VX] = MAX(position_end[VX],ff.position_end[VX]);
			combined.position_end[VY] = MAX(position_end[VY],ff.position_end[VY]);
			combined.position_end[VZ] = MAX(position_end[VZ],ff.position_end[VZ]);
			return combined;
		}

		inline frame& operator = (const frame& ff)
		{ 
			position_start = ff.position_start;
			position_end = ff.position_end;
			return *this;
		}
		
		vec3* getEdgePoints()
		{
			vec3* pts = new vec3[8];
			pts[0] = position_start;
			pts[1] = position_start+vec3(position_end[VX]-position_start[VX],0,0);
			pts[2] = position_start+vec3(0,position_end[VY]-position_start[VY],0);
			pts[3] = position_start+vec3(position_end[VX]-position_start[VX],position_end[VY]-position_start[VY],0);
			pts[4] = position_start+vec3(0,0,position_end[VZ]-position_start[VZ]);
			pts[5] = position_start+vec3(position_end[VX]-position_start[VX],0,position_end[VZ]-position_start[VZ]);
			pts[6] = position_start+vec3(0,position_end[VY]-position_start[VY],position_end[VZ]-position_start[VZ]);
			pts[7] = position_end;
			return pts;
		}
		 
	public:
		inline void calcDimensions()
		{
			vec3 size(position_end-position_start);
			mindimension = MIN(MIN(size[VX],size[VY]),size[VZ]); 	
			maxdimension = MAX(MAX(size[VX],size[VY]),size[VZ]); 	
		}	
};

class pointframe: public frame
{
	public:
		vec3* m_vPoints;
		int m_iPointCount;
		vec3 movingPoint;
		
	public:
		pointframe()
		{	
			m_iPointCount = 0;
		}
		virtual ~pointframe()
		{
			if (m_iPointCount > 0)
				delete [] m_vPoints;
			m_iPointCount = 0;
		}
		pointframe(vec3& pstart, vec3& pend, int divX = 1, int divY = 1, int divZ = 1, bool facePoints = true)
		{
			m_iPointCount = 0;
			setupPointframe(pstart, pend, divX, divY, divZ, facePoints);
			
		}
		pointframe(frame& ff,int divX = 1, int divY = 1, int divZ = 1, bool facePoints = true)
		{
			m_iPointCount = 0;
			setupPointframe(ff.position_start, ff.position_end, divX, divY, divZ, facePoints);
		}
		
		void setupPointframe(vec3& pstart, vec3& pend, int divX = 1, int divY = 1, int divZ = 1, bool facePoints = true)
		{
			position_start = pstart;
			position_end = pend;
			movingPoint = 0.5*(position_start + position_end); 
			if (divX == 0 || divY == 0 || divZ == 0) {
				vec3 midPt;
				midPt[VX] = (position_end[VX]+position_start[VX])/2;
				midPt[VY] = (position_end[VY]+position_start[VY])/2;
				midPt[VZ] = (position_end[VZ]+position_start[VZ])/2;
				m_iPointCount = 0;
				movingPoint = midPt;
				calcDimensions();
				return;
			}
			double spaceX = (position_end[VX]-position_start[VX])/divX;
			double spaceY = (position_end[VY]-position_start[VY])/divY;
			double spaceZ = (position_end[VZ]-position_start[VZ])/divZ;
			if (m_iPointCount > 0)
				delete [] m_vPoints;
			m_iPointCount = 2*(divX+1)*(divY+1)+2*(divX+divY)*(divZ-1);
			if (facePoints)
				m_iPointCount += 2*divX*divY + 2*divX*divZ + 2*divY*divZ;
			vec3* pts = new vec3[m_iPointCount];
			int iP = 0;
			for (int iX = 0; iX <= divX; iX++) {
				for (int iY = 0; iY <= divY; iY++) {
					pts[iP++] = position_start+vec3(spaceX*iX,spaceY*iY,0);
					pts[iP++] = position_start+vec3(spaceX*iX,spaceY*iY,spaceZ*divZ);
					if (facePoints && iX < divX && iY < divY) {
						pts[iP++] = position_start+vec3(spaceX*(iX+0.5),spaceY*(iY+0.5),0);
						pts[iP++] = position_start+vec3(spaceX*(iX+0.5),spaceY*(iY+0.5),spaceZ*divZ);
					}
				}
				for (int iZ = 0; iZ < divZ; iZ++) {
					if (iZ > 0) {
						pts[iP++] = position_start+vec3(spaceX*iX,0,spaceZ*iZ);
						pts[iP++] = position_start+vec3(spaceX*iX,spaceY*divY,spaceZ*iZ);
					}
					if (facePoints && iX < divX) {
						pts[iP++] = position_start+vec3(spaceX*(iX+0.5),0,spaceZ*(iZ+0.5));
						pts[iP++] = position_start+vec3(spaceX*(iX+0.5),spaceY*divY,spaceZ*(iZ+0.5));
					}
				}
			}
			for (int iY = 0; iY < divY; iY++) {
				for (int iZ = 0; iZ < divZ; iZ++) {
					if (iY > 0 && iZ > 0) {
						pts[iP++] = position_start+vec3(0,spaceY*iY,spaceZ*iZ);
						pts[iP++] = position_start+vec3(spaceX*divX,spaceY*iY,spaceZ*iZ);
					}
					if (facePoints) {
						pts[iP++] = position_start+vec3(0,spaceY*(iY+0.5),spaceZ*(iZ+0.5));
						pts[iP++] = position_start+vec3(spaceX*divX,spaceY*(iY+0.5),spaceZ*(iZ+0.5));
					}
				}
			}
			m_vPoints = pts;
			calcDimensions();
		}
	
		inline void transformByMatrix(mat4& mt, double scale = 1)
		{
			movingPoint = mt*movingPoint;
			position_start = mt*m_vPoints[0];
			position_end = mt*m_vPoints[0];
			for (int i=0; i<m_iPointCount; i++)
			{
				m_vPoints[i] = mt * m_vPoints[i]; 	
				if (m_vPoints[i][VX] > position_end[VX])
					position_end[VX] = m_vPoints[i][VX]; 
				if (m_vPoints[i][VX] < position_start[VX])
					position_start[VX] = m_vPoints[i][VX]; 
				if (m_vPoints[i][VY] > position_end[VY])
					position_end[VY] = m_vPoints[i][VY]; 
				if (m_vPoints[i][VY] < position_start[VY])
					position_start[VY] = m_vPoints[i][VY]; 
				if (m_vPoints[i][VZ] > position_end[VZ])
					position_end[VZ] = m_vPoints[i][VZ]; 
				if (m_vPoints[i][VZ] < position_start[VZ])
					position_start[VZ] = m_vPoints[i][VZ]; 
			}	
			mindimension *= scale;
			maxdimension *= scale;
		} 
		
		inline void readFrameFromPoints(double scale = 1)
		{
			position_start = m_vPoints[0];
			position_end = m_vPoints[0];
			movingPoint = 0.5*(position_start + position_end); 
			for (int i=0; i<m_iPointCount; i++)
			{
				if (m_vPoints[i][VX] > position_end[VX])
					position_end[VX] = m_vPoints[i][VX]; 
				if (m_vPoints[i][VX] < position_start[VX])
					position_start[VX] = m_vPoints[i][VX]; 
				if (m_vPoints[i][VY] > position_end[VY])
					position_end[VY] = m_vPoints[i][VY]; 
				if (m_vPoints[i][VY] < position_start[VY])
					position_start[VY] = m_vPoints[i][VY]; 
				if (m_vPoints[i][VZ] > position_end[VZ])
					position_end[VZ] = m_vPoints[i][VZ]; 
				if (m_vPoints[i][VZ] < position_start[VZ])
					position_start[VZ] = m_vPoints[i][VZ]; 
			}	
			mindimension *= scale;
			maxdimension *= scale;
		}
		 
		inline void moveX(double dist)
		{
			for (int i=0; i<m_iPointCount; i++)
			{
				m_vPoints[i][VX] += dist; 	
			}	
			position_start[VX] += dist; 	
			position_end[VX] += dist;
			movingPoint[VX] += dist; 	
		}
		 
		inline void moveZ(double dist)
		{
			for (int i=0; i<m_iPointCount; i++)
			{
				m_vPoints[i][VZ] += dist; 	
			}	
			position_start[VZ] += dist; 	
			position_end[VZ] += dist; 	
			movingPoint[VZ] += dist; 	
		}
		 
		inline void moveY(double dist)
		{
			for (int i=0; i<m_iPointCount; i++)
			{
				m_vPoints[i][VY] += dist; 	
			}	
			position_start[VY] += dist; 	
			position_end[VY] += dist; 	
			movingPoint[VY] += dist; 	
		}

		inline pointframe& operator = (const pointframe& ff)
		{ 
			position_start = ff.position_start;
			position_end = ff.position_end;
			mindimension = ff.mindimension;
			maxdimension = ff.maxdimension;
			movingPoint = ff.movingPoint; 	
			if (m_iPointCount != ff.m_iPointCount)
			{
				if (m_iPointCount > 0)
					delete [] m_vPoints;
				m_iPointCount = ff.m_iPointCount;
				m_vPoints = new vec3[m_iPointCount];
			}
			for (int i = 0; i < m_iPointCount; i++)
			{
				m_vPoints[i] = ff.m_vPoints[i];
			}
			return *this;
		}
};

class line2 {
	protected:
	public:
		vec2 normal;
		vec2 direction;
		double zerodistance;	
		line2() { }
		line2(vec2& p1,vec2& p2)
		{
			if (p1[VX] != p2[VX]) {
				normal[VX] = (p2[VY]-p1[VY])/(p1[VX]-p2[VX]);
				normal[VY] = 1;
				normal.normalize();	
			} else {
				normal[VX] = 1;
				normal[VY] = 0;	
			}
			// direction = normal.rotate(90° cw)
			direction[VX] = normal[VY];
			direction[VY] = -normal[VX];
			zerodistance = normal[VX]*p1[VX] + normal[VY]*p1[VY]; 
		}

		inline double signedDistance(vec2& point)
		{
			return normal*point - zerodistance;	
		}
		
		inline vec2* intersectionPoint(line2& line)
		{
			double ddet = direction*line.normal;
			if (ddet == 0)
				return NULL;
			return new vec2(zerodistance*normal + (line.zerodistance-normal*line.normal*zerodistance)/ddet*direction);
		}
		
		inline line2& operator = (line2& line)
		{
			normal = line.normal;
			direction = line.direction;
			zerodistance = line.zerodistance;
			return *this;
		}
		
		inline void reverse()
		{
			normal = -1*normal;
			direction = -1*direction;
			zerodistance = -1*zerodistance;
		}
};

#endif /*ALGEBRA3EXTENSION_H_*/
