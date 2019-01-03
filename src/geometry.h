/***************************************************************************
 *   Copyright (C) 2009-2018 by Veselin Georgiev, Slavomir Kaslev,         *
 *                              Deyan Hadzhiev, Ivan Vankov et al          *
 *   admin@raytracing-bg.net                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/**
 * @File geometry.h
 * @Brief Contains declarations of geometry primitives.
 */
#pragma once

#include "vector.h"
#include <functional>
#include "matrix.h"
#include "scene.h"
#include <vector>

class Geometry;
struct IntersectionInfo {
	double dist;
	Vector ip;
	Vector norm, dNdx, dNdy;
	double u, v;
	Geometry* geom;
};

class Intersectable {
public:
	virtual ~Intersectable() {}
	
	/// returns true if the ray intersects the geometry
	virtual bool intersect(const Ray& ray, IntersectionInfo& info) = 0;
};

class Geometry: public Intersectable, public SceneElement {
public:
	ElementType getElementType() const { return ELEM_GEOMETRY; }
};

class Plane: public Geometry {
public:
	double limit;
	double height;
	Plane(double limit = 128, double height = 0): limit(limit), height(height) {}

	void fillProperties(ParsedBlock& pb)
	{
		pb.getDoubleProp("y", &height);
		pb.getDoubleProp("limit", &limit);
	}
	
	bool intersect(const Ray& ray, IntersectionInfo& info) override;
};

class Sphere: public Geometry {
public:
	Vector O = Vector(0, 0, 0);
	double R = 1;
	
	Sphere() {}
	Sphere(const Vector& position, double radius)
	{
		O = position;
		R = radius;
	}
	void fillProperties(ParsedBlock& pb)
	{
		pb.getVectorProp("O", &O);
		pb.getDoubleProp("R", &R);
	}
	
	bool intersect(const Ray& ray, IntersectionInfo& info) override;
};

class Cone : public Geometry {
	void copyFrom(const Cone& other);
public:
	Vector O = Vector(0, 0, 0);
	// f(x)=x*k+n
	double k = 1, n = 0;
	double lowerBound = 0, higherBound = 100;

	Cone() {}

	Cone(const Vector& position, double k, double n, double lowerBound = 0, double higherBound = 100);
	Cone(const Cone& other);
	Cone& operator=(const Cone& other);

	void fillProperties(ParsedBlock& pb)
	{
		pb.getVectorProp("O", &O);
		pb.getDoubleProp("k", &k);
		pb.getDoubleProp("n", &n);
		pb.getDoubleProp("lowerBound", &lowerBound);
		pb.getDoubleProp("higherBound", &higherBound);
	}

	bool intersect(const Ray& ray, IntersectionInfo& info) override;
};

class SOR : public Geometry {
private:
	void conesFromCurve(const Vector& origin, const std::vector<std::pair<double, double>>& curve);

	void loadCones(const char* fileName);

	std::vector<Cone> cones;
public:
	SOR() {}

	SOR(const Vector& origin, const std::vector<std::pair<double, double>>& curve);

	void fillProperties(ParsedBlock& pb)
	{
		char fileName[_MAX_PATH];
		pb.requiredProp("file");
		pb.getFilenameProp("file", fileName);
		loadCones(fileName);
	}

	bool intersect(const Ray& ray, IntersectionInfo& info) override;
};

class Cube: public Geometry {
	void intersectCubeSide(const Ray& ray, double start, double dir, double target,
							const Vector& normal, IntersectionInfo& info,
							std::function<void(const Vector&)> uv_mapping);

public:
	Vector O = Vector(0, 0, 0);
	double halfSide = 1;
	
	Cube() {}
	Cube(const Vector& position, double halfSide)
	{
		O = position;
		this->halfSide = halfSide;
	}
	void fillProperties(ParsedBlock& pb)
	{
		pb.getVectorProp("O", &O);
		pb.getDoubleProp("halfSide", &halfSide);
	}
	
	bool intersect(const Ray& ray, IntersectionInfo& info) override;
	
};

class CsgOp: public Geometry {
public:
	Geometry* left;
	Geometry* right;
	
	virtual bool boolOp(bool inLeft, bool inRight) = 0;

	void fillProperties(ParsedBlock& pb)
	{
		pb.requiredProp("left");
		pb.requiredProp("right");
		pb.getGeometryProp("left", &left);
		pb.getGeometryProp("right", &right);
	}
	
	bool intersect(const Ray& ray, IntersectionInfo& info) override;
};

class CsgPlus: public CsgOp {
public:
	bool boolOp(bool inLeft, bool inRight) override
	{
		return inLeft || inRight;
	}
};

class CsgIntersect: public CsgOp {
public:
	bool boolOp(bool inLeft, bool inRight) override
	{
		return inLeft && inRight;
	}
};

class CsgMinus: public CsgOp {
public:
	bool boolOp(bool inLeft, bool inRight) override
	{
		return inLeft && !inRight;
	}
};

struct Shader;
struct BumpTexture;
struct Node: public Intersectable, public SceneElement {
	Geometry* geometry = nullptr;
	Shader* shader = nullptr;
	Transform T;
	Texture* bump = nullptr;

	// from Intersectable:
	bool intersect(const Ray& ray, IntersectionInfo& info) override;

	// from SceneElement:
	ElementType getElementType() const { return ELEM_NODE; }
	void fillProperties(ParsedBlock& pb)
	{
		pb.getGeometryProp("geometry", &geometry);
		pb.getShaderProp("shader", &shader);
		pb.getTransformProp(T);
		pb.getTextureProp("bump", &bump);
	}
};

