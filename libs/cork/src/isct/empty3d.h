// +-------------------------------------------------------------------------
// | empty3d.h
// | 
// | Author: Gilbert Bernstein
// +-------------------------------------------------------------------------
// | COPYRIGHT:
// |    Copyright Gilbert Bernstein 2013
// |    See the included COPYRIGHT file for further details.
// |    
// |    This file is part of the Cork library.
// |
// |    Cork is free software: you can redistribute it and/or modify
// |    it under the terms of the GNU Lesser General Public License as
// |    published by the Free Software Foundation, either version 3 of
// |    the License, or (at your option) any later version.
// |
// |    Cork is distributed in the hope that it will be useful,
// |    but WITHOUT ANY WARRANTY; without even the implied warranty of
// |    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// |    GNU Lesser General Public License for more details.
// |
// |    You should have received a copy 
// |    of the GNU Lesser General Public License
// |    along with Cork.  If not, see <http://www.gnu.org/licenses/>.
// +-------------------------------------------------------------------------
#pragma once

#include "../math/vec.h"

class Quantization;

class Empty3d {

public:

	Empty3d(Quantization* _quantizer)
		: degeneracy_count(0)
		, exact_count(0)
		, callcount(0)
		, quantizer(_quantizer)
	{
	}

	struct TriIn
	{
		Vec3d p[3];
	};

	struct EdgeIn
	{
		Vec3d p[2];
	};

	struct TriEdgeIn
	{
		TriIn   tri;
		EdgeIn  edge;
	};
	bool isEmpty(const TriEdgeIn &input);
	Vec3d coords(const TriEdgeIn &input) const;
	bool emptyExact(const TriEdgeIn &input);
	Vec3d coordsExact(const TriEdgeIn &input) const;

	// Determines if an edge intersects a triangle
	//    Input:  an edge and a triangle
	//    Return: false if the edge intersects the triangle
	bool emptyApprox( const TriEdgeIn &input );
	Vec3d coordsApprox(const TriEdgeIn &input) const;

	struct TriTriTriIn
	{
		TriIn tri[3];
	};
	bool isEmpty(const TriTriTriIn &input);
	Vec3d coords(const TriTriTriIn &input) const;
	bool emptyExact(const TriTriTriIn &input);
	Vec3d coordsExact(const TriTriTriIn &input) const;

	inline bool hasDegeneracies() const { return degeneracy_count > 0; }

	void resetCounts()
	{
		degeneracy_count = 0;
		exact_count = 0;
		callcount = 0;
	}

	/*
	// exact versions


	bool emptyExact(const Cell3d0 &c0,
	const Complex3d2 &complex,
	const Metric3d2 &metric);

	void cell3d0toPointExact(SmVector3 &point,
	const Cell3d0 &c0,
	const Complex3d2 &complex,
	const Metric3d2 &metric);
	*/

protected:

	int emptyFilter(const TriEdgeIn &input) const;
	bool exactFallback(const TriEdgeIn &input);

	int emptyFilter(const TriTriTriIn &input) const;
	bool exactFallback(const TriTriTriIn &input);

	int degeneracy_count; // count degeneracies encountered
	int exact_count; // count of filter calls failed
	int callcount; // total call count

	Quantization* quantizer;

}; // end class Empty3d
