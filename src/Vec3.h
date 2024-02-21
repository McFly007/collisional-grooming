#ifndef Vec3_h_
#define Vec3_h_
#include <cmath>
// Three dimensional vector struct - can be more streamlined if time permits

struct Vec3;
Vec3 operator*(float r, const Vec3 &v);

struct Vec3
{
	union
	{
		struct
		{
			float x, y, z;
		};
		float D[3];
	};

	// Start with the empty vector
	Vec3() {}
	Vec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z){};

	float &operator[](unsigned int i) { return D[i]; }

	const float &operator[](unsigned int i) const { return D[i]; }

	// Finds the maximum value in the vector
	float maxComponent() const
	{
		float r = x;
		if (y > r)
			r = y;
		if (z > r)
			r = z;
		return r;
	}

	// Finds the minimum value in the vector
	float minComponent() const
	{
		float r = x;
		if (y < r)
			r = y;
		if (z < r)
			r = z;
		return r;
	}

	// Basic operators
	Vec3 operator+(const Vec3 &r) const { return Vec3(x + r.x, y + r.y, z + r.z); }
	Vec3 operator-(const Vec3 &r) const { return Vec3(x - r.x, y - r.y, z - r.z); }
	Vec3 operator*(float r) const { return Vec3(x * r, y * r, z * r); }
	Vec3 operator/(float r) const { return Vec3(x / r, y / r, z / r); }


	Vec3 cmul(const Vec3 &r) const
	{
		return Vec3(x * r.x, y * r.y, z * r.z);
	}

	Vec3 cdiv(const Vec3 &r) const
	{
		return Vec3(x / r.x, y / r.y, z / r.z);
	}

	Vec3 &operator+=(const Vec3 &r)
	{
		x += r.x;
		y += r.y;
		z += r.z;
		return *this;
	}

	Vec3 &operator-=(const Vec3 &r)
	{
		x -= r.x;
		y -= r.y;
		z -= r.z;
		return *this;
	}

	Vec3 &operator*=(float r)
	{
		x *= r;
		y *= r;
		z *= r;
		return *this;
	}

	// Inner/dot product
	float operator*(const Vec3 &r) const
	{
		return x * r.x + y * r.y + z * r.z;
	}

	// Vector Euclidean norm
	float norm() const { return sqrtf(x * x + y * y + z * z); }

	// Vector Eclidean norm squared
	float normSquared() const { return x * x + y * y + z * z; }

	// Cross product
	Vec3 operator^(const Vec3 &r) const
	{
		return Vec3(
			y * r.z - z * r.y,
			z * r.x - x * r.z,
			x * r.y - y * r.x);
	}

	// Normalize vector
	Vec3 normalized() const { return *this / norm(); }
};

inline Vec3 operator*(float r, const Vec3 &v) { return Vec3(v.x * r, v.y * r, v.z * r); }

#endif
