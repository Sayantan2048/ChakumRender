#if !defined(mathPrimitives_h__)
#define mathPrimitives_h__
#define PI 3.14159265

struct Vec {
  double x, y, z;
  //Constructor with default values
  Vec(double x_ = 0, double y_ = 0, double z_ = 0);
  //Operator overload: Add two vectors
  Vec operator+(const Vec &b) const;
  //Operator overload: Subtract two vectors
  Vec operator-(const Vec &b) const;
  //Operator overload: Scalar multiplication
  Vec operator*(double b) const;
  //Operator overload: Scalar division
  Vec operator/(double b) const;
  //Operator overload: cross product
  Vec operator%(Vec&b);
  //Dot product
  double dot(const Vec &b) const;
  //Component wise multiplication
  Vec mult(const Vec &b) const;
  // Normalize a vector
  Vec& norm();
  // Find the length of a vector.
  double length() const;
};

struct Ray {
  // Origin, Direction, as in o + td
  Vec o, d;
  //Constructor, weird initialization!!
  Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};

#endif
