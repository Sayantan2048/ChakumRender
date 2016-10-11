#if !defined(mathPrimitives_h__)
#define mathPrimitives_h__
#define PI 3.14159265359
#define INF 1e20

#define maX(a, b) ((a) > (b) ? (a) : (b))
#define miN(a, b) ((a) < (b) ? (a) : (b))

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

#define dV Vec(1., 1., 1.)
struct mat3 {
    double a[9];
    mat3 mul(mat3 m) {
      mat3 retm;
      retm.a[0] = a[0] * m.a[0] + a[1] * m.a[3] + a[2] * m.a[6];
      retm.a[1] = a[0] * m.a[1] + a[1] * m.a[4] + a[2] * m.a[7];
      retm.a[2] = a[0] * m.a[2] + a[1] * m.a[5] + a[2] * m.a[8];
      retm.a[3] = a[3] * m.a[0] + a[4] * m.a[3] + a[5] * m.a[6];
      retm.a[4] = a[3] * m.a[1] + a[4] * m.a[4] + a[5] * m.a[7];
      retm.a[5] = a[3] * m.a[2] + a[4] * m.a[5] + a[5] * m.a[8];
      retm.a[6] = a[6] * m.a[0] + a[7] * m.a[3] + a[8] * m.a[6];
      retm.a[7] = a[6] * m.a[1] + a[7] * m.a[4] + a[8] * m.a[7];
      retm.a[8] = a[6] * m.a[2] + a[7] * m.a[5] + a[8] * m.a[8];

      return retm;
    }

    Vec mul(Vec v) {
      return Vec(a[0] * v.x + a[1] * v.y + a[2] * v.z,
                 a[3] * v.x + a[4] * v.y + a[5] * v.z,
                 a[6] * v.x + a[7] * v.y + a[8] * v.z);
    }

    mat3(Vec col0 = dV, Vec col1 = dV, Vec col2 = dV) {
      a[0] = col0.x;
      a[3] = col0.y;
      a[6] = col0.z;

      a[1] = col1.x;
      a[4] = col1.y;
      a[7] = col1.z;

      a[2] = col2.x;
      a[5] = col2.y;
      a[8] = col2.z;
    }
};
#undef dV

#endif
