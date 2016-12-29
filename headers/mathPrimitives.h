#if !defined(mathPrimitives_h__)
#define mathPrimitives_h__
#define SQRT2 1.41421356
#define PI 3.14159265359
#define INF 1e20

#define maX(a, b) ((a) > (b) ? (a) : (b))
#define miN(a, b) ((a) < (b) ? (a) : (b))
#define isNotClose(a, b, eps) (fabs(a - b) > eps)
#define to_greyScale(L) ((L.x * 0.3 + L.y * 0.59 + L.z * 0.11))
#define clamp(x, min, max) (( (x) < (min) ) ? (min) : ( ((x) > (max)) ? (max) : (x) ))

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
  Vec operator%(const Vec &b) const;
  //Dot product
  double dot(const Vec &b) const;
  //Component wise multiplication
  Vec mult(const Vec &b) const;
  // Normalize a vector
  Vec& norm();
  // Find the length of a vector.
  double length() const;
  Vec& maxnorm();
  void show() const;
};

struct Ray {
  // Origin, Direction, as in o + td
  Vec o, d;
  //Constructor, weird initialization!!
  Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};

struct mat3 {
    double a[9];
    double det;
    mat3 mul(const mat3 &m) const {
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

      retm.det = det * m.det;
      return retm;
    }

    mat3 mul(const double s) const {
      mat3 retm;
      retm.a[0] = a[0] * s;
      retm.a[1] = a[1] * s;
      retm.a[2] = a[2] * s;
      retm.a[3] = a[3] * s;
      retm.a[4] = a[4] * s;
      retm.a[5] = a[5] * s;
      retm.a[6] = a[6] * s;
      retm.a[7] = a[7] * s;
      retm.a[8] = a[8] * s;

      retm.det = det * s * s * s;

      return retm;
    }

    Vec mul(const Vec &v) const {
      return Vec(a[0] * v.x + a[1] * v.y + a[2] * v.z,
                 a[3] * v.x + a[4] * v.y + a[5] * v.z,
                 a[6] * v.x + a[7] * v.y + a[8] * v.z);
    }

    mat3 inv() const {
      mat3 retm;

      retm.det = 1.0 / det;
      retm.a[0] = retm.det * (a[4] * a[8] - a[5] * a[7]);
      retm.a[1] = retm.det * (a[2] * a[7] - a[1] * a[8]);
      retm.a[2] = retm.det * (a[1] * a[5] - a[2] * a[4]);

      retm.a[3] = retm.det * (a[5] * a[6] - a[3] * a[8]);
      retm.a[4] = retm.det * (a[0] * a[8] - a[2] * a[6]);
      retm.a[5] = retm.det * (a[2] * a[3] - a[0] * a[5]);

      retm.a[6] = retm.det * (a[3] * a[7] - a[4] * a[6]);
      retm.a[7] = retm.det * (a[1] * a[6] - a[0] * a[7]);
      retm.a[8] = retm.det * (a[0] * a[4] - a[1] * a[3]);

      return retm;
    }

    mat3 transpose() const {
      mat3 retm;

      retm.a[0] = a[0];
      retm.a[1] = a[3];
      retm.a[2] = a[6];

      retm.a[3] = a[1];
      retm.a[4] = a[4];
      retm.a[5] = a[7];

      retm.a[6] = a[2];
      retm.a[7] = a[5];
      retm.a[8] = a[8];

      retm.det = det;

      return retm;
    }

    mat3(){
      a[1] = a[2] = a[3] = a[5] = a[6] = a[7] = 0;
      a[0] = a[4] = a[8] = 1;
      det = 1;;
    }

    mat3(const Vec &col0, const Vec &col1, const Vec &col2) {
      a[0] = col0.x;
      a[3] = col0.y;
      a[6] = col0.z;

      a[1] = col1.x;
      a[4] = col1.y;
      a[7] = col1.z;

      a[2] = col2.x;
      a[5] = col2.y;
      a[8] = col2.z;

      calcDet();
    }

    mat3(double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8) {
      a[0] = a0;
      a[3] = a1;
      a[6] = a2;

      a[1] = a3;
      a[4] = a4;
      a[7] = a5;

      a[2] = a6;
      a[5] = a7;
      a[8] = a8;

      calcDet();
    }

    void calcDet() {
      det = a[0] * (a[4] * a[8] - a[5] * a[7]) -
	     a[1] * (a[3] * a[8] - a[5] * a[6]) +
	     a[2] * (a[3] * a[7] - a[4] * a[6]);
    }
};

// Convert temperature in range of 1000K to 40000K to RGB.
extern Vec tempToColor(double temperature);
#endif
