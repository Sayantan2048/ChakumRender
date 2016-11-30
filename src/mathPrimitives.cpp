#include <cmath>
#include <iostream>
#include "mathPrimitives.h"

//Constructor with default values
Vec::Vec(double x_, double y_, double z_){ x=x_; y=y_; z=z_; }
//Operator overload: Add two vectors
Vec Vec::operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
//Operator overload: subtract two vectors
Vec Vec::operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
//Operator overload: scalar multiplication
Vec Vec::operator*(double b) const { return Vec(x * b, y * b, z * b); }
//Operator overload: scalar division
Vec Vec::operator/(double b) const { return Vec(x / b, y / b, z / b); }
//Operator overload: vector cross product. Why no const argument?? Why no const after funct declaration??
Vec Vec::operator%(Vec&b){return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);}
//dot product
double Vec::dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }
// Component wise multiply.
Vec Vec::mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
//Calulate norm, return reference? probably more efficient than calling copy constructor?
Vec& Vec::norm(){ return *this = *this * (1.0 / sqrt(x * x + y * y + z * z)); }
double Vec::length() const { return sqrt(x * x + y * y + z * z); }
void Vec::show() const {std::cout<<"("<<x<<","<<y<<","<<z<<")";}
