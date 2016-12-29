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
Vec Vec::operator%(const Vec &b) const {return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);}
//dot product
double Vec::dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }
// Component wise multiply.
Vec Vec::mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
//Calulate norm, return reference? probably more efficient than calling copy constructor?
Vec& Vec::norm(){ return *this = *this * (1.0 / sqrt(x * x + y * y + z * z)); }
double Vec::length() const { return sqrt(x * x + y * y + z * z); }
Vec& Vec::maxnorm() {
  double m = x;
  (m < y) && (m = y);
  (m < z) && (m = z);
  return *this = *this / m;
}
void Vec::show() const {std::cout<<"("<<x<<","<<y<<","<<z<<")";}

Vec tempToColor(double temperature) {
    temperature = clamp(temperature, 1000.0, 40000.0);
    temperature = temperature / 100.0;
    double red, green, blue;

    if( temperature <= 66 ) { 
        red = 255;
        
        green = temperature;
        green = 99.4708025861 * log(green) - 161.1195681661;
        
        if (temperature <= 19)
            blue = 0;
        else {
            blue = temperature - 10.0;
            blue = 138.5177312231 * log(blue) - 305.0447927307;
	}
    } else {
        red = temperature - 60.0;
        red = 329.698727446 * pow(red, -0.1332047592);
        
        green = temperature - 60;
        green = 288.1221695283 * pow(green, -0.0755148492 );

        blue = 255;
    }
    
    return Vec(clamp(red, 0.0, 255.0)/255.0, clamp(green, 0.0, 255.0)/255.0, clamp(blue, 0.0, 255.0)/255.0);
}
