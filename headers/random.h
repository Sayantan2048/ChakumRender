#if !defined (random_h__) 
#define random_h__ 

typedef unsigned int mt_uint32;

// Get an unsigned long
extern mt_uint32 randomMT(void);
// Get a double FP random number between a and b.
// The double FP number is as precise as the 32bit unsigned int.
extern double randomMTD(double a, double b);
// Seed the generator with a value
extern void seedMT(mt_uint32 seed); 


#endif