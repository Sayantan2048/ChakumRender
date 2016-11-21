#if !defined (random_h__)
#define random_h__

typedef unsigned int mt_uint32;

class Random {
  const static mt_uint32 N = 624; // length of state vector
  const static mt_uint32 M = 397; // a period parameter
  const static mt_uint32 K = 0x9908B0DFU; // a magic constant
  mt_uint32   state[N+1];     // state vector + 1 extra to not violate ANSI C
  mt_uint32   *next;          // next random value is computed from here
  int      left; 	      // can *next++ this many times before reloading
  mt_uint32 reloadMT(void);
  // Seed the generator with a value
  void seedMT(mt_uint32 seed);

public:
  Random(mt_uint32 seed) {
    left = -1;
    seedMT(seed);
  };
  // Get an unsigned long
  mt_uint32 randomMT(void);
  // Get a double FP random number between a and b.
  // The double FP number is as precise as the 32bit unsigned int.
  double randomMTD(double a, double b);

};

#endif