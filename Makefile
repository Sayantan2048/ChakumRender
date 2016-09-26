CC = g++
CFLAGS = -Wall -g
ROOT=$(realpath $(dir $(lastword $(MAKEFILE_LIST))))

INCLUDES = $(ROOT)/headers
SRC = $(ROOT)/src
OBJ = $(ROOT)/src

#make will look for files in these path apart from current location of makefile
VPATH = $(SRC) $(INCLUDES) $(OBJ)

main: domainSampler.o mathPrimitives.o random.o main.o materialTypes.o geometryPrimitives.o objects.o shader.o lightSources.o
	$(CC) $(OBJ)/domainSampler.o $(OBJ)/mathPrimitives.o $(OBJ)/random.o $(OBJ)/materialTypes.o $(OBJ)/geometryPrimitives.o $(OBJ)/objects.o $(OBJ)/shader.o $(OBJ)/lightSources.o $(OBJ)/main.o

main.o: main.cpp
	$(CC) $(CFLAGS) -I$(INCLUDES) -O -c $(SRC)/main.cpp -o $(OBJ)/main.o

random.o: random.h random.cpp
	$(CC) $(CFLAGS) -I$(INCLUDES) -O -c $(SRC)/random.cpp -o $(OBJ)/random.o

mathPrimitives.o: mathPrimitives.h mathPrimitives.cpp
	$(CC) $(CFLAGS) -I$(INCLUDES) -O -c $(SRC)/mathPrimitives.cpp -o $(OBJ)/mathPrimitives.o

domainSampler.o: domainSampler.cpp domainSampler.h
	$(CC) $(CFLAGS) -I$(INCLUDES) -O -c $(SRC)/domainSampler.cpp -o $(OBJ)/domainSampler.o

materialTypes.o: materialTypes.cpp materialTypes.h
	$(CC) $(CFLAGS) -I$(INCLUDES) -O -c $(SRC)/materialTypes.cpp -o $(OBJ)/materialTypes.o

geometryPrimitives.o: geometryPrimitives.cpp geometryPrimitives.h
	$(CC) $(CFLAGS) -I$(INCLUDES) -O -c $(SRC)/geometryPrimitives.cpp -o $(OBJ)/geometryPrimitives.o

objects.o: objects.h objects.cpp
	$(CC) $(CFLAGS) -I$(INCLUDES) -O -c $(SRC)/objects.cpp -o $(OBJ)/objects.o

shader.o: shader.h shader.cpp
	$(CC) $(CFLAGS) -I$(INCLUDES) -O -c $(SRC)/shader.cpp -o $(OBJ)/shader.o

lightSources.o: lightSources.h lightSources.cpp
	$(CC) $(CFLAGS) -I$(INCLUDES) -O -c $(SRC)/lightSources.cpp -o $(OBJ)/lightSources.o



.PHONY: clean

clean:
	cd $(OBJ) && rm *.o


