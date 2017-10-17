CC = g++
CC_OPT = -std=c++11 -O3 -Wall -L./lib -L/lib -fpermissive -fopenmp
CC_INC = -I./src -I./include
LD_LIB = -lSDL2 -lSDL2_image -lSDL2_test -lSOIL -l3ds -lGL -ltbb -ltbbmalloc -ltbbmalloc_proxy  -lomp
OUTPUT_DIR = ./bin
CFLAGS = -Wall
PROG = $(OUTPUT_DIR)/tp
TEST_COMMON_OBJ = 
TEST_COMMON_OBJ_VALGRIND = 
OBJ= $(OUTPUT_DIR)/Loader3ds.o $(OUTPUT_DIR)/Scene.o $(OUTPUT_DIR)/main.o $(OUTPUT_DIR)/sobol.o

$(PROG): $(OBJ)
	$(CC) $(CC_OPT) -o $@ $(OBJ) $(LD_LIB)
	ln -f -s $(PROG) .

clean:
	rm $(OBJ)
	rm $(PROG)

$(OUTPUT_DIR)/Loader3ds.o: src/Geometry/src/Loader3ds.cpp src/Geometry/Loader3ds.h \
 src/Geometry/Material.h src/Geometry/RGBColor.h src/Geometry/Texture.h \
 src/Math/Vectorf.h src/Math/Vector.h src/Geometry/Geometry.h \
 src/Geometry/CastedRay.h src/Geometry/Ray.h src/Geometry/Triangle.h \
 src/Geometry/RayTriangleIntersection.h src/Spy/Spy.h \
 src/Geometry/ComputeVertexNormals.h src/Math/Quaternion.h \
 src/System/aligned_allocator.h src/Math/Constant.h
	$(CC) -c $(CC_OPT) $(CC_INC) -o $@ $< $(LD_LIB)

$(OUTPUT_DIR)/Scene.o: src/Geometry/Scene.cpp
	$(CC) -c $(CC_OPT) $(CC_INC) -o $@ $< $(LD_LIB)

$(OUTPUT_DIR)/main.o: src/main.cpp src/Geometry/Texture.h src/Geometry/RGBColor.h \
 src/Math/Vectorf.h src/Math/Vector.h src/Geometry/Ray.h \
 src/Geometry/Triangle.h src/Geometry/Material.h src/Geometry/CastedRay.h \
 src/Geometry/RayTriangleIntersection.h src/Spy/Spy.h \
 src/Geometry/PointLight.h src/Geometry/Camera.h src/Math/Quaternion.h \
 src/Geometry/Cube.h src/Geometry/Geometry.h \
 src/Geometry/ComputeVertexNormals.h src/System/aligned_allocator.h \
 src/Math/Constant.h src/Geometry/Square.h src/Geometry/Disk.h \
 src/Geometry/Cylinder.h src/Geometry/Cone.h src/Visualizer/Visualizer.h \
 src/Geometry/Scene.h src/Geometry/BoundingBox.h \
 src/Math/RandomDirection.h src/Math/sobol.h src/Geometry/LightSampler.h \
 src/Geometry/Cornel.h src/Geometry/Loader3ds.h
	$(CC) -c $(CC_OPT) $(CC_INC) -o $@ $< $(LD_LIB)

$(OUTPUT_DIR)/sobol.o: src/Math/src/sobol.cpp src/Math/sobol.h
	$(CC) -c $(CC_OPT) $(CC_INC) -o $@ $< $(LD_LIB)

