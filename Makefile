default:
	g++ ./src/main.cpp -O3 -o octree -msse4.2 -ffast-math

debug:
	g++ ./src/main.cpp -o octree -g
clean:
	@rm -f octree 
