dir_maker=mkdir -p
output_directory=Outputs

default:
	g++ ./src/FreeParameters.cpp ./src/Elementary_Functions.cpp ./src/Outputs.cpp ./src/Records.cpp ./src/main.cpp -O3 -o coll_groom -msse4.2 -ffast-math -Wall
	$(dir_maker) $(output_directory)

debug:
	g++ ./src/Outputs.cpp ./src/main.cpp -o coll_groom -g
clean:
	@rm -f coll_groom 
