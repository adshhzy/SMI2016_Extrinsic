
INCLU_PATH = ./include
EIGEN_PATH = ./
CPP = g++
CARGS = -std=c++11 -m64 -O3 -Wall -g -I$(INCLU_PATH) -I$(EIGEN_PATH)
TARGET = extrinsic
OBJS = obj/main.o obj/ESDF_Core.o obj/geocurv_display.o obj/geocurv.o obj/geocurv_interface.o obj/geosur_display.o obj/geosur_interface.o obj/geosur.o obj/geovol_display.o obj/geovol_iso.o obj/geovol.o obj/my_mesh.o obj/readers.o

SOURCES = $(OBJS:.o=.cpp,obj=src)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CPP) $(CARGS) $(OBJS) -o $(TARGET)

.cpp.o:
	$(CPP) $(CARGS) $< -o $@

obj/main.o: src/main.cpp
	$(CPP) $(CARGS) -c src/main.cpp -o obj/main.o 

obj/ESDF_Core.o: src/ESDF_Core.cpp include/ESDF_Core.h
	$(CPP) $(CARGS) -c src/ESDF_Core.cpp -o obj/ESDF_Core.o

obj/geocurv_display.o: src/geocurv_display.cpp include/geo_curv.h
	$(CPP) $(CARGS) -c src/geocurv_display.cpp -o obj/geocurv_display.o

obj/geocurv.o: src/geocurv.cpp include/geo_curv.h
	$(CPP) $(CARGS) -c src/geocurv.cpp -o obj/geocurv.o

obj/geocurv_interface.o: src/geocurv_interface.cpp include/geo_curv.h
	$(CPP) $(CARGS) -c src/geocurv_interface.cpp -o obj/geocurv_interface.o

obj/geosur_display.o: src/geosur_display.cpp include/geo_sur.h
	$(CPP) $(CARGS) -c src/geosur_display.cpp -o obj/geosur_display.o

obj/geosur_interface.o: src/geosur_interface.cpp include/geo_sur.h
	$(CPP) $(CARGS) -c src/geosur_interface.cpp -o obj/geosur_interface.o

obj/geosur.o: src/geosur.cpp include/geo_sur.h
	$(CPP) $(CARGS) -c src/geosur.cpp -o obj/geosur.o

obj/geovol_display.o: src/geovol_display.cpp include/geo_vol.h
	$(CPP) $(CARGS) -c src/geovol_display.cpp -o obj/geovol_display.o

obj/geovol_iso.o: src/geovol_iso.cpp include/geo_vol.h
	$(CPP) $(CARGS) -c src/geovol_iso.cpp -o obj/geovol_iso.o

obj/geovol.o: src/geovol.cpp include/geo_vol.h
	$(CPP) $(CARGS) -c src/geovol.cpp -o obj/geovol.o

obj/my_mesh.o: src/geovol.cpp include/my_mesh.h
	$(CPP) $(CARGS) -c src/my_mesh.cpp -o obj/my_mesh.o

obj/readers.o:  src/readers.cpp include/readers.h
	$(CPP) $(CARGS) -c src/readers.cpp -o obj/readers.o

clean: 
	rm -f $(OBJS)
	rm -f $(TARGET)
