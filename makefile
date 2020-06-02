CC=/usr/bin/g++
### Directories
DEBUG = ./Debug/
Src=./src/
### Patterns 
C_SRC= $(wildcard *.cpp $(Src)/*.cpp )
C_OBJ= $(patsubst %.cpp,%.o,$(C_SRC))
### TMP OBJECT
TMP_OBJ= $(patsubst %.cpp,%.o,$(notdir $(C_SRC)))
F_OBJ= $(patsubst %.o,$(DEBUG)%.o,$(TMP_OBJ))

### Target
TARGET=ephconverter #bamboo #libpppclk.so
.PHONY: all clean
all:$(TARGET)
##sudo apt-get install  libcurl4-openssl-dev
### Compiler Flags
LCFLAGS=   -g -c   -std=c++11  
LLFLAGS=  -lpthread -lm
#######################################################################################
$(TARGET):$(C_OBJ)
	$(CC) -o  $(TARGET) $(TMP_OBJ)  $(LLFLAGS)
	mv $(TMP_OBJ) $(DEBUG) 
%.o:%.cpp
	$(CC) $(LCFLAGS) $<
clean:
	rm $(F_OBJ) $(TARGET) $(TMP_OBJ)
