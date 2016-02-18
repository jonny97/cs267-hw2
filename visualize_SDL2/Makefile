CC = g++

OSTYPE = $(shell uname -s)
ifeq ($(OSTYPE), darwin15)
  ### Old configurations for the old SDL for Mac. Change for SDL2 here ###
  INCLUDES =  -I/usr/local/include -L/usr/local/lib 
  LIBS = -lSDLmain -lSDL -lSDL2
  FRAMEWORKS = -framework SDL -framework OpenGL -framework Cocoa
else
  ### Configurations for SDL2 for Ubuntu -- Tested with Ubuntu 14.04 ###
  INCLUDES = -I/usr/local/include -I/usr/local/include/SDL2 -L/usr/local/lib 
  LIBS = -lSDL2main -lSDL2 -framework OpenGL
  FRAMEWORKS = 
endif

visualize: visualize.cpp
	$(CC) -o $@ $< $(INCLUDES) $(FRAMEWORKS) $(LIBS)

clean:
	rm visualize
