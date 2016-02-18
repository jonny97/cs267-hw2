# CC = g++

# OSTYPE = $(shell uname -s)
# ifeq ($(OSTYPE), Darwin)
#   INCLUDES = -I/Library/Frameworks/SDL.framework/Versions/A/Headers/
#   LIBS = -lSDLmain -lSDL
#   FRAMEWORKS = -framework SDL -framework OpenGL -framework Cocoa
# else
#   INCLUDES = -I/usr/include/SDL
#   LIBS = -lSDLmain -lSDL -lGL -lGLU
#   FRAMEWORKS =

# endif

# visualize: visualize.cpp
# 	$(CC) -o $@ $< $(INCLUDES) $(FRAMEWORKS) $(LIBS)

visualize: visualize.cpp
	g++ visualize.cpp -I/usr/local/include -L/usr/local/lib -lSDL2

clean:
	rm visualize
