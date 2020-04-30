all : TDSPGenerator

CCC = g++
CCFLAGS = -O3 -std=c++11
LIBS= -lm
TARGETDIR=.
CPPFLAGS += \
         -I.

OBJS2 = \
        $(TARGETDIR)/Client.o \
        $(TARGETDIR)/main.o \
        $(TARGETDIR)/Params.o \
        $(TARGETDIR)/ArcProfile.o \
        $(TARGETDIR)/PLFunction.o \
        $(TARGETDIR)/ShortestPath.o \

$(TARGETDIR)/TDSPGenerator: $(OBJS2)
	$(CCC)  $(CCFLAGS) -o $(TARGETDIR)/TDSPGenerator $(OBJS2)

$(TARGETDIR)/Client.o: Client.h Client.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Client.cpp -o $(TARGETDIR)/Client.o

$(TARGETDIR)/main.o: main.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c main.cpp -o $(TARGETDIR)/main.o
	
$(TARGETDIR)/Params.o: Params.h Params.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c Params.cpp -o $(TARGETDIR)/Params.o

$(TARGETDIR)/ArcProfile.o: ArcProfile.h ArcProfile.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c ArcProfile.cpp -o $(TARGETDIR)/ArcProfile.o

$(TARGETDIR)/PLFunction.o: PLFunction.h PLFunction.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c PLFunction.cpp -o $(TARGETDIR)/PLFunction.o

$(TARGETDIR)/ShortestPath.o: Params.h ShortestPath.cpp
	$(CCC) $(CCFLAGS) $(CPPFLAGS) -c ShortestPath.cpp -o $(TARGETDIR)/ShortestPath.o

clean:
	 $(RM) \
	$(TARGETDIR)/TDSPGenerator \
	$(TARGETDIR)/Client.o \
       $(TARGETDIR)/main.o \
     $(TARGETDIR)/Params.o \
     $(TARGETDIR)/ArcProfile.o \
     $(TARGETDIR)/PLFunction.o \
     $(TARGETDIR)/ShortestPath.o



