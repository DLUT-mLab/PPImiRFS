appname := PPImiRFS

CXX := g++
CXXFLAGS := -std=c++11

# srcfiles := $(shell find . -name "*.cpp")
srcfiles := $(wildcard *.cpp)
objects  := $(patsubst %.cpp, %.o, $(srcfiles))

all: $(appname)

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) -o $(appname) $(objects) $(LDLIBS)

clean:
	rm $(objects)
	rm $(appname)
	rm result/*

dist-clean: clean
	rm -f *~ .depend

# depend: .depend

# .depend: $(srcfiles)
    # rm -f ./.depend
    # $(CXX) $(CXXFLAGS) -MM $^>>./.depend;

# include .depend
run:
	./PPImiRFS -i ../inputfile/test_miRNA_pairs.txt -p ../inputfile/wppin.txt -o ./result/output.txt
