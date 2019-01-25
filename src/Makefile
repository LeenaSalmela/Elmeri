CFLAGS	= -O3 -g -Wall -fopenmp --std=c++11 -lboost_system -lboost_thread
LDFLAGS	= -fopenmp -lm -m64
CC	= g++

CFILES	= lmerindex.cpp elmeri-main.cpp multi-align.cpp main.cpp
HFILES  = lmerindex.hpp

EXTRA = Makefile README.md COPYING

OBJS	= $(CFILES:.cpp=.o)

# Files to compile Elmeri
OBJSE	= lmerindex.o elmeri-main.o multi-align.o

PROG	= elmeri
VERSION = 0.1

# Files to compile the index testing code
OBJSEI	= lmerindex.o main.o

PROGEI	= elmeri-index

all:		$(PROG)

$(PROG):	$(OBJSE)
		$(CC) $(LDFLAGS) -o $@ $(OBJSE)

$(OBJS):	%.o: %.cpp
		$(CC) $(CFLAGS) -c $<

$(PROGEI):	$(OBJSEI)
		$(CC) $(LDFLAGS) -o $@ $(OBJSEI)

clean:
		rm $(OBJS)

dist:
		mkdir $(PROG)-$(VERSION)
		cp $(CFILES) $(PROG)-$(VERSION)/
		cp $(HFILES) $(PROG)-$(VERSION)/
		cp $(EXTRA) $(PROG)-$(VERSION)/
#		sed 's/VERSION/$(VERSION)/g' < README > $(PROG)-$(VERSION)/README
		tar zcvf $(PROG)-$(VERSION).tar.gz $(PROG)-$(VERSION)

dist-clean:
		rm -r $(PROG)-$(VERSION)/
