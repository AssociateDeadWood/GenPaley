THISDIR=$(shell basename `pwd`)
DATE=$(shell date +"%m-%d-%y.%H-%M")


all : paley ;

paley : paley.c ffq.c ffq.h common.c common.h ;
	gcc -O3 -o paley paley.c ffq.c common.c -lm

clean : ;
	rm -f *.o paley

squeaky : clean ;

backup : clean ;
	tar cvf $(THISDIR)-$(DATE).tar --exclude='doc/*' --exclude='old/*' --exclude='paper/*' -C ../ ./$(THISDIR)/
	gzip $(THISDIR)-$(DATE).tar


