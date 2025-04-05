THISDIR=$(shell basename `pwd`)
DATE=$(shell date +"%m-%d-%y.%H-%M")


all : payley ;

payley : payley.c ffq.c ffq.h common.c common.h ;
	gcc -O3 -o payley payley.c ffq.c common.c -lm

clean : ;
	rm -f *.o payley

squeaky : clean ;

backup : clean ;
	tar cvf $(THISDIR)-$(DATE).tar --exclude='doc/*' --exclude='old/*' --exclude='paper/*' -C ../ ./$(THISDIR)/
	gzip $(THISDIR)-$(DATE).tar


