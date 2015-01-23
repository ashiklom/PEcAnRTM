INCL=-I/share/pkg/r/3.1.1/install/lib64/R/include
LIBS=-L/share/pkg/r/3.1.1/install/lib64/R/lib -lR

project_bayes.so: sample_alpha_n.o 
	gcc -shared -Xlinker $(LIBS) -o prospect_bayes.so sample_alpha_n.o 

sample_alpha_n.o: sample_alpha_n.c
	gcc -fpic -O2 -std=gnu99 $(INCL) -c sample_alpha_n.c -o sample_alpha_n.o

.PHONY: clean

clean: 
	rm -f prospect_bayes.so sample_alpha_n.o 

