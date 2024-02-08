CFLAGS= -Wall -m64 -g -w
#CFLAGS= -Wall -m64 -o2 
CXX=g++
#ILOG= /opt/ibm/ILOG/CPLEX_Studio1263
ILOG= /opt/ibm/ILOG/CPLEX_Studio1210
#ILOG= /opt/ibm/ILOG/CPLEX_Studio129
#ILOG= /opt/ibm/ILOG/CPLEX_Studio1271
CPPFLAGS= -DIL_STD -I$(ILOG)/cplex/include -I$(ILOG)/concert/include
CPLEXLIB=-L$(ILOG)/cplex/lib/x86-64_linux/static_pic -lilocplex -lcplex -L$(ILOG)/concert/lib/x86-64_linux/static_pic -lconcert -lm -lpthread

comp-bd1:  
	$(CXX) $(CFLAGS) $(CPPFLAGS) -o bd1  BD.cpp   $(CPLEXLIB) -ldl
comp-bd2n:  
	$(CXX) $(CFLAGS) $(CPPFLAGS) -o bd2n  BD-n.cpp   $(CPLEXLIB) -ldl
comp-bd2n2:  
	$(CXX) $(CFLAGS) $(CPPFLAGS) -o bd2n2  BD-n2.cpp   $(CPLEXLIB) -ldl
comp-bdmwn2:  
	$(CXX) $(CFLAGS) $(CPPFLAGS) -o bdmwn2  BDMW-n2.cpp   $(CPLEXLIB) -ldl

clean:
	rm -f  *.out *.aux *.log *.nav *.snm *.out *.toc 
