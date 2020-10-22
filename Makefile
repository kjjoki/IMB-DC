# Makefile for the interactive multibundle method IMB-DC for constrained multiobjective DC optimization

FF = gfortran
FFLAGS = -g -fbounds-check -Wall
OPEN =  


all: timbdc

timbdc: constants.o bundle2.o functions.o bundle1.o imbdc.o timbdc.o plqdf1.o 
	$(FF) -o timbdc $(FFLAGS) $(OPEN) constants.o bundle2.o functions.o bundle1.o imbdc.o timbdc.o plqdf1.o 

constants.mod: constants.o constants.f95
	$(FF) -c $(FFLAGS) $(OPEN) constants.f95
	
constants.o: constants.f95
	$(FF) -c $(FFLAGS) $(OPEN) constants.f95
		
bundle2.mod: constants.mod bundle2.o bundle2.f95
	$(FF) -c $(FFLAGS) $(OPEN) bundle2.f95 
	
bundle2.o: constants.mod bundle2.f95
	$(FF) -c $(FFLAGS) $(OPEN) bundle2.f95 	

functions.mod: constants.mod functions.o functions.f95
	$(FF) -c $(FFLAGS) $(OPEN) functions.f95 
	
functions.o: constants.mod functions.f95
	$(FF) -c $(FFLAGS) $(OPEN) functions.f95 
	
bundle1.mod: constants.mod functions.mod bundle1.o bundle1.f95 
	$(FF) -c $(FFLAGS) $(OPEN) bundle1.f95
	
bundle1.o: constants.mod functions.mod bundle1.f95
	$(FF) -c $(FFLAGS) $(OPEN) bundle1.f95 

imbdc.mod: constants.mod bundle2.mod functions.mod bundle1.mod imbdc.o imbdc.f95 
	$(FF) -c $(FFLAGS) $(OPEN) imbdc.f95	 
	
imbdc.o: constants.mod bundle2.mod functions.mod bundle1.mod imbdc.f95
	$(FF) -c $(FFLAGS) $(OPEN) imbdc.f95 
	
timbdc.o: constants.mod bundle2.mod functions.mod bundle1.mod imbdc.mod timbdc.f95
	$(FF) -c $(FFLAGS) $(OPEN) timbdc.f95 

plqdf1.o: plqdf1.f
	$(FF) -c $(FFLAGS) $(OPEN) plqdf1.f 

clean:	
	rm timbdc constants.mod constants.o bundle1.mod bundle1.o bundle2.mod bundle2.o functions.mod functions.o imbdc.mod imbdc.o timbdc.o plqdf1.o  
	echo Clean done