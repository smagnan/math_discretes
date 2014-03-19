#VARIABLES:
#----------
COMP=gcc
EDL=gcc
EXE=voyCom_a_distribuer
#Flags:
CPPFLAGS=-g
EDLFLAGS=
#Paths:
INCPATH=#-I /home/samuel/Documents/TpMultiTache
LIBPATH=#-L /home/samuel/Documents/TpMultiTache
#Affichage:
ECHO=@echo #pour le mode silencieux
#Suppression:
RM=@rm
RMFLAGS=-f
#
INT=voyCom_a_distribuer.h#Mere.h Clavier.h Entree.h Sortie.h
REAL=$(INT:.h=.cpp)
OBJ=$(INT:.h=.o)
LIBS=-lm -lrt#-ltp -lncurses -ltcl
#Clean:
CLEAN=clean
.PHONY:$(CLEAN)


#Règles:
#-------

$(EXE): $(OBJ)
	$(ECHO) + EDL de $(EXE)
	$(EDL) -o $(EXE) $(OBJ) $(EDLFLAGS) $(LIBPATH) $(LIBS)
voyCom_a_distribuer.o : voyCom_a_distribuer.c

#Patterns:
#---------

%.o : %.c
	$(ECHO) + Compliation de $<
	$(COMP) -c $(CPPFLAGS) $(INCPATH) $<


#Nettoyage:
#----------

$(CLEAN):
	$(ECHO) --- NETTOYAGE ---
	$(RM) $(RMFLAGS) $(OBJ) $(EXE) #Core 
