CPP = $(GCC) -g -std=c++0x
D2P = $(PWD)
APP = MCM

ifeq ($(PARALLEL),1)
CFLAGS+=-DPARALLEL 
endif

ifeq ($(DEBUG),1)
INCLUDES+=-I$(D2P)/common -I$(D2P)/$(APP)
CFLAGS+= -gdwarf-2 -DDEBUG
else
INCLUDES+=-I$(D2P)/common -I$(D2P)/$(APP)
CFLAGS += -$(OPT)
endif

$(D2P)/$(APP):./RecursiveFunctions.cpp ./$(APP)_Aux.cpp ./HelperFunctions.cpp ./main.cpp
	echo $(CPP)
	$(CPP) $(CFLAGS) $(INCLUDES) -o $@ $^

clean:
	rm -f $(APP)
	
