# Compiler
CXX         := g++
#CXX         := CC

# Linker
LD          := mpic++
#LD          := CC

# Library paths
LOCAL_PTH   := /usr/local
EIGEN_PTH   := /system/software/linux-x86_64/lib/eigen/3.2.0
#EIGEN_PTH   :=/home/e280/e280/mtortora
MPICC_PTH   := /system/software/arcus-b/lib/mpi/mvapich2/2.1.0/gcc-4.9.2

# Append mpic++ directories to the PATH environment variable
export PATH := $(LOCAL_PTH)/bin:$(PATH)
export PATH := $(MPICC_PTH)/bin:$(PATH)

# Target binary program
TARGET      := chiraldft

# Version control
VSN_MAJ     := 3
VSN_MIN     := 2

# Directories - sources, includes, objects, binaries, data and resources
SRCDIR      := src
INCDIR      := include
LIBDIR      := lib
BUILDDIR    := build
TARGETDIR   := bin
DATDIR      := data
RESDIR      := resources

# File extensions
SRCEXT      := cpp
HDREXT      := hpp
DATEXT	    := out
PRCEXT      := res
DEPEXT      := d
OBJEXT      := o

# Flags, libraries and includes
CXXFLAGS    := -Wno-logical-op-parentheses -O3
CFLAGS      := -std=c++0x -Werror -Wshadow -Wall -Wextra -msse4 -fno-common -fomit-frame-pointer -O3
FPATHS      := -D__DPATH__=$(CURDIR)/$(DATDIR)
FEXTRA      := -D__VERSION_MAJOR__=$(VSN_MAJ) -D__VERSION_MINOR__=$(VSN_MIN) -DNDEBUG
INC         := -I$(INCDIR) -I$(LIBDIR) -isystem $(LOCAL_PTH)/include -isystem $(EIGEN_PTH)/include -isystem $(MPICC_PTH)/include
LIB         := -L$(LIBDIR) -lRAPID

# Object files lists
SOURCES     := $(shell find $(SRCDIR) -type f -iname "*.$(SRCEXT)")
OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))
RLOBJ       := $(LIBDIR)/rapid2/RAPID.o $(LIBDIR)/rapid2/build.o $(LIBDIR)/rapid2/collide.o $(LIBDIR)/rapid2/overlap.o
RLCLEAN     := $(LIBDIR)/libRAPID.a $(RLOBJ) $(LIBDIR)/rapid2/*~

# Default make
all: resources $(TARGET)
	@echo -e "\033[1;31mCompiled $(TARGET) v$(VSN_MAJ).$(VSN_MIN) with MPI support\033[0m"
	@echo -e "\033[1;37m*****************\033[0m"

# Local run
run: cleandat all
	@mpirun -np 8 $(TARGETDIR)/$(TARGET)

# Background run
nrun: cleandat all
	@nohup mpirun -np 8 $(TARGETDIR)/$(TARGET) > $(TARGETDIR)/log.$(DATEXT)&

# Distributed run
submit: cleandat all
	@qsub $(TARGETDIR)/arcus_b.sh
#	@qsub $(TARGETDIR)/archer.sh

# Copy scripts from resources directory to target directory, substituting file paths
resources: directories
	@find $(RESDIR) -maxdepth 1 -type f -iname "*.sh" -exec basename {} \; | xargs -I{} sh -c 'sed -e "s|TARGETPATH|$(CURDIR)/$(TARGETDIR)|;s|DATPATH|$(CURDIR)/$(DATDIR)|;s|TARGET|$(TARGET)|" < "$(RESDIR)/{}" > "$(TARGETDIR)/{}"'

# Make directories
directories: FORCE
	@echo -e "\033[1;37m*****************\033[0m"
	@mkdir -p $(DATDIR)
	@mkdir -p $(BUILDDIR)
	@mkdir -p $(TARGETDIR)

# Clean objects
clean:
	@echo -e "\033[1;31mCleaning object files...\033[0m"
	@$(RM) -rf $(BUILDDIR)

# Clean data
cleandat:
	@echo -e "\033[1;31mCleaning data files...\033[0m"
	@$(RM) -f *.$(DATEXT) $(DATDIR)/*.$(DATEXT) $(DATDIR)/*.$(PRCEXT) $(TARGETDIR)/*.$(DATEXT)

# Full clean, data, objects and binaries
cleaner: cleandat clean cleanlib
	@echo -e "\033[1;31mCleaning binaries and library files...\033[0m"
	@$(RM) -rf $(TARGETDIR)

# Pull in dependency info for *existing* .o files
-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

# Link
$(TARGET): $(OBJECTS)
	@echo -e "\033[1;31mLinking object files...\033[0m"
	@$(LD) -o $(TARGETDIR)/$(TARGET) $^ $(LIB)

# Compile
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@echo -e "\033[1;34mBuilding: \033[0m" $<
	@mkdir -p $(dir $@)
	@$(CXX) $(CFLAGS) $(FPATHS) $(FEXTRA) $(INC) -c -o $@ $<
	@$(CXX) $(CFLAGS) $(FPATHS) $(FEXTRA) $(INC) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@$(RM) -f $(BUILDDIR)/$*.$(DEPEXT).tmp

# Build RAPID library
librapid: cleanlib $(RLOBJ)
	@ar ruv $(LIBDIR)/libRAPID.a $(RLOBJ)

# Clean library
cleanlib:
	@$(RM) -f $(RLCLEAN)

# Force dependency - include for rules to be run at every call
FORCE:

# Non-file targets
.PHONY: all remake run nrun submit resources directories clean cleandat cleaner librapid cleanlib
