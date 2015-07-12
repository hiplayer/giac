SOURCES := $(LOCAL_SRC_FILES)

OUTDIR := $(TOP)/bin
BUILDDIR := $(TOP)/obj/$(MODULE)_intermediates
INCLUDE := $(LOCAL_C_INCLUDES)
SHARED_LIBRARIES := $(LOCAL_SHARED_LIBRARIES)
STATIC_LIBRARIES := $(LOCAL_STATIC_LIBRARIES)
RELEASE_TYPE := $(LOCAL_RELEASE_TYPE)

ifeq ($(RELEASE_TYPE), debug)
CC_FLAGS += \
	-O0 \
	-g -rdynamic

endif

ifeq ($(RELEASE_TYPE), release)
CC_FLAGS += \
	-O3 \
	-funroll-loops \
	-fomit-frame-pointer

endif

CC_FLAGS += \
	-MMD \
	$(INCLUDE)

LD_FLAGS += \
	-L$(OUTDIR) \
	$(SHARED_LIBRARIES)

C_FILE := $(filter %.c, $(SOURCES))
CC_FILE += $(filter %.cc, $(SOURCES))
CPP_FILE := $(filter %.cpp, $(SOURCES))

#$(warning $(CPP_FILE))
OBJECTS :=
OBJECTS += $(addprefix $(BUILDDIR)/,$(C_FILE:.c=.o))
OBJECTS += $(addprefix $(BUILDDIR)/,$(CC_FILE:.cc=.o))
OBJECTS += $(addprefix $(BUILDDIR)/,$(CPP_FILE:.cpp=.o))
OBJECTS += $(addprefix $(OUTDIR)/,$(STATIC_LIBRARIES))

DEPENDS += $(addprefix $(BUILDDIR)/,$(C_FILE:.c=.d))
DEPENDS += $(addprefix $(BUILDDIR)/,$(CC_FILE:.cc=.d))
DEPENDS += $(addprefix $(BUILDDIR)/,$(CPP_FILE:.cpp=.d))

ALL_OBJECT += $(OUTDIR)/$(MODULE)
all:$(ALL_OBJECT)
$(OUTDIR)/$(MODULE): CC_FLAGS := $(CC_FLAGS)
$(OUTDIR)/$(MODULE): LD_FLAGS := $(LD_FLAGS)

ifeq ($(TARGET_TYPE), static)
$(OUTDIR)/$(MODULE): $(OBJECTS)
	mkdir -p $(@D)
	${AR} qc $@ $^
	${RANLIB} $@
endif

ifeq ($(TARGET_TYPE), shared)
$(OUTDIR)/$(MODULE): $(OBJECTS)
	mkdir -p $(@D)
	$(CPP) $^ -o $@ $(LD_FLAGS)
endif

ifeq ($(TARGET_TYPE), executable)
$(OUTDIR)/$(MODULE): $(OBJECTS)
	mkdir -p $(@D)
	$(CPP) $^ -o $@ $(LD_FLAGS)
endif

$(BUILDDIR)/%.o: %.c
	mkdir -p $(@D)
	$(CC) $(CC_FLAGS) -c $< -o $@

$(BUILDDIR)/%.o: %.cpp
	mkdir -p $(@D)
	$(CPP) $(CC_FLAGS) -c $< -o $@

$(BUILDDIR)/%.o: %.cc
	mkdir -p $(@D)
	$(CPP) $(CC_FLAGS) -c $< -o $@

-include $(DEPENDS)

my_CLEAN=$(MODULE)_clean
ALL_CLEAN+=$(my_CLEAN)
clean:$(ALL_CLEAN)
$(my_CLEAN): BUILDDIR := $(BUILDDIR)
$(my_CLEAN): OUTDIR := $(OUTDIR)
$(my_CLEAN): MODULE := $(MODULE)
$(my_CLEAN):
	rm -rf $(BUILDDIR)
	rm -rf $(OUTDIR)/$(MODULE)
