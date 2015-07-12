MODULE := $(LOCAL_MODULE).so
CC_FLAGS :=
LD_FLAGS :=
CC_FLAGS += -shared -fPIC
LD_FLAGS += -shared -fPIC

TARGET_TYPE := shared
include $(BUILD_INTERNAL)
