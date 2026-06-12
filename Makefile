BUILD_DIR ?= build
BUILD_TYPE ?= Debug
CMAKE ?= cmake

.PHONY: all clean configure

all: configure
	$(CMAKE) --build $(BUILD_DIR)

configure:
	$(CMAKE) -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=$(BUILD_TYPE)

clean:
	$(CMAKE) --build $(BUILD_DIR) --target clean

distclean:
	rm -rf $(BUILD_DIR)
