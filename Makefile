BUILD_DIR ?= build
COVERAGE_DIR ?= build-coverage
BUILD_TYPE ?= Debug
CMAKE ?= cmake

.PHONY: all clean configure test coverage coverage-configure

all: configure
	$(CMAKE) --build $(BUILD_DIR)

configure:
	$(CMAKE) -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=$(BUILD_TYPE)

test: configure
	$(CMAKE) --build $(BUILD_DIR)
	ctest --test-dir $(BUILD_DIR) --output-on-failure

coverage-configure:
	$(CMAKE) -B $(COVERAGE_DIR) -DCMAKE_BUILD_TYPE=Debug -DGIAC_ENABLE_COVERAGE=ON

coverage: coverage-configure
	$(CMAKE) --build $(COVERAGE_DIR)
	$(CMAKE) --build $(COVERAGE_DIR) --target coverage

clean:
	$(CMAKE) --build $(BUILD_DIR) --target clean

distclean:
	rm -rf $(BUILD_DIR) $(COVERAGE_DIR)
