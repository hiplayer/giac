all:
	$(MAKE) -C libtommath
	$(MAKE) -C giac
clean:
	$(MAKE) -C libtommath           clean
	$(MAKE) -C giac                 clean
