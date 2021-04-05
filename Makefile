
all:  core
 
core:
	python setup.py install
	git submodule update --init
	cd WPG; make clean
	cd WPG; make 
	cd WPG/; mv build/lib/*.so wpg/srw/
	export PYTHONPATH="${PYTHONPATH}:$(CURDIR)"

