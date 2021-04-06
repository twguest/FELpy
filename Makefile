
all:  core
 
core:
	conda env create -f felpy.yml
	source activate felpy
	pip3 install -r requirements.txt
	git submodule update --init
	cd WPG; make 
	cd WPG/; cp build/lib/*.so wpg/srw/
	export PYTHONPATH="${PYTHONPATH}:$(CURDIR)"

