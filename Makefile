.ONESHELL:

all:  core
 
core:
	conda env create -f felpy.yml
	conda activate felpy
	pip3 install -r requirements.txt
	git submodule update --init
	pwd
	cd WPG; make 
	pwd
	cd ../; cp WPG/build/lib/*.so WPG/wpg/srw/
	pwd	
	export PYTHONPATH="${PYTHONPATH}:$(CURDIR)"

