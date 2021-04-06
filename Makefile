SHELL := /bin/bash

all:  core
 
core:
	#pip3 install -r requirements.txt
	git submodule update --init
	pwd
	cd WPG; make 
	pwd
	cp WPG/build/lib/*.so WPG/wpg/srw/
	pwd	
	# find directory
	SITEDIR=$(python -m site --user-site)
	# create if it doesn't exist
	mkdir -p "$SITEDIR"
	# create new .pth file with our path
	echo "$(CURDIR)/" > "$SITEDIR/felpy.pth"
	echo "$(CURDIR)/WPG/" > "$SITEDIR/wpg.pth"
	echo "$(CURDIR)/WPG/srw/" > "$SITEDIR/srw.pth"
