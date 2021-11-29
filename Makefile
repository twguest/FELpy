SHELL := /bin/bash

all:  core

core:
	#pip3 install -r requirements.txt
	git submodule update --init
	cd WPG; make
	# find directory
	SITEDIR=$(python -m site --user-site)
	# create if it doesn't exist
	mkdir -p "$SITEDIR"
	CURDIR = $(pwd)
	# create new .pth file with our path
	echo "$(CURDIR)/" > "$SITEDIR/felpy.pth"
	echo "$(CURDIR)/WPG/" > "$SITEDIR/wpg.pth"
	echo "$(CURDIR)/WPG/wpg/" > "$SITEDIR/srw.pth"
