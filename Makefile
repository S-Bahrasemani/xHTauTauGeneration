.PHONY: dump

default: clean

clean-prod:
	rm -rf ./prod/*

clean-pyc:                                                                      
	find cluster gen -name "*.pyc" -exec rm {} \;

clean-jo:
	rm joboptions/*py

clean: clean-pyc

tarball:
	tar cvzf joboptions.tar.gz joboptions/*.py

