# USAGE: make [install|clean]
# 	 make RELEASE=yes



# -------- CAN BE MODIFIED --------

#PACKAGES=package1.tgz package2.tgz
PACKAGES=probmetab.tgz

# ------------------------

all: $(PACKAGES)
ifeq ($(RELEASE),yes)
	@echo "--Release version `date "+%Y%m%d"`--"
endif

# -------- CAN BE MODIFIED --------

probmetab.tgz: probmetab.r lib.r ProbMetab.xml static README.txt export.class.table-color-graph.R
	$(SEDVERSION) 
	tar --exclude=".svn" -zchf $@ $^

# ------------------------

RELEASE=no
SEDVERSION=
ifeq ($(RELEASE),yes)
	SEDVERSION=sed -i "s/version=\"[0-9]\+\"/version=\"`date "+%Y%m%d"`\"/" *.xml; sed -i "s/version=\"[0-9]\+\"/version=\"`date "+%Y%m%d"`\"/" *.r ; sed -i "s/version=\"[0-9]\+\"/version=\"`date "+%Y%m%d"`\"/" *.txt 
endif

# ------------------------

install: $(PACKAGES)
	cp $(PACKAGES) ~

clean:
	rm $(PACKAGES)
