EXEC=Hbinitio.x
all:
	cd src ; make
install:
	cd src ; make install
clean:
	cd src ; make clean
dav1d:
	cd test/davidson1D ; ../../bin/$(EXEC) inp_davidson1D
doc: Hbinitio.pdf
	cd doc ; pdflatex Hbinitio.tex
