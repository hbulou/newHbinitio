EXEC=Hbinitio.x
all:
	cd src ; make
install:
	mkdir -p bin ; cd src ; make install
clean:
	cd src ; make clean
clean_test:
	cd test/davidson1D ; find ./ -name "*.png" -exec rm {} \; ;  rm  ./davidson1D/*
dav1d:
	cd test/davidson1D ; ../../bin/$(EXEC) inp_davidson1D
experiment:
	cd test/Experiment ; ../../bin/$(EXEC) inp_Experiment
numerov:
	cd test/Numerov ; ../../bin/$(EXEC) inp_Numerov
doc: Hbinitio.pdf
	cd doc ; pdflatex Hbinitio.tex
movie:
	#cd test/davidson1D ; convert -delay 1 output*.png tunnel.gif ; convert tunnel.gif tunnel.avi ; rm *.gif
	cd test/davidson1D ; ffmpeg -r:v 10 -pattern_type glob -i "*.png" -codec:v libx264 -pix_fmt yuv420p -crf 28 -an tunnel.mp4
