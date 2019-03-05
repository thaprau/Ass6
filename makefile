galsim: galsim2.c
	gcc -g -O3 galsim2.c graphics/graphics.h graphics/graphics.c -o galsim -fopenmp -lX11 -lm
run:
	./galsim 2000 input_data/ellipse_N_02000.gal 200 0.00001 0.252 0 1
clean:
	rm galsim