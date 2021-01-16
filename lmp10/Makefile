aprox: main.o splines.o points.o aproks_czebyszewa.o gaus/libge.a
	$(CC) -o aprox  main.o splines.o points.o aproks_czebyszewa.o -L gaus -l ge

intrp: main.o splines.o points.o interpolator.o gaus/libge.a
	$(CC) -o intrp  main.o splines.o points.o interpolator.o -L gaus -l ge

prosta: main.o splines.o points.o prosta.o
	$(CC) -o prosta  main.o splines.o points.o prosta.o	

aproks_czebyszewa.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c aproks_czebyszewa.c

aproksymator_na_bazie.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c aproksymator_na_bazie.c

interpolator.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c interpolator.c

.PHONY: clean

clean:
	-rm *.o aprox intrp prosta

tests:
	./aprox -s result/spl_1 -p test/dane_1 -g result/myplot_1 -f 5.3 -t 5.8 -n 30
	./aprox -s result/spl_3 -p test/dane_3 -g result/myplot_3 -f -5 -t 5 -n 100
	./aprox -s result/spl_5 -p test/dane_5 -g result/myplot_5 -f -100 -t 5000 -n 1000
	./aprox -s result/spl_6 -p test/dane_6 -g result/myplot_6 -f -10 -t 10 -n 300	
