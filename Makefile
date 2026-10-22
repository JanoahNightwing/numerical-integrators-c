CC     = gcc
CFLAGS = -O2 -Wall -Wextra -std=c99
LDLIBS = -lm

PROGS = euler energy energy_long lorenz timing_O0 timing_O2

all: $(PROGS)

euler energy energy_long lorenz: %: %.c
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

timing_O0: timing.c
	$(CC) -O0 -Wall -Wextra -std=c99 $< -o $@ $(LDLIBS)

timing_O2: timing.c
	$(CC) -O2 -Wall -Wextra -std=c99 $< -o $@ $(LDLIBS)

run: all
	./euler
	./energy
	./energy_long
	./lorenz
	./timing_O0 timing_O0.csv
	./timing_O2 timing_O2.csv

plots: run
	python3 plot.py
	python3 plot_energy.py
	python3 plot_long.py
	python3 plot_lorenz.py
	python3 plot_timing.py

clean:
	rm -f $(PROGS) *.csv

.PHONY: all run plots clean