FC = gfortran
FFLAGS ?= -O2 -Wall -Wextra -std=f2008

SRC := src/coulomb_two_particles.f90
BIN := bin/coulomb_two_particles
MOD_DIR := build

.PHONY: all clean run

all: $(BIN)

$(BIN): $(SRC)
	@mkdir -p $(dir $@) $(MOD_DIR)
	$(FC) $(FFLAGS) -J $(MOD_DIR) -o $@ $<

run: $(BIN)
	./$(BIN)

clean:
	rm -rf $(BIN) $(MOD_DIR)/*.mod $(MOD_DIR)/*.o

