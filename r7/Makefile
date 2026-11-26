# FEA-R7: one command rebuilds every number, simulation and figure from the committed DFT results.
#   make            numbers + simulations + figures + LaTeX macros   (a few minutes, needs requirements.txt)
#   make check      the above, then compare key numbers with reference_r7.json
#   make extras     lever ranking (levers_rank.json) and side-flip test (side_flip.json); ~20 min, results are committed
#   make dft        rerun ALL DFT from scratch                          (days; needs environment-dft.yml)
PY ?= python3
.PHONY: all numbers sims figs macros check dft clean extras
all: macros
numbers:
	mkdir -p paper/figures
	cd opt && $(PY) bch_erasure.py 3000 > bch_erasure_output.txt
	cd paper && $(PY) fea_r5.py --figs > /dev/null && $(PY) fea_r6.py > /dev/null && $(PY) fea_r7.py
sims: numbers
	cd opt && L=$$($(PY) -c "import json; print(round(json.load(open('../paper/numbers_r7.json'))['lambda_eV']['minus'], 3))") && \
	          $(PY) fea_chip_sim.py $$L 0.5 1.0 1.5 > chip_sim_output.txt && \
	          $(PY) fea_chip_sim_adaptive.py $$L 1.0 1.5 > chip_sim_adaptive_output.txt
figs: sims
	cd paper && $(PY) fea_r7_figs.py
macros: figs
	cd paper && $(PY) make_macros_r7.py
check: all
	$(PY) check_r7.py
extras:
	cd opt && $(PY) levers_rank.py > levers_rank_output.txt && $(PY) side_flip.py > side_flip_output.txt
dft:
	./dft/run_all_dft.sh
clean:
	rm -f paper/numbers_r5.json paper/numbers_r6.json paper/numbers_r7.json paper/numbers_r7.tex \
	      opt/chip_sim_results.json opt/chip_sim_adaptive_results.json opt/*_output.txt opt/retention_lam*.json opt/bch_erasure.json
	rm -rf paper/figures
