CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall

SIMDIR   := simulations
# Version directories. Comments sit on their own lines: Make keeps whitespace
# that precedes a '#', and that whitespace would split a prerequisite list.
V1DIR    := $(SIMDIR)/v1
V2DIR    := $(SIMDIR)/v2
V3DIR    := $(SIMDIR)/v3
V1_SRC   := $(V1DIR)/FEA_sim_v1.cpp
V2_SRC   := $(V2DIR)/FEA_sim_v2.cpp
FZC_SRC  := $(V3DIR)/FEA_fzc_v3.cpp
SLING_SRC:= $(V3DIR)/FEA_slingshot_v3.cpp
REL_SRC  := $(V3DIR)/FEA_fzc_reliability_v3.cpp
RECOG_SRC:= $(V3DIR)/FEA_fzc_recognition_v3.cpp
OPEN_SRC := $(V3DIR)/FEA_fzc_opensystem_v3.cpp
SEL_SRC  := $(V3DIR)/FEA_fzc_selector_v3.cpp
SCR_SRC  := $(V3DIR)/FEA_fzc_screened_v3.cpp
PLAN_SRC := $(V3DIR)/FEA_fzc_floorplan_v3.cpp
E2E_SRC  := $(V3DIR)/FEA_fzc_e2e_v3.cpp
BUD_SRC  := $(V3DIR)/FEA_budget_v3.cpp
DIE_SRC  := $(V3DIR)/FEA_floorplan_v3.cpp
GAM_SRC  := $(V3DIR)/FEA_gamma_v3.cpp
RET_SRC  := $(V3DIR)/FEA_retention_v3.cpp
MF_SRC   := $(V3DIR)/FEA_multifire_v3.cpp
SC_SRC   := $(V3DIR)/FEA_secded_v3.cpp
CT_SRC   := $(V3DIR)/FEA_crosstalk_v3.cpp
RST_SRC  := $(V3DIR)/FEA_restoration_v3.cpp
REF_SRC  := $(V3DIR)/FEA_refresh_v3.cpp
CLK_SRC  := $(V3DIR)/FEA_clock_v3.cpp
BW_SRC   := $(V3DIR)/FEA_bandwidth_v3.cpp
CMP_SRC  := $(V3DIR)/FEA_compare_v3.cpp
FAB_SRC  := $(V3DIR)/FEA_fabrication_v3.cpp
LAY_SRC  := $(V3DIR)/FEA_layout_v3.cpp
PRG_SRC  := $(V3DIR)/FEA_fzc_program_v3.cpp
REC_SRC  := $(V3DIR)/FEA_recovery_v3.cpp
RESC_SRC := $(V3DIR)/FEA_fzc_rescue_v3.cpp
THERM_SRC := $(V3DIR)/FEA_thermal_v3.cpp
V1_BIN   := FEA_sim_v1
V2_BIN   := FEA_sim_v2
FZC_BIN  := FEA_fzc_v3
SLING_BIN:= FEA_slingshot_v3
REL_BIN  := FEA_fzc_reliability_v3
RECOG_BIN:= FEA_fzc_recognition_v3
OPEN_BIN := FEA_fzc_opensystem_v3
SEL_BIN  := FEA_fzc_selector_v3
SCR_BIN  := FEA_fzc_screened_v3
PLAN_BIN := FEA_fzc_floorplan_v3
E2E_BIN  := FEA_fzc_e2e_v3
BUD_BIN  := FEA_budget_v3
DIE_BIN  := FEA_floorplan_v3
GAM_BIN  := FEA_gamma_v3
RET_BIN  := FEA_retention_v3
MF_BIN   := FEA_multifire_v3
SC_BIN   := FEA_secded_v3
CT_BIN   := FEA_crosstalk_v3
RST_BIN  := FEA_restoration_v3
REF_BIN  := FEA_refresh_v3
CLK_BIN  := FEA_clock_v3
BW_BIN   := FEA_bandwidth_v3
CMP_BIN  := FEA_compare_v3
FAB_BIN  := FEA_fabrication_v3
LAY_BIN  := FEA_layout_v3
PRG_BIN  := FEA_fzc_program_v3
REC_BIN  := FEA_recovery_v3
RESC_BIN := FEA_fzc_rescue_v3
THERM_BIN := FEA_thermal_v3

.PHONY: all clean run run-v1 run-v2 run-fzc run-slingshot run-fzc-reliability run-fzc-recognition run-fzc-opensystem run-fzc-selector run-fzc-screened run-fzc-floorplan run-fzc-e2e run-budget run-floorplan run-gamma run-retention run-multifire run-secded run-crosstalk run-restoration run-refresh run-clock run-bandwidth run-compare run-fabrication run-layout run-program v1 v2 fzc slingshot fzc-reliability fzc-recognition fzc-opensystem fzc-selector fzc-screened fzc-floorplan fzc-e2e budget floorplan gamma retention multifire secded crosstalk restoration refresh clock bandwidth compare fabrication layout program recovery run-recovery rescue run-rescue thermal run-thermal

all: $(V2_BIN)

v1: $(V1_BIN)
v2: $(V2_BIN)
fzc: $(FZC_BIN)
slingshot: $(SLING_BIN)
fzc-reliability: $(REL_BIN)
fzc-recognition: $(RECOG_BIN)
fzc-opensystem: $(OPEN_BIN)
fzc-selector: $(SEL_BIN)
fzc-screened: $(SCR_BIN)
fzc-floorplan: $(PLAN_BIN)
fzc-e2e: $(E2E_BIN)
budget: $(BUD_BIN)
floorplan: $(DIE_BIN)
gamma: $(GAM_BIN)
retention: $(RET_BIN)
multifire: $(MF_BIN)
secded: $(SC_BIN)
crosstalk: $(CT_BIN)
restoration: $(RST_BIN)
refresh: $(REF_BIN)
clock: $(CLK_BIN)
bandwidth: $(BW_BIN)
compare: $(CMP_BIN)
fabrication: $(FAB_BIN)
layout: $(LAY_BIN)
program: $(PRG_BIN)
recovery: $(REC_BIN)
rescue: $(RESC_BIN)
thermal: $(THERM_BIN)

$(V1_BIN): $(V1_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(V2_BIN): $(V2_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(FZC_BIN): $(FZC_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(SLING_BIN): $(SLING_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(REL_BIN): $(REL_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(RECOG_BIN): $(RECOG_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(OPEN_BIN): $(OPEN_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(SEL_BIN): $(SEL_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(SCR_BIN): $(SCR_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(PLAN_BIN): $(PLAN_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(E2E_BIN): $(E2E_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BUD_BIN): $(BUD_SRC)
	$(CXX) $(CXXFLAGS) -I$(V3DIR) -o $@ $<

$(DIE_BIN): $(DIE_SRC)
	$(CXX) $(CXXFLAGS) -I$(V3DIR) -o $@ $<

$(GAM_BIN): $(GAM_SRC)
	$(CXX) $(CXXFLAGS) -I$(V3DIR) -o $@ $<

$(RET_BIN): $(RET_SRC)
	$(CXX) $(CXXFLAGS) -I$(V3DIR) -o $@ $<

$(MF_BIN): $(MF_SRC)
	$(CXX) $(CXXFLAGS) -I$(V3DIR) -o $@ $<

$(SC_BIN): $(SC_SRC)
	$(CXX) $(CXXFLAGS) -I$(V3DIR) -o $@ $<

$(CT_BIN): $(CT_SRC)
	$(CXX) $(CXXFLAGS) -I$(V3DIR) -o $@ $<

$(RST_BIN): $(RST_SRC)
	$(CXX) $(CXXFLAGS) -I$(V3DIR) -o $@ $<

$(REF_BIN): $(REF_SRC)
	$(CXX) $(CXXFLAGS) -I$(V3DIR) -o $@ $<

$(CLK_BIN): $(CLK_SRC)
	$(CXX) $(CXXFLAGS) -I$(V3DIR) -o $@ $<

$(BW_BIN): $(BW_SRC)
	$(CXX) $(CXXFLAGS) -I$(V3DIR) -o $@ $<

$(CMP_BIN): $(CMP_SRC)
	$(CXX) $(CXXFLAGS) -I$(V3DIR) -o $@ $<

$(FAB_BIN): $(FAB_SRC)
	$(CXX) $(CXXFLAGS) -I$(V3DIR) -o $@ $<

$(LAY_BIN): $(LAY_SRC)
	$(CXX) $(CXXFLAGS) -I$(V3DIR) -o $@ $<

run: run-v2

run-v1: $(V1_BIN)
	./$(V1_BIN)

run-v2: $(V2_BIN)
	./$(V2_BIN)

run-fzc: $(FZC_BIN)
	./$(FZC_BIN)

run-slingshot: $(SLING_BIN)
	./$(SLING_BIN)

run-fzc-reliability: $(REL_BIN)
	./$(REL_BIN)

run-fzc-recognition: $(RECOG_BIN)
	./$(RECOG_BIN)

run-fzc-opensystem: $(OPEN_BIN)
	./$(OPEN_BIN)

run-fzc-selector: $(SEL_BIN)
	./$(SEL_BIN)

run-fzc-screened: $(SCR_BIN)
	./$(SCR_BIN)

run-fzc-floorplan: $(PLAN_BIN)
	./$(PLAN_BIN)

run-fzc-e2e: $(E2E_BIN)
	./$(E2E_BIN)

run-budget: $(BUD_BIN)
	./$(BUD_BIN)

run-floorplan: $(DIE_BIN)
	./$(DIE_BIN)

run-gamma: $(GAM_BIN)
	./$(GAM_BIN)

run-retention: $(RET_BIN)
	./$(RET_BIN)

run-multifire: $(MF_BIN)
	./$(MF_BIN)

run-secded: $(SC_BIN)
	./$(SC_BIN)

run-crosstalk: $(CT_BIN)
	./$(CT_BIN)

run-restoration: $(RST_BIN)
	./$(RST_BIN)

run-refresh: $(REF_BIN)
	./$(REF_BIN)

run-clock: $(CLK_BIN)
	./$(CLK_BIN)

run-bandwidth: $(BW_BIN)
	./$(BW_BIN)

run-compare: $(CMP_BIN)
	./$(CMP_BIN)

run-fabrication: $(FAB_BIN)
	./$(FAB_BIN)

run-layout: $(LAY_BIN)
	./$(LAY_BIN)

$(PRG_BIN): $(PRG_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(REC_BIN): $(REC_SRC)
	$(CXX) $(CXXFLAGS) -I$(V3DIR) -o $@ $<

$(RESC_BIN): $(RESC_SRC)
	$(CXX) $(CXXFLAGS) -I$(V3DIR) -o $@ $<

$(THERM_BIN): $(THERM_SRC)
	$(CXX) $(CXXFLAGS) -I$(V3DIR) -o $@ $<

run-program: $(PRG_BIN)
	./$(PRG_BIN)

run-recovery: $(REC_BIN)
	./$(REC_BIN)

run-rescue: $(RESC_BIN)
	./$(RESC_BIN)

run-thermal: $(THERM_BIN)
	./$(THERM_BIN)

clean:
	rm -f $(V1_BIN) $(V2_BIN) $(FZC_BIN) $(SLING_BIN) $(REL_BIN) $(RECOG_BIN) $(OPEN_BIN) $(SEL_BIN) $(SCR_BIN) $(PLAN_BIN) $(E2E_BIN) $(BUD_BIN) $(DIE_BIN) $(GAM_BIN) $(RET_BIN) $(MF_BIN) $(SC_BIN) $(CT_BIN) $(RST_BIN) $(REF_BIN) $(CLK_BIN) $(BW_BIN) $(CMP_BIN) $(FAB_BIN) $(LAY_BIN) $(PRG_BIN) $(REC_BIN) $(RESC_BIN) $(THERM_BIN)

# ---------------------------------------------------------------------------
# check: build and run EVERY target with one command, and fail if any fails.
#
# Reproducibility entry point: `make check` rebuilds the whole
# suite and returns non-zero if anything breaks. Each target either prints PASS
# itself or throws and exits non-zero, so only the exit status is inspected;
# output is shown only for a failing target (nothing is written to a temp file,
# because /tmp is not assumed to be writable).
# ---------------------------------------------------------------------------
V3_CHECK_TARGETS := run-v1 run-v2 run-fzc run-slingshot run-fzc-reliability run-fzc-recognition run-fzc-opensystem run-fzc-selector run-fzc-screened run-fzc-floorplan run-fzc-e2e run-budget run-floorplan run-gamma run-retention run-multifire run-secded run-crosstalk run-restoration run-refresh run-clock run-bandwidth run-compare run-fabrication run-layout run-program run-recovery run-rescue run-thermal

check:
	@fail=0; \
	for t in $(V3_CHECK_TARGETS); do \
		printf '  %-24s ' "$$t"; \
		if $(MAKE) --no-print-directory "$$t" >/dev/null 2>&1; then \
			echo "PASS"; \
		else \
			echo "FAIL"; \
			$(MAKE) --no-print-directory "$$t" 2>&1 | tail -n 25 | sed 's/^/        /'; \
			fail=1; \
		fi; \
	done; \
	if [ $$fail -eq 0 ]; then echo "  ALL 29 TARGETS PASS"; \
	else echo "  SOME TARGETS FAILED"; fi; \
	exit $$fail

.PHONY: check
