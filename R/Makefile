# Top-level dirs
PROJECT_ROOT := $(realpath $(dir $(lastword $(MAKEFILE_LIST)))/..)
R := $(PROJECT_ROOT)/R
PYTHON := $(PROJECT_ROOT)/python
PERPLEX := $(PROJECT_ROOT)/perplex
SIMULATION := $(PROJECT_ROOT)/simulation
FIGS := $(PROJECT_ROOT)/figs
MARKDOWN := $(PROJECT_ROOT)/markdown

# Perplex models
PERPLEX_MODELS := $(shell find $(PERPLEX) -type f -name '*.dat' -exec basename {} .dat \;)

# Material parameters
KINETIC_FACTOR := 3.0e0 7.0e0 1.6e1 3.7e1 8.7e1 2.0e2 4.7e2 1.1e3 2.6e3 6.0e3 1.4e4 3.3e4 7.8e4 1.8e5 4.3e5 1.0e6 2.4e6 5.6e6 1.3e7 3.0e7 7.0e7

# Generate prm files for each kinetic factor
PLUME_PRMS := $(addprefix $(SIMULATION)/configs/plume-, $(addsuffix .prm, $(KINETIC_FACTOR)))
SLAB_PRMS := $(addprefix $(SIMULATION)/configs/slab-, $(addsuffix .prm, $(KINETIC_FACTOR)))
PRMS := $(PLUME_PRMS) $(SLAB_PRMS)

# Filepaths
DEPTH_AVERAGE_PATHS := $(shell find $(SIMULATION) -type f -path "*/results/*/depth_average.txt")

# Cleanup targets
CLEAN :=
DEEP_CLEAN :=

# Targets
.PHONY: all visualize adiabat simulation clean deep-clean help

all: visualize

visualize: adiabat simulation

adiabat: $(R)/adiabat/main.R
	@for model in $(PERPLEX_MODELS); do \
		IN_DIR=$(SIMULATION)/data; \
		OUT_FIG_DIR=$(FIGS); \
		Rscript $(R)/adiabat/main.R $$model $$IN_DIR $$OUT_FIG_DIR; \
	done

simulation: $(R)/simulation/main.R
	@Rscript $(R)/simulation/main.R $(PYTHON)/simulation/data $(FIGS)

clean:
	@for item in $(CLEAN); do \
		safe_rm() { \
			if [ -e "$$1" ]; then \
				ABS_PATH=$$(realpath "$$1"); \
				case "$$ABS_PATH" in \
					$(PROJECT_ROOT)*) \
						echo " -> Safely removing\n   $$ABS_PATH"; \
						rm -rf "$$ABS_PATH";; \
					*) \
						echo " -- Skipping (outside project root)\n   $$ABS_PATH";; \
				esac; \
			else \
				echo " -- Skipping (not found)\n   $$1"; \
			fi; \
		}; \
		safe_rm "$$item"; \
	done

deep-clean:
	@for item in $(DEEP_CLEAN); do \
		safe_rm() { \
			if [ -e "$$1" ]; then \
				ABS_PATH=$$(realpath "$$1"); \
				case "$$ABS_PATH" in \
					$(PROJECT_ROOT)*) \
						echo " -> Safely removing\n   $$ABS_PATH"; \
						rm -rf "$$ABS_PATH";; \
					*) \
						echo " -- Skipping (outside project root)\n   $$ABS_PATH";; \
				esac; \
			else \
				echo " -- Skipping (not found)\n   $$1"; \
			fi; \
		}; \
		safe_rm "$$item"; \
	done

help:
	@echo "    --------------------------------------------------"
	@echo "    Available targets:"
	@echo "    --------------------------------------------------"
	@echo "    visualize      Visualize all results"
	@echo "    adiabat        Visualize adiabatic profiles"
	@echo "    simulation     Visualize simulation results"
	@echo "    clean          Cleanup unnecessary files and directories (safe)"
	@echo "    deep-clean     Deep clean figures and results (use with caution!)"
	@echo "    help           Show this help message"
