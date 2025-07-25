# Python config
CONDA_ENV_BURNMAN := burnman
CONDA_ENV_PYVISTA := pyvista
CONDA_PYTHON_BURNMAN := $$(conda run -n $(CONDA_ENV_BURNMAN) which python)
CONDA_PYTHON_PYVISTA := $$(conda run -n $(CONDA_ENV_PYVISTA) which python)

# Top-level dirs
PROJECT_ROOT := $(realpath $(dir $(lastword $(MAKEFILE_LIST)))/..)
PERPLEX := $(PROJECT_ROOT)/perplex
PYTHON := $(PROJECT_ROOT)/python
SIMULATION := $(PROJECT_ROOT)/simulation
BASH := $(PROJECT_ROOT)/bash
FIGS := $(PROJECT_ROOT)/figs

# Perplex models
PERPLEX_MODELS := $(shell find $(PERPLEX) -type f -name '*.dat' -exec basename {} .dat \;)

# Simulations
2D_BOX_CONFIGS := $(patsubst $(SIMULATION)/%,%,$(shell find $(SIMULATION)/2d_box/configs -type f -name '*.prm'))
# 2D_BOX_CONFIGS := 2d_box/configs/simple-slab-Q9.prm 2d_box/configs/simple-plume-Q9.prm
2D_BOX_TIMESTEPS := 00 01 02 03 04 05 06 07 08 09 10 20 30 40 50 60 70 80 90
2D_SHELL_CONFIGS := $(patsubst $(SIMULATION)/%,%,$(shell find $(SIMULATION)/2d_shell/configs -type f -name '*.prm'))
2D_SHELL_TIMESTEPS := 00 01 02 03 04 05 06 07 08 09 10 20 30 40 50 60 70 80 90

# Script paths
MAIN_ADIABAT := $(PYTHON)/adiabat/main.py
MAIN_SIMULATION := $(PYTHON)/simulation/main.py

# Cleanup targets
CLEAN := $(PYTHON)/**/__pycache__
DEEP_CLEAN :=

# Targets
.PHONY: all visualize 2d-box-models 2d-shell-models adiabat clean deep-clean help

all: adiabat

visualize: adiabat 2d-box-models 2d-shell-models

sync-barkla:
	@source $(BASH)/simulation/sync-barkla2.sh "$(2D_BOX_CONFIGS)" "$(2D_BOX_TIMESTEPS)"


2d-box-models: $(MAIN_SIMULATION)
	@for prm in $(2D_BOX_CONFIGS); do \
		MODEL_DIR=$$(basename $$prm .prm | tr '-' '_'); \
		IN_DIR=$(SIMULATION)/2d_box/results/$$MODEL_DIR; \
		TSTEPS="$(2D_BOX_TIMESTEPS)"; \
		OUT_FIG_DIR=$(FIGS)/simulation/$$MODEL_DIR; \
		VIS_TYPE=2d_box; \
		$(CONDA_PYTHON_PYVISTA) -u $(MAIN_SIMULATION) \
		    --model-id "$$MODEL_DIR" \
		    --timesteps "$$TSTEPS" \
		    --in-dir "$$IN_DIR" \
		    --out-fig-dir "$$OUT_FIG_DIR" \
		    --visualization-type "$$VIS_TYPE"; \
	done

2d-shell-models: $(MAIN_SIMULATION)
	@for prm in $(2D_SHELL_CONFIGS); do \
		MODEL_DIR=$$(basename $$prm .prm | tr '-' '_'); \
		IN_DIR=$(SIMULATION)/2d_shell/results/$$MODEL_DIR; \
		TSTEPS="$(2D_SHELL_TIMESTEPS)"; \
		OUT_FIG_DIR="$(FIGS)/simulation/$$MODEL_DIR"; \
		VIS_TYPE="2d_shell"; \
		$(CONDA_PYTHON_PYVISTA) -u $(MAIN_SIMULATION) \
		    --model-id "$$MODEL_DIR" \
		    --timesteps "$$TSTEPS" \
		    --in-dir "$$IN_DIR" \
		    --out-fig-dir "$$OUT_FIG_DIR" \
		    --visualization-type "$$VIS_TYPE"; \
	done

adiabat: $(MAIN_ADIABAT)
	@for model in $(PERPLEX_MODELS); do \
		IN_DIR="$(PERPLEX)/$$model"; \
		OUT_DIR="$(SIMULATION)/data"; \
		OUT_FIG_DIR="$(FIGS)/adiabat"; \
		(cd $$IN_DIR && \
		    $(CONDA_PYTHON_BURNMAN) -u $(MAIN_ADIABAT) \
		        --model-id "$$model" \
		        --out-dir "$$OUT_DIR" \
		        --out-fig-dir "$$OUT_FIG_DIR" \
		); \
	done

clean:
	@for item in $(CLEAN); do \
		safe_rm() { \
			if [ -e "$$1" ]; then \
				ABS_PATH=$$(realpath "$$1"); \
				case "$$ABS_PATH" in \
					$(PROJECT_ROOT)*) \
						echo "-> Safely removing\n   $$ABS_PATH"; \
						rm -rf "$$ABS_PATH";; \
					*) \
						echo "!! Skipping (outside project root)\n   $$ABS_PATH";; \
				esac; \
			else \
				echo "-- Skipping (not found)\n   $$1"; \
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
						echo "-> Safely removing\n   $$ABS_PATH"; \
						rm -rf "$$ABS_PATH";; \
					*) \
						echo "!! Skipping (outside project root)\n   $$ABS_PATH";; \
				esac; \
			else \
				echo "-- Skipping (not found)\n   $$1"; \
			fi; \
		}; \
		safe_rm "$$item"; \
	done

help:
	@echo "Available targets:"
	@echo "  visualize        Visualize all results"
	@echo "  2d-shell-models  Visualize ASPECT shell models"
	@echo "  2d-box-models    Visualize ASPECT box models"
	@echo "  adiabat          Generate adiabatic profiles using BurnMan"
	@echo "  clean            Cleanup unnecessary files and directories (safe)"
	@echo "  deep-clean       Deep clean figures and results (use with caution!)"
	@echo "  help             Show this help message"
