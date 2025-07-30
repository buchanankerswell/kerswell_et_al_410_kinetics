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
# 2D_BOX_CONFIGS := 2d_box/configs/simple-slab-Q08.prm 2d_box/configs/simple-slab-Q09.prm
2D_BOX_TIMESTEPS_SYNC := $(shell seq -f "%03g" 0 1 150)
2D_BOX_TIMESTEPS_MESH := $(shell seq -f "%03g" 0 10 50)
2D_BOX_TIMESTEPS_PROFILE := $(shell seq -f "%03g" 0 10 50)
2D_BOX_TIMESTEPS_TOPOGRAPHY := $(shell seq -f "%03g" 0 1 150)
2D_SHELL_CONFIGS := $(patsubst $(SIMULATION)/%,%,$(shell find $(SIMULATION)/2d_shell/configs -type f -name '*.prm'))
2D_SHELL_TIMESTEPS_SYNC := $(shell seq -f "%03g" 0 1 150)
2D_SHELL_TIMESTEPS_MESH := $(shell seq -f "%03g" 0 10 50)
2D_SHELL_TIMESTEPS_PROFILE := $(shell seq -f "%03g" 0 10 50)
2D_SHELL_TIMESTEPS_TOPOGRAPHY := $(shell seq -f "%03g" 0 1 150)

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
	@source $(BASH)/simulation/sync-barkla2.sh "$(2D_BOX_CONFIGS)" "$(2D_BOX_TIMESTEPS_SYNC)"


2d-box-models: $(MAIN_SIMULATION)
	@MODEL_IDS=""; \
	IN_DIRS=""; \
	for prm in $(2D_BOX_CONFIGS); do \
		MODEL_DIR=$$(basename $$prm .prm | tr '-' '_'); \
		IN_DIR=$(SIMULATION)/2d_box/results/$$MODEL_DIR; \
		MODEL_IDS="$$MODEL_IDS $$MODEL_DIR"; \
		IN_DIRS="$$IN_DIRS $$IN_DIR"; \
	done; \
	TSTEPS_MESH="$(2D_BOX_TIMESTEPS_MESH)"; \
	TSTEPS_PROFILE="$(2D_BOX_TIMESTEPS_PROFILE)"; \
	TSTEPS_TOPOGRAPHY="$(2D_BOX_TIMESTEPS_TOPOGRAPHY)"; \
	OUT_FIG_DIR=$(FIGS)/simulation/2d_box; \
	VIS_TYPE=2d_box; \
	$(CONDA_PYTHON_PYVISTA) -u $(MAIN_SIMULATION) \
	    --model-ids $$MODEL_IDS \
	    --timesteps-mesh $$TSTEPS_MESH \
	    --timesteps-profile $$TSTEPS_PROFILE \
	    --timesteps-topography $$TSTEPS_TOPOGRAPHY \
	    --in-dirs $$IN_DIRS \
	    --out-fig-dir $$OUT_FIG_DIR \
	    --visualization-type $$VIS_TYPE

2d-shell-models: $(MAIN_SIMULATION)
	@MODEL_IDS=""; \
	IN_DIRS=""; \
	for prm in $(2D_SHELL_CONFIGS); do \
		MODEL_DIR=$$(basename $$prm .prm | tr '-' '_'); \
		IN_DIR=$(SIMULATION)/2d_shell/results/$$MODEL_DIR; \
		MODEL_IDS="$$MODEL_IDS $$MODEL_DIR"; \
		IN_DIRS="$$IN_DIRS $$IN_DIR"; \
	done; \
	TSTEPS_MESH="$(2D_SHELL_TIMESTEPS_MESH)"; \
	TSTEPS_PROFILE="$(2D_SHELL_TIMESTEPS_PROFILE)"; \
	TSTEPS_TOPOGRAPHY="$(2D_SHELL_TIMESTEPS_TOPOGRAPHY)"; \
	OUT_FIG_DIR=$(FIGS)/simulation/2d_shell; \
	VIS_TYPE=2d_shell; \
	$(CONDA_PYTHON_PYVISTA) -u $(MAIN_SIMULATION) \
	    --model-ids $$MODEL_IDS \
	    --timesteps-mesh $$TSTEPS_MESH \
	    --timesteps-profile $$TSTEPS_PROFILE \
	    --timesteps-topography $$TSTEPS_TOPOGRAPHY \
	    --in-dirs $$IN_DIRS \
	    --out-fig-dir $$OUT_FIG_DIR \
	    --visualization-type $$VIS_TYPE

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
