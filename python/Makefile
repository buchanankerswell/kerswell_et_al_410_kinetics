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

# Remote sim dir
SIMULATION_BARKLA := /users/kersweb/scratch/kerswell_et_al_dynp/simulation

# Perplex models
PERPLEX_MODELS := $(shell find $(PERPLEX) -type f -name '*.dat' -exec basename {} .dat \;)

# Material parameters
Z_FACTOR := 3.0e0 7.0e0 1.6e1 3.7e1 8.7e1 2.0e2 4.7e2 1.1e3 2.6e3 6.0e3 1.4e4 3.3e4 7.8e4 1.8e5 4.3e5 1.0e6 2.4e6 5.6e6 1.3e7 3.0e7 7.0e7
B_FACTOR := 2 4 6 8 10
ZB_COMBO := $(foreach z, $(Z_FACTOR), $(foreach b, $(B_FACTOR), Z$(z)-B$(b)))

# Generate prm files for each kinetic factor
PLUME_PRMS := $(addprefix $(SIMULATION)/configs/plume-, $(addsuffix .prm, $(ZB_COMBO)))
SLAB_PRMS := $(addprefix $(SIMULATION)/configs/slab-, $(addsuffix .prm, $(ZB_COMBO)))
PRMS := $(PLUME_PRMS) $(SLAB_PRMS)

# Timesteps
# TIMESTEPS_SYNC := $(shell seq -f "%03g" 0 10 1)
TIMESTEPS_SYNC := 10
TIMESTEPS_MESH := 10
TIMESTEPS_PROFILE := 10
TIMESTEPS_TOPOGRAPHY := 10

# Script paths
MAIN_ADIABAT := $(PYTHON)/adiabat/main.py
MAIN_SIMULATION := $(PYTHON)/simulation/main.py
INITIAL_SETUP := $(PYTHON)/simulation/initial-setup.py
WRITE_TABLE := $(PYTHON)/simulation/write-depth-profile-summary-table.py

# Cleanup targets
CLEAN := $(PYTHON)/**/__pycache__
DEEP_CLEAN := $(PYTHON)/*/data

# Targets
.PHONY: all visualize sync-barkla compositions initial-setup meshes adiabat clean deep-clean help

all: adiabat visualize

visualize: initial-setup meshes compositions

sync-barkla:
	@source $(BASH)/simulation/sync-barkla2.sh "$(PRMS)" "$(TIMESTEPS_SYNC)" "$(SIMULATION)" "$(SIMULATION_BARKLA)"

write-markdown-table:
	@$(CONDA_PYTHON_PYVISTA) -u $(WRITE_TABLE) --in-dir $(PYTHON)/simulation/data

compositions:
	@(cd $(BASH)/simulation && ./make-compositions.sh "$(FIGS)/simulation/meshes" "$(FIGS)/simulation/compositions" "0010" "Z3.0e0-B2" "Z4.7e2-B2" "Z7.0e7-B2")
	@(cd $(BASH)/simulation && ./make-compositions.sh "$(FIGS)/simulation/meshes" "$(FIGS)/simulation/compositions" "0010" "Z3.0e0-B4" "Z4.7e2-B4" "Z7.0e7-B4")
	@(cd $(BASH)/simulation && ./make-compositions.sh "$(FIGS)/simulation/meshes" "$(FIGS)/simulation/compositions" "0010" "Z3.0e0-B6" "Z4.7e2-B6" "Z7.0e7-B6")
	@(cd $(BASH)/simulation && ./make-compositions.sh "$(FIGS)/simulation/meshes" "$(FIGS)/simulation/compositions" "0010" "Z3.0e0-B8" "Z4.7e2-B8" "Z7.0e7-B8")
	@(cd $(BASH)/simulation && ./make-compositions.sh "$(FIGS)/simulation/meshes" "$(FIGS)/simulation/compositions" "0010" "Z3.0e0-B10" "Z4.7e2-B10" "Z7.0e7-B10")
	@(cd $(BASH)/simulation && ./make-compositions.sh "$(FIGS)/simulation/meshes" "$(FIGS)/simulation/compositions" "0010" "Z3.0e0-B2" "Z3.0e0-B6" "Z3.0e0-B10")
	@(cd $(BASH)/simulation && ./make-compositions.sh "$(FIGS)/simulation/meshes" "$(FIGS)/simulation/compositions" "0010" "Z4.7e2-B2" "Z4.7e2-B6" "Z4.7e2-B10")
	@(cd $(BASH)/simulation && ./make-compositions.sh "$(FIGS)/simulation/meshes" "$(FIGS)/simulation/compositions" "0010" "Z7.0e7-B2" "Z7.0e7-B6" "Z7.0e7-B10")

initial-setup: $(INITIAL_SETUP)
	@$(CONDA_PYTHON_PYVISTA) -u $(INITIAL_SETUP) --out-fig-dir $(FIGS)

meshes: $(MAIN_SIMULATION)
	@MODEL_IDS=""; \
	IN_DIRS=""; \
	for prm in $(PRMS); do \
		MODEL_DIR=$$(basename $$prm .prm | sed -E 's/-/_/g; s/(lnk)_([+-]?[0-9]+)/\1-\2/g'); \
		IN_DIR=$(SIMULATION)/results/$$MODEL_DIR; \
		MODEL_IDS="$$MODEL_IDS $$MODEL_DIR"; \
		IN_DIRS="$$IN_DIRS $$IN_DIR"; \
	done; \
	TSTEPS_MESH="$(TIMESTEPS_MESH)"; \
	TSTEPS_PROFILE="$(TIMESTEPS_PROFILE)"; \
	TSTEPS_TOPOGRAPHY="$(TIMESTEPS_TOPOGRAPHY)"; \
	OUT_FIG_DIR=$(FIGS)/simulation; \
	$(CONDA_PYTHON_PYVISTA) -u $(MAIN_SIMULATION) \
	    --model-ids $$MODEL_IDS \
	    --timesteps-mesh $$TSTEPS_MESH \
	    --timesteps-profile $$TSTEPS_PROFILE \
	    --timesteps-topography $$TSTEPS_TOPOGRAPHY \
	    --in-dirs $$IN_DIRS \
	    --out-fig-dir $$OUT_FIG_DIR

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
				if [[ "$$ABS_PATH" == $(PROJECT_ROOT)* ]]; then \
					echo " -> Safely removing: $$ABS_PATH"; \
					rm -rf "$$ABS_PATH"; \
				fi; \
			fi; \
		}; \
		safe_rm "$$item"; \
	done

deep-clean:
	@for item in $(DEEP_CLEAN); do \
		safe_rm() { \
			if [ -e "$$1" ]; then \
				ABS_PATH=$$(realpath "$$1"); \
				if [[ "$$ABS_PATH" == $(PROJECT_ROOT)* ]]; then \
					echo " -> Safely removing: $$ABS_PATH"; \
					rm -rf "$$ABS_PATH"; \
				fi; \
			fi; \
		}; \
		safe_rm "$$item"; \
	done

help:
	@echo "    --------------------------------------------------"
	@echo "    Available targets:"
	@echo "    --------------------------------------------------"
	@echo "    visualize         Visualize all results"
	@echo "    meshes            Visualize ASPECT models"
	@echo "    compositions      Compose individual simulation plots"
	@echo "    initial-setup     Visualize initial numerical setup"
	@echo "    sync-barkla       Fetch simulation results from barkla"
	@echo "    adiabat           Generate adiabatic material profiles using BurnMan"
	@echo "    clean             Cleanup unnecessary files and directories (safe)"
	@echo "    deep-clean        Deep clean figures and results (use with caution!)"
	@echo "    help              Show this help message"
