# Top-level dirs
PROJECT_ROOT := $(realpath $(dir $(lastword $(MAKEFILE_LIST)))/..)
BASH := $(PROJECT_ROOT)/bash
SIMULATION := $(PROJECT_ROOT)/simulation

# Check OS
UNAME_S := $(shell uname)
OS_ID := $(shell . /etc/os-release 2>/dev/null && echo $$ID)

# ASPECT config
ifeq ($(UNAME_S), Darwin)
  RUN_ASPECT := $(BASH)/simulation/run-aspect-macos.sh
  NTASKS := 6
else ifeq ($(OS_ID), rocky)
  RUN_ASPECT := $(BASH)/simulation/run-aspect-barkla2.sh
  GCC := gcc/14.2.0
  OPENMPI := openmpi/5.0.8-gcc14.2.0
  OPENBLAS := openblas/0.3.29/gcc-14.2.0
  CMAKE := cmake/3.30.5-gcc14.2.0
  NTASKS := 128
  NTASKSPERNODE := 168
  TIME := 2-12:00:00
  MAIL_ADDRESS := b.kerswell@liverpool.ac.uk
  MAIL_TYPE := ALL
else
  $(error Unsupported OS: Only macOS and Rocky Linux are supported!)
endif

# ASPECT executable
ASPECT_EXEC := $(PROJECT_ROOT)/aspect/build/aspect-release

# Material parameters
Z_FACTOR := 3.0e0 7.0e0 1.6e1 3.7e1 8.7e1 2.0e2 4.7e2 1.1e3 2.6e3 6.0e3 1.4e4 3.3e4 7.8e4 1.8e5 4.3e5 1.0e6 2.4e6 5.6e6 1.3e7 3.0e7 7.0e7
B_FACTOR := 2 4 6 8 10
ZB_COMBO := $(foreach z, $(Z_FACTOR), $(foreach b, $(B_FACTOR), Z$(z)-B$(b)))

# Templates prms
MATERIAL_TEMPLATE := $(SIMULATION)/parts/template-material.prm
PLUME_TEMPLATE := $(SIMULATION)/parts/template-plume.prm
SLAB_TEMPLATE := $(SIMULATION)/parts/template-slab.prm

# Generate prm files for each kinetic factor
MATERIAL_PRMS := $(addprefix $(SIMULATION)/parts/material-, $(addsuffix .prm, $(ZB_COMBO)))
PLUME_PRMS := $(addprefix $(SIMULATION)/configs/plume-, $(addsuffix .prm, $(ZB_COMBO)))
SLAB_PRMS := $(addprefix $(SIMULATION)/configs/slab-, $(addsuffix .prm, $(ZB_COMBO)))

# Simulation data
DATA := $(SIMULATION)/data/reaction-410-profile.txt

# Simulation plugins
PLUGINS := plugins/libphase-transition-kinetics.release.so

# Cleanup targets
CLEAN := \
	$(SIMULATION)/**/plugins/CMakeFiles \
	$(SIMULATION)/**/plugins/Makefile \
	$(SIMULATION)/**/plugins/CMakeCache.txt \
	$(SIMULATION)/**/plugins/cmake_install.cmake \
	$(SIMULATION)/**/plugins/.cache \
	$(MATERIAL_PRMS) \
	$(PLUME_PRMS) \
	$(SLAB_PRMS)
DEEP_CLEAN := \
	$(SIMULATION)/**/configs/logs

# Targets
.PHONY: all models prms clean deep-clean help

# Keep intermediate targets
.PRECIOUS: $(MATERIAL_PRMS) $(PLUME_PRMS) $(SLAB_PRMS)

all: models

models: plume-models slab-models

plume-models: $(ASPECT_EXEC) $(PLUGINS) $(DATA) $(MATERIAL_PRMS) $(PLUME_PRMS:.prm=.run)

slab-models: $(ASPECT_EXEC) $(PLUGINS) $(DATA) $(MATERIAL_PRMS) $(SLAB_PRMS:.prm=.run)

%.run: %.prm
	@ROOT_DIR=$(SIMULATION); \
	PRM=$<; \
	RUN_OUT=$@; \
	PRM_BASENAME=$$(basename $$PRM .prm); \
	MODEL_TYPE=$$(echo $$PRM_BASENAME | cut -d- -f1); \
	MODEL_DIR=$$(echo $$PRM_BASENAME | sed -E 's/-/_/g'); \
	RESULTS_DIR="$$ROOT_DIR/results"; \
	JOB_NAME=$$(echo $$PRM_BASENAME | tr '/_' '-'); \
	LOG_OUT="$$RESULTS_DIR/$$MODEL_DIR/$$JOB_NAME.out"; \
	if [ -f "$$RUN_OUT" ]; then \
		echo " -- $$RUN_OUT found!"; \
	else \
		if [ "$(OS_ID)" = "rocky" ]; then \
			sbatch --job-name="$$JOB_NAME" \
			       --output="$$LOG_OUT" \
			       --partition="nodes" \
			       --ntasks=$(NTASKS) \
			       --ntasks-per-node=$(NTASKSPERNODE) \
			       --time=$(TIME) \
			       --mail-user=$(MAIL_ADDRESS) \
			       --mail-type=$(MAIL_TYPE) \
			       $(RUN_ASPECT) $(ASPECT_EXEC) "$$PRM"; \
		else \
			mkdir -p $$(dirname "$$LOG_OUT"); \
			$(RUN_ASPECT) $(NTASKS) $(ASPECT_EXEC) "$$PRM" \
			2>&1 | tee -a "$$LOG_OUT"; \
		fi \
	fi

$(MATERIAL_PRMS): $(SIMULATION)/parts/material-%.prm: $(MATERIAL_TEMPLATE)
	@mkdir -p $(dir $@)
	@z=$$(echo $* | cut -d- -f1 | sed 's/^Z//'); \
	b=$$(echo $* | cut -d- -f2 | sed 's/^B//'); \
	sed "s/{z-factor}/$$z/g; s/{b-factor}/$$b/g" $< > $@

$(PLUME_PRMS): $(SIMULATION)/configs/plume-%.prm: $(PLUME_TEMPLATE)
	@mkdir -p $(dir $@)
	@z=$$(echo $* | cut -d- -f1 | sed 's/^Z//'); \
	b=$$(echo $* | cut -d- -f2 | sed 's/^B//'); \
	sed "s/{z-factor}/$$z/g; s/{b-factor}/$$b/g" $< > $@

$(SLAB_PRMS): $(SIMULATION)/configs/slab-%.prm: $(SLAB_TEMPLATE)
	@mkdir -p $(dir $@)
	@z=$$(echo $* | cut -d- -f1 | sed 's/^Z//'); \
	b=$$(echo $* | cut -d- -f2 | sed 's/^B//'); \
	sed "s/{z-factor}/$$z/g; s/{b-factor}/$$b/g" $< > $@

$(PLUGINS):
	@if [ "$(OS_ID)" = "rocky" ]; then \
		module load $(GCC) $(OPENMPI) $(OPENBLAS) $(CMAKE); \
		cd $(SIMULATION)/plugins; \
		cmake . -D CMAKE_EXPORT_COMPILE_COMMANDS=ON -D Aspect_DIR=$(dir $(ASPECT_EXEC)) -D CMAKE_CXX_COMPILER=mpicxx -W no-dev; \
		make; \
	else \
		cd $(SIMULATION)/plugins; \
		cmake . -D CMAKE_EXPORT_COMPILE_COMMANDS=ON -D Aspect_DIR=$(dir $(ASPECT_EXEC)) -W no-dev; \
		make; \
	fi

$(SIMULATION)/data/%.txt:
	@$(MAKE) --no-print-directory -C $(PYTHON) adiabat

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

	@find $(PROJECT_ROOT)/simulation -type l -path '*/plugins/aspect*' | while read -r item; do \
		case "$$item" in \
			$(PROJECT_ROOT)*) \
				echo " -> Safely removing symlink\n   $$item"; \
				rm -v "$$item";; \
			*) \
				echo " -- Skipping (outside project root)\n   $$item";; \
		esac; \
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
	@echo "    models        Run ASPECT 2d box models"
	@echo "    plume-models  Run ASPECT 2d box plume models"
	@echo "    slab-models   Run ASPECT 2d box slab models"
	@echo "    clean         Cleanup unnecessary files and directories (safe)"
	@echo "    deep-clean    Deep clean figures and results (use with caution!)"
	@echo "    help          Show this help message"
