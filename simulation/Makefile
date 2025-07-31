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

# Simulation data
2D_BOX_DATA := \
	$(SIMULATION)/data/forward-reaction-410-profile.txt \
	$(SIMULATION)/data/reverse-reaction-410-profile.txt
2D_SHELL_DATA := \
	$(SIMULATION)/data/PYR-adiabatic-profile.txt \
	$(SIMULATION)/data/PYR-material-table.txt \
	$(SIMULATION)/data/forward-reaction-410-profile.txt \
	$(SIMULATION)/data/reverse-reaction-410-profile.txt

# Simulation configs
2D_BOX_CONFIGS := $(patsubst $(SIMULATION)/%,%,$(shell find $(SIMULATION)/2d_box/configs -type f -name '*.prm'))
2D_SHELL_CONFIGS := $(patsubst $(SIMULATION)/%,%,$(shell find $(SIMULATION)/2d_shell/configs -type f -name '*.prm'))

# Simulation plugins
2D_BOX_PLUGINS := 2d_box/plugins/libphase-transition-kinetics.release.so
2D_SHELL_PLUGINS :=

# Cleanup targets
CLEAN := \
	$(SIMULATION)/**/plugins/CMakeFiles \
	$(SIMULATION)/**/plugins/Makefile \
	$(SIMULATION)/**/plugins/CMakeCache.txt \
	$(SIMULATION)/**/plugins/cmake_install.cmake \
	$(SIMULATION)/**/plugins/.cache
DEEP_CLEAN := \
	$(SIMULATION)/**/configs/logs

# Targets
.PHONY: all 2d-box-models 2d-shell-models clean deep-clean help

all: 2d-box-models 2d-shell-models

2d-shell-models: $(ASPECT_EXEC) $(2D_SHELL_PLUGINS) $(2D_SHELL_DATA) $(2D_SHELL_CONFIGS:.prm=.run)

2d-box-models: $(ASPECT_EXEC) $(2D_BOX_PLUGINS) $(2D_BOX_DATA) $(2D_BOX_CONFIGS:.prm=.run)

%.run: %.prm
	@ROOT_DIR=$(SIMULATION); \
	PRM=$<; \
	RUN_OUT=$@; \
	MODEL_TYPE=$$(echo $$PRM | cut -d/ -f1); \
	MODEL_DIR=$$(basename $$PRM .prm | sed -E 's/-/_/g; s/(lnk)_([+-]?[0-9]+)/\1-\2/g'); \
	PRM_DIR="$$ROOT_DIR/$$MODEL_TYPE/configs"; \
	RESULTS_DIR="$$ROOT_DIR/$$MODEL_TYPE/results"; \
	JOB_NAME=$$(echo $* | tr '/_' '-'); \
	LOG_OUT="$$RESULTS_DIR/$$MODEL_DIR/aspect.$$JOB_NAME.out"; \
	if [ -f "$$ROOT_DIR/$$RUN_OUT" ]; then \
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
			       $(RUN_ASPECT) $(ASPECT_EXEC) "$$ROOT_DIR/$$PRM"; \
		else \
			mkdir -p $$(dirname "$$LOG_OUT"); \
			$(RUN_ASPECT) $(NTASKS) $(ASPECT_EXEC) "$$ROOT_DIR/$$PRM" \
			2>&1 | tee -a "$$LOG_OUT"; \
		fi \
	fi

$(2D_BOX_PLUGINS):
	@if [ "$(OS_ID)" = "rocky" ]; then \
		module load $(GCC) $(OPENMPI) $(OPENBLAS) $(CMAKE); \
		cd $(SIMULATION)/2d_box/plugins; \
		cmake . -D CMAKE_EXPORT_COMPILE_COMMANDS=ON -D Aspect_DIR=$(dir $(ASPECT_EXEC)) -D CMAKE_CXX_COMPILER=mpicxx -W no-dev; \
		make; \
	else \
		cd $(SIMULATION)/2d_box/plugins; \
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
	@echo "    2d-shell-models  Run ASPECT shell models"
	@echo "    2d-box-models    Run ASPECT box models"
	@echo "    clean            Cleanup unnecessary files and directories (safe)"
	@echo "    deep-clean       Deep clean figures and results (use with caution!)"
	@echo "    help             Show this help message"

