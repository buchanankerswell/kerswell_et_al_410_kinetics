# Top-level dirs
PROJECT_ROOT := $(realpath $(dir $(lastword $(MAKEFILE_LIST)))/..)
R := $(PROJECT_ROOT)/R
PERPLEX := $(PROJECT_ROOT)/perplex
SIMULATION := $(PROJECT_ROOT)/simulation
FIGS := $(PROJECT_ROOT)/figs
MARKDOWN := $(PROJECT_ROOT)/markdown

# Perplex models
PERPLEX_MODELS := $(shell find $(PERPLEX) -type f -name '*.dat' -exec basename {} .dat \;)

# Simulations
2D_BOX_MODELS := $(shell find $(SIMULATION)/2d_box/configs -type f -name '*.prm' -exec basename {} .prm \;)
2D_SHELL_MODELS := $(shell find $(SIMULATION)/2d_shell/configs -type f -name '*.prm' -exec basename {} .prm \;)

# Filepaths
DEPTH_AVERAGE_PATHS := $(shell find $(SIMULATION) -type f -path "*/results/*/depth_average.txt")

# Cleanup targets
CLEAN :=
DEEP_CLEAN :=

# Targets
.PHONY: all visualize adiabat 2d-box-models 2d-shell-models misc manuscript clean deep-clean help

all: visualize

visualize: adiabat 2d-box-models 2d-shell-models misc

adiabat: $(R)/adiabat/main.R
	@for model in $(PERPLEX_MODELS); do \
		IN_DIR=$(SIMULATION)/data; \
		OUT_FIG_DIR=$(FIGS)/perplex/$$model; \
		Rscript $(R)/adiabat/main.R $$model $$IN_DIR $$OUT_FIG_DIR; \
	done

2d-box-models: $(R)/simulation/main.R
	@for model in $(2D_BOX_MODELS); do \
		DIR_NAME=$$(echo $$model | tr '/-' '_'); \
		IN_DIR=$(SIMULATION)/2d_box/results/$$DIR_NAME; \
		OUT_FIG_DIR=$(FIGS)/simulation/$$DIR_NAME; \
		VIS_TYPE=2d_box; \
		if [ -d "$$IN_DIR" ]; then \
			Rscript $(R)/simulation/main.R $$IN_DIR $$OUT_FIG_DIR; \
		fi \
	done

2d-shell-models: $(R)/simulation/main.R
	@for model in $(2D_SHELL_MODELS); do \
		DIR_NAME=$$(echo $$model | tr '/-' '_'); \
		IN_DIR=$(SIMULATION)/2d_shell/results/$$DIR_NAME; \
		OUT_FIG_DIR=$(FIGS)/simulation/$$DIR_NAME; \
		VIS_TYPE=2d_shell; \
		if [ -d "$$IN_DIR" ]; then \
			Rscript $(R)/simulation/main.R $$IN_DIR $$OUT_FIG_DIR; \
		fi \
	done

misc: $(R)/misc/main.R
	@Rscript $(R)/misc/main.R $(DEPTH_AVERAGE_PATHS) $(FIGS)/simulation $(FIGS)/misc

manuscript: $(R)/manuscript/main.R
	@Rscript $(R)/manuscript/main.R $(R)/manuscript/data/rheological-parameters.csv $(MARKDOWN)/manuscript/rheological-parameters.md

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
	@echo "  visualize       Visualize all results"
	@echo "  adiabat         Visualize adiabatic profiles"
	@echo "  2d-box-models   Visualize ASPECT box model statistics"
	@echo "  2d-shell-models Visualize ASPECT shell model statistics"
	@echo "  misc            Visualize miscelaneous plots"
	@echo "  clean           Cleanup unnecessary files and directories (safe)"
	@echo "  deep-clean      Deep clean figures and results (use with caution!)"
	@echo "  help            Show this help message"
