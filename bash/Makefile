# OS Detection
UNAME_S := $(shell uname)
OS_ID := $(shell [ -f /etc/os-release ] && . /etc/os-release && echo $$ID || echo unknown)

# Top-level dirs
PROJECT_ROOT := $(realpath $(dir $(lastword $(MAKEFILE_LIST)))/..)
DEALII := $(PROJECT_ROOT)/dealii
ASPECT := $(PROJECT_ROOT)/aspect
INSTALL := $(PROJECT_ROOT)/bash/install

# Platform scripts
DEALII_MACOS := $(INSTALL)/dealii-macos.sh
ASPECT_MACOS := $(INSTALL)/aspect-macos.sh
DEALII_BARKLA2 := $(INSTALL)/dealii-barkla2.sh
ASPECT_BARKLA2 := $(INSTALL)/aspect-barkla2.sh

# ASPECT executable
ASPECT_EXE := $(ASPECT)/aspect-release

# Install config
DEALII_VERSION := v9.6.1

ifeq ($(UNAME_S), Darwin)
  DEALII_LIB := $(DEALII)/deal.II-$(DEALII_VERSION)/lib/libdeal_II.dylib
  INSTALL_DEALII := $(DEALII_MACOS)
  INSTALL_ASPECT := $(ASPECT_MACOS)
  TRILINOS_VERSION := 16
else ifeq ($(OS_ID), rocky)
  DEALII_LIB := $(DEALII)/deal.II-$(DEALII_VERSION)/lib/libdeal_II.so
  INSTALL_DEALII := $(DEALII_BARKLA2)
  INSTALL_ASPECT := $(ASPECT_BARKLA2)
  TRILINOS_VERSION := 14
  GCC := gcc/14.2.0
  OPENMPI := openmpi/5.0.8-gcc14.2.0
  OPENBLAS := openblas/0.3.29/gcc-14.2.0
  CMAKE := cmake/3.30.5-gcc14.2.0
else
  $(error Unsupported OS: Only macOS and Rocky Linux are supported!)
endif

# Logging
DATE := $(shell date +"%d-%m-%Y")
LOG_FILE := $(INSTALL)/log-$(DATE).log
LOGGER := 2>&1 | tee -a $(LOG_FILE)

# Cleanup targets
CLEAN := $(INSTALL)/log-*.log $(DEALII)/tmp
DEEP_CLEAN :=

# Targets
.PHONY: all install uninstall install-aspect install-dealii uninstall-dealii uninstall-aspect clean help

all: install

install: install-dealii install-aspect

install-aspect:$(LOG_FILE) $(ASPECT_EXE)

install-dealii:$(LOG_FILE) $(DEALII_LIB)

uninstall: uninstall-aspect uninstall-dealii

uninstall-aspect:
	@echo "==> Uninstalling ASPECT ..."
	@for item in $(ASPECT); do \
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

uninstall-dealii:
	@echo "==> Uninstalling deal.II ..."
	@for item in $(DEALII); do \
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

$(ASPECT_EXE): $(LOG_FILE)
	@echo "==> Installing ASPECT ..." $(LOGGER)
	@if [ "$(OS_ID)" = "rocky" ]; then \
		if echo "$(shell hostname)" | grep -q "login"; then \
			echo "This script cannot be run on the BARKLA login node!" $(LOGGER); \
			echo "Start an interactive session with 'srun --ntasks=<N> --pty /bin/bash'" $(LOGGER); \
			exit 1; \
		fi; \
		bash -i -c "source $(INSTALL_ASPECT) $(DEALII_VERSION) $(GCC) $(OPENMPI) $(OPENBLAS) $(CMAKE)" $(LOGGER) || exit 1; \
	else \
		$(INSTALL_ASPECT) $(DEALII_VERSION) $(LOGGER); \
	fi

$(DEALII_LIB): $(LOG_FILE)
	@echo "==> Installing deal.II ..." $(LOGGER)
	@if [ "$(OS_ID)" = "rocky" ]; then \
		if echo "$(shell hostname)" | grep -q "login"; then \
			echo "This script cannot be run on the BARKLA login node!" $(LOGGER); \
			echo "Start an interactive session with 'srun --ntasks=<N> --pty /bin/bash'" $(LOGGER); \
			exit 1; \
		fi; \
		bash -i -c "source $(INSTALL_DEALII) $(DEALII_VERSION) $(TRILINOS_VERSION) $(GCC) $(OPENMPI) $(OPENBLAS) $(CMAKE)" $(LOGGER) || exit 1; \
	else \
		$(INSTALL_DEALII) $(DEALII_VERSION) $(LOGGER); \
	fi

$(LOG_FILE):
	@mkdir -p log
	@touch $(LOG_FILE)

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
	@echo "  install           Build deal.II and ASPECT"
	@echo "  install-aspect    Build ASPECT"
	@echo "  install-dealii    Build deal.II"
	@echo "  uninstall         Uninstall deal.II and ASPECT"
	@echo "  uninstall-aspect  UninstallASPECT"
	@echo "  uninstall-dealii  Uninstall deal.II"
	@echo "  clean             Cleanup unnecessary files and directories (safe)"
	@echo "  deep-clean        Deep clean figures and results (use with caution!)"
	@echo "  help              Show this help message"
