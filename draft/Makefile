# Top-level dirs
PROJECT_ROOT := $(realpath $(dir $(lastword $(MAKEFILE_LIST)))/..)
R            := $(PROJECT_ROOT)/R
PYTHON       := $(PROJECT_ROOT)/python
DRAFT        := $(PROJECT_ROOT)/draft

# Templates
TEMPLATE      := $(DRAFT)/agu.tex
TEMPLATE_NAME := $(notdir $(TEMPLATE))
BIB           := $(DRAFT)/main.bib

ifeq ($(TEMPLATE_NAME),agu.tex)
    CFG_YAML      := $(DRAFT)/agu.yaml
    CFG_SI_YAML   := $(DRAFT)/agu-si.yaml
    CFG_SI_TPL    := $(DRAFT)/agu-si.tex
else ifeq ($(TEMPLATE_NAME),eisvogel.tex)
    CFG_YAML      := $(DRAFT)/eisvogel.yaml
    CFG_SI_YAML   := $(DRAFT)/eisvogel.yaml
    CFG_SI_TPL    := $(TEMPLATE)
endif

# Configure pandoc
BASE_PANDOC_FLAGS := --from markdown --filter pandoc-crossref $(PANDOC_FLAGS)

$(DRAFT)/manuscript.%: TPL                      := $(TEMPLATE)
$(DRAFT)/manuscript.%: META                     := $(CFG_YAML)
$(DRAFT)/manuscript.%: BIB_FLAG                 := --bibliography=$(BIB)
$(DRAFT)/manuscript.%: CFG_EXTRA                := --citeproc --csl=$(DRAFT)/agu.csl --pdf-engine=pdflatex

$(DRAFT)/supplementary-information.%: TPL       := $(CFG_SI_TPL)
$(DRAFT)/supplementary-information.%: META      := $(CFG_SI_YAML)
$(DRAFT)/supplementary-information.%: BIB_FLAG  :=

$(DRAFT)/diff-manuscript.%: META                := $(CFG_YAML)
$(DRAFT)/diff-supplementary-information.%: META := $(CFG_SI_YAML)

$(DRAFT)/response.%: TPL                        :=
$(DRAFT)/response.%: META                       := $(DRAFT)/response.yaml
$(DRAFT)/response.%: BIB_FLAG                   :=

$(DRAFT)/erratum.%: CFG_EXTRA                   := --citeproc --csl=$(DRAFT)/agu.csl --pdf-engine=pdflatex
$(DRAFT)/erratum.%: BIB_FLAG                    := --bibliography=$(BIB)

TPL_FLAG  = $(if $(TPL),--template=$(TPL))
META_FLAG = $(if $(META),--metadata-file=$(META))
PANDOC_FLAGS = --standalone --from=markdown --filter=pandoc-crossref $(CFG_EXTRA) $(TPL_FLAG) $(BIB_FLAG) $(META_FLAG)

# Rendered files
DOCS := manuscript supplementary-information
PDFS := $(addprefix $(DRAFT)/, $(addsuffix .pdf, $(DOCS)))

# Erratum
ERRATUM     := erratum
ERRATUM_PDF := $(addprefix $(DRAFT)/, $(addsuffix .pdf, $(ERRATUM)))

# Configure pdflatex
PDFLATEX = pdflatex -interaction=nonstopmode -shell-escape -output-directory $(DRAFT)

# Configure latexdiff
BACKUP_DATE ?= 07-jul-2026
DIFF_FLAGS = --flatten --math-markup=0 --type=CFONT --graphics-markup=none --exclude-textcmd="ref,eqref,label" --add-to-config PICTUREENV=longtable --config="ARRENV=thebibliography|APACrefauthors"

# Cleanup targets
EXTS        := aux bbl blg log toc out tex pdf
LETTER_EXTS := $(foreach letter,$(ERRATUM),$(addprefix $(DRAFT)/$(letter).,$(EXTS)))
CLEAN       := $(foreach doc,$(DOCS) $(addprefix diff-,$(DOCS)),$(addprefix $(DRAFT)/$(doc).,$(EXTS))) $(addprefix $(DRAFT)/response.,$(EXTS)) $(LETTER_EXTS)
DEEP_CLEAN  := $(PROJECT_ROOT)/figs

# Prevent make from deleting .tex files on failure
.PRECIOUS: $(DRAFT)/%.tex

# Targets
.PHONY: all diff clean deep-clean help

all: $(PDFS)
erratum: $(ERRATUM_PDF)
diff: all $(DRAFT)/diff-manuscript.pdf $(DRAFT)/diff-supplementary-information.pdf $(DRAFT)/response.pdf

$(DRAFT)/%.pdf: $(DRAFT)/%.tex
	@echo " -> $@"
	@pandoc $(subst .tex,.md,$<) -o $@ $(PANDOC_FLAGS)

$(DRAFT)/bak/%.tex: ;

$(DRAFT)/%.tex: $(DRAFT)/%.md
	@echo " -> $@"
	@$(if $(META),sed -i '' -E "s/^date: \".*\"/date: \"$(shell date '+%d %B %Y')\"/" $(META))
	@pandoc $< -o $@ $(PANDOC_FLAGS)

$(DRAFT)/diff-%.pdf: $(DRAFT)/bak/%-$(BACKUP_DATE).tex $(DRAFT)/%.tex
	@echo " -> $@ [bak: $(BACKUP_DATE)]"
	@if [ ! -f "$<" ]; then \
	  echo "Error: Backup file $< not found!"; \
	  exit 1; \
	fi
	@latexdiff $(DIFF_FLAGS) $< $(word 2,$^) > $(DRAFT)/diff-$*.tex
	@$(PDFLATEX) $(DRAFT)/diff-$*.tex > /dev/null 2>&1
	@$(PDFLATEX) $(DRAFT)/diff-$*.tex > /dev/null 2>&1

clean:
	@for item in $(CLEAN); do \
	  safe_rm() { \
	    if [ -e "$$1" ]; then \
	      ABS_PATH=$$(realpath "$$1"); \
	      if [[ "$$ABS_PATH" == $(PROJECT_ROOT)* ]]; then \
	        echo " xx $$ABS_PATH"; \
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
	        echo " xx $$ABS_PATH"; \
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
	@echo "    all                                Render manuscript and SI"
	@echo "    erratum                            Render the draft erratum"
	@echo "    diff-manuscript.pdf                Generate manuscript diff"
	@echo "    diff-supplementary-information.pdf Generate SI diff"
	@echo "    clean                              Cleanup unnecessary files"
	@echo "    --------------------------------------------------"
	@echo "    Note: Use BACKUP_DATE=DD-MMM-YYYY to change diff source."
