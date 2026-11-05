# paper/Makefile: the paper's reproduction (README.md). SET = a key of sets.tsv, or main (D1-D8), supp, all.
#   make data SET=...     rebuild the sets' POD5s from their public sources and check them (data/extract.sh)
#   make squash SET=...   run SQUASH on them as the paper did (squash.sh)
#   make check SET=...    data + squash, then compare SQUASH's outputs with the paper's (check.sh, expected/)
#   make threads SET=...  run SQUASH again at THREADS2 threads (config.sh) and compare with the paper's outputs
#   make curate SET=...   remake the sets' read lists from their sources and check them (curation/remake.sh)
#   make tools            install the other tools at the paper's versions into TOOLS_DIR (tools/install.sh)
#   make basecall SET=... the sets' sup reads with move tables and HAC reads (basecall/run.sh; a GPU)
#   make timing SET=...   the paper's timed basecalling routes (basecall/run.sh with TIMED_ONLY=1)
#   make assemble SET=... the other assemblers on the sets' reads, timed (assemble/run.sh; TOOLS=... for some only)
#   make eval SET=...     the evaluation of SQUASH's and the other assemblers' runs (eval/eval.sh; SQUASH's run with
#                         KEEP_WORK=1: KEEP_WORK=1 make check SET=...); make eval-refs SET=...: only its references
#   make tables           the paper's tables and figure from the runs and evaluations under OUT_DIR (tables/)
#   make test             the tests of the POD5 tools, the selection, the tables, the CPU-timing subset and the
#                         evaluation (PYTHON: requirements.txt's packages)
# Settings: config.sh.
SHELL := /bin/bash
SET ?= cjk60_both60
SETS_TSV := sets.tsv
KEYS := $(shell awk -F'\t' -v s='$(SET)' 'NR > 1 && (s == "all" || $$1 == s || (s == "main" && $$2 != "supp") || (s == "supp" && $$2 == "supp")) {print $$1}' $(SETS_TSV))
CONF := source config.sh

.PHONY: data squash check threads curate tools basecall timing assemble eval eval-refs tables test
ifeq ($(strip $(KEYS)),)
$(error SET=$(SET): no such set in $(SETS_TSV) (a key, main, supp or all))
endif

data:
	@for k in $(KEYS); do bash data/extract.sh $$k || exit 1; done

squash:
	@for k in $(KEYS); do bash squash.sh $$k || exit 1; done

check: data squash
	@for k in $(KEYS); do bash check.sh $$k || exit 1; done

threads: data
	@$(CONF); for k in $(KEYS); do RUN_THREADS=$$THREADS2 RUN_OUT=$$OUT_DIR/$$k.t$$THREADS2 bash squash.sh $$k && \
	  bash check.sh $$k $$OUT_DIR/$$k.t$$THREADS2 || exit 1; done

curate:
	@for k in $(KEYS); do bash curation/remake.sh $$k || exit 1; done

tools:
	@bash tools/install.sh

basecall:
	@for k in $(KEYS); do bash basecall/run.sh $$k || exit 1; done

timing:
	@for k in $(KEYS); do TIMED_ONLY=1 bash basecall/run.sh $$k || exit 1; done

assemble:
	@for k in $(KEYS); do bash assemble/run.sh $$k $(TOOLS) || exit 1; done

eval:
	@for k in $(KEYS); do bash eval/eval.sh $$k || exit 1; done

eval-refs:
	@bash eval/refs.sh $(KEYS)

tables:
	@bash tables/make.sh

test:
	@$(CONF); $$PYTHON data/test_pod5_tools.py && $$PYTHON curation/test_select.py && $$PYTHON tables/test_tables.py && \
	  $$PYTHON basecall/test_cpu_subset.py && $$PYTHON eval/test_eval.py
