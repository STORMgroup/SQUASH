# SQUASH: `make` builds target/release/squash; `make test` runs the tests; `make install PREFIX=...` copies the binary.
PREFIX ?= /usr/local
CARGO ?= cargo
GIT_COMMIT := $(shell git rev-parse --short HEAD 2>/dev/null)

.PHONY: all build test install clean

all: build

build:
	SQUASH_GIT_COMMIT=$(GIT_COMMIT) $(CARGO) build --release --locked
	@echo "built target/release/squash"

test:
	$(CARGO) test --release --locked

install: build
	install -d $(PREFIX)/bin
	install -m 755 target/release/squash $(PREFIX)/bin/squash

clean:
	$(CARGO) clean
