# SQUASH only (no evaluation tools). Build: docker build -t squash:1.0.0 .
# Apptainer: apptainer build squash.sif docker-daemon://squash:1.0.0
FROM rust:1.95.0-slim-bookworm AS build
WORKDIR /src
COPY . .
RUN cargo build --release --locked

FROM debian:bookworm-slim
# pigz: SQUASH compresses its large intermediate tables with it when present (the same outputs without it)
RUN apt-get update && apt-get install -y --no-install-recommends pigz && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/target/release/squash /usr/local/bin/squash
ENTRYPOINT ["squash"]
