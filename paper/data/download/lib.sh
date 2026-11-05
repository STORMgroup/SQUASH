# paper/data/download/lib.sh: sourced by the download scripts (bash, curl, md5sum, tar, unzip).
#
# Every source goes to SOURCES_DIR/SOURCE/ (default paper/data/sources/SOURCE; paper/config.sh may set SOURCES_DIR). A
# script is resumable: each file is fetched with `curl -C -` into FILE.part, its size and MD5 checked, then renamed; a
# finished source has SOURCES_DIR/SOURCE/.complete and is skipped. An existing local copy is used instead of a download
# when the variable SRC_<SOURCE> (upper case, e.g. SRC_NVB_EC) names a folder holding the source's POD5 files as the
# download leaves them: SOURCES_DIR/SOURCE becomes a link to it (no checksum: it is your copy).
# Other settings: CURL (download command), PYTHON (a python with pod5 0.3.39: FAST5 conversion, the human cut),
# STREAM=1 (tar sources: unpack while downloading, MD5 on the stream: half the disk, no resume).

set -euo pipefail
DL_HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PAPER_DATA=$(cd "$DL_HERE/.." && pwd)
[[ -f $PAPER_DATA/../config.sh ]] && source "$PAPER_DATA/../config.sh"
SOURCES_DIR=${SOURCES_DIR:-$PAPER_DATA/sources}
CURL=${CURL:-curl -sSfL --retry 5 --retry-delay 30 --connect-timeout 60}
PYTHON=${PYTHON:-python3}
MANIFESTS=$DL_HERE/manifests

log() { echo "[$(date '+%F %T')] $*" >&2; }
die() { echo "download: $*" >&2; exit 1; }

# source_dir SOURCE: SOURCES_DIR/SOURCE, or the link to SRC_<SOURCE>; returns 1 (prints the dir) when already complete
source_dir() {
  local s=$1 var d
  var=SRC_$(echo "$s" | tr '[:lower:]-' '[:upper:]_')
  d=$SOURCES_DIR/$s
  mkdir -p "$SOURCES_DIR"
  if [[ -n ${!var:-} ]]; then
    [[ -d ${!var} ]] || die "$var=${!var} is not a folder"
    if [[ ! -L $d ]]; then
      [[ ! -e $d ]] || die "$d exists and is not a link: remove it to use $var"
      ln -s "$(cd "${!var}" && pwd)" "$d"
    fi
    log "$s: using the local copy ${!var} ($var)"
    echo "$d"; return 1
  fi
  mkdir -p "$d"
  echo "$d"
  [[ ! -e $d/.complete ]]
}

# fetch URL OUT BYTES [MD5]: one file, resumable; size (and MD5 if given) checked; OUT.part renamed to OUT when right
fetch() {
  local url=$1 out=$2 bytes=$3 md5=${4:-} try have
  [[ -s $out && $(stat -c %s "$out") == "$bytes" ]] && return 0
  mkdir -p "$(dirname "$out")"
  for try in 1 2 3 4 5 6 7 8; do
    have=0; [[ -f $out.part ]] && have=$(stat -c %s "$out.part")
    (( have > bytes )) && { rm -f -- "$out.part"; have=0; }
    (( have == bytes )) && break
    $CURL -C - -o "$out.part" "$url" || log "retry $try: $url"
    have=0; [[ -f $out.part ]] && have=$(stat -c %s "$out.part")
    (( have == bytes )) && break
    sleep $((try * 30))
  done
  have=0; [[ -f $out.part ]] && have=$(stat -c %s "$out.part")
  (( have == bytes )) || die "$url: $have bytes, expected $bytes"
  if [[ -n $md5 ]]; then
    local got; got=$(md5sum "$out.part" | cut -d' ' -f1)
    [[ $got == "$md5" ]] || { mv "$out.part" "$out.bad"; die "$url: MD5 $got, expected $md5 (kept as $out.bad)"; }
  fi
  mv "$out.part" "$out"
}

# tar_source DIR URL BYTES MD5 [TAR_OPTS]: a tar(.gz) unpacked into DIR (downloaded first, or streamed with STREAM=1)
tar_source() {
  local d=$1 url=$2 bytes=$3 md5=$4 opts=${5:-} name t
  name=$(basename "${url%%\?*}")
  mkdir -p "$d/.unpack"
  if [[ ${STREAM:-0} == 1 ]]; then
    log "streaming $name ($bytes bytes) into $d"
    rm -f -- "$d/.stream.$name.md5" "$d/.stream.$name.bytes"
    $CURL "$url" | tee >(md5sum | cut -d' ' -f1 > "$d/.stream.$name.md5") >(wc -c > "$d/.stream.$name.bytes") \
      | { tar -x $opts -C "$d/.unpack"; cat > /dev/null; }
    for t in $(seq 1 120); do [[ -s $d/.stream.$name.md5 && -s $d/.stream.$name.bytes ]] && break; sleep 2; done
    [[ $(tr -dc 0-9 < "$d/.stream.$name.bytes") == "$bytes" ]] || die "$name: size differs on the stream"
    [[ $(cat "$d/.stream.$name.md5") == "$md5" ]] || die "$name: MD5 differs on the stream"
    rm -f -- "$d/.stream.$name.md5" "$d/.stream.$name.bytes"
  else
    fetch "$url" "$d/.unpack/$name" "$bytes" "$md5"
    log "unpacking $name"
    tar -x $opts -f "$d/.unpack/$name" -C "$d/.unpack"
    rm -f -- "$d/.unpack/$name"
  fi
}

# mark_complete DIR SOURCE: count the POD5 files and write DIR/.complete
mark_complete() {
  local d=$1 s=$2 n
  n=$(find -L "$d" -name '*.pod5' -type f | wc -l)
  (( n > 0 )) || die "$s: no POD5 file in $d"
  printf '%s\t%s POD5 files\t%s\n' "$s" "$n" "$(date '+%F %T')" > "$d/.complete"
  log "$s: complete ($n POD5 files in $d)"
}

# nvb_source SOURCE FIGSHARE_FILE_ID BYTES MD5: a NanoVarBench 202309 sample (one figshare tar of POD5s); its POD5 files
# flattened into SOURCES_DIR/SOURCE/ (as the paper's download did)
nvb_source() {
  local s=$1 id=$2 bytes=$3 md5=$4 d
  d=$(source_dir "$s") || { log "$s: done before ($d)"; return 0; }
  tar_source "$d" "https://ndownloader.figshare.com/files/$id" "$bytes" "$md5"
  find "$d/.unpack" -name '*.pod5' -type f -exec mv -t "$d/" {} +
  rm -rf -- "$d/.unpack"
  mark_complete "$d" "$s"
}
