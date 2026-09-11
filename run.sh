#!/usr/bin/env bash
# run.sh -- launch Epitope-Evaluator 2 and report the port to tunnel to.
#
#   ./run.sh            pick the first free port from the list below
#   ./run.sh 8123       use a specific port
#
# On a shared login node the usual ports are often taken by other users, so the
# port actually used is printed before the app starts.
set -euo pipefail
cd "$(dirname "$0")"

port_free() { ! ss -ltn 2>/dev/null | awk '{print $4}' | grep -q ":$1\$"; }

if [ $# -ge 1 ]; then
  PORT=$1
  port_free "$PORT" || { echo "Port $PORT is already in use."; exit 1; }
else
  PORT=""
  for p in 8080 8081 8082 8090 8123 8456 9080; do
    if port_free "$p"; then PORT=$p; break; fi
  done
  [ -n "$PORT" ] || { echo "No free port found; pass one explicitly: ./run.sh 9999"; exit 1; }
fi

echo
echo "  Epitope-Evaluator 2 starting on $(hostname -s), port $PORT"
echo
echo "  From your laptop:"
echo "      ssh -N -L ${PORT}:127.0.0.1:${PORT} ${USER}@$(hostname -f)"
echo "  then open:"
echo "      http://localhost:${PORT}"
echo
echo "  Ctrl-C here to stop."
echo

exec Rscript -e "shiny::runApp('.', port = ${PORT}, host = '127.0.0.1', launch.browser = FALSE)"
