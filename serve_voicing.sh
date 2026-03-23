#!/bin/bash
# Launch Harp Voicing Helper in browser
cd "$(dirname "$0")"
PORT=8000
echo "Starting Harp Voicing Helper on http://localhost:$PORT/harp_voicing.html"
python3 -m http.server $PORT &
PID=$!
sleep 1
xdg-open "http://localhost:$PORT/harp_voicing.html" 2>/dev/null || \
  open "http://localhost:$PORT/harp_voicing.html" 2>/dev/null || \
  echo "Open http://localhost:$PORT/harp_voicing.html in your browser"
echo "Press Ctrl+C to stop"
trap "kill $PID 2>/dev/null" EXIT
wait $PID
