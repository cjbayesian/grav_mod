#!/bin/bash
## Usage:
##   ./push_to_server.sh <remote folder>
rsync -rva --delete --exclude-from="server_excludes" ./ darwin:~/CAISN/$1
